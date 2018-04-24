/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PCCTMC3Encoder_h
#define PCCTMC3Encoder_h

#include <assert.h>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include "ringbuf.h"

#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCPointSetProcessing.h"
#include "PCCTMC3Common.h"

#include "ArithmeticCodec.h"

namespace pcc {
struct PCCAttributeEncodeParamaters {
  size_t numberOfNearestNeighborsInPrediction;
  size_t levelOfDetailCount;
  size_t searchRange;
  std::vector<size_t> dist2;
  std::vector<size_t> quantizationSteps;
  std::vector<size_t> quantizationDeadZoneSizes;
};

struct PCCTMC3Encoder3Parameters {
  double positionQuantizationScale;
  bool mergeDuplicatedPoints;
  std::map<std::string, PCCAttributeEncodeParamaters> attributeEncodeParameters;
};

class PCCTMC3Encoder3 {
 public:
  PCCTMC3Encoder3() { init(); }
  PCCTMC3Encoder3(const PCCTMC3Encoder3 &) = default;
  PCCTMC3Encoder3 &operator=(const PCCTMC3Encoder3 &rhs) = default;
  ~PCCTMC3Encoder3() = default;
  void init() {
    minPositions = 0.0;
    boundingBox.min = uint32_t(0);
    boundingBox.max = uint32_t(0);
    pointCloud.clear();
  }

  size_t estimateBitstreamSize(const PCCPointSet3 &pointCloud,
                               const PCCTMC3Encoder3Parameters &params) {
    size_t bitstreamSize = pointCloud.getPointCount() * 3 * 4 + 1024;  // positions
    if (pointCloud.hasColors()) {                                      // colors
      bitstreamSize += pointCloud.getPointCount() * 3 * 2 + 1024;
    }
    if (pointCloud.hasReflectances()) {  // reflectance
      bitstreamSize += pointCloud.getPointCount() * 2 + 1024;
    }
    return bitstreamSize;
  }

  int compressWithLosslessGeometry(const PCCPointSet3 &inputPointCloud,
                                   const PCCTMC3Encoder3Parameters &params, PCCBitstream &bitstream,
                                   PCCPointSet3 *reconstructedCloud = nullptr) {
    init();
    pointCloud = inputPointCloud;
    const bool compressColors =
        inputPointCloud.hasColors() &&
        params.attributeEncodeParameters.find("color") != params.attributeEncodeParameters.end();
    if (!compressColors) {
      pointCloud.removeColors();
    }

    const bool compressReflectances = inputPointCloud.hasReflectances() &&
                                      params.attributeEncodeParameters.find("reflectance") !=
                                          params.attributeEncodeParameters.end();
    if (!compressReflectances) {
      pointCloud.removeReflectances();
    }

    PCCWriteToBuffer<uint32_t>(PCCTMC3MagicNumber, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint32_t>(PCCTMC3FormatVersion, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasColors()), bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasReflectances()), bitstream.buffer,
                              bitstream.size);

    if (pointCloud.hasColors()) {
      const auto &colorParams = params.attributeEncodeParameters.find("color")->second;
      uint64_t colorsSize = bitstream.size;
      if (int ret = encodeAttributeHeader(colorParams, "color", bitstream)) {
        return ret;
      }
      buildPredictors(colorParams);
      computeColorPredictionWeights(colorParams);
      if (int ret = encodeColors(colorParams, bitstream)) {
        return ret;
      }
      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B ("
                << (8.0 * colorsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;
    }

    if (pointCloud.hasReflectances()) {
      const auto &reflectanceParams = params.attributeEncodeParameters.find("reflectance")->second;
      uint64_t reflectancesSize = bitstream.size;
      if (int ret = encodeAttributeHeader(reflectanceParams, "reflectance", bitstream)) {
        return ret;
      }
      buildPredictors(reflectanceParams);
      computeReflectancePredictionWeights(reflectanceParams);
      if (int ret = encodeReflectances(reflectanceParams, bitstream)) {
        return ret;
      }
      reflectancesSize = bitstream.size - reflectancesSize;
      std::cout << "reflectances bitstream size " << reflectancesSize << " B ("
                << (8.0 * reflectancesSize) / inputPointCloud.getPointCount() << " bpp)"
                << std::endl;
    }

    if (reconstructedCloud) {
      (*reconstructedCloud) = pointCloud;
    }
    return 0;
  }
  int compress(const PCCPointSet3 &inputPointCloud, const PCCTMC3Encoder3Parameters &params,
               PCCBitstream &bitstream, PCCPointSet3 *reconstructedCloud = nullptr) {
    init();
    PCCWriteToBuffer<uint32_t>(PCCTMC3MagicNumber, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint32_t>(PCCTMC3FormatVersion, bitstream.buffer, bitstream.size);

    if (int ret = quantization(inputPointCloud, params)) {
      return ret;
    }
    uint64_t positionsSize = bitstream.size;
    if (int ret = encodePositionsHeader(params, bitstream)) {
      return ret;
    }
    if (int ret = encodePositions(params, bitstream)) {
      return ret;
    }
    positionsSize = bitstream.size - positionsSize;
    std::cout << "positions bitstream size " << positionsSize << " B ("
              << (8.0 * positionsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;

    if (pointCloud.hasColors()) {
      const auto &colorParams = params.attributeEncodeParameters.find("color")->second;
      uint64_t colorsSize = bitstream.size;
      if (int ret = encodeAttributeHeader(colorParams, "color", bitstream)) {
        return ret;
      }
      buildPredictors(colorParams);
      computeColorPredictionWeights(colorParams);
      if (int ret = encodeColors(colorParams, bitstream)) {
        return ret;
      }
      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B ("
                << (8.0 * colorsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;
    }

    if (pointCloud.hasReflectances()) {
      const auto &reflectanceParams = params.attributeEncodeParameters.find("reflectance")->second;
      uint64_t reflectancesSize = bitstream.size;
      if (int ret = encodeAttributeHeader(reflectanceParams, "reflectance", bitstream)) {
        return ret;
      }
      buildPredictors(reflectanceParams);
      computeReflectancePredictionWeights(reflectanceParams);
      if (int ret = encodeReflectances(reflectanceParams, bitstream)) {
        return ret;
      }
      reflectancesSize = bitstream.size - reflectancesSize;
      std::cout << "reflectances bitstream size " << reflectancesSize << " B ("
                << (8.0 * reflectancesSize) / inputPointCloud.getPointCount() << " bpp)"
                << std::endl;
    }

    reconstructedPointCloud(params, reconstructedCloud);
    return 0;
  }

 private:
  static void encodeAbsUInt32(const uint32_t value, const uint32_t bitCount,
                              o3dgc::Arithmetic_Codec &arithmeticEncoder,
                              o3dgc::Static_Bit_Model &binaryModel0) {
    uint32_t valueToEncode = PCCToLittleEndian<uint32_t>(value);
    for (uint32_t i = 0; i < bitCount; ++i) {
      arithmeticEncoder.encode(valueToEncode & 1, binaryModel0);
      valueToEncode >>= 1;
    }
  }
  static void encodeDiff0UInt32(const uint32_t value, const uint32_t maxAttributeValueDiff0,
                                const uint32_t adaptiveDataModelAlphabetSizeDiff0,
                                o3dgc::Arithmetic_Codec &arithmeticEncoder,
                                o3dgc::Adaptive_Data_Model &multiSymbolModelDiff0,
                                o3dgc::Adaptive_Bit_Model &binaryModelDiff0,
                                o3dgc::Static_Bit_Model &binaryModel0) {
    if (!maxAttributeValueDiff0) {
      return;
    } else if (value < adaptiveDataModelAlphabetSizeDiff0) {
      arithmeticEncoder.encode(value, multiSymbolModelDiff0);
    } else {
      arithmeticEncoder.encode(adaptiveDataModelAlphabetSizeDiff0, multiSymbolModelDiff0);
      arithmeticEncoder.ExpGolombEncode(value - adaptiveDataModelAlphabetSizeDiff0, 0, binaryModel0,
                                        binaryModelDiff0);
    }
  }
  static void encodeDiff1UInt32(const uint32_t value, const uint32_t modelIndex,
                                const uint32_t maxAttributeValueDiff1,
                                const uint32_t adaptiveDataModelAlphabetSizeDiff1,
                                o3dgc::Arithmetic_Codec &arithmeticEncoder,
                                std::vector<o3dgc::Adaptive_Data_Model> &multiSymbolModelDiff1,
                                std::vector<o3dgc::Adaptive_Bit_Model> &binaryModelDiff1,
                                o3dgc::Static_Bit_Model &binaryModel0) {
    if (!maxAttributeValueDiff1) {
      return;
    } else if (value < adaptiveDataModelAlphabetSizeDiff1) {
      arithmeticEncoder.encode(value, multiSymbolModelDiff1[modelIndex]);
    } else {
      arithmeticEncoder.encode(adaptiveDataModelAlphabetSizeDiff1,
                               multiSymbolModelDiff1[modelIndex]);
      arithmeticEncoder.ExpGolombEncode(value - adaptiveDataModelAlphabetSizeDiff1, 0, binaryModel0,
                                        binaryModelDiff1[modelIndex]);
    }
  }
  int encodeReflectances(const PCCAttributeEncodeParamaters &reflectanceParams,
                         PCCBitstream &bitstream) {
    const size_t pointCount = pointCloud.getPointCount();
    assert(pointCount == predictors.size());
    uint32_t maxAttributeValueDiff0 = 0;
    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      const auto &predictor = predictors[predictorIndex];
      if (predictor.maxNeighborCount) {
        const int64_t quantAttValue = pointCloud.getReflectance(predictor.index);
        const int64_t quantPredAttValue = predictor.predictReflectance(pointCloud);
        const uint32_t diffAttValue =
            uint32_t(o3dgc::IntToUInt(long(quantAttValue - quantPredAttValue)));
        if (maxAttributeValueDiff0 < diffAttValue) {
          maxAttributeValueDiff0 = diffAttValue;
        }
      }
    }

    uint32_t adaptiveDataModelAlphabetSizeDiff0 = 0;
    o3dgc::Adaptive_Bit_Model binaryModelDiff0;
    o3dgc::Adaptive_Data_Model multiSymbolModelDiff0;
    if (maxAttributeValueDiff0) {
      adaptiveDataModelAlphabetSizeDiff0 = 1 + maxAttributeValueDiff0;
      if (adaptiveDataModelAlphabetSizeDiff0 > PCCTMC3AdaptiveDataModelAlphabetMaxSize) {
        adaptiveDataModelAlphabetSizeDiff0 = PCCTMC3AdaptiveDataModelAlphabetMaxSize;
      }
      multiSymbolModelDiff0.set_alphabet(adaptiveDataModelAlphabetSizeDiff0 + 1);
      binaryModelDiff0.reset();
    }

    uint64_t startSize = bitstream.size;
    bitstream.size += 4;  // placehoder for bitstream size
    PCCWriteToBuffer<uint32_t>(maxAttributeValueDiff0, bitstream.buffer, bitstream.size);

    o3dgc::Arithmetic_Codec arithmeticEncoder;
    arithmeticEncoder.set_buffer(static_cast<uint32_t>(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticEncoder.start_encoder();
    o3dgc::Static_Bit_Model binaryModel0;
    o3dgc::Adaptive_Data_Model neighborCountModel;
    neighborCountModel.set_alphabet(
        uint32_t(reflectanceParams.numberOfNearestNeighborsInPrediction));
    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      const auto &predictor = predictors[predictorIndex];
      if (!predictor.neighborCount) {  // no prediction
        encodeAbsUInt32(pointCloud.getReflectance(predictor.index), 16, arithmeticEncoder,
                        binaryModel0);
      } else {
        const size_t lodIndex = predictor.levelOfDetailIndex;
        const int64_t qs = reflectanceParams.quantizationSteps[lodIndex];
        const int64_t dz = reflectanceParams.quantizationDeadZoneSizes[lodIndex];
        arithmeticEncoder.encode(uint32_t(predictor.neighborCount - 1), neighborCountModel);

        const int64_t quantAttValue = pointCloud.getReflectance(predictor.index);
        const int64_t quantPredAttValue = predictor.predictReflectance(pointCloud);
        const int64_t delta = PCCQuantization(quantAttValue - quantPredAttValue, qs, dz);
        const uint32_t attValue0 = uint32_t(o3dgc::IntToUInt(long(delta)));
        encodeDiff0UInt32(attValue0, maxAttributeValueDiff0, adaptiveDataModelAlphabetSizeDiff0,
                          arithmeticEncoder, multiSymbolModelDiff0, binaryModelDiff0, binaryModel0);

        const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs, dz);
        const int64_t reconstructedQuantAttValue = quantPredAttValue + reconstructedDelta;
        const uint16_t reconstructedReflectance = uint16_t(PCCClip(
            reconstructedQuantAttValue, int64_t(0), int64_t(std::numeric_limits<uint16_t>::max())));
        pointCloud.setReflectance(predictor.index, reconstructedReflectance);
      }
    }
    uint32_t compressedBitstreamSize = arithmeticEncoder.stop_encoder();
    bitstream.size += compressedBitstreamSize;
    PCCWriteToBuffer<uint32_t>(compressedBitstreamSize, bitstream.buffer, startSize);
    return 0;
  }

  int encodeAttributeHeader(const PCCAttributeEncodeParamaters &attributeParams,
                            const std::string &attributeName, PCCBitstream &bitstream) const {
    PCCWriteToBuffer<uint8_t>(uint8_t(attributeParams.numberOfNearestNeighborsInPrediction),
                              bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(attributeParams.levelOfDetailCount), bitstream.buffer,
                              bitstream.size);
    for (size_t lodIndex = 0; lodIndex < attributeParams.levelOfDetailCount; ++lodIndex) {
      const uint32_t d2 = uint32_t(attributeParams.dist2[lodIndex]);
      PCCWriteToBuffer<uint32_t>(d2, bitstream.buffer, bitstream.size);
    }
    for (size_t lodIndex = 0; lodIndex < attributeParams.levelOfDetailCount; ++lodIndex) {
      const uint16_t qs = uint16_t(attributeParams.quantizationSteps[lodIndex]);
      PCCWriteToBuffer<uint16_t>(qs, bitstream.buffer, bitstream.size);
    }
    for (size_t lodIndex = 0; lodIndex < attributeParams.levelOfDetailCount; ++lodIndex) {
      const uint16_t dz = uint16_t(attributeParams.quantizationDeadZoneSizes[lodIndex]);
      PCCWriteToBuffer<uint16_t>(dz, bitstream.buffer, bitstream.size);
    }
    return 0;
  }

  int encodeColors(const PCCAttributeEncodeParamaters &colorParams, PCCBitstream &bitstream) {
    const size_t pointCount = pointCloud.getPointCount();
    assert(pointCount == predictors.size());
    uint32_t maxAttributeValueDiff0 = 0;
    uint32_t maxAttributeValueDiff1 = 0;
    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      const auto &predictor = predictors[predictorIndex];
      if (predictor.maxNeighborCount) {
        const PCCColor3B color = pointCloud.getColor(predictor.index);
        const PCCColor3B predictedColor = predictor.predictColor(pointCloud);
        const int64_t quantAttValue = color[0];
        const int64_t quantPredAttValue = predictedColor[0];
        const uint32_t diffAttValue =
            uint32_t(o3dgc::IntToUInt(long(quantAttValue - quantPredAttValue)));
        if (maxAttributeValueDiff0 < diffAttValue) {
          maxAttributeValueDiff0 = diffAttValue;
        }
        for (size_t k = 1; k < 3; ++k) {
          const int64_t quantAttValue = color[k];
          const int64_t quantPredAttValue = predictedColor[k];
          const uint32_t diffAttValue =
              uint32_t(o3dgc::IntToUInt(long(quantAttValue - quantPredAttValue)));
          if (maxAttributeValueDiff1 < diffAttValue) {
            maxAttributeValueDiff1 = diffAttValue;
          }
        }
      }
    }

    uint32_t adaptiveDataModelAlphabetSizeDiff0 = 0;
    o3dgc::Adaptive_Bit_Model binaryModelDiff0;
    o3dgc::Adaptive_Data_Model multiSymbolModelDiff0;
    if (maxAttributeValueDiff0) {
      adaptiveDataModelAlphabetSizeDiff0 = 1 + maxAttributeValueDiff0;
      if (adaptiveDataModelAlphabetSizeDiff0 > PCCTMC3AdaptiveDataModelAlphabetMaxSize) {
        adaptiveDataModelAlphabetSizeDiff0 = PCCTMC3AdaptiveDataModelAlphabetMaxSize;
      }
      multiSymbolModelDiff0.set_alphabet(adaptiveDataModelAlphabetSizeDiff0 + 1);
      binaryModelDiff0.reset();
    }

    uint32_t adaptiveDataModelAlphabetSizeDiff1 = 0;
    std::vector<o3dgc::Adaptive_Bit_Model> binaryModelDiff1;
    std::vector<o3dgc::Adaptive_Data_Model> multiSymbolModelDiff1;
    if (maxAttributeValueDiff1) {
      adaptiveDataModelAlphabetSizeDiff1 = 1 + maxAttributeValueDiff1;
      if (adaptiveDataModelAlphabetSizeDiff1 > PCCTMC3AdaptiveDataModelAlphabetMaxSize) {
        adaptiveDataModelAlphabetSizeDiff1 = PCCTMC3AdaptiveDataModelAlphabetMaxSize;
      }
      multiSymbolModelDiff1.resize(PCCTMC3Diff1AdaptiveDataModelCount);
      binaryModelDiff1.resize(PCCTMC3Diff1AdaptiveDataModelCount);
      for (size_t m = 0; m < PCCTMC3Diff1AdaptiveDataModelCount; ++m) {
        multiSymbolModelDiff1[m].set_alphabet(adaptiveDataModelAlphabetSizeDiff1 + 1);
        binaryModelDiff1[m].reset();
      }
    }

    uint64_t startSize = bitstream.size;
    bitstream.size += 4;  // placehoder for bitstream size
    PCCWriteToBuffer<uint32_t>(maxAttributeValueDiff0, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint32_t>(maxAttributeValueDiff1, bitstream.buffer, bitstream.size);

    o3dgc::Arithmetic_Codec arithmeticEncoder;
    arithmeticEncoder.set_buffer(static_cast<uint32_t>(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticEncoder.start_encoder();
    o3dgc::Static_Bit_Model binaryModel0;
    o3dgc::Adaptive_Data_Model neighborCountModel;
    neighborCountModel.set_alphabet(uint32_t(colorParams.numberOfNearestNeighborsInPrediction));
    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      const auto &predictor = predictors[predictorIndex];
      const PCCColor3B color = pointCloud.getColor(predictor.index);

      if (!predictor.neighborCount) {  // no prediction
        for (size_t k = 0; k < 3; ++k) {
          encodeAbsUInt32(color[k], 8, arithmeticEncoder, binaryModel0);
        }
      } else {
        const PCCColor3B predictedColor = predictor.predictColor(pointCloud);
        const size_t lodIndex = predictor.levelOfDetailIndex;
        const int64_t qs = colorParams.quantizationSteps[lodIndex];
        const int64_t dz = colorParams.quantizationDeadZoneSizes[lodIndex];
        arithmeticEncoder.encode(uint32_t(predictor.neighborCount - 1), neighborCountModel);

        const int64_t quantAttValue = color[0];
        const int64_t quantPredAttValue = predictedColor[0];
        const int64_t delta = PCCQuantization(quantAttValue - quantPredAttValue, qs, dz);
        const uint32_t attValue0 = uint32_t(o3dgc::IntToUInt(long(delta)));
        const uint32_t modelIndex = (attValue0 < PCCTMC3Diff1AdaptiveDataModelCount)
                                        ? attValue0
                                        : PCCTMC3Diff1AdaptiveDataModelCount - 1;
        encodeDiff0UInt32(attValue0, maxAttributeValueDiff0, adaptiveDataModelAlphabetSizeDiff0,
                          arithmeticEncoder, multiSymbolModelDiff0, binaryModelDiff0, binaryModel0);

        const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs, dz);
        const int64_t reconstructedQuantAttValue = quantPredAttValue + reconstructedDelta;
        PCCColor3B reconstructedColor;
        reconstructedColor[0] =
            uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));
        for (size_t k = 1; k < 3; ++k) {
          const int64_t quantAttValue = color[k];
          const int64_t quantPredAttValue = predictedColor[k];
          const int64_t delta = PCCQuantization(quantAttValue - quantPredAttValue, qs, dz);
          const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs, dz);
          const int64_t reconstructedQuantAttValue = quantPredAttValue + reconstructedDelta;
          const uint32_t attValue1 = uint32_t(o3dgc::IntToUInt(long(delta)));
          encodeDiff1UInt32(attValue1, modelIndex, maxAttributeValueDiff1,
                            adaptiveDataModelAlphabetSizeDiff1, arithmeticEncoder,
                            multiSymbolModelDiff1, binaryModelDiff1, binaryModel0);
          reconstructedColor[k] =
              uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));
        }
        pointCloud.setColor(predictor.index, reconstructedColor);
      }
    }
    uint32_t compressedBitstreamSize = arithmeticEncoder.stop_encoder();
    bitstream.size += compressedBitstreamSize;
    PCCWriteToBuffer<uint32_t>(compressedBitstreamSize, bitstream.buffer, startSize);
    return 0;
  }
  void computeColorPredictionWeights(const PCCAttributeEncodeParamaters &attributeParams) {
    const size_t pointCount = predictors.size();
    double bestCost = std::numeric_limits<double>::max();
    size_t bestNeighborCount = 1;
    for (size_t neighborCount = 1;
         neighborCount <= attributeParams.numberOfNearestNeighborsInPrediction; ++neighborCount) {
      double cost = 0.0;
      for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
        auto &predictor = predictors[predictorIndex];
        if (!predictor.maxNeighborCount) {
          continue;
        }
        predictor.computeWeights(neighborCount);
        const PCCColor3B color = pointCloud.getColor(predictor.index);
        const PCCColor3B predictedColor = predictor.predictColor(pointCloud);

        for (int k = 0; k < 3; ++k) {
          const size_t lodIndex = predictor.levelOfDetailIndex;
          const int64_t quantAttValue = color[k];
          const int64_t quantPredAttValue = predictedColor[k];
          const int64_t delta0 =
              PCCQuantization(quantAttValue - quantPredAttValue,
                              int64_t(attributeParams.quantizationSteps[lodIndex]),
                              int64_t(attributeParams.quantizationDeadZoneSizes[lodIndex]));
          const uint32_t delta1 = uint32_t(o3dgc::IntToUInt(long(delta0)));
          cost += delta1;
        }
      }

      if (cost < bestCost) {
        bestCost = cost;
        bestNeighborCount = neighborCount;
      }
    }

    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      auto &predictor = predictors[predictorIndex];
      if (!predictor.maxNeighborCount) {
        continue;
      }
      predictor.computeWeights(bestNeighborCount);
    }
  }
  void computeReflectancePredictionWeights(const PCCAttributeEncodeParamaters &attributeParams) {
    const size_t pointCount = predictors.size();
    double bestCost = std::numeric_limits<double>::max();
    size_t bestNeighborCount = 1;
    for (size_t neighborCount = 1;
         neighborCount <= attributeParams.numberOfNearestNeighborsInPrediction; ++neighborCount) {
      double cost = 0.0;
      for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
        auto &predictor = predictors[predictorIndex];
        if (!predictor.maxNeighborCount) {
          continue;
        }
        predictor.computeWeights(neighborCount);
        const size_t lodIndex = predictor.levelOfDetailIndex;
        const int64_t quantAttValue = pointCloud.getReflectance(predictor.index);
        const int64_t quantPredAttValue = predictor.predictReflectance(pointCloud);
        const int64_t delta0 = PCCQuantization(
            quantAttValue - quantPredAttValue, int64_t(attributeParams.quantizationSteps[lodIndex]),
            int64_t(attributeParams.quantizationDeadZoneSizes[lodIndex]));
        const uint32_t delta1 = uint32_t(o3dgc::IntToUInt(long(delta0)));
        cost += delta1;
      }
      if (cost < bestCost) {
        bestCost = cost;
        bestNeighborCount = neighborCount;
      }
    }
    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      auto &predictor = predictors[predictorIndex];
      if (!predictor.maxNeighborCount) {
        continue;
      }
      predictor.computeWeights(bestNeighborCount);
    }
  }

  void buildPredictors(const PCCAttributeEncodeParamaters &attributeParams) {
    std::vector<uint32_t> numberOfPointsPerLOD;
    std::vector<uint32_t> indexes;
    PCCBuildPredictors(pointCloud, attributeParams.numberOfNearestNeighborsInPrediction,
                       attributeParams.levelOfDetailCount, attributeParams.dist2,
                       predictors, numberOfPointsPerLOD, indexes);
  }
  void reconstructedPointCloud(const PCCTMC3Encoder3Parameters &params,
                               PCCPointSet3 *reconstructedCloud) {
    if (reconstructedCloud == nullptr) {
      return;
    }
    const size_t pointCount = pointCloud.getPointCount();

    if (pointCloud.hasColors()) {
      reconstructedCloud->addColors();
    } else {
      pointCloud.removeColors();
    }
    if (pointCloud.hasReflectances()) {
      reconstructedCloud->addReflectances();
    } else {
      pointCloud.removeReflectances();
    }
    reconstructedCloud->resize(pointCount);

    const double minPositionQuantizationScale = 0.0000000001;
    const double invScale = fabs(params.positionQuantizationScale) > minPositionQuantizationScale
                                ? 1.0 / params.positionQuantizationScale
                                : 1.0;
    for (size_t i = 0; i < pointCount; ++i) {
      const auto quantizedPoint = pointCloud[i];
      auto &point = (*reconstructedCloud)[i];
      for (size_t k = 0; k < 3; ++k) {
         point[k] = quantizedPoint[k] * invScale + minPositions[k];
      }
      if (pointCloud.hasColors()) {
        reconstructedCloud->setColor(i, pointCloud.getColor(i));
      }
      if (pointCloud.hasReflectances()) {
        reconstructedCloud->setReflectance(i, pointCloud.getReflectance(i));
      }
    }
  }

  //-------------------------------------------------------------------------
  // Encode the number of points in a leaf node of the octree.

  int encodePositionLeafNumPoints(
    int count,
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    o3dgc::Adaptive_Bit_Model& ctxSinglePointPerBlock,
    o3dgc::Static_Bit_Model& ctxEquiProb,
    o3dgc::Adaptive_Bit_Model& ctxPointCountPerBlock
  ) {
    if (count == 1) {
      arithmeticEncoder->encode(1, ctxSinglePointPerBlock);
    } else {
      arithmeticEncoder->encode(0, ctxSinglePointPerBlock);
      arithmeticEncoder->ExpGolombEncode(
        uint32_t(count - 1), 0, ctxEquiProb, ctxPointCountPerBlock);
    }

    return count;
  }

  //-------------------------------------------------------------------------

  int encodePositions(
    const PCCTMC3Encoder3Parameters& params,
    PCCBitstream& bitstream
  ) {
    uint64_t startSize = bitstream.size;
    bitstream.size += 4;  // placehoder for bitstream size
    o3dgc::Arithmetic_Codec arithmeticEncoder;
    arithmeticEncoder.set_buffer(uint32_t(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticEncoder.start_encoder();
    o3dgc::Adaptive_Data_Model multiSymbolOccupancyModel0(257);

    o3dgc::Static_Bit_Model ctxEquiProb;
    o3dgc::Adaptive_Bit_Model ctxSinglePointPerBlock;
    o3dgc::Adaptive_Bit_Model ctxPointCountPerBlock;

    // init main fifo
    //  -- worst case size is the last level containing every input poit
    //     and each point being isolated in the previous level.
    pcc::ringbuf<PCCOctree3Node> fifo(pointCloud.getPointCount()+1);

    uint32_t maxBB = std::max({
      // todo(df): confirm minimum of 1 isn't needed
      1u, boundingBox.max[0], boundingBox.max[1], boundingBox.max[2]
    });

    // the current node dimension (log2) encompasing maxBB
    int nodeSizeLog2 = ceillog2(maxBB + 1);

    // push the first node
    PCCOctree3Node node00;
    node00.start = uint32_t(0);
    node00.end = uint32_t(pointCloud.getPointCount());
    node00.pos = uint32_t(0);
    fifo.push_back(node00);

    size_t processedPointCount = 0;
    std::vector<uint32_t> values;

    auto fifoCurrLvlEnd = fifo.end();

    for (; !fifo.empty(); fifo.pop_front()) {
      if (fifo.begin() == fifoCurrLvlEnd) {
        // transition to the next level
        fifoCurrLvlEnd = fifo.end();
        nodeSizeLog2--;
      }

      PCCOctree3Node& node0 = fifo.front();

      // split the current node into 8 children
      //  - perform an 8-way counting sort of the current node's points
      //  - (later) map to child nodes
      int childSizeLog2 = nodeSizeLog2 - 1;
      std::array<int,8> childCounts = {};
      countingSort(
        PCCPointSet3::iterator(&pointCloud, node0.start),
        PCCPointSet3::iterator(&pointCloud, node0.end),
        childCounts,
        [=](const PCCPointSet3::Proxy& proxy){
          const auto& point = *proxy;
          int bitpos = 1 << childSizeLog2;
          return
               !!(int(point[2]) & bitpos)
            | (!!(int(point[1]) & bitpos) << 1)
            | (!!(int(point[0]) & bitpos) << 2);
      });

      // determine occupancy
      int occupancy = 0;
      for (int i = 0; i < 8; i++) {
        if (childCounts[i]) {
          occupancy |= 1 << i;
        }
      }

      assert(occupancy > 0);
      arithmeticEncoder.encode(occupancy, multiSymbolOccupancyModel0);

      // when nodeSizeLog2 == 1, children are indivisible (ie leaf nodes)
      // and are immediately coded.  No further splitting occurs.
      if (nodeSizeLog2 == 1) {
        for (int i = 0; i < 8; i++) {
          if (!childCounts[i]) {
            // child is empty: skip
            continue;
          }

          // if the bitstream is configured to represent unique points,
          // no point count is sent.
          if (params.mergeDuplicatedPoints) {
            assert(childCounts[i] == 1);
            continue;
          }

          processedPointCount += encodePositionLeafNumPoints(
            childCounts[i], &arithmeticEncoder,
            ctxSinglePointPerBlock, ctxEquiProb, ctxPointCountPerBlock
          );
        }

        // leaf nodes do not get split
        continue;
      }

      // nodeSizeLog2 > 1: insert split children into fifo
      int childPointsStartIdx = node0.start;
      for (int i = 0; i < 8; i++) {
        if (!childCounts[i]) {
          // child is empty: skip
          continue;
        }

        // create new child
        fifo.emplace_back();
        auto& child = fifo.back();

        int x = !!(i & 4);
        int y = !!(i & 2);
        int z = !!(i & 1);

        child.pos[0] = node0.pos[0] + (x << childSizeLog2);
        child.pos[1] = node0.pos[1] + (y << childSizeLog2);
        child.pos[2] = node0.pos[2] + (z << childSizeLog2);

        child.start = childPointsStartIdx;
        childPointsStartIdx += childCounts[i];
        child.end = childPointsStartIdx;
      }
    }

    const uint32_t compressedBitstreamSize = arithmeticEncoder.stop_encoder();
    bitstream.size += compressedBitstreamSize;
    PCCWriteToBuffer<uint32_t>(compressedBitstreamSize, bitstream.buffer, startSize);
    return 0;
  }

  int encodePositionsHeader(const PCCTMC3Encoder3Parameters &params,
                            PCCBitstream &bitstream) const {
    PCCWriteToBuffer<uint32_t>(uint32_t(pointCloud.getPointCount()), bitstream.buffer,
                               bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasColors()), bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasReflectances()), bitstream.buffer,
                              bitstream.size);

    // todo(df): syntax element name: geometryPointsAreUnique?
    PCCWriteToBuffer<uint8_t>(
      uint8_t(params.mergeDuplicatedPoints), bitstream.buffer, bitstream.size);

    for (int k = 0; k < 3; ++k) {
      PCCWriteToBuffer<double>(minPositions[k], bitstream.buffer, bitstream.size);
    }
    for (int k = 0; k < 3; ++k) {
      PCCWriteToBuffer<uint32_t>(boundingBox.max[k], bitstream.buffer, bitstream.size);
    }
    PCCWriteToBuffer<double>(params.positionQuantizationScale, bitstream.buffer, bitstream.size);
    return 0;
  }
  void computeMinPositions(const PCCPointSet3 &inputPointCloud) {
    const size_t inputPointCount = inputPointCloud.getPointCount();
    minPositions = inputPointCloud[0];
    for (size_t i = 1; i < inputPointCount; ++i) {
      const auto point = inputPointCloud[i];
      for (int k = 0; k < 3; ++k) {
        if (minPositions[k] > point[k]) {
          minPositions[k] = point[k];
        }
      }
    }
  }
  int quantization(const PCCPointSet3 &inputPointCloud, const PCCTMC3Encoder3Parameters &params) {
    computeMinPositions(inputPointCloud);

    pointCloud.resize(0);
    const bool compressColors =
        inputPointCloud.hasColors() &&
        params.attributeEncodeParameters.find("color") != params.attributeEncodeParameters.end();
    if (compressColors) {
      pointCloud.addColors();
    } else {
      pointCloud.removeColors();
    }

    const bool compressReflectances = inputPointCloud.hasReflectances() &&
                                      params.attributeEncodeParameters.find("reflectance") !=
                                          params.attributeEncodeParameters.end();
    if (compressReflectances) {
      pointCloud.addReflectances();
    } else {
      pointCloud.removeReflectances();
    }

    const size_t inputPointCount = inputPointCloud.getPointCount();
    if (params.mergeDuplicatedPoints) {  // retain unique quantized points
      // filter points based on quantized positions
      std::set<PCCVector3<int64_t>> retainedPoints;
      for (size_t i = 0; i < inputPointCount; ++i) {
        const PCCVector3D point =
            (inputPointCloud[i] - minPositions) * params.positionQuantizationScale;
        PCCVector3<int64_t> quantizedPoint;
        quantizedPoint[0] = int64_t(std::round(point[0]));
        quantizedPoint[1] = int64_t(std::round(point[1]));
        quantizedPoint[2] = int64_t(std::round(point[2]));
        const auto it = retainedPoints.find(quantizedPoint);
        if (it == retainedPoints.end()) {
          retainedPoints.insert(quantizedPoint);
        }
      }
      // compute reconstructed point cloud
      pointCloud.resize(retainedPoints.size());
      const double invScale = 1.0 / params.positionQuantizationScale;
      size_t pointCounter = 0;
      for (const auto &quantizedPoint : retainedPoints) {
        auto &point = pointCloud[pointCounter++];
        for (size_t k = 0; k < 3; ++k) {
          point[k] = double(quantizedPoint[k]) * invScale + minPositions[k];
        }
      }

      if (compressColors) {  // transfer colors
        const auto &attributeParams = params.attributeEncodeParameters.find("color")->second;
        if (!PCCTransfertColors(inputPointCloud, int32_t(attributeParams.searchRange),
                                pointCloud)) {
          std::cout << "Error: can't transfer colors!" << std::endl;
          return -1;
        }
      }

      if (compressReflectances) {  // transfer reflectances
        const auto &attributeParams = params.attributeEncodeParameters.find("reflectance")->second;
        if (!PCCTransfertReflectances(inputPointCloud, int32_t(attributeParams.searchRange),
                                      pointCloud)) {
          std::cout << "Error: can't transfer reflectances!" << std::endl;
          return -1;
        }
      }

      // compute quantized coordinates
      pointCounter = 0;
      for (const auto &quantizedPoint : retainedPoints) {
        auto &point = pointCloud[pointCounter++];
        for (size_t k = 0; k < 3; ++k) {
          point[k] = double(quantizedPoint[k]);
        }
      }
    } else {
      // quantized point cloud
      pointCloud.resize(inputPointCount);
      for (size_t i = 0; i < inputPointCount; ++i) {
        const PCCVector3D point =
            (inputPointCloud[i] - minPositions) * params.positionQuantizationScale;
        auto &quantizedPoint = pointCloud[i];
        quantizedPoint[0] = std::round(point[0]);
        quantizedPoint[1] = std::round(point[1]);
        quantizedPoint[2] = std::round(point[2]);
        if (compressColors) {
          pointCloud.setColor(i, inputPointCloud.getColor(i));
        }
        if (compressReflectances) {
          pointCloud.setReflectance(i, inputPointCloud.getReflectance(i));
        }
      }
    }

    const size_t pointCount = pointCloud.getPointCount();
    boundingBox.min = uint32_t(0);
    boundingBox.max = uint32_t(0);
    for (size_t i = 0; i < pointCount; ++i) {
      const PCCVector3D point = pointCloud[i];
      for (int k = 0; k < 3; ++k) {
        const uint32_t coord = uint32_t(point[k]);
        if (boundingBox.max[k] < coord) {
          boundingBox.max[k] = coord;
        }
      }
    }
    return 0;
  }

 private:
  std::vector<PCCPredictor> predictors;
  PCCVector3D minPositions;
  PCCBox3<uint32_t> boundingBox;
  PCCPointSet3 pointCloud;
};
}

#endif /* PCCTMC3Encoder_h */
