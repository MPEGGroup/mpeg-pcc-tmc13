/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * <OWNER> = Apple Inc.
 * <ORGANIZATION> = Apple Inc.
 * <YEAR> = 2017
 *
 * Copyright (c) 2017, Apple Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PCCTMC3Decoder_h
#define PCCTMC3Decoder_h

#include <assert.h>
#include <queue>
#include <string>

#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCTMC3Common.h"

#include "ArithmeticCodec.h"

//#define TMCV3_DECODER_VERBOSE

namespace pcc {

class PCCTMC3Decoder3 {
 public:
  PCCTMC3Decoder3() { init(); }
  PCCTMC3Decoder3(const PCCTMC3Decoder3 &) = default;
  PCCTMC3Decoder3 &operator=(const PCCTMC3Decoder3 &rhs) = default;
  ~PCCTMC3Decoder3() = default;

  void init() {
    positionQuantizationScale = 1.0;
    minPositions = 0.0;
    boundingBox.min = uint32_t(0);
    boundingBox.max = uint32_t(0);
    numberOfNearestNeighborsInPrediction = 0;
    levelOfDetailCount = 0;
    dist2.clear();
    quantizationSteps.clear();
    quantizationDeadZoneSizes.clear();
    predictors.clear();
  }

  int decompress(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    init();
    uint32_t magicNumber = 0;
    uint32_t formatVersion = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, magicNumber, bitstream.size);
    if (magicNumber != PCCTMC3MagicNumber) {
      std::cout << "Error: corrupted bistream!" << std::endl;
      return -1;
    }
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, formatVersion, bitstream.size);
    if (formatVersion != PCCTMC3FormatVersion) {
      std::cout << "Error: bistream version not supported!" << std::endl;
      return -1;
    }

    uint64_t positionsSize = bitstream.size;
    if (int ret = decodePositionsHeader(bitstream, pointCloud)) {
      return ret;
    }
    if (int ret = decodePositions(bitstream, pointCloud)) {
      return ret;
    }
    positionsSize = bitstream.size - positionsSize;
    std::cout << "positions bitstream size " << positionsSize << " B" << std::endl;

    if (pointCloud.hasColors()) {
      uint64_t colorsSize = bitstream.size;
      if (int ret = decodeAttributeHeader("color", bitstream)) {
        return ret;
      }
      buildPredictors(pointCloud);
      if (int ret = decodeColors(bitstream, pointCloud)) {
        return ret;
      }
      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B" << std::endl;
      std::cout << std::endl;
    }

    if (pointCloud.hasReflectances()) {
      uint64_t reflectancesSize = bitstream.size;
      if (int ret = decodeAttributeHeader("reflectance", bitstream)) {
        return ret;
      }
      buildPredictors(pointCloud);
      if (int ret = decodeReflectances(bitstream, pointCloud)) {
        return ret;
      }
      reflectancesSize = bitstream.size - reflectancesSize;
      std::cout << "reflectances bitstream size " << reflectancesSize << " B" << std::endl;
    }

    inverseQuantization(pointCloud);

    return 0;
  }

 private:
  static uint32_t decodeAbsUInt32(const uint32_t bitCount,
                                  o3dgc::Arithmetic_Codec &arithmeticDecoder,
                                  o3dgc::Static_Bit_Model &binaryModel0) {
    uint32_t decodedValue = 0;
    for (uint32_t i = 0; i < bitCount; ++i) {
      decodedValue += (arithmeticDecoder.decode(binaryModel0) << i);
    }
    return PCCFromLittleEndian<uint32_t>(decodedValue);
  }
  static uint32_t decodeDiff0UInt32(const uint32_t maxAttributeValueDiff0,
                                    const uint32_t adaptiveDataModelAlphabetSizeDiff0,
                                    o3dgc::Arithmetic_Codec &arithmeticDecoder,
                                    o3dgc::Adaptive_Data_Model &multiSymbolModelDiff0,
                                    o3dgc::Adaptive_Bit_Model &binaryModelDiff0,
                                    o3dgc::Static_Bit_Model &binaryModel0) {
    if (!maxAttributeValueDiff0) {
      return 0;
    } else {
      uint32_t value = arithmeticDecoder.decode(multiSymbolModelDiff0);
      if (value == adaptiveDataModelAlphabetSizeDiff0) {
        value += arithmeticDecoder.ExpGolombDecode(0, binaryModel0, binaryModelDiff0);
      }
      return value;
    }
  }
  static uint32_t decodeDiff1UInt32(const uint32_t modelIndex,
                                    const uint32_t maxAttributeValueDiff1,
                                    const uint32_t adaptiveDataModelAlphabetSizeDiff1,
                                    o3dgc::Arithmetic_Codec &arithmeticDecoder,
                                    std::vector<o3dgc::Adaptive_Data_Model> &multiSymbolModelDiff1,
                                    std::vector<o3dgc::Adaptive_Bit_Model> &binaryModelDiff1,
                                    o3dgc::Static_Bit_Model &binaryModel0) {
    if (!maxAttributeValueDiff1) {
      return 0;
    } else {
      uint32_t value = arithmeticDecoder.decode(multiSymbolModelDiff1[modelIndex]);
      if (value == adaptiveDataModelAlphabetSizeDiff1) {
        value += arithmeticDecoder.ExpGolombDecode(0, binaryModel0, binaryModelDiff1[modelIndex]);
      }
      return value;
    }
  }
  int decodeReflectances(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    const size_t pointCount = pointCloud.getPointCount();
    assert(pointCount == predictors.size());

    uint32_t compressedBitstreamSize = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, compressedBitstreamSize, bitstream.size);
    uint32_t maxAttributeValueDiff0 = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, maxAttributeValueDiff0, bitstream.size);

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

    o3dgc::Arithmetic_Codec arithmeticDecoder;
    arithmeticDecoder.set_buffer(uint32_t(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticDecoder.start_decoder();
    o3dgc::Static_Bit_Model binaryModel0;
    o3dgc::Adaptive_Data_Model neighborCountModel;
    neighborCountModel.set_alphabet(uint32_t(numberOfNearestNeighborsInPrediction));

    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      auto &predictor = predictors[predictorIndex];
      uint16_t &reflectance = pointCloud.getReflectance(predictor.index);
      if (!predictor.neighborCount) {
        reflectance = decodeAbsUInt32(16, arithmeticDecoder, binaryModel0);
      } else {
        const size_t neighborCount = arithmeticDecoder.decode(neighborCountModel) + 1;
        predictor.computeWeights(neighborCount);
        const size_t lodIndex = predictor.levelOfDetailIndex;
        const int64_t qs = quantizationSteps[lodIndex];
        const int64_t dz = quantizationDeadZoneSizes[lodIndex];
        const uint32_t attValue0 = decodeDiff0UInt32(
            maxAttributeValueDiff0, adaptiveDataModelAlphabetSizeDiff0, arithmeticDecoder,
            multiSymbolModelDiff0, binaryModelDiff0, binaryModel0);
        const int64_t quantPredAttValue = predictor.predictReflectance(pointCloud);
        const int64_t delta = PCCInverseQuantization(o3dgc::UIntToInt(attValue0), qs, dz);
        const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
        reflectance = uint16_t(PCCClip(reconstructedQuantAttValue, int64_t(0),
                                       int64_t(std::numeric_limits<uint16_t>::max())));
      }
    }
    arithmeticDecoder.stop_decoder();
    bitstream.size += compressedBitstreamSize;
    return 0;
  }

  int decodeColors(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    const size_t pointCount = pointCloud.getPointCount();
    assert(pointCount == predictors.size());

    uint32_t compressedBitstreamSize = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, compressedBitstreamSize, bitstream.size);
    uint32_t maxAttributeValueDiff0 = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, maxAttributeValueDiff0, bitstream.size);
    uint32_t maxAttributeValueDiff1 = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, maxAttributeValueDiff1, bitstream.size);

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

    o3dgc::Arithmetic_Codec arithmeticDecoder;
    arithmeticDecoder.set_buffer(uint32_t(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticDecoder.start_decoder();
    o3dgc::Static_Bit_Model binaryModel0;
    o3dgc::Adaptive_Data_Model neighborCountModel;
    neighborCountModel.set_alphabet(uint32_t(numberOfNearestNeighborsInPrediction));

    for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
      auto &predictor = predictors[predictorIndex];
      PCCColor3B &color = pointCloud.getColor(predictor.index);
      if (!predictor.neighborCount) {
        for (size_t k = 0; k < 3; ++k) {
          color[k] = decodeAbsUInt32(8, arithmeticDecoder, binaryModel0);
        }
      } else {
        const size_t neighborCount = arithmeticDecoder.decode(neighborCountModel) + 1;
        predictor.computeWeights(neighborCount);
        const PCCColor3B predictedColor = predictor.predictColor(pointCloud);
        const size_t lodIndex = predictor.levelOfDetailIndex;
        const int64_t qs = quantizationSteps[lodIndex];
        const int64_t dz = quantizationDeadZoneSizes[lodIndex];
        const uint32_t attValue0 = decodeDiff0UInt32(
            maxAttributeValueDiff0, adaptiveDataModelAlphabetSizeDiff0, arithmeticDecoder,
            multiSymbolModelDiff0, binaryModelDiff0, binaryModel0);
        const uint32_t modelIndex = (attValue0 < PCCTMC3Diff1AdaptiveDataModelCount)
                                        ? attValue0
                                        : PCCTMC3Diff1AdaptiveDataModelCount - 1;
        const int64_t quantPredAttValue = predictedColor[0];
        const int64_t delta = PCCInverseQuantization(o3dgc::UIntToInt(attValue0), qs, dz);
        const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
        color[0] = uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));

        for (size_t k = 1; k < 3; ++k) {
          const uint32_t attValue = decodeDiff1UInt32(
              modelIndex, maxAttributeValueDiff1, adaptiveDataModelAlphabetSizeDiff1,
              arithmeticDecoder, multiSymbolModelDiff1, binaryModelDiff1, binaryModel0);
          const int64_t quantPredAttValue = predictedColor[k];
          const int64_t delta = PCCInverseQuantization(o3dgc::UIntToInt(attValue), qs, dz);
          const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
          color[k] = uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));
        }
      }
    }
    arithmeticDecoder.stop_decoder();
    bitstream.size += compressedBitstreamSize;
    return 0;
  }
  void buildPredictors(const PCCPointSet3 &pointCloud) {
    std::vector<uint32_t> numberOfPointsPerLOD;
    std::vector<uint32_t> indexes;
    PCCBuildPredictors(pointCloud, numberOfNearestNeighborsInPrediction, levelOfDetailCount, dist2,
                       predictors, numberOfPointsPerLOD, indexes);
#ifdef TMCV3_DECODER_VERBOSE
    std::cout << "building predictors" << std::endl;
    for (size_t lodIndex = 0; lodIndex < numberOfPointsPerLOD.size(); ++lodIndex) {
      std::cout << "\t LOD " << lodIndex << " -> " << numberOfPointsPerLOD[lodIndex] << std::endl;
    }
    for (size_t lodIndex = 0; lodIndex < numberOfPointsPerLOD.size(); ++lodIndex) {
      const size_t pointCount = numberOfPointsPerLOD[lodIndex];
      PCCPointSet3 pointSet;
      pointSet.resize(pointCount);
      for (size_t i = 0; i < pointCount; ++i) {
        pointSet[i] = pointCloud[indexes[i]];
      }
      std::stringstream fileName;
      fileName << "LOD_dec" << lodIndex << ".ply";
      pointSet.write(fileName.str().c_str());
    }
    std::cout << std::endl;
#endif  // TMCV3_DECODER_VERBOSE
  }

  int decodeAttributeHeader(const std::string &attributeName, PCCBitstream &bitstream) {
    uint8_t numberOfNearestNeighborsCount = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, numberOfNearestNeighborsCount, bitstream.size);
    numberOfNearestNeighborsInPrediction = numberOfNearestNeighborsCount;

    uint8_t lodCount = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, lodCount, bitstream.size);
    levelOfDetailCount = lodCount;

    dist2.resize(levelOfDetailCount);
    for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
      uint32_t d2 = 0;
      PCCReadFromBuffer<uint32_t>(bitstream.buffer, d2, bitstream.size);
      dist2[lodIndex] = d2;
    }

    quantizationSteps.resize(levelOfDetailCount);
    for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
      uint16_t qs = 0;
      PCCReadFromBuffer<uint16_t>(bitstream.buffer, qs, bitstream.size);
      quantizationSteps[lodIndex] = qs;
    }

    quantizationDeadZoneSizes.resize(levelOfDetailCount);
    for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
      uint16_t dz = 0;
      PCCReadFromBuffer<uint16_t>(bitstream.buffer, dz, bitstream.size);
      quantizationDeadZoneSizes[lodIndex] = dz;
    }

#ifdef TMCV3_DECODER_VERBOSE
    std::cout << attributeName << " header" << std::endl;
    std::cout << "\t numberOfNearestNeighborsInPrediction " << numberOfNearestNeighborsInPrediction
              << std::endl;
    std::cout << "\t levelOfDetailCount                   " << levelOfDetailCount << std::endl;
    std::cout << "\t dist2                                ";
    for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
      std::cout << dist2[lodIndex] << " ";
    }
    std::cout << std::endl;
    std::cout << "\t quantizationSteps                    ";
    for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
      std::cout << quantizationSteps[lodIndex] << " ";
    }
    std::cout << std::endl;
    std::cout << "\t quantizationDeadZoneSizes            ";
    for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
      std::cout << quantizationDeadZoneSizes[lodIndex] << " ";
    }
    std::cout << std::endl << std::endl;
#endif  // TMCV3_DECODER_VERBOSE
    return 0;
  }

  int decodePositionsHeader(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    uint32_t pointCount = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, pointCount, bitstream.size);
    uint8_t hasColors = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, hasColors, bitstream.size);
    uint8_t hasReflectances = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, hasReflectances, bitstream.size);
    if (hasColors) {
      pointCloud.addColors();
    } else {
      pointCloud.removeColors();
    }
    if (hasReflectances) {
      pointCloud.addReflectances();
    } else {
      pointCloud.removeReflectances();
    }
    pointCloud.resize(pointCount);

    for (int k = 0; k < 3; ++k) {
      PCCReadFromBuffer<double>(bitstream.buffer, minPositions[k], bitstream.size);
    }
    for (int k = 0; k < 3; ++k) {
      PCCReadFromBuffer<uint32_t>(bitstream.buffer, boundingBox.max[k], bitstream.size);
    }
    PCCReadFromBuffer<double>(bitstream.buffer, positionQuantizationScale, bitstream.size);
    const double minPositionQuantizationScale = 0.0000000001;
    if (fabs(positionQuantizationScale) < minPositionQuantizationScale) {
      positionQuantizationScale = 1.0;
    }

#ifdef TMCV3_DECODER_VERBOSE
    std::cout << "Positions Header" << std::endl;
    std::cout << "\t pointCount                " << pointCloud.getPointCount() << std::endl;
    std::cout << "\t hasColors                 " << pointCloud.hasColors() << std::endl;
    std::cout << "\t hasReflectances           " << pointCloud.hasReflectances() << std::endl;
    std::cout << "\t minPositions              " << minPositions[0] << ", " << minPositions[1]
              << ", " << minPositions[2] << std::endl;
    std::cout << "\t boundingBox               (" << boundingBox.min[0] << ", "
              << boundingBox.min[1] << ", " << boundingBox.min[2] << ") (" << boundingBox.max[0]
              << ", " << boundingBox.max[1] << ", " << boundingBox.max[2] << ")" << std::endl;
    std::cout << "\t positionQuantizationScale " << positionQuantizationScale << std::endl;
    std::cout << std::endl;
#endif  // TMCV3_DECODER_VERBOSE
    return 0;
  }
  int decodePositions(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    uint32_t compressedBitstreamSize;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, compressedBitstreamSize, bitstream.size);
    o3dgc::Arithmetic_Codec arithmeticDecoder;
    arithmeticDecoder.set_buffer(uint32_t(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticDecoder.start_decoder();
    o3dgc::Adaptive_Data_Model multiSymbolOccupacyModel0(257);
    std::queue<PCCOctree3Node> fifo;
    uint32_t maxBB = std::max(std::max(boundingBox.max[0], boundingBox.max[1]),
                              std::max(uint32_t(1), boundingBox.max[2]));
    // round up to the next highest power of 2
    maxBB--;
    maxBB |= maxBB >> 1;
    maxBB |= maxBB >> 2;
    maxBB |= maxBB >> 4;
    maxBB |= maxBB >> 8;
    maxBB |= maxBB >> 16;

    // push the first node
    PCCOctree3Node node0;
    node0.start = uint32_t(0);
    node0.end = uint32_t(pointCloud.getPointCount());
    node0.boundingBox.min = uint32_t(0);
    node0.boundingBox.max = maxBB;
    fifo.push(node0);

    size_t processedPointCount = 0;
    std::vector<uint32_t> values;
    o3dgc::Adaptive_Bit_Model singlePointPerBlock;
    o3dgc::Static_Bit_Model bModel0;
    o3dgc::Adaptive_Bit_Model bModelPointCountPerBlock;
    while (!fifo.empty()) {
      const PCCOctree3Node node0 = fifo.front();
      fifo.pop();
      const PCCVector3<uint32_t> range = node0.boundingBox.max - node0.boundingBox.min;
      if (range[0] == 0 && range[1] == 0 && range[2] == 0) {
        const bool isSinglePoint = arithmeticDecoder.decode(singlePointPerBlock) != 0;
        size_t count = 1;
        if (!isSinglePoint) {
          count = 1 + arithmeticDecoder.ExpGolombDecode(0, bModel0, bModelPointCountPerBlock);
        }
        const PCCVector3D point(node0.boundingBox.min[0], node0.boundingBox.min[1],
                                node0.boundingBox.min[2]);
        for (size_t i = 0; i < count; ++i) {
          pointCloud[processedPointCount++] = point;
        }
      } else {
        std::queue<PCCOctree3Node> tmpFifo;
        tmpFifo.push(node0);
        for (size_t splitAxis = 0; splitAxis < 3; ++splitAxis) {
          const size_t nodeCount = tmpFifo.size();
          for (size_t n = 0; n < nodeCount; ++n) {
            const PCCOctree3Node node = tmpFifo.front();
            tmpFifo.pop();
            const uint32_t splitValue =
                (node.boundingBox.max[splitAxis] + node.boundingBox.min[splitAxis]) / 2;
            PCCOctree3Node nodeLeft = node;
            nodeLeft.boundingBox.max[splitAxis] = splitValue;
            tmpFifo.push(nodeLeft);
            PCCOctree3Node nodeRight = node;
            nodeRight.boundingBox.min[splitAxis] = splitValue + 1;
            tmpFifo.push(nodeRight);
          }
        }

        assert(tmpFifo.size() == 8);
        const uint32_t occupacy = arithmeticDecoder.decode(multiSymbolOccupacyModel0);
        assert(occupacy > 0);

        uint32_t mask = 1;
        while (!tmpFifo.empty()) {
          const PCCOctree3Node node = tmpFifo.front();

          tmpFifo.pop();
          if (occupacy & mask) {
            fifo.push(node);
          }
          mask <<= 1;
        }
      }
    }
    arithmeticDecoder.stop_decoder();
    bitstream.size += compressedBitstreamSize;
    return 0;
  }
  void inverseQuantization(PCCPointSet3 &pointCloud) {
    const size_t pointCount = pointCloud.getPointCount();
    const double invScale = 1.0 / positionQuantizationScale;
    for (size_t i = 0; i < pointCount; ++i) {
      auto &point = pointCloud[i];
      for (size_t k = 0; k < 3; ++k) {
        point[k] = point[k] * invScale + minPositions[k];
        //          point[k] = std::round(point[k] * invScale + minPositions[k]);
      }
    }
  }

 private:
  std::vector<PCCPredictor> predictors;
  std::vector<size_t> dist2;
  std::vector<int64_t> quantizationSteps;
  std::vector<int64_t> quantizationDeadZoneSizes;
  PCCVector3D minPositions;
  PCCBox3<uint32_t> boundingBox;
  double positionQuantizationScale;
  size_t numberOfNearestNeighborsInPrediction;
  size_t levelOfDetailCount;
};
}

#endif /* PCCTMC3Decoder_h */
