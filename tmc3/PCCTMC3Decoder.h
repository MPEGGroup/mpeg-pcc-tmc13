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

#ifndef PCCTMC3Decoder_h
#define PCCTMC3Decoder_h

#include <assert.h>
#include <queue>
#include <string>

#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCTMC3Common.h"

#include "ArithmeticCodec.h"
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
  int decompressWithLosslessGeometry(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
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
    return 0;
  }
  int decompress(PCCBitstream &bitstream, PCCPointSet3 &pointCloud, const bool roundOutputPositions) {
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

    inverseQuantization(pointCloud, roundOutputPositions);
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
    return 0;
  }

  //-------------------------------------------------------------------------
  // Decode the number of points in a leaf node of the octree.

  int decodePositionLeafNumPoints(
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    o3dgc::Adaptive_Bit_Model& ctxSinglePointPerBlock,
    o3dgc::Static_Bit_Model& ctxEquiProb,
    o3dgc::Adaptive_Bit_Model& ctxPointCountPerBlock
  ) {
    const bool isSinglePoint =
        arithmeticDecoder->decode(ctxSinglePointPerBlock) != 0;

    int count = 1;
    if (!isSinglePoint) {
      count = 1 + arithmeticDecoder->ExpGolombDecode(
          0, ctxEquiProb, ctxPointCountPerBlock);
    }

    return count;
  }

  //-------------------------------------------------------------------------

  int decodePositions(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    uint32_t compressedBitstreamSize;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, compressedBitstreamSize, bitstream.size);
    o3dgc::Arithmetic_Codec arithmeticDecoder;
    arithmeticDecoder.set_buffer(uint32_t(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticDecoder.start_decoder();
    o3dgc::Adaptive_Data_Model multiSymbolOccupancyModel0(257);

    o3dgc::Static_Bit_Model ctxEquiProb;
    o3dgc::Adaptive_Bit_Model ctxSinglePointPerBlock;
    o3dgc::Adaptive_Bit_Model ctxPointCountPerBlock;

    // init main fifo
    std::deque<PCCOctree3Node> fifo;

    uint32_t maxBB = std::max({
      // todo(df): confirm minimum of 1 isn't needed
      1u, boundingBox.max[0], boundingBox.max[1], boundingBox.max[2]
    });
    // round up to the next highest power of 2 minus 1 >= maxBB.
    maxBB = ceilpow2(maxBB + 1) - 1;

    // push the first node
    PCCOctree3Node node00;
    node00.start = uint32_t(0);
    node00.end = uint32_t(pointCloud.getPointCount());
    node00.boundingBox.min = uint32_t(0);
    node00.boundingBox.max = maxBB;
    fifo.push_back(node00);

    size_t processedPointCount = 0;
    std::vector<uint32_t> values;

    for (; !fifo.empty(); fifo.pop_front()) {
      PCCOctree3Node& node0 = fifo.front();

      const PCCVector3<uint32_t> range = node0.boundingBox.max - node0.boundingBox.min;

      // decode the points from a leaf node at maximal depth
      if (range[0] == 0 && range[1] == 0 && range[2] == 0) {
        int numPoints = decodePositionLeafNumPoints(
          &arithmeticDecoder,
          ctxSinglePointPerBlock, ctxEquiProb, ctxPointCountPerBlock
        );

        const PCCVector3D point(
          node0.boundingBox.min[0], node0.boundingBox.min[1],
          node0.boundingBox.min[2]);

        for (int i = 0; i < numPoints; ++i)
          pointCloud[processedPointCount++] = point;

        // leaf nodes do not get split
        continue;
      }

      // split the current node based on occupancy.
      const uint32_t occupancy = arithmeticDecoder.decode(multiSymbolOccupancyModel0);
      assert(occupancy > 0);

      // split the current node
      const auto& bbox = node0.boundingBox;
      uint32_t splitX = (bbox.max[0] + bbox.min[0]) >> 1;
      uint32_t splitY = (bbox.max[1] + bbox.min[1]) >> 1;
      uint32_t splitZ = (bbox.max[2] + bbox.min[2]) >> 1;

      struct {
        uint32_t min;
        uint32_t max;
      } splitBounds[][2] = {
        {{bbox.min[0], splitX}, {splitX+1, bbox.max[0]}}, // x
        {{bbox.min[1], splitY}, {splitY+1, bbox.max[1]}}, // y
        {{bbox.min[2], splitZ}, {splitZ+1, bbox.max[2]}}  // z
      };

      for (int i = 0; i < 8; i++) {
        uint32_t mask = 1 << i;
        if (!(occupancy & mask)) {
          // child is empty: skip
          continue;
        }

        // create new child and set bounding box.
        fifo.emplace_back();
        auto& child = fifo.back();

        int x = !!(i & 4);
        int y = !!(i & 2);
        int z = !!(i & 1);

        child.boundingBox.min[0] = splitBounds[0][x].min;
        child.boundingBox.max[0] = splitBounds[0][x].max;
        child.boundingBox.min[1] = splitBounds[1][y].min;
        child.boundingBox.max[1] = splitBounds[1][y].max;
        child.boundingBox.min[2] = splitBounds[2][z].min;
        child.boundingBox.max[2] = splitBounds[2][z].max;
      }
    }

    arithmeticDecoder.stop_decoder();
    bitstream.size += compressedBitstreamSize;
    return 0;
  }
  void inverseQuantization(PCCPointSet3 &pointCloud, const bool roundOutputPositions) {
    const size_t pointCount = pointCloud.getPointCount();
    const double invScale = 1.0 / positionQuantizationScale;
      
      if (roundOutputPositions) {
          for (size_t i = 0; i < pointCount; ++i) {
              auto &point = pointCloud[i];
              for (size_t k = 0; k < 3; ++k) {
                  point[k] = std::round(point[k] * invScale + minPositions[k]);
              }
          }
      } else {
          for (size_t i = 0; i < pointCount; ++i) {
              auto &point = pointCloud[i];
              for (size_t k = 0; k < 3; ++k) {
                  point[k] = point[k] * invScale + minPositions[k];
              }
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
