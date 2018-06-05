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
#include "ringbuf.h"

#include "AttributeDecoder.h"
#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCTMC3Common.h"
#include "osspecific.h"

#include "ArithmeticCodec.h"
#include "tables.h"

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

    AttributeDecoder attrDecoder;

    if (pointCloud.hasColors()) {
      uint64_t colorsSize = bitstream.size;

      attrDecoder.decodeHeader("color", bitstream);
      attrDecoder.buildPredictors(pointCloud);
      attrDecoder.decodeColors(bitstream, pointCloud);

      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B" << std::endl;
      std::cout << std::endl;
    }

    if (pointCloud.hasReflectances()) {
      uint64_t reflectancesSize = bitstream.size;

      attrDecoder.decodeHeader("reflectance", bitstream);
      attrDecoder.buildPredictors(pointCloud);
      attrDecoder.decodeReflectances(bitstream, pointCloud);

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

    decodePositionsHeader(bitstream, pointCloud);
    decodePositions(bitstream, pointCloud);

    positionsSize = bitstream.size - positionsSize;
    std::cout << "positions bitstream size " << positionsSize << " B" << std::endl;

    if (pointCloud.hasColors()) {
      AttributeDecoder attrDecoder;
      uint64_t colorsSize = bitstream.size;

      attrDecoder.decodeHeader("color", bitstream);
      attrDecoder.buildPredictors(pointCloud);
      attrDecoder.decodeColors(bitstream, pointCloud);

      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B" << std::endl;
      std::cout << std::endl;
    }

    if (pointCloud.hasReflectances()) {
      AttributeDecoder attrDecoder;
      uint64_t reflectancesSize = bitstream.size;

      attrDecoder.decodeHeader("reflectance", bitstream);
      attrDecoder.buildPredictors(pointCloud);
      attrDecoder.decodeReflectances(bitstream, pointCloud);

      reflectancesSize = bitstream.size - reflectancesSize;
      std::cout << "reflectances bitstream size " << reflectancesSize << " B" << std::endl;
    }

    inverseQuantization(pointCloud, roundOutputPositions);
    return 0;
  }

  int decompressWithTrisoupGeometry(PCCBitstream &bitstream, PCCPointSet3 &pointCloud,
                                    const bool roundOutputPositions) {
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
    size_t depth;
    size_t level;
    double intToOrigScale;
    PCCPoint3D intToOrigTranslation;

    decodeTrisoupHeader(
      bitstream, pointCloud,
      depth, level, intToOrigScale, intToOrigTranslation);

    // Write out TMC1 geometry bitstream from the converged TMC13 bitstream.
    uint32_t binSize = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, binSize, bitstream.size);
    std::ofstream fout("outbin/outbin_0000.bin", std::ios::binary);
    if (!fout.is_open()) {
      // maybe because directory does not exist; try again...
      pcc::mkdir("outbin");
      fout.open("outbin/outbin_0000.bin", std::ios::binary);
      if (!fout.is_open()) return -1;
    }
    fout.write(reinterpret_cast<char *>(&bitstream.buffer[bitstream.size]), binSize);
    if (!fout) {
      return -1;
    }
    bitstream.size += binSize;
    fout.close();

    positionsSize = bitstream.size - positionsSize;
    std::cout << "positions bitstream size " << positionsSize << " B" << std::endl;

    // Decompress geometry with TMC1.
    std::string cmd;
    cmd = "7z e -aoa -bso0 -ooutbin outbin/outbin_0000.bin";
    system(cmd.c_str());
    cmd = "TMC1_geometryDecode -inbin outbin\\outbin -outref refinedVerticesDecoded.ply";
    cmd += " -depth " + std::to_string(depth) + " -level " + std::to_string(level) +
           " -format binary_little_endian";
    system(cmd.c_str());
    cmd = "TMC1_voxelize -inref refinedVerticesDecoded.ply -outvox quantizedPointCloudDecoded.ply";
    cmd += " -depth " + std::to_string(depth) + " -format binary_little_endian";
    system(cmd.c_str());
    pointCloud.read("quantizedPointCloudDecoded.ply");
    pointCloud.addColors();
    pointCloud.removeReflectances();

    if (pointCloud.hasColors()) {
      AttributeDecoder attrDecoder;
      uint64_t colorsSize = bitstream.size;

      attrDecoder.decodeHeader("color", bitstream);
      attrDecoder.buildPredictors(pointCloud);
      attrDecoder.decodeColors(bitstream, pointCloud);

      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B" << std::endl;
      std::cout << std::endl;
    }

    if (pointCloud.hasReflectances()) {
      AttributeDecoder attrDecoder;
      uint64_t reflectancesSize = bitstream.size;

      attrDecoder.decodeHeader("reflectance", bitstream);
      attrDecoder.buildPredictors(pointCloud);
      attrDecoder.decodeReflectances(bitstream, pointCloud);

      reflectancesSize = bitstream.size - reflectancesSize;
      std::cout << "reflectances bitstream size " << reflectancesSize << " B" << std::endl;
    }

    // Optionally write out point cloud with lossy geometry and lossy color, for later comparison.
    pcc::PCCPointSet3 tempPointCloud(pointCloud);
    tempPointCloud.convertYUVToRGB();
    tempPointCloud.write("decodedVoxelsAndDecodedColors.ply");

    inverseQuantization(pointCloud, roundOutputPositions);
    return 0;
  }

 private:
  void decodePositionsHeader(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    uint32_t pointCount = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, pointCount, bitstream.size);
    uint8_t hasColors = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, hasColors, bitstream.size);
    uint8_t hasReflectances = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, hasReflectances, bitstream.size);

    uint8_t u8value;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, u8value, bitstream.size);
    geometryPointsAreUnique = bool(u8value);

    PCCReadFromBuffer<uint8_t>(bitstream.buffer, u8value, bitstream.size);
    neighbourContextsEnabled = bool(u8value);

    PCCReadFromBuffer<uint8_t>(bitstream.buffer, u8value, bitstream.size);
    inferredDirectCodingModeEnabled = bool(u8value);

    pointCloud.addRemoveAttributes(hasColors, hasReflectances);
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
  }

  void decodeTrisoupHeader(PCCBitstream &bitstream, PCCPointSet3 &pointCloud, size_t &depth,
                          size_t &level, double &intToOrigScale, PCCPoint3D &intToOrigTranslation) {
    uint32_t pointCount = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, pointCount, bitstream.size);
    uint8_t hasColors = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, hasColors, bitstream.size);
    uint8_t hasReflectances = 0;
    PCCReadFromBuffer<uint8_t>(bitstream.buffer, hasReflectances, bitstream.size);

    pointCloud.addRemoveAttributes(hasColors, hasReflectances);
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
    PCCReadFromBuffer<size_t>(bitstream.buffer, depth, bitstream.size);
    PCCReadFromBuffer<size_t>(bitstream.buffer, level, bitstream.size);
    PCCReadFromBuffer<double>(bitstream.buffer, intToOrigScale, bitstream.size);
    PCCReadFromBuffer<PCCPoint3D>(bitstream.buffer, intToOrigTranslation, bitstream.size);
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
  // map the @occupancy pattern bits to take into account symmetries in the
  // neighbour configuration @neighPattern.
  //
  uint8_t
  mapGeometryOccupancyInv(uint8_t occupancy, uint8_t neighPattern)
  {
    switch (kOccMapRotateXIdFromPattern[neighPattern]) {
      case 1: occupancy = kOccMapRotateX270[occupancy]; break;
      case 2: occupancy = kOccMapRotateX270Y180[occupancy]; break;
      case 3: occupancy = kOccMapRotateX090Y180[occupancy]; break;
    }

    if (kOccMapRotateYIdFromPattern[neighPattern]) {
      occupancy = kOccMapRotateY090[occupancy];
    }

    bool flag_ud = (neighPattern & 16) && !(neighPattern & 32);
    if (flag_ud) {
      occupancy = kOccMapMirrorXY[occupancy];
    }

    switch (kOccMapRotateZIdFromPatternXY[neighPattern & 15]) {
      case 1: occupancy = kOccMapRotateZ270[occupancy]; break;
      case 2: occupancy = kOccMapRotateZ180[occupancy]; break;
      case 3: occupancy = kOccMapRotateZ090[occupancy]; break;
    }

    return occupancy;
  }

  //-------------------------------------------------------------------------
  // decode node occupancy bits
  //
  uint32_t decodeGeometryOccupancy(
    bool neighbourContextsEnabled,
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    o3dgc::Adaptive_Bit_Model& ctxSingleChild,
    o3dgc::Static_Bit_Model& ctxEquiProb,
    o3dgc::Adaptive_Data_Model (&ctxOccupancy)[10],
    const PCCOctree3Node& node0
  ) {
    if (!neighbourContextsEnabled) {
      return arithmeticDecoder->decode(ctxOccupancy[0]);
    }

    // neighbouring configuration with reduction from 64 to 10
    int neighPattern = node0.neighPattern;
    int neighPattern10 = kNeighPattern64to10[neighPattern];

    // decode occupancy pattern
    uint32_t occupancy;
    if (neighPattern10 == 0) {
      // neighbour empty and only one point => decode index, not pattern
      if (arithmeticDecoder->decode(ctxSingleChild)) {
        uint32_t cnt = arithmeticDecoder->decode(ctxEquiProb);
        cnt |= arithmeticDecoder->decode(ctxEquiProb) << 1;
        cnt |= arithmeticDecoder->decode(ctxEquiProb) << 2;
        occupancy = 1 << cnt;
      }
      else {
        occupancy = arithmeticDecoder->decode(ctxOccupancy[0]);
      }
    }
    else {
      occupancy = arithmeticDecoder->decode(ctxOccupancy[neighPattern10]);
      occupancy = mapGeometryOccupancyInv(occupancy, neighPattern);
    }

    return occupancy;
  }

  //-------------------------------------------------------------------------
  // Decode a position of a point in a given volume.

  PCCVector3<uint32_t>
  decodePointPosition(
    int nodeSizeLog2,
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    o3dgc::Static_Bit_Model& ctxPointPosBlock
  ){
    PCCVector3<uint32_t> delta{};
    for (int i = nodeSizeLog2; i > 0; i--) {
      delta <<= 1;
      delta[0] |= arithmeticDecoder->decode(ctxPointPosBlock);
      delta[1] |= arithmeticDecoder->decode(ctxPointPosBlock);
      delta[2] |= arithmeticDecoder->decode(ctxPointPosBlock);
    }

    return delta;
  }

  //-------------------------------------------------------------------------
  // Direct coding of position of points in node (early tree termination).
  // Decoded points are written to @outputPoints
  // Returns the number of points emitted.

  template<class OutputIt>
  int decodeDirectPosition(
    int nodeSizeLog2,
    const PCCOctree3Node& node,
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    o3dgc::Adaptive_Bit_Model& ctxBlockSkipTh,
    o3dgc::Adaptive_Bit_Model& ctxNumIdcmPointsEq1,
    o3dgc::Static_Bit_Model& ctxPointPosBlock,
    OutputIt outputPoints
  ) {
      bool isDirectMode = arithmeticDecoder->decode(ctxBlockSkipTh);
      if (!isDirectMode) {
        return 0;
      }

      int numPoints = 1;
      if (arithmeticDecoder->decode(ctxNumIdcmPointsEq1))
        numPoints++;

      for (int i = 0; i < numPoints; i++) {
        // convert node-relative position to world position
        PCCVector3<uint32_t> pos = node.pos + decodePointPosition(
          nodeSizeLog2, arithmeticDecoder, ctxPointPosBlock);

        *(outputPoints++) = {double(pos[0]), double(pos[1]), double(pos[2])};
      }

      return numPoints;
  }

  //-------------------------------------------------------------------------

  void decodePositions(PCCBitstream &bitstream, PCCPointSet3 &pointCloud) {
    uint32_t compressedBitstreamSize;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, compressedBitstreamSize, bitstream.size);
    o3dgc::Arithmetic_Codec arithmeticDecoder;
    arithmeticDecoder.set_buffer(uint32_t(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticDecoder.start_decoder();

    o3dgc::Static_Bit_Model ctxEquiProb;
    o3dgc::Adaptive_Bit_Model ctxSingleChild;
    o3dgc::Adaptive_Bit_Model ctxSinglePointPerBlock;
    o3dgc::Adaptive_Bit_Model ctxPointCountPerBlock;
    o3dgc::Adaptive_Bit_Model ctxBlockSkipTh;
    o3dgc::Adaptive_Bit_Model ctxNumIdcmPointsEq1;

    // pattern model using ten 6-neighbour configurations
    o3dgc::Adaptive_Data_Model ctxOccupancy[10];
    for (int i = 0; i < 10; i++) {
        ctxOccupancy[i].set_alphabet(256);
        if (neighbourContextsEnabled)
          ctxOccupancy[i].reset(kInitCtxOccupancy + 256 * i, true);
    }


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
    node00.neighPattern = 0;
    node00.numSiblingsPlus1 = 8;
    fifo.push_back(node00);

    size_t processedPointCount = 0;
    std::vector<uint32_t> values;

    auto fifoCurrLvlEnd = fifo.end();

    // this counter represents fifo.end() - fifoCurrLvlEnd().
    // ie, the number of nodes added to the next level of the tree
    int numNodesNextLvl = 0;

    for (; !fifo.empty(); fifo.pop_front()) {
      if (fifo.begin() == fifoCurrLvlEnd) {
        // transition to the next level
        fifoCurrLvlEnd = fifo.end();
        nodeSizeLog2--;
        numNodesNextLvl = 0;
      }
      PCCOctree3Node& node0 = fifo.front();

      // decode occupancy pattern
      uint8_t occupancy = decodeGeometryOccupancy(
          neighbourContextsEnabled, &arithmeticDecoder, ctxSingleChild,
          ctxEquiProb, ctxOccupancy, node0
      );
      assert(occupancy > 0);

      // population count of occupancy for IDCM
      int numOccupied = popcnt(occupancy);

      // split the current node
      for (int i = 0; i < 8; i++) {
        uint32_t mask = 1 << i;
        if (!(occupancy & mask)) {
          // child is empty: skip
          continue;
        }

        int x = !!(i & 4);
        int y = !!(i & 2);
        int z = !!(i & 1);

        int childSizeLog2 = nodeSizeLog2 - 1;

        // point counts for leaf nodes are coded immediately upon
        // encountering the leaf node.
        if (childSizeLog2 == 0) {
          int numPoints = 1;

          if (!geometryPointsAreUnique) {
            numPoints = decodePositionLeafNumPoints(
              &arithmeticDecoder,
              ctxSinglePointPerBlock, ctxEquiProb, ctxPointCountPerBlock
            );
          }

          const PCCVector3D point(
            node0.pos[0] + (x << childSizeLog2),
            node0.pos[1] + (y << childSizeLog2),
            node0.pos[2] + (z << childSizeLog2)
          );

          for (int i = 0; i < numPoints; ++i)
            pointCloud[processedPointCount++] = point;

          // do not recurse into leaf nodes
          continue;
        }

        // create & enqueue new child.
        fifo.emplace_back();
        auto& child = fifo.back();

        child.pos[0] = node0.pos[0] + (x << childSizeLog2);
        child.pos[1] = node0.pos[1] + (y << childSizeLog2);
        child.pos[2] = node0.pos[2] + (z << childSizeLog2);
        child.numSiblingsPlus1 = numOccupied;

        bool idcmEnabled = inferredDirectCodingModeEnabled;
        if (isDirectModeEligible(idcmEnabled, nodeSizeLog2, node0, child)) {
          int numPoints = decodeDirectPosition(
            childSizeLog2, child, &arithmeticDecoder,
            ctxBlockSkipTh, ctxNumIdcmPointsEq1, ctxEquiProb,
            &pointCloud[processedPointCount]
          );
          processedPointCount += numPoints;

          if (numPoints > 0) {
            // node fully decoded, do not split: discard child
            fifo.pop_back();

            // NB: no further siblings to decode by definition of IDCM
            assert(child.numSiblingsPlus1 == 1);
            break;
          }
        }

        numNodesNextLvl++;
        if (neighbourContextsEnabled) {
          updateGeometryNeighState(
            fifo.end(), numNodesNextLvl, childSizeLog2, child, i,
            node0.neighPattern, occupancy
          );
        }
      }
    }

    arithmeticDecoder.stop_decoder();
    bitstream.size += compressedBitstreamSize;
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
  PCCVector3D minPositions;
  PCCBox3<uint32_t> boundingBox;
  double positionQuantizationScale;

  // When not equal to zero, indicates that a count of points is present in
  // each geometry octree leaf node.
  bool geometryPointsAreUnique;

  // Controls the use of neighbour based contextualisation of octree
  // occupancy during geometry coding.
  bool neighbourContextsEnabled;

  // Controls the use of early termination of the geometry tree
  // by directly coding the position of isolated points.
  bool inferredDirectCodingModeEnabled;
};
}

#endif /* PCCTMC3Decoder_h */
