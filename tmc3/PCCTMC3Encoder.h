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
#include <deque>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include "ringbuf.h"

#include "AttributeEncoder.h"
#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCPointSetProcessing.h"
#include "PCCTMC3Common.h"

#include "ArithmeticCodec.h"
#include "tables.h"

namespace pcc {
struct PCCTMC3Encoder3Parameters {
  double positionQuantizationScale;
  bool mergeDuplicatedPoints;

  // Controls the use of neighbour based contextualisation of octree
  // occupancy during geometry coding.
  bool neighbourContextsEnabled;

  // Controls the use of early termination of the geometry tree
  // by directly coding the position of isolated points.
  bool inferredDirectCodingModeEnabled;

  struct TriSoup {
    // depth of voxels (reconstructed points) in trisoup geometry
    int depth;

    // level of triangles (reconstructed surface) in trisoup geometry
    int level;

    // orig_coords = integer_coords * intToOrigScale + intToOrigTranslation
    double intToOrigScale;

    // orig_coords = integer_coords * intToOrigScale + intToOrigTranslation
    // todo(df): convert to PCCVector3
    std::vector<double> intToOrigTranslation;
  } triSoup;

  std::map<std::string, PCCAttributeEncodeParamaters> attributeEncodeParameters;

  //--------------------------------------------------------------------------
  // Retrieve the a set of parameters, or nullptr if they either do not exist
  // or guard is false.

  const PCCAttributeEncodeParamaters*
  getAttrParams(bool guard, const std::string& what) const
  {
    if (!guard)
      return nullptr;
    auto it = attributeEncodeParameters.find(what);
    if (it == attributeEncodeParameters.end())
      return nullptr;
    return &it->second;
  }
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

    auto paramsColor =
      params.getAttrParams(inputPointCloud.hasColors(), "color");

    if (!paramsColor) {
      pointCloud.removeColors();
    }

    auto paramsReflectance =
      params.getAttrParams(inputPointCloud.hasReflectances(), "reflectance");

    if (!paramsReflectance) {
      pointCloud.removeReflectances();
    }

    PCCWriteToBuffer<uint32_t>(PCCTMC3MagicNumber, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint32_t>(PCCTMC3FormatVersion, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasColors()), bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasReflectances()), bitstream.buffer,
                              bitstream.size);

    if (paramsColor) {
      AttributeEncoder attrEncoder;
      uint64_t colorsSize = bitstream.size;

      attrEncoder.encodeHeader(*paramsColor, "color", bitstream);
      attrEncoder.buildPredictors(*paramsColor, pointCloud);
      attrEncoder.encodeColors(*paramsColor, pointCloud, bitstream);

      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B ("
                << (8.0 * colorsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;
    }

    if (paramsReflectance) {
      AttributeEncoder attrEncoder;
      uint64_t reflectancesSize = bitstream.size;

      attrEncoder.encodeHeader(*paramsReflectance, "reflectance", bitstream);
      attrEncoder.buildPredictors(*paramsReflectance, pointCloud);
      attrEncoder.encodeReflectances(*paramsReflectance, pointCloud, bitstream);

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
    encodePositionsHeader(params, bitstream);
    encodePositions(params, bitstream);

    positionsSize = bitstream.size - positionsSize;
    std::cout << "positions bitstream size " << positionsSize << " B ("
              << (8.0 * positionsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;

    auto paramsColor =
      params.getAttrParams(pointCloud.hasColors(), "color");

    if (paramsColor) {
      AttributeEncoder attrEncoder;
      uint64_t colorsSize = bitstream.size;

      attrEncoder.encodeHeader(*paramsColor, "color", bitstream);
      attrEncoder.buildPredictors(*paramsColor, pointCloud);
      attrEncoder.encodeColors(*paramsColor, pointCloud, bitstream);

      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B ("
                << (8.0 * colorsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;
    }

    auto paramsReflectance =
      params.getAttrParams(pointCloud.hasReflectances(), "reflectance");

    if (paramsReflectance) {
      AttributeEncoder attrEncoder;
      uint64_t reflectancesSize = bitstream.size;

      attrEncoder.encodeHeader(*paramsReflectance, "reflectance", bitstream);
      attrEncoder.buildPredictors(*paramsReflectance, pointCloud);
      attrEncoder.encodeReflectances(*paramsReflectance, pointCloud, bitstream);

      reflectancesSize = bitstream.size - reflectancesSize;
      std::cout << "reflectances bitstream size " << reflectancesSize << " B ("
                << (8.0 * reflectancesSize) / inputPointCloud.getPointCount() << " bpp)"
                << std::endl;
    }

    reconstructedPointCloud(params, reconstructedCloud);
    return 0;
  }
  int compressWithTrisoupGeometry(const PCCPointSet3 &inputPointCloud,
                                  const PCCTMC3Encoder3Parameters &params, PCCBitstream &bitstream,
                                  PCCPointSet3 *reconstructedCloud = nullptr) {
    init();
    PCCWriteToBuffer<uint32_t>(PCCTMC3MagicNumber, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint32_t>(PCCTMC3FormatVersion, bitstream.buffer, bitstream.size);

    // Write out inputPointCloud file, just so that we can get it into the TMC1 pipeline.
    pcc::PCCPointSet3 tempPointCloud(inputPointCloud);
    tempPointCloud.write("verticesInWorldCoords.ply", true);

    // Compress geometry with TMC1.
    pcc::mkdir("outbin");
    std::string cmd;
    cmd =
        "TMC1_coordinateTransform -inworld verticesInWorldCoords.ply -outframe "
        "verticesInFrameCoords.ply";
    cmd += " -depth " + std::to_string(params.triSoup.depth) + " -scale " +
           std::to_string(params.triSoup.intToOrigScale) + " -format binary_little_endian";
    system(cmd.c_str());
    cmd =
        "TMC1_geometryEncode -inframe verticesInFrameCoords.ply -outbin outbin -outref "
        "refinedVertices.ply";
    cmd += " -depth " + std::to_string(params.triSoup.depth) + " -level " + std::to_string(params.triSoup.level) +
           " -format binary_little_endian";
    system(cmd.c_str());
    cmd = "TMC1_voxelize -inref refinedVertices.ply -outvox voxelizedVertices.ply";
    cmd += " -depth " + std::to_string(params.triSoup.depth) + " -format binary_little_endian";
    system(cmd.c_str());

    // Read in quantizedPointCloud file, so that we can recolor it.
    pointCloud.read("voxelizedVertices.ply");

    // Transform from integer to original coordinates.
    const size_t pointCount = pointCloud.getPointCount();
    for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
      pointCloud[pointIndex] *= params.triSoup.intToOrigScale;
      // Should add intToOrigTranslation here.
    }

    auto paramsColor =
      params.getAttrParams(inputPointCloud.hasColors(), "color");

    auto paramsReflectance =
      params.getAttrParams(inputPointCloud.hasReflectances(), "reflectance");

    if (paramsColor) {
      pointCloud.addColors();
      int32_t searchRange = paramsColor->searchRange;
      if (!PCCTransfertColors(inputPointCloud, searchRange, pointCloud)) {
        std::cout << "Error: can't transfer colors!" << std::endl;
        return -1;
      }
    }

    // Transform from original to integer coordinates.
    for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
      // Should subtract intToOrigTranslation here.
      pointCloud[pointIndex] /= params.triSoup.intToOrigScale;
    }

    // Optionally write out point cloud with lossy geometry and lossless color, for later
    // comparison.
    tempPointCloud = pointCloud;  // deep copy
    tempPointCloud.convertYUVToRGB();
    tempPointCloud.write("decodedVoxelsAndTransferredColors.ply");

    // Encode positions.
    uint64_t positionsSize = bitstream.size;
    encodeTrisoupHeader(params, bitstream);

    // Read in TMC1 geometry bitstream and put it into the converged TMC13 bitstream.
    std::ifstream fin("outbin/outbin_0000.bin", std::ios::binary);
    if (!fin.is_open()) {
      return -1;
    }
    fin.seekg(0, std::ios::end);
    uint64_t binSize = fin.tellg();
    fin.seekg(0, std::ios::beg);
    std::unique_ptr<uint8_t[]> tempbuf(new uint8_t[binSize]);
    PCCWriteToBuffer<uint32_t>(uint32_t(binSize), bitstream.buffer, bitstream.size);
    fin.read(reinterpret_cast<char *>(&bitstream.buffer[bitstream.size]), binSize);
    if (!fin) {
      return -1;
    }
    bitstream.size += binSize;
    fin.close();

    positionsSize = bitstream.size - positionsSize;
    std::cout << "positions bitstream size " << positionsSize << " B ("
              << (8.0 * positionsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;

    if (paramsColor) {
      AttributeEncoder attrEncoder;
      uint64_t colorsSize = bitstream.size;

      attrEncoder.buildPredictors(*paramsColor, pointCloud);
      attrEncoder.encodeHeader(*paramsColor, "color", bitstream);
      attrEncoder.encodeColors(*paramsColor, pointCloud, bitstream);

      colorsSize = bitstream.size - colorsSize;
      std::cout << "colors bitstream size " << colorsSize << " B ("
                << (8.0 * colorsSize) / inputPointCloud.getPointCount() << " bpp)" << std::endl;
    }

    if (paramsReflectance) {
      AttributeEncoder attrEncoder;
      uint64_t reflectancesSize = bitstream.size;

      attrEncoder.buildPredictors(*paramsReflectance, pointCloud);
      attrEncoder.encodeHeader(*paramsReflectance, "reflectance", bitstream);
      attrEncoder.encodeReflectances(*paramsReflectance, pointCloud, bitstream);

      reflectancesSize = bitstream.size - reflectancesSize;
      std::cout << "reflectances bitstream size " << reflectancesSize << " B ("
                << (8.0 * reflectancesSize) / inputPointCloud.getPointCount() << " bpp)"
                << std::endl;
    }

    reconstructedPointCloud(params, reconstructedCloud);
    return 0;
  }

 private:
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
  // map the @occupancy pattern bits to take into account symmetries in the
  // neighbour configuration @neighPattern.
  //
  uint8_t
  mapGeometryOccupancy(uint8_t occupancy, uint8_t neighPattern)
  {
    switch (kOccMapRotateZIdFromPatternXY[neighPattern & 15]) {
      case 1: occupancy = kOccMapRotateZ090[occupancy]; break;
      case 2: occupancy = kOccMapRotateZ180[occupancy]; break;
      case 3: occupancy = kOccMapRotateZ270[occupancy]; break;
    }

    bool flag_ud = (neighPattern & 16) && !(neighPattern & 32);
    if (flag_ud) {
      occupancy = kOccMapMirrorXY[occupancy];
    }

    if (kOccMapRotateYIdFromPattern[neighPattern]) {
      occupancy = kOccMapRotateY270[occupancy];
    }

    switch (kOccMapRotateXIdFromPattern[neighPattern]) {
      case 1: occupancy = kOccMapRotateX090[occupancy]; break;
      case 2: occupancy = kOccMapRotateX270Y180[occupancy]; break;
      case 3: occupancy = kOccMapRotateX090Y180[occupancy]; break;
    }

    return occupancy;
  }

  //-------------------------------------------------------------------------
  // enecode a node's occupancy bits
  //
  void encodeGeometryOccupancy(
    bool neighbourContextsEnabled,
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    o3dgc::Adaptive_Bit_Model& ctxSingleChild,
    o3dgc::Static_Bit_Model& ctxEquiProb,
    o3dgc::Adaptive_Data_Model (&ctxOccupancy)[10],
    const PCCOctree3Node& node0,
    int occupancy
  ) {
    if (!neighbourContextsEnabled) {
      arithmeticEncoder->encode(occupancy, ctxOccupancy[0]);
      return;
    }

    // code occupancy using the neighbour configuration context
    // with reduction from 64 states to 10.
    int neighPattern = node0.neighPattern;
    int neighPattern10 = kNeighPattern64to10[neighPattern];

    if (neighPattern10 == 0) {
      bool singleChild = !popcntGt1(occupancy);
      arithmeticEncoder->encode(singleChild, ctxSingleChild);

      if (singleChild) {
        // no siblings => encode index = (z,y,x) not 8bit pattern
        arithmeticEncoder->encode(!!(occupancy & 0xaa), ctxEquiProb); // z
        arithmeticEncoder->encode(!!(occupancy & 0xcc), ctxEquiProb); // y
        arithmeticEncoder->encode(!!(occupancy & 0xf0), ctxEquiProb); // x
      } else {
        arithmeticEncoder->encode(occupancy, ctxOccupancy[0]);
      }
    }
    else {
      uint32_t mappedOccupancy = mapGeometryOccupancy(occupancy, neighPattern);
      arithmeticEncoder->encode(mappedOccupancy, ctxOccupancy[neighPattern10]);
    }
  }

  //-------------------------------------------------------------------------
  // Encode a position of a point in a given volume.

  void
  encodePointPosition(
    int nodeSizeLog2,
    const PCCVector3<uint32_t>& pos,
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    o3dgc::Static_Bit_Model& ctxPointPosBlock
  ){
    for (int mask = 1 << (nodeSizeLog2 - 1); mask; mask >>= 1) {
      arithmeticEncoder->encode(!!(pos[0] & mask), ctxPointPosBlock);
      arithmeticEncoder->encode(!!(pos[1] & mask), ctxPointPosBlock);
      arithmeticEncoder->encode(!!(pos[2] & mask), ctxPointPosBlock);
    }
  }

  //-------------------------------------------------------------------------
  // Direct coding of position of points in node (early tree termination).

  bool
  encodeDirectPosition(
    int nodeSizeLog2,
    const PCCOctree3Node& node,
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    o3dgc::Adaptive_Bit_Model& ctxBlockSkipTh,
    o3dgc::Adaptive_Bit_Model& ctxNumIdcmPointsEq1,
    o3dgc::Static_Bit_Model& ctxPointPosBlock
  ) {
      int numPoints = node.end - node.start;
      if (numPoints > MAX_NUM_DM_LEAF_POINTS) {
        arithmeticEncoder->encode(0, ctxBlockSkipTh);
        return false;
      }

      arithmeticEncoder->encode(1, ctxBlockSkipTh);
      arithmeticEncoder->encode(numPoints > 1, ctxNumIdcmPointsEq1);

      for (auto idx = node.start; idx < node.end; idx++) {
        // determine the point position relative to box edge
        encodePointPosition(nodeSizeLog2, PCCVector3<uint32_t>{
            int(pointCloud[idx][0]) - node.pos[0],
            int(pointCloud[idx][1]) - node.pos[1],
            int(pointCloud[idx][2]) - node.pos[2]},
          arithmeticEncoder, ctxPointPosBlock);
      }

      return true;
  }

  //-------------------------------------------------------------------------

  void encodePositions(
    const PCCTMC3Encoder3Parameters& params,
    PCCBitstream& bitstream
  ) {
    uint64_t startSize = bitstream.size;
    bitstream.size += 4;  // placehoder for bitstream size
    o3dgc::Arithmetic_Codec arithmeticEncoder;
    arithmeticEncoder.set_buffer(uint32_t(bitstream.capacity - bitstream.size),
                                 bitstream.buffer + bitstream.size);
    arithmeticEncoder.start_encoder();

    o3dgc::Static_Bit_Model ctxEquiProb;
    o3dgc::Adaptive_Bit_Model ctxSingleChild;
    o3dgc::Adaptive_Bit_Model ctxSinglePointPerBlock;
    o3dgc::Adaptive_Bit_Model ctxPointCountPerBlock;
    o3dgc::Adaptive_Bit_Model ctxBlockSkipTh;
    o3dgc::Adaptive_Bit_Model ctxNumIdcmPointsEq1;

    // occupancy map model using ten 6-neighbour configurations
    o3dgc::Adaptive_Data_Model ctxOccupancy[10];
    for (int i = 0; i < 10; i++) {
        ctxOccupancy[i].set_alphabet(256);
        if (params.neighbourContextsEnabled)
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

    // map of pointCloud idx to DM idx, used to reorder the points
    // after coding.
    std::vector<int> pointIdxToDmIdx(int(pointCloud.getPointCount()), -1);
    int nextDmIdx = 0;

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
        numNodesNextLvl = 0;
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

      // generate the bitmap of child occupancy and count
      // the number of occupied children in node0.
      int occupancy = 0;
      int numSiblings = 0;
      for (int i = 0; i < 8; i++) {
        if (childCounts[i]) {
          occupancy |= 1 << i;
          numSiblings++;
        }
      }

      // encode child occupancy map
      assert(occupancy > 0);
      encodeGeometryOccupancy(
        params.neighbourContextsEnabled, &arithmeticEncoder, ctxSingleChild,
        ctxEquiProb, ctxOccupancy, node0, occupancy
      );

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

      // nodeSizeLog2 > 1: for each child:
      //  - determine elegibility for IDCM
      //  - directly code point positions if IDCM allowed and selected
      //  - otherwise, insert split children into fifo while updating neighbour state
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
        child.numSiblingsPlus1 = numSiblings;

        bool idcmEnabled = params.inferredDirectCodingModeEnabled;
        if (isDirectModeEligible(idcmEnabled, nodeSizeLog2, node0, child)) {
          bool directModeUsed = encodeDirectPosition(
            childSizeLog2, child, &arithmeticEncoder,
            ctxBlockSkipTh, ctxNumIdcmPointsEq1, ctxEquiProb
          );

          if (directModeUsed) {
            // point reordering to match decoder's order
            for (auto idx = child.start; idx < child.end; idx++)
              pointIdxToDmIdx[idx] = nextDmIdx++;

            // NB: by definition, this is the only child node present
            assert(child.numSiblingsPlus1 == 1);

            // remove leaf node from fifo: it has been consumed and will
            // not be further split.
            fifo.pop_back();
            break;
          }
        }

        numNodesNextLvl++;
        if (params.neighbourContextsEnabled) {
          updateGeometryNeighState(
            fifo.end(), numNodesNextLvl, childSizeLog2, child, i,
            node0.neighPattern, occupancy
          );
        }
      }
    }

    ////
    // The following is to re-order the points according in the decoding
    // order since IDCM causes leaves to be coded earlier than they
    // otherwise would.
    PCCPointSet3 pointCloud2;
    if (pointCloud.hasColors())
      pointCloud2.addColors();
    if (pointCloud.hasReflectances())
      pointCloud2.addReflectances();
    pointCloud2.resize(pointCloud.getPointCount());

    // copy points with DM points first, the rest second
    int outIdx = nextDmIdx;
    for (int i = 0; i < pointIdxToDmIdx.size(); i++) {
      int dstIdx = pointIdxToDmIdx[i];
      if (dstIdx == -1) {
        dstIdx = outIdx++;
      }

      pointCloud2[dstIdx] = pointCloud[i];
      if (pointCloud.hasColors())
        pointCloud2.setColor(dstIdx, pointCloud.getColor(i));
      if (pointCloud.hasReflectances())
        pointCloud2.setReflectance(dstIdx, pointCloud.getReflectance(i));
    }

    swap(pointCloud, pointCloud2);

    const uint32_t compressedBitstreamSize = arithmeticEncoder.stop_encoder();
    bitstream.size += compressedBitstreamSize;
    PCCWriteToBuffer<uint32_t>(compressedBitstreamSize, bitstream.buffer, startSize);
  }

  void encodePositionsHeader(const PCCTMC3Encoder3Parameters &params,
                            PCCBitstream &bitstream) const {
    PCCWriteToBuffer<uint32_t>(uint32_t(pointCloud.getPointCount()), bitstream.buffer,
                               bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasColors()), bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasReflectances()), bitstream.buffer,
                              bitstream.size);

    // todo(df): syntax element name: geometryPointsAreUnique?
    PCCWriteToBuffer<uint8_t>(
      uint8_t(params.mergeDuplicatedPoints), bitstream.buffer, bitstream.size);

    PCCWriteToBuffer<uint8_t>(
      uint8_t(params.neighbourContextsEnabled), bitstream.buffer, bitstream.size);

    PCCWriteToBuffer<uint8_t>(
      uint8_t(params.inferredDirectCodingModeEnabled), bitstream.buffer, bitstream.size);

    for (int k = 0; k < 3; ++k) {
      PCCWriteToBuffer<double>(minPositions[k], bitstream.buffer, bitstream.size);
    }
    for (int k = 0; k < 3; ++k) {
      PCCWriteToBuffer<uint32_t>(boundingBox.max[k], bitstream.buffer, bitstream.size);
    }
    PCCWriteToBuffer<double>(params.positionQuantizationScale, bitstream.buffer, bitstream.size);
  }

  void encodeTrisoupHeader(const PCCTMC3Encoder3Parameters &params, PCCBitstream &bitstream) const {
    PCCWriteToBuffer<uint32_t>(uint32_t(pointCloud.getPointCount()), bitstream.buffer,
                               bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasColors()), bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<uint8_t>(uint8_t(pointCloud.hasReflectances()), bitstream.buffer,
                              bitstream.size);
    for (int k = 0; k < 3; ++k) {
      PCCWriteToBuffer<double>(minPositions[k], bitstream.buffer, bitstream.size);
    }
    for (int k = 0; k < 3; ++k) {
      PCCWriteToBuffer<uint32_t>(boundingBox.max[k], bitstream.buffer, bitstream.size);
    }
    PCCWriteToBuffer<double>(params.positionQuantizationScale, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<size_t>(params.triSoup.depth, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<size_t>(params.triSoup.level, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<double>(params.triSoup.intToOrigScale, bitstream.buffer, bitstream.size);
    PCCWriteToBuffer<PCCPoint3D>(
      PCCPoint3D(
        params.triSoup.intToOrigTranslation[0],
        params.triSoup.intToOrigTranslation[1],
        params.triSoup.intToOrigTranslation[2]
      ), bitstream.buffer, bitstream.size);
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

    auto paramsColor =
      params.getAttrParams(inputPointCloud.hasColors(), "color");

    if (paramsColor) {
      pointCloud.addColors();
    } else {
      pointCloud.removeColors();
    }

    auto paramsReflectance =
      params.getAttrParams(inputPointCloud.hasReflectances(), "reflectance");

    if (paramsReflectance) {
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

      if (paramsColor) {  // transfer colors
        int32_t searchRange = paramsColor->searchRange;
        if (!PCCTransfertColors(inputPointCloud, searchRange, pointCloud)) {
          std::cout << "Error: can't transfer colors!" << std::endl;
          return -1;
        }
      }

      if (paramsReflectance) {  // transfer reflectances
        int32_t searchRange = paramsReflectance->searchRange;
        if (!PCCTransfertReflectances(inputPointCloud, searchRange, pointCloud)) {
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
        if (paramsColor) {
          pointCloud.setColor(i, inputPointCloud.getColor(i));
        }
        if (paramsReflectance) {
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
  PCCVector3D minPositions;
  PCCBox3<uint32_t> boundingBox;
  PCCPointSet3 pointCloud;
};
}

#endif /* PCCTMC3Encoder_h */
