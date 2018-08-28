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
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include "ringbuf.h"

#include "AttributeEncoder.h"
#include "OctreeNeighMap.h"
#include "PayloadBuffer.h"
#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCPointSetProcessing.h"
#include "PCCTMC3Common.h"
#include "hls.h"
#include "io_hls.h"
#include "pcc_chrono.h"

#include "ArithmeticCodec.h"
#include "tables.h"

namespace pcc {

struct EncoderParams {
  SequenceParameterSet sps;
  GeometryParameterSet gps;

  // NB: information about attributes is split between the SPS and the APS.
  //  => The SPS enumerates the attributes, the APS controls coding params.
  std::vector<AttributeParameterSet> aps;

  // todo(df): this should go away
  std::map<std::string, int> attributeIdxMap;

  // Filename for saving recoloured point cloud.
  std::string postRecolorPath;
};

//============================================================================

class PCCTMC3Encoder3 {
public:
  PCCTMC3Encoder3() { init(); }
  PCCTMC3Encoder3(const PCCTMC3Encoder3&) = default;
  PCCTMC3Encoder3& operator=(const PCCTMC3Encoder3& rhs) = default;
  ~PCCTMC3Encoder3() = default;
  void init()
  {
    minPositions = 0.0;
    boundingBox.min = uint32_t(0);
    boundingBox.max = uint32_t(0);
    pointCloud.clear();
  }

  int compress(
    const PCCPointSet3& inputPointCloud,
    EncoderParams* params,
    std::function<void(const PayloadBuffer&)> outputFn,
    PCCPointSet3* reconstructedCloud = nullptr)
  {
    init();

    // fixup parameter set IDs
    params->sps.sps_seq_parameter_set_id = 0;
    params->gps.gps_seq_parameter_set_id = 0;
    params->gps.gps_geom_parameter_set_id = 0;
    for (int i = 0; i < params->aps.size(); i++) {
      params->aps[i].aps_seq_parameter_set_id = 0;
      params->aps[i].aps_attr_parameter_set_id = i;
    }

    // development level / header
    params->sps.profileCompatibility.profile_compatibility_flags = 0;
    params->sps.level = 0;

    // the encoder writes out the minPositions in the GBH:
    params->gps.geom_box_present_flag = true;

    // fixup attribute parameters
    for (auto it : params->attributeIdxMap) {
      auto& attr_sps = params->sps.attributeSets[it.second];
      attr_sps.attr_instance_id = it.second;

      attr_sps.cicp_colour_primaries_idx = 2;
      attr_sps.cicp_transfer_characteristics_idx = 2;
      attr_sps.cicp_matrix_coefficients_idx = 2;
      attr_sps.cicp_video_full_range_flag = true;

      if (it.first == "color") {
        attr_sps.attr_num_dimensions = 3;
        attr_sps.attributeLabel = KnownAttributeLabel::kColour;
      }

      if (it.first == "reflectance") {
        attr_sps.attr_num_dimensions = 1;
        attr_sps.attributeLabel = KnownAttributeLabel::kReflectance;
      }
    }

    // placeholder to "activate" the parameter sets
    _sps = &params->sps;
    _gps = &params->gps;
    _aps.clear();
    for (const auto& aps : params->aps) {
      _aps.push_back(&aps);
    }

    // write out all parameter sets prior to encoding
    outputFn(write(*_sps));
    outputFn(write(*_gps));
    for (const auto aps : _aps) {
      outputFn(write(*aps));
    }

    // geometry compression consists of the following stages:
    //  - prefilter/quantize geometry (non-normative)
    //  - recolour
    //  - encode geometry
    if (_gps->geom_codec_type == GeometryCodecType::kOctree) {
      quantization(inputPointCloud);
    }
    if (_gps->geom_codec_type == GeometryCodecType::kTriSoup) {
      // todo(?): this ought to be replaced with quantization(...)
      if (!inputPointCloud.write("verticesInWorldCoords.ply", true)) {
        std::cerr << "Failed to write verticesInWorldCoords.ply\n";
        return -1;
      }
      std::string cmd =
        "TMC1_coordinateTransform"
        " -inworld verticesInWorldCoords.ply"
        " -outframe verticesInFrameCoords.ply"
        " -depth "
        + std::to_string(_gps->trisoup_depth) + " -scale "
        + std::to_string(_sps->donotuse_trisoup_int_to_orig_scale)
        + " -format binary_little_endian";
      if (int err = system(cmd.c_str())) {
        std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
        return -1;
      }
    }

    // geometry encoding

    if (1) {
      PayloadBuffer payload(PayloadType::kGeometryBrick);

      pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
      clock_user.start();

      if (_gps->geom_codec_type == GeometryCodecType::kOctree) {
        encodePositions(&payload);
      }
      if (_gps->geom_codec_type == GeometryCodecType::kTriSoup) {
        bool hasColor = inputPointCloud.hasColors();
        bool hasReflectance = inputPointCloud.hasReflectances();
        if (encodeTrisoup(hasColor, hasReflectance, &payload))
          return -1;
        if (recolorTrisoup(inputPointCloud))
          return -1;
      }

      clock_user.stop();

      double bpp =
        double(8 * payload.size()) / inputPointCloud.getPointCount();
      std::cout << "positions bitstream size " << payload.size() << " B ("
                << bpp << " bpp)\n";

      auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
        clock_user.count());
      std::cout << "positions processing time (user): "
                << total_user.count() / 1000.0 << " s" << std::endl;

      outputFn(payload);
    }

    // attributeCoding

    // dump recoloured point cloud
    if (!params->postRecolorPath.empty()) {
      PCCPointSet3 tempPointCloud(pointCloud);
      tempPointCloud.convertYUVToRGB();
      tempPointCloud.write(params->postRecolorPath);
    }

    // for each attribute
    for (const auto& it : params->attributeIdxMap) {
      int attrIdx = it.second;
      const auto& attr_sps = _sps->attributeSets[attrIdx];
      const auto& attr_aps = *_aps[attrIdx];
      const auto& label = attr_sps.attributeLabel;

      PayloadBuffer payload(PayloadType::kAttributeBrick);

      pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
      clock_user.start();

      // todo(df): move elsewhere?
      AttributeBrickHeader abh;
      abh.attr_attr_parameter_set_id = attr_aps.aps_attr_parameter_set_id;
      abh.attr_sps_attr_idx = attrIdx;
      abh.attr_geom_brick_id = 0;
      write(abh, &payload);

      AttributeEncoder attrEncoder;
      attrEncoder.encode(attr_sps, attr_aps, pointCloud, &payload);
      clock_user.stop();

      int coded_size = int(payload.size());
      double bpp = double(8 * coded_size) / inputPointCloud.getPointCount();
      std::cout << label << "s bitstream size " << coded_size << " B (" << bpp
                << " bpp)\n";

      auto time_user = std::chrono::duration_cast<std::chrono::milliseconds>(
        clock_user.count());
      std::cout << label
                << "s processing time (user): " << time_user.count() / 1000.0
                << " s" << std::endl;

      outputFn(payload);
    }

    reconstructedPointCloud(reconstructedCloud);

    return 0;
  }

private:
  void reconstructedPointCloud(PCCPointSet3* reconstructedCloud)
  {
    if (reconstructedCloud == nullptr) {
      return;
    }
    const size_t pointCount = pointCloud.getPointCount();

    reconstructedCloud->addRemoveAttributes(
      pointCloud.hasColors(), pointCloud.hasReflectances());
    reconstructedCloud->resize(pointCount);

    const double minPositionQuantizationScale = 0.0000000001;
    const double invScale =
      fabs(_sps->seq_source_geom_scale_factor) > minPositionQuantizationScale
      ? 1.0 / _sps->seq_source_geom_scale_factor
      : 1.0;
    for (size_t i = 0; i < pointCount; ++i) {
      const auto quantizedPoint = pointCloud[i];
      auto& point = (*reconstructedCloud)[i];
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
    o3dgc::Adaptive_Bit_Model& ctxPointCountPerBlock)
  {
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
  uint8_t mapGeometryOccupancy(uint8_t occupancy, uint8_t neighPattern)
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
  // encode occupancy bits (neighPattern10 == 0 case)

  void encodeOccupancyNeighZ(
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    CtxModelOctreeOccupancy& ctxOccupancy,
    int mappedOccupancy)
  {
    int occupancyB0 = !!(mappedOccupancy & 2);
    int occupancyB1 = !!(mappedOccupancy & 128);
    int occupancyB2 = !!(mappedOccupancy & 32);
    int occupancyB3 = !!(mappedOccupancy & 8);
    int occupancyB4 = !!(mappedOccupancy & 4);
    int occupancyB5 = !!(mappedOccupancy & 16);
    int occupancyB6 = !!(mappedOccupancy & 64);
    int occupancyB7 = !!(mappedOccupancy & 1);

    int numOccupiedAcc = 0;

    arithmeticEncoder->encode(occupancyB0, ctxOccupancy.b0[numOccupiedAcc]);
    numOccupiedAcc += occupancyB0;

    arithmeticEncoder->encode(occupancyB1, ctxOccupancy.b1[numOccupiedAcc]);
    numOccupiedAcc += occupancyB1;

    arithmeticEncoder->encode(occupancyB2, ctxOccupancy.b2[numOccupiedAcc]);
    numOccupiedAcc += occupancyB2;

    arithmeticEncoder->encode(occupancyB3, ctxOccupancy.b3[numOccupiedAcc]);
    numOccupiedAcc += occupancyB3;

    arithmeticEncoder->encode(occupancyB4, ctxOccupancy.b4[numOccupiedAcc]);
    numOccupiedAcc += occupancyB4;

    arithmeticEncoder->encode(occupancyB5, ctxOccupancy.b5[numOccupiedAcc]);
    numOccupiedAcc += occupancyB5;

    // NB: There must be at least two occupied child nodes
    //  -- avoid coding the occupancyB if it is implied.
    if (numOccupiedAcc >= 1)
      arithmeticEncoder->encode(occupancyB6, ctxOccupancy.b6[numOccupiedAcc]);
    numOccupiedAcc += occupancyB6;

    if (numOccupiedAcc >= 2)
      arithmeticEncoder->encode(occupancyB7, ctxOccupancy.b7[numOccupiedAcc]);
  }

  //-------------------------------------------------------------------------
  // encode occupancy bits (neighPattern10 != 0 case)

  void encodeOccupancyNeighNZ(
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    CtxModelOctreeOccupancy& ctxOccupancy,
    int mappedOccupancy,
    int neighPattern10)
  {
    int occupancyB0 = !!(mappedOccupancy & 2);
    int occupancyB1 = !!(mappedOccupancy & 128);
    int occupancyB2 = !!(mappedOccupancy & 32);
    int occupancyB3 = !!(mappedOccupancy & 8);
    int occupancyB4 = !!(mappedOccupancy & 4);
    int occupancyB5 = !!(mappedOccupancy & 16);
    int occupancyB6 = !!(mappedOccupancy & 64);
    int occupancyB7 = !!(mappedOccupancy & 1);

    uint32_t partialOccupancy = 0;
    int idx;

    idx = neighPattern10;
    arithmeticEncoder->encode(occupancyB0, ctxOccupancy.b0[idx]);
    partialOccupancy |= occupancyB0;

    idx = (neighPattern10 << 1) + partialOccupancy;
    arithmeticEncoder->encode(occupancyB1, ctxOccupancy.b1[idx]);
    partialOccupancy |= occupancyB1 << 1;

    idx = (neighPattern10 << 2) + partialOccupancy;
    arithmeticEncoder->encode(occupancyB2, ctxOccupancy.b2[idx]);
    partialOccupancy |= occupancyB2 << 2;

    idx = (neighPattern10 << 3) + partialOccupancy;
    arithmeticEncoder->encode(occupancyB3, ctxOccupancy.b3[idx]);
    partialOccupancy |= occupancyB3 << 3;

    // todo(df): merge constants into lut.
    idx = ((neighPattern10 - 1) << 4) + partialOccupancy;
    idx = kOccMapBit4CtxIdx[idx] - 1 + 5;
    arithmeticEncoder->encode(occupancyB4, ctxOccupancy.b4[idx]);
    partialOccupancy |= occupancyB4 << 4;

    idx = ((neighPattern10 - 1) << 5) + partialOccupancy;
    idx = kOccMapBit5CtxIdx[idx] - 1 + 6;
    arithmeticEncoder->encode(occupancyB5, ctxOccupancy.b5[idx]);
    partialOccupancy |= occupancyB5 << 5;

    int neighPattern7 = kNeighPattern10to7[neighPattern10];
    idx = ((neighPattern7 - 1) << 6) + partialOccupancy;
    idx = kOccMapBit6CtxIdx[idx] - 1 + 7;
    arithmeticEncoder->encode(occupancyB6, ctxOccupancy.b6[idx]);
    partialOccupancy |= occupancyB6 << 6;

    int neighPattern5 = kNeighPattern7to5[neighPattern7];
    idx = ((neighPattern5 - 1) << 7) + partialOccupancy;
    idx = kOccMapBit7CtxIdx[idx] - 1 + 8;
    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    if (partialOccupancy)
      arithmeticEncoder->encode(occupancyB7, ctxOccupancy.b7[idx]);
  }

  //-------------------------------------------------------------------------
  // decode node occupancy bits
  //
  void encodeGeometryOccupancy(
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    o3dgc::Adaptive_Bit_Model& ctxSingleChild,
    o3dgc::Static_Bit_Model& ctxEquiProb,
    CtxModelOctreeOccupancy& ctxOccupancy,
    const PCCOctree3Node& node0,
    int occupancy)
  {
    // code occupancy using the neighbour configuration context
    // with reduction from 64 states to 10.
    int neighPattern = node0.neighPattern;
    int neighPattern10 = kNeighPattern64to10[neighPattern];

    uint32_t mappedOccupancy = mapGeometryOccupancy(occupancy, neighPattern);

    if (neighPattern10 == 0) {
      bool singleChild = !popcntGt1(occupancy);
      arithmeticEncoder->encode(singleChild, ctxSingleChild);

      if (singleChild) {
        // no siblings => encode index = (z,y,x) not 8bit pattern
        arithmeticEncoder->encode(!!(occupancy & 0xaa), ctxEquiProb);  // z
        arithmeticEncoder->encode(!!(occupancy & 0xcc), ctxEquiProb);  // y
        arithmeticEncoder->encode(!!(occupancy & 0xf0), ctxEquiProb);  // x
      } else {
        encodeOccupancyNeighZ(
          arithmeticEncoder, ctxOccupancy, mappedOccupancy);
      }
    } else {
      encodeOccupancyNeighNZ(
        arithmeticEncoder, ctxOccupancy, mappedOccupancy, neighPattern10);
    }
  }

  //-------------------------------------------------------------------------
  // Encode a position of a point in a given volume.

  void encodePointPosition(
    int nodeSizeLog2,
    const PCCVector3<uint32_t>& pos,
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    o3dgc::Static_Bit_Model& ctxPointPosBlock)
  {
    for (int mask = 1 << (nodeSizeLog2 - 1); mask; mask >>= 1) {
      arithmeticEncoder->encode(!!(pos[0] & mask), ctxPointPosBlock);
      arithmeticEncoder->encode(!!(pos[1] & mask), ctxPointPosBlock);
      arithmeticEncoder->encode(!!(pos[2] & mask), ctxPointPosBlock);
    }
  }

  //-------------------------------------------------------------------------
  // Direct coding of position of points in node (early tree termination).

  bool encodeDirectPosition(
    int nodeSizeLog2,
    const PCCOctree3Node& node,
    o3dgc::Arithmetic_Codec* arithmeticEncoder,
    o3dgc::Adaptive_Bit_Model& ctxBlockSkipTh,
    o3dgc::Adaptive_Bit_Model& ctxNumIdcmPointsEq1,
    o3dgc::Static_Bit_Model& ctxPointPosBlock)
  {
    int numPoints = node.end - node.start;
    if (numPoints > MAX_NUM_DM_LEAF_POINTS) {
      arithmeticEncoder->encode(0, ctxBlockSkipTh);
      return false;
    }

    arithmeticEncoder->encode(1, ctxBlockSkipTh);
    arithmeticEncoder->encode(numPoints > 1, ctxNumIdcmPointsEq1);

    for (auto idx = node.start; idx < node.end; idx++) {
      // determine the point position relative to box edge
      encodePointPosition(
        nodeSizeLog2,
        PCCVector3<uint32_t>{int(pointCloud[idx][0]) - node.pos[0],
                             int(pointCloud[idx][1]) - node.pos[1],
                             int(pointCloud[idx][2]) - node.pos[2]},
        arithmeticEncoder, ctxPointPosBlock);
    }

    return true;
  }

  //-------------------------------------------------------------------------

  void encodePositions(PayloadBuffer* buf)
  {
    uint32_t maxBB = std::max(
      {// todo(df): confirm minimum of 1 isn't needed
       1u, boundingBox.max[0], boundingBox.max[1], boundingBox.max[2]});

    // the current node dimension (log2) encompasing maxBB
    int nodeSizeLog2 = ceillog2(maxBB + 1);

    GeometryBrickHeader gbh;
    gbh.geom_geom_parameter_set_id = _gps->gps_geom_parameter_set_id;
    gbh.geomBoxOrigin.x() = int(minPositions.x());
    gbh.geomBoxOrigin.y() = int(minPositions.y());
    gbh.geomBoxOrigin.z() = int(minPositions.z());
    gbh.geom_box_log2_scale = 0;
    gbh.geom_max_node_size_log2 = nodeSizeLog2;
    gbh.geom_num_points = int(pointCloud.getPointCount());
    write(*_gps, gbh, buf);

    // todo(df): remove estimate when arithmetic codec is replaced
    int maxAcBufLen = int(pointCloud.getPointCount()) * 3 * 4 + 1024;
    o3dgc::Arithmetic_Codec arithmeticEncoder(maxAcBufLen, nullptr);
    arithmeticEncoder.start_encoder();

    o3dgc::Static_Bit_Model ctxEquiProb;
    o3dgc::Adaptive_Bit_Model ctxSingleChild;
    o3dgc::Adaptive_Bit_Model ctxSinglePointPerBlock;
    o3dgc::Adaptive_Bit_Model ctxPointCountPerBlock;
    o3dgc::Adaptive_Bit_Model ctxBlockSkipTh;
    o3dgc::Adaptive_Bit_Model ctxNumIdcmPointsEq1;
    CtxModelOctreeOccupancy ctxOccupancy;

    // init main fifo
    //  -- worst case size is the last level containing every input poit
    //     and each point being isolated in the previous level.
    pcc::ringbuf<PCCOctree3Node> fifo(pointCloud.getPointCount() + 1);

    // push the first node
    PCCOctree3Node node00;
    node00.start = uint32_t(0);
    node00.end = uint32_t(pointCloud.getPointCount());
    node00.pos = uint32_t(0);
    node00.neighPattern = 0;
    node00.numSiblingsPlus1 = 8;
    node00.siblingOccupancy = 0;
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

    MortonMap3D occupancyAtlas;
    if (_gps->neighbour_avail_boundary_log2) {
      occupancyAtlas.resize(_gps->neighbour_avail_boundary_log2);
      occupancyAtlas.clear();
    }
    PCCVector3<uint32_t> occupancyAtlasOrigin(0xffffffff);

    for (; !fifo.empty(); fifo.pop_front()) {
      if (fifo.begin() == fifoCurrLvlEnd) {
        // transition to the next level
        fifoCurrLvlEnd = fifo.end();
        numNodesNextLvl = 0;
        occupancyAtlasOrigin = 0xffffffff;
        nodeSizeLog2--;
      }

      PCCOctree3Node& node0 = fifo.front();

      // split the current node into 8 children
      //  - perform an 8-way counting sort of the current node's points
      //  - (later) map to child nodes
      int childSizeLog2 = nodeSizeLog2 - 1;
      std::array<int, 8> childCounts = {};
      countingSort(
        PCCPointSet3::iterator(&pointCloud, node0.start),
        PCCPointSet3::iterator(&pointCloud, node0.end), childCounts,
        [=](const PCCPointSet3::Proxy& proxy) {
          const auto& point = *proxy;
          int bitpos = 1 << childSizeLog2;
          return !!(int(point[2]) & bitpos) | (!!(int(point[1]) & bitpos) << 1)
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

      if (_gps->neighbour_avail_boundary_log2) {
        updateGeometryOccupancyAtlas(
          node0.pos, nodeSizeLog2, fifo, fifoCurrLvlEnd, &occupancyAtlas,
          &occupancyAtlasOrigin);

        node0.neighPattern =
          makeGeometryNeighPattern(node0.pos, nodeSizeLog2, occupancyAtlas);
      }

      // encode child occupancy map
      assert(occupancy > 0);
      encodeGeometryOccupancy(
        &arithmeticEncoder, ctxSingleChild, ctxEquiProb, ctxOccupancy, node0,
        occupancy);

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
          if (_gps->geom_unique_points_flag) {
            assert(childCounts[i] == 1);
            continue;
          }

          processedPointCount += encodePositionLeafNumPoints(
            childCounts[i], &arithmeticEncoder, ctxSinglePointPerBlock,
            ctxEquiProb, ctxPointCountPerBlock);
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
        child.siblingOccupancy = occupancy;

        bool idcmEnabled = _gps->inferred_direct_coding_mode_enabled_flag;
        if (isDirectModeEligible(idcmEnabled, nodeSizeLog2, node0, child)) {
          bool directModeUsed = encodeDirectPosition(
            childSizeLog2, child, &arithmeticEncoder, ctxBlockSkipTh,
            ctxNumIdcmPointsEq1, ctxEquiProb);

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

        // NB: when neighbourAvailBoundaryLog2 is set, an alternative
        //     implementation is used to calculate neighPattern.
        if (!_gps->neighbour_avail_boundary_log2) {
          updateGeometryNeighState(
            _gps->neighbour_context_restriction_flag, fifo.end(),
            numNodesNextLvl, childSizeLog2, child, i, node0.neighPattern,
            occupancy);
        }
      }
    }

    ////
    // The following is to re-order the points according in the decoding
    // order since IDCM causes leaves to be coded earlier than they
    // otherwise would.
    PCCPointSet3 pointCloud2;
    pointCloud2.addRemoveAttributes(
      pointCloud.hasColors(), pointCloud.hasReflectances());
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

    uint32_t acDataLen = arithmeticEncoder.stop_encoder();
    std::copy_n(
      arithmeticEncoder.buffer(), acDataLen, std::back_inserter(*buf));
  }

  //--------------------------------------------------------------------------

  int encodeTrisoup(bool hasColor, bool hasReflectance, PayloadBuffer* buf)
  {
    // todo(?): port TMC1 geometry encoder
    // todo(?): fix TMC1 to output where it is told not requiring an extra
    //          level of guesswork to find the output file.
    pcc::mkdir("outbin");
    std::string cmd =
      "TMC1_geometryEncode"
      " -inframe verticesInFrameCoords.ply"
      " -outbin outbin"
      " -outref refinedVertices.ply"
      " -depth "
      + std::to_string(_gps->trisoup_depth) + " -level "
      + std::to_string(_gps->trisoup_triangle_level)
      + " -format binary_little_endian";
    if (int err = system(cmd.c_str())) {
      std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
      return -1;
    }

    // todo(?): port interpolator
    cmd =
      "TMC1_voxelize"
      " -inref refinedVertices.ply"
      " -outvox voxelizedVertices.ply"
      " -depth "
      + std::to_string(_gps->trisoup_depth) + " -format binary_little_endian";
    if (int err = system(cmd.c_str())) {
      std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
      return -1;
    }

    // Read in quantizedPointCloud file.
    if (!pointCloud.read("voxelizedVertices.ply")) {
      std::cerr << "Failed to open voxelizedVertices.ply\n";
      return -1;
    }
    pointCloud.addRemoveAttributes(hasColor, hasReflectance);

    // encapsulate the TMC1 "bitstream"
    std::ifstream fin("outbin/outbin_0000.bin", std::ios::binary);
    if (!fin.is_open()) {
      std::cerr << "failed to open tmc1 bitstream\n";
      return -1;
    }

    std::copy(
      std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>(),
      std::back_inserter(*buf));
    return 0;
  }

  //--------------------------------------------------------------------------

  void computeMinPositions(const PCCPointSet3& inputPointCloud)
  {
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

  int quantization(const PCCPointSet3& inputPointCloud)
  {
    computeMinPositions(inputPointCloud);

    bool hasColour = inputPointCloud.hasColors();
    bool hasReflectance = inputPointCloud.hasReflectances();

    pointCloud.resize(0);
    pointCloud.addRemoveAttributes(hasColour, hasReflectance);

    const size_t inputPointCount = inputPointCloud.getPointCount();
    if (_gps->geom_unique_points_flag) {  // retain unique quantized points
      // filter points based on quantized positions
      std::set<PCCVector3<int64_t>> retainedPoints;
      for (size_t i = 0; i < inputPointCount; ++i) {
        const PCCVector3D point = (inputPointCloud[i] - minPositions)
          * _sps->seq_source_geom_scale_factor;
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
      const double invScale = 1.0 / _sps->seq_source_geom_scale_factor;
      size_t pointCounter = 0;
      for (const auto& quantizedPoint : retainedPoints) {
        auto& point = pointCloud[pointCounter++];
        for (size_t k = 0; k < 3; ++k) {
          point[k] = double(quantizedPoint[k]) * invScale + minPositions[k];
        }
      }

      if (hasColour) {  // transfer colors
        if (!PCCTransfertColors(inputPointCloud, pointCloud)) {
          std::cout << "Error: can't transfer colors!" << std::endl;
          return -1;
        }
      }

      if (hasReflectance) {  // transfer reflectances
        if (!PCCTransfertReflectances(inputPointCloud, pointCloud)) {
          std::cout << "Error: can't transfer reflectances!" << std::endl;
          return -1;
        }
      }

      // compute quantized coordinates
      pointCounter = 0;
      for (const auto& quantizedPoint : retainedPoints) {
        auto& point = pointCloud[pointCounter++];
        for (size_t k = 0; k < 3; ++k) {
          point[k] = double(quantizedPoint[k]);
        }
      }
    } else {
      // quantized point cloud
      pointCloud.resize(inputPointCount);
      for (size_t i = 0; i < inputPointCount; ++i) {
        const PCCVector3D point = (inputPointCloud[i] - minPositions)
          * _sps->seq_source_geom_scale_factor;
        auto& quantizedPoint = pointCloud[i];
        quantizedPoint[0] = std::round(point[0]);
        quantizedPoint[1] = std::round(point[1]);
        quantizedPoint[2] = std::round(point[2]);
        if (hasColour) {
          pointCloud.setColor(i, inputPointCloud.getColor(i));
        }
        if (hasReflectance) {
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

  //--------------------------------------------------------------------------

  int recolorTrisoup(const PCCPointSet3& inputPointCloud)
  {
    // Transform from integer to original coordinates.
    const size_t pointCount = pointCloud.getPointCount();
    for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
      pointCloud[pointIndex] *= _sps->donotuse_trisoup_int_to_orig_scale;
      // todo(pchou): should add intToOrigTranslation here.
    }

    // Recolour attributes.
    if (inputPointCloud.hasColors()) {
      if (!PCCTransfertColors(inputPointCloud, pointCloud)) {
        std::cout << "Error: can't transfer colors!" << std::endl;
        return -1;
      }
    }
    if (inputPointCloud.hasReflectances()) {
      if (!PCCTransfertReflectances(inputPointCloud, pointCloud)) {
        std::cout << "Error: can't transfer reflectance!" << std::endl;
        return -1;
      }
    }

    // Transform from original to integer coordinates.
    for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
      // todo(pchou): should subtract intToOrigTranslation here.
      // todo(?): don't use division here
      pointCloud[pointIndex] /= _sps->donotuse_trisoup_int_to_orig_scale;
    }

    return 0;
  }

private:
  PCCVector3D minPositions;
  PCCBox3<uint32_t> boundingBox;
  PCCPointSet3 pointCloud;

  // The active parameter sets
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;
  std::vector<const AttributeParameterSet*> _aps;
};
}  // namespace pcc

#endif /* PCCTMC3Encoder_h */
