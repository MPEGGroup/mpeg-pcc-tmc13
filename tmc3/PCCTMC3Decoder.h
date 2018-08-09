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
#include <functional>
#include <queue>
#include <string>
#include "ringbuf.h"

#include "AttributeDecoder.h"
#include "OctreeNeighMap.h"
#include "PayloadBuffer.h"
#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCTMC3Common.h"
#include "hls.h"
#include "io_hls.h"
#include "io_tlv.h"
#include "pcc_chrono.h"
#include "osspecific.h"

#include "ArithmeticCodec.h"
#include "tables.h"

namespace pcc {

struct DecoderParams {
  // Do not delete this structure -- it is for passing options to the decoder.
};

class PCCTMC3Decoder3 {
public:
  PCCTMC3Decoder3() { init(); }
  PCCTMC3Decoder3(const PCCTMC3Decoder3&) = default;
  PCCTMC3Decoder3& operator=(const PCCTMC3Decoder3& rhs) = default;
  ~PCCTMC3Decoder3() = default;

  void init()
  {
    _sps = nullptr;
    _gps = nullptr;
    _spss.clear();
    _gpss.clear();
    _apss.clear();
  }

  int decompress(
    const DecoderParams& params,
    const PayloadBuffer* buf,
    std::function<void(const PCCPointSet3&)> onOutputCloud)
  {
    if (!buf) {
      // flush decoder, output pending cloud if any
      onOutputCloud(_currentPointCloud);
      return 0;
    }

    switch (buf->type) {
    case PayloadType::kSequenceParameterSet:
      storeSps(parseSps(*buf));
      return 0;

    case PayloadType::kGeometryParameterSet:
      storeGps(parseGps(*buf));
      return 0;

    case PayloadType::kAttributeParameterSet:
      storeAps(parseAps(*buf));
      return 0;

    case PayloadType::kGeometryBrick:
      // todo(df): call onOutputCloud when starting a new cloud
      return decodeGeometryBrick(*buf);

    case PayloadType::kAttributeBrick: decodeAttributeBrick(*buf); return 0;
    }

    // todo(df): error, unhandled payload type
    return 1;
  }

  //--------------------------------------------------------------------------

  void storeSps(SequenceParameterSet&& sps)
  {
    // todo(df): handle replacement semantics
    _spss.emplace(std::make_pair(sps.sps_seq_parameter_set_id, sps));
  }

  //--------------------------------------------------------------------------

  void storeGps(GeometryParameterSet&& gps)
  {
    // todo(df): handle replacement semantics
    _gpss.emplace(std::make_pair(gps.gps_geom_parameter_set_id, gps));
  }

  //--------------------------------------------------------------------------

  void storeAps(AttributeParameterSet&& aps)
  {
    // todo(df): handle replacement semantics
    _apss.emplace(std::make_pair(aps.aps_attr_parameter_set_id, aps));
  }

  //==========================================================================

private:
  //--------------------------------------------------------------------------
  // Initialise the point cloud storage and decode a single geometry brick.
  //
  // NB: there is an implicit assumption that there is only a single geometry
  //     brick per point cloud.
  // todo(??): fixme

  int decodeGeometryBrick(const PayloadBuffer& buf)
  {
    assert(buf.type == PayloadType::kGeometryBrick);
    std::cout << "positions bitstream size " << buf.size() << " B\n";

    // HACK: assume activation of the first SPS and GPS
    // todo(df): parse brick header here for propper sps & gps activation
    //  -- this is currently inconsistent between trisoup and octree
    assert(!_spss.empty());
    assert(!_gpss.empty());
    _sps = &_spss.cbegin()->second;
    _gps = &_gpss.cbegin()->second;

    // todo(df): replace with attribute mapping
    bool hasColour = std::any_of(
      _sps->attributeSets.begin(), _sps->attributeSets.end(),
      [](const AttributeDescription& desc) {
        return desc.attributeLabel == KnownAttributeLabel::kColour;
      });

    bool hasReflectance = std::any_of(
      _sps->attributeSets.begin(), _sps->attributeSets.end(),
      [](const AttributeDescription& desc) {
        return desc.attributeLabel == KnownAttributeLabel::kReflectance;
      });

    _currentPointCloud.clear();
    _currentPointCloud.addRemoveAttributes(hasColour, hasReflectance);

    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
    clock_user.start();

    if (_gps->geom_codec_type == GeometryCodecType::kOctree) {
      decodePositions(buf, _currentPointCloud);
    }
    if (_gps->geom_codec_type == GeometryCodecType::kTriSoup) {
      if (decodeTrisoup(buf, _currentPointCloud))
        return -1;
    }
    clock_user.stop();

    auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << "positions processing time (user): "
              << total_user.count() / 1000.0 << " s\n";
    std::cout << std::endl;

    return 0;
  }

  //--------------------------------------------------------------------------
  // Initialise the point cloud storage and decode a single geometry brick.
  //
  // NB: there is an implicit assumption that there is only a single geometry
  //     brick per point cloud.

  void decodeAttributeBrick(const PayloadBuffer& buf)
  {
    assert(buf.type == PayloadType::kAttributeBrick);
    // todo(df): replace assertions with error handling
    assert(_sps);
    assert(_gps);

    // todo(df): validate that sps activation is not changed via the APS
    AttributeBrickHeader abh = parseAbh(buf, nullptr);
    const auto it_attr_aps = _apss.find(abh.attr_attr_parameter_set_id);

    assert(it_attr_aps != _apss.cend());
    const auto& attr_aps = it_attr_aps->second;

    assert(abh.attr_sps_attr_idx < _sps->attributeSets.size());
    const auto& attr_sps = _sps->attributeSets[abh.attr_sps_attr_idx];
    const auto& label = attr_sps.attributeLabel;

    AttributeDecoder attrDecoder;
    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;

    clock_user.start();
    attrDecoder.decode(attr_sps, attr_aps, buf, _currentPointCloud);
    clock_user.stop();

    std::cout << label << "s bitstream size " << buf.size() << " B\n";

    auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << label
              << "s processing time (user): " << total_user.count() / 1000.0
              << " s\n";
    std::cout << std::endl;
  }

  //--------------------------------------------------------------------------

  int decodeTrisoup(const PayloadBuffer& buf, PCCPointSet3& pointCloud)
  {
    assert(buf.type == PayloadType::kGeometryBrick);

    // Write out TMC1 geometry bitstream from the converged TMC13 bitstream.
    std::ofstream fout("outbin/outbin_0000.bin", std::ios::binary);
    if (!fout.is_open()) {
      // maybe because directory does not exist; try again...
      pcc::mkdir("outbin");
      fout.open("outbin/outbin_0000.bin", std::ios::binary);
      if (!fout.is_open())
        return -1;
    }
    fout.write(buf.data(), buf.size());
    if (!fout) {
      return -1;
    }
    fout.close();

    // Decompress geometry with TMC1.
    std::string cmd;
    cmd = "7z e -aoa -bso0 -ooutbin outbin/outbin_0000.bin";
    if (int err = system(cmd.c_str())) {
      std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
      return -1;
    }

    cmd =
      "TMC1_geometryDecode"
      " -inbin outbin/outbin"
      " -outref refinedVerticesDecoded.ply"
      " -depth "
      + std::to_string(_gps->trisoup_depth) + " -level "
      + std::to_string(_gps->trisoup_triangle_level)
      + " -format binary_little_endian";
    if (int err = system(cmd.c_str())) {
      std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
      return -1;
    }

    cmd =
      "TMC1_voxelize"
      " -inref refinedVerticesDecoded.ply"
      " -outvox quantizedPointCloudDecoded.ply"
      " -depth "
      + std::to_string(_gps->trisoup_depth) + " -format binary_little_endian";
    if (int err = system(cmd.c_str())) {
      std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
      return -1;
    }

    bool hasColors = pointCloud.hasColors();
    bool hasReflectances = pointCloud.hasReflectances();
    if (!pointCloud.read("quantizedPointCloudDecoded.ply")) {
      std::cerr << "Failed to read quantizedPointCloudDecoded.ply\n";
      return -1;
    }
    pointCloud.addRemoveAttributes(hasColors, hasReflectances);
    return 0;
  }

  //-------------------------------------------------------------------------
  // Decode the number of points in a leaf node of the octree.

  int decodePositionLeafNumPoints(
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    o3dgc::Adaptive_Bit_Model& ctxSinglePointPerBlock,
    o3dgc::Static_Bit_Model& ctxEquiProb,
    o3dgc::Adaptive_Bit_Model& ctxPointCountPerBlock)
  {
    const bool isSinglePoint =
      arithmeticDecoder->decode(ctxSinglePointPerBlock) != 0;

    int count = 1;
    if (!isSinglePoint) {
      count = 1
        + arithmeticDecoder->ExpGolombDecode(
            0, ctxEquiProb, ctxPointCountPerBlock);
    }

    return count;
  }

  //-------------------------------------------------------------------------
  // map the @occupancy pattern bits to take into account symmetries in the
  // neighbour configuration @neighPattern.
  //
  uint8_t mapGeometryOccupancyInv(uint8_t occupancy, uint8_t neighPattern)
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

  //---------------------------------------------------------------------------
  // decode occupancy bits (neighPattern10 == 0 case)

  int decodeOccupancyNeighZ(
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    CtxModelOctreeOccupancy& ctxOccupancy)
  {
    int numOccupiedAcc = 0;
    int bit;
    int occupancy = 0;

    bit = arithmeticDecoder->decode(ctxOccupancy.b0[numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << 1;

    bit = arithmeticDecoder->decode(ctxOccupancy.b1[numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << 7;

    bit = arithmeticDecoder->decode(ctxOccupancy.b2[numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << 5;

    bit = arithmeticDecoder->decode(ctxOccupancy.b3[numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << 3;

    bit = arithmeticDecoder->decode(ctxOccupancy.b4[numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << 2;

    bit = arithmeticDecoder->decode(ctxOccupancy.b5[numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << 4;

    // NB: There must be at least two occupied child nodes
    //  -- avoid coding the occupancyB if it is implied.
    bit = 1;
    if (numOccupiedAcc >= 1)
      bit = arithmeticDecoder->decode(ctxOccupancy.b6[numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << 6;

    bit = 1;
    if (numOccupiedAcc >= 2)
      bit = arithmeticDecoder->decode(ctxOccupancy.b7[numOccupiedAcc]);
    occupancy |= bit << 0;

    return occupancy;
  }

  //---------------------------------------------------------------------------
  // decode occupancy bits (neighPattern10 != 0 case)

  int decodeOccupancyNeighNZ(
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    CtxModelOctreeOccupancy& ctxOccupancy,
    int neighPattern10)
  {
    int occupancy = 0;
    int partialOccupancy = 0;
    int idx;
    int bit;

    idx = neighPattern10;
    bit = arithmeticDecoder->decode(ctxOccupancy.b0[idx]);
    partialOccupancy |= bit << 0;
    occupancy |= bit << 1;

    idx = (neighPattern10 << 1) + partialOccupancy;
    bit = arithmeticDecoder->decode(ctxOccupancy.b1[idx]);
    partialOccupancy |= bit << 1;
    occupancy |= bit << 7;

    idx = (neighPattern10 << 2) + partialOccupancy;
    bit = arithmeticDecoder->decode(ctxOccupancy.b2[idx]);
    partialOccupancy |= bit << 2;
    occupancy |= bit << 5;

    idx = (neighPattern10 << 3) + partialOccupancy;
    bit = arithmeticDecoder->decode(ctxOccupancy.b3[idx]);
    partialOccupancy |= bit << 3;
    occupancy |= bit << 3;

    // todo(df): merge constants into lut.
    idx = ((neighPattern10 - 1) << 4) + partialOccupancy;
    idx = kOccMapBit4CtxIdx[idx] - 1 + 5;
    bit = arithmeticDecoder->decode(ctxOccupancy.b4[idx]);
    partialOccupancy |= bit << 4;
    occupancy |= bit << 2;

    idx = ((neighPattern10 - 1) << 5) + partialOccupancy;
    idx = kOccMapBit5CtxIdx[idx] - 1 + 6;
    bit = arithmeticDecoder->decode(ctxOccupancy.b5[idx]);
    partialOccupancy |= bit << 5;
    occupancy |= bit << 4;

    int neighPattern7 = kNeighPattern10to7[neighPattern10];
    idx = ((neighPattern7 - 1) << 6) + partialOccupancy;
    idx = kOccMapBit6CtxIdx[idx] - 1 + 7;
    bit = arithmeticDecoder->decode(ctxOccupancy.b6[idx]);
    partialOccupancy += bit << 6;
    occupancy |= bit << 6;

    int neighPattern5 = kNeighPattern7to5[neighPattern7];
    idx = ((neighPattern5 - 1) << 7) + partialOccupancy;
    idx = kOccMapBit7CtxIdx[idx] - 1 + 8;
    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    bit = 1;
    if (partialOccupancy)
      bit = arithmeticDecoder->decode(ctxOccupancy.b7[idx]);
    occupancy |= bit << 0;

    return occupancy;
  }

  //-------------------------------------------------------------------------
  // decode node occupancy bits
  //
  uint32_t decodeGeometryOccupancy(
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    o3dgc::Adaptive_Bit_Model& ctxSingleChild,
    o3dgc::Static_Bit_Model& ctxEquiProb,
    CtxModelOctreeOccupancy& ctxOccupancy,
    const PCCOctree3Node& node0)
  {
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
      } else {
        occupancy = decodeOccupancyNeighZ(arithmeticDecoder, ctxOccupancy);
      }
    } else {
      occupancy = decodeOccupancyNeighNZ(
        arithmeticDecoder, ctxOccupancy, neighPattern10);

      occupancy = mapGeometryOccupancyInv(occupancy, neighPattern);
    }
    return occupancy;
  }

  //-------------------------------------------------------------------------
  // Decode a position of a point in a given volume.

  PCCVector3<uint32_t> decodePointPosition(
    int nodeSizeLog2,
    o3dgc::Arithmetic_Codec* arithmeticDecoder,
    o3dgc::Static_Bit_Model& ctxPointPosBlock)
  {
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
    OutputIt outputPoints)
  {
    bool isDirectMode = arithmeticDecoder->decode(ctxBlockSkipTh);
    if (!isDirectMode) {
      return 0;
    }

    int numPoints = 1;
    if (arithmeticDecoder->decode(ctxNumIdcmPointsEq1))
      numPoints++;

    for (int i = 0; i < numPoints; i++) {
      // convert node-relative position to world position
      PCCVector3<uint32_t> pos = node.pos
        + decodePointPosition(nodeSizeLog2, arithmeticDecoder,
                              ctxPointPosBlock);

      *(outputPoints++) = {double(pos[0]), double(pos[1]), double(pos[2])};
    }

    return numPoints;
  }

  //-------------------------------------------------------------------------

  void decodePositions(const PayloadBuffer& buf, PCCPointSet3& pointCloud)
  {
    assert(buf.type == PayloadType::kGeometryBrick);

    int gbhSize;
    GeometryBrickHeader gbh = parseGbh(*_gps, buf, &gbhSize);
    minPositions.x() = gbh.geomBoxOrigin.x();
    minPositions.y() = gbh.geomBoxOrigin.y();
    minPositions.z() = gbh.geomBoxOrigin.z();

    pointCloud.resize(gbh.geom_num_points);

    o3dgc::Arithmetic_Codec arithmeticDecoder(
      int(buf.size()) - gbhSize,
      reinterpret_cast<uint8_t*>(const_cast<char*>(buf.data() + gbhSize)));
    arithmeticDecoder.start_decoder();

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

    // the current node dimension (log2)
    int nodeSizeLog2 = gbh.geom_max_node_size_log2;

    // push the first node
    fifo.emplace_back();
    PCCOctree3Node& node00 = fifo.back();
    node00.start = uint32_t(0);
    node00.end = uint32_t(pointCloud.getPointCount());
    node00.pos = uint32_t(0);
    node00.neighPattern = 0;
    node00.numSiblingsPlus1 = 8;
    node00.siblingOccupancy = 0;

    size_t processedPointCount = 0;
    std::vector<uint32_t> values;

    auto fifoCurrLvlEnd = fifo.end();

    // this counter represents fifo.end() - fifoCurrLvlEnd().
    // ie, the number of nodes added to the next level of the tree
    int numNodesNextLvl = 0;

    PCCVector3<uint32_t> occupancyAtlasOrigin(0xffffffff);
    MortonMap3D occupancyAtlas;
    if (_gps->neighbour_avail_boundary_log2) {
      occupancyAtlas.resize(_gps->neighbour_avail_boundary_log2);
      occupancyAtlas.clear();
    }

    for (; !fifo.empty(); fifo.pop_front()) {
      if (fifo.begin() == fifoCurrLvlEnd) {
        // transition to the next level
        fifoCurrLvlEnd = fifo.end();
        nodeSizeLog2--;
        numNodesNextLvl = 0;
        occupancyAtlasOrigin = 0xffffffff;
      }
      PCCOctree3Node& node0 = fifo.front();

      if (_gps->neighbour_avail_boundary_log2) {
        updateGeometryOccupancyAtlas(
          node0.pos, nodeSizeLog2, fifo, fifoCurrLvlEnd, &occupancyAtlas,
          &occupancyAtlasOrigin);

        node0.neighPattern =
          makeGeometryNeighPattern(node0.pos, nodeSizeLog2, occupancyAtlas);
      }

      // decode occupancy pattern
      uint8_t occupancy = decodeGeometryOccupancy(
        &arithmeticDecoder, ctxSingleChild, ctxEquiProb, ctxOccupancy, node0);

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

          if (!_gps->geom_unique_points_flag) {
            numPoints = decodePositionLeafNumPoints(
              &arithmeticDecoder, ctxSinglePointPerBlock, ctxEquiProb,
              ctxPointCountPerBlock);
          }

          const PCCVector3D point(
            node0.pos[0] + (x << childSizeLog2),
            node0.pos[1] + (y << childSizeLog2),
            node0.pos[2] + (z << childSizeLog2));

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
        child.siblingOccupancy = occupancy;

        bool idcmEnabled = _gps->inferred_direct_coding_mode_enabled_flag;
        if (isDirectModeEligible(idcmEnabled, nodeSizeLog2, node0, child)) {
          int numPoints = decodeDirectPosition(
            childSizeLog2, child, &arithmeticDecoder, ctxBlockSkipTh,
            ctxNumIdcmPointsEq1, ctxEquiProb,
            &pointCloud[processedPointCount]);
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

        if (!_gps->neighbour_avail_boundary_log2) {
          updateGeometryNeighState(
            _gps->neighbour_context_restriction_flag, fifo.end(),
            numNodesNextLvl, childSizeLog2, child, i, node0.neighPattern,
            occupancy);
        }
      }
    }

    arithmeticDecoder.stop_decoder();
  }

  //==========================================================================

public:
  void inverseQuantization(PCCPointSet3& pointCloud)
  {
    const size_t pointCount = pointCloud.getPointCount();
    const double invScale = 1.0 / _sps->seq_source_geom_scale_factor;

    for (size_t i = 0; i < pointCount; ++i) {
      auto& point = pointCloud[i];
      for (size_t k = 0; k < 3; ++k) {
        point[k] = point[k] * invScale + minPositions[k];
      }
    }
  }

  //==========================================================================

private:
  PCCVector3D minPositions;

  // The point cloud currently being decoded
  PCCPointSet3 _currentPointCloud;

  // Received parameter sets, mapping parameter set id -> parameterset
  std::map<int, SequenceParameterSet> _spss;
  std::map<int, GeometryParameterSet> _gpss;
  std::map<int, AttributeParameterSet> _apss;

  // The active SPS
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;
};
}  // namespace pcc

#endif /* PCCTMC3Decoder_h */
