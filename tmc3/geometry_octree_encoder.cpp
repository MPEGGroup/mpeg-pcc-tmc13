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

#include "geometry.h"

#include "DualLutCoder.h"
#include "OctreeNeighMap.h"
#include "geometry_octree.h"
#include "geometry_intra_pred.h"
#include "io_hls.h"
#include "tables.h"
#include "quantization.h"

#include <set>

namespace pcc {

//============================================================================

class GeometryOctreeEncoder {
public:
  GeometryOctreeEncoder(
    const GeometryParameterSet& gps, EntropyEncoder* arithmeticEncoder);

  void beginOctreeLevel();

  int encodePositionLeafNumPoints(int count);

  void encodeOccupancyNeighZ(
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int occupancySkip);

  void encodeOccupancyNeighNZ(
    int neighPattern,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int occupancySkip);

  void encodeOccupancyBitwise(
    int neighPattern,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int occupancySkip);

  void encodeOccupancyBytewise(int neighPattern, int mappedOccupancy);

  void encodeOccupancy(
    int neighPattern,
    int occupancy,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int occupancyAdjGt0,
    int occupancyAdjGt1,
    int occupancyAdjUnocc,
    int occupancySkip);

  void encodePointPosition(
    const Vec3<int>& nodeSizeLog2, const Vec3<uint32_t>& pos);

  void encodeQpOffset(int dqp);

  bool encodeDirectPosition(
    const Vec3<int>& nodeSizeLog2,
    const PCCOctree3Node& node,
    const PCCPointSet3& pointCloud);

private:
  // selects between the bitwise and bytewise occupancy coders
  const bool _useBitwiseOccupancyCoder;

  const uint8_t (&_neighPattern64toR1)[64];

  EntropyEncoder* _arithmeticEncoder;
  StaticBitModel _ctxEquiProb;
  AdaptiveBitModel _ctxSingleChild;
  AdaptiveBitModel _ctxSinglePointPerBlock;
  AdaptiveBitModel _ctxPointCountPerBlock;
  AdaptiveBitModel _ctxBlockSkipTh;
  AdaptiveBitModel _ctxNumIdcmPointsEq1;

  AdaptiveBitModel _ctxQpOffsetIsZero;
  AdaptiveBitModel _ctxQpOffsetSign;
  AdaptiveBitModel _ctxQpOffsetAbsEgl;

  // For bitwise occupancy coding
  CtxModelOctreeOccupancy _ctxOccupancy;
  CtxMapOctreeOccupancy _ctxIdxMaps[18];

  // For bytewise occupancy coding
  DualLutCoder<true> _bytewiseOccupancyCoder[10];
};

//============================================================================

GeometryOctreeEncoder::GeometryOctreeEncoder(
  const GeometryParameterSet& gps, EntropyEncoder* arithmeticEncoder)
  : _useBitwiseOccupancyCoder(gps.bitwise_occupancy_coding_flag)
  , _neighPattern64toR1(neighPattern64toR1(gps))
  , _arithmeticEncoder(arithmeticEncoder)
  , _ctxOccupancy(gps.geom_occupancy_ctx_reduction_factor)
{
  if (!_useBitwiseOccupancyCoder) {
    for (int i = 0; i < 10; i++)
      _bytewiseOccupancyCoder[i].init(kDualLutOccupancyCoderInit[i]);
  }
}

//============================================================================

void
GeometryOctreeEncoder::beginOctreeLevel()
{
  for (int i = 0; i < 10; i++) {
    _bytewiseOccupancyCoder[i].resetLut();
  }
}

//============================================================================
// Encode the number of points in a leaf node of the octree.

int
GeometryOctreeEncoder::encodePositionLeafNumPoints(int count)
{
  if (count == 1) {
    _arithmeticEncoder->encode(1, _ctxSinglePointPerBlock);
  } else {
    _arithmeticEncoder->encode(0, _ctxSinglePointPerBlock);
    _arithmeticEncoder->encodeExpGolomb(
      uint32_t(count - 2), 0, _ctxEquiProb, _ctxPointCountPerBlock);
  }

  return count;
}

//-------------------------------------------------------------------------
// encode occupancy bits (neighPattern10 == 0 case)

void
GeometryOctreeEncoder::encodeOccupancyNeighZ(
  int mappedOccupancy,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int occupancySkip)
{
  int minOccupied = 2;
  int numOccupiedAcc = 0;
  int numCoded = 0;
  const int maxOccupied = numMaxOccupided[occupancySkip];

  for (int i = 0; i < 8; i++) {
    if (occupancySkip != 0) {
      if (
        (occupancySkip & 1) && (kOccBitCodingOrder[i] & 1))  // skip when z = 1
        continue;
      if (
        (occupancySkip & 2) && (kOccBitCodingOrder[i] & 2))  // skip when y = 1
        continue;
      if (
        (occupancySkip & 4) && (kOccBitCodingOrder[i] & 4))  // skip when x = 1
        continue;
    }

    int bit = (mappedOccupancy >> kOccBitCodingOrder[i]) & 1;
    int bitIsPredicted = (mappedOccIsPredicted >> kOccBitCodingOrder[i]) & 1;
    int bitPrediction = (mappedOccPrediction >> kOccBitCodingOrder[i]) & 1;
    int bitAdjGt0 = (mappedOccAdjGt0 >> kOccBitCodingOrder[i]) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> kOccBitCodingOrder[i]) & 1;
    int bitAdjUnocc = (mappedOccAdjUnocc >> kOccBitCodingOrder[i]) & 1;

    int numAdj = bitAdjGt0 + bitAdjGt1;
    int idxAdj = bitAdjUnocc + 2 * numAdj;
    if (i > 4) {
      static const int8_t kCtxIdxAdjReduc567[6] = {0, 0, 1, 2, 3, 3};
      idxAdj = kCtxIdxAdjReduc567[idxAdj];
    }

    int ctxIdxMapIdx = 3 * idxAdj + bitIsPredicted + bitPrediction;
    auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];

    int idx = numOccupiedAcc;
    int ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);

    // NB: There must be at least minOccupied child nodes
    //  -- avoid coding the occupancyBit if it is implied.
    if (numOccupiedAcc > minOccupied + numCoded - maxOccupied) {
      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);
    }
    numCoded++;
    numOccupiedAcc += bit;
  }
}

//-------------------------------------------------------------------------
// encode occupancy bits (neighPattern10 != 0 case)

void
GeometryOctreeEncoder::encodeOccupancyNeighNZ(
  int neighPattern,
  int mappedOccupancy,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int occupancySkip)
{
  // code occupancy using the neighbour configuration context
  // with reduction from 64 states to 9 (or 6).
  int neighPatternR1 = _neighPattern64toR1[neighPattern];

  //  int neighPattern9 = kNeighPattern64to9[neighPattern];
  int neighPattern5 = kNeighPattern9to5[neighPatternR1];
  int neighPattern3 = kNeighPattern9to3[neighPatternR1];

  //  int neighPattern7 = kNeighPattern10to7[neighPattern10];
  //  int neighPattern5 = kNeighPattern7to5[neighPattern7];

  uint32_t partialOccupancy = 0;

  const uint8_t* map = kNeighPatternInvMap[neighPattern];
  const int maxOccupied = numMaxOccupided[occupancySkip];
  int numCoded = 0;

  // NB: it is impossible for pattern to be 0 (handled in Z case).
  for (int i = 0; i < 8; i++) {
    if (occupancySkip != 0) {
      if (
        (occupancySkip & 1)
        && (map[kOccBitCodingOrder[i]] & 1))  // skip when z = 1
        continue;
      if (
        (occupancySkip & 2)
        && (map[kOccBitCodingOrder[i]] & 2))  // skip when y = 1
        continue;
      if (
        (occupancySkip & 4)
        && (map[kOccBitCodingOrder[i]] & 4))  // skip when x = 1
        continue;
    }

    int idx;
    if (i < 4) {
      idx = ((neighPatternR1 - 1) << i) + partialOccupancy + i + 1;
    } else if (i < 6) {
      idx = ((neighPattern5 - 1) << i) + partialOccupancy + i + 1;
    } else if (i == 6) {
      idx = ((neighPattern3 - 1) << i) + partialOccupancy + i + 1;
    } else if (i == 7) {
      idx = partialOccupancy + i + 1;
    } else {
      // work around clang -Wsometimes-uninitialized fault
      break;
    }

    int bit = (mappedOccupancy >> kOccBitCodingOrder[i]) & 1;
    int bitIsPredicted = (mappedOccIsPredicted >> kOccBitCodingOrder[i]) & 1;
    int bitPrediction = (mappedOccPrediction >> kOccBitCodingOrder[i]) & 1;
    int bitAdjGt0 = (mappedOccAdjGt0 >> kOccBitCodingOrder[i]) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> kOccBitCodingOrder[i]) & 1;
    int bitAdjUnocc = (mappedOccAdjUnocc >> kOccBitCodingOrder[i]) & 1;

    int numAdj = bitAdjGt0 + bitAdjGt1;
    int idxAdj = bitAdjUnocc + 2 * numAdj;
    if (i > 4) {
      static const int8_t kCtxIdxAdjReduc567[6] = {0, 0, 1, 2, 3, 3};
      idxAdj = kCtxIdxAdjReduc567[idxAdj];
    }

    int ctxIdxMapIdx = 3 * idxAdj + bitIsPredicted + bitPrediction;
    auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];

    int ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);

    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    if (numCoded < maxOccupied - 1 || partialOccupancy) {
      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);
    }
    numCoded++;
    partialOccupancy |= bit << i;
    if (numCoded == maxOccupied)
      break;
  }
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeOccupancyBitwise(
  int neighPattern,
  int mappedOccupancy,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int occupancySkip)
{
  if (neighPattern == 0) {
    encodeOccupancyNeighZ(
      mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
      mappedOccAdjGt0, mappedOccAdjGt1, mappedOccAdjUnocc, occupancySkip);
    return;
  }

  encodeOccupancyNeighNZ(
    neighPattern, mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
    mappedOccAdjGt0, mappedOccAdjGt1, mappedOccAdjUnocc, occupancySkip);
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeOccupancyBytewise(
  int neighPattern, int mappedOccupancy)
{
  // code occupancy using the neighbour configuration context
  // with reduction from 64 states to 10 (or 6).
  int neighPatternR1 = _neighPattern64toR1[neighPattern];
  auto& bytewiseCoder = _bytewiseOccupancyCoder[neighPatternR1];
  bytewiseCoder.encode(mappedOccupancy, _arithmeticEncoder);
}

//-------------------------------------------------------------------------
// decode node occupancy bits
//

void
GeometryOctreeEncoder::encodeOccupancy(
  int neighPattern,
  int occupancy,
  int occupancyIsPred,
  int occupancyPred,
  int occupancyAdjGt0,
  int occupancyAdjGt1,
  int occupancyAdjUnocc,
  int occupancySkip)
{
  if (neighPattern == 0) {
    bool singleChild = !popcntGt1(occupancy);
    _arithmeticEncoder->encode(singleChild, _ctxSingleChild);

    if (singleChild) {
      // no siblings => encode index = (z,y,x) not 8bit pattern
      if (!(occupancySkip & 1))
        _arithmeticEncoder->encode(!!(occupancy & 0xaa), _ctxEquiProb);  // z
      if (!(occupancySkip & 2))
        _arithmeticEncoder->encode(!!(occupancy & 0xcc), _ctxEquiProb);  // y
      if (!(occupancySkip & 4))
        _arithmeticEncoder->encode(!!(occupancy & 0xf0), _ctxEquiProb);  // x
      return;
    }
  }

  uint32_t mapOcc = mapGeometryOccupancy(occupancy, neighPattern);
  uint32_t mapOccIsP = mapGeometryOccupancy(occupancyIsPred, neighPattern);
  uint32_t mapOccP = mapGeometryOccupancy(occupancyPred, neighPattern);
  uint32_t mapAdjGt0 = mapGeometryOccupancy(occupancyAdjGt0, neighPattern);
  uint32_t mapAdjGt1 = mapGeometryOccupancy(occupancyAdjGt1, neighPattern);
  uint32_t mapAdjUnocc = mapGeometryOccupancy(occupancyAdjUnocc, neighPattern);

  if (_useBitwiseOccupancyCoder)
    encodeOccupancyBitwise(
      neighPattern, mapOcc, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1,
      mapAdjUnocc, occupancySkip);
  else
    encodeOccupancyBytewise(neighPattern, mapOcc);
}

//-------------------------------------------------------------------------
// Encode a position of a point in a given volume.
void
GeometryOctreeEncoder::encodePointPosition(
  const Vec3<int>& nodeSizeLog2, const Vec3<uint32_t>& pos)
{
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2[k] <= 0)
      continue;
    for (int mask = 1 << (nodeSizeLog2[k] - 1); mask; mask >>= 1) {
      _arithmeticEncoder->encode(!!(pos[k] & mask), _ctxEquiProb);
    }
  }
}

//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).
void
GeometryOctreeEncoder::encodeQpOffset(int dqp)
{
  _arithmeticEncoder->encode(dqp == 0, _ctxQpOffsetIsZero);
  if (dqp == 0) {
    return;
  }
  _arithmeticEncoder->encode(dqp > 0, _ctxQpOffsetSign);
  _arithmeticEncoder->encodeExpGolomb(
    abs(dqp) - 1, 0, _ctxEquiProb, _ctxQpOffsetAbsEgl);
}

//-------------------------------------------------------------------------

template<typename It>
void
calculateNodeQps(int baseQp, It nodesBegin, It nodesEnd)
{
  // determine delta qp for each node based on the point density
  int lowQp = baseQp - 4 >= 4 ? baseQp - 4 : 4;
  int mediumQp = baseQp;
  int highQp = baseQp + 6;
  std::vector<int> numPointsInNode;
  std::vector<double> cum_prob;
  int32_t numPointsInLvl = 0;
  for (auto it = nodesBegin; it != nodesEnd; ++it) {
    numPointsInNode.push_back(it->end - it->start);
    numPointsInLvl += it->end - it->start;
  }
  std::sort(numPointsInNode.begin(), numPointsInNode.end());
  double cc = 0;
  for (auto num : numPointsInNode) {
    cc += num;
    cum_prob.push_back(cc / numPointsInLvl);
  }
  int th1 = -1, th2 = -1;
  for (int i = 0; i < cum_prob.size(); i++) {
    if (th1 == -1 && cum_prob[i] > 0.05) {
      th1 = numPointsInNode[i];
    } else if (th2 == -1 && cum_prob[i] > 0.6)
      th2 = numPointsInNode[i];
  }
  for (auto it = nodesBegin; it != nodesEnd; ++it) {
    if (it->end - it->start < th1) {
      it->qp = highQp;
    } else if (it->end - it->start < th2)
      it->qp = mediumQp;
    else
      it->qp = lowQp;
  }
}

//-------------------------------------------------------------------------

void
geometryQuantization(
  PCCPointSet3& pointCloud, PCCOctree3Node& node, Vec3<int> nodeSizeLog2)
{
  QuantizerGeom quantizer = QuantizerGeom(node.qp);
  int qpShift = (node.qp - 4) / 6;

  for (int k = 0; k < 3; k++) {
    int quantBitsMask = (1 << nodeSizeLog2[k]) - 1;
    uint32_t clipMax = ((1 << nodeSizeLog2[k]) >> qpShift) - 1;

    for (int i = node.start; i < node.end; i++) {
      uint32_t pos = uint32_t(pointCloud[i][k]);
      uint32_t quantPos = quantizer.quantize(pos & quantBitsMask);
      quantPos = PCCClip(quantPos, 0, clipMax);

      // NB: this representation is: |pppppp00qqq|, which isn't the
      // same used by the decoder:  (|ppppppqqq|00)
      pointCloud[i][k] = (pos & ~quantBitsMask) | quantPos;
    }
  }
}

//-------------------------------------------------------------------------

void
geometryScale(
  PCCPointSet3& pointCloud, PCCOctree3Node& node, Vec3<int> quantNodeSizeLog2)
{
  QuantizerGeom quantizer = QuantizerGeom(node.qp);
  int qpShift = (node.qp - 4) / 6;

  for (int k = 0; k < 3; k++) {
    int quantBitsMask = (1 << quantNodeSizeLog2[k]) - 1;
    for (int i = node.start; i < node.end; i++) {
      uint32_t pos = uint32_t(pointCloud[i][k]);
      uint32_t quantPos = pos & quantBitsMask;
      pointCloud[i][k] = (pos & ~quantBitsMask) | quantizer.scale(quantPos);
    }
  }
}

//-------------------------------------------------------------------------

void
checkDuplicatePoints(
  PCCPointSet3& pointCloud,
  PCCOctree3Node& node,
  std::vector<int>& pointIdxToDmIdx)
{
  auto first = PCCPointSet3::iterator(&pointCloud, node.start);
  auto last = PCCPointSet3::iterator(&pointCloud, node.end);

  std::set<Vec3<double>> uniquePointsSet;
  for (auto i = first; i != last;) {
    if (uniquePointsSet.find(**i) == uniquePointsSet.end()) {
      uniquePointsSet.insert(**i);
      i++;
    } else {
      std::iter_swap(i, last - 1);
      last--;
      pointIdxToDmIdx[--node.end] = -2;  // mark as duplicate
    }
  }
}

//-------------------------------------------------------------------------

bool
GeometryOctreeEncoder::encodeDirectPosition(
  const Vec3<int>& nodeSizeLog2,
  const PCCOctree3Node& node,
  const PCCPointSet3& pointCloud)
{
  int numPoints = node.end - node.start;
  if (numPoints > MAX_NUM_DM_LEAF_POINTS) {
    _arithmeticEncoder->encode(0, _ctxBlockSkipTh);
    return false;
  }

  _arithmeticEncoder->encode(1, _ctxBlockSkipTh);
  _arithmeticEncoder->encode(numPoints > 1, _ctxNumIdcmPointsEq1);

  for (auto idx = node.start; idx < node.end; idx++) {
    encodePointPosition(
      nodeSizeLog2,
      Vec3<uint32_t>{uint32_t(pointCloud[idx][0]),
                     uint32_t(pointCloud[idx][1]),
                     uint32_t(pointCloud[idx][2])});
  }

  return true;
}

//-------------------------------------------------------------------------

void
encodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  EntropyEncoder* arithmeticEncoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining)
{
  GeometryOctreeEncoder encoder(gps, arithmeticEncoder);

  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  pcc::ringbuf<PCCOctree3Node> fifo(pointCloud.getPointCount() + 1);

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(pointCloud.getPointCount());
  node00.pos = uint32_t(0);
  node00.neighPattern = 0;
  node00.numSiblingsPlus1 = 8;
  node00.siblingOccupancy = 0;
  node00.qp = 4;

  // map of pointCloud idx to DM idx, used to reorder the points
  // after coding.
  std::vector<int> pointIdxToDmIdx(int(pointCloud.getPointCount()), -1);
  int nextDmIdx = 0;

  size_t processedPointCount = 0;
  std::vector<uint32_t> values;

  auto fifoCurrLvlEnd = fifo.end();

  // represents the largest dimension of the current node
  // NB: this is equal to the total depth of the tree
  int nodeMaxDimLog2 = gbh.geomMaxNodeSizeLog2(gps);

  // size of the current node (each dimension can vary due to qtbt)
  Vec3<int> nodeSizeLog2 = gbh.geomMaxNodeSizeLog2Xyz(gps);

  // update implicit qtbt parameters
  int maxNumImplicitQtbtBeforeOt = gps.max_num_implicit_qtbt_before_ot;
  int minSizeImplicitQtbt = gps.min_implicit_qtbt_size_log2;
  updateImplicitQtBtParameters(
    nodeSizeLog2, gps.trisoup_node_size_log2, &maxNumImplicitQtbtBeforeOt,
    &minSizeImplicitQtbt);

  // implicit qtbt for child nodes
  Vec3<int> childSizeLog2 = implicitQtBtDecision(
    nodeSizeLog2, maxNumImplicitQtbtBeforeOt, minSizeImplicitQtbt);

  // prepare parameters for partition and occupancy coding
  Vec3<int> pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);
  int occupancySkip = nonSplitQtBtAxes(nodeSizeLog2, childSizeLog2);
  int atlasShift = 7;

  // this counter represents fifo.end() - fifoCurrLvlEnd().
  // ie, the number of nodes added to the next level of the tree
  int numNodesNextLvl = 0;

  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2) {
    occupancyAtlas.resize(gps.neighbour_avail_boundary_log2);
    occupancyAtlas.clear();
  }
  Vec3<uint32_t> occupancyAtlasOrigin(0xffffffff);

  // the node size where quantisation is performed
  Vec3<int> quantNodeSizeLog2 = 0;
  int numLvlsUntilQuantization = -1;
  if (gps.geom_scaling_enabled_flag) {
    numLvlsUntilQuantization = 0;
    if (gbh.geom_octree_qp_offset_enabled_flag)
      numLvlsUntilQuantization = gbh.geom_octree_qp_offset_depth;
  }

  int sliceQp = gps.geom_base_qp + gbh.geom_slice_qp_offset;

  // applied to the root node, just set the node qp
  if (numLvlsUntilQuantization == 0) {
    quantNodeSizeLog2 = nodeSizeLog2;
    node00.qp = sliceQp;
  }

  for (; !fifo.empty(); fifo.pop_front()) {
    if (fifo.begin() == fifoCurrLvlEnd) {
      // transition to the next level
      fifoCurrLvlEnd = fifo.end();
      numNodesNextLvl = 0;
      occupancyAtlasOrigin = 0xffffffff;

      Vec3<int> parentNodeSizeLog2 = nodeSizeLog2;
      // implicit qtbt for current node
      nodeSizeLog2 = implicitQtBtDecision(
        nodeSizeLog2, maxNumImplicitQtbtBeforeOt, minSizeImplicitQtbt);

      // if one dimension is not split, atlasShift[k] = 0
      atlasShift = 7 & ~nonSplitQtBtAxes(parentNodeSizeLog2, nodeSizeLog2);

      if (maxNumImplicitQtbtBeforeOt)
        maxNumImplicitQtbtBeforeOt--;

      // if all dimensions have same size, then performing octree for remaining nodes
      if (
        nodeSizeLog2[0] == nodeSizeLog2[1]
        && nodeSizeLog2[1] == nodeSizeLog2[2])
        minSizeImplicitQtbt = 0;

      // implicit qtbt for child nodes
      childSizeLog2 = implicitQtBtDecision(
        nodeSizeLog2, maxNumImplicitQtbtBeforeOt, minSizeImplicitQtbt);

      pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);

      nodeMaxDimLog2--;
      encoder.beginOctreeLevel();

      // allow partial tree encoding using trisoup
      if (nodeMaxDimLog2 == gps.trisoup_node_size_log2)
        break;

      // determing a per node QP at the appropriate level
      // NB: this has no effect here if geom_octree_qp_offset_depth=0
      if (--numLvlsUntilQuantization == 0) {
        quantNodeSizeLog2 = nodeSizeLog2;
        calculateNodeQps(gps.geom_base_qp, fifo.begin(), fifoCurrLvlEnd);
      }
    }

    PCCOctree3Node& node0 = fifo.front();

    // encode delta qp for each octree block
    if (
      numLvlsUntilQuantization == 0 && gbh.geom_octree_qp_offset_enabled_flag)
      encoder.encodeQpOffset(node0.qp - sliceQp);

    int shiftBits = (node0.qp - 4) / 6;
    int effectiveNodeMaxDimLog2 = nodeMaxDimLog2 - shiftBits;
    auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;
    auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;

    // todo(??): the following needs to be reviewed, it is added to make
    // quantisation work with qtbt.
    Vec3<int> actualNodeSizeLog2, actualChildSizeLog2;
    for (int k = 0; k < 3; k++) {
      actualNodeSizeLog2[k] = std::max(nodeSizeLog2[k], shiftBits);
      actualChildSizeLog2[k] = std::max(childSizeLog2[k], shiftBits);
    }
    // todo(??): atlasShift may be wrong too
    occupancySkip = nonSplitQtBtAxes(actualNodeSizeLog2, actualChildSizeLog2);

    if (numLvlsUntilQuantization == 0) {
      geometryQuantization(pointCloud, node0, quantNodeSizeLog2);
      if (gps.geom_unique_points_flag)
        checkDuplicatePoints(pointCloud, node0, pointIdxToDmIdx);
    }

    // split the current node into 8 children
    //  - perform an 8-way counting sort of the current node's points
    //  - (later) map to child nodes
    std::array<int, 8> childCounts = {};
    countingSort(
      PCCPointSet3::iterator(&pointCloud, node0.start),
      PCCPointSet3::iterator(&pointCloud, node0.end), childCounts,
      [=](const PCCPointSet3::Proxy& proxy) {
        const auto& point = *proxy;
        return !!(int(point[2]) & (pointSortMask[2] >> shiftBits))
          | (!!(int(point[1]) & (pointSortMask[1] >> shiftBits)) << 1)
          | (!!(int(point[0]) & (pointSortMask[0] >> shiftBits)) << 2);
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

    int occupancyAdjacencyGt0 = 0;
    int occupancyAdjacencyGt1 = 0;
    int occupancyAdjacencyUnocc = 0;

    if (gps.neighbour_avail_boundary_log2) {
      updateGeometryOccupancyAtlas(
        node0.pos, atlasShift, fifo, fifoCurrLvlEnd, &occupancyAtlas,
        &occupancyAtlasOrigin);

      GeometryNeighPattern gnp = makeGeometryNeighPattern(
        gps.adjacent_child_contextualization_enabled_flag, node0.pos,
        atlasShift, occupancyAtlas);

      node0.neighPattern = gnp.neighPattern;
      occupancyAdjacencyGt0 = gnp.adjacencyGt0;
      occupancyAdjacencyGt1 = gnp.adjacencyGt1;
      occupancyAdjacencyUnocc = gnp.adjacencyUnocc;
    }

    int occupancyIsPredicted = 0;
    int occupancyPrediction = 0;

    // generate intra prediction
    if (effectiveNodeMaxDimLog2 < gps.intra_pred_max_node_size_log2) {
      predictGeometryOccupancyIntra(
        occupancyAtlas, node0.pos, atlasShift, &occupancyIsPredicted,
        &occupancyPrediction);
    }

    // update atlas for advanced neighbours
    if (gps.neighbour_avail_boundary_log2) {
      updateGeometryOccupancyAtlasOccChild(
        node0.pos, occupancy, &occupancyAtlas);
    }

    // when all points are quantized to a single point
    if (!isLeafNode(effectiveNodeSizeLog2)) {
      // encode child occupancy map
      assert(occupancy > 0);
      assert(occupancySkip != 7);
      encoder.encodeOccupancy(
        node0.neighPattern, occupancy, occupancyIsPredicted,
        occupancyPrediction, occupancyAdjacencyGt0, occupancyAdjacencyGt1,
        occupancyAdjacencyUnocc, occupancySkip);
    }

    // Leaf nodes are immediately coded.  No further splitting occurs.
    if (isLeafNode(effectiveChildSizeLog2)) {
      int childStart = node0.start;

      // inverse quantise any quantised positions
      geometryScale(pointCloud, node0, quantNodeSizeLog2);

      for (int i = 0; i < 8; i++) {
        if (!childCounts[i]) {
          // child is empty: skip
          continue;
        }

        int childEnd = childStart + childCounts[i];
        for (auto idx = childStart; idx < childEnd; idx++)
          pointIdxToDmIdx[idx] = nextDmIdx++;

        processedPointCount += childCounts[i];
        childStart = childEnd;

        // if the bitstream is configured to represent unique points,
        // no point count is sent.
        if (gps.geom_unique_points_flag) {
          assert(childCounts[i] == 1);
          processedPointCount++;
          continue;
        }

        encoder.encodePositionLeafNumPoints(childCounts[i]);
        processedPointCount += childCounts[i];
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

      child.qp = node0.qp;
      // only shift position if an occupancy bit was coded for the axis
      child.pos[0] = (node0.pos[0] << !(occupancySkip & 4)) + x;
      child.pos[1] = (node0.pos[1] << !(occupancySkip & 2)) + y;
      child.pos[2] = (node0.pos[2] << !(occupancySkip & 1)) + z;

      child.start = childPointsStartIdx;
      childPointsStartIdx += childCounts[i];
      child.end = childPointsStartIdx;
      child.numSiblingsPlus1 = numSiblings;
      child.siblingOccupancy = occupancy;

      bool idcmEnabled = gps.inferred_direct_coding_mode_enabled_flag;
      if (isDirectModeEligible(
            idcmEnabled, effectiveNodeMaxDimLog2, node0, child)) {
        bool directModeUsed = encoder.encodeDirectPosition(
          effectiveChildSizeLog2, child, pointCloud);

        if (directModeUsed) {
          // inverse quantise any quantised positions
          geometryScale(pointCloud, node0, quantNodeSizeLog2);

          // point reordering to match decoder's order
          for (auto idx = child.start; idx < child.end; idx++)
            pointIdxToDmIdx[idx] = nextDmIdx++;
          processedPointCount += child.end - child.start;

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
      if (!gps.neighbour_avail_boundary_log2) {
        updateGeometryNeighState(
          gps.neighbour_context_restriction_flag, fifo.end(), numNodesNextLvl,
          child, i, node0.neighPattern, occupancy);
      }
    }
  }

  // return partial coding result
  //  - add missing levels to node positions
  //  - inverse quantise the point cloud
  // todo(df): this does not yet support inverse quantisation of node.pos
  if (nodesRemaining) {
    for (auto& node : fifo) {
      for (int k = 0; k < 3; k++)
        node.pos[k] <<= nodeSizeLog2[k];
      geometryScale(pointCloud, node, quantNodeSizeLog2);
    }
    *nodesRemaining = std::move(fifo);
    return;
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
    } else if (dstIdx == -2) {  // ignore duplicated points
      continue;
    }

    pointCloud2[dstIdx] = pointCloud[i];
    if (pointCloud.hasColors())
      pointCloud2.setColor(dstIdx, pointCloud.getColor(i));
    if (pointCloud.hasReflectances())
      pointCloud2.setReflectance(dstIdx, pointCloud.getReflectance(i));
  }
  pointCloud2.resize(outIdx);
  swap(pointCloud, pointCloud2);
}

//============================================================================

void
encodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  EntropyEncoder* arithmeticEncoder)
{
  encodeGeometryOctree(gps, gbh, pointCloud, arithmeticEncoder, nullptr);
}

//============================================================================

}  // namespace pcc
