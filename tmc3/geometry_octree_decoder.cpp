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

namespace pcc {

//============================================================================

class GeometryOctreeDecoder {
public:
  GeometryOctreeDecoder(
    const GeometryParameterSet& gps, EntropyDecoder* arithmeticDecoder);

  void beginOctreeLevel();

  int decodePositionLeafNumPoints();

  int decodeOccupancyNeighZ(
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int occupancySkip);

  int decodeOccupancyNeighNZ(
    int neighPattern,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int occupancySkip);

  int decodeOccupancyBitwise(
    int neighPattern,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int occupancySkip);

  int decodeOccupancyBytewise(int neighPattern);

  uint32_t decodeOccupancy(
    int neighPattern,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int occupancyAdjGt0,
    int occupancyAdjGt1,
    int occupancyAdjUncc,
    int occupancySkip);

  Vec3<uint32_t> decodePointPosition(const Vec3<int>& nodeSizeLog2);

  int decodeQpOffset();

  template<class OutputIt>
  int decodeDirectPosition(
    const Vec3<int>& nodeSizeLog2,
    const PCCOctree3Node& node,
    OutputIt outputPoints);

private:
  // selects between the bitwise and bytewise occupancy coders
  const bool _useBitwiseOccupancyCoder;

  const uint8_t (&_neighPattern64toR1)[64];

  EntropyDecoder* _arithmeticDecoder;
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

GeometryOctreeDecoder::GeometryOctreeDecoder(
  const GeometryParameterSet& gps, EntropyDecoder* arithmeticDecoder)
  : _useBitwiseOccupancyCoder(gps.bitwise_occupancy_coding_flag)
  , _neighPattern64toR1(neighPattern64toR1(gps))
  , _arithmeticDecoder(arithmeticDecoder)
  , _ctxOccupancy(gps.geom_occupancy_ctx_reduction_factor)
{
  if (!_useBitwiseOccupancyCoder) {
    for (int i = 0; i < 10; i++)
      _bytewiseOccupancyCoder[i].init(kDualLutOccupancyCoderInit[i]);
  }
}
//============================================================================

void
GeometryOctreeDecoder::beginOctreeLevel()
{
  for (int i = 0; i < 10; i++) {
    _bytewiseOccupancyCoder[i].resetLut();
  }
}

//============================================================================
// Decode the number of points in a leaf node of the octree.

int
GeometryOctreeDecoder::decodePositionLeafNumPoints()
{
  const bool isSinglePoint =
    _arithmeticDecoder->decode(_ctxSinglePointPerBlock) != 0;

  int count = 1;
  if (!isSinglePoint) {
    count += 1
      + _arithmeticDecoder->decodeExpGolomb(
          0, _ctxEquiProb, _ctxPointCountPerBlock);
  }

  return count;
}

//---------------------------------------------------------------------------
// decode occupancy bits (neighPattern10 == 0 case)

int
GeometryOctreeDecoder::decodeOccupancyNeighZ(
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int occupancySkip)
{
  int minOccupied = 2;
  int numOccupiedAcc = 0;
  int occupancy = 0;
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
    int bit = 1;
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

    // NB: There must be at least two occupied child nodes
    //  -- avoid coding the occupancy bit if it is implied.
    if (numOccupiedAcc > minOccupied + numCoded - maxOccupied) {
      int ctxIdx = ctxIdxMap[i][numOccupiedAcc];
      bit = _arithmeticDecoder->decode(_ctxOccupancy[ctxIdx]);
    }
    numCoded++;
    ctxIdxMap.evolve(bit, &ctxIdxMap[i][numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << kOccBitCodingOrder[i];
  }

  return occupancy;
}

//---------------------------------------------------------------------------
// decode occupancy bits (neighPattern10 != 0 case)

int
GeometryOctreeDecoder::decodeOccupancyNeighNZ(
  int neighPattern,
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

  int occupancy = 0;
  int partialOccupancy = 0;

  const uint8_t* map = kNeighPatternInvMap[neighPattern];
  const int maxOccupied = numMaxOccupided[occupancySkip];
  int numCoded = 0;

  // NB: it is impossible for pattern to be 0 (handled in Z case).
  // NB: offsets are added since ctxIdxMap is shared between Z and NZ cases.
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
    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    int bit = 1;
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

    if (numCoded < maxOccupied - 1 || partialOccupancy) {
      int ctxIdx = ctxIdxMap[i][idx];
      bit = _arithmeticDecoder->decode(_ctxOccupancy[ctxIdx]);
    }
    numCoded++;

    ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);
    partialOccupancy |= bit << i;
    occupancy |= bit << kOccBitCodingOrder[i];
    if (numCoded == maxOccupied)
      break;
  }

  return occupancy;
}

//-------------------------------------------------------------------------

int
GeometryOctreeDecoder::decodeOccupancyBitwise(
  int neighPattern,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int occupancySkip)
{
  if (neighPattern == 0) {
    return decodeOccupancyNeighZ(
      mappedOccIsPredicted, mappedOccPrediction, mappedOccAdjGt0,
      mappedOccAdjGt1, mappedOccAdjUnocc, occupancySkip);
  }

  return decodeOccupancyNeighNZ(
    neighPattern, mappedOccIsPredicted, mappedOccPrediction, mappedOccAdjGt0,
    mappedOccAdjGt1, mappedOccAdjUnocc, occupancySkip);
}

//-------------------------------------------------------------------------

int
GeometryOctreeDecoder::decodeOccupancyBytewise(int neighPattern)
{
  // code occupancy using the neighbour configuration context
  // with reduction from 64 states to 10 (or 6).
  int neighPatternR1 = _neighPattern64toR1[neighPattern];
  auto& bytewiseCoder = _bytewiseOccupancyCoder[neighPatternR1];
  return bytewiseCoder.decode(_arithmeticDecoder);
}

//-------------------------------------------------------------------------
// decode node occupancy bits
//

uint32_t
GeometryOctreeDecoder::decodeOccupancy(
  int neighPattern,
  int occupancyIsPred,
  int occupancyPred,
  int occupancyAdjGt0,
  int occupancyAdjGt1,
  int occupancyAdjUnocc,
  int occupancySkip)
{
  // decode occupancy pattern
  uint32_t occupancy;
  if (neighPattern == 0) {
    // neighbour empty and only one point => decode index, not pattern
    if (_arithmeticDecoder->decode(_ctxSingleChild)) {
      uint32_t cnt = 0;
      if (!(occupancySkip & 1))
        cnt = _arithmeticDecoder->decode(_ctxEquiProb);
      if (!(occupancySkip & 2))
        cnt |= _arithmeticDecoder->decode(_ctxEquiProb) << 1;
      if (!(occupancySkip & 4))
        cnt |= _arithmeticDecoder->decode(_ctxEquiProb) << 2;
      occupancy = 1 << cnt;
      return occupancy;
    }
  }

  uint32_t mapOccIsP = mapGeometryOccupancy(occupancyIsPred, neighPattern);
  uint32_t mapOccP = mapGeometryOccupancy(occupancyPred, neighPattern);
  uint32_t mapAdjGt0 = mapGeometryOccupancy(occupancyAdjGt0, neighPattern);
  uint32_t mapAdjGt1 = mapGeometryOccupancy(occupancyAdjGt1, neighPattern);
  uint32_t mapAdjUnocc = mapGeometryOccupancy(occupancyAdjUnocc, neighPattern);
  uint32_t mappedOccupancy;

  if (_useBitwiseOccupancyCoder)
    mappedOccupancy = decodeOccupancyBitwise(
      neighPattern, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1, mapAdjUnocc,
      occupancySkip);
  else
    mappedOccupancy = decodeOccupancyBytewise(neighPattern);

  return mapGeometryOccupancyInv(mappedOccupancy, neighPattern);
}

//-------------------------------------------------------------------------
// Decode a position of a point in a given volume.
Vec3<uint32_t>
GeometryOctreeDecoder::decodePointPosition(const Vec3<int>& nodeSizeLog2)
{
  Vec3<uint32_t> delta{};
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2[k] <= 0)
      continue;
    for (int i = nodeSizeLog2[k]; i > 0; i--) {
      delta[k] <<= 1;
      delta[k] |= _arithmeticDecoder->decode(_ctxEquiProb);
    }
  }

  return delta;
}

int
GeometryOctreeDecoder::decodeQpOffset()
{
  int dqp = 0;
  if (!_arithmeticDecoder->decode(_ctxQpOffsetIsZero)) {
    int dqp_sign = _arithmeticDecoder->decode(_ctxQpOffsetSign);
    dqp =
      _arithmeticDecoder->decodeExpGolomb(0, _ctxEquiProb, _ctxQpOffsetAbsEgl)
      + 1;
    dqp = dqp_sign ? dqp : -dqp;
  }
  return dqp;
}

//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).
// Decoded points are written to @outputPoints
// Returns the number of points emitted.

template<class OutputIt>
int
GeometryOctreeDecoder::decodeDirectPosition(
  const Vec3<int>& nodeSizeLog2,
  const PCCOctree3Node& node,
  OutputIt outputPoints)
{
  bool isDirectMode = _arithmeticDecoder->decode(_ctxBlockSkipTh);
  if (!isDirectMode) {
    return 0;
  }

  int numPoints = 1;
  if (_arithmeticDecoder->decode(_ctxNumIdcmPointsEq1))
    numPoints++;

  for (int i = 0; i < numPoints; i++)
    *(outputPoints++) = decodePointPosition(nodeSizeLog2);

  return numPoints;
}

//-------------------------------------------------------------------------
// Helper to inverse quantise positions

Vec3<uint32_t>
invQuantPosition(int qp, int quantBitMask, const Vec3<uint32_t>& pos)
{
  // pos represents the position within the coded tree as follows:
  //   |pppppqqqqqq|00
  //  - p = unquantised bit
  //  - q = quantised bit
  //  - 0 = bits that were not coded (MSBs of q)
  // The reconstruction is:
  //   |ppppp00qqqqqq|

  QuantizerGeom quantizer(qp);
  int shiftBits = (qp - 4) / 6;
  Vec3<uint32_t> recon;
  for (int k = 0; k < 3; k++) {
    int posQuant = pos[k] & (quantBitMask >> shiftBits);
    recon[k] = (pos[k] ^ posQuant) << shiftBits;
    recon[k] |= quantizer.scale(posQuant);
  }

  return recon;
}

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minNodeSizeLog2,
  PCCPointSet3& pointCloud,
  EntropyDecoder* arithmeticDecoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining)
{
  GeometryOctreeDecoder decoder(gps, arithmeticDecoder);

  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  pcc::ringbuf<PCCOctree3Node> fifo(gbh.geom_num_points + 1);

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

  // Size of each child of the current node.
  // NB: the child node sizes may be different due to transitions
  // in the tree type (octree vs quadtree).

  // implicit qtbt for child nodes
  Vec3<int> childSizeLog2 = implicitQtBtDecision(
    nodeSizeLog2, maxNumImplicitQtbtBeforeOt, minSizeImplicitQtbt);

  // prepare parameters for partition and occupancy coding
  int occupancySkip = nonSplitQtBtAxes(nodeSizeLog2, childSizeLog2);
  int atlasShift = 7;

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(0);
  node00.pos = uint32_t(0);
  node00.neighPattern = 0;
  node00.numSiblingsPlus1 = 8;
  node00.siblingOccupancy = 0;
  node00.qp = 4;

  size_t processedPointCount = 0;
  std::vector<uint32_t> values;

  auto fifoCurrLvlEnd = fifo.end();

  // this counter represents fifo.end() - fifoCurrLvlEnd().
  // ie, the number of nodes added to the next level of the tree
  int numNodesNextLvl = 0;

  Vec3<uint32_t> occupancyAtlasOrigin(0xffffffff);
  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2) {
    occupancyAtlas.resize(gps.neighbour_avail_boundary_log2);
    occupancyAtlas.clear();
  }

  int sliceQp = gps.geom_base_qp + gbh.geom_slice_qp_offset;
  int numLvlsUntilQpOffset = -1;
  int posQuantBits = 0;

  if (gbh.geom_octree_qp_offset_enabled_flag)
    numLvlsUntilQpOffset = gbh.geom_octree_qp_offset_depth;
  else if (gps.geom_scaling_enabled_flag) {
    node00.qp = sliceQp;
    // determine the mask of LSBs used in quantisation
    posQuantBits = (1 << (nodeMaxDimLog2 - numLvlsUntilQpOffset)) - 1;
  }

  for (; !fifo.empty(); fifo.pop_front()) {
    if (fifo.begin() == fifoCurrLvlEnd) {
      // transition to the next level
      fifoCurrLvlEnd = fifo.end();

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

      occupancySkip = nonSplitQtBtAxes(nodeSizeLog2, childSizeLog2);

      nodeMaxDimLog2--;
      numNodesNextLvl = 0;
      occupancyAtlasOrigin = 0xffffffff;

      decoder.beginOctreeLevel();

      // allow partial tree encoding using trisoup
      if (nodeMaxDimLog2 == gps.trisoup_node_size_log2)
        break;

      // allow partial tree decoding
      if (nodeMaxDimLog2 == minNodeSizeLog2)
        break;

      numLvlsUntilQpOffset--;
    }

    PCCOctree3Node& node0 = fifo.front();

    if (numLvlsUntilQpOffset == 0)
      node0.qp = decoder.decodeQpOffset() + sliceQp;

    int shiftBits = (node0.qp - 4) / 6;
    int effectiveNodeMaxDimLog2 = nodeMaxDimLog2 - shiftBits;
    auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;
    auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;

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

    uint8_t occupancy = 1;
    if (!isLeafNode(effectiveNodeSizeLog2)) {
      assert(occupancySkip != 7);
      occupancy = decoder.decodeOccupancy(
        node0.neighPattern, occupancyIsPredicted, occupancyPrediction,
        occupancyAdjacencyGt0, occupancyAdjacencyGt1, occupancyAdjacencyUnocc,
        occupancySkip);
    }

    assert(occupancy > 0);

    // update atlas for advanced neighbours
    if (gps.neighbour_avail_boundary_log2) {
      updateGeometryOccupancyAtlasOccChild(
        node0.pos, occupancy, &occupancyAtlas);
    }

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

      // point counts for leaf nodes are coded immediately upon
      // encountering the leaf node.
      if (isLeafNode(effectiveChildSizeLog2)) {
        int numPoints = 1;

        if (!gps.geom_unique_points_flag) {
          numPoints = decoder.decodePositionLeafNumPoints();
        }

        // the final bits from the leaf:
        Vec3<uint32_t> pos{(node0.pos[0] << !(occupancySkip & 4)) + x,
                           (node0.pos[1] << !(occupancySkip & 2)) + y,
                           (node0.pos[2] << !(occupancySkip & 1)) + z};

        pos = invQuantPosition(node0.qp, posQuantBits, pos);
        const Vec3<double> point(pos[0], pos[1], pos[2]);

        for (int i = 0; i < numPoints; ++i)
          pointCloud[processedPointCount++] = point;

        // do not recurse into leaf nodes
        continue;
      }

      // create & enqueue new child.
      fifo.emplace_back();
      auto& child = fifo.back();

      child.qp = node0.qp;
      // only shift position if an occupancy bit was coded for the axis
      child.pos[0] = (node0.pos[0] << !(occupancySkip & 4)) + x;
      child.pos[1] = (node0.pos[1] << !(occupancySkip & 2)) + y;
      child.pos[2] = (node0.pos[2] << !(occupancySkip & 1)) + z;
      child.numSiblingsPlus1 = numOccupied;
      child.siblingOccupancy = occupancy;

      bool idcmEnabled = gps.inferred_direct_coding_mode_enabled_flag;
      if (isDirectModeEligible(
            idcmEnabled, effectiveNodeMaxDimLog2, node0, child)) {
        // todo(df): this should go away when output is integer
        Vec3<uint32_t> points[2]{};
        int numPoints =
          decoder.decodeDirectPosition(effectiveChildSizeLog2, child, points);

        for (int j = 0; j < numPoints; j++) {
          auto& point = points[j];
          for (int k = 0; k < 3; k++) {
            int shift = std::max(0, effectiveChildSizeLog2[k]);
            point[k] += child.pos[k] << shift;
          }

          point = invQuantPosition(node0.qp, posQuantBits, point);
          pointCloud[processedPointCount++] =
            Vec3<double>(point[0], point[1], point[2]);
        }

        if (numPoints > 0) {
          // node fully decoded, do not split: discard child
          fifo.pop_back();

          // NB: no further siblings to decode by definition of IDCM
          assert(child.numSiblingsPlus1 == 1);
          break;
        }
      }

      numNodesNextLvl++;

      if (!gps.neighbour_avail_boundary_log2) {
        updateGeometryNeighState(
          gps.neighbour_context_restriction_flag, fifo.end(), numNodesNextLvl,
          child, i, node0.neighPattern, occupancy);
      }
    }
  }

  // NB: the point cloud needs to be resized if partially decoded
  // OR: if geometry quantisation has changed the number of points
  // todo(df): this breaks the current definition of geom_num_points
  pointCloud.resize(processedPointCount);

  // return partial coding result
  //  - add missing levels to node positions and inverse quantise
  if (nodesRemaining) {
    for (auto& node : fifo) {
      int quantRemovedBits = (node.qp - 4) / 6;
      for (int k = 0; k < 3; k++)
        node.pos[k] <<= nodeSizeLog2[k] - quantRemovedBits;
      node.pos = invQuantPosition(node.qp, posQuantBits, node.pos);
    }
    *nodesRemaining = std::move(fifo);
    return;
  }
}

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  EntropyDecoder* arithmeticDecoder)
{
  decodeGeometryOctree(gps, gbh, 0, pointCloud, arithmeticDecoder, nullptr);
}

//-------------------------------------------------------------------------

void
decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  EntropyDecoder* arithmeticDecoder)
{
  pcc::ringbuf<PCCOctree3Node> nodes;
  decodeGeometryOctree(
    gps, gbh, minGeomNodeSizeLog2, pointCloud, arithmeticDecoder, &nodes);

  if (minGeomNodeSizeLog2 > 0) {
    size_t size =
      pointCloud.removeDuplicatePointInQuantizedPoint(minGeomNodeSizeLog2);

    pointCloud.resize(size + nodes.size());
    size_t processedPointCount = size;

    for (auto node0 : nodes) {
      const Vec3<double> point(node0.pos[0], node0.pos[1], node0.pos[2]);
      pointCloud[processedPointCount++] = point;
    }
  }
}

//============================================================================

}  // namespace pcc
