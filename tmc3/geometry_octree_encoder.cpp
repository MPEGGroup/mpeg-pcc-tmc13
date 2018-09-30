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
    int mappedOccAdjGt1);

  void encodeOccupancyNeighNZ(
    int neighPattern10,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1);

  void encodeOccupancyBitwise(
    int neighPattern,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1);

  void encodeOccupancyBytewise(int neighPattern, int mappedOccupancy);

  void encodeOccupancy(
    int neighPattern,
    int occupancy,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int occupancyAdjGt0,
    int occupancyAdjGt1);

  void encodePointPosition(int nodeSizeLog2, const PCCVector3<uint32_t>& pos);

  bool encodeDirectPosition(
    int nodeSizeLog2,
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

  // For bitwise occupancy coding
  //   maps 0 + i: no adjacency
  //   maps 3 + i: adjacency > 0
  //   maps 6 + i: adjacency > 1
  //   map i = 0 = not predicted
  //   map i = 1 = predicted unnoccupied
  //   map i = 2 = predicted occupied
  CtxModelOctreeOccupancy _ctxOccupancy;
  CtxMapOctreeOccupancy _ctxIdxMaps[9];

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
  int mappedOccAdjGt1)
{
  static const int8_t bitCodingOrder[8]{1, 7, 5, 3, 2, 4, 6, 0};
  int minOccupied = 2;
  int numOccupiedAcc = 0;

  for (int i = 0; i < 8; i++) {
    int bit = (mappedOccupancy >> bitCodingOrder[i]) & 1;
    int bitIsPredicted = (mappedOccIsPredicted >> bitCodingOrder[i]) & 1;
    int bitPrediction = (mappedOccPrediction >> bitCodingOrder[i]) & 1;
    int bitAdjGt0 = (mappedOccAdjGt0 >> bitCodingOrder[i]) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitCodingOrder[i]) & 1;

    int ctxIdxMapIdx =
      3 * (bitAdjGt0 + bitAdjGt1) + bitIsPredicted + bitPrediction;
    auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];

    int idx = numOccupiedAcc;
    int ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);

    // NB: There must be at least minOccupied child nodes
    //  -- avoid coding the occupancyBit if it is implied.
    if (numOccupiedAcc >= minOccupied + i - 7) {
      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);
    }

    numOccupiedAcc += bit;
  }
}

//-------------------------------------------------------------------------
// encode occupancy bits (neighPattern10 != 0 case)

void
GeometryOctreeEncoder::encodeOccupancyNeighNZ(
  int neighPattern10,
  int mappedOccupancy,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1)
{
  static const int8_t bitCodingOrder[8]{1, 7, 5, 3, 2, 4, 6, 0};

  int neighPattern7 = kNeighPattern10to7[neighPattern10];
  int neighPattern5 = kNeighPattern7to5[neighPattern7];

  uint32_t partialOccupancy = 0;

  // NB: it is impossible for pattern to be 0 (handled in Z case).
  for (int i = 0; i < 8; i++) {
    int idx;
    if (i < 6) {
      idx = ((neighPattern10 - 1) << i) + partialOccupancy + i + 1;
    } else if (i == 6) {
      idx = ((neighPattern7 - 1) << i) + partialOccupancy + i + 1;
    } else if (i == 7) {
      idx = ((neighPattern5 - 1) << i) + partialOccupancy + i + 1;
    } else {
      // work around clang -Wsometimes-uninitialized fault
      break;
    }

    int bit = (mappedOccupancy >> bitCodingOrder[i]) & 1;
    int bitIsPredicted = (mappedOccIsPredicted >> bitCodingOrder[i]) & 1;
    int bitPrediction = (mappedOccPrediction >> bitCodingOrder[i]) & 1;
    int bitAdjGt0 = (mappedOccAdjGt0 >> bitCodingOrder[i]) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitCodingOrder[i]) & 1;

    int ctxIdxMapIdx =
      3 * (bitAdjGt0 + bitAdjGt1) + bitIsPredicted + bitPrediction;
    auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];

    int ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);

    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    if (i < 7 || partialOccupancy)
      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);

    partialOccupancy |= bit << i;
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
  int mappedOccAdjGt1)
{
  if (neighPattern == 0) {
    encodeOccupancyNeighZ(
      mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
      mappedOccAdjGt0, mappedOccAdjGt1);
    return;
  }

  // code occupancy using the neighbour configuration context
  // with reduction from 64 states to 10 (or 6).
  int neighPatternR1 = _neighPattern64toR1[neighPattern];
  encodeOccupancyNeighNZ(
    neighPatternR1, mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
    mappedOccAdjGt0, mappedOccAdjGt1);
}

//-------------------------------------------------------------------------

void
GeometryOctreeEncoder::encodeOccupancyBytewise(
  int mappedOccupancy, int neighPattern)
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
  int occupancyAdjGt1)
{
  if (neighPattern == 0) {
    bool singleChild = !popcntGt1(occupancy);
    _arithmeticEncoder->encode(singleChild, _ctxSingleChild);

    if (singleChild) {
      // no siblings => encode index = (z,y,x) not 8bit pattern
      _arithmeticEncoder->encode(!!(occupancy & 0xaa), _ctxEquiProb);  // z
      _arithmeticEncoder->encode(!!(occupancy & 0xcc), _ctxEquiProb);  // y
      _arithmeticEncoder->encode(!!(occupancy & 0xf0), _ctxEquiProb);  // x
      return;
    }
  }

  uint32_t mapOcc = mapGeometryOccupancy(occupancy, neighPattern);
  uint32_t mapOccIsP = mapGeometryOccupancy(occupancyIsPred, neighPattern);
  uint32_t mapOccP = mapGeometryOccupancy(occupancyPred, neighPattern);
  uint32_t mapAdjGt0 = mapGeometryOccupancy(occupancyAdjGt0, neighPattern);
  uint32_t mapAdjGt1 = mapGeometryOccupancy(occupancyAdjGt1, neighPattern);

  if (_useBitwiseOccupancyCoder)
    encodeOccupancyBitwise(
      neighPattern, mapOcc, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1);
  else
    encodeOccupancyBytewise(neighPattern, mapOcc);
}

//-------------------------------------------------------------------------
// Encode a position of a point in a given volume.

void
GeometryOctreeEncoder::encodePointPosition(
  int nodeSizeLog2, const PCCVector3<uint32_t>& pos)
{
  for (int mask = 1 << (nodeSizeLog2 - 1); mask; mask >>= 1) {
    _arithmeticEncoder->encode(!!(pos[0] & mask), _ctxEquiProb);
    _arithmeticEncoder->encode(!!(pos[1] & mask), _ctxEquiProb);
    _arithmeticEncoder->encode(!!(pos[2] & mask), _ctxEquiProb);
  }
}

//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).

bool
GeometryOctreeEncoder::encodeDirectPosition(
  int nodeSizeLog2, const PCCOctree3Node& node, const PCCPointSet3& pointCloud)
{
  int numPoints = node.end - node.start;
  if (numPoints > MAX_NUM_DM_LEAF_POINTS) {
    _arithmeticEncoder->encode(0, _ctxBlockSkipTh);
    return false;
  }

  _arithmeticEncoder->encode(1, _ctxBlockSkipTh);
  _arithmeticEncoder->encode(numPoints > 1, _ctxNumIdcmPointsEq1);

  for (auto idx = node.start; idx < node.end; idx++) {
    // determine the point position relative to box edge
    encodePointPosition(
      nodeSizeLog2,
      PCCVector3<uint32_t>{int(pointCloud[idx][0]) - node.pos[0],
                           int(pointCloud[idx][1]) - node.pos[1],
                           int(pointCloud[idx][2]) - node.pos[2]});
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

  // map of pointCloud idx to DM idx, used to reorder the points
  // after coding.
  std::vector<int> pointIdxToDmIdx(int(pointCloud.getPointCount()), -1);
  int nextDmIdx = 0;

  size_t processedPointCount = 0;
  std::vector<uint32_t> values;

  auto fifoCurrLvlEnd = fifo.end();

  // the initial node size is the root node's
  int nodeSizeLog2 = gbh.geom_max_node_size_log2;

  // this counter represents fifo.end() - fifoCurrLvlEnd().
  // ie, the number of nodes added to the next level of the tree
  int numNodesNextLvl = 0;

  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2) {
    occupancyAtlas.resize(gps.neighbour_avail_boundary_log2);
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

      encoder.beginOctreeLevel();

      // allow partial tree encoding using trisoup
      if (nodeSizeLog2 == gps.trisoup_node_size_log2)
        break;
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

    int occupancyAdjacencyGt0 = 0;
    int occupancyAdjacencyGt1 = 0;

    if (gps.neighbour_avail_boundary_log2) {
      updateGeometryOccupancyAtlas(
        node0.pos, nodeSizeLog2, fifo, fifoCurrLvlEnd, &occupancyAtlas,
        &occupancyAtlasOrigin);

      GeometryNeighPattern gnp =
        makeGeometryNeighPattern(node0.pos, nodeSizeLog2, occupancyAtlas);
      node0.neighPattern = gnp.neighPattern;
      occupancyAdjacencyGt0 = gnp.adjacencyGt0;
      occupancyAdjacencyGt1 = gnp.adjacencyGt1;
    }

    int occupancyIsPredicted = 0;
    int occupancyPrediction = 0;

    // generate intra prediction
    if (nodeSizeLog2 < gps.intra_pred_max_node_size_log2) {
      predictGeometryOccupancyIntra(
        occupancyAtlas, node0.pos, nodeSizeLog2, &occupancyIsPredicted,
        &occupancyPrediction);
    }

    // update atlas for advanced neighbours
    updateGeometryOccupancyAtlasOccChild(
      node0.pos, nodeSizeLog2, occupancy, &occupancyAtlas);

    // encode child occupancy map
    assert(occupancy > 0);
    encoder.encodeOccupancy(
      node0.neighPattern, occupancy, occupancyIsPredicted, occupancyPrediction,
      occupancyAdjacencyGt0, occupancyAdjacencyGt1);

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
        if (gps.geom_unique_points_flag) {
          assert(childCounts[i] == 1);
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

      child.pos[0] = node0.pos[0] + (x << childSizeLog2);
      child.pos[1] = node0.pos[1] + (y << childSizeLog2);
      child.pos[2] = node0.pos[2] + (z << childSizeLog2);

      child.start = childPointsStartIdx;
      childPointsStartIdx += childCounts[i];
      child.end = childPointsStartIdx;
      child.numSiblingsPlus1 = numSiblings;
      child.siblingOccupancy = occupancy;

      bool idcmEnabled = gps.inferred_direct_coding_mode_enabled_flag;
      if (isDirectModeEligible(idcmEnabled, nodeSizeLog2, node0, child)) {
        bool directModeUsed =
          encoder.encodeDirectPosition(childSizeLog2, child, pointCloud);

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
      if (!gps.neighbour_avail_boundary_log2) {
        updateGeometryNeighState(
          gps.neighbour_context_restriction_flag, fifo.end(), numNodesNextLvl,
          childSizeLog2, child, i, node0.neighPattern, occupancy);
      }
    }
  }

  // return partial coding result
  if (nodesRemaining) {
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
    }

    pointCloud2[dstIdx] = pointCloud[i];
    if (pointCloud.hasColors())
      pointCloud2.setColor(dstIdx, pointCloud.getColor(i));
    if (pointCloud.hasReflectances())
      pointCloud2.setReflectance(dstIdx, pointCloud.getReflectance(i));
  }

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
