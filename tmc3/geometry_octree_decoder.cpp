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
    int mappedOccAdjGt1);

  int decodeOccupancyNeighNZ(
    int neighPattern10,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1);

  int decodeOccupancyBitwise(
    int neighPattern,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1);

  int decodeOccupancyBytewise(int neighPattern);

  uint32_t decodeOccupancy(
    int neighPattern,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int occupancyAdjGt0,
    int occupancyAdjGt1);

  PCCVector3<uint32_t> decodePointPosition(int nodeSizeLog2);

  template<class OutputIt>
  int decodeDirectPosition(
    int nodeSizeLog2, const PCCOctree3Node& node, OutputIt outputPoints);

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
  int mappedOccAdjGt1)
{
  static const int8_t bitCodingOrder[8]{1, 7, 5, 3, 2, 4, 6, 0};
  int minOccupied = 2;
  int numOccupiedAcc = 0;
  int occupancy = 0;

  for (int i = 0; i < 8; i++) {
    int bit = 1;
    int bitIsPredicted = (mappedOccIsPredicted >> bitCodingOrder[i]) & 1;
    int bitPrediction = (mappedOccPrediction >> bitCodingOrder[i]) & 1;
    int bitAdjGt0 = (mappedOccAdjGt0 >> bitCodingOrder[i]) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitCodingOrder[i]) & 1;

    int ctxIdxMapIdx =
      3 * (bitAdjGt0 + bitAdjGt1) + bitIsPredicted + bitPrediction;
    auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];

    // NB: There must be at least two occupied child nodes
    //  -- avoid coding the occupancy bit if it is implied.
    if (numOccupiedAcc >= minOccupied + i - 7) {
      int ctxIdx = ctxIdxMap[i][numOccupiedAcc];
      bit = _arithmeticDecoder->decode(_ctxOccupancy[ctxIdx]);
    }
    ctxIdxMap.evolve(bit, &ctxIdxMap[i][numOccupiedAcc]);
    numOccupiedAcc += bit;
    occupancy |= bit << bitCodingOrder[i];
  }

  return occupancy;
}

//---------------------------------------------------------------------------
// decode occupancy bits (neighPattern10 != 0 case)

int
GeometryOctreeDecoder::decodeOccupancyNeighNZ(
  int neighPattern10,
  int mappedOccIsPredicted,
  int mappedOccPrediction,
  int mappedOccAdjGt0,
  int mappedOccAdjGt1)
{
  static const int8_t bitCodingOrder[8]{1, 7, 5, 3, 2, 4, 6, 0};

  int neighPattern7 = kNeighPattern10to7[neighPattern10];
  int neighPattern5 = kNeighPattern7to5[neighPattern7];

  int occupancy = 0;
  int partialOccupancy = 0;

  // NB: it is impossible for pattern to be 0 (handled in Z case).
  // NB: offsets are added since ctxIdxMap is shared between Z and NZ cases.
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

    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    int bit = 1;
    int bitIsPredicted = (mappedOccIsPredicted >> bitCodingOrder[i]) & 1;
    int bitPrediction = (mappedOccPrediction >> bitCodingOrder[i]) & 1;
    int bitAdjGt0 = (mappedOccAdjGt0 >> bitCodingOrder[i]) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitCodingOrder[i]) & 1;

    int ctxIdxMapIdx =
      3 * (bitAdjGt0 + bitAdjGt1) + bitIsPredicted + bitPrediction;
    auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];

    if (i < 7 || partialOccupancy) {
      int ctxIdx = ctxIdxMap[i][idx];
      bit = _arithmeticDecoder->decode(_ctxOccupancy[ctxIdx]);
    }

    ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);
    partialOccupancy |= bit << i;
    occupancy |= bit << bitCodingOrder[i];
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
  int mappedOccAdjGt1)
{
  if (neighPattern == 0) {
    return decodeOccupancyNeighZ(
      mappedOccIsPredicted, mappedOccPrediction, mappedOccAdjGt0,
      mappedOccAdjGt1);
  }

  // code occupancy using the neighbour configuration context
  // with reduction from 64 states to 10 (or 6).
  int neighPatternR1 = _neighPattern64toR1[neighPattern];
  return decodeOccupancyNeighNZ(
    neighPatternR1, mappedOccIsPredicted, mappedOccPrediction, mappedOccAdjGt0,
    mappedOccAdjGt1);
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
  int occupancyAdjGt1)
{
  // decode occupancy pattern
  uint32_t occupancy;
  if (neighPattern == 0) {
    // neighbour empty and only one point => decode index, not pattern
    if (_arithmeticDecoder->decode(_ctxSingleChild)) {
      uint32_t cnt = _arithmeticDecoder->decode(_ctxEquiProb);
      cnt |= _arithmeticDecoder->decode(_ctxEquiProb) << 1;
      cnt |= _arithmeticDecoder->decode(_ctxEquiProb) << 2;
      occupancy = 1 << cnt;
      return occupancy;
    }
  }

  uint32_t mapOccIsP = mapGeometryOccupancy(occupancyIsPred, neighPattern);
  uint32_t mapOccP = mapGeometryOccupancy(occupancyPred, neighPattern);
  uint32_t mapAdjGt0 = mapGeometryOccupancy(occupancyAdjGt0, neighPattern);
  uint32_t mapAdjGt1 = mapGeometryOccupancy(occupancyAdjGt1, neighPattern);
  uint32_t mappedOccupancy;

  if (_useBitwiseOccupancyCoder)
    mappedOccupancy = decodeOccupancyBitwise(
      neighPattern, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1);
  else
    mappedOccupancy = decodeOccupancyBytewise(neighPattern);

  return mapGeometryOccupancyInv(mappedOccupancy, neighPattern);
}

//-------------------------------------------------------------------------
// Decode a position of a point in a given volume.

PCCVector3<uint32_t>
GeometryOctreeDecoder::decodePointPosition(int nodeSizeLog2)
{
  PCCVector3<uint32_t> delta{};
  for (int i = nodeSizeLog2; i > 0; i--) {
    delta <<= 1;
    delta[0] |= _arithmeticDecoder->decode(_ctxEquiProb);
    delta[1] |= _arithmeticDecoder->decode(_ctxEquiProb);
    delta[2] |= _arithmeticDecoder->decode(_ctxEquiProb);
  }

  return delta;
}

//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).
// Decoded points are written to @outputPoints
// Returns the number of points emitted.

template<class OutputIt>
int
GeometryOctreeDecoder::decodeDirectPosition(
  int nodeSizeLog2, const PCCOctree3Node& node, OutputIt outputPoints)
{
  bool isDirectMode = _arithmeticDecoder->decode(_ctxBlockSkipTh);
  if (!isDirectMode) {
    return 0;
  }

  int numPoints = 1;
  if (_arithmeticDecoder->decode(_ctxNumIdcmPointsEq1))
    numPoints++;

  for (int i = 0; i < numPoints; i++) {
    // convert node-relative position to world position
    PCCVector3<uint32_t> pos = node.pos + decodePointPosition(nodeSizeLog2);
    *(outputPoints++) = {double(pos[0]), double(pos[1]), double(pos[2])};
  }

  return numPoints;
}

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  EntropyDecoder* arithmeticDecoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining)
{
  GeometryOctreeDecoder decoder(gps, arithmeticDecoder);

  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  pcc::ringbuf<PCCOctree3Node> fifo(gbh.geom_num_points + 1);

  // the current node dimension (log2)
  int nodeSizeLog2 = gbh.geom_max_node_size_log2;

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(0);
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
  if (gps.neighbour_avail_boundary_log2) {
    occupancyAtlas.resize(gps.neighbour_avail_boundary_log2);
    occupancyAtlas.clear();
  }

  for (; !fifo.empty(); fifo.pop_front()) {
    if (fifo.begin() == fifoCurrLvlEnd) {
      // transition to the next level
      fifoCurrLvlEnd = fifo.end();
      nodeSizeLog2--;
      numNodesNextLvl = 0;
      occupancyAtlasOrigin = 0xffffffff;

      decoder.beginOctreeLevel();

      // allow partial tree encoding using trisoup
      if (nodeSizeLog2 == gps.trisoup_node_size_log2)
        break;
    }

    PCCOctree3Node& node0 = fifo.front();

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

    // decode occupancy pattern
    uint8_t occupancy = decoder.decodeOccupancy(
      node0.neighPattern, occupancyIsPredicted, occupancyPrediction,
      occupancyAdjacencyGt0, occupancyAdjacencyGt1);

    assert(occupancy > 0);

    // update atlas for advanced neighbours
    updateGeometryOccupancyAtlasOccChild(
      node0.pos, nodeSizeLog2, occupancy, &occupancyAtlas);

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

        if (!gps.geom_unique_points_flag) {
          numPoints = decoder.decodePositionLeafNumPoints();
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

      bool idcmEnabled = gps.inferred_direct_coding_mode_enabled_flag;
      if (isDirectModeEligible(idcmEnabled, nodeSizeLog2, node0, child)) {
        int numPoints = decoder.decodeDirectPosition(
          childSizeLog2, child, &pointCloud[processedPointCount]);
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
  decodeGeometryOctree(gps, gbh, pointCloud, arithmeticDecoder, nullptr);
}

//============================================================================

}  // namespace pcc
