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
  int encodePlanarMode(
    PCCOctree3Node& node0,
    int plane,
    int posxyz,
    int dist,
    int neighb,
    int& h,
    int planeId);

  void determinePlanarMode(
    int index,
    const int kNumPlanarPlanes,
    PCCOctree3Node& child,
    uint8_t planarMode,
    uint8_t planePosBits,
    int* plane1,
    int* plane2,
    int* plane3,
    int coord1,
    int coord2,
    int coord3,
    uint8_t neighPattern,
    int pos,
    int planarProb[3],
    int planarRate[3]);

  void determinePlanarMode(
    PCCPointSet3& pointCloud,
    bool planarEligible[3],
    const Vec3<int>& childSizeLog2,
    const int kNumPlanarPlanes,
    PCCOctree3Node& child,
    int* planes[9],
    uint8_t neighPattern,
    int x,
    int y,
    int z,
    int planarProb[3],
    int planarRate[3]);

  void encodeOccupancyNeighZ(
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int mappedPlanarMaskX,
    int mappedFixedMaskX0,
    bool planarPossibleX,
    int mappedPlanarMaskY,
    int mappedFixedMaskY0,
    bool planarPossibleY,
    int mappedPlanarMaskZ,
    int mappedFixedMaskZ0,
    bool planarPossibleZ);

  void encodeOccupancyNeighNZ(
    int neighPattern,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int mappedPlanarMaskX,
    int mappedFixedMaskX0,
    bool planarPossibleX,
    int mappedPlanarMaskY,
    int mappedFixedMaskY0,
    bool planarPossibleY,
    int mappedPlanarMaskZ,
    int mappedFixedMaskZ0,
    bool planarPossibleZ);

  void encodeOccupancyBitwise(
    int neighPattern,
    int mappedOccupancy,
    int mappedOccIsPredicted,
    int mappedOccPrediction,
    int mappedOccAdjGt0,
    int mappedOccAdjGt1,
    int mappedOccAdjUnocc,
    int mappedPlanarMaskX,
    int mappedFixedMaskX0,
    bool planarPossibleX,
    int mappedPlanarMaskY,
    int mappedFixedMaskY0,
    bool planarPossibleY,
    int mappedPlanarMaskZ,
    int mappedFixedMaskZ0,
    bool planarPossibleZ);

  void encodeOccupancyBytewise(int neighPattern, int mappedOccupancy);

  void encodeOccupancy(
    int neighPattern,
    int occupancy,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int occupancyAdjGt0,
    int occupancyAdjGt1,
    int occupancyAdjUnocc,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ);

  void encodePointPosition(
    const Vec3<int>& nodeSizeLog2,
    const Vec3<int32_t>& pos,
    uint8_t planarMode);

  void encodeQpOffset(int dqp);

  bool encodeDirectPosition(
    bool geom_unique_points_flag,
    const Vec3<int>& nodeSizeLog2,
    int shiftBits,
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
  AdaptiveBitModel _ctxSingleIdcmDupPoint;
  AdaptiveBitModel _ctxPointCountPerBlock;
  AdaptiveBitModel _ctxBlockSkipTh;
  AdaptiveBitModel _ctxNumIdcmPointsGt1;
  AdaptiveBitModel _ctxSameZ;

  AdaptiveBitModel _ctxQpOffsetIsZero;
  AdaptiveBitModel _ctxQpOffsetSign;
  AdaptiveBitModel _ctxQpOffsetAbsEgl;

  // for planar mode xyz
  AdaptiveBitModel _ctxPlanarMode[3][2][2];
  AdaptiveBitModel _ctxPlanarPlaneLastIndex[3][2][6][2];
  AdaptiveBitModel _ctxPlanarPlaneLastIndexZ[3];

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

//============================================================================

int
GeometryOctreeEncoder::encodePlanarMode(
  PCCOctree3Node& node0,
  int plane,
  int posxyz,
  int dist,
  int neighb,
  int& h,
  int planeId)
{
  const int mask0 = (1 << planeId);
  const int mask1[3] = {6, 5, 3};

  bool isPlanar = node0.planarMode & mask0;
  int planeBit = (node0.planePosBits & mask0) == 0 ? 0 : 1;

  int discreteDist = (dist <= 2 ? 0 : 1);
  _arithmeticEncoder->encode(
    isPlanar, _ctxPlanarMode[planeId][neighb][discreteDist]);

  if (!isPlanar) {
    node0.planarPossible &= mask1[planeId];
    return -1;
  }

  // encode the plane index
  if (plane < 0) {
    _arithmeticEncoder->encode(planeBit, _ctxPlanarPlaneLastIndexZ[planeId]);
    h = approxSymbolProbability(planeBit, _ctxPlanarPlaneLastIndexZ[planeId]);
  } else {
    discreteDist += (dist <= 16 ? 0 : 1);
    int lastIndexPlane2d = plane + (discreteDist << 1);
    _arithmeticEncoder->encode(
      planeBit,
      _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d][posxyz]);
    h = approxSymbolProbability(
      planeBit,
      _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d][posxyz]);
  }
  return planeBit;
}

//============================================================================
// determine Planar mode for one direction

void
GeometryOctreeEncoder::determinePlanarMode(
  int planeId,
  const int kNumPlanarPlanes,
  PCCOctree3Node& child,
  uint8_t planarMode,
  uint8_t planePosBits,
  int* plane1,
  int* plane2,
  int* plane3,
  int coord1,
  int coord2,
  int coord3,
  uint8_t neighPattern,
  int pos,
  int planarProb[3],
  int planarRate[3])
{
  const int kPlanarChildThreshold = 63;
  const int kAdjNeighIdxFromPlanePos[3][2] = {1, 0, 2, 3, 4, 5};
  const int planeSelector = 1 << planeId;

  child.planarMode |= planarMode & planeSelector;
  child.planePosBits |= planePosBits & planeSelector;

  const int shift = coord3 * kNumPlanarPlanes;
  int* localPlane1 = plane1 + shift;
  int* localPlane2 = plane2 + shift;
  int* localPlane3 = plane3 + shift;

  int minDist = std::abs(coord1 - localPlane1[kNumPlanarPlanes - 1])
    + std::abs(coord2 - localPlane2[kNumPlanarPlanes - 1]);
  int idxMinDist = kNumPlanarPlanes - 1;

  for (int idxP = 0; idxP < kNumPlanarPlanes - 1; idxP++) {
    int dist0 = std::abs(coord1 - localPlane1[idxP])
      + std::abs(coord2 - localPlane2[idxP]);
    if (dist0 < minDist) {
      idxMinDist = idxP;
      minDist = dist0;
    }
  }

  // push closest point front
  localPlane1[kNumPlanarPlanes - 1] = localPlane1[idxMinDist];
  localPlane2[kNumPlanarPlanes - 1] = localPlane2[idxMinDist];
  localPlane3[kNumPlanarPlanes - 1] = localPlane3[idxMinDist];

  int adjNeigh = (neighPattern >> kAdjNeighIdxFromPlanePos[planeId][pos]) & 1;
  int planeBit = encodePlanarMode(
    child, localPlane3[kNumPlanarPlanes - 1], pos, minDist, adjNeigh,
    planarProb[planeId], planeId);

  bool isPlanar = (child.planarMode & planeSelector)
    && planarProb[planeId] > kPlanarChildThreshold;

  planarRate[planeId] =
    (255 * planarRate[planeId] + (isPlanar ? 256 * 8 : 0) + 128) >> 8;

  for (int idxP = 0; idxP < kNumPlanarPlanes - 1;
       idxP++, localPlane1++, localPlane2++, localPlane3++) {
    *localPlane1 = *(localPlane1 + 1);
    *localPlane2 = *(localPlane2 + 1);
    *localPlane3 = *(localPlane3 + 1);
  }

  *localPlane1 = coord1;
  *localPlane2 = coord2;
  *localPlane3 = planeBit;
}

//============================================================================
// determine Planar mode for all directions

void
GeometryOctreeEncoder::determinePlanarMode(
  PCCPointSet3& pointCloud,
  bool planarEligible[3],
  const Vec3<int>& childSizeLog2,
  const int kNumPlanarPlanes,
  PCCOctree3Node& child,
  int* planes[9],
  uint8_t neighPattern,
  int x,
  int y,
  int z,
  int planarProb[3],
  int planarRate[3])
{
  // planarity
  uint8_t planarMode, planePosBits;
  isPlanarNode(
    pointCloud, child, childSizeLog2 - 1, planarMode, planePosBits,
    planarEligible);

  int xx = child.pos[0];
  int yy = child.pos[1];
  int zz = child.pos[2];

  // planar x
  if (planarEligible[0]) {
    determinePlanarMode(
      0, kNumPlanarPlanes, child, planarMode, planePosBits, planes[1],
      planes[2], planes[0], yy, zz, xx, neighPattern, x, planarProb,
      planarRate);
  }
  // planar y
  if (planarEligible[1]) {
    determinePlanarMode(
      1, kNumPlanarPlanes, child, planarMode, planePosBits, planes[3],
      planes[5], planes[4], xx, zz, yy, neighPattern, y, planarProb,
      planarRate);
  }
  // planar z
  if (planarEligible[2]) {
    determinePlanarMode(
      2, kNumPlanarPlanes, child, planarMode, planePosBits, planes[6],
      planes[7], planes[8], xx, yy, zz, neighPattern, z, planarProb,
      planarRate);
  }
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
  int mappedPlanarMaskX,
  int mappedFixedMaskX0,
  bool planarPossibleX,
  int mappedPlanarMaskY,
  int mappedFixedMaskY0,
  bool planarPossibleY,
  int mappedPlanarMaskZ,
  int mappedFixedMaskZ0,
  bool planarPossibleZ)
{
  int numOccupiedAcc = 0;

  int maxPerPlaneX = 4 - (mappedPlanarMaskX ? 2 : 1);
  int maxPerPlaneY = 4 - (mappedPlanarMaskY ? 2 : 1);
  int maxPerPlaneZ = 4 - (mappedPlanarMaskZ ? 2 : 1);
  bool sure_planarityX = mappedPlanarMaskX || !planarPossibleX;
  bool sure_planarityY = mappedPlanarMaskY || !planarPossibleY;
  bool sure_planarityZ = mappedPlanarMaskZ || !planarPossibleZ;

  int maskedOccupancy =
    mappedPlanarMaskX | mappedPlanarMaskY | mappedPlanarMaskZ;

  int coded0X[2] = {0, 0};
  int coded0Y[2] = {0, 0};
  int coded0Z[2] = {0, 0};
  if (maskedOccupancy) {
    for (int i = 0; i < 8; i++) {
      if ((maskedOccupancy >> i) & 1) {
        coded0X[(mappedFixedMaskX0 >> i) & 1]++;
        coded0Y[(mappedFixedMaskY0 >> i) & 1]++;
        coded0Z[(mappedFixedMaskZ0 >> i) & 1]++;
      }
    }
  }

  for (int i = 0; i < 8; i++) {
    int bitIdx = kOccBitCodingOrder[i];
    if ((maskedOccupancy >> bitIdx) & 1)
      continue;

    int bitAdjGt0 = (mappedOccAdjGt0 >> bitIdx) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitIdx) & 1;
    int bitAdjUnocc = (mappedOccAdjUnocc >> bitIdx) & 1;
    int numAdj = bitAdjGt0 + bitAdjGt1;
    int idxAdj = bitAdjUnocc + 2 * numAdj;
    if (i > 4) {
      static const int8_t kCtxIdxAdjReduc567[6] = {0, 0, 1, 2, 3, 3};
      idxAdj = kCtxIdxAdjReduc567[idxAdj];
    }

    int ctxIdxMapIdx = idxAdj;
    if (!maskedOccupancy) {
      int bitIsPredicted = (mappedOccIsPredicted >> bitIdx) & 1;
      int bitPrediction = (mappedOccPrediction >> bitIdx) & 1;
      ctxIdxMapIdx = 3 * idxAdj + bitIsPredicted + bitPrediction;
    }

    // NB: There must be at least minOccupied child nodes
    //  -- avoid coding the occupancyBit if it is implied.
    int mask0X = (mappedFixedMaskX0 >> bitIdx) & 1;
    bool bitIsOneX = (sure_planarityX && coded0X[mask0X] >= maxPerPlaneX)
      || (coded0X[0] + coded0X[1] >= 6);

    int mask0Y = (mappedFixedMaskY0 >> bitIdx) & 1;
    bool bitIsOneY = (sure_planarityY && coded0Y[mask0Y] >= maxPerPlaneY)
      || (coded0Y[0] + coded0Y[1] >= 6);

    int mask0Z = (mappedFixedMaskZ0 >> bitIdx) & 1;
    bool bitIsOneZ = (sure_planarityZ && coded0Z[mask0Z] >= maxPerPlaneZ)
      || (coded0Z[0] + coded0Z[1] >= 6);

    // masking for planar is here
    int bit = (mappedOccupancy >> bitIdx) & 1;
    if (!(bitIsOneX || bitIsOneY || bitIsOneZ)) {
      int ctxIdx;
      auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];
      ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][numOccupiedAcc]);

      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);

      if (!bit) {
        coded0X[mask0X]++;
        coded0Y[mask0Y]++;
        coded0Z[mask0Z]++;
      }
    }

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
  int mappedPlanarMaskX,
  int mappedFixedMaskX0,
  bool planarPossibleX,
  int mappedPlanarMaskY,
  int mappedFixedMaskY0,
  bool planarPossibleY,
  int mappedPlanarMaskZ,
  int mappedFixedMaskZ0,
  bool planarPossibleZ)
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

  bool sure_planarityX = mappedPlanarMaskX || !planarPossibleX;
  bool sure_planarityY = mappedPlanarMaskY || !planarPossibleY;
  bool sure_planarityZ = mappedPlanarMaskZ || !planarPossibleZ;

  int maskedOccupancy =
    mappedPlanarMaskX | mappedPlanarMaskY | mappedPlanarMaskZ;

  int coded0X[2] = {0, 0};
  int coded0Y[2] = {0, 0};
  int coded0Z[2] = {0, 0};
  if (maskedOccupancy) {
    for (int i = 0; i < 8; i++) {
      if ((maskedOccupancy >> i) & 1) {
        coded0X[(mappedFixedMaskX0 >> i) & 1]++;
        coded0Y[(mappedFixedMaskY0 >> i) & 1]++;
        coded0Z[(mappedFixedMaskZ0 >> i) & 1]++;
      }
    }
  }

  // NB: it is impossible for pattern to be 0 (handled in Z case).
  for (int i = 0; i < 8; i++) {
    int bitIdx = kOccBitCodingOrder[i];
    if ((maskedOccupancy >> bitIdx) & 1)
      continue;

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

    int bitAdjGt0 = (mappedOccAdjGt0 >> bitIdx) & 1;
    int bitAdjGt1 = (mappedOccAdjGt1 >> bitIdx) & 1;
    int bitAdjUnocc = (mappedOccAdjUnocc >> bitIdx) & 1;

    int numAdj = bitAdjGt0 + bitAdjGt1;
    int idxAdj = bitAdjUnocc + 2 * numAdj;
    if (i > 4) {
      static const int8_t kCtxIdxAdjReduc567[6] = {0, 0, 1, 2, 3, 3};
      idxAdj = kCtxIdxAdjReduc567[idxAdj];
    }

    int ctxIdxMapIdx = bitAdjGt0 + bitAdjGt1;
    if (!maskedOccupancy) {  // planar
      int bitIsPredicted = (mappedOccIsPredicted >> bitIdx) & 1;
      int bitPrediction = (mappedOccPrediction >> bitIdx) & 1;
      ctxIdxMapIdx = 3 * idxAdj + bitIsPredicted + bitPrediction;
    }

    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    // masking for planar is here
    int mask0X = (mappedFixedMaskX0 >> bitIdx) & 1;
    bool bitIsOneX = (sure_planarityX && coded0X[mask0X] >= 3)
      || (coded0X[0] + coded0X[1] >= 7);

    int mask0Y = (mappedFixedMaskY0 >> bitIdx) & 1;
    bool bitIsOneY = (sure_planarityY && coded0Y[mask0Y] >= 3)
      || (coded0Y[0] + coded0Y[1] >= 7);

    int mask0Z = (mappedFixedMaskZ0 >> bitIdx) & 1;
    bool bitIsOneZ = (sure_planarityZ && coded0Z[mask0Z] >= 3)
      || (coded0Z[0] + coded0Z[1] >= 7);

    int bit = (mappedOccupancy >> bitIdx) & 1;
    if (!(bitIsOneX || bitIsOneY || bitIsOneZ)) {
      auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];
      int ctxIdx = ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);
      _arithmeticEncoder->encode(bit, _ctxOccupancy[ctxIdx]);

      if (!bit) {
        coded0X[mask0X]++;
        coded0Y[mask0Y]++;
        coded0Z[mask0Z]++;
      }
    }

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
  int mappedOccAdjGt1,
  int mappedOccAdjUnocc,
  int mappedPlanarMaskX,
  int mappedFixedMaskX0,
  bool planarPossibleX,
  int mappedPlanarMaskY,
  int mappedFixedMaskY0,
  bool planarPossibleY,
  int mappedPlanarMaskZ,
  int mappedFixedMaskZ0,
  bool planarPossibleZ)

{
  if (neighPattern == 0) {
    encodeOccupancyNeighZ(
      mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
      mappedOccAdjGt0, mappedOccAdjGt1, mappedOccAdjUnocc, mappedPlanarMaskX,
      mappedFixedMaskX0, planarPossibleX, mappedPlanarMaskY, mappedFixedMaskY0,
      planarPossibleY, mappedPlanarMaskZ, mappedFixedMaskZ0, planarPossibleZ);
    return;
  }

  encodeOccupancyNeighNZ(
    neighPattern, mappedOccupancy, mappedOccIsPredicted, mappedOccPrediction,
    mappedOccAdjGt0, mappedOccAdjGt1, mappedOccAdjUnocc, mappedPlanarMaskX,
    mappedFixedMaskX0, planarPossibleX, mappedPlanarMaskY, mappedFixedMaskY0,
    planarPossibleY, mappedPlanarMaskZ, mappedFixedMaskZ0, planarPossibleZ);
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
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  bool planarPossibleX,
  bool planarPossibleY,
  bool planarPossibleZ)
{
  // 3 planars => single child and we know its position
  if (planarMaskX && planarMaskY && planarMaskZ)
    return;

  if (neighPattern == 0) {
    bool singleChild = !popcntGt1(occupancy);
    if (planarPossibleX && planarPossibleY && planarPossibleZ) {
      _arithmeticEncoder->encode(singleChild, _ctxSingleChild);
    }

    if (singleChild) {
      // no siblings => encode index = (z,y,x) not 8bit pattern
      // if mask is not zero, then planar, then child z known from plane index
      if (!planarMaskZ)
        _arithmeticEncoder->encode(!!(occupancy & 0xaa), _ctxEquiProb);  // z

      if (!planarMaskY)
        _arithmeticEncoder->encode(!!(occupancy & 0xcc), _ctxEquiProb);  // y

      if (!planarMaskX)
        _arithmeticEncoder->encode(!!(occupancy & 0xf0), _ctxEquiProb);  // x

      return;
    }
  }

  // at least two child nodes occupied and two planars => we know the occupancy
  if (neighPattern == 0) {
    if (planarMaskX && planarMaskY)
      return;
    if (planarMaskY && planarMaskZ)
      return;
    if (planarMaskX && planarMaskZ)
      return;
  }

  uint32_t mapOcc = mapGeometryOccupancy(occupancy, neighPattern);
  uint32_t mapOccIsP = mapGeometryOccupancy(occupancyIsPred, neighPattern);
  uint32_t mapOccP = mapGeometryOccupancy(occupancyPred, neighPattern);
  uint32_t mapAdjGt0 = mapGeometryOccupancy(occupancyAdjGt0, neighPattern);
  uint32_t mapAdjGt1 = mapGeometryOccupancy(occupancyAdjGt1, neighPattern);
  uint32_t mapAdjUnocc = mapGeometryOccupancy(occupancyAdjUnocc, neighPattern);

  uint32_t mapPlanarMaskX = mapGeometryOccupancy(planarMaskX, neighPattern);
  uint32_t mapPlanarMaskY = mapGeometryOccupancy(planarMaskY, neighPattern);
  uint32_t mapPlanarMaskZ = mapGeometryOccupancy(planarMaskZ, neighPattern);

  uint32_t mapFixedMaskX0 = mapGeometryOccupancy(0xf0, neighPattern);
  uint32_t mapFixedMaskY0 = mapGeometryOccupancy(0xcc, neighPattern);
  uint32_t mapFixedMaskZ0 = mapGeometryOccupancy(0xaa, neighPattern);

  if (_useBitwiseOccupancyCoder)
    encodeOccupancyBitwise(
      neighPattern, mapOcc, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1,
      mapAdjUnocc, mapPlanarMaskX, mapFixedMaskX0, planarPossibleX,
      mapPlanarMaskY, mapFixedMaskY0, planarPossibleY, mapPlanarMaskZ,
      mapFixedMaskZ0, planarPossibleZ);
  else
    encodeOccupancyBytewise(neighPattern, mapOcc);
}

//-------------------------------------------------------------------------
// Encode a position of a point in a given volume.
void
GeometryOctreeEncoder::encodePointPosition(
  const Vec3<int>& nodeSizeLog2, const Vec3<int32_t>& pos, uint8_t planarMode)
{
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2[k] <= 0)
      continue;
    int mask = 1 << (nodeSizeLog2[k] - 1);

    if (!(planarMode & (1 << k)))
      _arithmeticEncoder->encode(!!(pos[k] & mask), _ctxEquiProb);
    mask >>= 1;

    for (; mask; mask >>= 1) {
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
    int32_t clipMax = ((1 << nodeSizeLog2[k]) >> qpShift) - 1;

    for (int i = node.start; i < node.end; i++) {
      int32_t pos = int32_t(pointCloud[i][k]);
      int32_t quantPos = quantizer.quantize(pos & quantBitsMask);
      quantPos = PCCClip(quantPos, 0, clipMax);

      // NB: this representation is: |ppppppqqq00|, which, except for
      // the zero padding, is the same as the decoder.
      pointCloud[i][k] = (pos & ~quantBitsMask) | (quantPos << qpShift);
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
      int pos = pointCloud[i][k];
      int lowPart = (pos & quantBitsMask) >> qpShift;
      int lowPartScaled = PCCClip(quantizer.scale(lowPart), 0, quantBitsMask);
      int highPartScaled = pos & ~quantBitsMask;
      pointCloud[i][k] = highPartScaled | lowPartScaled;
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

  std::set<Vec3<int32_t>> uniquePointsSet;
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
  bool geom_unique_points_flag,
  const Vec3<int>& nodeSizeLog2,
  int shiftBits,
  const PCCOctree3Node& node,
  const PCCPointSet3& pointCloud)
{
  int numPoints = node.end - node.start;
  // Check for duplicated points only if there are less than 10.
  // NB: this limit is rather arbitrary
  if (numPoints > 10) {
    _arithmeticEncoder->encode(0, _ctxBlockSkipTh);
    return false;
  }

  bool allPointsAreEqual = numPoints > 1 && !geom_unique_points_flag;
  for (auto idx = node.start + 1; allPointsAreEqual && idx < node.end; idx++) {
    allPointsAreEqual &= pointCloud[node.start] == pointCloud[idx];
  }

  if (!allPointsAreEqual) {
    if (numPoints > MAX_NUM_DM_LEAF_POINTS) {
      _arithmeticEncoder->encode(0, _ctxBlockSkipTh);
      return false;
    }
    _arithmeticEncoder->encode(1, _ctxBlockSkipTh);
    _arithmeticEncoder->encode(numPoints > 1, _ctxNumIdcmPointsGt1);
    if (!geom_unique_points_flag && numPoints == 1)
      _arithmeticEncoder->encode(numPoints == 1, _ctxSinglePointPerBlock);
  } else {
    _arithmeticEncoder->encode(1, _ctxBlockSkipTh);
    _arithmeticEncoder->encode(0, _ctxNumIdcmPointsGt1);
    _arithmeticEncoder->encode(0, _ctxSinglePointPerBlock);
    _arithmeticEncoder->encode(numPoints == 2, _ctxSingleIdcmDupPoint);
    if (numPoints > 2)
      _arithmeticEncoder->encodeExpGolomb(
        numPoints - 3, 0, _ctxEquiProb, _ctxPointCountPerBlock);

    // only one actual psoition to code
    numPoints = 1;
  }

  for (auto idx = node.start; idx < node.start + numPoints; idx++) {
    encodePointPosition(
      nodeSizeLog2, pointCloud[idx] >> shiftBits, node.planarMode);
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
  node00.pos = int32_t(0);
  node00.neighPattern = 0;
  node00.numSiblingsPlus1 = 8;
  node00.siblingOccupancy = 0;
  node00.qp = 4;
  node00.planarMode = 0;

  // map of pointCloud idx to DM idx, used to reorder the points
  // after coding.
  std::vector<int> pointIdxToDmIdx(int(pointCloud.getPointCount()), -1);
  int nextDmIdx = 0;

  // represents the largest dimension of the current node
  // NB: this is equal to the total depth of the tree
  int nodeMaxDimLog2 = gbh.geomMaxNodeSizeLog2(gps);

  // generate the list of the node size for each level in the tree
  auto lvlNodeSizeLog2 = mkQtBtNodeSizeList(gps, gbh);
  auto nodeSizeLog2 = lvlNodeSizeLog2[0];

  // planar mode initialazation
  const int idcmThreshold = gps.geom_planar_idcm_threshold * 127 * 127;
  const int kNumPlanarPlanes = 4;
  std::vector<int> planes3x3;
  int* planes[9];

  if (gps.geom_planar_mode_enabled_flag) {
    const int maxPlaneSize = kNumPlanarPlanes << nodeMaxDimLog2;
    planes3x3.resize(maxPlaneSize * 9);
  }

  int planarRate[3] = {128 * 8, 128 * 8, 128 * 8};
  int localDensity = 1024 * 4;
  const int planarRateThreshold[3] = {gps.geom_planar_threshold0 << 4,
                                      gps.geom_planar_threshold1 << 4,
                                      gps.geom_planar_threshold2 << 4};

  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2) {
    occupancyAtlas.resize(gps.neighbour_avail_boundary_log2);
    occupancyAtlas.clear();
  }

  // the node size where quantisation is performed
  Vec3<int> quantNodeSizeLog2 = 0;
  int sliceQp = gps.geom_base_qp + gbh.geom_slice_qp_offset;
  int numLvlsUntilQuantization = 0;
  if (gps.geom_scaling_enabled_flag)
    numLvlsUntilQuantization = gbh.geom_octree_qp_offset_depth + 1;

  // the termination depth of the octree phase
  int maxDepth = nodeMaxDimLog2 - gps.trisoup_node_size_log2;
  for (int depth = 0; depth < maxDepth; depth++) {
    // setyo at the start of each level
    auto fifoCurrLvlEnd = fifo.end();
    int numNodesNextLvl = 0;
    Vec3<int32_t> occupancyAtlasOrigin = 0xffffffff;

    // derive per-level node size related parameters
    auto parentNodeSizeLog2 = nodeSizeLog2;
    nodeSizeLog2 = lvlNodeSizeLog2[depth];
    auto childSizeLog2 = lvlNodeSizeLog2[depth + 1];
    auto grandchildSizeLog2 = lvlNodeSizeLog2[depth + 2];

    nodeMaxDimLog2 =
      std::max({nodeSizeLog2[0], nodeSizeLog2[1], nodeSizeLog2[2]});

    // if one dimension is not split, atlasShift[k] = 0
    int atlasShift = 7 & ~nonSplitQtBtAxes(parentNodeSizeLog2, nodeSizeLog2);
    int occupancySkipLevel = nonSplitQtBtAxes(nodeSizeLog2, childSizeLog2);
    int childOccupancySkipLevel =
      nonSplitQtBtAxes(childSizeLog2, grandchildSizeLog2);

    auto pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);

    // determing a per node QP at the appropriate level
    if (--numLvlsUntilQuantization == 0) {
      quantNodeSizeLog2 = nodeSizeLog2;
      if (!depth)
        fifo.front().qp = sliceQp;
      else
        calculateNodeQps(gps.geom_base_qp, fifo.begin(), fifoCurrLvlEnd);
    }

    if (gps.geom_planar_mode_enabled_flag) {
      int planarDepth = gbh.geomMaxNodeSizeLog2(gps)
        - std::min({childSizeLog2[0], childSizeLog2[1], childSizeLog2[2]});
      planarInitPlanes(
        kNumPlanarPlanes, planarDepth, planes3x3.data(), planes);
    }

    encoder.beginOctreeLevel();

    // process all nodes within a single level
    for (; fifo.begin() != fifoCurrLvlEnd; fifo.pop_front()) {
      PCCOctree3Node& node0 = fifo.front();

      // encode delta qp for each octree block
      if (numLvlsUntilQuantization == 0)
        encoder.encodeQpOffset(node0.qp - sliceQp);

      int shiftBits = (node0.qp - 4) / 6;
      auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;
      auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;

      // make quantisation work with qtbt and planar.
      int occupancySkip = occupancySkipLevel;
      int childOccupancySkip = childOccupancySkipLevel;
      if (shiftBits != 0) {
        for (int k = 0; k < 3; k++) {
          if (effectiveChildSizeLog2[k] < 0)
            occupancySkip |= (4 >> k);
          if (effectiveChildSizeLog2[k] < 1)
            childOccupancySkip |= (4 >> k);
        }
      }

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
          return !!(int(point[2]) & pointSortMask[2])
            | (!!(int(point[1]) & pointSortMask[1]) << 1)
            | (!!(int(point[0]) & pointSortMask[0]) << 2);
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
      if (nodeMaxDimLog2 < gps.intra_pred_max_node_size_log2) {
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

        // planar mode for current node
        // mask to be used for the occupancy coding
        // (bit =1 => occupancy bit not coded due to not belonging to the plane)
        int planarMask[3] = {0, 0, 0};
        maskPlanar(node0, planarMask, occupancySkip);

        encoder.encodeOccupancy(
          node0.neighPattern, occupancy, occupancyIsPredicted,
          occupancyPrediction, occupancyAdjacencyGt0, occupancyAdjacencyGt1,
          occupancyAdjacencyUnocc, planarMask[0], planarMask[1], planarMask[2],
          node0.planarPossible & 1, node0.planarPossible & 2,
          node0.planarPossible & 4);
      }

      // planar eligibility
      bool planarEligible[3] = {false, false, false};
      if (gps.geom_planar_mode_enabled_flag) {
        // update the plane rate depending on the occupancy and local density
        updateplanarRate(planarRate, occupancy, localDensity, numSiblings);
        eligilityPlanar(
          planarEligible, planarRate, planarRateThreshold, localDensity);
        if (childOccupancySkip & 4)
          planarEligible[0] = false;
        if (childOccupancySkip & 2)
          planarEligible[1] = false;
        if (childOccupancySkip & 1)
          planarEligible[2] = false;

        // avoid mismatch when the next level will apply quantization
        if (numLvlsUntilQuantization == 1) {
          planarEligible[0] = false;
          planarEligible[1] = false;
          planarEligible[2] = false;
        }
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

          childStart = childEnd;

          // if the bitstream is configured to represent unique points,
          // no point count is sent.
          if (gps.geom_unique_points_flag) {
            assert(childCounts[i] == 1);
            continue;
          }

          encoder.encodePositionLeafNumPoints(childCounts[i]);
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
        child.pos[0] = (node0.pos[0] << !(occupancySkipLevel & 4)) + x;
        child.pos[1] = (node0.pos[1] << !(occupancySkipLevel & 2)) + y;
        child.pos[2] = (node0.pos[2] << !(occupancySkipLevel & 1)) + z;

        child.start = childPointsStartIdx;
        childPointsStartIdx += childCounts[i];
        child.end = childPointsStartIdx;
        child.numSiblingsPlus1 = numSiblings;
        child.siblingOccupancy = occupancy;

        // determine planarity if eligible
        int planarProb[3] = {127, 127, 127};
        if (
          gps.geom_planar_mode_enabled_flag
          && (planarEligible[0] || planarEligible[1] || planarEligible[2]))
          encoder.determinePlanarMode(
            pointCloud, planarEligible, childSizeLog2, kNumPlanarPlanes, child,
            planes, node0.neighPattern, x, y, z, planarProb, planarRate);

        // IDCM
        bool idcmEnabled = gps.inferred_direct_coding_mode_enabled_flag
          && planarProb[0] * planarProb[1] * planarProb[2] <= idcmThreshold;
        if (isDirectModeEligible(idcmEnabled, nodeMaxDimLog2, node0, child)) {
          bool directModeUsed = encoder.encodeDirectPosition(
            gps.geom_unique_points_flag, effectiveChildSizeLog2, shiftBits,
            child, pointCloud);

          if (directModeUsed) {
            // inverse quantise any quantised positions
            geometryScale(pointCloud, node0, quantNodeSizeLog2);

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
            gps.neighbour_context_restriction_flag, fifo.end(),
            numNodesNextLvl, child, i, node0.neighPattern, occupancy);
        }
      }
    }
  }

  // return partial coding result
  //  - add missing levels to node positions
  //  - inverse quantise the point cloud
  // todo(df): this does not yet support inverse quantisation of node.pos
  if (nodesRemaining) {
    nodeSizeLog2 = lvlNodeSizeLog2[maxDepth];
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
