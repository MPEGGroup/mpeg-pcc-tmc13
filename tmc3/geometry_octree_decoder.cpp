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
    int mappedPlanarMaskX,
    int mappedFixedMaskX0,
    bool planarPossibleX,
    int mappedPlanarMaskY,
    int mappedFixedMaskY0,
    bool planarPossibleY,
    int mappedPlanarMaskZ,
    int mappedFixedMaskZ0,
    bool planarPossibleZ);

  int decodeOccupancyNeighNZ(
    int neighPattern,
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

  int decodeOccupancyBitwise(
    int neighPattern,
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

  int decodeOccupancyBytewise(int neighPattern);

  int decodePlanarMode(
    PCCOctree3Node& node0,
    int planeZ,
    int posz,
    int dist,
    int neighb,
    int& h,
    int planeId);

  void determinePlanarMode(
    int index,
    const int kNumPlanarPlanes,
    PCCOctree3Node& child,
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
    bool planarEligible[3],
    const int kNumPlanarPlanes,
    PCCOctree3Node& child,
    int* planes[9],
    uint8_t neighPattern,
    int x,
    int y,
    int z,
    int planarProb[3],
    int planarRate[3]);

  uint32_t decodeOccupancy(
    int neighPattern,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int occupancyAdjGt0,
    int occupancyAdjGt1,
    int occupancyAdjUncc,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ);

  Vec3<int32_t> decodePointPosition(
    const Vec3<int>& nodeSizeLog2, uint8_t planarMode, uint8_t planePosBits);

  int decodeQpOffset();

  template<class OutputIt>
  int decodeDirectPosition(
    bool geom_unique_points_flag,
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
  AdaptiveBitModel _ctxSingleIdcmDupPoint;
  AdaptiveBitModel _ctxPointCountPerBlock;
  AdaptiveBitModel _ctxBlockSkipTh;
  AdaptiveBitModel _ctxNumIdcmPointsGt1;
  AdaptiveBitModel _ctxQpOffsetIsZero;
  AdaptiveBitModel _ctxQpOffsetSign;
  AdaptiveBitModel _ctxQpOffsetAbsEgl;
  AdaptiveBitModel _ctxSameZ;

  // for planar mode xyz
  AdaptiveBitModel _ctxPlanarMode[3][2][2];
  AdaptiveBitModel _ctxPlanarPlaneLastIndex[3][2][8][2];
  AdaptiveBitModel _ctxPlanarPlaneLastIndexZ[3];

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

//============================================================================

int
GeometryOctreeDecoder::decodePlanarMode(
  PCCOctree3Node& node0,
  int planeZ,
  int posz,
  int dist,
  int neighb,
  int& h,
  int planeId)
{
  const int mask0 = (1 << planeId);
  const int mask1[3] = {6, 5, 3};

  // decode planar mode
  int discreteDist = dist <= 2 ? 0 : 1;
  bool isPlanar =
    _arithmeticDecoder->decode(_ctxPlanarMode[planeId][neighb][discreteDist]);
  node0.planarMode |= isPlanar ? mask0 : 0;

  if (!isPlanar) {
    node0.planarPossible &= mask1[planeId];
    return -1;
  }

  // decode the plane index // encode the plane index
  int planeBit;
  if (planeZ < 0) {
    planeBit = _arithmeticDecoder->decode(_ctxPlanarPlaneLastIndexZ[planeId]);
    h = approxSymbolProbability(planeBit, _ctxPlanarPlaneLastIndexZ[planeId]);
  } else {
    discreteDist += (dist <= 16 ? 0 : 1);
    int lastIndexPlane2d = planeZ + (discreteDist << 1);
    planeBit = _arithmeticDecoder->decode(
      _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d][posz]);
    h = approxSymbolProbability(
      planeBit,
      _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d][posz]);
  }

  node0.planePosBits |= (planeBit << planeId);
  return planeBit;
}

//============================================================================

void
GeometryOctreeDecoder::determinePlanarMode(
  int planeId,
  const int kNumPlanarPlanes,
  PCCOctree3Node& child,
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
  int planeBit = decodePlanarMode(
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

void
GeometryOctreeDecoder::determinePlanarMode(
  bool planarEligible[3],
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
  int xx = child.pos[0];
  int yy = child.pos[1];
  int zz = child.pos[2];

  // planar x
  if (planarEligible[0]) {
    determinePlanarMode(
      0, kNumPlanarPlanes, child, planes[1], planes[2], planes[0], yy, zz, xx,
      neighPattern, x, planarProb, planarRate);
  }
  // planar y
  if (planarEligible[1]) {
    determinePlanarMode(
      1, kNumPlanarPlanes, child, planes[3], planes[5], planes[4], xx, zz, yy,
      neighPattern, y, planarProb, planarRate);
  }
  // planar z
  if (planarEligible[2]) {
    determinePlanarMode(
      2, kNumPlanarPlanes, child, planes[6], planes[7], planes[8], xx, yy, zz,
      neighPattern, z, planarProb, planarRate);
  }
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
  int occupancy = 0;

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
    int bit = 0;
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

    // NB: There must be at least two occupied child nodes
    //  -- avoid coding the occupancy bit if it is implied.
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
    if (bitIsOneX || bitIsOneY || bitIsOneZ) {
      bit = 1;
    } else {
      auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];
      int ctxIdx = ctxIdxMap[i][numOccupiedAcc];
      bit = _arithmeticDecoder->decode(_ctxOccupancy[ctxIdx]);
      ctxIdxMap.evolve(bit, &ctxIdxMap[i][numOccupiedAcc]);

      if (!bit) {
        coded0X[mask0X]++;
        coded0Y[mask0Y]++;
        coded0Z[mask0Z]++;
      }
    }

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

  int occupancy = 0;
  int partialOccupancy = 0;

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
  // NB: offsets are added since ctxIdxMap is shared between Z and NZ cases.
  for (int i = 0; i < 8; i++) {
    // NB: if firt 7 bits are 0, then the last is implicitly 1.
    int bit = 0;
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
    if (!maskedOccupancy) {  // !planar
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

    if (bitIsOneX || bitIsOneY || bitIsOneZ) {
      bit = 1;
    } else {
      auto& ctxIdxMap = _ctxIdxMaps[ctxIdxMapIdx];
      int ctxIdx = ctxIdxMap[i][idx];
      bit = _arithmeticDecoder->decode(_ctxOccupancy[ctxIdx]);
      ctxIdxMap.evolve(bit, &ctxIdxMap[i][idx]);

      if (!bit) {
        coded0X[mask0X]++;
        coded0Y[mask0Y]++;
        coded0Z[mask0Z]++;
      }
    }

    partialOccupancy |= bit << i;
    occupancy |= bit << kOccBitCodingOrder[i];
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
    return decodeOccupancyNeighZ(
      mappedOccIsPredicted, mappedOccPrediction, mappedOccAdjGt0,
      mappedOccAdjGt1, mappedOccAdjUnocc, mappedPlanarMaskX, mappedFixedMaskX0,
      planarPossibleX, mappedPlanarMaskY, mappedFixedMaskY0, planarPossibleY,
      mappedPlanarMaskZ, mappedFixedMaskZ0, planarPossibleZ);
  }

  return decodeOccupancyNeighNZ(
    neighPattern, mappedOccIsPredicted, mappedOccPrediction, mappedOccAdjGt0,
    mappedOccAdjGt1, mappedOccAdjUnocc, mappedPlanarMaskX, mappedFixedMaskX0,
    planarPossibleX, mappedPlanarMaskY, mappedFixedMaskY0, planarPossibleY,
    mappedPlanarMaskZ, mappedFixedMaskZ0, planarPossibleZ);
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
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  bool planarPossibleX,
  bool planarPossibleY,
  bool planarPossibleZ)
{
  // decode occupancy pattern
  uint32_t occupancy;

  // single child and we know its position
  if (planarMaskX && planarMaskY && planarMaskZ) {
    uint32_t cnt = (planarMaskZ & 1);
    cnt |= (planarMaskY & 1) << 1;
    cnt |= (planarMaskX & 1) << 2;
    occupancy = 1 << cnt;
    return occupancy;
  }

  // neighbour empty and only one point => decode index, not pattern
  if (neighPattern == 0) {
    bool singleChild = false;
    if (planarPossibleX && planarPossibleY && planarPossibleZ) {
      singleChild = _arithmeticDecoder->decode(_ctxSingleChild) == 1;
    }

    if (singleChild) {
      uint32_t cnt;
      if (!planarMaskZ)
        cnt = _arithmeticDecoder->decode(_ctxEquiProb);
      else
        cnt = (planarMaskZ & 1);

      if (!planarMaskY)
        cnt |= _arithmeticDecoder->decode(_ctxEquiProb) << 1;
      else
        cnt |= (planarMaskY & 1) << 1;

      if (!planarMaskX)
        cnt |= _arithmeticDecoder->decode(_ctxEquiProb) << 2;
      else
        cnt |= (planarMaskX & 1) << 2;

      occupancy = 1 << cnt;
      return occupancy;
    }
  }

  // at least two child nodes occupied and two planars => we know the occupancy
  if (neighPattern == 0) {
    if (planarMaskX && planarMaskY) {
      uint32_t cnt = ((planarMaskX & 1) << 2) | ((planarMaskY & 1) << 1);
      occupancy = (1 << cnt) | (1 << (cnt + 1));
      return occupancy;
    }

    if (planarMaskY && planarMaskZ) {
      uint32_t cnt = ((planarMaskY & 1) << 1) | (planarMaskZ & 1);
      occupancy = (1 << cnt) | (1 << (cnt + 4));
      return occupancy;
    }

    if (planarMaskX && planarMaskZ) {
      uint32_t cnt = ((planarMaskX & 1) << 2) | (planarMaskZ & 1);
      occupancy = (1 << cnt) | (1 << (cnt + 2));
      return occupancy;
    }
  }

  uint32_t mapOccIsP =
    mapGeometryOccupancy(occupancyIsPred, neighPattern, planarMaskZ);
  uint32_t mapOccP =
    mapGeometryOccupancy(occupancyPred, neighPattern, planarMaskZ);
  uint32_t mapAdjGt0 =
    mapGeometryOccupancy(occupancyAdjGt0, neighPattern, planarMaskZ);
  uint32_t mapAdjGt1 =
    mapGeometryOccupancy(occupancyAdjGt1, neighPattern, planarMaskZ);
  uint32_t mapAdjUnocc =
    mapGeometryOccupancy(occupancyAdjUnocc, neighPattern, planarMaskZ);

  uint32_t mapPlanarMaskX =
    mapGeometryOccupancy(planarMaskX, neighPattern, planarMaskZ);
  uint32_t mapPlanarMaskY =
    mapGeometryOccupancy(planarMaskY, neighPattern, planarMaskZ);
  uint32_t mapPlanarMaskZ =
    mapGeometryOccupancy(planarMaskZ, neighPattern, planarMaskZ);

  uint32_t mapFixedMaskX0 =
    mapGeometryOccupancy(0xf0, neighPattern, planarMaskZ);
  uint32_t mapFixedMaskY0 =
    mapGeometryOccupancy(0xcc, neighPattern, planarMaskZ);
  uint32_t mapFixedMaskZ0 =
    mapGeometryOccupancy(0xaa, neighPattern, planarMaskZ);

  uint32_t mappedOccupancy;

  if (_useBitwiseOccupancyCoder)
    mappedOccupancy = decodeOccupancyBitwise(
      neighPattern, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1, mapAdjUnocc,
      mapPlanarMaskX, mapFixedMaskX0, planarPossibleX, mapPlanarMaskY,
      mapFixedMaskY0, planarPossibleY, mapPlanarMaskZ, mapFixedMaskZ0,
      planarPossibleZ);
  else
    mappedOccupancy = decodeOccupancyBytewise(neighPattern);

  return mapGeometryOccupancyInv(mappedOccupancy, neighPattern, planarMaskZ);
}

//-------------------------------------------------------------------------
// Decode a position of a point in a given volume.
Vec3<int32_t>
GeometryOctreeDecoder::decodePointPosition(
  const Vec3<int>& nodeSizeLog2, uint8_t planarMode, uint8_t planePosBits)
{
  Vec3<int32_t> delta{};
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2[k] <= 0)
      continue;

    int i = nodeSizeLog2[k];
    delta[k] <<= 1;
    if (!(planarMode & (1 << k)))
      delta[k] |= _arithmeticDecoder->decode(_ctxEquiProb);
    else
      delta[k] |= (planePosBits & (1 << k) ? 1 : 0);
    i--;

    for (; i > 0; i--) {
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
  bool geom_unique_points_flag,
  const Vec3<int>& nodeSizeLog2,
  const PCCOctree3Node& node,
  OutputIt outputPoints)
{
  bool isDirectMode = _arithmeticDecoder->decode(_ctxBlockSkipTh);
  if (!isDirectMode) {
    return 0;
  }

  int numPoints = 1;
  bool numPointsGt1 = _arithmeticDecoder->decode(_ctxNumIdcmPointsGt1);
  numPoints += numPointsGt1;

  int numDuplicatePoints = 0;
  if (!geom_unique_points_flag && !numPointsGt1) {
    numDuplicatePoints = !_arithmeticDecoder->decode(_ctxSinglePointPerBlock);
    if (numDuplicatePoints) {
      bool singleDup = _arithmeticDecoder->decode(_ctxSingleIdcmDupPoint);
      if (!singleDup)
        numDuplicatePoints += 1
          + _arithmeticDecoder->decodeExpGolomb(
              0, _ctxEquiProb, _ctxPointCountPerBlock);
    }
  }

  Vec3<int32_t> pos;
  for (int i = 0; i < numPoints; i++)
    *(outputPoints++) = pos =
      decodePointPosition(nodeSizeLog2, node.planarMode, node.planePosBits);

  for (int i = 0; i < numDuplicatePoints; i++)
    *(outputPoints++) = pos;

  return numPoints + numDuplicatePoints;
}

//-------------------------------------------------------------------------
// Helper to inverse quantise positions

Vec3<int32_t>
invQuantPosition(int qp, Vec3<uint32_t> quantMasks, const Vec3<int32_t>& pos)
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
  Vec3<int32_t> recon;
  for (int k = 0; k < 3; k++) {
    int posQuant = pos[k] & (quantMasks[k] >> shiftBits);
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

  // implicit qtbt for grand-child node
  Vec3<int> grandchildSizeLog2 = implicitQtBtDecision(
    childSizeLog2, maxNumImplicitQtbtBeforeOt, minSizeImplicitQtbt);
  int childOccupancySkip = nonSplitQtBtAxes(childSizeLog2, grandchildSizeLog2);

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(0);
  node00.pos = int32_t(0);
  node00.neighPattern = 0;
  node00.numSiblingsPlus1 = 8;
  node00.siblingOccupancy = 0;
  node00.qp = 4;
  node00.planarMode = 0;

  size_t processedPointCount = 0;
  std::vector<uint32_t> values;

  auto fifoCurrLvlEnd = fifo.end();

  // planar mode initialazation
  const int th_idcm = gps.geom_planar_idcm_threshold * 127 * 127;

  const int kNumPlanarPlanes = 4;
  std::vector<int> planes3x3;
  int* planes[9];
  int depth = 1;

  if (gps.geom_planar_mode_enabled_flag) {
    const int max_plane_size = kNumPlanarPlanes << nodeMaxDimLog2;
    planes3x3.resize(max_plane_size * 9);
    depth = nodeMaxDimLog2;
    planarInitPlanes(kNumPlanarPlanes, depth, planes3x3.data(), planes);
  }
  int planarRate[3] = {128 * 8, 128 * 8, 128 * 8};
  int localDensity = 1024 * 4;
  const int planarRateThreshold[3] = {gps.geom_planar_threshold0 << 4,
                                      gps.geom_planar_threshold1 << 4,
                                      gps.geom_planar_threshold2 << 4};

  // this counter represents fifo.end() - fifoCurrLvlEnd().
  // ie, the number of nodes added to the next level of the tree
  int numNodesNextLvl = 0;

  Vec3<int32_t> occupancyAtlasOrigin{-1};
  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2) {
    occupancyAtlas.resize(gps.neighbour_avail_boundary_log2);
    occupancyAtlas.clear();
  }

  int sliceQp = gps.geom_base_qp + gbh.geom_slice_qp_offset;
  int numLvlsUntilQpOffset = -1;
  Vec3<uint32_t> posQuantBitMasks = 0xffffffff;

  if (gbh.geom_octree_qp_offset_enabled_flag)
    numLvlsUntilQpOffset = gbh.geom_octree_qp_offset_depth;
  else if (gps.geom_scaling_enabled_flag)
    node00.qp = sliceQp;

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

      // if all dimensions have same size, then use octree for remaining nodes
      if (
        nodeSizeLog2[0] == nodeSizeLog2[1]
        && nodeSizeLog2[1] == nodeSizeLog2[2])
        minSizeImplicitQtbt = 0;

      // implicit qtbt for child nodes
      childSizeLog2 = implicitQtBtDecision(
        nodeSizeLog2, maxNumImplicitQtbtBeforeOt, minSizeImplicitQtbt);

      // implicit qtbt for grand-child nodes
      int minSizeImplicitQtbt2 = minSizeImplicitQtbt;
      if (
        childSizeLog2[0] == childSizeLog2[1]
        && childSizeLog2[1] == childSizeLog2[2])
        minSizeImplicitQtbt2 = 0;

      grandchildSizeLog2 = implicitQtBtDecision(
        childSizeLog2,
        maxNumImplicitQtbtBeforeOt ? maxNumImplicitQtbtBeforeOt - 1 : 0,
        minSizeImplicitQtbt2);

      childOccupancySkip = nonSplitQtBtAxes(childSizeLog2, grandchildSizeLog2);

      nodeMaxDimLog2--;
      numNodesNextLvl = 0;
      occupancyAtlasOrigin = 0xffffffff;

      decoder.beginOctreeLevel();
      if (gps.geom_planar_mode_enabled_flag) {
        depth = gbh.geomMaxNodeSizeLog2(gps)
          - std::min({childSizeLog2[0], childSizeLog2[1], childSizeLog2[2]});
        planarInitPlanes(kNumPlanarPlanes, depth, planes3x3.data(), planes);
      }

      // allow partial tree encoding using trisoup
      if (nodeMaxDimLog2 == gps.trisoup_node_size_log2)
        break;

      // allow partial tree decoding
      if (nodeMaxDimLog2 == minNodeSizeLog2)
        break;

      // record the node size when quantisation is signalled -- all subsequnt
      // coded occupancy bits are quantised
      numLvlsUntilQpOffset--;
      if (!numLvlsUntilQpOffset)
        for (int k = 0; k < 3; k++)
          posQuantBitMasks[k] = (1 << nodeSizeLog2[k]) - 1;
    }

    PCCOctree3Node& node0 = fifo.front();

    if (numLvlsUntilQpOffset == 0)
      node0.qp = decoder.decodeQpOffset() + sliceQp;

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

      // planar mode for current node
      // mask to be used for the occupancy coding
      // (bit =1 => occupancy bit not coded due to not belonging to the plane)
      int mask_planar[3] = {0, 0, 0};
      maskPlanar(node0, mask_planar, occupancySkip);

      occupancy = decoder.decodeOccupancy(
        node0.neighPattern, occupancyIsPredicted, occupancyPrediction,
        occupancyAdjacencyGt0, occupancyAdjacencyGt1, occupancyAdjacencyUnocc,
        mask_planar[0], mask_planar[1], mask_planar[2],
        node0.planarPossible & 1, node0.planarPossible & 2,
        node0.planarPossible & 4);
    }

    assert(occupancy > 0);

    // update atlas for advanced neighbours
    if (gps.neighbour_avail_boundary_log2) {
      updateGeometryOccupancyAtlasOccChild(
        node0.pos, occupancy, &occupancyAtlas);
    }

    // population count of occupancy for IDCM
    int numOccupied = popcnt(occupancy);

    // planar eligibility
    bool planarEligible[3] = {false, false, false};
    if (gps.geom_planar_mode_enabled_flag) {
      // update the plane rate depending on the occupancy and local density
      updateplanarRate(planarRate, occupancy, localDensity, numOccupied);
      eligilityPlanar(
        planarEligible, planarRate, planarRateThreshold, localDensity);
      if (childOccupancySkip & 4)
        planarEligible[0] = false;
      if (childOccupancySkip & 2)
        planarEligible[1] = false;
      if (childOccupancySkip & 1)
        planarEligible[2] = false;
    }

    // nodeSizeLog2 > 1: for each child:
    //  - determine elegibility for IDCM
    //  - directly decode point positions if IDCM allowed and selected
    //  - otherwise, insert split children into fifo while updating neighbour state
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
        Vec3<int32_t> point{(node0.pos[0] << !(occupancySkip & 4)) + x,
                            (node0.pos[1] << !(occupancySkip & 2)) + y,
                            (node0.pos[2] << !(occupancySkip & 1)) + z};

        point = invQuantPosition(node0.qp, posQuantBitMasks, point);

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

      // decode planarity if eligible
      int planarProb[3] = {127, 127, 127};
      if (
        gps.geom_planar_mode_enabled_flag
        && (planarEligible[0] || planarEligible[1] || planarEligible[2]))
        decoder.determinePlanarMode(
          planarEligible, kNumPlanarPlanes, child, planes, node0.neighPattern,
          x, y, z, planarProb, planarRate);

      bool idcmEnabled = gps.inferred_direct_coding_mode_enabled_flag
        && planarProb[0] * planarProb[1] * planarProb[2] <= th_idcm;
      if (isDirectModeEligible(
            idcmEnabled, effectiveNodeMaxDimLog2, node0, child)) {
        int numPoints = decoder.decodeDirectPosition(
          gps.geom_unique_points_flag, effectiveChildSizeLog2, child,
          &pointCloud[processedPointCount]);

        for (int j = 0; j < numPoints; j++) {
          auto& point = pointCloud[processedPointCount++];
          for (int k = 0; k < 3; k++) {
            int shift = std::max(0, effectiveChildSizeLog2[k]);
            point[k] += child.pos[k] << shift;
          }

          point = invQuantPosition(node0.qp, posQuantBitMasks, point);
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
      node.pos = invQuantPosition(node.qp, posQuantBitMasks, node.pos);
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

    for (auto node0 : nodes)
      pointCloud[processedPointCount++] = node0.pos;
  }
}

//============================================================================

}  // namespace pcc
