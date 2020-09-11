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

class GeometryOctreeDecoder : protected GeometryOctreeContexts {
public:
  GeometryOctreeDecoder(
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    const GeometryOctreeContexts& ctxMem,
    EntropyDecoder* arithmeticDecoder);

  GeometryOctreeDecoder(const GeometryOctreeDecoder&) = default;
  GeometryOctreeDecoder(GeometryOctreeDecoder&&) = default;
  GeometryOctreeDecoder& operator=(const GeometryOctreeDecoder&) = default;
  GeometryOctreeDecoder& operator=(GeometryOctreeDecoder&&) = default;

  void beginOctreeLevel(const Vec3<int>& planarDepth);

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
    OctreeNodePlanar& planar,
    int planeZ,
    int dist,
    int neighb,
    int& h,
    int planeId,
    int contextAngle);

  void determinePlanarMode(
    int planeId,
    OctreeNodePlanar& child,
    OctreePlanarBuffer::Row* planeBuffer,
    int coord1,
    int coord2,
    int coord3,
    uint8_t neighPattern,
    int planarProb[3],
    int planarRate[3],
    int contextAngle);

  void determinePlanarMode(
    const bool planarEligible[3],
    PCCOctree3Node& child,
    OctreeNodePlanar& planar,
    uint8_t neighPattern,
    int planarProb[3],
    int contextAngle,
    int contextAnglePhiX,
    int contextAnglePhiY);

  uint32_t decodeOccupancy(
    const GeometryNeighPattern& gnp,
    int occupancyIsPredicted,
    int occupancyPrediction,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ);

  Vec3<int32_t> decodePointPosition(
    const Vec3<int>& nodeSizeLog2, Vec3<int32_t>& deltaPlanar);

  void decodeOrdered2ptPrefix(
    Vec3<bool> directIdcm,
    Vec3<int>& nodeSizeLog2AfterUnordered,
    Vec3<int32_t> deltaUnordered[2]);

  Vec3<int32_t> decodePointPositionAngular(
    const Vec3<int>& nodeSizeLog2,
    const Vec3<int>& nodeSizeLog2AfterPlanar,
    const PCCOctree3Node& node,
    const OctreeNodePlanar& planar,
    const Vec3<int>& headPos,
    const int* zLaser,
    const int* thetaLaser,
    Vec3<int32_t>& deltaPlanar);

  int decodeQpOffset();

  bool decodeIsIdcm();

  template<class OutputIt>
  int decodeDirectPosition(
    bool geom_unique_points_flag,
    bool joint_2pt_idcm_enabled_flag,
    const Vec3<int>& nodeSizeLog2,
    PCCOctree3Node& node,
    OctreeNodePlanar& planar,
    OutputIt outputPoints,
    bool angularIdcm,
    const Vec3<int>& headPos,
    const int* zLaser,
    const int* thetaLaser,
    int numLasers);

  int decodeThetaRes();

  const GeometryOctreeContexts& getCtx() const { return *this; }

public:
  // selects between the bitwise and bytewise occupancy coders
  bool _useBitwiseOccupancyCoder;

  const uint8_t* _neighPattern64toR1;

  EntropyDecoder* _arithmeticDecoder;

  // Planar state
  OctreePlanarState _planar;

  // Azimuthal buffer
  std::vector<int> _phiBuffer;

  // azimuthal elementary shifts
  AzimuthalPhiZi _phiZi;
};

//============================================================================

GeometryOctreeDecoder::GeometryOctreeDecoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const GeometryOctreeContexts& ctxtMem,
  EntropyDecoder* arithmeticDecoder)
  : GeometryOctreeContexts(ctxtMem)
  , _useBitwiseOccupancyCoder(gps.bitwise_occupancy_coding_flag)
  , _neighPattern64toR1(neighPattern64toR1(gps))
  , _arithmeticDecoder(arithmeticDecoder)
  , _planar(gps)
  , _phiBuffer(gps.geom_angular_num_lidar_lasers(), 0x80000000)
  , _phiZi(
      gps.geom_angular_num_lidar_lasers(), gps.geom_angular_num_phi_per_turn)
{
  if (!_useBitwiseOccupancyCoder && !gbh.entropy_continuation_flag) {
    for (int i = 0; i < 10; i++)
      _bytewiseOccupancyCoder[i].init(kDualLutOccupancyCoderInit[i]);
  }
}
//============================================================================

void
GeometryOctreeDecoder::beginOctreeLevel(const Vec3<int>& planarDepth)
{
  for (int i = 0; i < 10; i++) {
    _bytewiseOccupancyCoder[i].resetLut();
  }

  _planar.initPlanes(planarDepth);
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
    count++;
    count += _arithmeticDecoder->decodeExpGolomb(0, _ctxPointCountPerBlock);
  }

  return count;
}

//============================================================================

int
GeometryOctreeDecoder::decodePlanarMode(
  OctreeNodePlanar& planar,
  int planeZ,
  int dist,
  int neighb,
  int& h,
  int planeId,
  int contextAngle)
{
  const int mask0 = (1 << planeId);
  const int mask1[3] = {6, 5, 3};

  // decode planar mode
  bool isPlanar = _arithmeticDecoder->decode(_ctxPlanarMode[planeId]);
  planar.planarMode |= isPlanar ? mask0 : 0;

  if (!isPlanar) {
    planar.planarPossible &= mask1[planeId];
    return -1;
  }

  // decode the plane index // encode the plane index
  int planeBit;
  if (contextAngle == -1) {  // angular mode off
    if (planeZ < 0) {
      planeBit =
        _arithmeticDecoder->decode(_ctxPlanarPlaneLastIndexZ[planeId]);
      h =
        approxSymbolProbability(planeBit, _ctxPlanarPlaneLastIndexZ[planeId]);
    } else {
      int discreteDist = (dist <= (2 >> OctreePlanarBuffer::shiftAb) ? 0 : 1);
      discreteDist += (dist <= (16 >> OctreePlanarBuffer::shiftAb) ? 0 : 1);
      int lastIndexPlane2d = planeZ + (discreteDist << 1);

      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d]);
      h = approxSymbolProbability(
        planeBit, _ctxPlanarPlaneLastIndex[planeId][neighb][lastIndexPlane2d]);
    }
  } else {               // angular mode on
    if (planeId == 2) {  // angular
      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndexAngular[contextAngle]);
      h = approxSymbolProbability(
        planeBit, _ctxPlanarPlaneLastIndexAngular[contextAngle]);
    } else {  // azimuthal
      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndexAngularPhi[contextAngle]);
      h = approxSymbolProbability(
        planeBit, _ctxPlanarPlaneLastIndexAngularPhi[contextAngle]);
    }
  }

  planar.planePosBits |= (planeBit << planeId);
  return planeBit;
}

//============================================================================

void
GeometryOctreeDecoder::determinePlanarMode(
  int planeId,
  OctreeNodePlanar& planar,
  OctreePlanarBuffer::Row* planeBuffer,
  int coord1,
  int coord2,
  int coord3,
  uint8_t neighPattern,
  int planarProb[3],
  int planarRate[3],
  int contextAngle)
{
  const int kPlanarChildThreshold = 63;
  const int kAdjNeighIdxFromPlanePos[3][2] = {1, 0, 2, 3, 4, 5};
  const int planeSelector = 1 << planeId;

  OctreePlanarBuffer::Elmt* row;
  int rowLen = OctreePlanarBuffer::rowSize;
  int closestPlanarFlag;
  int closestDist;

  if (!planeBuffer) {
    // angular: buffer disabled
    closestPlanarFlag = 0;
    closestDist = 0;
  } else {
    coord1 =
      (coord1 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
    coord2 =
      (coord2 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
    coord3 = coord3 & OctreePlanarBuffer::maskC;

    row = planeBuffer[coord3];

    int minDist = std::abs(coord1 - int(row[rowLen - 1].a))
      + std::abs(coord2 - int(row[rowLen - 1].b));
    int idxMinDist = rowLen - 1;

    for (int idxP = 0; idxP < rowLen - 1; idxP++) {
      int dist0 = std::abs(coord1 - int(row[idxP].a))
        + std::abs(coord2 - int(row[idxP].b));
      if (dist0 < minDist) {
        idxMinDist = idxP;
        minDist = dist0;
      }
    }

    // push closest point front
    row[rowLen - 1] = row[idxMinDist];

    closestPlanarFlag = row[idxMinDist].planeIdx;
    closestDist = minDist;

    for (int idxP = 0; idxP < rowLen - 1; idxP++) {
      row[idxP] = row[idxP + 1];
    }
  }
  const int kAdjNeighIdxFromPlaneMask[3] = {0, 2, 4};
  int adjNeigh = (neighPattern >> kAdjNeighIdxFromPlaneMask[planeId]) & 3;
  int planeBit = decodePlanarMode(
    planar, closestPlanarFlag, closestDist, adjNeigh, planarProb[planeId],
    planeId, contextAngle);

  bool isPlanar = (planar.planarMode & planeSelector)
    && planarProb[planeId] > kPlanarChildThreshold;

  planarRate[planeId] =
    (255 * planarRate[planeId] + (isPlanar ? 256 * 8 : 0) + 128) >> 8;

  if (planeBuffer) {
    row[rowLen - 1] = {unsigned(coord1), planeBit, unsigned(coord2)};
  }
}

//============================================================================

void
GeometryOctreeDecoder::determinePlanarMode(
  const bool planarEligible[3],
  PCCOctree3Node& child,
  OctreeNodePlanar& planar,
  uint8_t neighPattern,
  int planarProb[3],
  int contextAngle,
  int contextAnglePhiX,
  int contextAnglePhiY)
{
  int xx = child.pos[0];
  int yy = child.pos[1];
  int zz = child.pos[2];

  auto& planeBuffer = _planar._planarBuffer;

  // planar x
  if (planarEligible[0]) {
    determinePlanarMode(
      0, planar, planeBuffer.getBuffer(0), yy, zz, xx, neighPattern,
      planarProb, _planar._rate.data(), contextAnglePhiX);
  }
  // planar y
  if (planarEligible[1]) {
    determinePlanarMode(
      1, planar, planeBuffer.getBuffer(1), xx, zz, yy, neighPattern,
      planarProb, _planar._rate.data(), contextAnglePhiY);
  }
  // planar z
  if (planarEligible[2]) {
    determinePlanarMode(
      2, planar, planeBuffer.getBuffer(2), xx, yy, zz, neighPattern,
      planarProb, _planar._rate.data(), contextAngle);
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

    int ctxIdxMapIdx = 3 * idxAdj;
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

    int ctxIdxMapIdx = 3 * idxAdj;
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
  const GeometryNeighPattern& gnp,
  int occupancyIsPred,
  int occupancyPred,
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
  if (gnp.neighPattern == 0) {
    bool singleChild = false;
    if (planarPossibleX && planarPossibleY && planarPossibleZ) {
      singleChild = _arithmeticDecoder->decode(_ctxSingleChild) == 1;
    }

    if (singleChild) {
      uint32_t cnt;
      if (!planarMaskZ)
        cnt = _arithmeticDecoder->decode();
      else
        cnt = (planarMaskZ & 1);

      if (!planarMaskY)
        cnt |= _arithmeticDecoder->decode() << 1;
      else
        cnt |= (planarMaskY & 1) << 1;

      if (!planarMaskX)
        cnt |= _arithmeticDecoder->decode() << 2;
      else
        cnt |= (planarMaskX & 1) << 2;

      occupancy = 1 << cnt;
      return occupancy;
    }
  }

  // at least two child nodes occupied and two planars => we know the occupancy
  if (gnp.neighPattern == 0) {
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

  auto neighPattern = gnp.neighPattern;
  auto mapOccIsP = mapGeometryOccupancy(occupancyIsPred, neighPattern);
  auto mapOccP = mapGeometryOccupancy(occupancyPred, neighPattern);
  auto mapAdjGt0 = mapGeometryOccupancy(gnp.adjacencyGt0, neighPattern);
  auto mapAdjGt1 = mapGeometryOccupancy(gnp.adjacencyGt1, neighPattern);
  auto mapAdjUnocc = mapGeometryOccupancy(gnp.adjacencyUnocc, neighPattern);

  auto mapPlanarMaskX = mapGeometryOccupancy(planarMaskX, neighPattern);
  auto mapPlanarMaskY = mapGeometryOccupancy(planarMaskY, neighPattern);
  auto mapPlanarMaskZ = mapGeometryOccupancy(planarMaskZ, neighPattern);

  auto mapFixedMaskX0 = mapGeometryOccupancy(0xf0, neighPattern);
  auto mapFixedMaskY0 = mapGeometryOccupancy(0xcc, neighPattern);
  auto mapFixedMaskZ0 = mapGeometryOccupancy(0xaa, neighPattern);

  uint32_t mappedOccupancy;

  if (_useBitwiseOccupancyCoder)
    mappedOccupancy = decodeOccupancyBitwise(
      neighPattern, mapOccIsP, mapOccP, mapAdjGt0, mapAdjGt1, mapAdjUnocc,
      mapPlanarMaskX, mapFixedMaskX0, planarPossibleX, mapPlanarMaskY,
      mapFixedMaskY0, planarPossibleY, mapPlanarMaskZ, mapFixedMaskZ0,
      planarPossibleZ);
  else
    mappedOccupancy = decodeOccupancyBytewise(neighPattern);

  return mapGeometryOccupancyInv(mappedOccupancy, neighPattern);
}

//-------------------------------------------------------------------------
// Decode a position of a point in a given volume.
Vec3<int32_t>
GeometryOctreeDecoder::decodePointPosition(
  const Vec3<int>& nodeSizeLog2, Vec3<int32_t>& deltaPlanar)
{
  Vec3<int32_t> delta = deltaPlanar;
  for (int k = 0; k < 3; k++) {
    if (nodeSizeLog2[k] <= 0)
      continue;

    for (int i = nodeSizeLog2[k]; i > 0; i--) {
      delta[k] <<= 1;
      delta[k] |= _arithmeticDecoder->decode();
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
    dqp = _arithmeticDecoder->decodeExpGolomb(0, _ctxQpOffsetAbsEgl) + 1;
    dqp = dqp_sign ? dqp : -dqp;
  }
  return dqp;
}

//-------------------------------------------------------------------------
// Decode part of the position of two unordred points  point in a given volume.
void
GeometryOctreeDecoder::decodeOrdered2ptPrefix(
  Vec3<bool> directIdcm, Vec3<int>& nodeSizeLog2, Vec3<int32_t> pointPrefix[2])
{
  if (nodeSizeLog2[0] >= 1 && directIdcm[0]) {
    int ctxIdx = 0;
    bool sameBit = true;
    while (nodeSizeLog2[0] && sameBit) {
      pointPrefix[0][0] <<= 1;
      pointPrefix[1][0] <<= 1;
      nodeSizeLog2[0]--;

      sameBit = _arithmeticDecoder->decode(_ctxSameBitHighx[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      if (sameBit) {
        int bit = _arithmeticDecoder->decode();
        pointPrefix[0][0] |= bit;
        pointPrefix[1][0] |= bit;
      } else {
        pointPrefix[1][0] |= 1;
      }
    }
  }

  if (nodeSizeLog2[1] >= 1 && directIdcm[1]) {
    int ctxIdx = 0;
    bool sameBit = true;
    bool sameX = !directIdcm[0] || pointPrefix[0][0] == pointPrefix[1][0];

    while (nodeSizeLog2[1] && sameBit) {
      pointPrefix[0][1] <<= 1;
      pointPrefix[1][1] <<= 1;
      nodeSizeLog2[1]--;

      sameBit = _arithmeticDecoder->decode(_ctxSameBitHighy[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      int bit = 0;
      if (!(sameX && !sameBit))
        bit = _arithmeticDecoder->decode();
      pointPrefix[0][1] |= bit;
      pointPrefix[1][1] |= sameBit ? bit : !bit;
    }
  }

  if (nodeSizeLog2[2] >= 1 && directIdcm[2]) {
    int ctxIdx = 0;
    bool sameBit = true;
    bool sameXy = (!directIdcm[0] || pointPrefix[0][0] == pointPrefix[1][0])
      && (!directIdcm[1] || pointPrefix[0][1] == pointPrefix[1][1]);

    while (nodeSizeLog2[2] && sameBit) {
      pointPrefix[0][2] <<= 1;
      pointPrefix[1][2] <<= 1;
      nodeSizeLog2[2]--;

      sameBit = _arithmeticDecoder->decode(_ctxSameBitHighz[ctxIdx]);
      ctxIdx = std::min(4, ctxIdx + 1);
      int bit = 0;
      if (!(sameXy && !sameBit))
        bit = _arithmeticDecoder->decode();
      pointPrefix[0][2] |= bit;
      pointPrefix[1][2] |= sameBit ? bit : !bit;
    }
  }
}

//-------------------------------------------------------------------------
// Decode a position of a point in a given volume, using elevation angle prior

Vec3<int32_t>
GeometryOctreeDecoder::decodePointPositionAngular(
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& nodeSizeLog2AfterPlanar,
  const PCCOctree3Node& child,
  const OctreeNodePlanar& planar,
  const Vec3<int>& headPos,
  const int* zLaser,
  const int* thetaLaser,
  Vec3<int32_t>& deltaPlanar)
{
  Vec3<int32_t> delta = deltaPlanar;
  Vec3<int> posXyz = {(child.pos[0] << nodeSizeLog2[0]) - headPos[0],
                      (child.pos[1] << nodeSizeLog2[1]) - headPos[1],
                      (child.pos[2] << nodeSizeLog2[2]) - headPos[2]};

  // -- PHI --
  // code x or y directly and compute phi of node
  bool codeXorY = std::abs(posXyz[0]) <= std::abs(posXyz[1]);
  if (codeXorY) {  // direct code y
    if (nodeSizeLog2AfterPlanar[1])
      for (int i = nodeSizeLog2AfterPlanar[1]; i > 0; i--) {
        delta[1] <<= 1;
        delta[1] |= _arithmeticDecoder->decode();
      }
    posXyz[1] += delta[1];
    posXyz[0] += delta[0] << nodeSizeLog2AfterPlanar[0];
  } else {  //direct code x
    if (nodeSizeLog2AfterPlanar[0])
      for (int i = nodeSizeLog2AfterPlanar[0]; i > 0; i--) {
        delta[0] <<= 1;
        delta[0] |= _arithmeticDecoder->decode();
      }
    posXyz[0] += delta[0];
    posXyz[1] += delta[1] << nodeSizeLog2AfterPlanar[1];
  }

  // find predictor
  int phiNode = iatan2(posXyz[1], posXyz[0]);
  int laserNode = int(child.laserIndex);

  // laser residual
  int laserIndex = laserNode + decodeThetaRes();

  int predPhi = _phiBuffer[laserIndex];
  if (predPhi == 0x80000000)
    predPhi = phiNode;

  // elementary shift predictor
  int nShift =
    ((predPhi - phiNode) * _phiZi.invDelta(laserIndex) + 536870912) >> 30;
  predPhi -= _phiZi.delta(laserIndex) * nShift;

  // choose x or y
  int* posXY = codeXorY ? &posXyz[0] : &posXyz[1];
  int idx = codeXorY ? 0 : 1;

  // azimuthal code x or y
  int mask2 = codeXorY
    ? (nodeSizeLog2AfterPlanar[0] > 0 ? 1 << (nodeSizeLog2AfterPlanar[0] - 1)
                                      : 0)
    : (nodeSizeLog2AfterPlanar[1] > 0 ? 1 << (nodeSizeLog2AfterPlanar[1] - 1)
                                      : 0);
  for (; mask2; mask2 >>= 1) {
    // angles left and right
    int phiR = codeXorY ? iatan2(posXyz[1], posXyz[0] + mask2)
                        : iatan2(posXyz[1] + mask2, posXyz[0]);
    int phiL = phiNode;

    // ctx azimutal
    int angleL = phiL - predPhi;
    int angleR = phiR - predPhi;
    int contextAnglePhi =
      (angleL >= 0 && angleR >= 0) || (angleL < 0 && angleR < 0) ? 2 : 0;
    angleL = std::abs(angleL);
    angleR = std::abs(angleR);
    if (angleL > angleR) {
      contextAnglePhi++;
      int temp = angleL;
      angleL = angleR;
      angleR = temp;
    }
    if (angleR > (angleL << 1))
      contextAnglePhi += 4;

    // entropy coding
    bool bit = _arithmeticDecoder->decode(
      _ctxPlanarPlaneLastIndexAngularPhiIDCM[contextAnglePhi]);
    delta[idx] <<= 1;
    if (bit) {
      delta[idx] |= 1;
      *posXY += mask2;
      phiNode = phiR;
      predPhi = _phiBuffer[laserIndex];
      if (predPhi == 0x80000000)
        predPhi = phiNode;

      // elementary shift predictor
      int nShift =
        ((predPhi - phiNode) * _phiZi.invDelta(laserIndex) + 536870912) >> 30;
      predPhi -= _phiZi.delta(laserIndex) * nShift;
    }
  }

  // update buffer phi
  _phiBuffer[laserIndex] = phiNode;

  // -- THETA --
  int maskz =
    nodeSizeLog2AfterPlanar[2] > 0 ? 1 << (nodeSizeLog2AfterPlanar[2] - 1) : 0;
  if (!maskz)
    return delta;

  int posz0 = posXyz[2];
  posXyz[2] += delta[2] << nodeSizeLog2AfterPlanar[2];

  // Since x and y are known,
  // r is known too and does not depend on the bit for z
  uint64_t xLidar = (int64_t(posXyz[0]) << 8) - 128;
  uint64_t yLidar = (int64_t(posXyz[1]) << 8) - 128;
  uint64_t r2 = xLidar * xLidar + yLidar * yLidar;
  int64_t rInv = irsqrt(r2);

  // code bits for z using angular. Eligility is implicit. Laser is known.
  int64_t hr = zLaser[laserIndex] * rInv;
  int fixedThetaLaser =
    thetaLaser[laserIndex] + int(hr >= 0 ? -(hr >> 17) : ((-hr) >> 17));

  int zShift = (rInv << nodeSizeLog2AfterPlanar[2]) >> 18;
  for (int bitIdxZ = nodeSizeLog2AfterPlanar[2]; bitIdxZ > 0;
       bitIdxZ--, maskz >>= 1, zShift >>= 1) {
    // determine non-corrected theta
    int64_t zLidar = ((posXyz[2] + maskz) << 1) - 1;
    int64_t theta = zLidar * rInv;
    int theta32 = theta >= 0 ? theta >> 15 : -((-theta) >> 15);
    int thetaLaserDelta = fixedThetaLaser - theta32;

    int thetaLaserDeltaBot = thetaLaserDelta + zShift;
    int thetaLaserDeltaTop = thetaLaserDelta - zShift;
    int contextAngle = thetaLaserDelta >= 0 ? 0 : 1;
    if (thetaLaserDeltaTop >= 0)
      contextAngle += 2;
    else if (thetaLaserDeltaBot < 0)
      contextAngle += 2;

    delta[2] <<= 1;
    delta[2] |= _arithmeticDecoder->decode(
      _ctxPlanarPlaneLastIndexAngularIdcm[contextAngle]);
    posXyz[2] = posz0 + (delta[2] << (bitIdxZ - 1));
  }

  return delta;
}

//-------------------------------------------------------------------------

bool
GeometryOctreeDecoder::decodeIsIdcm()
{
  return _arithmeticDecoder->decode(_ctxBlockSkipTh);
}

//-------------------------------------------------------------------------
// Direct coding of position of points in node (early tree termination).
// Decoded points are written to @outputPoints
// Returns the number of points emitted.

template<class OutputIt>
int
GeometryOctreeDecoder::decodeDirectPosition(
  bool geom_unique_points_flag,
  bool joint_2pt_idcm_enabled_flag,
  const Vec3<int>& nodeSizeLog2,
  PCCOctree3Node& node,
  OctreeNodePlanar& planar,
  OutputIt outputPoints,
  bool angularIdcm,
  const Vec3<int>& headPos,
  const int* zLaser,
  const int* thetaLaser,
  int numLasers)
{
  int numPoints = 1;
  bool numPointsGt1 = _arithmeticDecoder->decode(_ctxNumIdcmPointsGt1);
  numPoints += numPointsGt1;

  int numDuplicatePoints = 0;
  if (!geom_unique_points_flag && !numPointsGt1) {
    numDuplicatePoints = !_arithmeticDecoder->decode(_ctxSinglePointPerBlock);
    if (numDuplicatePoints) {
      bool singleDup = _arithmeticDecoder->decode(_ctxSingleIdcmDupPoint);
      if (!singleDup)
        numDuplicatePoints +=
          1 + _arithmeticDecoder->decodeExpGolomb(0, _ctxPointCountPerBlock);
    }
  }

  // update node size after planar and determine upper part of position from planar
  Vec3<int32_t> deltaPlanar{0, 0, 0};
  Vec3<int> nodeSizeLog2Rem = nodeSizeLog2;
  for (int k = 0; k < 3; k++)
    if (nodeSizeLog2Rem[k] > 0 && (planar.planarMode & (1 << k))) {
      deltaPlanar[k] |= (planar.planePosBits & (1 << k) ? 1 : 0);
      nodeSizeLog2Rem[k]--;
    }

  // Indicates which components are directly coded, or coded using angular
  // contextualisation.
  Vec3<bool> directIdcm = !angularIdcm;
  point_t posNodeLidar;

  if (angularIdcm) {
    posNodeLidar =
      point_t(
        node.pos[0] << nodeSizeLog2[0], node.pos[1] << nodeSizeLog2[1],
        node.pos[2] << nodeSizeLog2[2])
      - headPos;
    bool codeXorY = std::abs(posNodeLidar[0]) <= std::abs(posNodeLidar[1]);
    directIdcm.x() = !codeXorY;
    directIdcm.y() = codeXorY;
  }

  // decode nonordred two points
  Vec3<int32_t> deltaPos[2];
  deltaPos[0] = deltaPlanar;
  deltaPos[1] = deltaPlanar;
  if (numPoints == 2 && joint_2pt_idcm_enabled_flag)
    decodeOrdered2ptPrefix(directIdcm, nodeSizeLog2Rem, deltaPos);

  if (angularIdcm) {
    for (int idx = 0; idx < 3; ++idx) {
      int N = nodeSizeLog2[idx] - nodeSizeLog2Rem[idx];
      for (int mask = N ? 1 << (N - 1) : 0; mask; mask >>= 1) {
        if (deltaPos[0][idx] & mask)
          posNodeLidar[idx] += mask << nodeSizeLog2Rem[idx];
      }
      if (nodeSizeLog2Rem[idx])
        posNodeLidar[idx] += 1 << (nodeSizeLog2Rem[idx] - 1);
    }
    node.laserIndex = findLaser(posNodeLidar, thetaLaser, numLasers);
  }

  Vec3<int32_t> pos;
  for (int i = 0; i < numPoints; i++) {
    if (angularIdcm) {
      *(outputPoints++) = pos = decodePointPositionAngular(
        nodeSizeLog2, nodeSizeLog2Rem, node, planar, headPos, zLaser,
        thetaLaser, deltaPos[i]);
    } else
      *(outputPoints++) = pos =
        decodePointPosition(nodeSizeLog2Rem, deltaPos[i]);
  }

  for (int i = 0; i < numDuplicatePoints; i++)
    *(outputPoints++) = pos;

  return numPoints + numDuplicatePoints;
}

//-------------------------------------------------------------------------

int
GeometryOctreeDecoder::decodeThetaRes()
{
  if (_arithmeticDecoder->decode(_ctxThetaResIsZero))
    return 0;

  bool sign = _arithmeticDecoder->decode(_ctxThetaResSign);
  int thetaRes = _arithmeticDecoder->decode(_ctxThetaResIsOne) ? 1 : 2;
  if (thetaRes == 2)
    thetaRes += _arithmeticDecoder->decode(_ctxThetaResIsTwo) ? 0 : 1;
  if (thetaRes == 3)
    thetaRes += _arithmeticDecoder->decodeExpGolomb(1, _ctxThetaResExp);
  return sign ? thetaRes : -thetaRes;
}

//-------------------------------------------------------------------------
// Helper to inverse quantise positions

Vec3<int32_t>
invQuantPosition(int qp, Vec3<uint32_t> quantMasks, const Vec3<int32_t>& pos)
{
  // pos represents the position within the coded tree as follows:
  //     |pppppqqqqqq|00
  //  - p = unquantised bit
  //  - q = quantised bit
  //  - 0 = bits that were not coded (MSBs of q)
  // The reconstruction is:
  //   |ppppp00qqqqqq| <- just prior to scaling
  //   |pppppssssssss| <  after scaling (s = scale(q))

  QuantizerGeom quantizer(qp);
  int shiftBits = QuantizerGeom::qpShift(qp);
  Vec3<int32_t> recon;
  for (int k = 0; k < 3; k++) {
    int lowPart = pos[k] & (quantMasks[k] >> shiftBits);
    int highPart = pos[k] ^ lowPart;
    int lowPartScaled = PCCClip(quantizer.scale(lowPart), 0, quantMasks[k]);
    recon[k] = (highPart << shiftBits) | lowPartScaled;
  }

  return recon;
}

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyDecoder>>& arithmeticDecoders,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining)
{
  // init main fifo
  //  -- worst case size is the last level containing every input poit
  //     and each point being isolated in the previous level.
  // NB: some trisoup configurations can generate fewer points than
  //     octree nodes.  Blindly trusting the number of points to guide
  //     the ringbuffer size is problematic.
  // todo(df): derive buffer size from level limit
  size_t ringBufferSize = gbh.footer.geom_num_points_minus1 + 1;
  if (gbh.trisoup_node_size_log2)
    ringBufferSize = 1100000;
  pcc::ringbuf<PCCOctree3Node> fifo(ringBufferSize + 1);

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(0);
  node00.pos = int32_t(0);
  node00.numSiblingsPlus1 = 8;
  node00.siblingOccupancy = 0;
  node00.qp = 0;
  node00.idcmEligible = 0;

  size_t processedPointCount = 0;
  std::vector<uint32_t> values;

  const int idcmThreshold = gps.geom_planar_mode_enabled_flag
    ? gps.geom_planar_idcm_threshold * 127 * 127
    : 127 * 127 * 127;

  // Lidar angles for planar prediction
  const int numLasers = gps.geom_angular_num_lidar_lasers();
  const int* thetaLaser = gps.geom_angular_theta_laser.data();
  const int* zLaser = gps.geom_angular_z_laser.data();
  const int* numPhi = gps.geom_angular_num_phi_per_turn.data();

  // Lidar position relative to slice origin
  auto headPos = gps.geomAngularOrigin - gbh.geomBoxOrigin;

  int deltaAngle = 128 << 18;
  for (int i = 0; i < numLasers - 1; i++) {
    int d = std::abs(thetaLaser[i] - thetaLaser[i + 1]);
    if (deltaAngle > d) {
      deltaAngle = d;
    }
  }

  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2) {
    occupancyAtlas.resize(
      gps.adjacent_child_contextualization_enabled_flag,
      gps.neighbour_avail_boundary_log2);
    occupancyAtlas.clear();
  }

  Vec3<uint32_t> posQuantBitMasks = 0xffffffff;
  int idcmQp = 0;
  int sliceQp = gbh.sliceQp(gps);
  int numLvlsUntilQpOffset = 0;
  if (gps.geom_scaling_enabled_flag)
    numLvlsUntilQpOffset = gbh.geom_octree_qp_offset_depth + 1;

  // generate the list of the node size for each level in the tree
  //  - starts with the smallest node and works up
  std::vector<Vec3<int>> lvlNodeSizeLog2{gbh.trisoup_node_size_log2};
  for (auto split : inReverse(gbh.tree_lvl_coded_axis_list)) {
    Vec3<int> splitStv = {!!(split & 4), !!(split & 2), !!(split & 1)};
    lvlNodeSizeLog2.push_back(lvlNodeSizeLog2.back() + splitStv);
  }
  std::reverse(lvlNodeSizeLog2.begin(), lvlNodeSizeLog2.end());

  // Derived parameter used by trisoup.
  gbh.maxRootNodeDimLog2 = lvlNodeSizeLog2[0].max();

  // the termination depth of the octree phase
  // NB: minNodeSizeLog2 is only non-zero for partial decoding (not trisoup)
  int maxDepth = lvlNodeSizeLog2.size() - skipLastLayers - 1;

  // append a dummy entry to the list so that depth+2 access is always valid
  lvlNodeSizeLog2.emplace_back(lvlNodeSizeLog2.back());

  // NB: this needs to be after the root node size is determined to
  //     allocate the planar buffer
  auto arithmeticDecoderIt = arithmeticDecoders.begin();
  GeometryOctreeDecoder decoder(gps, gbh, ctxtMem, arithmeticDecoderIt->get());

  // saved state for use with parallel bistream coding.
  // the saved state is restored at the start of each parallel octree level
  std::unique_ptr<GeometryOctreeDecoder> savedState;

  // The number of nodes to wait before updating the planar rate.
  // This is to match the prior behaviour where planar is updated once
  // per coded occupancy.
  int nodesBeforePlanarUpdate = 1;

  for (int depth = 0; depth < maxDepth; depth++) {
    // setup at the start of each level
    auto fifoCurrLvlEnd = fifo.end();
    int numNodesNextLvl = 0;
    Vec3<int32_t> occupancyAtlasOrigin = 0xffffffff;

    // derive per-level node size related parameters
    auto nodeSizeLog2 = lvlNodeSizeLog2[depth];
    auto childSizeLog2 = lvlNodeSizeLog2[depth + 1];
    // represents the largest dimension of the current node
    int nodeMaxDimLog2 = nodeSizeLog2.max();

    // if one dimension is not split, atlasShift[k] = 0
    int codedAxesPrevLvl = depth ? gbh.tree_lvl_coded_axis_list[depth - 1] : 7;
    int codedAxesCurLvl = gbh.tree_lvl_coded_axis_list[depth];

    // Idcm quantisation applies to child nodes before per node qps
    if (--numLvlsUntilQpOffset > 0) {
      // If planar is enabled, the planar bits are not quantised (since
      // the planar mode is determined before quantisation)
      auto quantNodeSizeLog2 = nodeSizeLog2;
      if (gps.geom_planar_mode_enabled_flag)
        quantNodeSizeLog2 -= 1;

      for (int k = 0; k < 3; k++)
        quantNodeSizeLog2[k] = std::max(0, quantNodeSizeLog2[k]);

      // limit the idcmQp such that it cannot overquantise the node
      auto minNs = quantNodeSizeLog2.min();
      idcmQp = gps.geom_base_qp + gps.geom_idcm_qp_offset;
      idcmQp <<= gps.geom_qp_multiplier_log2;
      idcmQp = std::min(idcmQp, minNs * 8);

      for (int k = 0; k < 3; k++)
        posQuantBitMasks[k] = (1 << quantNodeSizeLog2[k]) - 1;
    }

    // record the node size when quantisation is signalled -- all subsequnt
    // coded occupancy bits are quantised
    // after the qp offset, idcm nodes do not receive special treatment
    if (!numLvlsUntilQpOffset) {
      idcmQp = 0;
      for (int k = 0; k < 3; k++)
        posQuantBitMasks[k] = (1 << nodeSizeLog2[k]) - 1;
    }

    // save context state for parallel coding
    if (depth == maxDepth - 1 - gbh.geom_stream_cnt_minus1)
      if (gbh.geom_stream_cnt_minus1)
        savedState.reset(new GeometryOctreeDecoder(decoder));

    // load context state for parallel coding starting one level later
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      decoder = *savedState;
      decoder._arithmeticDecoder = (++arithmeticDecoderIt)->get();
    }

    auto planarDepth = lvlNodeSizeLog2[0] - nodeSizeLog2;
    decoder.beginOctreeLevel(planarDepth);

    // process all nodes within a single level
    for (; fifo.begin() != fifoCurrLvlEnd; fifo.pop_front()) {
      PCCOctree3Node& node0 = fifo.front();

      if (numLvlsUntilQpOffset == 0) {
        node0.qp = sliceQp;
        node0.qp += decoder.decodeQpOffset() << gps.geom_qp_multiplier_log2;
      }

      int shiftBits = QuantizerGeom::qpShift(node0.qp);
      auto effectiveNodeSizeLog2 = nodeSizeLog2 - shiftBits;
      auto effectiveChildSizeLog2 = childSizeLog2 - shiftBits;

      // make quantisation work with qtbt and planar.
      auto codedAxesCurNode = codedAxesCurLvl;
      if (shiftBits != 0) {
        for (int k = 0; k < 3; k++) {
          if (effectiveChildSizeLog2[k] < 0)
            codedAxesCurNode &= ~(4 >> k);
        }
      }

      GeometryNeighPattern gnp{};
      if (gps.neighbour_avail_boundary_log2) {
        updateGeometryOccupancyAtlas(
          node0.pos, codedAxesPrevLvl, fifo, fifoCurrLvlEnd, &occupancyAtlas,
          &occupancyAtlasOrigin);

        gnp = makeGeometryNeighPattern(
          gps.adjacent_child_contextualization_enabled_flag, node0.pos,
          codedAxesPrevLvl, codedAxesCurLvl, occupancyAtlas);
      } else {
        // The position of the node in the parent's occupancy map
        int posInParent = 0;
        posInParent |= (node0.pos[0] & 1) << 2;
        posInParent |= (node0.pos[1] & 1) << 1;
        posInParent |= (node0.pos[2] & 1) << 0;
        posInParent &= codedAxesPrevLvl;

        gnp.neighPattern =
          neighPatternFromOccupancy(posInParent, node0.siblingOccupancy);
      }

      int contextAngle = -1;
      int contextAnglePhiX = -1;
      int contextAnglePhiY = -1;
      if (gps.geom_angular_mode_enabled_flag) {
        contextAngle = determineContextAngleForPlanar(
          node0, headPos, nodeSizeLog2, zLaser, thetaLaser, numLasers,
          deltaAngle, decoder._phiZi, decoder._phiBuffer.data(),
          &contextAnglePhiX, &contextAnglePhiY);
      }

      if (gps.geom_planar_mode_enabled_flag) {
        // update the plane rate depending on the occupancy and local density
        auto occupancy = node0.siblingOccupancy;
        auto numOccupied = node0.numSiblingsPlus1;
        if (!nodesBeforePlanarUpdate--) {
          decoder._planar.updateRate(occupancy, numOccupied);
          nodesBeforePlanarUpdate = numOccupied - 1;
        }
      }

      OctreeNodePlanar planar;
      if (!isLeafNode(effectiveNodeSizeLog2)) {
        // planar eligibility
        bool planarEligible[3] = {false, false, false};
        if (gps.geom_planar_mode_enabled_flag) {
          decoder._planar.isEligible(planarEligible);
          if (gps.geom_angular_mode_enabled_flag) {
            if (contextAngle != -1)
              planarEligible[2] = true;
            planarEligible[0] = (contextAnglePhiX != -1);
            planarEligible[1] = (contextAnglePhiY != -1);
          }

          for (int k = 0; k < 3; k++)
            planarEligible[k] &= (codedAxesCurNode >> (2 - k)) & 1;
        }

        int planarProb[3] = {127, 127, 127};
        decoder.determinePlanarMode(
          planarEligible, node0, planar, gnp.neighPattern, planarProb,
          contextAngle, contextAnglePhiX, contextAnglePhiY);

        node0.idcmEligible &=
          planarProb[0] * planarProb[1] * planarProb[2] <= idcmThreshold;
      }

      // At the scaling depth, it is possible for a node that has previously
      // been marked as being eligible for idcm to be fully quantised due
      // to the choice of QP.  There is therefore nothing to code with idcm.
      if (isLeafNode(effectiveNodeSizeLog2))
        node0.idcmEligible = false;

      if (node0.idcmEligible) {
        bool isDirectMode = decoder.decodeIsIdcm();
        if (isDirectMode) {
          auto idcmSize = effectiveNodeSizeLog2;
          if (idcmQp) {
            node0.qp = idcmQp;
            idcmSize = nodeSizeLog2 - QuantizerGeom::qpShift(idcmQp);
          }

          int numPoints = decoder.decodeDirectPosition(
            gps.geom_unique_points_flag, gps.joint_2pt_idcm_enabled_flag,
            idcmSize, node0, planar, &pointCloud[processedPointCount],
            gps.geom_angular_mode_enabled_flag, headPos, zLaser, thetaLaser,
            numLasers);

          for (int j = 0; j < numPoints; j++) {
            auto& point = pointCloud[processedPointCount++];
            for (int k = 0; k < 3; k++)
              point[k] += rotateLeft(node0.pos[k], idcmSize[k]);

            point = invQuantPosition(node0.qp, posQuantBitMasks, point);
          }

          // NB: no further siblings to decode by definition of IDCM
          if (gps.inferred_direct_coding_mode <= 1)
            assert(node0.numSiblingsPlus1 == 1);

          // This node has no children, ensure that future nodes avoid
          // accessing stale child occupancy data.
          if (gps.adjacent_child_contextualization_enabled_flag)
            updateGeometryOccupancyAtlasOccChild(
              node0.pos, 0, &occupancyAtlas);

          continue;
        }
      }

      int occupancyIsPredicted = 0;
      int occupancyPrediction = 0;

      // generate intra prediction
      if (
        nodeMaxDimLog2 < gps.intra_pred_max_node_size_log2
        && gps.neighbour_avail_boundary_log2 > 0) {
        predictGeometryOccupancyIntra(
          occupancyAtlas, node0.pos, codedAxesPrevLvl, &occupancyIsPredicted,
          &occupancyPrediction);
      }

      uint8_t occupancy = 1;
      if (!isLeafNode(effectiveNodeSizeLog2)) {
        // planar mode for current node
        // mask to be used for the occupancy coding
        // (bit =1 => occupancy bit not coded due to not belonging to the plane)
        int mask_planar[3] = {0, 0, 0};
        maskPlanar(planar, mask_planar, codedAxesCurNode);

        occupancy = decoder.decodeOccupancy(
          gnp, occupancyIsPredicted, occupancyPrediction, mask_planar[0],
          mask_planar[1], mask_planar[2], planar.planarPossible & 1,
          planar.planarPossible & 2, planar.planarPossible & 4);
      }

      assert(occupancy > 0);

      // update atlas for child neighbours
      // NB: the child occupancy atlas must be updated even if the current
      //     node has no occupancy coded in order to clear any stale state in
      //     the atlas.
      if (gps.adjacent_child_contextualization_enabled_flag)
        updateGeometryOccupancyAtlasOccChild(
          node0.pos, occupancy, &occupancyAtlas);

      // population count of occupancy for IDCM
      int numOccupied = popcnt(occupancy);

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
          Vec3<int32_t> point{(node0.pos[0] << !!(codedAxesCurLvl & 4)) + x,
                              (node0.pos[1] << !!(codedAxesCurLvl & 2)) + y,
                              (node0.pos[2] << !!(codedAxesCurLvl & 1)) + z};

          // remove any padding bits that were not coded
          for (int k = 0; k < 3; k++)
            point[k] = rotateLeft(point[k], effectiveChildSizeLog2[k]);

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
        child.pos[0] = (node0.pos[0] << !!(codedAxesCurLvl & 4)) + x;
        child.pos[1] = (node0.pos[1] << !!(codedAxesCurLvl & 2)) + y;
        child.pos[2] = (node0.pos[2] << !!(codedAxesCurLvl & 1)) + z;
        child.numSiblingsPlus1 = numOccupied;
        child.siblingOccupancy = occupancy;
        child.laserIndex = node0.laserIndex;

        child.idcmEligible = isDirectModeEligible(
          gps.inferred_direct_coding_mode, nodeMaxDimLog2, gnp.neighPattern,
          node0, child);

        numNodesNextLvl++;
      }
    }

    // Check that one level hasn't produced too many nodes
    // todo(df): this check is too weak to spot overflowing the fifo
    assert(numNodesNextLvl <= ringBufferSize);
  }

  // save the context state for re-use by a future slice if required
  ctxtMem = decoder.getCtx();

  // NB: the point cloud needs to be resized if partially decoded
  // OR: if geometry quantisation has changed the number of points
  pointCloud.resize(processedPointCount);

  // return partial coding result
  //  - add missing levels to node positions and inverse quantise
  if (nodesRemaining) {
    auto nodeSizeLog2 = lvlNodeSizeLog2[maxDepth];
    for (auto& node : fifo) {
      for (int k = 0; k < 3; k++)
        node.pos[k] <<= nodeSizeLog2[k] - QuantizerGeom::qpShift(node.qp);
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
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyDecoder>>& arithmeticDecoders)
{
  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMem, arithmeticDecoders, nullptr);
}

//-------------------------------------------------------------------------

void
decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyDecoder>>& arithmeticDecoders)
{
  pcc::ringbuf<PCCOctree3Node> nodes;
  decodeGeometryOctree(
    gps, gbh, minGeomNodeSizeLog2, pointCloud, ctxtMem, arithmeticDecoders,
    &nodes);

  if (minGeomNodeSizeLog2 > 0) {
    size_t size =
      pointCloud.removeDuplicatePointInQuantizedPoint(minGeomNodeSizeLog2);

    pointCloud.resize(size + nodes.size());
    size_t processedPointCount = size;

    if (minGeomNodeSizeLog2 > 1) {
      uint32_t mask = uint32_t(-1) << minGeomNodeSizeLog2;
      for (auto node0 : nodes) {
        for (int k = 0; k < 3; k++)
          node0.pos[k] &= mask;
        node0.pos += 1 << (minGeomNodeSizeLog2 - 1);
        pointCloud[processedPointCount++] = node0.pos;
      }
    } else {
      for (auto node0 : nodes)
        pointCloud[processedPointCount++] = node0.pos;
    }
  }
}

//============================================================================

}  // namespace pcc
