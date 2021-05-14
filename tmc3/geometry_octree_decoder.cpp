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
    int adjPlanes,
    int planeId,
    int contextAngle);

  void determinePlanarMode(
    bool adjacent_child_contextualization_enabled_flag,
    int planeId,
    OctreeNodePlanar& child,
    OctreePlanarBuffer::Row* planeBuffer,
    int coord1,
    int coord2,
    int coord3,
    int posInParent,
    const GeometryNeighPattern& gnp,
    uint8_t siblingOccupancy,
    int planarRate[3],
    int contextAngle);

  void determinePlanarMode(
    bool adjacent_child_contextualization_enabled_flag,
    const bool planarEligible[3],
    int posInParent,
    const GeometryNeighPattern& gnp,
    PCCOctree3Node& child,
    OctreeNodePlanar& planar,
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
    const OctreeAngPosScaler& quant,
    const Vec3<int>& nodeSizeLog2Rem,
    const Vec3<int>& angularOrigin,
    const int* zLaser,
    const int* thetaLaser,
    int laserIdx,
    const Vec3<int>& nodePos,
    Vec3<int> posXyz,
    Vec3<int> delta);

  bool decodeNodeQpOffsetsPresent();
  int decodeQpOffset();

  bool decodeIsIdcm();

  template<class OutputIt>
  int decodeDirectPosition(
    bool geom_unique_points_flag,
    bool joint_2pt_idcm_enabled_flag,
    bool geom_angular_mode_enabled_flag,
    const Vec3<int>& nodeSizeLog2,
    const Vec3<int>& posQuantBitMasks,
    const PCCOctree3Node& node,
    const OctreeNodePlanar& planar,
    const Vec3<int>& headPos,
    const int* zLaser,
    const int* thetaLaser,
    int numLasers,
    OutputIt outputPoints);

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
  , _phiBuffer(gps.numLasers(), 0x80000000)
  , _phiZi(gps.numLasers(), gps.angularNumPhiPerTurn)
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
  int val = _arithmeticDecoder->decode(_ctxDupPointCntGt0);
  if (val)
    val += _arithmeticDecoder->decodeExpGolomb(0, _ctxDupPointCntEgl);

  return val + 1;
}

//============================================================================

int
GeometryOctreeDecoder::decodePlanarMode(
  OctreeNodePlanar& planar,
  int planeZ,
  int dist,
  int adjPlanes,
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
    static const int kAdjPlaneCtx[4] = {0, 1, 2, 0};
    int planePosCtx = kAdjPlaneCtx[adjPlanes];
    if (planeZ < 0) {
      planeBit =
        _arithmeticDecoder->decode(_ctxPlanarPlaneLastIndexZ[planePosCtx]);
    } else {
      int discreteDist = dist > (8 >> OctreePlanarBuffer::shiftAb);
      int lastIndexPlane2d = planeZ + (discreteDist << 1);

      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndex[planeId][planePosCtx][lastIndexPlane2d]);
    }
  } else {               // angular mode on
    if (planeId == 2) {  // angular
      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndexAngular[contextAngle]);
    } else {  // azimuthal
      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndexAngularPhi[contextAngle]);
    }
  }

  planar.planePosBits |= (planeBit << planeId);
  return planeBit;
}

//============================================================================

void
GeometryOctreeDecoder::determinePlanarMode(
  bool adjacent_child_contextualization_enabled_flag,
  int planeId,
  OctreeNodePlanar& planar,
  OctreePlanarBuffer::Row* planeBuffer,
  int coord1,
  int coord2,
  int coord3,
  int posInParent,
  const GeometryNeighPattern& gnp,
  uint8_t siblingOccupancy,
  int planarRate[3],
  int contextAngle)
{
  const int kPlanarChildThreshold = 63;
  const int kAdjNeighIdxFromPlanePos[3][2] = {1, 0, 2, 3, 4, 5};
  const int planeSelector = 1 << planeId;
  static const uint8_t KAdjNeighIdxMask[3][2] = {0x0f, 0xf0, 0x33,
                                                 0xcc, 0x55, 0xaa};
  OctreePlanarBuffer::Elmt* row;
  int rowLen = OctreePlanarBuffer::rowSize;
  int closestPlanarFlag;
  int closestDist;
  int maxCoord;

  if (!planeBuffer) {
    // angular: buffer disabled
    closestPlanarFlag = -1;
    closestDist = 0;
  } else {
    coord1 =
      (coord1 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
    coord2 =
      (coord2 & OctreePlanarBuffer::maskAb) >> OctreePlanarBuffer::shiftAb;
    coord3 = coord3 & OctreePlanarBuffer::maskC;

    row = planeBuffer[coord3];

    maxCoord = std::max(coord1, coord2);
    closestDist = std::abs(maxCoord - int(row[rowLen - 1].pos));
    int idxMinDist = rowLen - 1;

    // push closest point front
    row[rowLen - 1] = row[idxMinDist];

    closestPlanarFlag = row[idxMinDist].planeIdx;
  }

  // The relative plane position (0|1) along the planeId axis.
  int pos = !(KAdjNeighIdxMask[planeId][0] & (1 << posInParent));

  // Determine which adjacent planes are occupied
  // The low plane is at position axis - 1
  bool lowAdjPlaneOccupied = adjacent_child_contextualization_enabled_flag
    ? KAdjNeighIdxMask[planeId][1] & gnp.adjNeighOcc[planeId]
    : (gnp.neighPattern >> kAdjNeighIdxFromPlanePos[planeId][0]) & 1;

  // The high adjacent plane is at position axis + 1
  bool highAdjPlaneOccupied = !pos
    ? KAdjNeighIdxMask[planeId][1] & siblingOccupancy
    : (gnp.neighPattern >> kAdjNeighIdxFromPlanePos[planeId][1]) & 1;

  int adjPlanes = (highAdjPlaneOccupied << 1) | lowAdjPlaneOccupied;
  int planeBit = decodePlanarMode(
    planar, closestPlanarFlag, closestDist, adjPlanes, planeId, contextAngle);

  bool isPlanar = (planar.planarMode & planeSelector);

  planarRate[planeId] =
    (255 * planarRate[planeId] + (isPlanar ? 256 * 8 : 0) + 128) >> 8;

  if (planeBuffer)
    row[rowLen - 1] = {unsigned(maxCoord), planeBit};
}

//============================================================================

void
GeometryOctreeDecoder::determinePlanarMode(
  bool adjacent_child_contextualization_enabled_flag,
  const bool planarEligible[3],
  int posInParent,
  const GeometryNeighPattern& gnp,
  PCCOctree3Node& child,
  OctreeNodePlanar& planar,
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
      adjacent_child_contextualization_enabled_flag, 0, planar,
      planeBuffer.getBuffer(0), yy, zz, xx, posInParent, gnp,
      child.siblingOccupancy, _planar._rate.data(), contextAnglePhiX);
  }
  // planar y
  if (planarEligible[1]) {
    determinePlanarMode(
      adjacent_child_contextualization_enabled_flag, 1, planar,
      planeBuffer.getBuffer(1), xx, zz, yy, posInParent, gnp,
      child.siblingOccupancy, _planar._rate.data(), contextAnglePhiY);
  }
  // planar z
  if (planarEligible[2]) {
    determinePlanarMode(
      adjacent_child_contextualization_enabled_flag, 2, planar,
      planeBuffer.getBuffer(2), xx, yy, zz, posInParent, gnp,
      child.siblingOccupancy, _planar._rate.data(), contextAngle);
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

bool
GeometryOctreeDecoder::decodeNodeQpOffsetsPresent()
{
  return _arithmeticDecoder->decode();
}

//-------------------------------------------------------------------------

int
GeometryOctreeDecoder::decodeQpOffset()
{
  if (!_arithmeticDecoder->decode(_ctxQpOffsetAbsGt0))
    return 0;

  int dqp = _arithmeticDecoder->decodeExpGolomb(0, _ctxQpOffsetAbsEgl) + 1;
  int dqp_sign = _arithmeticDecoder->decode(_ctxQpOffsetSign);
  return dqp_sign ? -dqp : dqp;
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
  const OctreeAngPosScaler& quant,
  const Vec3<int>& nodeSizeLog2Rem,
  const Vec3<int>& angularOrigin,
  const int* zLaser,
  const int* thetaLaser,
  int laserIdx,
  const Vec3<int>& nodePos,
  Vec3<int> posXyz,
  Vec3<int> delta)
{
  // -- PHI --
  // code x or y directly and compute phi of node
  bool directAxis = std::abs(posXyz[0]) <= std::abs(posXyz[1]);
  for (int i = nodeSizeLog2Rem[directAxis]; i > 0; i--) {
    delta[directAxis] <<= 1;
    delta[directAxis] |= _arithmeticDecoder->decode();
  }

  posXyz += quant.scaleEns(delta << nodeSizeLog2Rem);
  posXyz[directAxis] =
    quant.scaleEns(directAxis, nodePos[directAxis] + delta[directAxis])
    - angularOrigin[directAxis];

  // laser residual
  laserIdx += decodeThetaRes();

  // find predictor
  int phiNode = iatan2(posXyz[1], posXyz[0]);
  int predPhi = _phiBuffer[laserIdx];
  if (predPhi == 0x80000000)
    predPhi = phiNode;

  // elementary shift predictor
  int nShift =
    ((predPhi - phiNode) * _phiZi.invDelta(laserIdx) + (1 << 29)) >> 30;
  predPhi -= _phiZi.delta(laserIdx) * nShift;

  // azimuthal code x or y
  const int phiAxis = !directAxis;
  for (int mask = (1 << nodeSizeLog2Rem[phiAxis]) >> 1; mask; mask >>= 1) {
    // angles left and right
    int scaledMask = quant.scaleEns(phiAxis, mask);
    int phiL = phiNode;
    int phiR = directAxis ? iatan2(posXyz[1], posXyz[0] + scaledMask)
                          : iatan2(posXyz[1] + scaledMask, posXyz[0]);

    // ctx azimutal
    int angleL = phiL - predPhi;
    int angleR = phiR - predPhi;
    int contextAnglePhi =
      (angleL >= 0 && angleR >= 0) || (angleL < 0 && angleR < 0) ? 2 : 0;
    angleL = std::abs(angleL);
    angleR = std::abs(angleR);
    if (angleL > angleR) {
      contextAnglePhi++;
      std::swap(angleL, angleR);
    }
    if (angleR > (angleL << 1))
      contextAnglePhi += 4;

    // entropy coding
    auto& ctx = _ctxPlanarPlaneLastIndexAngularPhiIDCM[contextAnglePhi];
    bool bit = _arithmeticDecoder->decode(ctx);
    delta[phiAxis] <<= 1;
    if (bit) {
      delta[phiAxis] |= 1;
      posXyz[phiAxis] += scaledMask;
      phiNode = phiR;
      predPhi = _phiBuffer[laserIdx];
      if (predPhi == 0x80000000)
        predPhi = phiNode;

      // elementary shift predictor
      int nShift =
        ((predPhi - phiNode) * _phiZi.invDelta(laserIdx) + (1 << 29)) >> 30;
      predPhi -= _phiZi.delta(laserIdx) * nShift;
    }
  }

  // update buffer phi
  _phiBuffer[laserIdx] = phiNode;

  // -- THETA --
  int maskz = (1 << nodeSizeLog2Rem[2]) >> 1;
  if (!maskz)
    return delta;

  // Since x and y are known,
  // r is known too and does not depend on the bit for z
  uint64_t xLidar = (int64_t(posXyz[0]) << 8) - 128;
  uint64_t yLidar = (int64_t(posXyz[1]) << 8) - 128;
  uint64_t r2 = xLidar * xLidar + yLidar * yLidar;
  int64_t rInv = irsqrt(r2);

  // code bits for z using angular. Eligility is implicit. Laser is known.
  int64_t hr = zLaser[laserIdx] * rInv;
  int fixedThetaLaser =
    thetaLaser[laserIdx] + int(hr >= 0 ? -(hr >> 17) : ((-hr) >> 17));

  int zShift = (rInv * quant.scaleEns(2, 1 << nodeSizeLog2Rem[2])) >> 18;
  for (int bitIdxZ = nodeSizeLog2Rem[2]; bitIdxZ > 0;
       bitIdxZ--, maskz >>= 1, zShift >>= 1) {
    // determine non-corrected theta
    int scaledMaskZ = quant.scaleEns(2, maskz);
    int64_t zLidar = ((posXyz[2] + scaledMaskZ) << 1) - 1;
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

    auto& ctx = _ctxPlanarPlaneLastIndexAngularIdcm[contextAngle];
    delta[2] <<= 1;
    delta[2] |= _arithmeticDecoder->decode(ctx);
    if (delta[2] & 1)
      posXyz[2] += scaledMaskZ;
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
  bool geom_angular_mode_enabled_flag,
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& posQuantBitMask,
  const PCCOctree3Node& node,
  const OctreeNodePlanar& planar,
  const Vec3<int>& angularOrigin,
  const int* zLaser,
  const int* thetaLaser,
  int numLasers,
  OutputIt outputPoints)
{
  int numPoints = 1;
  bool numPointsGt1 = _arithmeticDecoder->decode(_ctxNumIdcmPointsGt1);
  numPoints += numPointsGt1;

  int numDuplicatePoints = 0;
  if (!geom_unique_points_flag && !numPointsGt1) {
    numDuplicatePoints = _arithmeticDecoder->decode(_ctxDupPointCntGt0);
    if (numDuplicatePoints) {
      numDuplicatePoints += _arithmeticDecoder->decode(_ctxDupPointCntGt1);
      if (numDuplicatePoints == 2)
        numDuplicatePoints +=
          _arithmeticDecoder->decodeExpGolomb(0, _ctxDupPointCntEgl);
    }
  }

  // nodeSizeLog2Rem indicates the number of bits left to decode
  // the first bit may be inferred from the planar information
  Vec3<int32_t> deltaPlanar{0, 0, 0};
  Vec3<int> nodeSizeLog2Rem = nodeSizeLog2;
  for (int k = 0; k < 3; k++)
    if (nodeSizeLog2Rem[k] > 0 && (planar.planarMode & (1 << k))) {
      deltaPlanar[k] |= (planar.planePosBits & (1 << k) ? 1 : 0);
      nodeSizeLog2Rem[k]--;
    }

  // quantised partial positions must be scaled for angular coding
  // nb, the decoded position remains quantised.
  OctreeAngPosScaler quant(node.qp, posQuantBitMask);

  // Indicates which components are directly coded
  Vec3<bool> directIdcm = true;

  // Position of the node relative to the angular origin
  point_t posNodeLidar;
  if (geom_angular_mode_enabled_flag) {
    posNodeLidar = quant.scaleEns(node.pos << nodeSizeLog2) - angularOrigin;
    bool directAxis = std::abs(posNodeLidar[0]) <= std::abs(posNodeLidar[1]);
    directIdcm = false;
    directIdcm[directAxis] = true;
  }

  // decode (ordred) two points
  Vec3<int32_t> deltaPos[2] = {deltaPlanar, deltaPlanar};
  if (numPoints == 2 && joint_2pt_idcm_enabled_flag)
    decodeOrdered2ptPrefix(directIdcm, nodeSizeLog2Rem, deltaPos);

  int laserIdx;
  if (geom_angular_mode_enabled_flag) {
    auto delta = (deltaPos[0] << nodeSizeLog2Rem);
    delta += (1 << nodeSizeLog2Rem) >> 1;
    delta = quant.scaleEns(delta);
    laserIdx = findLaser(posNodeLidar + delta, thetaLaser, numLasers);
  }

  Vec3<int32_t> pos;
  for (int i = 0; i < numPoints; i++) {
    if (geom_angular_mode_enabled_flag)
      *(outputPoints++) = pos = decodePointPositionAngular(
        quant, nodeSizeLog2Rem, angularOrigin, zLaser, thetaLaser, laserIdx,
        node.pos << nodeSizeLog2, posNodeLidar, deltaPos[i]);
    else
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
  if (!_arithmeticDecoder->decode(_ctxThetaRes[0]))
    return 0;

  int absVal = 1;
  absVal += _arithmeticDecoder->decode(_ctxThetaRes[1]);
  if (absVal > 1)
    absVal += _arithmeticDecoder->decode(_ctxThetaRes[2]);
  if (absVal == 3)
    absVal += _arithmeticDecoder->decodeExpGolomb(1, _ctxThetaResExp);

  bool sign = _arithmeticDecoder->decode(_ctxThetaResSign);
  return sign ? -absVal : absVal;
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
  EntropyDecoder& arithmeticDecoder,
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
  if (gbh.trisoupNodeSizeLog2(gps))
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

  // rotating mask used to enable idcm
  uint32_t idcmEnableMaskInit = mkIdcmEnableMask(gps);

  // Lidar angles for planar prediction
  const int numLasers = gps.numLasers();
  const int* thetaLaser = gps.angularTheta.data();
  const int* zLaser = gps.angularZ.data();

  // Lidar position relative to slice origin
  auto angularOrigin = gbh.geomAngularOrigin(gps);

  int deltaAngle = 128 << 18;
  for (int i = 0; i < numLasers - 1; i++) {
    int d = std::abs(thetaLaser[i] - thetaLaser[i + 1]);
    if (deltaAngle > d) {
      deltaAngle = d;
    }
  }

  MortonMap3D occupancyAtlas;
  if (gps.neighbour_avail_boundary_log2_minus1) {
    occupancyAtlas.resize(
      gps.adjacent_child_contextualization_enabled_flag,
      gps.neighbour_avail_boundary_log2_minus1 + 1);
    occupancyAtlas.clear();
  }

  Vec3<uint32_t> posQuantBitMasks = 0xffffffff;
  int idcmQp = 0;
  int sliceQp = gbh.sliceQp(gps);
  int nodeQpOffsetsSignalled = !gps.geom_scaling_enabled_flag;

  // generate the list of the node size for each level in the tree
  //  - starts with the smallest node and works up
  std::vector<Vec3<int>> lvlNodeSizeLog2{gbh.trisoupNodeSizeLog2(gps)};
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
  GeometryOctreeDecoder decoder(gps, gbh, ctxtMem, &arithmeticDecoder);

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

    // Determine if this is the level where node QPs are sent
    bool nodeQpOffsetsPresent =
      !nodeQpOffsetsSignalled && decoder.decodeNodeQpOffsetsPresent();

    // record the node size when quantisation is signalled -- all subsequnt
    // coded occupancy bits are quantised
    // after the qp offset, idcm nodes do not receive special treatment
    if (nodeQpOffsetsPresent) {
      nodeQpOffsetsSignalled = true;
      idcmQp = 0;
      posQuantBitMasks = Vec3<uint32_t>((1 << nodeSizeLog2) - 1);
    }

    // Idcm quantisation applies to child nodes before per node qps
    if (!nodeQpOffsetsSignalled) {
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

      posQuantBitMasks = Vec3<uint32_t>((1 << quantNodeSizeLog2) - 1);
    }

    // save context state for parallel coding
    if (depth == maxDepth - 1 - gbh.geom_stream_cnt_minus1)
      if (gbh.geom_stream_cnt_minus1)
        savedState.reset(new GeometryOctreeDecoder(decoder));

    // a new entropy stream starts one level after the context state is saved.
    // restore the saved state and flush the arithmetic decoder
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      decoder = *savedState;
      arithmeticDecoder.flushAndRestart();
    }

    // reset the idcm eligibility mask at the start of each level to
    // support multiple streams
    auto idcmEnableMask = rotateRight(idcmEnableMaskInit, depth);

    auto planarDepth = lvlNodeSizeLog2[0] - nodeSizeLog2;
    decoder.beginOctreeLevel(planarDepth);

    // process all nodes within a single level
    for (; fifo.begin() != fifoCurrLvlEnd; fifo.pop_front()) {
      PCCOctree3Node& node0 = fifo.front();

      if (nodeQpOffsetsPresent) {
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
      // The position of the node in the parent's occupancy map
      int posInParent = 0;
      posInParent |= (node0.pos[0] & 1) << 2;
      posInParent |= (node0.pos[1] & 1) << 1;
      posInParent |= (node0.pos[2] & 1) << 0;
      posInParent &= codedAxesPrevLvl;

      if (gps.neighbour_avail_boundary_log2_minus1) {
        updateGeometryOccupancyAtlas(
          node0.pos, codedAxesPrevLvl, fifo, fifoCurrLvlEnd, &occupancyAtlas,
          &occupancyAtlasOrigin);

        gnp = makeGeometryNeighPattern(
          gps.adjacent_child_contextualization_enabled_flag, node0.pos,
          codedAxesPrevLvl, codedAxesCurLvl, occupancyAtlas);
      } else {
        gnp.neighPattern =
          neighPatternFromOccupancy(posInParent, node0.siblingOccupancy);
      }

      int contextAngle = -1;
      int contextAnglePhiX = -1;
      int contextAnglePhiY = -1;
      if (gps.geom_angular_mode_enabled_flag) {
        contextAngle = determineContextAngleForPlanar(
          node0, nodeSizeLog2, angularOrigin, zLaser, thetaLaser, numLasers,
          deltaAngle, decoder._phiZi, decoder._phiBuffer.data(),
          &contextAnglePhiX, &contextAnglePhiY, posQuantBitMasks);
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

        decoder.determinePlanarMode(
          gps.adjacent_child_contextualization_enabled_flag, planarEligible,
          posInParent, gnp, node0, planar, contextAngle, contextAnglePhiX,
          contextAnglePhiY);
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
            gps.geom_angular_mode_enabled_flag, idcmSize, posQuantBitMasks,
            node0, planar, angularOrigin, zLaser, thetaLaser, numLasers,
            &pointCloud[processedPointCount]);

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

      uint8_t occupancy = 1;
      if (!isLeafNode(effectiveNodeSizeLog2)) {
        // planar mode for current node
        // mask to be used for the occupancy coding
        // (bit =1 => occupancy bit not coded due to not belonging to the plane)
        int planarMask[3] = {0, 0, 0};
        maskPlanar(planar, planarMask, codedAxesCurNode);

        // generate intra prediction
        bool intraPredUsed = !(planarMask[0] | planarMask[1] | planarMask[2]);
        int occupancyIsPredicted = 0;
        int occupancyPrediction = 0;
        if (
          nodeMaxDimLog2 < gps.intra_pred_max_node_size_log2
          && gps.neighbour_avail_boundary_log2_minus1 > 0 && intraPredUsed) {
          predictGeometryOccupancyIntra(
            occupancyAtlas, node0.pos, codedAxesPrevLvl, &occupancyIsPredicted,
            &occupancyPrediction);
        }

        occupancy = decoder.decodeOccupancy(
          gnp, occupancyIsPredicted, occupancyPrediction, planarMask[0],
          planarMask[1], planarMask[2], planar.planarPossible & 1,
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

        if (child.idcmEligible) {
          child.idcmEligible &= idcmEnableMask & 1;
          idcmEnableMask = rotateRight(idcmEnableMask, 1);
        }

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
      node.pos <<= nodeSizeLog2 - QuantizerGeom::qpShift(node.qp);
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
  EntropyDecoder& arithmeticDecoder)
{
  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMem, arithmeticDecoder, nullptr);
}

//-------------------------------------------------------------------------

void
decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder)
{
  pcc::ringbuf<PCCOctree3Node> nodes;
  decodeGeometryOctree(
    gps, gbh, minGeomNodeSizeLog2, pointCloud, ctxtMem, arithmeticDecoder,
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
