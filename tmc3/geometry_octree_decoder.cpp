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
#include "motionWip.h"

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

  // dynamic OBUF
  void resetMap(const bool& enableInter, const bool& enablePlanar) { GeometryOctreeContexts::resetMap(enableInter, enablePlanar); }
  void clearMap(const bool& enableInter, const bool& enablePlanar) { GeometryOctreeContexts::clearMap(enableInter, enablePlanar); };

  int decodePositionLeafNumPoints();

  int decodePlanarMode(
    bool planarEligibilityDynamicOBUF,
    OctreeNodePlanar& planar,
    int planeZ,
    int dist,
    int adjPlanes,
    int planeId,
    int contextAngle,
    bool* multiPlanarFlag,
    bool* multiPlanarEligible, 
    OctreeNodePlanar& planarRef,
    const OctreeNodePlanar adjNeighPlanar[7],
    const bool& neighAvai,
    const uint32_t& neighOccu

  );

  void derivePlanarPCMContextBuffer(
    OctreeNodePlanar& planar,
    OctreeNodePlanar& planarRef,
    OctreePlanarBuffer& planeBuffer,
    int xx,
    int yy,
    int zz);

  void determinePlanarMode(
    bool planarEligibilityDynamicOBUF,
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
    int contextAngle,
    bool* multiPlanarFlag,
    bool* multiPlanarEligible,
    OctreeNodePlanar& planarRef,
    const OctreeNodePlanar adjNeighPlanar[7],
    const bool& neighAvai,
    const uint32_t& neighOccu

  );

  void determinePlanarMode(
    bool planarEligibilityDynamicOBUF,
    bool adjacent_child_contextualization_enabled_flag,
    const bool planarEligible[3],
    int posInParent,
    const GeometryNeighPattern& gnp,
    PCCOctree3Node& child,
    OctreeNodePlanar& planar,
    int contextAngle,
    int contextAnglePhiX,
    int contextAnglePhiY
    , OctreeNodePlanar& planarRef
  );

  uint32_t decodeOccupancyFullNeihbourgsNZ(
    const GeometryNeighPattern& gnp,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ,
    const MortonMap3D& occupancyAtlas,
    Vec3<int32_t> &pos,
    const int atlasShift,
    int predOcc,
    bool flagNoSingle,
    const bool& planarEligibleKOctreeDepth);

  uint32_t decodeOccupancyFullNeihbourgs(
    const GeometryNeighPattern& gnp,
    int planarMaskX,
    int planarMaskY,
    int planarMaskZ,
    bool planarPossibleX,
    bool planarPossibleY,
    bool planarPossibleZ,
    const MortonMap3D& occupancyAtlas,
    Vec3<int32_t> pos,
    const int atlasShift,
    bool flagWord4,
    bool adjacent_child_contextualization_enabled_flag
    , int predOcc,
    const bool& planarEligibleKOctreeDepth
  );

  Vec3<int32_t> decodePointPosition(
    const Vec3<int>& nodeSizeLog2, Vec3<int32_t>& deltaPlanar);

  void decodeOrdered2ptPrefix(
    Vec3<bool> directIdcm,
    Vec3<int>& nodeSizeLog2AfterUnordered,
    Vec3<int32_t> deltaUnordered[2]);

  Vec3<int32_t> decodePointPositionAngular(
    const bool& enbleInter,
    const int& predLaserIdx,
    const OctreeAngPosScaler& quant,
    const Vec3<int>& nodeSizeLog2Rem,
    const Vec3<int>& angularOrigin,
    const int* zLaser,
    const int* thetaLaser,
    int laserIdx,
    const Vec3<int>& nodePos,
    Vec3<int> posXyz,
    Vec3<int> delta);

  inline int32_t decodePointPositionZAngular(
    const OctreeAngPosScaler& quant,
    const Vec3<int>& nodeSizeLog2Rem,
    const int* zLaser,
    const int* thetaLaser,
    int laserIdx,
    Vec3<int>& posXyz,
    int& deltaZ);

  inline int32_t decodePointPositionZAngularExtension(
    const Vec3<int>& angularOrigin,
    const Vec3<int>& nodePos,
    const int* zLaser,
    const int* thetaLaser,
    int laserIdx,
    int maskz,
    Vec3<int>& posXyz);

  bool decodeNodeQpOffsetsPresent();
  int decodeQpOffset();

  bool decodeIsIdcm();

  template<class OutputIt>
  int decodeDirectPosition(
    bool geom_inter_idcm_enabled_flag,
    const DirectMode& predMode,
    bool geom_unique_points_flag,
    bool joint_2pt_idcm_enabled_flag,
    bool geom_angular_mode_enabled_flag,
    const Vec3<int>& nodeSizeLog2,
    int shiftBits,
    const Vec3<int>& posQuantBitMasks,
    const PCCOctree3Node& node,
    const OctreeNodePlanar& planar,
    const PCCPointSet3& pointCloudWorld,
    const Vec3<int>& headPos,
    const int* zLaser,
    const int* thetaLaser,
    int numLasers,
    OutputIt outputPoints);

  int decodeThetaRes(int prevThetaRes);
  int decodeZRes();

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

  // last previously coded laser index residual for a given node laser index
  std::vector<int> _prevLaserIndexResidual;
  std::vector<int> _prevLaserInterIndexResidual;

  // azimuthal elementary shifts
  AzimuthalPhiZi _phiZi;

  // Octree extensions
  bool _angularExtension;
  ///< IDCM Coding
  bool _angularDecideIDCMEligible;
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
  , _prevLaserIndexResidual(gps.numLasers(), 0x00000000)
  , _prevLaserInterIndexResidual(gps.numLasers(), 0x00000000)
  , _phiZi(gps.numLasers(), gps.angularNumPhiPerTurn)
  , _angularExtension(gps.octree_angular_extension_flag)
  , _angularDecideIDCMEligible(gps.one_point_alone_laser_beam_flag) 
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
  bool planarEligibilityDynamicOBUF,
  OctreeNodePlanar& planar,
  int planeZ,
  int dist,
  int adjPlanes,
  int planeId,
  int contextAngle,
  bool* multiPlanarFlag,
  bool* multiPlanarEligible,
  OctreeNodePlanar& planarRef,
  const OctreeNodePlanar adjNeighPlanar[7],
  const bool& neighAvai,
  const uint32_t& neighOccu)
{
  const int mask0 = (1 << planeId);
  const int mask1[3] = {6, 5, 3};

  // decode planar mode
  bool isPlanarRef = planarRef.planarMode & mask0;
  int planeBitRef = (planarRef.planePosBits & mask0) == 0 ? 0 : 1;

  int ctxIdx_Planar_flag = planeId;
  if (isPlanarRef)
    ctxIdx_Planar_flag += 3 * (planeBitRef + 1);

  bool isPlanar = isPlanarRef;
  

  if (!planar.isPCM) {
    if (_planar._geom_multiple_planar_mode_enable_flag) {
      bool multiPlanarFlagFalse = true;
      static const int planeId2Index[3][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}};
      for (int i = 0; i < 3; i++) {
        multiPlanarFlagFalse &= !(multiPlanarFlag[planeId2Index[planeId][i]]);
      }
      bool inferredPlanarFalse =
        multiPlanarFlagFalse;  // todo: consider renaming inferredPlaneFalse

      if (multiPlanarFlagFalse) {
        if (planeId == 2) {
          if (multiPlanarEligible[0])  //xyz
            inferredPlanarFalse =
              !((planar.planarMode & 2) && (planar.planarMode & 1));
          else if (multiPlanarEligible[2])  //xz
            inferredPlanarFalse = !(planar.planarMode & 1);
          else if (multiPlanarEligible[3])  //yz
            inferredPlanarFalse = !(planar.planarMode & 2);

        } else if (planeId == 1) {
          if (multiPlanarEligible[1])  //xy
            inferredPlanarFalse = !(planar.planarMode & 1);
        }
      }

      if (inferredPlanarFalse)
        isPlanar =
          _arithmeticDecoder->decode(_ctxPlanarMode[ctxIdx_Planar_flag]);
      else if (!multiPlanarFlagFalse)
        isPlanar = true;
      else
        isPlanar = false;
    } else {
      isPlanar =
        _arithmeticDecoder->decode(_ctxPlanarMode[ctxIdx_Planar_flag]);
    }
  }

  planar.planarMode |= isPlanar ? mask0 : 0;

  if (!isPlanar) {
    planar.planarPossible &= mask1[planeId];
    return -1;
  }

  // decode the plane index // encode the plane index
  int planeBit;

  if (planar.isPCM) {
    planeBit = planeBitRef;
    planar.planePosBits |= (planeBit << planeId);
    return planeBit;
  }
  // Not PCM and signal the plane position bit
  if (
    planeId == planar.lastDirIdx && planar.isPreDirMatch && planar.allowPCM
    && isPlanarRef) {
    planeBit = (planeBitRef == 1) ? 0 : 1;
    planar.planePosBits |= (planeBit << planeId);
    return planeBit;
  }

  if (contextAngle == -1) {  // angular mode off
    static const int kAdjPlaneCtx[4] = {0, 1, 2, 0};
    int planePosCtx = kAdjPlaneCtx[adjPlanes];
    int discreteDist = 0;
    int lastIndexPlane2d = 0;
    int refPlane = 0;
    if (isPlanarRef) {
      refPlane = 1 + planeBitRef;
    }
    if (planarEligibilityDynamicOBUF) {
      if (planeZ < 0) {
        discreteDist = 1;
        planeZ = 0;
      }
      else {
        discreteDist = dist > (8 >> OctreePlanarBuffer::shiftAb);
      }
      lastIndexPlane2d = planeZ + (discreteDist << 1);

      if (neighAvai) {
        ///<get the planar information of co-planar neighbors
        int coPlaneNeighPlaneBits =
          (!!(adjNeighPlanar[0].planePosBits & mask0) << 2)
          | (!!(adjNeighPlanar[1].planePosBits & mask0) << 1)
          | !!(adjNeighPlanar[2].planePosBits & mask0);
        int coPlaneNeighPlaneMode =
          (!!(adjNeighPlanar[0].planarMode & mask0) << 2)
          | (!!(adjNeighPlanar[1].planarMode & mask0) << 1)
          | !!(adjNeighPlanar[2].planarMode & mask0);
        int coPlaneNeighbor =
          (coPlaneNeighPlaneBits << 3) | coPlaneNeighPlaneMode;  ///<6 bits
        ///<get the planar information of co-edger neighbors
        int coEdgesNeighPlaneBits =
          (!!(adjNeighPlanar[3].planePosBits & mask0) << 2)
          | (!!(adjNeighPlanar[4].planePosBits & mask0) << 1)
          | !!(adjNeighPlanar[5].planePosBits & mask0);
        int coEdgesNeighPlaneMode =
          (!!(adjNeighPlanar[3].planarMode & mask0) << 2)
          | (!!(adjNeighPlanar[4].planarMode & mask0) << 1)
          | !!(adjNeighPlanar[5].planarMode & mask0);
        int coEdgeNeighbor =
          (coEdgesNeighPlaneBits << 3) | coEdgesNeighPlaneMode;  ///<6 bits
        int coVertexNeighbor = (!!(adjNeighPlanar[6].planePosBits & mask0) << 1)
          | !!(adjNeighPlanar[6].planarMode & mask0);
        int ctx1 = (lastIndexPlane2d << 6) | coPlaneNeighbor;  ///<8bits
        int ctx2 = (planePosCtx << 8) | (coEdgeNeighbor << 2)
          | coVertexNeighbor;  ///< 10 bits
       
        planeBit = _MapPlanarPosition[refPlane][planeId].decodeEvolve(_arithmeticDecoder, _ctxMapPlanarPosition[planeId], ctx2, ctx1, &_planarOBUFleafNumber, _planarBufferOBUFleaves);

      }
      else {
        int ctx1 = 0b1 << 7;
        ctx1 |= (lastIndexPlane2d << 5) | ((planePosCtx & 3) << 3)
          | ((neighOccu >> 9) & 7);  ///<1+2+2+3=8bits
        int ctx2 = 0b1 << 9;
        ctx2 |= neighOccu & ((1 << 9) - 1);  ///<1+9=10bits
        
        planeBit = _MapPlanarPosition[refPlane][planeId].decodeEvolve(_arithmeticDecoder, _ctxMapPlanarPosition[planeId], ctx2, ctx1, &_planarOBUFleafNumber, _planarBufferOBUFleaves);
        
      }
    }
    else {
      if (planeZ < 0) {
        int planePostCtxTmp = planePosCtx;
        if (isPlanarRef) {
          planePostCtxTmp += 3 * (planeBitRef + 1);
        }
        planeBit =
          _arithmeticDecoder->decode(_ctxPlanarPlaneLastIndexZ[planePostCtxTmp]);

      }
      else {
        discreteDist = dist > (8 >> OctreePlanarBuffer::shiftAb);
        lastIndexPlane2d = planeZ + (discreteDist << 1);

        planeBit = _arithmeticDecoder->decode(
          _ctxPlanarPlaneLastIndex[refPlane][planeId][planePosCtx]
          [lastIndexPlane2d]);
      }
    }
    
  } else {  // angular mode on
    int refPlane = isPlanarRef ? (1 + planeBitRef) : 0;

    if (planeId == 2) {  // angular
      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndexAngular[refPlane][contextAngle]);
    } else {  // azimuthal
      planeBit = _arithmeticDecoder->decode(
        _ctxPlanarPlaneLastIndexAngularPhi[refPlane][contextAngle]);
    }
  }

  planar.planePosBits |= (planeBit << planeId);
  return planeBit;
}

//============================================================================
void
GeometryOctreeDecoder::derivePlanarPCMContextBuffer(
  OctreeNodePlanar& planar,
  OctreeNodePlanar& planarRef,
  OctreePlanarBuffer& planeBuffer,
  int xx,
  int yy,
  int zz)
{
  int matchedDir = 0;

  planarRef.ctxBufPCM = 4
    * (int(planar.eligible[0]) + int(planar.eligible[1])
       + int(planar.eligible[2]) - 1);
  assert(planarRef.ctxBufPCM >= 0);

  for (int planeId = 0; planeId < 3; planeId++) {
    if (planar.eligible[planeId]) {
      const int mask0 = (1 << planeId);

      bool isPlanarRef = planarRef.planarMode & mask0;
      int planeBitRef = (planarRef.planePosBits & mask0) == 0 ? 0 : 1;

      // Get the buffer information
      OctreePlanarBuffer::Row* planeBufferDir = planeBuffer.getBuffer(planeId);

      int coord3 = (planeId == 2) ? zz : (planeId == 1 ? yy : xx);

      const int rowLen = OctreePlanarBuffer::rowSize;

      if (planeBufferDir) {
        coord3 &= OctreePlanarBuffer::maskC;

        const auto row = planeBufferDir[coord3];

        const int idxMinDist = rowLen - 1;
        const int closestPlanarFlag = row[idxMinDist].planeIdx;
        const bool closestPL = (closestPlanarFlag > -1) ? true : false;
        const int closestPlane = closestPL ? closestPlanarFlag : 0;

        matchedDir +=
          int(closestPL == isPlanarRef && closestPlane == planeBitRef);

      }
    }
  }
  planarRef.ctxBufPCM += matchedDir;
}


//============================================================================

void
GeometryOctreeDecoder::determinePlanarMode(
  bool planarEligibilityDynamicOBUF,
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
  int contextAngle,
  bool* multiPlanarFlag,
  bool* multiPlanarEligible,
  OctreeNodePlanar& planarRef,
  const OctreeNodePlanar adjNeighPlanar[7],
  const bool& neighAvai,
  const uint32_t& neighOccu)
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

  int planeBit = decodePlanarMode(planarEligibilityDynamicOBUF,
    planar, closestPlanarFlag, closestDist, adjPlanes, planeId, contextAngle,
    multiPlanarFlag, multiPlanarEligible, planarRef, adjNeighPlanar, neighAvai,
    neighOccu);
  bool isPlanar = (planar.planarMode & planeSelector);

  planarRate[planeId] =
    (255 * planarRate[planeId] + (isPlanar ? 256 * 8 : 0) + 128) >> 8;

  if (planeBuffer)
    row[rowLen - 1] = {unsigned(maxCoord), planeBit};

  bool isPlanarRef = (planarRef.planarMode & planeSelector);
  int planeBitRef = (planarRef.planePosBits & planeSelector) == 0 ? 0 : 1;

  if (!((isPlanar == isPlanarRef) && (planeBit == planeBitRef))) {
    planar.isPreDirMatch = false;
  }
}

//============================================================================

void
GeometryOctreeDecoder::determinePlanarMode(
  bool planarEligibilityDynamicOBUF,
  bool adjacent_child_contextualization_enabled_flag,
  const bool planarEligible[3],
  int posInParent,
  const GeometryNeighPattern& gnp,
  PCCOctree3Node& child,
  OctreeNodePlanar& planar,
  int contextAngle,
  int contextAnglePhiX,
  int contextAnglePhiY
  , OctreeNodePlanar& planarRef
)
{
  int xx = child.pos[0];
  int yy = child.pos[1];
  int zz = child.pos[2];

  auto& planeBuffer = _planar._planarBuffer;

  OctreeNodePlanar adjNeighPlanar[7];
  if (planarEligibilityDynamicOBUF && gnp.neighOccuValid) {
    for (int idx = 0; idx < 7; ++idx)
      if (gnp.adjNeighOcc[idx])
        setPlanesFromOccupancy(gnp.adjNeighOcc[idx], adjNeighPlanar[idx]);
  }

  uint8_t planarEligibleMask = 0;
  planarEligibleMask |= planarEligible[2] << 2;
  planarEligibleMask |= planarEligible[1] << 1;
  planarEligibleMask |= planarEligible[0] << 0;
  planarRef.planarMode &= planarEligibleMask;
  planarRef.planePosBits &= planarEligibleMask;

  if (planar.allowPCM) {
    derivePlanarPCMContextBuffer(planar, planarRef, planeBuffer, xx, yy, zz);
  }

  if (!planar.isRead && planar.allowPCM) {
    planar.isPCM = _arithmeticDecoder->decode(
      _ctxPlanarCopyMode[planarRef.ctxBufPCM][planarRef.planarMode]);
    planar.isRead = true;
  }
  bool multiPlanarFlag[4] = {false, false, false, false};
  bool multiPlanarEligible[4] = {false, false, false, false};
  if (_planar._geom_multiple_planar_mode_enable_flag) {
    if (!planar.isPCM) {
      if (planarEligible[2] && planarEligible[1] && planarEligible[0]) {
        multiPlanarEligible[0] = true;
        multiPlanarFlag[0] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
      } else if (
        (!planarEligible[2]) && planarEligible[1] && planarEligible[0]) {  //xy
        multiPlanarEligible[1] = true;
        multiPlanarFlag[1] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
      } else if (
        planarEligible[2] && (!planarEligible[1]) && planarEligible[0]) {  //xz
        multiPlanarEligible[2] = true;
        multiPlanarFlag[2] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
      } else if (
        planarEligible[2] && planarEligible[1] && (!planarEligible[0])) {  //yz
        multiPlanarEligible[3] = true;
        multiPlanarFlag[3] = _arithmeticDecoder->decode(_ctxMultiPlanarMode);
      }
    }
  }



  // planar x
  if (planarEligible[0]) {
    determinePlanarMode(planarEligibilityDynamicOBUF,
      adjacent_child_contextualization_enabled_flag, 0, planar,
      planeBuffer.getBuffer(0), yy, zz, xx, posInParent, gnp,
      child.siblingOccupancy, _planar._rate.data(), contextAnglePhiX,
      multiPlanarFlag, multiPlanarEligible, planarRef, adjNeighPlanar,
      gnp.neighOccuValid, gnp.neighborOccu);
  }
  // planar y
  if (planarEligible[1]) {
    determinePlanarMode(planarEligibilityDynamicOBUF,
      adjacent_child_contextualization_enabled_flag, 1, planar,
      planeBuffer.getBuffer(1), xx, zz, yy, posInParent, gnp,
      child.siblingOccupancy, _planar._rate.data(), contextAnglePhiY,
      multiPlanarFlag, multiPlanarEligible, planarRef, adjNeighPlanar,
      gnp.neighOccuValid, gnp.neighborOccu

    );
  }
  // planar z
  if (planarEligible[2]) {
    determinePlanarMode(planarEligibilityDynamicOBUF,
      adjacent_child_contextualization_enabled_flag, 2, planar,
      planeBuffer.getBuffer(2), xx, yy, zz, posInParent, gnp,
      child.siblingOccupancy, _planar._rate.data(), contextAngle,
      multiPlanarFlag, multiPlanarEligible, planarRef, adjNeighPlanar,
      gnp.neighOccuValid, gnp.neighborOccu);
  }
}

//-------------------------------------------------------------------------
// decode node occupancy bits
//

static void (*pointer2FunctionContext[8])(
  OctreeNeighours&, int, int&, int&, bool&) = {
  makeGeometryAdvancedNeighPattern0,
  makeGeometryAdvancedNeighPattern1,
  makeGeometryAdvancedNeighPattern2,
  makeGeometryAdvancedNeighPattern3,
  makeGeometryAdvancedNeighPattern4,
  makeGeometryAdvancedNeighPattern5,
  makeGeometryAdvancedNeighPattern6,
  makeGeometryAdvancedNeighPattern7};

static const int LUTinitCoded0[27][6] = {
  {0, 0, 0, 0, 0, 0}, {4, 0, 2, 2, 2, 2}, {0, 4, 2, 2, 2, 2},
  {2, 2, 4, 0, 2, 2}, {4, 2, 4, 2, 3, 3}, {2, 4, 4, 2, 3, 3},
  {2, 2, 0, 4, 2, 2}, {4, 2, 2, 4, 3, 3}, {2, 4, 2, 4, 3, 3},
  {2, 2, 2, 2, 4, 0}, {4, 2, 3, 3, 4, 2}, {2, 4, 3, 3, 4, 2},
  {3, 3, 4, 2, 4, 2}, {4, 3, 4, 3, 4, 3}, {3, 4, 4, 3, 4, 3},
  {3, 3, 2, 4, 4, 2}, {4, 3, 3, 4, 4, 3}, {3, 4, 3, 4, 4, 3},
  {2, 2, 2, 2, 0, 4}, {4, 2, 3, 3, 2, 4}, {2, 4, 3, 3, 2, 4},
  {3, 3, 4, 2, 2, 4}, {4, 3, 4, 3, 3, 4}, {3, 4, 4, 3, 3, 4},
  {3, 3, 2, 4, 2, 4}, {4, 3, 3, 4, 3, 4}, {3, 4, 3, 4, 3, 4}};

uint32_t
GeometryOctreeDecoder::decodeOccupancyFullNeihbourgsNZ(
  const GeometryNeighPattern& gnp,
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  bool planarPossibleX,
  bool planarPossibleY,
  bool planarPossibleZ,
  const MortonMap3D& occupancyAtlas,
  Vec3<int32_t> &pos,
  const int atlasShift,
  int predOcc,
  bool flagNoSingle,
  const bool& planarEligibleKOctreeDepth)
{
  bool sure_planarityX = planarMaskX || !planarPossibleX;
  bool sure_planarityY = planarMaskY || !planarPossibleY;
  bool sure_planarityZ = planarMaskZ || !planarPossibleZ;
  const int maxPerPlaneX = planarMaskX && flagNoSingle ? 2 : 3;
  const int maxPerPlaneY = planarMaskY && flagNoSingle ? 2 : 3;
  const int maxPerPlaneZ = planarMaskZ && flagNoSingle ? 2 : 3;
  const int maxAll = flagNoSingle ? 6 : 7;

  //int  MaskConfig = !planarMaskX ? 0 : planarMaskX == 15 ? 1 : 2;
  //MaskConfig += !planarMaskY ? 0 : planarMaskY == 51 ? 3 : 6;
  //MaskConfig += !planarMaskZ ? 0 : planarMaskZ == 85 ? 9 : 18;
  int MaskConfig = (!!planarMaskX) * (1 + (planarMaskX != 0x0F));
  MaskConfig += (!!planarMaskY) * 3 * (1 + (planarMaskY != 0x33));
  MaskConfig += (!!planarMaskZ) * 9 * (1 + (planarMaskZ != 0x55));

  int coded0[6] = {0, 0, 0, 0, 0, 0};  // mask x0 x1 y0 y1 z0 z1
  if (MaskConfig) {
    memcpy(coded0, LUTinitCoded0[MaskConfig], 6 * sizeof(int));
  }

  OctreeNeighours octreeNeighours;
  prepareGeometryAdvancedNeighPattern(
    octreeNeighours, gnp, pos, atlasShift, occupancyAtlas, planarEligibleKOctreeDepth);

  // loop on occupancy bits from occupancy map
  uint32_t partialOccupancy = 0;
  uint32_t occupancy = 0;
  int maskedOccupancy = planarMaskX | planarMaskY | planarMaskZ;
  for (int i = 0; i < 8; i++) {
    if ((maskedOccupancy >> i) & 1) {
      //  bit is 0 because masked by QTBT or planar
      partialOccupancy <<= 1;
      continue;
    }

    int mask0X = (0xf0 >> i) & 1;
    int mask0Y = 2 + ((0xcc >> i) & 1);
    int mask0Z = 4 + ((0xaa >> i) & 1);

    bool bitIsOne =
      (sure_planarityX && coded0[mask0X] >= maxPerPlaneX)
      || (coded0[0] + coded0[1] >= maxAll)
      || (sure_planarityY && coded0[mask0Y] >= maxPerPlaneY)
      || (coded0[2] + coded0[3] >= maxAll)
      || (sure_planarityZ && coded0[mask0Z] >= maxPerPlaneZ)
      || (coded0[4] + coded0[5] >= maxAll);

    if (bitIsOne) {
      occupancy += 1 << i;
      partialOccupancy <<= 1;
      partialOccupancy |= 1;
      continue;
    }

    int bitPred = (predOcc >> i) & 1;
    int interCtx = bitPred;

    // OBUF contexts
    int ctx1, ctx2;
    bool Sparse;
    (*pointer2FunctionContext[i])(
      octreeNeighours, occupancy,  ctx1, ctx2, Sparse);

    // decode

    int bit;
    if (Sparse) {
      bit = _MapOccupancySparse[interCtx][i].decodeEvolve(_arithmeticDecoder, _CtxMapDynamicOBUF, ctx2, ctx1, &_OBUFleafNumber, _BufferOBUFleaves);
    }
    else {
      bit = _MapOccupancy[interCtx][i].decodeEvolve(_arithmeticDecoder, _CtxMapDynamicOBUF, ctx2, ctx1, &_OBUFleafNumber, _BufferOBUFleaves);
    }

    // update partial occupancy of current node
    occupancy += bit << i;
    coded0[mask0X] += !bit;
    coded0[mask0Y] += !bit;
    coded0[mask0Z] += !bit;
    partialOccupancy <<= 1;
    partialOccupancy |= bit;
  }

  return occupancy;
}

//-------------------------------------------------------------------------
// decode node occupancy bits
//
uint32_t
GeometryOctreeDecoder::decodeOccupancyFullNeihbourgs(
  const GeometryNeighPattern& gnp,
  int planarMaskX,
  int planarMaskY,
  int planarMaskZ,
  bool planarPossibleX,
  bool planarPossibleY,
  bool planarPossibleZ,
  const MortonMap3D& occupancyAtlas,
  Vec3<int32_t> pos,
  const int atlasShift,
  bool flagWord4,
  bool adjacent_child_contextualization_enabled_flag,
  int predOcc,
  const bool& planarEligibleKOctreeDepth)
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
  //------ Z occupancy decoding from here ----------------

  bool flagNoSingle = false;

  if (
    gnp.neighPattern == 0
    && (!predOcc || (planarMaskX | planarMaskY | planarMaskZ))) {
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

    flagNoSingle = true;
    // at least two child nodes occupied and two planars => we know the occupancy
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

  return decodeOccupancyFullNeihbourgsNZ(
    gnp, planarMaskX, planarMaskY, planarMaskZ,
    planarPossibleX, planarPossibleY, planarPossibleZ,
    occupancyAtlas, pos, atlasShift, predOcc, flagNoSingle, planarEligibleKOctreeDepth);
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
  const bool& enableInter,
  const int& predLaserIdx,
  const OctreeAngPosScaler& quant,
  const Vec3<int>& nodeSizeLog2Rem,
  const Vec3<int>& angularOrigin,
  const int* zLaser,
  const int* thetaLaser,
  int nodeLaserIdx,
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
  int resLaser = decodeThetaRes(enableInter ? _prevLaserInterIndexResidual[nodeLaserIdx]
    : _prevLaserIndexResidual[nodeLaserIdx]);
  int laserIdx = predLaserIdx + resLaser;

  if (_angularExtension)
    if (enableInter)
      _prevLaserInterIndexResidual[nodeLaserIdx] = resLaser;
    else
      _prevLaserIndexResidual[nodeLaserIdx] = resLaser;

  // find predictor
  const int thInterp = 1 << 13;

  int phiNode = iatan2(posXyz[1], posXyz[0]);
  int phiTop = directAxis
    ? iatan2(posXyz[1], posXyz[0] + (1 << nodeSizeLog2Rem[!directAxis]))
    : iatan2(posXyz[1] + (1 << nodeSizeLog2Rem[!directAxis]), posXyz[0]);
  int phiMiddle = (phiNode + phiTop) >> 1;
  if (_angularExtension && !(std::abs(phiNode - phiTop) < thInterp))
    phiMiddle = directAxis
      ? iatan2(
          posXyz[1], posXyz[0] + ((1 << nodeSizeLog2Rem[!directAxis]) >> 1))
      : iatan2(
          posXyz[1] + ((1 << nodeSizeLog2Rem[!directAxis]) >> 1), posXyz[0]);

  int predPhi = _phiBuffer[laserIdx];
  int phiRef = _angularExtension ? phiMiddle : phiNode;
  if (predPhi == 0x80000000)
    predPhi = phiRef;

  // elementary shift predictor
  int nShift =
    ((predPhi - phiRef) * _phiZi.invDelta(laserIdx) + (1 << 29)) >> 30;
  predPhi -= _phiZi.delta(laserIdx) * nShift;

  // azimuthal code x or y
  const int phiAxis = !directAxis;
  for (int mask = (1 << nodeSizeLog2Rem[phiAxis]) >> 1,
           shiftBits = nodeSizeLog2Rem[phiAxis];
       mask; mask >>= 1, shiftBits--) {
    // angles left and right
    int scaledMask = quant.scaleEns(phiAxis, mask);

    int phiL, phiR;

    if (_angularExtension) {
      const int offset = scaledMask - 1;
      const int offset2 = shiftBits > 1 ? (shiftBits > 2 ? 0 : 1) : 2;

      phiL =
        phiNode + ((offset - offset2) * (phiMiddle - phiNode) >> (shiftBits));
      phiR = phiMiddle
        + ((offset + offset2) * (phiMiddle - phiNode) >> (shiftBits));
    } else {
      phiL = phiNode;
      phiR = directAxis ? iatan2(posXyz[1], posXyz[0] + scaledMask)
                        : iatan2(posXyz[1] + scaledMask, posXyz[0]);
    }

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
    int ctxIndex = 0;
    if (_angularExtension)
      ctxIndex = determineContextIndexForAngularPhiIDCM(
        _phiZi.delta(laserIdx), std::abs(phiL - phiR));
    auto& ctx =
      _ctxPlanarPlaneLastIndexAngularPhiIDCM[contextAnglePhi][ctxIndex];
    bool bit = _arithmeticDecoder->decode(ctx);
    delta[phiAxis] <<= 1;
    if (bit) {
      delta[phiAxis] |= 1;
      posXyz[phiAxis] += scaledMask;
      if (_angularExtension)
        phiNode = phiMiddle;
      else {
        phiNode = phiR;
        predPhi = _phiBuffer[laserIdx];
        if (predPhi == 0x80000000)
          predPhi = phiNode;

        // elementary shift predictor
        int nShift =
          ((predPhi - phiNode) * _phiZi.invDelta(laserIdx) + (1 << 29)) >> 30;
        predPhi -= _phiZi.delta(laserIdx) * nShift;
      }
    } else if (_angularExtension)
      phiTop = phiMiddle;

    if (_angularExtension) {
      // update Phi middle
      if (std::abs(phiNode - phiTop) < thInterp)
        phiMiddle = (phiNode + phiTop) >> 1;
      else
        phiMiddle = directAxis
          ? iatan2(posXyz[1], posXyz[0] + (scaledMask >> 1))
          : iatan2(posXyz[1] + (scaledMask >> 1), posXyz[0]);

      // update elementary shift predictor
      int nShift =
        ((predPhi - phiMiddle) * _phiZi.invDelta(laserIdx) + (1 << 29)) >> 30;
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
  if (_angularExtension)
    delta[2] = decodePointPositionZAngularExtension(
      angularOrigin, nodePos, zLaser, thetaLaser, laserIdx, maskz, posXyz);
  else
    delta[2] = decodePointPositionZAngular(
      quant, nodeSizeLog2Rem, zLaser, thetaLaser, laserIdx, posXyz, delta[2]);

  return delta;
}

int32_t
GeometryOctreeDecoder::decodePointPositionZAngular(
  const OctreeAngPosScaler& quant,
  const Vec3<int>& nodeSizeLog2Rem,
  const int* zLaser,
  const int* thetaLaser,
  int laserIdx,
  Vec3<int>& posXyz,
  int& deltaZ)
{
  uint64_t xLidar = (int64_t(posXyz[0]) << 8) - 128;
  uint64_t yLidar = (int64_t(posXyz[1]) << 8) - 128;
  uint64_t r2 = xLidar * xLidar + yLidar * yLidar;
  int64_t rInv = irsqrt(r2);

  // code bits for z using angular. Eligility is implicit. Laser is known.
  int64_t hr = zLaser[laserIdx] * rInv;
  int fixedThetaLaser =
    thetaLaser[laserIdx] + int(hr >= 0 ? -(hr >> 17) : ((-hr) >> 17));

  int maskz = (1 << nodeSizeLog2Rem[2]) >> 1;
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
    deltaZ <<= 1;
    deltaZ |= _arithmeticDecoder->decode(ctx);
    if (deltaZ & 1)
      deltaZ += scaledMaskZ;
  }

  return deltaZ;
}

int32_t
GeometryOctreeDecoder::decodePointPositionZAngularExtension(
  const Vec3<int>& angularOrigin,
  const Vec3<int>& nodePos,
  const int* zLaser,
  const int* thetaLaser,
  int laserIdx,
  int maskz,
  Vec3<int>& posXyz)
{
  uint64_t xLidar = (int64_t(posXyz[0]) << 8);
  uint64_t yLidar = (int64_t(posXyz[1]) << 8);
  uint64_t r2 = xLidar * xLidar + yLidar * yLidar;
  int64_t r = isqrt(r2);

  // decode z
  int64_t zRec26 = thetaLaser[laserIdx] * r;
  zRec26 -= int64_t(zLaser[laserIdx]) << 23;
  int32_t zRec = divExp2RoundHalfInf(zRec26, 26);

  zRec = std::max(zRec, posXyz[2]);
  zRec = std::min(zRec, posXyz[2] + (2 * maskz - 1));

  int32_t zRes = decodeZRes();

  return zRes + zRec + angularOrigin[2] - nodePos[2];
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
  bool geom_inter_idcm_enabled_flag,
  const DirectMode& predMode,
  bool geom_unique_points_flag,
  bool joint_2pt_idcm_enabled_flag,
  bool geom_angular_mode_enabled_flag,
  const Vec3<int>& nodeSizeLog2,
  int shiftBits,
  const Vec3<int>& posQuantBitMask,
  const PCCOctree3Node& node,
  const OctreeNodePlanar& planar,
  const PCCPointSet3& pointCloudWorld,
  const Vec3<int>& angularOrigin,
  const int* zLaser,
  const int* thetaLaser,
  int numLasers,
  OutputIt outputPoints)
{
  int numPoints = 1;
  int numPredFramePoints = node.predEnd - node.predStart;
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
  if (predMode == DirectMode::kAllPointSame)
    numPredFramePoints = 1;

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

  numPredFramePoints =
    numPredFramePoints < numPoints ? numPredFramePoints : numPoints;
  Vec3<int> predPoints[2];
  for (int i = 0; i < numPredFramePoints; i++)
    predPoints[i] = pointCloudWorld[node.predStart + i] >> shiftBits;

  // decode (ordred) two points
  Vec3<int32_t> deltaPos[2] = {deltaPlanar, deltaPlanar};
  if (numPoints == 2 && joint_2pt_idcm_enabled_flag)
    decodeOrdered2ptPrefix(directIdcm, nodeSizeLog2Rem, deltaPos);

  int laserIdx;
  int predLaserIdx;
  if (geom_angular_mode_enabled_flag) {
    auto delta = (deltaPos[0] << nodeSizeLog2Rem);
    delta += (1 << nodeSizeLog2Rem) >> 1;
    delta = quant.scaleEns(delta);
    if (_angularExtension)
      laserIdx =
        findLaserPrecise(posNodeLidar + delta, thetaLaser, zLaser, numLasers);
    else
      laserIdx = findLaser(posNodeLidar + delta, thetaLaser, numLasers);
    predLaserIdx = laserIdx;
  }

  Vec3<int32_t> pos;
  const bool& canInterPred = geom_inter_idcm_enabled_flag
    && predMode != DirectMode::kUnavailable && numPredFramePoints > 0;
  for (int i = 0; i < numPoints; i++) {
    int predIdx = numPredFramePoints == 2 ? i : 0;
    point_t predPoint = predPoints[predIdx];
    if (geom_angular_mode_enabled_flag) {
      if (canInterPred) {
        if (_angularExtension)
          predLaserIdx = findLaserPrecise(
            quant.scaleEns(predPoint) - angularOrigin, thetaLaser, zLaser,
            numLasers);
        else
          predLaserIdx = findLaser(
            quant.scaleEns(predPoint) - angularOrigin, thetaLaser, numLasers);
      }
      *(outputPoints++) = pos = decodePointPositionAngular(canInterPred, predLaserIdx,
        quant, nodeSizeLog2Rem, angularOrigin, zLaser, thetaLaser, laserIdx,
        node.pos << nodeSizeLog2, posNodeLidar, deltaPos[i]);
    }
      
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
GeometryOctreeDecoder::decodeThetaRes(int prevThetaRes)
{
  int ctx = prevThetaRes != 0;

  if (!_arithmeticDecoder->decode(_ctxThetaRes[ctx][0]))
    return 0;

  int absVal = 1;
  absVal += _arithmeticDecoder->decode(_ctxThetaRes[ctx][1]);
  if (absVal > 1)
    absVal += _arithmeticDecoder->decode(_ctxThetaRes[ctx][2]);
  if (absVal == 3)
    absVal += _arithmeticDecoder->decodeExpGolomb(1, _ctxThetaResExp);

  int ctxSign = (prevThetaRes > 0) + 2 * (prevThetaRes < 0);
  bool sign = _arithmeticDecoder->decode(_ctxThetaResSign[ctxSign]);
  return sign ? -absVal : absVal;
}

//-------------------------------------------------------------------------

int
GeometryOctreeDecoder::decodeZRes()
{
  if (!_arithmeticDecoder->decode(_ctxZRes[0]))
    return 0;

  int absVal = 1;
  absVal += _arithmeticDecoder->decode(_ctxZRes[1]);
  if (absVal > 1)
    absVal += _arithmeticDecoder->decode(_ctxZRes[2]);
  if (absVal == 3)
    absVal += _arithmeticDecoder->decodeExpGolomb(1, _ctxZResExp);

  bool sign = _arithmeticDecoder->decode(_ctxZResSign);
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
compensateGlobalMotion(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& predPointCloud,
  PCCPointSet3& pointPredictorWorld,
  EntropyDecoder* arithmeticDecoder,
  const Vec3<int> minimum_position,
  bool predDir = 0)
{
  switch (gbh.lpu_type) {
  case 0:  // object and road classification
    if (predDir)
      compensateWithRoadObjClassfication(
        pointPredictorWorld, gbh.gm_matrix2, gbh.gm_trans2, gbh.gm_thresh2,
        minimum_position);
    else
      compensateWithRoadObjClassfication(
        pointPredictorWorld, gbh.gm_matrix, gbh.gm_trans, gbh.gm_thresh,
        minimum_position);
    break;
  case 1:  // cuboid partition
    decodeCompensateWithCuboidPartition(
      predPointCloud, pointPredictorWorld, gbh, minimum_position,
      arithmeticDecoder, predDir);
    break;
  }
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
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  PCCPointSet3& predPointCloud2,
  const Vec3<int> minimum_position
)
{
  const bool isInter = gbh.interPredictionEnabledFlag;

  PCCPointSet3 predPointCloud;

  if (isInter)
    predPointCloud = refFrame->cloud;

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
  node00.predStart = uint32_t(0);
  node00.predEnd = uint32_t(predPointCloud.getPointCount());
  node00.numSiblingsMispredicted = 0;
  bool enableBiPred = gbh.biPredictionEnabledFlag;
  node00.predStart2 = uint32_t(0);
  node00.predEnd2 =
    enableBiPred ? uint32_t(predPointCloud2.getPointCount()) : 0;
  node00.predDir = 0;
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

  // global motion is here
  auto updatePredictorWorld = [&](
                                PCCPointSet3& predictorCloud,
                                PCCPointSet3& cloud, const bool dir,
                                uint32_t& predEnd) {
    predictorCloud = cloud;
    if (gps.globalMotionEnabled) {
      const auto minPos = (gbh.min_zero_origin_flag) ? 0 : minimum_position;
      compensateGlobalMotion(
        gps, gbh, cloud, predictorCloud, &arithmeticDecoder, minPos, dir);
    }
    predEnd = uint32_t(predictorCloud.getPointCount());
    for (int i = 0; i < predEnd; i++) {
      predictorCloud[i] -= gbh.geomBoxOrigin;
    }
  };
  PCCPointSet3 pointPredictorWorld;
  if (isInter)
    updatePredictorWorld(
      pointPredictorWorld, predPointCloud, false, node00.predEnd);

  PCCPointSet3 pointPredictorWorld2;
  if (isInter && enableBiPred) {
    updatePredictorWorld(
      pointPredictorWorld2, predPointCloud2, true, node00.predEnd2);

    if (gps.frameMergeEnabledFlag) {
      pointPredictorWorld.append(pointPredictorWorld2);
      pointPredictorWorld2.clear();
      enableBiPred = false;

      node00.predEnd = pointPredictorWorld.getPointCount();
      node00.predEnd2 = 0;
    }
  }

  // The number of nodes to wait before updating the planar rate.
  // This is to match the prior behaviour where planar is updated once
  // per coded occupancy.
  int nodesBeforePlanarUpdate = 1;

  const bool planarEligibilityDynamicOBUF =
    gps.geom_planar_mode_enabled_flag
    && gps.geom_octree_planar_dynamic_obuf_eligibiity_enabled_flag
    && !gps.geom_angular_mode_enabled_flag;

  decoder.resetMap(isInter, planarEligibilityDynamicOBUF);

  bool planarEligibleKOctreeDepth = 0;
  int numPointsCodedByIdcm = 0;
  const bool checkPlanarEligibilityBasedOnOctreeDepth =
    gps.geom_planar_mode_enabled_flag
    && gps.geom_octree_depth_planar_eligibiity_enabled_flag
    && !gps.geom_angular_mode_enabled_flag;

  for (int depth = 0; depth < maxDepth; depth++) {
    int numSubnodes = 0;
    // setup at the start of each level
    auto fifoCurrLvlEnd = fifo.end();
    int numNodesNextLvl = 0;
    Vec3<int32_t> occupancyAtlasOrigin = 0xffffffff;

    // derive per-level node size related parameters
    auto nodeSizeLog2 = lvlNodeSizeLog2[depth];
    auto childSizeLog2 = lvlNodeSizeLog2[depth + 1];
    // represents the largest dimension of the current node
    int nodeMaxDimLog2 = nodeSizeLog2.max();

    enableBiPred &= ((1 << nodeSizeLog2[0]) >= 2);

    auto pointSortMask = qtBtChildSize(nodeSizeLog2, childSizeLog2);

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

      // inter prediction from the first reference frame is enabled
      bool enabledPred = isInter && (enableBiPred || !node0.predDir);
      // inter prediction from the second reference frame is enabled
      bool enabledPred2 = isInter && (enableBiPred || node0.predDir);

      // sort the predictor into eight child partitions
      //  - perform an 8-way counting sort of the current node's points
      //  - (later) map to child nodes

      auto sortPredicate = [=](const PCCPointSet3::Proxy& proxy) {
        const auto& point = *proxy;
        return !!(int(point[2]) & pointSortMask[2])
          | (!!(int(point[1]) & pointSortMask[1]) << 1)
          | (!!(int(point[0]) & pointSortMask[0]) << 2);
      };

      // sort and partition the predictor
      std::array<int, 8> predCounts = {};
      std::array<int, 8> predCounts2 = {};
      if (isInter) {
        if (enabledPred)
          countingSort(
            PCCPointSet3::iterator(&pointPredictorWorld, node0.predStart),
            PCCPointSet3::iterator(&pointPredictorWorld, node0.predEnd),
            predCounts, sortPredicate);
        if (enabledPred2)
          countingSort(
            PCCPointSet3::iterator(&pointPredictorWorld2, node0.predStart2),
            PCCPointSet3::iterator(&pointPredictorWorld2, node0.predEnd2),
            predCounts2, sortPredicate);     
      }
      // generate the bitmap of child occupancy and count
      // the number of occupied children in node0.

     auto determineOccupancy = [&](std::array<int, 8>& predCnt) -> int {
        int occ = 0;
        for (int i = 0; i < 8; i++)
          if (predCnt[i])
            occ |= 1 << i;
        return occ;
      };
      const int predOccupancy1 =
        enabledPred ? determineOccupancy(predCounts) : 0;
      const int predOccupancy2 =
        enabledPred2 ? determineOccupancy(predCounts2) : 0;

      int predOccupancy = node0.predDir ? predOccupancy2 : predOccupancy1;

      bool occupancyIsPredictable =
        predOccupancy && node0.numSiblingsMispredicted <= 5;
      // The predictor may be cleared for the purpose of context
      // selection if the prediction is unlikely to be good.
      // NB: any other tests should use the original prediction.
      int predOccupancyReal = predOccupancy;
      if (!occupancyIsPredictable) {
        predOccupancy = 0;
      }

      int occupancyIsPredicted = 0;
      int occupancyPrediction = 0;

      if (nodeQpOffsetsPresent) {
        node0.qp = sliceQp;
        node0.qp += decoder.decodeQpOffset() << gps.geom_qp_multiplier_log2;
      }

      OctreeNodePlanar planarRef;
      if (isInter)
        setPlanesFromOccupancy(predOccupancy, planarRef);

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
          gps.adjacent_child_contextualization_enabled_flag,
          node0.pos, codedAxesPrevLvl, occupancyAtlas, planarEligibilityDynamicOBUF && planarEligibleKOctreeDepth);

      } else {
        gnp.neighPattern =
          neighPatternFromOccupancy(posInParent, node0.siblingOccupancy);
      }

      DirectMode predMode = DirectMode::kUnavailable;
      if (
        gps.geom_inter_idcm_enabled_flag && !isLeafNode(effectiveNodeSizeLog2))
        predMode = canInterFrameEncodeDirectPosition(
          decoder._angularDecideIDCMEligible, gps.geom_unique_points_flag,
          node0, pointPredictorWorld, nodeSizeLog2, angularOrigin, thetaLaser,
          numLasers, deltaAngle, decoder._phiZi, posQuantBitMasks);
     
      bool isDirectMode = false;
      // At the scaling depth, it is possible for a node that has previously
      // been marked as being eligible for idcm to be fully quantised due
      // to the choice of QP.  There is therefore nothing to code with idcm.
      if (isLeafNode(effectiveNodeSizeLog2))
        node0.idcmEligible = false;
      bool planar_eligibility_idcm_angular = true;

      if (node0.idcmEligible) {
        if (gps.geom_planar_disabled_idcm_angular_flag) {
          isDirectMode = decoder.decodeIsIdcm();
          if (isDirectMode && gps.geom_angular_mode_enabled_flag)
            planar_eligibility_idcm_angular = false;
        }
      }

      int contextAngle = -1;
      int contextAnglePhiX = -1;
      int contextAnglePhiY = -1;
      if (
        gps.geom_angular_mode_enabled_flag
        && planar_eligibility_idcm_angular) {
        contextAngle = determineContextAngleForPlanar(
          node0, nodeSizeLog2, angularOrigin, zLaser, thetaLaser, numLasers,
          deltaAngle, decoder._phiZi, decoder._phiBuffer.data(),
          &contextAnglePhiX, &contextAnglePhiY, posQuantBitMasks);

      }

      if (
        gps.geom_planar_mode_enabled_flag && planar_eligibility_idcm_angular
        && !gps.geom_octree_depth_planar_eligibiity_enabled_flag) {
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
        if (
          gps.geom_planar_mode_enabled_flag
          && planar_eligibility_idcm_angular) {
          if (gps.geom_octree_depth_planar_eligibiity_enabled_flag) {
            if (gps.geom_angular_mode_enabled_flag) {
              if (contextAngle != -1)
                planarEligible[2] = true;
              planarEligible[0] = (contextAnglePhiX != -1);
              planarEligible[1] = (contextAnglePhiY != -1);
            } else if (planarEligibleKOctreeDepth) {
              planarEligible[0] = true;
              planarEligible[1] = true;
              planarEligible[2] = true;
            }
          } else {
            decoder._planar.isEligible(planarEligible);
            if (gps.geom_angular_mode_enabled_flag) {
              if (contextAngle != -1)
                planarEligible[2] = true;
              planarEligible[0] = (contextAnglePhiX != -1);
              planarEligible[1] = (contextAnglePhiY != -1);
            }
          }

          for (int k = 0; k < 3; k++)
            planarEligible[k] &= (codedAxesCurNode >> (2 - k)) & 1;
        }

        planar.allowPCM = (isInter) && (occupancyIsPredictable)
          && (planarEligible[0] || planarEligible[1] || planarEligible[2]);
        planar.isPreDirMatch = true;
        planar.eligible[0] = planarEligible[0];
        planar.eligible[1] = planarEligible[1];
        planar.eligible[2] = planarEligible[2];
        planar.lastDirIdx =
          planarEligible[2] ? 2 : (planarEligible[1] ? 1 : 0);

        if (planarEligible[0] || planarEligible[1] || planarEligible[2])
          decoder.determinePlanarMode(planarEligibilityDynamicOBUF,
            gps.adjacent_child_contextualization_enabled_flag, planarEligible,
            posInParent, gnp, node0, planar, contextAngle, contextAnglePhiX,
            contextAnglePhiY, planarRef);
      }

      if (node0.idcmEligible) {
        if (!gps.geom_planar_disabled_idcm_angular_flag)
          isDirectMode = decoder.decodeIsIdcm();
        int idcmShiftBits = 0;
        if (isDirectMode) {
          auto idcmSize = effectiveNodeSizeLog2;
          if (idcmQp) {
            node0.qp = idcmQp;
            idcmShiftBits = QuantizerGeom::qpShift(idcmQp);
            idcmSize = nodeSizeLog2 - idcmShiftBits;
          }

          int numPoints = decoder.decodeDirectPosition(
            gps.geom_inter_idcm_enabled_flag, predMode,
            gps.geom_unique_points_flag, gps.joint_2pt_idcm_enabled_flag,
            gps.geom_angular_mode_enabled_flag, idcmSize, idcmShiftBits,
            posQuantBitMasks, node0, planar, pointPredictorWorld,
            angularOrigin, zLaser, thetaLaser, numLasers,
            &pointCloud[processedPointCount]);

          // calculating the number of points coded by IDCM for determining eligibility of planar mode
          if (checkPlanarEligibilityBasedOnOctreeDepth)
            numPointsCodedByIdcm += numPoints;

          for (int j = 0; j < numPoints; j++) {
            auto& point = pointCloud[processedPointCount++];
            for (int k = 0; k < 3; k++)
              point[k] += rotateLeft(node0.pos[k], idcmSize[k]);

            point = invQuantPosition(node0.qp, posQuantBitMasks, point);
          }

          // NB: no further siblings to decode by definition of IDCM
          if (!gps.one_point_alone_laser_beam_flag && gps.inferred_direct_coding_mode <= 1)
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

        bool flagWord4 =
          gps.neighbour_avail_boundary_log2_minus1 > 0;  //&& intraPredUsed;
        occupancy = decoder.decodeOccupancyFullNeihbourgs(
          gnp, planarMask[0], planarMask[1], planarMask[2],
          planar.planarPossible & 1, planar.planarPossible & 2,
          planar.planarPossible & 4, occupancyAtlas, node0.pos,
          codedAxesPrevLvl, flagWord4,
          gps.adjacent_child_contextualization_enabled_flag, predOccupancy, planarEligibilityDynamicOBUF
          && planarEligibleKOctreeDepth);
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

      // calculating the number of subnodes for determining eligibility of planar mode
      if (checkPlanarEligibilityBasedOnOctreeDepth)
        numSubnodes += numOccupied;

      int predFailureCount =
        enabledPred ? popcnt(uint8_t(occupancy ^ predOccupancy1)) : 0;
      int predFailureCount2 =
        enabledPred2 ? popcnt(uint8_t(occupancy ^ predOccupancy2)) : 0;
      int predPointsStartIdx2 = node0.predStart2;
      int predPointsStartIdx = node0.predStart;
      // nodeSizeLog2 > 1: for each child:
      //  - determine elegibility for IDCM
      //  - directly decode point positions if IDCM allowed and selected
      //  - otherwise, insert split children into fifo while updating neighbour state
      for (int i = 0; i < 8; i++) {
        uint32_t mask = 1 << i;
        if (!(occupancy & mask)) {
          predPointsStartIdx += predCounts[i];
          predPointsStartIdx2 += predCounts2[i];
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
          Vec3<int32_t> point{
            (node0.pos[0] << !!(codedAxesCurLvl & 4)) + x,
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
        child.predStart = predPointsStartIdx;
        child.predEnd = predPointsStartIdx + predCounts[i];
        predPointsStartIdx = child.predEnd;

        child.predStart2 = predPointsStartIdx2;
        predPointsStartIdx2 += predCounts2[i];
        child.predEnd2 = predPointsStartIdx2;
        child.predDir = node0.predDir;
        if (enableBiPred) {
          if (!predCounts2[i])
            child.predDir = 0;
          else if (!predCounts[i])
            child.predDir = 1;
          else
            child.predDir = predFailureCount != predFailureCount2
              ? (predFailureCount >= predFailureCount2)
              : node0.predDir;
        }
        predFailureCount =
          node0.predDir ? predFailureCount2 : predFailureCount;
        child.numSiblingsMispredicted = predFailureCount;

        if (isInter && !gps.geom_angular_mode_enabled_flag)
          child.idcmEligible = isDirectModeEligible_Inter(
            gps.inferred_direct_coding_mode, nodeMaxDimLog2, gnp.neighPattern,
            node0, child, occupancyIsPredictable);
        else
          child.idcmEligible = isDirectModeEligible(
            gps.inferred_direct_coding_mode, nodeMaxDimLog2, gnp.neighPattern,
            node0, child, occupancyIsPredictable,
            gps.geom_angular_mode_enabled_flag);

        if (child.idcmEligible) {
          child.idcmEligible &= idcmEnableMask & 1;
          idcmEnableMask = rotateRight(idcmEnableMask, 1);
        }

        numNodesNextLvl++;
      }
    }
    if (checkPlanarEligibilityBasedOnOctreeDepth)
      planarEligibleKOctreeDepth =
        (ringBufferSize - numPointsCodedByIdcm) * 10 < numSubnodes * 13;

    // Check that one level hasn't produced too many nodes
    // todo(df): this check is too weak to spot overflowing the fifo
    assert(numNodesNextLvl <= ringBufferSize);
  }

  decoder.clearMap(isInter,planarEligibilityDynamicOBUF);

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
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  PCCPointSet3& predPointCloud2,
  const Vec3<int> minimum_position
)
{
  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMem, arithmeticDecoder, nullptr, refFrame,
    predPointCloud2, minimum_position);
}

//-------------------------------------------------------------------------

void
decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  PCCPointSet3& predPointCloud2)
{
  pcc::ringbuf<PCCOctree3Node> nodes;
  decodeGeometryOctree(
    gps, gbh, minGeomNodeSizeLog2, pointCloud, ctxtMem, arithmeticDecoder,
    &nodes, refFrame, predPointCloud2, {0, 0, 0});

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
