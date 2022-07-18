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

#pragma once

#include <cstdint>
#  include <cstring>
#include "DualLutCoder.h"
#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "geometry_params.h"
#include "hls.h"
#include "quantization.h"
#include "ringbuf.h"
#include "tables.h"
#include "TMC3.h"
namespace pcc {

//============================================================================

const int MAX_NUM_DM_LEAF_POINTS = 2;

//============================================================================

struct PCCOctree3Node {
  // 3D position of the current node's origin (local x,y,z = 0).
  Vec3<int32_t> pos;

  // Range of point indexes spanned by node
  uint32_t start;
  uint32_t end;

  // The current node's number of siblings plus one.
  // ie, the number of child nodes present in this node's parent.
  uint8_t numSiblingsPlus1;

  // Range of prediction's point indexes spanned by node
  uint32_t predStart;
  uint32_t predEnd;
  
  // The number of mispredictions in determining the occupancy
  // map of the child nodes in this node's parent.
  int8_t numSiblingsMispredicted;

  // The occupancy map used describing the current node and its siblings.
  uint8_t siblingOccupancy;

  // Indicatest hat the current node qualifies for IDCM
  bool idcmEligible;

  // The qp used for geometry quantisation.
  // NB: this qp value always uses a step size doubling interval of 8 qps
  int8_t qp;

  // angular
  uint8_t laserIndex = 255;
};

//============================================================================

struct OctreeNodePlanar {
  // planar; first bit for x, second bit for y, third bit for z
  uint8_t planarPossible = 7;
  uint8_t planePosBits = 0;
  uint8_t planarMode = 0;

  bool isPCM = false;
  bool isSignaled = false;
  bool isRead = false;

  bool allowPCM = true;
  bool isPreDirMatch = true;
  int lastDirIdx = 0;

  bool eligible[3] = {false, false, false};
  int ctxBufPCM = 0;

};

//---------------------------------------------------------------------------

int neighPatternFromOccupancy(int pos, int occupancy);

//---------------------------------------------------------------------------
// Determine if a node is a leaf node based on size.
// A node with all dimension = 0 is a leaf node.
// NB: some dimensions may be less than zero if coding of that dimension
// has already terminated.

inline bool
isLeafNode(const Vec3<int>& sizeLog2)
{
  return sizeLog2[0] <= 0 && sizeLog2[1] <= 0 && sizeLog2[2] <= 0;
}

//---------------------------------------------------------------------------
// Generates an idcm enable mask

uint32_t mkIdcmEnableMask(const GeometryParameterSet& gps);

//---------------------------------------------------------------------------
// Determine if direct coding is permitted.
// If tool is enabled:
//   - Block must not be near the bottom of the tree
//   - The parent / grandparent are sparsely occupied



inline bool
isDirectModeEligible(
  int intensity,
  int nodeSizeLog2,
  int nodeNeighPattern,
  const PCCOctree3Node& node,
  const PCCOctree3Node& child
  ,bool occupancyIsPredictable,
  bool isAngularModeEnabled

)
{
  if (!intensity)
    return false;
  if (occupancyIsPredictable && !isAngularModeEnabled)
    return false;

  if (intensity == 1)
    return (nodeSizeLog2 >= 2) && (nodeNeighPattern == 0)
      && (child.numSiblingsPlus1 == 1) && (node.numSiblingsPlus1 <= 2);

  if (intensity == 2)
    return (nodeSizeLog2 >= 2) && (nodeNeighPattern == 0);

  // This is basically unconditionally enabled.
  // If a node is that is IDCM-eligible is not coded with IDCM and has only
  // one child, then it is likely that the child would also not be able to
  // be coded with IDCM (eg, it still contains > 2 unique points).
  if (intensity == 3)
    return (nodeSizeLog2 >= 2) && (child.numSiblingsPlus1 > 1);

  return false;
}

inline bool
isDirectModeEligible_Inter(
  int intensity,
  int nodeSizeLog2,
  int nodeNeighPattern,
  const PCCOctree3Node& node,
  const PCCOctree3Node& child,
  bool occupancyIsPredictable
)
{

  if (!intensity)
    return false;

  if (occupancyIsPredictable)
    return false;

  return (nodeSizeLog2 >= 2) && (nodeNeighPattern == 0)
     && (child.numSiblingsPlus1 == 1) && (node.numSiblingsPlus1 <= 2);
}

//---------------------------------------------------------------------------
// Select the neighbour pattern reduction table according to GPS config.

inline const uint8_t*
neighPattern64toR1(const GeometryParameterSet& gps)
{
  if (gps.neighbour_avail_boundary_log2_minus1 > 0)
    return kNeighPattern64to9;
  return kNeighPattern64to6;
}

//---------------------------------------------------------------------------
// Encapsulates the derivation of ctxIdx for occupancy coding.

class CtxMapOctreeOccupancy {
public:
  struct CtxIdxMap {
    uint8_t b0[9];
    uint8_t b1[18];
    uint8_t b2[35];
    uint8_t b3[68];
    uint8_t b4[69];
    uint8_t b5[134];
    uint8_t b6[135];
    uint8_t b7[136];
  };

  CtxMapOctreeOccupancy();
  CtxMapOctreeOccupancy(const CtxMapOctreeOccupancy&);
  CtxMapOctreeOccupancy(CtxMapOctreeOccupancy&&);
  CtxMapOctreeOccupancy& operator=(const CtxMapOctreeOccupancy&);
  CtxMapOctreeOccupancy& operator=(CtxMapOctreeOccupancy&&);

  const uint8_t* operator[](int bit) const { return b[bit]; }

  uint8_t* operator[](int bit) { return b[bit]; }

  // return *ctxIdx and update *ctxIdx according to bit
  static uint8_t evolve(bool bit, uint8_t* ctxIdx);

private:
  std::unique_ptr<CtxIdxMap> map;
  std::array<uint8_t*, 8> b;
};

//----------------------------------------------------------------------------

inline uint8_t
CtxMapOctreeOccupancy::evolve(bool bit, uint8_t* ctxIdx)
{
  uint8_t retval = *ctxIdx;

  if (bit)
    *ctxIdx += kCtxMapOctreeOccupancyDelta[(255 - *ctxIdx) >> 4];
  else
    *ctxIdx -= kCtxMapOctreeOccupancyDelta[*ctxIdx >> 4];

  return retval;
}

//============================================================================
class CtxMapOccupancy {
public:
  int S1 = 0;  // 16;
  int S2 = 0;  // 128 * 2 * 8;

  uint8_t* CtxIdxMap;  // S1*S2
  uint8_t* kDown;      // S2
  uint8_t* Nseen;      // S2

  //  allocate and reset CtxIdxMap to 127
  void reset(int userBitS1, int userBitS2);

  //  deallocate CtxIdxMap
  void clear();

  //  return ctxIdx according to bit
  uint8_t get(int i, int j);

  //  update *ctxIdx according to bit
  void evolve(bool bit, int i, int j);
  //  update *ctxIdx according to bit
  uint8_t getEvolve(bool bit, int i, int j);

private:
  int maxTreeDepth = 3;
  int minkTree = 0;
  int LUTth[16];

  //  update kDown
  void decreaseKdown(int iP, int j, int iTree, int kDown0, int kTree);
  int idx(int i, int j);
};

inline void
CtxMapOccupancy::reset(int userBitS1, int userBitS2)
{
  S1 = 1 << userBitS1;
  S2 = 1 << userBitS2;
  maxTreeDepth = userBitS1 - 3;

  minkTree = userBitS1 - maxTreeDepth;

  for (int k = 0; k < 16; k++)
    LUTth[k] = std::min(64, 8 << std::max(0, minkTree - k - 1));

  kDown = new uint8_t[(1 << maxTreeDepth) * S2];
  memset(kDown, userBitS1, sizeof *kDown * (1 << maxTreeDepth) * S2);
  Nseen = new uint8_t[(1 << maxTreeDepth) * S2];
  memset(Nseen, 0, sizeof *Nseen * (1 << maxTreeDepth) * S2);
  CtxIdxMap = new uint8_t[S1 * S2];
  memset(CtxIdxMap, 127, sizeof *CtxIdxMap * S1 * S2);
}

inline void
CtxMapOccupancy::clear()
{
  if (!S1 || !S2)
    return;

  delete[] kDown;
  delete[] Nseen;
  delete[] CtxIdxMap;
}

inline uint8_t
CtxMapOccupancy::get(int i, int j)
{
  int kDown0 = kDown[idx(i >> minkTree, j)];
  int iP = (i >> kDown0) << kDown0;
  return CtxIdxMap[idx(iP, j)];
}

inline void
CtxMapOccupancy::evolve(bool bit, int i, int j)
{
  int iTree = i >> minkTree;
  int kDown0 = kDown[idx(iTree, j)];
  int iP = (i >> kDown0) << kDown0;
  uint8_t* ctxIdx = &(CtxIdxMap[idx(iP, j)]);

  if (bit)
    *ctxIdx += kCtxMapOctreeOccupancyDelta[(255 - *ctxIdx) >> 4];
  else
    *ctxIdx -= kCtxMapOctreeOccupancyDelta[*ctxIdx >> 4];

  if (kDown0) {
    int kTree = std::max(0, kDown0 - minkTree);
    iTree = (iTree >> kTree) << kTree;
    if (
      ++Nseen[idx(iTree, j)]
      >= LUTth[kDown0])  // if more than th stats per pack
      decreaseKdown(iP, j, iTree, kDown0, kTree);
  }
}

inline uint8_t
CtxMapOccupancy::getEvolve(bool bit, int i, int j)
{
  int iTree = i >> minkTree;
  int kDown0 = kDown[idx(iTree, j)];
  int iP = (i >> kDown0) << kDown0;

  uint8_t* ctxIdx = &(CtxIdxMap[idx(iP, j)]);
  uint8_t out = *ctxIdx;

  if (bit)
    *ctxIdx += kCtxMapOctreeOccupancyDelta[(255 - *ctxIdx) >> 4];
  else
    *ctxIdx -= kCtxMapOctreeOccupancyDelta[*ctxIdx >> 4];

  if (kDown0) {
    int kTree = std::max(0, kDown0 - minkTree);
    iTree = (iTree >> kTree) << kTree;
    if (
      ++Nseen[idx(iTree, j)]
      >= LUTth[kDown0])  // if more than th stats per pack
      decreaseKdown(iP, j, iTree, kDown0, kTree);
  }

  return out;
}

inline void
CtxMapOccupancy::decreaseKdown(int iP, int j, int iTree, int kDown0, int kTree)
{
  int idxTree = idx(iTree, j);
  Nseen[idxTree] = 0;  // setting other Nseen unneeded because initialized to 0
  if (kTree) {         // binary-tree based dynamic OBUF
    int iEnd = S2 << kTree;
    for (int ii = 0; ii < iEnd; ii += S2)
      kDown[idxTree + ii]--;

    auto* p = &CtxIdxMap[idx(iP, j)];
    p[S2 << kDown0 - 1] = *p;
  } else {  // simple dynamic OBUF on a subblock of states attached to a leaf of binary-tree
    kDown[idxTree]--;

    int S1Simple = S1 >> maxTreeDepth;
    int SS1 = (S1Simple >> kDown0);
    int stride = S2 << kDown0;

    auto* p = &CtxIdxMap[idx(iTree << minkTree, j)];
    for (int iL = 0; iL < SS1; iL++, p += stride)  // pack
      p[stride >> 1] = *p;
  }
}

inline int
CtxMapOccupancy::idx(int i, int j)
{
  return i * S2 + j;
}

//============================================================================
struct CtxModelOctreeOccupancy {
  static const int kCtxFactorShift = 4;
  AdaptiveBitModelFast contexts[256 >> kCtxFactorShift];

  AdaptiveBitModelFast& operator[](int idx)
  {
    return contexts[idx >> kCtxFactorShift];
  }
};

//---------------------------------------------------------------------------

struct CtxModelOctreeOccupancy2 {
  static const int kCtxFactorShift = 2;
  AdaptiveBitModelFast contexts[256 >> kCtxFactorShift];

  AdaptiveBitModelFast& operator[](int idx)
  {
    return contexts[idx >> kCtxFactorShift];
  }
};

//---------------------------------------------------------------------------
// generate an array of node sizes according to subsequent qtbt decisions

std::vector<Vec3<int>> mkQtBtNodeSizeList(
  const GeometryParameterSet& gps,
  const QtBtParameters& qtbt,
  const GeometryBrickHeader& gbh);

//---------------------------------------------------------------------------

inline Vec3<int>
qtBtChildSize(const Vec3<int>& nodeSizeLog2, const Vec3<int>& childSizeLog2)
{
  Vec3<int> bitpos = 0;
  for (int k = 0; k < 3; k++) {
    if (childSizeLog2[k] != nodeSizeLog2[k])
      bitpos[k] = 1 << childSizeLog2[k];
  }
  return bitpos;
}

//---------------------------------------------------------------------------

inline int
nonSplitQtBtAxes(const Vec3<int>& nodeSizeLog2, const Vec3<int>& childSizeLog2)
{
  int indicator = 0;
  for (int k = 0; k < 3; k++) {
    indicator <<= 1;
    indicator |= nodeSizeLog2[k] == childSizeLog2[k];
  }
  return indicator;
}

//============================================================================
// Scales quantized positions used internally in angular coding.
//
// NB: this is not used to scale output positions since generated positions
//     are not clipped to node boundaries.
//
// NB: there are two different position representations used in the codec:
//        ppppppssssss = original position
//        ppppppqqqq00 = pos, (quantisation) node size aligned -> use scaleNs()
//        00ppppppqqqq = pos, effective node size aligned -> use scaleEns()
//     where p are unquantised bits, q are quantised bits, and 0 are zero bits.

class OctreeAngPosScaler {
  QuantizerGeom _quant;
  Vec3<uint32_t> _mask;
  int _qp;

public:
  OctreeAngPosScaler(int qp, const Vec3<uint32_t>& quantMaskBits)
    : _quant(qp), _qp(qp), _mask(quantMaskBits)
  {}

  // Scale an effectiveNodeSize aligned position as the k-th position component.
  int scaleEns(int k, int pos) const;

  // Scale an effectiveNodeSize aligned position.
  Vec3<int> scaleEns(Vec3<int> pos) const;

  // Scale a NodeSize aligned position.
  Vec3<int> scaleNs(Vec3<int> pos) const;
};

//----------------------------------------------------------------------------

inline int
OctreeAngPosScaler::scaleEns(int k, int pos) const
{
  if (!_qp)
    return pos;

  int shiftBits = QuantizerGeom::qpShift(_qp);
  int lowPart = pos & (_mask[k] >> shiftBits);
  int highPart = pos ^ lowPart;
  int lowPartScaled = _quant.scale(lowPart);

  return (highPart << shiftBits) + lowPartScaled;
}

//----------------------------------------------------------------------------

inline Vec3<int32_t>
OctreeAngPosScaler::scaleEns(Vec3<int32_t> pos) const
{
  if (!_qp)
    return pos;

  for (int k = 0; k < 3; k++)
    pos[k] = scaleEns(k, pos[k]);

  return pos;
}
//----------------------------------------------------------------------------

inline Vec3<int32_t>
OctreeAngPosScaler::scaleNs(Vec3<int32_t> pos) const
{
  if (!_qp)
    return pos;

  // convert pos to effectiveNodeSize form
  return scaleEns(pos >> QuantizerGeom::qpShift(_qp));
}

//============================================================================

struct OctreePlanarBuffer {
  static constexpr unsigned numBitsC = 14;
  static constexpr unsigned numBitsAb = 5;
  static constexpr unsigned rowSize = 1;
  static_assert(numBitsC >= 0 && numBitsC <= 32, "0 <= numBitsC <= 32");
  static_assert(numBitsAb >= 0 && numBitsAb <= 32, "0 <= numBitsAb <= 32");
  static_assert(rowSize > 0, "rowSize must be greater than 0");
  static constexpr unsigned shiftAb = 3;
  static constexpr int maskAb = ((1 << numBitsAb) - 1) << shiftAb;
  static constexpr int maskC = (1 << numBitsC) - 1;

#pragma pack(push)
#pragma pack(1)
  struct Elmt {
    // maximum of two position components
    unsigned int pos : numBitsAb;

    // -2: not used, -1: not planar, 0: plane 0, 1: plane 1
    int planeIdx : 2;
  };
#pragma pack(pop)

  typedef Elmt Row[rowSize];

  OctreePlanarBuffer();
  OctreePlanarBuffer(const OctreePlanarBuffer& rhs);
  OctreePlanarBuffer(OctreePlanarBuffer&& rhs);
  ~OctreePlanarBuffer();

  OctreePlanarBuffer& operator=(const OctreePlanarBuffer& rhs);
  OctreePlanarBuffer& operator=(OctreePlanarBuffer&& rhs);

  void resize(Vec3<int> numBufferRows);
  void clear();

  // Access to a particular buffer column (dimension)
  Row* getBuffer(int dim) { return _col[dim]; }

private:
  // Backing storage for the underlying buffer
  std::vector<Elmt> _buf;

  // Base pointers for the first, second and third position components.
  std::array<Row*, 3> _col = {{nullptr, nullptr, nullptr}};
};

//============================================================================

struct OctreePlanarState {
  OctreePlanarState(const GeometryParameterSet&);

  OctreePlanarState(const OctreePlanarState&);
  OctreePlanarState(OctreePlanarState&&);
  OctreePlanarState& operator=(const OctreePlanarState&);
  OctreePlanarState& operator=(OctreePlanarState&&);

  bool _planarBufferEnabled;
  OctreePlanarBuffer _planarBuffer;

  std::array<int, 3> _rate{{128 * 8, 128 * 8, 128 * 8}};
  int _localDensity = 1024 * 4;

  std::array<int, 3> _rateThreshold;

  void initPlanes(const Vec3<int>& planarDepth);
  void updateRate(int occupancy, int numSiblings);
  void isEligible(bool eligible[3]);
};

// determine if a 222 block is planar
void setPlanesFromOccupancy(int occupancy, OctreeNodePlanar& planar);

int maskPlanarX(const OctreeNodePlanar& planar);
int maskPlanarY(const OctreeNodePlanar& planar);
int maskPlanarZ(const OctreeNodePlanar& planar);

void maskPlanar(OctreeNodePlanar& planar, int mask[3], int codedAxes);

int determineContextAngleForPlanar(
  PCCOctree3Node& node,
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& angularOrigin,
  const int* zLaser,
  const int* thetaLaser,
  const int numLasers,
  int deltaAngle,
  const AzimuthalPhiZi& phiZi,
  int* phiBuffer,
  int* contextAnglePhiX,
  int* contextAnglePhiY,
  Vec3<uint32_t> quantMasks);

//---------------------------------------------------------------------------

inline int determineContextIndexForAngularPhiIDCM(int deltaPhi, int phiLRDiff)
{
  return (3*deltaPhi < phiLRDiff<<2) + (deltaPhi < phiLRDiff<<1);
}

//----------------------------------------------------------------------------

int findLaser(point_t point, const int* thetaList, const int numTheta);

//============================================================================

class GeometryOctreeContexts {
public:
  void reset();

protected:
  AdaptiveBitModel _ctxSingleChild;
  AdaptiveBitModel _ctxZ[8][7][4];

  AdaptiveBitModel _ctxDupPointCntGt0;
  AdaptiveBitModel _ctxDupPointCntGt1;
  AdaptiveBitModel _ctxDupPointCntEgl;

  AdaptiveBitModel _ctxBlockSkipTh;
  AdaptiveBitModel _ctxNumIdcmPointsGt1;
  AdaptiveBitModel _ctxSameZ;

  // IDCM unordered
  AdaptiveBitModel _ctxSameBitHighx[5];
  AdaptiveBitModel _ctxSameBitHighy[5];
  AdaptiveBitModel _ctxSameBitHighz[5];

  // residual laser index
  AdaptiveBitModel _ctxThetaRes[2][3];
  AdaptiveBitModel _ctxThetaResSign[3];
  AdaptiveBitModel _ctxThetaResExp;

  // residual z
  AdaptiveBitModel _ctxZRes[3];
  AdaptiveBitModel _ctxZResSign;
  AdaptiveBitModel _ctxZResExp;

  AdaptiveBitModel _ctxQpOffsetAbsGt0;
  AdaptiveBitModel _ctxQpOffsetSign;
  AdaptiveBitModel _ctxQpOffsetAbsEgl;

  // for planar mode xyz

  AdaptiveBitModel _ctxPlanarCopyMode[16][8];
  AdaptiveBitModel _ctxPlanarMode[9];
  AdaptiveBitModel _ctxPlanarPlaneLastIndex[3][3][3][4];
  AdaptiveBitModel _ctxPlanarPlaneLastIndexZ[9];

  AdaptiveBitModel _ctxPlanarPlaneLastIndexAngular[3][4];
  AdaptiveBitModel _ctxPlanarPlaneLastIndexAngularIdcm[4];

  AdaptiveBitModel _ctxPlanarPlaneLastIndexAngularPhi[3][8];
 
  AdaptiveBitModel _ctxPlanarPlaneLastIndexAngularPhiIDCM[8][3];

  // For bitwise occupancy coding
  CtxModelOctreeOccupancy _ctxOccupancy;
  CtxMapOctreeOccupancy _ctxIdxMaps[24];


  // For bytewise occupancy coding
  DualLutCoder<true> _bytewiseOccupancyCoder[10];
  // OBUF somplified
  CtxMapOccupancy _MapOccupancy[4][8];
  CtxMapOccupancy _MapOccupancySparse[4][8];
  CtxModelOctreeOccupancy2 _ctxMapOccupancy;
};

//----------------------------------------------------------------------------

inline void
GeometryOctreeContexts::reset()
{
  this->~GeometryOctreeContexts();
  new (this) GeometryOctreeContexts;
}

//============================================================================
// :: octree encoder exposing internal ringbuffer

void encodeGeometryOctree(
  const OctreeEncOpts& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining,
  PCCPointSet3& predPointCloud,
  const Vec3<int> minimum_position);

void decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder& arithmeticDecoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining
  , PCCPointSet3 &predPointCloud
  , const Vec3<int> minimum_position

);

void applyGlobalMotion_with_shift(
  PCCPointSet3& PC,
  const std::vector<int> gm_matrix,
  const Vec3<int> gm_trans,
  const std::pair<int, int> thresh,
  const Vec3<int> minimum_position);
//============================================================================

}  // namespace pcc
