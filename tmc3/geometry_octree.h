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
#include "geometry.h"
namespace pcc {

//============================================================================
struct AngularInformation {
  AngularInformation()
  {
    thetaEligible = false;
    phiEligible = false;
    valid = false;
    contextAngle = -1;
    contextAnglePhiX = -1;
    contextAnglePhiY = -1;
  }
  bool valid;
  bool thetaEligible;
  bool phiEligible;
  int laserIndex;
  int phiNode;
  int phiNode0;
  int theta32;
  uint64_t rInv;
  int contextAngle;
  int contextAnglePhiX;
  int contextAnglePhiY;
};
  enum class DirectMode
  {
    kUnavailable,
    kAllPointSame,
    kTwoPoints
  };
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

  // Range of second predictor's point indexes spanned by node
  uint32_t predStart2;
  uint32_t predEnd2;

  // Indicates the prediction direction
  bool predDir;

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
  const PCCOctree3Node& child,
  bool occupancyIsPredictable,
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
  bool occupancyIsPredictable)
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
    *ctxIdx += kCtxMapDynamicOBUFDelta[(255 - *ctxIdx) >> 4];
  else
    *ctxIdx -= kCtxMapDynamicOBUFDelta[*ctxIdx >> 4];

  return retval;
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

struct CtxModelDynamicOBUF {
  static const int kCtxFactorShift = 3;
  static const int kNumContexts = 256 >> kCtxFactorShift;
  AdaptiveBitModelFast contexts[kNumContexts];
  static const int kContextsInitProbability[kNumContexts];
  uint16_t obufSingleBound[33];

  CtxModelDynamicOBUF()
  {
    for (int i = 0; i < kNumContexts; i++)
      contexts[i].probability = kContextsInitProbability[i];
    for (int i = 0; i <= 32; i++) {
      obufSingleBound[i] = obufSingleBoundOrigin[i];
    }
  }

  AdaptiveBitModelFast& operator[](int idx)
  {
    return contexts[idx >> kCtxFactorShift];
  }
};

//============================================================================

class CtxMapDynamicOBUF {
public:
  static constexpr int kLeafDepth = 4;
  static constexpr int kLeafBufferSize = 20000;

  int S1 = 0; // 16;
  int S2 = 0; // 128 * 2 * 8;

  std::vector<uint8_t> CtxIdxMap; // S1*S2
  std::vector<uint8_t> kDown; // S2
  std::vector<uint8_t> Nseen; // S2

  ~CtxMapDynamicOBUF() { clear(); }

  //  allocate and reset CtxIdxMap to 127
  void reset(int userBitS1, int userBitS2);

  // initialize coder LUT
  void init(const uint8_t* initValue);

  //  deallocate CtxIdxMap
  void clear();

  //  decode bit  and update *ctxIdx according to bit
  int decodeEvolve(
    EntropyDecoder* _arithmeticDecoder,
    CtxModelDynamicOBUF& _ctxMapOccupancy,
    int i,
    int j,
    int* OBUFleafNumber,
    uint8_t* BufferOBUFleaves);

  //  get and update *ctxIdx according to bit
  uint8_t getEvolve(bool bit, int i, int j, int* OBUFleafNumber, uint8_t* BufferOBUFleaves);

private:
  int maxTreeDepth = 0;
  int minkTree = 0;

  //  update kDown
  void decreaseKdown(int idxTree, int kDownTree);
  void createLeaf(int idxTree, int kDownTree, int* OBUFleafNumber, uint8_t * BufferOBUFleaves, int ctx, int i);
  bool createLeafElement(int leafPos, uint8_t * BufferOBUFleaves, uint8_t ctx);
  uint8_t getEvolveLeaf(int leafPos, uint8_t * BufferOBUFleaves, bool bit, int i);
  int decodeEvolveLeaf(EntropyDecoder * _arithmeticDecoder, CtxModelDynamicOBUF & _ctxMapOccupancy, int leafPos, uint8_t * BufferOBUFleaves, int i);
  int idx(int i, int j);
};

//----------------------------------------------------------------------------

inline void
CtxMapDynamicOBUF::reset(int userBitS1, int userBitS2)
{
  S1 = 1 << userBitS1;
  S2 = 1 << userBitS2;

  maxTreeDepth = userBitS1 - kLeafDepth;
  minkTree = kLeafDepth;

  // tree of size (1 << maxTreeDepth) * S2
  const int treeSize = (1 << maxTreeDepth) * S2;
  kDown.resize(treeSize);
  Nseen.resize(treeSize);
  CtxIdxMap.resize(treeSize);

  std::fill_n(kDown.begin(), treeSize, userBitS1);
  std::fill_n(Nseen.begin(), S2, 0); // only needed for the S2 root nodes
  std::fill_n(CtxIdxMap.begin(), S2, 127); // only needed for the S2 root nodes
}

//----------------------------------------------------------------------------

inline void
CtxMapDynamicOBUF::init(const uint8_t* initValue) {
  for (int j = 0; j < S2; j++)
    CtxIdxMap[j] = initValue[j];
}

//----------------------------------------------------------------------------

inline void
CtxMapDynamicOBUF::clear()
{
  if (!S1 || !S2)
    return;

  kDown.resize(0);
  Nseen.resize(0);
  CtxIdxMap.resize(0);

  S1 = S2 = 0;
}

//----------------------------------------------------------------------------
inline bool
CtxMapDynamicOBUF::createLeafElement(int leafPos, uint8_t* BufferOBUFleaves, uint8_t ctx)
{
  int firstCtxIdx = leafPos * (1 << kLeafDepth);
  if (!BufferOBUFleaves[firstCtxIdx]) {
    std::fill_n(&BufferOBUFleaves[firstCtxIdx], (1 << kLeafDepth), ctx);
    return true;
  }
  return false;
}

//----------------------------------------------------------------------------
inline uint8_t
CtxMapDynamicOBUF::getEvolveLeaf(int leafPos, uint8_t* BufferOBUFleaves, bool bit, int i)
{
  int maskI = (1 << kLeafDepth) - 1;
  uint8_t* ctxIdx = &BufferOBUFleaves[leafPos * (1 << kLeafDepth) + (i & maskI)];
  uint8_t out = *ctxIdx;

  // coder index evolves
  if (bit)
    *ctxIdx += kCtxMapDynamicOBUFDelta[(255 - *ctxIdx) >> 4];
  else
    *ctxIdx -= kCtxMapDynamicOBUFDelta[*ctxIdx >> 4];

  return out;
}


//----------------------------------------------------------------------------
inline int
CtxMapDynamicOBUF::decodeEvolveLeaf(EntropyDecoder* _arithmeticDecoder, CtxModelDynamicOBUF& _ctxMapOccupancy, int leafPos, uint8_t* BufferOBUFleaves, int i) {
  int maskI = (1 << kLeafDepth) - 1;
  uint8_t* ctxIdx = &BufferOBUFleaves[leafPos * (1 << kLeafDepth) + (i & maskI)];
  int bit = _arithmeticDecoder->decode(
    (*ctxIdx) >> 3, _ctxMapOccupancy[*ctxIdx],
    _ctxMapOccupancy.obufSingleBound);

  // coder index evolves
  if (bit)
    *ctxIdx += kCtxMapDynamicOBUFDelta[(255 - *ctxIdx) >> 4];
  else
    *ctxIdx -= kCtxMapDynamicOBUFDelta[*ctxIdx >> 4];

  return bit;
}


//----------------------------------------------------------------------------
inline int
CtxMapDynamicOBUF::decodeEvolve(EntropyDecoder* _arithmeticDecoder, CtxModelDynamicOBUF& _ctxMapOccupancy, int i, int j, int* OBUFleafNumber, uint8_t* BufferOBUFleaves)
{
  int iTree = i >> kLeafDepth; // drop the bits that are in OBUF leaf
  int kDown0 = kDown[idx(iTree, j)];
  int bit;

  // ------------------ in Tree ---------------------
  if (kDown0 >= kLeafDepth) { // still in tree , not in OBUF leaf
    int kDownTree = kDown0 - kLeafDepth; // kdown in the tree part >=0
    int iP = (iTree >> kDownTree) << kDownTree; // erase bits
    int idxTree = idx(iP, j);  // index ofelements in the tree tables

    uint8_t* ctxIdx = &(CtxIdxMap[idxTree]); // get coder index
    bit = _arithmeticDecoder->decode(
      (*ctxIdx) >> 3, _ctxMapOccupancy[*ctxIdx],
      _ctxMapOccupancy.obufSingleBound);

    // coder index evolves
    if (bit)
      *ctxIdx += kCtxMapDynamicOBUFDelta[(255 - *ctxIdx) >> 4];
    else
      *ctxIdx -= kCtxMapDynamicOBUFDelta[*ctxIdx >> 4];

    // decrease number if erased bits if seens >= th
    int th = 3 + (std::abs(int(*ctxIdx) - 127) >> 4);
    if (++Nseen[idxTree] >= th) {
      if (kDownTree > 0) // we'll stay in tree
        decreaseKdown(idxTree, kDownTree); // kDownTree >0
      else  // we'll go to a leaf to be created int othe buffer
        createLeaf(idxTree, kDownTree, OBUFleafNumber, BufferOBUFleaves, *ctxIdx, i);

    }
  }
  // ------------------ in Leaf  ---------------------
  else { // in OBUF leaf
    int leafIdx = (CtxIdxMap[idx(iTree, j)] << 8) + Nseen[idx(iTree, j)]; // 16bit pointer hidden in CtxIdx and Nseen
    bit = decodeEvolveLeaf(_arithmeticDecoder, _ctxMapOccupancy, leafIdx, BufferOBUFleaves, i);
  }

  return bit;
}


//----------------------------------------------------------------------------
inline uint8_t
CtxMapDynamicOBUF::getEvolve(bool bit, int i, int j, int* OBUFleafNumber, uint8_t* BufferOBUFleaves)
{
  int iTree = i >> kLeafDepth; // drop the bits that are in OBUF leaf
  int kDown0 = kDown[idx(iTree, j)];
  uint8_t out;

  // ------------------ in Tree ---------------------
  if (kDown0 >= kLeafDepth) { // still in tree , not in OBUF leaf
    int kDownTree = kDown0 - kLeafDepth; // kdown in the tree part >=0
    int iP = (iTree >> kDownTree) << kDownTree; // erase bits
    int idxTree = idx(iP, j);  // index ofelements in the tree tables

    uint8_t* ctxIdx = &(CtxIdxMap[idxTree]); // get coder index
    out = *ctxIdx;

    // coder index evolves
    if (bit)
      *ctxIdx += kCtxMapDynamicOBUFDelta[(255 - *ctxIdx) >> 4];
    else
      *ctxIdx -= kCtxMapDynamicOBUFDelta[*ctxIdx >> 4];

    // decrease number if erased bits if seens >= th
    int th = 3 + (std::abs(int(*ctxIdx) - 127) >> 4);
    if (++Nseen[idxTree] >= th) {
      if (kDownTree > 0) // we'll stay in tree
        decreaseKdown(idxTree, kDownTree); // kDownTree >0
      else  // we'll go to a leaf to be created int othe buffer
        createLeaf(idxTree, kDownTree, OBUFleafNumber, BufferOBUFleaves, *ctxIdx, i);

    }

  }
  // ------------------ in Leaf  ---------------------
  else { // in OBUF leaf
    int leafIdx = (CtxIdxMap[idx(iTree, j)] << 8) + Nseen[idx(iTree, j)]; // 16bit pointer hidden in CtxIdx and Nseen
    out = getEvolveLeaf(leafIdx, BufferOBUFleaves, bit, i);
  }

  return out;
}


//----------------------------------------------------------------------------
inline void
CtxMapDynamicOBUF::decreaseKdown(int idxTree, int kDownTree)
{
  Nseen[idxTree] = 0;  // reintitlaize number of seen
  Nseen[idxTree + (S2 << kDownTree - 1)] = 0;
  int iEnd = S2 << kDownTree;
  for (int ii = 0; ii < iEnd; ii += S2)
    kDown[idxTree + ii]--; // decrease number of erased bits for all possible i involved (there are 2^kDownTree)

  auto* p = &CtxIdxMap[idxTree]; // coder index of first leaf in tree is here
  p[S2 << kDownTree - 1] = *p; // copy coder index to second leaf in tree
}


//----------------------------------------------------------------------------
inline void
CtxMapDynamicOBUF::createLeaf(int idxTree, int kDownTree, int* OBUFleafNumber, uint8_t* BufferOBUFleaves, int ctx, int i)
{
  bool  bufferAvailable = createLeafElement(*OBUFleafNumber, BufferOBUFleaves, ctx);
  if (bufferAvailable) {
    Nseen[idxTree] = (*OBUFleafNumber) & 255;// lower 8 bits
    CtxIdxMap[idxTree] = (*OBUFleafNumber) >> 8; // upper 8 bits
    *OBUFleafNumber += 1;
  }
  else {
    int dmin = 256;
    int bmin = *OBUFleafNumber;
    const int maskI = (1 << kLeafDepth) - 1;

    for (int b = *OBUFleafNumber; b < *OBUFleafNumber + 20 && b < kLeafBufferSize; b++) {
      int d = std::abs(ctx - BufferOBUFleaves[b * (1 << kLeafDepth) + (i & maskI)]);
      if (d < dmin) {
        dmin = d;
        bmin = b;
      }
    }
    Nseen[idxTree] = bmin & 255;// lower 8 bits
    CtxIdxMap[idxTree] = bmin >> 8; // upper 8 bits
    *OBUFleafNumber = bmin + 1;

  }

  if (*OBUFleafNumber >= kLeafBufferSize) // buffer not full
    *OBUFleafNumber = 0;
  kDown[idxTree]--; // same as  kDown[idx(iTree, j)]--;  kdown should be equal to kLeafDepth - 1 now
}


//----------------------------------------------------------------------------
inline int
CtxMapDynamicOBUF::idx(int i, int j)
{
  return i * S2 + j;
}


//============================================================================
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
  bool _geom_multiple_planar_mode_enable_flag;
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

void fracReduction(Rational& scale);

void compensateZCoordinate(
  PCCPointSet3& pointCloud,
  const GeometryParameterSet* gps,
  Rational scale,
  const Vec3<int>& angularOrigin,
  double& outputUnitLength,
  Vec3<int>& outputOrigin);

//---------------------------------------------------------------------------

inline int
determineContextIndexForAngularPhiIDCM(int deltaPhi, int phiLRDiff)
{
  return (3 * deltaPhi < phiLRDiff << 2) + (deltaPhi < phiLRDiff << 1);
}

//----------------------------------------------------------------------------

int findLaser(point_t point, const int* thetaList, const int numTheta);

//============================================================================

class GeometryOctreeContexts {
public:
  void reset();

  // dynamic OBUF
  void resetMap(const bool& enableInter, const bool& enablePlanar);
  void clearMap(const bool& enableInter, const bool& enablePlanar);

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
  AdaptiveBitModel _ctxMultiPlanarMode;
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
  CtxMapDynamicOBUF _MapOccupancy[4][8];
  CtxMapDynamicOBUF _MapOccupancySparse[4][8];
  CtxModelDynamicOBUF _CtxMapDynamicOBUF;

  uint8_t _BufferOBUFleaves[CtxMapDynamicOBUF::kLeafBufferSize * (1 << CtxMapDynamicOBUF::kLeafDepth)];
  int _OBUFleafNumber = 0;

  CtxMapDynamicOBUF _MapPlanarPosition[3][3];
  CtxModelDynamicOBUF _ctxMapPlanarPosition[3];
  uint8_t _planarBufferOBUFleaves[CtxMapDynamicOBUF::kLeafBufferSize * (1 << CtxMapDynamicOBUF::kLeafDepth)];
  int _planarOBUFleafNumber = 0;
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
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  const InterGeomEncOpts& interParams,
  const BiPredictionEncodeParams& biPredEncodeParams);

void decodeGeometryOctree(
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

);
//----------------------------------------------------------------------------------
void IsThetaPhiEligible(
  PCCOctree3Node& node,
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& angularOrigin,
  const int* thetaLaser,
  const int numLasers,
  int deltaAngle,
  const AzimuthalPhiZi& phiZi,
  Vec3<uint32_t> quantMasks,
  AngularInformation& angularInf);
//============================================================================================
inline DirectMode
canInterFrameEncodeDirectPosition(
  bool angularDecideIDCMEligible,
  bool geom_unique_points_flag,
  PCCOctree3Node& node,
  const PCCPointSet3& pointCloud,
  const Vec3<int>& nodeSizeLog2,
  const Vec3<int>& angularOrigin,
  const int* thetaLaser,
  const int numLasers,
  int deltaAngle,
  const AzimuthalPhiZi& phiZi,
  Vec3<uint32_t> posQuantBitMasks)
{
  AngularInformation angularInf;
  if (angularDecideIDCMEligible) {
    IsThetaPhiEligible(
      node, nodeSizeLog2, angularOrigin, thetaLaser, numLasers, deltaAngle,
      phiZi, posQuantBitMasks, angularInf);
    if (!geom_unique_points_flag)
      node.idcmEligible = angularInf.thetaEligible || angularInf.phiEligible;
    else
      node.idcmEligible = angularInf.thetaEligible && angularInf.phiEligible;
  }

  int numPoints = node.predEnd - node.predStart;
  // Check for duplicated points only if there are less than 10.
  // NB: this limit is rather arbitrary
  if (numPoints > 10)
    return DirectMode::kUnavailable;

  bool allPointsAreEqual = numPoints > 1 && !geom_unique_points_flag;
  for (auto idx = node.predStart + 1; allPointsAreEqual && idx < node.predEnd;
    idx++)
    allPointsAreEqual &= pointCloud[node.predStart] == pointCloud[idx];

  if (allPointsAreEqual)
    return DirectMode::kAllPointSame;

  if (numPoints > MAX_NUM_DM_LEAF_POINTS)
    return DirectMode::kUnavailable;

  return DirectMode::kTwoPoints;
}


}  // namespace pcc
