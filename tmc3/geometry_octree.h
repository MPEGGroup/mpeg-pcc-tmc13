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

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "hls.h"
#include "ringbuf.h"
#include "tables.h"

namespace pcc {

//============================================================================

const int MAX_NUM_DM_LEAF_POINTS = 2;

//============================================================================

struct PCCOctree3Node {
  // 3D position of the current node's origin (local x,y,z = 0).
  Vec3<int32_t> pos;
  // 3D position of the current node's origin (local x,y,z = 0) for decoder
  // reconstruction with in-tree geometry quantization.
  Vec3<int32_t> posQ;

  // Range of point indexes spanned by node
  uint32_t start;
  uint32_t end;

  // address of the current node in 3D morton order.
  int64_t mortonIdx;

  // pattern denoting occupied neighbour nodes.
  //    32 8 (y)
  //     |/
  //  2--n--1 (x)
  //    /|
  //   4 16 (z)
  uint8_t neighPattern = 0;

  // The current node's number of siblings plus one.
  // ie, the number of child nodes present in this node's parent.
  uint8_t numSiblingsPlus1;

  // The occupancy map used describing the current node and its siblings.
  uint8_t siblingOccupancy;

  // The qp used for geometry quantisation
  int qp;

  // planar; first bit for x, second bit for y, third bit for z
  uint8_t planarPossible = 7;
  uint8_t planePosBits = 0;
  uint8_t planarMode = 0;
};

//---------------------------------------------------------------------------
uint8_t mapGeometryOccupancy(uint8_t occupancy, uint8_t neighPattern);
uint8_t mapGeometryOccupancyInv(uint8_t occupancy, uint8_t neighPattern);

void updateGeometryNeighState(
  bool siblingRestriction,
  const ringbuf<PCCOctree3Node>::iterator& bufEnd,
  int64_t numNodesNextLvl,
  PCCOctree3Node& child,
  int childIdx,
  uint8_t neighPattern,
  uint8_t parentOccupancy);

//---------------------------------------------------------------------------
// :: octree encoder exposing internal ringbuffer

void encodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  EntropyEncoder* arithmeticEncoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining);

void decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minNodeSizeLog2,
  PCCPointSet3& pointCloud,
  EntropyDecoder* arithmeticDecoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining);

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
// Determine if direct coding is permitted.
// If tool is enabled:
//   - Block must not be near the bottom of the tree
//   - The parent / grandparent are sparsely occupied

inline bool
isDirectModeEligible(
  bool featureEnabled,
  int nodeSizeLog2,
  const PCCOctree3Node& node,
  const PCCOctree3Node& child)
{
  return featureEnabled && (nodeSizeLog2 >= 2) && (node.neighPattern == 0)
    && (child.numSiblingsPlus1 == 1) && (node.numSiblingsPlus1 <= 2);
}

//---------------------------------------------------------------------------
// Select the neighbour pattern reduction table according to GPS config.

using NeighPattern64toX = const uint8_t[64];

inline NeighPattern64toX&
neighPattern64toR1(const GeometryParameterSet& gps)
{
  if (!gps.neighbour_context_restriction_flag)
    return kNeighPattern64to9;
  return kNeighPattern64to6;
}

//---------------------------------------------------------------------------

struct CtxModelOctreeOccupancy {
  AdaptiveBitModelFast contexts[256];
  int ctxFactorShift;

  CtxModelOctreeOccupancy(int ctxFactorShift) : ctxFactorShift(ctxFactorShift)
  {}

  AdaptiveBitModelFast& operator[](int idx)
  {
    return contexts[idx >> ctxFactorShift];
  }
};

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

  const uint8_t* operator[](int bit) const { return b[bit]; }

  uint8_t* operator[](int bit) { return b[bit]; }

  // return *ctxIdx and update *ctxIdx according to bit
  static uint8_t evolve(bool bit, uint8_t* ctxIdx);

private:
  std::unique_ptr<CtxIdxMap> map;
  uint8_t* b[8];
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

//---------------------------------------------------------------------------

inline Vec3<int>
implicitQtBtDecision(
  Vec3<int> nodeSizeLog2,
  int maxNumImplicitQtbtBeforeOt,
  int minDepthImplicitQtbt)
{
  int nodeMinDimLog2 =
    std::min({nodeSizeLog2[0], nodeSizeLog2[1], nodeSizeLog2[2]});

  if (maxNumImplicitQtbtBeforeOt || nodeMinDimLog2 == minDepthImplicitQtbt) {
    // implicit qt bt
    int nodeMaxDimLog2 =
      std::max({nodeSizeLog2[0], nodeSizeLog2[1], nodeSizeLog2[2]});
    for (int k = 0; k < 3; k++) {
      if (nodeSizeLog2[k] == nodeMaxDimLog2)
        nodeSizeLog2[k]--;
    }
  } else  // octree partition
    nodeSizeLog2 = nodeSizeLog2 - 1;

  return nodeSizeLog2;
}

//---------------------------------------------------------------------------

inline void
updateImplicitQtBtParameters(
  const Vec3<int>& nodeSizeLog2,
  int trisoup_node_size_log2,
  int* maxNumImplicitQtbtBeforeOt,
  int* minSizeImplicitQtbt)
{
  int nodeMinDimLog2 =
    std::min({nodeSizeLog2[0], nodeSizeLog2[1], nodeSizeLog2[2]});
  int nodeMaxDimLog2 =
    std::max({nodeSizeLog2[0], nodeSizeLog2[1], nodeSizeLog2[2]});

  // max number of implicit qtbt before ot is bounded by difference between
  // max and min node size
  if (*maxNumImplicitQtbtBeforeOt > (nodeMaxDimLog2 - nodeMinDimLog2))
    *maxNumImplicitQtbtBeforeOt = nodeMaxDimLog2 - nodeMinDimLog2;
  // min depth of implicit qtbt is bounded by min node size
  if (*minSizeImplicitQtbt > nodeMinDimLog2)
    *minSizeImplicitQtbt = nodeMinDimLog2;
  // if all dimensions have same size, min depth of implicit qtbt should be 0
  if (nodeMaxDimLog2 == nodeMinDimLog2) {
    *minSizeImplicitQtbt = 0;
  }

  // if trisoup is enabled, perform qtbt first before ot
  if (trisoup_node_size_log2 != 0) {
    *maxNumImplicitQtbtBeforeOt = nodeMaxDimLog2 - nodeMinDimLog2;
    *minSizeImplicitQtbt = 0;
  }
}

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

//---------------------------------------------------------------------------

// determine if a 222 block is planar
void isPlanarNode(
  PCCPointSet3& pointCloud,
  PCCOctree3Node& node0,
  const Vec3<int>& sizeLog2,
  uint8_t& planarMode,
  uint8_t& planePosBits,
  bool planarEligible[3]);

void planarInitPlanes(
  const int kNumPlanarPlanes, int depth, int* planes3x3, int* planes[9]);
void updateplanarRate(
  int planarRate[3], int occupancy, int& localDensity, int numSiblings);
void eligilityPlanar(
  bool planarEligible[3],
  int plane_rate[3],
  const int threshold[3],
  int localDensity);

int maskPlanarX(PCCOctree3Node& node0, bool activatable);
int maskPlanarY(PCCOctree3Node& node0, bool activatable);
int maskPlanarZ(PCCOctree3Node& node0, bool activatable);

void maskPlanar(PCCOctree3Node& node0, int mask[3], const int occupancySkip);

//---------------------------------------------------------------------------

}  // namespace pcc
