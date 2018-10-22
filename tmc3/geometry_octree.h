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

#include "ArithmeticCodec.h"
#include "PCCMath.h"
#include "PCCPointSet.h"
#include "hls.h"
#include "ringbuf.h"

namespace pcc {

//============================================================================

const int MAX_NUM_DM_LEAF_POINTS = 2;

//============================================================================

struct PCCOctree3Node {
  // 3D position of the current node's origin (local x,y,z = 0).
  PCCVector3<uint32_t> pos;

  // Range of point indexes spanned by node
  uint32_t start;
  uint32_t end;

  // address of the current node in 3D morton order.
  uint64_t mortonIdx;

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
};

//---------------------------------------------------------------------------

uint8_t mapGeometryOccupancy(uint8_t occupancy, uint8_t neighPattern);
uint8_t mapGeometryOccupancyInv(uint8_t occupancy, uint8_t neighPattern);

void updateGeometryNeighState(
  bool siblingRestriction,
  const ringbuf<PCCOctree3Node>::iterator& bufEnd,
  int64_t numNodesNextLvl,
  int childSizeLog2,
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
  o3dgc::Arithmetic_Codec* arithmeticEncoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining);

void decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  o3dgc::Arithmetic_Codec* arithmeticDecoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining);

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

struct CtxModelOctreeOccupancy {
  o3dgc::Adaptive_Bit_Model_Fast b0[10];
  o3dgc::Adaptive_Bit_Model_Fast b1[2 * 10];
  o3dgc::Adaptive_Bit_Model_Fast b2[4 * 10];
  o3dgc::Adaptive_Bit_Model_Fast b3[8 * 10];
  o3dgc::Adaptive_Bit_Model_Fast b4[75];
  o3dgc::Adaptive_Bit_Model_Fast b5[112];
  o3dgc::Adaptive_Bit_Model_Fast b6[117];
  o3dgc::Adaptive_Bit_Model_Fast b7[127];

  o3dgc::Adaptive_Bit_Model* operator[](int bit_pos)
  {
    assert(unsigned(bit_pos) < 8);
    switch (bit_pos) {
    case 0: return b0;
    case 1: return b1;
    case 2: return b2;
    case 3: return b3;
    case 4: return b4;
    case 5: return b5;
    case 6: return b6;
    case 7: return b7;
    default: return nullptr;
    }
  }
};

//---------------------------------------------------------------------------

}  // namespace pcc
