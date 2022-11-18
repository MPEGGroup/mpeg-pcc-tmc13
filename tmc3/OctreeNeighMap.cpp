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

#include "OctreeNeighMap.h"

#include <iostream>

#if WITH_MEMCHECK
#  include <valgrind/memcheck.h>
#else
#  define VALGRIND_MAKE_MEM_UNDEFINED(a, b) (void)0
#endif

namespace pcc {

//============================================================================

void
MortonMap3D::clearUpdates()
{
  for (const auto byteIndex : _updates) {
    _buffer[byteIndex] = uint8_t(0);
  }
  _updates.resize(0);

  // Access to the child array is supposed to be guarded by checking the main
  // map first.  It is therefore not necessary to clear the array between
  // updates.  Setting it to zero just hides the issue -- ie, tools like
  // valgrind don't complain, but the logic error is still present.
  //
  // The following undoes the effect of any writes that have cleared the
  // undefined state.
  VALGRIND_MAKE_MEM_UNDEFINED(_childOccupancy.get(), _bufferSizeInBytes << 3);
}

//----------------------------------------------------------------------------

void
MortonMap3D::clear()
{
  memset(_buffer.get(), 0, _bufferSizeInBytes);
  _updates.resize(0);

  // See clearUpdates()
  VALGRIND_MAKE_MEM_UNDEFINED(_childOccupancy.get(), _bufferSizeInBytes << 3);
}

//============================================================================

void
updateGeometryOccupancyAtlas(
  const Vec3<int32_t>& currentPosition,
  const int atlasShift,
  const pcc::ringbuf<PCCOctree3Node>& fifo,
  const pcc::ringbuf<PCCOctree3Node>::iterator& fifoCurrLvlEnd,
  MortonMap3D* occupancyAtlas,
  Vec3<int32_t>* atlasOrigin)
{
  const uint32_t mask = (1 << occupancyAtlas->cubeSizeLog2()) - 1;
  const int shift = occupancyAtlas->cubeSizeLog2();
  const int shiftX = (atlasShift & 4 ? 1 : 0);
  const int shiftY = (atlasShift & 2 ? 1 : 0);
  const int shiftZ = (atlasShift & 1 ? 1 : 0);

  const auto currentOrigin = currentPosition >> shift;

  // only refresh the atlas if the current position lies outside the
  // the current atlas.
  if (*atlasOrigin == currentOrigin) {
    return;
  }

  *atlasOrigin = currentOrigin;
  occupancyAtlas->clearUpdates();

  for (auto it = fifo.begin(); it != fifoCurrLvlEnd; ++it) {
    if (currentOrigin != it->pos >> shift)
      break;
    const uint32_t x = (it->pos[0] & mask) >> shiftX;
    const uint32_t y = (it->pos[1] & mask) >> shiftY;
    const uint32_t z = (it->pos[2] & mask) >> shiftZ;
    occupancyAtlas->setByte(x, y, z, it->siblingOccupancy);
  }
}

//----------------------------------------------------------------------------

void
updateGeometryOccupancyAtlasOccChild(
  const Vec3<int32_t>& pos,
  uint8_t childOccupancy,
  MortonMap3D* occupancyAtlas)
{
  uint32_t mask = (1 << occupancyAtlas->cubeSizeLog2()) - 1;
  uint32_t x = pos[0] & mask;
  uint32_t y = pos[1] & mask;
  uint32_t z = pos[2] & mask;

  occupancyAtlas->setChildOcc(x, y, z, childOccupancy);
}

//----------------------------------------------------------------------------
// neighIdx: 0 => (x-1), 1 => (y-1), 2 => (z-1)
//
static GeometryNeighPattern
updatePatternFromNeighOccupancy(
  const MortonMap3D& occupancyAtlas,
  int x,
  int y,
  int z,
  GeometryNeighPattern gnp,
  int neighIdx,
  bool codedAxisCurLvl)
{
  static const uint8_t childMasks[] = {
    0xf0 /* x-1 */, 0xcc /* y-1 */, 0xaa /* z-1 */
  };

  uint32_t patternBit = 1 << (1 << neighIdx);
  uint8_t childMask = childMasks[neighIdx];

  // conversions between neighbour occupancy and adjacency:
  //  x: >> 4, y: >> 2, z: >> 1
  // Always inspect the adjacent children, taking into account that their
  // position changes depending upon whther the current axis is coded or not.
  if (!codedAxisCurLvl) {
    childMask ^= 0xff;
  }
  uint8_t child_occ = occupancyAtlas.getChildOcc(x, y, z);
  gnp.adjNeighOcc[neighIdx] = child_occ;

  return gnp;
}

//----------------------------------------------------------------------------

GeometryNeighPattern
makeGeometryNeighPattern(
  bool adjacent_child_contextualization_enabled_flag,
  const Vec3<int32_t>& position,
  int codedAxesPrevLvl,
  const MortonMap3D& occupancyAtlas)
{
  const int mask = occupancyAtlas.cubeSize() - 1;
  const int cubeSizeMinusOne = mask;
  const int32_t x = position[0] & mask;
  const int32_t y = position[1] & mask;
  const int32_t z = position[2] & mask;
  uint8_t neighPattern;

  const int sx = codedAxesPrevLvl & 4 ? 1 : 0;
  const int sy = codedAxesPrevLvl & 2 ? 1 : 0;
  const int sz = codedAxesPrevLvl & 1 ? 1 : 0;

  if (
    x > 0 && x < cubeSizeMinusOne && y > 0 && y < cubeSizeMinusOne && z > 0
    && z < cubeSizeMinusOne) {
    neighPattern = occupancyAtlas.get(x + 1, y, z, sx, sy, sz);
    neighPattern |= occupancyAtlas.get(x - 1, y, z, sx, sy, sz) << 1;
    neighPattern |= occupancyAtlas.get(x, y - 1, z, sx, sy, sz) << 2;
    neighPattern |= occupancyAtlas.get(x, y + 1, z, sx, sy, sz) << 3;
    neighPattern |= occupancyAtlas.get(x, y, z - 1, sx, sy, sz) << 4;
    neighPattern |= occupancyAtlas.get(x, y, z + 1, sx, sy, sz) << 5;
  } else {
    neighPattern = occupancyAtlas.getWithCheck(x + 1, y, z, sx, sy, sz);
    neighPattern |= occupancyAtlas.getWithCheck(x - 1, y, z, sx, sy, sz) << 1;
    neighPattern |= occupancyAtlas.getWithCheck(x, y - 1, z, sx, sy, sz) << 2;
    neighPattern |= occupancyAtlas.getWithCheck(x, y + 1, z, sx, sy, sz) << 3;
    neighPattern |= occupancyAtlas.getWithCheck(x, y, z - 1, sx, sy, sz) << 4;
    neighPattern |= occupancyAtlas.getWithCheck(x, y, z + 1, sx, sy, sz) << 5;
  }

  // Above, the neighbour pattern corresponds directly to the six same
  // sized neighbours of the given node.
  // The patten is then refined by examining the available children
  // of the same neighbours.

  // NB: the process of updating neighpattern below also derives
  // the occupancy contextualisation bits.
  GeometryNeighPattern gnp = {neighPattern, 0, 0, 0, {0, 0, 0}};

  if (!gnp.neighPattern || !adjacent_child_contextualization_enabled_flag)
    return gnp;

  if (x > 0 && gnp.neighPattern & 2)
    gnp.adjNeighOcc[0] = occupancyAtlas.getChildOcc(x - 1, y, z);

  if (y > 0 && gnp.neighPattern & 4)
    gnp.adjNeighOcc[1] = occupancyAtlas.getChildOcc(x, y - 1, z);

  if (z > 0 && gnp.neighPattern & 16)
    gnp.adjNeighOcc[2] = occupancyAtlas.getChildOcc(x, y, z - 1);

  return gnp;
}


//----------------------------------------------------------------------------

void
makeGeometryAdvancedNeighPattern(
  int neighPattern,
  const Vec3<int32_t>& position,
  int atlasShift,
  const MortonMap3D& occupancyAtlas,
  int Word7Adj[8],
  bool Sparse[8])
{
  const int mask = occupancyAtlas.cubeSize() - 1;
  const int32_t x = position[0] & mask;
  const int32_t y = position[1] & mask;
  const int32_t z = position[2] & mask;

  // ----- neighbours  FLB -----
  int occLeft = 0;
  if (neighPattern & 2)  //Neighbour L
    occLeft = occupancyAtlas.getChildOcc(x - 1, y, z);

  int occFront = 0;
  if (neighPattern & 4) //Neighbour F
    occFront = occupancyAtlas.getChildOcc(x, y - 1, z);

  int occBottom = 0;
  if (neighPattern & 16) //Neighbour B
    occBottom = occupancyAtlas.getChildOcc(x, y, z - 1);

  int occL = occLeft >> 4;
  int occF = ((occFront >> 2) & 3) | ((occFront >> 4) & 12);
  int occB = ((occBottom >> 1) & 1) | ((occBottom >> 2) & 2) | ((occBottom >> 3) & 4) | ((occBottom >> 4) & 8);
  int occOrLFBfb = occLeft | occFront | occBottom;

  // ----- neighbours  LB, FB, LF -----
  const int sx = atlasShift & 4 ? 1 : 0;
  const int sy = atlasShift & 2 ? 1 : 0;
  const int sz = atlasShift & 1 ? 1 : 0;

  int edgeBits = 0;
  if (occupancyAtlas.getWithCheck(x - 1, y, z - 1, sx, sy, sz)) {
    int occLB = occupancyAtlas.getChildOcc(x - 1, y, z - 1);
    edgeBits = ((occLB & 32) >> 5) | ((occLB & 128) >> 6);
  }

  if (occupancyAtlas.getWithCheck(x, y - 1, z - 1, sx, sy, sz)) {
    int occFB = occupancyAtlas.getChildOcc(x, y - 1, z - 1);
    edgeBits |= ((occFB & 8) >> 1) | ((occFB & 128) >> 4);
  }

  if (occupancyAtlas.getWithCheck(x - 1, y - 1, z, sx, sy, sz)) {
    int occLF = occupancyAtlas.getChildOcc(x - 1, y - 1, z);
    edgeBits |= (occLF & 0b11000000) >> 2;
  }

  int N3 = ((neighPattern >> 3) & 4) | ((neighPattern >> 2) & 2) | (neighPattern & 1);// gnp.neighPattern & 0b101001
  int N2 = N3 & 3;

  // ------ bit 0  ------ 12+0 = 12 bits
  int NN = !!(occL & 1) + !!(occL & 2) + !!(occL & 4) + !!(occL & 8);
  NN += !!(occF & 1) + !!(occF & 2) + !!(occF & 4) + !!(occF & 8);
  NN += !!(occB & 1) + !!(occB & 2) + !!(occB & 4) + !!(occB & 8);

  int NLFB = !!occL + !!occF + !!occB;
  if (NN>1) {
    if (NLFB == 3) {
      Word7Adj[0] = 0b1000 << 8; // put tag at the head as it is the most important info
      Word7Adj[0] |= (occL & 7);
      Word7Adj[0] |= (occF & 7) << 3;
      Word7Adj[0] |= (occB & 7) << 6;
    }
    else if (NLFB == 2) {
      if (occL && occB) {
        Word7Adj[0] = 0b0100 << 8;
        Word7Adj[0] |= occL;
        Word7Adj[0] |= occB << 4;
      }
      if (occF && occB) {
        Word7Adj[0] = 0b0101 << 8;
        Word7Adj[0] |= occF;
        Word7Adj[0] |= occB << 4;
      }
      if (occL && occF) {
        Word7Adj[0] = 0b0110 << 8;
        Word7Adj[0] |= occL;
        Word7Adj[0] |= occF << 4;
      }
    }
    else //NLFB == 1
    {
      if (occL) {
        Word7Adj[0] = 0b0000 << 8;
        Word7Adj[0] |= occL << 4;
        Word7Adj[0] |= (edgeBits & 0b001100);
      }
      if (occF) {
        Word7Adj[0] = 0b0001 << 8;
        Word7Adj[0] |= occF << 4;
        Word7Adj[0] |= (edgeBits & 0b000011) << 2;

      }
      if (occB) {
        Word7Adj[0] = 0b0010 << 8;
        Word7Adj[0] |= occB << 4;
        Word7Adj[0] |= (edgeBits & 0b110000) >> 2;
      }
    }
  }
  else {
    int neighPatternLFB = ((neighPattern & 0b110) >> 1) | ((neighPattern & 16) >> 2);
    Word7Adj[0] = neighPatternLFB << 4; // 9 bits
    if (NN) {
      if (occL)
        Word7Adj[0] |= 1 << 7;
      if (occF)
        Word7Adj[0] |= 2 << 7;
      if (occB)
        Word7Adj[0] |= 3 << 7;
    }
    if (neighPatternLFB) {
      if (occOrLFBfb & 1) {
        Word7Adj[0] |= 8;
        Word7Adj[0] |= (occLeft & 1);
        Word7Adj[0] |= (occFront & 1) << 1;
        Word7Adj[0] |= (occBottom & 1) << 2;
      }
      else
      {
        Word7Adj[0] |= ((occLeft & 2) || (occFront & 16) || (occBottom & 16));
        Word7Adj[0] |= ((occLeft & 4) || (occFront & 2) || (occBottom & 4)) << 1;
        Word7Adj[0] |= !edgeBits << 2;
      }
    }
    else {
      Word7Adj[0] |= !(edgeBits & 0b000011) << 1;
      Word7Adj[0] |= !(edgeBits & 0b001100) << 2;
      Word7Adj[0] |= !(edgeBits & 0b110000) << 3;
    }

    Sparse[0] = true;
  }


  // ------ bit 1  ------ 8+3+1 = 12 bits
  if (occL || occF) {
    Word7Adj[1] = occL << 3;
    Word7Adj[1] |= occF << 7;
    Word7Adj[1] |= N3;
  }
  else { // 5+1+1 = 7 bits
    Word7Adj[1] = (N3 >> 2) << 5; //Top
    if (occOrLFBfb & 2) {
      Word7Adj[1] |= 16;
      Word7Adj[1] |= (occLeft & 2) << 0;
      Word7Adj[1] |= (occFront & 2) << 1;
      Word7Adj[1] |= (occBottom & 2) << 2;
    }
    else {
      Word7Adj[1] |= ((occLeft & 1) || (occFront & 1)) << 1;
      Word7Adj[1] |= ((occLeft & 8) || (occFront & 32)) << 2;
      Word7Adj[1] |= !(edgeBits & 0b110101) << 3;
    }

    Word7Adj[1] |= !occB;

    Sparse[1] = true;
  }

  // ------ bit 2  ------ 8+2+2 = 12 bits
  if (occL || occB) {
    Word7Adj[2] = occL << 2;
    Word7Adj[2] |= occB << 6;
    Word7Adj[2] |= N2;
  }
  else { // 5+1+2 = 8 bits
    Word7Adj[2] = !occF;
    if (occOrLFBfb & 4) {
      Word7Adj[2] |= 16;
      Word7Adj[2] |= (occLeft & 4) << 1;
      Word7Adj[2] |= (occBottom & 4);
      Word7Adj[2] |= (occFront & 4) >> 1;
    }
    else {
      Word7Adj[2] |= ((occLeft & 1) || (occBottom & 1)) << 3;
      Word7Adj[2] |= ((occLeft & 8) || (occBottom & 64)) << 2;
      Word7Adj[2] |= !(edgeBits & 0b000011) << 1;
    }

    Sparse[2] = true;
    Word7Adj[2] |= ((N2 >> 1) & 1) << 5; // Back
  }

  // ------ bit 3  ------ 4+3+3 = 10 bits
  if (occL) {
    Word7Adj[3] = occL << 3;
    Word7Adj[3] |= N3;
  }
  else { // 6+2+3 = 11 bits
    Word7Adj[3] = !occF;
    Word7Adj[3] |= !occB << 1;
    if (occOrLFBfb & 8) {
      Word7Adj[3] |= 32;
      Word7Adj[3] |= (occLeft & 8) >> 1;
      Word7Adj[3] |= (occFront & 8);
      Word7Adj[3] |= (occBottom & 8) << 1;
    }
    else {
      Word7Adj[3] |= (occLeft & 0b110) << 2;
      Word7Adj[3] |= !(edgeBits & 0b110010) << 2;
    }

    Sparse[3] = true;
    Word7Adj[3] |= (N3 >> 1) << 6; // Back + Top
  }

  // ------ bit 4  ------ 6+2+4 = 12 bits
  if (occF || occB) {
    Word7Adj[4] = ((occF & 1) | ((occF >> 1) & 6)) << 2; // &0b1101;
    Word7Adj[4] |= ((occB & 1) | ((occB >> 1) & 6)) << 5; // &0b1101;
    Word7Adj[4] |= N2;
  }
  else { // 5+1+4 = 10 bits
    Word7Adj[4] = !occL;
    if (occOrLFBfb & 16) {
      Word7Adj[4] |= 16;
      Word7Adj[4] |= (occFront & 16) >> 1;
      Word7Adj[4] |= (occBottom & 16) >> 2;
      Word7Adj[4] |= (occLeft & 16) >> 3;
    }
    else {
      Word7Adj[4] |= ((occFront & 32) || (occBottom & 1)) << 3;
      Word7Adj[4] |= ((occFront & 64) || (occBottom & 64)) << 2;
      Word7Adj[4] |= !(edgeBits & 0b001100) << 1;
    }

    Sparse[4] = true;
    Word7Adj[4] |= (N2 & 1) << 5; // Right
  }

  // ------ bit 5  ------ 3+3+5 = 11 bits
  if (occF) {
    Word7Adj[5] = (occF & 0b1110) << 2;
    Word7Adj[5] |= N3;
  }
  else { // 4+2+5 = 11 bits
    Word7Adj[5] = !occL;
    Word7Adj[5] |= !occB << 1;
    Word7Adj[5] |= !(edgeBits & 0b111100) << 2;
    Word7Adj[5] |= !(occOrLFBfb & 32) << 3;
    Sparse[5] = true;
    Word7Adj[5] |= (N3 & 1) << 4; // Right
    Word7Adj[5] |= (N3 >> 2) << 5; //  and Top
  }

  // ------ bit 6  ------   3+2+6 = 11 bits
  if (occB) {
    Word7Adj[6] = (occB & 0b1110) << 1;
    Word7Adj[6] |= N2;
  }
  else { // 4+2+6 = 12 bits
    Word7Adj[6] = !occL;
    Word7Adj[6] |= !occF << 1;
    Word7Adj[6] |= !(edgeBits & 0b001010) << 2;
    Word7Adj[6] |= !(occOrLFBfb & 64) << 3;
    Sparse[6] = true;
    Word7Adj[6] |= N2 << 4; // Right and Back
  }

  // ------ bit 7  ------  3+7 = 10 bits, never sparse
  Word7Adj[7] = N3; // Right and Back and Top
}

//============================================================================

}  // namespace pcc