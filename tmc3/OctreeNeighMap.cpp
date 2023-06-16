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
static const int LUTLineardx[9] = { 1, 1,  1,  1, 0,  0,  0, -1, -1 };
static const int LUTLineardy[9] = { 1, 0,  0, -1, 1,  1,  -1, 1,  0 };
static const int LUTLineardz[9] = { 0, 1, -1,  0, 1, -1,  1,  0,  1 };
GeometryNeighPattern
makeGeometryNeighPattern(
  bool adjacent_child_contextualization_enabled_flag,
  const Vec3<int32_t>& position,
  int codedAxesPrevLvl,
  const MortonMap3D& occupancyAtlas,
  const bool& planarEligibleKOctreeDepth)
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
  GeometryNeighPattern gnp = { neighPattern, 0, 0, 0, {0, 0, 0, 0, 0, 0, 0}, 0, false };

  if ((!gnp.neighPattern && !planarEligibleKOctreeDepth) || !adjacent_child_contextualization_enabled_flag)
    return gnp;

  if (gnp.neighPattern) {
    if (gnp.neighPattern & 2)
      gnp.adjNeighOcc[0] = occupancyAtlas.getChildOcc(x - 1, y, z);

    if (gnp.neighPattern & 4)
      gnp.adjNeighOcc[1] = occupancyAtlas.getChildOcc(x, y - 1, z);

    if (gnp.neighPattern & 16)
      gnp.adjNeighOcc[2] = occupancyAtlas.getChildOcc(x, y, z - 1);
  }


  if (planarEligibleKOctreeDepth) {
    if (
      x > 0 && x < cubeSizeMinusOne && y > 0 && y < cubeSizeMinusOne && z > 0
      && z < cubeSizeMinusOne) {
      if (occupancyAtlas.get(x - 1, y - 1, z, sx, sy, sz))
        gnp.adjNeighOcc[3] = occupancyAtlas.getChildOcc(x - 1, y - 1, z);
      if (occupancyAtlas.get(x - 1, y, z - 1, sx, sy, sz))
        gnp.adjNeighOcc[4] = occupancyAtlas.getChildOcc(x - 1, y, z - 1);
      if (occupancyAtlas.get(x, y - 1, z - 1, sx, sy, sz))
        gnp.adjNeighOcc[5] = occupancyAtlas.getChildOcc(x, y - 1, z - 1);
      if (occupancyAtlas.get(x - 1, y - 1, z - 1, sx, sy, sz))
        gnp.adjNeighOcc[6] = occupancyAtlas.getChildOcc(x - 1, y - 1, z - 1);
    }
    else {
      if (occupancyAtlas.getWithCheck(x - 1, y - 1, z, sx, sy, sz))
        gnp.adjNeighOcc[3] = occupancyAtlas.getChildOcc(x - 1, y - 1, z);
      if (occupancyAtlas.getWithCheck(x - 1, y, z - 1, sx, sy, sz))
        gnp.adjNeighOcc[4] = occupancyAtlas.getChildOcc(x - 1, y, z - 1);
      if (occupancyAtlas.getWithCheck(x, y - 1, z - 1, sx, sy, sz))
        gnp.adjNeighOcc[5] = occupancyAtlas.getChildOcc(x, y - 1, z - 1);
      if (occupancyAtlas.getWithCheck(x - 1, y - 1, z - 1, sx, sy, sz))
        gnp.adjNeighOcc[6] = occupancyAtlas.getChildOcc(x - 1, y - 1, z - 1);
    }
    bool& neighAvailable = gnp.neighOccuValid;
    neighAvailable = false;
    for (int idx = 0; idx < 7 && !neighAvailable; ++idx)
      neighAvailable |= gnp.adjNeighOcc[idx];

    //< If the advanced neighbor occupancy is not available
    if (!neighAvailable) {
      auto& neighOccu = gnp.neighborOccu;
      neighOccu = (!!(gnp.neighPattern & 1) << 11)
        | (!!(gnp.neighPattern & 8) << 10) | (!!(gnp.neighPattern & 32) << 9);


      if (x > 0 && x < cubeSizeMinusOne && y > 0 && y < cubeSizeMinusOne && z > 0
        && z < cubeSizeMinusOne)
        for (int n = 0; n < 9; n++) {
          neighOccu |= occupancyAtlas.get(
            x + LUTLineardx[n], y + LUTLineardy[n], z + LUTLineardz[n], sx, sy, sz) << n;
        }

      else
        for (int n = 0; n < 9; n++) {
          neighOccu |= occupancyAtlas.getWithCheck(
            x + LUTLineardx[n], y + LUTLineardy[n], z + LUTLineardz[n], sx, sy, sz) << n;
        }
    }
  }

  return gnp;
}

//----------------------------------------------------------------------------

static const int LUTdx[20] =
  { -1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 };
static const int LUTdy[20] =
  { -1,-1,-1, 0, 0, 1, 1, 1,-1,-1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1 };
static const int LUTdz[20] =
  { -1, 0, 1,-1, 1,-1, 0, 1,-1, 1,-1, 1,-1, 0, 1,-1, 1,-1, 0, 1 };

void prepareGeometryAdvancedNeighPattern(
  OctreeNeighours& octreeNeighours,
  const GeometryNeighPattern& gnp,
  const Vec3<int32_t>& position,
  int atlasShift,
  const MortonMap3D& occupancyAtlas,
  const bool& planarEligibleKOctreeDepth)
{
  // prepare neighbours
  const int neighPattern = gnp.neighPattern;
  const int mask = occupancyAtlas.cubeSize() - 1;
  const int32_t x = position[0] & mask;
  const int32_t y = position[1] & mask;
  const int32_t z = position[2] & mask;

  const int sx = atlasShift & 4 ? 1 : 0;
  const int sy = atlasShift & 2 ? 1 : 0;
  const int sz = atlasShift & 1 ? 1 : 0;

  const int cubeSizeMinusOne = mask;

  // prepare 20 neighbours
  int neighb20 = 0;
  if (x > 0 && x < cubeSizeMinusOne && y > 0 && y < cubeSizeMinusOne && z > 0
      && z < cubeSizeMinusOne)
    for (int n = 0; n < 20; n++)
      neighb20 |= occupancyAtlas.get(
        x + LUTdx[n], y + LUTdy[n], z + LUTdz[n], sx, sy, sz) << n;
  else
    for (int n = 0; n < 20; n++)
      neighb20 |= occupancyAtlas.getWithCheck(
        x + LUTdx[n], y + LUTdy[n], z + LUTdz[n], sx, sy, sz) << n;

  octreeNeighours.neighb20 = neighb20;

  // ----- neighbours  FLB -----
  octreeNeighours.occLeft = gnp.adjNeighOcc[0];
  octreeNeighours.occFront = gnp.adjNeighOcc[1];
  octreeNeighours.occBottom = gnp.adjNeighOcc[2];

  octreeNeighours.occL = octreeNeighours.occLeft >> 4;
  octreeNeighours.occF =
    ((octreeNeighours.occFront >> 2) & 3)
    | ((octreeNeighours.occFront >> 4) & 12);
  octreeNeighours.occB =
    ((octreeNeighours.occBottom >> 1) & 1)
    | ((octreeNeighours.occBottom >> 2) & 2)
    | ((octreeNeighours.occBottom >> 3) & 4)
    | ((octreeNeighours.occBottom >> 4) & 8);
  octreeNeighours.occOrLFBfb =
    octreeNeighours.occLeft
    | octreeNeighours.occFront
    | octreeNeighours.occBottom;

  // ----- neighbours  LB, FB, LF -----
  octreeNeighours.edgeBits = 0;
  if ((neighb20 >> 3) & 1) {
    int occLB = planarEligibleKOctreeDepth ? gnp.adjNeighOcc[4] : occupancyAtlas.getChildOcc(x - 1, y, z - 1);
    octreeNeighours.edgeBits = ((occLB & 32) >> 5) | ((occLB & 128) >> 6);
  }

  if ((neighb20 >> 8) & 1) {
    int occFB = planarEligibleKOctreeDepth ? gnp.adjNeighOcc[5] : occupancyAtlas.getChildOcc(x, y - 1, z - 1);
    octreeNeighours.edgeBits |= ((occFB & 8) >> 1) | ((occFB & 128) >> 4);
  }

  if ((neighb20 >> 1) & 1) {
    int occLF = planarEligibleKOctreeDepth ? gnp.adjNeighOcc[3] : occupancyAtlas.getChildOcc(x - 1, y - 1, z);
    octreeNeighours.edgeBits |= (occLF & 0b11000000) >> 2;
  }

  // neighPattern & 0b101001  -> right  back top
  octreeNeighours.N3 =
    ((neighPattern >> 3) & 4)
    | ((neighPattern >> 2) & 2)
    | (neighPattern & 1);
  // -> right  back
  octreeNeighours.N2 = octreeNeighours.N3 & 3;
  // LFB pattern
  octreeNeighours.neighPatternLFB =
    ((neighPattern & 0b110) >> 1)
    | ((neighPattern & 16) >> 2);
}

// ---------------------------------------------------------------------------

inline int getBit(int w, const int n)
{
  return (w >> n) & 1;
}

inline int getBit(int w, const int n1, const int n2)
{
  return ((w >> n1 - 1) & 2) | ((w >> n2) & 1);
}

inline int getBit(int w, const int n1, const int n2, const int n3)
{
  return ((w >> n1 - 2) & 4) | ((w >> n2 - 1) & 2) | ((w >> n3) & 1);
}

inline int getBit(
  int w, const int n1, const int n2, const int n3, const int n4)
{
  return ((w >> n1 - 3) & 8) | ((w >> n2 - 2) & 4) | ((w >> n3 - 1) & 2)
    | ((w >> n4) & 1);
}

// ---------------------------------------------------------------------------

static const int LUTNN[16] = { 0,1,1,2 , 1,2,2,3 , 1,2,2,3 , 2,3,3,4};

// ------------------------------------ bit 0 --------------------------------

void
makeGeometryAdvancedNeighPattern0(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;

  int NN =
    LUTNN[octreeNeighours.occL]
    + LUTNN[octreeNeighours.occF]
    + LUTNN[octreeNeighours.occB];

  if (NN > 1) {
    int NLFB =
      !!octreeNeighours.occL
      + !!octreeNeighours.occF
      + !!octreeNeighours.occB;

    if (NLFB == 3) {
      // put tag at the head as it is the most important info
      infoNeigh = 0b100 << 16;
      // face
      infoNeigh |= (octreeNeighours.occB & 1) << 15;
      infoNeigh |= (octreeNeighours.occF & 1) << 14;
      infoNeigh |= (octreeNeighours.occL & 1) << 13;
      //edge
      infoNeigh |= (octreeNeighours.occB & 0b110) << 11 - 1;
      infoNeigh |= (octreeNeighours.occF & 0b110) << 9 - 1;
      infoNeigh |= (octreeNeighours.occL & 0b110) << 7 - 1;

      infoNeigh |= octreeNeighours.N3 << 4;
      infoNeigh |= getBit(N20, 8, 3, 1, 0);
    }
    else {
      if (NLFB == 2) {
        if (octreeNeighours.occL && octreeNeighours.occB) {
          infoNeigh = 0b101 << 16;
          //face
          infoNeigh |= (octreeNeighours.occB & 1) << 15;
          infoNeigh |= (octreeNeighours.occL & 1) << 14;
          //edge
          infoNeigh |= (octreeNeighours.occB & 0b110) << 12 - 1;
          infoNeigh |= (octreeNeighours.occL & 0b110) << 10 - 1;
          //vertex
          infoNeigh |= !(octreeNeighours.occB & 8) << 9;
          infoNeigh |= !(octreeNeighours.occL & 8) << 8;

          infoNeigh |= !(octreeNeighours.N3 & 2) << 7; //Back
        }
        if (octreeNeighours.occF && octreeNeighours.occB) {
          infoNeigh = 0b110 << 16;
          //face
          infoNeigh |= (octreeNeighours.occB & 1) << 15;
          infoNeigh |= (octreeNeighours.occF & 1) << 14;
          //edge
          infoNeigh |= (octreeNeighours.occB & 0b110) << 12 - 1;
          infoNeigh |= (octreeNeighours.occF & 0b110) << 10 - 1;
          //vertex
          infoNeigh |= !(octreeNeighours.occB & 8) << 9;
          infoNeigh |= !(octreeNeighours.occF & 8) << 8;

          infoNeigh |= !(octreeNeighours.N3 & 1) << 7; // Right
        }
        if (octreeNeighours.occL && octreeNeighours.occF) {
          infoNeigh = 0b111 << 16;
          //face
          infoNeigh |= (octreeNeighours.occF & 1) << 15;
          infoNeigh |= (octreeNeighours.occL & 1) << 14;
          // edge
          infoNeigh |= (octreeNeighours.occF & 0b110) << 12 - 1;
          infoNeigh |= (octreeNeighours.occL & 0b110) << 10 - 1;
          //vertex
          infoNeigh |= !(octreeNeighours.occF & 8) << 9;
          infoNeigh |= !(octreeNeighours.occL & 8) << 8;

          infoNeigh |= !(octreeNeighours.N3 & 4) << 7; //Top
        }
      }
      else //NLFB == 1
      {
        if (octreeNeighours.occL) {
          infoNeigh = 0b000 << 16;
          //face
          infoNeigh |= (octreeNeighours.occL & 1) << 15;
          //edge
          infoNeigh |= (octreeNeighours.occL & 0b110) << 13 - 1;
          //vertex
          infoNeigh |= !(octreeNeighours.occL & 8) << 12;
          infoNeigh |= (octreeNeighours.edgeBits & 0b001100) << 10 - 2;
        }
        else if (octreeNeighours.occF) {
          infoNeigh = 0b001 << 16;
          //face
          infoNeigh |= (octreeNeighours.occF & 1) << 15;
          //edge
          infoNeigh |= (octreeNeighours.occF & 0b110) << 13 - 1;
          //vertex
          infoNeigh |= !(octreeNeighours.occF & 8) << 12;
          infoNeigh |= (octreeNeighours.edgeBits & 0b000011) << 10;
        }
        else { //if (octreeNeighours.occB) {
          infoNeigh = 0b010 << 16;
          //face
          infoNeigh |= (octreeNeighours.occB & 1) << 15;
          //edge
          infoNeigh |= (octreeNeighours.occB & 0b110) << 13 - 1;
          //vertex
          infoNeigh |= !(octreeNeighours.occB & 8) << 12;
          infoNeigh |= (octreeNeighours.edgeBits & 0b110000) << 10 - 4;
        }
        infoNeigh |= octreeNeighours.N3 << 7;
      }

      infoNeigh |= getBit(N20, 8, 3, 1, 0) << 3;
      infoNeigh |= getBit(N20, 18, 19, 11);
    }

    Sparse = false;
    ctx1 = infoNeigh >> 13;
    ctx2 = infoNeigh & 0x1FFF;
  }
  else {
    int neighPatternLFB = octreeNeighours.neighPatternLFB;
    if (NN) {
      if (octreeNeighours.occL) {
        infoNeigh = 1 << 14;
        infoNeigh |= !(octreeNeighours.occL & 1) << 13;  // face
        infoNeigh |= !(neighPatternLFB & 4) << 12;  //edge
        infoNeigh |= !(neighPatternLFB & 2) << 11;  //edge
      }
      else if (octreeNeighours.occF) {
        infoNeigh = 2 << 14;
        infoNeigh |= !(octreeNeighours.occF & 1) << 13;  // face
        infoNeigh |= !(neighPatternLFB & 4) << 12; //edge
        infoNeigh |= !(neighPatternLFB & 1) << 11;  //edge
      }
      else { //if (octreeNeighours.occB) {
        infoNeigh = 3 << 14;
        infoNeigh |= !(octreeNeighours.occB & 1) << 13; // face
        infoNeigh |= !(neighPatternLFB & 2) << 12; //edge
        infoNeigh |= !(neighPatternLFB & 1) << 11; //edge
      }
    }
    else {
      infoNeigh = 0 << 14;
      infoNeigh |= neighPatternLFB << 11;
    }

    infoNeigh |= getBit(N20, 1, 3) << 9;
    infoNeigh |= getBit(N20, 8, 0) << 7;

    if (neighPatternLFB) {
      if (octreeNeighours.occOrLFBfb & 1) {
        infoNeigh |= 1 << 6;
        infoNeigh |= (octreeNeighours.occBottom & 1) << 5;
        infoNeigh |= (octreeNeighours.occFront & 1) << 4;
        infoNeigh |= (octreeNeighours.occLeft & 1) << 3;
      }
      else
      {
        infoNeigh |= !octreeNeighours.edgeBits << 5;
        infoNeigh |=
          ((octreeNeighours.occLeft & 4)
          || (octreeNeighours.occFront & 2)
          || (octreeNeighours.occBottom & 4)) << 4;
        infoNeigh |=
          ((octreeNeighours.occLeft & 2)
          || (octreeNeighours.occFront & 16)
          || (octreeNeighours.occBottom & 16)) << 3;
      }
    }
    else {
      infoNeigh |= !(octreeNeighours.edgeBits & 0b110000) << 6;
      infoNeigh |= !(octreeNeighours.edgeBits & 0b001100) << 5;
      infoNeigh |= !(octreeNeighours.edgeBits & 0b000011) << 4;
      // bit 3 is unused
    }

    infoNeigh |= getBit(N20, 18, 19, 11);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

// ------------------------------------ bit 1 --------------------------------

void
makeGeometryAdvancedNeighPattern1(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;

  if (octreeNeighours.occF) {
    //face
    infoNeigh = (occupancy & 1) << 18;
    infoNeigh |= !(octreeNeighours.occF & 0b0010) << 17;
    infoNeigh |= !octreeNeighours.occL << 16;

    if (octreeNeighours.occL) {
      //face
      infoNeigh |= !(octreeNeighours.occL & 0b0010) << 15;
      infoNeigh |= !(octreeNeighours.N3 & 4) << 14; //Top
      //edge
      infoNeigh |= !(octreeNeighours.occF & 0b0001) << 13;
      infoNeigh |= !(octreeNeighours.occF & 0b1000) << 12;
      infoNeigh |= !(octreeNeighours.occL & 0b0001) << 11;
      infoNeigh |= !(octreeNeighours.occL & 0b1000) << 10;
      //vertex
      infoNeigh |= !(octreeNeighours.occF & 0b0100) << 9;
      infoNeigh |= !(octreeNeighours.occL & 0b0100) << 8;

      infoNeigh |= (octreeNeighours.N3 & 1) << 7; // Right

      infoNeigh |= getBit(N20, 9, 4, 1, 2) << 3;
    }
    else {
      infoNeigh |= !(octreeNeighours.N3 & 4) << 15; //Top
      // edge
      infoNeigh |= !(octreeNeighours.occF & 0b0001) << 14;
      infoNeigh |= !(octreeNeighours.occF & 0b1000) << 13;
      //vertex
      infoNeigh |= !(octreeNeighours.occF & 0b0100) << 12;

      infoNeigh |= getBit(N20, 9, 4, 1, 2) << 8;

      infoNeigh |= !(octreeNeighours.occBottom & 2) << 7;
      infoNeigh |= !(octreeNeighours.occFront & 2) << 6;
      infoNeigh |= !(octreeNeighours.occLeft & 2) << 5;
      infoNeigh |= (octreeNeighours.N3 & 3) << 3; // Right Back
    }

    infoNeigh |= getBit(N20, 11, 16, 19);

    Sparse = false;
    ctx1 = infoNeigh >> 13;
    ctx2 = infoNeigh & 0x1FFF;
  }
  else {
    //face
    infoNeigh = (occupancy & 1) << 18;
    infoNeigh |= !(octreeNeighours.occL & 0b0010) << 17;
    infoNeigh |= !(octreeNeighours.N3 & 4) << 16; //Top
    //edge
    infoNeigh |= !(octreeNeighours.occL & 0b0001) << 15;
    infoNeigh |= !(octreeNeighours.occL & 0b1000) << 14;
    //vertex
    infoNeigh |= !(octreeNeighours.occL & 0b0100) << 13;

    infoNeigh |= (octreeNeighours.N3 & 1) << 12; // Right

    infoNeigh |= getBit(N20, 1, 4) << 10;
    infoNeigh |= getBit(N20, 9, 2) << 8;

    if (octreeNeighours.occOrLFBfb & 2) {
      infoNeigh |= 1 << 7;
      infoNeigh |= !(octreeNeighours.occBottom & 2) << 6;
      infoNeigh |= !(octreeNeighours.occFront & 2) << 5;
      infoNeigh |= !(octreeNeighours.occLeft & 2) << 4;
    }
    else {
      infoNeigh |= !(octreeNeighours.edgeBits & 0b110101) << 6;
      infoNeigh |=
        ((octreeNeighours.occLeft & 8)
        || (octreeNeighours.occFront & 32)) << 5;
      infoNeigh |=
        ((octreeNeighours.occLeft & 1)
        || (octreeNeighours.occFront & 1)) << 4;
    }

    infoNeigh |= !octreeNeighours.occB << 3;
    infoNeigh |= getBit(N20, 11, 16, 19);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

// ------------------------------------ bit 2 --------------------------------

void
makeGeometryAdvancedNeighPattern2(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;

  if (octreeNeighours.occB) {
    // face
    infoNeigh = (occupancy & 1) << 18;
    infoNeigh |= !(octreeNeighours.occB & 0b0010) << 17;
    infoNeigh |= !octreeNeighours.occL << 16;

    if (octreeNeighours.occL) {
      //face
      infoNeigh |= !(octreeNeighours.occL & 0b0100) << 15;
      infoNeigh |= !(octreeNeighours.N3 & 2) << 14; //back
      //edge
      infoNeigh |= !(occupancy & 2) << 13;
      infoNeigh |= !(octreeNeighours.occB & 0b1000) << 12;
      infoNeigh |= !(octreeNeighours.occL & 0b1000) << 11;
      infoNeigh |= !(octreeNeighours.occL & 0b0001) << 10;
      infoNeigh |= !(octreeNeighours.occB & 0b0001) << 9;

      infoNeigh |= getBit(N20, 10, 6, 3) << 6;

      //vertex
      infoNeigh |= !(octreeNeighours.occB & 0b0100) << 5;
      infoNeigh |= !(octreeNeighours.occL & 0b0010) << 4;
    }
    else {
      //face
      infoNeigh |= !(octreeNeighours.N3 & 2) << 15; //back
      //edge
      infoNeigh |= !(occupancy & 2) << 14;
      infoNeigh |= !(octreeNeighours.occB & 0b0001) << 13;
      infoNeigh |= !(octreeNeighours.occB & 0b1000) << 12;
      //vertex
      infoNeigh |= !(octreeNeighours.occB & 0b0100) << 11;

      infoNeigh |= getBit(N20, 10, 6, 3) << 8;

      //others
      infoNeigh |= !(octreeNeighours.N3 & 4)  << 7;
      infoNeigh |= !(octreeNeighours.occLeft & 4) << 6;
      infoNeigh |= !(octreeNeighours.occBottom & 4) << 5;
      infoNeigh |= !(octreeNeighours.occFront & 4) << 4;
    }

    infoNeigh |= getBit(N20, 0) << 3;
    infoNeigh |= getBit(N20, 18, 19, 11);

    Sparse = false;
    ctx1 = infoNeigh >> 13;
    ctx2 = infoNeigh & 0x1FFF;
  }
  else {
    // face
    infoNeigh = (occupancy & 1) << 18;
    infoNeigh |= !(octreeNeighours.occL & 0b0100) << 17;
    infoNeigh |= !(octreeNeighours.N3 & 2) << 16; //back
    //edge
    infoNeigh |= !(occupancy & 2) << 15;
    infoNeigh |= !(octreeNeighours.occL & 0b1000) << 14;
    infoNeigh |= !(octreeNeighours.occL & 0b0001) << 13;
    //vertex
    infoNeigh |= !(octreeNeighours.occL & 0b0010) << 12;

    infoNeigh |= getBit(N20, 3, 6, 10, 5) << 8;

    if (octreeNeighours.occOrLFBfb & 4) {
      infoNeigh |= 1 << 7;
      infoNeigh |= !(octreeNeighours.occLeft & 4) << 6;
      infoNeigh |= !(octreeNeighours.occBottom & 4) << 5;
      infoNeigh |= !(octreeNeighours.occFront & 4) << 4;
    }
    else {
      infoNeigh |=
        ((octreeNeighours.occLeft & 1)
        || (octreeNeighours.occBottom & 1)) << 6;
      infoNeigh |=
        ((octreeNeighours.occLeft & 8)
        || (octreeNeighours.occBottom & 64)) << 5;
      infoNeigh |= !(octreeNeighours.edgeBits & 0b000011) << 4;
    }

    infoNeigh |= !octreeNeighours.occF << 3;
    infoNeigh |= getBit(N20, 18, 19, 11);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

// ------------------------------------ bit 3 --------------------------------

void
makeGeometryAdvancedNeighPattern3(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;

  int NN = LUTNN[octreeNeighours.occL] + LUTNN[occupancy & 7];

  if (NN > 1) {
    // face
    infoNeigh = !(occupancy & 4) << 16;
    infoNeigh |= !(occupancy & 2) << 15;
    infoNeigh |= !(octreeNeighours.occL & 8 ) << 14;
    infoNeigh |= octreeNeighours.N3 << 11;
    // edge
    infoNeigh |= !(occupancy & 1) << 10;
    infoNeigh |= !(octreeNeighours.occL & 4 ) << 9;
    infoNeigh |= !(octreeNeighours.occL & 2 ) << 8;
    // vertex
    infoNeigh |= (octreeNeighours.occL & 1) << 7;

    infoNeigh |= getBit(N20, 11, 6, 4, 0) << 3;
    infoNeigh |= getBit(N20, 16, 19, 18);

    Sparse = false;
    ctx1 = infoNeigh >> 11;
    ctx2 = infoNeigh & 0x07FF;
  }
  else {
    // NN <= 1
    int occup = occupancy & 7;
    infoNeigh = !(occup) << 17;
    if (occup) { // occL is empty, at most 1 in occupancy & 7
      infoNeigh |= (!!occup + !!(occup >> 1) + !!(occup >> 2)) << 15;
    }
    else {
      // occL may not be empty;
      // empty and furthest neigh occupied together in 0
      infoNeigh |=
        (!!(octreeNeighours.occL >> 1)
        + !!(octreeNeighours.occL >> 2)
        + !!(octreeNeighours.occL >> 3)) << 15;
    }

    infoNeigh |= (octreeNeighours.N3 >> 1) << 13; // Back + Top
    infoNeigh |= getBit(N20, 4, 6, 11, 7) << 9;

    if (octreeNeighours.occOrLFBfb & 8) {
      infoNeigh |= 1 << 8;
      infoNeigh |= !(octreeNeighours.occBottom & 8) << 7;
      infoNeigh |= !(octreeNeighours.occFront & 8) << 6;
      infoNeigh |= !(octreeNeighours.occLeft & 8) << 5;
    }
    else {
      infoNeigh |= (octreeNeighours.occLeft & 0b110) << 5;
      infoNeigh |= !(octreeNeighours.edgeBits & 0b110010) << 5;
    }
    infoNeigh |= !octreeNeighours.occB << 4;
    infoNeigh |= !octreeNeighours.occF << 3;

    infoNeigh |= getBit(N20, 18, 19, 16);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

// ------------------------------------ bit 4 --------------------------------

void
makeGeometryAdvancedNeighPattern4(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;
  const int occupancyLeft = occupancy & 15;

  int NN =
    LUTNN[occupancyLeft]
    + LUTNN[octreeNeighours.occF]
    + LUTNN[octreeNeighours.occB];

  if (NN > 1) {
    int NLFB =
      !!occupancyLeft
      + !!octreeNeighours.occF
      + !!octreeNeighours.occB;

    if (NLFB == 3) {
      // put tag at the head as it is the most important info
      infoNeigh = 0b1000 << 15;
      //face
      infoNeigh |= !(octreeNeighours.occB & 4) << 17;
      infoNeigh |= !(octreeNeighours.occF & 4) << 16;
      infoNeigh |= (occupancyLeft & 1) << 15;
      infoNeigh |= !(octreeNeighours.N3 & 1) << 14; //right
      //edge
      infoNeigh |= !(octreeNeighours.occB & 1) << 13;
      infoNeigh |= !(octreeNeighours.occB & 8) << 12;
      infoNeigh |= !(octreeNeighours.occF & 1) << 11;
      infoNeigh |= !(octreeNeighours.occF & 8) << 10;
      infoNeigh |= ! (occupancyLeft & 2) << 9;
      infoNeigh |= ! (occupancyLeft & 4) << 8;
      //vertex
      infoNeigh |= !(octreeNeighours.occB & 2) << 7;
      infoNeigh |= !(octreeNeighours.occF & 2) << 6;

      infoNeigh |= (octreeNeighours.N3 >> 1) << 4; // back top
      infoNeigh |= getBit(N20, 15, 13, 8, 12);
    }
    else if (NLFB == 2) {
      if (occupancyLeft && octreeNeighours.occB) {
        infoNeigh = 0b0100 << 15;
        //face
        infoNeigh |= !(octreeNeighours.occB & 4) << 14;
        infoNeigh |= !(occupancyLeft & 1) << 13;
        infoNeigh |= !(octreeNeighours.N3 & 1) << 12; //right
        //edge
        infoNeigh |= !(octreeNeighours.occB & 1) << 11;
        infoNeigh |= !(octreeNeighours.occB & 8) << 10;
        infoNeigh |= !(occupancyLeft & 2) << 9;
        infoNeigh |= !(occupancyLeft & 4) << 8;
        //vertex
        infoNeigh |= !(octreeNeighours.occB & 2) << 7;
        infoNeigh |= !(occupancyLeft & 8) << 6;
      }
      else if (octreeNeighours.occF && octreeNeighours.occB) {
        infoNeigh = 0b0101 << 15;
        //face
        infoNeigh |= !(octreeNeighours.occB & 4) << 14;
        infoNeigh |= !(octreeNeighours.occF & 4) << 13;
        infoNeigh |= !(octreeNeighours.N3 & 1) << 12; //right
        //edge
        infoNeigh |= !(octreeNeighours.occB & 1) << 11;
        infoNeigh |= !(octreeNeighours.occB & 8) << 10;
        infoNeigh |= !(octreeNeighours.occF & 1) << 9;
        infoNeigh |= !(octreeNeighours.occF & 8) << 8;
        //vertex
        infoNeigh |= !(octreeNeighours.occB & 2) << 7;
        infoNeigh |= !(octreeNeighours.occF & 2) << 6;
      }
      else { //if (occupancyLeft && octreeNeighours.occF) {
        infoNeigh = 0b0110 << 15;
        //face
        infoNeigh |= !(octreeNeighours.occF & 4) << 14;
        infoNeigh |= !(occupancyLeft & 1) << 13;
        infoNeigh |= !(octreeNeighours.N3 & 1) << 12; //right
        //edge
        infoNeigh |= !(octreeNeighours.occF & 1) << 11;
        infoNeigh |= !(octreeNeighours.occF & 8) << 10;
        infoNeigh |= ! (occupancyLeft & 2) << 9;
        infoNeigh |= ! (occupancyLeft & 4) << 8;
        //vertex
        infoNeigh |= !(octreeNeighours.occF & 2) << 7;
        infoNeigh |= !(occupancyLeft & 8) << 6;
      }

      infoNeigh |= getBit(N20, 15, 13, 8) << 3;
      infoNeigh |= getBit(N20, 12, 16, 18);
    }
    else //NLFB == 1
    {
      if (occupancyLeft) {
        infoNeigh = 0b0000 << 15;
        infoNeigh |= (occupancyLeft & 1) << 14; //face
        infoNeigh |= !(octreeNeighours.N3 & 1) << 13; //right
        infoNeigh |= (occupancyLeft & 0b110) << 11-1; //edge
        infoNeigh |= !(occupancyLeft & 8) << 10; //vertex
        infoNeigh |= (octreeNeighours.edgeBits & 0b001100) << 8-2;
      }
      else if (octreeNeighours.occF) {
        infoNeigh = 0b0001 << 15;
        infoNeigh |= !(octreeNeighours.occF & 0b0100) << 14; //face
        infoNeigh |= !(octreeNeighours.N3 & 1) << 13; //right
        infoNeigh |= !(octreeNeighours.occF & 0b0001) << 12; //edge
        infoNeigh |= !(octreeNeighours.occF & 0b1000) << 11; //edge
        infoNeigh |= !(octreeNeighours.occF & 0b0010) << 10; //vertex
        infoNeigh |= (octreeNeighours.edgeBits & 0b000011) << 8;
      }
      else { //if (octreeNeighours.occB) {
        infoNeigh = 0b0010 << 15;
        infoNeigh |= !(octreeNeighours.occB & 0b0100) << 14;  //face
        infoNeigh |= !(octreeNeighours.N3 & 1) << 12; //right
        infoNeigh |= !(octreeNeighours.occB & 0b0001) << 12; //edge
        infoNeigh |= !(octreeNeighours.occB & 0b1000) << 11; //edge
        infoNeigh |= !(octreeNeighours.occB & 0b0010) << 10; //vertex
        infoNeigh |= (octreeNeighours.edgeBits & 0b110000) << 8-4;
      }

      infoNeigh |= (octreeNeighours.N3 >> 1 ) << 6;

      infoNeigh |= getBit(N20, 15, 13, 8) << 3;
      infoNeigh |= getBit(N20, 12, 16, 18);
    }

    Sparse = false;
    ctx1 = infoNeigh >> 13;
    ctx2 = infoNeigh & 0x1FFF;
  }
  else {
    int neighPatternLFB = octreeNeighours.neighPatternLFB;
    if (NN) {
      if (occupancyLeft) {
        infoNeigh = 1 << 14;
        infoNeigh |= !(occupancyLeft & 1) << 13; //face
        infoNeigh |= !(neighPatternLFB & 4) << 12; //edge
        infoNeigh |= !(neighPatternLFB & 2) << 11; //edge
      }
      else if (octreeNeighours.occF) {
        infoNeigh = 2 << 14;
        infoNeigh |= !(octreeNeighours.occF & 1) << 13; //face
        infoNeigh |= !(neighPatternLFB & 4) << 12; //edge
        infoNeigh |= !(neighPatternLFB & 1) << 11; //edge
      }
      else { //if (octreeNeighours.occB) {
        infoNeigh = 3 << 14;
        infoNeigh |= !(octreeNeighours.occB & 1) << 13; //face
        infoNeigh |= !(neighPatternLFB & 2) << 12; //edge
        infoNeigh |= !(neighPatternLFB & 1) << 11; //edge
      }
    }
    else {
      infoNeigh = 0 << 14;
      infoNeigh |= neighPatternLFB << 11;
    }

    infoNeigh |= getBit(N20, 8, 13, 15, 12) << 7;

    if (neighPatternLFB) {
      if (octreeNeighours.occOrLFBfb & 16) {
        infoNeigh |= 1 << 6;
        infoNeigh |= !(octreeNeighours.occBottom & 16) << 5;
        infoNeigh |= !(octreeNeighours.occFront & 16) << 4;
        infoNeigh |= !(octreeNeighours.occLeft & 16) << 3;
      }
      else
      {
        infoNeigh |= !octreeNeighours.edgeBits << 5;
        infoNeigh |=
          ((octreeNeighours.occLeft & 64)
          || (octreeNeighours.occFront & 8)
          || (octreeNeighours.occBottom & 8)) << 4;
        infoNeigh |= ((octreeNeighours.occLeft & 32)
          || (octreeNeighours.occFront & 64)
          || (octreeNeighours.occBottom & 32)) << 3;
      }
    }
    else {
      infoNeigh |= !(octreeNeighours.edgeBits & 0b110000) << 6;
      infoNeigh |= !(octreeNeighours.edgeBits & 0b001100) << 5;
      infoNeigh |= !(octreeNeighours.edgeBits & 0b000011) << 4;
      // bit3 is unused
    }

    infoNeigh |= getBit(N20, 16, 18, 19);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

// ------------------------------------ bit 5 --------------------------------

void
makeGeometryAdvancedNeighPattern5(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;
  const int occupancyLeft = occupancy & 15;

  if (octreeNeighours.occF) {
    //face
    infoNeigh = ((occupancy >> 4) & 1) << 18;
    infoNeigh |= !(octreeNeighours.occF & 0b1000) << 17;
    infoNeigh |= !occupancyLeft << 16;

    if (occupancyLeft) {
      //face
      infoNeigh |= !(occupancyLeft & 0b0010) << 15;
      infoNeigh |= !(octreeNeighours.N3 & 4) << 14; //Top
      infoNeigh |= !(octreeNeighours.N3 & 1) << 13; //Right
      //edge
      infoNeigh |= !(octreeNeighours.occF & 0b0010) << 12;
      infoNeigh |= !(octreeNeighours.occF & 0b0100) << 11;
      infoNeigh |= !(occupancyLeft & 0b0001) << 10;
      infoNeigh |= !(occupancyLeft & 0b1000) << 9;
      //vertex
      infoNeigh |= !(octreeNeighours.occF & 0b0001) << 8;
      infoNeigh |= !(occupancyLeft & 0b0100) << 7;

      infoNeigh |= getBit(N20, 16, 13, 9, 14) << 3;
    }
    else {
      //face
      infoNeigh |= !(octreeNeighours.N3 & 4) << 15; //Top
      infoNeigh |= !(octreeNeighours.N3 & 1) << 14; //Right
      //edge
      infoNeigh |= !(octreeNeighours.occF & 0b0010) << 13;
      infoNeigh |= !(octreeNeighours.occF & 0b0100) << 12;
      //vertex
      infoNeigh |= !(octreeNeighours.occF & 0b0001) << 11;

      infoNeigh |= getBit(N20, 16, 13, 9, 14) << 7;

      infoNeigh |= !(octreeNeighours.occBottom & 32) << 6;
      infoNeigh |= !(octreeNeighours.occFront & 32) << 5;
      infoNeigh |= !(octreeNeighours.occLeft & 32) << 4;
      infoNeigh |= !(octreeNeighours.N3 & 2) << 3; //back
    }

    infoNeigh |= getBit(N20, 18, 19, 11);

    Sparse = false;
    ctx1 = infoNeigh >> 13;
    ctx2 = infoNeigh & 0x1FFF;
  }
  else {
    // face
    infoNeigh = !((occupancy >> 4) & 1) << 18;
    infoNeigh |= !(occupancyLeft & 0b0010) << 17;
    infoNeigh |= !(octreeNeighours.N3 & 4) << 16; //Top
    infoNeigh |= !(octreeNeighours.N3 & 1) << 15; //Right
    //edge
    infoNeigh |= !(occupancyLeft & 0b0001) << 14;
    infoNeigh |= !(occupancyLeft & 0b1000) << 13;
    //vertex
    infoNeigh |= !(octreeNeighours.occL & 0b0100) << 12;

    infoNeigh |= getBit(N20, 9, 13, 16, 14) << 8;

    if (octreeNeighours.occOrLFBfb & 32) {
      infoNeigh |= 1 << 7;
      infoNeigh |= !(octreeNeighours.occBottom & 32) << 6;
      infoNeigh |= !(octreeNeighours.occFront & 32) << 5;
      infoNeigh |= !(octreeNeighours.occLeft & 32) << 4;
    }
    else {
      infoNeigh |= !(octreeNeighours.edgeBits & 0b111100) << 6;
      infoNeigh |=
        ((octreeNeighours.occLeft & 128)
        || (octreeNeighours.occFront & 2)) << 5;
      infoNeigh |=
        ((octreeNeighours.occLeft & 16)
        || (octreeNeighours.occFront & 16)) << 4;
    }

    infoNeigh |= !octreeNeighours.occB << 3;
    infoNeigh |= getBit(N20, 18, 19, 11);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

// ------------------------------------ bit 6 --------------------------------

void
makeGeometryAdvancedNeighPattern6(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;
  const int occupancyLeft = occupancy & 15;

  if (octreeNeighours.occB) {
    //face
    infoNeigh = !((occupancy >> 4) & 1) << 18;
    infoNeigh |= !(octreeNeighours.occB & 0b1000) << 17;
    infoNeigh |= !occupancyLeft << 16;

    if (occupancyLeft) {
      //face
      infoNeigh |= !(occupancyLeft & 0b0100) << 15;
      infoNeigh |= !(octreeNeighours.N3 & 1) << 14; //right
      infoNeigh |= !(octreeNeighours.N3 & 2) << 13; //back
      //edge
      infoNeigh |= !((occupancy >> 4) & 2) << 12;
      infoNeigh |= !(octreeNeighours.occB & 0b0010) << 11;
      infoNeigh |= !(occupancyLeft & 0b0001) << 10;
      infoNeigh |= !(occupancyLeft & 0b1000) << 9;
      infoNeigh |= !(octreeNeighours.occB & 0b0100) << 8;

      infoNeigh |= getBit(N20, 18, 15, 10) << 5;

      //vertex
      infoNeigh |= !(octreeNeighours.occB & 0b0001) << 4;
      infoNeigh |= !(occupancyLeft & 0b0010) << 3;

      infoNeigh |= getBit(N20, 17) << 2;
      infoNeigh |= getBit(N20, 0) << 1;
      infoNeigh |= getBit(N20, 11) << 0;
    }
    else {
      //face
      infoNeigh |= !(octreeNeighours.N3 & 2) << 15; //back
      infoNeigh |= !(octreeNeighours.N3 & 1) << 14; //right
      //edge
      infoNeigh |= !((occupancy >> 4) & 2) << 13;
      infoNeigh |= !(octreeNeighours.occB & 0b0010) << 12;
      infoNeigh |= !(octreeNeighours.occB & 0b0100) << 11;
      //vertex
      infoNeigh |= !(octreeNeighours.occB & 0b0001) << 10;

      infoNeigh |= !(octreeNeighours.occLeft & 64) << 9;
      infoNeigh |= !(octreeNeighours.occBottom & 64) << 8;
      infoNeigh |= !(octreeNeighours.occFront & 64) << 7;

      infoNeigh |= getBit(N20, 18, 15, 10, 17) << 3;
      infoNeigh |= getBit(N20, 0) << 2;
      infoNeigh |= getBit(N20, 11, 19);
    }

    Sparse = false;
    ctx1 = infoNeigh >> 13;
    ctx2 = infoNeigh & 0x1FFF;
  }
  else {
    //face
    infoNeigh = !((occupancy >> 4) & 1) << 18;
    infoNeigh |= !(occupancyLeft & 0b0100) << 17;
    infoNeigh |= !(octreeNeighours.N3 & 1) << 16; //right
    //edge
    infoNeigh |= !((occupancy >> 4) & 2) << 15;
    infoNeigh |= !(occupancyLeft & 0b1000) << 14;
    infoNeigh |= !(occupancyLeft & 0b0001) << 13;
    //vertex
    infoNeigh |= !(occupancyLeft & 0b0010) << 12;

    infoNeigh |= getBit(N20, 17, 18, 15, 10) << 8;

    if (octreeNeighours.occOrLFBfb & 64) {
      infoNeigh |= 1<<7;
      infoNeigh |= !(octreeNeighours.occLeft & 64) << 6;
      infoNeigh |= !(octreeNeighours.occBottom & 64) << 5;
      infoNeigh |= !(octreeNeighours.occFront & 64) << 4;
    }
    else {
      infoNeigh |=
        ((octreeNeighours.occLeft & 1)
        || (octreeNeighours.occBottom & 1)) << 6;
      infoNeigh |=
        ((octreeNeighours.occLeft & 8)
        || (octreeNeighours.occBottom & 64)) << 5;
      infoNeigh |= !(octreeNeighours.edgeBits & 0b000011) << 4;
    }

    infoNeigh |= !octreeNeighours.occF << 3;
    infoNeigh |= getBit(N20, 19, 16, 11);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

// ------------------------------------ bit 7 --------------------------------

void
makeGeometryAdvancedNeighPattern7(
  OctreeNeighours& octreeNeighours,
  int occupancy,
  int& ctx1,
  int& ctx2,
  bool& Sparse)
{
  int infoNeigh = 0;
  const int N20 = octreeNeighours.neighb20;
  const int occupancyLeft = occupancy & 15;

  int NN = LUTNN[occupancyLeft] + LUTNN[(occupancy >> 4) & 7];

  // NN cannot be 0! otherwise 7 bit empty end b7 is inferred to 1
  if (NN > 1) {
    // face opposition
    infoNeigh = !((occupancy >> 4) & 4) << 16;
    infoNeigh |= !((occupancy >> 4) & 2) << 15;
    infoNeigh |= !(occupancyLeft & 8 ) << 14;
    infoNeigh |= octreeNeighours.N3 << 11;
    // edge opposition
    infoNeigh |= !((occupancy >> 4) & 1) << 10;
    infoNeigh |= getBit(N20, 11) << 9;
    infoNeigh |= !(occupancyLeft & 4 ) << 8;
    infoNeigh |= getBit(N20, 16) << 7;
    infoNeigh |= !(occupancyLeft & 2 ) << 6;
    infoNeigh |= getBit(N20, 18) << 5;
    // vertex opposion
    infoNeigh |= (occupancyLeft & 1 ) << 4;
    infoNeigh |= getBit(N20, 19) << 3;

    // more neighbours
    infoNeigh |= getBit(N20, 0) << 2;
    //infoNeigh |= getBit(N20, 17) << 1;
    infoNeigh |= getBit(N20, 17, 10);

    Sparse = false;
    ctx1 = infoNeigh >> 11;
    ctx2 = infoNeigh & 0x07FF;
  }
  else {
    //NN==1  (can't be zero otherwise bit 7 is infered upstream to 1 )
    int occup = (occupancy >> 4) & 7;
    infoNeigh = !occup << 17;
    if (occup) { // occupancyLeft is empty,  1 in occup & 7
      infoNeigh |= (!!occup + !!(occup >> 1) + !!(occup >> 2)) << 15;
      infoNeigh |= !(octreeNeighours.N3 & 2) << 14; // Back
    }
    else { // occupancyLeft is not empty!!!: otherwise 7 bit empty end b7 is 1
      infoNeigh |=
        (!!(occupancyLeft >> 1)
        + !!(occupancyLeft >> 2)
        + !!(occupancyLeft >> 3)) << 15;
      infoNeigh |= !(octreeNeighours.N3 & 1) << 14; // Right
    }
    infoNeigh |= !(octreeNeighours.N3 & 4) << 13; // Top
    infoNeigh |= getBit(N20, 11, 16, 18, 19) << 9;

    if (octreeNeighours.occOrLFBfb & 128) {
      infoNeigh |= 1 << 8;
      infoNeigh |= !(octreeNeighours.occLeft & 128) << 7;
      infoNeigh |= !(octreeNeighours.occFront & 128) << 6;
      infoNeigh |= !(octreeNeighours.occBottom & 128) << 5;
    }
    else {
      infoNeigh |= (octreeNeighours.occLeft & 0b01100000) << 1;
      infoNeigh |=
        (octreeNeighours.occF & 0b0011
        || octreeNeighours.occB & 0b0110) << 5;  // LB FB LF
    }
    infoNeigh |= !octreeNeighours.occB << 4;
    infoNeigh |= !octreeNeighours.occF << 3;

    infoNeigh |= getBit(N20, 7, 17, 10);

    Sparse = true;
    ctx1 = infoNeigh >> 12;
    ctx2 = infoNeigh & 0x0FFF;
  }
}

//============================================================================

}  // namespace pcc