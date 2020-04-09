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

namespace pcc {

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
  bool occupancySkip)
{
  static const uint8_t childMasks[] = {
    0xf0 /* x-1 */, 0xcc /* y-1 */, 0xaa /* z-1 */
  };

  uint32_t patternBit = 1 << (1 << neighIdx);
  uint8_t childMask = childMasks[neighIdx];

  // conversions between neighbour occupancy and adjacency:
  //  x: >> 4, y: >> 2, z: >> 1
  int adjacencyShift = 4 >> neighIdx;

  if (occupancySkip) {
    childMask ^= 0xff;
    adjacencyShift = 0;
  }

  if (gnp.neighPattern & patternBit) {
    uint8_t child_occ = occupancyAtlas.getChildOcc(x, y, z);
    uint8_t child_unocc = ~child_occ;
    child_occ &= childMask;
    if (!child_occ) {
      /* neighbour is falsely occupied */
      gnp.neighPattern ^= patternBit;
    } else {
      child_occ >>= adjacencyShift;
      gnp.adjacencyGt1 |= gnp.adjacencyGt0 & child_occ;
      gnp.adjacencyGt0 |= child_occ;
    }

    // map of children with any unoccupied adjacent child
    gnp.adjacencyUnocc |= (child_unocc & childMask) >> adjacencyShift;
  }

  return gnp;
}

//----------------------------------------------------------------------------

GeometryNeighPattern
makeGeometryNeighPattern(
  bool adjacent_child_contextualization_enabled_flag,
  const Vec3<int32_t>& position,
  const int atlasShift,
  const MortonMap3D& occupancyAtlas)
{
  const int mask = occupancyAtlas.cubeSize() - 1;
  const int cubeSizeMinusOne = mask;
  const int32_t x = position[0] & mask;
  const int32_t y = position[1] & mask;
  const int32_t z = position[2] & mask;
  uint8_t neighPattern;
  if (atlasShift == 0) {
    if (
      x > 0 && x < cubeSizeMinusOne && y > 0 && y < cubeSizeMinusOne && z > 0
      && z < cubeSizeMinusOne) {
      neighPattern = occupancyAtlas.get(x + 1, y, z);
      neighPattern |= occupancyAtlas.get(x - 1, y, z) << 1;
      neighPattern |= occupancyAtlas.get(x, y - 1, z) << 2;
      neighPattern |= occupancyAtlas.get(x, y + 1, z) << 3;
      neighPattern |= occupancyAtlas.get(x, y, z - 1) << 4;
      neighPattern |= occupancyAtlas.get(x, y, z + 1) << 5;
    } else {
      neighPattern = occupancyAtlas.getWithCheck(x + 1, y, z);
      neighPattern |= occupancyAtlas.getWithCheck(x - 1, y, z) << 1;
      neighPattern |= occupancyAtlas.getWithCheck(x, y - 1, z) << 2;
      neighPattern |= occupancyAtlas.getWithCheck(x, y + 1, z) << 3;
      neighPattern |= occupancyAtlas.getWithCheck(x, y, z - 1) << 4;
      neighPattern |= occupancyAtlas.getWithCheck(x, y, z + 1) << 5;
    }
  } else {
    const int shiftX = (atlasShift & 4 ? 1 : 0);
    const int shiftY = (atlasShift & 2 ? 1 : 0);
    const int shiftZ = (atlasShift & 1 ? 1 : 0);

    if (
      x > 0 && x < cubeSizeMinusOne && y > 0 && y < cubeSizeMinusOne && z > 0
      && z < cubeSizeMinusOne) {
      neighPattern = occupancyAtlas.get(x + 1, y, z, shiftX, shiftY, shiftZ);
      neighPattern |= occupancyAtlas.get(x - 1, y, z, shiftX, shiftY, shiftZ)
        << 1;
      neighPattern |= occupancyAtlas.get(x, y - 1, z, shiftX, shiftY, shiftZ)
        << 2;
      neighPattern |= occupancyAtlas.get(x, y + 1, z, shiftX, shiftY, shiftZ)
        << 3;
      neighPattern |= occupancyAtlas.get(x, y, z - 1, shiftX, shiftY, shiftZ)
        << 4;
      neighPattern |= occupancyAtlas.get(x, y, z + 1, shiftX, shiftY, shiftZ)
        << 5;
    } else {
      neighPattern =
        occupancyAtlas.getWithCheck(x + 1, y, z, shiftX, shiftY, shiftZ);
      neighPattern |=
        occupancyAtlas.getWithCheck(x - 1, y, z, shiftX, shiftY, shiftZ) << 1;
      neighPattern |=
        occupancyAtlas.getWithCheck(x, y - 1, z, shiftX, shiftY, shiftZ) << 2;
      neighPattern |=
        occupancyAtlas.getWithCheck(x, y + 1, z, shiftX, shiftY, shiftZ) << 3;
      neighPattern |=
        occupancyAtlas.getWithCheck(x, y, z - 1, shiftX, shiftY, shiftZ) << 4;
      neighPattern |=
        occupancyAtlas.getWithCheck(x, y, z + 1, shiftX, shiftY, shiftZ) << 5;
    }
  }

  // Above, the neighbour pattern corresponds directly to the six same
  // sized neighbours of the given node.
  // The patten is then refined by examining the available children
  // of the same neighbours.

  // NB: the process of updating neighpattern below also derives
  // the occupancy contextualisation bits.
  GeometryNeighPattern gnp = {neighPattern, 0, 0, 0};

  if (!adjacent_child_contextualization_enabled_flag)
    return gnp;

  if (x > 0)
    gnp = updatePatternFromNeighOccupancy(
      occupancyAtlas, x - 1, y, z, gnp, 0, !(atlasShift & 4));

  if (y > 0)
    gnp = updatePatternFromNeighOccupancy(
      occupancyAtlas, x, y - 1, z, gnp, 1, !(atlasShift & 2));

  if (z > 0)
    gnp = updatePatternFromNeighOccupancy(
      occupancyAtlas, x, y, z - 1, gnp, 2, !(atlasShift & 1));

  return gnp;
}

//============================================================================

}  // namespace pcc