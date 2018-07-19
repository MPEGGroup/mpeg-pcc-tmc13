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
  const PCCVector3<uint32_t>& currentPosition,
  const int nodeSizeLog2,
  const pcc::ringbuf<PCCOctree3Node>& fifo,
  const pcc::ringbuf<PCCOctree3Node>::iterator& fifoCurrLvlEnd,
  MortonMap3D* occupancyAtlas,
  PCCVector3<uint32_t>* atlasOrigin)
{
  const uint32_t mask = (1 << occupancyAtlas->cubeSizeLog2()) - 1;
  const int shift = occupancyAtlas->cubeSizeLog2() + nodeSizeLog2;
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

    const uint32_t x = (it->pos[0] >> nodeSizeLog2) & mask;
    const uint32_t y = (it->pos[1] >> nodeSizeLog2) & mask;
    const uint32_t z = (it->pos[2] >> nodeSizeLog2) & mask;
    occupancyAtlas->setByte(x, y, z, it->siblingOccupancy);
  }
}

//----------------------------------------------------------------------------

uint32_t
makeGeometryNeighPattern(
  const PCCVector3<uint32_t>& position,
  const int nodeSizeLog2,
  const MortonMap3D& occupancyAtlas)
{
  const int mask = occupancyAtlas.cubeSize() - 1;
  const int cubeSizeMinusOne = mask;
  const int32_t x = (position[0] >> nodeSizeLog2) & mask;
  const int32_t y = (position[1] >> nodeSizeLog2) & mask;
  const int32_t z = (position[2] >> nodeSizeLog2) & mask;
  uint32_t neighPattern;
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
  return neighPattern;
}

//============================================================================

}  // namespace pcc