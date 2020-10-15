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

#include "PCCMath.h"
#include "geometry_octree.h"
#include "ringbuf.h"
#include "tables.h"

#include <cassert>
#include <memory>
#include <vector>

namespace pcc {

//============================================================================
// Provides a mapping of (x,y,z) co-ordinate to a bit flag.
//
// Internal representation uses a morton code to access an array of flags.
//
// Updates to the array are made byte-wise, allowing 8 flags (in morton order)
// to be stored in a single operation.

class MortonMap3D {
public:
  void resize(bool childOccupancyEnabled, const uint32_t cubeSizeLog2)
  {
    assert(cubeSizeLog2 < 10);
    const uint32_t halfCubeSizeLog2 = cubeSizeLog2 ? cubeSizeLog2 - 1 : 0;
    _cubeSizeLog2 = cubeSizeLog2;
    _cubeSize = 1 << cubeSizeLog2;
    _bufferSizeInBytes = 1 << (3 * cubeSizeLog2);
    _buffer.reset(new uint8_t[_bufferSizeInBytes]);
    if (childOccupancyEnabled)
      _childOccupancy.reset(new uint8_t[_bufferSizeInBytes << 3]);
    _updates.reserve(1 << 16);
  }

  int cubeSize() const { return _cubeSize; }
  int cubeSizeLog2() const { return _cubeSizeLog2; }

  // Removes all updates, zeros all map entries.  Does not affect child map.
  void clear();

  // Reverts all updates, zeroing affected map entries.
  // Only modified bytes are touched.
  void clearUpdates();

  void setByte(
    const int32_t x, const int32_t y, const int32_t z, const uint8_t value)
  {
    assert(
      x >= 0 && y >= 0 && z >= 0 && x < _cubeSize && y < _cubeSize
      && z < _cubeSize);
    if (value) {
      const uint32_t byteIndex = getByteIndex(x, y, z);
      _buffer[byteIndex] = value;
      _updates.push_back(byteIndex);
    }
  }

  uint32_t get(
    const int32_t x,
    const int32_t y,
    const int32_t z,
    const int shiftX,
    const int shiftY,
    const int shiftZ) const
  {
    assert(
      x >= 0 && y >= 0 && z >= 0 && x < _cubeSize && y < _cubeSize
      && z < _cubeSize);
    return (_buffer[getByteIndex(x >> shiftX, y >> shiftY, z >> shiftZ)]
            >> getBitIndex(shiftX ? x : 0, shiftY ? y : 0, shiftZ ? z : 0))
      & 1;
  }

  uint32_t getWithCheck(
    const int32_t x,
    const int32_t y,
    const int32_t z,
    const int shiftX,
    const int shiftY,
    const int shiftZ) const
  {
    if (
      x < 0 || x >= _cubeSize || y < 0 || y >= _cubeSize || z < 0
      || z >= _cubeSize) {
      return false;
    }
    return get(x, y, z, shiftX, shiftY, shiftZ);
  }

  uint32_t get(const int32_t x, const int32_t y, const int32_t z) const
  {
    assert(
      x >= 0 && y >= 0 && z >= 0 && x < _cubeSize && y < _cubeSize
      && z < _cubeSize);

    return get(x, y, z, 1, 1, 1);
  }

  uint32_t
  getWithCheck(const int32_t x, const int32_t y, const int32_t z) const
  {
    if (
      x < 0 || x >= _cubeSize || y < 0 || y >= _cubeSize || z < 0
      || z >= _cubeSize) {
      return false;
    }
    return getWithCheck(x, y, z, 1, 1, 1);
  }

  void setChildOcc(int32_t x, int32_t y, int32_t z, uint8_t childOccupancy)
  {
    _childOccupancy[getByteIndex(x, y, z)] = childOccupancy;
  }

  uint8_t getChildOcc(int32_t x, int32_t y, int32_t z) const
  {
    uint8_t childOccupancy = _childOccupancy[getByteIndex(x, y, z)];
    return childOccupancy;
  }

private:
  int32_t getBitIndex(const int32_t x, const int32_t y, const int32_t z) const
  {
    return (z & 1) + ((y & 1) << 1) + ((x & 1) << 2);
  }

  uint32_t
  getByteIndex(const int32_t x, const int32_t y, const int32_t z) const
  {
    return kMortonCode256X[x] | kMortonCode256Y[y] | kMortonCode256Z[z];
  }

  int _cubeSize = 0;
  int _cubeSizeLog2 = 0;

  uint32_t _bufferSizeInBytes = 0;
  std::unique_ptr<uint8_t[]> _buffer;

  // A list of indexes in _buffer that are dirty
  std::vector<uint32_t> _updates;

  // Child occupancy values
  std::unique_ptr<uint8_t[]> _childOccupancy;
};

//============================================================================

struct GeometryNeighPattern {
  // Mask indicating presence of neigbours of the corresponding tree node
  //    32 8 (y)
  //     |/
  //  2--n--1 (x)
  //    /|
  //   4 16 (z)
  uint8_t neighPattern;

  // mask indicating the number of external child neighbours
  uint8_t adjacencyGt0;
  uint8_t adjacencyGt1;

  // mask indicating unoccupied external child neighbours
  uint8_t adjacencyUnocc;

  // occupancy map of {x-1, y-1, z-1} neighbours.
  uint8_t adjNeighOcc[3];
};

//============================================================================
// determine the occupancy pattern of the six neighbours of the node at
// @position.  If @adjacent_child_contextualization_enabled_flag is true,
// the occupancy state of previously coded neighbours is used to refine
// the neighbour pattern and derive external adjacency counts for each child.
GeometryNeighPattern makeGeometryNeighPattern(
  bool adjacent_child_contextualization_enabled_flag,
  const Vec3<int32_t>& currentPosition,
  int codedAxesPrevLvl,
  int codedAxesCurLvl,
  const MortonMap3D& occupancyAtlas);

// populate (if necessary) the occupancy atlas with occupancy information
// from @fifo.
void updateGeometryOccupancyAtlas(
  const Vec3<int32_t>& position,
  const int atlasShift,
  const ringbuf<PCCOctree3Node>& fifo,
  const ringbuf<PCCOctree3Node>::iterator& fifoCurrLvlEnd,
  MortonMap3D* occupancyAtlas,
  Vec3<int32_t>* atlasOrigin);

void updateGeometryOccupancyAtlasOccChild(
  const Vec3<int32_t>& pos,
  uint8_t childOccupancy,
  MortonMap3D* occupancyAtlas);

}  // namespace pcc
