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
#include <vector>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "hls.h"

namespace pcc {

//============================================================================

enum class PartitionMethod
{
  // Don't partition input
  kNone = 0,

  // Partition according to uniform geometry
  kUniformGeom = 2,

  // Partition according to the depth of octree
  kOctreeUniform = 3,

  // TBD
  kUniformSquare = 4,

  // Paritition into n-point slices
  kNpoints = 5,
};

//============================================================================

struct PartitionParams {
  // Method for partitioning a point cloud
  PartitionMethod method;

  // Depth of octree used in partitioning
  int octreeDepth;

  // Maximum number of reconstructed points per slice with Trisoup
  int sliceMaxPointsTrisoup;

  // Maximum number of points per slice
  int sliceMaxPoints;

  // Minimum number of points per slice
  int sliceMinPoints;

  // Baseline tile width. (0 => disabled)
  int tileSize;

  // use safe partitionning for Trisoup, aligning nodes to a global grid
  bool safeTrisoupPartionning;
};

//============================================================================

struct Partition {
  // The *_slice_id to encode this partition with
  int sliceId;

  // The identifier of the tileset that describes the bounding box of all
  // slices with the same tile_id.  A value of -1 indicates no tile_id
  // mapping.
  int tileId;

  // The value of geom_box_origin for this partition, using the
  // translated+scaled co-ordinate system.
  Vec3<int> origin;

  // Some metadata used by the partitioning process
  Vec3<int> location;

  // Point indexes of the source point cloud that form this partition.
  std::vector<int32_t> pointIndexes;
};

//----------------------------------------------------------------------------

struct PartitionSet {
  TileInventory tileInventory;
  std::vector<Partition> slices;
};

//============================================================================

std::vector<Partition> partitionByUniformGeom(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  int tileID,
  int paritionBoundaryLog2);

std::vector<Partition> partitionByOctreeDepth(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  int tileID,
  bool splitByDepth = false);

std::vector<Partition> partitionByUniformSquare(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  int tileID,
  int paritionBoundaryLog2);

std::vector<Partition> partitionByNpts(
  const PartitionParams& params, const PCCPointSet3& cloud, int tileID);

std::vector<Partition> partitionNone(
  const PartitionParams& params, const PCCPointSet3& cloud, int tileID);

//============================================================================

std::vector<std::vector<int32_t>>
tilePartition(const PartitionParams& params, const PCCPointSet3& cloud);

//============================================================================

void refineSlicesByAdjacentInfo(
  const PartitionParams& params,
  const PCCPointSet3& inputPointCloud,
  Vec3<int> sliceArrNum,
  std::vector<Partition>& slices,
  int partitionBoundary = 0);

//============================================================================

}  // namespace pcc
