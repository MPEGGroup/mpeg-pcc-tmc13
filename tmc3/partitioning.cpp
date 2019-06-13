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

#include "partitioning.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <stdint.h>

namespace pcc {

//============================================================================
// Determine whether half of slices are smaller than maxPoints

bool
halfQualified(const std::vector<Partition> slices, int maxPoints)
{
  int Qualified = 0;
  for (int i = 0; i < slices.size(); i++) {
    if (slices[i].pointIndexes.size() < maxPoints)
      Qualified++;
  }

  return ((double)Qualified / (double)slices.size()) > 0.5;
}

//============================================================================

template<typename T>
static int
longestAxis(const Box3<T>& curBox)
{
  int edgeAxis = 0;

  for (int i = 1; i < 3; i++) {
    T axisLength = curBox.max[i] - curBox.min[i];
    T longestLength = curBox.max[edgeAxis] - curBox.min[edgeAxis];
    if (axisLength > longestLength)
      edgeAxis = i;
  }

  return edgeAxis;
}

//----------------------------------------------------------------------------

template<typename T>
static int
shortestAxis(const Box3<T>& curBox)
{
  int edgeAxis = 0;

  for (int i = 1; i < 3; i++) {
    T axisLength = curBox.max[i] - curBox.min[i];
    T shortestLength = curBox.max[edgeAxis] - curBox.min[edgeAxis];
    if (axisLength < shortestLength)
      edgeAxis = i;
  }

  return edgeAxis;
}

//----------------------------------------------------------------------------
// Split point cloud into slices along the longest axis.
// numPartitions describes the number of slices to produce (if zero, the
// ratio of longest:shortest axis is used).
// No tile metadata is generated.

std::vector<Partition>
partitionByUniformGeom(
  const PartitionParams& params, const PCCPointSet3& cloud, int tileID)
{
  std::vector<Partition> slices;

  Box3<double> bbox = cloud.computeBoundingBox();

  int maxEdgeAxis = longestAxis(bbox);
  int maxEdge = bbox.max[maxEdgeAxis] - bbox.min[maxEdgeAxis];

  int minEdgeAxis = shortestAxis(bbox);
  int minEdge = bbox.max[minEdgeAxis] - bbox.min[minEdgeAxis];

  int sliceNum = minEdge ? (maxEdge / minEdge) : 1;
  int sliceSize = minEdge ? minEdge : maxEdge;

  while (1) {
    slices.clear();
    slices.resize(sliceNum);

    for (int i = 0; i < sliceNum; i++) {
      auto& slice = slices[i];
      slice.sliceId = i;
      slice.tileId = tileID;
      slice.origin = Vec3<int>{0};
    }

    for (int n = 0; n < cloud.getPointCount(); n++) {
      for (int p = sliceNum - 1; p >= 0; p--) {
        if (
          cloud[n][maxEdgeAxis]
          >= int(p * sliceSize + bbox.min[maxEdgeAxis])) {
          auto& slice = slices[p];

          slice.pointIndexes.push_back(n);
          break;
        }
      }
    }

    if (halfQualified(slices, params.sliceMaxPoints))
      break;

    sliceNum *= 2;
    sliceSize = maxEdge / sliceNum;
  }

  // Delete the slice that with no points
  slices.erase(
    std::remove_if(
      slices.begin(), slices.end(),
      [](const Partition& p) { return p.pointIndexes.empty(); }),
    slices.end());

  return slices;
}

//----------------------------------------------------------------------------
// Split point cloud into several parts according to octree depth.
// No tile metadata is generated.

std::vector<Partition>
partitionByOctreeDepth(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  int tileID,
  bool splitByDepth)
{
  std::vector<Partition> slices;

  // noting that there is a correspondence between point position
  // and octree node, calculate the position mask and shift required
  // to determine the node address for a point.
  Box3<double> bbox = cloud.computeBoundingBox();
  int maxBb = (int)std::max({bbox.max[0], bbox.max[1], bbox.max[2]});

  int cloudSizeLog2 = ceillog2(maxBb + 1);
  int depOctree = splitByDepth ? params.octreeDepth : 1;

  do {
    slices.clear();
    int posShift = cloudSizeLog2 - depOctree;
    int posMask = (1 << depOctree) - 1;

    // initially: number of points in each partition
    // then: mapping of partId to sliceId
    std::vector<int> partMap(1 << (3 * depOctree));

    // per-point indexes used for assigning to a partition
    std::vector<int> pointToPartId(cloud.getPointCount());

    // for each point, determine a partition based upon the position
    for (int i = 0, last = cloud.getPointCount(); i < last; i++) {
      int x = ((int(cloud[i].x()) >> posShift) & posMask) << (2 * depOctree);
      int y = ((int(cloud[i].y()) >> posShift) & posMask) << depOctree;
      int z = (int(cloud[i].z()) >> posShift) & posMask;
      int partId = x | y | z;
      partMap[partId]++;
      pointToPartId[i] = partId;
    }

    // generate slice mapping
    //  - allocate slice map storage and determine contiguous sliceIds
    //    NB: the sliceIds replace partPointCount.
    //  - map points to each slice.

    int numSlices =
      partMap.size() - std::count(partMap.begin(), partMap.end(), 0);
    slices.resize(numSlices);

    int sliceId = 0;
    for (auto& part : partMap) {
      if (!part)
        continue;

      auto& slice = slices[sliceId];
      slice.sliceId = sliceId;
      slice.tileId = tileID;
      slice.origin = Vec3<int>{0};
      slice.pointIndexes.reserve(part);
      part = sliceId++;
    }

    for (int i = 0, last = cloud.getPointCount(); i < last; i++) {
      int partId = pointToPartId[i];
      int sliceId = partMap[partId];
      slices[sliceId].pointIndexes.push_back(i);
    }

    if (halfQualified(slices, params.sliceMaxPoints))
      break;

    depOctree++;
  } while (!splitByDepth);

  return slices;
}

//=============================================================================
// Split point cloud into several tiles according to tileSize

std::vector<std::vector<int32_t>>
tilePartition(const PartitionParams& params, const PCCPointSet3& cloud)
{
  std::vector<std::vector<int32_t>> tilePartition;
  int tileSize = params.tileSize;

  // for each point determine the tile to which it belongs
  // let tile_origin = floor(pos / tile_size)
  // append pointIdx to tileMap[tile_origin]
  Box3<double> bbox = cloud.computeBoundingBox();
  int maxtileNum =
    std::ceil(std::max({bbox.max[0], bbox.max[1], bbox.max[2]}) / tileSize);
  int tileNumlog2 = ceillog2(maxtileNum);
  std::vector<int> partMap(1 << (3 * tileNumlog2));

  // per-point indexes used for assigning to a partition
  std::vector<uint64_t> pointToPartId(cloud.getPointCount());

  // for each point, determine a partition based upon the position
  std::vector<uint8_t> tilePos(3);
  uint64_t mortonTileID;
  for (int32_t i = 0, last = cloud.getPointCount(); i < last; i++) {
    mortonTileID = 0;
    for (int k = 0; k < 3; k++) {
      tilePos[k] = std::floor(cloud[i][k] / tileSize);
    }

    for (int p = 0; p < 8; p++) {
      mortonTileID |= ((tilePos[0] >> p) & 1) << (3 * p + 2);
      mortonTileID |= ((tilePos[1] >> p) & 1) << (3 * p + 1);
      mortonTileID |= ((tilePos[2] >> p) & 1) << (3 * p);
    }

    partMap[mortonTileID]++;
    pointToPartId[i] = mortonTileID;
  }

  // generate tile mapping
  // - allocate tile map storage and determine contiguous tileIds
  //   NB: the tileIds replace partPointCount.
  // - map points to each tile.

  int numTiles =
    partMap.size() - std::count(partMap.begin(), partMap.end(), 0);
  tilePartition.resize(numTiles);

  int tileId = 0;
  for (auto& part : partMap) {
    if (!part)
      continue;

    tilePartition[tileId].reserve(part);
    part = tileId++;
  }

  for (int i = 0, last = cloud.getPointCount(); i < last; i++) {
    int partId = pointToPartId[i];
    int tileId = partMap[partId];
    tilePartition[tileId].push_back(i);
  }

  return tilePartition;
}

//============================================================================

}  // namespace pcc
