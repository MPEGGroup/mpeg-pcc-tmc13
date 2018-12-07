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

namespace pcc {

//============================================================================

template<typename T>
static int
longestAxis(const PCCBox3<T>& curBox)
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
shortestAxis(const PCCBox3<T>& curBox)
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

PartitionSet
partitionByUniformGeom(const PCCPointSet3& cloud, int numPartitions)
{
  PartitionSet partitions;

  PCCBox3D bbox = cloud.computeBoundingBox();

  int maxEdgeAxis = longestAxis(bbox);
  int maxEdge = bbox.max[maxEdgeAxis] - bbox.min[maxEdgeAxis];

  int minEdgeAxis = shortestAxis(bbox);
  int minEdge = bbox.max[minEdgeAxis] - bbox.min[minEdgeAxis];
  int sliceNum;
  int sliceSize;

  if (!numPartitions) {
    sliceNum = maxEdge / minEdge;
    sliceSize = minEdge;
  } else {
    sliceNum = numPartitions;
    sliceSize = maxEdge / sliceNum;
  }

  partitions.slices.resize(sliceNum);

  for (int i = 0; i < sliceNum; i++) {
    auto& slice = partitions.slices[i];
    slice.sliceId = i;
    slice.tileId = -1;
    slice.origin = PCCVector3<int>{0};
  }

  for (int n = 0; n < cloud.getPointCount(); n++) {
    for (int p = sliceNum - 1; p >= 0; p--) {
      if (
        cloud[n][maxEdgeAxis] >= int(p * sliceSize + bbox.min[maxEdgeAxis])) {
        auto& slice = partitions.slices[p];

        slice.pointIndexes.push_back(n);
        break;
      }
    }
  }

  // Delete the slice that with no points
  partitions.slices.erase(
    std::remove_if(
      partitions.slices.begin(), partitions.slices.end(),
      [](const Partition& p) { return p.pointIndexes.empty(); }),
    partitions.slices.end());

  return partitions;
}

//----------------------------------------------------------------------------
// Split point cloud into several parts according to octree depth.
// No tile metadata is generated.

PartitionSet
partitionByOctreeDepth(const PCCPointSet3& cloud, int depOctree)
{
  PartitionSet partitions;

  // noting that there is a correspondence between point position
  // and octree node, calculate the position mask and shift required
  // to determine the node address for a point.
  PCCBox3<double> bbox = cloud.computeBoundingBox();
  int maxBb = (int)std::max({bbox.max[0], bbox.max[1], bbox.max[2]});

  int cloudSizeLog2 = ceillog2(maxBb + 1);
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
  partitions.slices.resize(numSlices);

  int sliceId = 0;
  for (auto& part : partMap) {
    if (!part)
      continue;

    auto& slice = partitions.slices[sliceId];
    slice.sliceId = sliceId;
    slice.tileId = -1;
    slice.origin = PCCVector3<int>{0};
    slice.pointIndexes.reserve(part);
    part = sliceId++;
  }

  for (int i = 0, last = cloud.getPointCount(); i < last; i++) {
    int partId = pointToPartId[i];
    int sliceId = partMap[partId];
    partitions.slices[sliceId].pointIndexes.push_back(i);
  }

  return partitions;
}

//============================================================================

}  // namespace pcc
