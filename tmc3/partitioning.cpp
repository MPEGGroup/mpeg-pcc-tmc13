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

struct CM_Node {  //node == partition
  int cnt;
  int idx;
  std::vector<int32_t> pointCloudIndex;
};

struct CM_Nodes {  //nodes == merged partition
  int total;
  double xEvg;
  double yEvg;
  std::vector<CM_Node> nodes;
};

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

std::vector<Partition>
partitionByNpts(const PartitionParams& params, const PCCPointSet3& cloud)
{
  std::vector<Partition> slices;
  int sliceMaxPts = params.sliceMaxPoints;
  int numSlices = (cloud.getPointCount() + sliceMaxPts - 1) / sliceMaxPts;

  for (int sliceId = 0; sliceId < numSlices; sliceId++) {
    slices.emplace_back();
    auto& slice = slices.back();
    slice.sliceId = sliceId;
    slice.tileId = 0;

    // generate the point range
    int firstPointIdx = sliceMaxPts * sliceId;
    int numPoints =
      std::min(int(cloud.getPointCount()) - firstPointIdx, sliceMaxPts);
    slice.pointIndexes.resize(numPoints);
    for (int i = 0; i < numPoints; i++)
      slice.pointIndexes[i] = firstPointIdx + i;

    // work out the origin
    auto bbox = cloud.computeBoundingBox(
      slice.pointIndexes.begin(), slice.pointIndexes.end());

    for (int k = 0; k < 3; k++)
      slice.origin[k] = int(bbox.min[k]);
  }

  return slices;
}

//----------------------------------------------------------------------------
// Split point cloud into slices along the longest axis.
// numPartitions describes the number of slices to produce (if zero, the
// ratio of longest:shortest axis is used).
// No tile metadata is generated.

std::vector<Partition>
partitionByUniformGeom(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  int tileID,
  int partitionBoundaryLog2)
{
  std::vector<Partition> slices;

  Box3<int32_t> bbox = cloud.computeBoundingBox();

  int maxEdgeAxis = longestAxis(bbox);
  int maxEdge = bbox.max[maxEdgeAxis] - bbox.min[maxEdgeAxis];

  int minEdgeAxis = shortestAxis(bbox);
  int minEdge = bbox.max[minEdgeAxis] - bbox.min[minEdgeAxis];

  int sliceNum = minEdge ? (maxEdge / minEdge) : 1;
  int sliceSize = minEdge ? minEdge : maxEdge;

  // In order to avoid issues with trisoup, don't partition points within
  // a trisoup node, otherwise there will be issues fitting triangles.
  int partitionBoundary = 1 << partitionBoundaryLog2;
  if (sliceSize % partitionBoundary) {
    sliceSize = (1 + sliceSize / partitionBoundary) * partitionBoundary;
  }

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
    if (sliceSize % partitionBoundary) {
      sliceSize = (1 + sliceSize / partitionBoundary) * partitionBoundary;
    }
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

std::vector<Partition>
partitionByUniformSquare(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  int tileID,
  int partitionBoundaryLog2)
{
  std::vector<Partition> slices;

  Box3<int32_t> bbox = cloud.computeBoundingBox();

  int maxEdgeAxis = longestAxis(bbox);
  int maxEdge = bbox.max[maxEdgeAxis] - bbox.min[maxEdgeAxis];

  int minEdgeAxis = shortestAxis(bbox);
  int minEdge = (bbox.max[minEdgeAxis] - bbox.min[minEdgeAxis]);

  int midEdgeAxis = 3 - maxEdgeAxis - minEdgeAxis;
  int midEdge = bbox.max[midEdgeAxis] - bbox.min[midEdgeAxis];

  int firstSliceNum = maxEdge / minEdge + 1;
  int secondSliceNum = midEdge / minEdge + 1;
  int sliceNum = firstSliceNum * secondSliceNum;
  int sliceSize = minEdge ? minEdge : maxEdge;

  int partitionBoundary = 1 << partitionBoundaryLog2;
  if (sliceSize % partitionBoundary) {
    sliceSize = (1 + sliceSize / partitionBoundary) * partitionBoundary;
  }

  while (1) {
    slices.clear();
    slices.resize(sliceNum);

    std::vector<Partition> slicesFirstAxis;
    slicesFirstAxis.clear();
    slicesFirstAxis.resize(firstSliceNum);

    int count = 0;
    for (int i = 0; i < firstSliceNum; i++) {
      for (int j = 0; j < secondSliceNum; j++) {
        auto& slice = slices[i * secondSliceNum + j];
        slice.sliceId = count;
        slice.tileId = tileID;
        slice.origin = Vec3<int>{0};
        count++;
      }
    }
    int p = 0;
    for (int n = 0; n < cloud.getPointCount(); n++) {
      p = int(cloud[n][maxEdgeAxis]) / sliceSize;
      auto& slice = slicesFirstAxis[p];
      slice.pointIndexes.push_back(n);
    }

    int q = 0;
    for (int s = 0; s < slicesFirstAxis.size(); s++) {
      auto& sliceRough = slicesFirstAxis[s];
      for (int n = 0; n < sliceRough.pointIndexes.size(); n++) {
        q = int(cloud[sliceRough.pointIndexes[n]][midEdgeAxis]) / sliceSize;
        auto& slice = slices[s * secondSliceNum + q];
        slice.pointIndexes.push_back(sliceRough.pointIndexes[n]);
      }
    }
    slicesFirstAxis.clear();

    break;
  }

  int count = 0;
  for (int i = 0; i < slices.size(); i++) {
    auto& slice = slices[i];
    slice.sliceId = count;
    count++;
  }
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
  Box3<int32_t> bbox = cloud.computeBoundingBox();
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
  Box3<int32_t> bbox = cloud.computeBoundingBox();
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

//=============================================================================

Vec3<int>
minOrigin(const Vec3<int> a, const Vec3<int> b)
{
  Vec3<int> newOrigin;
  for (int i = 0; i < 3; i++) {
    newOrigin[i] = (a[i] < b[i]) ? a[i] : b[i];
  }
  return newOrigin;
}

//----------------------------------------------------------------------------
// get the axis of the longest edge

int
maxEdgeAxis(const PCCPointSet3& cloud, std::vector<int32_t>& sliceIndexes)
{
  int maxEdge = 0;
  int maxAxis = 0;
  for (int i = 0; i < 3; i++) {
    int maxpoint, minpoint;
    maxpoint = minpoint = cloud[sliceIndexes[0]][i];

    for (int k = 1; k < sliceIndexes.size(); k++) {
      if (cloud[sliceIndexes[k]][i] < minpoint)
        minpoint = cloud[sliceIndexes[k]][i];
      else if (cloud[sliceIndexes[k]][i] > maxpoint)
        maxpoint = cloud[sliceIndexes[k]][i];
    }

    if (maxEdge < (maxpoint - minpoint)) {
      maxEdge = maxpoint - minpoint;
      maxAxis = i;
    }
  }

  return maxAxis;
}

//============================================================================
// evenly split slice into several partitions no larger than maxPoints

void
splitSlice(
  CM_Nodes& splitingSlice,
  std::vector<CM_Nodes>& newlist,
  const PCCPointSet3& cloud,
  int maxPoints)
{
  auto& aIndexes = splitingSlice.nodes[0].pointCloudIndex;

  // Split along the longest edge at the median point
  int splitAxis = maxEdgeAxis(cloud, aIndexes);
  std::stable_sort(
    aIndexes.begin(), aIndexes.end(), [&](int32_t a, int32_t b) {
      return cloud[a][splitAxis] < cloud[b][splitAxis];
    });

  int numSplit = std::ceil((double)aIndexes.size() / (double)maxPoints);
  int splitsize = aIndexes.size() / numSplit;
  std::vector<Partition> splitPartitions;
  splitPartitions.resize(numSplit);

  for (int i = 0; i < numSplit - 1; i++) {
    auto& indexes = splitPartitions[i].pointIndexes;
    indexes.insert(
      indexes.begin(), aIndexes.begin() + i * splitsize,
      aIndexes.begin() + (i + 1) * splitsize);
  }
  auto& indexes = splitPartitions[numSplit - 1].pointIndexes;
  indexes.insert(
    indexes.begin(), aIndexes.begin() + (numSplit - 1) * splitsize,
    aIndexes.end());

  newlist.resize(numSplit);
  for (int i = 0; i < numSplit; i++) {
    auto& Indexes = splitPartitions[i].pointIndexes;
    newlist[i].xEvg = -1;
    newlist[i].yEvg = -1;
    newlist[i].total = Indexes.size();
    newlist[i].nodes.push_back(CM_Node());
    newlist[i].nodes[0].pointCloudIndex.resize(Indexes.size());
    newlist[i].nodes[0].pointCloudIndex = Indexes;
    newlist[i].nodes[0].cnt = Indexes.size();
  }
}

//----------------------------------------------------------------------------
// combine the two slices into one

void
mergePartition(Partition& a, const Partition& b)
{
  auto& aIndexes = a.pointIndexes;
  auto& bIndexes = b.pointIndexes;
  aIndexes.insert(aIndexes.end(), bIndexes.begin(), bIndexes.end());
  a.origin = minOrigin(a.origin, b.origin);
}

//=============================================================================
// first split slices still having too many points(more than maxPoints)
// then merge slices with too few points(less than minPoints)

void
refineSlicesByAdjacentInfo(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  std::vector<Partition>& slices)
{
  int maxPoints = params.sliceMaxPoints;
  int minPoints = params.sliceMinPoints;

  //initialize
  const auto& bBox = cloud.computeBoundingBox();
  int minEdgeAxis = shortestAxis(bBox);
  int minlength = (bBox.max[minEdgeAxis] - bBox.min[minEdgeAxis]);

  int maxEdgeAxis = longestAxis(bBox);
  int maxEdge = bBox.max[maxEdgeAxis] - bBox.min[maxEdgeAxis];

  int midEdgeAxis = 3 - maxEdgeAxis - minEdgeAxis;
  int midEdge = bBox.max[midEdgeAxis] - bBox.min[midEdgeAxis];

  int xNum = floor(maxEdge / minlength) + 1;
  int yNum = floor(midEdge / minlength) + 1;

  std::vector<std::vector<Partition>> slicePartition;
  std::vector<std::vector<int>> listCnt;
  slicePartition.clear();
  slicePartition.resize(xNum);
  listCnt.resize(xNum);
  for (int i = 0; i < xNum; i++) {
    slicePartition[i].resize(yNum);
    listCnt[i].resize(yNum);
  }

  for (int i = 0; i < xNum; i++) {
    for (int j = 0; j < yNum; j++) {
      int idx = i * yNum + j;
      if (idx > (slices.size() - 1))
        continue;
      slicePartition[i][j] = slices[i * yNum + j];
    }
  }

  std::vector<CM_Nodes> list;
  list.clear();
  list.resize(xNum * yNum);
  for (int i = 0; i < xNum * yNum; i++) {
    list[i].nodes.resize(1);
    list[i].xEvg = i / yNum;
    list[i].yEvg = i % yNum;
    list[i].total = slicePartition[i / yNum][i % yNum].pointIndexes.size();
    list[i].nodes[0].cnt =
      slicePartition[i / yNum][i % yNum].pointIndexes.size();
    list[i].nodes[0].idx = i;
  }

  for (int i = 0; i < list.size(); i++) {
    for (int n = 0; n < list[i].nodes.size(); n++) {
      list[i].nodes[n].pointCloudIndex.resize(list[i].total);
      for (int idx = 0; idx < list[i].total; idx++) {
        list[i].nodes[n].pointCloudIndex[idx] =
          slicePartition[i / yNum][i % yNum].pointIndexes[idx];
      }
    }
  }

  //list erase
  for (int i = 0; i < list.size();) {
    if (list[i].total == 0) {
      list.erase(list.begin() + i);
    } else {
      i++;
    }
  }

  //initial sort
  int min_idx;
  for (int i = 0; i < list.size() - 1; i++) {
    min_idx = i;
    for (int j = i + 1; j < list.size(); j++) {
      if (list[min_idx].total > list[j].total)
        min_idx = j;
    }
    std::swap(list[min_idx], list[i]);
  }

  std::vector<CM_Nodes> newlist;
  std::vector<CM_Nodes> newSlice;
  std::vector<CM_Nodes> tmplist = list;
  int listNum = list.size();
  for (int i = 0; i < listNum; i++) {
    if (tmplist[i].total > maxPoints) {
      splitSlice(list[i], newlist, cloud, maxPoints);
      for (int j = 0; j < newlist.size(); j++) {
        newSlice.push_back(newlist[j]);
      }
    }
  }

  tmplist.erase(
    std::remove_if(
      tmplist.begin(), tmplist.end(),
      [=](const CM_Nodes& p) { return p.total > maxPoints; }),
    tmplist.end());

  list = tmplist;
  tmplist.clear();
  for (int i = 0; i < list.size(); i++) {
    for (int n = 0; n < list[i].nodes.size(); n++) {
      int xIdx = list[i].nodes[n].idx / yNum;  //partition index x
      int yIdx = list[i].nodes[n].idx % yNum;  //partition index y
      listCnt[xIdx][yIdx] = i + 1;
    }
  }

  //find adjacent slice
  for (int i = 0; i < list.size();) {
    int minLidx = -1;
    double mindist = maxPoints;
    for (int n = 0; n < list[i].nodes.size(); n++) {
      int xIdxP = list[i].nodes[n].idx / yNum;
      int yIdxP = list[i].nodes[n].idx % yNum;
      for (int f = 0; f < 4; f++) {
        int xIdxPTmp = 0;
        int yIdxPTmp = 0;
        if (f == 0) {
          xIdxPTmp = xIdxP - 1;  //left
          yIdxPTmp = yIdxP;
          if (xIdxPTmp < 0)
            continue;
        } else if (f == 1) {
          xIdxPTmp = xIdxP + 1;  //right
          yIdxPTmp = yIdxP;
          if (xIdxPTmp >= xNum)
            continue;
        } else if (f == 2) {
          xIdxPTmp = xIdxP;
          yIdxPTmp = yIdxP - 1;  //down
          if (yIdxPTmp < 0)
            continue;
        } else if (f == 3) {
          xIdxPTmp = xIdxP;
          yIdxPTmp = yIdxP + 1;  //up
          if (yIdxPTmp >= yNum)
            continue;
        }

        if (
          (listCnt[xIdxPTmp][yIdxPTmp] == 0)
          || (listCnt[xIdxPTmp][yIdxPTmp] == i + 1)) {
          continue;
        }

        int lIdxTmp = listCnt[xIdxPTmp][yIdxPTmp] - 1;  //list index
        double distTmp = sqrt(
          pow((list[i].xEvg - xIdxPTmp), 2)
          + pow((list[i].yEvg - yIdxPTmp), 2));
        if (mindist >= distTmp) {
          if ((list[lIdxTmp].total + list[i].total) < maxPoints) {
            mindist = distTmp;
            minLidx = lIdxTmp;
          }
        }
      }
    }

    if (minLidx <= 0) {
      i++;
    } else {
      int pre_num = list[i].nodes.size();
      list[i].nodes.resize(pre_num + list[minLidx].nodes.size());

      for (int m = pre_num; m < list[i].nodes.size(); m++) {
        list[i].nodes[m].cnt = list[minLidx].nodes[m - pre_num].cnt;
        list[i].nodes[m].idx = list[minLidx].nodes[m - pre_num].idx;
      }
      list[i].total += list[minLidx].total;
      list.erase(list.begin() + minLidx);
      for (int j = i; j < list.size() - 1; j++) {
        if (list[j].total > list[j + 1].total)
          std::swap(list[j], list[j + 1]);
      }

      for (int j = i; j < list.size(); j++) {
        double xSum = 0;
        double ySum = 0;
        for (int n = 0; n < list[j].nodes.size(); n++) {
          int xIdx = list[j].nodes[n].idx / yNum;
          int yIdx = list[j].nodes[n].idx % yNum;
          xSum += xIdx;
          ySum += yIdx;
          listCnt[xIdx][yIdx] = j + 1;
        }
        list[j].xEvg = xSum / list[j].nodes.size();
        list[j].yEvg = ySum / list[j].nodes.size();
      }
    }
  }

  //split
  std::vector<Partition> refinedSlice;
  refinedSlice.clear();
  refinedSlice.resize(list.size() + newSlice.size());
  for (int i = 0; i < list.size(); i++) {
    for (int n = 0; n < list[i].nodes.size(); n++) {
      int xIdxP = list[i].nodes[n].idx / yNum;
      int yIdxP = list[i].nodes[n].idx % yNum;
      mergePartition(refinedSlice[i], slicePartition[xIdxP][yIdxP]);
    }
  }

  for (int i = list.size(); i < list.size() + newSlice.size(); i++) {
    refinedSlice[i].pointIndexes.resize(0);
    refinedSlice[i].pointIndexes.insert(
      refinedSlice[i].pointIndexes.end(),
      newSlice[i - list.size()].nodes[0].pointCloudIndex.begin(),
      newSlice[i - list.size()].nodes[0].pointCloudIndex.end());
  }
  slices = refinedSlice;
  for (int i = 0; i < slices.size(); i++) {
    auto& slice = slices[i];
    slice.sliceId = i;
    slice.tileId = -1;
  }
}

//=============================================================================

std::vector<Partition>::iterator
splitSlice(
  const PCCPointSet3& cloud,
  std::vector<Partition>& slices,
  std::vector<Partition>::iterator toBeSplit,
  int maxPoints)
{
  auto& sliceA = (*toBeSplit);
  auto& aIndexes = sliceA.pointIndexes;

  // Split along the longest edge at the median point
  int splitAxis = maxEdgeAxis(cloud, aIndexes);
  std::stable_sort(
    aIndexes.begin(), aIndexes.end(), [&](int32_t a, int32_t b) {
      return cloud[a][splitAxis] < cloud[b][splitAxis];
    });

  int numSplit = std::ceil((double)aIndexes.size() / (double)maxPoints);
  int splitsize = aIndexes.size() / numSplit;
  std::vector<Partition> splitPartitions;
  splitPartitions.resize(numSplit);

  // The 2nd to the penultimate partitions
  for (int i = 1; i < numSplit - 1; i++) {
    splitPartitions[i].sliceId = sliceA.sliceId + i;
    splitPartitions[i].tileId = sliceA.tileId;
    splitPartitions[i].origin = Vec3<int>{0};

    auto& Indexes = splitPartitions[i].pointIndexes;
    Indexes.insert(
      Indexes.begin(), aIndexes.begin() + i * splitsize,
      aIndexes.begin() + (i + 1) * splitsize);
  }
  // The last split partition
  auto& Indexes = splitPartitions[numSplit - 1].pointIndexes;
  Indexes.insert(
    Indexes.begin(), aIndexes.begin() + (numSplit - 1) * splitsize,
    aIndexes.end());

  aIndexes.erase(aIndexes.begin() + splitsize, aIndexes.end());

  toBeSplit = slices.insert(
    toBeSplit + 1, splitPartitions.begin() + 1, splitPartitions.end());

  return toBeSplit + (numSplit - 1);
}

//----------------------------------------------------------------------------
// combine the two slices into one

std::vector<Partition>::iterator
mergeSlice(
  std::vector<Partition>& slices,
  std::vector<Partition>::iterator a,
  std::vector<Partition>::iterator b)
{
  auto& aIndexes = (*a).pointIndexes;
  auto& bIndexes = (*b).pointIndexes;

  aIndexes.insert(aIndexes.end(), bIndexes.begin(), bIndexes.end());
  (*a).origin = minOrigin((*a).origin, (*b).origin);

  return slices.erase(b);
}

//=============================================================================
// first split slices still having too many points(more than maxPoints)
// then merge slices with too few points(less than minPoints)

void
refineSlices(
  const PartitionParams& params,
  const PCCPointSet3& cloud,
  std::vector<Partition>& slices)
{
  int maxPoints = params.sliceMaxPoints;
  int minPoints = params.sliceMinPoints;

  std::vector<Partition>::iterator it = slices.begin();
  while (it != slices.end()) {
    if ((*it).pointIndexes.size() > maxPoints) {
      it = splitSlice(cloud, slices, it, maxPoints);
    } else {
      it++;
    }
  }

  it = slices.begin();
  while (it != slices.end() && slices.size() > 1) {
    if ((*it).pointIndexes.size() < minPoints) {
      std::vector<Partition>::iterator toBeMerge;
      bool isFront = 0;

      // - a slice could only merge with the one before or after it
      // - the first/last slice could only merge with the next/front one
      // - let mergerfront = point number of slice after merging with the front one
      //   mergenext = point number of slice after merging with the subsquent one
      // - if both mergefront and mergenext < maxPoints, choose the larger one;
      // - if one of them < maxPoints and another > maxPoints,
      //     choose the small one to meet the requirement
      // - if both mergefront and mergenext > maxPoints, choose the larger one
      //     and do one more split after merge. We deem the slice after split
      //     as acceptable whether they are larger than minPoints or not
      // NB: if the merger slice is still smaller than minPoints,
      //     go on merging the same slice
      if (it == slices.begin()) {
        assert(it + 1 != slices.end());
        toBeMerge = it + 1;
        isFront = 0;
      } else if (it == slices.end() - 1) {
        assert(it - 1 != slices.end());
        toBeMerge = it - 1;
        isFront = 1;
      } else {
        int mergefront =
          (*it).pointIndexes.size() + (*(it - 1)).pointIndexes.size();
        int mergenext =
          (*it).pointIndexes.size() + (*(it + 1)).pointIndexes.size();

        if (
          (mergefront > maxPoints && mergenext > maxPoints)
          || (mergefront < maxPoints && mergenext < maxPoints)) {
          toBeMerge = mergefront > mergenext ? (it - 1) : (it + 1);
          isFront = mergefront > mergenext ? 1 : 0;
        } else {
          toBeMerge = mergefront < mergenext ? (it - 1) : (it + 1);
          isFront = mergefront < mergenext ? 1 : 0;
        }
      }

      it = isFront ? mergeSlice(slices, toBeMerge, it)
                   : mergeSlice(slices, it, toBeMerge);

      if ((*(it - 1)).pointIndexes.size() > maxPoints)
        it = splitSlice(cloud, slices, it - 1, maxPoints);
      else if ((*(it - 1)).pointIndexes.size() < minPoints)
        it--;

    } else {
      it++;
    }
  }

  for (int i = 0; i < slices.size(); i++) {
    auto& slice = slices[i];
    slice.sliceId = i;
    slice.tileId = -1;
  }
}

//============================================================================

}  // namespace pcc
