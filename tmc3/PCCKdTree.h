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

#ifndef PCCKdTree_h
#define PCCKdTree_h

#include <algorithm>
#include <cstddef>

#include "PCCPointSet.h"

namespace pcc {
enum PCCAxis3
{
  PCC_AXIS3_UNDEFINED = -1,
  PCC_AXIS3_X = 0,
  PCC_AXIS3_Y = 1,
  PCC_AXIS3_Z = 2
};

struct PCCPointDistInfo {
  double dist2;
  uint32_t index;
  bool operator<(const PCCPointDistInfo& rhs) const
  {
    return dist2 < rhs.dist2;
  }
};

struct PCCNNResult {
  PCCPointDistInfo* neighbors;
  size_t resultCount;
};

struct PCCRangeResult {
  PCCPointDistInfo* neighbors;
  size_t resultCount;
};

struct PCCNNQuery3 {
  PCCPoint3D point;
  double radius;
  size_t nearestNeighborCount;
};

struct PCCRangeQuery3 {
  PCCPoint3D point;
  double radius;
  size_t maxResultCount;
};

template<typename T>
class PCCHeap {
public:
  PCCHeap(void) = default;
  PCCHeap(const PCCHeap&) = delete;
  PCCHeap& operator=(const PCCHeap&) = delete;
  PCCHeap(T* array, size_t capacity) : array(array), cap(capacity), sz(0) {}
  T& top()
  {
    assert(sz > 0);
    return array[0];
  }
  const T& top() const
  {
    assert(sz > 0);
    return array[0];
  }
  void push(const T& elem)
  {
    if (sz < cap) {
      size_t i = sz;
      array[i] = elem;
      size_t p = parent(i);
      while (i > 0 && array[p] < elem) {
        std::swap(array[i], array[p]);
        i = p;
        p = parent(i);
      }
      ++sz;
    } else if (elem < array[0]) {
      array[0] = elem;
      maxHeapify(0);
    }
  }
  void pop()
  {
    assert(sz > 0);
    array[0] = array[--sz];
    maxHeapify(0);
  }
  void clear() { sz = 0; }
  size_t size() const { return sz; }
  size_t capacity() const { return cap; }

private:
  void maxHeapify(size_t i)
  {
    size_t l = left(i);
    size_t r = right(i);
    size_t g = (l < sz && array[i] < array[l]) ? l : i;
    if (r < sz && array[g] < array[r]) {
      g = r;
    }
    if (g != i) {
      std::swap(array[i], array[g]);
      maxHeapify(g);
    }
  }
  size_t parent(size_t index) const { return (index + (index & 1) - 1) / 2; }
  size_t left(size_t index) const { return 2 * index + 1; }
  size_t right(size_t index) const { return 2 * index + 2; }

private:
  T* array;
  size_t cap;
  size_t sz;
};

class PCCIncrementalKdTree3 {
  struct PCCIncrementalKdTree3Node {
    PCCPoint3D pos;
    uint32_t id;
    uint32_t left;
    uint32_t right;
    PCCAxis3 axis;
  };

public:
  PCCIncrementalKdTree3(void) { root = PCC_UNDEFINED_INDEX; }
  PCCIncrementalKdTree3(const PCCIncrementalKdTree3&) = default;
  PCCIncrementalKdTree3& operator=(const PCCIncrementalKdTree3&) = default;
  ~PCCIncrementalKdTree3(void) = default;

  void build(const PCCPointSet3& pointCloud)
  {
    clear();
    insert(pointCloud);
  }

  size_t size() const { return nodes.size(); }
  size_t capacity() const { return nodes.capacity(); }
  void clear()
  {
    nodes.resize(0);
    root = PCC_UNDEFINED_INDEX;
  }
  void reserve(const size_t pointCount) { nodes.reserve(pointCount); }
  uint32_t insert(const PCCPoint3D point)
  {
    const uint32_t id = static_cast<uint32_t>(nodes.size());
    nodes.resize(id + 1);
    PCCIncrementalKdTree3Node& node = nodes[id];
    node.pos = point;
    node.id = id;
    node.left = PCC_UNDEFINED_INDEX;
    node.right = PCC_UNDEFINED_INDEX;
    if (root == PCC_UNDEFINED_INDEX) {
      root = id;
      node.axis = PCC_AXIS3_X;
    } else {
      insert(node, root);
    }
    return id;
  }
  void insert(const PCCPointSet3& pointCloud)
  {
    if (pointCloud.getPointCount()) {
      const uint32_t start = static_cast<uint32_t>(nodes.size());
      const uint32_t end =
        start + static_cast<uint32_t>(pointCloud.getPointCount());
      nodes.resize(end);
      for (uint32_t index = start; index < end; ++index) {
        PCCIncrementalKdTree3Node& node = nodes[index];
        node.pos = pointCloud[index - start];
        node.id = index;
      }
      root = balance(0, uint32_t(nodes.size()));
    }
  }
  void balance()
  {
    if (!nodes.empty()) {
      root = balance(0, uint32_t(nodes.size()));
    }
  }
  void
  findNearestNeighbors(const PCCNNQuery3& query, PCCNNResult& result) const
  {
    assert(result.neighbors && query.nearestNeighborCount > 0);
    result.resultCount = 0;
    const PCCNNQuery3 query2 = {query.point, query.radius * query.radius,
                                query.nearestNeighborCount};
    PCCHeap<PCCPointDistInfo> neighbors(
      result.neighbors, query.nearestNeighborCount);
    findNearestNeighbors(root, query2, neighbors);
    std::sort(result.neighbors, result.neighbors + neighbors.size());
    result.resultCount = neighbors.size();
  }
  void
  findNearestNeighbors2(const PCCNNQuery3& query2, PCCNNResult& result) const
  {
    assert(result.neighbors && query2.nearestNeighborCount > 0);
    result.resultCount = 0;
    PCCHeap<PCCPointDistInfo> neighbors(
      result.neighbors, query2.nearestNeighborCount);
    findNearestNeighbors(root, query2, neighbors);
    std::sort(result.neighbors, result.neighbors + neighbors.size());
    result.resultCount = neighbors.size();
  }
  void findNeighbors(const PCCRangeQuery3& query, PCCRangeResult& result) const
  {
    assert(result.neighbors && query.maxResultCount > 0);
    result.resultCount = 0;
    const PCCRangeQuery3 query2 = {query.point, query.radius * query.radius,
                                   query.maxResultCount};
    findNeighbors(root, query2, result);
  }

  uint32_t
  hasNeighborWithinRange(const PCCVector3D& point, const double radius2) const
  {
    return hasNeighborWithinRange(root, point, radius2);
  }

private:
  static PCCAxis3 nextAxis(const PCCAxis3 axis)
  {
    switch (axis) {
    case PCC_AXIS3_X: return PCC_AXIS3_Y;
    case PCC_AXIS3_Y: return PCC_AXIS3_Z;
    case PCC_AXIS3_Z: return PCC_AXIS3_X;
    case PCC_AXIS3_UNDEFINED: return PCC_AXIS3_X;
    default: return PCC_AXIS3_X;
    }
  }
  void insert(PCCIncrementalKdTree3Node& node, const uint32_t parent)
  {
    const PCCAxis3 axis = nodes[parent].axis;
    uint32_t& index = (node.pos[axis] < nodes[parent].pos[axis])
      ? nodes[parent].left
      : nodes[parent].right;
    if (index == PCC_UNDEFINED_INDEX) {
      index = node.id;
      node.axis = nextAxis(axis);
    } else {
      insert(node, index);
    }
  }
  PCCAxis3 computeSplitAxis(const uint32_t start, const uint32_t end) const
  {
    PCCPoint3D minBB = nodes[start].pos;
    PCCPoint3D maxBB = nodes[start].pos;
    for (size_t i = start + 1; i < end; ++i) {
      const PCCPoint3D& pt = nodes[i].pos;
      for (int32_t k = 0; k < 3; ++k) {
        if (minBB[k] > pt[k]) {
          minBB[k] = pt[k];
        } else if (maxBB[k] < pt[k]) {
          maxBB[k] = pt[k];
        }
      }
    }
    PCCPoint3D d = maxBB - minBB;
    if (d.x() > d.y() && d.x() > d.z()) {
      return PCC_AXIS3_X;
    } else if (d.y() > d.z()) {
      return PCC_AXIS3_Y;
    } else {
      return PCC_AXIS3_Z;
    }
  }
  uint32_t findMedian(uint32_t start, uint32_t end, const PCCAxis3 splitAxis)
  {
    assert(start < end);
    if (end == start + 1) {
      return start;
    }
    const uint32_t medianIndex = start + (end - start) / 2;
    while (1) {
      double pivot = nodes[medianIndex].pos[splitAxis];
      std::swap(nodes[medianIndex], nodes[end - 1]);
      uint32_t store, p;
      for (store = p = start; p < end; p++) {
        if (nodes[p].pos[splitAxis] < pivot) {
          if (p != store) {
            std::swap(nodes[p], nodes[store]);
          }
          ++store;
        }
      }
      std::swap(nodes[store], nodes[end - 1]);

      while (store < medianIndex
             && nodes[store].pos[splitAxis]
               == nodes[store + 1].pos[splitAxis]) {
        // optimization in case of duplicated values
        ++store;
      }

      if (store == medianIndex) {
        return medianIndex;
      } else if (store > medianIndex) {
        end = store;
      } else {
        start = store + 1;
      }
    }
  }
  uint32_t balance(const uint32_t start, const uint32_t end)
  {
    assert(start < end);
    uint32_t index = start;
    PCCAxis3 splitAxis = computeSplitAxis(start, end);
    index = findMedian(start, end, splitAxis);
    PCCIncrementalKdTree3Node& node = nodes[index];
    node.axis = splitAxis;
    if (start < index) {
      node.left = balance(start, index);
    } else {
      node.left = PCC_UNDEFINED_INDEX;
    }
    if (index + 1 < end) {
      node.right = balance(index + 1, end);
    } else {
      node.right = PCC_UNDEFINED_INDEX;
    }
    return index;
  }
  void findNearestNeighbors(
    const uint32_t current,
    const PCCNNQuery3& query2,
    PCCHeap<PCCPointDistInfo>& neighbors) const
  {
    if (current == PCC_UNDEFINED_INDEX) {
      return;
    }
    const PCCIncrementalKdTree3Node& node = nodes[current];
    const double dist2 = (query2.point - node.pos).getNorm2();
    const PCCAxis3 splitAxis = node.axis;
    const double coord = node.pos[splitAxis];
    const double dx = query2.point[splitAxis] - coord;
    uint32_t first, second;
    if (dx < 0.0) {
      first = node.left;
      second = node.right;
    } else {
      first = node.right;
      second = node.left;
    }
    findNearestNeighbors(first, query2, neighbors);
    if (dist2 < query2.radius) {
      PCCPointDistInfo distInfo = {dist2, node.id};
      neighbors.push(distInfo);
    }
    const double dx2 = dx * dx;
    if (
      dx2 > query2.radius
      || (neighbors.size() == query2.nearestNeighborCount
          && dx2 > neighbors.top().dist2)) {
      return;
    }
    findNearestNeighbors(second, query2, neighbors);
  }
  void findNeighbors(
    const uint32_t current,
    const PCCRangeQuery3& query2,
    PCCRangeResult& result) const
  {
    if (current == PCC_UNDEFINED_INDEX) {
      return;
    }
    const PCCIncrementalKdTree3Node& node = nodes[current];
    const double dist2 = (query2.point - node.pos).getNorm2();
    const PCCAxis3 splitAxis = node.axis;
    const double coord = node.pos[splitAxis];
    const double dx = query2.point[splitAxis] - coord;
    uint32_t first, second;
    if (dx < 0.0) {
      first = node.left;
      second = node.right;
    } else {
      first = node.right;
      second = node.left;
    }
    findNeighbors(first, query2, result);
    if (result.resultCount == query2.maxResultCount) {
      return;
    }
    if (dist2 < query2.radius) {
      result.neighbors[result.resultCount].dist2 = dist2;
      result.neighbors[result.resultCount].index = node.id;
      result.resultCount++;
    }
    if (
      (dx * dx) > query2.radius
      || result.resultCount == query2.maxResultCount) {
      return;
    }
    findNeighbors(second, query2, result);
  }

  uint32_t hasNeighborWithinRange(
    const uint32_t current,
    const PCCVector3D& point,
    const double radius2) const
  {
    if (current == PCC_UNDEFINED_INDEX) {
      return PCC_UNDEFINED_INDEX;
    }
    const PCCIncrementalKdTree3Node& node = nodes[current];
    const double dist2 = (point - node.pos).getNorm2();
    if (dist2 < radius2) {
      return node.id;
    }
    const PCCAxis3 splitAxis = node.axis;
    const int32_t coord = node.pos[splitAxis];
    const int32_t dx = point[splitAxis] - coord;
    uint32_t first, second;
    if (dx < 0) {
      first = node.left;
      second = node.right;
    } else {
      first = node.right;
      second = node.left;
    }
    const uint32_t index = hasNeighborWithinRange(first, point, radius2);
    if (index != PCC_UNDEFINED_INDEX) {
      return index;
    }
    if ((dx * dx) > radius2) {
      return PCC_UNDEFINED_INDEX;
    }
    return hasNeighborWithinRange(second, point, radius2);
  }

private:
  std::vector<PCCIncrementalKdTree3Node> nodes;
  uint32_t root;
};

//---------------------------------------
class PCCKdTree3 {
  struct PCCKdTree3Node {
    PCCBox3D BB;
    PCCPoint3D centd;
    uint32_t id;
    uint32_t start;
    uint32_t end;
    PCCAxis3 axis;
    uint32_t median;
    uint32_t medianIdx;
  };
  struct PointIDNode {
    PCCPoint3D pos;
    uint32_t id;
    bool isVisisted;
  };

public:
  PCCKdTree3(const PCCPointSet3& pointCloud, uint8_t depth)
  {
    init(pointCloud, depth);
  }
  PCCKdTree3(const PCCKdTree3&) = default;
  PCCKdTree3& operator=(const PCCKdTree3&) = default;
  ~PCCKdTree3(void) = default;

  std::vector<PCCKdTree3Node> nodes;
  void build()
  {
    uint32_t nodeCount =
      (1 << (kdDepth + 1)) - 1;  // std::pow(2, kdDepth + 1) - 1
    for (size_t nodeIt = 1; nodeIt < nodeCount; nodeIt++) {
      bool isLeftNode = (nodeIt & 1) != 0;
      uint32_t parentNodeIdx =
        isLeftNode ? (nodeIt - 1) >> 1 : (nodeIt - 2) >> 1;
      const uint32_t start = isLeftNode
        ? static_cast<uint32_t>(nodes[parentNodeIdx].start)
        : static_cast<uint32_t>(nodes[parentNodeIdx].medianIdx + 1);
      const uint32_t end = isLeftNode
        ? static_cast<uint32_t>(nodes[parentNodeIdx].medianIdx)
        : static_cast<uint32_t>(nodes[parentNodeIdx].end);
      PCCBox3D BB = nodes[parentNodeIdx].BB;
      isLeftNode
        ? BB.max[nodes[parentNodeIdx].axis] = nodes[parentNodeIdx].median
        : BB.min[nodes[parentNodeIdx].axis] = nodes[parentNodeIdx].median;

      PCCPoint3D nodeMean = computePCCMean(start, end);
      PCCAxis3 axis = computeSplitAxisVar(start, end, nodeMean);
      uint32_t medianIdx = findMedian(start, end, axis);

      PCCKdTree3Node& node = nodes[nodeIt];
      node.BB = BB;
      node.centd = nodeMean;
      node.id = nodeIt;
      node.start = start;
      node.end = end;
      node.axis = axis;
      node.median = pointCloudTemp[medianIdx].pos[axis];
      node.medianIdx = medianIdx;
    }
  }

  void init(const PCCPointSet3& pointCloud, uint8_t depth)
  {
    kdDepth = depth;
    nodes.resize(0);
    nodes.resize(std::pow(2, kdDepth + 1) - 1);
    uint32_t pointCount = pointCloud.getPointCount();
    pointCloudTemp.resize(pointCount);
    for (size_t i = 0; i < pointCount; i++) {
      pointCloudTemp[i].pos = pointCloud[i];
      pointCloudTemp[i].id = i;
      pointCloudTemp[i].isVisisted = false;
    }

    PCCBox3D BB = computeBoundingBox(0, pointCount);
    PCCPoint3D nodeMean = computePCCMean(0, pointCount);
    PCCAxis3 axis = computeSplitAxisVar(0, pointCount, nodeMean);
    uint32_t medianIdx = findMedian(0, pointCount, axis);

    PCCKdTree3Node& rootNode = nodes[0];
    rootNode.BB = BB;
    rootNode.centd = nodeMean;
    rootNode.id = 0;
    rootNode.start = static_cast<uint32_t>(0);
    rootNode.end = static_cast<uint32_t>(pointCount - 1);
    ;
    rootNode.axis = axis;
    rootNode.median = pointCloudTemp[medianIdx].pos[axis];
    rootNode.medianIdx = medianIdx;
  }

  uint32_t searchClosestAvailablePoint(PCCPoint3D queryPoint)
  {
    uint32_t idToClosestPoint = -1;
    uint32_t id = 0;
    for (int8_t d = 0; d < kdDepth; d++) {
      id = (queryPoint[nodes[id].axis] <= nodes[id].median) ? 2 * id + 1
                                                            : 2 * id + 2;
    }
    const uint32_t start = nodes[id].start;
    const uint32_t end = nodes[id].end;
    const uint32_t closestDistThr = 1;
    uint32_t smallestDist = -1;
    uint32_t closestID = 0;
    for (size_t i = start; i <= end; ++i) {
      if (!pointCloudTemp[i].isVisisted) {
        PCCPoint3D diff = pointCloudTemp[i].pos - queryPoint;
        uint32_t dist =
          std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
        if (dist <= closestDistThr) {
          return pointCloudTemp[i].id;
        }
        if (dist < smallestDist) {
          smallestDist = dist;
          idToClosestPoint = pointCloudTemp[i].id;
          closestID = i;
        }
      }
    }
    if (smallestDist != -1) {
      pointCloudTemp[closestID].isVisisted = true;
    }
    return idToClosestPoint;
  }

private:
  PCCBox3D computeBoundingBox(const uint32_t start, const uint32_t end) const
  {
    PCCPoint3D minBB = pointCloudTemp[start].pos;
    PCCPoint3D maxBB = pointCloudTemp[start].pos;
    for (size_t i = start + 1; i < end; ++i) {
      const PCCPoint3D& pt = pointCloudTemp[i].pos;
      for (int32_t k = 0; k < 3; ++k) {
        if (minBB[k] > pt[k]) {
          minBB[k] = pt[k];
        } else if (maxBB[k] < pt[k]) {
          maxBB[k] = pt[k];
        }
      }
    }
    PCCBox3D BB;
    {
      BB.min = minBB;
      BB.max = maxBB;
    }
    return BB;
  }
  PCCAxis3 computeSplitAxis(const uint32_t start, const uint32_t end) const
  {
    PCCBox3D BB = computeBoundingBox(start, end);
    PCCPoint3D d = BB.max - BB.min;
    if (d.x() > d.y() && d.x() > d.z()) {
      return PCC_AXIS3_X;
    } else if (d.y() > d.z()) {
      return PCC_AXIS3_Y;
    } else {
      return PCC_AXIS3_Z;
    }
  }
  PCCAxis3 computeSplitAxisVar(
    const uint32_t start, const uint32_t end, PCCPoint3D nodeMean) const
  {
    double nodeVar[3] = {0, 0, 0};
    for (size_t axis = 0; axis < 3; ++axis) {
      double acc = 0, diff = 0;
      for (size_t i = start; i < end; i++) {
        diff = pointCloudTemp[i].pos[axis] - nodeMean[axis];
        acc += diff * diff;
      }
      nodeVar[axis] = acc / (end - start);
    }
    if (nodeVar[0] > nodeVar[1] && nodeVar[0] > nodeVar[2]) {
      return PCC_AXIS3_X;
    } else if (nodeVar[1] > nodeVar[2]) {
      return PCC_AXIS3_Y;
    } else {
      return PCC_AXIS3_Z;
    }
  }
  PCCPoint3D computePCCMean(const uint32_t start, const uint32_t end)
  {
    assert(end >= start);
    PCCPoint3D nodeMean;
    for (size_t axis = 0; axis < 3; ++axis) {
      uint32_t acc = 0;
      for (size_t i = start; i < end; i++) {
        acc += pointCloudTemp[i].pos[axis];
      }
      nodeMean[axis] = round(acc / (end - start));
    }
    return nodeMean;
  }

  uint32_t findMedian(uint32_t start, uint32_t end, const PCCAxis3 splitAxis)
  {
    assert(start < end);
    if (end == start + 1) {
      return start;
    }
    const uint32_t medianIndex = start + (end - start) / 2;
    while (1) {
      double pivot = pointCloudTemp[medianIndex].pos[splitAxis];
      std::swap(pointCloudTemp[medianIndex], pointCloudTemp[end - 1]);
      uint32_t store, p;
      for (store = p = start; p < end; p++) {
        if (pointCloudTemp[p].pos[splitAxis] < pivot) {
          if (p != store) {
            std::swap(pointCloudTemp[p], pointCloudTemp[store]);
          }
          ++store;
        }
      }
      std::swap(pointCloudTemp[store], pointCloudTemp[end - 1]);

      while (store < medianIndex
             && pointCloudTemp[store].pos[splitAxis]
               == pointCloudTemp[store + 1].pos[splitAxis]) {
        // optimization in case of duplicated values
        ++store;
      }

      if (store == medianIndex) {
        return medianIndex;
      } else if (store > medianIndex) {
        end = store;
      } else {
        start = store + 1;
      }
    }
  }

private:
  std::vector<PointIDNode> pointCloudTemp;
  uint8_t kdDepth;
};
//---------------------------------------
}  // namespace pcc
#endif /* PCCKdTree_h */
