/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * <OWNER> = Apple Inc.
 * <ORGANIZATION> = Apple Inc.
 * <YEAR> = 2017
 *
 * Copyright (c) 2017, Apple Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PCCKdTree_h
#define PCCKdTree_h

#include <algorithm>

#include "PCCPointSet.h"

namespace pcc {
enum PCCAxis3 { PCC_AXIS3_UNDEFINED = -1, PCC_AXIS3_X = 0, PCC_AXIS3_Y = 1, PCC_AXIS3_Z = 2 };

struct PCCPointDistInfo {
  double dist2;
  uint32_t index;
  bool operator<(const PCCPointDistInfo &rhs) const { return dist2 < rhs.dist2; }
};

struct PCCNNResult {
  PCCPointDistInfo *neighbors;
  size_t resultCount;
};

struct PCCRangeResult {
  PCCPointDistInfo *neighbors;
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

template <typename T>
class PCCHeap {
 public:
  PCCHeap(void) = default;
  PCCHeap(const PCCHeap &) = delete;
  PCCHeap &operator=(const PCCHeap &) = delete;
  PCCHeap(T *array, size_t capacity) : array(array), cap(capacity), sz(0) {}
  T &top() {
    assert(sz > 0);
    return array[0];
  }
  const T &top() const {
    assert(sz > 0);
    return array[0];
  }
  void push(const T &elem) {
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
  void pop() {
    assert(sz > 0);
    array[0] = array[--sz];
    maxHeapify(0);
  }
  void clear() { sz = 0; }
  size_t size() const { return sz; }
  size_t capacity() const { return cap; }

 private:
  void maxHeapify(size_t i) {
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
  T *array;
  size_t cap;
  size_t sz;
};

class PCCStaticKdTree3 {
  struct PCCKdTree3Node {
    PCCPoint3D pos;
    uint32_t id;
    PCCAxis3 axis;
  };

 public:
  PCCStaticKdTree3(void) = default;
  PCCStaticKdTree3(const PCCStaticKdTree3 &) = default;
  PCCStaticKdTree3 &operator=(const PCCStaticKdTree3 &) = default;
  ~PCCStaticKdTree3(void) = default;
  size_t size() const { return nodes.size(); }
  size_t capacity() const { return nodes.capacity(); }
  void clear() { nodes.resize(0); }
  void reserve(const size_t pointCount) { nodes.reserve(pointCount); }
  void build(const PCCPointSet3 &pointCloud) {
    const size_t pointCount = pointCloud.getPointCount();
    init(pointCloud);
    if (pointCount) {
      build(0, uint32_t(pointCount));
    }
  }
  void findNearestNeighbors(const PCCNNQuery3 &query, PCCNNResult &result) const {
    assert(result.neighbors && query.nearestNeighborCount > 0);
    result.resultCount = 0;
    const PCCNNQuery3 query2 = {query.point, query.radius * query.radius,
                                query.nearestNeighborCount};
    PCCHeap<PCCPointDistInfo> neighbors(result.neighbors, query.nearestNeighborCount);
    findNearestNeighbors(0, uint32_t(nodes.size()), query2, neighbors);
    std::sort(result.neighbors, result.neighbors + neighbors.size());
    result.resultCount = neighbors.size();
  }
  void findNeighbors(const PCCRangeQuery3 &query, PCCRangeResult &result) const {
    assert(result.neighbors && query.maxResultCount > 0);
    result.resultCount = 0;
    const PCCRangeQuery3 query2 = {query.point, query.radius * query.radius, query.maxResultCount};
    findNeighbors(0, uint32_t(nodes.size()), query2, result);
  }

 private:
  void init(const PCCPointSet3 &pointCloud) {
    const size_t size = pointCloud.getPointCount();
    nodes.resize(size);
    for (size_t i = 0; i < size; ++i) {
      nodes[i].id = uint32_t(i);
      nodes[i].axis = PCC_AXIS3_UNDEFINED;
      nodes[i].pos = pointCloud[i];
    }
  }
  void build(const uint32_t start, const uint32_t end) {
    assert(start < end);
    uint32_t index = start;
    if (end == start + 1) {
      nodes[start].axis = PCC_AXIS3_X;
    } else {
      PCCAxis3 splitAxis = computeSplitAxis(start, end);
      index = findMedian(start, end, splitAxis);
      PCCKdTree3Node &current = nodes[index];
      current.axis = splitAxis;
      if (start < index) {
        build(start, index);
      }
      if (index + 1 < end) {
        build(index + 1, end);
      }
    }
  }
  PCCAxis3 computeSplitAxis(const uint32_t start, const uint32_t end) {
    PCCPoint3D minBB = nodes[start].pos;
    PCCPoint3D maxBB = nodes[start].pos;
    for (size_t i = start + 1; i < end; ++i) {
      const PCCPoint3D &pt = nodes[i].pos;
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
  uint32_t findMedian(uint32_t start, uint32_t end, const PCCAxis3 splitAxis) {
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

      while (store < medianIndex &&
             PCCApproximatelyEqual(
                 nodes[store].pos[splitAxis],
                 nodes[store + 1].pos[splitAxis])) {  // optimization in case of duplicated values
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
  void findNearestNeighbors(const uint32_t start, const uint32_t end, const PCCNNQuery3 &query2,
                            PCCHeap<PCCPointDistInfo> &neighbors) const {
    if (start >= end) {
      return;
    } else if (start + 1 == end) {
      const PCCKdTree3Node &medianNode = nodes[start];
      const double dist2 = (query2.point - medianNode.pos).getNorm2();

      if (dist2 < query2.radius) {
        PCCPointDistInfo distInfo = {dist2, medianNode.id};
        neighbors.push(distInfo);
      }
    } else {
      const uint32_t medianIndex = start + (end - start) / 2;
      const PCCKdTree3Node &medianNode = nodes[medianIndex];
      const double dist2 = (query2.point - medianNode.pos).getNorm2();
      const PCCAxis3 splitAxis = medianNode.axis;
      const double coord = medianNode.pos[splitAxis];
      const double dx = query2.point[splitAxis] - coord;
      const double dx2 = dx * dx;

      if (dx < 0.0) {
        findNearestNeighbors(start, medianIndex, query2, neighbors);
        if (dist2 < query2.radius) {
          PCCPointDistInfo distInfo = {dist2, medianNode.id};
          neighbors.push(distInfo);
        }
        if (dx2 > query2.radius ||
            (neighbors.size() == query2.nearestNeighborCount && dx2 > neighbors.top().dist2)) {
          return;
        }
        findNearestNeighbors(medianIndex + 1, end, query2, neighbors);
      } else {
        findNearestNeighbors(medianIndex + 1, end, query2, neighbors);
        if (dist2 < query2.radius) {
          PCCPointDistInfo distInfo = {dist2, medianNode.id};
          neighbors.push(distInfo);
        }
        if (dx2 > query2.radius ||
            (neighbors.size() == query2.nearestNeighborCount && dx2 > neighbors.top().dist2)) {
          return;
        }
        findNearestNeighbors(start, medianIndex, query2, neighbors);
      }
    }
  }
  void findNeighbors(const uint32_t start, const uint32_t end, const PCCRangeQuery3 &query2,
                     PCCRangeResult &result) const {
    if (start >= end || result.resultCount == query2.maxResultCount) {
      return;
    } else if (start + 1 == end) {
      const PCCKdTree3Node &medianNode = nodes[start];
      const double dist2 = (query2.point - medianNode.pos).getNorm2();

      if (dist2 < query2.radius) {
        result.neighbors[result.resultCount].dist2 = dist2;
        result.neighbors[result.resultCount].index = medianNode.id;
        result.resultCount++;
      }
    } else {
      const uint32_t medianIndex = start + (end - start) / 2;
      const PCCKdTree3Node &medianNode = nodes[medianIndex];
      const double dist2 = (query2.point - medianNode.pos).getNorm2();
      const PCCAxis3 splitAxis = medianNode.axis;
      const double coord = medianNode.pos[splitAxis];
      const double dx = query2.point[splitAxis] - coord;
      const double dx2 = dx * dx;

      if (dx < 0.0) {
        findNeighbors(start, medianIndex, query2, result);
        if (result.resultCount == query2.maxResultCount) {
          return;
        }
        if (dist2 < query2.radius) {
          result.neighbors[result.resultCount].dist2 = dist2;
          result.neighbors[result.resultCount].index = medianNode.id;
          result.resultCount++;
        }
        if (dx2 > query2.radius || result.resultCount == query2.maxResultCount) {
          return;
        }
        findNeighbors(medianIndex + 1, end, query2, result);
      } else {
        findNeighbors(medianIndex + 1, end, query2, result);
        if (result.resultCount == query2.maxResultCount) {
          return;
        }
        if (dist2 < query2.radius) {
          result.neighbors[result.resultCount].dist2 = dist2;
          result.neighbors[result.resultCount].index = medianNode.id;
          result.resultCount++;
        }
        if (dx2 > query2.radius || result.resultCount == query2.maxResultCount) {
          return;
        }
        findNeighbors(start, medianIndex, query2, result);
      }
    }
  }

 private:
  std::vector<PCCKdTree3Node> nodes;
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
  PCCIncrementalKdTree3(const PCCIncrementalKdTree3 &) = default;
  PCCIncrementalKdTree3 &operator=(const PCCIncrementalKdTree3 &) = default;
  ~PCCIncrementalKdTree3(void) = default;

  void build(const PCCPointSet3 &pointCloud) {
    clear();
    insert(pointCloud);
  }

  size_t size() const { return nodes.size(); }
  size_t capacity() const { return nodes.capacity(); }
  void clear() {
    nodes.resize(0);
    root = PCC_UNDEFINED_INDEX;
  }
  void reserve(const size_t pointCount) { nodes.reserve(pointCount); }
  uint32_t insert(const PCCPoint3D point) {
    const uint32_t id = static_cast<uint32_t>(nodes.size());
    nodes.resize(id + 1);
    PCCIncrementalKdTree3Node &node = nodes[id];
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
  void insert(const PCCPointSet3 &pointCloud) {
    if (pointCloud.getPointCount()) {
      const uint32_t start = static_cast<uint32_t>(nodes.size());
      const uint32_t end = start + static_cast<uint32_t>(pointCloud.getPointCount());
      nodes.resize(end);
      for (uint32_t index = start; index < end; ++index) {
        PCCIncrementalKdTree3Node &node = nodes[index];
        node.pos = pointCloud[index - start];
        node.id = index;
      }
      root = balance(0, uint32_t(nodes.size()));
    }
  }
  void balance() {
    if (!nodes.empty()) {
      root = balance(0, uint32_t(nodes.size()));
    }
  }
  void findNearestNeighbors(const PCCNNQuery3 &query, PCCNNResult &result) const {
    assert(result.neighbors && query.nearestNeighborCount > 0);
    result.resultCount = 0;
    const PCCNNQuery3 query2 = {query.point, query.radius * query.radius,
                                query.nearestNeighborCount};
    PCCHeap<PCCPointDistInfo> neighbors(result.neighbors, query.nearestNeighborCount);
    findNearestNeighbors(root, query2, neighbors);
    std::sort(result.neighbors, result.neighbors + neighbors.size());
    result.resultCount = neighbors.size();
  }
  void findNearestNeighbors2(const PCCNNQuery3 &query2, PCCNNResult &result) const {
    assert(result.neighbors && query2.nearestNeighborCount > 0);
    result.resultCount = 0;
    PCCHeap<PCCPointDistInfo> neighbors(result.neighbors, query2.nearestNeighborCount);
    findNearestNeighbors(root, query2, neighbors);
    std::sort(result.neighbors, result.neighbors + neighbors.size());
    result.resultCount = neighbors.size();
  }
  void findNeighbors(const PCCRangeQuery3 &query, PCCRangeResult &result) const {
    assert(result.neighbors && query.maxResultCount > 0);
    result.resultCount = 0;
    const PCCRangeQuery3 query2 = {query.point, query.radius * query.radius, query.maxResultCount};
    findNeighbors(root, query2, result);
  }
  void append(const PCCPoint3D point) {
    const uint32_t id = static_cast<uint32_t>(nodes.size());
    nodes.resize(id + 1);
    PCCIncrementalKdTree3Node &node = nodes[id];
    node.pos = point;
    node.id = id;
    node.left = PCC_UNDEFINED_INDEX;
    node.right = PCC_UNDEFINED_INDEX;
  }
  void append(const PCCPointSet3 &pointCloud) {
    if (pointCloud.getPointCount()) {
      const uint32_t start = static_cast<uint32_t>(nodes.size());
      const uint32_t end = start + static_cast<uint32_t>(pointCloud.getPointCount());
      nodes.resize(end);
      for (uint32_t index = start; index < end; ++index) {
        PCCIncrementalKdTree3Node &node = nodes[index];
        node.pos = pointCloud[index - start];
        node.id = index;
        node.left = PCC_UNDEFINED_INDEX;
        node.right = PCC_UNDEFINED_INDEX;
      }
    }
  }

 private:
  static PCCAxis3 nextAxis(const PCCAxis3 axis) {
    switch (axis) {
      case PCC_AXIS3_X:
        return PCC_AXIS3_Y;
      case PCC_AXIS3_Y:
        return PCC_AXIS3_Z;
      case PCC_AXIS3_Z:
        return PCC_AXIS3_X;
      case PCC_AXIS3_UNDEFINED:
        return PCC_AXIS3_X;
      default:
        return PCC_AXIS3_X;
    }
  }
  void insert(PCCIncrementalKdTree3Node &node, const uint32_t parent) {
    const PCCAxis3 axis = nodes[parent].axis;
    uint32_t &index =
        (node.pos[axis] < nodes[parent].pos[axis]) ? nodes[parent].left : nodes[parent].right;
    if (index == PCC_UNDEFINED_INDEX) {
      index = node.id;
      node.axis = nextAxis(axis);
    } else {
      insert(node, index);
    }
  }
  PCCAxis3 computeSplitAxis(const uint32_t start, const uint32_t end) const {
    PCCPoint3D minBB = nodes[start].pos;
    PCCPoint3D maxBB = nodes[start].pos;
    for (size_t i = start + 1; i < end; ++i) {
      const PCCPoint3D &pt = nodes[i].pos;
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
  uint32_t findMedian(uint32_t start, uint32_t end, const PCCAxis3 splitAxis) {
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

      while (store < medianIndex &&
             PCCApproximatelyEqual(
                 nodes[store].pos[splitAxis],
                 nodes[store + 1].pos[splitAxis])) {  // optimization in case of duplicated values
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
  uint32_t balance(const uint32_t start, const uint32_t end) {
    assert(start < end);
    uint32_t index = start;
    PCCAxis3 splitAxis = computeSplitAxis(start, end);
    index = findMedian(start, end, splitAxis);
    PCCIncrementalKdTree3Node &node = nodes[index];
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
  void findNearestNeighbors(const uint32_t current, const PCCNNQuery3 &query2,
                            PCCHeap<PCCPointDistInfo> &neighbors) const {
    if (current == PCC_UNDEFINED_INDEX) {
      return;
    }
    const PCCIncrementalKdTree3Node &node = nodes[current];
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
    if (dx2 > query2.radius ||
        (neighbors.size() == query2.nearestNeighborCount && dx2 > neighbors.top().dist2)) {
      return;
    }
    findNearestNeighbors(second, query2, neighbors);
  }
  void findNeighbors(const uint32_t current, const PCCRangeQuery3 &query2,
                     PCCRangeResult &result) const {
    if (current == PCC_UNDEFINED_INDEX) {
      return;
    }
    const PCCIncrementalKdTree3Node &node = nodes[current];
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
    if ((dx * dx) > query2.radius || result.resultCount == query2.maxResultCount) {
      return;
    }
    findNeighbors(second, query2, result);
  }

 private:
  std::vector<PCCIncrementalKdTree3Node> nodes;
  uint32_t root;
};
}
#endif /* PCCKdTree_h */
