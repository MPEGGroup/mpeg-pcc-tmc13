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

//---------------------------------------
class PCCKdTree3 {
  struct PCCKdTree3Node {
    Box3<double> BB;
    Vec3<double> centd;
    uint32_t id;
    uint32_t start;
    uint32_t end;
    PCCAxis3 axis;
    uint32_t median;
    uint32_t medianIdx;
  };
  struct PointIDNode {
    Vec3<double> pos;
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
      Box3<double> BB = nodes[parentNodeIdx].BB;
      isLeftNode
        ? BB.max[nodes[parentNodeIdx].axis] = nodes[parentNodeIdx].median
        : BB.min[nodes[parentNodeIdx].axis] = nodes[parentNodeIdx].median;

      Vec3<double> nodeMean = computePCCMean(start, end);
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

    Box3<double> BB = computeBoundingBox(0, pointCount);
    Vec3<double> nodeMean = computePCCMean(0, pointCount);
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

  uint32_t searchClosestAvailablePoint(Vec3<double> queryPoint)
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
        Vec3<double> diff = pointCloudTemp[i].pos - queryPoint;
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
  Box3<double>
  computeBoundingBox(const uint32_t start, const uint32_t end) const
  {
    Vec3<double> minBB = pointCloudTemp[start].pos;
    Vec3<double> maxBB = pointCloudTemp[start].pos;
    for (size_t i = start + 1; i < end; ++i) {
      const Vec3<double>& pt = pointCloudTemp[i].pos;
      for (int32_t k = 0; k < 3; ++k) {
        if (minBB[k] > pt[k]) {
          minBB[k] = pt[k];
        } else if (maxBB[k] < pt[k]) {
          maxBB[k] = pt[k];
        }
      }
    }
    Box3<double> BB;
    {
      BB.min = minBB;
      BB.max = maxBB;
    }
    return BB;
  }
  PCCAxis3 computeSplitAxis(const uint32_t start, const uint32_t end) const
  {
    Box3<double> BB = computeBoundingBox(start, end);
    Vec3<double> d = BB.max - BB.min;
    if (d.x() > d.y() && d.x() > d.z()) {
      return PCC_AXIS3_X;
    } else if (d.y() > d.z()) {
      return PCC_AXIS3_Y;
    } else {
      return PCC_AXIS3_Z;
    }
  }
  PCCAxis3 computeSplitAxisVar(
    const uint32_t start, const uint32_t end, Vec3<double> nodeMean) const
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
  Vec3<double> computePCCMean(const uint32_t start, const uint32_t end)
  {
    assert(end >= start);
    Vec3<double> nodeMean;
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
