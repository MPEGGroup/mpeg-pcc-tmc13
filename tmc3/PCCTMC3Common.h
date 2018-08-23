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

#ifndef PCCTMC3Common_h
#define PCCTMC3Common_h

#include "ArithmeticCodec.h"
#include "PCCKdTree.h"
#include "PCCMath.h"
#include "ringbuf.h"

#include "nanoflann.hpp"

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pcc {
const uint32_t PCCTMC3MaxPredictionNearestNeighborCount = 3;
const uint32_t PCCTMC3SymbolCount = 64;

const int MAX_NUM_DM_LEAF_POINTS = 2;

struct PCCOctree3Node {
  // 3D position of the current node's origin (local x,y,z = 0).
  PCCVector3<uint32_t> pos;

  // Range of point indexes spanned by node
  uint32_t start;
  uint32_t end;

  // address of the current node in 3D morton order.
  uint64_t mortonIdx;

  // pattern denoting occupied neighbour nodes.
  //    32 8 (y)
  //     |/
  //  2--n--1 (x)
  //    /|
  //   4 16 (z)
  uint8_t neighPattern = 0;

  // The current node's number of siblings plus one.
  // ie, the number of child nodes present in this node's parent.
  uint8_t numSiblingsPlus1;

  // The occupancy map used describing the current node and its siblings.
  uint8_t siblingOccupancy;
};

// Structure for sorting weights.
struct WeightWithIndex {
  float weight;
  int index;

  WeightWithIndex() = default;

  WeightWithIndex(const int index, const float weight)
    : weight(weight), index(index)
  {}

  // NB: this definition ranks larger weights before smaller values.
  bool operator<(const WeightWithIndex& rhs) const
  {
    // NB: index used to maintain stable sort
    if (weight == rhs.weight)
      return index < rhs.index;
    return weight > rhs.weight;
  }
};

struct PCCNeighborInfo {
  double weight;
  uint32_t index;
  uint32_t predictorIndex;
};

struct PCCPredictor {
  size_t neighborCount;
  uint32_t index;
  uint32_t levelOfDetailIndex;
  PCCNeighborInfo neighbors[PCCTMC3MaxPredictionNearestNeighborCount];

  PCCColor3B predictColor(const PCCPointSet3& pointCloud) const
  {
    PCCVector3D predicted(0.0);
    for (size_t i = 0; i < neighborCount; ++i) {
      const PCCColor3B color = pointCloud.getColor(neighbors[i].index);
      const double w = neighbors[i].weight;
      for (size_t k = 0; k < 3; ++k) {
        predicted[k] += w * color[k];
      }
    }
    return PCCColor3B(
      uint8_t(std::round(predicted[0])), uint8_t(std::round(predicted[1])),
      uint8_t(std::round(predicted[2])));
  }
  uint16_t predictReflectance(const PCCPointSet3& pointCloud) const
  {
    double predicted(0.0);
    for (size_t i = 0; i < neighborCount; ++i) {
      predicted +=
        neighbors[i].weight * pointCloud.getReflectance(neighbors[i].index);
    }
    return uint16_t(std::round(predicted));
  }

  void computeWeights()
  {
    if (neighborCount <= 1) {
      neighbors[0].weight = 1.0;
      return;
    }
    double sum = 0.0;
    for (size_t n = 0; n < neighborCount; ++n) {
      sum += neighbors[n].weight;
    }
    assert(sum > 0.0);
    for (size_t n = 0; n < neighborCount; ++n) {
      neighbors[n].weight /= sum;
    }
  }

  void init(
    const uint32_t current,
    const uint32_t reference,
    const uint32_t predictorIndex,
    const uint32_t lodIndex)
  {
    index = current;
    neighborCount = (reference != PCC_UNDEFINED_INDEX) ? 1 : 0;
    neighbors[0].index = reference;
    neighbors[0].predictorIndex = predictorIndex;
    neighbors[0].weight = 1.0;
    levelOfDetailIndex = lodIndex;
  }
};
inline int64_t
PCCQuantization(const int64_t value, const int64_t qs)
{
  const int64_t shift = (qs / 3);
  if (!qs) {
    return value;
  }
  if (value >= 0) {
    return (value + shift) / qs;
  }
  return -((shift - value) / qs);
}

inline int64_t
PCCInverseQuantization(const int64_t value, const int64_t qs)
{
  return qs == 0 ? value : (value * qs);
}

inline int64_t
PCCQuantization(const double value, const int64_t qs)
{
  return int64_t(
    value >= 0.0 ? std::floor(value / qs + 1.0 / 3.0)
                 : -std::floor(-value / qs + 1.0 / 3.0));
}

//---------------------------------------------------------------------------

template<typename T>
void
PCCLiftPredict(
  const std::vector<PCCPredictor>& predictors,
  const size_t startIndex,
  const size_t endIndex,
  const bool direct,
  std::vector<T>& attributes)
{
  const size_t predictorCount = endIndex - startIndex;
  for (size_t index = 0; index < predictorCount; ++index) {
    const size_t predictorIndex = predictorCount - index - 1 + startIndex;
    const auto& predictor = predictors[predictorIndex];
    T predicted(0.0);
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const size_t neighborPredIndex = predictor.neighbors[i].predictorIndex;
      const double weight = predictor.neighbors[i].weight;
      assert(neighborPredIndex < startIndex);
      predicted += weight * attributes[neighborPredIndex];
    }
    if (direct) {
      attributes[predictorIndex] -= predicted;
    } else {
      attributes[predictorIndex] += predicted;
    }
  }
}

//---------------------------------------------------------------------------

template<typename T>
void
PCCLiftUpdate(
  const std::vector<PCCPredictor>& predictors,
  const std::vector<double>& quantizationWeights,
  const size_t startIndex,
  const size_t endIndex,
  const bool direct,
  std::vector<T>& attributes)
{
  std::vector<double> updateWeights;
  updateWeights.resize(startIndex, 0.0);
  std::vector<T> updates;
  updates.resize(startIndex);
  for (size_t index = 0; index < startIndex; ++index) {
    updates[index] = 0.0;
  }
  const size_t predictorCount = endIndex - startIndex;
  for (size_t index = 0; index < predictorCount; ++index) {
    const size_t predictorIndex = predictorCount - index - 1 + startIndex;
    const auto& predictor = predictors[predictorIndex];
    const double currentQuantWeight = quantizationWeights[predictorIndex];
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const size_t neighborPredIndex = predictor.neighbors[i].predictorIndex;
      const double weight = predictor.neighbors[i].weight * currentQuantWeight;
      assert(neighborPredIndex < startIndex);
      updateWeights[neighborPredIndex] += weight;
      updates[neighborPredIndex] += weight * attributes[predictorIndex];
    }
  }
  for (size_t predictorIndex = 0; predictorIndex < startIndex;
       ++predictorIndex) {
    const double sumWeights = updateWeights[predictorIndex];
    if (sumWeights > 0.0) {
      const double alpha = 1.0 / sumWeights;
      if (direct) {
        attributes[predictorIndex] += alpha * updates[predictorIndex];
      } else {
        attributes[predictorIndex] -= alpha * updates[predictorIndex];
      }
    }
  }
}

//---------------------------------------------------------------------------

inline void
PCCComputeQuantizationWeights(
  const std::vector<PCCPredictor>& predictors,
  std::vector<double>& quantizationWeights)
{
  const size_t pointCount = predictors.size();
  quantizationWeights.resize(pointCount);
  for (size_t i = 0; i < pointCount; ++i) {
    quantizationWeights[i] = 1.0;
  }
  for (size_t i = 0; i < pointCount; ++i) {
    const size_t predictorIndex = pointCount - i - 1;
    const auto& predictor = predictors[predictorIndex];
    const double currentQuantWeight = quantizationWeights[predictorIndex];
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const size_t neighborPredIndex = predictor.neighbors[i].predictorIndex;
      const double weight = predictor.neighbors[i].weight;
      quantizationWeights[neighborPredIndex] += weight * currentQuantWeight;
    }
  }
}

//---------------------------------------------------------------------------

inline void
CheckPoint(
  const PCCPointSet3& pointCloud,
  const int32_t index,
  const double radius2,
  const int32_t searchRange1,
  const int32_t searchRange2,
  uint32_t& lastNNIndex,
  std::vector<bool>& visited,
  std::vector<uint32_t>& retained,
  PCCIncrementalKdTree3& kdtree)
{
  const auto& point = pointCloud[index];
  if (retained.empty()) {
    kdtree.insert(point);
    visited[index] = true;
    retained.push_back(index);
    return;
  }

  const int32_t retainedPointCount = retained.size();
  for (int32_t k = 0, j = retainedPointCount - 1; j >= 0 && k < searchRange1;
       --j, ++k) {
    if ((point - pointCloud[retained[j]]).getNorm2() < radius2) {
      return;
    }
  }

  if (lastNNIndex != PCC_UNDEFINED_INDEX) {
    for (int32_t k = 0; k < searchRange2; ++k) {
      const int32_t j1 = lastNNIndex - k;
      if (j1 > 0 && (point - pointCloud[retained[j1]]).getNorm2() < radius2) {
        lastNNIndex = j1;
        return;
      }
      const int32_t j2 = lastNNIndex + k;
      if (
        j2 < retainedPointCount
        && (point - pointCloud[retained[j2]]).getNorm2() < radius2) {
        lastNNIndex = j2;
        return;
      }
    }
  }
  const uint32_t nnIndex = kdtree.hasNeighborWithinRange(point, radius2);
  if (nnIndex != PCC_UNDEFINED_INDEX) {
    lastNNIndex = nnIndex;
    return;
  }
  kdtree.insert(point);
  visited[index] = true;
  retained.push_back(index);
}

//---------------------------------------------------------------------------

inline void
PCCBuildLevelOfDetail(
  const PCCPointSet3& pointCloud,
  const int levelOfDetailCount,
  const std::vector<int64_t>& dist2,
  std::vector<uint32_t>& numberOfPointsPerLOD,
  std::vector<uint32_t>& indexes)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<bool> visited(pointCount, false);
  indexes.resize(0);
  indexes.reserve(pointCount);
  numberOfPointsPerLOD.resize(0);
  numberOfPointsPerLOD.reserve(levelOfDetailCount);
  size_t prevIndexesSize = 0;
  const int32_t searchRange1 = 16;
  const int32_t searchRange2 = 8;
  PCCIncrementalKdTree3 kdtree;
  kdtree.reserve(pointCount / 2);
  for (size_t lodIndex = 0;
       lodIndex < levelOfDetailCount && indexes.size() < pointCount;
       ++lodIndex) {
    if ((lodIndex + 1) != levelOfDetailCount) {
      const double radius2 = dist2[lodIndex];
      uint32_t lastNNIndex = PCC_UNDEFINED_INDEX;
      kdtree.balance();
      for (uint32_t current = 0; current < pointCount; ++current) {
        if (!visited[current]) {
          CheckPoint(
            pointCloud, current, radius2, searchRange1, searchRange2,
            lastNNIndex, visited, indexes, kdtree);
        }
      }
    } else {
      for (uint32_t current = 0; current < pointCount; ++current) {
        if (!visited[current]) {
          indexes.push_back(current);
        }
      }
    }
    numberOfPointsPerLOD.push_back(uint32_t(indexes.size()));
    prevIndexesSize = indexes.size();
  }
}

struct PointCloudWrapper {
  PointCloudWrapper(
    const PCCPointSet3& pointCloud, const std::vector<uint32_t>& indexes)
    : _pointCloud(pointCloud), _indexes(indexes)
  {}
  inline size_t kdtree_get_point_count() const { return _pointCount; }
  inline void kdtree_set_point_count(const size_t pointCount)
  {
    assert(pointCount < _indexes.size());
    assert(pointCount < _pointCloud.getPointCount());
    _pointCount = pointCount;
  }
  inline double kdtree_get_pt(const size_t idx, int dim) const
  {
    assert(idx < _pointCount && dim < 3);
    return _pointCloud[_indexes[idx]][dim];
  }
  template<class BBOX>
  bool kdtree_get_bbox(BBOX& /* bb */) const
  {
    return false;
  }

  const PCCPointSet3& _pointCloud;
  const std::vector<uint32_t>& _indexes;
  size_t _pointCount = 0;
};

//---------------------------------------------------------------------------

inline void
PCCComputePredictors2(
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& numberOfPointsPerLOD,
  const std::vector<uint32_t>& indexes,
  const size_t numberOfNearestNeighborsInPrediction,
  std::vector<PCCPredictor>& predictors)
{
  const size_t pointCount = pointCloud.getPointCount();
  const size_t lodCount = numberOfPointsPerLOD.size();
  predictors.resize(pointCount);
  PCCPointDistInfo nearestNeighbors[PCCTMC3MaxPredictionNearestNeighborCount];
  PCCNNQuery3 nNQuery = {PCCVector3D(0.0), std::numeric_limits<double>::max(),
                         numberOfNearestNeighborsInPrediction};
  PCCNNResult nNResult = {nearestNeighbors, 0};
  PCCIncrementalKdTree3 kdtree;
  uint32_t i0 = 0;
  size_t balanceTargetPointCount = 16;
  for (uint32_t lodIndex = 0; lodIndex < lodCount; ++lodIndex) {
    const uint32_t i1 = numberOfPointsPerLOD[lodIndex];
    for (uint32_t i = i0; i < i1; ++i) {
      const uint32_t pointIndex = indexes[i];
      const auto& point = pointCloud[pointIndex];
      auto& predictor = predictors[i];
      nNQuery.point = point;
      kdtree.findNearestNeighbors2(nNQuery, nNResult);
      if (!nNResult.resultCount) {
        predictor.init(
          pointIndex, PCC_UNDEFINED_INDEX, PCC_UNDEFINED_INDEX,
          uint32_t(lodIndex));
      } else if (
        nNResult.neighbors[0].dist2 == 0.0 || nNResult.resultCount == 1) {
        const uint32_t reference =
          static_cast<uint32_t>(indexes[nNResult.neighbors[0].index]);
        assert(nNResult.neighbors[0].index < i);
        predictor.init(
          pointIndex, reference, nNResult.neighbors[0].index,
          uint32_t(lodIndex));
      } else {
        predictor.levelOfDetailIndex = uint32_t(lodIndex);
        predictor.index = pointIndex;
        predictor.neighborCount = nNResult.resultCount;
        for (size_t n = 0; n < nNResult.resultCount; ++n) {
          const uint32_t neighborIndex = nNResult.neighbors[n].index;
          assert(neighborIndex < i);
          assert(predictors[neighborIndex].levelOfDetailIndex <= lodIndex);
          predictor.neighbors[n].predictorIndex = neighborIndex;
          predictor.neighbors[n].index = indexes[neighborIndex];
          predictor.neighbors[n].weight = 1.0 / nNResult.neighbors[n].dist2;
        }
      }
      kdtree.insert(point);
      if (balanceTargetPointCount <= kdtree.size()) {
        kdtree.balance();
        balanceTargetPointCount = 2 * kdtree.size();
      }
    }
    i0 = i1;
  }
}

//---------------------------------------------------------------------------

inline void
PCCComputePredictors(
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& numberOfPointsPerLOD,
  const std::vector<uint32_t>& indexes,
  const size_t numberOfNearestNeighborsInPrediction,
  std::vector<PCCPredictor>& predictors)
{
  const size_t pointCount = pointCloud.getPointCount();
  const size_t lodCount = numberOfPointsPerLOD.size();
  assert(lodCount);
  predictors.resize(pointCount);

  // delta prediction for LOD0
  uint32_t i0 = numberOfPointsPerLOD[0];
  for (uint32_t i = 0; i < i0; ++i) {
    const uint32_t pointIndex = indexes[i];
    auto& predictor = predictors[i];
    if (i == 0) {
      predictor.init(pointIndex, PCC_UNDEFINED_INDEX, PCC_UNDEFINED_INDEX, 0);
    } else {
      predictor.init(pointIndex, predictors[i - 1].index, i - 1, 0);
    }
  }
  PointCloudWrapper pointCloudWrapper(pointCloud, indexes);
  const nanoflann::SearchParams params(10, 0.0f, true);
  size_t indices[PCCTMC3MaxPredictionNearestNeighborCount];
  double sqrDist[PCCTMC3MaxPredictionNearestNeighborCount];
  nanoflann::KNNResultSet<double> resultSet(
    numberOfNearestNeighborsInPrediction);
  for (uint32_t lodIndex = 1; lodIndex < lodCount; ++lodIndex) {
    pointCloudWrapper.kdtree_set_point_count(i0);
    nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, PointCloudWrapper>,
      PointCloudWrapper, 3>
      kdtree(
        3, pointCloudWrapper, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    kdtree.buildIndex();
    const uint32_t i1 = numberOfPointsPerLOD[lodIndex];

    for (uint32_t i = i0; i < i1; ++i) {
      const uint32_t pointIndex = indexes[i];
      const auto& point = pointCloud[pointIndex];
      auto& predictor = predictors[i];
      resultSet.init(indices, sqrDist);
      kdtree.findNeighbors(resultSet, &point[0], params);
      const uint32_t resultCount = resultSet.size();
      if (sqrDist[0] == 0.0 || resultCount == 1) {
        const uint32_t predIndex = indices[0];
        predictor.init(
          pointIndex, predictors[predIndex].index, predIndex, lodIndex);
      } else {
        predictor.levelOfDetailIndex = uint32_t(lodIndex);
        predictor.index = pointIndex;
        predictor.neighborCount = resultCount;
        for (size_t n = 0; n < resultCount; ++n) {
          const uint32_t predIndex = indices[n];
          assert(predIndex < i);
          assert(predictors[predIndex].levelOfDetailIndex <= lodIndex);
          predictor.neighbors[n].predictorIndex = predIndex;
          const uint32_t pointIndex1 = indexes[predIndex];
          predictor.neighbors[n].index = pointIndex1;
          const auto& point1 = pointCloud[pointIndex1];
          predictor.neighbors[n].weight = 1.0 / (point - point1).getNorm2();
        }
      }
    }
    i0 = i1;
  }
}

//---------------------------------------------------------------------------
// Update the neighbour pattern flags for a node and the 'left' neighbour on
// each axis.  This update should be applied to each newly inserted node.
//
// @param siblingRestriction limits neighbours to direct siblings of child

inline void
updateGeometryNeighState(
  bool siblingRestriction,
  const ringbuf<PCCOctree3Node>::iterator& bufEnd,
  int64_t numNodesNextLvl,
  int childSizeLog2,
  PCCOctree3Node& child,
  int childIdx,
  uint8_t neighPattern,
  uint8_t parentOccupancy)
{
  uint64_t midx;
  if (!siblingRestriction) {
    midx = child.mortonIdx = mortonAddr(child.pos, childSizeLog2);
  }

  static const struct {
    int childIdxBitPos;
    int axis;
    int patternFlagUs;
    int patternFlagThem;
  } neighParamMap[] = {
    {4, 2, 1 << 1, 1 << 0},  // x
    {2, 1, 1 << 2, 1 << 3},  // y
    {1, 0, 1 << 4, 1 << 5},  // z
  };

  for (const auto& param : neighParamMap) {
    // skip expensive check if parent's flags indicate adjacent neighbour
    // is not present.
    if ((childIdx & param.childIdxBitPos) == 0) {
      // $axis co-ordinate = 0
      if (parentOccupancy & (1 << (childIdx + param.childIdxBitPos)))
        child.neighPattern |= param.patternFlagThem;

      if (!(neighPattern & param.patternFlagUs))
        continue;
    } else {
      if (parentOccupancy & (1 << (childIdx - param.childIdxBitPos)))
        child.neighPattern |= param.patternFlagUs;

      // no external search is required for $axis co-ordinate = 1
      continue;
    }

    if (siblingRestriction)
      continue;

    // calculate the morton address of the 'left' neighbour,
    // the delta is then used as the starting position for a search
    int64_t mortonIdxNeigh =
      morton3dAxisDec(midx, param.axis) & ~0x8000000000000000ull;
    int64_t mortonDelta = midx - mortonIdxNeigh;

    if (mortonDelta < 0) {
      // no neighbour due to being in zero-th col/row/plane
      continue;
    }

    // NB: fifo already contains current node, no point searching it
    auto posEnd = bufEnd;
    std::advance(posEnd, -1);

    auto posStart = bufEnd;
    std::advance(posStart, -std::min(numNodesNextLvl, mortonDelta + 2));

    auto found = std::lower_bound(
      posStart, posEnd, mortonIdxNeigh,
      [](const PCCOctree3Node& node, uint64_t mortonIdx) {
        return node.mortonIdx < mortonIdx;
      });

    // NB: found is always valid (see posEnd) => can skip check.
    if (found->mortonIdx != mortonIdxNeigh) {
      // neighbour isn't present => must have been empty
      continue;
    }

    // update both node's neighbour pattern
    // NB: neighours being present implies occupancy
    child.neighPattern |= param.patternFlagUs;
    found->neighPattern |= param.patternFlagThem;
  }
}

//---------------------------------------------------------------------------
// Determine if direct coding is permitted.
// If tool is enabled:
//   - Block must not be near the bottom of the tree
//   - The parent / grandparent are sparsely occupied

inline bool
isDirectModeEligible(
  bool featureEnabled,
  int nodeSizeLog2,
  const PCCOctree3Node& node,
  const PCCOctree3Node& child)
{
  return featureEnabled && (nodeSizeLog2 >= 2) && (node.neighPattern == 0)
    && (child.numSiblingsPlus1 == 1) && (node.numSiblingsPlus1 <= 2);
}

//---------------------------------------------------------------------------

struct CtxModelOctreeOccupancy {
  o3dgc::Adaptive_Bit_Model b0[10];
  o3dgc::Adaptive_Bit_Model b1[2 * 10];
  o3dgc::Adaptive_Bit_Model b2[4 * 10];
  o3dgc::Adaptive_Bit_Model b3[8 * 10];
  o3dgc::Adaptive_Bit_Model b4[75];
  o3dgc::Adaptive_Bit_Model b5[112];
  o3dgc::Adaptive_Bit_Model b6[117];
  o3dgc::Adaptive_Bit_Model b7[127];

  CtxModelOctreeOccupancy()
  {
    for (auto& ctx : b0)
      ctx.reset(true);
    for (auto& ctx : b1)
      ctx.reset(true);
    for (auto& ctx : b2)
      ctx.reset(true);
    for (auto& ctx : b3)
      ctx.reset(true);
    for (auto& ctx : b4)
      ctx.reset(true);
    for (auto& ctx : b5)
      ctx.reset(true);
    for (auto& ctx : b6)
      ctx.reset(true);
    for (auto& ctx : b7)
      ctx.reset(true);
  }
};

//---------------------------------------------------------------------------

}  // namespace pcc

#endif /* PCCTMC3Common_h */
