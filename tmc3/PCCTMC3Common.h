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

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "constants.h"
#include "hls.h"

#include "nanoflann.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>

namespace pcc {

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

//---------------------------------------------------------------------------

struct MortonCodeWithIndex {
  int64_t mortonCode;
  int32_t index;
  bool operator<(const MortonCodeWithIndex& rhs) const
  {
    // NB: index used to maintain stable sort
    if (mortonCode == rhs.mortonCode)
      return index < rhs.index;
    return mortonCode < rhs.mortonCode;
  }
};

//---------------------------------------------------------------------------

struct PCCNeighborInfo {
  uint64_t weight;
  uint32_t predictorIndex;
  uint32_t insertIndex;
  bool operator<(const PCCNeighborInfo& rhs) const
  {
    return (weight == rhs.weight) ? insertIndex < rhs.insertIndex
                                  : weight < rhs.weight;
  }
};

//---------------------------------------------------------------------------

struct PCCPredictor {
  uint32_t neighborCount;
  PCCNeighborInfo neighbors[kAttributePredictionMaxNeighbourCount];
  int8_t predMode;
  int64_t maxDiff;

  Vec3<attr_t> predictColor(
    const PCCPointSet3& pointCloud, const std::vector<uint32_t>& indexes) const
  {
    Vec3<int64_t> predicted(0);
    if (predMode > neighborCount) {
      /* nop */
    } else if (predMode > 0) {
      const Vec3<attr_t> color =
        pointCloud.getColor(indexes[neighbors[predMode - 1].predictorIndex]);
      for (size_t k = 0; k < 3; ++k) {
        predicted[k] += color[k];
      }
    } else {
      for (size_t i = 0; i < neighborCount; ++i) {
        const Vec3<attr_t> color =
          pointCloud.getColor(indexes[neighbors[i].predictorIndex]);
        const uint32_t w = neighbors[i].weight;
        for (size_t k = 0; k < 3; ++k) {
          predicted[k] += w * color[k];
        }
      }
      for (uint32_t k = 0; k < 3; ++k) {
        predicted[k] =
          divExp2RoundHalfInf(predicted[k], kFixedPointWeightShift);
      }
    }
    return Vec3<attr_t>(predicted[0], predicted[1], predicted[2]);
  }

  int64_t predictReflectance(
    const PCCPointSet3& pointCloud, const std::vector<uint32_t>& indexes) const
  {
    int64_t predicted(0);
    if (predMode > neighborCount) {
      /* nop */
    } else if (predMode > 0) {
      predicted = pointCloud.getReflectance(
        indexes[neighbors[predMode - 1].predictorIndex]);
    } else {
      for (size_t i = 0; i < neighborCount; ++i) {
        predicted += neighbors[i].weight
          * pointCloud.getReflectance(indexes[neighbors[i].predictorIndex]);
      }
      predicted = divExp2RoundHalfInf(predicted, kFixedPointWeightShift);
    }
    return predicted;
  }

  void computeWeights()
  {
    const uint32_t shift = (1 << kFixedPointWeightShift);
    int32_t n = 0;
    while ((neighbors[0].weight >> n) >= shift) {
      ++n;
    }
    if (n > 0) {
      for (size_t i = 0; i < neighborCount; ++i) {
        neighbors[i].weight = (neighbors[i].weight + (1 << (n - 1))) >> n;
      }
    }
    while (neighborCount > 1) {
      if (
        neighbors[neighborCount - 1].weight
        >= (neighbors[neighborCount - 2].weight << kFixedPointWeightShift)) {
        --neighborCount;
      } else {
        break;
      }
    }
    if (neighborCount <= 1) {
      neighbors[0].weight = shift;
    } else if (neighborCount == 2) {
      const uint64_t d0 = neighbors[0].weight;
      const uint64_t d1 = neighbors[1].weight;
      const uint64_t sum = d1 + d0;
      const uint64_t w1 = (d0 << kFixedPointWeightShift) / sum;
      const uint64_t w0 = shift - w1;
      neighbors[0].weight = uint32_t(w0);
      neighbors[1].weight = uint32_t(w1);
    } else {
      neighborCount = 3;
      const uint64_t d0 = neighbors[0].weight;
      const uint64_t d1 = neighbors[1].weight;
      const uint64_t d2 = neighbors[2].weight;
      const uint64_t sum = d1 * d2 + d0 * d2 + d0 * d1;
      const uint64_t w2 = ((d0 * d1) << kFixedPointWeightShift) / sum;
      const uint64_t w1 = ((d0 * d2) << kFixedPointWeightShift) / sum;
      const uint64_t w0 = shift - (w1 + w2);
      neighbors[0].weight = uint32_t(w0);
      neighbors[1].weight = uint32_t(w1);
      neighbors[2].weight = uint32_t(w2);
    }
  }

  void init(const uint32_t predictorIndex)
  {
    neighborCount = (predictorIndex != PCC_UNDEFINED_INDEX) ? 1 : 0;
    neighbors[0].predictorIndex = predictorIndex;
    neighbors[0].weight = 1;
    predMode = 0;
  }

  void init()
  {
    neighborCount = 0;
    memset(
      neighbors, 0,
      sizeof(PCCNeighborInfo) * kAttributePredictionMaxNeighbourCount);
  }

  void insertNeighbor(
    const uint32_t reference,
    const uint64_t weight,
    const uint32_t maxNeighborCount,
    const uint32_t insertIndex)
  {
    bool sort = false;
    assert(
      maxNeighborCount > 0
      && maxNeighborCount <= kAttributePredictionMaxNeighbourCount);
    if (neighborCount < maxNeighborCount) {
      PCCNeighborInfo& neighborInfo = neighbors[neighborCount];
      neighborInfo.weight = weight;
      neighborInfo.predictorIndex = reference;
      neighborInfo.insertIndex = insertIndex;
      ++neighborCount;
      sort = true;
    } else {
      PCCNeighborInfo& neighborInfo = neighbors[maxNeighborCount - 1];
      if (
        weight < neighborInfo.weight
        || (weight == neighborInfo.weight
            && insertIndex < neighborInfo.insertIndex)) {
        neighborInfo.weight = weight;
        neighborInfo.predictorIndex = reference;
        neighborInfo.insertIndex = insertIndex;
        sort = true;
      }
    }
    for (int32_t k = neighborCount - 1; k > 0 && sort; --k) {
      if (neighbors[k] < neighbors[k - 1])
        std::swap(neighbors[k], neighbors[k - 1]);
      else
        return;
    }
  }
};

//---------------------------------------------------------------------------

inline void
PCCLiftPredict(
  const std::vector<PCCPredictor>& predictors,
  const size_t startIndex,
  const size_t endIndex,
  const bool direct,
  std::vector<Vec3<int64_t>>& attributes)
{
  const size_t predictorCount = endIndex - startIndex;
  for (size_t index = 0; index < predictorCount; ++index) {
    const size_t predictorIndex = predictorCount - index - 1 + startIndex;
    const auto& predictor = predictors[predictorIndex];
    auto& attribute = attributes[predictorIndex];
    Vec3<int64_t> predicted(int64_t(0));
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const size_t neighborPredIndex = predictor.neighbors[i].predictorIndex;
      const uint32_t weight = predictor.neighbors[i].weight;
      assert(neighborPredIndex < startIndex);
      predicted += weight * attributes[neighborPredIndex];
    }
    predicted = divExp2RoundHalfInf(predicted, kFixedPointWeightShift);
    if (direct) {
      attribute -= predicted;
    } else {
      attribute += predicted;
    }
  }
}

inline void
PCCLiftPredict(
  const std::vector<PCCPredictor>& predictors,
  const size_t startIndex,
  const size_t endIndex,
  const bool direct,
  std::vector<int64_t>& attributes)
{
  const size_t predictorCount = endIndex - startIndex;
  for (size_t index = 0; index < predictorCount; ++index) {
    const size_t predictorIndex = predictorCount - index - 1 + startIndex;
    const auto& predictor = predictors[predictorIndex];
    auto& attribute = attributes[predictorIndex];
    int64_t predicted(int64_t(0));
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const size_t neighborPredIndex = predictor.neighbors[i].predictorIndex;
      const uint32_t weight = predictor.neighbors[i].weight;
      assert(neighborPredIndex < startIndex);
      predicted += weight * attributes[neighborPredIndex];
    }
    predicted = divExp2RoundHalfInf(predicted, kFixedPointWeightShift);
    if (direct) {
      attribute -= predicted;
    } else {
      attribute += predicted;
    }
  }
}

//---------------------------------------------------------------------------

inline void
PCCLiftUpdate(
  const std::vector<PCCPredictor>& predictors,
  const std::vector<uint64_t>& quantizationWeights,
  const size_t startIndex,
  const size_t endIndex,
  const bool direct,
  std::vector<Vec3<int64_t>>& attributes)
{
  std::vector<uint64_t> updateWeights;
  updateWeights.resize(startIndex, uint64_t(0));
  std::vector<Vec3<int64_t>> updates;
  updates.resize(startIndex);
  for (size_t index = 0; index < startIndex; ++index) {
    updates[index] = int64_t(0);
  }
  const size_t predictorCount = endIndex - startIndex;
  for (size_t index = 0; index < predictorCount; ++index) {
    const size_t predictorIndex = predictorCount - index - 1 + startIndex;
    const auto& predictor = predictors[predictorIndex];
    const auto currentQuantWeight = quantizationWeights[predictorIndex];
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const size_t neighborPredIndex = predictor.neighbors[i].predictorIndex;
      const uint64_t weight =
        predictor.neighbors[i].weight * currentQuantWeight;
      assert(neighborPredIndex < startIndex);
      updateWeights[neighborPredIndex] += weight;
      updates[neighborPredIndex] += weight * attributes[predictorIndex];
    }
  }
  for (size_t predictorIndex = 0; predictorIndex < startIndex;
       ++predictorIndex) {
    const uint32_t sumWeights = updateWeights[predictorIndex];
    if (sumWeights) {
      auto& update = updates[predictorIndex];
      update = (update + sumWeights / 2) / sumWeights;
      auto& attribute = attributes[predictorIndex];
      if (direct) {
        attribute += update;
      } else {
        attribute -= update;
      }
    }
  }
}

inline void
PCCLiftUpdate(
  const std::vector<PCCPredictor>& predictors,
  const std::vector<uint64_t>& quantizationWeights,
  const size_t startIndex,
  const size_t endIndex,
  const bool direct,
  std::vector<int64_t>& attributes)
{
  std::vector<uint64_t> updateWeights;
  updateWeights.resize(startIndex, uint64_t(0));
  std::vector<int64_t> updates;
  updates.resize(startIndex);
  for (size_t index = 0; index < startIndex; ++index) {
    updates[index] = int64_t(0);
  }
  const size_t predictorCount = endIndex - startIndex;
  for (size_t index = 0; index < predictorCount; ++index) {
    const size_t predictorIndex = predictorCount - index - 1 + startIndex;
    const auto& predictor = predictors[predictorIndex];
    const auto currentQuantWeight = quantizationWeights[predictorIndex];
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const size_t neighborPredIndex = predictor.neighbors[i].predictorIndex;
      const uint64_t weight =
        predictor.neighbors[i].weight * currentQuantWeight;
      assert(neighborPredIndex < startIndex);
      updateWeights[neighborPredIndex] += weight;
      updates[neighborPredIndex] += weight * attributes[predictorIndex];
    }
  }
  for (size_t predictorIndex = 0; predictorIndex < startIndex;
       ++predictorIndex) {
    const uint32_t sumWeights = updateWeights[predictorIndex];
    if (sumWeights > 0.0) {
      auto& update = updates[predictorIndex];
      update = (update + sumWeights / 2) / sumWeights;
      auto& attribute = attributes[predictorIndex];
      if (direct) {
        attribute += update;
      } else {
        attribute -= update;
      }
    }
  }
}

//---------------------------------------------------------------------------

inline void
PCCComputeQuantizationWeights(
  const std::vector<PCCPredictor>& predictors,
  std::vector<uint64_t>& quantizationWeights)
{
  const size_t pointCount = predictors.size();
  quantizationWeights.resize(pointCount);
  for (size_t i = 0; i < pointCount; ++i) {
    quantizationWeights[i] = (1 << kFixedPointWeightShift);
  }
  for (size_t i = 0; i < pointCount; ++i) {
    const size_t predictorIndex = pointCount - i - 1;
    const auto& predictor = predictors[predictorIndex];
    const auto currentQuantWeight = quantizationWeights[predictorIndex];
    for (size_t j = 0; j < predictor.neighborCount; ++j) {
      const size_t neighborPredIndex = predictor.neighbors[j].predictorIndex;
      const auto weight = predictor.neighbors[j].weight;
      auto& neighborQuantWeight = quantizationWeights[neighborPredIndex];
      neighborQuantWeight += divExp2RoundHalfInf(
        weight * currentQuantWeight, kFixedPointWeightShift);
    }
  }
  for (auto& w : quantizationWeights) {
    w = isqrt(w);
  }
}

//---------------------------------------------------------------------------

inline void
computeQuantizationWeightsScalable(
  const std::vector<PCCPredictor>& predictors,
  const std::vector<uint32_t>& numberOfPointsPerLOD,
  size_t geom_num_points,
  int32_t minGeomNodeSizeLog2,
  std::vector<uint64_t>& quantizationWeights)
{
  const size_t pointCount = predictors.size();
  quantizationWeights.resize(pointCount);
  for (size_t i = 0; i < pointCount; ++i) {
    quantizationWeights[i] = (1 << kFixedPointWeightShift);
  }

  const size_t lodCount = numberOfPointsPerLOD.size();
  for (size_t lodIndex = 0; lodIndex < lodCount; ++lodIndex) {
    const size_t startIndex =
      (lodIndex == 0) ? 0 : numberOfPointsPerLOD[lodIndex - 1];
    const size_t endIndex = numberOfPointsPerLOD[lodIndex];

    const size_t predictorCount = endIndex - startIndex;
    for (size_t index = 0; index < predictorCount; ++index) {
      const size_t predictorIndex = index + startIndex;
      const double currentQuantWeight = (geom_num_points / predictorCount);

      if (!minGeomNodeSizeLog2 && (lodIndex == lodCount - 1)) {
        quantizationWeights[predictorIndex] = (1 << kFixedPointWeightShift);
      } else {
        quantizationWeights[predictorIndex] =
          currentQuantWeight * (1 << kFixedPointWeightShift);
      }
    }
  }

  for (auto& w : quantizationWeights) {
    w = isqrt(w);
  }
}

//---------------------------------------------------------------------------

inline uint32_t
FindNeighborWithinDistance(
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const int32_t index,
  const double radius2,
  const int32_t searchRange,
  std::vector<uint32_t>& retained)
{
  const auto& point = pointCloud[packedVoxel[index].index];
  const int32_t retainedSize = retained.size();
  int32_t j = retainedSize - 2;
  int32_t k = 0;
  while (j >= 0 && ++k < searchRange) {
    const int32_t index1 = retained[j];
    const int32_t pointIndex1 = packedVoxel[index1].index;
    const auto& point1 = pointCloud[pointIndex1];
    const auto d2 = (point1 - point).getNorm2<double>();
    if (d2 <= radius2) {
      return index1;
    }
    --j;
  }
  return PCC_UNDEFINED_INDEX;
}

//---------------------------------------------------------------------------

inline bool
checkDistance(
  const point_t& point, const point_t& refpoint, int32_t nodeSizeLog2)
{
  uint32_t mask = uint32_t(-1) << nodeSizeLog2;

  int32_t minX = (refpoint.x() & mask) * 2 - 1;
  int32_t minY = (refpoint.y() & mask) * 2 - 1;
  int32_t minZ = (refpoint.z() & mask) * 2 - 1;

  int32_t maxX = minX + (2 << nodeSizeLog2);
  int32_t maxY = minY + (2 << nodeSizeLog2);
  int32_t maxZ = minZ + (2 << nodeSizeLog2);

  int32_t x = point.x() * 2;
  int32_t y = point.y() * 2;
  int32_t z = point.z() * 2;

  return minX < x && maxX > x && minY < y && maxY > y && minZ < z && maxZ > z;
}

//---------------------------------------------------------------------------

inline uint32_t
findNeighborWithinVoxel(
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  int32_t index,
  int32_t octreeNodeSizeLog2,
  int32_t searchRange,
  std::vector<uint32_t>& retained)
{
  const auto& point = pointCloud[packedVoxel[index].index];
  int32_t retainedSize = retained.size();
  int32_t j = retainedSize - 1;
  int32_t k = 0;
  while (j >= 0 && ++k < searchRange) {
    int32_t index1 = retained[j];
    int32_t pointIndex1 = packedVoxel[index1].index;
    const auto& point1 = pointCloud[pointIndex1];

    if (checkDistance(point, point1, octreeNodeSizeLog2))
      return index1;
    --j;
  }
  return PCC_UNDEFINED_INDEX;
}

//---------------------------------------------------------------------------

inline point_t
clacIntermediatePosition(
  bool enabled, int32_t nodeSizeLog2, const point_t& point)
{
  if (!enabled || !nodeSizeLog2)
    return point;

  uint32_t mask = (uint32_t(-1)) << nodeSizeLog2;
  int32_t centerX = (point.x() & mask) + (1 << (nodeSizeLog2 - 1));
  int32_t centerY = (point.y() & mask) + (1 << (nodeSizeLog2 - 1));
  int32_t centerZ = (point.z() & mask) + (1 << (nodeSizeLog2 - 1));

  point_t newPoint{centerX, centerY, centerZ};

  return newPoint;
}

//---------------------------------------------------------------------------

inline int
indexTieBreaker(int a, int b)
{
  return a > b ? ((a - b) << 1) - 1 : ((b - a) << 1);
}

//---------------------------------------------------------------------------

inline void
updateNearestNeighbor(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const int32_t nodeSizeLog2,
  const int32_t predictorIndex,
  const point_t& point,
  PCCPredictor& predictor)
{
  const int32_t pointIndex1 = packedVoxel[predictorIndex].index;
  const auto point1 = clacIntermediatePosition(
    aps.scalable_lifting_enabled_flag, nodeSizeLog2, pointCloud[pointIndex1]);

  double norm2 = times(point - point1, aps.lod_neigh_bias).getNorm2<double>();

  if (nodeSizeLog2 > 0 && point == point1) {
    norm2 = double(1 << (nodeSizeLog2 - 1));
    norm2 = norm2 * norm2;
  }
  predictor.insertNeighbor(
    pointIndex1, norm2, aps.num_pred_nearest_neighbours, 0);
}

//---------------------------------------------------------------------------

inline void
computeNearestNeighbors(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& retained,
  const int32_t startIndex,
  const int32_t endIndex,
  const int32_t nodeSizeLog2,
  std::vector<uint32_t>& indexes,
  std::vector<PCCPredictor>& predictors,
  std::vector<uint32_t>& pointIndexToPredictorIndex,
  int32_t& predIndex,
  std::vector<Box3<int32_t>>& bBoxes)
{
  constexpr auto searchRangeNear = 2;
  const int32_t retainedSize = retained.size();
  const int32_t bucketSize = 8;
  bBoxes.resize((retainedSize + bucketSize - 1) / bucketSize);
  for (int32_t i = 0, b = 0; i < retainedSize; ++b) {
    auto& bBox = bBoxes[b];
    bBox.min = bBox.max = clacIntermediatePosition(
      aps.scalable_lifting_enabled_flag, nodeSizeLog2,
      pointCloud[packedVoxel[retained[i++]].index]);
    for (int32_t k = 1; k < bucketSize && i < retainedSize; ++k, ++i) {
      const int32_t pointIndex = packedVoxel[retained[i]].index;
      const auto point = clacIntermediatePosition(
        aps.scalable_lifting_enabled_flag, nodeSizeLog2,
        pointCloud[pointIndex]);
      for (int32_t p = 0; p < 3; ++p) {
        bBox.min[p] = std::min(bBox.min[p], point[p]);
        bBox.max[p] = std::max(bBox.max[p], point[p]);
      }
    }
  }

  std::vector<Box3<int32_t>> bBoxesI;
  const int32_t indexesSize = endIndex - startIndex;
  if (aps.intra_lod_prediction_enabled_flag) {
    bBoxesI.resize((indexesSize + bucketSize - 1) / bucketSize);
    for (int32_t i = startIndex, b = 0; i < endIndex; ++b) {
      auto& bBox = bBoxesI[b];
      bBox.min = bBox.max = pointCloud[packedVoxel[indexes[i++]].index];
      for (int32_t k = 1; k < bucketSize && i < endIndex; ++k, ++i) {
        const int32_t pointIndex = packedVoxel[indexes[i]].index;
        const auto& point = pointCloud[pointIndex];
        for (int32_t p = 0; p < 3; ++p) {
          bBox.min[p] = std::min(bBox.min[p], point[p]);
          bBox.max[p] = std::max(bBox.max[p], point[p]);
        }
      }
    }
  }

  const int32_t index0 = aps.num_pred_nearest_neighbours - 1;

  for (int32_t i = startIndex, j = 0; i < endIndex; ++i) {
    const int32_t index = indexes[i];
    const int64_t mortonCode = packedVoxel[index].mortonCode;
    const int32_t pointIndex = packedVoxel[index].index;
    const auto point = clacIntermediatePosition(
      aps.scalable_lifting_enabled_flag, nodeSizeLog2, pointCloud[pointIndex]);
    indexes[i] = pointIndex;
    while (j < retainedSize - 1
           && mortonCode >= packedVoxel[retained[j]].mortonCode)
      ++j;
    auto& predictor = predictors[--predIndex];
    pointIndexToPredictorIndex[pointIndex] = predIndex;

    predictor.init();

    const int32_t k0 = std::max(0, j - aps.search_range);
    const int32_t k1 = std::min(retainedSize - 1, j + aps.search_range);

    if (retainedSize)
      updateNearestNeighbor(
        aps, pointCloud, packedVoxel, nodeSizeLog2, retained[j], point,
        predictor);

    for (int32_t n = 1; n <= searchRangeNear; ++n) {
      const int32_t kp = j + n;
      if (kp <= k1) {
        updateNearestNeighbor(
          aps, pointCloud, packedVoxel, nodeSizeLog2, retained[kp], point,
          predictor);
      }
      const int32_t kn = j - n;
      if (kn >= k0) {
        updateNearestNeighbor(
          aps, pointCloud, packedVoxel, nodeSizeLog2, retained[kn], point,
          predictor);
      }
    }

    const int32_t p0 = j - searchRangeNear - 1;
    const int32_t p1 = j + searchRangeNear + 1;
    // process p1..k1
    const int32_t bucketIndex1 = (k1 + bucketSize - 1) / bucketSize;
    for (int32_t bucketIndex = p1 / bucketSize; bucketIndex < bucketIndex1;
         ++bucketIndex) {
      if (
        predictor.neighborCount < aps.num_pred_nearest_neighbours
        || bBoxes[bucketIndex].getDist2<int64_t>(point)
          <= predictor.neighbors[index0].weight) {
        const int32_t indexBucketAligned = bucketIndex * bucketSize;
        const int32_t h0 = std::max(p1, indexBucketAligned);
        const int32_t h1 = std::min(k1, indexBucketAligned + bucketSize - 1);
        for (int32_t k = h0; k <= h1; ++k) {
          updateNearestNeighbor(
            aps, pointCloud, packedVoxel, nodeSizeLog2, retained[k], point,
            predictor);
        }
      }
    }

    // process k0..p0
    const int32_t bucketIndex0 = k0 / bucketSize;
    for (int32_t bucketIndex = p0 / bucketSize; bucketIndex >= bucketIndex0;
         --bucketIndex) {
      if (
        predictor.neighborCount < aps.num_pred_nearest_neighbours
        || bBoxes[bucketIndex].getDist2<int64_t>(point)
          <= predictor.neighbors[index0].weight) {
        const int32_t indexBucketAligned = bucketIndex * bucketSize;
        const int32_t h0 = std::max(k0, indexBucketAligned);
        const int32_t h1 = std::min(p0, indexBucketAligned + bucketSize - 1);
        for (int32_t k = h1; k >= h0; --k) {
          updateNearestNeighbor(
            aps, pointCloud, packedVoxel, nodeSizeLog2, retained[k], point,
            predictor);
        }
      }
    }

    if (aps.intra_lod_prediction_enabled_flag) {
      const int32_t k00 = i + 1;
      const int32_t k01 = std::min(endIndex - 1, k00 + searchRangeNear);
      for (int32_t k = k00; k <= k01; ++k) {
        updateNearestNeighbor(
          aps, pointCloud, packedVoxel, nodeSizeLog2, indexes[k], point,
          predictor);
      }

      const int32_t k0 = k01 + 1 - startIndex;
      const int32_t k1 =
        std::min(endIndex - 1, k00 + aps.search_range) - startIndex;
      const int32_t bucketIndex0 = k0 / bucketSize;
      const int32_t bucketIndex1 = (k1 + bucketSize - 1) / bucketSize;
      for (int32_t bucketIndex = bucketIndex0; bucketIndex < bucketIndex1;
           ++bucketIndex) {
        if (
          predictor.neighborCount < aps.num_pred_nearest_neighbours
          || bBoxesI[bucketIndex].getDist2<int64_t>(point)
            <= predictor.neighbors[index0].weight) {
          const int32_t indexBucketAligned = bucketIndex * bucketSize;
          const int32_t h0 = std::max(k0, indexBucketAligned);
          const int32_t h1 = std::min(k1, indexBucketAligned + bucketSize);
          for (int32_t h = h0; h < h1; ++h) {
            const int32_t k = startIndex + h;
            updateNearestNeighbor(
              aps, pointCloud, packedVoxel, nodeSizeLog2, indexes[k], point,
              predictor);
          }
        }
      }
    }
    assert(predictor.neighborCount <= aps.num_pred_nearest_neighbours);
  }
}

//---------------------------------------------------------------------------

inline void
subsampleByDistance(
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  const double radius2,
  const int32_t searchRange,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes)
{
  if (input.size() == 1) {
    indexes.push_back(input[0]);
  } else {
    for (const auto index : input) {
      if (retained.empty()) {
        retained.push_back(index);
        continue;
      }
      const auto& point = pointCloud[packedVoxel[index].index];
      const auto& lastRetained =
        pointCloud[packedVoxel[retained.back()].index];
      auto distanceToLastRetained = (lastRetained - point).getNorm2<double>();
      if (
        distanceToLastRetained <= radius2
        || FindNeighborWithinDistance(
             pointCloud, packedVoxel, index, radius2, searchRange, retained)
          != PCC_UNDEFINED_INDEX) {
        indexes.push_back(index);
      } else {
        retained.push_back(index);
      }
    }
  }
}

//---------------------------------------------------------------------------

inline void
subsampleByOctree(
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  int32_t octreeNodeSizeLog2,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes)
{
  int32_t nodeSizeLog2OfUpperLayer = octreeNodeSizeLog2 + 1;
  int32_t searchRange = 1;

  if (input.size() == 1) {
    indexes.push_back(input[0]);
  } else {
    for (const auto index : input) {
      if (retained.empty()) {
        retained.push_back(index);
        continue;
      }
      const auto& point = pointCloud[packedVoxel[index].index];
      const auto& retainedPoint =
        pointCloud[packedVoxel[retained.back()].index];

      if (
        checkDistance(point, retainedPoint, nodeSizeLog2OfUpperLayer)
        || findNeighborWithinVoxel(
             pointCloud, packedVoxel, index, nodeSizeLog2OfUpperLayer,
             searchRange, retained)
          != PCC_UNDEFINED_INDEX) {
        indexes.push_back(index);
      } else {
        retained.push_back(index);
      }
    }
  }
}

//---------------------------------------------------------------------------

inline void
subsampleByDecimation(
  const std::vector<uint32_t>& input,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes)
{
  static const int kLodUniformQuant = 4;
  const int indexCount = int(input.size());
  for (int i = 0; i < indexCount; ++i) {
    if (i % kLodUniformQuant == 0)
      retained.push_back(input[i]);
    else
      indexes.push_back(input[i]);
  }
}

//---------------------------------------------------------------------------

inline void
subsample(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  const int32_t lodIndex,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes)
{
  if (aps.scalable_lifting_enabled_flag) {
    int32_t octreeNodeSizeLog2 = lodIndex;
    subsampleByOctree(
      pointCloud, packedVoxel, input, octreeNodeSizeLog2, retained, indexes);
  } else if (aps.lod_decimation_enabled_flag) {
    subsampleByDecimation(input, retained, indexes);
  } else {
    double radius2 = aps.dist2[lodIndex];
    subsampleByDistance(
      pointCloud, packedVoxel, input, radius2, aps.search_range, retained,
      indexes);
  }
}

//---------------------------------------------------------------------------

inline void
computeMortonCodes(
  const PCCPointSet3& pointCloud,
  std::vector<MortonCodeWithIndex>& packedVoxel)
{
  const int32_t pointCount = int32_t(pointCloud.getPointCount());
  packedVoxel.resize(pointCount);
  for (int n = 0; n < pointCount; n++) {
    const auto& position = pointCloud[n];
    packedVoxel[n].mortonCode = mortonAddr(
      int32_t(position[0]), int32_t(position[1]), int32_t(position[2]));
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());
}

//---------------------------------------------------------------------------

inline void
updatePredictors(
  const std::vector<uint32_t>& pointIndexToPredictorIndex,
  std::vector<PCCPredictor>& predictors)
{
  for (auto& predictor : predictors) {
    if (predictor.neighborCount < 2) {
      predictor.neighbors[0].weight = 1;
    } else if (predictor.neighbors[0].weight == 0) {
      predictor.neighborCount = 1;
      predictor.neighbors[0].weight = 1;
    }
    for (int32_t k = 0; k < predictor.neighborCount; ++k) {
      auto& neighbor = predictor.neighbors[k];
      neighbor.predictorIndex =
        pointIndexToPredictorIndex[neighbor.predictorIndex];
    }
  }
}

//---------------------------------------------------------------------------

inline void
buildPredictorsFast(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  int32_t minGeomNodeSizeLog2,
  std::vector<PCCPredictor>& predictors,
  std::vector<uint32_t>& numberOfPointsPerLevelOfDetail,
  std::vector<uint32_t>& indexes)
{
  const int32_t pointCount = int32_t(pointCloud.getPointCount());
  assert(pointCount);

  std::vector<MortonCodeWithIndex> packedVoxel;
  computeMortonCodes(pointCloud, packedVoxel);

  std::vector<uint32_t> retained, input, pointIndexToPredictorIndex;
  pointIndexToPredictorIndex.resize(pointCount);
  retained.reserve(pointCount);
  input.resize(pointCount);
  for (uint32_t i = 0; i < pointCount; ++i) {
    input[i] = i;
  }

  // prepare output buffers
  predictors.resize(pointCount);
  numberOfPointsPerLevelOfDetail.resize(0);
  indexes.resize(0);
  indexes.reserve(pointCount);
  numberOfPointsPerLevelOfDetail.reserve(21);
  numberOfPointsPerLevelOfDetail.push_back(pointCount);

  int32_t num_detail_levels = aps.num_detail_levels;
  if (aps.scalable_lifting_enabled_flag == 1) {
    // NB: when partial decoding is enabled, LoDs correspond to octree levels
    num_detail_levels = std::numeric_limits<int>::max();
    ;
  }

  // prepare temporal buffer
  std::vector<uint32_t> indexesOfSbsample;
  indexesOfSbsample.resize(0);
  if (aps.scalable_lifting_enabled_flag) {
    indexesOfSbsample.reserve(pointCount);
  }

  std::vector<Box3<int32_t>> bBoxes;
  int32_t predIndex = int32_t(pointCount);
  for (auto lodIndex = minGeomNodeSizeLog2;
       !input.empty() && lodIndex <= num_detail_levels; ++lodIndex) {
    const int32_t startIndex = indexes.size();
    if (lodIndex == num_detail_levels) {
      for (const auto index : input) {
        indexes.push_back(index);
      }
    } else {
      subsample(
        aps, pointCloud, packedVoxel, input, lodIndex, retained, indexes);
    }
    const int32_t endIndex = indexes.size();
    if (aps.scalable_lifting_enabled_flag) {
      indexesOfSbsample.resize(endIndex);
      if (lodIndex == 0) {
        for (int32_t i = startIndex; i < endIndex; i++) {
          indexesOfSbsample[i] = indexes[i];
        }
      }
      if (lodIndex == 1) {
        if (endIndex - startIndex > startIndex) {
          for (int32_t i = 0; i < startIndex; i++) {
            indexes[i] = indexesOfSbsample[i];
          }
          predIndex = int32_t(pointCount);
          computeNearestNeighbors(
            aps, pointCloud, packedVoxel, retained, 0, startIndex,
            lodIndex - 1, indexes, predictors, pointIndexToPredictorIndex,
            predIndex, bBoxes);
        }
      }
    }

    computeNearestNeighbors(
      aps, pointCloud, packedVoxel, retained, startIndex, endIndex, lodIndex,
      indexes, predictors, pointIndexToPredictorIndex, predIndex, bBoxes);

    if (!retained.empty()) {
      numberOfPointsPerLevelOfDetail.push_back(retained.size());
    }
    input.resize(0);
    std::swap(retained, input);
  }
  std::reverse(indexes.begin(), indexes.end());
  updatePredictors(pointIndexToPredictorIndex, predictors);
  std::reverse(
    numberOfPointsPerLevelOfDetail.begin(),
    numberOfPointsPerLevelOfDetail.end());
}

//---------------------------------------------------------------------------

}  // namespace pcc

#endif /* PCCTMC3Common_h */
