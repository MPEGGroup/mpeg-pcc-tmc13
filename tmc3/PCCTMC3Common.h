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
  int regionQPOffset;
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
      const uint64_t w1 = divApprox(d0, sum, kFixedPointWeightShift);
      const uint64_t w0 = shift - w1;
      neighbors[0].weight = uint32_t(w0);
      neighbors[1].weight = uint32_t(w1);
    } else {
      neighborCount = 3;
      const uint64_t d0 = neighbors[0].weight;
      const uint64_t d1 = neighbors[1].weight;
      const uint64_t d2 = neighbors[2].weight;
      const uint64_t sum = d1 * d2 + d0 * d2 + d0 * d1;
      const uint64_t w2 = divApprox(d0 * d1, sum, kFixedPointWeightShift);
      const uint64_t w1 = divApprox(d0 * d2, sum, kFixedPointWeightShift);
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
};

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
    auto& attribute = attributes[predictorIndex];
    T predicted(T(0));
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

template<typename T>
void
PCCLiftUpdate(
  const std::vector<PCCPredictor>& predictors,
  const std::vector<uint64_t>& quantizationWeights,
  const size_t startIndex,
  const size_t endIndex,
  const bool direct,
  std::vector<T>& attributes)
{
  std::vector<uint64_t> updateWeights;
  updateWeights.resize(startIndex, uint64_t(0));
  std::vector<T> updates;
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
      update = divApprox(update, sumWeights, 0);
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
  size_t numPoints,
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
      const double currentQuantWeight =
        ((numPoints - startIndex) / predictorCount);

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

inline void
insertNeighbour(
  const uint32_t reference,
  const uint64_t weight,
  const uint32_t maxNeighborCountMinus1,
  const uint32_t insertIndex,
  uint32_t& neighborCount,
  PCCNeighborInfo* neighbors)
{
  bool sort = false;
  assert(
    maxNeighborCountMinus1 >= 0
    && maxNeighborCountMinus1 < kAttributePredictionMaxNeighbourCount);
  if (neighborCount <= maxNeighborCountMinus1) {
    PCCNeighborInfo& neighborInfo = neighbors[neighborCount];
    neighborInfo.weight = weight;
    neighborInfo.predictorIndex = reference;
    neighborInfo.insertIndex = insertIndex;
    ++neighborCount;
    sort = true;
  } else {
    PCCNeighborInfo& neighborInfo = neighbors[maxNeighborCountMinus1];
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

//---------------------------------------------------------------------------

inline void
updateNearestNeighbor(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const int32_t nodeSizeLog2,
  const int32_t predictorIndex,
  const point_t& point,
  uint32_t& neighborCount,
  PCCNeighborInfo* neighbors)
{
  const int32_t pointIndex1 = packedVoxel[predictorIndex].index;
  const auto point1 = clacIntermediatePosition(
    aps.scalable_lifting_enabled_flag, nodeSizeLog2, pointCloud[pointIndex1]);

  double norm1 = times(point - point1, aps.lodNeighBias).getNorm1();

  if (nodeSizeLog2 > 0 && point == point1) {
    norm1 = double(1 << (nodeSizeLog2 - 1));
  }

  insertNeighbour(
    pointIndex1, norm1, aps.num_pred_nearest_neighbours_minus1, 0,
    neighborCount, neighbors);
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

  const int32_t index0 = aps.num_pred_nearest_neighbours_minus1;

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

    PCCNeighborInfo localNeighbors[kAttributePredictionMaxNeighbourCount];
    uint32_t localNeighborCount = 0;

    const int32_t k0 = std::max(0, j - aps.search_range);
    const int32_t k1 = std::min(retainedSize - 1, j + aps.search_range);

    if (retainedSize)
      updateNearestNeighbor(
        aps, pointCloud, packedVoxel, nodeSizeLog2, retained[j], point,
        localNeighborCount, localNeighbors);

    for (int32_t n = 1; n <= searchRangeNear; ++n) {
      const int32_t kp = j + n;
      if (kp <= k1) {
        updateNearestNeighbor(
          aps, pointCloud, packedVoxel, nodeSizeLog2, retained[kp], point,
          localNeighborCount, localNeighbors);
      }
      const int32_t kn = j - n;
      if (kn >= k0) {
        updateNearestNeighbor(
          aps, pointCloud, packedVoxel, nodeSizeLog2, retained[kn], point,
          localNeighborCount, localNeighbors);
      }
    }

    const int32_t p0 = j - searchRangeNear - 1;
    const int32_t p1 = j + searchRangeNear + 1;
    // process p1..k1
    const int32_t bucketIndex1 = (k1 + bucketSize - 1) / bucketSize;
    for (int32_t bucketIndex = p1 / bucketSize; bucketIndex < bucketIndex1;
         ++bucketIndex) {
      if (
        localNeighborCount <= aps.num_pred_nearest_neighbours_minus1
        || bBoxes[bucketIndex].getDist1(point)
          <= localNeighbors[index0].weight) {
        const int32_t indexBucketAligned = bucketIndex * bucketSize;
        const int32_t h0 = std::max(p1, indexBucketAligned);
        const int32_t h1 = std::min(k1, indexBucketAligned + bucketSize - 1);
        for (int32_t k = h0; k <= h1; ++k) {
          updateNearestNeighbor(
            aps, pointCloud, packedVoxel, nodeSizeLog2, retained[k], point,
            localNeighborCount, localNeighbors);
        }
      }
    }

    // process k0..p0
    const int32_t bucketIndex0 = k0 / bucketSize;
    for (int32_t bucketIndex = p0 / bucketSize; bucketIndex >= bucketIndex0;
         --bucketIndex) {
      if (
        localNeighborCount <= aps.num_pred_nearest_neighbours_minus1
        || bBoxes[bucketIndex].getDist1(point)
          <= localNeighbors[index0].weight) {
        const int32_t indexBucketAligned = bucketIndex * bucketSize;
        const int32_t h0 = std::max(k0, indexBucketAligned);
        const int32_t h1 = std::min(p0, indexBucketAligned + bucketSize - 1);
        for (int32_t k = h1; k >= h0; --k) {
          updateNearestNeighbor(
            aps, pointCloud, packedVoxel, nodeSizeLog2, retained[k], point,
            localNeighborCount, localNeighbors);
        }
      }
    }

    if (aps.intra_lod_prediction_enabled_flag) {
      const int32_t k00 = i + 1;
      const int32_t k01 = std::min(endIndex - 1, k00 + searchRangeNear);
      for (int32_t k = k00; k <= k01; ++k) {
        updateNearestNeighbor(
          aps, pointCloud, packedVoxel, nodeSizeLog2, indexes[k], point,
          localNeighborCount, localNeighbors);
      }

      const int32_t k0 = k01 + 1 - startIndex;
      const int32_t k1 =
        std::min(endIndex - 1, k00 + aps.search_range) - startIndex;
      const int32_t bucketIndex0 = k0 / bucketSize;
      const int32_t bucketIndex1 = (k1 + bucketSize - 1) / bucketSize;
      for (int32_t bucketIndex = bucketIndex0; bucketIndex < bucketIndex1;
           ++bucketIndex) {
        if (
          localNeighborCount < aps.num_pred_nearest_neighbours_minus1
          || bBoxesI[bucketIndex].getDist1(point)
            <= localNeighbors[index0].weight) {
          const int32_t indexBucketAligned = bucketIndex * bucketSize;
          const int32_t h0 = std::max(k0, indexBucketAligned);
          const int32_t h1 = std::min(k1, indexBucketAligned + bucketSize);
          for (int32_t h = h0; h < h1; ++h) {
            const int32_t k = startIndex + h;
            updateNearestNeighbor(
              aps, pointCloud, packedVoxel, nodeSizeLog2, indexes[k], point,
              localNeighborCount, localNeighbors);
          }
        }
      }
    }

    assert(
      predictor.neighborCount <= aps.num_pred_nearest_neighbours_minus1 + 1);

    predictor.neighborCount = localNeighborCount;
    for (int i = 0; i < predictor.neighborCount; ++i) {
      predictor.neighbors[i] = localNeighbors[i];

      // use L2 norm for the final weight
      const int32_t pointIndex1 = localNeighbors[i].predictorIndex;
      const auto point1 = clacIntermediatePosition(
        aps.scalable_lifting_enabled_flag, nodeSizeLog2,
        pointCloud[pointIndex1]);

      double norm2 =
        times(point - point1, aps.lodNeighBias).getNorm2<double>();

      if (nodeSizeLog2 > 0 && point == point1) {
        norm2 = (double)(1 << (nodeSizeLog2 - 1));
        norm2 = norm2 * norm2;
      }
      predictor.neighbors[i].weight = norm2;
    }
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
partitionRange(
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes,
  int start,
  int end,
  bool retainFirst)
{
  int range = end - start;

  if (retainFirst) {
    retained.push_back(input[start]);
    for (int i = 1; i < range; ++i)
      indexes.push_back(input[start + i]);
  } else {
    // last index
    range--;
    for (int i = 0; i < range; ++i)
      indexes.push_back(input[start + i]);
    retained.push_back(input[start + range]);
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
  if (input.size() == 1) {
    indexes.push_back(input[0]);
    return;
  }

  bool retainFirst = !(octreeNodeSizeLog2 % 2);
  uint64_t kLodUniformQuant = pow(8, octreeNodeSizeLog2 + 1);
  const int indexCount = int(input.size());
  uint64_t endIdx = 0;
  uint64_t parentIdx = 0;
  int startIdx = 0;
  if (indexCount > 0) {
    auto mortonCode = packedVoxel[input[0]].mortonCode;
    endIdx = mortonCode + (kLodUniformQuant - (mortonCode % kLodUniformQuant));
  }

  for (int i = 1; i < indexCount; ++i) {
    int index = input[i];
    auto mortonCode = packedVoxel[index].mortonCode;
    if (mortonCode >= endIdx) {
      partitionRange(
        packedVoxel, input, retained, indexes, startIdx, i, retainFirst);
      startIdx = i;
      endIdx =
        mortonCode + (kLodUniformQuant - (mortonCode % kLodUniformQuant));
    }
  }

  if (startIdx < indexCount) {
    partitionRange(
      packedVoxel, input, retained, indexes, startIdx, indexCount,
      retainFirst);
  }
}

//---------------------------------------------------------------------------

inline void
subsampleByDecimation(
  const std::vector<uint32_t>& input,
  int lodSamplingPeriod,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes)
{
  const int indexCount = int(input.size());
  for (int i = 0, j = 1; i < indexCount; ++i) {
    if (--j)
      indexes.push_back(input[i]);
    else {
      retained.push_back(input[i]);
      j = lodSamplingPeriod;
    }
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
    auto samplingPeriod = aps.lodSamplingPeriod[lodIndex];
    subsampleByDecimation(input, samplingPeriod, retained, indexes);
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
  int geom_num_points_minus1,
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
  if (aps.scalable_lifting_enabled_flag) {
    // NB: when partial decoding is enabled, LoDs correspond to octree levels
    num_detail_levels = std::numeric_limits<int>::max();
  }

  bool concatenateLayers = aps.scalable_lifting_enabled_flag;
  std::vector<uint32_t> indexesOfSubsample;
  if (concatenateLayers)
    indexesOfSubsample.reserve(pointCount);

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

    if (concatenateLayers) {
      indexesOfSubsample.resize(endIndex);
      if (startIndex != endIndex) {
        for (int32_t i = startIndex; i < endIndex; i++)
          indexesOfSubsample[i] = indexes[i];

        int32_t numOfPointInSkipped = geom_num_points_minus1 + 1 - pointCount;
        if (endIndex - startIndex <= startIndex + numOfPointInSkipped) {
          concatenateLayers = false;
        } else {
          for (int32_t i = 0; i < startIndex; i++)
            indexes[i] = indexesOfSubsample[i];

          // reset predIndex
          predIndex = pointCount;
          for (int lod = 0; lod < lodIndex - minGeomNodeSizeLog2; lod++) {
            int divided_startIndex =
              pointCount - numberOfPointsPerLevelOfDetail[lod];
            int divided_endIndex =
              pointCount - numberOfPointsPerLevelOfDetail[lod + 1];

            computeNearestNeighbors(
              aps, pointCloud, packedVoxel, retained, divided_startIndex,
              divided_endIndex, lod + minGeomNodeSizeLog2, indexes, predictors,
              pointIndexToPredictorIndex, predIndex, bBoxes);
          }
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
