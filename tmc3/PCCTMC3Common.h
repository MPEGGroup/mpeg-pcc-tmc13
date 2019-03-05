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

#include "PCCKdTree.h"
#include "PCCMath.h"
#include "constants.h"

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
  uint64_t mortonCode;
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
  bool operator<(const PCCNeighborInfo& rhs) const
  {
    return (weight == rhs.weight) ? predictorIndex < rhs.predictorIndex
                                  : weight < rhs.weight;
  }
};

//---------------------------------------------------------------------------

struct PCCPredictor {
  uint32_t neighborCount;
  PCCNeighborInfo neighbors[kAttributePredictionMaxNeighbourCount];
  int8_t predMode;
  int64_t maxDiff;

  Vec3<uint8_t> predictColor(
    const PCCPointSet3& pointCloud, const std::vector<uint32_t>& indexes) const
  {
    Vec3<int64_t> predicted(0);
    if (predMode > neighborCount) {
      /* nop */
    } else if (predMode > 0) {
      const Vec3<uint8_t> color =
        pointCloud.getColor(indexes[neighbors[predMode - 1].predictorIndex]);
      for (size_t k = 0; k < 3; ++k) {
        predicted[k] += color[k];
      }
    } else {
      for (size_t i = 0; i < neighborCount; ++i) {
        const Vec3<uint8_t> color =
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
    return Vec3<uint8_t>(predicted[0], predicted[1], predicted[2]);
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

  void init() { neighborCount = 0; }

  void insertNeighbor(
    const uint32_t reference,
    const uint64_t weight,
    const uint32_t maxNeighborCount)
  {
    bool sort = false;
    assert(
      maxNeighborCount > 0
      && maxNeighborCount <= kAttributePredictionMaxNeighbourCount);
    if (neighborCount < maxNeighborCount) {
      PCCNeighborInfo& neighborInfo = neighbors[neighborCount];
      neighborInfo.weight = weight;
      neighborInfo.predictorIndex = reference;
      ++neighborCount;
      sort = true;
    } else {
      PCCNeighborInfo& neighborInfo = neighbors[maxNeighborCount - 1];
      if (weight < neighborInfo.weight) {
        neighborInfo.weight = weight;
        neighborInfo.predictorIndex = reference;
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

inline int64_t
PCCQuantization(const int64_t value, const int64_t qs, bool isup = false)
{
  const int64_t shift = (qs / 3);
  if (!qs) {
    return value;
  }
  if (isup) {
    if (value >= 0) {
      return ((value << kFixedPointAttributeShift) + shift) / qs;
    }
    return -((shift - (value << kFixedPointAttributeShift)) / qs);
  } else if (value >= 0) {
    return (value + shift) / qs;
  }
  return -((shift - value) / qs);
}

inline int64_t
PCCInverseQuantization(
  const int64_t value, const int64_t qs, bool isdown = false)
{
  if (!qs)
    return value;

  if (isdown) {
    const int offset = 1 << (kFixedPointAttributeShift - 1);
    return (value * qs + offset) >> kFixedPointAttributeShift;
  }
  return value * qs;
}

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
    const auto d2 = (point1 - point).getNorm2();
    if (d2 <= radius2) {
      return index1;
    }
    --j;
  }
  return PCC_UNDEFINED_INDEX;
}

//---------------------------------------------------------------------------

inline void
computeNearestNeighbors(
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& retained,
  const int32_t startIndex,
  const int32_t endIndex,
  const int32_t searchRange,
  const int32_t numberOfNearestNeighborsInPrediction,
  std::vector<uint32_t>& indexes,
  std::vector<PCCPredictor>& predictors,
  std::vector<uint32_t>& pointIndexToPredictorIndex,
  int32_t& predIndex,
  std::vector<Box3<double>>& bBoxes)
{
  const int32_t retainedSize = retained.size();
  const int32_t bucketSize = 8;
  bBoxes.resize((retainedSize + bucketSize - 1) / bucketSize);
  for (int32_t i = 0, b = 0; i < retainedSize; ++b) {
    auto& bBox = bBoxes[b];
    bBox.min = bBox.max = pointCloud[packedVoxel[retained[i++]].index];
    for (int32_t k = 1; k < bucketSize && i < retainedSize; ++k, ++i) {
      const int32_t pointIndex = packedVoxel[retained[i]].index;
      const auto& point = pointCloud[pointIndex];
      for (int32_t p = 0; p < 3; ++p) {
        bBox.min[p] = std::min(bBox.min[p], point[p]);
        bBox.max[p] = std::max(bBox.max[p], point[p]);
      }
    }
  }

  const int32_t index0 = numberOfNearestNeighborsInPrediction - 1;

  for (int32_t i = startIndex, j = 0; i < endIndex; ++i) {
    const int32_t index = indexes[i];
    const uint64_t mortonCode = packedVoxel[index].mortonCode;
    const int32_t pointIndex = packedVoxel[index].index;
    const auto& point = pointCloud[pointIndex];
    indexes[i] = pointIndex;
    while (j < retainedSize
           && mortonCode >= packedVoxel[retained[j]].mortonCode)
      ++j;
    j = std::min(retainedSize - 1, j);
    auto& predictor = predictors[--predIndex];
    pointIndexToPredictorIndex[pointIndex] = predIndex;

    predictor.init();

    const int32_t j0 = std::max(0, j - searchRange);
    const int32_t j1 = std::min(retainedSize, j + searchRange + 1);

    const int32_t bucketIndex0 = j / bucketSize;
    int32_t k0 = std::max(bucketIndex0 * bucketSize, j0);
    int32_t k1 = std::min((bucketIndex0 + 1) * bucketSize, j1);

    for (int32_t k = k0; k < k1; ++k) {
      const int32_t pointIndex1 = packedVoxel[retained[k]].index;
      const auto& point1 = pointCloud[pointIndex1];
      predictor.insertNeighbor(
        pointIndex1, (point - point1).getNorm2(),
        numberOfNearestNeighborsInPrediction);
    }

    for (int32_t s0 = 1, sr = (1 + searchRange / bucketSize); s0 < sr; ++s0) {
      for (int32_t s1 = 0; s1 < 2; ++s1) {
        const int32_t bucketIndex1 =
          s1 == 0 ? bucketIndex0 + s0 : bucketIndex0 - s0;
        if (bucketIndex1 < 0 || bucketIndex1 >= bBoxes.size()) {
          continue;
        }
        if (
          predictor.neighborCount < numberOfNearestNeighborsInPrediction
          || bBoxes[bucketIndex1].getDist2(point)
            < predictor.neighbors[index0].weight) {
          const int32_t k0 = std::max(bucketIndex1 * bucketSize, j0);
          const int32_t k1 = std::min((bucketIndex1 + 1) * bucketSize, j1);
          for (int32_t k = k0; k < k1; ++k) {
            const int32_t pointIndex1 = packedVoxel[retained[k]].index;
            const auto& point1 = pointCloud[pointIndex1];
            predictor.insertNeighbor(
              pointIndex1, (point - point1).getNorm2(),
              numberOfNearestNeighborsInPrediction);
          }
        }
      }
    }
    assert(predictor.neighborCount <= numberOfNearestNeighborsInPrediction);
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
      if (
        (pointCloud[packedVoxel[retained.back()].index] - point).getNorm2()
          <= radius2
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
  bool useDecimation,
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  const double radius2,
  const int32_t searchRange,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes)
{
  if (useDecimation)
    subsampleByDecimation(input, retained, indexes);
  else
    subsampleByDistance(
      pointCloud, packedVoxel, input, radius2, searchRange, retained, indexes);
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
// LoD generation using Binary-tree
inline void
buildLevelOfDetailBinaryTree(
  const PCCPointSet3& pointCloud,
  std::vector<uint32_t>& numberOfPointsPerLOD,
  std::vector<uint32_t>& indexes)
{
  const uint32_t pointCount = pointCloud.getPointCount();
  uint8_t const btDepth = std::log2(round(pointCount / 2)) - 1;
  PCCKdTree3 kdtree(pointCloud, btDepth);
  kdtree.build();

  bool skipLayer = true;  //if true, skips alternate layer of BT
  std::vector<bool> skipDepth(btDepth + 1, false);
  if (skipLayer) {
    for (int i = skipDepth.size() - 2; i >= 0; i--) {
      skipDepth[i] = !skipDepth[i + 1];  //if true, that layer is skipped
    }
  }

  indexes.resize(pointCount);
  std::vector<bool> visited(pointCount, false);
  uint32_t lod = 0;
  uint32_t start = 0;
  uint32_t end = 0;
  for (size_t i = 0; i < btDepth + 1; i++) {
    start = end;
    end = start + (1 << i);
    if (!skipDepth[i]) {
      for (int j = start; j < end; j++) {
        auto indx = kdtree.searchClosestAvailablePoint(kdtree.nodes[j].centd);
        if (indx != PCC_UNDEFINED_INDEX && !visited[indx]) {
          indexes[lod++] = indx;
          visited[indx] = true;
        }
      }
      numberOfPointsPerLOD.push_back(lod);
    }
  }
  for (size_t i = 0; i < pointCount; i++) {
    if (!visited[i]) {
      indexes[lod++] = i;
    }
  }
  numberOfPointsPerLOD.push_back(lod);
}

//---------------------------------------------------------------------------
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
computePredictors(
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& numberOfPointsPerLOD,
  const std::vector<uint32_t>& indexes,
  const size_t numberOfNearestNeighborsInPrediction,
  std::vector<PCCPredictor>& predictors)
{
  const uint32_t PCCTMC3MaxPredictionNearestNeighborCount = 3;
  const size_t pointCount = pointCloud.getPointCount();
  const size_t lodCount = numberOfPointsPerLOD.size();
  assert(lodCount);
  predictors.resize(pointCount);

  // delta prediction for LOD0
  uint32_t i0 = numberOfPointsPerLOD[0];
  for (uint32_t i = 0; i < i0; ++i) {
    auto& predictor = predictors[i];
    if (i == 0) {
      predictor.init(PCC_UNDEFINED_INDEX);
    } else {
      predictor.init(PCC_UNDEFINED_INDEX);
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
        predictor.init(predIndex);
      } else {
        predictor.neighborCount = resultCount;
        for (size_t n = 0; n < resultCount; ++n) {
          const uint32_t predIndex = indices[n];
          assert(predIndex < i);
          predictor.neighbors[n].predictorIndex = predIndex;
          const uint32_t pointIndex1 = indexes[predIndex];
          const auto& point1 = pointCloud[pointIndex1];
          predictor.neighbors[n].weight = 1.0 / (point - point1).getNorm2();
        }
      }
    }
    i0 = i1;
  }
}

//---------------------------------------------------------------------------
inline void
buildPredictorsFast(
  const PCCPointSet3& pointCloud,
  bool lod_decimation_enabled_flag,
  const std::vector<int64_t>& dist2,
  const int32_t levelOfDetailCount,
  const int32_t numberOfNearestNeighborsInPrediction,
  const int32_t searchRange1,
  const int32_t searchRange2,
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

  std::vector<Box3<double>> bBoxes;
  int32_t predIndex = int32_t(pointCount);
  for (uint32_t lodIndex = 0; !input.empty() && lodIndex <= levelOfDetailCount;
       ++lodIndex) {
    const int32_t startIndex = indexes.size();
    if (lodIndex == levelOfDetailCount) {
      std::reverse(input.begin(), input.end());
      for (const auto index : input) {
        indexes.push_back(index);
      }
    } else {
      const double radius2 = dist2[lodIndex];
      subsample(
        lod_decimation_enabled_flag, pointCloud, packedVoxel, input, radius2,
        searchRange1, retained, indexes);
    }
    const int32_t endIndex = indexes.size();

    if (retained.empty()) {
      for (int32_t i = startIndex; i < endIndex; ++i) {
        const int32_t index = indexes[i];
        const int32_t pointIndex = packedVoxel[index].index;
        indexes[i] = pointIndex;
        auto& predictor = predictors[--predIndex];
        assert(predIndex >= 0);
        pointIndexToPredictorIndex[pointIndex] = predIndex;
        predictor.init(PCC_UNDEFINED_INDEX);
      }
      break;
    } else {
      computeNearestNeighbors(
        pointCloud, packedVoxel, retained, startIndex, endIndex, searchRange2,
        numberOfNearestNeighborsInPrediction, indexes, predictors,
        pointIndexToPredictorIndex, predIndex, bBoxes);
    }

    numberOfPointsPerLevelOfDetail.push_back(retained.size());
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

inline void
buildPredictorsFastNoLod(
  const PCCPointSet3& pointCloud,
  const int32_t numberOfNearestNeighborsInPrediction,
  const int32_t searchRange,
  std::vector<PCCPredictor>& predictors,
  std::vector<uint32_t>& indexes)
{
  const int32_t pointCount = int32_t(pointCloud.getPointCount());
  assert(pointCount);

  indexes.resize(pointCount);
  // re-order points
  {
    std::vector<MortonCodeWithIndex> packedVoxel;
    computeMortonCodes(pointCloud, packedVoxel);
    for (int32_t i = 0; i < pointCount; ++i) {
      indexes[i] = packedVoxel[i].index;
    }
  }
  predictors.resize(pointCount);
  for (int32_t i = 0; i < pointCount; ++i) {
    const int32_t index = indexes[i];
    const auto& point = pointCloud[index];
    auto& predictor = predictors[i];
    predictor.init();
    const int32_t k0 = std::max(0, i - searchRange);
    const int32_t k1 = i - 1;
    for (int32_t k = k1; k >= k0; --k) {
      const int32_t index1 = indexes[k];
      const auto& point1 = pointCloud[index1];
      predictor.insertNeighbor(
        k, (point - point1).getNorm2(), numberOfNearestNeighborsInPrediction);
    }
    assert(predictor.neighborCount <= numberOfNearestNeighborsInPrediction);
    if (predictor.neighborCount < 2) {
      predictor.neighbors[0].weight = 1;
    } else if (predictor.neighbors[0].weight == 0) {
      predictor.neighborCount = 1;
      predictor.neighbors[0].weight = 1;
    }
  }
}

//---------------------------------------------------------------------------

}  // namespace pcc

#endif /* PCCTMC3Common_h */
