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
#include <memory>
#include <vector>

namespace pcc {

//============================================================================
// Hierachichal bounding boxes.
// Insert points (into the base layer), then generate the hierarchy via update.

template<int32_t BucketSizeLog2, int32_t LevelCount>
class BoxHierarchy {
public:
  void resize(const int32_t pointCount)
  {
    constexpr auto BucketSize = 1 << BucketSizeLog2;
    constexpr auto BucketSizeMinus1 = BucketSize - 1;
    int32_t count = pointCount;
    for (int i = 0; i < LevelCount; ++i) {
      count = (count + BucketSizeMinus1) >> BucketSizeLog2;
      _bBoxes[i].clear();
      _bBoxes[i].resize(count, Box3<int32_t>(INT32_MAX, INT32_MIN));
    }
  }

  void insert(const Vec3<int32_t>& point, const int32_t index)
  {
    const auto bindex = (index >> BucketSizeLog2);
    assert(bindex >= 0 && bindex < _bBoxes[0].size());
    _bBoxes[0][bindex].insert(point);
  }

  void update()
  {
    constexpr auto LevelCountMinus1 = LevelCount - 1;
    for (int i = 0; i < LevelCountMinus1; ++i) {
      for (int32_t j = 0, count = int32_t(_bBoxes[i].size()); j < count; ++j) {
        _bBoxes[i + 1][j >> BucketSizeLog2].merge(_bBoxes[i][j]);
      }
    }
  }

  const Box3<int32_t>& bBox(int32_t bindex, int32_t level) const
  {
    return _bBoxes[level][bindex];
  }

  int32_t bucketSizeLog2(int32_t level = 0) const
  {
    return BucketSizeLog2 * (1 + level);
  }

  int32_t bucketSize(int32_t level = 0) const
  {
    return 1 << bucketSizeLog2(level);
  }

private:
  std::vector<Box3<int32_t>> _bBoxes[LevelCount];
};

//============================================================================

class MortonIndexMap3d {
public:
  struct Range {
    int32_t start;
    int32_t end;
  };

  void resize(const int32_t cubeSizeLog2)
  {
    _cubeSizeLog2 = cubeSizeLog2;
    _cubeSize = 1 << cubeSizeLog2;
    _bufferSize = 1 << (3 * cubeSizeLog2);
    _mask = _bufferSize - 1;
    _buffer.reset(new Range[_bufferSize]);
  }

  void reserve(const uint32_t sz) { _updates.reserve(sz); }
  int cubeSize() const { return _cubeSize; }
  int cubeSizeLog2() const { return _cubeSizeLog2; }

  void init()
  {
    for (int32_t i = 0; i < _bufferSize; ++i) {
      _buffer[i] = {-1, -1};
    }
    _updates.resize(0);
  }

  void clearUpdates()
  {
    for (const auto index : _updates) {
      _buffer[index] = {-1, -1};
    }
    _updates.resize(0);
  }

  void set(const int64_t mortonCode, const int32_t index)
  {
    const int64_t mortonAddr = mortonCode & _mask;
    auto& unit = _buffer[mortonAddr];
    if (unit.start == -1) {
      unit.start = index;
    }
    unit.end = index + 1;
    _updates.push_back(mortonAddr);
  }

  Range get(const int64_t mortonCode) const
  {
    return _buffer[mortonCode & _mask];
  }

private:
  int32_t _cubeSize = 0;
  int32_t _cubeSizeLog2 = 0;
  int32_t _bufferSize = 0;
  int64_t _mask = 0;
  std::unique_ptr<Range[]> _buffer;

  // A list of indexes in _buffer that are dirty
  std::vector<int32_t> _updates;
};

//============================================================================

struct MortonCodeWithIndex {
  int64_t mortonCode;

  // The position used to generate the mortonCode
  Vec3<int32_t> position;

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
    return weight < rhs.weight;
  }
};

//---------------------------------------------------------------------------

struct PCCPredictor {
  uint32_t neighborCount;
  PCCNeighborInfo neighbors[kAttributePredictionMaxNeighbourCount];
  int8_t predMode;

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
        neighbors[i].weight = (neighbors[i].weight + (1ull << (n - 1))) >> n;
      }
    }
    while (neighborCount > 1) {
      if (
        neighbors[neighborCount - 1].weight
        >= (neighbors[0].weight << kFixedPointWeightShift)) {
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

  void
  blendWeights(const PCCPointSet3& cloud, const std::vector<uint32_t>& indexes)
  {
    int w0 = neighbors[0].weight;
    int w1 = neighbors[1].weight;
    int w2 = neighbors[2].weight;

    if (neighborCount != 3)
      return;

    auto neigh0Pos = cloud[indexes[neighbors[0].predictorIndex]];
    auto neigh1Pos = cloud[indexes[neighbors[1].predictorIndex]];
    auto neigh2Pos = cloud[indexes[neighbors[2].predictorIndex]];

    constexpr bool variant = 1;
    const auto d = variant ? 10 : 8;
    const auto bb = variant ? 1 : 4;
    const auto cc = variant ? 5 : 4;

    auto dist01 = (neigh0Pos - neigh1Pos).getNorm2<int64_t>();
    auto dist02 = (neigh0Pos - neigh2Pos).getNorm2<int64_t>();
    auto& dist10 = dist01;
    auto dist12 = (neigh1Pos - neigh2Pos).getNorm2<int64_t>();
    auto& dist20 = dist02;
    auto& dist21 = dist12;

    auto b1 = dist01 <= dist02 ? bb : cc;
    auto b2 = dist10 <= dist12 ? cc : bb;
    auto b3 = dist20 <= dist21 ? bb : cc;

    Vec3<int> w;
    w[0] = (w0 * d + w1 * (16 - d - b2) + w2 * b3) >> 4;
    w[1] = (w0 * b1 + w1 * d + w2 * (16 - d - b3)) >> 4;
    w[2] = 256 - w[0] - w[1];

    for (int i = 0; i < 3; i++)
      neighbors[i].weight = w[i];
  }

  void pruneDistanceGt(uint64_t maxDistance)
  {
    for (int i = 1; i < neighborCount; i++) {
      if (neighbors[i].weight > maxDistance) {
        neighborCount = i;
        break;
      }
    }
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
      const auto weight = divExp2RoundHalfInf(
        predictor.neighbors[i].weight * currentQuantWeight,
        kFixedPointWeightShift);
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
    const uint64_t currentQuantWeight =
      (numPoints / numberOfPointsPerLOD[lodIndex]) << kFixedPointWeightShift;

    const size_t predictorCount = endIndex - startIndex;
    for (size_t index = 0; index < predictorCount; ++index) {
      const size_t predictorIndex = index + startIndex;

      if (!minGeomNodeSizeLog2 && (lodIndex == lodCount - 1)) {
        quantizationWeights[predictorIndex] = (1 << kFixedPointWeightShift);
      } else {
        quantizationWeights[predictorIndex] = currentQuantWeight;
      }
    }
  }
}

//---------------------------------------------------------------------------

inline void
computeQuantizationWeights(
  const std::vector<PCCPredictor>& predictors,
  std::vector<int64_t>& quantizationWeights,
  Vec3<int32_t> neighWeight)
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
      auto& neighborQuantWeight = quantizationWeights[neighborPredIndex];
      neighborQuantWeight += divExp2RoundHalfInf(
        neighWeight[j] * currentQuantWeight, kFixedPointWeightShift);
    }
  }
}

//---------------------------------------------------------------------------

inline point_t
clacIntermediatePosition(
  bool enabled, int32_t nodeSizeLog2, const point_t& point)
{
  if (!enabled || !nodeSizeLog2)
    return point;

  uint32_t mask = (uint32_t(-1)) << nodeSizeLog2;
  int32_t centerX = point.x() & mask;
  int32_t centerY = point.y() & mask;
  int32_t centerZ = point.z() & mask;

  point_t newPoint{centerX, centerY, centerZ};

  return newPoint;
}

//---------------------------------------------------------------------------

inline void
updateNearestNeigh(
  const Vec3<int32_t>& point0,
  const Vec3<int32_t>& point1,
  int32_t index,
  int32_t (&localIndexes)[3],
  int64_t (&minDistances)[3])
{
  auto d = (point0 - point1).getNorm1();

  if (d >= minDistances[2]) {
    // do nothing
  } else if (d < minDistances[0]) {
    minDistances[2] = minDistances[1];
    minDistances[1] = minDistances[0];
    minDistances[0] = d;

    localIndexes[2] = localIndexes[1];
    localIndexes[1] = localIndexes[0];
    localIndexes[0] = index;
  } else if (d < minDistances[1]) {
    minDistances[2] = minDistances[1];
    minDistances[1] = d;
    localIndexes[2] = localIndexes[1];
    localIndexes[1] = index;
  } else {
    minDistances[2] = d;
    localIndexes[2] = index;
  }
}

//---------------------------------------------------------------------------

inline void
updateNearestNeighWithCheck(
  const Vec3<int32_t>& point0,
  const Vec3<int32_t>& point1,
  const int32_t index,
  int32_t (&localIndexes)[3],
  int64_t (&minDistances)[3])
{
  if (
    index == localIndexes[0] || index == localIndexes[1]
    || index == localIndexes[2])
    return;

  updateNearestNeigh(point0, point1, index, localIndexes, minDistances);
}

//---------------------------------------------------------------------------

inline void
computeNearestNeighbors(
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& retained,
  int32_t startIndex,
  int32_t endIndex,
  int32_t lodIndex,
  std::vector<uint32_t>& indexes,
  std::vector<PCCPredictor>& predictors,
  std::vector<uint32_t>& pointIndexToPredictorIndex,
  int32_t& predIndex,
  MortonIndexMap3d& atlas)
{
  constexpr auto searchRangeNear = 2;
  constexpr auto bucketSizeLog2 = 5;
  constexpr auto bucketSize = 1 << bucketSizeLog2;
  constexpr auto bucketSizeMinus1 = bucketSize - 1;
  constexpr auto levelCount = 3;

  const int32_t shiftBits = aps.scalable_lifting_enabled_flag
    ? 1 + lodIndex
    : 1 + aps.dist2 + abh.attr_dist2_delta + lodIndex;
  const int32_t shiftBits3 = 3 * shiftBits;
  const int32_t log2CubeSize = atlas.cubeSizeLog2();
  const int32_t atlasBits = 3 * log2CubeSize;
  // NB: when the atlas boundary is greater than 2^63, all points belong
  //     to a single atlas.  The clipping is necessary to avoid undefined
  //     behaviour of shifts greater than or equal to the word size.
  const int32_t atlasBoundaryBit = std::min(63, shiftBits3 + atlasBits);

  const int32_t retainedSize = retained.size();
  const int32_t indexesSize = endIndex - startIndex;
  const auto rangeInterLod = aps.inter_lod_search_range;
  const auto rangeIntraLod = aps.intra_lod_search_range;

  static const uint8_t kNeighOffset[27] = {
    7,   // { 0,  0,  0} 0
    3,   // {-1,  0,  0} 1
    5,   // { 0, -1,  0} 2
    6,   // { 0,  0, -1} 3
    35,  // { 1,  0,  0} 4
    21,  // { 0,  1,  0} 5
    14,  // { 0,  0,  1} 6
    28,  // { 0,  1,  1} 7
    42,  // { 1,  0,  1} 8
    49,  // { 1,  1,  0} 9
    12,  // { 0, -1,  1} 10
    10,  // {-1,  0,  1} 11
    17,  // {-1,  1,  0} 12
    20,  // { 0,  1, -1} 13
    34,  // { 1,  0, -1} 14
    33,  // { 1, -1,  0} 15
    4,   // { 0, -1, -1} 16
    2,   // {-1,  0, -1} 17
    1,   // {-1, -1,  0} 18
    56,  // { 1,  1,  1} 19
    24,  // {-1,  1,  1} 20
    40,  // { 1, -1,  1} 21
    48,  // { 1,  1, -1} 22
    32,  // { 1, -1, -1} 23
    16,  // {-1,  1, -1} 24
    8,   // {-1, -1,  1} 25
    0    // {-1, -1, -1} 26
  };

  // The point positions biased by lodNieghBias
  // todo(df): preserve this
  std::vector<point_t> biasedPos;
  biasedPos.reserve(packedVoxel.size());
  for (const auto& src : packedVoxel) {
    auto point = clacIntermediatePosition(
      aps.scalable_lifting_enabled_flag, lodIndex, src.position);
    biasedPos.push_back(times(point, aps.lodNeighBias));
  }

  atlas.reserve(retainedSize);
  std::vector<int32_t> neighborIndexes;
  neighborIndexes.reserve(64);

  BoxHierarchy<bucketSizeLog2, levelCount> hBBoxes;
  hBBoxes.resize(retainedSize);
  for (int32_t i = 0, b = 0; i < retainedSize; ++b) {
    hBBoxes.insert(biasedPos[retained[i]], i);
    ++i;
    for (int32_t k = 1; k < bucketSize && i < retainedSize; ++k, ++i) {
      hBBoxes.insert(biasedPos[retained[i]], i);
    }
  }
  hBBoxes.update();

  BoxHierarchy<bucketSizeLog2, levelCount> hIntraBBoxes;
  if (lodIndex >= aps.intra_lod_prediction_skip_layers) {
    hIntraBBoxes.resize(indexesSize);
    for (int32_t i = startIndex, b = 0; i < endIndex; ++b) {
      hIntraBBoxes.insert(biasedPos[indexes[i]], i - startIndex);
      ++i;
      for (int32_t k = 1; k < bucketSize && i < endIndex; ++k, ++i) {
        hIntraBBoxes.insert(biasedPos[indexes[i]], i - startIndex);
      }
    }
    hIntraBBoxes.update();
  }

  const auto bucketSize0Log2 = hBBoxes.bucketSizeLog2(0);
  const auto bucketSize1Log2 = hBBoxes.bucketSizeLog2(1);
  const auto bucketSize2Log2 = hBBoxes.bucketSizeLog2(2);

  int64_t curAtlasId = -1;
  int64_t lastMortonCodeShift3 = -1;
  int64_t cubeIndex = 0;
  for (int32_t i = startIndex, j = 0; i < endIndex; ++i) {
    int32_t localIndexes[3] = {-1, -1, -1};
    int64_t minDistances[3] = {std::numeric_limits<int64_t>::max(),
                               std::numeric_limits<int64_t>::max(),
                               std::numeric_limits<int64_t>::max()};

    const int32_t index = indexes[i];
    const auto& pv = packedVoxel[index];
    const int64_t mortonCode = pv.mortonCode;
    const int64_t pointAtlasId = mortonCode >> atlasBoundaryBit;
    const int64_t mortonCodeShiftBits3 = mortonCode >> shiftBits3;
    const int32_t pointIndex = pv.index;
    const auto bpoint = biasedPos[index];
    indexes[i] = pointIndex;
    auto& predictor = predictors[--predIndex];
    pointIndexToPredictorIndex[pointIndex] = predIndex;

    if (retainedSize) {
      while (j < retainedSize - 1
             && mortonCode >= packedVoxel[retained[j]].mortonCode) {
        ++j;
      }

      if (curAtlasId != pointAtlasId) {
        atlas.clearUpdates();
        curAtlasId = pointAtlasId;
        while (
          cubeIndex < retainedSize
          && (packedVoxel[retained[cubeIndex]].mortonCode >> atlasBoundaryBit)
            == curAtlasId) {
          atlas.set(
            packedVoxel[retained[cubeIndex]].mortonCode >> shiftBits3,
            cubeIndex);
          ++cubeIndex;
        }
      }

      if (lastMortonCodeShift3 != mortonCodeShiftBits3) {
        lastMortonCodeShift3 = mortonCodeShiftBits3;
        const auto basePosition = morton3dAdd(mortonCodeShiftBits3, -1ll);
        neighborIndexes.resize(0);
        for (int32_t n = 0; n < 27; ++n) {
          const auto neighbMortonCode =
            morton3dAdd(basePosition, kNeighOffset[n]);
          if ((neighbMortonCode >> atlasBits) != curAtlasId) {
            continue;
          }
          const auto range = atlas.get(neighbMortonCode);
          for (int32_t k = range.start; k < range.end; ++k) {
            neighborIndexes.push_back(k);
          }
        }
      }

      for (const auto k : neighborIndexes) {
        updateNearestNeigh(
          bpoint, biasedPos[retained[k]], k, localIndexes, minDistances);
      }

      if (localIndexes[2] == -1) {
        const auto center = localIndexes[0] == -1 ? j : localIndexes[0];
        const auto k0 = std::max(0, center - rangeInterLod);
        const auto k1 = std::min(retainedSize - 1, center + rangeInterLod);
        updateNearestNeighWithCheck(
          bpoint, biasedPos[retained[center]], center, localIndexes,
          minDistances);
        for (int32_t n = 1; n <= searchRangeNear; ++n) {
          const int32_t kp = center + n;
          if (kp <= k1) {
            updateNearestNeighWithCheck(
              bpoint, biasedPos[retained[kp]], kp, localIndexes, minDistances);
          }
          const int32_t kn = center - n;
          if (kn >= k0) {
            updateNearestNeighWithCheck(
              bpoint, biasedPos[retained[kn]], kn, localIndexes, minDistances);
          }
        }

        const int32_t p1 =
          std::min(retainedSize - 1, center + searchRangeNear + 1);
        const int32_t p0 = std::max(0, center - searchRangeNear - 1);

        // search p1...k1
        const int32_t b21 = k1 >> bucketSize2Log2;
        const int32_t b20 = p1 >> bucketSize2Log2;
        const int32_t b11 = k1 >> bucketSize1Log2;
        const int32_t b10 = p1 >> bucketSize1Log2;
        const int32_t b01 = k1 >> bucketSize0Log2;
        const int32_t b00 = p1 >> bucketSize0Log2;
        for (int32_t b2 = b20; b2 <= b21; ++b2) {
          if (
            localIndexes[2] != -1
            && hBBoxes.bBox(b2, 2).getDist1(bpoint) >= minDistances[2])
            continue;

          const auto alignedIndex1 = b2 << bucketSizeLog2;
          const auto start1 = std::max(b10, alignedIndex1);
          const auto end1 = std::min(b11, alignedIndex1 + bucketSizeMinus1);
          for (int32_t b1 = start1; b1 <= end1; ++b1) {
            if (
              localIndexes[2] != -1
              && hBBoxes.bBox(b1, 1).getDist1(bpoint) >= minDistances[2])
              continue;

            const auto alignedIndex0 = b1 << bucketSizeLog2;
            const auto start0 = std::max(b00, alignedIndex0);
            const auto end0 = std::min(b01, alignedIndex0 + bucketSizeMinus1);
            for (int32_t b0 = start0; b0 <= end0; ++b0) {
              if (
                localIndexes[2] != -1
                && hBBoxes.bBox(b0, 0).getDist1(bpoint) >= minDistances[2])
                continue;

              const int32_t alignedIndex = b0 << bucketSizeLog2;
              const int32_t h0 = std::max(p1, alignedIndex);
              const int32_t h1 = std::min(k1, alignedIndex + bucketSizeMinus1);
              for (int32_t k = h0; k <= h1; ++k) {
                updateNearestNeighWithCheck(
                  bpoint, biasedPos[retained[k]], k, localIndexes,
                  minDistances);
              }
            }
          }
        }

        // search k0...p1
        const int32_t c21 = p0 >> bucketSize2Log2;
        const int32_t c20 = k0 >> bucketSize2Log2;
        const int32_t c11 = p0 >> bucketSize1Log2;
        const int32_t c10 = k0 >> bucketSize1Log2;
        const int32_t c01 = p0 >> bucketSize0Log2;
        const int32_t c00 = k0 >> bucketSize0Log2;
        for (int32_t c2 = c21; c2 >= c20; --c2) {
          if (
            localIndexes[2] != -1
            && hBBoxes.bBox(c2, 2).getDist1(bpoint) >= minDistances[2])
            continue;

          const auto alignedIndex1 = c2 << bucketSizeLog2;
          const auto start1 = std::max(c10, alignedIndex1);
          const auto end1 = std::min(c11, alignedIndex1 + bucketSizeMinus1);
          for (int32_t c1 = end1; c1 >= start1; --c1) {
            if (
              localIndexes[2] != -1
              && hBBoxes.bBox(c1, 1).getDist1(bpoint) >= minDistances[2])
              continue;

            const auto alignedIndex0 = c1 << bucketSizeLog2;
            const auto start0 = std::max(c00, alignedIndex0);
            const auto end0 = std::min(c01, alignedIndex0 + bucketSizeMinus1);
            for (int32_t c0 = end0; c0 >= start0; --c0) {
              if (
                localIndexes[2] != -1
                && hBBoxes.bBox(c0, 0).getDist1(bpoint) >= minDistances[2])
                continue;

              const int32_t alignedIndex = c0 << bucketSizeLog2;
              const int32_t h0 = std::max(k0, alignedIndex);
              const int32_t h1 = std::min(p0, alignedIndex + bucketSizeMinus1);
              for (int32_t k = h1; k >= h0; --k) {
                updateNearestNeighWithCheck(
                  bpoint, biasedPos[retained[k]], k, localIndexes,
                  minDistances);
              }
            }
          }
        }
      }

      predictor.neighborCount = (localIndexes[0] != -1)
        + (localIndexes[1] != -1) + (localIndexes[2] != -1);

      for (int32_t h = 0; h < predictor.neighborCount; ++h)
        localIndexes[h] = retained[localIndexes[h]];
    }

    if (lodIndex >= aps.intra_lod_prediction_skip_layers) {
      const int32_t k00 = i + 1;
      const int32_t k01 = std::min(endIndex - 1, k00 + searchRangeNear);
      for (int32_t k = k00; k <= k01; ++k) {
        updateNearestNeigh(
          bpoint, biasedPos[indexes[k]], indexes[k], localIndexes,
          minDistances);
      }
      const int32_t k0 = k01 + 1 - startIndex;
      const int32_t k1 =
        std::min(endIndex - 1, k00 + rangeIntraLod) - startIndex;

      // search k0...k1
      const int32_t b21 = k1 >> bucketSize2Log2;
      const int32_t b20 = k0 >> bucketSize2Log2;
      const int32_t b11 = k1 >> bucketSize1Log2;
      const int32_t b10 = k0 >> bucketSize1Log2;
      const int32_t b01 = k1 >> bucketSize0Log2;
      const int32_t b00 = k0 >> bucketSize0Log2;
      for (int32_t b2 = b20; b2 <= b21; ++b2) {
        if (
          localIndexes[2] != -1
          && hIntraBBoxes.bBox(b2, 2).getDist1(bpoint) >= minDistances[2])
          continue;

        const auto alignedIndex1 = b2 << bucketSizeLog2;
        const auto start1 = std::max(b10, alignedIndex1);
        const auto end1 = std::min(b11, alignedIndex1 + bucketSizeMinus1);
        for (int32_t b1 = start1; b1 <= end1; ++b1) {
          if (
            localIndexes[2] != -1
            && hIntraBBoxes.bBox(b1, 1).getDist1(bpoint) >= minDistances[2])
            continue;

          const auto alignedIndex0 = b1 << bucketSizeLog2;
          const auto start0 = std::max(b00, alignedIndex0);
          const auto end0 = std::min(b01, alignedIndex0 + bucketSizeMinus1);
          for (int32_t b0 = start0; b0 <= end0; ++b0) {
            if (
              localIndexes[2] != -1
              && hIntraBBoxes.bBox(b0, 0).getDist1(bpoint) >= minDistances[2])
              continue;

            const int32_t alignedIndex = b0 << bucketSizeLog2;
            const int32_t h0 = std::max(k0, alignedIndex);
            const int32_t h1 = std::min(k1, alignedIndex + bucketSizeMinus1);
            for (int32_t h = h0; h <= h1; ++h) {
              const int32_t k = startIndex + h;
              updateNearestNeigh(
                bpoint, biasedPos[indexes[k]], indexes[k], localIndexes,
                minDistances);
            }
          }
        }
      }
    }

    predictor.neighborCount = std::min(
      aps.num_pred_nearest_neighbours_minus1 + 1,
      (localIndexes[0] != -1) + (localIndexes[1] != -1)
        + (localIndexes[2] != -1));
    for (int32_t h = 0; h < predictor.neighborCount; ++h) {
      auto& neigh = predictor.neighbors[h];
      neigh.predictorIndex = packedVoxel[localIndexes[h]].index;
      neigh.weight = (biasedPos[localIndexes[h]] - bpoint).getNorm2<int64_t>();
    }

    // Prune neighbours based upon max neigh range.
    if (aps.scalable_lifting_enabled_flag) {
      int maxNeighRange = aps.max_neigh_range_minus1 + 1;
      int64_t maxDistance = 3ll * maxNeighRange << 2 * lodIndex;
      if (aps.lodNeighBias == 1) {
        predictor.pruneDistanceGt(maxDistance);
      } else {
        auto curPt = clacIntermediatePosition(true, lodIndex, pv.position);

        for (int h = 1; h < predictor.neighborCount; h++) {
          auto neighPt = clacIntermediatePosition(
            true, lodIndex, packedVoxel[localIndexes[h]].position);

          // Discard this and subsequent points if distance limit exceeded
          auto norm2 = (curPt - neighPt).getNorm2<int64_t>();
          if (norm2 > maxDistance) {
            predictor.neighborCount = h;
            break;
          }
        }
      }
    }

    if (predictor.neighborCount > 1) {
      if (predictor.neighbors[0].weight > predictor.neighbors[1].weight)
        std::swap(predictor.neighbors[1], predictor.neighbors[0]);
      if (predictor.neighborCount == 3) {
        if (predictor.neighbors[1].weight > predictor.neighbors[2].weight) {
          std::swap(predictor.neighbors[2], predictor.neighbors[1]);
          if (predictor.neighbors[0].weight > predictor.neighbors[1].weight)
            std::swap(predictor.neighbors[1], predictor.neighbors[0]);
        }
      }
    }
  }
}

//---------------------------------------------------------------------------

inline void
subsampleByDistance(
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  const int32_t shiftBits0,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes,
  MortonIndexMap3d& atlas)
{
  assert(retained.empty());
  if (input.size() == 1) {
    indexes.push_back(input[0]);
    return;
  }

  const int64_t radius2 = 3ll << (shiftBits0 << 1);
  const int32_t shiftBits = shiftBits0 + 1;
  const int32_t shiftBits3 = 3 * shiftBits;
  const int32_t atlasBits = 3 * atlas.cubeSizeLog2();
  // NB: when the atlas boundary is greater than 2^63, all points belong
  //     to a single atlas.  The clipping is necessary to avoid undefined
  //     behaviour of shifts greater than or equal to the word size.
  const int32_t atlasBoundaryBit = std::min(63, shiftBits3 + atlasBits);

  // these neighbour offsets are relative to basePosition
  static const uint8_t kNeighOffset[20] = {
    7,   // { 0,  0,  0}
    3,   // {-1,  0,  0}
    5,   // { 0, -1,  0}
    6,   // { 0,  0, -1}
    12,  // { 0, -1,  1}
    10,  // {-1,  0,  1}
    17,  // {-1,  1,  0}
    20,  // { 0,  1, -1}
    34,  // { 1,  0, -1}
    33,  // { 1, -1,  0}
    4,   // { 0, -1, -1}
    2,   // {-1,  0, -1}
    1,   // {-1, -1,  0}
    24,  // {-1,  1,  1}
    40,  // { 1, -1,  1}
    48,  // { 1,  1, -1}
    32,  // { 1, -1, -1}
    16,  // {-1,  1, -1}
    8,   // {-1, -1,  1}
    0,   // {-1, -1, -1}
  };

  atlas.reserve(indexes.size() >> 1);
  int64_t curAtlasId = -1;
  int64_t lastRetainedMortonCode = -1;

  for (const auto index : input) {
    const auto& point = packedVoxel[index].position;
    const int64_t mortonCode = packedVoxel[index].mortonCode;
    const int64_t pointAtlasId = mortonCode >> atlasBoundaryBit;
    const int64_t mortonCodeShiftBits3 = mortonCode >> shiftBits3;

    if (curAtlasId != pointAtlasId) {
      atlas.clearUpdates();
      curAtlasId = pointAtlasId;
    }

    if (retained.empty()) {
      retained.push_back(index);
      lastRetainedMortonCode = mortonCodeShiftBits3;
      atlas.set(lastRetainedMortonCode, int32_t(retained.size()) - 1);
      continue;
    }

    if (lastRetainedMortonCode == mortonCodeShiftBits3) {
      indexes.push_back(index);
      continue;
    }

    // the position of the parent, offset by (-1,-1,-1)
    const auto basePosition = morton3dAdd(mortonCodeShiftBits3, -1ll);
    bool found = false;
    for (int32_t n = 0; n < 20 && !found; ++n) {
      const auto neighbMortonCode = morton3dAdd(basePosition, kNeighOffset[n]);
      if ((neighbMortonCode >> atlasBits) != curAtlasId)
        continue;

      const auto unit = atlas.get(neighbMortonCode);
      for (int32_t k = unit.start; k < unit.end; ++k) {
        const auto delta = (packedVoxel[retained[k]].position - point);
        if (delta.getNorm2<int64_t>() <= radius2) {
          found = true;
          break;
        }
      }
    }

    if (found) {
      indexes.push_back(index);
    } else {
      retained.push_back(index);
      lastRetainedMortonCode = mortonCodeShiftBits3;
      atlas.set(lastRetainedMortonCode, int32_t(retained.size()) - 1);
    }
  }
}

//---------------------------------------------------------------------------

inline int32_t
subsampleByOctreeWithCentroid(
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  int32_t octreeNodeSizeLog2,
  const bool backward,
  const std::vector<uint32_t>& voxels)
{
  point_t centroid(0);
  int count = 0;
  for (const auto t : voxels) {
    // forward direction
    point_t pos = clacIntermediatePosition(
      true, octreeNodeSizeLog2, pointCloud[packedVoxel[t].index]);

    centroid += pos;
    count++;
  }

  int32_t nnIndex = backward ? voxels.size() - 1 : 0;
  int64_t minNorm2 = std::numeric_limits<int64_t>::max();

  if (backward) {
    int num = voxels.size() - 1;
    for (auto t = voxels.rbegin(), e = voxels.rend(); t != e; t++) {
      // backward direction
      point_t pos = clacIntermediatePosition(
        true, octreeNodeSizeLog2, pointCloud[packedVoxel[*t].index]);
      pos *= count;
      int64_t m = (pos - centroid).getNorm1();
      if (minNorm2 > m) {
        minNorm2 = m;
        nnIndex = num;
      }
      num--;
    }
  } else {
    int num = 0;
    for (const auto t : voxels) {
      // forward direction
      point_t pos = clacIntermediatePosition(
        true, octreeNodeSizeLog2, pointCloud[packedVoxel[t].index]);
      pos *= count;
      int64_t m = (pos - centroid).getNorm1();
      if (minNorm2 > m) {
        minNorm2 = m;
        nnIndex = num;
      }
      num++;
    }
  }

  return voxels[nnIndex];
}

//---------------------------------------------------------------------------

inline void
subsampleByOctree(
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  int32_t octreeNodeSizeLog2,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes,
  bool direction,
  int lodSamplingPeriod = 0)
{
  const int indexCount = int(input.size());
  if (indexCount == 1) {
    indexes.push_back(input[0]);
    return;
  }

  uint64_t lodUniformQuant = 3 * (octreeNodeSizeLog2 + 1);
  uint64_t currVoxelPos;

  std::vector<uint32_t> voxels;
  voxels.reserve(8);

  for (int i = 0; i < indexCount; ++i) {
    uint64_t nextVoxelPos = currVoxelPos =
      (packedVoxel[input[i]].mortonCode >> lodUniformQuant);

    if (i < indexCount - 1)
      nextVoxelPos = (packedVoxel[input[i + 1]].mortonCode >> lodUniformQuant);

    voxels.push_back(input[i]);

    if (i == (indexCount - 1) || currVoxelPos < nextVoxelPos) {
      if ((voxels.size() < lodSamplingPeriod) && (i != (indexCount - 1)))
        continue;

      uint32_t picked = subsampleByOctreeWithCentroid(
        pointCloud, packedVoxel, octreeNodeSizeLog2, direction, voxels);

      for (const auto idx : voxels) {
        if (picked == idx)
          retained.push_back(idx);
        else
          indexes.push_back(idx);
      }
      voxels.clear();
    }
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
  const AttributeBrickHeader& abh,
  const PCCPointSet3& pointCloud,
  const std::vector<MortonCodeWithIndex>& packedVoxel,
  const std::vector<uint32_t>& input,
  const int32_t lodIndex,
  std::vector<uint32_t>& retained,
  std::vector<uint32_t>& indexes,
  MortonIndexMap3d& atlas)
{
  if (aps.scalable_lifting_enabled_flag) {
    int32_t octreeNodeSizeLog2 = lodIndex;
    bool direction = octreeNodeSizeLog2 & 1;
    subsampleByOctree(
      pointCloud, packedVoxel, input, octreeNodeSizeLog2, retained, indexes,
      direction);
  } else if (aps.lod_decimation_type == LodDecimationMethod::kPeriodic) {
    auto samplingPeriod = aps.lodSamplingPeriod[lodIndex];
    subsampleByDecimation(input, samplingPeriod, retained, indexes);
  } else if (aps.lod_decimation_type == LodDecimationMethod::kCentroid) {
    auto samplingPeriod = aps.lodSamplingPeriod[lodIndex];
    int32_t octreeNodeSizeLog2 = aps.dist2 + abh.attr_dist2_delta + lodIndex;
    subsampleByOctree(
      pointCloud, packedVoxel, input, octreeNodeSizeLog2, retained, indexes,
      true, samplingPeriod);
  } else {
    const auto shiftBits = aps.dist2 + abh.attr_dist2_delta + lodIndex;
    subsampleByDistance(
      packedVoxel, input, shiftBits, retained, indexes, atlas);
  }
}

//---------------------------------------------------------------------------

inline void
computeMortonCodesUnsorted(
  const PCCPointSet3& pointCloud,
  const Vec3<int32_t> lodNeighBias,
  std::vector<MortonCodeWithIndex>& packedVoxel)
{
  const int32_t pointCount = int32_t(pointCloud.getPointCount());
  packedVoxel.resize(pointCount);

  for (int n = 0; n < pointCount; n++) {
    auto& pv = packedVoxel[n];
    pv.position = pointCloud[n];
    pv.mortonCode = mortonAddr(pv.position);
    pv.index = n;
  }
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
  const AttributeBrickHeader& abh,
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
  computeMortonCodesUnsorted(pointCloud, aps.lodNeighBias, packedVoxel);

  if (!aps.canonical_point_order_flag)
    std::sort(packedVoxel.begin(), packedVoxel.end());

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

  bool concatenateLayers = aps.scalable_lifting_enabled_flag;
  std::vector<uint32_t> indexesOfSubsample;
  if (concatenateLayers)
    indexesOfSubsample.reserve(pointCount);

  std::vector<Box3<int32_t>> bBoxes;

  const int32_t log2CubeSize = 7;
  MortonIndexMap3d atlas;
  atlas.resize(log2CubeSize);
  atlas.init();

  auto maxNumDetailLevels = aps.maxNumDetailLevels();
  int32_t predIndex = int32_t(pointCount);
  for (auto lodIndex = minGeomNodeSizeLog2;
       !input.empty() && lodIndex < maxNumDetailLevels; ++lodIndex) {
    const int32_t startIndex = indexes.size();
    if (lodIndex == maxNumDetailLevels - 1) {
      for (const auto index : input) {
        indexes.push_back(index);
      }
    } else {
      subsample(
        aps, abh, pointCloud, packedVoxel, input, lodIndex, retained, indexes,
        atlas);
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
              aps, abh, packedVoxel, retained, divided_startIndex,
              divided_endIndex, lod + minGeomNodeSizeLog2, indexes, predictors,
              pointIndexToPredictorIndex, predIndex, atlas);
          }
        }
      }
    }

    computeNearestNeighbors(
      aps, abh, packedVoxel, retained, startIndex, endIndex, lodIndex, indexes,
      predictors, pointIndexToPredictorIndex, predIndex, atlas);

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
