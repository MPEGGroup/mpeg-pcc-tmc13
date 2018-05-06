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

#include <vector>

namespace pcc {
const uint32_t PCCTMC3MagicNumber = 20110904;
const uint32_t PCCTMC3FormatVersion = 1;
const uint32_t PCCTMC3MaxPredictionNearestNeighborCount = 8;
const uint32_t PCCTMC3Diff1AdaptiveDataModelCount = 512;
const uint32_t PCCTMC3AdaptiveDataModelAlphabetMaxSize = 2047;

struct PCCOctree3Node {
  PCCBox3<uint32_t> boundingBox;
  size_t start;
  size_t end;
};

struct PCCNeighborInfo {
  double invDistance;
  double weight;
  uint32_t index;
};

struct PCCPredictor {
  size_t maxNeighborCount;
  size_t neighborCount;
  uint32_t index;
  uint32_t levelOfDetailIndex;
  PCCNeighborInfo neighbors[PCCTMC3MaxPredictionNearestNeighborCount];

  PCCColor3B predictColor(const PCCPointSet3 &pointCloud) const {
    PCCVector3D predicted(0.0);
    for (size_t i = 0; i < neighborCount; ++i) {
      const PCCColor3B color = pointCloud.getColor(neighbors[i].index);
      const double w = neighbors[i].weight;
      for (size_t k = 0; k < 3; ++k) {
        predicted[k] += w * color[k];
      }
    }
    return PCCColor3B(uint8_t(std::round(predicted[0])), uint8_t(std::round(predicted[1])),
                      uint8_t(std::round(predicted[2])));
  }
  uint16_t predictReflectance(const PCCPointSet3 &pointCloud) const {
    double predicted(0.0);
    for (size_t i = 0; i < neighborCount; ++i) {
      predicted += neighbors[i].weight * pointCloud.getReflectance(neighbors[i].index);
    }
    return uint16_t(std::round(predicted));
  }

  void computeWeights(const size_t neighborCount0) {
    neighborCount = std::min(neighborCount0, maxNeighborCount);
    double sum = 0.0;
    for (size_t n = 0; n < neighborCount; ++n) {
      neighbors[n].weight = neighbors[n].invDistance;
      sum += neighbors[n].weight;
    }
    assert(sum > 0.0);
    for (size_t n = 0; n < neighborCount; ++n) {
      neighbors[n].weight /= sum;
    }
  }

  void init(const uint32_t current, const uint32_t reference, const uint32_t lodIndex) {
    index = current;
    neighborCount = (reference != PCC_UNDEFINED_INDEX);
    maxNeighborCount = neighborCount;
    neighbors[0].index = reference;
    neighbors[0].invDistance = 1.0;
    neighbors[0].weight = 1.0;
    levelOfDetailIndex = lodIndex;
  }
};
inline int64_t PCCQuantization(const int64_t value, const int64_t qs, const int64_t deadZone) {
  if (!qs) {
    return value;
  } else {
    const int64_t qs2 = qs / 2;
    if (value > deadZone) {
      return (value - deadZone + qs2) / qs + 1;
    } else if (value < -deadZone) {
      return (value + deadZone - qs2) / qs - 1;
    } else {
      return int64_t(0);
    }
  }
}
inline int64_t PCCInverseQuantization(const int64_t value, const int64_t qs, int64_t deadZone) {
  if (!qs) {
    return value;
  } else {
    if (value > 0) {
      return (value - 1) * qs + deadZone;
    } else if (value < 0) {
      return (value + 1) * qs - deadZone;
    } else {
      return int64_t(0);
    }
  }
}
inline void PCCBuildPredictors(const PCCPointSet3 &pointCloud,
                               const size_t numberOfNearestNeighborsInPrediction,
                               const size_t levelOfDetailCount, const std::vector<size_t> &dist2,
                               std::vector<PCCPredictor> &predictors,
                               std::vector<uint32_t> &numberOfPointsPerLOD,
                               std::vector<uint32_t> &indexes) {
  const size_t pointCount = pointCloud.getPointCount();
  predictors.resize(pointCount);
  const size_t initialBalanceTargetPointCount = 16;
  std::vector<bool> visited;
  visited.resize(pointCount, false);
  indexes.resize(0);
  indexes.reserve(pointCount);
  numberOfPointsPerLOD.resize(0);
  numberOfPointsPerLOD.reserve(levelOfDetailCount);
  size_t prevIndexesSize = 0;
  PCCPointDistInfo nearestNeighbors[PCCTMC3MaxPredictionNearestNeighborCount];
  PCCNNQuery3 nNQuery = {PCCVector3D(0.0), 0.0, 1};
  PCCNNQuery3 nNQuery2 = {PCCVector3D(0.0), std::numeric_limits<double>::max(),
                          numberOfNearestNeighborsInPrediction};
  PCCNNResult nNResult = {nearestNeighbors, 0};
  PCCIncrementalKdTree3 kdtree;

  for (size_t lodIndex = 0; lodIndex < levelOfDetailCount && indexes.size() < pointCount;
       ++lodIndex) {
    const bool filterByDistance = (lodIndex + 1) != levelOfDetailCount;
    const double minDist2 = dist2[lodIndex];
    nNQuery.radius = minDist2;
    size_t balanceTargetPointCount = initialBalanceTargetPointCount;
    for (uint32_t current = 0; current < pointCount; ++current) {
      if (!visited[current]) {
        if (filterByDistance) {
          nNQuery.point = pointCloud[current];
          kdtree.findNearestNeighbors2(nNQuery, nNResult);
        }
        if (!filterByDistance || !nNResult.resultCount || nNResult.neighbors[0].dist2 >= minDist2) {
          nNQuery2.point = pointCloud[current];
          kdtree.findNearestNeighbors2(nNQuery2, nNResult);
          auto &predictor = predictors[indexes.size()];
          if (!nNResult.resultCount) {
            predictor.init(current, PCC_UNDEFINED_INDEX, uint32_t(lodIndex));
          } else if (nNResult.neighbors[0].dist2 == 0.0 || nNResult.resultCount == 1) {
            const uint32_t reference = static_cast<uint32_t>(indexes[nNResult.neighbors[0].index]);
            predictor.init(current, reference, uint32_t(lodIndex));
          } else {
            predictor.levelOfDetailIndex = uint32_t(lodIndex);
            predictor.index = current;
            predictor.neighborCount = predictor.maxNeighborCount = nNResult.resultCount;
            for (size_t n = 0; n < nNResult.resultCount; ++n) {
              predictor.neighbors[n].index = indexes[nNResult.neighbors[n].index];
              predictor.neighbors[n].weight = predictor.neighbors[n].invDistance =
                  1.0 / nNResult.neighbors[n].dist2;
            }
          }
          indexes.push_back(current);
          visited[current] = true;
          kdtree.insert(pointCloud[current]);
          if (balanceTargetPointCount <= kdtree.size()) {
            kdtree.balance();
            balanceTargetPointCount = 2 * kdtree.size();
          }
        }
      }
    }
    numberOfPointsPerLOD.push_back(uint32_t(indexes.size()));
    prevIndexesSize = indexes.size();
  }
}
}

#endif /* PCCTMC3Common_h */
