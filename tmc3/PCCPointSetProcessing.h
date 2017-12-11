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

#ifndef PCCPointSetProcessing_h
#define PCCPointSetProcessing_h

#include <vector>

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "PCCKdTree.h"
#include "PCCPointSet.h"

namespace pcc {

inline bool PCCTransfertColors(const PCCPointSet3 &source, const int32_t searchRange,
                               PCCPointSet3 &target) {
  const size_t pointCountSource = source.getPointCount();
  const size_t pointCountTarget = target.getPointCount();
  if (!pointCountSource || !pointCountTarget || !source.hasColors()) {
    return false;
  }
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeTarget(3, target, 10);
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeSource(3, source, 10);

  target.addColors();
  std::vector<PCCColor3B> refinedColors1;
  std::vector<std::vector<PCCColor3B>> refinedColors2;
  refinedColors1.resize(pointCountTarget);
  refinedColors2.resize(pointCountTarget);
  const size_t num_results = 1;
  std::vector<size_t> indices(num_results);
  std::vector<double> sqrDist(num_results);
  nanoflann::KNNResultSet<double> resultSet(num_results);
  for (size_t index = 0; index < pointCountTarget; ++index) {
    resultSet.init(&indices[0], &sqrDist[0]);
    kdtreeSource.index->findNeighbors(resultSet, &target[index][0], nanoflann::SearchParams(10));
    refinedColors1[index] = source.getColor(indices[0]);
  }
  for (size_t index = 0; index < pointCountSource; ++index) {
    const PCCColor3B color = source.getColor(index);
    resultSet.init(&indices[0], &sqrDist[0]);
    kdtreeTarget.index->findNeighbors(resultSet, &source[index][0], nanoflann::SearchParams(10));
    refinedColors2[indices[0]].push_back(color);
  }
  for (size_t index = 0; index < pointCountTarget; ++index) {
    const PCCColor3B color1 = refinedColors1[index];
    const std::vector<PCCColor3B> &colors2 = refinedColors2[index];
    if (colors2.empty()) {
      target.setColor(index, color1);
    } else {
      const double H = double(colors2.size());
      const PCCVector3D centroid1(color1[0], color1[1], color1[2]);
      PCCVector3D centroid2(0.0);
      for (const auto color2 : colors2) {
        for (size_t k = 0; k < 3; ++k) {
          centroid2[k] += color2[k];
        }
      }
      centroid2 /= H;

      double D2 = 0.0;
      for (const auto color2 : colors2) {
        for (size_t k = 0; k < 3; ++k) {
          const double d2 = centroid2[k] - color2[k];
          D2 += d2 * d2;
        }
      }
      const double r = double(pointCountTarget) / double(pointCountSource);
      const double delta2 = (centroid2 - centroid1).getNorm2();
      const double eps = 0.000001;
      if (delta2 > eps) {  // centroid2 != centroid1
        double w = 0.0;
        const double alpha = D2 / delta2;
        const double a = H * r - 1.0;
        const double c = alpha * r - 1.0;
        if (fabs(a) < eps) {
          w = -0.5 * c;
        } else {
          const double delta = 1.0 - a * c;
          if (delta >= 0.0) {
            w = (-1.0 + sqrt(delta)) / a;
          }
        }
        const double oneMinusW = 1.0 - w;
        PCCVector3D color0;
        for (size_t k = 0; k < 3; ++k) {
          color0[k] = PCCClip(round(w * centroid1[k] + oneMinusW * centroid2[k]), 0.0, 255.0);
        }
        const double rSource = 1.0 / double(pointCountSource);
        const double rTarget = 1.0 / double(pointCountTarget);
        const double maxValue = std::numeric_limits<uint8_t>::max();
        double minError = std::numeric_limits<double>::max();
        PCCVector3D bestColor(color0);
        PCCVector3D color;
        for (int32_t s1 = -searchRange; s1 <= searchRange; ++s1) {
          color[0] = PCCClip(color0[0] + s1, 0.0, maxValue);
          for (int32_t s2 = -searchRange; s2 <= searchRange; ++s2) {
            color[1] = PCCClip(color0[1] + s2, 0.0, maxValue);
            for (int32_t s3 = -searchRange; s3 <= searchRange; ++s3) {
              color[2] = PCCClip(color0[2] + s3, 0.0, maxValue);

              double e1 = 0.0;
              for (size_t k = 0; k < 3; ++k) {
                const double d = color[k] - color1[k];
                e1 += d * d;
              }
              e1 *= rTarget;

              double e2 = 0.0;
              for (const auto color2 : colors2) {
                for (size_t k = 0; k < 3; ++k) {
                  const double d = color[k] - color2[k];
                  e2 += d * d;
                }
              }
              e2 *= rSource;

              const double error = std::max(e1, e2);
              if (error < minError) {
                minError = error;
                bestColor = color;
              }
            }
          }
        }
        target.setColor(
            index, PCCColor3B(uint8_t(bestColor[0]), uint8_t(bestColor[1]), uint8_t(bestColor[2])));
      } else {  // centroid2 == centroid1
        target.setColor(index, color1);
      }
    }
  }
  return true;
}

inline bool PCCTransfertReflectances(const PCCPointSet3 &source, const int32_t searchRange,
                                     PCCPointSet3 &target) {
  const size_t pointCountSource = source.getPointCount();
  const size_t pointCountTarget = target.getPointCount();
  if (!pointCountSource || !pointCountTarget || !source.hasReflectances()) {
    return false;
  }
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeTarget(3, target, 10);
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeSource(3, source, 10);

  target.addReflectances();
  std::vector<uint16_t> refined1;
  std::vector<std::vector<uint16_t>> refined2;
  refined1.resize(pointCountTarget);
  refined2.resize(pointCountTarget);
  const size_t num_results = 1;
  std::vector<size_t> indices(num_results);
  std::vector<double> sqrDist(num_results);
  nanoflann::KNNResultSet<double> resultSet(num_results);
  for (size_t index = 0; index < pointCountTarget; ++index) {
    resultSet.init(&indices[0], &sqrDist[0]);
    kdtreeSource.index->findNeighbors(resultSet, &target[index][0], nanoflann::SearchParams(10));
    refined1[index] = source.getReflectance(indices[0]);
  }
  for (size_t index = 0; index < pointCountSource; ++index) {
    const uint16_t reflectance = source.getReflectance(index);
    resultSet.init(&indices[0], &sqrDist[0]);
    kdtreeTarget.index->findNeighbors(resultSet, &source[index][0], nanoflann::SearchParams(10));
    refined2[indices[0]].push_back(reflectance);
  }
  for (size_t index = 0; index < pointCountTarget; ++index) {
    const uint16_t reflectance1 = refined1[index];
    const std::vector<uint16_t> &reflectances2 = refined2[index];
    if (reflectances2.empty()) {
      target.setReflectance(index, reflectance1);
    } else {
      const double H = double(reflectances2.size());
      const double centroid1 = reflectance1;
      double centroid2 = 0.0;
      for (const auto reflectance2 : reflectances2) {
        centroid2 += reflectance2;
      }
      centroid2 /= H;

      double D2 = 0.0;
      for (const auto reflectance2 : reflectances2) {
        const double d2 = centroid2 - reflectance2;
        D2 += d2 * d2;
      }

      const double delta2 = pow(centroid2 - centroid1, 2.0);
      const double eps = 0.000001;

      if (delta2 > eps) {  // centroid2 != centroid1
        double w = 0.0;
        const double alpha = D2 / delta2;
        const double r = double(pointCountTarget) / double(pointCountSource);
        const double a = H * r - 1.0;
        const double c = alpha * r - 1.0;
        if (fabs(a) < eps) {
          w = -0.5 * c;
        } else {
          const double delta = 1.0 - a * c;
          if (delta >= 0.0) {
            w = (-1.0 + sqrt(delta)) / a;
          }
        }
        const double oneMinusW = 1.0 - w;
        const double reflectance0 =
            PCCClip(round(w * centroid1 + oneMinusW * centroid2), 0.0, 255.0);
        const double rSource = 1.0 / double(pointCountSource);
        const double rTarget = 1.0 / double(pointCountTarget);
        double minError = std::numeric_limits<double>::max();
        double bestReflectance = reflectance0;
        double reflectance;
        const double maxValue = std::numeric_limits<uint16_t>::max();
        for (int32_t s = -searchRange; s <= searchRange; ++s) {
          reflectance = PCCClip(reflectance0 + s, 0.0, maxValue);
          const double d = reflectance - reflectance1;
          const double e1 = (d * d) * rTarget;

          double e2 = 0.0;
          for (const auto reflectance2 : reflectances2) {
            const double d = reflectance - reflectance2;
            e2 += d * d;
          }
          e2 *= rSource;

          const double error = std::max(e1, e2);
          if (error < minError) {
            minError = error;
            bestReflectance = reflectance;
          }
        }
        target.setReflectance(index, uint16_t(bestReflectance));
      } else {  // centroid2 == centroid1
        target.setReflectance(index, reflectance1);
      }
    }
  }
  return true;
}
};

#endif /* PCCPointSetProcessing_h */
