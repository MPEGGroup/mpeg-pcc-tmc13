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

#ifndef PCCPointSetProcessing_h
#define PCCPointSetProcessing_h

#include <vector>

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "PCCKdTree.h"
#include "PCCPointSet.h"

namespace pcc {

inline bool
PCCTransfertColors(const PCCPointSet3& source, PCCPointSet3& target)
{
  const size_t pointCountSource = source.getPointCount();
  const size_t pointCountTarget = target.getPointCount();
  if (!pointCountSource || !pointCountTarget || !source.hasColors()) {
    return false;
  }
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeTarget(
    3, target, 10);
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeSource(
    3, source, 10);

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
    kdtreeSource.index->findNeighbors(
      resultSet, &target[index][0], nanoflann::SearchParams(10));
    refinedColors1[index] = source.getColor(indices[0]);
  }
  for (size_t index = 0; index < pointCountSource; ++index) {
    const PCCColor3B color = source.getColor(index);
    resultSet.init(&indices[0], &sqrDist[0]);
    kdtreeTarget.index->findNeighbors(
      resultSet, &source[index][0], nanoflann::SearchParams(10));
    refinedColors2[indices[0]].push_back(color);
  }
  for (size_t index = 0; index < pointCountTarget; ++index) {
    const PCCColor3B color1 = refinedColors1[index];
    const std::vector<PCCColor3B>& colors2 = refinedColors2[index];
    if (colors2.empty()) {
      target.setColor(index, color1);
    } else {
      PCCVector3D centroid2(0.0);
      for (const auto& color2 : colors2) {
        for (size_t k = 0; k < 3; ++k) {
          centroid2[k] += color2[k];
        }
      }
      centroid2 /= colors2.size();
      PCCColor3B color0;
      for (size_t k = 0; k < 3; ++k) {
        color0[k] = uint8_t(PCCClip(round(centroid2[k]), 0.0, 255.0));
      }
      target.setColor(index, color0);
    }
  }
  return true;
}

inline bool
PCCTransfertReflectances(const PCCPointSet3& source, PCCPointSet3& target)
{
  const size_t pointCountSource = source.getPointCount();
  const size_t pointCountTarget = target.getPointCount();
  if (!pointCountSource || !pointCountTarget || !source.hasReflectances()) {
    return false;
  }
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeTarget(
    3, target, 10);
  KDTreeVectorOfVectorsAdaptor<PCCPointSet3, double> kdtreeSource(
    3, source, 10);

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
    kdtreeSource.index->findNeighbors(
      resultSet, &target[index][0], nanoflann::SearchParams(10));
    refined1[index] = source.getReflectance(indices[0]);
  }
  for (size_t index = 0; index < pointCountSource; ++index) {
    const uint16_t reflectance = source.getReflectance(index);
    resultSet.init(&indices[0], &sqrDist[0]);
    kdtreeTarget.index->findNeighbors(
      resultSet, &source[index][0], nanoflann::SearchParams(10));
    refined2[indices[0]].push_back(reflectance);
  }
  for (size_t index = 0; index < pointCountTarget; ++index) {
    const uint16_t reflectance1 = refined1[index];
    const std::vector<uint16_t>& reflectances2 = refined2[index];
    if (reflectances2.empty()) {
      target.setReflectance(index, reflectance1);
    } else {
      double centroid2 = 0.0;
      for (const auto reflectance2 : reflectances2) {
        centroid2 += reflectance2;
      }
      centroid2 = PCCClip(
        std::round(centroid2 / reflectances2.size()), 0.0,
        double(std::numeric_limits<uint16_t>::max()));
      target.setReflectance(index, uint16_t(centroid2));
    }
  }
  return true;
}
};  // namespace pcc

#endif /* PCCPointSetProcessing_h */
