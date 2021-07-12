/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2020, ISO/IEC
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

#pragma once

#include <cmath>
#include <cstdint>

#include "entropy.h"
#include "PCCMath.h"
#include "hls.h"

namespace pcc {

//============================================================================

struct GPredicter {
  enum Mode
  {
    None,
    Delta,
    Linear2,
    Linear3
  };

  int32_t index[3];
  int pgeom_min_radius;

  bool isValid(Mode mode);

  Vec3<int32_t> predict(const Vec3<int32_t>* points, Mode mode, bool angular);
};

//============================================================================

struct GNode {
  static const int32_t MaxChildrenCount = 3;

  int numDups;
  int32_t parent;
  int32_t childrenCount;
  int32_t children[MaxChildrenCount];
};

//============================================================================

class PredGeomContexts {
public:
  void reset();

protected:
  AdaptiveBitModel _ctxNumChildren[3];
  AdaptiveBitModel _ctxPredMode[3];
  AdaptiveBitModel _ctxResGt0[3];
  AdaptiveBitModel _ctxSign[3];
  AdaptiveBitModel _ctxNumBits[5][3][31];
  AdaptiveBitModel _ctxNumDupPointsGt0;
  AdaptiveBitModel _ctxNumDupPoints;

  AdaptiveBitModel _ctxResidual2GtN[2][3];
  AdaptiveBitModel _ctxSign2[3];
  AdaptiveBitModel _ctxEG2Prefix[3][5];
  AdaptiveBitModel _ctxEG2Suffix[3][4];

  AdaptiveBitModel _ctxQpOffsetAbsGt0;
  AdaptiveBitModel _ctxQpOffsetSign;
  AdaptiveBitModel _ctxQpOffsetAbsEgl;

  AdaptiveBitModel _ctxPhiGtN[2];
  AdaptiveBitModel _ctxSignPhi;
  AdaptiveBitModel _ctxEGPhi;
  AdaptiveBitModel _ctxResidualPhi[7];

  AdaptiveBitModel _ctxEndOfTrees;
};

//----------------------------------------------------------------------------

inline void
PredGeomContexts::reset()
{
  this->~PredGeomContexts();
  new (this) PredGeomContexts;
}

//============================================================================

template<typename LookupFn>
GPredicter
makePredicter(
  int32_t curNodeIdx,
  GPredicter::Mode mode,
  int minRadius,
  LookupFn nodeIdxToParentIdx)
{
  if (mode == GPredicter::None)
    mode = GPredicter::Delta;

  GPredicter predIdx;
  predIdx.pgeom_min_radius = minRadius;

  switch (mode) {
  default:
  case GPredicter::None:
  case GPredicter::Delta:
  case GPredicter::Linear2:
  case GPredicter::Linear3:
    for (int i = 0; i < int(mode); i++) {
      if (curNodeIdx < 0)
        break;
      predIdx.index[i] = curNodeIdx = nodeIdxToParentIdx(curNodeIdx);
    }
    break;
  }
  return predIdx;
}

//============================================================================

inline bool
GPredicter::isValid(GPredicter::Mode mode)
{
  int numPredictors = int(mode);
  for (int i = 0; i < numPredictors; i++) {
    if (this->index[i] < 0)
      return false;
  }
  return true;
}

//============================================================================

inline Vec3<int32_t>
GPredicter::predict(
  const Vec3<int32_t>* points, GPredicter::Mode mode, bool angular)
{
  Vec3<int32_t> pred;
  switch (mode) {
  case GPredicter::None: {
    pred = 0;

    if (angular)
      pred[0] = pgeom_min_radius;

    if (this->index[0] >= 0 && angular) {
      const auto& p0 = points[this->index[0]];
      pred[1] = p0[1];
      pred[2] = p0[2];
    }
    break;
  }

  case GPredicter::Delta: {
    pred = points[this->index[0]];
    break;
  }

  case GPredicter::Linear2: {
    const auto& p0 = points[this->index[0]];
    const auto& p1 = points[this->index[1]];
    pred = 2 * p0 - p1;
    break;
  }

  default:
  case GPredicter::Linear3: {
    const auto& p0 = points[this->index[0]];
    const auto& p1 = points[this->index[1]];
    const auto& p2 = points[this->index[2]];
    pred = p0 + p1 - p2;
    break;
  }
  }
  return pred;
}

//============================================================================

class SphericalToCartesian {
public:
  SphericalToCartesian(const GeometryParameterSet& gps)
    : log2ScaleRadius(gps.geom_angular_radius_inv_scale_log2)
    , log2ScalePhi(gps.geom_angular_azimuth_scale_log2_minus11 + 12)
    , tanThetaLaser(gps.angularTheta.data())
    , zLaser(gps.angularZ.data())
  {}

  Vec3<int32_t> operator()(Vec3<int32_t> sph)
  {
    int64_t r = sph[0] << log2ScaleRadius;
    int64_t z = divExp2RoundHalfInf(
      tanThetaLaser[sph[2]] * r << 2, log2ScaleTheta - log2ScaleZ);

    return Vec3<int32_t>(Vec3<int64_t>{
      divExp2RoundHalfInf(r * icos(sph[1], log2ScalePhi), kLog2ISineScale),
      divExp2RoundHalfInf(r * isin(sph[1], log2ScalePhi), kLog2ISineScale),
      divExp2RoundHalfInf(z - zLaser[sph[2]], log2ScaleZ)});
  }

private:
  static constexpr int log2ScaleZ = 3;
  static constexpr int log2ScaleTheta = 20;
  int log2ScaleRadius;
  int log2ScalePhi;
  const int* tanThetaLaser;
  const int* zLaser;
};

//============================================================================

class CartesianToSpherical {
public:
  CartesianToSpherical(const GeometryParameterSet& gps)
    : sphToCartesian(gps)
    , log2ScaleRadius(gps.geom_angular_radius_inv_scale_log2)
    , scalePhi(1 << (gps.geom_angular_azimuth_scale_log2_minus11 + 12))
    , numLasers(gps.angularTheta.size())
    , tanThetaLaser(gps.angularTheta.data())
    , zLaser(gps.angularZ.data())
  {}

  Vec3<int32_t> operator()(Vec3<int32_t> xyz)
  {
    int64_t r0 = int64_t(std::round(hypot(xyz[0], xyz[1])));
    int32_t thetaIdx = 0;
    int32_t minError = std::numeric_limits<int32_t>::max();
    for (int idx = 0; idx < numLasers; ++idx) {
      int64_t z = divExp2RoundHalfInf(
        tanThetaLaser[idx] * r0 << 2, log2ScaleTheta - log2ScaleZ);
      int64_t z1 = divExp2RoundHalfInf(z - zLaser[idx], log2ScaleZ);
      int32_t err = abs(z1 - xyz[2]);
      if (err < minError) {
        thetaIdx = idx;
        minError = err;
      }
    }

    auto phi0 = std::round((atan2(xyz[1], xyz[0]) / (2.0 * M_PI)) * scalePhi);

    Vec3<int32_t> sphPos{int32_t(divExp2RoundHalfUp(r0, log2ScaleRadius)),
                         int32_t(phi0), thetaIdx};

    // local optmization
    auto minErr = (sphToCartesian(sphPos) - xyz).getNorm1();
    int32_t dt0 = 0;
    int32_t dr0 = 0;
    for (int32_t dt = -2; dt <= 2 && minErr; ++dt) {
      for (int32_t dr = -2; dr <= 2; ++dr) {
        auto sphPosCand = sphPos + Vec3<int32_t>{dr, dt, 0};
        auto err = (sphToCartesian(sphPosCand) - xyz).getNorm1();
        if (err < minErr) {
          minErr = err;
          dt0 = dt;
          dr0 = dr;
        }
      }
    }
    sphPos[0] += dr0;
    sphPos[1] += dt0;

    return sphPos;
  }

private:
  SphericalToCartesian sphToCartesian;
  static constexpr int32_t log2ScaleZ = 3;
  static constexpr int32_t log2ScaleTheta = 20;
  int32_t log2ScaleRadius;
  int32_t scalePhi;
  int numLasers;
  const int* tanThetaLaser;
  const int* zLaser;
};

//============================================================================

}  // namespace pcc
