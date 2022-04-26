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
#include "PCCPointSet.h"
#include "geometry_octree.h"
#include <map>


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
  static const uint8_t interFlagBufferMask = 0x1F;

protected:
  static const int NPredDelta = 4;
  AdaptiveBitModel _ctxNumChildren[3];
  AdaptiveBitModel _ctxPredMode[3];
  AdaptiveBitModel _ctxPredIdx[NPredDelta];
  AdaptiveBitModel _ctxResGt0[2][3];
  AdaptiveBitModel _ctxSign[2][3];
  AdaptiveBitModel _ctxNumBits[2][5][3][31];
  AdaptiveBitModel _ctxNumDupPointsGt0;
  AdaptiveBitModel _ctxNumDupPoints;
  AdaptiveBitModel _ctxInterFlag[32];
  AdaptiveBitModel _ctxRefNodeFlag;

  AdaptiveBitModel _ctxResidual2GtN[2][3];
  AdaptiveBitModel _ctxSign2[3];
  AdaptiveBitModel _ctxEG2Prefix[3][5];
  AdaptiveBitModel _ctxEG2Suffix[3][4];

  AdaptiveBitModel _ctxQpOffsetAbsGt0;
  AdaptiveBitModel _ctxQpOffsetSign;
  AdaptiveBitModel _ctxQpOffsetAbsEgl;

  AdaptiveBitModel _ctxPhiGtN[2][2];
  AdaptiveBitModel _ctxSignPhi[2];
  AdaptiveBitModel _ctxEGPhi[2];
  AdaptiveBitModel _ctxResidualPhi[2][7];

  AdaptiveBitModel _ctxEndOfTrees;
  AdaptiveBitModel _ctxResRIsZero[2][4];
  AdaptiveBitModel _ctxResRIsOne[2][4];
  AdaptiveBitModel _ctxResRIsTwo[2][4];
  AdaptiveBitModel _ctxResRExpGolombPre[2][4][10];
  AdaptiveBitModel _ctxResRExpGolombSuf[2][4][10];

  AdaptiveBitModel _ctxResPhiIsZero[2][2];
  int _resPhiOldSign = 2;
  AdaptiveBitModel _ctxResPhiSign[2][2+2];
  AdaptiveBitModel _ctxResPhiIsOne[2][2];
  AdaptiveBitModel _ctxResPhiExpGolombPre[2][2][4];
  AdaptiveBitModel _ctxResPhiExpGolombSuf[2][2][4];

  bool prevInterFlag = false;
  int              _precAzimuthStepDelta = 0;
  int              _precSignR = 0;

  AdaptiveBitModel _ctxResRSign[3][2][8];
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
    pred = 0;

    pred[0] = pgeom_min_radius;

    if (this->index[0] >= 0) {
      pred = points[this->index[0]];
    }
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
class CartesianToSphericalSimple {
public:
  CartesianToSphericalSimple(const GeometryParameterSet& gps)
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

    Vec3<int32_t> sphPos{
      int32_t(divExp2RoundHalfUp(r0, log2ScaleRadius)), int32_t(phi0),
      thetaIdx};    

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
typedef std::vector<Vec3<int>> VecVec3;
typedef std::vector<VecVec3> VecVecVec3;
typedef std::vector<VecVecVec3> VecVecVecVec3;
inline int64_t
divExp2RoundHalfInfPositiveShift(
  const int64_t scalar, const unsigned int shift, const unsigned int s0)
{
  return scalar >= 0 ? (s0 + scalar) >> shift : -((s0 - scalar) >> shift);
}
class MotionParameters {
private:
  int numFrames;
  std::vector<std::vector<int>> motionMatrix;
  std::vector<Vec3<int>> transVec;
  std::vector<std::pair<int, int>> threshVec;
  int frameCounter;

public:
  MotionParameters() : numFrames(0), frameCounter(-1) {}
  bool EndOfFrames() const { return frameCounter >= numFrames; }
  void AdvanceFrame() { frameCounter++; }
  int getFrameCtr() { return frameCounter; }
  float quantizeParameter(float x)
  {
    int scaleFactor = 65536;
    float y = x * scaleFactor;
    if (y < 0)
      return int(y - 0.5) / float(scaleFactor);
    else
      return int(y + 0.5) / float(scaleFactor);
  }
  void parseFile(std::string fileName)
  {
    if (numFrames)
      return;

    std::ifstream fin(fileName);
    float tmp;
    std::vector<float> out;
    std::copy(
      std::istream_iterator<float>(fin), std::istream_iterator<float>(),
      std::back_inserter(out));
    numFrames = out.size() / 14;
    motionMatrix.resize(numFrames);
    transVec.resize(numFrames);
    threshVec.resize(numFrames);
    const int scaleFactor = 65536;
    auto it = out.begin();
    for (auto i = 0; i < numFrames; i++) {
      motionMatrix[i].resize(9);
      for (int j = 0; j < 9; j++)
        if (j % 3 == j / 3)
          motionMatrix[i][j] = std::round((*(it++) - 1) * scaleFactor) + 65536;
        else
          motionMatrix[i][j] = std::round((*(it++)) * scaleFactor);
      for (int j = 0; j < 3; j++)
        transVec[i][j] = std::round(*(it++));
      threshVec[i].first = std::round(*(it++));  // quantizeParameter(*(it++));
      threshVec[i].second = std::round(*(it++)); // quantizeParameter(*(it++));
    }
  }
  void
    setMotionParams(std::pair<int, int> thresh, std::vector<int> matrix, Vec3<int> trans)
  {
    motionMatrix.push_back(matrix);
    threshVec.push_back(thresh);
    transVec.push_back(trans);
  }
  template <typename T>
  void 
  getMotionParams(
    std::pair<int, int> &th, std::vector<T> &mat, Vec3<int> &tr) const
  {
    if (frameCounter > threshVec.size()) {
      std::cout << "Accessing unassigned values\n";
    }
    th = threshVec[frameCounter];
    mat.resize(9);
    for (auto i = 0; i < 9; i++)
      mat[i] = motionMatrix[frameCounter][i];
    tr = transVec[frameCounter];
  }
};

class PredGeomPredictor {
public:
  PredGeomPredictor() : numLasers(0), azimScaleLog2(0), interEnabled(false) {}

  void init(int azim_scale_log2, int nLasers, const bool globMotionEnabled)
  {
    if (!numLasers) {
      azimScaleLog2 = azim_scale_log2;
      numLasers = nLasers;
      globalMotionEnabled = globMotionEnabled;
      refPointVals.resize(nLasers);
      refPointValsGlob.resize(nLasers);
      refPointValsCur.resize(nLasers);
    } else {
    }  // Already initialized
  }

  void insert(const std::vector<Vec3<int32_t>>& refFramePosSph)
  {
    assert(numLasers > 0);
    for (auto const& pt : refFramePosSph) {
      refPointValsCur[pt[2]].insert({computePhiQuantized(pt[1]), pt});
    }
  }

  std::pair<bool, point_t> getClosestPred(int currAzim, int currLaserId)
  {
    const auto& refPoints = refPointVals;
    auto quantizedPhi = computePhiQuantized(currAzim);
    if (refPoints[currLaserId].size()) {
      auto idx = refPoints[currLaserId].upper_bound(quantizedPhi);
      if (idx == refPoints[currLaserId].end())
        return std::pair<bool, point_t>(false, 0);
      else
        return std::pair<bool, point_t>(true, idx->second);
    }
    return std::pair<bool, point_t>(false, 0);
  }

  std::pair<bool, point_t> getNextClosestPred(int currAzim, int currLaserId)
  {
    const auto& refPoints = refPointVals;
    auto quantizedPhi = computePhiQuantized(currAzim);
    if (refPoints[currLaserId].size()) {
      auto idx = refPoints[currLaserId].upper_bound(quantizedPhi);
      if (idx == refPoints[currLaserId].end())
        return std::pair<bool, point_t>(false, 0);
      else
        quantizedPhi = computePhiQuantized(idx->second[1]);
      idx = refPoints[currLaserId].upper_bound(quantizedPhi);
      if (idx == refPoints[currLaserId].end())
        return std::pair<bool, point_t>(false, 0);
      else
        return std::pair<bool, point_t>(true, idx->second);
    }
    return std::pair<bool, point_t>(false, 0);
  }

  int computePhiQuantized(const int val)
  {
    int offset = azimScaleLog2 ? (1 << (azimScaleLog2 - 1)) : 0;
    return val >= 0 ? (val + offset) >> azimScaleLog2
                    : -((-val + offset) >> azimScaleLog2);
  }

  bool isInterEnabled() { return interEnabled; }
  void setInterEnabled(bool enableFlag) { interEnabled = enableFlag; }
  void advanceFrame() { motionParams.AdvanceFrame(); }
  void parseMotionParams(std::string fileName)
  {
    motionParams.parseFile(fileName);
  }
  void updateFrame(const GeometryParameterSet& gps)
  {
    motionParams.AdvanceFrame();
    if (globalMotionEnabled) {
      CartesianToSphericalSimple cartToSpherical(gps);
      SphericalToCartesian sphericalToCart(gps);
      std::pair<int, int> currThresh;
      std::vector<int64_t> currMatrix;
      Vec3<int> currVector;

      for (auto& laserPts : refPointValsGlob)
        laserPts.clear();
      motionParams.getMotionParams(currThresh, currMatrix, currVector);
      for (auto laserId = 0; laserId < numLasers; laserId++) {
        for (auto ptIter : refPointValsCur[laserId]) {
          auto pt = sphericalToCart(ptIter.second);
          if (pt[2] > currThresh.first || pt[2] < currThresh.second) {
            auto p = pt;
            for (auto k = 0; k < 3; k++) {
              int64_t x =
                divExp2RoundHalfInfPositiveShift(
                  currMatrix[3 * k + 0] * p[0] + currMatrix[3 * k + 1] * p[1]
                    + currMatrix[3 * k + 2] * p[2],
                  16, 1 << 15)
                + currVector[k];
              pt[k] = x;
            }

            pt = cartToSpherical(pt);
          } else
            pt = ptIter.second;
          const auto phiQ = computePhiQuantized(pt[1]);
          if (!refPointValsGlob[pt[2]].count(phiQ))
            refPointValsGlob[pt[2]].insert({phiQ, pt});
          else if (refPointValsGlob[pt[2]].at(phiQ)[0] > pt[0])
            refPointValsGlob[pt[2]].at(phiQ) = pt;
        }
      }
      for (auto laserId = 0; laserId < numLasers; laserId++) {
        auto& ptsZero = refPointValsCur[laserId];
        auto& ptsGlob = refPointValsGlob[laserId];
        for (auto& ptIter : ptsZero) {
          pcc::point_t ptA = 0, ptB = 0;
          auto& pt = ptIter.second;
          auto phiQ = computePhiQuantized(pt[1]);
          if (ptsGlob.count(phiQ)) {
            const auto& colPt = ptsGlob[phiQ];
            ptA = colPt;
            if (colPt[1] < pt[1]) {
              auto idx = ptsGlob.upper_bound(phiQ);
              ptB = (idx == ptsGlob.end()) ? ptA : idx->second;
            } else if (colPt[1] > pt[1]) {
              auto idx = ptsGlob.lower_bound(phiQ);
              ptB = (idx == ptsGlob.begin()) ? ptA : std::prev(idx)->second;
            } else
              ptB = ptA;
          } else {
            auto idx = ptsGlob.upper_bound(phiQ);
            auto idx1 = idx;
            if (idx != ptsGlob.begin())
              idx1 = std::prev(idx);
            if (idx == ptsGlob.end())
              idx = idx1;
            ptA = idx->second;
            ptB = idx1->second;
          }
          int64_t delAzim = ptA[1] - ptB[1];
          int64_t delRad = (ptA[0] - ptB[0]);
          if (!delAzim || !delRad)
            pt[0] = ptA[0];
          else {
            const auto nr = delRad * (pt[1] - ptA[1]);
            const auto dr = delAzim;
            const bool sign =
              ((nr > 0 && dr > 0) || (nr < 0 && dr < 0)) ? 0 : 1;
            pt[0] =
              ptA[0] + (1 - 2 * sign) * (abs(nr) + (abs(dr) >> 1)) / abs(dr);
          }
        }
      }
      for (auto laserId = 0; laserId < numLasers; laserId++)
        refPointVals[laserId] = std::move(refPointValsCur[laserId]);
    } else {
      for (auto laserId = 0; laserId < numLasers; laserId++)
        refPointVals[laserId] = std::move(refPointValsCur[laserId]);
    }
  }
  template <typename T>
  void getMotionParams(std::pair<int, int> &th, std::vector<T> &mat, Vec3<int> &tr) const
  {
    motionParams.getMotionParams(th, mat, tr);
  }
  void setMotionParams(std::pair<int, int> &th, std::vector<int> &mat, Vec3<int> &tr) 
  {
    motionParams.setMotionParams(th, mat, tr);
  }


private:
  std::vector<std::map<int, pcc::point_t>> refPointVals;
  std::vector<std::map<int, pcc::point_t>> refPointValsGlob;
  std::vector<std::map<int, pcc::point_t>> refPointValsCur;
  MotionParameters motionParams;
  int numLasers;
  int azimScaleLog2;
  bool interEnabled;
  bool globalMotionEnabled;
};

//============================================================================

}  // namespace pcc
