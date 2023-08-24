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
#include "TMC3.h"

#include "motionWip.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <set>
#include <vector>
#include <unordered_map>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================
struct MotionEntropy {
  AdaptiveBitModelFast _ctxIsWorld;
  MotionEntropy();
};

//----------------------------------------------------------------------------
class MotionEntropyEncoder : public MotionEntropy {
public:
  MotionEntropyEncoder(EntropyEncoder* arithmeticEncoder)
    : _arithmeticEncoder(arithmeticEncoder)
  {}

  void encodeIsWorld(bool hasMotion);

private:
  EntropyEncoder* _arithmeticEncoder;
};

//----------------------------------------------------------------------------
class MotionEntropyDecoder : public MotionEntropy {
public:
  MotionEntropyDecoder(EntropyDecoder* arithmeticDecoder)
    : _arithmeticDecoder(arithmeticDecoder)
  {}

  bool decodeIsWorld();

private:
  EntropyDecoder* _arithmeticDecoder;
};

//============================================================================
MotionEntropy::MotionEntropy()
{}

//----------------------------------------------------------------------------
// Encode the presence of an horizontal segment
inline void
MotionEntropyEncoder::encodeIsWorld(bool present)
{
  _arithmeticEncoder->encode(present, _ctxIsWorld);
}

//----------------------------------------------------------------------------
// Decode horizontal segment presence flag
inline bool
MotionEntropyDecoder::decodeIsWorld()
{
  return _arithmeticDecoder->decode(_ctxIsWorld);
}

//============================================================================
static const int LUT_LOG2[64]{
  INT_MIN, 0,  16, 25, 32, 37, 41, 45, 48, 51, 53, 55, 57, 59, 61, 63,
  64,      65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 79,
  80,      81, 81, 82, 83, 83, 84, 85, 85, 86, 86, 87, 87, 88, 88, 89,
  89,      90, 90, 91, 91, 92, 92, 93, 93, 93, 94, 94, 95, 95, 95, 96};

//----------------------------------------------------------------------------
int
plus1log2shifted4(int x)
{
  x++;
  int result = 0;
  while (x >= 64) {
    x >>= 1;
    result += 16;
  }

  return result + LUT_LOG2[x];
}

//----------------------------------------------------------------------------
static double
calcCostOfGlobalMotion(
  std::vector<Vec3<int>>& Window,
  std::vector<Vec3<int>>& Block0,
  int* bufferPoints,
  int wSize)
{
  if (!Window.size())
    return DBL_MAX;

  double cost = DBL_MAX;
  const int samples = 4;
  const int decimate = 6;

  if (Window.size() > samples * std::max(int(Block0.size()), 16))
    wSize >>= 1;
  int maxDistance = wSize << 1;

  int dist = 0;
  int* pBuffer = bufferPoints;

  auto* itB = Block0.data();
  int jumpBlock = 1
    + (Block0.size()
       >> decimate);  // (kind of) random sampling of the original block to code
  int nSamples = 0;
  for (int Nb = 0; Nb < Block0.size();
       Nb += jumpBlock, itB += jumpBlock, nSamples++) {
    Vec3<int> b = *itB;
    int min_d = maxDistance;
    for (const auto w : Window) {
      int a[3];
      for (auto i = 0; i < 3; i++) {
        pBuffer[i] = b[i] - w[i];
        a[i] = std::abs(pBuffer[i]);
      }
      if (a[0] <= wSize && a[1] <= wSize && a[2] <= wSize)
        pBuffer += 3;

      a[0] += a[1] + a[2];  // replaces int local_d
      if (a[0] < min_d)
        min_d = a[0];
    }
    dist += plus1log2shifted4(min_d);  // 1/0.0625 = 16 times log
  }
  cost = jumpBlock * dist;
  return cost;
}

//============================================================================

void
populateWindowList(
  VecVecVec3& windowList,
  PCCPointSet3& pointPred,
  std::vector<int> motion_block_size,
  int LPUnumInAxis[3],
  const Vec3<int> min)
{
  for (size_t i = 0; i < pointPred.getPointCount(); ++i) {
    const Vec3<int> point = pointPred[i];
    const int xidx =
      motion_block_size[0] ? (point[0] - min.x()) / motion_block_size[0] : 0;
    if (xidx < 0 || xidx >= LPUnumInAxis[0])
      continue;
    const int yidx =
      motion_block_size[1] ? (point[1] - min.y()) / motion_block_size[1] : 0;
    if (yidx < 0 || yidx >= LPUnumInAxis[1])
      continue;
    const int zidx =
      motion_block_size[2] ? (point[2] - min.z()) / motion_block_size[2] : 0;
    if (zidx < 0 || zidx >= LPUnumInAxis[2])
      continue;
    const int idx = (xidx * LPUnumInAxis[1] + yidx) * LPUnumInAxis[2] + zidx;
    windowList[idx].push_back(point);
  }
}

void
compensateCuboidGlobalMotion(
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictor,
  PCCPointSet3& pointPredictorWorld,
  std::vector<bool>& isWorldList,
  std::vector<int> motion_block_size,
  int LPUnumInAxis[3],
  const Vec3<int> min)
{
  VecVecVec3 Window_W_List, Window_V_List;
  const int blocksize = LPUnumInAxis[0] * LPUnumInAxis[1] * LPUnumInAxis[2];
  Window_W_List.resize(blocksize);
  Window_V_List.resize(blocksize);
  populateWindowList(
    Window_W_List, pointPredictorWorld, motion_block_size, LPUnumInAxis, min);
  populateWindowList(
    Window_V_List, pointPredictor, motion_block_size, LPUnumInAxis, min);

  std::vector<Vec3<double>> pointPredictorMC;
  for (int idx = 0; idx < blocksize; idx++) {
    const auto& windowList =
      isWorldList[idx] ? Window_W_List[idx] : Window_V_List[idx];
    for (const auto& point : windowList)
      pointPredictorMC.push_back(point);
  }

  compensatedPointCloud->resize(pointPredictorMC.size());
  int counter = 0;
  for (const auto& p : pointPredictorMC) {
    auto& predPoint = (*compensatedPointCloud)[counter++];
    predPoint[0] = p[0];
    predPoint[1] = p[1];
    predPoint[2] = p[2];
  }
}
void
populateCuboidBlocks(
  VecVecVec3& windowList,
  const PCCPointSet3& pointCloud,
  std::vector<int> motion_block_size,
  std::vector<int> th_dists,
  Box3<int32_t> bbox,
  int LPUnumInAxis[3])
{
  const int samples = 4;
  for (size_t i = 0; i < pointCloud.getPointCount(); i += samples) {
    std::unordered_map<int, bool> lpuToAdd;
    const Vec3<int> point = pointCloud[i];
    for (size_t m = 0; m < th_dists.size(); m++) {
      const int xidx = motion_block_size[0]
        ? (point[0] + th_dists[m] - bbox.min.x()) / motion_block_size[0]
        : 0;
      if (xidx < 0 || xidx >= LPUnumInAxis[0])
        continue;
      for (size_t n = 0; n < th_dists.size(); n++) {
        const int yidx = motion_block_size[1]
          ? (point[1] + th_dists[n] - bbox.min.y()) / motion_block_size[1]
          : 0;
        if (yidx < 0 || yidx >= LPUnumInAxis[1])
          continue;
        for (size_t k = 0; k < th_dists.size(); k++) {
          const int zidx = motion_block_size[2]
            ? (point[2] + th_dists[k] - bbox.min.z()) / motion_block_size[2]
            : 0;
          if (zidx < 0 || zidx >= LPUnumInAxis[2])
            continue;
          const int idx =
            (xidx * LPUnumInAxis[1] + yidx) * LPUnumInAxis[2] + zidx;
          lpuToAdd[idx] = true;
        }
      }
    }
    for (auto idx : lpuToAdd)
      windowList[idx.first].push_back(point);
  }
}
void
encodeCuboidGlobalMotion(
  const PCCPointSet3& pointCloud,
  std::vector<int> motion_block_size,
  int motion_window_size,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3& pointPredictor,
  int* bufferPoints,
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorWorld)
{
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  const auto bbox = pointPredictor.computeBoundingBox();

  int LPUnumInAxis[3] = {1, 1, 1};
  for (int i = 0; i < 3; i++) {
    if (motion_block_size[i])
      LPUnumInAxis[i] = (bbox.max[i] - bbox.min[i] + motion_block_size[i] - 1)
        / motion_block_size[i];
  }
  const int extendedWindow = motion_window_size;
  std::vector<int> th_dists;
  th_dists.push_back(extendedWindow);
  if (motion_window_size)
    th_dists.push_back(-extendedWindow);

  const int blockSize = LPUnumInAxis[0] * LPUnumInAxis[1] * LPUnumInAxis[2];

  VecVecVec3 Block0_List, Window_W_List, Window_V_List;
  Block0_List.resize(blockSize);
  Window_W_List.resize(blockSize);
  Window_V_List.resize(blockSize);
  populateCuboidBlocks(
    Block0_List, pointCloud, motion_block_size, th_dists, bbox, LPUnumInAxis);
  populateCuboidBlocks(
    Window_W_List, pointPredictorWorld, motion_block_size, th_dists, bbox,
    LPUnumInAxis);
  populateCuboidBlocks(
    Window_V_List, pointPredictor, motion_block_size, th_dists, bbox,
    LPUnumInAxis);

  std::vector<bool> isWorldList(blockSize, true);
  for (int ith = 0; ith < blockSize; ith++) {
    if (
      !Block0_List[ith].size()
      || (!Window_W_List[ith].size() && !Window_V_List[ith].size())) {
      continue;
    }

    double costWorld = calcCostOfGlobalMotion(
      Window_W_List[ith], Block0_List[ith], bufferPoints, motion_window_size);

    double costVehicle = calcCostOfGlobalMotion(
      Window_V_List[ith], Block0_List[ith], bufferPoints, motion_window_size);

    if (!Window_W_List[ith].size() || costWorld >= costVehicle)
      isWorldList[ith] = false;
  }

  Window_V_List.clear();
  Window_W_List.clear();
  Block0_List.clear();

  // encoding
  for (int idx = 0; idx < blockSize; idx++)
    motionEncoder.encodeIsWorld(isWorldList[idx]);

  // compensation
  compensateCuboidGlobalMotion(
    compensatedPointCloud, pointPredictor, pointPredictorWorld, isWorldList,
    motion_block_size, LPUnumInAxis, bbox.min);
}

//----------------------------------------------------------------------------
void
decodeCuboidGlobalMotion(
  std::vector<int> motion_block_size,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3& pointPredictor,
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorWorld)
{
  MotionEntropyDecoder motionDecoder(arithmeticDecoder);

  const auto bbox = pointPredictor.computeBoundingBox();
  int LPUnumInAxis[3] = {1, 1, 1};
  for (int i = 0; i < 3; i++) {
    if (motion_block_size[i])
      LPUnumInAxis[i] = (bbox.max[i] - bbox.min[i] + motion_block_size[i] - 1)
        / motion_block_size[i];
  }

  const int blockSize = LPUnumInAxis[0] * LPUnumInAxis[1] * LPUnumInAxis[2];
  std::vector<bool> isWorldList(blockSize, false);

  // decoding
  for (int idx = 0; idx < blockSize; idx++)
    isWorldList[idx] = motionDecoder.decodeIsWorld();

  // compensation
  compensateCuboidGlobalMotion(
    compensatedPointCloud, pointPredictor, pointPredictorWorld, isWorldList,
    motion_block_size, LPUnumInAxis, bbox.min);
}

//============================================================================
void
quantizeGlobalMotion(double Mat_GM[4][3], int32_t Mat_GM_Q[4][3])
{
  double scale = motionParamScale;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++) {
      if (l == c)
        Mat_GM_Q[l][c] =
          roundIntegerHalfInf((Mat_GM[l][c] - 1.) * scale) + scale;
      else if (l < 3)
        Mat_GM_Q[l][c] = roundIntegerHalfInf(Mat_GM[l][c] * scale);
      else
        Mat_GM_Q[l][c] = roundIntegerHalfInf(Mat_GM[l][c]);
    }
}

//----------------------------------------------------------------------------
void
applyGlobalMotion(std::vector<Vec3<int>>& listPoints, double Mat_GM[4][3])
{
  for (auto& b : listPoints) {
    Vec3<int> point;
    for (auto i = 0; i < 3; i++)
      point[i] = roundIntegerHalfInf(
        Mat_GM[0][i] * b[0] + Mat_GM[1][i] * b[1] + Mat_GM[2][i] * b[2]
        + Mat_GM[3][i]);
    b = point;
  }
}

//----------------------------------------------------------------------------
void
applyGlobalMotion(
  PCCPointSet3& PC,
  const int32_t Mat_GM_Q[4][3],
  Vec3<double> vehicle_position,
  const int32_t global_thresholds[2])
{
  // unquantize
  double Mat_GM[4][3];
  double scale = motionParamScale;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++) {
      if (l == c)
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]) / scale + 1.;
      else if (l < 3)
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]) / scale;
      else
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]);
    }

  int topZ = global_thresholds[0];
  int bottomZ = global_thresholds[1];

  // apply motion
  for (int n = 0; n < PC.getPointCount(); n++) {
    Vec3<double> b = PC[n] + vehicle_position;

    if (b[2] < bottomZ || b[2] > topZ) {
      for (auto i = 0; i < 3; i++)
        PC[n][i] = roundIntegerHalfInf(
          Mat_GM[0][i] * b[0] + Mat_GM[1][i] * b[1] + Mat_GM[2][i] * b[2]
          + Mat_GM[3][i] - vehicle_position[i]);
    }
  }
}

//----------------------------------------------------------------------------
int
roundIntegerHalfInf(const double x)
{
  return (x >= 0) ? int(x + 0.5) : -int(-x + 0.5);
}

//----------------------------------------------------------------------------
double
map_reference(
  VecVec3& pcWorldTarget,
  const VecVec3& pointPredictorCentered,
  VecVec3& pcWorldRef)
{
  std::vector<int> accu_m;
  int meanM = 0;
  for (const auto& b : pcWorldTarget) {
    int dmin = 1 << 30;
    Vec3<int> closest;
    for (const auto& w : pointPredictorCentered) {
      const int L =
        std::abs(w[0] - b[0]) + std::abs(w[1] - b[1]) + std::abs(w[2] - b[2]);
      if (L < dmin) {
        dmin = L;
        closest = w;
      }
    }
    pcWorldRef.push_back(closest);

    accu_m.push_back(dmin);
    meanM += dmin;
  }

  double err = double(meanM) / pcWorldTarget.size() / 3.;

  // eliminate outliers
  int count = 0;
  auto it_accu = accu_m.begin();
  auto it_target = pcWorldTarget.begin();
  auto it_ref = pcWorldRef.begin();
  auto it_target_fill = pcWorldTarget.begin();
  auto it_ref_fill = pcWorldRef.begin();
  for (; it_accu != accu_m.end(); it_accu++, it_target++, it_ref++) {
    if (*it_accu * accu_m.size() <= 2 * meanM) {
      *it_target_fill++ = *it_target;
      *it_ref_fill++ = *it_ref;
      count++;
    }
  }

  pcWorldTarget.resize(count);
  pcWorldRef.resize(count);

  return err;
}

//----------------------------------------------------------------------------
void
LMS3D(
  std::vector<Vec3<int>>& P1,
  std::vector<Vec3<int>>& P2,
  std::vector<Vec3<int>>& pointPredictor_centered,
  uint32_t maxBB,
  double Mat_GM[4][3])
{
  // determine correlation matrix M in (X,Y,Z,MV_unity)
  const int MV_unity = maxBB >> 4;  //  // for better matrix conditioning
  double M[4][4] = {
    {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}};

  for (auto it1 = P1.begin(); it1 != P1.end(); it1++) {
    Vec3<double> Pref = *it1;
    //{(*it1)[0], (*it1)[1], (*it1)[2]};

    // X
    M[0][0] += Pref[0] * Pref[0];
    M[0][1] += Pref[0] * Pref[1];
    M[0][2] += Pref[0] * Pref[2];
    M[0][3] += Pref[0] * MV_unity;
    // Y
    M[1][1] += Pref[1] * Pref[1];
    M[1][2] += Pref[1] * Pref[2];
    M[1][3] += Pref[1] * MV_unity;
    // Z
    M[2][2] += Pref[2] * Pref[2];
    M[2][3] += Pref[2] * MV_unity;
    // 1
    M[3][3] += MV_unity * MV_unity;
  }
  M[1][0] = M[0][1];
  M[2][0] = M[0][2];
  M[2][1] = M[1][2];
  M[3][0] = M[0][3];
  M[3][1] = M[1][3];
  M[3][2] = M[2][3];

  // inverse M by Gauss pivoting
  double invM[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  for (int pivot = 0; pivot < 3; pivot++)  // descente
  {
    double value_pivot = M[pivot][pivot];
    for (int l = pivot + 1; l < 4; l++) {
      double factor = -M[l][pivot] / value_pivot;
      for (int c = 0; c < 4; c++) {
        M[l][c] += M[pivot][c] * factor;
        invM[l][c] += invM[pivot][c] * factor;
      }
    }
  }

  for (int pivot = 3; pivot > 0; pivot--)  // mont�e
  {
    double value_pivot = M[pivot][pivot];
    for (int l = pivot - 1; l >= 0; l--) {
      double factor = -M[l][pivot] / value_pivot;
      for (int c = 0; c < 4; c++) {
        M[l][c] += M[pivot][c] * factor;
        invM[l][c] += invM[pivot][c] * factor;
      }
    }
  }

  for (int pivot = 0; pivot < 4; pivot++)  // normalisation
  {
    double factor = 1 / M[pivot][pivot];
    for (int c = 0; c < 4; c++)
      invM[pivot][c] *= factor;
  }

  // determine rhs matrix R
  double R[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  for (auto it1 = P1.begin(), it2 = P2.begin(); it1 != P1.end();
       it1++, it2++) {
    Vec3<double> Pref = *it1;
    Vec3<double> Ptarget = *it2;

    // X
    R[0][0] += Ptarget[0] * Pref[0];
    R[1][0] += Ptarget[0] * Pref[1];
    R[2][0] += Ptarget[0] * Pref[2];
    R[3][0] += Ptarget[0] * MV_unity;
    // Y
    R[0][1] += Ptarget[1] * Pref[0];
    R[1][1] += Ptarget[1] * Pref[1];
    R[2][1] += Ptarget[1] * Pref[2];
    R[3][1] += Ptarget[1] * MV_unity;
    // Z
    R[0][2] += Ptarget[2] * Pref[0];
    R[1][2] += Ptarget[2] * Pref[1];
    R[2][2] += Ptarget[2] * Pref[2];
    R[3][2] += Ptarget[2] * MV_unity;
  }

  // apply inv M to R to get the transformation matrix T
  double T[4][3];
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      T[l][c] = invM[l][0] * R[0][c] + invM[l][1] * R[1][c]
        + invM[l][2] * R[2][c] + invM[l][3] * R[3][c];

  // deconditioning of 1 <-> MV_unity
  for (int c = 0; c < 3; c++)
    T[3][c] *= double(MV_unity);

  // penalization
  double lambda = 1.0;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      T[l][c] *= lambda;
  T[0][0] += 1 - lambda;
  T[1][1] += 1 - lambda;
  T[2][2] += 1 - lambda;

  // apply T to global motion matrix  Mat_GM
  double Mat_GM1[4][3];  //copy old GM matrix
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      Mat_GM1[l][c] = Mat_GM[l][c];

  for (int l = 0; l < 3; l++)  // deformation part
    for (int c = 0; c < 3; c++)
      Mat_GM[l][c] = Mat_GM1[l][0] * T[0][c] + Mat_GM1[l][1] * T[1][c]
        + Mat_GM1[l][2] * T[2][c];

  for (int c = 0; c < 3; c++)  // translation part
    Mat_GM[3][c] = Mat_GM1[3][0] * T[0][c] + Mat_GM1[3][1] * T[1][c]
      + Mat_GM1[3][2] * T[2][c] + T[3][c];
}

//----------------------------------------------------------------------------
void
PopulatePCLikelyWorld(
  const int blocknumInAxis,
  const int th_dist,
  const int bsize,
  const int top_z_boundary,
  const int bottom_z_boundary,
  const bool useCuboidalRegionsInGMEstimation,
  const PCCPointSet3& pointCloud,
  const PCCPointSet3& pointPredictor,
  std::vector<Vec3<int>>& pc_likely_world)
{
  int th_dists[2] = {th_dist, -th_dist};
  std::vector<bool> ref_planar_region;

  if (!useCuboidalRegionsInGMEstimation) {
    ref_planar_region.resize(
      blocknumInAxis * blocknumInAxis * blocknumInAxis, false);

    // scan each point in reference frame
    // use cubic block as the process unit
    for (size_t i = 0; i < pointPredictor.getPointCount(); i++) {
      const Vec3<double> point = pointPredictor[i];
      double x = point[0], y = point[1], z = point[2];
      for (size_t m = 0; m < 2; m++) {
        int xidx = (x + th_dists[m]) / bsize;
        if (xidx < 0 || xidx >= blocknumInAxis)
          continue;
        for (size_t n = 0; n < 2; n++) {
          int yidx = (y + th_dists[n]) / bsize;
          if (yidx < 0 || yidx >= blocknumInAxis)
            continue;
          for (size_t k = 0; k < 2; k++) {
            int zidx = (z + th_dists[k]) / bsize;
            if (zidx < 0 || zidx >= blocknumInAxis)
              continue;
            int idx = (xidx * blocknumInAxis + yidx) * blocknumInAxis + zidx;
            ref_planar_region[idx] = true;
          }
        }
      }
    }
    //scan each point in current frame
    for (size_t i = 0; i < pointCloud.getPointCount(); i++) {
      const Vec3<double> point = pointCloud[i];
      int xidx = point[0] / bsize;
      int yidx = point[1] / bsize;
      int zidx = point[2] / bsize;
      int idx = (xidx * blocknumInAxis + yidx) * blocknumInAxis + zidx;
      if (idx >= ref_planar_region.size() || !ref_planar_region[idx])
        continue;
      const Vec3<int> b = {int(point[0]), int(point[1]), int(point[2])};
      if ((b[2] < bottom_z_boundary || b[2] > top_z_boundary)) {
        pc_likely_world.push_back(b);
      }
    }
  } else {  // Use cuboidal regions
    ref_planar_region.resize(blocknumInAxis * blocknumInAxis, false);

    // scan each point in reference frame
    // use cubic block as the process unit
    for (size_t i = 0; i < pointPredictor.getPointCount(); i++) {
      const Vec3<double> point = pointPredictor[i];
      double x = point[0], y = point[1], z = point[2];
      for (size_t m = 0; m < 2; m++) {
        int xidx = (x + th_dists[m]) / bsize;
        if (xidx < 0 || xidx >= blocknumInAxis)
          continue;
        for (size_t n = 0; n < 2; n++) {
          int yidx = (y + th_dists[n]) / bsize;
          if (yidx < 0 || yidx >= blocknumInAxis)
            continue;
          int idx = xidx * blocknumInAxis + yidx;
          ref_planar_region[idx] = true;
        }
      }
    }
    //scan each point in current frame
    for (size_t i = 0; i < pointCloud.getPointCount(); i++) {
      const Vec3<double> point = pointCloud[i];
      int xidx = point[0] / bsize;
      int yidx = point[1] / bsize;
      int idx = xidx * blocknumInAxis + yidx;
      if (idx >= ref_planar_region.size() || !ref_planar_region[idx])
        continue;
      const Vec3<int> b = {int(point[0]), int(point[1]), int(point[2])};
      if ((b[2] < bottom_z_boundary || b[2] > top_z_boundary)) {
        pc_likely_world.push_back(b);
      }
    }
  }
}

//----------------------------------------------------------------------------
void
SearchGlobalMotion(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  double QS,
  int bsize,
  int th_dist,
  uint32_t maxBB,
  const bool useCuboidalRegionsInGMEstimation,
  std::vector<int>& gm_matrix,
  Vec3<int>& gm_trans,
  const std::pair<int, int> thresh)
{
  // ------------- first pass: find world-referential-likely LPU ----
  // number of LCU
  uint32_t maxBB_Scalled = (maxBB);

  // loop on LCU
  std::vector<Vec3<int>> pcLikelyWorld;

  int topZ = thresh.first;
  int bottomZ = thresh.second;
  int blocknumInAxis = (maxBB_Scalled % bsize) ? (maxBB_Scalled / bsize + 1)
                                               : (maxBB_Scalled / bsize);

  PopulatePCLikelyWorld(
    blocknumInAxis, th_dist, bsize, topZ, bottomZ,
    useCuboidalRegionsInGMEstimation, pointCloud, pointPredictor,
    pcLikelyWorld);

  Vec3<double> minPositions = {0., 0., 0.};
  Vec3<double> minPositions_local = minPositions * QS;
  Vec3<int> minPositions_int = Vec3<int>(
    int(minPositions_local[0]), int(minPositions_local[1]),
    int(minPositions_local[2]));
  for (auto& b : pcLikelyWorld)
    b += minPositions_int;

  std::vector<Vec3<int>> pointPredictorCentered;
  for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
    const Vec3<double> point = pointPredictor[i] + minPositions_local;
    ;
    pointPredictorCentered.push_back(
      Vec3<int>(int(point[0]), int(point[1]), int(point[2])));
  }
  std::vector<Vec3<int>> pointPredictorCentered0 = pointPredictorCentered;

  // global motion found iteratively using LMS

  int NLMS = 1;
  int nb_points = 100;

  int jump = 1 + (pcLikelyWorld.size() / nb_points);

  double Mat_GM[4][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {0, 0, 0}};  // global motion is identity as a start

  //std::cout << "(points/err) = ";

  for (int i = 0; i < NLMS; i++) {
    // sample pc_likely_world
    std::vector<Vec3<int>> pcWorldTarget;
    for (int N = i; N < pcLikelyWorld.size(); N += jump) {
      auto it = pcLikelyWorld.begin() + N;
      pcWorldTarget.push_back(*it);
    }

    // map reference to pc_world
    std::vector<Vec3<int>> pcWorldRef;
    map_reference(pcWorldTarget, pointPredictorCentered, pcWorldRef);
    // Least Mean Square 3D
    LMS3D(
      pcWorldRef, pcWorldTarget, pointPredictorCentered, maxBB_Scalled,
      Mat_GM);

    // apply global motion
    if (NLMS > 1) {  // Unnecessary when NLMS = 1
      pointPredictorCentered = pointPredictorCentered0;
      applyGlobalMotion(pointPredictorCentered, Mat_GM);
    }
  }
  //std::cout << std::endl;
  int Mat_GM_Q[4][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {0, 0, 0}};  // global motion is identity as a start
  quantizeGlobalMotion(Mat_GM, Mat_GM_Q);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gm_matrix[3 * i + j] = Mat_GM_Q[j][i];
    }
    gm_trans[i] = Mat_GM_Q[3][i];
  }
}

//----------------------------------------------------------------------------
void
SearchGlobalMotionPerTile(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  double QS,
  GeometryBrickHeader& gbh,
  int th_dist,
  const bool useCuboidalRegionsInGMEstimation,
  const bool predDir)
{
  int maxBB = (1 << gbh.maxRootNodeDimLog2) - 1;

  if (predDir)
    SearchGlobalMotion(
      pointCloud, pointPredictor, QS, gbh.motion_block_size[2], th_dist, maxBB,
      useCuboidalRegionsInGMEstimation, gbh.gm_matrix2, gbh.gm_trans2,
      gbh.gm_thresh2);
  else
    SearchGlobalMotion(
      pointCloud, pointPredictor, QS, gbh.motion_block_size[2], th_dist, maxBB,
      useCuboidalRegionsInGMEstimation, gbh.gm_matrix, gbh.gm_trans,
      gbh.gm_thresh);
}

//----------------------------------------------------------------------------
void
applyGlobalMotion_with_shift(
  PCCPointSet3& pointPredictorWorld,
  const std::vector<int> gm_matrix,
  const Vec3<int> gm_trans,
  const Vec3<int> minimum_position)
{
  static const unsigned int motionParamPrec = 16;
  static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);

  // apply
  Vec3<int> b = 0;
  for (int n = 0; n < pointPredictorWorld.getPointCount(); n++) {
    b = pointPredictorWorld[n];

    b[0] = b[0] + minimum_position[0];
    b[1] = b[1] + minimum_position[1];
    b[2] = b[2] + minimum_position[2];

    for (int i = 0; i < 3; i++) {
      pointPredictorWorld[n][i] =
        divExp2RoundHalfInfPositiveShift(
          gm_matrix[3 * i] * b[0] + gm_matrix[3 * i + 1] * b[1]
            + gm_matrix[3 * i + 2] * b[2],
          motionParamPrec, motionParamOffset)
        + gm_trans[i] - minimum_position[i];  // translation and offset
    }
  }
}

//----------------------------------------------------------------------------
void
compensateWithRoadObjClassfication(
  PCCPointSet3& pointPredictorWorld,
  const std::vector<int> gm_matrix,
  const Vec3<int> gm_trans,
  const std::pair<int, int> thresh,
  const Vec3<int> minimum_position)
{
  static const unsigned int motionParamPrec = 16;
  static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);

  // apply
  Vec3<int> b = 0;
  for (int n = 0; n < pointPredictorWorld.getPointCount(); n++) {
    b = pointPredictorWorld[n];

    b[0] = b[0] + minimum_position[0];
    b[1] = b[1] + minimum_position[1];
    b[2] = b[2] + minimum_position[2];

    if ((b[2] < thresh.second) || (b[2] > thresh.first)) {
      for (int i = 0; i < 3; i++) {
        pointPredictorWorld[n][i] =
          divExp2RoundHalfInfPositiveShift(
            gm_matrix[3 * i] * b[0] + gm_matrix[3 * i + 1] * b[1]
              + gm_matrix[3 * i + 2] * b[2],
            motionParamPrec, motionParamOffset)
          + gm_trans[i] - minimum_position[i];  // translation and offset
      }
    }
  }
}

//----------------------------------------------------------------------------
// output: pointPredictorWorld : predPointCloud + applied global motion
//         compensatedPointCloud : predPointCloud + Octree-based partition
void
compensateWithCuboidPartition(
  PCCPointSet3& pointCloud,
  PCCPointSet3& predPointCloud,
  PCCPointSet3& pointPredictorWorld,
  const GeometryBrickHeader& gbh,
  int motion_window_size,
  const Vec3<int> minimum_position,
  EntropyEncoder* arithmeticEncoder,
  const bool predDir)
{
  if (predDir)
    applyGlobalMotion_with_shift(
      pointPredictorWorld, gbh.gm_matrix2, gbh.gm_trans2, minimum_position);
  else
    applyGlobalMotion_with_shift(
      pointPredictorWorld, gbh.gm_matrix, gbh.gm_trans, minimum_position);

  std::unique_ptr<int[]> bufferPoints;
  bufferPoints.reset(new int[3 * 3000 * 10000]);
  PCCPointSet3 compensatedPointCloud;

  encodeCuboidGlobalMotion(
    pointCloud, gbh.motion_block_size, motion_window_size, arithmeticEncoder,
    predPointCloud, bufferPoints.get(), &compensatedPointCloud,
    pointPredictorWorld);

  bufferPoints.reset();
  pointPredictorWorld.clear();
  pointPredictorWorld = compensatedPointCloud;
}

//----------------------------------------------------------------------------
void
decodeCompensateWithCuboidPartition(
  PCCPointSet3& predPointCloud,
  PCCPointSet3& pointPredictorWorld,
  const GeometryBrickHeader& gbh,
  const Vec3<int> minimum_position,
  EntropyDecoder* arithmeticDecoder,
  const bool predDir)
{
  if (predDir)
    applyGlobalMotion_with_shift(
      pointPredictorWorld, gbh.gm_matrix2, gbh.gm_trans2, minimum_position);
  else
    applyGlobalMotion_with_shift(
      pointPredictorWorld, gbh.gm_matrix, gbh.gm_trans, minimum_position);

  PCCPointSet3 compensatedPointCloud;

  decodeCuboidGlobalMotion(
    gbh.motion_block_size, arithmeticDecoder, predPointCloud,
    &compensatedPointCloud, pointPredictorWorld);

  pointPredictorWorld = compensatedPointCloud;
}

//----------------------------------------------------------------------------
}  // namespace pcc
