/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2020, ISO/IEC
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

namespace pcc {

//=============================================================================
typedef std::vector<int> VecInt;
typedef std::vector<double> VecDouble;
typedef std::vector<std::vector<double>> VecVecDouble;

class MotionParameters {
private:
  int numFrames;
  std::vector<std::vector<int>> motionMatrix;
  VecVecDouble motionMatrix_double;
  std::vector<Vec3<int>> transVec;
  std::vector<std::pair<int, int>> threshVec;
  int frameCounter;
  int refFrameCounter;
  std::vector<bool> perFrameMovingStatus; 

public:
  MotionParameters() : numFrames(0), frameCounter(-1) {}
  bool EndOfFrames() const { return frameCounter >= numFrames; }
  void AdvanceFrame() { frameCounter++; }
  void setFrameCtr(const int frameCtr) { frameCounter = frameCtr; }
  void setRefFrameCtr(const int refFrameCtr) { refFrameCounter = refFrameCtr; }
  int getRefFrameCtr() { return refFrameCounter; }
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
  void parseFile(std::string fileName, double qs)
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
    motionMatrix_double.resize(numFrames);
    transVec.resize(numFrames);
    threshVec.resize(numFrames);
    const int scaleFactor = 65536;
    perFrameMovingStatus.resize(numFrames);

    auto it = out.begin();
    for (auto i = 0; i < numFrames; i++) {
      motionMatrix[i].resize(9);
      motionMatrix_double[i].resize(9);
      for (int j = 0; j < 9; j++) {
        motionMatrix_double[i][j] = *(it);
        if (j % 3 == j / 3)
          motionMatrix[i][j] = std::round((*(it++) - 1) * scaleFactor) + 65536;
        else
          motionMatrix[i][j] = std::round((*(it++)) * scaleFactor);
      }
      for (int j = 0; j < 3; j++)
        transVec[i][j] = std::round((*(it++)) * qs);
      threshVec[i].first =
        std::round((*(it++)) * qs);  // quantizeParameter(*(it++));
      threshVec[i].second =
        std::round((*(it++)) * qs);  // quantizeParameter(*(it++));
      
      perFrameMovingStatus[i] = checkMovingStatus(motionMatrix[i], transVec[i]);
    }
  }

  bool checkMovingStatus(std::vector<int> matrix, Vec3<int> vec){
    const int frameDistance = 1;
    const double scale = 65536.;
    const double thr1_pre = 0.1 / frameDistance;
    const double thr2_pre = 250.;
    const double thr1_tan_pre = tan(M_PI * thr1_pre / 180);
    const double thr1_sin_pre = sin(M_PI * thr1_pre / 180);
    const double Rx = std::abs(
      (matrix[5] / scale) / (1. + matrix[8] / scale));
    const double Ry = std::abs(matrix[2] / scale);
    const double Rz = std::abs(
      (matrix[1] / scale) / (1. + matrix[0] / scale));
    const double Sx = std::abs(vec[0]);
    const double Sy = std::abs(vec[1]);
    const double Sz = std::abs(vec[2]);
    return !(Rx < thr1_tan_pre && Ry < thr1_sin_pre && Rz < thr1_tan_pre
        && Sx < thr2_pre && Sy < thr2_pre && Sz < thr2_pre);    
  }

  void setMotionParams(
    std::pair<int, int> thresh, std::vector<int> matrix, Vec3<int> trans)
  {
    motionMatrix.push_back(matrix);
    threshVec.push_back(thresh);
    transVec.push_back(trans);
  }
  bool getMovingStatus() 
  {
    if (frameCounter > perFrameMovingStatus.size()) {
      std::cout << "Accessing unassigned values\n";
    }
    return perFrameMovingStatus[frameCounter];
  }
  void setMovingStatus(const bool movsts)  // decoder use
  {
    perFrameMovingStatus.push_back(movsts);
  }
  void updateNextMovingStatus(const bool movsts)
  {
    if (frameCounter + 1 > perFrameMovingStatus.size())
      return;
    perFrameMovingStatus[frameCounter + 1] = movsts;
  }
  void updateMovingStatus(const bool status, const int Ctr){
    if (Ctr >= perFrameMovingStatus.size())
      return;
    perFrameMovingStatus[Ctr] = status;
  }
  template<typename T>
  void getMotionParams(
    std::pair<int, int>& th, std::vector<T>& mat, Vec3<int>& tr) const
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
  void updateThresholds(
    const int frameCounter, const int leftThresh, const int rightTresh)
  {
    // todo: use variable instead of array for thresholds
    if (threshVec.size() > frameCounter)
      threshVec[frameCounter] = {rightTresh, leftThresh};
    else if (threshVec.size() == frameCounter) {
      /* Do nothing - last frame*/
    } else {
      throw std::runtime_error(
        "Check: frameCounter > threshVec.size()");
    }
  }
  void populateMatrix(
    VecVecDouble& mat, const VecDouble& rotMat, const Vec3<int> tran) const
  {
    for (size_t j = 0; j < 9; j++)
      mat[j / 3][j % 3] = rotMat[j];

    mat[3][0] = tran[0];
    mat[3][1] = tran[1];
    mat[3][2] = tran[2];
  }
  void getMotionParamsMultiple(
    VecVecDouble& mat, const int _curPicIndex, const int _prePicIndex) const
  {
    if (_curPicIndex > threshVec.size() || _prePicIndex > threshVec.size()) {
      std::cout << "Accessing unassigned values\n";
    }
    VecVecDouble mat_tmp = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
    VecVecDouble mat_accum = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

    populateMatrix(
      mat, motionMatrix_double[_prePicIndex], transVec[_prePicIndex]);

    for (size_t i = _prePicIndex + 1; i < _curPicIndex; i++) {
      populateMatrix(mat_tmp, motionMatrix_double[i], transVec[i]);

      for (int l = 0; l < 3; l++)
        for (int c = 0; c < 3; c++)
          mat_accum[l][c] = mat[l][0] * mat_tmp[0][c]
            + mat[l][1] * mat_tmp[1][c] + mat[l][2] * mat_tmp[2][c];
      for (int i = 0; i < 3; i++)
        mat_accum[3][i] = mat[3][0] * mat_tmp[0][i] + mat[3][1] * mat_tmp[1][i]
          + mat[3][2] * mat_tmp[2][i] + mat_tmp[3][i];
      mat = mat_accum;
    }
  }

  void ReverseGlobalMotion(const VecVecDouble& M1, VecVecDouble& re) const
  {
    for (int i = 0; i < 3; i++) {
      re[3][i] = -M1[3][i];
    }
    double rota_M1[3][3];
    for (int l = 0; l < 3; l++) {
      for (int c = 0; c < 3; c++) {
        rota_M1[l][c] = M1[l][c];
      }
    }

    double flag = GetA(rota_M1, 3);
    double t[3][3];
    double tmp[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 2; k++) {
          for (int t = 0; t < 2; t++) {
            tmp[k][t] = rota_M1[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
          }
        }

        t[j][i] = GetA(tmp, 2);
        if ((i + j) % 2 == 1) {
          t[j][i] = -t[j][i];
        }
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        re[i][j] = t[i][j] / flag;
      }
    }
  }

  double GetA(double arcs[3][3], int n) const
  {
    if (n == 1)
      return arcs[0][0];
    double ans = 0;
    double tmp[3][3] = {0.};
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n - 1; j++) {
        for (int k = 0; k < n - 1; k++) {
          tmp[j][k] = arcs[j + 1][(k >= i) ? k + 1 : k];
        }
      }
      double t = GetA(tmp, n - 1);

      ans += (i % 2 ? -1 : 1) * arcs[0][i] * t;
    }
    return ans;
  }

  template<typename T>
  void getMotionParamsMultiple2(
    std::pair<int, int>& th,
    std::vector<T>& mat,
    Vec3<int>& tr,
    const int _curPicIndex,
    const int _prePicIndex) const
  {
    if (_curPicIndex > threshVec.size() || _prePicIndex > threshVec.size()) {
      std::cout << "Accessing unassigned values\n";
    }

    if (motionMatrix_double.size()){
      th = threshVec
        [(_prePicIndex < threshVec.size()) ? _prePicIndex : (threshVec.size() - 1)];

      VecVecDouble mat_tmp = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

      if (_prePicIndex < _curPicIndex) {
        getMotionParamsMultiple(mat_tmp, _curPicIndex, _prePicIndex);
      } else {
        getMotionParamsMultiple(mat_tmp, _prePicIndex, _curPicIndex);
        VecVecDouble mat_tmp2 = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        ReverseGlobalMotion(mat_tmp, mat_tmp2);
        mat_tmp = std::move(mat_tmp2);
      }

      const float scaleFactor = 65536.;

      mat.resize(9);
      for (int i = 0; i < 9; i++) {
        mat[i] = (i / 3 == i % 3)
          ? (std::round((mat_tmp[i / 3][i % 3] - 1) * scaleFactor) + scaleFactor)
          : (std::round(mat_tmp[i / 3][i % 3] * scaleFactor));
      }

      for (int i = 0; i < 3; i++) {
        tr[i] = mat_tmp[3][i];
      }      
    }
    else{
      getMotionParams(th, mat, tr);
    }
  }
};
//=============================================================================

struct QtBtParameters {
  // maximum number of qtbt partitions before performing octree partitioning.
  int maxNumQtBtBeforeOt;

  // minimum size of qtbt partitions.
  int minQtbtSizeLog2;

  bool trisoupEnabled;

  bool angularTweakEnabled;
  int angularMaxNodeMinDimLog2ToSplitV;
  int angularMaxDiffToSplitZ;
};

//----------------------------------------------------------------------------

struct OctreeEncOpts {
  QtBtParameters qtbt;

  // Method used to derive in-tree quantisation parameters
  enum class QpMethod
  {
    kUniform = 0,
    kRandom = 1,
    kByDensity = 2,
  } qpMethod;

  // Tree depth at which to apply geometry quantisation
  int qpOffsetDepth;

  // Node size (rather than depth) at which to apply geometry quantisation
  int qpOffsetNodeSizeLog2;
};

//----------------------------------------------------------------------------

struct TrisoupEncOpts {
  bool improvedVertexDetermination;
  bool nodeUniqueDSE;
};

//----------------------------------------------------------------------------

struct PredGeomEncOpts {
  enum SortMode
  {
    kNoSort,
    kSortMorton,
    kSortAzimuth,
    kSortRadius,
    kSortLaserAngle
  } sortMode;

  // limit on number of points per tree
  int maxPtsPerTree;

  // Reciprocal bin width used in azimuthal sorting.
  //  0 => full precision
  float azimuthSortRecipBinWidth;

  int maxPredIdxTested;
  int radiusThresholdForNewPred;

  // Enable partition for Predgeom
  bool enablePartition;
  // height for partition
  int splitter;
};


//----------------------------------------------------------------------------

struct InterGeomEncOpts {
  enum MotionSource
  {
    kExternalGMSrc = 0,
    kInternalLMSGMSrc = 1,
    kInternalICPGMSrc = 2,
  } motionSrc;

  enum LPUType
  {
    kRoadObjClassfication = 0,
    kCuboidPartition = 1,
  } lpuType;

  MotionParameters motionParams;
	bool useCuboidalRegionsInGMEstimation;
  bool deriveGMThreshold;
  float gmThresholdHistScale;
  int gmThresholdMinZ;
  int gmThresholdMaxZ;
  float gmThresholdLeftScale;
  float gmThresholdRightScale;

  std::vector<int> motion_block_size;
  int motion_window_size;
  int th_dist;

};

//=============================================================================

}  // namespace pcc
