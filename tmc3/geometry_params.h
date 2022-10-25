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
        transVec[i][j] = std::round((*(it++)) * qs);
      threshVec[i].first = std::round((*(it++)) * qs);  // quantizeParameter(*(it++));
      threshVec[i].second = std::round((*(it++)) * qs); // quantizeParameter(*(it++));
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

  std::vector<int> motion_block_size;
  int motion_window_size;
  int th_dist;

};

//=============================================================================

}  // namespace pcc
