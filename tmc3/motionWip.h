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

#pragma once
#include "TMC3.h"
#include <vector>

#include "PCCPointSet.h"
#include "PCCTMC3Encoder.h"
#include "entropy.h"
#include "hls.h"

namespace pcc {

//============================================================================
static const unsigned int motionParamPrec = 16;
static const unsigned int motionParamScale = 1 << motionParamPrec;
static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);
int plus1log2shifted4(int x);

//============================================================================

void quantizeGlobalMotion(double Mat_GM[4][3], int32_t Mat_GM_Q[4][3]);

void applyGlobalMotion(
  std::vector<Vec3<int>>& listPoints, double Mat_GM[4][3]);

void applyGlobalMotion(
  PCCPointSet3& PC,
  const int32_t Mat_GM_Q[4][3],
  Vec3<double> vehicle_position,
  const int32_t global_thresholds[2]
);

double map_reference(
  std::vector<Vec3<int>>& pc_world_target,
  const std::vector<Vec3<int>>& pointPredictor_centered,
  std::vector<Vec3<int>>& pc_world_ref);

void LMS3D(
  std::vector<Vec3<int>>& P1,
  std::vector<Vec3<int>>& P2,
  std::vector<Vec3<int>>& pointPredictor_centered,
  uint32_t maxBB,
  double Mat_GM[4][3]);

void SearchGlobalMotion(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  double QS,
  int bsize,
	int th_dist,
  uint32_t maxBB,
  const bool useCuboidalRegionsInGMEstimation,
  std::vector<int>& gm_matrix,
  Vec3<int>& gm_trans,
  const std::pair<int, int> thresh);

void SearchGlobalMotionPerTile(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  double QS,
  GeometryBrickHeader& gbh,
	int th_dist,
  const bool useCuboidalRegionsInGMEstimation,
  const bool predDir);

void applyGlobalMotion_with_shift(
	PCCPointSet3& PC,
	const int32_t Mat_GM_Q[4][3],
	const Vec3<int> minimum_position
);

int roundIntegerHalfInf(const double x);

void
compensateWithRoadObjClassfication(
	PCCPointSet3& pointPredictorWorld,
  const std::vector<int> gm_matrix,
  const Vec3<int> gm_trans,
  const std::pair<int, int> thresh,
  const Vec3<int> minimum_position
);

void
compensateWithCuboidPartition(
	PCCPointSet3& pointCloud,
	PCCPointSet3& predPointCloud,
	PCCPointSet3& pointPredictorWorld,
  const GeometryBrickHeader& gbh,
  int motion_window_size,
  const Vec3<int> minimum_position,
  EntropyEncoder* arithmeticEncoder,
  const bool predDir);

void
decodeCompensateWithCuboidPartition(
	PCCPointSet3& predPointCloud,
	PCCPointSet3& pointPredictorWorld,
  const GeometryBrickHeader& gbh,
	const Vec3<int> minimum_position,
  EntropyDecoder* arithmeticDecoder,
  const bool predDir);

  //============================================================================

}  // namespace pcc
