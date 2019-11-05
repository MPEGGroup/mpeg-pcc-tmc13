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

#include <map>

#include "PCCPointSet.h"

namespace pcc {

//============================================================================

struct RecolourParams {
  double distOffsetFwd;
  double distOffsetBwd;
  double maxGeometryDist2Fwd;
  double maxGeometryDist2Bwd;
  double maxAttributeDist2Fwd;
  double maxAttributeDist2Bwd;

  int searchRange;
  int numNeighboursFwd;
  int numNeighboursBwd;

  bool useDistWeightedAvgFwd;
  bool useDistWeightedAvgBwd;
  bool skipAvgIfIdenticalSourcePointPresentFwd;
  bool skipAvgIfIdenticalSourcePointPresentBwd;
};

//============================================================================
// Quantise the geometry of a point cloud, retaining unique points only.
// Points in the @src point cloud are translated by -@offset, quantised by a
// multiplicitive @scaleFactor with rounding, then clamped to @clamp.
//
// The destination and source point clouds may be the same object.
//
// NB: attributes are not processed.

void quantizePositionsUniq(
  const float scaleFactor,
  const Vec3<int> offset,
  const Box3<int> clamp,
  const PCCPointSet3& src,
  PCCPointSet3* dst,
  std::multimap<Vec3<double>, int32_t>& doubleQuantizedToOrigin);

//============================================================================
// Quantise the geometry of a point cloud, retaining duplicate points.
// Points in the @src point cloud are translated by -@offset, then quantised
// by a multiplicitive @scaleFactor with rounding.
//
// The destination and source point clouds may be the same object.
//
// NB: attributes are preserved

void quantizePositions(
  const float scaleFactor,
  const Vec3<int> offset,
  const Box3<int> clamp,
  const PCCPointSet3& src,
  PCCPointSet3* dst);

//============================================================================
// Clamp point co-ordinates in @cloud to @bbox, preserving attributes.

void clampVolume(Box3<double> bbox, PCCPointSet3* cloud);

//============================================================================
// Determine colour attribute values from a reference/source point cloud.
// For each point of the target p_t:
//  - Find the N_1 (1 < N_1) nearest neighbours in source to p_t and create
//    a set of points denoted by Ψ_1.
//  - Find the set of source points that p_t belongs to their set of N_2
//    nearest neighbours. Denote this set of points by Ψ_2.
//  - Compute the distance-weighted average of points in Ψ_1 and Ψ_2 by:
//        \bar{Ψ}_k = ∑_{q∈Ψ_k} c(q)/Δ(q,p_t)
//                    ----------------------- ,
//                    ∑_{q∈Ψ_k} 1/Δ(q,p_t)
//
// where Δ(a,b) denotes the Euclidian distance between the points a and b,
// and c(q) denotes the colour of point q.  Compute the average (or the
// weighted average with the number of points of each set as the weights)
// of \bar{Ψ}̅_1 and \bar{Ψ}̅_2 and transfer it to p_t.
//
// Differences in the scale and translation of the target and source point
// clouds, is handled according to:
//    posInTgt = (posInSrc - targetToSourceOffset) * sourceToTargetScaleFactor

bool recolourColour(
  const RecolourParams& params,
  const PCCPointSet3& source,
  double sourceToTargetScaleFactor,
  Vec3<double> targetToSourceOffset,
  PCCPointSet3& target);

//============================================================================
// Determine reflectance attribute values from a reference/source point cloud.
// For each point of the target p_t:
//  - Find the N_1 (1 < N_1) nearest neighbours in source to p_t and create
//    a set of points denoted by Ψ_1.
//  - Find the set of source points that p_t belongs to their set of N_2
//    nearest neighbours. Denote this set of points by Ψ_2.
//  - Compute the distance-weighted average of points in Ψ_1 and Ψ_2 by:
//        \bar{Ψ}_k = ∑_{q∈Ψ_k} c(q)/Δ(q,p_t)
//                    ----------------------- ,
//                    ∑_{q∈Ψ_k} 1/Δ(q,p_t)
//
// where Δ(a,b) denotes the Euclidian distance between the points a and b,
// and c(q) denotes the colour of point q.  Compute the average (or the
// weighted average with the number of points of each set as the weights)
// of \bar{Ψ}̅_1 and \bar{Ψ}̅_2 and transfer it to p_t.
//
// Differences in the scale and translation of the target and source point
// clouds, is handled according to:
//    posInTgt = (posInSrc - targetToSourceOffset) * sourceToTargetScaleFactor

bool recolourReflectance(
  const RecolourParams& cfg,
  const PCCPointSet3& source,
  double sourceToTargetScaleFactor,
  Vec3<double> targetToSourceOffset,
  PCCPointSet3& target);

//============================================================================
// Recolour attributes based on a source/reference point cloud.
//
// Differences in the scale and translation of the target and source point
// clouds, is handled according to:
//   posInTgt =
//     (posInSrc - targetToSourceOffset) * sourceToTargetScaleFactor - offset

int recolour(
  const RecolourParams& cfg,
  const PCCPointSet3& source,
  float sourceToTargetScaleFactor,
  Vec3<int> targetToSourceOffset,
  Vec3<int> offset,
  PCCPointSet3* target);

//============================================================================

void convertGbrToYCbCrBt709(PCCPointSet3&);
void convertYCbCrBt709ToGbr(PCCPointSet3&);

//============================================================================

}  // namespace pcc
