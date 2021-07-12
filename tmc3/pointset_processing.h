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
#include "hls.h"

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
// Represents a quatized point cloud with index mappings to source positions.

struct SrcMappedPointSet {
  // A subsampled or quantised pointcloud generated from a source.
  PCCPointSet3 cloud;

  // Maps cloud indexes to the generating source index.
  std::vector<int> idxToSrcIdx;

  // Linked lists of source indexes that map to the same generated position.
  // Each element is the index of the next element in the chain.
  // The end of a chain is indicated by: srcIdxDupList[i] == i.
  std::vector<int> srcIdxDupList;
};

//============================================================================
// Subsample a point cloud, retaining unique points only.
// Uniqueness is assessed by quantising each position by a multiplicative
// @sampleScale.  Output points are quantised by @quantScale with rounding,
// and translated by -@offset.
//
// NB: attributes are not processed.

SrcMappedPointSet samplePositionsUniq(
  float sampleScale,
  float quantScale,
  Vec3<int> offset,
  const PCCPointSet3& src);

//============================================================================
// Quantise the geometry of a point cloud, retaining unique points only.
// Points in the @src point cloud are translated by -@offset, quantised by a
// multiplicitive @scaleFactor with rounding, then clamped to @clamp.
//
// NB: attributes are not processed.

SrcMappedPointSet quantizePositionsUniq(
  const float scaleFactor,
  const Vec3<int> offset,
  const Box3<int> clamp,
  const PCCPointSet3& src);

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

void clampVolume(Box3<int32_t> bbox, PCCPointSet3* cloud);

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
//    posInTgt = posInSrc * sourceToTargetScaleFactor - targetToSourceOffset

bool recolourColour(
  const AttributeDescription& desc,
  const RecolourParams& params,
  const PCCPointSet3& source,
  double sourceToTargetScaleFactor,
  point_t targetToSourceOffset,
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
//    posInTgt = posInSrc * sourceToTargetScaleFactor - targetToSourceOffset

bool recolourReflectance(
  const AttributeDescription& desc,
  const RecolourParams& cfg,
  const PCCPointSet3& source,
  double sourceToTargetScaleFactor,
  point_t targetToSourceOffset,
  PCCPointSet3& target);

//============================================================================
// Recolour attributes based on a source/reference point cloud.
//
// Differences in the scale and translation of the target and source point
// clouds, is handled according to:
//   posInTgt = posInSrc * sourceToTargetScaleFactor - tgtToSrcOffset

int recolour(
  const AttributeDescription& desc,
  const RecolourParams& cfg,
  const PCCPointSet3& source,
  float sourceToTargetScaleFactor,
  point_t tgtToSrcOffset,
  PCCPointSet3* target);

//============================================================================

void convertGbrToYCbCrBt709(PCCPointSet3&);
void convertYCbCrBt709ToGbr(PCCPointSet3&);

void convertGbrToYCgCoR(int bitDepth, PCCPointSet3&);
void convertYCgCoRToGbr(int bitDepth, PCCPointSet3&);

//============================================================================

// Generate index order sorted by azimuth angle.
std::vector<int32_t> orderByAzimuth(
  PCCPointSet3&, int start, int end, int numDigit, Vec3<int32_t> origin);

std::vector<int32_t>
orderByRadius(PCCPointSet3&, int start, int end, Vec3<int32_t> origin);

std::vector<int32_t> orderByLaserAngle(
  PCCPointSet3&, int start, int end, int numDigit, Vec3<int32_t> origin);

// Sorts points according to azimuth angle.
void sortByAzimuth(
  PCCPointSet3&,
  int start,
  int end,
  double recipBinWidth,
  Vec3<int32_t> origin = 0);

void sortByRadius(PCCPointSet3&, int start, int end, Vec3<int32_t> origin = 0);

void sortByLaserAngle(
  PCCPointSet3&,
  int start,
  int end,
  double recipBinWidth,
  Vec3<int32_t> origin = 0);

//============================================================================

}  // namespace pcc
