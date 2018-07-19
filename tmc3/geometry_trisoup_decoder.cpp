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

#include <cstdio>

#include "geometry_trisoup.h"

#include "PCCPointSetProcessing.h"
#include "geometry.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================

bool
operator<(const TrisoupSegment& s1, const TrisoupSegment& s2)
{
  // assert all quantities are at most 21 bits
  uint64_t s1startpos = (uint64_t(s1.startpos[0]) << 42)
    | (uint64_t(s1.startpos[1]) << 21) | s1.startpos[2];

  uint64_t s2startpos = (uint64_t(s2.startpos[0]) << 42)
    | (uint64_t(s2.startpos[1]) << 21) | s2.startpos[2];

  if (s1startpos < s2startpos)
    return true;

  if (s1startpos != s2startpos)
    return false;

  uint64_t s1endpos = (uint64_t(s1.endpos[0]) << 42)
    | (uint64_t(s1.endpos[1]) << 21) | s1.endpos[2];

  uint64_t s2endpos = (uint64_t(s2.endpos[0]) << 42)
    | (uint64_t(s2.endpos[1]) << 21) | s2.endpos[2];

  if (s1endpos < s2endpos)
    return true;

  if (s1endpos == s2endpos)
    if (s1.index < s2.index)  // stable sort
      return true;

  return false;
}

//============================================================================

void
decodeGeometryTrisoup(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  o3dgc::Arithmetic_Codec* arithmeticDecoder)
{
  // trisoup uses octree coding until reaching the triangulation level.
  pcc::ringbuf<PCCOctree3Node> nodes;
  decodeGeometryOctree(gps, gbh, pointCloud, arithmeticDecoder, &nodes);

  int blockWidth = 1 << (gps.trisoup_depth - gps.trisoup_triangle_level);

  uint32_t symbolCount;
  o3dgc::Adaptive_Bit_Model ctxTemp;
  o3dgc::Static_Bit_Model ctxBypass;

  // Decode segind from bitstream.
  symbolCount = arithmeticDecoder->ExpGolombDecode(0, ctxBypass, ctxTemp);

  o3dgc::Adaptive_Data_Model multiSymbolSegindModel0(256);
  std::vector<bool> segind;
  for (uint32_t i = 0; i < symbolCount; i++) {
    const uint8_t c = arithmeticDecoder->decode(multiSymbolSegindModel0);
    for (int b = 7; b >= 0; b--) {
      segind.push_back(!!(c & (1 << b)));
    }
  }

  // Decode vertices from bitstream.
  symbolCount = arithmeticDecoder->ExpGolombDecode(0, ctxBypass, ctxTemp);
  o3dgc::Adaptive_Data_Model multiSymbolVerticesModel0(blockWidth);
  std::vector<uint8_t> vertices;
  for (uint32_t i = 0; i < symbolCount; i++) {
    const uint8_t c = arithmeticDecoder->decode(multiSymbolVerticesModel0);
    vertices.push_back(c);
  }

  // Compute refinedVertices.
  decodeTrisoupCommon(nodes, segind, vertices, pointCloud, blockWidth);

  // Voxelize.
  int32_t maxval = (1 << gps.trisoup_depth) - 1;
  PCCBox3<int32_t> clampBox{{0, 0, 0}, {maxval, maxval, maxval}};
  quantizePositionsUniq(1.0, {0.0}, clampBox, pointCloud, &pointCloud);
}

//---------------------------------------------------------------------------
// Trisoup geometry decoding, at both encoder and decoder.
// Compute from leaves, segment indicators, and vertices
// a set of triangles, refine the triangles, and output their vertices.
//
// @param leaves  list of blocks containing the surface
// @param segind, indicators for edges of blocks if they intersect the surface
// @param vertices, locations of intersections

void
decodeTrisoupCommon(
  const ringbuf<PCCOctree3Node>& leaves,
  const std::vector<bool>& segind,
  const std::vector<uint8_t>& vertices,
  PCCPointSet3& pointCloud,
  const int defaultBlockWidth)
{
  // Put all leaves' sgements into a list.
  std::vector<TrisoupSegment> segments;
  for (int i = 0; i < leaves.size(); i++) {
    auto leaf = leaves[i];

    // Width of block.
    const uint32_t blockWidth =
      defaultBlockWidth;  // in future, may override with leaf blockWidth

    // Eight corners of block.
    const PCCVector3<uint32_t> pos000({0, 0, 0});
    const PCCVector3<uint32_t> posW00({blockWidth, 0, 0});
    const PCCVector3<uint32_t> pos0W0({0, blockWidth, 0});
    const PCCVector3<uint32_t> posWW0({blockWidth, blockWidth, 0});
    const PCCVector3<uint32_t> pos00W({0, 0, blockWidth});
    const PCCVector3<uint32_t> posW0W({blockWidth, 0, blockWidth});
    const PCCVector3<uint32_t> pos0WW({0, blockWidth, blockWidth});
    const PCCVector3<uint32_t> posWWW({blockWidth, blockWidth, blockWidth});

    // x: left to right; y: bottom to top; z: far to near
    segments.push_back(  // far bottom edge
      {leaf.pos + pos000, leaf.pos + posW00, 12 * i + 0, -1, -1});
    segments.push_back(  // far left edge
      {leaf.pos + pos000, leaf.pos + pos0W0, 12 * i + 1, -1, -1});
    segments.push_back(  // far top edge
      {leaf.pos + pos0W0, leaf.pos + posWW0, 12 * i + 2, -1, -1});
    segments.push_back(  // far right edge
      {leaf.pos + posW00, leaf.pos + posWW0, 12 * i + 3, -1, -1});
    segments.push_back(  // bottom left edge
      {leaf.pos + pos000, leaf.pos + pos00W, 12 * i + 4, -1, -1});
    segments.push_back(  // top left edge
      {leaf.pos + pos0W0, leaf.pos + pos0WW, 12 * i + 5, -1, -1});
    segments.push_back(  // top right edge
      {leaf.pos + posWW0, leaf.pos + posWWW, 12 * i + 6, -1, -1});
    segments.push_back(  // bottom right edge
      {leaf.pos + posW00, leaf.pos + posW0W, 12 * i + 7, -1, -1});
    segments.push_back(  // near bottom edge
      {leaf.pos + pos00W, leaf.pos + posW0W, 12 * i + 8, -1, -1});
    segments.push_back(  // near left edge
      {leaf.pos + pos00W, leaf.pos + pos0WW, 12 * i + 9, -1, -1});
    segments.push_back(  // near top edge
      {leaf.pos + pos0WW, leaf.pos + posWWW, 12 * i + 10, -1, -1});
    segments.push_back(  // near right edge
      {leaf.pos + posW0W, leaf.pos + posWWW, 12 * i + 11, -1, -1});
  }

  // Copy list of segments to another list to be sorted.
  std::vector<TrisoupSegment> sortedSegments;
  for (int i = 0; i < segments.size(); i++)
    sortedSegments.push_back(segments[i]);

  // Sort the list and find unique segments.
  std::sort(sortedSegments.begin(), sortedSegments.end());
  std::vector<TrisoupSegment> uniqueSegments;
  uniqueSegments.push_back(sortedSegments[0]);
  segments[sortedSegments[0].index].uniqueIndex = 0;
  for (int i = 1; i < sortedSegments.size(); i++) {
    if (
      uniqueSegments.back().startpos != sortedSegments[i].startpos
      || uniqueSegments.back().endpos != sortedSegments[i].endpos) {
      // sortedSegment[i] is different from uniqueSegments.back().
      // Start a new uniqueSegment.
      uniqueSegments.push_back(sortedSegments[i]);
    }
    segments[sortedSegments[i].index].uniqueIndex = uniqueSegments.size() - 1;
  }

  // Get vertex for each unique segment that intersects the surface.
  int vertexCount = 0;
  for (int i = 0; i < uniqueSegments.size(); i++) {
    if (segind[i]) {  // intersects the surface
      uniqueSegments[i].vertex = vertices[vertexCount++];
    } else {  // does not intersect the surface
      uniqueSegments[i].vertex = -1;
    }
  }

  // Copy vertices back to original (non-unique, non-sorted) segments.
  for (int i = 0; i < segments.size(); i++) {
    segments[i].vertex = uniqueSegments[segments[i].uniqueIndex].vertex;
  }

  // Create list of refined vertices, one leaf at a time.
  std::vector<PCCPoint3D> refinedVertices;
  for (int i = 0; i < leaves.size(); i++) {
    uint32_t blockWidth = 0;

    // Representation for a vertex in preparation for sorting.
    struct Vertex {
      PCCPoint3D pos;     // position of vertex
      double theta;       // angle of vertex when projected along dominant axis
      double tiebreaker;  // coordinate of vertex along dominant axis
      bool operator()(Vertex v1, Vertex v2)
      {
        if (v1.theta > v2.theta)
          return true;  // sort in decreasing order of theta
        if (v1.theta == v2.theta && v1.tiebreaker < v2.tiebreaker)
          return true;
        return false;
      }
    } vertex;

    // Find up to 12 vertices for this leaf.
    std::vector<Vertex> leafVertices;
    for (int j = 0; j < 12; j++) {
      TrisoupSegment& segment = segments[i * 12 + j];
      if (segment.vertex < 0)
        continue;  // skip segments that do not intersect the surface

      // Get distance along edge of vertex.
      // Vertex code is the index of the voxel along the edge of the block
      // of surface intersection./ Put decoded vertex at center of voxel,
      // unless voxel is first or last along the edge, in which case put the
      // decoded vertex at the start or endpoint of the segment.
      PCCVector3<uint32_t> direction = segment.endpos - segment.startpos;
      blockWidth =
        std::max(direction[0], std::max(direction[1], direction[2]));
      double distance;
      if (segment.vertex == 0)
        distance = 0.0;
      else if (segment.vertex == blockWidth - 1)
        distance = blockWidth;
      else  // 0 < segment.vertex < blockWidth-1
        distance = segment.vertex + 0.5;

      // Get 3D position of point of intersection.
      PCCVector3D point(
        segment.startpos[0], segment.startpos[1], segment.startpos[2]);
      if (direction[0] > 0)
        point[0] += distance;
      else if (direction[1] > 0)
        point[1] += distance;
      else  // direction[2] > 0
        point[2] += distance;

      // Add vertex to list of points.
      leafVertices.push_back({point, 0.0, 0.0});
    }

    // Skip leaves that have fewer than 3 vertices.
    if (leafVertices.size() < 3)
      continue;

    // Compute mean of leaf vertices.
    PCCPoint3D blockCentroid = 0;
    for (int j = 0; j < leafVertices.size(); j++) {
      blockCentroid += leafVertices[j].pos;
    }
    blockCentroid /= leafVertices.size();

    // Compute variance of each component of leaf vertices.
    PCCVector3D SS = 0;
    for (int j = 0; j < leafVertices.size(); j++) {
      PCCVector3D S = leafVertices[j].pos - blockCentroid;
      SS += {S[0] * S[0], S[1] * S[1], S[2] * S[2]};
    }

    // Dominant axis is the coordinate minimizing the variance.
    double minSS = SS[0];
    int dominantAxis = 0;
    for (int j = 1; j < 3; j++) {
      if (minSS > SS[j]) {
        minSS = SS[j];
        dominantAxis = j;
      }
    }

    // Project vertices along dominant axis (i.e., into YZ, XZ, or XY plane).
    // Sort projected vertices by decreasing angle in [-pi,+pi] around center
    // of block (i.e., clockwise from 9:00) breaking ties in angle by
    // increasing distance along the dominant axis.
    PCCVector3<uint32_t> bc = leaves[i].pos + (blockWidth / 2);
    PCCVector3D blockCenter = {double(bc[0]), double(bc[1]), double(bc[2])};
    for (int j = 0; j < leafVertices.size(); j++) {
      PCCVector3D S = leafVertices[j].pos - blockCenter;
      switch (dominantAxis) {
      case 0:  // dominant axis is X so project into YZ plane
        leafVertices[j].theta = atan2(S[2], S[1]);
        leafVertices[j].tiebreaker = S[0];
        break;
      case 1:  // dominant axis is Y so project into XZ plane
        leafVertices[j].theta = atan2(S[2], S[0]);
        leafVertices[j].tiebreaker = S[1];
        break;
      case 2:  // dominant axis is Z so project into XY plane
        leafVertices[j].theta = atan2(S[1], S[0]);
        leafVertices[j].tiebreaker = S[2];
        break;
      }
    }
    std::sort(leafVertices.begin(), leafVertices.end(), vertex);

    // Table of triangles that make up an n-gon.
    const int polyTriangles[][3] = {
      {0, 1, 2},                                                  // n = 3
      {0, 1, 2},   {2, 3, 0},                                     // n = 4
      {0, 1, 2},   {2, 3, 4}, {4, 0, 2},                          // n = 5
      {0, 1, 2},   {2, 3, 4}, {4, 5, 0},  {0, 2, 4},              // n = 6
      {0, 1, 2},   {2, 3, 4}, {4, 5, 6},  {6, 0, 2},  {2, 4, 6},  // n = 7
      {0, 1, 2},   {2, 3, 4}, {4, 5, 6},  {6, 7, 0},  {0, 2, 4},
      {4, 6, 0},  // n = 8
      {0, 1, 2},   {2, 3, 4}, {4, 5, 6},  {6, 7, 8},  {8, 0, 2},
      {2, 4, 6},   {6, 8, 2},  // n = 9
      {0, 1, 2},   {2, 3, 4}, {4, 5, 6},  {6, 7, 8},  {8, 9, 0},
      {0, 2, 4},   {4, 6, 8}, {8, 0, 4},  // n = 10
      {0, 1, 2},   {2, 3, 4}, {4, 5, 6},  {6, 7, 8},  {8, 9, 10},
      {10, 0, 2},  {2, 4, 6}, {6, 8, 10}, {10, 2, 6},  // n = 11
      {0, 1, 2},   {2, 3, 4}, {4, 5, 6},  {6, 7, 8},  {8, 9, 10},
      {10, 11, 0}, {0, 2, 4}, {4, 6, 8},  {8, 10, 0}, {0, 4, 8}  // n = 12
    };

    // Divide vertices into triangles according to table
    // and upsample each triangle by an upsamplingFactor.
    int triCount = leafVertices.size() - 2;
    int triStart = (triCount - 1) * triCount / 2;
    for (int triIndex = 0; triIndex < triCount; triIndex++) {
      int j0 = polyTriangles[triStart + triIndex][0];
      int j1 = polyTriangles[triStart + triIndex][1];
      int j2 = polyTriangles[triStart + triIndex][2];
      PCCVector3D v0 = leafVertices[j0].pos;
      PCCVector3D v1 = leafVertices[j1].pos;
      PCCVector3D v2 = leafVertices[j2].pos;
      PCCVector3D dv1 = v1 - v0;
      PCCVector3D dv2 = v2 - v0;
      int upsamplingFactor = blockWidth * 2;
      for (int k1 = 0; k1 <= upsamplingFactor; k1++) {
        for (int k2 = 0; k2 <= upsamplingFactor - k1; k2++) {
          PCCPoint3D v = v0 + dv1 * (double(k1) / upsamplingFactor)
            + dv2 * (double(k2) / upsamplingFactor);
          //// Clip to leaf boundaries.
          //for (int m = 0; m < 3; m++)
          //  v[m] = PCCClip<double>(v[m], leaves[i].pos[m], leaves[i].pos[m] + blockWidth - 0.5);
          // todo(pac?): why does subtracting 0.5 improve performance slightly?
          refinedVertices.push_back(v - 0.5);
        }
      }
    }
  }

  // Move list of points to pointCloud.
  pointCloud.resize(refinedVertices.size());
  for (int i = 0; i < refinedVertices.size(); i++)
    pointCloud[i] = refinedVertices[i];
}

//============================================================================

}  // namespace pcc
