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

#include "geometry_trisoup.h"

#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================

void
encodeGeometryTrisoup(
  const OctreeEncOpts& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders)
{
  // trisoup uses octree coding until reaching the triangulation level.
  pcc::ringbuf<PCCOctree3Node> nodes;
  encodeGeometryOctree(
    opt, gps, gbh, pointCloud, ctxtMemOctree, arithmeticEncoders, &nodes);

  // resume encoding with the last encoder
  auto arithmeticEncoder = arithmeticEncoders.back().get();

  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);

  // Determine segind and vertices.
  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  determineTrisoupVertices(nodes, segind, vertices, pointCloud, blockWidth);

  // Decode refinedVertices from segind and vertices.
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;

  // Decode vertices with certain sampling value
  int subsample = 1;
  if (gps.trisoup_sampling_value > 0) {
    subsample = gps.trisoup_sampling_value;
    decodeTrisoupCommon(
      nodes, segind, vertices, pointCloud, blockWidth, maxval, subsample);
  } else {
    int maxSubsample = 1 << gbh.trisoupNodeSizeLog2(gps);
    for (subsample = 1; subsample <= maxSubsample; subsample++) {
      decodeTrisoupCommon(
        nodes, segind, vertices, pointCloud, blockWidth, maxval, subsample);
      if (pointCloud.getPointCount() <= gbh.footer.geom_num_points_minus1 + 1)
        break;
    }
  }

  gbh.trisoup_sampling_value_minus1 = subsample - 1;
  gbh.num_unique_segments_minus1 = segind.size() - 1;
  gbh.num_unique_segments_bits_minus1 =
    numBits(gbh.num_unique_segments_minus1) - 1;
  assert(segind.size() > 0);

  // Encode segind to bitstream.
  int numVertices = 0;
  AdaptiveBitModel ctxTempSeg;
  for (int i = 0; i <= gbh.num_unique_segments_minus1; i++) {
    arithmeticEncoder->encode((int)segind[i], ctxTempSeg);
    numVertices += segind[i];
  }

  // Encode vertices to bitstream.
  assert(numVertices == vertices.size());
  for (int i = 0; i < numVertices; i++) {
    auto vertex = vertices[i];
    for (int b = gbh.trisoupNodeSizeLog2(gps) - 1; b >= 0; b--)
      arithmeticEncoder->encode((vertex >> b) & 1);
  }
}

//---------------------------------------------------------------------------
// Determine where the surface crosses each leaf
// (i.e., determine the segment indicators and vertices)
// from the set of leaves and the points in each leaf.
//
// @param leaves, list of blocks containing the surface
// @param segind, indicators for edges of blocks if they intersect the surface
// @param vertices, locations of intersections

void
determineTrisoupVertices(
  const ringbuf<PCCOctree3Node>& leaves,
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  const PCCPointSet3& pointCloud,
  const int defaultBlockWidth)
{
  // Put all leaves' edges into a list.
  std::vector<TrisoupSegmentEnc> segments;
  for (int i = 0; i < leaves.size(); i++) {
    const auto& leaf = leaves[i];

    // Width of block.
    // in future, may override with leaf blockWidth
    const int32_t blockWidth = defaultBlockWidth;

    // Eight corners of block.
    const Vec3<int32_t> pos000({0, 0, 0});
    const Vec3<int32_t> posW00({blockWidth, 0, 0});
    const Vec3<int32_t> pos0W0({0, blockWidth, 0});
    const Vec3<int32_t> posWW0({blockWidth, blockWidth, 0});
    const Vec3<int32_t> pos00W({0, 0, blockWidth});
    const Vec3<int32_t> posW0W({blockWidth, 0, blockWidth});
    const Vec3<int32_t> pos0WW({0, blockWidth, blockWidth});
    const Vec3<int32_t> posWWW({blockWidth, blockWidth, blockWidth});

    // x: left to right; y: bottom to top; z: far to near
    TrisoupSegmentEnc seg000W00 =  // far bottom edge
      {leaf.pos + pos000, leaf.pos + posW00, 12 * i + 0, -1, -1, 0, 0};
    TrisoupSegmentEnc seg0000W0 =  // far left edge
      {leaf.pos + pos000, leaf.pos + pos0W0, 12 * i + 1, -1, -1, 0, 0};
    TrisoupSegmentEnc seg0W0WW0 =  // far top edge
      {leaf.pos + pos0W0, leaf.pos + posWW0, 12 * i + 2, -1, -1, 0, 0};
    TrisoupSegmentEnc segW00WW0 =  // far right edge
      {leaf.pos + posW00, leaf.pos + posWW0, 12 * i + 3, -1, -1, 0, 0};
    TrisoupSegmentEnc seg00000W =  // bottom left edge
      {leaf.pos + pos000, leaf.pos + pos00W, 12 * i + 4, -1, -1, 0, 0};
    TrisoupSegmentEnc seg0W00WW =  // top left edge
      {leaf.pos + pos0W0, leaf.pos + pos0WW, 12 * i + 5, -1, -1, 0, 0};
    TrisoupSegmentEnc segWW0WWW =  // top right edge
      {leaf.pos + posWW0, leaf.pos + posWWW, 12 * i + 6, -1, -1, 0, 0};
    TrisoupSegmentEnc segW00W0W =  // bottom right edge
      {leaf.pos + posW00, leaf.pos + posW0W, 12 * i + 7, -1, -1, 0, 0};
    TrisoupSegmentEnc seg00WW0W =  // near bottom edge
      {leaf.pos + pos00W, leaf.pos + posW0W, 12 * i + 8, -1, -1, 0, 0};
    TrisoupSegmentEnc seg00W0WW =  // near left edge
      {leaf.pos + pos00W, leaf.pos + pos0WW, 12 * i + 9, -1, -1, 0, 0};
    TrisoupSegmentEnc seg0WWWWW =  // near top edge
      {leaf.pos + pos0WW, leaf.pos + posWWW, 12 * i + 10, -1, -1, 0, 0};
    TrisoupSegmentEnc segW0WWWW =  // near right edge
      {leaf.pos + posW0W, leaf.pos + posWWW, 12 * i + 11, -1, -1, 0, 0};

    // Each voxel votes for a position along each edge it is close to
    for (int j = leaf.start; j < leaf.end; j++) {
      Vec3<int> voxel;
      voxel.x() = int(pointCloud[j].x()) - leaf.pos.x();
      voxel.y() = int(pointCloud[j].y()) - leaf.pos.y();
      voxel.z() = int(pointCloud[j].z()) - leaf.pos.z();
      // parameter indicating threshold of how close voxels must be to edge
      // to be relevant
      int tmin = 1;
      int tmax = blockWidth - tmin - 1;
      if (voxel.y() < tmin && voxel.z() < tmin) {
        seg000W00.count++;
        seg000W00.distanceSum += voxel.x();
      }  // far bottom edge
      if (voxel.x() < tmin && voxel.z() < tmin) {
        seg0000W0.count++;
        seg0000W0.distanceSum += voxel.y();
      }  // far left edge
      if (voxel.y() > tmax && voxel.z() < tmin) {
        seg0W0WW0.count++;
        seg0W0WW0.distanceSum += voxel.x();
      }  // far top edge
      if (voxel.x() > tmax && voxel.z() < tmin) {
        segW00WW0.count++;
        segW00WW0.distanceSum += voxel.y();
      }  // far right edge
      if (voxel.x() < tmin && voxel.y() < tmin) {
        seg00000W.count++;
        seg00000W.distanceSum += voxel.z();
      }  // bottom left edge
      if (voxel.x() < tmin && voxel.y() > tmax) {
        seg0W00WW.count++;
        seg0W00WW.distanceSum += voxel.z();
      }  // top left edge
      if (voxel.x() > tmax && voxel.y() > tmax) {
        segWW0WWW.count++;
        segWW0WWW.distanceSum += voxel.z();
      }  // top right edge
      if (voxel.x() > tmax && voxel.y() < tmin) {
        segW00W0W.count++;
        segW00W0W.distanceSum += voxel.z();
      }  // bottom right edge
      if (voxel.y() < tmin && voxel.z() > tmax) {
        seg00WW0W.count++;
        seg00WW0W.distanceSum += voxel.x();
      }  // near bottom edge
      if (voxel.x() < tmin && voxel.z() > tmax) {
        seg00W0WW.count++;
        seg00W0WW.distanceSum += voxel.y();
      }  // near left edge
      if (voxel.y() > tmax && voxel.z() > tmax) {
        seg0WWWWW.count++;
        seg0WWWWW.distanceSum += voxel.x();
      }  // near top edge
      if (voxel.x() > tmax && voxel.z() > tmax) {
        segW0WWWW.count++;
        segW0WWWW.distanceSum += voxel.y();
      }  // near right edge
    }

    // Push segments onto list.
    segments.push_back(seg000W00);  // far bottom edge
    segments.push_back(seg0000W0);  // far left edge
    segments.push_back(seg0W0WW0);  // far top edge
    segments.push_back(segW00WW0);  // far right edge
    segments.push_back(seg00000W);  // bottom left edge
    segments.push_back(seg0W00WW);  // top left edge
    segments.push_back(segWW0WWW);  // top right edge
    segments.push_back(segW00W0W);  // bottom right edge
    segments.push_back(seg00WW0W);  // near bottom edge
    segments.push_back(seg00W0WW);  // near left edge
    segments.push_back(seg0WWWWW);  // near top edge
    segments.push_back(segW0WWWW);  // near right edge
  }

  // Copy list of segments to another list to be sorted.
  std::vector<TrisoupSegmentEnc> sortedSegments;
  for (int i = 0; i < segments.size(); i++)
    sortedSegments.push_back(segments[i]);

  // Sort the list and find unique segments.
  std::sort(sortedSegments.begin(), sortedSegments.end());
  std::vector<TrisoupSegmentEnc> uniqueSegments;
  uniqueSegments.push_back(sortedSegments[0]);
  segments[sortedSegments[0].index].uniqueIndex = 0;
  for (int i = 1; i < sortedSegments.size(); i++) {
    if (
      uniqueSegments.back().startpos != sortedSegments[i].startpos
      || uniqueSegments.back().endpos != sortedSegments[i].endpos) {
      // sortedSegment[i] is different from uniqueSegments.back().
      // Start a new uniqueSegment.
      uniqueSegments.push_back(sortedSegments[i]);  // unique segment
    } else {
      // sortedSegment[i] is the same as uniqueSegments.back().
      // Accumulate into uniqueSegment.back().
      uniqueSegments.back().count += sortedSegments[i].count;
      uniqueSegments.back().distanceSum += sortedSegments[i].distanceSum;
    }
    segments[sortedSegments[i].index].uniqueIndex = uniqueSegments.size() - 1;
  }

  // Compute vertex for each unique segment that intersects the surface.
  for (int i = 0; i < uniqueSegments.size(); i++) {
    if (uniqueSegments[i].count > 0) {  // intersects the surface
      segind.push_back(true);
      uint8_t vertex =
        (uniqueSegments[i].distanceSum + uniqueSegments[i].count / 2)
        / uniqueSegments[i].count;
      vertices.push_back(vertex);
      uniqueSegments[i].vertex = vertex;
    } else {  // does not intersect the surface
      segind.push_back(false);
      uniqueSegments[i].vertex = -1;
    }
  }

  // Copy vertices back to original (non-unique, non-sorted) segments.
  for (int i = 0; i < segments.size(); i++) {
    segments[i].vertex = uniqueSegments[segments[i].uniqueIndex].vertex;
  }
}

//============================================================================

}  // namespace pcc
