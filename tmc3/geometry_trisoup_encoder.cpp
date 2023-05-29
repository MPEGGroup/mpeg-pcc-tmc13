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

#include <numeric>
#include <algorithm>
#include <unordered_set>
#include "geometry_trisoup.h"

#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================

void
encodeGeometryTrisoup(
  const TrisoupEncOpts& opt,
  const OctreeEncOpts& optOctree,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointCloudPadding,
  GeometryOctreeContexts& ctxtMemOctree,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  const InterGeomEncOpts& interParams)
{
  // trisoup uses octree coding until reaching the triangulation level.
  pcc::ringbuf<PCCOctree3Node> nodes;

  BiPredictionEncodeParams biPredEncodeParams;
  
  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);

  // Partition padding into nodes
  Box3<int32_t> box = pointCloud.computeBoundingBox();
  Box3<int32_t> boxP = pointCloudPadding.computeBoundingBox();

  Vec3<int32_t> min = {0};
  Vec3<int32_t> max = {0};

  max[0] += (box.max[0] % blockWidth == 0 ? box.max[0] + 16 : box.max[0] + (blockWidth - box.max[0] % blockWidth));
  max[1] += (box.max[1] % blockWidth == 0 ? box.max[1] + 16 : box.max[1] + (blockWidth - box.max[1] % blockWidth));
  max[2] += (box.max[2] % blockWidth == 0 ? box.max[2] + 16 : box.max[2] + (blockWidth - box.max[2] % blockWidth));

  Box3<int32_t> originalBox = Box3<int32_t>(min, max);

  std::vector<PCCOctree3Node> nodesPadded;
  std::vector<int> indices(pointCloudPadding.getPointCount());
  if (pointCloudPadding.getPointCount() != 0) {
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<Vec3<int>> mapping(pointCloudPadding.getPointCount());
    for (int i = 0; i < pointCloudPadding.getPointCount(); i++) {
      mapping[i][0] = (pointCloudPadding[i][0] < 0 ? pointCloudPadding[i][0]/blockWidth - 1 : pointCloudPadding[i][0]/blockWidth) * blockWidth;
      mapping[i][1] = (pointCloudPadding[i][1] < 0 ? pointCloudPadding[i][1]/blockWidth - 1 : pointCloudPadding[i][1]/blockWidth) * blockWidth;
      mapping[i][2] = (pointCloudPadding[i][2] < 0 ? pointCloudPadding[i][2]/blockWidth - 1 : pointCloudPadding[i][2]/blockWidth) * blockWidth;
    }
    std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool {return mapping[A] < mapping[B];});

    Vec3<int> v = mapping[indices[0]];
    PCCOctree3Node n;
    n.pos = mapping[indices[0]];
    n.start = 0;
    for (int i = 1; i < pointCloudPadding.getPointCount(); i++){
      if (v != mapping[indices[i]]){
        v = mapping[indices[i]];
        n.end = i; // -1 ?
        nodesPadded.push_back(n);
        n.pos = mapping[indices[i]];
        n.start = i;
      }
      if (i == pointCloudPadding.getPointCount() - 1){
        n.end = i;
        nodesPadded.push_back(n);
      }
    }
  }
  // End: Partition padding into nodes

  encodeGeometryOctree(
    optOctree, gps, gbh, pointCloud, ctxtMemOctree, arithmeticEncoders, &nodes,
    refFrame, sps, interParams, biPredEncodeParams);

  // resume encoding with the last encoder
  pcc::EntropyEncoder* arithmeticEncoder = arithmeticEncoders.back().get();

  const int maxVertexPrecisionLog2 = gbh.trisoup_vertex_quantization_bits
    ? gbh.trisoup_vertex_quantization_bits
    : gbh.trisoupNodeSizeLog2(gps);
  const int bitDropped =
    std::max(0, gbh.trisoupNodeSizeLog2(gps) - maxVertexPrecisionLog2);
  const bool isCentroidDriftActivated =
    gbh.trisoup_centroid_vertex_residual_flag;

  // Determine vertices
  std::cout << "Number of points = " << pointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes = " << nodes.size() << "\n";
  int distanceSearchEncoder = 1;
  if (opt.improvedVertexDetermination) {
    float estimatedSampling = float(nodes.size());
    estimatedSampling /= pointCloud.getPointCount();
    estimatedSampling = std::sqrt(estimatedSampling);
    estimatedSampling *= blockWidth;
    estimatedSampling = std::max(1.f, estimatedSampling);
    std::cout << "Estimation of sampling = " << estimatedSampling << "\n";

    distanceSearchEncoder = (1 << std::max(0, bitDropped - 2)) - 1;
    distanceSearchEncoder += int(std::round(estimatedSampling + 0.1f));
    distanceSearchEncoder = std::max(1, std::min(8, distanceSearchEncoder));
    std::cout << "distanceSearchEncoder = " << distanceSearchEncoder << "\n";
  }

  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  determineTrisoupVertices(
    nodes, segind, vertices, pointCloud,
    gps, gbh,
    blockWidth, bitDropped,
    distanceSearchEncoder, nodesPadded, pointCloudPadding, indices, originalBox);

  // Determine neighbours
  std::vector<uint16_t> neighbNodes;
  std::vector<std::array<int, 18>> edgePattern;
  determineTrisoupNeighbours(nodes, neighbNodes, edgePattern, blockWidth);

  gbh.num_unique_segments_minus1 = segind.size() - 1;
  gbh.num_unique_segments_bits_minus1 =
    numBits(gbh.num_unique_segments_minus1) - 1;

  // Encode vertex presence and position into bitstream
  assert(segind.size() > 0);
  encodeTrisoupVertices(segind, vertices, neighbNodes, edgePattern, bitDropped, gps, gbh, arithmeticEncoder);

  // Decode vertices with certain sampling value
  bool haloFlag = gbh.trisoup_halo_flag;
  bool adaptiveHaloFlag = gbh.trisoup_adaptive_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;

  PCCPointSet3 recPointCloud;
  recPointCloud.addRemoveAttributes(pointCloud);

  std::vector<CentroidDrift> drifts;
  int subsample = 1;
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;
  std::cout << "Loop on sampling for max "
            << (gbh.footer.geom_num_points_minus1 + 1) << " points \n";
  if (gps.trisoup_sampling_value > 0) {
    subsample = gps.trisoup_sampling_value;
    decodeTrisoupCommon(
      nodes, segind, vertices, drifts, pointCloud, recPointCloud,
      gps, gbh, blockWidth,
      maxval, subsample, bitDropped, isCentroidDriftActivated, false,
      haloFlag, adaptiveHaloFlag, fineRayFlag, NULL);
    std::cout << "Sub-sampling " << subsample << " gives "
              << recPointCloud.getPointCount() << " points \n";
  } else {
    int maxSubsample = 1 << gbh.trisoupNodeSizeLog2(gps);
    for (subsample = 1; subsample <= maxSubsample; subsample++) {
      decodeTrisoupCommon(
        nodes, segind, vertices, drifts, pointCloud, recPointCloud,
        gps, gbh, blockWidth,
        maxval, subsample, bitDropped, isCentroidDriftActivated, false,
        haloFlag, adaptiveHaloFlag, fineRayFlag, NULL);
      std::cout << "Sub-sampling " << subsample << " gives "
                << recPointCloud.getPointCount() << " points \n";
      if (recPointCloud.getPointCount() <= gbh.footer.geom_num_points_minus1 + 1)
        break;
    }
  }

  pointCloud.resize(0);
  pointCloud = std::move(recPointCloud);

  gbh.trisoup_sampling_value_minus1 = subsample - 1;

  // encoder centroid residua into bitstream
  if (isCentroidDriftActivated)
    encodeTrisoupCentroidResidue(drifts, arithmeticEncoder);
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
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const int defaultBlockWidth,
  const int bitDropped,
  int distanceSearchEncoder,
  const std::vector<PCCOctree3Node>& nodesPadded,
  const PCCPointSet3& pointCloudPadding,
  std::vector<int> indices,
  Box3<int32_t> originalBox)
{
  Box3<int32_t> sliceBB;
  sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
  sliceBB.max = sliceBB.min + ( gbh.slice_bb_width << gbh.slice_bb_width_log2_scale );
  std::cout << "SliceBB: " << sliceBB << std::endl;

  // Put all leaves' edges into a list.
  std::vector<TrisoupSegmentEnc> segments;
  std::unordered_set<TrisoupSegmentEnc, SegmentHash> lookup;
  segments.reserve(12 * (leaves.size()) + 12 * (nodesPadded.size()));
  for (int i = 0; i < leaves.size() + nodesPadded.size(); i++) {
    PCCOctree3Node leaf =
      i < leaves.size() ? leaves[i] : nodesPadded[i - leaves.size()];

    // Width of block.
    // in future, may override with leaf blockWidth
    const int32_t blockWidth = defaultBlockWidth;

    Vec3<int32_t> newP, newW, corner[8];
    nonCubicNode( gps, gbh, leaf.pos, defaultBlockWidth, sliceBB, newP, newW, corner );

    // x: left to right; y: bottom to top; z: far to near
    TrisoupSegmentEnc seg000W00 =  // far bottom edge
      {newP + corner[POS_000], newP + corner[POS_W00], 12 * i + 0, -1, -1, 0, 0, 0, 0}; // direction x
    TrisoupSegmentEnc seg0000W0 =  // far left edge
      {newP + corner[POS_000], newP + corner[POS_0W0], 12 * i + 1, -1, -1, 0, 0, 0, 0}; // direction y
    TrisoupSegmentEnc seg0W0WW0 =  // far top edge
      {newP + corner[POS_0W0], newP + corner[POS_WW0], 12 * i + 2, -1, -1, 0, 0, 0, 0}; // direction x
    TrisoupSegmentEnc segW00WW0 =  // far right edge
      {newP + corner[POS_W00], newP + corner[POS_WW0], 12 * i + 3, -1, -1, 0, 0, 0, 0}; // direction y
    TrisoupSegmentEnc seg00000W =  // bottom left edge
      {newP + corner[POS_000], newP + corner[POS_00W], 12 * i + 4, -1, -1, 0, 0, 0, 0}; // direction z
    TrisoupSegmentEnc seg0W00WW =  // top left edge
      {newP + corner[POS_0W0], newP + corner[POS_0WW], 12 * i + 5, -1, -1, 0, 0, 0, 0}; // direction z
    TrisoupSegmentEnc segWW0WWW =  // top right edge
      {newP + corner[POS_WW0], newP + corner[POS_WWW], 12 * i + 6, -1, -1, 0, 0, 0, 0}; // direction z
    TrisoupSegmentEnc segW00W0W =  // bottom right edge
      {newP + corner[POS_W00], newP + corner[POS_W0W], 12 * i + 7, -1, -1, 0, 0, 0, 0}; // direction z
    TrisoupSegmentEnc seg00WW0W =  // near bottom edge
      {newP + corner[POS_00W], newP + corner[POS_W0W], 12 * i + 8, -1, -1, 0, 0, 0, 0}; // direction x
    TrisoupSegmentEnc seg00W0WW =  // near left edge
      {newP + corner[POS_00W], newP + corner[POS_0WW], 12 * i + 9, -1, -1, 0, 0, 0, 0}; // direction y
    TrisoupSegmentEnc seg0WWWWW =  // near top edge
      {newP + corner[POS_0WW], newP + corner[POS_WWW], 12 * i + 10, -1, -1, 0, 0, 0, 0}; // direction x
    TrisoupSegmentEnc segW0WWWW =  // near right edge
      {newP + corner[POS_W0W], newP + corner[POS_WWW], 12 * i + 11, -1, -1, 0, 0, 0, 0}; // direction y

    // Each voxel votes for a position along each edge it is close to
    const int tmin = 1;
    const Vec3<int> tmax( newW.x() - tmin - 1,
                          newW.y() - tmin - 1,
                          newW.z() - tmin - 1 );
    const int tmin2 = distanceSearchEncoder;
    const Vec3<int> tmax2( newW.x() - tmin2 - 1,
                           newW.y() - tmin2 - 1,
                           newW.z() - tmin2 - 1 );
    Vec3<int> voxel;
    for (int j = leaf.start; j < leaf.end; j++) {
      if (i < leaves.size()){
        voxel = pointCloud[j] - newP;
      } else {
        voxel = pointCloudPadding[indices[j]] - newP;
      }
      // parameter indicating threshold of how close voxels must be to edge
      // ----------- 1 -------------------
      // to be relevant
      if (voxel[1] < tmin && voxel[2] < tmin) {
        seg000W00.count++;
        seg000W00.distanceSum += voxel[0];
      }  // far bottom edge
      if (voxel[0] < tmin && voxel[2] < tmin) {
        seg0000W0.count++;
        seg0000W0.distanceSum += voxel[1];
      }  // far left edge
      if (voxel[1] > tmax.y() && voxel[2] < tmin) {
        seg0W0WW0.count++;
        seg0W0WW0.distanceSum += voxel[0];
      }  // far top edge
      if (voxel[0] > tmax.x() && voxel[2] < tmin) {
        segW00WW0.count++;
        segW00WW0.distanceSum += voxel[1];
      }  // far right edge
      if (voxel[0] < tmin && voxel[1] < tmin) {
        seg00000W.count++;
        seg00000W.distanceSum += voxel[2];
      }  // bottom left edge
      if (voxel[0] < tmin && voxel[1] > tmax.y()) {
        seg0W00WW.count++;
        seg0W00WW.distanceSum += voxel[2];
      }  // top left edge
      if (voxel[0] > tmax.x() && voxel[1] > tmax.y()) {
        segWW0WWW.count++;
        segWW0WWW.distanceSum += voxel[2];
      }  // top right edge
      if (voxel[0] > tmax.x() && voxel[1] < tmin) {
        segW00W0W.count++;
        segW00W0W.distanceSum += voxel[2];
      }  // bottom right edge
      if (voxel[1] < tmin && voxel[2] > tmax.z()) {
        seg00WW0W.count++;
        seg00WW0W.distanceSum += voxel[0];
      }  // near bottom edge
      if (voxel[0] < tmin && voxel[2] > tmax.z()) {
        seg00W0WW.count++;
        seg00W0WW.distanceSum += voxel[1];
      }  // near left edge
      if (voxel[1] > tmax.y() && voxel[2] > tmax.z()) {
        seg0WWWWW.count++;
        seg0WWWWW.distanceSum += voxel[0];
      }  // near top edge
      if (voxel[0] > tmax.x() && voxel[2] > tmax.z()) {
        segW0WWWW.count++;
        segW0WWWW.distanceSum += voxel[1];
      }  // near right edge

      // parameter indicating threshold of how close voxels must be to edge
      // ----------- 2 -------------------
      // to be relevant
      if (voxel[1] < tmin2 && voxel[2] < tmin2) {
        seg000W00.count2++;
        seg000W00.distanceSum2 += voxel[0];
      }  // far bottom edge
      if (voxel[0] < tmin2 && voxel[2] < tmin2) {
        seg0000W0.count2++;
        seg0000W0.distanceSum2 += voxel[1];
      }  // far left edge
      if (voxel[1] > tmax2.y() && voxel[2] < tmin2) {
        seg0W0WW0.count2++;
        seg0W0WW0.distanceSum2 += voxel[0];
      }  // far top edge
      if (voxel[0] > tmax2.x() && voxel[2] < tmin2) {
        segW00WW0.count2++;
        segW00WW0.distanceSum2 += voxel[1];
      }  // far right edge
      if (voxel[0] < tmin2 && voxel[1] < tmin2) {
        seg00000W.count2++;
        seg00000W.distanceSum2 += voxel[2];
      }  // bottom left edge
      if (voxel[0] < tmin2 && voxel[1] > tmax2.y()) {
        seg0W00WW.count2++;
        seg0W00WW.distanceSum2 += voxel[2];
      }  // top left edge
      if (voxel[0] > tmax2.x() && voxel[1] > tmax2.y()) {
        segWW0WWW.count2++;
        segWW0WWW.distanceSum2 += voxel[2];
      }  // top right edge
      if (voxel[0] > tmax2.x() && voxel[1] < tmin2) {
        segW00W0W.count2++;
        segW00W0W.distanceSum2 += voxel[2];
      }  // bottom right edge
      if (voxel[1] < tmin2 && voxel[2] > tmax2.z()) {
        seg00WW0W.count2++;
        seg00WW0W.distanceSum2 += voxel[0];
      }  // near bottom edge
      if (voxel[0] < tmin2 && voxel[2] > tmax2.z()) {
        seg00W0WW.count2++;
        seg00W0WW.distanceSum2 += voxel[1];
      }  // near left edge
      if (voxel[1] > tmax2.y() && voxel[2] > tmax2.z()) {
        seg0WWWWW.count2++;
        seg0WWWWW.distanceSum2 += voxel[0];
      }  // near top edge
      if (voxel[0] > tmax2.x() && voxel[2] > tmax2.z()) {
        segW0WWWW.count2++;
        segW0WWWW.distanceSum2 += voxel[1];
      }  // near right edge
    }

    // Push segments onto list.
    if (i < leaves.size()){
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

      lookup.insert(seg000W00);  // far bottom edge
      lookup.insert(seg0000W0);  // far left edge
      lookup.insert(seg0W0WW0);  // far top edge
      lookup.insert(segW00WW0);  // far right edge
      lookup.insert(seg00000W);  // bottom left edge
      lookup.insert(seg0W00WW);  // top left edge
      lookup.insert(segWW0WWW);  // top right edge
      lookup.insert(segW00W0W);  // bottom right edge
      lookup.insert(seg00WW0W);  // near bottom edge
      lookup.insert(seg00W0WW);  // near left edge
      lookup.insert(seg0WWWWW);  // near top edge
      lookup.insert(segW0WWWW);  // near right edge
    } else {
      if (lookup.find(seg000W00) != lookup.end())
        segments.push_back(seg000W00);  // far bottom edge
      if (lookup.find(seg0000W0) != lookup.end())
        segments.push_back(seg0000W0);  // far left edge
      if (lookup.find(seg0W0WW0) != lookup.end())
        segments.push_back(seg0W0WW0);  // far top edge
      if (lookup.find(segW00WW0) != lookup.end())
        segments.push_back(segW00WW0);  // far right edge
      if (lookup.find(seg00000W) != lookup.end())
        segments.push_back(seg00000W);  // bottom left edge
      if (lookup.find(seg0W00WW) != lookup.end())
        segments.push_back(seg0W00WW);  // top left edge
      if (lookup.find(segWW0WWW) != lookup.end())
        segments.push_back(segWW0WWW);  // top right edge
      if (lookup.find(segW00W0W) != lookup.end())
        segments.push_back(segW00W0W);  // bottom right edge
      if (lookup.find(seg00WW0W) != lookup.end())
        segments.push_back(seg00WW0W);  // near bottom edge
      if (lookup.find(seg00W0WW) != lookup.end())
        segments.push_back(seg00W0WW);  // near left edge
      if (lookup.find(seg0WWWWW) != lookup.end())
        segments.push_back(seg0WWWWW);  // near top edge
      if (lookup.find(segW0WWWW) != lookup.end())
        segments.push_back(segW0WWWW);  // near right edge
    }
  }

  // Sort the list and find unique segments.
  std::sort(segments.begin(), segments.end());

  TrisoupSegmentEnc localSegment = segments[0];
  auto it = segments.begin() + 1;
  int i = 0;
  for (; it != segments.end(); it++) {
    if (
      localSegment.startpos != it->startpos
      || localSegment.endpos != it->endpos) {
      // Segment[i] is different from localSegment
      // Start a new uniqueSegment.
      segind.push_back(localSegment.count > 0 || localSegment.count2 > 1);
      if (segind.back()) {  // intersects the surface
        int temp = ((2 * localSegment.distanceSum + localSegment.distanceSum2)
                    << (10 - bitDropped))
          / (2 * localSegment.count + localSegment.count2);
        int8_t vertex = (temp + (1 << 9 - bitDropped)) >> 10;
        vertices.push_back(vertex);
      }
      localSegment = *it;  // unique segment
    } else {
      // Segment[i] is the same as localSegment
      // Accumulate
      localSegment.count += it->count;
      localSegment.distanceSum += it->distanceSum;
      localSegment.count2 += it->count2;
      localSegment.distanceSum2 += it->distanceSum2;
    }
  }
  segind.push_back(localSegment.count > 0 || localSegment.count2 > 1);
  if (segind.back()) {  // intersects the surface
    int temp = ((2 * localSegment.distanceSum + localSegment.distanceSum2)
                << (10 - bitDropped))
      / (2 * localSegment.count + localSegment.count2);
    int8_t vertex = (temp + (1 << 9 - bitDropped)) >> 10;
    vertices.push_back(vertex);
  }
}

//------------------------------------------------------------------------------------
void
encodeTrisoupVertices(
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  int bitDropped,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder)
{
  const int nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;
  const int max2bits = nbitsVertices > 1 ? 3 : 1;
  const int mid2bits = nbitsVertices > 1 ? 2 : 1;

  int iV = 0;
  std::vector<int> correspondanceSegment2V(segind.size(), -1);

  AdaptiveBitModel ctxTempV2[120];

  CtxModelDynamicOBUF ctxTriSoup;
  CtxMapDynamicOBUF MapOBUFTriSoup[3];
  MapOBUFTriSoup[0].reset(14+1, 7); // flag
  MapOBUFTriSoup[1].reset(10+1, 6);  // first bit position
  MapOBUFTriSoup[2].reset(10+1+3+1, 6 + 1);  // second bit position

  const uint8_t initValue0[128] = {
     15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  42,  96,  71,  37,  15,
     15,  22,  51,  15,  15,  30,  27,  15,  15,  64,  15,  48,  15, 224, 171,
    127,  24, 127,  34,  80,  46, 141,  44,  66,  49, 127, 116, 140, 116, 105,
     39, 127, 116, 114,  46, 172, 109,  60,  73, 181, 161, 112,  65, 240, 159,
    127, 127, 127,  87, 183, 127, 116, 116, 195,  88, 152, 141, 228, 141, 127,
     80, 127, 127, 160,  92, 224, 167, 129, 135, 240, 183, 240, 184, 240, 240,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127
  };
  const uint8_t initValue1[64] = {
    116, 127, 118,  15, 104,  56,  97,  15,  96,  15,  29,  15,  95,  15,  46,
     15, 196, 116, 182,  53, 210, 104, 163,  69, 169,  15, 114,  15, 121,  15,
    167,  63, 240, 127, 184,  92, 240, 163, 197,  77, 239,  73, 179,  59, 213,
     48, 185, 108, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127
  };
  const uint8_t initValue2[128] = {
    141, 127, 127, 127, 189,  81,  36, 127, 143, 105, 103, 116, 201,  60,  38,
    116, 116, 127,  15, 127, 153,  59,  15, 116,  69, 105,  15, 127, 158,  93,
     36,  79, 141, 161, 116, 127, 197, 102,  53, 127, 177, 125,  88,  79, 209,
     75, 102,  28,  95,  74,  72,  56, 189,  62,  78,  18,  88, 116,  28,  45,
    237, 100, 152,  35, 141, 240, 127, 127, 208, 133, 101, 141, 186, 210, 168,
     98, 201, 124, 138,  15, 195, 194, 103,  94, 229,  82, 167,  23,  92, 197,
    112,  59, 185,  87, 156,  79, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127
  };
  MapOBUFTriSoup[0].init(initValue0);
  MapOBUFTriSoup[1].init(initValue1);
  MapOBUFTriSoup[2].init(initValue2);

  uint8_t _BufferOBUFleaves[CtxMapDynamicOBUF::kLeafBufferSize * (1 << CtxMapDynamicOBUF::kLeafDepth)];
  memset(_BufferOBUFleaves, 0, sizeof(uint8_t) * CtxMapDynamicOBUF::kLeafBufferSize * (1 << CtxMapDynamicOBUF::kLeafDepth));
  int _OBUFleafNumber = 0;

  for (int i = 0; i <= gbh.num_unique_segments_minus1; i++) {
    // reduced neighbour contexts
    int ctxE = (!!(neighbNodes[i] & 1)) + (!!(neighbNodes[i] & 2)) + (!!(neighbNodes[i] & 4)) + (!!(neighbNodes[i] & 8)) - 1; // at least one node is occupied
    int ctx0 = (!!(neighbNodes[i] & 16)) + (!!(neighbNodes[i] & 32)) + (!!(neighbNodes[i] & 64)) + (!!(neighbNodes[i] & 128));
    int ctx1 = (!!(neighbNodes[i] & 256)) + (!!(neighbNodes[i] & 512)) + (!!(neighbNodes[i] & 1024)) + (!!(neighbNodes[i] & 2048));
    int direction = neighbNodes[i] >> 13; // 0=x, 1=y, 2=z

    // construct pattern
    auto patternIdx = edgePattern[i];
    int pattern = 0;
    int patternClose  = 0;
    int patternClosest  = 0;
    int nclosestPattern = 0;

    int towardOrAway[18] = { // 0 = toward; 1 = away
      0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    int mapping18to9[3][9] = {
      { 0, 1, 2, 3,  4, 15, 14, 5,  7},
      { 0, 1, 2, 3,  9, 15, 14, 7, 12},
      { 0, 1, 2, 9, 10, 15, 14, 7, 12}
    };

    for (int v = 0; v < 9; v++) {
      int v18 = mapping18to9[direction][v];

      if (patternIdx[v18] != -1) {
        int idxEdge = patternIdx[v18];
        if (segind[idxEdge]) {
          pattern |= 1 << v;
          int vertexPos2bits = vertices[correspondanceSegment2V[idxEdge]] >> std::max(0, nbitsVertices - 2);
          if (towardOrAway[v18])
            vertexPos2bits = max2bits - vertexPos2bits; // reverses for away
          if (vertexPos2bits >= mid2bits)
            patternClose |= 1 << v;
          if (vertexPos2bits >= max2bits)
            patternClosest |= 1 << v;
          nclosestPattern += vertexPos2bits >= max2bits && v <= 4;
        }
      }
    }

    int missedCloseStart = /*!(pattern & 1)*/ + !(pattern & 2) + !(pattern & 4);
    int nclosestStart = !!(patternClosest & 1) + !!(patternClosest & 2) + !!(patternClosest & 4);
    if (direction == 0) {
      missedCloseStart +=  !(pattern & 8) + !(pattern & 16);
      nclosestStart +=  !!(patternClosest & 8) + !!(patternClosest & 16);
    }
    if (direction == 1) {
      missedCloseStart +=  !(pattern & 8);
      nclosestStart +=  !!(patternClosest & 8) - !!(patternClosest & 16) ;
    }
    if (direction == 2) {
      nclosestStart +=  - !!(patternClosest & 8) - !!(patternClosest & 16) ;
    }

    // reorganize neighbours of vertex /edge (endpoint) independently on xyz
    int neighbEdge = (neighbNodes[i] >> 0) & 15;
    int neighbEnd = (neighbNodes[i] >> 4) & 15;
    int neighbStart = (neighbNodes[i] >> 8) & 15;
    if (direction == 2) {
      neighbEdge = ((neighbNodes[i] >> 0 + 0) & 1);
      neighbEdge += ((neighbNodes[i] >> 0 + 3) & 1) << 1;
      neighbEdge += ((neighbNodes[i] >> 0 + 1) & 1) << 2;
      neighbEdge += ((neighbNodes[i] >> 0 + 2) & 1) << 3;

      neighbEnd = ((neighbNodes[i] >> 4 + 0) & 1);
      neighbEnd += ((neighbNodes[i] >> 4 + 3) & 1) << 1;
      neighbEnd += ((neighbNodes[i] >> 4 + 1) & 1) << 2;
      neighbEnd += ((neighbNodes[i] >> 4 + 2) & 1) << 3;

      neighbStart = ((neighbNodes[i] >> 8 + 0) & 1);
      neighbStart += ((neighbNodes[i] >> 8 + 3) & 1) << 1;
      neighbStart += ((neighbNodes[i] >> 8 + 1) & 1) << 2;
      neighbStart += ((neighbNodes[i] >> 8 + 2) & 1) << 3;
    }


    // encode flag vertex
    int ctxMap1 = std::min(nclosestPattern, 2) * 15 * 2 +  (neighbEdge-1) * 2 + ((ctx1 == 4));    // 2* 15 *3 = 90 -> 7 bits

    int ctxMap2 = neighbEnd << 11;
    ctxMap2 |= (patternClose & (0b00000110)) << 9 - 1 ; // perp that do not depend on direction = to start
    ctxMap2 |= direction << 7;
    ctxMap2 |= (patternClose & (0b00011000))<< 5-3; // perp that  depend on direction = to start or to end
    ctxMap2 |= (patternClose & (0b00000001))<< 4;  // before
    int orderedPclosePar = (((pattern >> 5) & 3) << 2) + (!!(pattern & 128) << 1) + !!(pattern & 256);
    ctxMap2 |= orderedPclosePar;
    auto index0 = MapOBUFTriSoup[0].getEvolve(
      segind[i], ctxMap2, ctxMap1, &_OBUFleafNumber, _BufferOBUFleaves);
    arithmeticEncoder->encode(
      (int)segind[i], index0 >> 3, ctxTriSoup[index0],
      ctxTriSoup.obufSingleBound);

    // encode position vertex
    if (segind[i]) {
      int v = 0;
      auto vertex = vertices[iV];
      correspondanceSegment2V[i] = iV;

      int ctxFullNbounds = (4 * (ctx0 <= 1 ? 0 : (ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctx1) - 1)) * 2 + (ctxE == 3);
      int b = nbitsVertices - 1;

      // first bit
      ctxMap1 = ctxFullNbounds * 2 + (nclosestStart > 0);
      ctxMap2 = missedCloseStart << 8;
      ctxMap2 |= (patternClosest & 1) << 7;
      ctxMap2 |= direction << 5;
      ctxMap2 |= patternClose & (0b00011111);
      int orderedPclosePar = (((patternClose >> 5) & 3) << 2) + (!!(patternClose & 128) << 1) + !!(patternClose & 256);

      int bit = (vertex >> b--) & 1;

      auto index1 = MapOBUFTriSoup[1].getEvolve(
        bit, ctxMap2, ctxMap1, &_OBUFleafNumber, _BufferOBUFleaves);
      arithmeticEncoder->encode(
        bit, index1 >> 3, ctxTriSoup[index1], ctxTriSoup.obufSingleBound);
      v = bit;

      // second bit
      if (b >= 0) {
        ctxMap1 = ctxFullNbounds * 2 + (nclosestStart > 0);
        ctxMap2 = missedCloseStart << 8;
        ctxMap2 |= (patternClose & 1) << 7;
        ctxMap2 |= (patternClosest & 1) << 6;
        ctxMap2 |= direction << 4;
        ctxMap2 |= (patternClose & (0b00011111)) >> 1;
        ctxMap2 = (ctxMap2 << 4) + orderedPclosePar;

        bit = (vertex >> b--) & 1;
        auto index2 = MapOBUFTriSoup[2].getEvolve(
          bit, ctxMap2, (ctxMap1 << 1) + v, &_OBUFleafNumber,
          _BufferOBUFleaves);
        arithmeticEncoder->encode(
          bit, index2 >> 3, ctxTriSoup[index2], ctxTriSoup.obufSingleBound);
        v = (v << 1) | bit;
      }

      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 = (5 * (ctx0 >> 1) + missedCloseStart) * 2 + (ctxE == 3);
        bit = (vertex >> b--) & 1;
        arithmeticEncoder->encode(
          bit, ctxTempV2[4 * ctxFullNboundsReduced1 + v]);
        v = (v << 1) | bit;
      }

      // remaining bits are bypassed
      for (; b >= 0; b--)
        arithmeticEncoder->encode((vertex >> b) & 1);
      iV++;
    }
  }

}

//-------------------------------------------------------------------------------------
void
encodeTrisoupCentroidResidue(
  std::vector<CentroidDrift>& drifts, pcc::EntropyEncoder* arithmeticEncoder)
{
  AdaptiveBitModel ctxDrift0[9];
  AdaptiveBitModel ctxDriftSign[3][8][8];
  AdaptiveBitModel ctxDriftMag[4];
  for (int i = 0; i < drifts.size(); i++) {
    int driftQ = drifts[i].driftQ;
    arithmeticEncoder->encode(driftQ == 0, ctxDrift0[drifts[i].ctxMinMax]);

    // if not 0
    // drift in [-lowBound; highBound]
    if (driftQ) {
      int lowBound = drifts[i].lowBound;
      int highBound = drifts[i].highBound;
      // code sign
      int lowS = std::min(7, drifts[i].lowBoundSurface);
      int highS = std::min(7, drifts[i].highBoundSurface);
      if (highBound && lowBound) {  // otherwise sign is known
        arithmeticEncoder->encode(
          driftQ > 0,
          ctxDriftSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)]
                      [lowS][highS]);
      }

      // code remaining bits 1 to 7 at most
      int magBound = (driftQ > 0 ? highBound : lowBound) - 1;

      int magDrift = std::abs(driftQ) - 1;
      int ctx = 0;
      while (magBound > 0 && magDrift >= 0) {
        if (ctx < 4)
          arithmeticEncoder->encode(magDrift == 0, ctxDriftMag[ctx]);
        else
          arithmeticEncoder->encode(magDrift == 0);
        magDrift--;
        magBound--;
        ctx++;
      }
    }  // end if not 0

  }  // end loop on drifts
}

//============================================================================

}  // namespace pcc
