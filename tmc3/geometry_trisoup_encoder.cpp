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
  const bool isFaceVertexActivated = gbh.trisoup_face_vertex_flag;

  // Determine vertices
  std::cout << "Number of points = " << pointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes = " << nodes.size() << "\n";
  float estimatedSampling = 1;
  int distanceSearchEncoder = 1;
  if (opt.improvedVertexDetermination) {
    estimatedSampling = float(nodes.size());
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

  std::vector<node6nei> nodes6nei;
  if (isFaceVertexActivated) {
    determineTrisoupNodeNeighbours(nodes, nodes6nei, blockWidth);
  }

  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  std::vector<TrisoupNodeEdgeVertex> eVerts;
  determineTrisoupVertices(
    nodes, segind, vertices, pointCloud,
    gps, gbh,
    blockWidth, bitDropped, eVerts,
    distanceSearchEncoder, nodesPadded, pointCloudPadding, indices,
    originalBox, estimatedSampling, opt.nodeUniqueDSE);

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


  std::vector<TrisoupCentroidVertex> cVerts;
  std::vector<CentroidDrift> drifts;
  std::vector<Vec3<int32_t>> normVs;
  std::vector<Vec3<int32_t>> gravityCenter;
  determineTrisoupCentroids(
    pointCloud, nodes, gps, gbh, blockWidth, bitDropped,
    isCentroidDriftActivated, eVerts, gravityCenter, drifts, cVerts, normVs);

  std::vector<TrisoupFace> faces;
  std::vector<TrisoupFace> limited_faces;
  std::vector<TrisoupNodeFaceVertex> fVerts;
  fVerts.resize(nodes.size());
  if (isFaceVertexActivated) {
    determineTrisoupFaceVertices(
      pointCloud, nodes, gps, gbh, nodes6nei, blockWidth,
      distanceSearchEncoder, eVerts, gravityCenter, cVerts, fVerts, normVs,
      limited_faces, faces);
  }

  // Decode vertices with certain sampling value
  bool haloFlag = gbh.trisoup_halo_flag;
  bool adaptiveHaloFlag = gbh.trisoup_adaptive_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;

  PCCPointSet3 recPointCloud;
  recPointCloud.addRemoveAttributes(pointCloud);

  int subsample = 1;
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;
  std::cout << "Loop on sampling for max "
            << (gbh.footer.geom_num_points_minus1 + 1) << " points \n";
  if (gps.trisoup_sampling_value > 0) {
    subsample = gps.trisoup_sampling_value;
    decodeTrisoupCommon(
      nodes, nodes6nei, eVerts, cVerts, gravityCenter, normVs, faces, fVerts,
      recPointCloud, gps, gbh, blockWidth, maxval, subsample, bitDropped,
      isCentroidDriftActivated, isFaceVertexActivated, haloFlag,
      adaptiveHaloFlag, fineRayFlag);
    std::cout << "Sub-sampling " << subsample << " gives "
              << recPointCloud.getPointCount() << " points \n";
  } else {
    int maxSubsample = 1 << gbh.trisoupNodeSizeLog2(gps);
    for (subsample = 1; subsample <= maxSubsample; subsample++) {
      decodeTrisoupCommon(
        nodes, nodes6nei, eVerts, cVerts, gravityCenter, normVs, faces,
        fVerts, recPointCloud, gps, gbh, blockWidth, maxval, subsample,
        bitDropped, isCentroidDriftActivated, isFaceVertexActivated, haloFlag,
        adaptiveHaloFlag, fineRayFlag);
      std::cout << "Sub-sampling " << subsample << " gives "
                << recPointCloud.getPointCount() << " points \n";
      if (recPointCloud.getPointCount() <= gbh.footer.geom_num_points_minus1 + 1)
        break;
    }
  }

  pointCloud.resize(0);
  pointCloud = std::move(recPointCloud);

  gbh.trisoup_sampling_value_minus1 = subsample - 1;


  // encode centroid residual and node-faces into bitstream
  if (isCentroidDriftActivated){
    encodeTrisoupCentroidResidue(drifts, arithmeticEncoder);
    if (isFaceVertexActivated) {
      pcc::EntropyEncoder* _aec = arithmeticEncoders.back().get();
      encodeTrisoupFaceList(limited_faces, gps, gbh, arithmeticEncoder);
    }
  }
}

//----------------------------------------------------------------------------

float
estimatedSampling1(PCCOctree3Node leaf, Vec3<int32_t> newW)
{
  return
    sqrt(float(newW.max() * newW.mid())) / sqrt(float(leaf.end - leaf.start));
}

//----------------------------------------------------------------------------

float
estimatedSampling2(PCCOctree3Node leaf, const PCCPointSet3& pointCloud)
{
  Vec3<int> min = pointCloud[leaf.start] - leaf.pos;
  Vec3<int> max = pointCloud[leaf.start] - leaf.pos;
  for (int j = leaf.start; j < leaf.end; j++) {
    Vec3<int> currentVoxel = pointCloud[j] - leaf.pos;
    min[0] = currentVoxel[0] < min[0] ? currentVoxel[0] : min[0];
    min[1] = currentVoxel[1] < min[1] ? currentVoxel[1] : min[1];
    min[2] = currentVoxel[2] < min[2] ? currentVoxel[2] : min[2];
    max[0] = currentVoxel[0] > max[0] ? currentVoxel[0] : max[0];
    max[1] = currentVoxel[1] > max[1] ? currentVoxel[1] : max[1];
    max[2] = currentVoxel[2] > max[2] ? currentVoxel[2] : max[2];
  }
  Vec3<int> dim = max - min;
  return sqrt(float(dim.max() * dim.mid()) / float(leaf.end - leaf.start));
}

//----------------------------------------------------------------------------

float
estimatedSampling3(PCCOctree3Node leaf, const PCCPointSet3& pointCloud)
{
  std::vector<std::vector<float>> vec_nn(leaf.end - leaf.start);
  std::vector<int> vec_one(leaf.end - leaf.start, 0);
  int N = 4;
  int cnt1 = 0;
  float es = 0;
  for (int j = leaf.start; j < leaf.end; j++) {
    Vec3<int> currentVoxel = pointCloud[j] - leaf.pos;
    int cnt2 = cnt1 + 1;
    for (int i = leaf.start + cnt2; i < leaf.end; i++) {
      float distance = sqrt(
        float(((currentVoxel - (pointCloud[i] - leaf.pos))).getNorm2()));
      if (vec_nn[cnt1].size() < N) {
        vec_nn[cnt1].push_back(distance);
        std::sort(vec_nn[cnt1].begin(), vec_nn[cnt1].end());
      } else if (distance < vec_nn[cnt1].back() && vec_one[cnt1] < N) {
        vec_nn[cnt1][N-1] = distance;
        std::sort(vec_nn[cnt1].begin(), vec_nn[cnt1].end());
      }
      if (vec_nn[cnt2].size() < N) {
        vec_nn[cnt2].push_back(distance);
        std::sort(vec_nn[cnt2].begin(), vec_nn[cnt2].end());
      } else if (distance < vec_nn[cnt2].back() && vec_one[cnt2] < N) {
        vec_nn[cnt2][N-1] = distance;
        std::sort(vec_nn[cnt2].begin(), vec_nn[cnt2].end());
      }
      if (distance <= 1.0f) {
        ++ vec_one[cnt1];
        ++ vec_one[cnt2];
      }
      ++ cnt2;
    }
    float nn = 0;
    int n = vec_nn[cnt1].size();
    for (int k = 0; k < n; k++) {
      nn += vec_nn[cnt1][k];
    }
    es += (nn / float(n));
    ++ cnt1;
  }
  return es / float(leaf.end - leaf.start);
}

//----------------------------------------------------------------------------
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
  std::vector<TrisoupNodeEdgeVertex>& eVerts,
  int distanceSearchEncoder,
  const std::vector<PCCOctree3Node>& nodesPadded,
  const PCCPointSet3& pointCloudPadding,
  std::vector<int> indices,
  Box3<int32_t> originalBox,
  float estimatedSampling,
  bool nodeUniqueDSE)
{
  // not use
  std::vector<uint16_t> neighbNodes;
  std::vector<std::array<int, 18>> edgePattern;
  pcc::EntropyDecoder arithmeticDecoder;

  processTrisoupVertices(
    gps, gbh, leaves, defaultBlockWidth, bitDropped, false, pointCloud,
    distanceSearchEncoder, neighbNodes, edgePattern, arithmeticDecoder,
    eVerts, segind, vertices, nodesPadded, pointCloudPadding, indices,
    originalBox, estimatedSampling, nodeUniqueDSE);
}

// ---------------------------------------------------------------------------

void processTrisoupVertices(
  // common input
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const ringbuf<PCCOctree3Node>& leaves,
  const int defaultBlockWidth,
  const int bitDropped,

  bool isDecoder,

  // input only encoder
  const PCCPointSet3& pointCloud,
  const int distanceSearchEncoder,

  // input only decoder
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  pcc::EntropyDecoder& arithmeticDecoder,

  // common output
  std::vector<TrisoupNodeEdgeVertex>& eVerts,

  // output only encoder
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,

  const std::vector<PCCOctree3Node>& nodesPadded,
  const PCCPointSet3& pointCloudPadding,
  std::vector<int> indices,
  Box3<int32_t> originalBox,
  float estimatedSampling,
  bool nodeUniqueDSE)
{
  if (isDecoder) {
    decodeTrisoupVerticesSub(
      segind, vertices, neighbNodes, edgePattern, bitDropped, gps, gbh,
      arithmeticDecoder);
  }

  Box3<int32_t> sliceBB;
  sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
  sliceBB.max = sliceBB.min + ( gbh.slice_bb_width << gbh.slice_bb_width_log2_scale );

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
    nonCubicNode(
      gps, gbh, leaf.pos, defaultBlockWidth, sliceBB, newP, newW, corner);

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

    if (!isDecoder) {
      // Each voxel votes for a position along each edge it is close to
      const int tmin = 1;
      const Vec3<int> tmax( newW.x() - tmin - 1,
                            newW.y() - tmin - 1,
                            newW.z() - tmin - 1 );

      int localDistanceSearchEncoder = -1;
      if (nodeUniqueDSE) {
        // Desicion tree
        float es = estimatedSampling;
        if(estimatedSampling > 1.0f) {
          es = estimatedSampling1(leaf, newW);
          if (abs(estimatedSampling - es) > 0.5f){
            es = estimatedSampling2(
              leaf, i < leaves.size() ? pointCloud : pointCloudPadding);
            if (abs(estimatedSampling - es) > 0.5f){
              if (leaf.end - leaf.start > 1) {
                es = estimatedSampling3(
                  leaf, i < leaves.size() ? pointCloud : pointCloudPadding);
              } else {es = estimatedSampling;}
              es = std::min(es, estimatedSampling + 1);
            }
          } else {es = estimatedSampling;}
        } else {es = estimatedSampling;}
        es = std::min(es, float(blockWidth/4));
        localDistanceSearchEncoder = (1 << std::max(0, bitDropped - 2)) - 1;
        localDistanceSearchEncoder += int(std::round(es + 0.1f));
        localDistanceSearchEncoder =
          std::max(1, std::min(8, localDistanceSearchEncoder));
        // End of decision tree
      }

      const int tmin2 =
        nodeUniqueDSE ? localDistanceSearchEncoder : distanceSearchEncoder;
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
    }

    // Push segments onto list.
    if (i < leaves.size()) {
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

      if (!isDecoder) {
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
      }
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

  std::vector<TrisoupSegmentEnc> segmentsPerNode;
  std::copy(
    segments.begin(), segments.end(), std::back_inserter(segmentsPerNode));

  // Sort the list and find unique segments.
  std::sort(segments.begin(), segments.end());

  if (!isDecoder) {
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
          int temp =
            ((2 * localSegment.distanceSum + localSegment.distanceSum2)
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

  // Sort the list and find unique segments.
  std::vector<TrisoupSegment> uniqueSegments;
  uniqueSegments.push_back(segments[0]);
  segmentsPerNode[segments[0].index].uniqueIndex = 0;
  for (int i = 1; i < segments.size(); i++) {
    if (uniqueSegments.back().startpos != segments[i].startpos
        || uniqueSegments.back().endpos != segments[i].endpos) {
      // sortedSegment[i] is different from uniqueSegments.back().
      // Start a new uniqueSegment.
      uniqueSegments.push_back(segments[i]);
    }
    if (segments[i].index < leaves.size() * 12) {
      segmentsPerNode[segments[i].index].uniqueIndex =
        uniqueSegments.size() - 1;
    }
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
  for (int i = 0; i < leaves.size()*12; i++) {
    segmentsPerNode[i].vertex =
      uniqueSegments[segmentsPerNode[i].uniqueIndex].vertex;
  }


  eVerts.clear();
  // ----------- loop on leaf nodes ----------------------
  for (int i = 0; i < leaves.size(); i++) {
    TrisoupNodeEdgeVertex neVertex;
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode(
      gps, gbh, leaves[i].pos, defaultBlockWidth, sliceBB, nodepos, nodew,
      corner);

    // Find up to 12 vertices for this leaf.
    for (int j = 0; j < 12; j++) {
      TrisoupSegment& segment = segmentsPerNode[i * 12 + j];
      if (segment.vertex < 0)
        continue;  // skip segments that do not intersect the surface

      // Get distance along edge of vertex.
      // Vertex code is the index of the voxel along the edge of the block
      // of surface intersection./ Put decoded vertex at center of voxel,
      // unless voxel is first or last along the edge, in which case put the
      // decoded vertex at the start or endpoint of the segment.
      Vec3<int32_t> direction = segment.endpos - segment.startpos;
      uint32_t segment_len = direction.max();

      // Get 3D position of point of intersection.
      Vec3<int32_t> point = (segment.startpos - nodepos) << kTrisoupFpBits;
      point -= kTrisoupFpHalf; // the volume is [-0.5; B-0.5]^3

      // points on edges are located at integer values
      int halfDropped = 0; // bitDropped ? 1 << bitDropped - 1 : 0;
      int32_t distance =
        (segment.vertex << (kTrisoupFpBits + bitDropped)) +
        (kTrisoupFpHalf << bitDropped);
      if (direction[0]) {
        point[0] += distance; // in {0,1,...,B-1}
      } else if (direction[1]) {
        point[1] += distance;
      } else {  // direction[2]
        point[2] += distance;
      }
      // Add vertex to list of points.
      neVertex.vertices.push_back({ point, 0, 0 });
    }

    // compute centroid (gravity center)
    int vtxCount = (int)neVertex.vertices.size();
    Vec3<int32_t> gCenter = 0;
    for (int j = 0; j < vtxCount; j++) {
      gCenter += neVertex.vertices[j].pos;
    }
    if (vtxCount) {
      gCenter /= vtxCount;
    }

    // order vertices along a dominant axis only if more than three (otherwise only one triangle, whatever...)
    neVertex.dominantAxis = findDominantAxis(neVertex.vertices, nodew, gCenter);
    eVerts.push_back(neVertex);
  }
}

void determineTrisoupCentroids(
  PCCPointSet3& pointCloud,
  const ringbuf<PCCOctree3Node>& leaves,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int defaultBlockWidth,
  const int bitDropped,
  const bool isCentroidDriftActivated,
  const std::vector<TrisoupNodeEdgeVertex> eVerts,
  std::vector<Vec3<int32_t>>& gravityCenter,
  std::vector<CentroidDrift>& drifts,
  std::vector<TrisoupCentroidVertex>& cVerts,
  std::vector<Vec3<int32_t>>& normVs)
{
  Box3<int32_t> sliceBB;
  sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
  sliceBB.max = sliceBB.min + ( gbh.slice_bb_width << gbh.slice_bb_width_log2_scale );

  // ----------- loop on leaf nodes ----------------------
  for (int i = 0; i < leaves.size(); i++) {
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode(
      gps, gbh, leaves[i].pos, defaultBlockWidth, sliceBB, nodepos, nodew,
      corner);

    // Skip leaves that have fewer than 3 vertices.
    if (eVerts[i].vertices.size() < 3) {
      cVerts.push_back({ false, { 0, 0, 0 }, 0, true });
      normVs.push_back({ 0, 0, 0 });
      gravityCenter.push_back({ 0, 0, 0 });
      continue;
    }

    Vec3<int32_t> gCenter={ 0 }, normalV={ 0 };
    TrisoupCentroidContext cctx={ 0 };

    // judgment for refinement of the centroid along the domiannt axis
    bool driftCondition =
      determineNormVandCentroidContexts(
        nodew, eVerts[i], bitDropped, gCenter, normalV, cctx);

    if(!(driftCondition && isCentroidDriftActivated)) {
      cVerts.push_back({ false, gCenter, 0, true });
      normVs.push_back(normalV);
      gravityCenter.push_back(gCenter);
      continue;
    }

    // Refinement of the centroid along the domiannt axis
    // determine qauntized drift
    Vec3<int32_t> blockCentroid = gCenter;
    int counter = 0;
    int driftQ = 0, drift = 0;
    int bitDropped2 = bitDropped;
    int maxD = std::max(int(3), bitDropped2);

    // determine quantized drift
    for (int p = leaves[i].start; p < leaves[i].end; p++) {
      Vec3<int32_t> point = (pointCloud[p] - nodepos) << kTrisoupFpBits;
      Vec3<int64_t> CP =
        crossProduct(normalV, point - blockCentroid) >> kTrisoupFpBits;
      int64_t dist = isqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);
      dist >>= kTrisoupFpBits;
      if ((dist << 10) <= 1774 * maxD) {
        int32_t w = (1 << 10) + 4 * (1774 * maxD - ((1 << 10) * dist));
        counter += w >> 10;
        drift +=
          (w >> 10) * ((normalV * (point - blockCentroid)) >> kTrisoupFpBits);
      }
    }

    if (counter) { // drift is shift by kTrisoupFpBits
      drift = (drift >> kTrisoupFpBits - 6) / counter; // drift is shift by 6 bits
    }

    int half = 1 << 5 + bitDropped2;
    int DZ = 2*half/3;

    if (abs(drift) >= DZ) {
      driftQ = (abs(drift) - DZ + 2*half+2*half/3) >> 6 + bitDropped2;
      if (drift < 0)
        driftQ = -driftQ;
    }
    // drift in [-lowBound; highBound]
    driftQ = std::min(std::max(driftQ, -cctx.lowBound), cctx.highBound);
    // push quantized drift  to buffeer for encoding
    drifts.push_back({ driftQ, cctx.lowBound, cctx.highBound, cctx.ctxMinMax,
        cctx.lowBoundSurface, cctx.highBoundSurface });

    // dequantize and apply drift
    int driftDQ = 0;
    if (driftQ) {
      driftDQ = std::abs(driftQ) << bitDropped2 + 6;
      int half = 1 << 5 + bitDropped2;
      int DZ = 2*half/3;
      driftDQ += DZ - half;
      if (driftQ < 0)
        driftDQ = -driftDQ;
    }

    blockCentroid += (driftDQ * normalV) >> 6;
    blockCentroid[0] = std::max(-kTrisoupFpHalf, blockCentroid[0]);
    blockCentroid[1] = std::max(-kTrisoupFpHalf, blockCentroid[1]);
    blockCentroid[2] = std::max(-kTrisoupFpHalf, blockCentroid[2]);
    int blockWidth = defaultBlockWidth;
    blockCentroid[0] = std::min(
      ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1,
      blockCentroid[0]);
    blockCentroid[1] = std::min(
      ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1,
      blockCentroid[1]);
    blockCentroid[2] = std::min(
      ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1,
      blockCentroid[2]);

    bool boundaryInside = true;
    if (!nodeBoundaryInsideCheck( nodew<<kTrisoupFpBits, blockCentroid)) {
      boundaryInside = false;
    }
    cVerts.push_back({ true, blockCentroid, driftDQ, boundaryInside });
    normVs.push_back(normalV);
    gravityCenter.push_back(gCenter);

    // end refinement of the centroid
  }
}


//----------------------------------------------------------------------------

void
determineTrisoupFaceVertices(
  PCCPointSet3& pointCloud,
  const ringbuf<PCCOctree3Node>& leaves,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const std::vector<node6nei> nodes6nei,
  int defaultBlockWidth,
  const int distanceSearchEncoder,
  const std::vector<TrisoupNodeEdgeVertex> eVerts,
  const std::vector<Vec3<int32_t>> gravityCenter,
  const std::vector<TrisoupCentroidVertex> cVerts,
  std::vector<TrisoupNodeFaceVertex>& fVerts,
  std::vector<Vec3<int32_t>> normVs, // currently not used
  std::vector<TrisoupFace>& limited_faces,
  std::vector<TrisoupFace>& faces)
{
  int _ones = 0, _zeros = 0;
  limited_faces.clear();
  faces.clear();
  const int32_t tmin1 = 2*4;
  const int32_t tmin2 = distanceSearchEncoder*4;
  Box3<int32_t> sliceBB;
  sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
  sliceBB.max =
    sliceBB.min + (gbh.slice_bb_width << gbh.slice_bb_width_log2_scale);
  int32_t w = defaultBlockWidth;
  Vec3<int32_t> nodePosOfst[6] = {
    { 0, 0, -w }, { 0, 0, w }, { 0, -w, 0 },
    { 0, w, 0 }, { -w, 0, 0 }, { w, 0, 0 } };
  // For 6-neighbour nodes of three which have bigger coordinates than current
  // node, if current node and its neighbour node has refined centroid vertex
  // each other, and when the centroids are connected, if there is no bite on
  // the current surface, then the intersection of centroid connection segment
  // and node boundary face is defined as the temporary face vertex.
  // And if original points are distributed around the temporary face vertex,
  // it is defined as a true face vertex and determine to connect these
  // centroids, and then face-flag becomes true.
  for (int i = 0; i < leaves.size(); i++) {
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode(
      gps, gbh, leaves[i].pos, defaultBlockWidth, sliceBB, nodepos, nodew,
      corner);
    // to generate 3 faces per node, neighbour direction loop must be placed
    // at the most outer loop.
    // z,y,x-axis order : nei=0,1,2  6nei-idx(j)=1,3,5
    for (int j = 1, nei = 0; j < 6; j += 2, nei++) {

      TrisoupFace face( false );

      if (cVerts[i].valid && cVerts[i].boundaryInside) {
        int ii = nodes6nei[i].idx[j];
        // neighbour-node exists on this direction
        if (-1 != ii) {
          // centroid of the neighbour-node exists and inside of the boundary
          if (cVerts[ii].valid && cVerts[ii].boundaryInside) {
            int eIdx[2][2] = { -1 };
            int axis = 2 - nei; // z,y,x-axis order 2,1,0
            Vec3<int32_t> nodeW = nodew << kTrisoupFpBits;
            Vec3<int32_t> zeroW =     0 << kTrisoupFpBits;
            int neVtxBoundaryFace =
              countTrisoupEdgeVerticesOnFace(eVerts[ i], nodeW, axis);
            if (2 == neVtxBoundaryFace || 3 == neVtxBoundaryFace) {

              Vertex fVert[2];
              findTrisoupFaceVertex(
                i, nei, nodes6nei[i], cVerts, nodew, fVert);

              determineTrisoupEdgeBoundaryLine(
                i, eVerts[ i], cVerts[ i], gravityCenter[ i], nodeW, axis,
                fVert[0], eIdx[0]);
              determineTrisoupEdgeBoundaryLine(
                ii, eVerts[ii], cVerts[ii], gravityCenter[ii], zeroW, axis,
                fVert[1], eIdx[1]);

              if ((-1 != eIdx[0][0]) && (-1 != eIdx[0][1])) {
                bool judge = determineTrisoupDirectionOfCentroidsAndFvert(
                  eVerts[i], cVerts, gravityCenter, nodew, i, nei,
                  nodes6nei[i], ii, w, eIdx[0][0], eIdx[0][1], fVert);
                if (judge) {
                  // c0, c1 and face vertex on the same side of the surface
                  int32_t weight1 = 0, weight2 = 0;
                  uint32_t st[2] = { leaves[i].start, leaves[ii].start };
                  uint32_t ed[2] = { leaves[i].end,   leaves[ii].end   };
                  // order - 0:z, 1:y, 2:x
                  Vec3<int32_t> neiOfst[2][3] = {
                    { { 0, 0, 0 },
                      { 0, 0, 0 },
                      { 0, 0, 0 } },
                    { { 0, 0, nodew[2] },
                      { 0, nodew[1], 0 },
                      { nodew[0], 0, 0 } }
                  };
                  // 0:current-node  1:nei-node
                  for (int n = 0; n < 2; n++) {
                    Vec3<uint32_t> ofst[2] = { {0}, nodePosOfst[j] };
                    for (int k = st[n]; k < ed[n]; k++) {
                      Vec3<int32_t> dist =
                        fVert[n].pos
                          - ((pointCloud[k] - nodepos - neiOfst[n][nei])
                            << kTrisoupFpBits);
                      int32_t d =
                        (dist.abs().max() + kTrisoupFpHalf) >> kTrisoupFpBits;
                      if(d < tmin1) { weight1++; }
                      if(d < tmin2) { weight2++; }
                    }
                  }
                  if (weight1 > 0 || weight2 > 1) {
                    face.connect = true;
                    fVerts[ i].formerEdgeVertexIdx.push_back(eIdx[0][0]);
                    fVerts[ i].vertices.push_back(fVert[0]);
                    fVerts[ii].formerEdgeVertexIdx.push_back(eIdx[1][0]);
                    fVerts[ii].vertices.push_back(fVert[1]);
                  }
                  if (face.connect) { _ones++;  }
                  else {              _zeros++; }
                  limited_faces.push_back(face);
                }
              }
            } // if( 2 or 3 == neVtxBoundaryFace )

          }
        }
      }
      faces.push_back( face );
    }
  }
  return;
}

//----------------------------------------------------------------------------

void encodeTrisoupFaceList(
  std::vector<TrisoupFace>& faces,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder)
{
  AdaptiveBitModel ctxFaces;

  for (int i = 0; i < faces.size(); i++) {
    arithmeticEncoder->encode((int)faces[i].connect, ctxFaces);
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
