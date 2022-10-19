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
  const TrisoupEncOpts& opt,
  const OctreeEncOpts& optOctree,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  PCCPointSet3& predPointCloud,
  const SequenceParameterSet& sps,
  const InterGeomEncOpts& interParams)
{
  // trisoup uses octree coding until reaching the triangulation level.
  pcc::ringbuf<PCCOctree3Node> nodes;
  encodeGeometryOctree(
    optOctree, gps, gbh, pointCloud, ctxtMemOctree, arithmeticEncoders, &nodes,
    predPointCloud, sps, interParams);

  // resume encoding with the last encoder
  pcc::EntropyEncoder* arithmeticEncoder = arithmeticEncoders.back().get();

  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);
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
    nodes, segind, vertices, pointCloud, blockWidth, bitDropped,
    distanceSearchEncoder);

  // Determine neighbours
  std::vector<uint16_t> neighbNodes;
  std::vector<int> indexBefore;
  std::vector<std::vector<int>> perpVertexStart;
  determineTrisoupNeighbours(
    nodes, neighbNodes, indexBefore, perpVertexStart, blockWidth);

  gbh.num_unique_segments_minus1 = segind.size() - 1;
  gbh.num_unique_segments_bits_minus1 =
    numBits(gbh.num_unique_segments_minus1) - 1;

  // Encode vertex presence and position into bitstream
  assert(segind.size() > 0);
  encodeTrisoupVertices(
    segind, vertices, neighbNodes, indexBefore, perpVertexStart, bitDropped,
    gps, gbh, arithmeticEncoder);

  // Decode vertices with certain sampling value
  bool haloFlag = gbh.trisoup_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;

  std::vector<CentroidDrift> drifts;
  int subsample = 1;
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;
  std::cout << "Loop on sampling for max "
            << (gbh.footer.geom_num_points_minus1 + 1) << " points \n";
  if (gps.trisoup_sampling_value > 0) {
    subsample = gps.trisoup_sampling_value;
    decodeTrisoupCommon(
      nodes, segind, vertices, drifts, pointCloud, blockWidth, maxval,
      subsample, bitDropped, isCentroidDriftActivated, false, haloFlag,
      fineRayFlag, NULL);
    std::cout << "Sub-sampling " << subsample << " gives "
              << pointCloud.getPointCount() << " points \n";
  } else {
    int maxSubsample = 1 << gbh.trisoupNodeSizeLog2(gps);
    for (subsample = 1; subsample <= maxSubsample; subsample++) {
      decodeTrisoupCommon(
        nodes, segind, vertices, drifts, pointCloud, blockWidth, maxval,
        subsample, bitDropped, isCentroidDriftActivated, false, haloFlag,
        fineRayFlag, NULL);

      std::cout << "Sub-sampling " << subsample << " gives "
                << pointCloud.getPointCount() << " points \n";
      if (pointCloud.getPointCount() <= gbh.footer.geom_num_points_minus1 + 1)
        break;
    }
  }

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
  const int defaultBlockWidth,
  const int bitDropped,
  int distanceSearchEncoder)
{
  // Put all leaves' edges into a list.
  std::vector<TrisoupSegmentEnc> segments;
  segments.reserve(12 * leaves.size());
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
      {leaf.pos + pos000,
       leaf.pos + posW00,
       12 * i + 0,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc seg0000W0 =  // far left edge
      {leaf.pos + pos000,
       leaf.pos + pos0W0,
       12 * i + 1,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction y
    TrisoupSegmentEnc seg0W0WW0 =  // far top edge
      {leaf.pos + pos0W0,
       leaf.pos + posWW0,
       12 * i + 2,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc segW00WW0 =  // far right edge
      {leaf.pos + posW00,
       leaf.pos + posWW0,
       12 * i + 3,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction y
    TrisoupSegmentEnc seg00000W =  // bottom left edge
      {leaf.pos + pos000,
       leaf.pos + pos00W,
       12 * i + 4,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc seg0W00WW =  // top left edge
      {leaf.pos + pos0W0,
       leaf.pos + pos0WW,
       12 * i + 5,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc segWW0WWW =  // top right edge
      {leaf.pos + posWW0,
       leaf.pos + posWWW,
       12 * i + 6,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc segW00W0W =  // bottom right edge
      {leaf.pos + posW00,
       leaf.pos + posW0W,
       12 * i + 7,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction z
    TrisoupSegmentEnc seg00WW0W =  // near bottom edge
      {leaf.pos + pos00W,
       leaf.pos + posW0W,
       12 * i + 8,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc seg00W0WW =  // near left edge
      {leaf.pos + pos00W,
       leaf.pos + pos0WW,
       12 * i + 9,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction y
    TrisoupSegmentEnc seg0WWWWW =  // near top edge
      {leaf.pos + pos0WW,
       leaf.pos + posWWW,
       12 * i + 10,
       -1,
       -1,
       0,
       0,
       0,
       0};                         // direction x
    TrisoupSegmentEnc segW0WWWW =  // near right edge
      {leaf.pos + posW0W,
       leaf.pos + posWWW,
       12 * i + 11,
       -1,
       -1,
       0,
       0,
       0,
       0};  // // direction y

    // Each voxel votes for a position along each edge it is close to
    const int tmin = 1;
    const int tmax = blockWidth - tmin - 1;
    const int tmin2 = distanceSearchEncoder;
    const int tmax2 = blockWidth - tmin2 - 1;
    for (int j = leaf.start; j < leaf.end; j++) {
      Vec3<int> voxel = pointCloud[j] - leaf.pos;

      // parameter indicating threshold of how close voxels must be to edge ----------- 1 -------------------
      // to be relevant
      if (voxel[1] < tmin && voxel[2] < tmin) {
        seg000W00.count++;
        seg000W00.distanceSum += voxel[0];
      }  // far bottom edge
      if (voxel[0] < tmin && voxel[2] < tmin) {
        seg0000W0.count++;
        seg0000W0.distanceSum += voxel[1];
      }  // far left edge
      if (voxel[1] > tmax && voxel[2] < tmin) {
        seg0W0WW0.count++;
        seg0W0WW0.distanceSum += voxel[0];
      }  // far top edge
      if (voxel[0] > tmax && voxel[2] < tmin) {
        segW00WW0.count++;
        segW00WW0.distanceSum += voxel[1];
      }  // far right edge
      if (voxel[0] < tmin && voxel[1] < tmin) {
        seg00000W.count++;
        seg00000W.distanceSum += voxel[2];
      }  // bottom left edge
      if (voxel[0] < tmin && voxel[1] > tmax) {
        seg0W00WW.count++;
        seg0W00WW.distanceSum += voxel[2];
      }  // top left edge
      if (voxel[0] > tmax && voxel[1] > tmax) {
        segWW0WWW.count++;
        segWW0WWW.distanceSum += voxel[2];
      }  // top right edge
      if (voxel[0] > tmax && voxel[1] < tmin) {
        segW00W0W.count++;
        segW00W0W.distanceSum += voxel[2];
      }  // bottom right edge
      if (voxel[1] < tmin && voxel[2] > tmax) {
        seg00WW0W.count++;
        seg00WW0W.distanceSum += voxel[0];
      }  // near bottom edge
      if (voxel[0] < tmin && voxel[2] > tmax) {
        seg00W0WW.count++;
        seg00W0WW.distanceSum += voxel[1];
      }  // near left edge
      if (voxel[1] > tmax && voxel[2] > tmax) {
        seg0WWWWW.count++;
        seg0WWWWW.distanceSum += voxel[0];
      }  // near top edge
      if (voxel[0] > tmax && voxel[2] > tmax) {
        segW0WWWW.count++;
        segW0WWWW.distanceSum += voxel[1];
      }  // near right edge

      // parameter indicating threshold of how close voxels must be to edge ----------- 2 -------------------
      // to be relevant
      if (voxel[1] < tmin2 && voxel[2] < tmin2) {
        seg000W00.count2++;
        seg000W00.distanceSum2 += voxel[0];
      }  // far bottom edge
      if (voxel[0] < tmin2 && voxel[2] < tmin2) {
        seg0000W0.count2++;
        seg0000W0.distanceSum2 += voxel[1];
      }  // far left edge
      if (voxel[1] > tmax2 && voxel[2] < tmin2) {
        seg0W0WW0.count2++;
        seg0W0WW0.distanceSum2 += voxel[0];
      }  // far top edge
      if (voxel[0] > tmax2 && voxel[2] < tmin2) {
        segW00WW0.count2++;
        segW00WW0.distanceSum2 += voxel[1];
      }  // far right edge
      if (voxel[0] < tmin2 && voxel[1] < tmin2) {
        seg00000W.count2++;
        seg00000W.distanceSum2 += voxel[2];
      }  // bottom left edge
      if (voxel[0] < tmin2 && voxel[1] > tmax2) {
        seg0W00WW.count2++;
        seg0W00WW.distanceSum2 += voxel[2];
      }  // top left edge
      if (voxel[0] > tmax2 && voxel[1] > tmax2) {
        segWW0WWW.count2++;
        segWW0WWW.distanceSum2 += voxel[2];
      }  // top right edge
      if (voxel[0] > tmax2 && voxel[1] < tmin2) {
        segW00W0W.count2++;
        segW00W0W.distanceSum2 += voxel[2];
      }  // bottom right edge
      if (voxel[1] < tmin2 && voxel[2] > tmax2) {
        seg00WW0W.count2++;
        seg00WW0W.distanceSum2 += voxel[0];
      }  // near bottom edge
      if (voxel[0] < tmin2 && voxel[2] > tmax2) {
        seg00W0WW.count2++;
        seg00W0WW.distanceSum2 += voxel[1];
      }  // near left edge
      if (voxel[1] > tmax2 && voxel[2] > tmax2) {
        seg0WWWWW.count2++;
        seg0WWWWW.distanceSum2 += voxel[0];
      }  // near top edge
      if (voxel[0] > tmax2 && voxel[2] > tmax2) {
        segW0WWWW.count2++;
        segW0WWWW.distanceSum2 += voxel[1];
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
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  int bitDropped,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder)
{
  const int nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;
  const int max2bits = nbitsVertices > 1 ? 3 : 1;

  int iV = 0;
  std::vector<int> correspondanceSegment2V(segind.size(), -1);

  AdaptiveBitModel ctxTempV2[144];

  CtxModelTrisoup ctxTriSoup;
  CtxMapTriSoup MapOBUFTriSoup[3];
  MapOBUFTriSoup[0].reset(10, 7);      // flag
  MapOBUFTriSoup[1].reset(10, 6);      // first bit position
  MapOBUFTriSoup[2].reset(10, 6 + 1);  // second bit position
  for (int i = 0; i <= gbh.num_unique_segments_minus1; i++) {
    // reduced neighbour contexts
    int ctxE = (!!(neighbNodes[i] & 1)) + (!!(neighbNodes[i] & 2))
      + (!!(neighbNodes[i] & 4)) + (!!(neighbNodes[i] & 8))
      - 1;  // at least one node is occupied
    int ctx0 = (!!(neighbNodes[i] & 16)) + (!!(neighbNodes[i] & 32))
      + (!!(neighbNodes[i] & 64)) + (!!(neighbNodes[i] & 128));
    int ctx1 = (!!(neighbNodes[i] & 256)) + (!!(neighbNodes[i] & 512))
      + (!!(neighbNodes[i] & 1024)) + (!!(neighbNodes[i] & 2048));
    int direction = neighbNodes[i] >> 13;

    int beforeCtx = 0;
    int Vbefore = 0;
    int nclose = 0;
    int nfar = 0;

    // compute data relative to preceding (before) vertex along edge direction (presence and position)
    // vertoces are on gbh.trisoupNodeSizeLog2(gps) - bitDropped bits
    if (indexBefore[i] != -1 && gbh.trisoupNodeSizeLog2(gps) <= 4) {
      beforeCtx = segind[indexBefore[i]];
      if (correspondanceSegment2V[indexBefore[i]] != -1) {
        Vbefore = 1 + vertices[correspondanceSegment2V[indexBefore[i]]]
          >> std::max(0, nbitsVertices - 2);  // on 2 bits

        int v2bits = max2bits
          - (vertices[correspondanceSegment2V[indexBefore[i]]]
             >> std::max(0, nbitsVertices - 2));  // on 2 bits
        if (v2bits <= 0)
          nclose++;
        if (v2bits >= max2bits)
          nfar++;
      }
    }

    // count number of perp vertices close to current edge; find their position
    int occupPerp = 0;
    auto perpVS = perpVertexStart[i];
    int maxVal = (1 << nbitsVertices) - 1;

    for (int k = 0; k < perpVS.size(); k++) {
      const int maskIdx = (1 << 30) - 1;
      int idxEdge = perpVS[k] & maskIdx;

      if (segind[idxEdge]) {
        occupPerp++;
        int idxVertex = correspondanceSegment2V[idxEdge];
        int vertexPos = vertices[idxVertex] >> std::max(0, nbitsVertices - 2);  // on 2 bits

        int orientation = perpVS[k] >> 30;
        if (orientation) {            // if toward then reverse to away
          vertexPos = max2bits - vertexPos;  // 0 is closest, 3 is farthest
        }

        if (vertexPos <= 0)
          nclose++;
        if (vertexPos >= max2bits)
          nfar++;
      }
    }

    int perpStartCtx =
      (4 - direction) - occupPerp + 1 - beforeCtx;  // in [0,5]
    int nmiddle = occupPerp + beforeCtx - nfar - nclose;
    bool flagTouch = nclose > 0;

    // encode flag vertex
    int ctxMap1 = (Vbefore * 4 + ctxE) * 4 + std::min(nclose, 3);
    bool full01 = ((ctx0 == 4) || (ctx1 == 4));
    int ctxMap2 = full01 << 8;
    if (!full01) {  // none is full
      ctxMap2 |= (neighbNodes[i] >> 4) & 255;
    } else {  // one is full, very rarely both
      ctxMap2 |= (ctx0 == 4) << 7;
      if (ctx0 == 4)
        ctxMap2 |= ((neighbNodes[i] >> 8) & 15) << 3;
      else
        ctxMap2 |= ((neighbNodes[i] >> 4) & 15) << 3;
      // three bits to spare
      ctxMap2 |= perpStartCtx;
    }

    arithmeticEncoder->encode(
      (int)segind[i],
      ctxTriSoup[MapOBUFTriSoup[0].getEvolve(segind[i], ctxMap2, ctxMap1)]);

    // encode position vertex
    if (segind[i]) {
      int v = 0;
      auto vertex = vertices[iV];
      correspondanceSegment2V[i] = iV;

      int ctxFullNbounds =
        (4 * (std::max(1, ctx0) - 1) + (std::max(1, ctx1) - 1)) * 2
        + (ctxE == 3);
      int b = nbitsVertices - 1;

      // first bit
      ctxMap1 = ctxFullNbounds * 2 + flagTouch;
      ctxMap2 = perpStartCtx << 7;
      ctxMap2 |= beforeCtx << 6;
      ctxMap2 |= std::min(3, nfar) << 4;
      ctxMap2 |= std::min(3, nclose) << 2;
      ctxMap2 |= std::min(3, nmiddle);
      int bit = (vertex >> b--) & 1;
      arithmeticEncoder->encode(
        bit, ctxTriSoup[MapOBUFTriSoup[1].getEvolve(bit, ctxMap2, ctxMap1)]);
      v = bit;

      // second bit
      if (b >= 0) {
        bit = (vertex >> b--) & 1;
        arithmeticEncoder->encode(
          bit,
          ctxTriSoup[MapOBUFTriSoup[2].getEvolve(
            bit, ctxMap2, (ctxMap1 << 1) + v)]);
        v = (v << 1) | bit;
      }

      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 =
          (6 * (ctx0 >> 1) + perpStartCtx) * 2 + (ctxE == 3);
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
