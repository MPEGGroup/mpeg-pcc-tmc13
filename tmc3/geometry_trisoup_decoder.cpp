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

#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================

// The number of fractional bits used in trisoup triangle voxelisation
const int kTrisoupFpBits = 8;

// The value 1 in fixed-point representation
const int kTrisoupFpOne = 1 << (kTrisoupFpBits);
const int kTrisoupFpHalf = 1 << (kTrisoupFpBits - 1);
const int truncateValue = kTrisoupFpHalf;

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
bool
operator<(const TrisoupSegmentNeighbours& s1, const TrisoupSegmentNeighbours& s2)
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
  GeometryOctreeContexts& ctxtMemOctree,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const Vec3<int> minimum_position)
{
  // trisoup uses octree coding until reaching the triangulation level.
  // todo(df): pass trisoup node size rather than 0?
  pcc::ringbuf<PCCOctree3Node> nodes;

  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMemOctree, arithmeticDecoder, &nodes,
    refFrame, minimum_position);

  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);
  const int maxVertexPrecisionLog2 = gbh.trisoup_vertex_quantization_bits ? gbh.trisoup_vertex_quantization_bits : gbh.trisoupNodeSizeLog2(gps);
  const int bitDropped =  std::max(0, gbh.trisoupNodeSizeLog2(gps) - maxVertexPrecisionLog2);
  const bool isCentroidDriftActivated = gbh.trisoup_centroid_vertex_residual_flag;


  // Determine neighbours
  std::vector<uint16_t> neighbNodes;
  std::vector<std::array<int, 18>> edgePattern;
  determineTrisoupNeighbours(nodes, neighbNodes, edgePattern, blockWidth);;


  // Decode vertex presence and position into bitstream  
  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  decodeTrisoupVertices(segind, vertices, neighbNodes, edgePattern, bitDropped, gps, gbh, arithmeticDecoder);


  PCCPointSet3 recPointCloud;
  recPointCloud.addRemoveAttributes(pointCloud);

  // Compute refinedVertices.
  std::vector<CentroidDrift> drifts;
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;
  bool haloFlag = gbh.trisoup_halo_flag;
  bool adaptiveHaloFlag = gbh.trisoup_adaptive_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;
  decodeTrisoupCommon(
    nodes, segind, vertices, drifts, pointCloud, recPointCloud, gps, gbh,
    blockWidth, maxval, gbh.trisoup_sampling_value_minus1 + 1, bitDropped,
    isCentroidDriftActivated, true, haloFlag, adaptiveHaloFlag, fineRayFlag,
    &arithmeticDecoder);

  pointCloud.resize(0);
  pointCloud = std::move(recPointCloud);
}


//============================================================================
void determineTrisoupNeighbours(
  const ringbuf<PCCOctree3Node>& leaves,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  const int defaultBlockWidth) {

  // Width of block.
  // in future, may override with leaf blockWidth
  const int32_t blockWidth = defaultBlockWidth;

  // Eight corners of block.
  const Vec3<int32_t> pos000({ 0, 0, 0 });
  const Vec3<int32_t> posW00({ blockWidth, 0, 0 });
  const Vec3<int32_t> pos0W0({ 0, blockWidth, 0 });
  const Vec3<int32_t> posWW0({ blockWidth, blockWidth, 0 });
  const Vec3<int32_t> pos00W({ 0, 0, blockWidth });
  const Vec3<int32_t> posW0W({ blockWidth, 0, blockWidth });
  const Vec3<int32_t> pos0WW({ 0, blockWidth, blockWidth });
  const Vec3<int32_t> posWWW({ blockWidth, blockWidth, blockWidth });

  // Put all leaves' edges into a list.
  std::vector<TrisoupSegmentNeighbours> segments;
  segments.reserve(36 * leaves.size());
  for (int i = 0; i < leaves.size(); i++) {

    const auto& leaf = leaves[i];
    int ii = 36 * i;
    int ii2 = ii + 12;
    int ii3 = ii + 24;
    // x: left to right; y: bottom to top; z: far to near
    auto posNode = leaf.pos + blockWidth;

    // ------------ edges along x
    // in node
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos000, posNode + posW00, ii + 0, 1 })); // far bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos0W0, posNode + posWW0, ii + 2,  2 })); // far top edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos00W, posNode + posW0W, ii + 8,  4 })); // near bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos0WW, posNode + posWWW, ii + 10,  8 })); // near top edge
    // left 
    auto posLeft = posNode - posW00;
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos000, posLeft + posW00, ii2 + 0, 16 })); // far bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos0W0, posLeft + posWW0, ii2 + 2,  32 })); // far top edge
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos00W, posLeft + posW0W, ii2 + 8,  64 })); // near bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos0WW, posLeft + posWWW, ii2 + 10,  128 })); // near top edge
    //right
    auto posRight = posNode + posW00;
    segments.push_back(TrisoupSegmentNeighbours({ posRight + pos000, posRight + posW00, ii3 + 0, 256 })); // far bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posRight + pos0W0, posRight + posWW0, ii3 + 2,  512 })); // far top edge
    segments.push_back(TrisoupSegmentNeighbours({ posRight + pos00W, posRight + posW0W, ii3 + 8,  1024 })); // near bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posRight + pos0WW, posRight + posWWW, ii3 + 10,  2048 })); // near top edge                                                                                                       

    // ------------ edges along y 
    // in node
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos000, posNode + pos0W0, ii + 1, 1 + (1 << 13) })); // far left edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + posW00, posNode + posWW0, ii + 3,  2 + (1 << 13) })); // far right edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos00W, posNode + pos0WW, ii + 9,  4 + (1 << 13) })); // near left edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + posW0W, posNode + posWWW, ii + 11,  8 + (1 << 13) })); // near right edge   
    // bottom
    auto posBottom = posNode - pos0W0;
    segments.push_back(TrisoupSegmentNeighbours({ posBottom + pos000, posBottom + pos0W0, ii2 + 1, 16 })); // far left edge
    segments.push_back(TrisoupSegmentNeighbours({ posBottom + posW00, posBottom + posWW0, ii2 + 3,  32 })); // far right edge
    segments.push_back(TrisoupSegmentNeighbours({ posBottom + pos00W, posBottom + pos0WW, ii2 + 9,  64 })); // near left edge
    segments.push_back(TrisoupSegmentNeighbours({ posBottom + posW0W, posBottom + posWWW, ii2 + 11,  128 })); // near right edge  
    // top
    auto posTop = posNode + pos0W0;
    segments.push_back(TrisoupSegmentNeighbours({ posTop + pos000, posTop + pos0W0, ii3 + 1, 256 })); // far left edge
    segments.push_back(TrisoupSegmentNeighbours({ posTop + posW00, posTop + posWW0, ii3 + 3,  512 })); // far right edge
    segments.push_back(TrisoupSegmentNeighbours({ posTop + pos00W, posTop + pos0WW, ii3 + 9,  1024 })); // near left edge
    segments.push_back(TrisoupSegmentNeighbours({ posTop + posW0W, posTop + posWWW, ii3 + 11,  2048 })); // near right edge 

    // ------------ edges along z 
    // in node
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos000, posNode + pos00W, ii + 4, 1 + (1 << 14) })); // bottom left edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos0W0, posNode + pos0WW, ii + 5,  2 + (1 << 14) })); // top left edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + posWW0, posNode + posWWW, ii + 6,  4 + (1 << 14) })); // top right edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + posW00, posNode + posW0W, ii + 7,  8 + (1 << 14) })); // bottom right edge    
    // near 
    auto posNear = posNode - pos00W;
    segments.push_back(TrisoupSegmentNeighbours({ posNear + pos000, posNear + pos00W, ii2 + 4, 16 })); // bottom left edge
    segments.push_back(TrisoupSegmentNeighbours({ posNear + pos0W0, posNear + pos0WW, ii2 + 5,  32 })); // top left edge
    segments.push_back(TrisoupSegmentNeighbours({ posNear + posWW0, posNear + posWWW, ii2 + 6,  64 })); // top right edge
    segments.push_back(TrisoupSegmentNeighbours({ posNear + posW00, posNear + posW0W, ii2 + 7,  128 })); // bottom right edge
    // far 
    auto posFar = posNode + pos00W;
    segments.push_back(TrisoupSegmentNeighbours({ posFar + pos000, posFar + pos00W, ii3 + 4, 256 })); // bottom left edge
    segments.push_back(TrisoupSegmentNeighbours({ posFar + pos0W0, posFar + pos0WW, ii3 + 5,  512 })); // top left edge
    segments.push_back(TrisoupSegmentNeighbours({ posFar + posWW0, posFar + posWWW, ii3 + 6,  1024 })); // top right edge
    segments.push_back(TrisoupSegmentNeighbours({ posFar + posW00, posFar + posW0W, ii3 + 7,  2048 })); // bottom right edge
  }

  // Sort the list and find unique segments.
  std::sort(segments.begin(), segments.end());

  // find neighbourgs for unique segments
  TrisoupSegmentNeighbours localSegment = segments[0];


  auto it = segments.begin() + 1;
  neighbNodes.clear();
  std::vector<int> correspondanceUnique(segments.size(), -1);

  int uniqueIndex = 0;
  std::array<int, 18> pattern = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

  for (; it != segments.end(); it++) {
    if (localSegment.startpos != it->startpos || localSegment.endpos != it->endpos) {

      if (localSegment.neighboursMask & 15) {
        // the segment is a true segment and not only accumulation of copies; then push and jump to next segment
        neighbNodes.push_back(localSegment.neighboursMask);
        edgePattern.push_back(pattern);

        uniqueIndex++;
        for (int v = 0; v < 18; v++)
          pattern[v] = -1;
      }
      localSegment = *it;
    }
    else {
      // segment[i] is the same as localSegment
      // Accumulate into localSegment       
      localSegment.neighboursMask |= it->neighboursMask;
    }
    correspondanceUnique[it->index] = uniqueIndex;

    // ---------- neighbouring vertex parallel before
    if (it->neighboursMask >= 256 && it->neighboursMask <= 2048) { // lookinbg for vertex before (x left or y front or z below)
      int indexBefore = it->index - 24;

      if (correspondanceUnique[indexBefore] != -1) {
        pattern[0] = correspondanceUnique[indexBefore];
      }
    }

    // ---------    8-bit pattern = 0 before, 1-4 perp, 5-12 others
    static const int localEdgeindex[12][11] = {
      { 4,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 0
      { 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 1
      { 1,  5,  4,  9,  0,  8, -1, -1, -1, -1, -1}, // vertex 2
      { 0,  7,  4,  8,  2, 10,  1,  9, -1, -1, -1}, // vertex 3
      {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 4
      { 1,  0,  9,  4, -1, -1, -1, -1, -1, -1, -1}, // vertex 5
      { 3,  2,  0, 10, 11,  9,  8,  7,  5,  4, -1}, // vertex 6
      { 0,  1,  2,  8, 10,  4,  5, -1, -1, -1, -1}, // vertex 7
      { 4,  9,  1,  0, -1, -1, -1, -1, -1, -1, -1}, // vertex 8
      { 4,  0,  1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 9
      { 5,  9,  1,  2,  8,  0, -1, -1, -1, -1, -1}, // vertex 10
      { 7,  8,  0, 10,  5,  2,  3,  9,  1, -1, -1}  // vertex 11
    };
    static const int patternIndex[12][11] = {
      { 3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 0
      { 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 1
      { 2,  3,  5,  8, 15, 17, -1, -1, -1, -1, -1}, // vertex 2
      { 2,  3,  5,  8,  9, 12, 15, 17, -1, -1, -1}, // vertex 3
      {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 4
      { 1,  7, 10, 14, -1, -1, -1, -1, -1, -1, -1}, // vertex 5
      { 1,  2,  6,  9, 10, 11, 13, 14, 15, 16, -1}, // vertex 6
      { 2,  5,  8,  9, 12, 15, 17, -1, -1, -1, -1}, // vertex 7
      { 1,  4,  7, 14, -1, -1, -1, -1, -1, -1, -1}, // vertex 8
      { 1,  7, 14, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 9
      { 1,  2,  6, 14, 15, 16, -1, -1, -1, -1, -1}, // vertex 10
      { 1,  2,  6,  9, 11, 13, 14, 15, 16, -1, -1}  // vertex 11
    };

    if ((it->neighboursMask & 4095) <= 8) { // true edge, not a copy; so done as many times a nodes for the  true edge
      int indexLow = it->index % 12;   // true edge index within node
      for (int v = 0; v < 11; v++) {
        if (localEdgeindex[indexLow][v] == -1)
          break;

        int indexV = it->index - indexLow + localEdgeindex[indexLow][v]; // index of segment
        int Vidx = correspondanceUnique[indexV];
        if (Vidx != -1)  // check if already coded
          pattern[patternIndex[indexLow][v]] = Vidx;
        else
          std::cout << "#" << indexLow << "/" << v << " ";

      }
    }

  }
  if (localSegment.neighboursMask & 15) {
    neighbNodes.push_back(localSegment.neighboursMask);
    edgePattern.push_back(pattern);
  }
}
//============================================================================

template<typename T>
Vec3<T>
crossProduct(const Vec3<T> a, const Vec3<T> b)
{
  Vec3<T> ret;
  ret[0] = a[1] * b[2] - a[2] * b[1];
  ret[1] = a[2] * b[0] - a[0] * b[2];
  ret[2] = a[0] * b[1] - a[1] * b[0];
  return ret;
}

//---------------------------------------------------------------------------

Vec3<int32_t>
truncate(const Vec3<int32_t> in, const int32_t offset)
{
  Vec3<int32_t> out = in + offset;  
  if (out[0] < 0)
    out[0] = 0;
  if (out[1] < 0)
    out[1] = 0;
  if (out[2] < 0)
    out[2] = 0;

  return out;
}



//---------------------------------------------------------------------------

int32_t
trisoupVertexArc(int32_t x, int32_t y, int32_t Width_x, int32_t Width_y)
{
  int32_t score;

  if (x >= Width_x) {
    score = y;
  } else if (y >= Width_y) {
    score = Width_y + Width_x - x;
  } else if (x <= 0) {
    score = Width_y*2 + Width_x - y;
  } else {
    score = Width_y*2 + Width_x + x;
  }

  return score;
}

//---------------------------------------------------------------------------
bool
boundaryinsidecheck(const Vec3<int32_t> a, const int bbsize)
{
  return a[0] >= 0 && a[0] <= bbsize && a[1] >= 0 && a[1] <= bbsize && a[2] >= 0 && a[2] <= bbsize;
}

//---------------------------------------------------------------------------

bool
rayIntersectsTriangle(
  const Vec3<int32_t>& rayOrigin,  
  const Vec3<int32_t>& TriangleVertex0,
  const Vec3<int32_t>& edge1,
  const Vec3<int32_t>& edge2,
  const Vec3<int32_t>& h,
  int32_t a,
  Vec3<int32_t>& outIntersectionPoint,
  Vec3<int32_t>& outIntersectionPointUp,
  Vec3<int32_t>& outIntersectionPointDown,
  int direction,
  int haloTriangle,
  int thickness)
{  
  Vec3<int32_t> s = rayOrigin - TriangleVertex0;
  int32_t u = (s * h) / a;  

  Vec3<int32_t> q = crossProduct(s, edge1) ;
  //int32_t v = (rayVector * q) / a;
  int32_t v =  q[direction]  / a;

  int w = kTrisoupFpOne - u - v;

  int32_t t = (edge2 * (q >> kTrisoupFpBits)) / a;  
  // outIntersectionPoint = rayOrigin + ((rayVector * t) >> kTrisoupFpBits);  
  outIntersectionPoint[direction] +=  t;

  outIntersectionPointUp = outIntersectionPoint;
  outIntersectionPointUp[direction] += thickness;

  outIntersectionPointDown = outIntersectionPoint;
  outIntersectionPointDown[direction] -= thickness;

  return u >= -haloTriangle && v >= -haloTriangle && w >= -haloTriangle;
}

//---------------------------------------------------------------------------

void nonCubicNode
(
 const GeometryParameterSet& gps,
 const GeometryBrickHeader& gbh,
 const Vec3<int32_t>& leafpos,
 const int32_t blockWidth,
 const Box3<int32_t>& bbox,
 Vec3<int32_t>& newp,
 Vec3<int32_t>& neww,
 Vec3<int32_t>* corner )
{
  bool flag_n = gps.non_cubic_node_start_edge && ( gbh.slice_bb_pos_bits   > 0 );
  bool flag_f = gps.non_cubic_node_end_edge   && ( gbh.slice_bb_width_bits > 0 );
  for( int k=0; k<3; k++ ) {
    newp[k] = ( ( flag_n ) && ( leafpos[k] < bbox.min[k] ) ) ? bbox.min[k] : leafpos[k];
    neww[k] = ( ( flag_n ) && ( leafpos[k] < bbox.min[k] ) ) ?
      (blockWidth-(bbox.min[k]-leafpos[k])) :
      ( flag_f ) ? std::min(bbox.max[k]-leafpos[k]+1, blockWidth) : blockWidth;
  }
  corner[POS_000] = {       0,       0,       0 };
  corner[POS_W00] = { neww[0],       0,       0 };
  corner[POS_0W0] = {       0, neww[1],       0 };
  corner[POS_WW0] = { neww[0], neww[1],       0 };
  corner[POS_00W] = {       0,       0, neww[2] };
  corner[POS_W0W] = { neww[0],       0, neww[2] };
  corner[POS_0WW] = {       0, neww[1], neww[2] };
  corner[POS_WWW] = { neww[0], neww[1], neww[2] };
  return;
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
  std::vector<CentroidDrift>& drifts,
  PCCPointSet3& pointCloud,
  PCCPointSet3& recPointCloud,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int defaultBlockWidth,
  int poistionClipValue,
  uint32_t samplingValue,
  const int bitDropped,
  const bool isCentroidDriftActivated,
  bool isDecoder,
  bool haloFlag,
  bool adaptiveHaloFlag,
  bool fineRayflag,
  pcc::EntropyDecoder* arithmeticDecoder)
{
  // clear drifst vecause of encoder multi pass 
  drifts.clear();

  Box3<int32_t> sliceBB;
  sliceBB.min = gbh.slice_bb_pos << gbh.slice_bb_pos_log2_scale;
  sliceBB.max = sliceBB.min + ( gbh.slice_bb_width << gbh.slice_bb_width_log2_scale );

  // Put all leaves' sgements into a list.
  std::vector<TrisoupSegment> segments;
  segments.resize(12*leaves.size());

  // Width of block.
  // in future, may override with leaf blockWidth
  const int32_t blockWidth = defaultBlockWidth;

  for (int i = 0; i < leaves.size(); i++) {
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode( gps, gbh, leaves[i].pos, defaultBlockWidth, sliceBB, nodepos, nodew, corner );

    int ii = 12 * i;
    for (int j = 0; j < 12; j++) {
      int iii = ii + j;
      segments[iii] = { nodepos, nodepos, iii,-1,-1 };
    }
    segments[ii + 0].startpos += corner[POS_000]; // far bottom edge
    segments[ii + 0].endpos += corner[POS_W00];
    segments[ii + 1].startpos += corner[POS_000]; // far left edge
    segments[ii + 1].endpos += corner[POS_0W0];
    segments[ii + 2].startpos += corner[POS_0W0]; // far top edge
    segments[ii + 2].endpos += corner[POS_WW0];
    segments[ii + 3].startpos += corner[POS_W00]; // far right edge
    segments[ii + 3].endpos += corner[POS_WW0];
    segments[ii + 4].startpos += corner[POS_000]; //bottom left edge
    segments[ii + 4].endpos += corner[POS_00W];
    segments[ii + 5].startpos += corner[POS_0W0]; // top left edge
    segments[ii + 5].endpos += corner[POS_0WW];
    segments[ii + 6].startpos += corner[POS_WW0]; //top right edge
    segments[ii + 6].endpos += corner[POS_WWW];
    segments[ii + 7].startpos += corner[POS_W00]; // bottom right edge
    segments[ii + 7].endpos += corner[POS_W0W];
    segments[ii + 8].startpos += corner[POS_00W]; // near bottom edge
    segments[ii + 8].endpos += corner[POS_W0W];
    segments[ii + 9].startpos += corner[POS_00W]; // near left edge
    segments[ii + 9].endpos += corner[POS_0WW];
    segments[ii + 10].startpos += corner[POS_0WW]; // near top edge
    segments[ii + 10].endpos += corner[POS_WWW];
    segments[ii + 11].startpos += corner[POS_W0W]; // near right edge
    segments[ii + 11].endpos += corner[POS_WWW];
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

  // contexts for drift centroids 
  AdaptiveBitModel ctxDrift0[9];
  AdaptiveBitModel ctxDriftSign[3][8][8];
  AdaptiveBitModel ctxDriftMag[4];

  // Create list of refined vertices, one leaf at a time.
  std::vector<Vec3<int32_t>> refinedVertices;


  // ----------- loop on leaf nodes ----------------------
  for (int i = 0; i < leaves.size(); i++) {
    Vec3<int32_t> nodepos, nodew, corner[8];
    nonCubicNode( gps, gbh, leaves[i].pos, defaultBlockWidth, sliceBB, nodepos, nodew, corner );

    // Find up to 12 vertices for this leaf.
    std::vector<Vertex> leafVertices;
    std::vector<Vec3<int32_t>> refinedVerticesBlock;
    refinedVerticesBlock.reserve(blockWidth * blockWidth * 4);

    for (int j = 0; j < 12; j++) {
      TrisoupSegment& segment = segments[i * 12 + j];
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
      int32_t distance = (segment.vertex << (kTrisoupFpBits + bitDropped)) + (kTrisoupFpHalf << bitDropped);
      if (direction[0])
        point[0] += distance; // in {0,1,...,B-1}
      else if (direction[1])
        point[1] += distance;
      else  // direction[2]
        point[2] += distance;

      leafVertices.push_back({ point, 0, 0 });

      // Add vertex to list of points.
      if (bitDropped || samplingValue > 1) {
        Vec3<int32_t> foundvoxel = (point + truncateValue) >> kTrisoupFpBits;
        if (boundaryinsidecheck(foundvoxel, blockWidth - 1))
          refinedVerticesBlock.push_back(nodepos + foundvoxel);
      }

    }

    // Skip leaves that have fewer than 3 vertices.
    int triCount = (int)leafVertices.size();
    if (triCount < 3) {
      std::sort(refinedVerticesBlock.begin(), refinedVerticesBlock.end());
      refinedVerticesBlock.erase(std::unique(refinedVerticesBlock.begin(), refinedVerticesBlock.end()), refinedVerticesBlock.end());
      refinedVertices.insert(refinedVertices.end(), refinedVerticesBlock.begin(), refinedVerticesBlock.end());
      continue;
    }

    // compute centroid
    Vec3<int32_t> blockCentroid = 0;
    for (int j = 0; j < triCount; j++) {
      blockCentroid += leafVertices[j].pos;
    }
    blockCentroid /= triCount; 

    // order vertices along a dominant axis only if more than three (otherwise only one triangle, whatever...)
    int dominantAxis = findDominantAxis(leafVertices, nodew, blockCentroid);

    // compute centroid
    std::vector<int> Weigths(leafVertices.size(), 0);
    int Wtotal = 0;
    for (int k = 0; k < triCount; k++) {
      int k2 = k + 1;
      if (k2 >= triCount)
        k2 -= triCount;
      Vec3<int32_t> segment = (leafVertices[k].pos - leafVertices[k2].pos).abs();
      int weight = segment[0] + segment[1] + segment[2];

      Weigths[k] += weight;
      Weigths[k2] += weight;
      Wtotal += 2 * weight;
    }

    Vec3<int64_t> blockCentroid2 = 0;
    for (int j = 0; j < triCount; j++) {
      blockCentroid2 += int64_t(Weigths[j]) * leafVertices[j].pos;
    }
    blockCentroid2 /= int64_t(Wtotal);
    blockCentroid = { int(blockCentroid2[0]),int(blockCentroid2[1]), int(blockCentroid2[2]) };

    // Refinement of the centroid along the domiannt axis
    // deactivated if sampling is too big
    if (triCount > 3 && isCentroidDriftActivated && samplingValue <= 4) {

      int bitDropped2 = bitDropped;
      int halfDropped2 = bitDropped2 == 0 ? 0 : 1 << bitDropped2 - 1;

      // contextual information  for drift coding 
      int minPos = leafVertices[0].pos[dominantAxis];
      int maxPos = leafVertices[0].pos[dominantAxis];
      for (int k = 1; k < triCount; k++) {
        if (leafVertices[k].pos[dominantAxis] < minPos)
          minPos = leafVertices[k].pos[dominantAxis];
        if (leafVertices[k].pos[dominantAxis] > maxPos)
          maxPos = leafVertices[k].pos[dominantAxis];
      }

      // find normal vector 
      Vec3<int64_t> accuNormal = 0;
      for (int k = 0; k < triCount; k++) {
        int k2 = k + 1;
        if (k2 >= triCount)
          k2 -= triCount;        
        accuNormal += crossProduct(leafVertices[k].pos - blockCentroid, leafVertices[k2].pos - blockCentroid);
      }
      int64_t normN = isqrt(accuNormal[0]* accuNormal[0] + accuNormal[1] * accuNormal[1] + accuNormal[2] * accuNormal[2]);       
      Vec3<int32_t> normalV = (accuNormal<< kTrisoupFpBits) / normN;
     
      //  drift bounds     
      int ctxMinMax = std::min(8, (maxPos - minPos) >> (kTrisoupFpBits + bitDropped));
      int bound = (int(nodew[dominantAxis]) - 1) << kTrisoupFpBits;
      int m = 1;
      for (; m < nodew[dominantAxis]; m++) {
        Vec3<int32_t> temp = blockCentroid + m * normalV;
        if (temp[0]<0 || temp[1]<0 || temp[2]<0 || temp[0]>bound || temp[1]>bound || temp[2]> bound)
          break;
      }
      int highBound = (m - 1) + halfDropped2 >> bitDropped2;

      m = 1;
      for (; m < nodew[dominantAxis]; m++) {
        Vec3<int32_t> temp = blockCentroid - m * normalV;
        if (temp[0]<0 || temp[1]<0 || temp[2]<0 || temp[0]>bound || temp[1]>bound || temp[2]> bound)
          break;
      }
      int lowBound = (m - 1) + halfDropped2 >> bitDropped2;
      int lowBoundSurface =  ((blockCentroid[dominantAxis] - minPos) + kTrisoupFpHalf >> kTrisoupFpBits)   + halfDropped2 >> bitDropped2;
      int highBoundSurface =  ((maxPos - blockCentroid[dominantAxis]) + kTrisoupFpHalf >> kTrisoupFpBits) + halfDropped2 >> bitDropped2;

      int driftQ = 0;
      if (!isDecoder) { // encoder 
        // determine qauntized drift 
        int counter = 0;
        int drift = 0;
        int maxD = std::max(int(samplingValue), bitDropped2);

        // determine quantized drift         
        for (int p = leaves[i].start; p < leaves[i].end; p++) {
          auto point = (pointCloud[p] - nodepos) << kTrisoupFpBits;

         Vec3<int32_t> CP = crossProduct(normalV, point - blockCentroid) >> kTrisoupFpBits;
         int dist = std::max(std::max(std::abs(CP[0]) , std::abs(CP[1])) , std::abs(CP[2]));
         dist >>= kTrisoupFpBits;
          
          if (dist <= maxD) {
            int w = 1 + 4 * (maxD - dist);
            counter += w;
            drift += w * ( (normalV * (point - blockCentroid)) >> kTrisoupFpBits );
          }
        }

        if (counter) { // drift is shift by kTrisoupFpBits                    
          drift = (drift >> kTrisoupFpBits - 6) / counter; // drift is shift by 6 bits 
        }

        int half = 1 << 5 + bitDropped2;
        int DZ = 2*half/3; 

        if (abs(drift) >= DZ) {
          driftQ = (abs(drift) - DZ + 2*half) >> 6 + bitDropped2;
          if (drift < 0)
            driftQ = -driftQ;
        }
        driftQ = std::min(std::max(driftQ, -lowBound), highBound);  // drift in [-lowBound; highBound]
        // push quantized drift  to buffeer for encoding 
        drifts.push_back({ driftQ, lowBound, highBound, ctxMinMax, lowBoundSurface, highBoundSurface });       

      } // end encoder 

      else { // decode drift 

        driftQ = arithmeticDecoder->decode(ctxDrift0[ctxMinMax]) ? 0 : 1;

        // if not 0, drift in [-lowBound; highBound]
        if (driftQ) {
          // code sign
          int lowS = std::min(7,lowBoundSurface);
          int highS = std::min(7,highBoundSurface);

          int sign = 1;
          if (highBound && lowBound) // otherwise sign is knwow 
            sign = arithmeticDecoder->decode(ctxDriftSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)][lowS][highS]);
          else if (!highBound) // highbound is 0 , so sign is negative; otherwise sign is already set to positive 
            sign = 0;

          // code remaining bits 1 to 7 at most 
          int magBound = (sign ? highBound : lowBound) - 1;

          int ctx = 0;
          while (magBound > 0) {
            int bit;
            if (ctx < 4)
              bit = arithmeticDecoder->decode(ctxDriftMag[ctx]);
            else
              bit = arithmeticDecoder->decode();

            if (bit) // magDrift==0 and magnitude coding is finished 
              break;

            driftQ++;
            magBound--;
            ctx++;
          }

          if (!sign)
            driftQ = -driftQ;

        } // end decoder 
      }

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
      blockCentroid[0] = std::min(((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1, blockCentroid[0]);
      blockCentroid[1] = std::min(((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1, blockCentroid[1]);
      blockCentroid[2] = std::min(((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1, blockCentroid[2]);
    } // end refinement of the centroid


    // Divide vertices into triangles around centroid
    // and upsample each triangle by an upsamplingFactor.    
    Vec3<int32_t> v2 = triCount == 3 ? leafVertices[2].pos : blockCentroid;
    Vec3<int32_t> v1 = leafVertices[0].pos;

    Vec3<int32_t> posNode = nodepos << kTrisoupFpBits;

    if (triCount > 3) {
      Vec3<int32_t> foundvoxel = (blockCentroid + truncateValue) >> kTrisoupFpBits;
      if (boundaryinsidecheck(foundvoxel, blockWidth - 1))
        refinedVerticesBlock.push_back(nodepos + foundvoxel);
    }

    for (int triIndex = 0; triIndex < (triCount == 3 ? 1 : triCount); triIndex++) {
      int j0 = triIndex;
      int j1 = triIndex + 1;
      if (j1 >= triCount)
        j1 -= triCount;

      Vec3<int32_t> v0 = v1;
      v1 = leafVertices[j1].pos;

      // range      
      int minRange[3];
      int maxRange[3];
      for (int k = 0; k < 3; k++) {
        minRange[k] = std::max(0, std::min(std::min(v0[k], v1[k]), v2[k]) + truncateValue >> kTrisoupFpBits);
        maxRange[k] = std::min(127, std::max(std::max(v0[k], v1[k]), v2[k]) + truncateValue >> kTrisoupFpBits);
      }

      // precompute for rays 
      Vec3<int32_t> edge1 = v1 - v0;
      Vec3<int32_t> edge2 = v2 - v0;
      int minDir = 1 << 28;
      int directionExcluded = 0;
      for (int k = 0; k <= 2; k++) {
        Vec3<int32_t> rayVector = 0;
        rayVector[k] = 1 << kTrisoupFpBits;
        Vec3<int32_t> h = crossProduct(edge1, edge2) >> kTrisoupFpBits;
        int32_t a = (rayVector * h) >> kTrisoupFpBits;
        if (std::abs(a) < minDir) {
          minDir = std::abs(a);
          directionExcluded = k;
        }
      }

      // applying ray tracing along direction
      for (int direction = 0; direction < 3; direction++) {
        if (directionExcluded == direction) // exclude most parallel direction
          continue;
        rayTracingAlongdirection(
          refinedVerticesBlock, direction, samplingValue, bitDropped, blockWidth, nodepos , minRange,
          maxRange, edge1, edge2, v0, haloFlag,
          adaptiveHaloFlag, fineRayflag);
      }

    }  // end loop on triangles

    std::sort(refinedVerticesBlock.begin(), refinedVerticesBlock.end());
    refinedVerticesBlock.erase(std::unique(refinedVerticesBlock.begin(), refinedVerticesBlock.end()), refinedVerticesBlock.end());
    refinedVertices.insert(refinedVertices.end(), refinedVerticesBlock.begin(), refinedVerticesBlock.end());

  }// end loop on leaves

  // Move list of points to pointCloud.
  recPointCloud.resize(refinedVertices.size());
  for (int i = 0; i < refinedVertices.size(); i++) {
    recPointCloud[i] = refinedVertices[i];
  }
}


// ---------------------------------------------------------------------------
void decodeTrisoupVertices(
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  int bitDropped,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyDecoder& arithmeticDecoder) {
  const int nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;
  const int max2bits = nbitsVertices > 1 ? 3 : 1;
  const int mid2bits = nbitsVertices > 1 ? 2 : 1;

  int iV = 0;
  std::vector<int> correspondanceSegment2V;

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
    bool c = MapOBUFTriSoup[0].decodeEvolve(&arithmeticDecoder, ctxTriSoup, ctxMap2, ctxMap1, &_OBUFleafNumber, _BufferOBUFleaves);

    segind.push_back(c);
    correspondanceSegment2V.push_back(-1);

    // encode position vertex 
    if (c) {
      correspondanceSegment2V.back() = iV;

      uint8_t v = 0;
      int ctxFullNbounds = (4 * (ctx0 <= 1 ? 0 : (ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctx1) - 1)) * 2 + (ctxE == 3);
      int b = nbitsVertices - 1;

      // first bit
      ctxMap1 = ctxFullNbounds * 2 + (nclosestStart > 0);
      ctxMap2 = missedCloseStart << 8;
      ctxMap2 |= (patternClosest & 1) << 7;
      ctxMap2 |= direction << 5;
      ctxMap2 |= patternClose & (0b00011111);
      int orderedPclosePar = (((patternClose >> 5) & 3) << 2) + (!!(patternClose & 128) << 1) + !!(patternClose & 256);

      int bit = MapOBUFTriSoup[1].decodeEvolve(&arithmeticDecoder, ctxTriSoup, ctxMap2, ctxMap1, &_OBUFleafNumber, _BufferOBUFleaves);
      v = (v << 1) | bit;
      b--;

      // second bit
      if (b >= 0) {
        ctxMap1 = ctxFullNbounds * 2 + (nclosestStart > 0);
        ctxMap2 = missedCloseStart << 8;
        ctxMap2 |= (patternClose & 1) << 7;
        ctxMap2 |= (patternClosest & 1) << 6;
        ctxMap2 |= direction << 4;
        ctxMap2 |= (patternClose & (0b00011111)) >> 1;
        ctxMap2 = (ctxMap2 << 4) + orderedPclosePar;

        bit = MapOBUFTriSoup[2].decodeEvolve(&arithmeticDecoder, ctxTriSoup, ctxMap2, (ctxMap1 << 1) + v, &_OBUFleafNumber, _BufferOBUFleaves);
        v = (v << 1) | bit;
        b--;
      }


      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 = (5 * (ctx0 >> 1) + missedCloseStart) * 2 + (ctxE == 3);
        v = (v << 1) | arithmeticDecoder.decode(ctxTempV2[4 * ctxFullNboundsReduced1 + v]);
        b--;
      }

      // remaining bits are bypassed
      for (; b >= 0; b--)
        v = (v << 1) | arithmeticDecoder.decode();
      vertices.push_back(v);
      iV++;
    }

  }

}


//-----------------------
// Project vertices along dominant axis (i.e., into YZ, XZ, or XY plane).
// Sort projected vertices by decreasing angle in [-pi,+pi] around center
// of block (i.e., clockwise) breaking ties in angle by
// increasing distance along the dominant axis.   

int findDominantAxis(
  std::vector<Vertex>& leafVertices,
  Vec3<uint32_t> blockWidth,
  Vec3<int32_t> blockCentroid ) {

  int dominantAxis = 0;
  int triCount = leafVertices.size();
  if (triCount > 3) {
    Vertex vertex;
    Vec3<int32_t> Width = blockWidth << kTrisoupFpBits;

    const int sIdx1[3] = { 2,2,1 };
    const int sIdx2[3] = { 1,0,0 };

    int maxNormTri = 0;
    for (int axis = 0; axis <= 2; axis++) {
      // order along axis
      for (int j = 0; j < triCount; j++) {
        Vec3<int32_t> s = leafVertices[j].pos + kTrisoupFpHalf; // back to [0,B]^3 for ordering
        leafVertices[j].theta = trisoupVertexArc(s[sIdx1[axis]], s[sIdx2[axis]],
                                                 Width[sIdx1[axis]], Width[sIdx2[axis]]);
        leafVertices[j].tiebreaker = s[axis];
      }
      std::sort(leafVertices.begin(), leafVertices.end(), vertex);

      // compute sum normal
      int32_t accuN = 0;
      for (int k = 0; k < triCount; k++) {
        int k2 = k + 1;
        if (k2 >= triCount)
          k2 -= triCount;
        Vec3<int32_t> h = crossProduct(leafVertices[k].pos - blockCentroid, leafVertices[k2].pos - blockCentroid);
        accuN += std::abs(h[axis]);
      }

      // if sumnormal is bigger , this is dominantAxis
      if (accuN > maxNormTri) {
        maxNormTri = accuN;
        dominantAxis = axis;
      }
    }

    for (int j = 0; j < leafVertices.size(); j++) {
      Vec3<int32_t> s = leafVertices[j].pos + kTrisoupFpHalf; // back to [0,B]^3 for ordering
      leafVertices[j].theta = trisoupVertexArc(s[sIdx1[dominantAxis]], s[sIdx2[dominantAxis]],
                                               Width[sIdx1[dominantAxis]], Width[sIdx2[dominantAxis]]);
      leafVertices[j].tiebreaker = s[dominantAxis];
    }
    std::sort(leafVertices.begin(), leafVertices.end(), vertex);
  } // end find dominant axis 

  return dominantAxis;
}


// -------------------------------------------
void rayTracingAlongdirection(
  std::vector<Vec3<int32_t>>& refinedVerticesBlock,
  int direction,
  uint32_t samplingValue,
  int bitDropped,
  int blockWidth,
  Vec3<int32_t> nodepos,
  int minRange[3],
  int maxRange[3],
  Vec3<int32_t> edge1,
  Vec3<int32_t> edge2,
  Vec3<int32_t> v0,
  bool haloFlag,
  bool adaptiveHaloFlag,
  bool fineRayflag) {

  // check if ray tracing is valid; if not skip the direction
  Vec3<int32_t> rayVector = 0;
  rayVector[direction] = 1 << kTrisoupFpBits;
  Vec3<int32_t> h = crossProduct(rayVector, edge2) >> kTrisoupFpBits;
  int32_t a = (edge1 * h) >> kTrisoupFpBits;
  if (std::abs(a) <= kTrisoupFpOne)
    return;

  //bounds
  const int g1pos[3] = { 1, 0, 0 };
  const int g2pos[3] = { 2, 2, 1 };
  const int32_t startposG1 = minRange[g1pos[direction]];
  const int32_t startposG2 = minRange[g2pos[direction]];
  const int32_t endposG1 = maxRange[g1pos[direction]];
  const int32_t endposG2 = maxRange[g2pos[direction]];
  const int32_t rayStart = minRange[direction] << kTrisoupFpBits;
  Vec3<int32_t>  rayOrigin = rayStart;


  // ray tracing
  int haloTriangle = 0;
  int haloBit = (((1 << bitDropped) - 1) << kTrisoupFpBits) / blockWidth; // 28
  haloBit = (haloBit * 24) / 32;
  haloBit = haloBit > 40 ? 40 : haloBit;

  if (haloFlag) {
    if (samplingValue > 1) {
      haloTriangle = haloFlag ? (adaptiveHaloFlag ? 50 * samplingValue : 50) : 0;
      haloTriangle = haloTriangle > 100 ? 100 : haloTriangle;
    }
    else {
      haloTriangle = haloBit;
    }
  }

  int thickness = samplingValue > 1 ? 16 : 32;

  for (int32_t g1 = startposG1; g1 <= endposG1; g1 += samplingValue) {
    rayOrigin[g1pos[direction]] = g1 << kTrisoupFpBits;


    for (int32_t g2 = startposG2; g2 <= endposG2; g2 += samplingValue) {
      rayOrigin[g2pos[direction]] = g2 << kTrisoupFpBits;

      // middle ray at integer position 
      Vec3<int32_t>  intersection = rayOrigin;
      Vec3<int32_t> intersectionUp = rayOrigin;
      Vec3<int32_t> intersectionDown = rayOrigin;
      bool foundIntersection = rayIntersectsTriangle(rayOrigin, v0, edge1, edge2, h, a, intersection, intersectionUp, intersectionDown, direction, haloTriangle, thickness);
      if (foundIntersection) {
        Vec3<int32_t>foundvoxel;

        if (true || samplingValue == 1) {
          foundvoxel = (intersectionUp + truncateValue) >> kTrisoupFpBits;
          if (boundaryinsidecheck(foundvoxel, blockWidth-1)) {
            refinedVerticesBlock.push_back(nodepos + foundvoxel);
          }
          foundvoxel = (intersectionDown + truncateValue) >> kTrisoupFpBits;
          if (boundaryinsidecheck(foundvoxel, blockWidth-1)) {
            refinedVerticesBlock.push_back(nodepos + foundvoxel);
          }
        }

        foundvoxel = (intersection + truncateValue) >> kTrisoupFpBits;
        if (boundaryinsidecheck(foundvoxel, blockWidth-1)) {
          refinedVerticesBlock.push_back(nodepos + foundvoxel);
          continue; // ray interected , no need to launch other rays
        }

      }

      // if ray not interected then  augment +- offset
      if (samplingValue == 1 && fineRayflag) {
        const int Offset1[8] = { 0,  0, -1, +1, -1, -1, +1, +1 };
        const int Offset2[8] = { -1, +1,  0,  0, -1, +1, -1, +1 };
        const int offset = kTrisoupFpHalf >> 2;

        for (int pos = 0; pos < 8; pos++) {

          Vec3<int32_t> rayOrigin2 = rayOrigin;
          rayOrigin2[g1pos[direction]] += Offset1[pos] * offset;
          rayOrigin2[g2pos[direction]] += Offset2[pos] * offset;

          Vec3<int32_t> intersection = rayOrigin2;
          if (rayIntersectsTriangle(rayOrigin2, v0, edge1, edge2, h, a, intersection, intersectionUp, intersectionDown, direction, haloTriangle, thickness)) {
            Vec3<int32_t> foundvoxel = (intersection + truncateValue) >> kTrisoupFpBits;
            if (boundaryinsidecheck(foundvoxel, blockWidth-1)) {
              refinedVerticesBlock.push_back(nodepos + foundvoxel);
              break; // ray interected , no need to launch other rays
            }
          }

        } //pos

      } // augment


    }// loop g2 
  }//loop g1

}



//============================================================================

}  // namespace pcc
