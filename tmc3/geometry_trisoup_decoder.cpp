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
  PCCPointSet3& predPointCloud,
  const Vec3<int> minimum_position)
{
  // trisoup uses octree coding until reaching the triangulation level.
  // todo(df): pass trisoup node size rather than 0?
  pcc::ringbuf<PCCOctree3Node> nodes;

  decodeGeometryOctree(
    gps, gbh, 0, pointCloud, ctxtMemOctree, arithmeticDecoder, &nodes,
    predPointCloud, minimum_position);

  int blockWidth = 1 << gbh.trisoupNodeSizeLog2(gps);
  const int maxVertexPrecisionLog2 = gbh.trisoup_vertex_quantization_bits ? gbh.trisoup_vertex_quantization_bits : gbh.trisoupNodeSizeLog2(gps);
  const int bitDropped =  std::max(0, gbh.trisoupNodeSizeLog2(gps) - maxVertexPrecisionLog2);
  const bool isCentroidDriftActivated = gbh.trisoup_centroid_vertex_residual_flag;


  // Determine neighbours
  std::vector<uint16_t> neighbNodes;
  std::vector<int> indexBefore;
  std::vector<std::vector<int>> perpVertexStart;
  determineTrisoupNeighbours(nodes, neighbNodes, indexBefore, perpVertexStart, blockWidth);


  // Decode vertex presence and position into bitstream  
  std::vector<bool> segind;
  std::vector<uint8_t> vertices;
  decodeTrisoupVertices(segind, vertices, neighbNodes, indexBefore, perpVertexStart, bitDropped, gps, gbh, arithmeticDecoder); 


  // Compute refinedVertices.
  std::vector<CentroidDrift> drifts;
  int32_t maxval = (1 << gbh.maxRootNodeDimLog2) - 1;
  bool haloFlag = gbh.trisoup_halo_flag;
  bool fineRayFlag = gbh.trisoup_fine_ray_tracing_flag;
  decodeTrisoupCommon(
    nodes, segind, vertices, drifts, pointCloud, blockWidth, maxval,
    gbh.trisoup_sampling_value_minus1 + 1, bitDropped, isCentroidDriftActivated, true, haloFlag, fineRayFlag, &arithmeticDecoder);
}


//============================================================================
void determineTrisoupNeighbours(
  const ringbuf<PCCOctree3Node>& leaves,
  std::vector<uint16_t>& neighbNodes,
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
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
    auto posNode = leaf.pos;

    // ------------ edges along x
    // in node
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos000, posNode + posW00, ii + 0, 1})); // far bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos0W0, posNode + posWW0, ii + 2,  2})); // far top edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos00W, posNode + posW0W, ii + 8,  4})); // near bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posNode + pos0WW, posNode + posWWW, ii + 10,  8})); // near top edge 
    // left 
    auto posLeft = posNode - posW00;
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos000, posLeft + posW00, ii2 + 0, 16 })); // far bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos0W0, posLeft + posWW0, ii2 + 2,  32})); // far top edge
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos00W, posLeft + posW0W, ii2 + 8,  64})); // near bottom edge
    segments.push_back(TrisoupSegmentNeighbours({ posLeft + pos0WW, posLeft + posWWW, ii2 + 10,  128})); // near top edge
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
  indexBefore.clear();
  std::vector<int> correspondanceUnique(segments.size(), -1);

  int uniqueIndex = 0; 
  int uniqueIndexBefore = -1;
  std::vector<int> uniquePerpVertexStart;  
  
  for (; it != segments.end() ; it++) {
    if (localSegment.startpos != it->startpos || localSegment.endpos != it->endpos) {      
      
      if (localSegment.neighboursMask & 15) {
        neighbNodes.push_back(localSegment.neighboursMask);
        indexBefore.push_back(uniqueIndexBefore); 
        perpVertexStart.push_back(uniquePerpVertexStart);       
        
        uniqueIndex++;
        uniqueIndexBefore = -1;
        uniquePerpVertexStart.clear();
      }
      localSegment =*it;     
    } 
    else {
      // segment[i] is the same as localSegment
      // Accumulate into localSegment       
      localSegment.neighboursMask |= it->neighboursMask;
    }   
    correspondanceUnique[it->index] = uniqueIndex;

    
    if (it->neighboursMask >=256 && it->neighboursMask <= 2048) { // lookinbg for vertex before (x left or y front or z below) 
      int indexBefore = it->index - 24;
      
      if (correspondanceUnique[indexBefore] != -1) {
        uniqueIndexBefore = correspondanceUnique[indexBefore];       
      }      
    }
    static const int perpVertexStartIndex[12][2] =       { {1,4}, {0,4}, {1,5}, {0,7}, {0,1},{1,2},{2,3},{0,3},{4,9},{4,8},{5,9},{7,8} };
    static const int perpVertexStartOrientation[12][2] = { {0,0}, {0,0}, {1,0}, {1,0}, {0,0},{1,0},{1,1},{1,0},{1,0},{1,0},{1,1},{1,1} }; // 0=away from, 1=toward  StartPos

    if (( it->neighboursMask  & 4095)<=8) { // looking for perpendicalar vertices from startPoint 
      const int maskIdx = (1 << 30) - 1;
      int indexLow = it->index % 12;   
      
      int indexPerp = it->index  - indexLow + perpVertexStartIndex[indexLow][0];
      int Vidx = correspondanceUnique[indexPerp];
      if (Vidx != -1) {
        bool flag = false;
        for (int k = 0; k < uniquePerpVertexStart.size(); k++) {
          flag = flag || Vidx == (uniquePerpVertexStart[k] & maskIdx);
        }
        if (!flag)
          uniquePerpVertexStart.push_back(Vidx + (perpVertexStartOrientation[indexLow][0] << 30));
      }
     
       
      indexPerp = it->index  - indexLow + perpVertexStartIndex[indexLow][1];
      Vidx = correspondanceUnique[indexPerp];
      if (Vidx != -1) {
        bool flag = false;
        for (int k = 0; k < uniquePerpVertexStart.size(); k++) {
          flag = flag || Vidx == (uniquePerpVertexStart[k] & maskIdx);
        }
        if (!flag)
          uniquePerpVertexStart.push_back(Vidx + (perpVertexStartOrientation[indexLow][1] << 30));
      }
      
    }

  }
  if (localSegment.neighboursMask & 15) {
    neighbNodes.push_back(localSegment.neighboursMask);
    indexBefore.push_back(uniqueIndexBefore);
    perpVertexStart.push_back(uniquePerpVertexStart);
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
trisoupVertexArc(int32_t x, int32_t y, int32_t Width)
{
  int32_t score;

  if (x >= Width) {
    score = y;   
  } else if (y >= Width) {
    score = 2*Width - x ;
  } else if (x <= 0) {
    score = 3*Width -y ;
  }
  else {
    score = 3 * Width + x;
  }

  return score;
}

//---------------------------------------------------------------------------
bool
boundaryinsidecheck(const Vec3<int32_t> a, const int bbsize)
{
  if (a[0] < 0 || a[0] > bbsize)
    return false;
  if (a[1] < 0 || a[1] > bbsize)
    return false;
  if (a[2] < 0 || a[2] > bbsize)
    return false;
  return true;
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
  int direction,
  int haloTriangle)
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

  return u >= -haloTriangle && v >= -haloTriangle && w >= -haloTriangle;
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
  int defaultBlockWidth,
  int poistionClipValue,
  uint32_t samplingValue,
  const int bitDropped,
  const bool isCentroidDriftActivated,
  bool isDecoder,
  bool haloFlag,
  bool fineRayflag,
  pcc::EntropyDecoder* arithmeticDecoder)
{
  // clear drifst vecause of encoder multi pass 
  drifts.clear();
  
  // Put all leaves' sgements into a list.
  std::vector<TrisoupSegment> segments;
  segments.resize(12*leaves.size());

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

  for (int i = 0; i < leaves.size(); i++) {
   
    int ii = 12 * i;
    auto pos = leaves[i].pos;
    for (int j = 0; j < 12; j++) {      
      int iii = ii + j;
      segments[iii] = { pos, pos, iii,-1,-1 };      
    }
   
    segments[ii + 0].startpos += pos000; // far bottom edge
    segments[ii + 0].endpos += posW00;
    segments[ii + 1].startpos += pos000; // far left edge
    segments[ii + 1].endpos += pos0W0;
    segments[ii + 2].startpos += pos0W0; // far top edge
    segments[ii + 2].endpos += posWW0;
    segments[ii + 3].startpos += posW00; // far right edge
    segments[ii + 3].endpos += posWW0;
    segments[ii + 4].startpos += pos000; //bottom left edge
    segments[ii + 4].endpos += pos00W;
    segments[ii + 5].startpos += pos0W0; // top left edge
    segments[ii + 5].endpos += pos0WW;       
    segments[ii + 6].startpos += posWW0; //top right edge
    segments[ii + 6].endpos += posWWW;
    segments[ii + 7].startpos += posW00; // bottom right edge
    segments[ii + 7].endpos += posW0W;   
    segments[ii + 8].startpos += pos00W; // near bottom edge
    segments[ii + 8].endpos += posW0W;
    segments[ii + 9].startpos += pos00W; // near left edge
    segments[ii + 9].endpos += pos0WW;
    segments[ii + 10].startpos += pos0WW; // near top edge
    segments[ii + 10].endpos += posWWW;
    segments[ii + 11].startpos += posW0W; // near right edge
    segments[ii + 11].endpos += posWWW;    
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
    uint32_t blockWidth = 0;

    // Find up to 12 vertices for this leaf.
    std::vector<Vertex> leafVertices;
    std::vector<Vec3<int32_t>> refinedVerticesBlock;

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
      blockWidth = direction.max();
      
      // vertex to list of points 
      Vec3<int32_t> foundvoxel = segment.startpos;
      for (int k = 0; k <= 2; k++) {
        if (direction[k])
          foundvoxel[k] += segment.vertex == (blockWidth >> bitDropped) - 1 ? blockWidth - 1 : segment.vertex << bitDropped;
        if (segment.startpos[k] - leaves[i].pos[k] > 0) // back to B-1 if B 
          foundvoxel[k]--;
      }
      if (boundaryinsidecheck(foundvoxel, poistionClipValue))
        refinedVerticesBlock.push_back(foundvoxel);

      // Get 3D position of point of intersection.      
      Vec3<int32_t> point = (segment.startpos - leaves[i].pos) << kTrisoupFpBits;
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

      // Add vertex to list of points.     
      leafVertices.push_back({ point, 0, 0 });
    }

    // Skip leaves that have fewer than 3 vertices.
    if (leafVertices.size() < 3) {
      refinedVertices.insert(refinedVertices.end(), refinedVerticesBlock.begin(), refinedVerticesBlock.end());
      continue;
    }

    // compute centroid 
    int triCount = (int)leafVertices.size();
    Vec3<int32_t> blockCentroid = 0;
    for (int j = 0; j < triCount; j++) {
      blockCentroid += leafVertices[j].pos;
    }
    blockCentroid /= triCount; 


    // order vertices along a dominant axis only if more than three (otherwise only one triangle, whatever...)
    int dominantAxis = findDominantAxis(leafVertices, blockWidth, blockCentroid);    


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
      int bound = (int(blockWidth) - 1) << kTrisoupFpBits;
      int m = 1;
      for (; m < blockWidth; m++) {
        Vec3<int32_t> temp = blockCentroid + m * normalV;
        if (temp[0]<0 || temp[1]<0 || temp[2]<0 || temp[0]>bound || temp[1]>bound || temp[2]> bound)
          break;
      }
      int highBound = (m - 1) + halfDropped2 >> bitDropped2;

      m = 1;
      for (; m < blockWidth; m++) {
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
          auto point = (pointCloud[p] - leaves[i].pos) << kTrisoupFpBits;
          
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
    } // end refinement of the centroid 


    // Divide vertices into triangles around centroid
    // and upsample each triangle by an upsamplingFactor.    
    Vec3<int32_t> v2 = triCount == 3 ? leafVertices[2].pos : blockCentroid;
    Vec3<int32_t> v1 = leafVertices[0].pos;

    Vec3<int32_t> posNode = leaves[i].pos << kTrisoupFpBits;       

    if (triCount > 3) {
      Vec3<int32_t> foundvoxel = (posNode + blockCentroid + truncateValue) >> kTrisoupFpBits;
      if (boundaryinsidecheck(foundvoxel, poistionClipValue)) 
        refinedVerticesBlock.push_back(foundvoxel);
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
        maxRange[k] = std::min(33, std::max(std::max(v0[k], v1[k]), v2[k]) + truncateValue >> kTrisoupFpBits);
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
        rayTracingAlongdirection(refinedVerticesBlock, direction, samplingValue, posNode, minRange, maxRange, edge1, edge2, v0, poistionClipValue, haloFlag, fineRayflag);   
      }    

    }  // end loop on triangles

    std::sort(refinedVerticesBlock.begin(), refinedVerticesBlock.end());
    refinedVerticesBlock.erase(std::unique(refinedVerticesBlock.begin(), refinedVerticesBlock.end()), refinedVerticesBlock.end());
    refinedVertices.insert(refinedVertices.end(), refinedVerticesBlock.begin(), refinedVerticesBlock.end());

  }// end loop on leaves

 // remove points present twice or more 
  std::sort(refinedVertices.begin(), refinedVertices.end());
  refinedVertices.erase( std::unique(refinedVertices.begin(), refinedVertices.end()), refinedVertices.end());
  
  // Move list of points to pointCloud.
  pointCloud.resize(refinedVertices.size());
  for (int i = 0; i < refinedVertices.size(); i++) {
    pointCloud[i] = refinedVertices[i];
  }
}


// ---------------------------------------------------------------------------
void decodeTrisoupVertices(
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  int bitDropped,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyDecoder& arithmeticDecoder) {

  int iV = 0;
  std::vector<int> correspondanceSegment2V;

  AdaptiveBitModel ctxTempV2[144];

  CtxModelTrisoup ctxTriSoup;
  CtxMapTriSoup MapOBUFTriSoup[3];
  MapOBUFTriSoup[0].reset(10, 7); // flag
  MapOBUFTriSoup[1].reset(10, 6);  // first bit position
  MapOBUFTriSoup[2].reset(10, 6 + 1);  // second bit position
  for (int i = 0; i <= gbh.num_unique_segments_minus1; i++) {
    int ctxE = (!!(neighbNodes[i] & 1)) + (!!(neighbNodes[i] & 2)) + (!!(neighbNodes[i] & 4)) + (!!(neighbNodes[i] & 8)) - 1; // at least one node is occupied 
    int ctx0 = (!!(neighbNodes[i] & 16)) + (!!(neighbNodes[i] & 32)) + (!!(neighbNodes[i] & 64)) + (!!(neighbNodes[i] & 128));
    int ctx1 = (!!(neighbNodes[i] & 256)) + (!!(neighbNodes[i] & 512)) + (!!(neighbNodes[i] & 1024)) + (!!(neighbNodes[i] & 2048));

    int direction = neighbNodes[i] >> 13;

    int beforeCtx = 0;
    int Vbefore = 0;
    int nclose = 0;
    int nfar = 0;

    // compute data relative to preceding (before) vertex along edge direction (presence and position) 
    // vertices are on gbh.trisoupNodeSizeLog2(gps) - bitDropped bits
    int nbitsVertices = gbh.trisoupNodeSizeLog2(gps) - bitDropped;

    if (indexBefore[i] != -1 && gbh.trisoupNodeSizeLog2(gps) <= 4) {
      beforeCtx = segind[indexBefore[i]];
      if (correspondanceSegment2V[indexBefore[i]] != -1) {
        Vbefore = 1 + vertices[correspondanceSegment2V[indexBefore[i]]] >> std::max(0, nbitsVertices - 2);

        int v2bits = 3 - (vertices[correspondanceSegment2V[indexBefore[i]]] >> nbitsVertices - 2);
        if (v2bits <= 0)
          nclose++;
        if (v2bits >= 3)
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
        int vertexPos = vertices[idxVertex] >> gbh.trisoupNodeSizeLog2(gps) - bitDropped - 2; // on 2 bits

        int orientation = perpVS[k] >> 30;
        if (orientation) {// if toward then reverse to away 
          vertexPos = 3 - vertexPos; // 0 is closest, 3 is farthest         
        }

        if (vertexPos <= 0)
          nclose++;
        if (vertexPos >= 3)
          nfar++;
      }
    }
    int perpStartCtx = (4 - direction) - occupPerp + 1 - beforeCtx;   // in [0,5]
    int nmiddle = occupPerp + beforeCtx - nfar - nclose;
    bool flagTouch = nclose > 0;

    // encode flag vertex
    bool c;


    int ctxMap1 = (Vbefore * 4 + ctxE) * 4 + std::min(nclose, 3);
    bool full01 = ((ctx0 == 4) || (ctx1 == 4));
    int ctxMap2 = full01 << 8;
    if (!full01) { // none is full 
      ctxMap2 |= (neighbNodes[i] >> 4) & 255;
    }
    else { // one is full, very rarely both 
      ctxMap2 |= (ctx0 == 4) << 7;
      if (ctx0 == 4)
        ctxMap2 |= ((neighbNodes[i] >> 8) & 15) << 3;
      else
        ctxMap2 |= ((neighbNodes[i] >> 4) & 15) << 3;
      // three bits to spare 
      ctxMap2 |= perpStartCtx;
    }
    c = arithmeticDecoder.decode(ctxTriSoup[MapOBUFTriSoup[0].get(ctxMap2, ctxMap1)]);
    MapOBUFTriSoup[0].evolve(c, ctxMap2, ctxMap1);


    segind.push_back(c);
    correspondanceSegment2V.push_back(-1);

    // encode position vertex 
    if (c) {
      correspondanceSegment2V.back() = iV;

      uint8_t v = 0;
      int ctxFullNbounds = (4 * (std::max(1, ctx0) - 1) + (std::max(1, ctx1) - 1)) * 2 + (ctxE == 3);
      int b = gbh.trisoupNodeSizeLog2(gps) - bitDropped - 1;

      // first bit
      ctxMap1 = ctxFullNbounds * 2 + flagTouch;
      ctxMap2 = perpStartCtx << 7;
      ctxMap2 |= beforeCtx << 6;
      ctxMap2 |= std::min(3, nfar) << 4;
      ctxMap2 |= std::min(3, nclose) << 2;
      ctxMap2 |= std::min(3, nmiddle);
      int bit = arithmeticDecoder.decode(ctxTriSoup[MapOBUFTriSoup[1].get(ctxMap2, ctxMap1)]);
      MapOBUFTriSoup[1].evolve(bit, ctxMap2, ctxMap1);
      v = (v << 1) | bit;
      b--;

      // second bit       
      bit = arithmeticDecoder.decode(ctxTriSoup[MapOBUFTriSoup[2].get(ctxMap2, (ctxMap1 << 1) + v)]);
      MapOBUFTriSoup[2].evolve(bit, ctxMap2, (ctxMap1 << 1) + v);
      v = (v << 1) | bit;
      b--;

      // third bit
      if (b >= 0) {
        int ctxFullNboundsReduced1 = (6 * (ctx0 >> 1) + perpStartCtx) * 2 + (ctxE == 3);
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
  uint32_t blockWidth, 
  Vec3<int32_t> blockCentroid ) {

  int dominantAxis = 0;
  int triCount = leafVertices.size();
  if (triCount > 3) {
    Vertex vertex;
    int32_t Width = blockWidth << kTrisoupFpBits;

    const int sIdx1[3] = { 2,2,1 };
    const int sIdx2[3] = { 1,0,0 };

    int maxNormTri = 0;
    for (int axis = 0; axis <= 2; axis++) {
      // order along axis
      for (int j = 0; j < triCount; j++) {
        Vec3<int32_t> s = leafVertices[j].pos + kTrisoupFpHalf; // back to [0,B]^3 for ordering 
        leafVertices[j].theta = trisoupVertexArc(s[sIdx1[axis]], s[sIdx2[axis]], Width);
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
      leafVertices[j].theta = trisoupVertexArc(s[sIdx1[dominantAxis]], s[sIdx2[dominantAxis]], Width);
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
  Vec3<int32_t> posNode,
  int minRange[3],
  int maxRange[3],
  Vec3<int32_t> edge1,
  Vec3<int32_t> edge2,
  Vec3<int32_t> v0,
  int poistionClipValue,
  bool haloFlag,
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
  const int haloTriangle = haloFlag ? 32 : 0;
  for (int32_t g1 = startposG1; g1 <= endposG1; g1 += samplingValue) {
    rayOrigin[g1pos[direction]] = g1 << kTrisoupFpBits;


    for (int32_t g2 = startposG2; g2 <= endposG2; g2 += samplingValue) {
      rayOrigin[g2pos[direction]] = g2 << kTrisoupFpBits;

      // middle ray at integer position 
      Vec3<int32_t>  intersection = rayOrigin;
      bool foundIntersection = rayIntersectsTriangle(rayOrigin, v0, edge1, edge2, h, a, intersection, direction, haloTriangle);
      if (foundIntersection) {
        Vec3<int32_t> foundvoxel = (posNode + intersection + truncateValue) >> kTrisoupFpBits;
        if (boundaryinsidecheck(foundvoxel, poistionClipValue)) {
          refinedVerticesBlock.push_back(foundvoxel);
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
          if (rayIntersectsTriangle(rayOrigin2, v0, edge1, edge2, h, a, intersection, direction, haloTriangle)) {
            Vec3<int32_t> foundvoxel = (posNode + intersection + truncateValue) >> kTrisoupFpBits;
            if (boundaryinsidecheck(foundvoxel, poistionClipValue)) {
              refinedVerticesBlock.push_back(foundvoxel);
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
