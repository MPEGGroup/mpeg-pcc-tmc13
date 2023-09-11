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

#include <cstdint>
#include <vector>
#include <cstring>

#include "PCCPointSet.h"
#include "geometry_octree.h"
#include "ringbuf.h"

namespace pcc {

//============================================================================

// The number of fractional bits used in trisoup triangle voxelisation
const int kTrisoupFpBits = 8;

// The value 1 in fixed-point representation
const int kTrisoupFpOne = 1 << (kTrisoupFpBits);
const int kTrisoupFpHalf = 1 << (kTrisoupFpBits - 1);
const int truncateValue = kTrisoupFpHalf;

//============================================================================

  struct CentroidDrift {
    int driftQ;
    int lowBound;
    int highBound;
    int ctxMinMax;
    int lowBoundSurface;
    int  highBoundSurface;

  };


  // Representation for a vertex in preparation for sorting.
  struct Vertex {
    Vec3<int32_t> pos;  // position of vertex
    int32_t theta;      // angle of vertex when projected along dominant axis
    int32_t tiebreaker;  // coordinate of vertex along dominant axis
    bool operator()(const Vertex& v1, const Vertex& v2) const
    {
      if (v1.theta > v2.theta)
        return true;  // sort in decreasing order of theta
      if (v1.theta == v2.theta && v1.tiebreaker < v2.tiebreaker)
        return true;
      return false;
    }

    bool operator==(const Vertex& v1) const
    {
      return this->pos == v1.pos;
    }

    bool operator<(const Vertex& v1) const
    {
      return this->pos < v1.pos;
    }

  } ;

//============================================================================

  struct TrisoupNodeEdgeVertex {
    int dominantAxis;
    // ( x - 0.5 ) << kTrisoupFpBits ( [-0.5, W-0.5] x256 )
    std::vector<Vertex> vertices;
  };


  struct TrisoupCentroidVertex {
    bool valid;  // this represents centroid existence
    Vec3<int32_t> pos;
    int32_t drift;
    bool boundaryInside;  // true if pos is inside of node boundary
  };


  struct TrisoupNodeFaceVertex {
    std::vector<Vertex> vertices;
    std::vector<int> formerEdgeVertexIdx;
  };


  struct duplicateNodes {
    Vec3<int32_t> pos;
    int idx;  // duplicateListIdx = (origIdx<<3)+dupNodeIdx
    bool operator<( const duplicateNodes& n ) const {
      return this->pos < n.pos;
    }
  };


  struct node6nei {
    Vec3<int32_t> pos;
    // indices of octree-nodes
    int idx[7]; // {-z,z,-y,y,-x,x,O}

    node6nei(){ this->clear(); }
    void clear(void) {
      pos = {-1,-1,-1};
      for( int i=0; i<7; i++){
        idx[i]=-1;
      }
    }

    bool operator<( const node6nei& n ) const {
      return this->idx[6] < n.idx[6];
    }
  };


  struct TrisoupFace {
    bool connect;

    TrisoupFace( bool cn ) : connect(cn) {}

    TrisoupFace(){ this->clear(); }

    void clear(void) {
      connect = false;
    }
  };

  struct TrisoupCentroidContext {
    int lowBound;
    int highBound;
    int ctxMinMax;
    int lowBoundSurface;
    int highBoundSurface;
  };

//============================================================================

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
  float estimatedSampling = 0.f,
  bool nodeUniqueDSE = false);


void determineTrisoupVertices(
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
  bool nodeUniqueDSE);

void determineTrisoupNeighbours(
  const ringbuf<PCCOctree3Node>& leaves, 
  std::vector<uint16_t>& neighbNodes, 
  std::vector<std::array<int, 18>>& edgePattern,
  const int defaultBlockWidth);

bool determineNormVandCentroidContexts
(
 const Vec3<int32_t>& nodeWidth,
 const TrisoupNodeEdgeVertex& eVert,
 const int bitDropped,
 Vec3<int32_t>& gravityCenter,
 Vec3<int32_t>& normalV,
 TrisoupCentroidContext& cctx
);

void determineTrisoupNodeNeighbours
(
  const ringbuf<PCCOctree3Node>& leaves,
  std::vector<node6nei>& nodes6nei,
  const int defaultBlockWidth
  );

void determineTrisoupCentroids
(
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
  std::vector<Vec3<int32_t>>& normVs
 );

void determineTrisoupFaceVertices
(
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
  std::vector<Vec3<int32_t>> normVs,
  std::vector<TrisoupFace>& limited_faces,
  std::vector<TrisoupFace>& faces
  );

void encodeTrisoupFaceList(
  std::vector<TrisoupFace>& faceList,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder );

void decodeTrisoupFaceList(
  const ringbuf<PCCOctree3Node>& leaves,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int defaultBlockWidth,
  const int distanceSearchEncoder,
  const std::vector<node6nei> nodes6nei,
  const std::vector<TrisoupNodeEdgeVertex> eVerts,
  const std::vector<Vec3<int32_t>> gravityCenter,
  const std::vector<TrisoupCentroidVertex> cVerts,
  std::vector<TrisoupNodeFaceVertex>& fVerts,
  std::vector<Vec3<int32_t>> normVs,
  std::vector<TrisoupFace>& faces,
  pcc::EntropyDecoder& arithmeticDecoder );

void decodeTrisoupCentroids(
  const ringbuf<PCCOctree3Node>& leaves,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int defaultBlockWidth,
  const int bitDropped,
  const std::vector<TrisoupNodeEdgeVertex>& eVerts,
  std::vector<Vec3<int32_t>>& gravityCenter,
  std::vector<TrisoupCentroidVertex>& cVerts,
  std::vector<Vec3<int32_t>>& normVs,
  pcc::EntropyDecoder* arithmeticDecoder );

bool nodeBoundaryInsideCheck( Vec3<int32_t> bw, Vec3<int32_t> pt );

void findTrisoupFaceVertex
( const int nodeIdx,
  const int neiOrderIdx,
  const node6nei& nodes6nei,
  const std::vector<TrisoupCentroidVertex>& cVerts,
  const Vec3<int32_t>& nodew,
  Vertex *fVert );

void encodeTrisoupVertices(
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  int bitDropped,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder);

void decodeTrisoupVertices(
  const ringbuf<PCCOctree3Node>& leaves,
  const int defaultBlockWidth,
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  int bitDropped,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  std::vector<TrisoupNodeEdgeVertex>& eVerts,
  pcc::EntropyDecoder& arithmeticDecoder);

void decodeTrisoupVerticesSub(
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<std::array<int, 18>>& edgePattern,
  int bitDropped,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyDecoder& arithmeticDecoder);


void encodeTrisoupCentroidResidue(
  std::vector<CentroidDrift>& drifts,
  pcc::EntropyEncoder* arithmeticEncoder);

void decodeTrisoupCommon(
  const ringbuf<PCCOctree3Node>& leaves,
  const std::vector<node6nei> nodes6nei,
  const std::vector<TrisoupNodeEdgeVertex>& eVerts,
  const std::vector<TrisoupCentroidVertex>& cVerts,
  const std::vector<Vec3<int32_t>>& gravityCenter,
  std::vector<Vec3<int32_t>> normVs,
  std::vector<TrisoupFace>& faces,
  std::vector<TrisoupNodeFaceVertex>& fVerts,
  PCCPointSet3& recPointCloud,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int defaultBlockWidth,
  int positionClipValue,
  uint32_t samplingValue,
  const int bitDropped,
  const bool isCentroidDriftActivated,
  const bool isFaceVertexActivated,
  bool haloFlag,
  bool adaptiveHaloFlag,
  bool fineRayflag );

int findDominantAxis(
  std::vector<Vertex>& leafVertices,
  Vec3<uint32_t> blockWidth,
  Vec3<int32_t> blockCentroid);

int countTrisoupEdgeVerticesOnFace
( const TrisoupNodeEdgeVertex& eVerts, Vec3<int32_t>& nodeW, int axis );

void determineTrisoupEdgeBoundaryLine
( int i,
  const TrisoupNodeEdgeVertex& eVerts,
  const TrisoupCentroidVertex& cVerts,
  const Vec3<int32_t>& gravityCenter,
  Vec3<int32_t>& nodeW, int axis, Vertex& fvert, int *eIdx );

bool determineTrisoupDirectionOfCentroidsAndFvert
( const TrisoupNodeEdgeVertex& eVerts,
  const std::vector<TrisoupCentroidVertex>& cVerts,
  const std::vector<Vec3<int32_t>>& gravityCenter,
  Vec3<int32_t>& nodew, int i, int nei,
  const node6nei& nodes6nei, int neiNodeIdx,
  int w, int e0, int e1, Vertex *fVert );


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
  bool fineRayflag);

//============================================================================

struct TrisoupSegment {
  Vec3<int32_t> startpos;  // start point of edge segment
  Vec3<int32_t> endpos;    // end point of edge segment

  int index;        // index of segment, to reorder after sorting
  int uniqueIndex;  // index of uniqueSegment
  int vertex;       // distance along segment for intersection (else -1)
};

struct TrisoupSegmentEnc : public TrisoupSegment {
  TrisoupSegmentEnc(
    const Vec3<int32_t>& startpos,
    const Vec3<int32_t>& endpos,
    int index,
    int uniqueIndex,
    int vertex,
    int count,
    int distanceSum,
    int count2,
    int distanceSum2)
    : TrisoupSegment{startpos, endpos, index, uniqueIndex, vertex}
    , count(count)
    , distanceSum(distanceSum)
    , count2(count2)
    , distanceSum2(distanceSum2)
  {}

  int count;        // count of voxels adjacent to this segment
  int distanceSum;  // sum of distances (along segment) of adjacent voxels  
  int count2;        // count of voxels adjacent to this segment
  int distanceSum2;  // sum of distances (along segment) of adjacent voxels
};

struct SegmentHash {
  std::size_t operator()(const TrisoupSegmentEnc& s) const
  {
    std::size_t h1 = std::hash<int>()(s.startpos[0]);
    std::size_t h2 = std::hash<int>()(s.startpos[1]);
    std::size_t h3 = std::hash<int>()(s.startpos[2]);
    std::size_t h4 = std::hash<int>()(s.endpos[0]);
    std::size_t h5 = std::hash<int>()(s.endpos[1]);
    std::size_t h6 = std::hash<int>()(s.endpos[2]);
    return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4) ^ (h6 << 5); // Combine the six hash values
  }
};

struct TrisoupSegmentNeighbours {
  Vec3<int32_t> startpos;  // start point of edge segment
  Vec3<int32_t> endpos;    // end point of edge segment

  int index;        // index of segment, to reorder after sorting 
  uint16_t neighboursMask;   
};


enum{
  POS_000 = 0,
  POS_W00 = 1,
  POS_0W0 = 2,
  POS_WW0 = 3,
  POS_00W = 4,
  POS_W0W = 5,
  POS_0WW = 6,
  POS_WWW = 7
};

void nonCubicNode
(
 const GeometryParameterSet& gps,
 const GeometryBrickHeader& gbh,
 const Vec3<int32_t>& leafpos,
 const int32_t blockWidth,
 const Box3<int32_t>& bbox,
 Vec3<int32_t>& newp,
 Vec3<int32_t>& neww,
 Vec3<int32_t>* corner
 );


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


//----------------------------------------------------------------------------
// comparison for sorting

bool operator<(const TrisoupSegment& s1, const TrisoupSegment& s2);

bool operator<(const TrisoupSegmentNeighbours& s1, const TrisoupSegmentNeighbours& s2);

bool operator==(const TrisoupSegment& s1, const TrisoupSegment& s2);

//============================================================================

}  // namespace pcc
