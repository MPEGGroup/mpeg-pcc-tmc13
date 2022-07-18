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
// ----------- dynamic OBUF for TriSoup is here -----------------------------
//============================================================================
  class CtxMapTriSoup {
  public:
    int S1 = 0; // 16;
    int S2 = 0;// 128 * 2 * 8;


    uint8_t* CtxIdxMap; // S1*S2  
    uint8_t* kDown; // S2 
    uint8_t* Nseen; // S2

                    //  allocate and reset CtxIdxMap to 127
    void reset(int userBitS1, int userBitS2);

    //  deallocate CtxIdxMap
    void clear();

    //  return ctxIdx according to bit
    uint8_t get(int i, int j);

    //  update *ctxIdx according to bit
    void evolve(bool bit, int i, int j);
    //  update *ctxIdx according to bit
    uint8_t getEvolve(bool bit, int i, int j);

  private:
    int maxTreeDepth = 3;
    int minkTree = 0;
    int LUTth[16];

    //  update kDown
    void  decreaseKdown(int iP, int j, int iTree, int kDown0, int  kTree);
    int idx(int i, int j);
  };

  inline void
    CtxMapTriSoup::reset(int userBitS1, int userBitS2)
  {
    S1 = 1 << userBitS1;
    S2 = 1 << userBitS2;
    maxTreeDepth = userBitS1 - 0;

    minkTree = userBitS1 - maxTreeDepth;

    for (int k =0; k <16;k++)
      LUTth[k] = std::min(32, 8 << std::max(0, minkTree - k - 1))/2;

    kDown = new uint8_t[(1 << maxTreeDepth) * S2];
    std::memset(kDown, userBitS1, sizeof * kDown * (1 << maxTreeDepth) * S2);
    Nseen = new uint8_t[(1 << maxTreeDepth) * S2];
    std::memset(Nseen, 0, sizeof * Nseen * (1 << maxTreeDepth) * S2);
    CtxIdxMap = new uint8_t[S1 * S2];
    std::memset(CtxIdxMap, 127, sizeof * CtxIdxMap * S1 * S2);
  }

  inline void
    CtxMapTriSoup::clear()
  {
    if (!S1 || !S2)
      return;

    delete[] kDown ;
    delete[] Nseen;
    delete[] CtxIdxMap;
  }

  inline uint8_t
    CtxMapTriSoup::get(int i, int j)
  {
    int kDown0 = kDown[idx(i >> minkTree, j)];
    int iP = (i >> kDown0) << kDown0;
    return  CtxIdxMap[idx(iP, j)];
  }

  inline void
    CtxMapTriSoup::evolve(bool bit, int i, int j)
  {  
    int iTree = i >> minkTree;
    int kDown0 = kDown[idx(iTree, j)];
    int iP = (i >> kDown0) << kDown0;
    uint8_t* ctxIdx = &(CtxIdxMap[idx(iP, j)]);

    if (bit)
      *ctxIdx += kCtxMapTriSoupDelta[(255 - *ctxIdx) >> 4];
    else
      *ctxIdx -= kCtxMapTriSoupDelta[*ctxIdx >> 4];  

    if (kDown0) { 
      int kTree = std::max(0, kDown0 - minkTree);
      iTree = (iTree >> kTree) << kTree;
      if (++Nseen[idx(iTree, j)] >= LUTth[kDown0])     // if more than th stats per pack    
        decreaseKdown(iP, j, iTree, kDown0, kTree);
    }
  }

  inline uint8_t
    CtxMapTriSoup::getEvolve(bool bit, int i, int j)
  {  
    int iTree = i >> minkTree;
    int kDown0 = kDown[idx(iTree, j)];
    int iP = (i >> kDown0) << kDown0;

    uint8_t* ctxIdx = &(CtxIdxMap[idx(iP, j)]);
    uint8_t out = *ctxIdx;

    if (bit)
      *ctxIdx += kCtxMapTriSoupDelta[(255 - *ctxIdx) >> 4];
    else
      *ctxIdx -= kCtxMapTriSoupDelta[*ctxIdx >> 4];

    if (kDown0) { 
      int kTree = std::max(0, kDown0 - minkTree);
      iTree = (iTree >> kTree) << kTree;
      if (++Nseen[idx(iTree, j)] >= LUTth[kDown0])     // if more than th stats per pack     
        decreaseKdown(iP, j, iTree, kDown0, kTree);
    }

    return out;
  }

  inline void
    CtxMapTriSoup::decreaseKdown(int iP, int j, int iTree, int kDown0, int  kTree)
  {
    int idxTree = idx(iTree, j);
    Nseen[idxTree] = 0; // setting other Nseen unneeded because initialized to 0
    if (kTree) { // binary-tree based dynamic OBUF  
      int iEnd = S2 << kTree;
      for (int ii = 0; ii < iEnd; ii += S2)
        kDown[idxTree + ii]--;

      auto *p = &CtxIdxMap[idx(iP, j)];
      p[S2 << kDown0 - 1] = *p;
    }
    else {// simple dynamic OBUF on a subblock of states attached to a leaf of binary-tree
      kDown[idxTree]--;

      int S1Simple = S1 >> maxTreeDepth;
      int SS1 = (S1Simple >> kDown0) ;    
      int stride = S2 << kDown0;

      auto *p = &CtxIdxMap[idx(iTree << minkTree, j)];
      for (int iL = 0; iL < SS1; iL++, p += stride)  // pack       
        p[stride >> 1] = *p;
    }
  }

  inline int
    CtxMapTriSoup::idx(int i, int j)
  {  
    return i * S2 + j;
  }

  //---------------------------------------------------------------------------

  struct CtxModelTrisoup {
    static const int kCtxFactorShift = 4;
    AdaptiveBitModelFast contexts[256>>kCtxFactorShift];  

    AdaptiveBitModelFast& operator[](int idx)
    {
      return contexts[idx >> kCtxFactorShift];
    }
  };



// ----------- end dynamic OBUF -------
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
    bool operator()(Vertex v1, Vertex v2)
    {
      if (v1.theta > v2.theta)
        return true;  // sort in decreasing order of theta
      if (v1.theta == v2.theta && v1.tiebreaker < v2.tiebreaker)
        return true;
      return false;
    }

    bool operator==(Vertex v1)
    {
      return this->pos == v1.pos;
    }

    bool operator<(Vertex v1)
    {
      return this->pos < v1.pos;
    }

  } ;


 //============================================================================
void determineTrisoupVertices(
  const ringbuf<PCCOctree3Node>& leaves,
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,  
  const PCCPointSet3& pointCloud,
  const int defaultBlockWidth,
  const int bitDropped,
  int distanceSearchEncoder);

void determineTrisoupNeighbours(
  const ringbuf<PCCOctree3Node>& leaves, 
  std::vector<uint16_t>& neighbNodes, 
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  const int defaultBlockWidth);

void encodeTrisoupVertices(  
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  int bitDropped,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder);

void decodeTrisoupVertices(  
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  int bitDropped,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyDecoder& arithmeticDecoder);


void encodeTrisoupCentroidResidue(
  std::vector<CentroidDrift>& drifts,
  pcc::EntropyEncoder* arithmeticEncoder);

void decodeTrisoupCommon(
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
  pcc::EntropyDecoder* arithmeticDecoder);

int findDominantAxis(
  std::vector<Vertex>& leafVertices,
  uint32_t blockWidth, 
  Vec3<int32_t> blockCentroid ); 

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

struct TrisoupSegmentNeighbours {
  Vec3<int32_t> startpos;  // start point of edge segment
  Vec3<int32_t> endpos;    // end point of edge segment

  int index;        // index of segment, to reorder after sorting 
  uint16_t neighboursMask;   
};



//----------------------------------------------------------------------------
// comparison for sorting

bool operator<(const TrisoupSegment& s1, const TrisoupSegment& s2);

bool operator<(const TrisoupSegmentNeighbours& s1, const TrisoupSegmentNeighbours& s2);




//============================================================================

}  // namespace pcc
