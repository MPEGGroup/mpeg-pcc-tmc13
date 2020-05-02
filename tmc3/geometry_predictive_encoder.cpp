/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2020, ISO/IEC
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

#include "geometry_predictive.h"
#include "geometry.h"
#include "pointset_processing.h"

#include "PCCMisc.h"

#include "nanoflann.hpp"

namespace pcc {

//============================================================================

namespace {
  struct NanoflannCloud {
    std::vector<Vec3<int32_t>> pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    int32_t kdtree_get_pt(const size_t idx, const size_t dim) const
    {
      return pts[idx][dim];
    }

    template<class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
      return false;
    }
  };
}  // namespace

//============================================================================

namespace {
  float estimate(bool bit, AdaptiveBitModel model)
  {
    return -log2(dirac::approxSymbolProbability(bit, model) / 128.);
  }
}  // namespace

//============================================================================

class PredGeomEncoder : public PredGeomCodec {
public:
  PredGeomEncoder(const PredGeomEncoder&) = delete;
  PredGeomEncoder& operator=(const PredGeomEncoder&) = delete;

  PredGeomEncoder(const GeometryParameterSet&, EntropyEncoder* aec);

  void encode(
    const Vec3<int32_t>* cloud,
    const GNode* nodes,
    int numNodes,
    int* codedOrder);

  int encodeTree(
    const Vec3<int32_t>* cloud,
    const GNode* nodes,
    int numNodes,
    int rootIdx,
    int* codedOrder);

  void encodeNumDuplicatePoints(int numDupPoints);
  void encodeNumChildren(int numChildren);
  void encodePredMode(GPredicter::Mode mode);
  void encodeResidual(const Vec3<int32_t>& residual, GPredicter::Mode mode);

  float estimateBits(GPredicter::Mode mode, const Vec3<int32_t>& residual);

private:
  EntropyEncoder* _aec;
  std::vector<int32_t> _stack;
  bool _geom_unique_points_flag;
};

//============================================================================

PredGeomEncoder::PredGeomEncoder(
  const GeometryParameterSet& gps, EntropyEncoder* aec)
  : _aec(aec), _geom_unique_points_flag(gps.geom_unique_points_flag)
{
  _stack.reserve(1024);
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeNumDuplicatePoints(int numDupPoints)
{
  _aec->encode(numDupPoints > 0, _ctxNumDupPointsGt0);
  if (numDupPoints)
    _aec->encodeExpGolomb(numDupPoints - 1, 0, _ctxBypass, _ctxNumDupPoints);
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeNumChildren(int numChildren)
{
  _aec->encode(numChildren & 1, _ctxNumChildren[0]);
  _aec->encode((numChildren >> 1) & 1, _ctxNumChildren[1 + (numChildren & 1)]);
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodePredMode(GPredicter::Mode mode)
{
  int iMode = int(mode);
  _aec->encode(iMode & 1, _ctxPredMode[0]);
  _aec->encode((iMode >> 1) & 1, _ctxPredMode[1 + (iMode & 1)]);
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeResidual(
  const Vec3<int32_t>& residual, GPredicter::Mode mode)
{
  int iMode = int(mode);
  for (int k = 0, ctxIdx = 0; k < 3; k++) {
    const auto res = residual[k];
    const bool isZero = res == 0;
    _aec->encode(isZero, _ctxIsZero[k]);
    if (isZero)
      continue;

    if (iMode > 0) {
      _aec->encode(res > 0, _ctxSign[k]);
    }

    int32_t value = abs(res) - 1;
    int32_t numBits = 1 + ilog2(uint32_t(value));

    AdaptiveBitModel* ctxs = _ctxNumBits[ctxIdx][k];
    _aec->encode(numBits & 1, ctxs[0]);
    _aec->encode((numBits >> 1) & 1, ctxs[1 + (numBits & 1)]);
    _aec->encode((numBits >> 2) & 1, ctxs[3 + (numBits & 3)]);
    _aec->encode((numBits >> 3) & 1, ctxs[7 + (numBits & 7)]);
    _aec->encode((numBits >> 4) & 1, ctxs[15 + (numBits & 15)]);

    if (!k)
      ctxIdx = (numBits + 1) >> 1;

    --numBits;
    for (int32_t i = 0; i < numBits; ++i)
      _aec->encode((value >> i) & 1, _ctxBypass);
  }
}

//----------------------------------------------------------------------------

float
PredGeomEncoder::estimateBits(
  GPredicter::Mode mode, const Vec3<int32_t>& residual)
{
  int iMode = int(mode);
  float bits = 0.;
  bits += estimate(iMode & 1, _ctxPredMode[0]);
  bits += estimate((iMode >> 1) & 1, _ctxPredMode[1 + (iMode & 1)]);

  for (int k = 0, ctxIdx = 0; k < 3; k++) {
    const auto res = residual[k];
    const bool isZero = res == 0;
    bits += estimate(isZero, _ctxIsZero[k]);
    if (isZero)
      continue;

    if (iMode > 0) {
      bits += estimate(res > 0, _ctxSign[k]);
    }

    int32_t value = abs(res) - 1;
    int32_t numBits = 1 + ilog2(uint32_t(value));

    AdaptiveBitModel* ctxs = _ctxNumBits[ctxIdx][k];
    bits += estimate(numBits & 1, ctxs[0]);
    bits += estimate((numBits >> 1) & 1, ctxs[1 + (numBits & 1)]);
    bits += estimate((numBits >> 2) & 1, ctxs[3 + (numBits & 3)]);
    bits += estimate((numBits >> 3) & 1, ctxs[7 + (numBits & 7)]);
    bits += estimate((numBits >> 4) & 1, ctxs[15 + (numBits & 15)]);

    if (!k)
      ctxIdx = (numBits + 1) >> 1;

    bits += numBits - 1;
  }

  return bits;
}

//----------------------------------------------------------------------------

int
PredGeomEncoder::encodeTree(
  const Vec3<int32_t>* cloud,
  const GNode* nodes,
  int numNodes,
  int rootIdx,
  int* codedOrder)
{
  int processedNodes = 0;

  _stack.push_back(rootIdx);
  while (!_stack.empty()) {
    const auto nodeIdx = _stack.back();
    _stack.pop_back();

    const auto& node = nodes[nodeIdx];
    const auto point = cloud[nodeIdx];

    struct {
      float bits;
      GPredicter::Mode mode;
      Vec3<int32_t> residual;
    } best;

    // mode decision to pick best prediction from available set
    for (int iMode = 0; iMode < 4; iMode++) {
      GPredicter::Mode mode = GPredicter::Mode(iMode);
      GPredicter predicter = makePredicter(
        nodeIdx, mode, [=](int idx) { return nodes[idx].parent; });

      if (!predicter.isValid(mode))
        continue;

      auto pred = predicter.predict(&cloud[0], mode);
      auto residual = point - pred;
      auto bits = estimateBits(mode, residual);

      if (iMode == 0 || bits < best.bits) {
        best.residual = residual;
        best.mode = mode;
        best.bits = bits;
      }
    }

    assert(node.childrenCount <= GNode::MaxChildrenCount);
    if (!_geom_unique_points_flag)
      encodeNumDuplicatePoints(node.numDups);
    encodeNumChildren(node.childrenCount);
    encodePredMode(best.mode);
    encodeResidual(best.residual, best.mode);

    // NB: the coded order of duplicate points assumes that the duplicates
    // are consecutive -- in order that the correct attributes are coded.
    codedOrder[processedNodes++] = nodeIdx;
    for (int i = 1; i <= node.numDups; i++)
      codedOrder[processedNodes++] = nodeIdx + i;

    for (int i = 0; i < node.childrenCount; i++) {
      _stack.push_back(node.children[i]);
    }
  }

  return processedNodes;
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encode(
  const Vec3<int32_t>* cloud,
  const GNode* nodes,
  int numNodes,
  int32_t* codedOrder)
{
  int32_t processedNodes = 0;
  for (int32_t rootIdx = 0; rootIdx < numNodes; rootIdx++) {
    // find the root node(s)
    if (nodes[rootIdx].parent >= 0)
      continue;

    int numSubtreeNodes =
      encodeTree(cloud, nodes, numNodes, rootIdx, codedOrder + processedNodes);
    processedNodes += numSubtreeNodes;
  }
  assert(processedNodes == numNodes);
}

//============================================================================

std::vector<GNode>
generateGeomPredictionTree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const Vec3<int32_t>* begin,
  const Vec3<int32_t>* end)
{
  const int32_t pointCount = std::distance(begin, end);

  // Build the prediction tree, roughly as follows:
  //  - For each point, find the node with the nearest prediction
  //    with empty children slots:
  //     - generate new predictions based upon the current point and
  //       insert into the search tree
  using NanoflannKdTreeT = nanoflann::KDTreeSingleIndexDynamicAdaptor<
    nanoflann::L2_Simple_Adaptor<int32_t, NanoflannCloud, int64_t>,
    NanoflannCloud, 3, int32_t>;

  using NanoflannResultT = nanoflann::KNNResultSet<int64_t, int32_t>;

  // the predicted point positions, used for searching.
  // each node will generate up to three predicted positions
  NanoflannCloud predictedPoints;
  NanoflannKdTreeT predictedPointsTree(3, predictedPoints);
  predictedPoints.pts.reserve(3 * pointCount);

  // mapping of predictedPoints indicies to prediction tree nodes
  std::vector<int32_t> predictedPointIdxToNodeIdx;
  predictedPointIdxToNodeIdx.reserve(3 * pointCount);

  // the prediction tree, one node for each point
  std::vector<GNode> nodes(pointCount);

  for (int nodeIdx = 0, nodeIdxN; nodeIdx < pointCount; nodeIdx = nodeIdxN) {
    auto& node = nodes[nodeIdx];
    auto queryPoint = begin[nodeIdx];

    // scan for duplicate points
    // NB: the tree coder assumes that duplicate point indices are consecutive
    // If this guarantee is changed, then the tree coder must be modified too.
    node.numDups = 0;
    for (nodeIdxN = nodeIdx + 1; nodeIdxN < pointCount; nodeIdxN++) {
      if (queryPoint != begin[nodeIdxN])
        break;
      node.numDups++;
    }

    int32_t nnPredictedPointIdx[GNode::MaxChildrenCount];
    int64_t nnPredictedDist[GNode::MaxChildrenCount];
    NanoflannResultT resultSet(GNode::MaxChildrenCount);

    resultSet.init(nnPredictedPointIdx, nnPredictedDist);
    predictedPointsTree.findNeighbors(resultSet, &queryPoint[0], {});

    // find a suitable parent.  default case: node has no parent
    node.parent = -1;
    node.childrenCount = 0;
    for (size_t r = 0; r < resultSet.size(); ++r) {
      auto parentIdx = predictedPointIdxToNodeIdx[nnPredictedPointIdx[r]];
      auto& pnode = nodes[parentIdx];
      if (pnode.childrenCount < GNode::MaxChildrenCount) {
        node.parent = parentIdx;
        pnode.children[pnode.childrenCount++] = nodeIdx;
        break;
      }
    }

    // update the search tree with new predictions from this point
    const auto size0 = predictedPoints.pts.size();

    // set the indicies for prediction
    GPredicter predicter;
    predicter.index[0] = nodeIdx;
    predicter.index[1] = nodes[predicter.index[0]].parent;
    if (predicter.index[1] >= 0)
      predicter.index[2] = nodes[predicter.index[1]].parent;

    for (int iMode = 0; iMode < 4; iMode++) {
      GPredicter::Mode mode = GPredicter::Mode(iMode);

      // don't add zero prediction
      if (mode == GPredicter::None)
        continue;

      if (!predicter.isValid(mode))
        continue;

      auto prediction = predicter.predict(begin, mode);
      predictedPointIdxToNodeIdx.push_back(nodeIdx);
      predictedPoints.pts.push_back(prediction);
    }

    const auto size1 = predictedPoints.pts.size();
    if (size0 != size1)
      predictedPointsTree.addPoints(size0, size1 - 1);
  }

  return nodes;
}

//============================================================================

static void
mortonSort(PCCPointSet3& cloud, int begin, int end, int depth)
{
  radixSort8(
    depth, PCCPointSet3::iterator(&cloud, begin),
    PCCPointSet3::iterator(&cloud, end),
    [=](int depth, const PCCPointSet3::Proxy& proxy) {
      const auto& point = *proxy;
      int mask = 1 << depth;
      return !!(point[2] & mask) | (!!(point[1] & mask) << 1)
        | (!!(point[0] & mask) << 2);
    });
}

//============================================================================

void
encodePredictiveGeometry(
  const PredGeomEncOpts& opt,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& cloud,
  EntropyEncoder* arithmeticEncoder)
{
  auto numPoints = cloud.getPointCount();

  // Origin relative to slice origin
  auto origin = gps.geomAngularOrigin - gbh.geomBoxOrigin;

  // storage for reordering the output point cloud
  PCCPointSet3 outCloud;
  outCloud.addRemoveAttributes(cloud.hasColors(), cloud.hasReflectances());
  outCloud.resize(numPoints);

  // src indexes in coded order
  std::vector<int32_t> codedOrder(numPoints, -1);

  // determine each geometry tree, and encode.  Size of trees is limited
  // by maxPtsPerTree.
  PredGeomEncoder enc(gps, arithmeticEncoder);
  int maxPtsPerTree = std::min(opt.maxPtsPerTree, int(numPoints));
  ;
  for (int i = 0; i < numPoints;) {
    int iEnd = std::min(i + maxPtsPerTree, int(numPoints));
    auto* begin = &cloud[i];
    auto* end = &cloud[0] + iEnd;

    // first, put the points in this tree into a sorted order
    // this can significantly improve the constructed tree
    if (opt.sortMode == PredGeomEncOpts::kSortMorton)
      mortonSort(cloud, i, iEnd, gbh.maxRootNodeDimLog2);
    else if (opt.sortMode == PredGeomEncOpts::kSortAzimuth)
      sortByAzimuth(cloud, i, iEnd, origin);
    else if (opt.sortMode == PredGeomEncOpts::kSortRadius)
      sortByRadius(cloud, i, iEnd, origin);

    // then build and encode the tree
    auto nodes = generateGeomPredictionTree(gps, gbh, begin, end);
    enc.encode(begin, nodes.data(), nodes.size(), codedOrder.data() + i);

    // put points in output cloud in decoded order
    for (auto iBegin = i; i < iEnd; i++) {
      auto srcIdx = iBegin + codedOrder[i];
      assert(srcIdx >= 0);
      outCloud[i] = cloud[srcIdx];
      if (cloud.hasColors())
        outCloud.setColor(i, cloud.getColor(srcIdx));
      if (cloud.hasReflectances())
        outCloud.setReflectance(i, cloud.getReflectance(srcIdx));
    }
  }

  swap(cloud, outCloud);
}

//============================================================================

}  // namespace pcc
