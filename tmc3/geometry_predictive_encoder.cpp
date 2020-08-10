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
#include "quantization.h"

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

class PredGeomEncoder : protected PredGeomContexts {
public:
  PredGeomEncoder(const PredGeomEncoder&) = delete;
  PredGeomEncoder& operator=(const PredGeomEncoder&) = delete;

  PredGeomEncoder(
    const GeometryParameterSet&,
    const GeometryBrickHeader&,
    const PredGeomContexts& ctxtMem,
    EntropyEncoder* aec);

  int qpSelector(const GNode& node) const { return _sliceQp; }

  void encode(
    const Vec3<int32_t>* cloudA,
    Vec3<int32_t>* cloudB,
    const GNode* nodes,
    int numNodes,
    int* codedOrder);

  int encodeTree(
    const Vec3<int32_t>* cloudA,
    Vec3<int32_t>* cloudB,
    const GNode* nodes,
    int numNodes,
    int rootIdx,
    int* codedOrder);

  void encodeNumDuplicatePoints(int numDupPoints);
  void encodeNumChildren(int numChildren);
  void encodePredMode(GPredicter::Mode mode);
  void encodeResidual(const Vec3<int32_t>& residual);
  void encodeResidual2(const Vec3<int32_t>& residual);
  void encodePhiMultiplier(const int32_t multiplier);
  void encodeQpOffset(int dqp);

  float estimateBits(GPredicter::Mode mode, const Vec3<int32_t>& residual);

  const PredGeomContexts& getCtx() const { return *this; }

private:
  EntropyEncoder* _aec;
  std::vector<int32_t> _stack;
  bool _geom_unique_points_flag;

  bool _geom_angular_mode_enabled_flag;
  Vec3<int32_t> origin;
  SphericalToCartesian _sphToCartesian;
  int _geom_angular_azimuth_speed;

  bool _geom_scaling_enabled_flag;
  int _sliceQp;
  int _qpOffsetInterval;

  Vec3<int> _maxAbsResidualMinus1Log2;
  Vec3<int> _pgeom_resid_abs_log2_bits;
};

//============================================================================

PredGeomEncoder::PredGeomEncoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const PredGeomContexts& ctxtMem,
  EntropyEncoder* aec)
  : PredGeomContexts(ctxtMem)
  , _aec(aec)
  , _geom_unique_points_flag(gps.geom_unique_points_flag)
  , _geom_angular_mode_enabled_flag(gps.geom_angular_mode_enabled_flag)
  , origin()
  , _sphToCartesian(gps)
  , _geom_angular_azimuth_speed(gps.geom_angular_azimuth_speed)
  , _geom_scaling_enabled_flag(gps.geom_scaling_enabled_flag)
  , _sliceQp(0)
  , _pgeom_resid_abs_log2_bits(gbh.pgeom_resid_abs_log2_bits)
{
  if (gps.geom_scaling_enabled_flag) {
    _sliceQp = gbh.sliceQp(gps);
    int qpIntervalLog2 =
      gps.geom_qp_offset_intvl_log2 + gbh.geom_qp_offset_intvl_log2_delta;
    _qpOffsetInterval = (1 << qpIntervalLog2) - 1;
  }

  if (gps.geom_angular_mode_enabled_flag)
    origin = gps.geomAngularOrigin - gbh.geomBoxOrigin;

  _stack.reserve(1024);

  for (int k = 0; k < 3; k++)
    _maxAbsResidualMinus1Log2[k] = (1 << gbh.pgeom_resid_abs_log2_bits[k]) - 1;
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeNumDuplicatePoints(int numDupPoints)
{
  _aec->encode(numDupPoints > 0, _ctxNumDupPointsGt0);
  if (numDupPoints)
    _aec->encodeExpGolomb(numDupPoints - 1, 0, _ctxNumDupPoints);
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
PredGeomEncoder::encodeResidual(const Vec3<int32_t>& residual)
{
  for (int k = 0, ctxIdx = 0; k < 3; k++) {
    const auto res = residual[k];
    const bool isZero = res == 0;
    _aec->encode(isZero, _ctxIsZero[k]);
    if (isZero)
      continue;

    _aec->encode(res > 0, _ctxSign[k]);

    int32_t value = abs(res) - 1;
    int32_t numBits = 1 + ilog2(uint32_t(value));

    AdaptiveBitModel* ctxs = _ctxNumBits[ctxIdx][k] - 1;
    for (int ctxIdx = 1, n = _pgeom_resid_abs_log2_bits[k] - 1; n >= 0; n--) {
      auto bin = (numBits >> n) & 1;
      _aec->encode(bin, ctxs[ctxIdx]);
      ctxIdx = (ctxIdx << 1) | bin;
    }

    if (!k && !_geom_angular_mode_enabled_flag)
      ctxIdx = (numBits + 1) >> 1;

    --numBits;
    for (int32_t i = 0; i < numBits; ++i)
      _aec->encode((value >> i) & 1);
  }
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeResidual2(const Vec3<int32_t>& residual)
{
  for (int k = 0; k < 3; k++) {
    const auto res = residual[k];
    const bool isZero = res == 0;
    _aec->encode(isZero, _ctxIsZero2[k]);
    if (isZero)
      continue;

    _aec->encode(res > 0, _ctxSign2[k]);

    const bool isOne = res == 1 || res == -1;
    _aec->encode(isOne, _ctxIsOne2[k]);
    if (isOne)
      continue;

    int32_t value = abs(res) - 2;
    auto& ctxs = _ctxResidual2[k];
    if (value < 15) {
      _aec->encode(value & 1, ctxs[0]);
      _aec->encode((value >> 1) & 1, ctxs[1 + (value & 1)]);
      _aec->encode((value >> 2) & 1, ctxs[3 + (value & 3)]);
      _aec->encode((value >> 3) & 1, ctxs[7 + (value & 7)]);
    } else {
      _aec->encode(1, ctxs[0]);
      _aec->encode(1, ctxs[2]);
      _aec->encode(1, ctxs[6]);
      _aec->encode(1, ctxs[14]);
      _aec->encodeExpGolomb(value - 15, 0, _ctxEG2[k]);
    }
  }
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodePhiMultiplier(int32_t multiplier)
{
  bool isZero = multiplier == 0;
  _aec->encode(isZero, _ctxIsZeroPhi);
  if (isZero)
    return;

  _aec->encode(multiplier > 0, _ctxSignPhi);

  int32_t value = abs(multiplier) - 1;
  bool isOne = !value;
  _aec->encode(isOne, _ctxIsOnePhi);
  if (isOne)
    return;

  value--;
  auto& ctxs = _ctxResidualPhi;
  if (value < 15) {
    _aec->encode(value & 1, ctxs[0]);
    _aec->encode((value >> 1) & 1, ctxs[1 + (value & 1)]);
    _aec->encode((value >> 2) & 1, ctxs[3 + (value & 3)]);
    _aec->encode((value >> 3) & 1, ctxs[7 + (value & 7)]);
  } else {
    _aec->encode(1, ctxs[0]);
    _aec->encode(1, ctxs[2]);
    _aec->encode(1, ctxs[6]);
    _aec->encode(1, ctxs[14]);
    _aec->encodeExpGolomb(value - 15, 0, _ctxEGPhi);
  }
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeQpOffset(int dqp)
{
  _aec->encode(dqp == 0, _ctxQpOffsetIsZero);
  if (dqp == 0) {
    return;
  }
  _aec->encode(dqp > 0, _ctxQpOffsetSign);
  _aec->encodeExpGolomb(abs(dqp) - 1, 0, _ctxQpOffsetAbsEgl);
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

    AdaptiveBitModel* ctxs = _ctxNumBits[ctxIdx][k] - 1;
    for (int ctxIdx = 1, n = _pgeom_resid_abs_log2_bits[k] - 1; n >= 0; n--) {
      auto bin = (numBits >> n) & 1;
      bits += estimate(bin, ctxs[ctxIdx]);
      ctxIdx = (ctxIdx << 1) | bin;
    }

    if (!k && !_geom_angular_mode_enabled_flag)
      ctxIdx = (numBits + 1) >> 1;

    bits += std::max(0, numBits - 1);
  }

  return bits;
}

//----------------------------------------------------------------------------

int
PredGeomEncoder::encodeTree(
  const Vec3<int32_t>* srcPts,
  Vec3<int32_t>* reconPts,
  const GNode* nodes,
  int numNodes,
  int rootIdx,
  int* codedOrder)
{
  QuantizerGeom quantizer(_sliceQp);
  int nodesUntilQpOffset = 0;
  int processedNodes = 0;

  _stack.push_back(rootIdx);
  while (!_stack.empty()) {
    const auto nodeIdx = _stack.back();
    _stack.pop_back();

    const auto& node = nodes[nodeIdx];
    const auto point = srcPts[nodeIdx];

    struct {
      float bits;
      GPredicter::Mode mode;
      Vec3<int32_t> residual;
      Vec3<int32_t> prediction;
    } best;

    if (_geom_scaling_enabled_flag && !nodesUntilQpOffset--) {
      int qp = qpSelector(node);
      quantizer = QuantizerGeom(qp);
      encodeQpOffset(qp - _sliceQp);
      nodesUntilQpOffset = _qpOffsetInterval;
    }

    // mode decision to pick best prediction from available set
    int qphi;
    for (int iMode = 0; iMode < 4; iMode++) {
      GPredicter::Mode mode = GPredicter::Mode(iMode);
      GPredicter predicter = makePredicter(
        nodeIdx, mode, [=](int idx) { return nodes[idx].parent; });

      if (!predicter.isValid(mode))
        continue;

      auto pred = predicter.predict(&srcPts[0], mode);
      if (_geom_angular_mode_enabled_flag) {
        if (iMode == GPredicter::Mode::Delta) {
          int32_t phi0 = srcPts[predicter.index[0]][1];
          int32_t phi1 = point[1];
          int32_t deltaPhi = phi1 - phi0;
          qphi = deltaPhi >= 0
            ? (deltaPhi + (_geom_angular_azimuth_speed >> 1))
              / _geom_angular_azimuth_speed
            : -(-deltaPhi + (_geom_angular_azimuth_speed >> 1))
              / _geom_angular_azimuth_speed;
          pred[1] += qphi * _geom_angular_azimuth_speed;
        }
      }

      // The residual in the spherical domain is loesslessly coded
      auto residual = point - pred;
      if (!_geom_angular_mode_enabled_flag)
        for (int k = 0; k < 3; k++)
          residual[k] = int32_t(quantizer.quantize(residual[k]));

      // Check if the prediction residual can be represented with the
      // current configuration.  If it can't, don't use this mode.
      bool isOverflow = false;
      for (int k = 0; k < 3; k++) {
        if (residual[k])
          if ((abs(residual[k]) - 1) >> _maxAbsResidualMinus1Log2[k])
            isOverflow = true;
      }
      if (isOverflow)
        continue;

      auto bits = estimateBits(mode, residual);

      if (iMode == 0 || bits < best.bits) {
        best.prediction = pred;
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

    if (
      _geom_angular_mode_enabled_flag && best.mode == GPredicter::Mode::Delta)
      encodePhiMultiplier(qphi);

    encodeResidual(best.residual);

    // convert spherical prediction to cartesian and re-calculate residual
    if (_geom_angular_mode_enabled_flag) {
      best.prediction = origin + _sphToCartesian(point);
      best.residual = reconPts[nodeIdx] - best.prediction;
      for (int k = 0; k < 3; k++)
        best.residual[k] = int32_t(quantizer.quantize(best.residual[k]));

      encodeResidual2(best.residual);
    }

    // write the reconstructed position back to the point cloud
    for (int k = 0; k < 3; k++)
      best.residual[k] = int32_t(quantizer.scale(best.residual[k]));
    reconPts[nodeIdx] = best.prediction + best.residual;

    for (int k = 0; k < 3; k++)
      reconPts[nodeIdx][k] = std::max(0, reconPts[nodeIdx][k]);

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
  const Vec3<int32_t>* cloudA,
  Vec3<int32_t>* cloudB,
  const GNode* nodes,
  int numNodes,
  int32_t* codedOrder)
{
  int32_t processedNodes = 0;
  for (int32_t rootIdx = 0; rootIdx < numNodes; rootIdx++) {
    // find the root node(s)
    if (nodes[rootIdx].parent >= 0)
      continue;

    int numSubtreeNodes = encodeTree(
      cloudA, cloudB, nodes, numNodes, rootIdx, codedOrder + processedNodes);
    processedNodes += numSubtreeNodes;
  }
  assert(processedNodes == numNodes);
}

//============================================================================

std::vector<GNode>
generateGeomPredictionTree(
  const GeometryParameterSet& gps,
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

//----------------------------------------------------------------------------

std::vector<GNode>
generateGeomPredictionTreeAngular(
  const GeometryParameterSet& gps,
  const Vec3<int32_t> origin,
  const Vec3<int32_t>* begin,
  const Vec3<int32_t>* end,
  Vec3<int32_t>* beginSph)
{
  int32_t pointCount = std::distance(begin, end);
  int32_t numLasers = gps.geom_angular_num_lidar_lasers();

  // the prediction tree, one node for each point
  std::vector<GNode> nodes(pointCount);
  std::vector<int32_t> prevNodes(numLasers, -1);
  std::vector<int32_t> firstNodes(numLasers, -1);

  CartesianToSpherical cartToSpherical(gps);

  for (int nodeIdx = 0, nodeIdxN; nodeIdx < pointCount; nodeIdx = nodeIdxN) {
    auto& node = nodes[nodeIdx];
    auto curPoint = begin[nodeIdx];
    node.childrenCount = 0;

    // scan for duplicate points
    // NB: the tree coder assumes that duplicate point indices are consecutive
    // If this guarantee is changed, then the tree coder must be modified too.
    node.numDups = 0;
    for (nodeIdxN = nodeIdx + 1; nodeIdxN < pointCount; nodeIdxN++) {
      if (curPoint != begin[nodeIdxN])
        break;
      node.numDups++;
    }

    // cartesian to spherical coordinates
    const auto carPos = curPoint - origin;
    auto& sphPos = beginSph[nodeIdx] = cartToSpherical(carPos);
    auto thetaIdx = sphPos[2];

    node.parent = prevNodes[thetaIdx];
    if (node.parent != -1) {
      auto& pnode = nodes[prevNodes[thetaIdx]];
      pnode.children[pnode.childrenCount++] = nodeIdx;
    } else
      firstNodes[thetaIdx] = nodeIdx;

    prevNodes[thetaIdx] = nodeIdx;
  }

  int32_t n0 = 0;
  while (firstNodes[n0] == -1)
    ++n0;

  for (int32_t n = n0 + 1, parentIdx = firstNodes[n0]; n < numLasers; ++n) {
    auto nodeIdx = firstNodes[n];
    if (nodeIdx < 0)
      continue;

    auto& pnode = nodes[parentIdx];
    if (pnode.childrenCount < GNode::MaxChildrenCount) {
      nodes[nodeIdx].parent = parentIdx;
      pnode.children[pnode.childrenCount++] = nodeIdx;
    }
    parentIdx = nodeIdx;
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
  GeometryBrickHeader& gbh,
  PCCPointSet3& cloud,
  PredGeomContexts& ctxtMem,
  EntropyEncoder* arithmeticEncoder)
{
  auto numPoints = cloud.getPointCount();

  // Origin relative to slice origin
  auto origin = gps.geomAngularOrigin - gbh.geomBoxOrigin;

  // storage for reordering the output point cloud
  PCCPointSet3 outCloud;
  outCloud.addRemoveAttributes(cloud.hasColors(), cloud.hasReflectances());
  outCloud.resize(numPoints);

  // storage for spherical point co-ordinates determined in angular mode
  std::vector<Vec3<int32_t>> sphericalPos;
  if (gps.geom_angular_mode_enabled_flag)
    sphericalPos.resize(numPoints);

  // src indexes in coded order
  std::vector<int32_t> codedOrder(numPoints, -1);

  // Assume that the number of bits required for residuals is equal to the
  // root node size.  This allows every position to be coded using PCM.
  for (int k = 0; k < 3; k++)
    gbh.pgeom_resid_abs_log2_bits[k] =
      ilog2(uint32_t(gbh.rootNodeSizeLog2[k])) + 1;

  // Number of residual bits bits for angular mode.  This is slightly
  // pessimistic in the calculation of r.
  if (gps.geom_angular_mode_enabled_flag) {
    auto xyzBboxLog2 = gbh.rootNodeSizeLog2;
    auto rDivLog2 = gps.geom_angular_radius_inv_scale_log2;
    auto azimuthBits = gps.geom_angular_azimuth_scale_log2;

    // first work out the maximum number of bits for the residual
    // NB: the slice coordinate system is used here: ie, minX|minY = 0
    int maxX = (1 << xyzBboxLog2[0]) - 1;
    int maxY = (1 << xyzBboxLog2[1]) - 1;
    int maxAbsDx = std::max(std::abs(origin[0]), std::abs(maxX - origin[0]));
    int maxAbsDy = std::max(std::abs(origin[1]), std::abs(maxY - origin[1]));
    auto r = std::round(std::hypot(maxAbsDx, maxAbsDy));

    Vec3<int> residualBits;
    residualBits[0] = ceillog2(divExp2RoundHalfUp(int64_t(r), rDivLog2));
    residualBits[1] = gps.geom_angular_azimuth_scale_log2;
    residualBits[2] = ceillog2(gps.geom_angular_theta_laser.size() - 1);

    // the number of prefix bits required
    for (int k = 0; k < 3; k++)
      gbh.pgeom_resid_abs_log2_bits[k] = ilog2(uint32_t(residualBits[k])) + 1;
  }

  // determine each geometry tree, and encode.  Size of trees is limited
  // by maxPtsPerTree.
  PredGeomEncoder enc(gps, gbh, ctxtMem, arithmeticEncoder);
  int maxPtsPerTree = std::min(opt.maxPtsPerTree, int(numPoints));

  for (int i = 0; i < numPoints;) {
    int iEnd = std::min(i + maxPtsPerTree, int(numPoints));
    auto* begin = &cloud[i];
    auto* beginSph = &sphericalPos[i];
    auto* end = &cloud[0] + iEnd;

    // first, put the points in this tree into a sorted order
    // this can significantly improve the constructed tree
    if (opt.sortMode == PredGeomEncOpts::kSortMorton)
      mortonSort(cloud, i, iEnd, gbh.maxRootNodeDimLog2);
    else if (opt.sortMode == PredGeomEncOpts::kSortAzimuth)
      sortByAzimuth(cloud, i, iEnd, opt.azimuthSortRecipBinWidth, origin);
    else if (opt.sortMode == PredGeomEncOpts::kSortRadius)
      sortByRadius(cloud, i, iEnd, origin);

    // then build and encode the tree
    auto nodes = gps.geom_angular_mode_enabled_flag
      ? generateGeomPredictionTreeAngular(gps, origin, begin, end, beginSph)
      : generateGeomPredictionTree(gps, begin, end);

    auto* a = gps.geom_angular_mode_enabled_flag ? beginSph : begin;
    auto* b = begin;

    enc.encode(a, b, nodes.data(), nodes.size(), codedOrder.data() + i);

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

  // save the context state for re-use by a future slice if required
  ctxtMem = enc.getCtx();

  swap(cloud, outCloud);
}

//============================================================================

}  // namespace pcc
