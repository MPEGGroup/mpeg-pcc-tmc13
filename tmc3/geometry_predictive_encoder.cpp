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
#include <algorithm>
#include <bitset>

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

  // NB: The following constraint must be honoured:
  //      qp % (1 << geom_qp_multiplier_log2) == 0
  int qpSelector(const GNode& node) const { return _sliceQp; }

  void encode(
    Vec3<int32_t>* cloudA,
    Vec3<int32_t>* cloudB,
    const GNode* nodes,
    int numNodes,
    int* codedOrder);

  int encodeTree(
    Vec3<int32_t>* cloudA,
    Vec3<int32_t>* cloudB,
    const GNode* nodes,
    int numNodes,
    int rootIdx,
    int* codedOrder);

  void encodeNumDuplicatePoints(int numDupPoints);
  void encodeNumChildren(int numChildren);
  void encodePredMode(GPredicter::Mode mode);
  void encodeResidual(const Vec3<int32_t>& residual, int iMode);
  void encodeResidual2(const Vec3<int32_t>& residual);
  void encodePhiMultiplier(const int32_t multiplier);
  void encodeQpOffset(int dqp);
  void encodeEndOfTreesFlag(int endFlag);

  float
  estimateBits(GPredicter::Mode mode, const Vec3<int32_t>& residual, int qphi);

  const PredGeomContexts& getCtx() const { return *this; }

  void setMinRadius(int value) { _pgeom_min_radius = value; }

private:
  EntropyEncoder* _aec;
  std::vector<int32_t> _stack;
  bool _geom_unique_points_flag;

  bool _geom_angular_mode_enabled_flag;
  Vec3<int32_t> origin;
  int _numLasers;
  SphericalToCartesian _sphToCartesian;
  bool _azimuth_scaling_enabled_flag;
  int _geomAngularAzimuthSpeed;

  bool _geom_scaling_enabled_flag;
  int _geom_qp_multiplier_log2;
  int _sliceQp;
  int _qpOffsetInterval;

  int _azimuthTwoPiLog2;

  Vec3<int> _maxAbsResidualMinus1Log2;
  Vec3<int> _pgeom_resid_abs_log2_bits;

  // Minimum radius used for prediction in angular coding
  int _pgeom_min_radius;
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
  , _numLasers(gps.numLasers())
  , _sphToCartesian(gps)
  , _azimuth_scaling_enabled_flag(gps.azimuth_scaling_enabled_flag)
  , _geomAngularAzimuthSpeed(gps.geom_angular_azimuth_speed_minus1 + 1)
  , _geom_scaling_enabled_flag(gps.geom_scaling_enabled_flag)
  , _geom_qp_multiplier_log2(gps.geom_qp_multiplier_log2)
  , _sliceQp(0)
  , _maxAbsResidualMinus1Log2((1 << gbh.pgeom_resid_abs_log2_bits) - 1)
  , _pgeom_resid_abs_log2_bits(gbh.pgeom_resid_abs_log2_bits)
  , _azimuthTwoPiLog2(gps.geom_angular_azimuth_scale_log2_minus11 + 12)
  , _pgeom_min_radius(gbh.pgeom_min_radius)
{
  if (gps.geom_scaling_enabled_flag) {
    _sliceQp = gbh.sliceQp(gps);
    int qpIntervalLog2 =
      gps.geom_qp_offset_intvl_log2 + gbh.geom_qp_offset_intvl_log2_delta;
    _qpOffsetInterval = (1 << qpIntervalLog2) - 1;
  }

  if (gps.geom_angular_mode_enabled_flag)
    origin = gbh.geomAngularOrigin(gps);
  ;

  _stack.reserve(1024);
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
  // Mapping order: 0, 1, 3, 2
  int val = numChildren ^ 1;

  _aec->encode(val > 0, _ctxNumChildren[0]);
  if (val > 0) {
    _aec->encode(val > 1, _ctxNumChildren[1]);
    if (val > 1)
      _aec->encode(val - 2, _ctxNumChildren[2]);
  }
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodePredMode(GPredicter::Mode mode)
{
  int iMode = int(mode);
  _aec->encode((iMode >> 1) & 1, _ctxPredMode[0]);
  _aec->encode(iMode & 1, _ctxPredMode[1 + (iMode >> 1)]);
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeResidual(const Vec3<int32_t>& residual, int iMode)
{
  for (int k = 0, ctxIdx = 0; k < 3; k++) {
    // The last component (delta laseridx) isn't coded if there is one laser
    if (_geom_angular_mode_enabled_flag && _numLasers == 1 && k == 2)
      continue;

    const auto res = residual[k];
    _aec->encode(res != 0, _ctxResGt0[k]);
    if (!res)
      continue;

    int32_t value = abs(res) - 1;
    int32_t numBits = 1 + ilog2(uint32_t(value));

    AdaptiveBitModel* ctxs = &_ctxNumBits[ctxIdx][k][0] - 1;
    for (int ctxIdx = 1, n = _pgeom_resid_abs_log2_bits[k] - 1; n >= 0; n--) {
      auto bin = (numBits >> n) & 1;
      _aec->encode(bin, ctxs[ctxIdx]);
      ctxIdx = (ctxIdx << 1) | bin;
    }

    if (!k && !_geom_angular_mode_enabled_flag)
      ctxIdx = std::min(4, (numBits + 1) >> 1);

    --numBits;
    for (int32_t i = 0; i < numBits; ++i)
      _aec->encode((value >> i) & 1);

    if (iMode || k)
      _aec->encode(res < 0, _ctxSign[k]);
  }
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeResidual2(const Vec3<int32_t>& residual)
{
  for (int k = 0; k < 3; k++) {
    const auto res = residual[k];
    _aec->encode(res != 0, _ctxResidual2GtN[0][k]);
    if (!res)
      continue;

    int value = abs(res) - 1;
    _aec->encode(value > 0, _ctxResidual2GtN[1][k]);
    if (value)
      _aec->encodeExpGolomb(value - 1, 0, _ctxEG2Prefix[k], _ctxEG2Suffix[k]);

    _aec->encode(res < 0, _ctxSign2[k]);
  }
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodePhiMultiplier(int32_t multiplier)
{
  _aec->encode(multiplier != 0, _ctxPhiGtN[0]);
  if (!multiplier)
    return;

  int32_t value = abs(multiplier) - 1;
  _aec->encode(value > 0, _ctxPhiGtN[1]);
  if (!value) {
    _aec->encode(multiplier < 0, _ctxSignPhi);
    return;
  }

  value--;
  int valueMinus7 = value - 7;
  value = std::min(value, 7);
  _aec->encode((value >> 2) & 1, _ctxResidualPhi[0]);
  _aec->encode((value >> 1) & 1, _ctxResidualPhi[1 + (value >> 2)]);
  _aec->encode((value >> 0) & 1, _ctxResidualPhi[3 + (value >> 1)]);

  if (valueMinus7 >= 0)
    _aec->encodeExpGolomb(valueMinus7, 0, _ctxEGPhi);

  _aec->encode(multiplier < 0, _ctxSignPhi);
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeQpOffset(int dqp)
{
  _aec->encode(dqp != 0, _ctxQpOffsetAbsGt0);
  if (dqp == 0)
    return;

  _aec->encode(dqp < 0, _ctxQpOffsetSign);
  _aec->encodeExpGolomb(abs(dqp) - 1, 0, _ctxQpOffsetAbsEgl);
}

//----------------------------------------------------------------------------

void
PredGeomEncoder::encodeEndOfTreesFlag(int end_of_trees_flag)
{
  _aec->encode(end_of_trees_flag, _ctxEndOfTrees);
}

//----------------------------------------------------------------------------

float
PredGeomEncoder::estimateBits(
  GPredicter::Mode mode, const Vec3<int32_t>& residual, int multiplier)
{
  int iMode = int(mode);
  float bits = 0.;
  bits += estimate((iMode >> 1) & 1, _ctxPredMode[0]);
  bits += estimate(iMode & 1, _ctxPredMode[1 + (iMode >> 1)]);

  if (_geom_angular_mode_enabled_flag) {
    bits += estimate(multiplier != 0, _ctxPhiGtN[0]);

    if (multiplier) {
      int32_t value = abs(multiplier) - 1;
      bits += estimate(value > 0, _ctxPhiGtN[1]);
      bits += estimate(multiplier < 0, _ctxSignPhi);
      if (value) {
        value--;

        int valueMinus7 = value - 7;
        value = std::min(value, 7);
        bits += estimate((value >> 2) & 1, _ctxResidualPhi[0]);
        bits += estimate((value >> 1) & 1, _ctxResidualPhi[1 + (value >> 2)]);
        bits += estimate((value >> 0) & 1, _ctxResidualPhi[3 + (value >> 1)]);

        if (valueMinus7 >= 0)
          bits += (1 + 2.0 * log2(valueMinus7 + 1));
      }
    }
  }

  for (int k = 0, ctxIdx = 0; k < 3; k++) {
    // The last component (delta laseridx) isn't coded if there is one laser
    if (_geom_angular_mode_enabled_flag && _numLasers == 1 && k == 2)
      continue;

    const auto res = residual[k];
    bits += estimate(res != 0, _ctxResGt0[k]);
    if (res == 0)
      continue;

    if (iMode > 0 || k)
      bits += estimate(res < 0, _ctxSign[k]);

    int32_t value = abs(res) - 1;
    int32_t numBits = 1 + ilog2(uint32_t(value));

    AdaptiveBitModel* ctxs = &_ctxNumBits[ctxIdx][k][0] - 1;
    for (int ctxIdx = 1, n = _pgeom_resid_abs_log2_bits[k] - 1; n >= 0; n--) {
      auto bin = (numBits >> n) & 1;
      bits += estimate(bin, ctxs[ctxIdx]);
      ctxIdx = (ctxIdx << 1) | bin;
    }

    if (!k && !_geom_angular_mode_enabled_flag)
      ctxIdx = std::min(4, (numBits + 1) >> 1);

    bits += std::max(0, numBits - 1);
  }

  return bits;
}

//----------------------------------------------------------------------------

int
PredGeomEncoder::encodeTree(
  Vec3<int32_t>* srcPts,
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
    const auto& point = srcPts[nodeIdx];

    struct {
      float bits;
      GPredicter::Mode mode;
      Vec3<int32_t> residual;
      Vec3<int32_t> prediction;
      int qphi;
    } best;

    if (_geom_scaling_enabled_flag && !nodesUntilQpOffset--) {
      int qp = qpSelector(node);
      quantizer = QuantizerGeom(qp);
      encodeQpOffset((qp - _sliceQp) >> _geom_qp_multiplier_log2);
      nodesUntilQpOffset = _qpOffsetInterval;
    }

    // mode decision to pick best prediction from available set
    int qphi = 0;
    std::bitset<4> unusable;
    for (int iMode = 0; iMode < 4; iMode++) {
      GPredicter::Mode mode = GPredicter::Mode(iMode);
      GPredicter predicter =
        makePredicter(nodeIdx, mode, _pgeom_min_radius, [=](int idx) {
          return nodes[idx].parent;
        });

      if (!predicter.isValid(mode))
        continue;

      auto pred =
        predicter.predict(&srcPts[0], mode, _geom_angular_mode_enabled_flag);
      if (_geom_angular_mode_enabled_flag) {
        int32_t phi0 = pred[1];
        int32_t phi1 = point[1];
        int32_t deltaPhi = phi1 - phi0;
        qphi = deltaPhi >= 0 ? (deltaPhi + (_geomAngularAzimuthSpeed >> 1))
            / _geomAngularAzimuthSpeed
                             : -(-deltaPhi + (_geomAngularAzimuthSpeed >> 1))
            / _geomAngularAzimuthSpeed;
        pred[1] += qphi * _geomAngularAzimuthSpeed;
      }

      auto residual = point - pred;
      if (!_geom_angular_mode_enabled_flag)
        for (int k = 0; k < 3; k++)
          residual[k] = int32_t(quantizer.quantize(residual[k]));
      else {
        if (_azimuth_scaling_enabled_flag) {
          // Quantize phi by 8r/2^n
          auto arc = int64_t(residual[1]) * (pred[0] + residual[0]) << 3;
          residual[1] = int32_t(divExp2RoundHalfInf(arc, _azimuthTwoPiLog2));
        }
      }

      // Check if the prediction residual can be represented with the
      // current configuration.  If it can't, don't use this mode.
      for (int k = 0; k < 3; k++) {
        if (residual[k])
          if ((abs(residual[k]) - 1) >> _maxAbsResidualMinus1Log2[k])
            unusable[iMode] = true;
      }

      if (unusable[iMode]) {
        if (iMode == 3 && unusable.all())
          throw std::runtime_error(
            "predgeom: can't represent residual in any mode");
        if (iMode > 0)
          continue;
      }

      // penalize bit cost so that it is only used if all modes overfow
      auto bits = estimateBits(mode, residual, qphi);
      if (unusable[iMode])
        bits = std::numeric_limits<decltype(bits)>::max();

      if (iMode == 0 || bits < best.bits) {
        best.prediction = pred;
        best.residual = residual;
        best.mode = mode;
        best.bits = bits;
        best.qphi = qphi;
      }
    }

    assert(node.childrenCount <= GNode::MaxChildrenCount);
    if (!_geom_unique_points_flag)
      encodeNumDuplicatePoints(node.numDups);
    encodeNumChildren(node.childrenCount);
    encodePredMode(best.mode);

    if (_geom_angular_mode_enabled_flag)
      encodePhiMultiplier(best.qphi);

    encodeResidual(best.residual, best.mode);

    // convert spherical prediction to cartesian and re-calculate residual
    if (_geom_angular_mode_enabled_flag) {
      if (_azimuth_scaling_enabled_flag) {
        auto r = (best.prediction[0] + best.residual[0]) << 3;
        if (!r)
          r = 1;

        // scale the azimuth residual
        int32_t rInvScaleLog2;
        int64_t rInv = recipApprox(r, rInvScaleLog2);
        best.residual[1] =
          divExp2(best.residual[1] * rInv, rInvScaleLog2 - _azimuthTwoPiLog2);

        // update original spherical position with reconstructed spherical
        // position, for use in, for instance, attribute coding.
        srcPts[nodeIdx] /* == point */ = best.prediction + best.residual;
        for (int i = 1; i <= node.numDups; i++)
          srcPts[nodeIdx + i] = srcPts[nodeIdx];
      }

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
  Vec3<int32_t>* cloudA,
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

    // NB: this is just in case this call to encode needs to encode an
    // additional tree.  The actual end_of_trees_flag = true is signalled
    // at the end.
    if (processedNodes != numNodes)
      encodeEndOfTreesFlag(false);
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

      auto prediction = predicter.predict(begin, mode, false);
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
  int32_t numLasers = gps.numLasers();

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

    // propagate converted coordinates over duplicate points
    for (int i = nodeIdx + 1; i < nodeIdxN; i++)
      beginSph[i] = sphPos;

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

Vec3<int>
originFromLaserAngle(const PCCPointSet3& cloud)
{
  Vec3<int> origin;
  auto numPoints = cloud.getPointCount();
  int i;

  for (i = 0; i < numPoints; i++) {
    if (cloud.getLaserAngle(i) == 90) {
      origin = cloud[i];
      break;
    }
  }

  for (; i < numPoints; i++) {
    if (cloud.getLaserAngle(i) == 90) {
      // the most left point
      if (origin[0] > cloud[i][0])
        origin = cloud[i];
    }
  }

  return origin;
}

//============================================================================

void
encodePredictiveGeometry(
  const PredGeomEncOpts& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& cloud,
  std::vector<point_t>* reconPosSph,
  PredGeomContexts& ctxtMem,
  EntropyEncoder* arithmeticEncoder)
{
  auto numPoints = cloud.getPointCount();

  // Origin relative to slice origin
  //  - if angular disabled, try to find a presort origin using laser angles
  Vec3<int> origin = gbh.geomAngularOrigin(gps);
  if (!gps.geom_angular_mode_enabled_flag && cloud.hasLaserAngles())
    origin = originFromLaserAngle(cloud);

  // storage for reordering the output point cloud
  PCCPointSet3 outCloud;
  outCloud.addRemoveAttributes(cloud.hasColors(), cloud.hasReflectances());
  outCloud.resize(numPoints);

  // storage for spherical point co-ordinates determined in angular mode
  std::vector<Vec3<int32_t>> sphericalPos;
  if (!gps.geom_angular_mode_enabled_flag)
    reconPosSph = nullptr;
  else {
    sphericalPos.resize(numPoints);
    if (reconPosSph)
      reconPosSph->resize(numPoints);
  }

  // src indexes in coded order
  std::vector<int32_t> codedOrder(numPoints, -1);

  // Assume that the number of bits required for residuals is equal to the
  // quantised root node size.  This allows every position to be coded with PCM
  if (!gps.geom_angular_mode_enabled_flag) {
    QuantizerGeom quant(gbh.sliceQp(gps));
    for (int k = 0; k < 3; k++) {
      int max = quant.quantize((1 << gbh.rootNodeSizeLog2[k]) - 1);
      gbh.pgeom_resid_abs_log2_bits[k] = numBits(ceillog2(std::max(1, max)));
    }
  }

  // Number of residual bits bits for angular mode.  This is slightly
  // pessimistic in the calculation of r.
  if (gps.geom_angular_mode_enabled_flag) {
    auto xyzBboxLog2 = gbh.rootNodeSizeLog2;
    auto rDivLog2 = gps.geom_angular_radius_inv_scale_log2;

    // first work out the maximum number of bits for the residual
    // NB: the slice coordinate system is used here: ie, minX|minY = 0
    int maxX = (1 << xyzBboxLog2[0]) - 1;
    int maxY = (1 << xyzBboxLog2[1]) - 1;
    int maxAbsDx = std::max(std::abs(origin[0]), std::abs(maxX - origin[0]));
    int maxAbsDy = std::max(std::abs(origin[1]), std::abs(maxY - origin[1]));
    auto r = (int64_t)std::round(std::hypot(maxAbsDx, maxAbsDy));

    Vec3<int> residualBits;
    residualBits[0] = ceillog2(divExp2RoundHalfUp(r, rDivLog2));
    residualBits[2] = ceillog2(gps.numLasers() - 1);
    if (!gps.azimuth_scaling_enabled_flag)
      residualBits[1] =
        ceillog2((gps.geom_angular_azimuth_speed_minus1 + 1) >> 1);
    else {
      int maxError = ((gps.geom_angular_azimuth_speed_minus1 + 1) >> 1) + 1;
      int twoPi = gps.geom_angular_azimuth_scale_log2_minus11 + 12;
      residualBits[1] = ceillog2(divExp2RoundHalfInf(
        maxError * divExp2RoundHalfUp(r << 3, rDivLog2), twoPi));
    }

    // the number of prefix bits required
    for (int k = 0; k < 3; k++)
      gbh.pgeom_resid_abs_log2_bits[k] = ilog2(uint32_t(residualBits[k])) + 1;
  }

  // set the minimum value to zero, if encoder don't use this info.
  gbh.pgeom_min_radius = 0;

  // determine each geometry tree, and encode.  Size of trees is limited
  // by maxPtsPerTree.
  PredGeomEncoder enc(gps, gbh, ctxtMem, arithmeticEncoder);
  int maxPtsPerTree = std::min(opt.maxPtsPerTree, int(numPoints));

  for (int i = 0; i < numPoints;) {
    int iEnd = std::min(i + maxPtsPerTree, int(numPoints));
    auto* begin = &cloud[i];
    auto* beginSph = sphericalPos.data() + i;
    auto* end = &cloud[0] + iEnd;

    // first, put the points in this tree into a sorted order
    // this can significantly improve the constructed tree
    if (opt.sortMode == PredGeomEncOpts::kSortMorton)
      mortonSort(cloud, i, iEnd, gbh.maxRootNodeDimLog2);
    else if (opt.sortMode == PredGeomEncOpts::kSortAzimuth)
      sortByAzimuth(cloud, i, iEnd, opt.azimuthSortRecipBinWidth, origin);
    else if (opt.sortMode == PredGeomEncOpts::kSortRadius)
      sortByRadius(cloud, i, iEnd, origin);
    else if (opt.sortMode == PredGeomEncOpts::kSortLaserAngle)
      sortByLaserAngle(cloud, i, iEnd, opt.azimuthSortRecipBinWidth, origin);

    // then build and encode the tree
    auto nodes = gps.geom_angular_mode_enabled_flag
      ? generateGeomPredictionTreeAngular(gps, origin, begin, end, beginSph)
      : generateGeomPredictionTree(gps, begin, end);

    // Determining minimum radius for prediction.
    // NB: this value is per-slice, but if the slice generates multiple
    // trees due to the maxPtsPerTree limit, the radius of all points in
    // cthe slice is not known and the feature disabled.
    if (gps.geom_angular_mode_enabled_flag && numPoints <= maxPtsPerTree) {
      int min = beginSph[i][0];
      for (int j = i + 1; j < iEnd; j++)
        min = std::min(min, beginSph[j][0]);
      gbh.pgeom_min_radius = min;
      enc.setMinRadius(gbh.pgeom_min_radius);
    }

    auto* a = gps.geom_angular_mode_enabled_flag ? beginSph : begin;
    auto* b = begin;

    if (i > 0)
      enc.encodeEndOfTreesFlag(false);
    enc.encode(a, b, nodes.data(), nodes.size(), codedOrder.data() + i);

    // put points in output cloud in decoded order
    for (auto iBegin = i; i < iEnd; i++) {
      auto srcIdx = iBegin + codedOrder[i];
      assert(srcIdx >= 0);
      outCloud[i] = cloud[srcIdx];
      if (reconPosSph)
        (*reconPosSph)[i] = sphericalPos[srcIdx];
      if (cloud.hasColors())
        outCloud.setColor(i, cloud.getColor(srcIdx));
      if (cloud.hasReflectances())
        outCloud.setReflectance(i, cloud.getReflectance(srcIdx));
    }
  }

  // the slice has finished
  enc.encodeEndOfTreesFlag(true);

  // save the context state for re-use by a future slice if required
  ctxtMem = enc.getCtx();

  swap(cloud, outCloud);
}

//============================================================================

}  // namespace pcc
