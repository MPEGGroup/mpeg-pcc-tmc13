
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
#include "hls.h"
#include "quantization.h"

#include <vector>

namespace pcc {

//============================================================================

class PredGeomDecoder : protected PredGeomContexts {
public:
  PredGeomDecoder(const PredGeomDecoder&) = delete;
  PredGeomDecoder& operator=(const PredGeomDecoder&) = delete;

  PredGeomDecoder(
    const GeometryParameterSet&,
    const GeometryBrickHeader& gbh,
    const PredGeomContexts& ctxtMem,
    EntropyDecoder* aed);

  /**
   * decodes a sequence of decoded geometry trees.
   * @returns the number of points decoded.
   */
  int decode(int numPoints, Vec3<int32_t>* outputPoints);

  /**
   * decodes a single predictive geometry tree.
   * @returns the number of points decoded.
   */
  int decodeTree(Vec3<int32_t>* outA, Vec3<int32_t>* outB);

  const PredGeomContexts& getCtx() const { return *this; }

private:
  int decodeNumDuplicatePoints();
  int decodeNumChildren();
  GPredicter::Mode decodePredMode();
  Vec3<int32_t> decodeResidual();
  Vec3<int32_t> decodeResidual2();
  int32_t decodePhiMultiplier(GPredicter::Mode mode);
  int32_t decodeQpOffset();

private:
  EntropyDecoder* _aed;
  std::vector<int32_t> _stack;
  std::vector<int32_t> _nodeIdxToParentIdx;
  bool _geom_unique_points_flag;

  bool _geom_angular_mode_enabled_flag;
  Vec3<int32_t> origin;
  int numLasers;
  SphericalToCartesian _sphToCartesian;
  int _geom_angular_azimuth_speed;

  bool _geom_scaling_enabled_flag;
  int _sliceQp;
  int _qpOffsetInterval;

  Vec3<int> _pgeom_resid_abs_log2_bits;
};

//============================================================================

PredGeomDecoder::PredGeomDecoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  const PredGeomContexts& ctxtMem,
  EntropyDecoder* aed)
  : PredGeomContexts(ctxtMem)
  , _aed(aed)
  , _geom_unique_points_flag(gps.geom_unique_points_flag)
  , _geom_angular_mode_enabled_flag(gps.geom_angular_mode_enabled_flag)
  , origin()
  , numLasers(gps.geom_angular_theta_laser.size())
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
}

//----------------------------------------------------------------------------

int
PredGeomDecoder::decodeNumDuplicatePoints()
{
  bool num_dup_points_gt0 = _aed->decode(_ctxNumDupPointsGt0);
  if (!num_dup_points_gt0)
    return 0;
  return 1 + _aed->decodeExpGolomb(0, _ctxNumDupPoints);
}

//----------------------------------------------------------------------------

int
PredGeomDecoder::decodeNumChildren()
{
  int numChildren = _aed->decode(_ctxNumChildren[0]);
  numChildren += _aed->decode(_ctxNumChildren[1 + numChildren]) << 1;
  return numChildren;
}

//----------------------------------------------------------------------------

GPredicter::Mode
PredGeomDecoder::decodePredMode()
{
  int mode = _aed->decode(_ctxPredMode[0]);
  mode += _aed->decode(_ctxPredMode[1 + mode]) << 1;
  return GPredicter::Mode(mode);
}

//----------------------------------------------------------------------------

Vec3<int32_t>
PredGeomDecoder::decodeResidual2()
{
  Vec3<int32_t> residual;
  for (int k = 0; k < 3; ++k) {
    if (_aed->decode(_ctxIsZero2[k])) {
      residual[k] = 0;
      continue;
    }

    auto sign = _aed->decode(_ctxSign2[k]);

    if (_aed->decode(_ctxIsOne2[k])) {
      residual[k] = sign ? 1 : -1;
      continue;
    }

    auto& ctxs = _ctxResidual2[k];
    int32_t value = _aed->decode(ctxs[0]);
    value += _aed->decode(ctxs[1 + (value & 1)]) << 1;
    value += _aed->decode(ctxs[3 + (value & 3)]) << 2;
    value += _aed->decode(ctxs[7 + (value & 7)]) << 3;
    if (value == 15)
      value += _aed->decodeExpGolomb(0, _ctxEG2[k]);

    residual[k] = sign ? (value + 2) : -(value + 2);
  }
  return residual;
}

//----------------------------------------------------------------------------

int32_t
PredGeomDecoder::decodePhiMultiplier(GPredicter::Mode mode)
{
  if (!_geom_angular_mode_enabled_flag || mode != GPredicter::Mode::Delta)
    return 0;

  if (_aed->decode(_ctxIsZeroPhi))
    return 0;

  const auto sign = _aed->decode(_ctxSignPhi);
  if (_aed->decode(_ctxIsOnePhi))
    return sign ? 1 : -1;

  auto& ctxs = _ctxResidualPhi;
  int32_t value = _aed->decode(ctxs[0]);
  value += _aed->decode(ctxs[1 + (value & 1)]) << 1;
  value += _aed->decode(ctxs[3 + (value & 3)]) << 2;
  value += _aed->decode(ctxs[7 + (value & 7)]) << 3;
  if (value == 15)
    value += _aed->decodeExpGolomb(0, _ctxEGPhi);

  return sign ? (value + 2) : -(value + 2);
}

//----------------------------------------------------------------------------

int32_t
PredGeomDecoder::decodeQpOffset()
{
  int dqp = 0;
  if (!_aed->decode(_ctxQpOffsetIsZero)) {
    int dqp_sign = _aed->decode(_ctxQpOffsetSign);
    dqp = _aed->decodeExpGolomb(0, _ctxQpOffsetAbsEgl) + 1;
    dqp = dqp_sign ? dqp : -dqp;
  }
  return dqp;
}

//----------------------------------------------------------------------------

Vec3<int32_t>
PredGeomDecoder::decodeResidual()
{
  Vec3<int32_t> residual;
  for (int k = 0, ctxIdx = 0; k < 3; ++k) {
    if (_aed->decode(_ctxIsZero[k])) {
      residual[k] = 0;
      continue;
    }

    auto sign = _aed->decode(_ctxSign[k]);

    AdaptiveBitModel* ctxs = _ctxNumBits[ctxIdx][k] - 1;
    int32_t numBits = 1;
    for (int n = 0; n < _pgeom_resid_abs_log2_bits[k]; n++)
      numBits = (numBits << 1) | _aed->decode(ctxs[numBits]);
    numBits ^= 1 << _pgeom_resid_abs_log2_bits[k];

    if (!k && !_geom_angular_mode_enabled_flag)
      ctxIdx = (numBits + 1) >> 1;

    int32_t res = 0;
    --numBits;
    if (numBits <= 0) {
      res = 2 + numBits;
    } else {
      res = 1 + (1 << numBits);
      for (int i = 0; i < numBits; ++i) {
        res += _aed->decode() << i;
      }
    }
    residual[k] = sign ? res : -res;
  }

  return residual;
}

//----------------------------------------------------------------------------

int
PredGeomDecoder::decodeTree(Vec3<int32_t>* outA, Vec3<int32_t>* outB)
{
  QuantizerGeom quantizer(_sliceQp);
  int nodesUntilQpOffset = 0;
  int nodeCount = 0;
  _stack.push_back(-1);

  while (!_stack.empty()) {
    auto parentNodeIdx = _stack.back();
    _stack.pop_back();

    if (_geom_scaling_enabled_flag && !nodesUntilQpOffset--) {
      int qpOffset = decodeQpOffset();
      int qp = _sliceQp + qpOffset;
      quantizer = QuantizerGeom(qp);
      nodesUntilQpOffset = _qpOffsetInterval;
    }

    // allocate point in traversal order (depth first)
    auto curNodeIdx = nodeCount++;
    _nodeIdxToParentIdx[curNodeIdx] = parentNodeIdx;

    int numDuplicatePoints = 0;
    if (!_geom_unique_points_flag)
      numDuplicatePoints = decodeNumDuplicatePoints();
    int numChildren = decodeNumChildren();
    auto mode = decodePredMode();
    int qphi = decodePhiMultiplier(mode);

    auto residual = decodeResidual();
    if (!_geom_angular_mode_enabled_flag)
      for (int k = 0; k < 3; k++)
        residual[k] = int32_t(quantizer.scale(residual[k]));

    auto predicter = makePredicter(
      curNodeIdx, mode, [&](int idx) { return _nodeIdxToParentIdx[idx]; });

    auto pred = predicter.predict(outA, mode);
    if (_geom_angular_mode_enabled_flag)
      if (mode == GPredicter::Mode::Delta)
        pred[1] += qphi * _geom_angular_azimuth_speed;

    auto pos = pred + residual;
    if (!_geom_angular_mode_enabled_flag)
      for (int k = 0; k < 3; k++)
        pos[k] = std::max(0, pos[k]);
    outA[curNodeIdx] = pos;

    // convert pos from spherical to cartesian, add secondary residual
    if (_geom_angular_mode_enabled_flag) {
      residual = decodeResidual2();
      for (int k = 0; k < 3; k++)
        residual[k] = int32_t(quantizer.scale(residual[k]));

      assert(pos[2] < numLasers && pos[2] >= 0);
      pred = origin + _sphToCartesian(pos);
      outB[curNodeIdx] = pred + residual;
      for (int k = 0; k < 3; k++)
        outB[curNodeIdx][k] = std::max(0, outB[curNodeIdx][k]);
    }

    // copy duplicate point output
    for (int i = 0; i < numDuplicatePoints; i++, nodeCount++) {
      outA[nodeCount] = outA[curNodeIdx];
      outB[nodeCount] = outB[curNodeIdx];
    }

    for (int i = 0; i < numChildren; i++)
      _stack.push_back(curNodeIdx);
  }

  return nodeCount;
}

//----------------------------------------------------------------------------

int
PredGeomDecoder::decode(int numPoints, Vec3<int32_t>* outputPoints)
{
  _nodeIdxToParentIdx.resize(numPoints);

  // An intermediate buffer used for reconstruction of the spherical
  // co-ordinates.
  auto* reconA = outputPoints;
  std::vector<Vec3<int32_t>> sphericalPos;
  if (_geom_angular_mode_enabled_flag) {
    sphericalPos.resize(numPoints);
    reconA = sphericalPos.data();
  }

  int32_t pointCount = 0;
  while (pointCount < numPoints) {
    auto numSubtreePoints = decodeTree(reconA, outputPoints);
    outputPoints += numSubtreePoints;
    reconA += numSubtreePoints;
    pointCount += numSubtreePoints;
  }

  return pointCount;
}

//============================================================================

void
decodePredictiveGeometry(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  PredGeomContexts& ctxtMem,
  EntropyDecoder* aed)
{
  PredGeomDecoder dec(gps, gbh, ctxtMem, aed);
  dec.decode(gbh.footer.geom_num_points_minus1 + 1, &pointCloud[0]);
  ctxtMem = dec.getCtx();
}

//============================================================================

}  // namespace pcc
