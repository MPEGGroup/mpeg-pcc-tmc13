
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

#include <vector>

namespace pcc {

//============================================================================

class PredGeomDecoder : public PredGeomCodec {
public:
  PredGeomDecoder(const PredGeomDecoder&) = delete;
  PredGeomDecoder& operator=(const PredGeomDecoder&) = delete;

  PredGeomDecoder(EntropyDecoder* aed) : _aed(aed) { _stack.reserve(1024); };

  /**
   * decodes a sequence of decoded geometry trees.
   * @returns the number of points decoded.
   */
  int decode(int numPoints, Vec3<int32_t>* outputPoints);

  /**
   * decodes a single predictive geometry tree.
   * @returns the number of points decoded.
   */
  int decodeTree(Vec3<int32_t>* outputPoints);

private:
  int decodeNumChildren();
  GPredicter::Mode decodePredMode();
  Vec3<int32_t> decodeResidual();

private:
  EntropyDecoder* _aed;
  std::vector<int32_t> _stack;
  std::vector<int32_t> _nodeIdxToParentIdx;
};

//============================================================================

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
PredGeomDecoder::decodeResidual()
{
  Vec3<int32_t> residual;
  for (int k = 0, ctxIdx = 0; k < 3; ++k) {
    if (_aed->decode(_ctxIsZero[k])) {
      residual[k] = 0;
      continue;
    }

    auto sign = _aed->decode(_ctxSign[k]);

    AdaptiveBitModel* ctxs = _ctxNumBits[ctxIdx][k];
    int32_t numBits;
    numBits = _aed->decode(ctxs[0]);
    numBits += _aed->decode(ctxs[1 + numBits]) << 1;
    numBits += _aed->decode(ctxs[3 + numBits]) << 2;
    numBits += _aed->decode(ctxs[7 + numBits]) << 3;
    numBits += _aed->decode(ctxs[15 + numBits]) << 4;

    if (!k)
      ctxIdx = (numBits + 1) >> 1;

    int32_t res = 0;
    --numBits;
    if (numBits <= 0) {
      res = 2 + numBits;
    } else {
      res = 1 + (1 << numBits);
      for (int i = 0; i < numBits; ++i) {
        res += _aed->decode(_ctxBypass) << i;
      }
    }
    residual[k] = sign ? res : -res;
  }

  return residual;
}

//----------------------------------------------------------------------------

int
PredGeomDecoder::decodeTree(Vec3<int32_t>* outputPoints)
{
  int nodeCount = 0;
  _stack.push_back(-1);

  while (!_stack.empty()) {
    auto parentNodeIdx = _stack.back();
    _stack.pop_back();

    // allocate point in traversal order (depth first)
    auto curNodeIdx = nodeCount++;
    _nodeIdxToParentIdx[curNodeIdx] = parentNodeIdx;

    int numChildren = decodeNumChildren();
    auto mode = decodePredMode();
    auto residual = decodeResidual();

    auto predicter = makePredicter(
      curNodeIdx, mode, [&](int idx) { return _nodeIdxToParentIdx[idx]; });

    auto pred = predicter.predict(outputPoints, mode);
    outputPoints[curNodeIdx] = pred + residual;

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

  int32_t pointCount = 0;
  while (pointCount < numPoints) {
    auto numSubtreePoints = decodeTree(outputPoints);
    outputPoints += numSubtreePoints;
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
  EntropyDecoder* aed)
{
  PredGeomDecoder dec(aed);
  dec.decode(gbh.geom_num_points_minus1 + 1, &pointCloud[0]);
}

//============================================================================

}  // namespace pcc
