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

#pragma once

#include <cstdint>

#include "entropy.h"
#include "PCCMath.h"

namespace pcc {

//============================================================================

struct GPredicter {
  enum Mode
  {
    None,
    Delta,
    Linear2,
    Linear3
  };

  int32_t index[3];

  bool isValid(Mode mode);
  Vec3<int32_t> predict(const Vec3<int32_t>* points, Mode mode);
};

//============================================================================

struct GNode {
  static const int32_t MaxChildrenCount = 3;

  int32_t parent;
  int32_t childrenCount;
  int32_t children[MaxChildrenCount];
};

//============================================================================

struct PredGeomCodec {
  StaticBitModel _ctxBypass;
  AdaptiveBitModel _ctxNumChildren[3];
  AdaptiveBitModel _ctxPredMode[3];
  AdaptiveBitModel _ctxIsZero[3];
  AdaptiveBitModel _ctxSign[3];
  AdaptiveBitModel _ctxNumBits[12][3][31];
};

//============================================================================

template<typename LookupFn>
GPredicter
makePredicter(
  int32_t curNodeIdx, GPredicter::Mode mode, LookupFn nodeIdxToParentIdx)
{
  GPredicter predIdx;
  switch (mode) {
  default:
  case GPredicter::None:
  case GPredicter::Delta:
  case GPredicter::Linear2:
  case GPredicter::Linear3:
    for (int i = 0; i < int(mode); i++) {
      if (curNodeIdx < 0)
        break;
      predIdx.index[i] = curNodeIdx = nodeIdxToParentIdx(curNodeIdx);
    }
    break;
  }
  return predIdx;
}

//============================================================================

inline bool
GPredicter::isValid(GPredicter::Mode mode)
{
  int numPredictors = int(mode);
  for (int i = 0; i < numPredictors; i++) {
    if (this->index[i] < 0)
      return false;
  }
  return true;
}

//============================================================================

inline Vec3<int32_t>
GPredicter::predict(const Vec3<int32_t>* points, GPredicter::Mode mode)
{
  Vec3<int32_t> pred;
  switch (mode) {
  case GPredicter::None: pred = 0; break;

  case GPredicter::Delta: {
    pred = points[this->index[0]];
    break;
  }

  case GPredicter::Linear2: {
    const auto& p0 = points[this->index[0]];
    const auto& p1 = points[this->index[1]];
    pred = 2 * p0 - p1;
    break;
  }

  default:
  case GPredicter::Linear3: {
    const auto& p0 = points[this->index[0]];
    const auto& p1 = points[this->index[1]];
    const auto& p2 = points[this->index[2]];
    pred = p0 + p1 - p2;
    break;
  }
  }
  return pred;
}

//============================================================================

}  // namespace pcc
