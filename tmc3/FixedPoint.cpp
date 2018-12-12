/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2019, ISO/IEC
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

#include "FixedPoint.h"

namespace pcc {

//============================================================================

void
FixedPoint::operator/=(const FixedPoint& that)
{
  if (this->val < 0) {
    if (that.val < 0)
      this->val =
        -(((-that.val) >> 1) + ((-this->val) << kFracBits)) / that.val;
    else
      this->val =
        -(((+that.val) >> 1) + ((-this->val) << kFracBits)) / that.val;
  } else {
    if (that.val < 0)
      this->val =
        +(((-that.val) >> 1) + ((+this->val) << kFracBits)) / that.val;
    else
      this->val =
        +(((+that.val) >> 1) + ((+this->val) << kFracBits)) / that.val;
  }
}

static const uint32_t kSqrtLut[256] = {
  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,
  4,  4,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,
  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,
  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  10, 10, 10, 10,
  10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
  13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14,
  14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
  14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
  15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16, 16};

uint32_t
sqrtFixedpoint(uint64_t x)
{
  uint32_t a;

  // Initial guess
  if (x >= ((int64_t)0x1 << 32)) {
    if (x >= ((int64_t)0x1 << 48)) {
      if (x >= ((int64_t)0x1 << 56))
        a = (kSqrtLut[x >> 56] << 28) - 1;
      else
        a = kSqrtLut[x >> 48] << 24;
    } else {
      if (x >= ((int64_t)0x1 << 40))
        a = kSqrtLut[x >> 40] << 20;
      else
        a = kSqrtLut[x >> 32] << 16;
    }
  } else {
    if (x >= ((int64_t)0x1 << 16)) {
      if (x >= ((int64_t)0x1 << 24))
        a = kSqrtLut[x >> 24] << 12;
      else
        a = kSqrtLut[x >> 16] << 8;
    } else {
      if (x >= ((int64_t)0x1 << 8))
        a = kSqrtLut[x >> 8] << 4;
      else
        return kSqrtLut[x];
    }
  }

  a = (a + x / a) >> 1;
  return (a + x / a + 1) >> 1;
}

//============================================================================

}  // namespace pcc
