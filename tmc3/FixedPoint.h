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

#pragma once

#include <cstdint>

namespace pcc {

//============================================================================

class FixedPoint {
public:
  // Number of fractional bits in fixed point representation
  static const int kFracBits = 15;
  static const int kOneHalf = 1 << (kFracBits - 1);

  // Fixed point value
  int64_t val;

  FixedPoint() = default;
  FixedPoint(const FixedPoint&) = default;
  FixedPoint(FixedPoint&&) = default;
  FixedPoint& operator=(const FixedPoint&) = default;
  FixedPoint& operator=(FixedPoint&&) = default;

  FixedPoint(int val) { this->operator=(int64_t(val)); }
  FixedPoint(int64_t val) { this->operator=(val); }
  FixedPoint(double val) { this->val = int64_t(val * (1 << kFracBits)); }

  // return the rounded integer value
  int64_t round();

  void operator=(const int64_t val);

  void operator+=(const FixedPoint& that);
  void operator-=(const FixedPoint& that);
  void operator*=(const FixedPoint& that);
  void operator/=(const FixedPoint& that);
};

//============================================================================

inline int64_t
FixedPoint::round()
{
  if (this->val > 0)
    return (kOneHalf + this->val) >> kFracBits;
  return -((kOneHalf - this->val) >> kFracBits);
}

//----------------------------------------------------------------------------

inline void
FixedPoint::operator=(int64_t val)
{
  if (val > 0)
    this->val = val << kFracBits;
  else
    this->val = -((-val) << kFracBits);
}

//----------------------------------------------------------------------------

inline void
FixedPoint::operator+=(const FixedPoint& that)
{
  this->val += that.val;
}

//----------------------------------------------------------------------------

inline void
FixedPoint::operator-=(const FixedPoint& that)
{
  this->val -= that.val;
}

//----------------------------------------------------------------------------

inline void
FixedPoint::operator*=(const FixedPoint& that)
{
  this->val *= that.val;

  if (this->val < 0)
    this->val = -((kOneHalf - this->val) >> kFracBits);
  else
    this->val = +((kOneHalf + this->val) >> kFracBits);
}

//============================================================================

}  // namespace pcc
