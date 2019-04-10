/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2018-2019, ISO/IEC
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

#include "PCCMisc.h"

#include <cstdint>
#include <cstdio>
#include <string>

namespace pcc {

//============================================================================

std::string
expandNum(const std::string& src, int num)
{
  int idx = 0;
  auto nextChar = [&]() { return idx + 1 < src.size() ? src[idx + 1] : '\0'; };

  std::string out;
  out.reserve(src.size());

  // Search for each occurrence of a %d format pattern,
  //  - validate it, and pass off to snprintf for formatting
  for (; nextChar() != '\0'; idx++) {
    // Find start of pattern:
    int prev = idx;
    int start = idx = src.find('%', prev);

    // copy intermediate bytes
    if (start == std::string::npos) {
      out.append(src, prev, std::string::npos);
      break;
    }
    out.append(src, prev, start - prev);

    char c;
    while ((c = nextChar()) && (c == '#' || c == '0' || c == ' '))
      idx++; /* flags */

    if ((c = nextChar()) && (c >= '1' && c <= '9')) {
      idx++; /* width[0] */

      while ((c = nextChar()) && (c >= '0' && c <= '9'))
        idx++; /* width[>0] */
    }

    if ((c = nextChar()) && (c == '.')) {
      idx++; /* precision[0] */

      if ((c = nextChar()) && (c == '-'))
        idx++; /* precision[1] */

      while ((c = nextChar()) && (c >= '0' && c <= '9'))
        idx++; /* precision[>=1] */
    }

    c = nextChar();
    // Permit a %% escape to handle a literal %d
    if (c == '%' && idx == start) {
      out.push_back('%');
      idx++;
    } else if (c == 'd' && idx - start < 30) {
      char fmt[32];
      char buf[32];
      int fmtlen = src.copy(fmt, idx - start + 2, start);
      fmt[fmtlen] = '\0';
      int len = snprintf(buf, 32, fmt, num);
      if (len < 32) {
        out.append(buf);
      } else {
        out.append(fmt);
      }
      idx++;
    } else {
      out.append(src, start, 2);
      idx++;
    }
  }

  return out;
}

//============================================================================

static const int8_t kSqrtLut[256] = {
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

//============================================================================

uint32_t
isqrt(uint64_t x)
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
