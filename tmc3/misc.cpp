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
  1,  1,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  5,  5,
  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  7,
  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,
  8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
  9,  9,  9,  9,  9,  9,  10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
  10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13,
  13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14,
  14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
  14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
  15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
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
// Fixed point inverse square root

namespace rsqrt {
  static const uint64_t k3timesR[96] = {
    3196059648, 3145728000, 3107979264, 3057647616, 3019898880, 2969567232,
    2931818496, 2894069760, 2868903936, 2831155200, 2793406464, 2768240640,
    2730491904, 2705326080, 2667577344, 2642411520, 2617245696, 2592079872,
    2566914048, 2541748224, 2516582400, 2491416576, 2466250752, 2441084928,
    2428502016, 2403336192, 2378170368, 2365587456, 2340421632, 2327838720,
    2302672896, 2290089984, 2264924160, 2252341248, 2239758336, 2214592512,
    2202009600, 2189426688, 2164260864, 2151677952, 2139095040, 2126512128,
    2113929216, 2101346304, 2088763392, 2076180480, 2051014656, 2038431744,
    2025848832, 2013265920, 2000683008, 2000683008, 1988100096, 1962934272,
    1962934272, 1950351360, 1937768448, 1925185536, 1912602624, 1900019712,
    1900019712, 1887436800, 1874853888, 1862270976, 1849688064, 1849688064,
    1837105152, 1824522240, 1811939328, 1811939328, 1799356416, 1786773504,
    1786773504, 1774190592, 1761607680, 1761607680, 1749024768, 1736441856,
    1736441856, 1723858944, 1723858944, 1711276032, 1698693120, 1698693120,
    1686110208, 1686110208, 1673527296, 1660944384, 1660944384, 1648361472,
    1648361472, 1635778560, 1635778560, 1623195648, 1623195648, 1610612736};

  static const uint64_t kRcubed[96] = {
    4195081216, 3999986688, 3857709056, 3673323520, 3538940928, 3364924416,
    3238224896, 3114735616, 3034196992, 2915990528, 2800922624, 2725880832,
    2615890944, 2544223232, 2439185408, 2370818048, 2303728640, 2237913088,
    2173355008, 2110061568, 2048008192, 1987165184, 1927563264, 1869150208,
    1840392192, 1783783424, 1728321536, 1701024768, 1647311872, 1620883456,
    1568898048, 1543306240, 1492993024, 1468236800, 1443762176, 1395656704,
    1372007424, 1348605952, 1302626304, 1280060416, 1257736192, 1235650560,
    1213861888, 1192294400, 1171008512, 1149979648, 1108673536, 1088379904,
    1068352512, 1048567808, 1029031936, 1029036032, 1009729536, 971888640,
    971882496,  953319424,  934993920,  916897792,  899011584,  881389568,
    881392640,  864009216,  846846976,  829900800,  813182976,  813201408,
    796721152,  780459008,  764412928,  764417024,  748601344,  732995584,
    733017088,  717624320,  702468096,  702466048,  687520768,  672786432,
    672787456,  658258944,  658256896,  643947520,  629854208,  629862400,
    615976960,  615952384,  602276864,  588779520,  588804096,  575512576,
    575526912,  562433024,  562439168,  549556224,  549564416,  536876032};
}  // namespace rsqrt

uint64_t
irsqrt(uint64_t a64)
{
  using namespace rsqrt;

  if (!a64)
    return 0;

  const uint64_t th32 = (uint64_t(1) << 31);
  int shift = -3;
  while (a64 & 0xffffffff00000000) {
    a64 >>= 2;
    shift--;
  }

  uint32_t a = a64;
  while (!(a & 0xc0000000)) {
    a <<= 2;
    shift++;
  }

  // initial approximation and first fixed point iteration:
  //   r' = (3r-ar^3)/2 (divide by 2 will by handled by shifts)
  int idx = (a >> 25) - 32;
  uint64_t r = k3timesR[idx] - ((kRcubed[idx] * a) >> 32);

  // second fixed point iteration
  uint64_t ar = (r * a) >> 32;
  uint64_t s = 0x30000000 - ((r * ar) >> 32);
  r = (r * s) >> 32;

  // denormalize
  if (shift > 0)
    return r << shift;
  else
    return r >> -shift;
}

//============================================================================

}  // namespace pcc
