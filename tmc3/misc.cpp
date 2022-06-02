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
#include "PCCMath.h"

#include <cmath>
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
  if (x <= (uint64_t(1) << 46))
    return 1 + ((x * irsqrt(x)) >> 40);
  else {
    uint64_t x0 = (x + 65536) >> 16;
    return 1 + ((x0 * irsqrt(x0)) >> 32);
  }
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
// Fixed point arc tangent

namespace atan2 {
  // LUT with 20 bit precision on angle
  const int kAsin[364] = {
    0,      2048,   4096,   6144,   8192,   10240,  12288,  14336,  16385,
    18433,  20481,  22530,  24578,  26627,  28676,  30724,  32773,  34822,
    36872,  38921,  40970,  43020,  45070,  47120,  49170,  51220,  53271,
    55322,  57373,  59424,  61475,  63527,  65579,  67631,  69683,  71736,
    73789,  75842,  77896,  79949,  82004,  84058,  86113,  88168,  90223,
    92279,  94335,  96392,  98449,  100506, 102563, 104621, 106680, 108739,
    110798, 112858, 114918, 116978, 119040, 121101, 123163, 125225, 127288,
    129352, 131416, 133480, 135545, 137611, 139677, 141743, 143810, 145878,
    147946, 150015, 152085, 154155, 156225, 158297, 160368, 162441, 164514,
    166588, 168662, 170737, 172813, 174890, 176967, 179045, 181123, 183203,
    185283, 187363, 189445, 191527, 193610, 195694, 197779, 199864, 201950,
    204037, 206125, 208214, 210303, 212393, 214485, 216577, 218669, 220763,
    222858, 224954, 227050, 229148, 231246, 233345, 235445, 237547, 239649,
    241752, 243856, 245961, 248068, 250175, 252283, 254392, 256502, 258614,
    260726, 262840, 264954, 267070, 269187, 271305, 273424, 275544, 277666,
    279788, 281912, 284037, 286163, 288290, 290419, 292549, 294680, 296812,
    298945, 301080, 303216, 305354, 307492, 309632, 311773, 313916, 316060,
    318206, 320352, 322500, 324650, 326801, 328953, 331107, 333262, 335419,
    337577, 339737, 341898, 344061, 346225, 348391, 350558, 352727, 354897,
    357069, 359243, 361418, 363595, 365773, 367953, 370135, 372318, 374503,
    376690, 378879, 381069, 383261, 385455, 387650, 389847, 392046, 394247,
    396450, 398655, 400861, 403069, 405279, 407491, 409705, 411921, 414139,
    416359, 418581, 420804, 423030, 425258, 427488, 429720, 431954, 434190,
    436428, 438668, 440910, 443155, 445401, 447650, 449901, 452155, 454410,
    456668, 458928, 461190, 463455, 465722, 467991, 470262, 472536, 474813,
    477091, 479373, 481656, 483942, 486231, 488522, 490815, 493111, 495410,
    497711, 500015, 502322, 504631, 506943, 509257, 511574, 513894, 516217,
    518542, 520870, 523201, 525535, 527872, 530211, 532553, 534899, 537247,
    539598, 541952, 544310, 546670, 549033, 551399, 553769, 556142, 558517,
    560896, 563278, 565664, 568052, 570444, 572839, 575238, 577640, 580045,
    582454, 584866, 587282, 589701, 592123, 594549, 596979, 599412, 601849,
    604290, 606734, 609183, 611634, 614090, 616549, 619013, 621480, 623951,
    626426, 628905, 631388, 633875, 636366, 638862, 641361, 643865, 646373,
    648885, 651401, 653922, 656447, 658976, 661510, 664049, 666592, 669139,
    671691, 674248, 676809, 679375, 681946, 684522, 687103, 689688, 692278,
    694874, 697474, 700080, 702690, 705306, 707927, 710553, 713184, 715821,
    718463, 721111, 723764, 726423, 729087, 731757, 734433, 737115, 739802,
    742495, 745194, 747899, 750611, 753328, 756051, 758781, 761517, 764259,
    767008, 769763, 772525, 775294, 778069, 780850, 783639, 786435, 789237,
    792047, 794863, 797687, 800518, 803357, 806202, 809056, 811917, 814785,
    817662, 820546, 823438, 823438};
}  // namespace atan2

//----------------------------------------------------------------------------

int
iatan2Core(int y, int x)
{
  using namespace atan2;

  // ratio y/x knowing that x,y >0 and y<=x
  if (x == 0)
    return 0;

  uint64_t rinv =
    irsqrt(uint64_t(x) * uint64_t(x) + uint64_t(y) * uint64_t(y));
  int r = (y * rinv) >> 20;  // 40 - 20 = 20 bits precision
  int idx = r >> 11;
  int lambda = r - (idx << 11);
  return kAsin[idx] + (lambda * (kAsin[idx + 1] - kAsin[idx]) >> 11);
}

//----------------------------------------------------------------------------

int
iatan2(int y, int x)
{
  int xa = std::abs(x);
  int ya = std::abs(y);

  // atan or pi/2 - atan
  int t = ya <= xa ? iatan2Core(ya, xa) : 1647099 - iatan2Core(xa, ya);
  if (x < 0)
    t = 3294199 - t;  // pi -atan

  return y < 0 ? -t : t;
}

//============================================================================

const uint16_t kDivApproxDivisor[256] = {
  65535, 32767, 21844, 16383, 13106, 10922, 9361, 8191, 7281, 6553, 5957, 5460,
  5040,  4680,  4368,  4095,  3854,  3640,  3448, 3276, 3120, 2978, 2848, 2730,
  2620,  2520,  2426,  2340,  2259,  2184,  2113, 2047, 1985, 1927, 1871, 1819,
  1770,  1724,  1679,  1637,  1597,  1559,  1523, 1488, 1455, 1424, 1393, 1364,
  1336,  1310,  1284,  1259,  1236,  1213,  1191, 1169, 1149, 1129, 1110, 1091,
  1073,  1056,  1039,  1023,  1007,  992,   977,  963,  949,  935,  922,  909,
  897,   885,   873,   861,   850,   839,   829,  818,  808,  798,  789,  779,
  770,   761,   752,   744,   735,   727,   719,  711,  704,  696,  689,  682,
  675,   668,   661,   654,   648,   642,   635,  629,  623,  617,  611,  606,
  600,   595,   589,   584,   579,   574,   569,  564,  559,  554,  550,  545,
  541,   536,   532,   528,   523,   519,   515,  511,  507,  503,  499,  495,
  492,   488,   484,   481,   477,   474,   470,  467,  464,  461,  457,  454,
  451,   448,   445,   442,   439,   436,   433,  430,  427,  425,  422,  419,
  416,   414,   411,   409,   406,   404,   401,  399,  396,  394,  391,  389,
  387,   385,   382,   380,   378,   376,   373,  371,  369,  367,  365,  363,
  361,   359,   357,   355,   353,   351,   349,  348,  346,  344,  342,  340,
  339,   337,   335,   333,   332,   330,   328,  327,  325,  323,  322,  320,
  319,   317,   316,   314,   313,   311,   310,  308,  307,  305,  304,  302,
  301,   300,   298,   297,   296,   294,   293,  292,  290,  289,  288,  286,
  285,   284,   283,   281,   280,   279,   278,  277,  276,  274,  273,  272,
  271,   270,   269,   268,   266,   265,   264,  263,  262,  261,  260,  259,
  258,   257,   256,   255};

//============================================================================
// Convert a floating point value to a rational representation

template<typename T>
static Rational
fromReal(T val, int maxQ)
{
  if (val == T())
    return Rational(0, 1);

  // Find a best rational approximation in the interval ndL <= val <= ndH.
  T ndL[2] = {std::nextafter(val, -std::numeric_limits<T>::infinity()), 1};
  T ndH[2] = {std::nextafter(val, +std::numeric_limits<T>::infinity()), 1};

  // pq is the state used in the calculation of the convergent P/Q using the
  // recurrence formula applied to successive terms of a simple continued
  // fraction.
  int pq[2][2] = {{1, 0}, {0, 1}};

  // Calculate the terms of the continued fractions of ndL and ndH.
  // Calculation proceeds until the terms (aiL, aiH) diverge.
  // The convergent P/Q is updated with each new term until the denominator
  // hits a limit.
  for (int i = 0; i < 10; i++) {
    int aiL = int(ndL[0] / ndL[1]);
    int aiH = int(ndH[0] / ndH[1]);

    int ai = aiL == aiH ? aiL : std::min(aiL, aiH) + 1;

    int p = ai * pq[0][0] + pq[1][0];
    int q = ai * pq[0][1] + pq[1][1];

    if (q > maxQ)
      break;

    pq[1][0] = pq[0][0];
    pq[1][1] = pq[0][1];
    pq[0][0] = p;
    pq[0][1] = q;

    if (aiL != aiH)
      break;

    auto remL = std::fmod(ndL[0], ndL[1]);
    auto remH = std::fmod(ndH[0], ndH[1]);

    ndL[0] = ndL[1];
    ndL[1] = remL;
    ndH[0] = ndH[1];
    ndH[1] = remH;
  }

  return Rational(pq[0][0], pq[0][1]);
}

//----------------------------------------------------------------------------

Rational::Rational(float val)
{
  *this = fromReal(val, 1 << 16);
}

//----------------------------------------------------------------------------

Rational::Rational(double val)
{
  *this = fromReal(val, 1 << 16);
}

//============================================================================

}  // namespace pcc
