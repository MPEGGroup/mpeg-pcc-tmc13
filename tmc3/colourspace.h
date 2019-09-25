/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
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

#include "PCCMath.h"

namespace pcc {

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformGbrToYCbCrBt709(T<Tv>& gbr)
{
  const Tv g = gbr[0];
  const Tv b = gbr[1];
  const Tv r = gbr[2];
  const double y =
    PCCClip(std::round(0.212600 * r + 0.715200 * g + 0.072200 * b), 0., 255.);
  const double u = PCCClip(
    std::round(-0.114572 * r - 0.385428 * g + 0.5 * b + 128.0), 0., 255.);
  const double v = PCCClip(
    std::round(0.5 * r - 0.454153 * g - 0.045847 * b + 128.0), 0., 255.);
  return {Tv(y), Tv(u), Tv(v)};
}

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformYCbCrBt709ToGbr(T<Tv>& ycbcr)
{
  const double y1 = ycbcr[0];
  const double u1 = ycbcr[1] - 128.0;
  const double v1 = ycbcr[2] - 128.0;
  const double r =
    PCCClip(round(y1 /*- 0.00000 * u1*/ + 1.57480 * v1), 0.0, 255.0);
  const double g =
    PCCClip(round(y1 - 0.18733 * u1 - 0.46813 * v1), 0.0, 255.0);
  const double b =
    PCCClip(round(y1 + 1.85563 * u1 /*+ 0.00000 * v1*/), 0.0, 255.0);
  return {Tv(g), Tv(b), Tv(r)};
}

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformGbrToYCgCoR(int bitDepth, T<Tv>& gbr)
{
  int g = gbr[0];
  int b = gbr[1];
  int r = gbr[2];

  int co = r - b;
  int t = b + (co >> 1);
  int cg = g - t;
  int y = t + (cg >> 1);

  int offset = 1 << bitDepth;

  // NB: YCgCoR needs extra 1-bit for chroma
  return {Tv(y), Tv(cg + offset), Tv(co + offset)};
}

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformYCgCoRToGbr(int bitDepth, T<Tv>& ycgco)
{
  int offset = 1 << bitDepth;
  int y0 = ycgco[0];
  int cg = ycgco[1] - offset;
  int co = ycgco[2] - offset;

  int t = y0 - (cg >> 1);

  int g = cg + t;
  int b = t - (co >> 1);
  int r = co + b;

  int maxVal = (1 << bitDepth) - 1;
  g = PCCClip(g, 0, maxVal);
  b = PCCClip(b, 0, maxVal);
  r = PCCClip(r, 0, maxVal);

  return {Tv(g), Tv(b), Tv(r)};
}

//============================================================================

}  // namespace pcc
