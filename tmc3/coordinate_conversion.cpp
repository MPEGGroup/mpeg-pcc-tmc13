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

#include "coordinate_conversion.h"

#include "geometry_octree.h"

namespace pcc {

//============================================================================

Box3<int>
convertXyzToRpl(
  Vec3<int> laserOrigin,
  const int* laserThetaList,
  int numTheta,
  const Vec3<int>* begin,
  const Vec3<int>* end,
  Vec3<int>* dst)
{
  Box3<int> bbox(INT32_MAX, INT32_MIN);

  for (auto it = begin; it != end; it++, dst++) {
    auto pos = *it - laserOrigin;
    auto laser = findLaser(pos, laserThetaList, numTheta);

    int64_t xLaser = pos[0] << 8;
    int64_t yLaser = pos[1] << 8;
    (*dst)[0] = isqrt(xLaser * xLaser + yLaser * yLaser) >> 8;
    (*dst)[1] = (iatan2(yLaser, xLaser) + 3294199) >> 8;
    (*dst)[2] = laser;

    bbox.insert(*dst);
  }

  return bbox;
}

//----------------------------------------------------------------------------

Vec3<int>
normalisedAxesWeights(Box3<int>& bbox, int forcedMaxLog2)
{
  auto width = bbox.max - bbox.min + 1;
  auto maxWidth = width.max();

  // Otherwise morton code would overflow
  assert(maxWidth < 1 << (21 + 8));

  bool underflow = false;
  if (forcedMaxLog2 > 0) {
    for (int k = 0; k < 3; k++)
      if (width[k] > 1 << (forcedMaxLog2 + 8)) {
        std::cerr << "Warning: normalizedAxesWeight[" << k << "] underflow\n";
        underflow = true;
      }

    while (maxWidth > 1 << (forcedMaxLog2 + 8))
      ++forcedMaxLog2;

    if (underflow)
      std::cerr << "Using " << forcedMaxLog2 << " scaling instead\n";
    maxWidth = 1 << forcedMaxLog2;
  }

  maxWidth = std::min(1 << 21, maxWidth);

  Vec3<int> axesWeight;
  for (int k = 0; k < 3; k++)
    axesWeight[k] = (maxWidth << 8) / width[k];

  return axesWeight;
}

//----------------------------------------------------------------------------

void
offsetAndScale(
  const Vec3<int>& minPos,
  const Vec3<int>& axisWeight,
  Vec3<int>* begin,
  Vec3<int>* end)
{
  for (auto it = begin; it != end; it++)
    *it = times((*it - minPos), axisWeight) + (1 << 7) >> 8;
}

//============================================================================

}  // namespace pcc
