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

#include "geometry_intra_pred.h"

namespace pcc {

//============================================================================

static const int8_t LUT_dist[26][8] = {
  {2, 4, 4, 6, 4, 6, 6, 7}, {1, 1, 3, 3, 3, 3, 5, 5}, {4, 2, 6, 4, 6, 4, 7, 6},
  {1, 3, 1, 3, 3, 5, 3, 5}, {0, 0, 0, 0, 2, 2, 2, 2}, {3, 1, 3, 1, 5, 3, 5, 3},
  {4, 6, 2, 4, 6, 7, 4, 6}, {3, 3, 1, 1, 5, 5, 3, 3}, {6, 4, 4, 2, 7, 6, 6, 4},
  {1, 3, 3, 5, 1, 3, 3, 5}, {0, 0, 2, 2, 0, 0, 2, 2}, {3, 1, 5, 3, 3, 1, 5, 3},
  {0, 2, 0, 2, 0, 2, 0, 2}, {2, 0, 2, 0, 2, 0, 2, 0}, {3, 5, 1, 3, 3, 5, 1, 3},
  {2, 2, 0, 0, 2, 2, 0, 0}, {5, 3, 3, 1, 5, 3, 3, 1}, {4, 6, 6, 7, 2, 4, 4, 6},
  {3, 3, 5, 5, 1, 1, 3, 3}, {6, 4, 7, 6, 4, 2, 6, 4}, {3, 5, 3, 5, 1, 3, 1, 3},
  {2, 2, 2, 2, 0, 0, 0, 0}, {5, 3, 5, 3, 3, 1, 3, 1}, {6, 7, 4, 6, 4, 6, 2, 4},
  {5, 5, 3, 3, 3, 3, 1, 1}, {7, 6, 6, 4, 6, 4, 4, 2}};

//----------------------------------------------------------------------------

static const int LUT1[8] = {27, 39, 20, 8, 18, 4, 11, 18};
static const int LUT0[8] = {-1, -6, 12, 20, 14, 28, 22, 12};
static const int LUT_th0[5] = {62, 60, 61, 59, 59};
static const int LUT_th1[5] = {67, 66, 65, 66, 64};

//============================================================================

void
predictGeometryOccupancyIntra(
  const MortonMap3D& occupancyAtlas,
  PCCVector3<uint32_t> pos,
  int nodeSizeLog2,
  int* occupancyIsPredicted,
  int* occupancyPrediction)
{
  uint32_t mask = occupancyAtlas.cubeSize() - 1;
  int32_t x = (pos[0] >> nodeSizeLog2) & mask;
  int32_t y = (pos[1] >> nodeSizeLog2) & mask;
  int32_t z = (pos[2] >> nodeSizeLog2) & mask;

  int score[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int numOccupied = 0;

  const int8_t* p = &LUT_dist[0][0];

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        if (dz == 0 && dy == 0 && dx == 0)
          continue;

        // todo(df): remove unnecessary checks
        bool occupied = occupancyAtlas.getWithCheck(x + dx, y + dy, z + dz);

        if (occupied) {
          for (int i = 0; i < 8; i++)
            score[i] += LUT1[p[i]];
          numOccupied++;
        } else {
          for (int i = 0; i < 8; i++)
            score[i] += LUT0[p[i]];
        }

        // next pattern
        p += 8;
      }
    }
  }

  if (numOccupied <= 8) {
    *occupancyIsPredicted = 0;
    *occupancyPrediction = 0;
    return;
  }

  numOccupied -= 9;
  if (numOccupied >= 4)
    numOccupied = 4;

  int th0 = LUT_th0[numOccupied] * 26;
  int th1 = LUT_th1[numOccupied] * 26;

  int occIsPredicted = 0;
  int occPrediction = 0;

  for (int i = 0; i < 8; i++) {
    int score_i = score[i] << 2;
    if (score_i <= th0) {
      occIsPredicted |= 1 << i;
    } else if (score_i >= th1) {
      occIsPredicted |= 1 << i;
      occPrediction |= 1 << i;
    }
  }

  *occupancyIsPredicted = occIsPredicted;
  *occupancyPrediction = occPrediction;
}

//============================================================================

}  // namespace pcc
