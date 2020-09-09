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
// Mapping of neighbour influence to affected (childIdx + 1).
// NB: zero indicates no influence.

static const int8_t kNeighToChildIdx[26][4] = {
  {1, 0, 0, 0}, {1, 2, 0, 0}, {2, 0, 0, 0}, {1, 3, 0, 0}, {1, 2, 3, 4},
  {2, 4, 0, 0}, {3, 0, 0, 0}, {3, 4, 0, 0}, {4, 0, 0, 0}, {1, 5, 0, 0},
  {1, 2, 5, 6}, {2, 6, 0, 0}, {1, 3, 5, 7}, {2, 4, 6, 8}, {3, 7, 0, 0},
  {3, 4, 7, 8}, {4, 8, 0, 0}, {5, 0, 0, 0}, {5, 6, 0, 0}, {6, 0, 0, 0},
  {5, 7, 0, 0}, {5, 6, 7, 8}, {6, 8, 0, 0}, {7, 0, 0, 0}, {7, 8, 0, 0},
  {8, 0, 0, 0}};

//============================================================================

void
predictGeometryOccupancyIntra(
  const MortonMap3D& occupancyAtlas,
  Vec3<int32_t> pos,
  const int atlasShift,
  int* occupancyIsPredicted,
  int* occupancyPrediction)
{
  uint32_t mask = occupancyAtlas.cubeSize() - 1;
  int32_t x = pos[0] & mask;
  int32_t y = pos[1] & mask;
  int32_t z = pos[2] & mask;

  int score[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  const int8_t* kNeigh = &kNeighToChildIdx[0][0];
  int numOccupied = 0;

  const int shiftX = (atlasShift & 4 ? 1 : 0);
  const int shiftY = (atlasShift & 2 ? 1 : 0);
  const int shiftZ = (atlasShift & 1 ? 1 : 0);

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        if (dz == 0 && dy == 0 && dx == 0)
          continue;

        // todo(df): remove unnecessary checks
        bool occupied = occupancyAtlas.getWithCheck(
          x + dx, y + dy, z + dz, shiftX, shiftY, shiftZ);

        if (occupied) {
          numOccupied++;
          for (int i = 0; i < 4; i++) {
            if (!kNeigh[i])
              break;
            score[kNeigh[i] - 1]++;
          }
        }
        // next pattern
        kNeigh += 4;
      }
    }
  }

  int th0 = 2;
  int th1 = numOccupied < 14 ? 4 : 5;
  int occIsPredicted = 0;
  int occPrediction = 0;

  for (int i = 0; i < 8; i++) {
    if (score[i] <= th0)
      occIsPredicted |= 1 << i;
    else if (score[i] >= th1) {
      occIsPredicted |= 1 << i;
      occPrediction |= 1 << i;
    }
  }

  *occupancyIsPredicted = occIsPredicted;
  *occupancyPrediction = occPrediction;
}

//============================================================================

}  // namespace pcc
