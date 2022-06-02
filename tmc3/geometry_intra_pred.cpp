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

void

construct26NeighbourWord(

  const MortonMap3D& occupancyAtlas,
  Vec3<int32_t> pos,
  const int atlasShift,
  int Word4[8]
)
{
  uint32_t mask = occupancyAtlas.cubeSize() - 1;
  int32_t x = pos[0] & mask;
  int32_t y = pos[1] & mask;
  int32_t z = pos[2] & mask;

  const int shiftX = (atlasShift & 4 ? 1 : 0);
  const int shiftY = (atlasShift & 2 ? 1 : 0);
  const int shiftZ = (atlasShift & 1 ? 1 : 0);
  int WordDiag[8] = { 0,0,0,0,0,0,0,0 };
  static const int LUTdx[20] = { -1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1 };
  static const int LUTdy[20] = { -1,-1,-1,0,0,1,1,1,-1,-1,1,1,-1,-1,-1,0,0,1,1,1 };
  static const int LUTdz[20] = { -1,0,1,-1,1,-1,0,1,-1,1,-1,1,-1,0,1,-1,1,-1,0,1 };
  static const bool LUTsumdydydzIs3[20] = { 1,0,1,0,0,1,0,1,0,0,0,0,1,0,1,0,0,1,0,1 };
  static const int LUTi2[12][2] = { {0,1}, {0,2}, {1,3}, {2,3}, {0,4}, {1,5}, {2,6}, {3,7}, {4,5}, {4,6}, {5,7}, {6,7} };

  int i3 = 0;
  const int *i2 = LUTi2[0];
  for (int n = 0; n < 20; n++) {
    int occupied = occupancyAtlas.getWithCheck(x + LUTdx[n], y +LUTdy[n], z + LUTdz[n], shiftX, shiftY, shiftZ);

    if (LUTsumdydydzIs3[n])
    {
      WordDiag[i3++] |= occupied;
    }
    else {
      Word4[*i2] <<= 1;
      Word4[*i2++] |= occupied;
      Word4[*i2] <<= 1;
      Word4[*i2++] |= occupied;
    }
  }

  for (int i = 0; i < 8; i++) {
    Word4[i] <<= 1;
    Word4[i] |= WordDiag[i];
  }

}

//============================================================================

}  // namespace pcc
