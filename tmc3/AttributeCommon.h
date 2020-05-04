/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2019, ISO/IEC
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

#include <stdint.h>
#include <vector>

#include "hls.h"
#include "PCCTMC3Common.h"

namespace pcc {

//============================================================================
// Constants for coefficient coding

constexpr int kCoeffIntervals = 16;
constexpr int kCoeffIntervalMax = 255;

// The starting coefficient value for each interval
extern const uint16_t kCoeffIntervalStart[kCoeffIntervals + 1];

// Map of coefficient values to intervals
extern const uint8_t kCoeffToIntervalIdx[kCoeffIntervalMax + 1];

//============================================================================

struct AttributeLods {
  // Indicates if the generated LoDs are compatible with the provided aps
  bool isReusable(const AttributeParameterSet& aps) const;

  bool empty() const { return numPointsInLod.empty(); };

  void generate(
    const AttributeParameterSet& aps,
    int geom_num_points_minus1,
    int minGeomNodeSizeLog2,
    const PCCPointSet3& cloud);

  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numPointsInLod;
  std::vector<uint32_t> indexes;

private:
  // This is the aps that was used to generate the LoDs.  It is used to check
  // if the generated LoDs are reusable.
  AttributeParameterSet _aps;
};

//============================================================================

}  // namespace pcc
