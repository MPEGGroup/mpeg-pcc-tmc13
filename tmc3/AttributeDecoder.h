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

#include <stdint.h>

#include "PayloadBuffer.h"
#include "PCCTMC3Common.h"

namespace pcc {

//============================================================================
// Opaque definitions (Internal detail)

struct PCCResidualsDecoder;

//============================================================================

class AttributeDecoder {
public:
  void decode(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    const PayloadBuffer&,
    PCCPointSet3& pointCloud);

protected:
  // todo(df): consider alternative encapsulation

  void decodeReflectancesLift(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud);

  void decodeColorsLift(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud);

  void decodeReflectancesPred(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud);

  void decodeColorsPred(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud);

  void decodeReflectancesRaht(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud);

  void decodeColorsRaht(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud);

  static void computeColorPredictionWeights(
    const AttributeParameterSet& aps,
    const PCCPointSet3& pointCloud,
    const std::vector<uint32_t>& indexes,
    PCCPredictor& predictor,
    PCCResidualsDecoder& decoder);

  static void computeReflectancePredictionWeights(
    const AttributeParameterSet& aps,
    const PCCPointSet3& pointCloud,
    const std::vector<uint32_t>& indexes,
    PCCPredictor& predictor,
    PCCResidualsDecoder& decoder);
};

//============================================================================

} /* namespace pcc */
