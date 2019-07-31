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

#include "quantization.h"

#include "constants.h"
#include "hls.h"
#include "tables.h"

namespace pcc {

//============================================================================

Quantizer::Quantizer(int qp)
{
  int qpShift = qp / 6;
  _stepSize = kQpStep[qp % 6] << qpShift;
  _stepSizeRecip = kQpStepRecip[qp % 6] >> qpShift;
}

//============================================================================

Quantizers
deriveQuantizers(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh)
{
  int sliceQpLuma = attr_aps.init_qp;
  int sliceQpChroma = attr_aps.init_qp + attr_aps.aps_chroma_qp_offset;

  if (attr_aps.aps_slice_qp_deltas_present_flag) {
    sliceQpLuma += abh.attr_qp_delta_luma;
    sliceQpChroma += abh.attr_qp_delta_chroma;
  }

  // the lifting transform has extra fractional bits that equate to
  // increasing the QP.
  if (attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform) {
    int fixedPointQpOffset = (kFixedPointWeightShift / 2) * 6;
    sliceQpLuma += fixedPointQpOffset;
    sliceQpChroma += fixedPointQpOffset;
  }

  return {Quantizer{sliceQpLuma}, Quantizer{sliceQpChroma}};
}

//============================================================================

}  // namespace pcc