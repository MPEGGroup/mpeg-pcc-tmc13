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
  qp = std::max(qp, 4);
  int qpShift = qp / 6;
  _stepSize = kQpStep[qp % 6] << qpShift;
  _stepSizeRecip = kQpStepRecip[qp % 6] >> qpShift;
}

//============================================================================

Qps
deriveQps(
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  int qpLayer)
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

  if (abh.attr_layer_qp_present_flag()) {
    sliceQpLuma += abh.attr_layer_qp_delta_luma[qpLayer];
    sliceQpChroma += abh.attr_layer_qp_delta_chroma[qpLayer];
  }

  return {sliceQpLuma, sliceQpChroma};
}

//============================================================================

QpLayers
deriveLayerQps(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh)
{
  QpLayers layers;

  layers.push_back(deriveQps(attr_aps, abh, 0));

  if (abh.attr_layer_qp_present_flag()) {
    int numLayers = abh.attr_num_qp_layers_minus1() + 1;
    for (int layer = 1; layer < numLayers; layer++) {
      layers.push_back(deriveQps(attr_aps, abh, layer));
    }
  }

  return layers;
}

//============================================================================

QpRegionOffset
deriveQpRegions(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh)
{
  QpRegionOffset qpRegionOffset;
  qpRegionOffset.valid = abh.attr_region_qp_present_flag;
  if (qpRegionOffset.valid) {
    qpRegionOffset.qpOffset = abh.attr_region_qp_delta;
    for (int i = 0; i < 3; i++) {
      qpRegionOffset.region.min[i] = abh.attr_region_qp_origin[i];
      qpRegionOffset.region.max[i] =
        abh.attr_region_qp_origin[i] + abh.attr_region_qp_whd[i];
    }
  }
  return qpRegionOffset;
}

//============================================================================

QpSet
deriveQpSet(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh)
{
  return {deriveLayerQps(attr_aps, abh), deriveQpRegions(attr_aps, abh)};
}

//============================================================================
//for PredLift
Quantizers
QpSet::quantizers(const Vec3<int32_t>& point, int qpLayer) const
{
  int qpRegionOffset = 0;
  if (regionOffset.valid && regionOffset.region.contains(point)) {
    qpRegionOffset = regionOffset.qpOffset;
  }
  return {Quantizer(layers[qpLayer][0] + qpRegionOffset),
          Quantizer(layers[qpLayer][1] + qpRegionOffset)};
}

//============================================================================
//for RAHT region QP Offset
int
QpSet::regionQpOffset(const Vec3<int32_t>& point) const
{
  int qpRegionOffset = 0;
  if (regionOffset.valid && regionOffset.region.contains(point)) {
    qpRegionOffset = regionOffset.qpOffset;
  }
  return qpRegionOffset;
}

//============================================================================

const int32_t QuantizerGeom::kQpStep[6] = {659445, 741374,  831472,
                                           933892, 1048576, 1175576};

const int32_t QuantizerGeom::kQpStepRecip[6] = {1667328, 1483072, 1322368,
                                                1177344, 1048576, 935296};

//============================================================================

}  // namespace pcc
