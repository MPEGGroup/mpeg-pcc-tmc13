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
  int sliceQpLuma = attr_aps.init_qp_minus4 + 4;
  int sliceQpChroma = attr_aps.aps_chroma_qp_offset;

  if (attr_aps.aps_slice_qp_deltas_present_flag) {
    sliceQpLuma += abh.attr_qp_delta_luma;
    sliceQpChroma += abh.attr_qp_delta_chroma;
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

QpRegionList
deriveQpRegions(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh)
{
  QpRegionList regions;
  regions.reserve(abh.qpRegions.size());

  for (int i = 0; i < abh.qpRegions.size(); i++) {
    regions.emplace_back();
    auto& region = regions.back();
    const auto& src = abh.qpRegions[i];

    region.qpOffset = src.attr_region_qp_offset;
    region.region.min = src.regionOrigin;
    region.region.max = src.regionOrigin + src.regionSize;
  }

  return regions;
}

//============================================================================

QpSet
deriveQpSet(
  const AttributeDescription& attrDesc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh)
{
  QpSet qpset;
  qpset.layers = deriveLayerQps(attr_aps, abh);
  qpset.regions = deriveQpRegions(attr_aps, abh);

  // The mimimum Qp = 4 is always lossless; the maximum varies according to
  // bitdepth.
  qpset.maxQp = 51 + 6 * (attrDesc.bitdepth - 8);

  // the lifting transform has extra fractional bits that equate to
  // increasing the QP.
  qpset.fixedPointQpOffset = 0;
  if (attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform)
    qpset.fixedPointQpOffset = (kFixedPointWeightShift / 2) * 6;

  return qpset;
}

//============================================================================
// Determines the quantizers at a given layer
Quantizers
QpSet::quantizers(int qpLayer, Qps qpOffset) const
{
  int qp0 = PCCClip(layers[qpLayer][0] + qpOffset[0], 4, maxQp);
  int qp1 = PCCClip(layers[qpLayer][1] + qpOffset[1] + qp0, 4, maxQp);
  qp0 = qp0 + fixedPointQpOffset;
  qp1 = qp1 + fixedPointQpOffset;

  return {Quantizer(qp0), Quantizer(qp1)};
}

//============================================================================
// Determines the quantizers for a point at a given layer
Quantizers
QpSet::quantizers(const Vec3<int32_t>& point, int qpLayer) const
{
  for (const auto& region : regions) {
    if (region.region.contains(point))
      return quantizers(qpLayer, region.qpOffset);
  }

  return quantizers(qpLayer, {0, 0});
}

//============================================================================
//for RAHT region QP Offset
Qps
QpSet::regionQpOffset(const Vec3<int32_t>& point) const
{
  for (const auto& region : regions) {
    if (region.region.contains(point))
      return region.qpOffset;
  }

  return {0, 0};
}

//============================================================================

const int32_t QuantizerGeom::kQpStep[8] = {8, 9, 10, 11, 12, 13, 14, 15};

const int32_t QuantizerGeom::kQpStepRecip[8] = {
  1 << 20, 932068, 838861, 762601, 699051, 645278, 599186, 559241};

//============================================================================

}  // namespace pcc
