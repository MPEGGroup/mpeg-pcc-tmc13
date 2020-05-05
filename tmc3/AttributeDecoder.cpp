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

#include "AttributeDecoder.h"

#include "DualLutCoder.h"
#include "constants.h"
#include "entropy.h"
#include "io_hls.h"
#include "RAHT.h"
#include "FixedPoint.h"

namespace pcc {

//============================================================================
// An encapsulation of the entropy decoding methods used in attribute coding

struct PCCResidualsDecoder {
  EntropyDecoder arithmeticDecoder;
  StaticBitModel binaryModel0;
  AdaptiveBitModel binaryModelDiff[7];
  AdaptiveBitModel binaryModelIsZero[7];
  AdaptiveBitModel ctxPredMode[2];
  AdaptiveBitModel ctxZeroCnt[3];
  AdaptiveBitModel binaryModelIsOne[7];
  AdaptiveBitModel ctxSymbolBit[2];
  AdaptiveBitModel ctxSetIdx[2][16];

  void start(const SequenceParameterSet& sps, const char* buf, int buf_len);
  void stop();
  int decodePredMode(int max);
  int decodeZeroCnt(int max);
  uint32_t decodeSymbol(int k1, int k2, int k3);
  int decodeInterval(int k3);
  void decode(int32_t values[3]);
  int32_t decode();
};

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::start(
  const SequenceParameterSet& sps, const char* buf, int buf_len)
{
  arithmeticDecoder.setBuffer(buf_len, buf);
  arithmeticDecoder.enableBypassStream(sps.cabac_bypass_stream_enabled_flag);
  arithmeticDecoder.start();
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::stop()
{
  arithmeticDecoder.stop();
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodePredMode(int maxMode)
{
  int mode = 0;

  if (maxMode == 0)
    return mode;

  int ctxIdx = 0;
  while (arithmeticDecoder.decode(ctxPredMode[ctxIdx])) {
    ctxIdx = 1;
    mode++;
    if (mode == maxMode)
      break;
  }

  return mode;
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeZeroCnt(int maxMode)
{
  int mode = 0;

  if (maxMode == 0)
    return mode;

  int ctxIdx = 0;
  while (arithmeticDecoder.decode(ctxZeroCnt[ctxIdx])) {
    ctxIdx = (ctxIdx == 0 ? 1 : 2);
    mode++;
    if (mode == maxMode)
      break;
  }
  return mode;
}

//----------------------------------------------------------------------------

uint32_t
PCCResidualsDecoder::decodeSymbol(int k1, int k2, int k3)
{
  if (arithmeticDecoder.decode(binaryModelIsZero[k1]))
    return 0u;

  if (arithmeticDecoder.decode(binaryModelIsOne[k2]))
    return 1u;

  uint32_t value = decodeInterval(k3);
  if (value == kAttributeResidualAlphabetSize) {
    value +=
      arithmeticDecoder.decodeExpGolomb(0, binaryModel0, binaryModelDiff[k1]);
  }

  return value + 2;
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeInterval(int k3)
{
  // Decoding of interval index
  auto& aed = arithmeticDecoder;
  int setIdx = 0;
  setIdx = (setIdx << 1) | aed.decode(ctxSetIdx[k3][0]);
  setIdx = (setIdx << 1) | aed.decode(ctxSetIdx[k3][1 + setIdx]);
  setIdx = (setIdx << 1) | aed.decode(ctxSetIdx[k3][3 + setIdx]);
  setIdx = (setIdx << 1) | aed.decode(ctxSetIdx[k3][7 + setIdx]);

  // Decode position within interval
  int intervalStart = kCoeffIntervalStart[setIdx];
  int intervalEnd = kCoeffIntervalStart[setIdx + 1];
  int intervalRange = intervalEnd - intervalStart;

  // Number of bits to code is log2(intervalEnd - intervalStart)
  // The following assumes that the range is a power of two
  int symbolIdx = 0;
  for (int mask = intervalRange - 1, i = 0; mask; mask >>= 1, ++i)
    symbolIdx |= aed.decode(ctxSymbolBit[k3]) << i;

  // Reconstruct
  return intervalStart + symbolIdx;
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::decode(int32_t value[3])
{
  value[0] = decodeSymbol(0, 0, 0);
  int b0 = value[0] == 0;
  int b1 = value[0] <= 1;
  value[1] = decodeSymbol(1 + b0, 1 + b1, 1);
  int b2 = value[1] == 0;
  int b3 = value[1] <= 1;
  value[2] = decodeSymbol(3 + (b0 << 1) + b2, 3 + (b1 << 1) + b3, 1);

  if (value[0] && arithmeticDecoder.decode(binaryModel0))
    value[0] = -value[0];
  if (value[1] && arithmeticDecoder.decode(binaryModel0))
    value[1] = -value[1];
  if (value[2] && arithmeticDecoder.decode(binaryModel0))
    value[2] = -value[2];
}

//----------------------------------------------------------------------------

int32_t
PCCResidualsDecoder::decode()
{
  auto mag = decodeSymbol(0, 0, 0) + 1;
  bool sign = arithmeticDecoder.decode(binaryModel0);
  return sign ? -mag : mag;
}

//============================================================================
// AttributeDecoderIntf

AttributeDecoderIntf::~AttributeDecoderIntf() = default;

//============================================================================
// AttributeDecoder factory

std::unique_ptr<AttributeDecoderIntf>
makeAttributeDecoder()
{
  return std::unique_ptr<AttributeDecoder>(new AttributeDecoder());
}

//============================================================================
// AttributeDecoder Members

void
AttributeDecoder::decode(
  const SequenceParameterSet& sps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  const PayloadBuffer& payload,
  PCCPointSet3& pointCloud)
{
  int abhSize;
  AttributeBrickHeader abh = parseAbh(sps, attr_aps, payload, &abhSize);

  QpSet qpSet = deriveQpSet(attr_desc, attr_aps, abh);

  PCCResidualsDecoder decoder;
  decoder.start(sps, payload.data() + abhSize, payload.size() - abhSize);

  // generate LoDs if necessary
  if (attr_aps.lodParametersPresent() && _lods.empty())
    _lods.generate(
      attr_aps, geom_num_points_minus1, minGeomNodeSizeLog2, pointCloud);

  if (attr_desc.attr_num_dimensions_minus1 == 0) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeReflectancesRaht(attr_desc, attr_aps, qpSet, decoder, pointCloud);
      break;

    case AttributeEncoding::kPredictingTransform:
      decodeReflectancesPred(attr_desc, attr_aps, qpSet, decoder, pointCloud);
      break;

    case AttributeEncoding::kLiftingTransform:
      decodeReflectancesLift(
        attr_desc, attr_aps, qpSet, geom_num_points_minus1,
        minGeomNodeSizeLog2, decoder, pointCloud);
      break;
    }
  } else if (attr_desc.attr_num_dimensions_minus1 == 2) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeColorsRaht(attr_desc, attr_aps, qpSet, decoder, pointCloud);
      break;

    case AttributeEncoding::kPredictingTransform:
      decodeColorsPred(attr_desc, attr_aps, qpSet, decoder, pointCloud);
      break;

    case AttributeEncoding::kLiftingTransform:
      decodeColorsLift(
        attr_desc, attr_aps, qpSet, geom_num_points_minus1,
        minGeomNodeSizeLog2, decoder, pointCloud);
      break;
    }
  } else {
    assert(
      attr_desc.attr_num_dimensions_minus1 == 0
      || attr_desc.attr_num_dimensions_minus1 == 2);
  }

  decoder.stop();
}

//----------------------------------------------------------------------------

bool
AttributeDecoder::isReusable(const AttributeParameterSet& aps) const
{
  return _lods.isReusable(aps);
}

//----------------------------------------------------------------------------

void
AttributeDecoder::computeReflectancePredictionWeights(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexes,
  PCCPredictor& predictor,
  PCCResidualsDecoder& decoder)
{
  predictor.predMode = 0;
  int64_t maxDiff = 0;

  if (predictor.neighborCount > 1 && aps.max_num_direct_predictors) {
    int64_t minValue = 0;
    int64_t maxValue = 0;
    for (int i = 0; i < predictor.neighborCount; ++i) {
      const attr_t reflectanceNeighbor = pointCloud.getReflectance(
        indexes[predictor.neighbors[i].predictorIndex]);
      if (i == 0 || reflectanceNeighbor < minValue) {
        minValue = reflectanceNeighbor;
      }
      if (i == 0 || reflectanceNeighbor > maxValue) {
        maxValue = reflectanceNeighbor;
      }
    }
    maxDiff = maxValue - minValue;
  }

  if (maxDiff >= aps.adaptive_prediction_threshold) {
    predictor.predMode = decoder.decodePredMode(aps.max_num_direct_predictors);
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesPred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  const int64_t maxReflectance = (1ll << desc.bitdepth) - 1;
  int zero_cnt = decoder.decodeZeroCnt(pointCount);
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }
    const uint32_t pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);
    auto& predictor = _lods.predictors[predictorIndex];

    computeReflectancePredictionWeights(
      aps, pointCloud, _lods.indexes, predictor, decoder);
    attr_t& reflectance = pointCloud.getReflectance(pointIndex);
    int32_t attValue0 = 0;
    if (zero_cnt > 0) {
      zero_cnt--;
    } else {
      attValue0 = decoder.decode();
      zero_cnt = decoder.decodeZeroCnt(pointCount);
    }
    const int64_t quantPredAttValue =
      predictor.predictReflectance(pointCloud, _lods.indexes);
    const int64_t delta =
      divExp2RoundHalfUp(quant[0].scale(attValue0), kFixedPointAttributeShift);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
    reflectance =
      attr_t(PCCClip(reconstructedQuantAttValue, int64_t(0), maxReflectance));
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::computeColorPredictionWeights(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexes,
  PCCPredictor& predictor,
  PCCResidualsDecoder& decoder)
{
  int64_t maxDiff = 0;

  if (predictor.neighborCount > 1 && aps.max_num_direct_predictors) {
    int64_t minValue[3] = {0, 0, 0};
    int64_t maxValue[3] = {0, 0, 0};
    for (int i = 0; i < predictor.neighborCount; ++i) {
      const Vec3<attr_t> colorNeighbor =
        pointCloud.getColor(indexes[predictor.neighbors[i].predictorIndex]);
      for (size_t k = 0; k < 3; ++k) {
        if (i == 0 || colorNeighbor[k] < minValue[k]) {
          minValue[k] = colorNeighbor[k];
        }
        if (i == 0 || colorNeighbor[k] > maxValue[k]) {
          maxValue[k] = colorNeighbor[k];
        }
      }
    }
    maxDiff = (std::max)(
      maxValue[2] - minValue[2],
      (std::max)(maxValue[0] - minValue[0], maxValue[1] - minValue[1]));
  }

  if (maxDiff >= aps.adaptive_prediction_threshold) {
    predictor.predMode = decoder.decodePredMode(aps.max_num_direct_predictors);
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeColorsPred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();

  Vec3<int64_t> clipMax{(1 << desc.bitdepth) - 1,
                        (1 << desc.bitdepthSecondary) - 1,
                        (1 << desc.bitdepthSecondary) - 1};

  int32_t values[3];
  int zero_cnt = decoder.decodeZeroCnt(pointCount);
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }
    const uint32_t pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);
    auto& predictor = _lods.predictors[predictorIndex];

    computeColorPredictionWeights(
      aps, pointCloud, _lods.indexes, predictor, decoder);
    if (zero_cnt > 0) {
      values[0] = values[1] = values[2] = 0;
      zero_cnt--;
    } else {
      decoder.decode(values);
      zero_cnt = decoder.decodeZeroCnt(pointCount);
    }
    Vec3<attr_t>& color = pointCloud.getColor(pointIndex);
    const Vec3<attr_t> predictedColor =
      predictor.predictColor(pointCloud, _lods.indexes);

    int64_t residual0 = 0;
    for (int k = 0; k < 3; ++k) {
      const auto& q = quant[std::min(k, 1)];
      const int64_t residual =
        divExp2RoundHalfUp(q.scale(values[k]), kFixedPointAttributeShift);
      const int64_t recon = predictedColor[k] + residual + residual0;
      color[k] = attr_t(PCCClip(recon, int64_t(0), clipMax[k]));

      if (!k && aps.inter_component_prediction_enabled_flag)
        residual0 = residual;
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n]);
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Morton codes
  int64_t* mortonCode = new int64_t[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
  }

  // Entropy decode
  const int attribCount = 1;
  int* coefficients = new int[attribCount * voxelCount];
  Qps* pointQpOffsets = new Qps[voxelCount];
  int zero_cnt = decoder.decodeZeroCnt(voxelCount);
  for (int n = 0; n < voxelCount; ++n) {
    uint32_t value = 0;
    if (zero_cnt > 0) {
      zero_cnt--;
    } else {
      value = decoder.decode();
      zero_cnt = decoder.decodeZeroCnt(voxelCount);
    }
    coefficients[n] = value;
    pointQpOffsets[n] = qpSet.regionQpOffset(pointCloud[packedVoxel[n].index]);
  }

  int* attributes = new int[attribCount * voxelCount];
  const int rahtPredThreshold[2] = {aps.raht_prediction_threshold0,
                                    aps.raht_prediction_threshold1};

  regionAdaptiveHierarchicalInverseTransform(
    aps.raht_prediction_enabled_flag, rahtPredThreshold, qpSet, pointQpOffsets,
    mortonCode, attributes, attribCount, voxelCount, coefficients);

  const int64_t maxReflectance = (1 << desc.bitdepth) - 1;
  const int64_t minReflectance = 0;
  for (int n = 0; n < voxelCount; n++) {
    int64_t val = attributes[attribCount * n];
    const attr_t reflectance =
      attr_t(PCCClip(val, minReflectance, maxReflectance));
    pointCloud.setReflectance(packedVoxel[n].index, reflectance);
  }

  // De-allocate arrays.
  delete[] mortonCode;
  delete[] attributes;
  delete[] coefficients;
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeColorsRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n]);
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Morton codes
  int64_t* mortonCode = new int64_t[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
  }

  // Entropy decode
  const int attribCount = 3;
  int zero_cnt = decoder.decodeZeroCnt(voxelCount);
  int* coefficients = new int[attribCount * voxelCount];
  Qps* pointQpOffsets = new Qps[voxelCount];

  for (int n = 0; n < voxelCount; ++n) {
    int32_t values[3];
    if (zero_cnt > 0) {
      values[0] = values[1] = values[2] = 0;
      zero_cnt--;
    } else {
      decoder.decode(values);
      zero_cnt = decoder.decodeZeroCnt(voxelCount);
    }
    for (int d = 0; d < attribCount; ++d) {
      coefficients[voxelCount * d + n] = values[d];
    }
    pointQpOffsets[n] = qpSet.regionQpOffset(pointCloud[packedVoxel[n].index]);
  }

  int* attributes = new int[attribCount * voxelCount];
  const int rahtPredThreshold[2] = {aps.raht_prediction_threshold0,
                                    aps.raht_prediction_threshold1};

  regionAdaptiveHierarchicalInverseTransform(
    aps.raht_prediction_enabled_flag, rahtPredThreshold, qpSet, pointQpOffsets,
    mortonCode, attributes, attribCount, voxelCount, coefficients);

  Vec3<int> clipMax{(1 << desc.bitdepth) - 1,
                    (1 << desc.bitdepthSecondary) - 1,
                    (1 << desc.bitdepthSecondary) - 1};

  for (int n = 0; n < voxelCount; n++) {
    const int r = attributes[attribCount * n];
    const int g = attributes[attribCount * n + 1];
    const int b = attributes[attribCount * n + 2];
    Vec3<attr_t> color;
    color[0] = attr_t(PCCClip(r, 0, clipMax[0]));
    color[1] = attr_t(PCCClip(g, 0, clipMax[1]));
    color[2] = attr_t(PCCClip(b, 0, clipMax[2]));
    pointCloud.setColor(packedVoxel[n].index, color);
  }

  // De-allocate arrays.
  delete[] mortonCode;
  delete[] attributes;
  delete[] coefficients;
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeColorsLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint64_t> weights;

  if (!aps.scalable_lifting_enabled_flag) {
    PCCComputeQuantizationWeights(_lods.predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      _lods.predictors, _lods.numPointsInLod, geom_num_points_minus1 + 1,
      minGeomNodeSizeLog2, weights);
  }

  const size_t lodCount = _lods.numPointsInLod.size();
  std::vector<Vec3<int64_t>> colors;
  colors.resize(pointCount);

  // NB: when partially decoding, the truncated unary limit for zero_run
  // must be the original value.  geom_num_points may be the case.  However,
  // the encoder does lie sometimes, and there are actually more points.
  int zeroCntLimit = std::max(geom_num_points_minus1 + 1, int(pointCount));

  // decompress
  int zero_cnt = decoder.decodeZeroCnt(zeroCntLimit);
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }
    const uint32_t pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);

    int32_t values[3];
    if (zero_cnt > 0) {
      values[0] = values[1] = values[2] = 0;
      zero_cnt--;
    } else {
      decoder.decode(values);
      zero_cnt = decoder.decodeZeroCnt(zeroCntLimit);
    }

    const int64_t iQuantWeight = irsqrt(weights[predictorIndex]);
    auto& color = colors[predictorIndex];
    const int64_t delta = values[0];
    const int64_t reconstructedDelta = quant[0].scale(delta);
    color[0] = divExp2RoundHalfInf(reconstructedDelta * iQuantWeight, 40);
    for (size_t d = 1; d < 3; ++d) {
      const int64_t delta = values[d];
      const int64_t reconstructedDelta = quant[1].scale(delta);
      color[d] = divExp2RoundHalfInf(reconstructedDelta * iQuantWeight, 40);
    }
  }

  // reconstruct
  for (size_t lodIndex = 1; lodIndex < lodCount; ++lodIndex) {
    const size_t startIndex = _lods.numPointsInLod[lodIndex - 1];
    const size_t endIndex = _lods.numPointsInLod[lodIndex];
    PCCLiftUpdate(
      _lods.predictors, weights, startIndex, endIndex, false, colors);
    PCCLiftPredict(_lods.predictors, startIndex, endIndex, false, colors);
  }

  Vec3<int64_t> clipMax{(1 << desc.bitdepth) - 1,
                        (1 << desc.bitdepthSecondary) - 1,
                        (1 << desc.bitdepthSecondary) - 1};

  for (size_t f = 0; f < pointCount; ++f) {
    const auto color0 =
      divExp2RoundHalfInf(colors[f], kFixedPointAttributeShift);
    Vec3<attr_t> color;
    for (size_t d = 0; d < 3; ++d) {
      color[d] = attr_t(PCCClip(color0[d], int64_t(0), clipMax[d]));
    }
    pointCloud.setColor(_lods.indexes[f], color);
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint64_t> weights;

  if (!aps.scalable_lifting_enabled_flag) {
    PCCComputeQuantizationWeights(_lods.predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      _lods.predictors, _lods.numPointsInLod, geom_num_points_minus1 + 1,
      minGeomNodeSizeLog2, weights);
  }

  const size_t lodCount = _lods.numPointsInLod.size();
  std::vector<int64_t> reflectances;
  reflectances.resize(pointCount);

  // NB: when partially decoding, the truncated unary limit for zero_run
  // must be the original value.  geom_num_points may be the case.  However,
  // the encoder does lie sometimes, and there are actually more points.
  int zeroCntLimit = std::max(geom_num_points_minus1 + 1, int(pointCount));

  // decompress
  int zero_cnt = decoder.decodeZeroCnt(zeroCntLimit);
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }
    const uint32_t pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);

    int64_t detail = 0;
    if (zero_cnt > 0) {
      zero_cnt--;
    } else {
      detail = decoder.decode();
      zero_cnt = decoder.decodeZeroCnt(zeroCntLimit);
    }
    const int64_t iQuantWeight = irsqrt(weights[predictorIndex]);
    auto& reflectance = reflectances[predictorIndex];
    const int64_t delta = detail;
    const int64_t reconstructedDelta = quant[0].scale(delta);
    reflectance = divExp2RoundHalfInf(reconstructedDelta * iQuantWeight, 40);
  }

  // reconstruct
  for (size_t lodIndex = 1; lodIndex < lodCount; ++lodIndex) {
    const size_t startIndex = _lods.numPointsInLod[lodIndex - 1];
    const size_t endIndex = _lods.numPointsInLod[lodIndex];
    PCCLiftUpdate(
      _lods.predictors, weights, startIndex, endIndex, false, reflectances);
    PCCLiftPredict(
      _lods.predictors, startIndex, endIndex, false, reflectances);
  }
  const int64_t maxReflectance = (1 << desc.bitdepth) - 1;
  for (size_t f = 0; f < pointCount; ++f) {
    const auto refl =
      divExp2RoundHalfInf(reflectances[f], kFixedPointAttributeShift);
    pointCloud.setReflectance(
      _lods.indexes[f], attr_t(PCCClip(refl, int64_t(0), maxReflectance)));
  }
}

//============================================================================

} /* namespace pcc */
