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
  DualLutCoder<false> symbolCoder[2];

  void start(const char* buf, int buf_len);
  void stop();
  int decodePredMode(int max);
  int decodeZeroCnt(int max);
  uint32_t decodeSymbol(int k1, int k2);
  void decode(uint32_t values[3]);
  uint32_t decode();
};

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::start(const char* buf, int buf_len)
{
  arithmeticDecoder.setBuffer(buf_len, buf);
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
PCCResidualsDecoder::decodeSymbol(int k1, int k2)
{
  if (arithmeticDecoder.decode(binaryModelIsZero[k1]))
    return 0u;

  if (arithmeticDecoder.decode(binaryModelIsOne[k1]))
    return 1u;

  uint32_t value = symbolCoder[k2].decode(&arithmeticDecoder);
  if (value == kAttributeResidualAlphabetSize) {
    value +=
      arithmeticDecoder.decodeExpGolomb(0, binaryModel0, binaryModelDiff[k1]);
  }

  return value + 2;
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::decode(uint32_t value[3])
{
  value[0] = decodeSymbol(0, 0);
  int b0 = value[0] == 0;
  value[1] = decodeSymbol(1 + b0, 1);
  int b1 = value[1] == 0;
  value[2] = decodeSymbol(3 + (b0 << 1) + b1, 1);

  int d = (value[0] == value[1] && value[0] == value[2]);
  for (int k = 0; k < 3; k++) {
    value[k] += d;
  }
}

//----------------------------------------------------------------------------

uint32_t
PCCResidualsDecoder::decode()
{
  return decodeSymbol(0, 0) + 1;
}

//============================================================================
// AttributeDecoder Members

void
AttributeDecoder::decode(
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const PayloadBuffer& payload,
  PCCPointSet3& pointCloud)
{
  int abhSize;
  AttributeBrickHeader abh = parseAbh(attr_aps, payload, &abhSize);
  Quantizers qstep = deriveQuantSteps(attr_aps, abh);

  PCCResidualsDecoder decoder;
  decoder.start(payload.data() + abhSize, payload.size() - abhSize);

  if (attr_desc.attr_num_dimensions == 1) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeReflectancesRaht(attr_desc, attr_aps, qstep, decoder, pointCloud);
      break;

    case AttributeEncoding::kPredictingTransform:
      decodeReflectancesPred(attr_desc, attr_aps, qstep, decoder, pointCloud);
      break;

    case AttributeEncoding::kLiftingTransform:
      decodeReflectancesLift(attr_desc, attr_aps, qstep, decoder, pointCloud);
      break;
    }
  } else if (attr_desc.attr_num_dimensions == 3) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeColorsRaht(attr_desc, attr_aps, qstep, decoder, pointCloud);
      break;

    case AttributeEncoding::kPredictingTransform:
      decodeColorsPred(attr_desc, attr_aps, qstep, decoder, pointCloud);
      break;

    case AttributeEncoding::kLiftingTransform:
      decodeColorsLift(attr_desc, attr_aps, qstep, decoder, pointCloud);
      break;
    }
  } else {
    assert(
      attr_desc.attr_num_dimensions == 1
      || attr_desc.attr_num_dimensions == 3);
  }

  decoder.stop();
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
  predictor.computeWeights();
  int64_t maxDiff = 0;

  if (predictor.neighborCount > 1) {
    int64_t minValue = 0;
    int64_t maxValue = 0;
    for (int i = 0; i < predictor.neighborCount; ++i) {
      const uint16_t reflectanceNeighbor = pointCloud.getReflectance(
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
  const Quantizers& qstep,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    pointCloud, aps.lod_decimation_enabled_flag,
    aps.intra_lod_prediction_enabled_flag, aps.dist2, aps.num_detail_levels,
    aps.num_pred_nearest_neighbours, aps.search_range, aps.search_range,
    predictors, numberOfPointsPerLOD, indexesLOD);

  const int64_t maxReflectance = (1ll << desc.attr_bitdepth) - 1;
  int zero_cnt = decoder.decodeZeroCnt(pointCount);
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    const int64_t qs = qstep[0];
    computeReflectancePredictionWeights(
      aps, pointCloud, indexesLOD, predictor, decoder);
    const uint32_t pointIndex = indexesLOD[predictorIndex];
    uint16_t& reflectance = pointCloud.getReflectance(pointIndex);
    uint32_t attValue0 = 0;
    if (zero_cnt > 0) {
      zero_cnt--;
    } else {
      attValue0 = decoder.decode();
      zero_cnt = decoder.decodeZeroCnt(pointCount);
    }
    const int64_t quantPredAttValue =
      predictor.predictReflectance(pointCloud, indexesLOD);
    const int64_t delta =
      PCCInverseQuantization(UIntToInt(attValue0), qs, true);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
    reflectance = uint16_t(
      PCCClip(reconstructedQuantAttValue, int64_t(0), maxReflectance));
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
  predictor.computeWeights();
  int64_t maxDiff = 0;

  if (predictor.neighborCount > 1) {
    int64_t minValue[3] = {0, 0, 0};
    int64_t maxValue[3] = {0, 0, 0};
    for (int i = 0; i < predictor.neighborCount; ++i) {
      const Vec3<uint8_t> colorNeighbor =
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
  const Quantizers& qstep,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    pointCloud, aps.lod_decimation_enabled_flag,
    aps.intra_lod_prediction_enabled_flag, aps.dist2, aps.num_detail_levels,
    aps.num_pred_nearest_neighbours, aps.search_range, aps.search_range,
    predictors, numberOfPointsPerLOD, indexesLOD);

  uint32_t values[3];
  int zero_cnt = decoder.decodeZeroCnt(pointCount);
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    const int64_t qs = qstep[0];
    const int64_t qs2 = qstep[1];
    computeColorPredictionWeights(
      aps, pointCloud, indexesLOD, predictor, decoder);
    if (zero_cnt > 0) {
      values[0] = values[1] = values[2] = 0;
      zero_cnt--;
    } else {
      decoder.decode(values);
      zero_cnt = decoder.decodeZeroCnt(pointCount);
    }
    const uint32_t pointIndex = indexesLOD[predictorIndex];
    Vec3<uint8_t>& color = pointCloud.getColor(pointIndex);
    const Vec3<uint8_t> predictedColor =
      predictor.predictColor(pointCloud, indexesLOD);
    const int64_t quantPredAttValue = predictedColor[0];
    const int64_t delta =
      PCCInverseQuantization(UIntToInt(values[0]), qs, true);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
    int64_t clipMax = (1 << desc.attr_bitdepth) - 1;
    color[0] =
      uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));
    for (size_t k = 1; k < 3; ++k) {
      const int64_t quantPredAttValue = predictedColor[k];
      const int64_t delta =
        PCCInverseQuantization(UIntToInt(values[k]), qs2, true);
      const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
      color[k] =
        uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n], 0);
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
  int zero_cnt = decoder.decodeZeroCnt(voxelCount);
  for (int n = 0; n < voxelCount; ++n) {
    uint32_t value = 0;
    if (zero_cnt > 0) {
      zero_cnt--;
    } else {
      value = decoder.decode();
      zero_cnt = decoder.decodeZeroCnt(voxelCount);
    }
    coefficients[n] = UIntToInt(value);
  }

  int* attributes = new int[attribCount * voxelCount];

  regionAdaptiveHierarchicalInverseTransform(
    aps.raht_prediction_enabled_flag, qstep, mortonCode, attributes,
    attribCount, voxelCount, coefficients);

  const int64_t maxReflectance = (1 << desc.attr_bitdepth) - 1;
  const int64_t minReflectance = 0;
  for (int n = 0; n < voxelCount; n++) {
    int64_t val = attributes[attribCount * n];
    const uint16_t reflectance =
      (uint16_t)PCCClip(val, minReflectance, maxReflectance);
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
  const Quantizers& qstep,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n], 0);
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

  for (int n = 0; n < voxelCount; ++n) {
    uint32_t values[3];
    if (zero_cnt > 0) {
      values[0] = values[1] = values[2] = 0;
      zero_cnt--;
    } else {
      decoder.decode(values);
      zero_cnt = decoder.decodeZeroCnt(voxelCount);
    }
    for (int d = 0; d < attribCount; ++d) {
      coefficients[voxelCount * d + n] = UIntToInt(values[d]);
    }
  }

  int* attributes = new int[attribCount * voxelCount];

  regionAdaptiveHierarchicalInverseTransform(
    aps.raht_prediction_enabled_flag, qstep, mortonCode, attributes,
    attribCount, voxelCount, coefficients);

  const int clipMax = (1 << desc.attr_bitdepth) - 1;
  for (int n = 0; n < voxelCount; n++) {
    const int r = attributes[attribCount * n];
    const int g = attributes[attribCount * n + 1];
    const int b = attributes[attribCount * n + 2];
    Vec3<uint8_t> color;
    color[0] = uint8_t(PCCClip(r, 0, clipMax));
    color[1] = uint8_t(PCCClip(g, 0, clipMax));
    color[2] = uint8_t(PCCClip(b, 0, clipMax));
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
  const Quantizers& qstep,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    pointCloud, aps.lod_decimation_enabled_flag,
    aps.intra_lod_prediction_enabled_flag, aps.dist2, aps.num_detail_levels,
    aps.num_pred_nearest_neighbours, aps.search_range, aps.search_range,
    predictors, numberOfPointsPerLOD, indexesLOD);

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    predictors[predictorIndex].computeWeights();
  }
  std::vector<uint64_t> weights;
  PCCComputeQuantizationWeights(predictors, weights);
  const size_t lodCount = numberOfPointsPerLOD.size();
  std::vector<Vec3<int64_t>> colors;
  colors.resize(pointCount);
  // decompress
  int zero_cnt = decoder.decodeZeroCnt(pointCount);
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    uint32_t values[3];
    if (zero_cnt > 0) {
      values[0] = values[1] = values[2] = 0;
      zero_cnt--;
    } else {
      decoder.decode(values);
      zero_cnt = decoder.decodeZeroCnt(pointCount);
    }
    const int64_t qs = qstep[0] << (kFixedPointWeightShift / 2);
    const int64_t qs2 = qstep[1] << (kFixedPointWeightShift / 2);
    // + kFixedPointAttributeShift ???
    const int64_t quantWeight = weights[predictorIndex];
    auto& color = colors[predictorIndex];
    const int64_t delta = UIntToInt(values[0]);
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs);
    color[0] = reconstructedDelta / quantWeight;
    for (size_t d = 1; d < 3; ++d) {
      const int64_t delta = UIntToInt(values[d]);
      const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs2);
      color[d] = reconstructedDelta / quantWeight;
    }
  }

  // reconstruct
  for (size_t lodIndex = 1; lodIndex < lodCount; ++lodIndex) {
    const size_t startIndex = numberOfPointsPerLOD[lodIndex - 1];
    const size_t endIndex = numberOfPointsPerLOD[lodIndex];
    PCCLiftUpdate(predictors, weights, startIndex, endIndex, false, colors);
    PCCLiftPredict(predictors, startIndex, endIndex, false, colors);
  }

  const int64_t clipMax = (1 << desc.attr_bitdepth) - 1;
  for (size_t f = 0; f < pointCount; ++f) {
    const auto color0 =
      divExp2RoundHalfInf(colors[f], kFixedPointAttributeShift);
    Vec3<uint8_t> color;
    for (size_t d = 0; d < 3; ++d) {
      color[d] = uint8_t(PCCClip(color0[d], int64_t(0), clipMax));
    }
    pointCloud.setColor(indexesLOD[f], color);
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    pointCloud, aps.lod_decimation_enabled_flag,
    aps.intra_lod_prediction_enabled_flag, aps.dist2, aps.num_detail_levels,
    aps.num_pred_nearest_neighbours, aps.search_range, aps.search_range,
    predictors, numberOfPointsPerLOD, indexesLOD);

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    predictors[predictorIndex].computeWeights();
  }

  std::vector<uint64_t> weights;
  PCCComputeQuantizationWeights(predictors, weights);
  const size_t lodCount = numberOfPointsPerLOD.size();
  std::vector<int64_t> reflectances;
  reflectances.resize(pointCount);

  // decompress
  int zero_cnt = decoder.decodeZeroCnt(pointCount);
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    int64_t detail = 0;
    if (zero_cnt > 0) {
      zero_cnt--;
    } else {
      detail = decoder.decode();
      zero_cnt = decoder.decodeZeroCnt(pointCount);
    }
    const int64_t qs = qstep[0] << (kFixedPointWeightShift / 2);
    const int64_t quantWeight = weights[predictorIndex];
    auto& reflectance = reflectances[predictorIndex];
    const int64_t delta = UIntToInt(detail);
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs);
    reflectance = reconstructedDelta / quantWeight;
  }

  // reconstruct
  for (size_t lodIndex = 1; lodIndex < lodCount; ++lodIndex) {
    const size_t startIndex = numberOfPointsPerLOD[lodIndex - 1];
    const size_t endIndex = numberOfPointsPerLOD[lodIndex];
    PCCLiftUpdate(
      predictors, weights, startIndex, endIndex, false, reflectances);
    PCCLiftPredict(predictors, startIndex, endIndex, false, reflectances);
  }
  const int64_t maxReflectance = (1 << desc.attr_bitdepth) - 1;
  for (size_t f = 0; f < pointCount; ++f) {
    const auto refl =
      divExp2RoundHalfInf(reflectances[f], kFixedPointAttributeShift);
    pointCloud.setReflectance(
      indexesLOD[f], uint16_t(PCCClip(refl, int64_t(0), maxReflectance)));
  }
}

//============================================================================

} /* namespace pcc */
