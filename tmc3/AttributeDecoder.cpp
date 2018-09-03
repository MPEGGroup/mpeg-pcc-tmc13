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

namespace pcc {

//============================================================================
// An encapsulation of the entropy decoding methods used in attribute coding

struct PCCResidualsDecoder {
  EntropyDecoder arithmeticDecoder;
  StaticBitModel binaryModel0;
  AdaptiveBitModel binaryModelPred;
  AdaptiveBitModel binaryModelDiff[7];
  AdaptiveBitModel binaryModelIsZero[7];
  DualLutCoder<false> symbolCoder[2];

  void start(const char* buf, int buf_len);
  void stop();
  bool decodePred();
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

bool
PCCResidualsDecoder::decodePred()
{
  return arithmeticDecoder.decode(binaryModelPred);
}

//----------------------------------------------------------------------------

uint32_t
PCCResidualsDecoder::decodeSymbol(int k1, int k2)
{
  if (arithmeticDecoder.decode(binaryModelIsZero[k1])) {
    return 0u;
  }

  uint32_t value = symbolCoder[k2].decode(&arithmeticDecoder);
  if (value == kAttributeResidualAlphabetSize) {
    value +=
      arithmeticDecoder.decodeExpGolomb(0, binaryModel0, binaryModelDiff[k1]);
  }
  ++value;

  return value;
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
}

//----------------------------------------------------------------------------

uint32_t
PCCResidualsDecoder::decode()
{
  return decodeSymbol(0, 0);
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
  /* AttributeBrickHeader abh = */ parseAbh(payload, &abhSize);

  PCCResidualsDecoder decoder;
  decoder.start(payload.data() + abhSize, payload.size() - abhSize);

  if (attr_desc.attr_num_dimensions == 1) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeReflectancesRaht(attr_desc, attr_aps, decoder, pointCloud);
      break;

    case AttributeEncoding::kPredictingTransform:
      decodeReflectancesPred(attr_desc, attr_aps, decoder, pointCloud);
      break;

    case AttributeEncoding::kLiftingTransform:
      decodeReflectancesLift(attr_desc, attr_aps, decoder, pointCloud);
      break;
    }
  } else if (attr_desc.attr_num_dimensions == 3) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      decodeColorsRaht(attr_desc, attr_aps, decoder, pointCloud);
      break;

    case AttributeEncoding::kPredictingTransform:
      decodeColorsPred(attr_desc, attr_aps, decoder, pointCloud);
      break;

    case AttributeEncoding::kLiftingTransform:
      decodeColorsLift(attr_desc, attr_aps, decoder, pointCloud);
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
  PCCPredictor& predictor,
  PCCResidualsDecoder& decoder)
{
  predictor.computeWeights();
  if (predictor.neighborCount > 1) {
    int64_t minValue = 0;
    int64_t maxValue = 0;
    for (int i = 0; i < aps.num_pred_nearest_neighbours; ++i) {
      const uint16_t reflectanceNeighbor =
        pointCloud.getReflectance(predictor.neighbors[i].index);
      if (i == 0 || reflectanceNeighbor < minValue) {
        minValue = reflectanceNeighbor;
      }
      if (i == 0 || reflectanceNeighbor > maxValue) {
        maxValue = reflectanceNeighbor;
      }
    }
    const int64_t maxDiff = maxValue - minValue;
    if (maxDiff > aps.adaptive_prediction_threshold) {
      const bool predIndex = decoder.decodePred();
      if (predIndex == 0) {
        predictor.neighborCount = 1;
        predictor.neighbors[0].weight = 1.0;
      }
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesPred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;
  PCCBuildLevelOfDetail(
    pointCloud, aps.numDetailLevels, aps.dist2, numberOfPointsPerLOD,
    indexesLOD);
  PCCComputePredictors2(
    pointCloud, numberOfPointsPerLOD, indexesLOD,
    aps.num_pred_nearest_neighbours, predictors);
  const size_t pointCount = pointCloud.getPointCount();
  const int64_t maxReflectance = (1ll << desc.attr_bitdepth) - 1;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    const int64_t qs = aps.quant_step_size_luma;
    computeReflectancePredictionWeights(aps, pointCloud, predictor, decoder);
    uint16_t& reflectance = pointCloud.getReflectance(predictor.index);

    const uint32_t attValue0 = decoder.decode();
    const int64_t quantPredAttValue = predictor.predictReflectance(pointCloud);
    const int64_t delta = PCCInverseQuantization(UIntToInt(attValue0), qs);
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
  PCCPredictor& predictor,
  PCCResidualsDecoder& decoder)
{
  predictor.computeWeights();
  if (predictor.neighborCount > 1) {
    int64_t minValue[3] = {0, 0, 0};
    int64_t maxValue[3] = {0, 0, 0};
    for (int i = 0; i < aps.num_pred_nearest_neighbours; ++i) {
      const PCCColor3B colorNeighbor =
        pointCloud.getColor(predictor.neighbors[i].index);
      for (size_t k = 0; k < 3; ++k) {
        if (i == 0 || colorNeighbor[k] < minValue[k]) {
          minValue[k] = colorNeighbor[k];
        }
        if (i == 0 || colorNeighbor[k] > maxValue[k]) {
          maxValue[k] = colorNeighbor[k];
        }
      }
    }
    const int64_t maxDiff = (std::max)(
      maxValue[2] - minValue[2],
      (std::max)(maxValue[0] - minValue[0], maxValue[1] - minValue[1]));
    if (maxDiff > aps.adaptive_prediction_threshold) {
      const bool predIndex = decoder.decodePred();
      if (predIndex == 0) {
        predictor.neighborCount = 1;
        predictor.neighbors[0].weight = 1.0;
      }
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeColorsPred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;
  PCCBuildLevelOfDetail(
    pointCloud, aps.numDetailLevels, aps.dist2, numberOfPointsPerLOD,
    indexesLOD);
  PCCComputePredictors2(
    pointCloud, numberOfPointsPerLOD, indexesLOD,
    aps.num_pred_nearest_neighbours, predictors);
  const size_t pointCount = pointCloud.getPointCount();
  uint32_t values[3];
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    const int64_t qs = aps.quant_step_size_luma;
    const int64_t qs2 = aps.quant_step_size_chroma;
    computeColorPredictionWeights(aps, pointCloud, predictor, decoder);
    decoder.decode(values);
    PCCColor3B& color = pointCloud.getColor(predictor.index);
    const PCCColor3B predictedColor = predictor.predictColor(pointCloud);
    const int64_t quantPredAttValue = predictedColor[0];
    const int64_t delta = PCCInverseQuantization(UIntToInt(values[0]), qs);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;

    int64_t clipMax = (1 << desc.attr_bitdepth) - 1;
    color[0] =
      uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));
    for (size_t k = 1; k < 3; ++k) {
      const int64_t quantPredAttValue = predictedColor[k];
      const int64_t delta = PCCInverseQuantization(UIntToInt(values[k]), qs2);
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
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    const auto position = pointCloud[n];
    int x = int(position[0]);
    int y = int(position[1]);
    int z = int(position[2]);
    long long mortonCode = 0;
    for (int b = 0; b < aps.raht_depth; b++) {
      mortonCode |= (long long)((x >> b) & 1) << (3 * b + 2);
      mortonCode |= (long long)((y >> b) & 1) << (3 * b + 1);
      mortonCode |= (long long)((z >> b) & 1) << (3 * b);
    }
    packedVoxel[n].mortonCode = mortonCode;
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Moroton codes
  long long* mortonCode = new long long[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
  }

  // Re-obtain weights at the decoder by calling RAHT without any attributes.
  float* weight = new float[voxelCount];
  int* binaryLayer = new int[voxelCount];
  regionAdaptiveHierarchicalTransform(
    mortonCode, nullptr, weight, binaryLayer, 0, voxelCount, aps.raht_depth);

  // Sort integerized attributes by weight
  std::vector<WeightWithIndex> sortedWeight(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    sortedWeight[n].weight = weight[n];
    sortedWeight[n].index = n;
  }
  sort(sortedWeight.begin(), sortedWeight.end());

  // Entropy decode
  int* sortedIntegerizedAttributes = new int[voxelCount];
  for (int n = 0; n < voxelCount; ++n) {
    const uint32_t attValue0 = decoder.decode();
    sortedIntegerizedAttributes[n] = UIntToInt(attValue0);
  }

  // Unsort integerized attributes by weight.
  int* integerizedAttributes = new int[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    integerizedAttributes[sortedWeight[n].index] =
      sortedIntegerizedAttributes[n];
  }

  // Inverse Quantize.
  float* attributes = new float[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    attributes[n] = integerizedAttributes[n] * aps.quant_step_size_luma;
  }

  regionAdaptiveHierarchicalInverseTransform(
    mortonCode, attributes, 1, voxelCount, aps.raht_depth);

  const int maxReflectance = (1 << desc.attr_bitdepth) - 1;
  const int minReflectance = 0;
  for (int n = 0; n < voxelCount; n++) {
    const int reflectance =
      PCCClip((int)round(attributes[n]), minReflectance, maxReflectance);
    pointCloud.setReflectance(packedVoxel[n].index, uint16_t(reflectance));
  }

  // De-allocate arrays.
  delete[] binaryLayer;
  delete[] mortonCode;
  delete[] attributes;
  delete[] integerizedAttributes;
  delete[] sortedIntegerizedAttributes;
  delete[] weight;
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeColorsRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    const auto position = pointCloud[n];
    int x = int(position[0]);
    int y = int(position[1]);
    int z = int(position[2]);
    long long mortonCode = 0;
    for (int b = 0; b < aps.raht_depth; b++) {
      mortonCode |= (long long)((x >> b) & 1) << (3 * b + 2);
      mortonCode |= (long long)((y >> b) & 1) << (3 * b + 1);
      mortonCode |= (long long)((z >> b) & 1) << (3 * b);
    }
    packedVoxel[n].mortonCode = mortonCode;
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Moroton codes
  long long* mortonCode = new long long[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
  }

  // Re-obtain weights at the decoder by calling RAHT without any attributes.
  float* weight = new float[voxelCount];
  int* binaryLayer = new int[voxelCount];
  regionAdaptiveHierarchicalTransform(
    mortonCode, nullptr, weight, binaryLayer, 0, voxelCount, aps.raht_depth);

  // Sort integerized attributes by weight
  std::vector<WeightWithIndex> sortedWeight(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    sortedWeight[n].weight = weight[n];
    sortedWeight[n].index = n;
  }
  sort(sortedWeight.begin(), sortedWeight.end());

  // Entropy decode
  const int attribCount = 3;
  uint32_t values[3];
  int* sortedIntegerizedAttributes = new int[attribCount * voxelCount];
  for (int n = 0; n < voxelCount; ++n) {
    if (
      binaryLayer[sortedWeight[n].index] >= aps.raht_binary_level_threshold) {
      decoder.decode(values);
      sortedIntegerizedAttributes[n] = UIntToInt(values[0]);
      for (int d = 1; d < 3; ++d) {
        sortedIntegerizedAttributes[voxelCount * d + n] = UIntToInt(values[d]);
      }
    } else {
      values[0] = decoder.decode();
      sortedIntegerizedAttributes[n] = UIntToInt(values[0]);
      for (int d = 1; d < 3; d++) {
        sortedIntegerizedAttributes[voxelCount * d + n] = 0;
      }
    }
  }

  // Unsort integerized attributes by weight.
  int* integerizedAttributes = new int[attribCount * voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    for (int k = 0; k < attribCount; k++) {
      // Pull sorted integerized attributes out of column-major order.
      integerizedAttributes[attribCount * sortedWeight[n].index + k] =
        sortedIntegerizedAttributes[voxelCount * k + n];
    }
  }

  // Inverse Quantize.
  float* attributes = new float[attribCount * voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    for (int k = 0; k < attribCount; k++) {
      attributes[attribCount * n + k] =
        integerizedAttributes[attribCount * n + k] * aps.quant_step_size_luma;
    }
  }

  regionAdaptiveHierarchicalInverseTransform(
    mortonCode, attributes, attribCount, voxelCount, aps.raht_depth);

  const int clipMax = (1 << desc.attr_bitdepth) - 1;
  for (int n = 0; n < voxelCount; n++) {
    const int r = (int)round(attributes[attribCount * n]);
    const int g = (int)round(attributes[attribCount * n + 1]);
    const int b = (int)round(attributes[attribCount * n + 2]);
    PCCColor3B color;
    color[0] = uint8_t(PCCClip(r, 0, clipMax));
    color[1] = uint8_t(PCCClip(g, 0, clipMax));
    color[2] = uint8_t(PCCClip(b, 0, clipMax));
    pointCloud.setColor(packedVoxel[n].index, color);
  }

  // De-allocate arrays.
  delete[] binaryLayer;
  delete[] mortonCode;
  delete[] attributes;
  delete[] integerizedAttributes;
  delete[] sortedIntegerizedAttributes;
  delete[] weight;
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeColorsLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;
  PCCBuildLevelOfDetail(
    pointCloud, aps.numDetailLevels, aps.dist2, numberOfPointsPerLOD,
    indexesLOD);
  const size_t lodCount = numberOfPointsPerLOD.size();
  std::vector<PCCPredictor> predictors;
  PCCComputePredictors(
    pointCloud, numberOfPointsPerLOD, indexesLOD,
    aps.num_pred_nearest_neighbours, predictors);
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    predictors[predictorIndex].computeWeights();
  }
  std::vector<PCCVector3D> colors;
  colors.resize(pointCount);
  std::vector<double> weights;
  PCCComputeQuantizationWeights(predictors, weights);
  // decompress
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    uint32_t values[3];
    decoder.decode(values);
    const int64_t qs = aps.quant_step_size_luma;
    const int64_t qs2 = aps.quant_step_size_chroma;
    const double quantWeight = sqrt(weights[predictorIndex]);
    auto& color = colors[predictorIndex];
    const int64_t delta = UIntToInt(values[0]);
    const double reconstructedDelta = PCCInverseQuantization(delta, qs);
    color[0] = reconstructedDelta / quantWeight;
    for (size_t d = 1; d < 3; ++d) {
      const int64_t delta = UIntToInt(values[d]);
      const double reconstructedDelta = PCCInverseQuantization(delta, qs2);
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

  const double clipMax = (1 << desc.attr_bitdepth) - 1;
  for (size_t f = 0; f < pointCount; ++f) {
    const auto& predictor = predictors[f];
    PCCColor3B color;
    for (size_t d = 0; d < 3; ++d) {
      color[d] = uint8_t(PCCClip(std::round(colors[f][d]), 0.0, clipMax));
    }
    pointCloud.setColor(predictor.index, color);
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;
  PCCBuildLevelOfDetail(
    pointCloud, aps.numDetailLevels, aps.dist2, numberOfPointsPerLOD,
    indexesLOD);
  const size_t lodCount = numberOfPointsPerLOD.size();
  std::vector<PCCPredictor> predictors;
  PCCComputePredictors(
    pointCloud, numberOfPointsPerLOD, indexesLOD,
    aps.num_pred_nearest_neighbours, predictors);
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    predictors[predictorIndex].computeWeights();
  }
  std::vector<double> reflectances;
  reflectances.resize(pointCount);
  std::vector<double> weights;
  PCCComputeQuantizationWeights(predictors, weights);

  // decompress
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    const int64_t detail = decoder.decode();
    const int64_t qs = aps.quant_step_size_luma;
    const double quantWeight = sqrt(weights[predictorIndex]);
    auto& reflectance = reflectances[predictorIndex];
    const int64_t delta = UIntToInt(detail);
    const double reconstructedDelta = PCCInverseQuantization(delta, qs);
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
  const double maxReflectance = (1 << desc.attr_bitdepth) - 1;
  for (size_t f = 0; f < pointCount; ++f) {
    pointCloud.setReflectance(
      predictors[f].index,
      uint16_t(PCCClip(std::round(reflectances[f]), 0.0, maxReflectance)));
  }
}

//============================================================================

} /* namespace pcc */
