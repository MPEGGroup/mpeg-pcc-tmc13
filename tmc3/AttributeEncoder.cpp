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

#include "AttributeEncoder.h"

#include "ArithmeticCodec.h"
#include "DualLutCoder.h"
#include "constants.h"
#include "entropy.h"
#include "quantization.h"
#include "RAHT.h"
#include "FixedPoint.h"

// todo(df): promote to per-attribute encoder parameter
static const double kAttrPredLambdaR = 0.01;
static const double kAttrPredLambdaC = 0.01;

namespace pcc {
//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

struct PCCResidualsEncoder {
  EntropyEncoder arithmeticEncoder;
  StaticBitModel binaryModel0;
  AdaptiveBitModel binaryModelDiff[7];
  AdaptiveBitModel binaryModelIsZero[7];
  AdaptiveBitModel ctxPredMode[2];
  AdaptiveBitModel ctxZeroCnt[3];
  AdaptiveBitModel binaryModelIsOne[7];
  DualLutCoder<false> symbolCoder[2];

  void start(int numPoints);
  int stop();
  void encodePredMode(int value, int max);
  void encodeZeroCnt(int value, int max);
  void encodeSymbol(uint32_t value, int k1, int k2);
  void encode(uint32_t value0, uint32_t value1, uint32_t value2);
  void encode(uint32_t value);
};

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::start(int pointCount)
{
  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = pointCount * 3 * 2 + 1024;
  arithmeticEncoder.setBuffer(maxAcBufLen, nullptr);
  arithmeticEncoder.start();
}

//----------------------------------------------------------------------------

int
PCCResidualsEncoder::stop()
{
  return arithmeticEncoder.stop();
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodePredMode(int mode, int maxMode)
{
  // max = 0 => no direct predictors are used
  if (maxMode == 0)
    return;

  int ctxIdx = 0;
  for (int i = 0; i < mode; i++) {
    arithmeticEncoder.encode(1, ctxPredMode[ctxIdx]);
    ctxIdx = 1;
  }

  // Truncated unary
  if (mode != maxMode)
    arithmeticEncoder.encode(0, ctxPredMode[ctxIdx]);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeZeroCnt(int mode, int maxMode)
{
  // max = 0 => no direct predictors are used
  if (maxMode == 0)
    return;

  int ctxIdx = 0;
  for (int i = 0; i < mode; i++) {
    arithmeticEncoder.encode(1, ctxZeroCnt[ctxIdx]);
    ctxIdx = (ctxIdx == 0 ? 1 : 2);
  }

  // Truncated unary
  if (mode != maxMode)
    arithmeticEncoder.encode(0, ctxZeroCnt[ctxIdx]);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeSymbol(uint32_t value, int k1, int k2)
{
  bool isZero = value == 0;
  arithmeticEncoder.encode(isZero, binaryModelIsZero[k1]);
  if (isZero) {
    return;
  }

  bool isOne = value == 1;
  arithmeticEncoder.encode(isOne, binaryModelIsOne[k1]);
  if (isOne) {
    return;
  }
  value -= 2;

  if (value < kAttributeResidualAlphabetSize) {
    symbolCoder[k2].encode(value, &arithmeticEncoder);
  } else {
    int alphabetSize = kAttributeResidualAlphabetSize;
    symbolCoder[k2].encode(alphabetSize, &arithmeticEncoder);
    arithmeticEncoder.encodeExpGolomb(
      value - alphabetSize, 0, binaryModel0, binaryModelDiff[k1]);
  }
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(uint32_t value0, uint32_t value1, uint32_t value2)
{
  if (value0 == value1 && value0 == value2) {
    value0--;
    value1--;
    value2--;
  }

  int b0 = value0 == 0;
  int b1 = value1 == 0;
  encodeSymbol(value0, 0, 0);
  encodeSymbol(value1, 1 + b0, 1);
  encodeSymbol(value2, 3 + (b0 << 1) + b1, 1);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(uint32_t value)
{
  encodeSymbol(value - 1, 0, 0);
}

//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

struct PCCResidualsEntropyEstimator {
  size_t freq0[kAttributeResidualAlphabetSize + 1];
  size_t freq1[kAttributeResidualAlphabetSize + 1];
  size_t symbolCount0;
  size_t symbolCount1;
  size_t isZero0Count;
  size_t isZero1Count;
  PCCResidualsEntropyEstimator() { init(); }
  void init();
  double bitsDetail(
    const uint32_t detail,
    const size_t symbolCount,
    const size_t* const freq) const;
  double bits(const uint32_t value0) const;
  void update(const uint32_t value0);
  double bits(
    const uint32_t value0, const uint32_t value1, const uint32_t value2) const;
  void
  update(const uint32_t value0, const uint32_t value1, const uint32_t value2);
};

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::init()
{
  for (size_t i = 0; i <= kAttributeResidualAlphabetSize; ++i) {
    freq0[i] = 1;
    freq1[i] = 1;
  }
  symbolCount0 = kAttributeResidualAlphabetSize + 1;
  symbolCount1 = kAttributeResidualAlphabetSize + 1;
  isZero1Count = isZero0Count = symbolCount0 / 2;
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bitsDetail(
  const uint32_t detail,
  const size_t symbolCount,
  const size_t* const freq) const
{
  const uint32_t detailClipped =
    std::min(detail, uint32_t(kAttributeResidualAlphabetSize));
  const double pDetail =
    PCCClip(double(freq[detailClipped]) / symbolCount, 0.001, 0.999);
  double bits = -log2(pDetail);
  if (detail >= kAttributeResidualAlphabetSize) {
    const double x = double(detail) - double(kAttributeResidualAlphabetSize);
    bits += 2.0 * std::floor(log2(x + 1.0)) + 1.0;
  }
  return bits;
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bits(const uint32_t value0) const
{
  const bool isZero0 = value0 == 0;
  const double pIsZero0 = isZero0
    ? double(isZero0Count) / symbolCount0
    : double(symbolCount0 - isZero0Count) / symbolCount0;
  double bits = -log2(PCCClip(pIsZero0, 0.001, 0.999));
  if (!isZero0) {
    bits += bitsDetail(value0 - 1, symbolCount0, freq0);
  }
  return bits;
}

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::update(const uint32_t value0)
{
  const bool isZero0 = value0 == 0;
  ++symbolCount0;
  if (!isZero0) {
    ++freq0[std::min(value0 - 1, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero0Count;
  }
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bits(
  const uint32_t value0, const uint32_t value1, const uint32_t value2) const
{
  const bool isZero0 = value0 == 0;
  const double pIsZero0 = isZero0
    ? double(isZero0Count) / symbolCount0
    : double(symbolCount0 - isZero0Count) / symbolCount0;
  double bits = -log2(PCCClip(pIsZero0, 0.001, 0.999));
  if (!isZero0) {
    bits += bitsDetail(value0 - 1, symbolCount0, freq0);
  }

  const bool isZero1 = value1 == 0 && value2 == 0;
  const double pIsZero1 = isZero1
    ? double(isZero1Count) / symbolCount0
    : double(symbolCount0 - isZero1Count) / symbolCount0;
  bits -= log2(PCCClip(pIsZero1, 0.001, 0.999));
  if (!isZero1) {
    bits += bitsDetail(value1, symbolCount1, freq1);
    bits += bitsDetail(value2, symbolCount1, freq1);
  }
  return bits;
}

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::update(
  const uint32_t value0, const uint32_t value1, const uint32_t value2)
{
  const bool isZero0 = value0 == 0;
  ++symbolCount0;
  if (!isZero0) {
    ++freq0[std::min(value0 - 1, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero0Count;
  }

  const bool isZero1 = value1 == 0 && value2 == 0;
  symbolCount1 += 2;
  if (!isZero1) {
    ++freq1[std::min(value1, uint32_t(kAttributeResidualAlphabetSize))];
    ++freq1[std::min(value2, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero1Count;
  }
}

//============================================================================
// AttributeEncoder Members

void
AttributeEncoder::encode(
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  PCCPointSet3& pointCloud,
  PayloadBuffer* payload)
{
  Quantizers qstep = deriveQuantSteps(attr_aps, abh);

  PCCResidualsEncoder encoder;
  encoder.start(int(pointCloud.getPointCount()));

  if (desc.attr_num_dimensions == 1) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      encodeReflectancesTransformRaht(
        desc, attr_aps, qstep, pointCloud, encoder);
      break;

    case AttributeEncoding::kPredictingTransform:
      encodeReflectancesPred(desc, attr_aps, qstep, pointCloud, encoder);
      break;

    case AttributeEncoding::kLiftingTransform:
      encodeReflectancesLift(desc, attr_aps, qstep, pointCloud, encoder);
      break;
    }
  } else if (desc.attr_num_dimensions == 3) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      encodeColorsTransformRaht(desc, attr_aps, qstep, pointCloud, encoder);
      break;

    case AttributeEncoding::kPredictingTransform:
      encodeColorsPred(desc, attr_aps, qstep, pointCloud, encoder);
      break;

    case AttributeEncoding::kLiftingTransform:
      encodeColorsLift(desc, attr_aps, qstep, pointCloud, encoder);
      break;
    }
  } else {
    assert(desc.attr_num_dimensions == 1 || desc.attr_num_dimensions == 3);
  }

  uint32_t acDataLen = encoder.stop();
  std::copy_n(
    encoder.arithmeticEncoder.buffer(), acDataLen,
    std::back_inserter(*payload));
}

//----------------------------------------------------------------------------

int64_t
AttributeEncoder::computeReflectanceResidual(
  const uint64_t reflectance,
  const uint64_t predictedReflectance,
  const int64_t qs)
{
  const int64_t quantAttValue = reflectance;
  const int64_t quantPredAttValue = predictedReflectance;
  const int64_t delta =
    PCCQuantization(quantAttValue - quantPredAttValue, qs, true);

  return IntToUInt(delta);
}

//----------------------------------------------------------------------------

void
AttributeEncoder::computeReflectancePredictionWeights(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexesLOD,
  const uint32_t predictorIndex,
  PCCPredictor& predictor,
  PCCResidualsEncoder& encoder,
  PCCResidualsEntropyEstimator& context,
  const int64_t qs)
{
  predictor.computeWeights();
  predictor.maxDiff = -1;
  if (predictor.neighborCount > 1) {
    int64_t minValue = 0;
    int64_t maxValue = 0;
    for (size_t i = 0; i < predictor.neighborCount; ++i) {
      const uint64_t reflectanceNeighbor = pointCloud.getReflectance(
        indexesLOD[predictor.neighbors[i].predictorIndex]);
      if (i == 0 || reflectanceNeighbor < minValue) {
        minValue = reflectanceNeighbor;
      }
      if (i == 0 || reflectanceNeighbor > maxValue) {
        maxValue = reflectanceNeighbor;
      }
    }
    const int64_t maxDiff = maxValue - minValue;
    predictor.maxDiff = maxDiff;
    if (maxDiff > aps.adaptive_prediction_threshold) {
      uint64_t attrValue =
        pointCloud.getReflectance(indexesLOD[predictorIndex]);

      // base case: weighted average of n neighbours
      predictor.predMode = 0;
      uint64_t attrPred = predictor.predictReflectance(pointCloud, indexesLOD);
      int64_t attrResidualQuant =
        computeReflectanceResidual(attrValue, attrPred, qs);

      double best_score = attrResidualQuant
        + kAttrPredLambdaR * (qs >> kFixedPointAttributeShift);

      for (int i = 0; i < predictor.neighborCount; i++) {
        if (i == aps.max_num_direct_predictors)
          break;

        attrPred = pointCloud.getReflectance(
          indexesLOD[predictor.neighbors[i].predictorIndex]);
        attrResidualQuant =
          computeReflectanceResidual(attrValue, attrPred, qs);

        double idxBits = i + (i == aps.max_num_direct_predictors - 1 ? 1 : 2);
        double score = attrResidualQuant
          + idxBits * kAttrPredLambdaR * (qs >> kFixedPointAttributeShift);

        if (score < best_score) {
          best_score = score;
          predictor.predMode = i + 1;
          // NB: setting predictor.neighborCount = 1 will cause issues
          // with reconstruction.
        }
      }
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeReflectancesPred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const uint32_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  if (aps.num_detail_levels <= 1) {
    buildPredictorsFastNoLod(
      pointCloud, aps.num_pred_nearest_neighbours, aps.search_range,
      predictors, indexesLOD);
  } else {
    buildPredictorsFast(
      pointCloud, aps.lod_decimation_enabled_flag,
      aps.intra_lod_prediction_enabled_flag, aps.dist2, aps.num_detail_levels,
      aps.num_pred_nearest_neighbours, aps.search_range, aps.search_range,
      predictors, numberOfPointsPerLOD, indexesLOD);
  }
  const int64_t clipMax = (1ll << desc.attr_bitdepth) - 1;
  PCCResidualsEntropyEstimator context;
  int zero_cnt = 0;
  std::vector<int> zerorun;
  zerorun.reserve(pointCount);
  std::vector<uint32_t> residual;
  residual.resize(pointCount);

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    const int64_t qs = qstep[0];
    computeReflectancePredictionWeights(
      aps, pointCloud, indexesLOD, predictorIndex, predictor, encoder, context,
      qs);

    const uint32_t pointIndex = indexesLOD[predictorIndex];
    const uint64_t reflectance = pointCloud.getReflectance(pointIndex);
    const uint16_t predictedReflectance =
      predictor.predictReflectance(pointCloud, indexesLOD);
    const int64_t quantAttValue = reflectance;
    const int64_t quantPredAttValue = predictedReflectance;
    const int64_t delta =
      PCCQuantization(quantAttValue - quantPredAttValue, qs, true);
    const uint32_t attValue0 = uint32_t(IntToUInt(long(delta)));
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs, true);
    const int64_t reconstructedQuantAttValue =
      quantPredAttValue + reconstructedDelta;
    const uint16_t reconstructedReflectance =
      uint16_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));

    if (!attValue0)
      ++zero_cnt;
    else {
      zerorun.push_back(zero_cnt);
      zero_cnt = 0;
    }
    residual[predictorIndex] = attValue0;
    pointCloud.setReflectance(pointIndex, reconstructedReflectance);
  }

  zerorun.push_back(zero_cnt);
  int run_index = 0;
  encoder.encodeZeroCnt(zerorun[run_index], pointCount);
  zero_cnt = zerorun[run_index++];

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    if (predictor.maxDiff > aps.adaptive_prediction_threshold) {
      encoder.encodePredMode(
        predictor.predMode, aps.max_num_direct_predictors);
    }
    if (zero_cnt > 0)
      zero_cnt--;
    else {
      encoder.encode(residual[predictorIndex]);
      if (predictorIndex != pointCount - 1)
        encoder.encodeZeroCnt(zerorun[run_index], pointCount);
      zero_cnt = zerorun[run_index++];
    }
  }
}

//----------------------------------------------------------------------------

Vec3<int64_t>
AttributeEncoder::computeColorResiduals(
  const Vec3<uint8_t> color,
  const Vec3<uint8_t> predictedColor,
  const int64_t qs,
  const int64_t qs2)
{
  Vec3<int64_t> residuals;
  const int64_t quantAttValue = color[0];
  const int64_t quantPredAttValue = predictedColor[0];
  const int64_t delta =
    PCCQuantization(quantAttValue - quantPredAttValue, qs, true);
  residuals[0] = IntToUInt(delta);
  for (size_t k = 1; k < 3; ++k) {
    const int64_t quantAttValue = color[k];
    const int64_t quantPredAttValue = predictedColor[k];
    const int64_t delta =
      PCCQuantization(quantAttValue - quantPredAttValue, qs2, true);
    residuals[k] = IntToUInt(delta);
  }
  return residuals;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::computeColorPredictionWeights(
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexesLOD,
  const uint32_t predictorIndex,
  PCCPredictor& predictor,
  PCCResidualsEncoder& encoder,
  PCCResidualsEntropyEstimator& context,
  const int64_t qs,
  const int64_t qs2)
{
  predictor.computeWeights();
  predictor.maxDiff = -1;
  if (predictor.neighborCount > 1) {
    int64_t minValue[3] = {0, 0, 0};
    int64_t maxValue[3] = {0, 0, 0};
    for (int i = 0; i < predictor.neighborCount; ++i) {
      const Vec3<uint8_t> colorNeighbor =
        pointCloud.getColor(indexesLOD[predictor.neighbors[i].predictorIndex]);
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
    predictor.maxDiff = maxDiff;

    if (maxDiff > aps.adaptive_prediction_threshold) {
      Vec3<uint8_t> attrValue =
        pointCloud.getColor(indexesLOD[predictorIndex]);

      // base case: weighted average of n neighbours
      predictor.predMode = 0;
      Vec3<uint8_t> attrPred = predictor.predictColor(pointCloud, indexesLOD);
      Vec3<int64_t> attrResidualQuant =
        computeColorResiduals(attrValue, attrPred, qs, qs2);

      double best_score = attrResidualQuant[0] + attrResidualQuant[1]
        + attrResidualQuant[2]
        + kAttrPredLambdaC * (double)(qs >> kFixedPointAttributeShift);

      for (int i = 0; i < predictor.neighborCount; i++) {
        if (i == aps.max_num_direct_predictors)
          break;

        attrPred = pointCloud.getColor(
          indexesLOD[predictor.neighbors[i].predictorIndex]);
        attrResidualQuant =
          computeColorResiduals(attrValue, attrPred, qs, qs2);

        double idxBits = i + (i == aps.max_num_direct_predictors - 1 ? 1 : 2);
        double score = attrResidualQuant[0] + attrResidualQuant[1]
          + attrResidualQuant[2]
          + idxBits * kAttrPredLambdaC * (qs >> kFixedPointAttributeShift);

        if (score < best_score) {
          best_score = score;
          predictor.predMode = i + 1;
          // NB: setting predictor.neighborCount = 1 will cause issues
          // with reconstruction.
        }
      }
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeColorsPred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  if (aps.num_detail_levels <= 1) {
    buildPredictorsFastNoLod(
      pointCloud, aps.num_pred_nearest_neighbours, aps.search_range,
      predictors, indexesLOD);
  } else {
    buildPredictorsFast(
      pointCloud, aps.lod_decimation_enabled_flag,
      aps.intra_lod_prediction_enabled_flag, aps.dist2, aps.num_detail_levels,
      aps.num_pred_nearest_neighbours, aps.search_range, aps.search_range,
      predictors, numberOfPointsPerLOD, indexesLOD);
  }
  const int64_t clipMax = (1ll << desc.attr_bitdepth) - 1;
  uint32_t values[3];
  PCCResidualsEntropyEstimator context;

  int zero_cnt = 0;
  std::vector<int> zerorun;
  std::vector<uint32_t> residual[3];
  for (int i = 0; i < 3; i++) {
    residual[i].resize(pointCount);
  }

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    const int64_t qs = qstep[0];
    const int64_t qs2 = qstep[1];
    computeColorPredictionWeights(
      aps, pointCloud, indexesLOD, predictorIndex, predictor, encoder, context,
      qs, qs2);
    const auto pointIndex = indexesLOD[predictorIndex];
    const Vec3<uint8_t> color = pointCloud.getColor(pointIndex);
    const Vec3<uint8_t> predictedColor =
      predictor.predictColor(pointCloud, indexesLOD);
    const int64_t quantAttValue = color[0];
    const int64_t quantPredAttValue = predictedColor[0];
    const int64_t delta =
      PCCQuantization(quantAttValue - quantPredAttValue, qs, true);
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs, true);
    const int64_t reconstructedQuantAttValue =
      quantPredAttValue + reconstructedDelta;
    values[0] = uint32_t(IntToUInt(long(delta)));
    Vec3<uint8_t> reconstructedColor;
    reconstructedColor[0] =
      uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));
    for (size_t k = 1; k < 3; ++k) {
      const int64_t quantAttValue = color[k];
      const int64_t quantPredAttValue = predictedColor[k];
      const int64_t delta =
        PCCQuantization(quantAttValue - quantPredAttValue, qs2, true);
      const int64_t reconstructedDelta =
        PCCInverseQuantization(delta, qs2, true);
      const int64_t reconstructedQuantAttValue =
        quantPredAttValue + reconstructedDelta;
      values[k] = uint32_t(IntToUInt(long(delta)));
      reconstructedColor[k] =
        uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));
    }
    pointCloud.setColor(pointIndex, reconstructedColor);

    if (!values[0] && !values[1] && !values[2]) {
      ++zero_cnt;
    } else {
      zerorun.push_back(zero_cnt);
      zero_cnt = 0;
    }

    for (int i = 0; i < 3; i++) {
      residual[i][predictorIndex] = values[i];
    }
  }

  zerorun.push_back(zero_cnt);
  int run_index = 0;
  encoder.encodeZeroCnt(zerorun[run_index], pointCount);
  zero_cnt = zerorun[run_index++];
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = predictors[predictorIndex];
    if (predictor.maxDiff > aps.adaptive_prediction_threshold) {
      encoder.encodePredMode(
        predictor.predMode, aps.max_num_direct_predictors);
    }
    if (zero_cnt > 0)
      zero_cnt--;
    else {
      for (size_t k = 0; k < 3; k++)
        values[k] = residual[k][predictorIndex];

      encoder.encode(values[0], values[1], values[2]);
      if (predictorIndex != pointCount - 1)
        encoder.encodeZeroCnt(zerorun[run_index], pointCount);
      zero_cnt = zerorun[run_index++];
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeReflectancesTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n], 0);
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Allocate arrays.
  int64_t* mortonCode = new int64_t[voxelCount];
  const int attribCount = 1;
  FixedPoint* attributes = new FixedPoint[attribCount * voxelCount];
  int* integerizedAttributes = new int[attribCount * voxelCount];
  uint64_t* weight = new uint64_t[voxelCount];
  int* binaryLayer = new int[voxelCount];

  // Populate input arrays.
  for (int n = 0; n < voxelCount; n++) {
    weight[n] = 1;
    mortonCode[n] = packedVoxel[n].mortonCode;
    const auto reflectance = pointCloud.getReflectance(packedVoxel[n].index);
    attributes[attribCount * n] = reflectance;
  }

  // Transform.
  regionAdaptiveHierarchicalTransform(
    FixedPoint(qstep[0]), mortonCode, attributes, weight, binaryLayer,
    attribCount, voxelCount, integerizedAttributes);

  // Entropy encode.
  int zero_cnt = 0;
  uint32_t value;
  for (int n = 0; n < voxelCount; ++n) {
    const int64_t detail = IntToUInt(integerizedAttributes[n]);
    assert(detail < std::numeric_limits<uint32_t>::max());
    value = uint32_t(detail);
    if (!value)
      ++zero_cnt;
    else {
      encoder.encodeZeroCnt(zero_cnt, voxelCount);
      encoder.encode(value);
      zero_cnt = 0;
    }
  }
  encoder.encodeZeroCnt(zero_cnt, voxelCount);

  // local decode
  std::fill_n(attributes, attribCount * voxelCount, FixedPoint(0));
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
    weight[n] = 1;
  }

  regionAdaptiveHierarchicalInverseTransform(
    FixedPoint(qstep[0]), mortonCode, attributes, weight, attribCount,
    voxelCount, integerizedAttributes);

  const int64_t maxReflectance = (1 << desc.attr_bitdepth) - 1;
  const int64_t minReflectance = 0;
  for (int n = 0; n < voxelCount; n++) {
    int64_t val = attributes[attribCount * n].round();
    const uint16_t reflectance =
      (uint16_t)PCCClip(val, minReflectance, maxReflectance);
    pointCloud.setReflectance(packedVoxel[n].index, reflectance);
  }

  // De-allocate arrays.
  delete[] binaryLayer;
  delete[] mortonCode;
  delete[] attributes;
  delete[] integerizedAttributes;
  delete[] weight;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeColorsTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n], 0);
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Allocate arrays.
  int64_t* mortonCode = new int64_t[voxelCount];
  const int attribCount = 3;
  FixedPoint* attributes = new FixedPoint[attribCount * voxelCount];
  int* integerizedAttributes = new int[attribCount * voxelCount];
  uint64_t* weight = new uint64_t[voxelCount];
  int* binaryLayer = new int[voxelCount];

  // Populate input arrays.
  for (int n = 0; n < voxelCount; n++) {
    weight[n] = 1;
    mortonCode[n] = packedVoxel[n].mortonCode;
    const auto color = pointCloud.getColor(packedVoxel[n].index);
    attributes[attribCount * n] = color[0];
    attributes[attribCount * n + 1] = color[1];
    attributes[attribCount * n + 2] = color[2];
  }

  // Transform.
  regionAdaptiveHierarchicalTransform(
    FixedPoint(qstep[0]), mortonCode, attributes, weight, binaryLayer,
    attribCount, voxelCount, integerizedAttributes);

  // Entropy encode.
  uint32_t values[attribCount];
  int zero_cnt = 0;
  for (int n = 0; n < voxelCount; ++n) {
    for (int d = 0; d < attribCount; ++d) {
      const int64_t detail =
        IntToUInt(integerizedAttributes[voxelCount * d + n]);
      assert(detail < std::numeric_limits<uint32_t>::max());
      values[d] = uint32_t(detail);
    }
    if (!values[0] && !values[1] && !values[2])
      ++zero_cnt;
    else {
      encoder.encodeZeroCnt(zero_cnt, voxelCount);
      encoder.encode(values[0], values[1], values[2]);
      zero_cnt = 0;
    }
  }

  encoder.encodeZeroCnt(zero_cnt, voxelCount);

  // local decode
  std::fill_n(attributes, attribCount * voxelCount, FixedPoint(0));
  for (int n = 0; n < voxelCount; n++) {
    weight[n] = 1;
    mortonCode[n] = packedVoxel[n].mortonCode;
  }

  regionAdaptiveHierarchicalInverseTransform(
    FixedPoint(qstep[0]), mortonCode, attributes, weight, attribCount,
    voxelCount, integerizedAttributes);

  const int clipMax = (1 << desc.attr_bitdepth) - 1;
  for (int n = 0; n < voxelCount; n++) {
    const int r = attributes[attribCount * n].round();
    const int g = attributes[attribCount * n + 1].round();
    const int b = attributes[attribCount * n + 2].round();
    Vec3<uint8_t> color;
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
  delete[] weight;
}
//----------------------------------------------------------------------------

void
AttributeEncoder::encodeColorsLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
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

  for (size_t index = 0; index < pointCount; ++index) {
    const auto& color = pointCloud.getColor(indexesLOD[index]);
    for (size_t d = 0; d < 3; ++d) {
      colors[index][d] = int32_t(color[d]) << kFixedPointAttributeShift;
    }
  }

  for (size_t i = 0; (i + 1) < lodCount; ++i) {
    const size_t lodIndex = lodCount - i - 1;
    const size_t startIndex = numberOfPointsPerLOD[lodIndex - 1];
    const size_t endIndex = numberOfPointsPerLOD[lodIndex];
    PCCLiftPredict(predictors, startIndex, endIndex, true, colors);
    PCCLiftUpdate(predictors, weights, startIndex, endIndex, true, colors);
  }

  // compress
  int zero_cnt = 0;
  const int64_t qs = qstep[0] << (kFixedPointWeightShift / 2);
  const size_t qs2 = qstep[1] << (kFixedPointWeightShift / 2);
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    const int64_t quantWeight = weights[predictorIndex];
    auto& color = colors[predictorIndex];
    const int64_t delta = PCCQuantization(color[0] * quantWeight, qs);
    const int64_t detail = IntToUInt(delta);
    assert(detail < std::numeric_limits<uint32_t>::max());
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs);
    color[0] = reconstructedDelta / quantWeight;
    uint32_t values[3];
    values[0] = uint32_t(detail);
    for (size_t d = 1; d < 3; ++d) {
      const int64_t delta = PCCQuantization(color[d] * quantWeight, qs2);
      const int64_t detail = IntToUInt(delta);
      assert(detail < std::numeric_limits<uint32_t>::max());
      const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs2);
      color[d] = reconstructedDelta / quantWeight;
      values[d] = uint32_t(detail);
    }
    if (!values[0] && !values[1] && !values[2])
      ++zero_cnt;
    else {
      encoder.encodeZeroCnt(zero_cnt, pointCount);
      encoder.encode(values[0], values[1], values[2]);
      zero_cnt = 0;
    }
  }
  encoder.encodeZeroCnt(zero_cnt, pointCount);

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
      color[d] = uint8_t(PCCClip(color0[d], 0.0, clipMax));
    }
    pointCloud.setColor(indexesLOD[f], color);
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeReflectancesLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Quantizers& qstep,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
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

  for (size_t index = 0; index < pointCount; ++index) {
    reflectances[index] = int32_t(pointCloud.getReflectance(indexesLOD[index]))
      << kFixedPointAttributeShift;
  }

  for (size_t i = 0; (i + 1) < lodCount; ++i) {
    const size_t lodIndex = lodCount - i - 1;
    const size_t startIndex = numberOfPointsPerLOD[lodIndex - 1];
    const size_t endIndex = numberOfPointsPerLOD[lodIndex];
    PCCLiftPredict(predictors, startIndex, endIndex, true, reflectances);
    PCCLiftUpdate(
      predictors, weights, startIndex, endIndex, true, reflectances);
  }

  // compress
  int zero_cnt = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    const int64_t qs = qstep[0] << (kFixedPointWeightShift / 2);
    const int64_t quantWeight = weights[predictorIndex];
    auto& reflectance = reflectances[predictorIndex];
    const int64_t delta = PCCQuantization(reflectance * quantWeight, qs);
    const int64_t detail = IntToUInt(delta);
    assert(detail < std::numeric_limits<uint32_t>::max());
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs);
    reflectance = reconstructedDelta / quantWeight;
    if (!detail)
      ++zero_cnt;
    else {
      encoder.encodeZeroCnt(zero_cnt, pointCount);
      encoder.encode(detail);
      zero_cnt = 0;
    }
  }
  encoder.encodeZeroCnt(zero_cnt, pointCount);

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
    const int64_t refl =
      divExp2RoundHalfInf(reflectances[f], kFixedPointAttributeShift);
    pointCloud.setReflectance(
      indexesLOD[f], uint16_t(PCCClip(refl, int64_t(0), maxReflectance)));
  }
}

//============================================================================

} /* namespace pcc */
