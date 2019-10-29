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

  void start(const SequenceParameterSet& sps, int numPoints);
  int stop();
  void encodePredMode(int value, int max);
  void encodeZeroCnt(int value, int max);
  void encodeSymbol(uint32_t value, int k1, int k2);
  void encode(uint32_t value0, uint32_t value1, uint32_t value2);
  void encode(uint32_t value);
};

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::start(const SequenceParameterSet& sps, int pointCount)
{
  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = pointCount * 3 * 2 + 1024;
  arithmeticEncoder.setBuffer(maxAcBufLen, nullptr);
  arithmeticEncoder.enableBypassStream(sps.cabac_bypass_stream_enabled_flag);
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
  const SequenceParameterSet& sps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  PCCPointSet3& pointCloud,
  PayloadBuffer* payload)
{
  std::vector<Quantizers> quantLayers = deriveQuantizerLayers(attr_aps, abh);

  PCCResidualsEncoder encoder;
  encoder.start(sps, int(pointCloud.getPointCount()));

  if (desc.attr_num_dimensions == 1) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      encodeReflectancesTransformRaht(
        desc, attr_aps, quantLayers, pointCloud, encoder);
      break;

    case AttributeEncoding::kPredictingTransform:
      encodeReflectancesPred(desc, attr_aps, quantLayers, pointCloud, encoder);
      break;

    case AttributeEncoding::kLiftingTransform:
      encodeReflectancesLift(desc, attr_aps, quantLayers, pointCloud, encoder);
      break;
    }
  } else if (desc.attr_num_dimensions == 3) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      encodeColorsTransformRaht(
        desc, attr_aps, quantLayers, pointCloud, encoder);
      break;

    case AttributeEncoding::kPredictingTransform:
      encodeColorsPred(desc, attr_aps, quantLayers, pointCloud, encoder);
      break;

    case AttributeEncoding::kLiftingTransform:
      encodeColorsLift(desc, attr_aps, quantLayers, pointCloud, encoder);
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
  const Quantizer& quant)
{
  const int64_t quantAttValue = reflectance;
  const int64_t quantPredAttValue = predictedReflectance;
  const int64_t delta = quant.quantize(
    (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);

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
  const Quantizer& quant)
{
  predictor.computeWeights();
  predictor.maxDiff = 0;
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
    if (maxDiff >= aps.adaptive_prediction_threshold) {
      uint64_t attrValue =
        pointCloud.getReflectance(indexesLOD[predictorIndex]);

      // base case: weighted average of n neighbours
      predictor.predMode = 0;
      uint64_t attrPred = predictor.predictReflectance(pointCloud, indexesLOD);
      int64_t attrResidualQuant =
        computeReflectanceResidual(attrValue, attrPred, quant);

      double best_score = attrResidualQuant
        + kAttrPredLambdaR * (quant.stepSize() >> kFixedPointAttributeShift);

      for (int i = 0; i < predictor.neighborCount; i++) {
        if (i == aps.max_num_direct_predictors)
          break;

        attrPred = pointCloud.getReflectance(
          indexesLOD[predictor.neighbors[i].predictorIndex]);
        attrResidualQuant =
          computeReflectanceResidual(attrValue, attrPred, quant);

        double idxBits = i + (i == aps.max_num_direct_predictors - 1 ? 1 : 2);
        double score = attrResidualQuant
          + idxBits * kAttrPredLambdaR
            * (quant.stepSize() >> kFixedPointAttributeShift);

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
  const std::vector<Quantizers>& quantLayers,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const uint32_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    aps, pointCloud, 0, predictors, numberOfPointsPerLOD, indexesLOD);

  const int64_t clipMax = (1ll << desc.attr_bitdepth) - 1;
  PCCResidualsEntropyEstimator context;
  int zero_cnt = 0;
  std::vector<int> zerorun;
  zerorun.reserve(pointCount);
  std::vector<uint32_t> residual;
  residual.resize(pointCount);

  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == numberOfPointsPerLOD[quantLayer]) {
      quantLayer = std::min(int(quantLayers.size()) - 1, quantLayer + 1);
    }
    auto& quant = quantLayers[quantLayer];
    auto& predictor = predictors[predictorIndex];

    computeReflectancePredictionWeights(
      aps, pointCloud, indexesLOD, predictorIndex, predictor, encoder, context,
      quant[0]);

    const uint32_t pointIndex = indexesLOD[predictorIndex];
    const uint64_t reflectance = pointCloud.getReflectance(pointIndex);
    const uint16_t predictedReflectance =
      predictor.predictReflectance(pointCloud, indexesLOD);
    const int64_t quantAttValue = reflectance;
    const int64_t quantPredAttValue = predictedReflectance;
    const int64_t delta = quant[0].quantize(
      (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
    const uint32_t attValue0 = uint32_t(IntToUInt(long(delta)));
    const int64_t reconstructedDelta =
      divExp2RoundHalfUp(quant[0].scale(delta), kFixedPointAttributeShift);
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
    if (predictor.maxDiff >= aps.adaptive_prediction_threshold) {
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
  const Quantizers& quant)
{
  Vec3<int64_t> residuals;
  const int64_t quantAttValue = color[0];
  const int64_t quantPredAttValue = predictedColor[0];
  const int64_t delta = quant[0].quantize(
    (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
  residuals[0] = IntToUInt(delta);
  for (size_t k = 1; k < 3; ++k) {
    const int64_t quantAttValue = color[k];
    const int64_t quantPredAttValue = predictedColor[k];
    const int64_t delta = quant[1].quantize(
      (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
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
  const Quantizers& quant)
{
  predictor.computeWeights();
  predictor.maxDiff = 0;
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

    if (maxDiff >= aps.adaptive_prediction_threshold) {
      Vec3<uint8_t> attrValue =
        pointCloud.getColor(indexesLOD[predictorIndex]);

      // base case: weighted average of n neighbours
      predictor.predMode = 0;
      Vec3<uint8_t> attrPred = predictor.predictColor(pointCloud, indexesLOD);
      Vec3<int64_t> attrResidualQuant =
        computeColorResiduals(attrValue, attrPred, quant);

      double best_score = attrResidualQuant[0] + attrResidualQuant[1]
        + attrResidualQuant[2]
        + kAttrPredLambdaC
          * (double)(quant[0].stepSize() >> kFixedPointAttributeShift);

      for (int i = 0; i < predictor.neighborCount; i++) {
        if (i == aps.max_num_direct_predictors)
          break;

        attrPred = pointCloud.getColor(
          indexesLOD[predictor.neighbors[i].predictorIndex]);
        attrResidualQuant = computeColorResiduals(attrValue, attrPred, quant);

        double idxBits = i + (i == aps.max_num_direct_predictors - 1 ? 1 : 2);
        double score = attrResidualQuant[0] + attrResidualQuant[1]
          + attrResidualQuant[2]
          + idxBits * kAttrPredLambdaC
            * (quant[0].stepSize() >> kFixedPointAttributeShift);

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
  const std::vector<Quantizers>& quantLayers,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    aps, pointCloud, 0, predictors, numberOfPointsPerLOD, indexesLOD);

  const int64_t clipMax = (1ll << desc.attr_bitdepth) - 1;
  uint32_t values[3];
  PCCResidualsEntropyEstimator context;
  int zero_cnt = 0;
  std::vector<int> zerorun;
  std::vector<uint32_t> residual[3];
  for (int i = 0; i < 3; i++) {
    residual[i].resize(pointCount);
  }

  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == numberOfPointsPerLOD[quantLayer]) {
      quantLayer = std::min(int(quantLayers.size()) - 1, quantLayer + 1);
    }
    auto& quant = quantLayers[quantLayer];
    auto& predictor = predictors[predictorIndex];

    computeColorPredictionWeights(
      aps, pointCloud, indexesLOD, predictorIndex, predictor, encoder, context,
      quant);
    const auto pointIndex = indexesLOD[predictorIndex];
    const Vec3<uint8_t> color = pointCloud.getColor(pointIndex);
    const Vec3<uint8_t> predictedColor =
      predictor.predictColor(pointCloud, indexesLOD);

    Vec3<uint8_t> reconstructedColor;
    for (int k = 0; k < 3; ++k) {
      const auto& q = quant[std::min(k, 1)];
      int64_t residual = color[k] - predictedColor[k];

      int64_t residualQ = q.quantize(residual << kFixedPointAttributeShift);
      int64_t residualR =
        divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);

      values[k] = uint32_t(IntToUInt(long(residualQ)));

      int64_t recon = predictedColor[k] + residualR;
      reconstructedColor[k] = uint8_t(PCCClip(recon, int64_t(0), clipMax));
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
    if (predictor.maxDiff >= aps.adaptive_prediction_threshold) {
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
  const std::vector<Quantizers>& quantLayers,
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
  int* attributes = new int[attribCount * voxelCount];
  int* coefficients = new int[attribCount * voxelCount];

  // Populate input arrays.
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
    const auto reflectance = pointCloud.getReflectance(packedVoxel[n].index);
    attributes[attribCount * n] = reflectance;
  }

  // Transform.
  regionAdaptiveHierarchicalTransform(
    aps.raht_prediction_enabled_flag, quantLayers, mortonCode, attributes,
    attribCount, voxelCount, coefficients);

  // Entropy encode.
  int zero_cnt = 0;
  uint32_t value;
  for (int n = 0; n < voxelCount; ++n) {
    const int64_t detail = IntToUInt(coefficients[n]);
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
AttributeEncoder::encodeColorsTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const std::vector<Quantizers>& quantLayers,
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
  int* attributes = new int[attribCount * voxelCount];
  int* coefficients = new int[attribCount * voxelCount];

  // Populate input arrays.
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
    const auto color = pointCloud.getColor(packedVoxel[n].index);
    attributes[attribCount * n] = color[0];
    attributes[attribCount * n + 1] = color[1];
    attributes[attribCount * n + 2] = color[2];
  }

  // Transform.
  regionAdaptiveHierarchicalTransform(
    aps.raht_prediction_enabled_flag, quantLayers, mortonCode, attributes,
    attribCount, voxelCount, coefficients);

  // Entropy encode.
  uint32_t values[attribCount];
  int zero_cnt = 0;
  for (int n = 0; n < voxelCount; ++n) {
    for (int d = 0; d < attribCount; ++d) {
      const int64_t detail = IntToUInt(coefficients[voxelCount * d + n]);
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
AttributeEncoder::encodeColorsLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const std::vector<Quantizers>& quantLayers,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    aps, pointCloud, 0, predictors, numberOfPointsPerLOD, indexesLOD);

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    predictors[predictorIndex].computeWeights();
  }
  std::vector<uint64_t> weights;
  if (!aps.scalable_lifting_enabled_flag) {
    PCCComputeQuantizationWeights(predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      predictors, numberOfPointsPerLOD, pointCount, 0, weights);
  }
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
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == numberOfPointsPerLOD[quantLayer]) {
      quantLayer = std::min(int(quantLayers.size()) - 1, quantLayer + 1);
    }
    auto& quant = quantLayers[quantLayer];

    const int64_t quantWeight = weights[predictorIndex];
    auto& color = colors[predictorIndex];
    const int64_t delta = quant[0].quantize(color[0] * quantWeight);
    const int64_t detail = IntToUInt(delta);
    assert(detail < std::numeric_limits<uint32_t>::max());
    const int64_t reconstructedDelta = quant[0].scale(delta);
    color[0] = reconstructedDelta / quantWeight;
    uint32_t values[3];
    values[0] = uint32_t(detail);
    for (size_t d = 1; d < 3; ++d) {
      const int64_t delta = quant[1].quantize(color[d] * quantWeight);
      const int64_t detail = IntToUInt(delta);
      assert(detail < std::numeric_limits<uint32_t>::max());
      const int64_t reconstructedDelta = quant[1].scale(delta);
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
  const std::vector<Quantizers>& quantLayers,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<PCCPredictor> predictors;
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexesLOD;

  buildPredictorsFast(
    aps, pointCloud, 0, predictors, numberOfPointsPerLOD, indexesLOD);

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    predictors[predictorIndex].computeWeights();
  }
  std::vector<uint64_t> weights;
  if (!aps.scalable_lifting_enabled_flag) {
    PCCComputeQuantizationWeights(predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      predictors, numberOfPointsPerLOD, pointCount, 0, weights);
  }

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
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == numberOfPointsPerLOD[quantLayer]) {
      quantLayer = std::min(int(quantLayers.size()) - 1, quantLayer + 1);
    }
    auto& quant = quantLayers[quantLayer];

    const int64_t quantWeight = weights[predictorIndex];
    auto& reflectance = reflectances[predictorIndex];
    const int64_t delta = quant[0].quantize(reflectance * quantWeight);
    const int64_t detail = IntToUInt(delta);
    assert(detail < std::numeric_limits<uint32_t>::max());
    const int64_t reconstructedDelta = quant[0].scale(delta);
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
