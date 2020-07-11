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
#include "io_hls.h"
#include "quantization.h"
#include "RAHT.h"
#include "FixedPoint.h"

#include <algorithm>

// todo(df): promote to per-attribute encoder parameter
static const double kAttrPredLambdaR = 0.01;
static const double kAttrPredLambdaC = 0.01;

namespace pcc {
//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

class PCCResidualsEncoder : protected AttributeContexts {
public:
  PCCResidualsEncoder(
    const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem);

  EntropyEncoder arithmeticEncoder;

  const AttributeContexts& getCtx() const { return *this; }

  void start(const SequenceParameterSet& sps, int numPoints);
  int stop();

  void encodeLastCompPredCoeffs(const std::vector<int8_t>& coeffs);
  void encodePredMode(int value, int max);
  void encodeRunLength(int runLength);
  void encodeSymbol(uint32_t value, int k1, int k2, int k3);
  void encode(int32_t value0, int32_t value1, int32_t value2);
  void encode(int32_t value);
};

//----------------------------------------------------------------------------

PCCResidualsEncoder::PCCResidualsEncoder(
  const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem)
  : AttributeContexts(ctxtMem)
{}

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
PCCResidualsEncoder::encodeLastCompPredCoeffs(
  const std::vector<int8_t>& coeffs)
{
  for (auto coeff : coeffs) {
    bool last_comp_pred_coeff_ne0 = coeff != 0;
    arithmeticEncoder.encode(last_comp_pred_coeff_ne0);
    if (last_comp_pred_coeff_ne0) {
      bool last_comp_pred_coeff_sign = coeff == -1;
      arithmeticEncoder.encode(last_comp_pred_coeff_sign);
    }
  }
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
PCCResidualsEncoder::encodeRunLength(int runLength)
{
  auto* ctx = ctxRunLen;
  for (int i = 0; i < std::min(3, runLength); i++, ctx++)
    arithmeticEncoder.encode(1, *ctx);

  if (runLength < 3) {
    arithmeticEncoder.encode(0, *ctx);
    return;
  }
  runLength -= 3;

  auto prefix = runLength >> 1;
  for (int i = 0; i < std::min(4, prefix); i++)
    arithmeticEncoder.encode(1, *ctx);

  if (runLength < 8) {
    arithmeticEncoder.encode(0, *ctx);
    arithmeticEncoder.encode(runLength & 1);
    return;
  }
  runLength -= 8;

  arithmeticEncoder.encodeExpGolomb(runLength, 2, *++ctx);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeSymbol(uint32_t value, int k1, int k2, int k3)
{
  bool isZero = value == 0;
  arithmeticEncoder.encode(isZero, ctxCoeffEqN[0][k1]);
  if (isZero)
    return;

  bool is1 = value == 1;
  arithmeticEncoder.encode(is1, ctxCoeffEqN[1][k2]);
  if (is1)
    return;

  arithmeticEncoder.encodeExpGolomb(
    value - 2, 1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(int32_t value0, int32_t value1, int32_t value2)
{
  int mag0 = abs(value0);
  int mag1 = abs(value1);
  int mag2 = abs(value2);

  int b0 = (mag1 == 0);
  int b1 = (mag1 <= 1);
  int b2 = (mag2 == 0);
  int b3 = (mag2 <= 1);
  encodeSymbol(mag1, 0, 0, 1);
  encodeSymbol(mag2, 1 + b0, 1 + b1, 1);

  auto mag0minusX = b0 && b2 ? mag0 - 1 : mag0;
  assert(mag0minusX >= 0);
  encodeSymbol(mag0minusX, 3 + (b0 << 1) + b2, 3 + (b1 << 1) + b3, 0);

  if (mag0)
    arithmeticEncoder.encode(value0 < 0);
  if (mag1)
    arithmeticEncoder.encode(value1 < 0);
  if (mag2)
    arithmeticEncoder.encode(value2 < 0);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(int32_t value)
{
  int mag = abs(value) - 1;
  encodeSymbol(mag, 0, 0, 0);
  arithmeticEncoder.encode(value < 0);
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
// AttributeEncoderIntf

AttributeEncoderIntf::~AttributeEncoderIntf() = default;

//============================================================================
// AttributeEncoder factory

std::unique_ptr<AttributeEncoderIntf>
makeAttributeEncoder()
{
  return std::unique_ptr<AttributeEncoder>(new AttributeEncoder());
}

//============================================================================
// AttributeEncoder Members

void
AttributeEncoder::encode(
  const SequenceParameterSet& sps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  AttributeBrickHeader& abh,
  AttributeContexts& ctxtMem,
  PCCPointSet3& pointCloud,
  PayloadBuffer* payload)
{
  QpSet qpSet = deriveQpSet(desc, attr_aps, abh);

  // generate LoDs if necessary
  if (attr_aps.lodParametersPresent() && _lods.empty())
    _lods.generate(
      attr_aps, abh, pointCloud.getPointCount() - 1, 0, pointCloud);

  // write abh
  write(sps, attr_aps, abh, payload);

  PCCResidualsEncoder encoder(abh, ctxtMem);
  encoder.start(sps, int(pointCloud.getPointCount()));

  if (desc.attr_num_dimensions_minus1 == 0) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      encodeReflectancesTransformRaht(
        desc, attr_aps, qpSet, pointCloud, encoder);
      break;

    case AttributeEncoding::kPredictingTransform:
      encodeReflectancesPred(desc, attr_aps, qpSet, pointCloud, encoder);
      break;

    case AttributeEncoding::kLiftingTransform:
      encodeReflectancesLift(desc, attr_aps, qpSet, pointCloud, encoder);
      break;
    }
  } else if (desc.attr_num_dimensions_minus1 == 2) {
    switch (attr_aps.attr_encoding) {
    case AttributeEncoding::kRAHTransform:
      encodeColorsTransformRaht(desc, attr_aps, qpSet, pointCloud, encoder);
      break;

    case AttributeEncoding::kPredictingTransform:
      encodeColorsPred(desc, attr_aps, qpSet, pointCloud, encoder);
      break;

    case AttributeEncoding::kLiftingTransform:
      encodeColorsLift(desc, attr_aps, qpSet, pointCloud, encoder);
      break;
    }
  } else {
    assert(
      desc.attr_num_dimensions_minus1 == 0
      || desc.attr_num_dimensions_minus1 == 2);
  }

  uint32_t acDataLen = encoder.stop();
  std::copy_n(
    encoder.arithmeticEncoder.buffer(), acDataLen,
    std::back_inserter(*payload));

  // save the context state for re-use by a future slice if required
  ctxtMem = encoder.getCtx();
}

//----------------------------------------------------------------------------

bool
AttributeEncoder::isReusable(
  const AttributeParameterSet& aps, const AttributeBrickHeader& abh) const
{
  return _lods.isReusable(aps, abh);
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
  predictor.predMode = 0;
  predictor.maxDiff = 0;
  if (predictor.neighborCount > 1 && aps.max_num_direct_predictors) {
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

      // base case: start with the first neighbour
      // NB: skip evaluation of mode 0 (weighted average of n neighbours)
      predictor.predMode = 1;
      uint64_t attrPred = predictor.predictReflectance(pointCloud, indexesLOD);
      int64_t attrResidualQuant =
        computeReflectanceResidual(attrValue, attrPred, quant);

      // NB: idxBits is not included in the score
      int64_t best_score = attrResidualQuant;

      for (int i = 1; i < predictor.neighborCount; i++) {
        if (i == aps.max_num_direct_predictors)
          break;

        attrPred = pointCloud.getReflectance(
          indexesLOD[predictor.neighbors[i].predictorIndex]);
        attrResidualQuant =
          computeReflectanceResidual(attrValue, attrPred, quant);

        if (attrResidualQuant < best_score) {
          best_score = attrResidualQuant;
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
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const uint32_t pointCount = pointCloud.getPointCount();
  const int64_t clipMax = (1ll << desc.bitdepth) - 1;
  PCCResidualsEntropyEstimator context;
  int zero_cnt = 0;
  std::vector<int> zerorun;
  zerorun.reserve(pointCount);
  std::vector<uint32_t> residual;
  residual.resize(pointCount);

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
      aps, pointCloud, _lods.indexes, predictorIndex, predictor, encoder,
      context, quant[0]);

    const uint64_t reflectance = pointCloud.getReflectance(pointIndex);
    const attr_t predictedReflectance =
      predictor.predictReflectance(pointCloud, _lods.indexes);
    const int64_t quantAttValue = reflectance;
    const int64_t quantPredAttValue = predictedReflectance;
    const int64_t delta = quant[0].quantize(
      (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
    const auto attValue0 = delta;
    const int64_t reconstructedDelta =
      divExp2RoundHalfUp(quant[0].scale(delta), kFixedPointAttributeShift);
    const int64_t reconstructedQuantAttValue =
      quantPredAttValue + reconstructedDelta;
    const attr_t reconstructedReflectance =
      attr_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));

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
  encoder.encodeRunLength(zerorun[run_index]);
  zero_cnt = zerorun[run_index++];

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = _lods.predictors[predictorIndex];
    if (predictor.maxDiff >= aps.adaptive_prediction_threshold) {
      encoder.encodePredMode(
        predictor.predMode, aps.max_num_direct_predictors);
    }
    if (zero_cnt > 0)
      zero_cnt--;
    else {
      encoder.encode(residual[predictorIndex]);
      encoder.encodeRunLength(zerorun[run_index]);
      zero_cnt = zerorun[run_index++];
    }
  }
}

//----------------------------------------------------------------------------

Vec3<int64_t>
AttributeEncoder::computeColorResiduals(
  const AttributeParameterSet& aps,
  const Vec3<attr_t> color,
  const Vec3<attr_t> predictedColor,
  const Quantizers& quant)
{
  Vec3<int64_t> residuals;
  const int64_t quantAttValue = color[0];
  const int64_t quantPredAttValue = predictedColor[0];
  const int64_t delta = quant[0].quantize(
    (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
  residuals[0] = IntToUInt(delta);
  const int64_t residual0 =
    divExp2RoundHalfUp(quant[0].scale(delta), kFixedPointAttributeShift);
  for (size_t k = 1; k < 3; ++k) {
    const int64_t quantAttValue = color[k];
    const int64_t quantPredAttValue = predictedColor[k];
    if (aps.inter_component_prediction_enabled_flag) {
      const int64_t delta = quant[1].quantize(
        (quantAttValue - quantPredAttValue - residual0)
        << kFixedPointAttributeShift);
      residuals[k] = IntToUInt(delta);
    } else {
      const int64_t delta = quant[1].quantize(
        (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
      residuals[k] = IntToUInt(delta);
    }
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
  predictor.maxDiff = 0;
  if (predictor.neighborCount > 1 && aps.max_num_direct_predictors) {
    int64_t minValue[3] = {0, 0, 0};
    int64_t maxValue[3] = {0, 0, 0};
    for (int i = 0; i < predictor.neighborCount; ++i) {
      const Vec3<attr_t> colorNeighbor =
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
      Vec3<attr_t> attrValue = pointCloud.getColor(indexesLOD[predictorIndex]);

      // base case: weighted average of n neighbours
      predictor.predMode = 0;
      Vec3<attr_t> attrPred = predictor.predictColor(pointCloud, indexesLOD);
      Vec3<int64_t> attrResidualQuant =
        computeColorResiduals(aps, attrValue, attrPred, quant);

      double best_score = attrResidualQuant[0] + attrResidualQuant[1]
        + attrResidualQuant[2]
        + kAttrPredLambdaC
          * (double)(quant[0].stepSize() >> kFixedPointAttributeShift);

      for (int i = 0; i < predictor.neighborCount; i++) {
        if (i == aps.max_num_direct_predictors)
          break;

        attrPred = pointCloud.getColor(
          indexesLOD[predictor.neighbors[i].predictorIndex]);
        attrResidualQuant =
          computeColorResiduals(aps, attrValue, attrPred, quant);

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
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const size_t pointCount = pointCloud.getPointCount();

  Vec3<int64_t> clipMax{(1 << desc.bitdepth) - 1,
                        (1 << desc.bitdepthSecondary) - 1,
                        (1 << desc.bitdepthSecondary) - 1};

  int32_t values[3];
  PCCResidualsEntropyEstimator context;
  int zero_cnt = 0;
  std::vector<int> zerorun;
  std::vector<int32_t> residual[3];
  for (int i = 0; i < 3; i++) {
    residual[i].resize(pointCount);
  }
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }
    const auto pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);
    auto& predictor = _lods.predictors[predictorIndex];

    computeColorPredictionWeights(
      aps, pointCloud, _lods.indexes, predictorIndex, predictor, encoder,
      context, quant);
    const Vec3<attr_t> color = pointCloud.getColor(pointIndex);
    const Vec3<attr_t> predictedColor =
      predictor.predictColor(pointCloud, _lods.indexes);

    Vec3<attr_t> reconstructedColor;
    int64_t residual0 = 0;
    for (int k = 0; k < 3; ++k) {
      const auto& q = quant[std::min(k, 1)];
      int64_t residual = color[k] - predictedColor[k];

      int64_t residualQ = q.quantize(residual << kFixedPointAttributeShift);
      int64_t residualR =
        divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);

      if (aps.inter_component_prediction_enabled_flag && k > 0) {
        residual = residual - residual0;
        residualQ = q.quantize(residual << kFixedPointAttributeShift);
        residualR = residual0
          + divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);
      }

      if (k == 0)
        residual0 = residualR;

      values[k] = residualQ;

      int64_t recon = predictedColor[k] + residualR;
      reconstructedColor[k] = attr_t(PCCClip(recon, int64_t(0), clipMax[k]));
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
  encoder.encodeRunLength(zerorun[run_index]);
  zero_cnt = zerorun[run_index++];
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    auto& predictor = _lods.predictors[predictorIndex];
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
      encoder.encodeRunLength(zerorun[run_index]);
      zero_cnt = zerorun[run_index++];
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeReflectancesTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n]);
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Allocate arrays.
  const int attribCount = 1;
  std::vector<int64_t> mortonCode(voxelCount);
  std::vector<int> attributes(attribCount * voxelCount);
  std::vector<int> coefficients(attribCount * voxelCount);
  std::vector<Qps> pointQpOffsets(voxelCount);

  // Populate input arrays.
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
    const auto reflectance = pointCloud.getReflectance(packedVoxel[n].index);
    attributes[attribCount * n] = reflectance;
    pointQpOffsets[n] = qpSet.regionQpOffset(pointCloud[packedVoxel[n].index]);
  }

  const int rahtPredThreshold[2] = {aps.raht_prediction_threshold0,
                                    aps.raht_prediction_threshold1};

  // Transform.
  regionAdaptiveHierarchicalTransform(
    aps.raht_prediction_enabled_flag, rahtPredThreshold, qpSet,
    pointQpOffsets.data(), mortonCode.data(), attributes.data(), attribCount,
    voxelCount, coefficients.data());

  // Entropy encode.
  int zero_cnt = 0;
  for (int n = 0; n < voxelCount; ++n) {
    auto value = coefficients[n];
    if (!value)
      ++zero_cnt;
    else {
      encoder.encodeRunLength(zero_cnt);
      encoder.encode(value);
      zero_cnt = 0;
    }
  }
  encoder.encodeRunLength(zero_cnt);

  const int64_t maxReflectance = (1 << desc.bitdepth) - 1;
  const int64_t minReflectance = 0;
  for (int n = 0; n < voxelCount; n++) {
    int64_t val = attributes[attribCount * n];
    const attr_t reflectance =
      attr_t(PCCClip(val, minReflectance, maxReflectance));
    pointCloud.setReflectance(packedVoxel[n].index, reflectance);
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeColorsTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    packedVoxel[n].mortonCode = mortonAddr(pointCloud[n]);
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Allocate arrays.
  const int attribCount = 3;
  std::vector<int64_t> mortonCode(voxelCount);
  std::vector<int> attributes(attribCount * voxelCount);
  std::vector<int> coefficients(attribCount * voxelCount);
  std::vector<Qps> pointQpOffsets(voxelCount);

  // Populate input arrays.
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
    const auto color = pointCloud.getColor(packedVoxel[n].index);
    attributes[attribCount * n] = color[0];
    attributes[attribCount * n + 1] = color[1];
    attributes[attribCount * n + 2] = color[2];
    pointQpOffsets[n] = qpSet.regionQpOffset(pointCloud[packedVoxel[n].index]);
  }

  const int rahtPredThreshold[2] = {aps.raht_prediction_threshold0,
                                    aps.raht_prediction_threshold1};

  // Transform.
  regionAdaptiveHierarchicalTransform(
    aps.raht_prediction_enabled_flag, rahtPredThreshold, qpSet,
    pointQpOffsets.data(), mortonCode.data(), attributes.data(), attribCount,
    voxelCount, coefficients.data());

  // Entropy encode.
  int values[attribCount];
  int zero_cnt = 0;
  for (int n = 0; n < voxelCount; ++n) {
    for (int d = 0; d < attribCount; ++d) {
      values[d] = coefficients[voxelCount * d + n];
    }
    if (!values[0] && !values[1] && !values[2])
      ++zero_cnt;
    else {
      encoder.encodeRunLength(zero_cnt);
      encoder.encode(values[0], values[1], values[2]);
      zero_cnt = 0;
    }
  }
  encoder.encodeRunLength(zero_cnt);

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
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeColorsLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint64_t> weights;

  if (!aps.scalable_lifting_enabled_flag) {
    PCCComputeQuantizationWeights(_lods.predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      _lods.predictors, _lods.numPointsInLod, pointCount, 0, weights);
  }

  const size_t lodCount = _lods.numPointsInLod.size();
  std::vector<Vec3<int64_t>> colors;
  colors.resize(pointCount);

  for (size_t index = 0; index < pointCount; ++index) {
    const auto& color = pointCloud.getColor(_lods.indexes[index]);
    for (size_t d = 0; d < 3; ++d) {
      colors[index][d] = int32_t(color[d]) << kFixedPointAttributeShift;
    }
  }

  for (size_t i = 0; (i + 1) < lodCount; ++i) {
    const size_t lodIndex = lodCount - i - 1;
    const size_t startIndex = _lods.numPointsInLod[lodIndex - 1];
    const size_t endIndex = _lods.numPointsInLod[lodIndex];
    PCCLiftPredict(_lods.predictors, startIndex, endIndex, true, colors);
    PCCLiftUpdate(
      _lods.predictors, weights, startIndex, endIndex, true, colors);
  }

  // Per level-of-detail coefficients {-1,0,1} for last component prediction
  int8_t lastCompPredCoeff = 0;
  std::vector<int8_t> lastCompPredCoeffs;
  if (aps.last_component_prediction_enabled_flag) {
    lastCompPredCoeffs = computeLastComponentPredictionCoeff(colors);
    encoder.encodeLastCompPredCoeffs(lastCompPredCoeffs);
    lastCompPredCoeff = lastCompPredCoeffs[0];
  }

  int zero_cnt = 0;
  int quantLayer = 0;
  int lod = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }

    if (predictorIndex == _lods.numPointsInLod[lod]) {
      lod++;
      if (aps.last_component_prediction_enabled_flag)
        lastCompPredCoeff = lastCompPredCoeffs[lod];
    }

    const auto pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);

    const int64_t iQuantWeight = irsqrt(weights[predictorIndex]);
    const int64_t quantWeight =
      (weights[predictorIndex] * iQuantWeight + (1ull << 39)) >> 40;

    auto& color = colors[predictorIndex];
    int values[3];
    values[0] = quant[0].quantize(color[0] * quantWeight);
    int64_t scaled = quant[0].scale(values[0]);
    color[0] = divExp2RoundHalfInf(scaled * iQuantWeight, 40);

    values[1] = quant[1].quantize(color[1] * quantWeight);
    scaled = quant[1].scale(values[1]);
    color[1] = divExp2RoundHalfInf(scaled * iQuantWeight, 40);

    color[2] -= lastCompPredCoeff * color[1];
    scaled *= lastCompPredCoeff;

    values[2] = quant[1].quantize(color[2] * quantWeight);
    scaled += quant[1].scale(values[2]);
    color[2] = divExp2RoundHalfInf(scaled * iQuantWeight, 40);

    if (!values[0] && !values[1] && !values[2])
      ++zero_cnt;
    else {
      encoder.encodeRunLength(zero_cnt);
      encoder.encode(values[0], values[1], values[2]);
      zero_cnt = 0;
    }
  }
  encoder.encodeRunLength(zero_cnt);

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
      color[d] = attr_t(PCCClip(color0[d], 0, clipMax[d]));
    }
    pointCloud.setColor(_lods.indexes[f], color);
  }
}

//----------------------------------------------------------------------------

std::vector<int8_t>
AttributeEncoder::computeLastComponentPredictionCoeff(
  const std::vector<Vec3<int64_t>>& coeffs)
{
  std::vector<int8_t> signs(_lods.numPointsInLod.size(), 0);

  int numGt0 = 0;
  int numLt0 = 0;
  int lod = 0;
  int lodPointIdx = 0;
  for (size_t coeffIdx = 0; coeffIdx < coeffs.size(); ++coeffIdx) {
    auto& color = coeffs[coeffIdx];
    lodPointIdx++;

    int mult = color[1] * color[2];
    if (mult > 0)
      numGt0++;
    else if (mult < 0)
      numLt0++;

    if (coeffIdx == _lods.numPointsInLod[lod] - 1) {
      constexpr double threshold = 0.7;
      if (numGt0 > threshold * lodPointIdx)
        signs[lod] = 1;
      else if (numLt0 > threshold * lodPointIdx)
        signs[lod] = -1;
      else
        signs[lod] = 0;

      lodPointIdx = 0;
      numGt0 = 0;
      numLt0 = 0;
      lod++;
    }
  }

  // check if residual prediction is effective
  lod = 0;
  int64_t sumPredCoeff = 0;
  int64_t sumOrigCoeff = 0;
  for (size_t coeffIdx = 0; coeffIdx < coeffs.size(); ++coeffIdx) {
    auto& coeff = coeffs[coeffIdx];
    if (signs[lod] == 1) {
      sumPredCoeff += abs(coeff[2] - coeff[1]);
      sumOrigCoeff += abs(coeff[2]);
    } else if (signs[lod] == -1) {
      sumPredCoeff += abs(coeff[2] + coeff[1]);
      sumOrigCoeff += abs(coeff[2]);
    }

    if (coeffIdx == _lods.numPointsInLod[lod] - 1) {
      constexpr double threshold2 = 0.9;
      if (signs[lod] != 0 && sumPredCoeff > threshold2 * sumOrigCoeff)
        signs[lod] = 0;

      sumPredCoeff = 0;
      sumOrigCoeff = 0;
      lod++;
    }
  }

  return signs;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeReflectancesLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint64_t> weights;

  if (!aps.scalable_lifting_enabled_flag) {
    PCCComputeQuantizationWeights(_lods.predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      _lods.predictors, _lods.numPointsInLod, pointCount, 0, weights);
  }

  const size_t lodCount = _lods.numPointsInLod.size();
  std::vector<int64_t> reflectances;
  reflectances.resize(pointCount);

  for (size_t index = 0; index < pointCount; ++index) {
    reflectances[index] =
      int32_t(pointCloud.getReflectance(_lods.indexes[index]))
      << kFixedPointAttributeShift;
  }

  for (size_t i = 0; (i + 1) < lodCount; ++i) {
    const size_t lodIndex = lodCount - i - 1;
    const size_t startIndex = _lods.numPointsInLod[lodIndex - 1];
    const size_t endIndex = _lods.numPointsInLod[lodIndex];
    PCCLiftPredict(_lods.predictors, startIndex, endIndex, true, reflectances);
    PCCLiftUpdate(
      _lods.predictors, weights, startIndex, endIndex, true, reflectances);
  }

  // compress
  int zero_cnt = 0;
  int quantLayer = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }
    const auto pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);

    const int64_t iQuantWeight = irsqrt(weights[predictorIndex]);
    const int64_t quantWeight =
      (weights[predictorIndex] * iQuantWeight + (1ull << 39)) >> 40;

    auto& reflectance = reflectances[predictorIndex];
    const int64_t delta = quant[0].quantize(reflectance * quantWeight);
    const auto detail = delta;
    const int64_t reconstructedDelta = quant[0].scale(delta);
    reflectance = divExp2RoundHalfInf(reconstructedDelta * iQuantWeight, 40);
    if (!detail)
      ++zero_cnt;
    else {
      encoder.encodeRunLength(zero_cnt);
      encoder.encode(detail);
      zero_cnt = 0;
    }
  }
  encoder.encodeRunLength(zero_cnt);

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
    const int64_t refl =
      divExp2RoundHalfInf(reflectances[f], kFixedPointAttributeShift);
    pointCloud.setReflectance(
      _lods.indexes[f], attr_t(PCCClip(refl, int64_t(0), maxReflectance)));
  }
}

//============================================================================
// estimation of dist2

int
estimateDist2(
  const PCCPointSet3& cloud,
  int32_t samplingPeriod,
  int32_t searchRange,
  float percentileEstimate)
{
  int32_t pointCount = cloud.getPointCount();
  if (pointCount < 2)
    return 0;

  std::vector<int64_t> dists;
  dists.reserve(pointCount / samplingPeriod + 1);

  for (int32_t index = 0; index < pointCount; index += samplingPeriod) {
    auto k0 = std::max(0, index - searchRange);
    auto k1 = std::min(pointCount - 1, index + searchRange);
    auto d2 = std::numeric_limits<int64_t>::max();
    for (auto k = k0; k <= k1; ++k) {
      if (k == index)
        continue;

      d2 = std::min(d2, (cloud[index] - cloud[k]).getNorm2<int64_t>());
    }
    dists.push_back(d2);
  }

  int p = int(std::floor(dists.size() * percentileEstimate));

  std::nth_element(dists.begin(), dists.begin() + p, dists.end());
  int64_t dist2 = dists[p];
  int shiftBits = 0;
  while ((int64_t(3) << (shiftBits << 1)) < dist2 && shiftBits < 20)
    ++shiftBits;

  return shiftBits;
}

//============================================================================

} /* namespace pcc */
