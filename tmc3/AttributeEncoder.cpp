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

#include "DualLutCoder.h"
#include "attribute_raw.h"
#include "constants.h"
#include "entropy.h"
#include "io_hls.h"
#include "quantization.h"
#include "RAHT.h"
#include "FixedPoint.h"

#include <algorithm>

// todo(df): promote to per-attribute encoder parameter
static const double kAttrPredLambdaR = 0.01;
static const double kAttrPredLambdaC = 0.14;

namespace pcc {
//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

class PCCResidualsEncoder : protected AttributeContexts {
public:
  PCCResidualsEncoder(
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    const AttributeContexts& ctxtMem);

  EntropyEncoder arithmeticEncoder;

  const AttributeContexts& getCtx() const { return *this; }

  void start(const SequenceParameterSet& sps, int numPoints);
  int stop();

  void encodeRunLength(int runLength);
  void encodeSymbol(uint32_t value, int k1, int k2, int k3);
  void encode(int32_t value0, int32_t value1, int32_t value2);
  void encode(int32_t value);

  int availPredModes;
  double bitsPtColor(Vec3<int32_t> value, int parity);
  double bitsPtRefl(int32_t value, int parity);

  // Encoder side residual cost calculation
  const int scaleRes = 1 << 20;
  const int windowLog2 = 6;
  int probResGt0[3];  //prob of residuals larger than 0: 1 for each component
  int probResGt1[3];  //prob of residuals larger than 1: 1 for each component
  void resStatUpdateColor(Vec3<int32_t> values);
  void resStatUpdateRefl(int32_t values);
  void resStatReset();
};

//----------------------------------------------------------------------------

PCCResidualsEncoder::PCCResidualsEncoder(
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  const AttributeContexts& ctxtMem)
  : AttributeContexts(ctxtMem)
{
  availPredModes =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;

  resStatReset();
}

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
PCCResidualsEncoder::resStatReset()
{
  for (int k = 0; k < 3; k++)
    probResGt0[k] = probResGt1[k] = (scaleRes >> 1);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatUpdateColor(Vec3<int32_t> value)
{
  for (int k = 0; k < 3; k++) {
    probResGt0[k] += value[k] ? (scaleRes - probResGt0[k]) >> windowLog2
                              : -((probResGt0[k]) >> windowLog2);
    if (value[k])
      probResGt1[k] += abs(value[k]) > 1
        ? (scaleRes - probResGt1[k]) >> windowLog2
        : -((probResGt1[k]) >> windowLog2);
  }
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatUpdateRefl(int32_t value)
{
  probResGt0[0] += value ? (scaleRes - probResGt0[0]) >> windowLog2
                         : -(probResGt0[0] >> windowLog2);
  if (value)
    probResGt1[0] += abs(value) > 1 ? (scaleRes - probResGt1[0]) >> windowLog2
                                    : -(probResGt1[0] >> windowLog2);
}

//----------------------------------------------------------------------------

double
PCCResidualsEncoder::bitsPtColor(Vec3<int32_t> value, int mode)
{
  if (availPredModes == 4) {
    value[1] = 2 * abs(value[1]) + (mode >> 1);
    value[2] = 2 * abs(value[2]) + (mode & 1);
  } else if (availPredModes == 3) {
    value[1] = 2 * abs(value[1]) + (mode > 0);
    if (mode > 0)
      value[2] = 2 * abs(value[2]) + (mode - 1);
  } else if (availPredModes == 2) {
    value[1] = 2 * abs(value[1]) + (mode & 1);
  }

  int log2scaleRes = ilog2(uint32_t(scaleRes));
  double bits = 0;
  for (int k = 0; k < 3; k++) {
    bits += value[k] ? log2scaleRes - log2(probResGt0[k])
                     : log2scaleRes - log2(scaleRes - probResGt0[k]);  //Gt0
    int mag = abs(value[k]);
    if (mag) {
      bits += mag > 1 ? log2scaleRes - log2(probResGt1[k])
                      : log2scaleRes - log2(scaleRes - probResGt1[k]);  //Gt1
      bits += 1;  //sign bit.
      if (mag > 1)
        bits += 2.0 * log2(mag - 1.0) + 1.0;  //EG0 approximation.
    }
  }

  return bits;
}

//----------------------------------------------------------------------------

double
PCCResidualsEncoder::bitsPtRefl(int32_t value, int mode)
{
  if (availPredModes == 4) {
    value = (abs(value) << 2) + (mode);
  } else if (availPredModes == 3) {
    if (mode > 0)
      value = (abs(value) << 1) + (mode - 1);
    value = (abs(value) << 1) + (mode > 0);
  } else if (availPredModes == 2) {
    value = (abs(value) << 1) + (mode & 1);
  }

  int log2scaleRes = ilog2((uint32_t)scaleRes);
  double bits = 0;
  bits += value ? log2scaleRes - log2(probResGt0[0])
                : log2scaleRes - log2(scaleRes - probResGt0[0]);  //Gt0
  int mag = abs(value);
  if (mag) {
    bits += mag > 1 ? log2scaleRes - log2(probResGt1[0])
                    : log2scaleRes - log2(scaleRes - probResGt1[0]);  //Gt1
    bits += 1;  //sign bit.
    if (mag > 1)
      bits += 2.0 * log2(mag - 1.0) + 1.0;  //EG0 approximation.
  }
  return bits;
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
  arithmeticEncoder.encode(value > 0, ctxCoeffGtN[0][k1]);
  if (!value)
    return;

  arithmeticEncoder.encode(--value > 0, ctxCoeffGtN[1][k2]);
  if (!value)
    return;

  arithmeticEncoder.encodeExpGolomb(
    --value, 1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);
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
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    AttrRawEncoder::encode(sps, desc, attr_aps, abh, pointCloud, payload);
    return;
  }

  // Encoders are able to modify the slice header:
  _abh = &abh;

  QpSet qpSet = deriveQpSet(desc, attr_aps, abh);

  // generate LoDs if necessary
  if (attr_aps.lodParametersPresent() && _lods.empty())
    _lods.generate(
      attr_aps, abh, pointCloud.getPointCount() - 1, 0, pointCloud);

  PCCResidualsEncoder encoder(attr_aps, abh, ctxtMem);
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

    case AttributeEncoding::kRaw:
      // Already handled
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

    case AttributeEncoding::kRaw:
      // Already handled
      break;
    }
  } else {
    assert(
      desc.attr_num_dimensions_minus1 == 0
      || desc.attr_num_dimensions_minus1 == 2);
  }

  uint32_t acDataLen = encoder.stop();

  // write abh
  write(sps, attr_aps, abh, payload);
  _abh = nullptr;

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
  return delta;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::decidePredModeRefl(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexesLOD,
  const uint32_t predictorIndex,
  PCCPredictor& predictor,
  PCCResidualsEncoder& encoder,
  PCCResidualsEntropyEstimator& context,
  const Quantizer& quant)
{
  uint64_t attrValue = pointCloud.getReflectance(indexesLOD[predictorIndex]);

  // base case: start with the first neighbour
  // NB: skip evaluation of mode 0 (weighted average of n neighbours)
  int startpredIndex = aps.direct_avg_predictor_disabled_flag;
  predictor.predMode = startpredIndex;
  uint64_t attrPred = predictor.predictReflectance(pointCloud, indexesLOD);
  int64_t attrResidualQuant =
    computeReflectanceResidual(attrValue, attrPred, quant);

  // NB: idxBits is not included in the score
  int mode = predictor.predMode - aps.direct_avg_predictor_disabled_flag;
  int64_t best_score = encoder.bitsPtRefl(attrResidualQuant, mode);

  for (int i = startpredIndex; i < predictor.neighborCount; i++) {
    if (i == aps.max_num_direct_predictors)
      break;

    attrPred = pointCloud.getReflectance(
      indexesLOD[predictor.neighbors[i].predictorIndex]);
    attrResidualQuant = computeReflectanceResidual(attrValue, attrPred, quant);

    mode = i + !aps.direct_avg_predictor_disabled_flag;
    int64_t score = encoder.bitsPtRefl(attrResidualQuant, mode);
    if (score < best_score) {
      best_score = score;
      predictor.predMode = i + 1;
      // NB: setting predictor.neighborCount = 1 will cause issues
      // with reconstruction.
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodePredModeRefl(
  const AttributeParameterSet& aps, int predMode, int32_t& coeff)
{
  int coeffSign = coeff < 0 ? -1 : 1;
  int coeffAbs = abs(coeff);

  int mode = predMode - aps.direct_avg_predictor_disabled_flag;
  int maxcand =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;
  switch (maxcand) {
  case 4: coeff = coeffSign * ((coeffAbs << 2) + mode); break;

  case 3:
    if (mode > 0)
      coeffAbs = ((coeffAbs << 1) + (mode - 1));
    coeffAbs = ((coeffAbs << 1) + (mode > 0));
    coeff = coeffSign * coeffAbs;
    break;

  case 2: coeff = coeffSign * ((coeffAbs << 1) + mode); break;

  default: assert(mode == 0);
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
  int zeroRunAcc = 0;
  std::vector<int> zerorun;
  zerorun.reserve(pointCount);
  std::vector<uint32_t> residual;
  residual.resize(pointCount);

  int quantLayer = 0;

  std::vector<int64_t> quantWeights;
  computeQuantizationWeights(
    _lods.predictors, quantWeights, aps.quant_neigh_weight);

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }
    const uint32_t pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);
    auto& predictor = _lods.predictors[predictorIndex];
    predictor.predMode = 0;

    bool predModeEligible =
      predModeEligibleRefl(desc, aps, pointCloud, _lods.indexes, predictor);
    if (predModeEligible)
      decidePredModeRefl(
        desc, aps, pointCloud, _lods.indexes, predictorIndex, predictor,
        encoder, context, quant[0]);

    const uint64_t reflectance = pointCloud.getReflectance(pointIndex);
    const attr_t predictedReflectance =
      predictor.predictReflectance(pointCloud, _lods.indexes);
    const int64_t quantAttValue = reflectance;
    const int64_t quantPredAttValue = predictedReflectance;

    int64_t qStep = quant[0].stepSize();
    int64_t weight =
      std::min(quantWeights[predictorIndex], qStep) >> kFixedPointWeightShift;
    const int64_t delta = quant[0].quantize(
      ((quantAttValue - quantPredAttValue) * weight)
      << kFixedPointAttributeShift);
    int32_t attValue0 = delta;
    int64_t reconstructedDelta =
      divExp2RoundHalfUp(quant[0].scale(delta), kFixedPointAttributeShift);
    reconstructedDelta /= weight;

    if (predModeEligible)
      encodePredModeRefl(aps, predictor.predMode, attValue0);

    const int64_t reconstructedQuantAttValue =
      quantPredAttValue + reconstructedDelta;
    const attr_t reconstructedReflectance =
      attr_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));

    if (!attValue0)
      ++zeroRunAcc;
    else {
      zerorun.push_back(zeroRunAcc);
      zeroRunAcc = 0;
    }
    residual[predictorIndex] = attValue0;
    pointCloud.setReflectance(pointIndex, reconstructedReflectance);
    encoder.resStatUpdateRefl(attValue0);
  }
  if (zeroRunAcc)
    zerorun.push_back(zeroRunAcc);

  int runIdx = 0;
  int zeroRunRem = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (--zeroRunRem < 0) {
      zeroRunRem = zerorun[runIdx++];
      encoder.encodeRunLength(zeroRunRem);
    }

    if (!zeroRunRem)
      encoder.encode(residual[predictorIndex]);
  }
}

//----------------------------------------------------------------------------

Vec3<int64_t>
AttributeEncoder::computeColorResiduals(
  const AttributeParameterSet& aps,
  const Vec3<attr_t> color,
  const Vec3<attr_t> predictedColor,
  const Vec3<int8_t> icpCoeff,
  const Quantizers& quant)
{
  Vec3<int64_t> residuals;
  const int64_t quantAttValue = color[0];
  const int64_t quantPredAttValue = predictedColor[0];
  const int64_t delta = quant[0].quantize(
    (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
  residuals[0] = delta;
  const int64_t residual0 =
    divExp2RoundHalfUp(quant[0].scale(delta), kFixedPointAttributeShift);
  for (size_t k = 1; k < 3; ++k) {
    const int64_t quantAttValue = color[k];
    const int64_t quantPredAttValue = predictedColor[k];

    if (aps.inter_component_prediction_enabled_flag) {
      auto err = quantAttValue - quantPredAttValue
        - ((icpCoeff[k] * residual0 + 2) >> 2);

      auto delta = quant[1].quantize(err << kFixedPointAttributeShift);
      residuals[k] = delta;
    } else {
      const int64_t delta = quant[1].quantize(
        (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
      residuals[k] = delta;
    }
  }
  return residuals;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::decidePredModeColor(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexesLOD,
  const uint32_t predictorIndex,
  PCCPredictor& predictor,
  PCCResidualsEncoder& encoder,
  PCCResidualsEntropyEstimator& context,
  const Vec3<int8_t>& icpCoeff,
  const Quantizers& quant)
{
  Vec3<attr_t> attrValue = pointCloud.getColor(indexesLOD[predictorIndex]);

  // base case: weighted average of n neighbours
  int startpredIndex = aps.direct_avg_predictor_disabled_flag;
  predictor.predMode = startpredIndex;
  Vec3<attr_t> attrPred = predictor.predictColor(pointCloud, indexesLOD);
  Vec3<int64_t> attrResidualQuant =
    computeColorResiduals(aps, attrValue, attrPred, icpCoeff, quant);
  auto attrDistortion =
    computeColorDistortions(desc, attrValue, attrPred, quant);

  double rate = encoder.bitsPtColor(attrResidualQuant, 0);
  double best_score = attrDistortion
    + rate * kAttrPredLambdaC
      * (quant[0].stepSize() >> kFixedPointAttributeShift);

  for (int i = startpredIndex; i < predictor.neighborCount; i++) {
    if (i == aps.max_num_direct_predictors)
      break;

    attrPred =
      pointCloud.getColor(indexesLOD[predictor.neighbors[i].predictorIndex]);
    attrResidualQuant =
      computeColorResiduals(aps, attrValue, attrPred, icpCoeff, quant);
    attrDistortion = computeColorDistortions(desc, attrValue, attrPred, quant);

    int sigIdx = i + !aps.direct_avg_predictor_disabled_flag;
    double rate = encoder.bitsPtColor(attrResidualQuant, sigIdx);
    double score = attrDistortion
      + rate * kAttrPredLambdaC
        * (quant[0].stepSize() >> kFixedPointAttributeShift);

    if (score < best_score) {
      best_score = score;
      predictor.predMode = i + 1;
      // NB: setting predictor.neighborCount = 1 will cause issues
      // with reconstruction.
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodePredModeColor(
  const AttributeParameterSet& aps, int predMode, Vec3<int32_t>& coeff)
{
  int signk1 = coeff[1] < 0 ? -1 : 1;
  int signk2 = coeff[2] < 0 ? -1 : 1;
  int coeffAbsk1 = abs(coeff[1]);
  int coeffAbsk2 = abs(coeff[2]);

  int mode = predMode - aps.direct_avg_predictor_disabled_flag;
  int maxcand =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;
  assert(mode < maxcand);
  switch (maxcand) {
    int parityk1, parityk2;
  case 4:
    parityk1 = mode >> 1;
    parityk2 = mode & 1;
    coeff[1] = signk1 * ((coeffAbsk1 << 1) + parityk1);
    coeff[2] = signk2 * ((coeffAbsk2 << 1) + parityk2);
    break;

  case 3:
    parityk1 = mode ? 1 : 0;
    coeff[1] = signk1 * ((coeffAbsk1 << 1) + parityk1);
    if (parityk1) {
      parityk2 = mode - parityk1;
      coeff[2] = signk2 * ((coeffAbsk2 << 1) + parityk2);
    }
    break;

  case 2:
    parityk1 = mode;
    coeff[1] = signk1 * ((coeffAbsk1 << 1) + parityk1);
    break;

  default: assert(mode == 0);
  }
}

//----------------------------------------------------------------------------

std::vector<Vec3<int8_t>>
AttributeEncoder::computeInterComponentPredictionCoeffs(
  const AttributeParameterSet& aps, const PCCPointSet3& pointCloud)
{
  int maxNumDetailLevels = aps.maxNumDetailLevels();
  assert(_lods.numPointsInLod.size() <= maxNumDetailLevels);

  // Two secondary colour components (positive sign set)
  // NB: k=0 is never used
  std::vector<Vec3<int8_t>> signs(maxNumDetailLevels, {0, 1, 1});

  // Estimate residual using original neighbour as predictor
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<Vec3<int32_t>> residual(pointCount);

  for (size_t predIdx = 0; predIdx < pointCount; ++predIdx) {
    const auto pointIdx = _lods.indexes[predIdx];
    auto& predictor = _lods.predictors[predIdx];

    // taking first neighbor for simplicity
    predictor.predMode = 1;
    auto predAttr = predictor.predictColor(pointCloud, _lods.indexes);
    auto srcAttr = pointCloud.getColor(pointIdx);
    residual[predIdx] = Vec3<int>(srcAttr) - Vec3<int>(predAttr);

    // reset is needed as RD would be done later.
    predictor.predMode = 0;
  }

  const int nWeights = 8;
  const int nShift = 2;  // from log2(nWeights >> 1)
  std::vector<Vec3<int64_t>> sumPredCoeff(nWeights, 0);
  Vec3<int64_t> sumOrigCoeff = 0;

  int lod = 0;
  for (size_t predIdx = 0; predIdx < pointCount; ++predIdx) {
    Vec3<int32_t> resid = residual[predIdx];

    for (int w = 0; w < nWeights; w++) {
      for (int k = 1; k < 3; k++)
        sumPredCoeff[w][k] +=
          abs(resid[k] - signs[lod][k] * (((w + 1) * resid[0] + 2) >> nShift));
    }

    for (int k = 1; k < 3; k++)
      sumOrigCoeff[k] += abs(resid[k]);

    // at LoD transition, determine the sign coeff
    if (predIdx != _lods.numPointsInLod[lod] - 1)
      continue;

    // find the best weight
    for (int k = 1; k < 3; k++) {
      auto best = std::min_element(
        sumPredCoeff.begin(), sumPredCoeff.end(),
        [=](Vec3<int64_t>& a, Vec3<int64_t>& b) { return a[k] < b[k]; });

      int coeff = 1 + std::distance(sumPredCoeff.begin(), best);
      signs[lod][k] *= coeff;

      assert(signs[lod][k] < nWeights + 1 && signs[lod][k] > -(nWeights + 1));

      if ((*best)[k] > sumOrigCoeff[k])
        signs[lod][k] = 0;
    }

    for (int w = 0; w < nWeights; w++)
      sumPredCoeff[w] = 0;
    sumOrigCoeff = 0;
    lod++;
  }

  // NB: there may be more coefficients than actual detail levels
  // Set any unused detail level coefficients to 0
  for (; lod < maxNumDetailLevels; lod++)
    signs[lod] = 0;

  return signs;
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

  int64_t clipMax = (1 << desc.bitdepth) - 1;
  Vec3<int32_t> values;
  PCCResidualsEntropyEstimator context;
  int zeroRunAcc = 0;
  std::vector<int> zerorun;
  std::vector<int32_t> residual[3];
  for (int i = 0; i < 3; i++) {
    residual[i].resize(pointCount);
  }

  bool icpPresent = _abh->icpPresent(desc, aps);
  if (icpPresent)
    _abh->icpCoeffs = computeInterComponentPredictionCoeffs(aps, pointCloud);
  auto icpCoeff = icpPresent ? _abh->icpCoeffs[0] : 0;

  int lod = 0;
  int quantLayer = 0;

  std::vector<int64_t> quantWeights;
  computeQuantizationWeights(
    _lods.predictors, quantWeights, aps.quant_neigh_weight);

  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (predictorIndex == _lods.numPointsInLod[quantLayer]) {
      quantLayer = std::min(int(qpSet.layers.size()) - 1, quantLayer + 1);
    }

    if (icpPresent && predictorIndex == _lods.numPointsInLod[lod])
      icpCoeff = _abh->icpCoeffs[++lod];

    const auto pointIndex = _lods.indexes[predictorIndex];
    auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);
    auto& predictor = _lods.predictors[predictorIndex];
    predictor.predMode = 0;

    bool predModeEligible =
      predModeEligibleColor(desc, aps, pointCloud, _lods.indexes, predictor);
    if (predModeEligible)
      decidePredModeColor(
        desc, aps, pointCloud, _lods.indexes, predictorIndex, predictor,
        encoder, context, icpCoeff, quant);

    const Vec3<attr_t> color = pointCloud.getColor(pointIndex);
    const Vec3<attr_t> predictedColor =
      predictor.predictColor(pointCloud, _lods.indexes);

    Vec3<attr_t> reconstructedColor;
    int64_t residual0 = 0;
    for (int k = 0; k < 3; ++k) {
      const auto& q = quant[std::min(k, 1)];
      int64_t residual = color[k] - predictedColor[k];

      int64_t qStep = q.stepSize();
      int64_t weight = std::min(quantWeights[predictorIndex], qStep)
        >> kFixedPointWeightShift;
      int64_t residualQ =
        q.quantize((residual * weight) << kFixedPointAttributeShift);
      int64_t residualR =
        divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);
      residualR /= weight;

      if (aps.inter_component_prediction_enabled_flag && k > 0) {
        residual = residual - ((icpCoeff[k] * residual0 + 2) >> 2);
        residualQ =
          q.quantize((residual * weight) << kFixedPointAttributeShift);
        residualR =
          divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);
        residualR /= weight;
        residualR += ((icpCoeff[k] * residual0 + 2) >> 2);
      }

      if (k == 0)
        residual0 = residualR;

      values[k] = residualQ;

      int64_t recon = predictedColor[k] + residualR;
      reconstructedColor[k] = attr_t(PCCClip(recon, int64_t(0), clipMax));
    }

    if (predModeEligible)
      encodePredModeColor(aps, predictor.predMode, values);

    pointCloud.setColor(pointIndex, reconstructedColor);
    encoder.resStatUpdateColor(values);

    if (!values[0] && !values[1] && !values[2]) {
      ++zeroRunAcc;
    } else {
      zerorun.push_back(zeroRunAcc);
      zeroRunAcc = 0;
    }

    for (int i = 0; i < 3; i++) {
      residual[i][predictorIndex] = values[i];
    }
  }
  if (zeroRunAcc)
    zerorun.push_back(zeroRunAcc);

  int runIdx = 0;
  int zeroRunRem = 0;
  for (size_t predictorIndex = 0; predictorIndex < pointCount;
       ++predictorIndex) {
    if (--zeroRunRem < 0) {
      zeroRunRem = zerorun[runIdx++];
      encoder.encodeRunLength(zeroRunRem);
    }

    if (!zeroRunRem) {
      for (size_t k = 0; k < 3; k++)
        values[k] = residual[k][predictorIndex];

      encoder.encode(values[0], values[1], values[2]);
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
  int zeroRun = 0;
  for (int n = 0; n < voxelCount; ++n) {
    auto value = coefficients[n];
    if (!value)
      ++zeroRun;
    else {
      encoder.encodeRunLength(zeroRun);
      encoder.encode(value);
      zeroRun = 0;
    }
  }
  if (zeroRun)
    encoder.encodeRunLength(zeroRun);

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
  int zeroRun = 0;
  for (int n = 0; n < voxelCount; ++n) {
    for (int d = 0; d < attribCount; ++d) {
      values[d] = coefficients[voxelCount * d + n];
    }
    if (!values[0] && !values[1] && !values[2])
      ++zeroRun;
    else {
      encoder.encodeRunLength(zeroRun);
      encoder.encode(values[0], values[1], values[2]);
      zeroRun = 0;
    }
  }
  if (zeroRun)
    encoder.encodeRunLength(zeroRun);

  int clipMax = (1 << desc.bitdepth) - 1;
  for (int n = 0; n < voxelCount; n++) {
    const int r = attributes[attribCount * n];
    const int g = attributes[attribCount * n + 1];
    const int b = attributes[attribCount * n + 2];
    Vec3<attr_t> color;
    color[0] = attr_t(PCCClip(r, 0, clipMax));
    color[1] = attr_t(PCCClip(g, 0, clipMax));
    color[2] = attr_t(PCCClip(b, 0, clipMax));
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

  // Per level-of-detail coefficients for last component prediction
  int8_t lastCompPredCoeff = 0;
  if (aps.last_component_prediction_enabled_flag) {
    _abh->attrLcpCoeffs = computeLastComponentPredictionCoeff(aps, colors);
    lastCompPredCoeff = _abh->attrLcpCoeffs[0];
  }

  int zeroRun = 0;
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
        lastCompPredCoeff = _abh->attrLcpCoeffs[lod];
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

    color[2] -= (lastCompPredCoeff * color[1]) >> 2;
    scaled *= lastCompPredCoeff;
    scaled >>= 2;

    values[2] = quant[1].quantize(color[2] * quantWeight);
    scaled += quant[1].scale(values[2]);
    color[2] = divExp2RoundHalfInf(scaled * iQuantWeight, 40);

    if (!values[0] && !values[1] && !values[2])
      ++zeroRun;
    else {
      encoder.encodeRunLength(zeroRun);
      encoder.encode(values[0], values[1], values[2]);
      zeroRun = 0;
    }
  }
  if (zeroRun)
    encoder.encodeRunLength(zeroRun);

  // reconstruct
  for (size_t lodIndex = 1; lodIndex < lodCount; ++lodIndex) {
    const size_t startIndex = _lods.numPointsInLod[lodIndex - 1];
    const size_t endIndex = _lods.numPointsInLod[lodIndex];
    PCCLiftUpdate(
      _lods.predictors, weights, startIndex, endIndex, false, colors);
    PCCLiftPredict(_lods.predictors, startIndex, endIndex, false, colors);
  }

  int64_t clipMax = (1 << desc.bitdepth) - 1;
  for (size_t f = 0; f < pointCount; ++f) {
    const auto color0 =
      divExp2RoundHalfInf(colors[f], kFixedPointAttributeShift);
    Vec3<attr_t> color;
    for (size_t d = 0; d < 3; ++d) {
      color[d] = attr_t(PCCClip(color0[d], 0, clipMax));
    }
    pointCloud.setColor(_lods.indexes[f], color);
  }
}

//----------------------------------------------------------------------------

std::vector<int8_t>
AttributeEncoder::computeLastComponentPredictionCoeff(
  const AttributeParameterSet& aps, const std::vector<Vec3<int64_t>>& coeffs)
{
  int maxNumDetailLevels = aps.maxNumDetailLevels();
  assert(_lods.numPointsInLod.size() <= maxNumDetailLevels);
  std::vector<int8_t> signs(maxNumDetailLevels, 0);

  int64_t sumk1k2 = 0;
  int64_t sumk1k1 = 0;
  int lod = 0;
  for (size_t coeffIdx = 0; coeffIdx < coeffs.size(); ++coeffIdx) {
    auto& attr = coeffs[coeffIdx];
    int mult = attr[1] * attr[2];
    int mult2 = attr[1] * attr[1];
    sumk1k2 += mult;
    sumk1k1 += mult2;

    // compute prediction coefficient at end of detail level
    if (coeffIdx != _lods.numPointsInLod[lod] - 1)
      continue;

    int scale = 0;
    if (sumk1k2 && sumk1k1) {
      // sign(sumk1k2) * sign(sumk1k1)
      int sign = (sumk1k2 < 0) ^ (sumk1k1 < 0) ? -1 : 1;
      scale = ((sumk1k2 << 2) + sign * (sumk1k1 >> 1)) / sumk1k1;
    }
    sumk1k2 = sumk1k1 = 0;

    // NB: coding range is limited to +-8
    signs[lod] = PCCClip(scale, -8, 8);
    lod++;
  }

  // NB: there may be more coefficients than actual detail levels
  // Propagate the last value to all unused levels to minimise useless cost
  for (; lod < maxNumDetailLevels; lod++)
    signs[lod] = signs[lod - 1];

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
  int zeroRun = 0;
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
      ++zeroRun;
    else {
      encoder.encodeRunLength(zeroRun);
      encoder.encode(detail);
      zeroRun = 0;
    }
  }
  if (zeroRun)
    encoder.encodeRunLength(zeroRun);

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

int
AttributeEncoder::computeColorDistortions(
  const AttributeDescription& desc,
  const Vec3<attr_t> color,
  const Vec3<attr_t> predictedColor,
  const Quantizers& quant)
{
  int64_t clipMax = (1 << desc.bitdepth) - 1;

  Vec3<attr_t> reconstructedColor;
  for (int k = 0; k < 3; ++k) {
    const auto& q = quant[std::min(k, 1)];
    int64_t residual = color[k] - predictedColor[k];

    int64_t residualQ = q.quantize(residual << kFixedPointAttributeShift);
    int64_t residualR =
      divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);

    int64_t recon = predictedColor[k] + residualR;
    reconstructedColor[k] = attr_t(PCCClip(recon, int64_t(0), clipMax));
  }

  int distortion = 0;
  for (int k = 0; k < 3; ++k)
    distortion += std::abs(color[k] - reconstructedColor[k]);

  return distortion;
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
