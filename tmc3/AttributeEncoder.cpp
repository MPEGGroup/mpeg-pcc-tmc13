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

namespace pcc {

//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

struct PCCResidualsEncoder {
  uint32_t alphabetSize;
  o3dgc::Arithmetic_Codec arithmeticEncoder;
  o3dgc::Static_Bit_Model binaryModel0;
  o3dgc::Adaptive_Bit_Model binaryModelDiff0;
  o3dgc::Adaptive_Data_Model multiSymbolModelDiff0;
  o3dgc::Adaptive_Bit_Model binaryModelDiff1;
  o3dgc::Adaptive_Data_Model multiSymbolModelDiff1;

  PCCResidualsEncoder() { alphabetSize = 0; }

  void start(PCCBitstream &bitstream, const uint32_t alphabetSize = 64);
  uint32_t stop();
  inline void encode0(const uint32_t value);
  inline void encode1(const uint32_t value);
};

//----------------------------------------------------------------------------

void PCCResidualsEncoder::start(
  PCCBitstream &bitstream, const uint32_t alphabetSize
) {
  this->alphabetSize = alphabetSize;
  multiSymbolModelDiff0.set_alphabet(alphabetSize + 1);
  binaryModelDiff0.reset();
  multiSymbolModelDiff1.set_alphabet(alphabetSize + 1);
  binaryModelDiff1.reset();
  arithmeticEncoder.set_buffer(
    static_cast<uint32_t>(bitstream.capacity - bitstream.size),
    bitstream.buffer + bitstream.size
  );
  arithmeticEncoder.start_encoder();
}

//----------------------------------------------------------------------------

uint32_t PCCResidualsEncoder::stop() {
  return arithmeticEncoder.stop_encoder();
}

//----------------------------------------------------------------------------

inline void PCCResidualsEncoder::encode0(const uint32_t value) {
  if (value < alphabetSize) {
    arithmeticEncoder.encode(value, multiSymbolModelDiff0);
  } else {
    arithmeticEncoder.encode(alphabetSize, multiSymbolModelDiff0);
    arithmeticEncoder.ExpGolombEncode(
      value - alphabetSize, 0, binaryModel0, binaryModelDiff0);
  }
}

//----------------------------------------------------------------------------

inline void PCCResidualsEncoder::encode1(const uint32_t value) {
  if (value < alphabetSize) {
    arithmeticEncoder.encode(value, multiSymbolModelDiff1);
  } else {
    arithmeticEncoder.encode(alphabetSize, multiSymbolModelDiff1);
    arithmeticEncoder.ExpGolombEncode(
      value - alphabetSize, 0, binaryModel0, binaryModelDiff1);
  }
}

//============================================================================
// AttributeEncoder Members

void AttributeEncoder::buildPredictors(
  const PCCAttributeEncodeParamaters &attributeParams,
  const PCCPointSet3 &pointCloud
) {
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexes;
  PCCBuildPredictors(
    pointCloud, attributeParams.numberOfNearestNeighborsInPrediction,
    attributeParams.levelOfDetailCount, attributeParams.dist2,
    predictors, numberOfPointsPerLOD, indexes);

  for (auto &predictor : predictors) {
    predictor.computeWeights(
      attributeParams.numberOfNearestNeighborsInPrediction);
  }
}

//----------------------------------------------------------------------------

int AttributeEncoder::encodeHeader(
  const PCCAttributeEncodeParamaters &attributeParams,
  const std::string &attributeName,
  PCCBitstream &bitstream
) const {
  PCCWriteToBuffer<uint8_t>(
    uint8_t(attributeParams.numberOfNearestNeighborsInPrediction),
    bitstream.buffer, bitstream.size);

  PCCWriteToBuffer<uint8_t>(
    uint8_t(attributeParams.levelOfDetailCount),
    bitstream.buffer, bitstream.size);

  for (size_t lodIndex = 0; lodIndex < attributeParams.levelOfDetailCount; ++lodIndex) {
    const uint32_t d2 = uint32_t(attributeParams.dist2[lodIndex]);
    PCCWriteToBuffer<uint32_t>(d2, bitstream.buffer, bitstream.size);
  }
  for (size_t lodIndex = 0; lodIndex < attributeParams.levelOfDetailCount; ++lodIndex) {
    const uint32_t qs = uint32_t(attributeParams.quantizationSteps[lodIndex]);
    PCCWriteToBuffer<uint32_t>(qs, bitstream.buffer, bitstream.size);
  }
  return 0;
}

//----------------------------------------------------------------------------

int AttributeEncoder::encodeReflectances(
  const PCCAttributeEncodeParamaters &reflectanceParams,
  PCCPointSet3 &pointCloud,
  PCCBitstream &bitstream
) {
  uint64_t startSize = bitstream.size;
  bitstream.size += 4;  // placehoder for bitstream size
  PCCResidualsEncoder encoder;
  const uint32_t alphabetSize = 64;
  encoder.start(bitstream, alphabetSize);

  const size_t pointCount = predictors.size();
  for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
    const auto &predictor = predictors[predictorIndex];
    const size_t lodIndex = predictor.levelOfDetailIndex;
    const int64_t qs = reflectanceParams.quantizationSteps[lodIndex];

    const int64_t quantAttValue = pointCloud.getReflectance(predictor.index);
    const int64_t quantPredAttValue = predictor.predictReflectance(pointCloud);
    const int64_t delta = PCCQuantization(quantAttValue - quantPredAttValue, qs);
    const uint32_t attValue0 = uint32_t(o3dgc::IntToUInt(long(delta)));
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + reconstructedDelta;
    const uint16_t reconstructedReflectance = uint16_t(PCCClip(
        reconstructedQuantAttValue, int64_t(0), int64_t(std::numeric_limits<uint16_t>::max())));
    encoder.encode0(attValue0);
    pointCloud.setReflectance(predictor.index, reconstructedReflectance);
  }
  uint32_t compressedBitstreamSize = encoder.stop();
  bitstream.size += compressedBitstreamSize;
  PCCWriteToBuffer<uint32_t>(compressedBitstreamSize, bitstream.buffer, startSize);
  return 0;
}

//----------------------------------------------------------------------------

int AttributeEncoder::encodeColors(
  const PCCAttributeEncodeParamaters &colorParams,
  PCCPointSet3 &pointCloud,
  PCCBitstream &bitstream
) {
  uint64_t startSize = bitstream.size;
  bitstream.size += 4;  // placehoder for bitstream size
  PCCResidualsEncoder encoder;
  const uint32_t alphabetSize = 64;
  encoder.start(bitstream, alphabetSize);

  const size_t pointCount = predictors.size();
  for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
    const auto &predictor = predictors[predictorIndex];
    const PCCColor3B color = pointCloud.getColor(predictor.index);
    const PCCColor3B predictedColor = predictor.predictColor(pointCloud);
    const size_t lodIndex = predictor.levelOfDetailIndex;
    const int64_t qs = colorParams.quantizationSteps[lodIndex];
    const int64_t quantAttValue = color[0];
    const int64_t quantPredAttValue = predictedColor[0];
    const int64_t delta = PCCQuantization(quantAttValue - quantPredAttValue, qs);
    const uint32_t attValue0 = uint32_t(o3dgc::IntToUInt(long(delta)));
    const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + reconstructedDelta;
    encoder.encode0(attValue0);
    PCCColor3B reconstructedColor;
    reconstructedColor[0] =
        uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));
    for (size_t k = 1; k < 3; ++k) {
      const int64_t quantAttValue = color[k];
      const int64_t quantPredAttValue = predictedColor[k];
      const int64_t delta = PCCQuantization(quantAttValue - quantPredAttValue, qs);
      const int64_t reconstructedDelta = PCCInverseQuantization(delta, qs);
      const int64_t reconstructedQuantAttValue = quantPredAttValue + reconstructedDelta;
      const uint32_t attValue1 = uint32_t(o3dgc::IntToUInt(long(delta)));
      encoder.encode1(attValue1);
      reconstructedColor[k] =
          uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));
    }
    pointCloud.setColor(predictor.index, reconstructedColor);
  }
  uint32_t compressedBitstreamSize = encoder.stop();
  bitstream.size += compressedBitstreamSize;
  PCCWriteToBuffer<uint32_t>(compressedBitstreamSize, bitstream.buffer, startSize);
  return 0;
}

//============================================================================

} /* namespace pcc */
