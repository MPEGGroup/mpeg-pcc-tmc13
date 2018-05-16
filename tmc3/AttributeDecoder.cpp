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

namespace pcc {

//============================================================================
// An encapsulation of the entropy decoding methods used in attribute coding

struct PCCResidualsDecoder {
  uint32_t alphabetSize;
  o3dgc::Arithmetic_Codec arithmeticDecoder;
  o3dgc::Static_Bit_Model binaryModel0;
  o3dgc::Adaptive_Bit_Model binaryModelDiff0;
  o3dgc::Adaptive_Data_Model multiSymbolModelDiff0;
  o3dgc::Adaptive_Bit_Model binaryModelDiff1;
  o3dgc::Adaptive_Data_Model multiSymbolModelDiff1;

  PCCResidualsDecoder() { alphabetSize = 0; }

  void start(PCCBitstream &bitstream, const uint32_t alphabetSize = 64);
  void stop();
  uint32_t decode0();
  uint32_t decode1();
};

//----------------------------------------------------------------------------

void PCCResidualsDecoder::start(
  PCCBitstream &bitstream, const uint32_t alphabetSize
) {
  this->alphabetSize = alphabetSize;
  multiSymbolModelDiff0.set_alphabet(alphabetSize + 1);
  binaryModelDiff0.reset();
  multiSymbolModelDiff1.set_alphabet(alphabetSize + 1);
  binaryModelDiff1.reset();
  arithmeticDecoder.set_buffer(
    static_cast<uint32_t>(bitstream.capacity - bitstream.size),
    bitstream.buffer + bitstream.size
  );
  arithmeticDecoder.start_decoder();
}

//----------------------------------------------------------------------------

void PCCResidualsDecoder::stop() {
  arithmeticDecoder.stop_decoder();
}

//----------------------------------------------------------------------------

uint32_t PCCResidualsDecoder::decode0() {
  uint32_t value = arithmeticDecoder.decode(multiSymbolModelDiff0);
  if (value == alphabetSize) {
    value += arithmeticDecoder.ExpGolombDecode(
      0, binaryModel0, binaryModelDiff0
    );
  }
  return value;
}

//----------------------------------------------------------------------------

uint32_t PCCResidualsDecoder::decode1() {
  uint32_t value = arithmeticDecoder.decode(multiSymbolModelDiff1);
  if (value == alphabetSize) {
    value += arithmeticDecoder.ExpGolombDecode(
      0, binaryModel0, binaryModelDiff1
    );
  }
  return value;
}

//============================================================================
// AttributeDecoder Members

void AttributeDecoder::buildPredictors(const PCCPointSet3 &pointCloud) {
  std::vector<uint32_t> numberOfPointsPerLOD;
  std::vector<uint32_t> indexes;
  PCCBuildPredictors(
    pointCloud, numberOfNearestNeighborsInPrediction, levelOfDetailCount,
    dist2, predictors, numberOfPointsPerLOD, indexes
  );

  for (auto &predictor : predictors) {
    predictor.computeWeights(numberOfNearestNeighborsInPrediction);
  }
}

//----------------------------------------------------------------------------

int AttributeDecoder::decodeHeader(
  const std::string &attributeName,
  PCCBitstream &bitstream
) {
  uint8_t numberOfNearestNeighborsCount = 0;
  PCCReadFromBuffer<uint8_t>(bitstream.buffer, numberOfNearestNeighborsCount, bitstream.size);
  numberOfNearestNeighborsInPrediction = numberOfNearestNeighborsCount;

  uint8_t lodCount = 0;
  PCCReadFromBuffer<uint8_t>(bitstream.buffer, lodCount, bitstream.size);
  levelOfDetailCount = lodCount;

  dist2.resize(levelOfDetailCount);
  for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
    uint32_t d2 = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, d2, bitstream.size);
    dist2[lodIndex] = d2;
  }

  quantizationSteps.resize(levelOfDetailCount);
  for (size_t lodIndex = 0; lodIndex < levelOfDetailCount; ++lodIndex) {
    uint32_t qs = 0;
    PCCReadFromBuffer<uint32_t>(bitstream.buffer, qs, bitstream.size);
    quantizationSteps[lodIndex] = qs;
  }

  return 0;
}

//----------------------------------------------------------------------------

int AttributeDecoder::decodeReflectances(
  PCCBitstream &bitstream,
  PCCPointSet3 &pointCloud
) {
  uint32_t compressedBitstreamSize = 0;
  PCCReadFromBuffer<uint32_t>(bitstream.buffer, compressedBitstreamSize, bitstream.size);
  PCCResidualsDecoder decoder;
  const uint32_t alphabetSize = 64;
  decoder.start(bitstream, alphabetSize);

  const size_t pointCount = predictors.size();
  for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
    auto &predictor = predictors[predictorIndex];
    uint16_t &reflectance = pointCloud.getReflectance(predictor.index);
    const size_t lodIndex = predictor.levelOfDetailIndex;
    const int64_t qs = quantizationSteps[lodIndex];
    const uint32_t attValue0 = decoder.decode0();
    const int64_t quantPredAttValue = predictor.predictReflectance(pointCloud);
    const int64_t delta = PCCInverseQuantization(o3dgc::UIntToInt(attValue0), qs);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
    reflectance = uint16_t(PCCClip(reconstructedQuantAttValue, int64_t(0),
                                   int64_t(std::numeric_limits<uint16_t>::max())));
  }
  decoder.stop();
  bitstream.size += compressedBitstreamSize;
  return 0;
}

//----------------------------------------------------------------------------

int AttributeDecoder::decodeColors(
  PCCBitstream &bitstream,
  PCCPointSet3 &pointCloud
) {
  uint32_t compressedBitstreamSize = 0;
  PCCReadFromBuffer<uint32_t>(bitstream.buffer, compressedBitstreamSize, bitstream.size);
  PCCResidualsDecoder decoder;
  const uint32_t alphabetSize = 64;
  decoder.start(bitstream, alphabetSize);

  const size_t pointCount = predictors.size();
  for (size_t predictorIndex = 0; predictorIndex < pointCount; ++predictorIndex) {
    auto &predictor = predictors[predictorIndex];
    PCCColor3B &color = pointCloud.getColor(predictor.index);
    const PCCColor3B predictedColor = predictor.predictColor(pointCloud);
    const size_t lodIndex = predictor.levelOfDetailIndex;
    const int64_t qs = quantizationSteps[lodIndex];
    const uint32_t attValue0 = decoder.decode0();
    const int64_t quantPredAttValue = predictedColor[0];
    const int64_t delta = PCCInverseQuantization(o3dgc::UIntToInt(attValue0), qs);
    const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
    color[0] = uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));

    for (size_t k = 1; k < 3; ++k) {
      const uint32_t attValue = decoder.decode1();
      const int64_t quantPredAttValue = predictedColor[k];
      const int64_t delta = PCCInverseQuantization(o3dgc::UIntToInt(attValue), qs);
      const int64_t reconstructedQuantAttValue = quantPredAttValue + delta;
      color[k] = uint8_t(PCCClip(reconstructedQuantAttValue, int64_t(0), int64_t(255)));
    }
  }
  decoder.stop();
  bitstream.size += compressedBitstreamSize;
  return 0;
}

//============================================================================

} /* namespace pcc */
