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
#include "RAHT.h"

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

  uint8_t transType;
  PCCReadFromBuffer<uint8_t>(bitstream.buffer, transType, bitstream.size);
  transformType = TransformType(transType);

  PCCReadFromBuffer<uint8_t>(bitstream.buffer, depthRaht, bitstream.size);
  PCCReadFromBuffer<uint8_t>(bitstream.buffer, binaryLevelThresholdRaht, bitstream.size);
  PCCReadFromBuffer<uint32_t>(bitstream.buffer, quantizationStepRaht, bitstream.size);
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

  switch (transformType) {
  case TransformType::kRAHT:
    decodeReflectancesRaht(decoder, pointCloud);
    break;

  case TransformType::kIntegerLift:
    decodeReflectancesIntegerLift(decoder, pointCloud);
    break;
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

  switch (transformType) {
  case TransformType::kRAHT:
    decodeColorsRaht(decoder, pointCloud);
    break;

  case TransformType::kIntegerLift:
    decodeColorsIntegerLift(decoder, pointCloud);
    break;
  }

  decoder.stop();
  bitstream.size += compressedBitstreamSize;
  return 0;
}

//----------------------------------------------------------------------------

void AttributeDecoder::decodeReflectancesIntegerLift(
  PCCResidualsDecoder &decoder,
  PCCPointSet3 &pointCloud
) {
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
}

//----------------------------------------------------------------------------

void AttributeDecoder::decodeColorsIntegerLift(
  PCCResidualsDecoder &decoder,
  PCCPointSet3 &pointCloud
) {
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
}

//----------------------------------------------------------------------------

void AttributeDecoder::decodeReflectancesRaht(
  PCCResidualsDecoder &decoder,
  PCCPointSet3 &pointCloud
) {
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    const auto position = pointCloud[n];
    int x = int(position[0]);
    int y = int(position[1]);
    int z = int(position[2]);
    long long mortonCode = 0;
    for (int b = 0; b < depthRaht; b++) {
      mortonCode |= (long long)((x >> b) & 1) << (3 * b + 2);
      mortonCode |= (long long)((y >> b) & 1) << (3 * b + 1);
      mortonCode |= (long long)((z >> b) & 1) << (3 * b);
    }
    packedVoxel[n].mortonCode = mortonCode;
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Moroton codes
  long long *mortonCode = new long long[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
  }

  // Re-obtain weights at the decoder by calling RAHT without any attributes.
  float *weight = new float[voxelCount];
  int *binaryLayer = new int[voxelCount];
  regionAdaptiveHierarchicalTransform(
      mortonCode, nullptr, weight, binaryLayer, 0, voxelCount, depthRaht);

  // Sort integerized attributes by weight
  std::vector<WeightWithIndex> sortedWeight(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    sortedWeight[n].weight = weight[n];
    sortedWeight[n].index = n;
  }
  sort(sortedWeight.begin(), sortedWeight.end());

  // Entropy decode
  int *sortedIntegerizedAttributes = new int[voxelCount];
  for (int n = 0; n < voxelCount; ++n) {
    const uint32_t attValue0 = decoder.decode0();
    sortedIntegerizedAttributes[n] = o3dgc::UIntToInt(attValue0);
  }

  // Unsort integerized attributes by weight.
  int *integerizedAttributes = new int[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    integerizedAttributes[sortedWeight[n].index] = sortedIntegerizedAttributes[n];
  }

  // Inverse Quantize.
  float *attributes = new float[voxelCount];
  const int qstep = int(quantizationStepRaht);
  for (int n = 0; n < voxelCount; n++) {
    attributes[n] = integerizedAttributes[n] * qstep;
  }

  regionAdaptiveHierarchicalInverseTransform(
      mortonCode, attributes, 1, voxelCount, depthRaht);

  const int maxReflectance = std::numeric_limits<uint16_t>::max();
  const int minReflectance = 0;
  for (int n = 0; n < voxelCount; n++) {
    const int reflectance = PCCClip((int)round(attributes[n]), minReflectance, maxReflectance);
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

void AttributeDecoder::decodeColorsRaht(
  PCCResidualsDecoder &decoder,
  PCCPointSet3 &pointCloud
) {
  const int voxelCount = int(pointCloud.getPointCount());
  std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    const auto position = pointCloud[n];
    int x = int(position[0]);
    int y = int(position[1]);
    int z = int(position[2]);
    long long mortonCode = 0;
    for (int b = 0; b < depthRaht; b++) {
      mortonCode |= (long long)((x >> b) & 1) << (3 * b + 2);
      mortonCode |= (long long)((y >> b) & 1) << (3 * b + 1);
      mortonCode |= (long long)((z >> b) & 1) << (3 * b);
    }
    packedVoxel[n].mortonCode = mortonCode;
    packedVoxel[n].index = n;
  }
  sort(packedVoxel.begin(), packedVoxel.end());

  // Moroton codes
  long long *mortonCode = new long long[voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    mortonCode[n] = packedVoxel[n].mortonCode;
  }

  // Re-obtain weights at the decoder by calling RAHT without any attributes.
  float *weight = new float[voxelCount];
  int *binaryLayer = new int[voxelCount];
  regionAdaptiveHierarchicalTransform(
      mortonCode, nullptr, weight, binaryLayer, 0, voxelCount, depthRaht);

  // Sort integerized attributes by weight
  std::vector<WeightWithIndex> sortedWeight(voxelCount);
  for (int n = 0; n < voxelCount; n++) {
    sortedWeight[n].weight = weight[n];
    sortedWeight[n].index = n;
  }
  sort(sortedWeight.begin(), sortedWeight.end());

  // Entropy decode
  const int attribCount = 3;
  int *sortedIntegerizedAttributes = new int[attribCount * voxelCount];
  for (int n = 0; n < voxelCount; ++n) {
    const uint32_t attValue0 = decoder.decode0();
    sortedIntegerizedAttributes[n] = o3dgc::UIntToInt(attValue0);
    if (binaryLayer[sortedWeight[n].index] >= binaryLevelThresholdRaht) {
      for (int d = 1; d < 3; ++d) {
        const uint32_t attValue1 = decoder.decode1();
        sortedIntegerizedAttributes[voxelCount * d + n] = o3dgc::UIntToInt(attValue1);
      }
    } else {
      for (int d = 1; d < 3; d++) {
        sortedIntegerizedAttributes[voxelCount * d + n] = 0;
      }
    }
  }

  // Unsort integerized attributes by weight.
  int *integerizedAttributes = new int[attribCount * voxelCount];
  for (int n = 0; n < voxelCount; n++) {
    for (int k = 0; k < attribCount; k++) {
      // Pull sorted integerized attributes out of column-major order.
      integerizedAttributes[attribCount * sortedWeight[n].index + k] =
          sortedIntegerizedAttributes[voxelCount * k + n];
    }
  }

  // Inverse Quantize.
  float *attributes = new float[attribCount * voxelCount];
  const int qstep = int(quantizationStepRaht);
  for (int n = 0; n < voxelCount; n++) {
    for (int k = 0; k < attribCount; k++) {
      attributes[attribCount * n + k] = integerizedAttributes[attribCount * n + k] * qstep;
    }
  }

  regionAdaptiveHierarchicalInverseTransform(
      mortonCode, attributes, attribCount, voxelCount, depthRaht);

  for (int n = 0; n < voxelCount; n++) {
    const int r = (int)round(attributes[attribCount * n]);
    const int g = (int)round(attributes[attribCount * n + 1]);
    const int b = (int)round(attributes[attribCount * n + 2]);
    PCCColor3B color;
    color[0] = uint8_t(PCCClip(r, 0, 255));
    color[1] = uint8_t(PCCClip(g, 0, 255));
    color[2] = uint8_t(PCCClip(b, 0, 255));
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

//============================================================================

} /* namespace pcc */
