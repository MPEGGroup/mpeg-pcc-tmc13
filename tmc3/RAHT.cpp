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

#include "RAHT.h"

#include <cstddef>
#include <utility>

namespace pcc {

//============================================================================

void
rahtFixedPointRotation(
  FixedPoint quantStepSizeLuma,
  uint64_t weightLeft,
  uint64_t weightRight,
  FixedPoint* attributeLeft,
  FixedPoint* attributeRight,
  FixedPoint* attributeTransformedLow,
  FixedPoint* attributeTransformedHigh,
  size_t attributeCount)
{
  FixedPoint b, adjustedQuantStepSize;

  b.val = (weightRight << b.kFracBits) / (weightLeft + weightRight);

  adjustedQuantStepSize.val = sqrtFixedpoint(
    ((quantStepSizeLuma.val * quantStepSizeLuma.val)
     * (weightLeft + weightRight))
    / (weightLeft * weightRight));

  while (attributeCount--) {
    *attributeTransformedHigh = *attributeRight;
    *attributeTransformedHigh -= *attributeLeft;
    *attributeTransformedLow = *attributeTransformedHigh;
    *attributeTransformedLow *= b;
    *attributeTransformedLow += *attributeLeft;
    *attributeTransformedHigh /= adjustedQuantStepSize;

    attributeLeft++;
    attributeRight++;
    attributeTransformedLow++;
    attributeTransformedHigh++;
  }
}

void
rahtFixedPointInverseRotation(
  FixedPoint quantStepSizeLuma,
  uint64_t weightLeft,
  uint64_t weightRight,
  FixedPoint* attributeLeft,
  FixedPoint* attributeRight,
  FixedPoint* attributeTransformedLow,
  FixedPoint* attributeTransformedHigh,
  size_t attributeCount)
{
  FixedPoint b, adjustedQuantStepSize;

  b.val = (weightRight << b.kFracBits) / (weightLeft + weightRight);

  adjustedQuantStepSize.val = sqrtFixedpoint(
    ((quantStepSizeLuma.val * quantStepSizeLuma.val)
     * (weightLeft + weightRight))
    / (weightLeft * weightRight));

  while (attributeCount--) {
    *attributeRight = *attributeTransformedHigh;
    *attributeRight *= adjustedQuantStepSize;
    *attributeLeft = *attributeRight;
    *attributeLeft *= b;
    *attributeLeft -= *attributeTransformedLow;
    *attributeRight -= *attributeLeft;

    attributeLeft->val = -attributeLeft->val;

    attributeLeft++;
    attributeRight++;
    attributeTransformedLow++;
    attributeTransformedHigh++;
  }
}

//============================================================================
/*
 * RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 *
 * Outputs:
 * weights = list of 'voxelCount' weights associated with each transform coefficient
 * integerizedAttributes = quantized transformed attributes array, in column-major order
 * binaryLayer = binary layer where each coefficient was generated
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalTransform(
  FixedPoint quantStepSizeLuma,
  long long* mortonCode,
  FixedPoint* attributes,
  uint64_t* weight,
  int* binaryLayer,
  const int attribCount,
  const int voxelCount,
  int* integerizedAttributes)
{
  size_t M, N;
  size_t S, d, i, j;

  N = voxelCount;

  FixedPoint* attributesTransformed = new FixedPoint[voxelCount * attribCount];
  uint64_t* weightTransformed = new uint64_t[voxelCount];
  long long* mortonCodeTransformed = new long long[voxelCount];

  d = 0;
  while (N > 1 || d % 2) {
    d++;
    i = 0;
    M = 0;
    S = N;

    while (i < S) {
      j = i + 1;
      mortonCodeTransformed[M] = mortonCode[i] >> 1;

      if (j < S && mortonCodeTransformed[M] == (mortonCode[j] >> 1)) {
        N--;

        weightTransformed[M] = weightTransformed[N] = weight[i] + weight[j];

        rahtFixedPointRotation(
          quantStepSizeLuma, weight[i], weight[j],
          attributes + i * attribCount, attributes + j * attribCount,
          attributesTransformed + M * attribCount,
          attributesTransformed + N * attribCount, attribCount);

        binaryLayer[N] = d;

        i += 2;
      } else {
        weightTransformed[M] = weight[i];
        for (size_t k = 0; k < attribCount; k++)
          attributesTransformed[M * attribCount + k] =
            attributes[i * attribCount + k];
        i += 1;
      }
      M++;
    }

    i = S * attribCount;
    j = N * attribCount;
    while (i-- > j)
      attributes[i] = attributesTransformed[i];

    std::swap(attributes, attributesTransformed);
    std::swap(mortonCode, mortonCodeTransformed);
    std::swap(weight, weightTransformed);
  }
  binaryLayer[0] = d;

  delete[] weightTransformed;
  delete[] mortonCodeTransformed;
  delete[] attributesTransformed;

  // Quantization of DC coefficients
  quantStepSizeLuma.val = sqrtFixedpoint(
    (quantStepSizeLuma.val * quantStepSizeLuma.val) / weight[0]);
  for (size_t k = 0; k < attribCount; k++)
    attributes[k] /= quantStepSizeLuma;

  for (i = 0; i < voxelCount; i++)
    for (size_t k = 0; k < attribCount; k++)
      integerizedAttributes[i + voxelCount * k] =
        attributes[i * attribCount + k].round();
}

//============================================================================
/*
 * inverse RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 * integerizedAttributes = quantized transformed attributes array, in column-major order
 *
 * Outputs:
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalInverseTransform(
  FixedPoint quantStepSizeLuma,
  long long* mortonCode,
  FixedPoint* attributes,
  uint64_t* weight,
  const int attribCount,
  const int voxelCount,
  int* integerizedAttributes)
{
  size_t M, N;
  size_t S, d, i, j;

  M = N = voxelCount;

  FixedPoint* attributesTransformed = new FixedPoint[voxelCount * attribCount];
  uint64_t* weightTransformed = new uint64_t[voxelCount];
  long long* mortonCodeTransformed = new long long[voxelCount];

  long long** mortonCodeBuffer;
  uint64_t** weightBuffer;
  size_t* M_buffer;
  {
    size_t depth = 0;
    long long maxMortonCode = mortonCode[voxelCount - 1];

    while (maxMortonCode--) {
      depth++;
      maxMortonCode >>= 1;
    }

    mortonCodeBuffer = new long long*[depth + 1];
    weightBuffer = new uint64_t*[depth + 1];
    M_buffer = new size_t[depth + 1];
  }

  for (i = 0; i < voxelCount; i++)
    for (size_t k = 0; k < attribCount; k++)
      attributesTransformed[i * attribCount + k] =
        (double)integerizedAttributes[i + k * voxelCount];

  // Dequantization of DC coefficients
  {
    FixedPoint quantStepSizeLuma2;
    quantStepSizeLuma2.val =
      sqrtFixedpoint((quantStepSizeLuma.val * quantStepSizeLuma.val) / N);
    for (size_t k = 0; k < attribCount; k++)
      attributesTransformed[k] *= quantStepSizeLuma2;
  }

  // Re-obtain weights at the decoder by partially executing the encoder
  d = 0;
  while (N > 1) {
    long long* _mortonCode = new long long[M];
    uint64_t* _weight = new uint64_t[M];

    M_buffer[d] = M;
    mortonCodeBuffer[d] = _mortonCode;
    weightBuffer[d] = _weight;

    for (i = 0; i < M; i++) {
      _mortonCode[i] = mortonCode[i] >> 1;
      _weight[i] = weight[i];
    }
    d++;

    i = 0;
    S = M;
    M = 0;

    while (i < S) {
      mortonCodeTransformed[M] = _mortonCode[i];
      j = i + 1;

      if (j < S && _mortonCode[i] == _mortonCode[j]) {
        N--;
        weightTransformed[M] = weight[i] + weight[j];
        i += 2;
      } else {
        weightTransformed[M] = weight[i];
        i += 1;
      }
      M++;
    }

    std::swap(mortonCode, mortonCodeTransformed);
    std::swap(weight, weightTransformed);
  }

  if (d % 2) {
    delete[] mortonCode;
    delete[] weight;
  } else {
    delete[] mortonCodeTransformed;
    delete[] weightTransformed;
  }

  // Inverse transform
  while (d--) {
    S = M_buffer[d];

    mortonCode = mortonCodeBuffer[d];
    weight = weightBuffer[d];

    M = 0;
    N = S;
    i = 0;

    while (i < S) {
      j = i + 1;

      if (j < S && mortonCode[i] == mortonCode[j]) {
        N--;

        rahtFixedPointInverseRotation(
          quantStepSizeLuma, weight[i], weight[j],
          attributes + i * attribCount, attributes + j * attribCount,
          attributesTransformed + M * attribCount,
          attributesTransformed + N * attribCount, attribCount);

        i += 2;
      } else {
        for (size_t k = 0; k < attribCount; k++)
          attributes[i * attribCount + k] =
            attributesTransformed[M * attribCount + k];
        i += 1;
      }
      M++;
    }

    delete[] mortonCode;
    delete[] weight;

    i = S * attribCount;
    while (i--)
      attributesTransformed[i] = attributes[i];
  }

  delete[] M_buffer;
  delete[] weightBuffer;
  delete[] mortonCodeBuffer;
  delete[] attributesTransformed;
}

//============================================================================

}  // namespace pcc
