/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2019, ISO/IEC
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

#pragma once

#include <array>
#include <vector>

#include "constants.h"
#include "PCCMath.h"

namespace pcc {

struct AttributeDescription;
struct AttributeParameterSet;
struct AttributeBrickHeader;

//============================================================================
// Quantisation methods

class Quantizer {
public:
  // Derives step sizes from qp
  Quantizer(int qp);
  Quantizer(const Quantizer&) = default;
  Quantizer& operator=(const Quantizer&) = default;

  // The quantizer's step size
  int stepSize() const { return _stepSize; }

  // Quantise a value
  int64_t quantize(int64_t x) const;

  // Scale (inverse quantise) a quantised value
  int64_t scale(int64_t x) const;

private:
  // Quantisation step size
  int _stepSize;

  // Reciprocal stepsize for forward quantisation optimisation
  int _stepSizeRecip;
};

//---------------------------------------------------------------------------

inline int64_t
Quantizer::quantize(int64_t x) const
{
  // Forward quantisation avoids division by using the multiplicative inverse
  // with 18 fractional bits.
  int64_t fracBits = 18 + kFixedPointAttributeShift;

  // NB, the folowing offsets quantizes with a different deadzone to the
  // reconstruction function.
  int64_t offset = (1ll << fracBits) / 3;

  if (x >= 0) {
    return (x * _stepSizeRecip + offset) >> fracBits;
  }
  return -((offset - x * _stepSizeRecip) >> fracBits);
}

//---------------------------------------------------------------------------

inline int64_t
Quantizer::scale(int64_t x) const
{
  return x * _stepSize;
}

//============================================================================
// Encapslation of multi-component attribute quantizer values.

typedef std::array<Quantizer, 2> Quantizers;
typedef std::array<int, 2> Qps;
typedef std::vector<Qps> QpLayers;

//============================================================================

struct QpRegionOffset {
  Qps qpOffset;
  Box3<int32_t> region;
};

typedef std::vector<QpRegionOffset> QpRegionList;

//============================================================================

struct QpSet {
  QpLayers layers;
  QpRegionList regions;
  int maxQp;

  // Derive the quantizers at a given layer after applying qpOffset
  Quantizers quantizers(int qpLayer, Qps qpOffset) const;

  // Derive the quantizer for a point at a particular layer
  Quantizers quantizers(const Vec3<int32_t>& point, int qpLayer) const;

  Qps regionQpOffset(const Vec3<int32_t>& point) const;
};

//============================================================================
// Determine the Qps for a particular layer in an attribute slice
Qps deriveQps(
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  int qpLayer);

// Determine the base layer QPs for an attribute slice
QpLayers deriveLayerQps(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh);

// Determine a list of Qp offsets per region
QpRegionList deriveQpRegions(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh);

// Determine the Qp configuration for an attribute slice
QpSet deriveQpSet(
  const AttributeDescription& attrDesc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh);

//============================================================================
// Quantisation methods for geometry

class QuantizerGeom {
public:
  // Derives step sizes from qp
  QuantizerGeom(int qp)
  {
    _stepSize = (8 + (qp % 8)) << qpShift(qp);
    _stepSizeRecip = kQpStepRecip[qp % 8] >> qpShift(qp);
  }

  QuantizerGeom(const QuantizerGeom&) = default;
  QuantizerGeom& operator=(const QuantizerGeom&) = default;

  // The quantizer's step size
  int stepSize() const { return _stepSize; }

  // Quantise a value
  int64_t quantize(int64_t x) const;

  // Scale (inverse quantise) a quantised value
  int64_t scale(int64_t x) const;

  // The number of bits eliminated by a given qp
  static int qpShift(int qp) { return qp >> 3; }

private:
  // Quantisation step size
  int _stepSize;

  // Reciprocal stepsize for forward quantisation optimisation
  int _stepSizeRecip;

  static const int _shift = 20;

  static const int32_t kQpStep[8];
  static const int32_t kQpStepRecip[8];
};

//---------------------------------------------------------------------------

inline int64_t
QuantizerGeom::quantize(int64_t x) const
{
  // NB: geometry positions are only ever positive
  return (x * _stepSizeRecip + (1 << 19)) >> _shift;
}

//---------------------------------------------------------------------------

inline int64_t
QuantizerGeom::scale(int64_t x) const
{
  return (x * _stepSize + 4) >> 3;
}

//============================================================================

}  // namespace pcc
