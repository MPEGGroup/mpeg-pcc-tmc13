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

#pragma once

#include <stddef.h>

namespace pcc {

//============================================================================
// :: Entropy codec interface (Encoder)
//
// The base class must implement the following methods:
//  - void setBuffer(size_t size, const uint8_t *buffer);
//  - void start();
//  - size_t stop();
//  - void encode(int symbol, StaticBitModel&);
//  - void encode(int symbol, StaticMAryModel&);
//  - void encode(int symbol, AdaptiveBitModel&);
//  - void encode(int symbol, AdaptiveBitModelFast&);
//  - void encode(int symbol, AdaptiveMAryModel&);

template<class Base>
class EntropyEncoderWrapper : protected Base {
public:
  using Base::Base;
  using Base::buffer;
  using Base::encode;
  using Base::setBuffer;
  using Base::start;
  using Base::stop;

  //--------------------------------------------------------------------------
  // :: encoding / common binarisation methods

  void encodeExpGolomb(
    unsigned int symbol,
    int k,
    StaticBitModel& bModel0,
    AdaptiveBitModel& bModel1);

  void encodeIntACEGC(
    long predResidual,
    AdaptiveMAryModel& mModelValues,
    StaticBitModel& bModel0,
    AdaptiveBitModel& bModel1,
    const unsigned long M);

  void encodeUIntACEGC(
    long predResidual,
    AdaptiveMAryModel& mModelValues,
    StaticBitModel& bModel0,
    AdaptiveBitModel& bModel1,
    const unsigned long M);
};

//============================================================================
// :: Entropy codec interface (Decoder)
//
// The base class must implement the following methods:
//  - void setBuffer(size_t size, const uint8_t *buffer);
//  - void start();
//  - void stop();
//  - int decode(StaticBitModel&);
//  - int decode(StaticMAryModel&);
//  - int decode(AdaptiveBitModel&);
//  - int decode(AdaptiveBitModelFast&);
//  - int decode(AdaptiveMAryModel&);

template<class Base>
class EntropyDecoderWrapper : protected Base {
public:
  EntropyDecoderWrapper() : Base() {}

  using Base::decode;
  using Base::setBuffer;
  using Base::start;
  using Base::stop;

  //--------------------------------------------------------------------------
  // :: encoding / common binarisation methods

  unsigned int
  decodeExpGolomb(int k, StaticBitModel& bModel0, AdaptiveBitModel& bModel1);

  long decodeIntACEGC(
    AdaptiveMAryModel& mModelValues,
    StaticBitModel& bModel0,
    AdaptiveBitModel& bModel1,
    unsigned long exp_k,
    unsigned long M);

  unsigned long decodeUIntACEGC(
    AdaptiveMAryModel& mModelValues,
    StaticBitModel& bModel0,
    AdaptiveBitModel& bModel1,
    unsigned long exp_k,
    unsigned long M);
};

//============================================================================
// :: Various binarisation forms

inline unsigned long
IntToUInt(long value)
{
  return (value < 0) ? static_cast<unsigned long>(-1 - (2 * value))
                     : static_cast<unsigned long>(2 * value);
}

//----------------------------------------------------------------------------

inline long
UIntToInt(unsigned long uiValue)
{
  return (uiValue & 1) ? -(static_cast<long>((uiValue + 1) >> 1))
                       : (static_cast<long>(uiValue >> 1));
}

//----------------------------------------------------------------------------

template<class Base>
inline void
EntropyEncoderWrapper<Base>::encodeExpGolomb(
  unsigned int symbol,
  int k,
  StaticBitModel& bModel0,
  AdaptiveBitModel& bModel1)
{
  while (1) {
    if (symbol >= static_cast<unsigned int>(1 << k)) {
      encode(1, bModel1);
      symbol = symbol - (1 << k);
      k++;
    } else {
      encode(0, bModel1);  // now terminated zero of unary part
      while (k--)          // next binary part
      {
        encode(static_cast<signed short>((symbol >> k) & 1), bModel0);
      }
      break;
    }
  }
}

//----------------------------------------------------------------------------

template<class Base>
inline unsigned int
EntropyDecoderWrapper<Base>::decodeExpGolomb(
  int k, StaticBitModel& bModel0, AdaptiveBitModel& bModel1)
{
  unsigned int l;
  int symbol = 0;
  int binary_symbol = 0;
  do {
    l = decode(bModel1);
    if (l == 1) {
      symbol += (1 << k);
      k++;
    }
  } while (l != 0);
  while (k--)  //next binary part
    if (decode(bModel0) == 1) {
      binary_symbol |= (1 << k);
    }
  return static_cast<unsigned int>(symbol + binary_symbol);
}

//----------------------------------------------------------------------------

template<class Base>
inline long
EntropyDecoderWrapper<Base>::decodeIntACEGC(
  AdaptiveMAryModel& mModelValues,
  StaticBitModel& bModel0,
  AdaptiveBitModel& bModel1,
  unsigned long exp_k,
  unsigned long M)
{
  unsigned long uiValue = decode(mModelValues);
  if (uiValue == M) {
    uiValue += decodeExpGolomb(exp_k, bModel0, bModel1);
  }
  return UIntToInt(uiValue);
}

//----------------------------------------------------------------------------

template<class Base>
inline unsigned long
EntropyDecoderWrapper<Base>::decodeUIntACEGC(
  AdaptiveMAryModel& mModelValues,
  StaticBitModel& bModel0,
  AdaptiveBitModel& bModel1,
  unsigned long exp_k,
  unsigned long M)
{
  unsigned long uiValue = decode(mModelValues);
  if (uiValue == M) {
    uiValue += decodeExpGolomb(exp_k, bModel0, bModel1);
  }
  return uiValue;
}

//----------------------------------------------------------------------------

template<class Base>
inline void
EntropyEncoderWrapper<Base>::encodeIntACEGC(
  long predResidual,
  AdaptiveMAryModel& mModelValues,
  StaticBitModel& bModel0,
  AdaptiveBitModel& bModel1,
  unsigned long M)
{
  unsigned long uiValue = IntToUInt(predResidual);
  if (uiValue < M) {
    encode(uiValue, mModelValues);
  } else {
    encode(M, mModelValues);
    encodeExpGolomb(uiValue - M, 0, bModel0, bModel1);
  }
}

//----------------------------------------------------------------------------

template<class Base>
inline void
EntropyEncoderWrapper<Base>::encodeUIntACEGC(
  long predResidual,
  AdaptiveMAryModel& mModelValues,
  StaticBitModel& bModel0,
  AdaptiveBitModel& bModel1,
  unsigned long M)
{
  unsigned long uiValue = static_cast<unsigned long>(predResidual);
  if (uiValue < M) {
    encode(uiValue, mModelValues);
  } else {
    encode(M, mModelValues);
    encodeExpGolomb(uiValue - M, 0, bModel0, bModel1);
  }
}

//============================================================================

}  // namespace pcc
