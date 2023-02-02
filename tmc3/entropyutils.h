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

#include <algorithm>
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
  using Base::enableBypassStream;
  using Base::encode;
  using Base::setBuffer;
  using Base::setBypassBinCodingWithoutProbUpdate;
  using Base::start;
  using Base::stop;

  //--------------------------------------------------------------------------
  // :: encoding / common binarisation methods

  void encodeExpGolomb(unsigned int symbol, int k, AdaptiveBitModel& bModel1);

  template<size_t NumPrefixCtx, size_t NumSuffixCtx>
  void encodeExpGolomb(
    unsigned int symbol,
    int k,
    AdaptiveBitModel (&ctxPrefix)[NumPrefixCtx],
    AdaptiveBitModel (&ctxSuffix)[NumSuffixCtx]);
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
  using Base::enableBypassStream;
  using Base::flushAndRestart;
  using Base::setBuffer;
  using Base::setBypassBinCodingWithoutProbUpdate;
  using Base::start;
  using Base::stop;

  //--------------------------------------------------------------------------
  // :: encoding / common binarisation methods

  unsigned int decodeExpGolomb(int k, AdaptiveBitModel& bModel1);

  template<size_t NumPrefixCtx, size_t NumSuffixCtx>
  unsigned int decodeExpGolomb(
    int k,
    AdaptiveBitModel (&ctxPrefix)[NumPrefixCtx],
    AdaptiveBitModel (&ctxSuffix)[NumSuffixCtx]);
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
  unsigned int symbol, int k, AdaptiveBitModel& ctxPrefix)
{
  while (1) {
    if (symbol >= (1u << k)) {
      encode(1, ctxPrefix);
      symbol -= (1u << k);
      k++;
    } else {
      encode(0, ctxPrefix);
      while (k--)
        encode((symbol >> k) & 1);
      break;
    }
  }
}

//----------------------------------------------------------------------------

template<class Base>
template<size_t NumPrefixCtx, size_t NumSuffixCtx>
inline void
EntropyEncoderWrapper<Base>::encodeExpGolomb(
  unsigned int symbol,
  int k,
  AdaptiveBitModel (&ctxPrefix)[NumPrefixCtx],
  AdaptiveBitModel (&ctxSuffix)[NumSuffixCtx])
{
  constexpr int maxPrefixIdx = NumPrefixCtx - 1;
  constexpr int maxSuffixIdx = NumSuffixCtx - 1;
  const int k0 = k;

  while (symbol >= (1u << k)) {
    encode(1, ctxPrefix[std::min(maxPrefixIdx, k - k0)]);
    symbol -= 1u << k;
    k++;
  }
  encode(0, ctxPrefix[std::min(maxPrefixIdx, k - k0)]);

  while (k--)
    encode((symbol >> k) & 1, ctxSuffix[std::min(maxSuffixIdx, k)]);
}

//----------------------------------------------------------------------------

template<class Base>
inline unsigned int
EntropyDecoderWrapper<Base>::decodeExpGolomb(
  int k, AdaptiveBitModel& ctxPrefix)
{
  unsigned int l;
  int symbol = 0;
  int binary_symbol = 0;
  do {
    l = decode(ctxPrefix);
    if (l == 1) {
      symbol += (1 << k);
      k++;
    }
  } while (l != 0);
  while (k--)  //next binary part
    if (decode() == 1) {
      binary_symbol |= (1 << k);
    }
  return static_cast<unsigned int>(symbol + binary_symbol);
}

//----------------------------------------------------------------------------

template<class Base>
template<size_t NumPrefixCtx, size_t NumSuffixCtx>
inline unsigned int
EntropyDecoderWrapper<Base>::decodeExpGolomb(
  int k,
  AdaptiveBitModel (&ctxPrefix)[NumPrefixCtx],
  AdaptiveBitModel (&ctxSuffix)[NumSuffixCtx])
{
  constexpr int maxPrefixIdx = NumPrefixCtx - 1;
  constexpr int maxSuffixIdx = NumSuffixCtx - 1;
  const int k0 = k;
  unsigned int l;
  int symbol = 0;
  int binary_symbol = 0;

  do {
    l = decode(ctxPrefix[std::min(maxPrefixIdx, k - k0)]);
    if (l == 1) {
      symbol += (1 << k);
      k++;
    }
  } while (l != 0);

  while (k--)
    binary_symbol |= decode(ctxSuffix[std::min(maxSuffixIdx, k)]) << k;

  return static_cast<unsigned int>(symbol + binary_symbol);
}

//============================================================================

}  // namespace pcc
