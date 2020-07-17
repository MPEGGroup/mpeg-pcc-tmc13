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

#include <cassert>
#include <cstdint>

#include "entropy.h"

namespace pcc {

//============================================================================
// A forward and inverse symbol to frequency-sorted index lookup table.
//
// Maintains mappings for the most @_lutSize frequent symbols from an
// @_alphabetSize alphabet.
//
// Mappings are recomputed using frequency statistics maintained for
// each symbol on a periodic basis following symbol insertion according
// to an exponential back-off and maximum update period.
//

template<int _lutSize, int _alphabetSize>
class FrequencySortingLut {
public:
  // Implementation detail assumes that:
  //    _alphabetSize * _maxHistogramCount <= UINT32_MAX
  static_assert(_alphabetSize <= 256, "Symbols must be no more than 8-bit");

  // It makes little sense for the LUT to be bigger than the alphabet
  // xxx
  static_assert(_lutSize <= _alphabetSize, "LUT is too large");

  static const int kMaxUpdatePeriod = 0x33333333;
  static const int kMaxHistogramCount = 1 << 24;

  FrequencySortingLut(
    int maxUpdatePeriod = kMaxUpdatePeriod,
    int _maxHistogramCount = kMaxHistogramCount,
    const uint8_t buffer[_lutSize] = nullptr);

  FrequencySortingLut(const FrequencySortingLut&) = default;
  FrequencySortingLut(FrequencySortingLut&&) = default;
  FrequencySortingLut& operator=(const FrequencySortingLut&) = default;
  FrequencySortingLut& operator=(FrequencySortingLut&&) = default;

  void init(const uint8_t buffer[_lutSize]);

  void reset() { _reset = true; }

  // Update the histogram entry for @symbol and recompute LUT if sufficient
  // symbols have been processed.
  void pushSymbol(int symbol);

  // The sorted index for @symbol
  int getIndex(int symbol) const
  {
    assert(unsigned(symbol) < _alphabetSize);
    return _toIndex[symbol];
  }

  // The @index -th symbol
  int getSymbol(int index) const
  {
    assert(unsigned(index) < _lutSize);
    return _toSymbol[index];
  }

  int end() const { return kUndefinedIndex; }

private:
  static const int kInitialUpdatePeriod = 16;
  static const int kUndefinedIndex = -1;

  // Re-compute the LUT based on symbol histogram.
  void update();

  // per-symbol occurence counts
  int _histogram[_alphabetSize];

  // mapping of symbol to LUT index
  // NB: type must allow storage of distinct kUndefinedIndex
  int8_t _toIndex[_alphabetSize];

  // mapping of LUT index to symbol
  uint8_t _toSymbol[_lutSize];

  int _maxHistogramCount;
  unsigned _maxUpdatePeriod;

  unsigned _updatePeriod = kInitialUpdatePeriod;
  unsigned _symbolsUntilUpdate = kInitialUpdatePeriod;

  bool _reset = false;
};

//============================================================================
// A forward and inverse cache of recently used symbols.
//
// A least recently used eviction policy is used to update the cache.

template<int _cacheSize, int _alphabetSize>
class FrequentSymbolCache {
public:
  // Implementation detail assumes 8-bit alphabet:
  static_assert(_alphabetSize <= 256, "Symbols must be no more than 8-bit");

  // It makes little sense for the LUT to be bigger than the alphabet
  static_assert(_cacheSize <= _alphabetSize, "LUT is larger than alphabet?");

  FrequentSymbolCache();
  FrequentSymbolCache(const FrequentSymbolCache&) = default;
  FrequentSymbolCache(FrequentSymbolCache&&) = default;
  FrequentSymbolCache& operator=(const FrequentSymbolCache&) = default;
  FrequentSymbolCache& operator=(FrequentSymbolCache&&) = default;

  void pushSymbol(int symbol);

  int getIndex(int symbol) const
  {
    assert(unsigned(symbol) < _alphabetSize);
    return _toIndex[symbol];
  }

  int getSymbol(int index) const
  {
    assert(unsigned(index) < _cacheSize);
    return _toSymbol[index];
  }

  int end() const { return kUndefinedIndex; }

private:
  static const int kUndefinedIndex = -1;

  // mapping of symbol to cached index
  int8_t _toIndex[_alphabetSize];

  // mapping of cached index to symbol
  uint8_t _toSymbol[_cacheSize];

  unsigned _last;
};

//============================================================================

template<bool _limitedContextMode>
class DualLutCoder {
  static const int kLog2CacheSize = 4;
  static const int kCacheSize = 1 << kLog2CacheSize;
  static const int kLog2LutSize = 5;
  static const int kLutSize = 1 << kLog2LutSize;
  static const int kNumLutContexts = _limitedContextMode ? 5 : 31;

public:
  DualLutCoder();
  DualLutCoder(const uint8_t initTable[kLutSize]);

  DualLutCoder(const DualLutCoder&) = default;
  DualLutCoder(DualLutCoder&&) = default;
  DualLutCoder& operator=(const DualLutCoder&) = default;
  DualLutCoder& operator=(DualLutCoder&&) = default;

  void init(const uint8_t initTable[32]);
  void resetLut();

  void encode(int symbol, EntropyEncoder* arithmeticEncoder);
  int decode(EntropyDecoder* arithmeticDecoder);

private:
  void encodeFrequencySortedLutIndex(int index, EntropyEncoder* entropy);

  int decodeFrequencySortedLutIndex(EntropyDecoder* entropy);

  //  bool _limitedContextMode;

  FrequentSymbolCache<kCacheSize, 256> _cache;
  FrequencySortingLut<kLutSize, 256> _adaptiveLut;

  AdaptiveBitModelFast _ctxLutHit;
  AdaptiveBitModelFast _ctxCacheHit;
  AdaptiveBitModelFast _ctxSymbolBit;
  AdaptiveBitModelFast _ctxLutIndex[kNumLutContexts];
};

//============================================================================

}  // namespace pcc
