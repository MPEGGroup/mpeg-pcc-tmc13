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

#include "DualLutCoder.h"

#include <algorithm>

namespace pcc {

//============================================================================
// :: FrequencySortingLut methods

template<int _lutSize, int _alphabetSize>
FrequencySortingLut<_lutSize, _alphabetSize>::FrequencySortingLut(
  int maxUpdatePeriod,
  int maxHistogramCount,
  const uint8_t initTable[_lutSize])
  : _maxHistogramCount(maxHistogramCount), _maxUpdatePeriod(maxUpdatePeriod)
{
  assert(maxHistogramCount <= kMaxHistogramCount);
  assert(maxUpdatePeriod <= kMaxUpdatePeriod);
  init(initTable);
}

//----------------------------------------------------------------------------

template<int _lutSize, int _alphabetSize>
void
FrequencySortingLut<_lutSize, _alphabetSize>::init(
  const uint8_t initTable[_lutSize])
{
  for (int k = 0; k < _alphabetSize; ++k) {
    _histogram[k] = 0;
    _toIndex[k] = kUndefinedIndex;
  }

  // Initialise LUT with a default mapping
  for (int k = 0; k < _lutSize; ++k) {
    int symbol = initTable ? initTable[k] : k;
    _toSymbol[k] = symbol;
    _toIndex[symbol] = k;
    _histogram[symbol] = 1;
  }
}

//----------------------------------------------------------------------------

template<int _lutSize, int _alphabetSize>
void
FrequencySortingLut<_lutSize, _alphabetSize>::update()
{
  // NB: This expression limits the value of _maxUpdatePeriod to 0x33333333
  _updatePeriod = std::min((5u * _updatePeriod) >> 2, _maxUpdatePeriod);
  _symbolsUntilUpdate = _updatePeriod;

  // Sort the symbols by occurrence.
  // NB: stability is guaranteed by including symbol value to break ties.
  uint32_t tmp[_alphabetSize];
  for (int symbol = 0; symbol < _alphabetSize; ++symbol)
    tmp[symbol] = ((~_histogram[symbol]) << 8) + symbol;

  std::nth_element(tmp, tmp + _lutSize, tmp + _alphabetSize);
  std::sort(tmp, tmp + _lutSize);

  // Remove any existing mappings
  for (int k = 0; k < _lutSize; ++k)
    _toIndex[_toSymbol[k]] = kUndefinedIndex;

  // Re-populate the LUTs
  for (int k = 0; k < _lutSize; ++k) {
    uint32_t symbol = (tmp[k] & 255);
    _toSymbol[k] = symbol;
    _toIndex[symbol] = k;
  }

  // Reset the symbol counts
  if (_reset) {
    _reset = false;

    for (int k = 0; k < _alphabetSize; ++k)
      _histogram[k] = 0;

    for (int k = 0; k < _lutSize; ++k)
      _histogram[_toSymbol[k]] = 1;
  }
}

//----------------------------------------------------------------------------

template<int _lutSize, int _alphabetSize>
void
FrequencySortingLut<_lutSize, _alphabetSize>::pushSymbol(int symbol)
{
  assert(unsigned(symbol) < _alphabetSize);

  if (++_histogram[symbol] == _maxHistogramCount) {
    for (int k = 0; k < _alphabetSize; ++k)
      _histogram[k] = _histogram[k] >> 1;
  }

  if (!(--_symbolsUntilUpdate))
    update();
}

//============================================================================
// :: LRU Cache methods

template<int _cacheSize, int _alphabetSize>
FrequentSymbolCache<_cacheSize, _alphabetSize>::FrequentSymbolCache()
{
  for (int k = 0; k < _cacheSize; ++k) {
    _toSymbol[k] = k;
    _toIndex[k] = k;
  }

  for (int k = _cacheSize; k < _alphabetSize; ++k)
    _toIndex[k] = kUndefinedIndex;

  _last = 0;
}

//----------------------------------------------------------------------------

template<int _cacheSize, int _alphabetSize>
void
FrequentSymbolCache<_cacheSize, _alphabetSize>::pushSymbol(int symbol)
{
  assert(symbol < _alphabetSize);
  const int index = _toIndex[symbol];
  const int index0 = (_last++) % _cacheSize;
  const int symbol0 = _toSymbol[index0];
  assert(symbol0 < _alphabetSize);

  std::swap(_toIndex[symbol], _toIndex[symbol0]);
  if (index == kUndefinedIndex) {
    _toSymbol[index0] = symbol;
  } else {
    std::swap(_toSymbol[index0], _toSymbol[index]);
  }
}

//============================================================================
// :: Entropy coding

template<bool _limitedContextMode>
DualLutCoder<_limitedContextMode>::DualLutCoder() : DualLutCoder(nullptr)
{}

template<>
DualLutCoder<true>::DualLutCoder(const uint8_t initTable[32])
  : _adaptiveLut(1024, 1024, initTable)
{}

template<>
DualLutCoder<false>::DualLutCoder(const uint8_t initTable[32])
  : _adaptiveLut(0x33333333, 1 << 24, initTable)
{}

//----------------------------------------------------------------------------

template<bool _limitedContextMode>
void
DualLutCoder<_limitedContextMode>::init(const uint8_t initTable[32])
{
  _adaptiveLut.init(initTable);
}

//----------------------------------------------------------------------------

template<bool _limitedContextMode>
void
DualLutCoder<_limitedContextMode>::resetLut()
{
  _adaptiveLut.reset();
}

//----------------------------------------------------------------------------

template<>
void
DualLutCoder<true>::encodeFrequencySortedLutIndex(
  int index, EntropyEncoder* entropy)
{
  bool b4 = index & 1;
  bool b3 = (index >> 1) & 1;
  bool b2 = (index >> 2) & 1;
  bool b1 = (index >> 3) & 1;
  bool b0 = (index >> 4) & 1;

  entropy->encode(b0, _ctxLutIndex[0]);
  if (b0) {
    entropy->encode(b1);
    entropy->encode(b2);
    entropy->encode(b3);
    entropy->encode(b4);
    return;
  }

  entropy->encode(b1, _ctxLutIndex[1]);
  if (b1) {
    entropy->encode(b2);
    entropy->encode(b3);
    entropy->encode(b4);
    return;
  }

  entropy->encode(b2, _ctxLutIndex[2]);
  if (b2) {
    entropy->encode(b3);
    entropy->encode(b4);
    return;
  }

  entropy->encode(b3, _ctxLutIndex[3]);
  entropy->encode(b4, _ctxLutIndex[4]);
}

//----------------------------------------------------------------------------

template<>
void
DualLutCoder<false>::encodeFrequencySortedLutIndex(
  int index, EntropyEncoder* entropy)
{
  entropy->encode((index >> 4) & 1, _ctxLutIndex[0]);
  entropy->encode((index >> 3) & 1, _ctxLutIndex[1 + (index >> 4)]);
  entropy->encode((index >> 2) & 1, _ctxLutIndex[3 + (index >> 3)]);
  entropy->encode((index >> 1) & 1, _ctxLutIndex[7 + (index >> 2)]);
  entropy->encode((index >> 0) & 1, _ctxLutIndex[15 + (index >> 1)]);
}

//----------------------------------------------------------------------------

template<bool _limitedContextMode>
void
DualLutCoder<_limitedContextMode>::encode(int value, EntropyEncoder* entropy)
{
  // One of three coding methods are used:
  //  - Encode position in LUT (if present)
  //  - Encode position in Cache (if present)
  //  - Encode symbol directly

  // LUT index coding
  int index = _adaptiveLut.getIndex(value);
  bool inLut = index != _adaptiveLut.end();
  _adaptiveLut.pushSymbol(value);

  entropy->encode(inLut, _ctxLutHit);
  if (inLut) {
    encodeFrequencySortedLutIndex(index, entropy);
    return;
  }

  // Cache index coding
  index = _cache.getIndex(value);
  bool inCache = index != _cache.end();
  _cache.pushSymbol(value);

  entropy->encode(inCache, _ctxCacheHit);
  if (inCache) {
    for (int i = 0; i < kLog2CacheSize; ++i) {
      entropy->encode(index & 1);
      index >>= 1;
    }
    return;
  }

  // Direct coding
  for (int i = 0; i < 8; ++i) {
    entropy->encode(value & 1, _ctxSymbolBit);
    value >>= 1;
  }
}

//----------------------------------------------------------------------------

template<>
int
DualLutCoder<true>::decodeFrequencySortedLutIndex(EntropyDecoder* entropy)
{
  bool b0, b1, b2, b3, b4;

  b0 = entropy->decode(_ctxLutIndex[0]);
  if (b0) {
    b1 = entropy->decode();
    b2 = entropy->decode();
    b3 = entropy->decode();
    b4 = entropy->decode();
  } else {
    b1 = entropy->decode(_ctxLutIndex[1]);
    if (b1) {
      b2 = entropy->decode();
      b3 = entropy->decode();
      b4 = entropy->decode();
    } else {
      b2 = entropy->decode(_ctxLutIndex[2]);
      if (b2) {
        b3 = entropy->decode();
        b4 = entropy->decode();
      } else {
        b3 = entropy->decode(_ctxLutIndex[3]);
        b4 = entropy->decode(_ctxLutIndex[4]);
      }
    }
  }

  return (b0 << 4) | (b1 << 3) | (b2 << 2) | (b3 << 1) | b4;
}

//----------------------------------------------------------------------------

template<>
int
DualLutCoder<false>::decodeFrequencySortedLutIndex(EntropyDecoder* entropy)
{
  int index = 0;
  index = (index << 1) | entropy->decode(_ctxLutIndex[0]);
  index = (index << 1) | entropy->decode(_ctxLutIndex[1 + index]);
  index = (index << 1) | entropy->decode(_ctxLutIndex[3 + index]);
  index = (index << 1) | entropy->decode(_ctxLutIndex[7 + index]);
  index = (index << 1) | entropy->decode(_ctxLutIndex[15 + index]);

  return index;
}

//----------------------------------------------------------------------------

template<bool _limitedContextMode>
int
DualLutCoder<_limitedContextMode>::decode(EntropyDecoder* entropy)
{
  int symbol;

  bool inLut = entropy->decode(_ctxLutHit);
  if (inLut) {
    int index = decodeFrequencySortedLutIndex(entropy);
    symbol = _adaptiveLut.getSymbol(index);
  }

  if (!inLut) {
    bool inCache = entropy->decode(_ctxCacheHit);
    if (inCache) {
      int index = 0;
      for (int i = 0; i < kLog2CacheSize; ++i) {
        index |= entropy->decode() << i;
      }
      symbol = _cache.getSymbol(index);
    } else {
      symbol = 0;
      for (int i = 0; i < 8; ++i) {
        symbol |= entropy->decode(_ctxSymbolBit) << i;
      }
    }

    _cache.pushSymbol(symbol);
  }

  _adaptiveLut.pushSymbol(symbol);

  return symbol;
}

//============================================================================

template class DualLutCoder<true>;
template class DualLutCoder<false>;

//============================================================================

}  // namespace pcc
