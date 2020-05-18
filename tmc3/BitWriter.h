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

#include "PCCMisc.h"

#include <cstdint>
#include <cstdlib>

namespace pcc {

//============================================================================

template<class OutputIt>
class BitWriter {
public:
  BitWriter(OutputIt bs) : _bs(bs), _num_bits(), _buffer() {}

  void write(bool bit);

  // Write value using a fixed-width (num_bits) literal encoding (big endian).
  void writeUn(int num_bits, uint64_t value);

  template<typename T>
  void writeUn(int num_bits, const T& value)
  {
    writeUn(num_bits, uint64_t(value));
  }

  // Write value using a fixed-width (num_bits+sign) literal encoding (big
  // endian).
  void writeSn(int num_bits, int64_t value);

  void writeUe(uint32_t value);

  // NB: probably best to think twice before using this method
  void writeUe64(uint64_t value);

  template<typename T>
  void writeUe(const T& value)
  {
    // Avoid accidental truncation
    constexpr bool use_writeUe64 =
      std::is_same<uint64_t, T>::value || std::is_same<int64_t, T>::value;
    static_assert(!use_writeUe64, "use explicit writeUe64()");

    writeUe(uint32_t(value));
  }

  void writeSe(int32_t value);

  void writeF(float value);

  void byteAlign();

private:
  OutputIt _bs;
  int _num_bits;
  uint8_t _buffer;
};

//============================================================================

class InsertionCounter {
public:
  InsertionCounter(int* counter) : _counter(counter) {}
  InsertionCounter& operator++(int) { return *this; }
  InsertionCounter& operator++() { return *this; }
  InsertionCounter& operator*() { return *this; }

  template<typename T>
  InsertionCounter& operator=(const T&)
  {
    (*_counter)++;
    return *this;
  }

private:
  int* _counter;
};

//============================================================================

template<typename T>
BitWriter<T>
makeBitWriter(T t)
{
  return BitWriter<T>(t);
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::write(bool bit)
{
  _buffer <<= 1;
  _buffer |= uint8_t(bit);
  _num_bits++;
  if (_num_bits == 8) {
    *_bs++ = static_cast<char>(_buffer);
    _buffer = 0;
    _num_bits = 0;
  }
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::byteAlign()
{
  if (!_num_bits)
    return;

  _buffer <<= 8 - _num_bits;
  *_bs++ = static_cast<char>(_buffer);
  _buffer = 0;
  _num_bits = 0;
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::writeUn(int num_bits, uint64_t value)
{
  if (!num_bits)
    return;
  for (uint64_t mask = uint64_t(1) << (num_bits - 1); mask; mask >>= 1) {
    write(!!(value & mask));
  }
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::writeSn(int num_bits, int64_t value)
{
  writeUn(num_bits, uint64_t(::llabs(value)));
  write(value < 0);
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::writeUe(uint32_t value)
{
  value++;
  int len = ilog2(value);
  writeUn(len, 0);
  writeUn(len + 1, value);
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::writeUe64(uint64_t value)
{
  value++;
  int len = ilog2(value);
  writeUn(len, 0);
  writeUn(len + 1, value);
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::writeSe(int32_t value)
{
  bool sign = value > 0;
  value = uint32_t(::llabs(value)) << 1;
  writeUe(value - sign);
}

//----------------------------------------------------------------------------

template<class OutputIt>
void
BitWriter<OutputIt>::writeF(float value)
{
  char* data = reinterpret_cast<char*>(&value);
  uint32_t val = *reinterpret_cast<uint32_t*>(data);
  writeUn(32, val);
}

//============================================================================

}  // namespace pcc