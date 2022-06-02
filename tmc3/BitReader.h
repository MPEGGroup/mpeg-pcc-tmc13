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

#include <cstdint>
#include <cstdlib>

namespace pcc {

//============================================================================

template<class ForwardIt>
class BitReader {
public:
  BitReader(ForwardIt bs, ForwardIt bs_end)
    : _bs(bs), _bsEnd(bs_end), _bitMask(), _buffer()
  {}

  // return the current bytestream read position;
  ForwardIt pos() { return _bs; }

  bool read();

  uint64_t readUn(int num_bits);

  int64_t readSn(int num_bits);

  uint64_t readUe();

  int64_t readSe();

  float readF();

  template<typename T>
  void read(T* value)
  {
    *value = T(read());
  }

  template<typename T>
  void readUn(int num_bits, T* value)
  {
    *value = T(readUn(num_bits));
  }

  template<typename T>
  void readSn(int num_bits, T* value)
  {
    *value = T(readSn(num_bits));
  }

  template<typename T>
  void readUe(T* value)
  {
    *value = T(readUe());
  }

  template<typename T>
  void readSe(T* value)
  {
    *value = T(readSe());
  }

  template<typename T>
  void readF(T* value)
  {
    *value = T(readF());
  }

  void byteAlign();

private:
  ForwardIt _bs;
  ForwardIt _bsEnd;
  int _bitMask;
  uint8_t _buffer;
};

//============================================================================

template<typename T>
BitReader<T>
makeBitReader(T start, T end)
{
  return BitReader<T>(start, end);
}

//============================================================================

template<class ForwardIt>
bool
BitReader<ForwardIt>::read()
{
  if (_bitMask == 0) {
    if (_bs == _bsEnd)
      return false;

    _buffer = *_bs;
    _bs++;
    _bitMask = 1 << 7;
  }

  bool value = _buffer & _bitMask;
  _bitMask >>= 1;
  return value;
}

//----------------------------------------------------------------------------

template<class ForwardIt>
void
BitReader<ForwardIt>::byteAlign()
{
  _bitMask = 0;
  return;
}

//----------------------------------------------------------------------------

template<class ForwardIt>
uint64_t
BitReader<ForwardIt>::readUn(int num_bits)
{
  uint64_t value = 0;
  for (int i = 0; i < num_bits; i++) {
    value <<= 1;
    value |= unsigned(read());
  }
  return value;
}

//----------------------------------------------------------------------------

template<class ForwardIt>
int64_t
BitReader<ForwardIt>::readSn(int num_bits)
{
  int64_t value = readUn(num_bits);
  bool sign = read();
  return sign ? -value : value;
}

//----------------------------------------------------------------------------

template<class ForwardIt>
uint64_t
BitReader<ForwardIt>::readUe()
{
  int len = 0;
  while (!read())
    len++;

  uint64_t value = (1ull << len) | readUn(len);
  return value - 1;
}

//----------------------------------------------------------------------------

template<class ForwardIt>
int64_t
BitReader<ForwardIt>::readSe()
{
  uint64_t value = readUe();
  bool sign = value & 1;
  value = (value + sign) >> 1;
  return sign ? value : -value;
}

//----------------------------------------------------------------------------

template<class ForwardIt>
float
BitReader<ForwardIt>::readF()
{
  uint32_t bits = uint32_t(readUn(32));
  char* data = reinterpret_cast<char*>(&bits);
  return *reinterpret_cast<float*>(data);
}

//============================================================================

}  // namespace pcc