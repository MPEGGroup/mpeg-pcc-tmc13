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

#include "dependencies/schroedinger/schroarith.h"
#include "entropychunk.h"

#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include <vector>

namespace pcc {
namespace dirac {

  //==========================================================================

  const uint16_t diraclut[256] = {
    //LUT corresponds to window = 16 @ p0=0.5 & 256 @ p=1.0
    0,    2,    5,    8,    11,   15,   20,   24,   29,   35,   41,   47,
    53,   60,   67,   74,   82,   89,   97,   106,  114,  123,  132,  141,
    150,  160,  170,  180,  190,  201,  211,  222,  233,  244,  256,  267,
    279,  291,  303,  315,  327,  340,  353,  366,  379,  392,  405,  419,
    433,  447,  461,  475,  489,  504,  518,  533,  548,  563,  578,  593,
    609,  624,  640,  656,  672,  688,  705,  721,  738,  754,  771,  788,
    805,  822,  840,  857,  875,  892,  910,  928,  946,  964,  983,  1001,
    1020, 1038, 1057, 1076, 1095, 1114, 1133, 1153, 1172, 1192, 1211, 1231,
    1251, 1271, 1291, 1311, 1332, 1352, 1373, 1393, 1414, 1435, 1456, 1477,
    1498, 1520, 1541, 1562, 1584, 1606, 1628, 1649, 1671, 1694, 1716, 1738,
    1760, 1783, 1806, 1828, 1851, 1874, 1897, 1920, 1935, 1942, 1949, 1955,
    1961, 1968, 1974, 1980, 1985, 1991, 1996, 2001, 2006, 2011, 2016, 2021,
    2025, 2029, 2033, 2037, 2040, 2044, 2047, 2050, 2053, 2056, 2058, 2061,
    2063, 2065, 2066, 2068, 2069, 2070, 2071, 2072, 2072, 2072, 2072, 2072,
    2072, 2071, 2070, 2069, 2068, 2066, 2065, 2063, 2060, 2058, 2055, 2052,
    2049, 2045, 2042, 2038, 2033, 2029, 2024, 2019, 2013, 2008, 2002, 1996,
    1989, 1982, 1975, 1968, 1960, 1952, 1943, 1934, 1925, 1916, 1906, 1896,
    1885, 1874, 1863, 1851, 1839, 1827, 1814, 1800, 1786, 1772, 1757, 1742,
    1727, 1710, 1694, 1676, 1659, 1640, 1622, 1602, 1582, 1561, 1540, 1518,
    1495, 1471, 1447, 1422, 1396, 1369, 1341, 1312, 1282, 1251, 1219, 1186,
    1151, 1114, 1077, 1037, 995,  952,  906,  857,  805,  750,  690,  625,
    553,  471,  376,  255};

  //==========================================================================

  struct SchroContext {
    uint16_t probability = 0x8000;  // p=0.5

    template<class... Args>
    void reset(Args...)
    {
      probability = 0x8000;
    }
  };

  //--------------------------------------------------------------------------
  // The approximate (7 bit) probability of a symbol being 1 or 0 according
  // to a context model.

  inline int approxSymbolProbability(int bit, SchroContext& model)
  {
    int p = std::max(1, model.probability >> 9);
    return bit ? 128 - p : p;
  }

  //==========================================================================

  struct SchroContextFixed {
    //uint16_t probability = 0x8000; // p=0.5
  };

  //==========================================================================
  // context definition that automatically binarises an m-ary symbol

  struct SchroMAryContext {
    SchroMAryContext() = default;

    SchroMAryContext(int numsyms) { set_alphabet(numsyms); }

    void set_alphabet(int numsyms);

    int numsyms;
    std::vector<uint16_t> probabilities;
  };

  //==========================================================================

  struct InvalidContext {
    void set_alphabet(int) {}
  };

  //==========================================================================

  class ArithmeticEncoder {
  public:
    ArithmeticEncoder() = default;
    ArithmeticEncoder(size_t bufferSize, std::nullptr_t)
    {
      setBuffer(bufferSize, nullptr);
    }

    //------------------------------------------------------------------------

    void setBuffer(size_t size, uint8_t* buffer)
    {
      _bufSize = size;
      if (buffer)
        _buf = _bufWr = buffer;
      else {
        allocatedBuffer.reset(new uint8_t[size]);
        _buf = _bufWr = allocatedBuffer.get();
      }
    }

    //------------------------------------------------------------------------

    void enableBypassStream(bool cabac_bypass_stream_enabled_flag)
    {
      _cabac_bypass_stream_enabled_flag = cabac_bypass_stream_enabled_flag;
    }

    //------------------------------------------------------------------------

    void setBypassBinCodingWithoutProbUpdate(bool bypass_bin_coding_without_prob_update)
    {
      _bypass_bin_coding_without_prob_update = bypass_bin_coding_without_prob_update;
    }

    //------------------------------------------------------------------------

    void start()
    {
      if (!_cabac_bypass_stream_enabled_flag)
        schro_arith_encode_init(&impl, &writeByteCallback, this);
      else {
        _chunkStream.reset(_bufWr, _bufSize);
        schro_arith_encode_init(&impl, &writeChunkCallback, &_chunkStream);
      }
    }

    //------------------------------------------------------------------------

    size_t stop()
    {
      schro_arith_flush(&impl);

      if (_cabac_bypass_stream_enabled_flag) {
        _chunkStream.flush();
        return _chunkStream.size();
      }

      return _bufWr - _buf;
    }

    //------------------------------------------------------------------------

    const char* buffer() { return reinterpret_cast<char*>(_buf); }

    //------------------------------------------------------------------------

    void encode(int bit, SchroContextFixed&) { encode(bit); }

    //------------------------------------------------------------------------

    void encode(int bit)
    {
      if (!_cabac_bypass_stream_enabled_flag) {
        if (_bypass_bin_coding_without_prob_update)
          schro_arith_encode_bypass_bit(&impl, bit);
        else {
          uint16_t probability = 0x8000;  // p=0.5
          schro_arith_encode_bit(&impl, &probability, bit);
        }
        return;
      }

      _chunkStream.writeBypassBit(bit);
    }

    //------------------------------------------------------------------------

    void encode(int data, SchroMAryContext& model);

    //------------------------------------------------------------------------

    void encode(int data, InvalidContext& model) { assert(0); }

    void encode(int bit, SchroContext& model)
    {
      schro_arith_encode_bit(&impl, &model.probability, bit);
    }

    //------------------------------------------------------------------------

    void encode(
      int bit,
      const int& offset,
      SchroContext& model,
      uint16_t* obufSingleBound)
    {
      uint16_t* probability = &model.probability;
      uint16_t& lowTh = obufSingleBound[offset + 1];
      uint16_t& HighTh = obufSingleBound[offset];
      if (*probability > HighTh) {
        *probability = HighTh;
        HighTh += (diraclut[255 - ((HighTh) >> 8)] >> 2);
        if (offset > 0 && HighTh > obufSingleBound[offset - 1]) {
          HighTh = obufSingleBound[offset - 1];
        }
      } else if (*probability < lowTh) {
        *probability = lowTh;
        lowTh -= (diraclut[(lowTh) >> 8] >> 2);
        if (offset < 31 && lowTh < obufSingleBound[offset + 2]) {
          lowTh = obufSingleBound[offset + 2];
        }
      }
      schro_arith_encode_bit(&impl, probability, bit);
      return;
    }

    //------------------------------------------------------------------------
  private:
    static void writeByteCallback(uint8_t byte, void* thisptr)
    {
      auto _this = reinterpret_cast<ArithmeticEncoder*>(thisptr);
      if (_this->_bufSize == 0)
        throw std::runtime_error("Aec stream overflow");
      _this->_bufSize--;
      *_this->_bufWr++ = byte;
    }

    //------------------------------------------------------------------------

    static void writeChunkCallback(uint8_t byte, void* thisptr)
    {
      auto _this = reinterpret_cast<ChunkStreamBuilder*>(thisptr);
      _this->writeAecByte(byte);
    }
    //------------------------------------------------------------------------

  private:
    ::SchroArith impl;
    uint8_t* _buf;
    uint8_t* _bufWr;
    size_t _bufSize;
    std::unique_ptr<uint8_t[]> allocatedBuffer;

    // Controls entropy coding method for bypass bins
    bool _cabac_bypass_stream_enabled_flag = false;

    // Controls separate coding of bypass bins
    bool _bypass_bin_coding_without_prob_update = false;

    ChunkStreamBuilder _chunkStream;
  };

  //==========================================================================

  class ArithmeticDecoder {
  public:
    void setBuffer(size_t size, const char* buffer)
    {
      _buffer = reinterpret_cast<const uint8_t*>(buffer);
      _bufferLen = size;
    }

    //------------------------------------------------------------------------

    void enableBypassStream(bool cabac_bypass_stream_enabled_flag)
    {
      _cabac_bypass_stream_enabled_flag = cabac_bypass_stream_enabled_flag;
    }

    //------------------------------------------------------------------------

    void setBypassBinCodingWithoutProbUpdate(bool bypass_bin_coding_without_prob_update)
    {
      _bypass_bin_coding_without_prob_update = bypass_bin_coding_without_prob_update;
    }

    //------------------------------------------------------------------------

    void start()
    {
      if (_cabac_bypass_stream_enabled_flag) {
        _chunkReader.reset(_buffer, _bufferLen);
        schro_arith_decode_init(&impl, &readChunkCallback, &_chunkReader);
      } else {
        schro_arith_decode_init(&impl, &readByteCallback, this);
      }
    }

    //------------------------------------------------------------------------

    void stop() { schro_arith_decode_flush(&impl); }

    //------------------------------------------------------------------------
    // Terminate the arithmetic decoder, and reinitialise to start decoding
    // the next entropy stream.

    void flushAndRestart()
    {
      stop();
      if (_cabac_bypass_stream_enabled_flag) {
        _chunkReader.nextStream();
        schro_arith_decode_init(&impl, &readChunkCallback, &_chunkReader);
      } else {
        schro_arith_decode_init(&impl, &readByteCallback, this);
      }
    }

    //------------------------------------------------------------------------

    int decode(SchroContextFixed&) { return decode(); }

    //------------------------------------------------------------------------

    int decode()
    {
      if (!_cabac_bypass_stream_enabled_flag) {
        if (_bypass_bin_coding_without_prob_update)
          return schro_arith_decode_bypass_bit(&impl);
        else {
          uint16_t probability = 0x8000;  // p=0.5
          return schro_arith_decode_bit(&impl, &probability);
        }
      }

      return _chunkReader.readBypassBit();
    }

    //------------------------------------------------------------------------

    int decode(SchroMAryContext& model);

    //------------------------------------------------------------------------

    int decode(InvalidContext& model)
    {
      assert(0);
      return 0;
    }

    //------------------------------------------------------------------------

    int decode(SchroContext& model)
    {
      return schro_arith_decode_bit(&impl, &model.probability);
    }

    //------------------------------------------------------------------------

    int
    decode(const int& offset, SchroContext& model, uint16_t* obufSingleBound)
    {
      uint16_t* probability = &model.probability;
      uint16_t& lowTh = obufSingleBound[offset + 1];
      uint16_t& HighTh = obufSingleBound[offset];
      if (*probability > HighTh) {
        *probability = HighTh;
        HighTh += (diraclut[255 - ((HighTh) >> 8)] >> 2);
        if (offset > 0 && HighTh > obufSingleBound[offset - 1]) {
          HighTh = obufSingleBound[offset - 1];
        }
      } else if (*probability < lowTh) {
        *probability = lowTh;
        lowTh -= (diraclut[(lowTh) >> 8] >> 2);
        if (offset < 31 && lowTh < obufSingleBound[offset + 2]) {
          lowTh = obufSingleBound[offset + 2];
        }
      }

      return schro_arith_decode_bit(&impl, probability);
    }

    //------------------------------------------------------------------------
  private:
    static uint8_t readByteCallback(void* thisptr)
    {
      auto _this = reinterpret_cast<ArithmeticDecoder*>(thisptr);
      if (!_this->_bufferLen)
        return 0xff;
      _this->_bufferLen--;
      return *_this->_buffer++;
    }

    //------------------------------------------------------------------------

    static uint8_t readChunkCallback(void* thisptr)
    {
      auto _this = reinterpret_cast<ChunkStreamReader*>(thisptr);
      return _this->readAecByte();
    }

    //------------------------------------------------------------------------

  private:
    ::SchroArith impl;

    // the user supplied buffer.
    const uint8_t* _buffer;

    // the length of the user supplied buffer
    size_t _bufferLen;

    // Controls entropy coding method for bypass bins
    bool _cabac_bypass_stream_enabled_flag = false;

    // Controls separate coding of bypass bins
    bool _bypass_bin_coding_without_prob_update = false;

    // Parser for chunked bypass stream representation
    ChunkStreamReader _chunkReader;
  };

  //==========================================================================

}  // namespace dirac

using StaticBitModel = dirac::SchroContextFixed;
using StaticMAryModel = dirac::InvalidContext;
using AdaptiveBitModel = dirac::SchroContext;
using AdaptiveBitModelFast = dirac::SchroContext;
using AdaptiveMAryModel = dirac::SchroMAryContext;

//============================================================================

}  // namespace pcc
