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

#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <memory>
#include <vector>

namespace pcc {
namespace dirac {

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
      aec_stream_len = 0;
      buf_size = int(size);

      if (buffer)
        buf = buffer;
      else {
        allocatedBuffer.reset(new uint8_t[size]);
        buf = allocatedBuffer.get();
      }

      // bypass data is written backwards, starting at the end of the buffer
      bypassPtr = buf + size;
      bypassCount = 8;
    }

    //------------------------------------------------------------------------

    void enableBypassStream(bool cabac_bypass_stream_enabled_flag)
    {
      _cabac_bypass_stream_enabled_flag = cabac_bypass_stream_enabled_flag;
    }

    //------------------------------------------------------------------------

    void start() { schro_arith_encode_init(&impl, &writeByteCallback, this); }

    //------------------------------------------------------------------------

    size_t stop()
    {
      schro_arith_flush(&impl);
      if (bypassCount != 8) {
        bypassAccum <<= bypassCount;
        *--bypassPtr = bypassAccum;
      }

      // make the bypass data contiguous with the arithmetic coded data
      // since they share an initial buffer
      auto end = std::move(bypassPtr, buf + buf_size, buf + aec_stream_len);

      return end - buf;
    }

    //------------------------------------------------------------------------

    const char* buffer() { return reinterpret_cast<char*>(buf); }

    //------------------------------------------------------------------------

    void encode(int bit, SchroContextFixed&)
    {
      if (!_cabac_bypass_stream_enabled_flag) {
        uint16_t probability = 0x8000;  // p=0.5
        schro_arith_encode_bit(&impl, &probability, bit);
        return;
      }

      bypassAccum <<= 1;
      bypassAccum |= bit;

      if (--bypassCount)
        return;

      bypassCount = 8;
      *--bypassPtr = bypassAccum;
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
  private:
    static void writeByteCallback(uint8_t byte, void* thisptr)
    {
      return reinterpret_cast<ArithmeticEncoder*>(thisptr)->writeByte(byte);
    }

    void writeByte(uint8_t byte) { buf[aec_stream_len++] = byte; }

    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

  private:
    ::SchroArith impl;
    uint8_t* buf;
    int buf_size;
    std::unique_ptr<uint8_t[]> allocatedBuffer;

    // Controls entropy coding method for bypass bins
    bool _cabac_bypass_stream_enabled_flag = false;

    int aec_stream_len;

    // State related to bypass stream coding.
    // The bypass stream is stored as a byte-reversed sequence starting at
    // the end of buf and growing downwards.

    // Pointer to the tail of the bypass stream
    uint8_t* bypassPtr;

    // Number of bins in the bypass accumulator
    int bypassCount;

    // Accumulator for bypass bins to construct
    uint8_t bypassAccum;
  };

  //==========================================================================

  class ArithmeticDecoder {
  public:
    void setBuffer(size_t size, const char* buffer)
    {
      buf = reinterpret_cast<const uint8_t*>(buffer);
      buf_end = buf + size;

      bypassPtr = buf_end - 1;
      bypassCount = 0;
    }

    //------------------------------------------------------------------------

    void enableBypassStream(bool cabac_bypass_stream_enabled_flag)
    {
      _cabac_bypass_stream_enabled_flag = cabac_bypass_stream_enabled_flag;
    }

    //------------------------------------------------------------------------

    void start() { schro_arith_decode_init(&impl, &readByteCallback, this); }

    //------------------------------------------------------------------------

    void stop() { schro_arith_decode_flush(&impl); }

    //------------------------------------------------------------------------

    int decode(SchroContextFixed&)
    {
      if (!_cabac_bypass_stream_enabled_flag) {
        uint16_t probability = 0x8000;  // p=0.5
        return schro_arith_decode_bit(&impl, &probability);
      }

      if (!bypassCount--) {
        bypassAccum = *bypassPtr--;
        bypassCount = 7;
      }
      int bit = !!(bypassAccum & 0x80);
      bypassAccum <<= 1;
      return bit;
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
  private:
    static uint8_t readByteCallback(void* thisptr)
    {
      return reinterpret_cast<ArithmeticDecoder*>(thisptr)->readByte();
    }

    uint8_t readByte()
    {
      if (buf < buf_end)
        return *buf++;
      else
        return 0xff;
    }

    //------------------------------------------------------------------------

  private:
    ::SchroArith impl;
    const uint8_t* buf;
    const uint8_t* buf_end;

    // Controls entropy coding method for bypass bins
    bool _cabac_bypass_stream_enabled_flag = false;

    // State related to bypass stream coding.
    // The bypass stream is stored as a byte-reversed sequence starting at
    // the end of buf and growing downwards.

    // Pointer to the tail of the bypass stream
    const uint8_t* bypassPtr;

    // Number of bins in the bypass accumulator
    int bypassCount;

    // Accumulator for bypass bins to construct
    uint8_t bypassAccum;
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
