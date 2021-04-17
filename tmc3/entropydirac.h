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
        uint16_t probability = 0x8000;  // p=0.5
        schro_arith_encode_bit(&impl, &probability, bit);
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
        uint16_t probability = 0x8000;  // p=0.5
        return schro_arith_decode_bit(&impl, &probability);
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
