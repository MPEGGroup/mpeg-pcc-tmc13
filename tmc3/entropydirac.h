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
      buf.length = int(size);

      if (buffer)
        buf.data = buffer;
      else {
        allocatedBuffer.reset(new uint8_t[size]);
        buf.data = allocatedBuffer.get();
      }
    }

    //------------------------------------------------------------------------

    void start() { schro_arith_encode_init(&impl, &buf); }

    //------------------------------------------------------------------------

    size_t stop()
    {
      schro_arith_flush(&impl);
      return impl.offset;
    }

    //------------------------------------------------------------------------

    const char* buffer() { return reinterpret_cast<char*>(buf.data); }

    //------------------------------------------------------------------------

    void encode(int bit, SchroContextFixed& model)
    {
      uint16_t probability = 0x8000;  // p=0.5
      schro_arith_encode_bit(&impl, &probability, bit);
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
    ::SchroArith impl;
    ::SchroBuffer buf;
    std::unique_ptr<uint8_t[]> allocatedBuffer;
  };

  //==========================================================================

  class ArithmeticDecoder {
  public:
    void setBuffer(size_t size, const char* buffer)
    {
      buf.data = reinterpret_cast<uint8_t*>(const_cast<char*>(buffer));
      buf.length = int(size);
    }

    //------------------------------------------------------------------------

    void start() { schro_arith_decode_init(&impl, &buf); }

    //------------------------------------------------------------------------

    void stop() { schro_arith_decode_flush(&impl); }

    //------------------------------------------------------------------------

    int decode(SchroContextFixed& model)
    {
      uint16_t probability = 0x8000;  // p=0.5
      return schro_arith_decode_bit(&impl, &probability);
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
    ::SchroArith impl;
    ::SchroBuffer buf;
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
