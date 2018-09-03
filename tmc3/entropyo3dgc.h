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

#include "ArithmeticCodec.h"

#include <stdint.h>
#include <stddef.h>

namespace pcc {
namespace o3dgc {

  //============================================================================

  class AdaptiveBitModelFast : public ::o3dgc::Adaptive_Bit_Model {
  public:
    AdaptiveBitModelFast() { ::o3dgc::Adaptive_Bit_Model::reset(true); }

    void reset() { ::o3dgc::Adaptive_Bit_Model::reset(true); }
  };

  //============================================================================

  class ArithmeticEncoder : public ::o3dgc::Arithmetic_Codec {
  public:
    using ::o3dgc::Arithmetic_Codec::Arithmetic_Codec;

    void setBuffer(size_t size, char* buffer)
    {
      set_buffer(size, reinterpret_cast<uint8_t*>(buffer));
    }

    void start() { start_encoder(); }

    size_t stop() { return stop_encoder(); }

    uint8_t* buffer() { return ::o3dgc::Arithmetic_Codec::buffer(); }

    void encode(int bit, ::o3dgc::Static_Bit_Model& model)
    {
      ::o3dgc::Arithmetic_Codec::encode(bit, model);
    }

    void encode(int data, ::o3dgc::Static_Data_Model& model)
    {
      ::o3dgc::Arithmetic_Codec::encode(data, model);
    }

    void encode(int bit, ::o3dgc::Adaptive_Bit_Model& model)
    {
      ::o3dgc::Arithmetic_Codec::encode(bit, model);
    }

    void encode(int bit, AdaptiveBitModelFast& model)
    {
      ::o3dgc::Arithmetic_Codec::encode(bit, model);
    }

    void encode(int data, ::o3dgc::Adaptive_Data_Model& model)
    {
      ::o3dgc::Arithmetic_Codec::encode(data, model);
    }
  };

  //============================================================================

  class ArithmeticDecoder : public ::o3dgc::Arithmetic_Codec {
  public:
    void setBuffer(size_t size, const char* buffer)
    {
      set_buffer(size, reinterpret_cast<uint8_t*>(const_cast<char*>(buffer)));
    }

    void start() { start_decoder(); }

    void stop() { stop_decoder(); }

    int decode(::o3dgc::Static_Bit_Model& model)
    {
      return ::o3dgc::Arithmetic_Codec::decode(model);
    }

    int decode(::o3dgc::Static_Data_Model& model)
    {
      return ::o3dgc::Arithmetic_Codec::decode(model);
    }

    int decode(::o3dgc::Adaptive_Bit_Model& model)
    {
      return ::o3dgc::Arithmetic_Codec::decode(model);
    }

    int decode(AdaptiveBitModelFast& model)
    {
      return ::o3dgc::Arithmetic_Codec::decode(model);
    }

    int decode(::o3dgc::Adaptive_Data_Model& model)
    {
      return ::o3dgc::Arithmetic_Codec::decode(model);
    }
  };

  //============================================================================

}  // namespace o3dgc
}  // namespace pcc

//============================================================================

namespace pcc {

using StaticBitModel = ::o3dgc::Static_Bit_Model;
using StaticMAryModel = ::o3dgc::Static_Data_Model;
using AdaptiveBitModel = ::o3dgc::Adaptive_Bit_Model;
using AdaptiveMAryModel = ::o3dgc::Adaptive_Data_Model;
using ::pcc::o3dgc::AdaptiveBitModelFast;

}  // namespace pcc
