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

#include "io_tlv.h"

#include <cstdint>

namespace pcc {

//============================================================================

std::ostream&
writeTlv(const PayloadBuffer& buf, std::ostream& os)
{
  uint32_t length = uint32_t(buf.size());

  os.put(char(buf.type));
  os.put(char(length >> 24));
  os.put(char(length >> 16));
  os.put(char(length >> 8));
  os.put(char(length >> 0));

  os.write(buf.data(), length);
  return os;
}

//============================================================================

std::istream&
readTlv(std::istream& is, PayloadBuffer* buf)
{
  buf->resize(0);
  buf->type = PayloadType(static_cast<unsigned>(is.get()));

  uint32_t length = 0;
  length = (length << 8) | static_cast<unsigned>(is.get());
  length = (length << 8) | static_cast<unsigned>(is.get());
  length = (length << 8) | static_cast<unsigned>(is.get());
  length = (length << 8) | static_cast<unsigned>(is.get());

  if (!is)
    return is;

  buf->resize(length);
  is.read(buf->data(), length);
  return is;
}

//============================================================================

}  // namespace pcc