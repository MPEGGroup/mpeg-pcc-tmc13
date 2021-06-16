/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2021, ISO/IEC
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

#include "attribute_raw.h"

#include "BitWriter.h"
#include "io_hls.h"

namespace pcc {

//============================================================================

void
AttrRawEncoder::encode(
  const SequenceParameterSet& sps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  AttributeBrickHeader& abh,
  PCCPointSet3& cloud,
  PayloadBuffer* payload)
{
  write(sps, attr_aps, abh, payload);
  auto bs = makeBitWriter(std::back_inserter(*payload));

  // todo(df): add support for variable length coding
  assert(!attr_aps.raw_attr_variable_len_flag);
  int valueBits = desc.bitdepth;

  // todo(df): update to correctly map attribute types
  if (desc.attr_num_dimensions_minus1 == 0) {
    for (size_t i = 0; i < cloud.getPointCount(); i++) {
      bs.writeUn(valueBits, cloud.getReflectance(i));
    }
  } else if (desc.attr_num_dimensions_minus1 == 2) {
    for (size_t i = 0; i < cloud.getPointCount(); i++) {
      auto& value = cloud.getColor(i);
      for (int c = 0; c < 3; c++)
        bs.writeUn(valueBits, value[c]);
    }
  } else {
    assert(
      desc.attr_num_dimensions_minus1 == 0
      || desc.attr_num_dimensions_minus1 == 2);
  }

  bs.byteAlign();
}

//============================================================================

} /* namespace pcc */
