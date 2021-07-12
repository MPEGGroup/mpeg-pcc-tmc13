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

namespace pcc {

//============================================================================
// A free-running frame counter, synchronised by updates

class FrameCtr {
public:
  operator int() const { return _frameCtr; }

  // Update the frame counter using the current value of frame_ctr_lsb.
  void update(int frame_ctr_lsb, int frame_ctr_lsb_bits);

  // Query whether frame_ctr_lsb does not match the current frame counter.
  bool isDifferentFrame(int frame_ctr_lsb, int frame_ctr_lsb_bits) const
  {
    return frame_ctr_lsb != (_frameCtr & ((1 << frame_ctr_lsb_bits) - 1));
  }

private:
  // The reconstructed frame counter value.
  int _frameCtr = 0;
};

//----------------------------------------------------------------------------

inline void
FrameCtr::update(int frame_ctr_lsb, int frame_ctr_lsb_bits)
{
  int window = (1 << frame_ctr_lsb_bits) >> 1;
  int curLsb = unsigned(_frameCtr) & ((1 << frame_ctr_lsb_bits) - 1);
  int curMsb = unsigned(_frameCtr) >> frame_ctr_lsb_bits;

  if (frame_ctr_lsb < curLsb && curLsb - frame_ctr_lsb >= window)
    curMsb++;
  else if (frame_ctr_lsb > curLsb && frame_ctr_lsb - curLsb > window)
    curMsb--;

  _frameCtr = (curMsb << frame_ctr_lsb_bits) + frame_ctr_lsb;
}

//============================================================================

}  // namespace pcc
