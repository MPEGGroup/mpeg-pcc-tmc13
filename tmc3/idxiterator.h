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

#pragma once

#include <cstddef>
#include <iterator>

namespace pcc {

//============================================================================
// An index iterator.
// Dereferencing the iterator gives the index, *IdxIterator(i) == i.

class IdxIterator {
public:
  using value_type = size_t;
  using difference_type = ptrdiff_t;
  using pointer = value_type;
  using reference = value_type;
  using iterator_category = std::random_access_iterator_tag;

  IdxIterator() : _idx(value_type(-1)){};
  IdxIterator(size_t idx) : _idx(idx){};

  IdxIterator& operator++()
  {
    ++_idx;
    return *this;
  }
  IdxIterator& operator--()
  {
    --_idx;
    return *this;
  }

  IdxIterator operator++(int) { return IdxIterator(_idx++); }
  IdxIterator operator--(int) { return IdxIterator(_idx--); }

  IdxIterator& operator+=(difference_type n)
  {
    _idx += n;
    return *this;
  }
  IdxIterator& operator-=(difference_type n)
  {
    _idx -= n;
    return *this;
  }

  IdxIterator operator+(difference_type n) const { return {_idx + n}; }
  IdxIterator operator-(difference_type n) const { return {_idx - n}; }

  friend IdxIterator operator+(difference_type n, IdxIterator rhs)
  {
    return rhs + n;
  }

  friend IdxIterator operator-(difference_type n, IdxIterator rhs)
  {
    return rhs - n;
  }

  friend difference_type operator-(IdxIterator lhs, IdxIterator rhs)
  {
    return lhs._idx - rhs._idx;
  }

  reference operator*() { return _idx; }
  reference operator[](difference_type n) { return _idx + n; }

  friend bool operator==(IdxIterator lhs, IdxIterator rhs)
  {
    return lhs._idx == rhs._idx;
  }

  friend bool operator!=(IdxIterator lhs, IdxIterator rhs)
  {
    return !(lhs == rhs);
  }

  friend bool operator<(IdxIterator lhs, IdxIterator rhs)
  {
    return lhs._idx < rhs._idx;
  }

  friend bool operator>(IdxIterator lhs, IdxIterator rhs) { return rhs < lhs; }

  friend bool operator>=(IdxIterator lhs, IdxIterator rhs)
  {
    return !(lhs < rhs);
  }

  friend bool operator<=(IdxIterator lhs, IdxIterator rhs)
  {
    return !(lhs > rhs);
  }

private:
  value_type _idx;
};

//============================================================================

}  // namespace pcc
