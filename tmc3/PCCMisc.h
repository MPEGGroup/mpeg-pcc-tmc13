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

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <utility>
#include <string>
#include <vector>

#if _MSC_VER
#  define DEPRECATED_MSVC __declspec(deprecated)
#  define DEPRECATED
#else
#  define DEPRECATED_MSVC
#  define DEPRECATED __attribute__((deprecated))
#endif

#if _MSC_VER && !defined(__attribute__)
#  define __attribute__(...)
#endif

namespace pcc {
const uint32_t PCC_UNDEFINED_INDEX = -1;

enum PCCEndianness
{
  PCC_BIG_ENDIAN = 0,
  PCC_LITTLE_ENDIAN = 1
};

inline PCCEndianness
PCCSystemEndianness()
{
  uint32_t num = 1;
  return (*(reinterpret_cast<char*>(&num)) == 1) ? PCC_LITTLE_ENDIAN
                                                 : PCC_BIG_ENDIAN;
}

//---------------------------------------------------------------------------
// Replace any occurence of %d with formatted number.  The %d format
// specifier may use the formatting conventions of snprintf().
std::string expandNum(const std::string& src, int num);

//---------------------------------------------------------------------------
// Population count -- return the number of bits set in @x.
//
inline int
popcnt(uint32_t x)
{
  x = x - ((x >> 1) & 0x55555555u);
  x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
  return ((x + (x >> 4) & 0xF0F0F0Fu) * 0x1010101u) >> 24;
}

//---------------------------------------------------------------------------
// Population count -- return the number of bits set in @x.
//
inline int
popcnt(uint8_t x)
{
  uint32_t val = x * 0x08040201u;
  val >>= 3;
  val &= 0x11111111u;
  val *= 0x11111111u;
  return val >> 28;
}

//---------------------------------------------------------------------------
// Test if population count is greater than 1.
// Returns non-zero if true.
//
inline uint32_t
popcntGt1(uint32_t x)
{
  return x & (x - 1);
}

//---------------------------------------------------------------------------
// Round @x up to next power of two.
//
inline uint32_t
ceilpow2(uint32_t x)
{
  x--;
  x = x | (x >> 1);
  x = x | (x >> 2);
  x = x | (x >> 4);
  x = x | (x >> 8);
  x = x | (x >> 16);
  return x + 1;
}

//---------------------------------------------------------------------------
// Round @x up to next power of two.
//
inline uint64_t
ceilpow2(uint64_t x)
{
  x--;
  x = x | (x >> 1);
  x = x | (x >> 2);
  x = x | (x >> 4);
  x = x | (x >> 8);
  x = x | (x >> 16);
  x = x | (x >> 32);
  return x + 1;
}

//---------------------------------------------------------------------------
// Compute \left\floor \text{log}_2(x) \right\floor.
// NB: ilog2(0) = -1.

inline int
ilog2(uint32_t x)
{
  x = ceilpow2(x + 1) - 1;
  return popcnt(x) - 1;
}

//---------------------------------------------------------------------------
// Compute \left\floor \text{log}_2(x) \right\floor.
// NB: ilog2(0) = -1.

inline int
ilog2(uint64_t x)
{
  x = ceilpow2(x + 1) - 1;
  return popcnt(uint32_t(x >> 32)) + popcnt(uint32_t(x)) - 1;
}

//---------------------------------------------------------------------------
// Compute \left\ceil \text{log}_2(x) \right\ceil.
// NB: ceillog2(0) = 32.

inline int
ceillog2(uint32_t x)
{
  return ilog2(x - 1) + 1;
}

//---------------------------------------------------------------------------
// The number of bits required to represent x.
// NB: x must be >= 0.
// NB: numBits(0) = 1.

inline int
numBits(int x)
{
  return std::max(0, ilog2(uint32_t(x))) + 1;
}

//---------------------------------------------------------------------------
// rotate left by n bits.
// If n is negative, the effect is to rotate right.
// NB: Signed types are treated as unsigned.

template<typename T>
T
rotateLeft(T val, int n)
{
  n &= sizeof(T) * 8 - 1;
  using unsigned_t = typename std::make_unsigned<T>::type;
  return (unsigned_t(val) << n) | (unsigned_t(val) >> (8 * sizeof(val) - n));
}

//---------------------------------------------------------------------------
// rotate right by n bits.
// If n is negative, the effect is to rotate left.
// NB: Signed types are treated as unsigned.

template<typename T>
T
rotateRight(T val, int n)
{
  n &= sizeof(T) * 8 - 1;
  using unsigned_t = typename std::make_unsigned<T>::type;
  return (unsigned_t(val) >> n) | (unsigned_t(val) << (8 * sizeof(val) - n));
}

//---------------------------------------------------------------------------
// Compute an approximation of \left\floor \sqrt{x} \right\floor

uint32_t isqrt(uint64_t x) __attribute__((const));

//---------------------------------------------------------------------------
// Compute an approximation of reciprocal sqrt

uint64_t irsqrt(uint64_t a64) __attribute__((const));

//---------------------------------------------------------------------------
// Compute an approximation of atan2

int iatan2(int y, int x);

//---------------------------------------------------------------------------
// Decrement the @axis-th dimension of 3D morton code @x.
//
inline int64_t
morton3dAxisDec(int64_t val, int axis)
{
  const int64_t mask0 = 0x9249249249249249llu << axis;
  return ((val & mask0) - 1 & mask0) | (val & ~mask0);
}

//---------------------------------------------------------------------------
// add the three dimentional addresses @a + @b;
//
inline uint64_t
morton3dAdd(uint64_t a, uint64_t b)
{
  uint64_t mask = 0x9249249249249249llu;
  uint64_t val = 0;

  for (int i = 0; i < 3; i++) {
    val |= (a | ~mask) + (b & mask) & mask;
    mask <<= 1;
  }

  return val;
}

//---------------------------------------------------------------------------
// Sort the elements in range [@first, @last) using a counting sort.
//
// The value of each element is determined by the function
// @value_of(RandomIt::value_type) and must be in the range [0, Radix).
//
// A supplied output array of @counts represents the histogram of values,
// and may be used to calculate the output position of each value span.
//
// NB: This is an in-place implementation and is not a stable sort.

template<class RandomIt, class ValueOp, std::size_t Radix>
void
countingSort(
  RandomIt first,
  RandomIt last,
  std::array<int, Radix>& counts,
  ValueOp value_of)
{
  // step 1: count each radix
  for (auto it = first; it != last; ++it) {
    counts[value_of(*it)]++;
  }

  // step 2: determine the output offsets
  std::array<RandomIt, Radix> ptrs = {{first}};
  for (int i = 1; i < Radix; i++) {
    ptrs[i] = std::next(ptrs[i - 1], counts[i - 1]);
  }

  // step 3: re-order, completing each radix in turn.
  RandomIt ptr_orig_last = first;
  for (int i = 0; i < Radix; i++) {
    std::advance(ptr_orig_last, counts[i]);
    while (ptrs[i] != ptr_orig_last) {
      int radix = value_of(*ptrs[i]);
      std::iter_swap(ptrs[i], ptrs[radix]);
      ++ptrs[radix];
    }
  }
}

//---------------------------------------------------------------------------

struct NoOp {
  template<typename... Args>
  void operator()(Args...)
  {}
};

//---------------------------------------------------------------------------

template<typename It, typename ValueOp, typename AccumOp>
void
radixSort8WithAccum(int maxValLog2, It begin, It end, ValueOp op, AccumOp acc)
{
  std::array<int, 8> counts = {};
  countingSort(begin, end, counts, [=](decltype(*begin)& it) {
    return op(maxValLog2, it);
  });

  acc(maxValLog2, counts);

  if (--maxValLog2 < 0)
    return;

  auto childBegin = begin;
  for (int i = 0; i < counts.size(); i++) {
    if (!counts[i])
      continue;
    auto childEnd = std::next(childBegin, counts[i]);
    radixSort8WithAccum(maxValLog2, childBegin, childEnd, op, acc);
    childBegin = childEnd;
  }
}

//---------------------------------------------------------------------------

template<typename It, typename ValueOp>
void
radixSort8(int maxValLog2, It begin, It end, ValueOp op)
{
  radixSort8WithAccum(maxValLog2, begin, end, op, NoOp());
}

//============================================================================
// A wrapper to reverse the iteration order of a range based for loop

template<typename T>
struct fwd_is_reverse_iterator {
  T& obj;
};

template<typename T>
auto
begin(fwd_is_reverse_iterator<T> w) -> decltype(w.obj.rbegin())
{
  return w.obj.rbegin();
}

template<typename T>
auto
end(fwd_is_reverse_iterator<T> w) -> decltype(w.obj.rend())
{
  return w.obj.rend();
}

template<typename T>
fwd_is_reverse_iterator<T>
inReverse(T&& obj)
{
  return {obj};
}

//----------------------------------------------------------------------------

}  // namespace pcc
