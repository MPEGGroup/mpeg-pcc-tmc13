/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * <OWNER> = Apple Inc.
 * <ORGANIZATION> = Apple Inc.
 * <YEAR> = 2017
 *
 * Copyright (c) 2017, Apple Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PCCMisc_h
#define PCCMisc_h

#include <vector>

//#define PCC_UNDEFINED_INDEX (uint32_t(-1))

namespace pcc {
const uint32_t PCC_UNDEFINED_INDEX = -1;

enum PCCEndianness { PCC_BIG_ENDIAN = 0, PCC_LITTLE_ENDIAN = 1 };

struct PCCBitstream {
  uint8_t *buffer;
  uint64_t size;
  uint64_t capacity;
};

inline PCCEndianness PCCSystemEndianness() {
  uint32_t num = 1;
  return (*(reinterpret_cast<char *>(&num)) == 1) ? PCC_LITTLE_ENDIAN : PCC_BIG_ENDIAN;
}
template <typename T>
const T PCCEndianSwap(const T u) {
  union {
    T u;
    uint8_t u8[sizeof(T)];
  } source, dest;

  source.u = u;

  for (size_t k = 0; k < sizeof(T); k++) dest.u8[k] = source.u8[sizeof(T) - k - 1];

  return dest.u;
}
template <typename T>
const T PCCToLittleEndian(const T u) {
  return (PCCSystemEndianness() == PCC_BIG_ENDIAN) ? PCCEndianSwap(u) : u;
}
template <typename T>
const T PCCFromLittleEndian(const T u) {
  return (PCCSystemEndianness() == PCC_BIG_ENDIAN) ? PCCEndianSwap(u) : u;
}
template <typename T>
void PCCWriteToBuffer(const T u, uint8_t *const buffer, uint64_t &size) {
  union {
    T u;
    uint8_t u8[sizeof(T)];
  } source;
  source.u = u;
  if (PCCSystemEndianness() == PCC_LITTLE_ENDIAN) {
    for (size_t k = 0; k < sizeof(T); k++) {
      buffer[size++] = source.u8[k];
    }
  } else {
    for (size_t k = 0; k < sizeof(T); k++) {
      buffer[size++] = source.u8[sizeof(T) - k - 1];
    }
  }
}
template <typename T>
void PCCReadFromBuffer(const uint8_t *const buffer, T &u, uint64_t &size) {
  union {
    T u;
    uint8_t u8[sizeof(T)];
  } dest;

  if (PCCSystemEndianness() == PCC_LITTLE_ENDIAN) {
    for (size_t k = 0; k < sizeof(T); k++) {
      dest.u8[k] = buffer[size++];
    }
  } else {
    for (size_t k = 0; k < sizeof(T); k++) {
      dest.u8[sizeof(T) - k - 1] = buffer[size++];
    }
  }
  u = dest.u;
}

void PCCDivideRange(const size_t start, const size_t end, const size_t chunckCount,
                    std::vector<size_t> &subRanges) {
  const size_t elementCount = end - start;
  if (elementCount <= chunckCount) {
    subRanges.resize(elementCount + 1);
    for (size_t i = start; i <= end; ++i) {
      subRanges[i - start] = i;
    }
  } else {
    subRanges.resize(chunckCount + 1);
    const double step = static_cast<double>(elementCount) / (chunckCount + 1);
    double pos = static_cast<double>(start);
    for (size_t i = 0; i < chunckCount; ++i) {
      subRanges[i] = static_cast<size_t>(pos);
      pos += step;
    }
    subRanges[chunckCount] = end;
  }
}
uint32_t PCCGetNumberOfBitsInFixedLengthRepresentation(uint32_t range) {
  if (range <= 1) {
    return 0;
  } else if (range == 2) {
    return 1;
  } else if (range <= 4) {
    return 2;
  } else if (range <= 8) {
    return 3;
  } else if (range <= 16) {
    return 4;
  } else if (range <= 32) {
    return 5;
  } else if (range <= 64) {
    return 6;
  } else if (range <= 128) {
    return 7;
  } else if (range <= 256) {
    return 8;
  } else if (range <= 512) {
    return 9;
  } else if (range <= 1024) {
    return 10;
  } else if (range <= 2048) {
    return 11;
  } else if (range <= 4096) {
    return 12;
  } else if (range <= 8192) {
    return 13;
  } else if (range <= 16384) {
    return 14;
  } else if (range <= 32768) {
    return 15;
  } else if (range <= 65536) {
    return 16;
  } else if (range <= 131072) {
    return 17;
  } else if (range <= 262144) {
    return 18;
  } else if (range <= 524288) {
    return 19;
  } else if (range <= 1048576) {
    return 20;
  } else if (range <= 2097152) {
    return 21;
  } else if (range <= 4194304) {
    return 22;
  } else if (range <= 8388608) {
    return 23;
  } else if (range <= 16777216) {
    return 24;
  } else if (range <= 33554432) {
    return 25;
  } else if (range <= 67108864) {
    return 26;
  } else if (range <= 134217728) {
    return 27;
  } else if (range <= 268435456) {
    return 28;
  } else if (range <= 536870912) {
    return 29;
  } else if (range <= 1073741824) {
    return 30;
  } else if (range <= 2147483648) {
    return 31;
  } else {
    return 32;
  }
}
}

#endif /* PCCMisc_h */
