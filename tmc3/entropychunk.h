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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace pcc {

//=============================================================================
// This multiplexer takes two input streams (one bytewise and one bitwise)
// and assembles them into chunks within an output buffer

class ChunkStreamBuilder {
public:
  ChunkStreamBuilder(const ChunkStreamBuilder&) = delete;
  ChunkStreamBuilder(ChunkStreamBuilder&&) = delete;
  ChunkStreamBuilder& operator=(const ChunkStreamBuilder&) = delete;
  ChunkStreamBuilder& operator=(ChunkStreamBuilder&&) = delete;

  ChunkStreamBuilder() : _chunkBase(nullptr), _outputSizeRemaining(0) {}

  ChunkStreamBuilder(uint8_t* buf, size_t size) { reset(buf, size); }

  void reset(uint8_t* buf = nullptr, size_t size = 0);
  size_t size() const;
  void writeAecByte(uint8_t byte);
  void writeBypassBit(bool bit);
  void flush();

private:
  void reserveChunkByte();
  void finaliseChunk();
  void startNextChunk();

private:
  static const int kChunkSize = 256;

  // number of bytes remaining in the output buffer
  size_t _outputSizeRemaining;

  // number of bytes written to the output
  size_t _outputLength{0};

  // start of the curernt chunk;
  uint8_t* _chunkBase;

  // number of bytes available in the current chunk
  int _chunkBytesRemaining;

  uint8_t* _aecPtr;

  uint8_t* _bypassPtr;

  int _bypassBitIdx;
  int _bypassByteAllocCounter;

  int _bypassBitCounterHack;
};

//=============================================================================

inline void
ChunkStreamBuilder::reset(uint8_t* buf, size_t size)
{
  _outputLength = 0;
  if (!buf)
    return;

  // allocate the first chunk, (fixup the start address first)
  _chunkBase = buf - kChunkSize;
  _outputSizeRemaining = size;
  startNextChunk();
}

//-----------------------------------------------------------------------------

inline size_t
ChunkStreamBuilder::size() const
{
  return _outputLength;
}

//-----------------------------------------------------------------------------

inline void
ChunkStreamBuilder::writeAecByte(uint8_t byte)
{
  reserveChunkByte();
  *_aecPtr++ = byte;
}

//-----------------------------------------------------------------------------

inline void
ChunkStreamBuilder::writeBypassBit(bool bit)
{
  if (_bypassByteAllocCounter < 1) {
    reserveChunkByte();
    _bypassByteAllocCounter += 8;
  }
  _bypassByteAllocCounter--;

  if (--_bypassBitIdx < 0) {
    _bypassPtr--;
    _bypassBitIdx = 7;
  }

  *_bypassPtr = (*_bypassPtr << 1) | bit;
}

//-----------------------------------------------------------------------------

inline void
ChunkStreamBuilder::flush()
{
  if (!_chunkBase)
    return;

  // if nothing has been written to the chunk, remove it
  if (_chunkBytesRemaining == kChunkSize - 1) {
    _outputLength -= kChunkSize;
    return;
  }

  finaliseChunk();

  // if it isn't a full chunk, truncate
  _outputLength -= _chunkBytesRemaining;
}

//-----------------------------------------------------------------------------
// Ensures that there is space to write a byte in the chunk.  If not,
// the current chunk is finalised and starts the next.

inline void
ChunkStreamBuilder::reserveChunkByte()
{
  if (--_chunkBytesRemaining >= 0)
    return;

  // _chunkBytesRemaning is negative: set to zero (since there are none left)
  _chunkBytesRemaining = 0;
  finaliseChunk();
  startNextChunk();

  // reserve a byte from the chunk
  _chunkBytesRemaining--;
}

//-----------------------------------------------------------------------------

inline void
ChunkStreamBuilder::finaliseChunk()
{
  int chunk_num_ae_bytes = _aecPtr - _chunkBase - 1;
  int bypassLen = kChunkSize - _chunkBytesRemaining - chunk_num_ae_bytes - 1;
  if (bypassLen) {
    // the number of padding bits (less the symtax element size)
    int chunk_bypass_num_flushed_bits = _bypassBitIdx - 3;

    // first, add padding bits to current partial byte
    *_bypassPtr <<= _bypassBitIdx;

    // there may be an extra byte at the end if last byte occupancy > 5
    if (chunk_bypass_num_flushed_bits < 0) {
      *--_bypassPtr = 0;
      chunk_bypass_num_flushed_bits += 8;
    }

    *_bypassPtr |= uint8_t(chunk_bypass_num_flushed_bits);
    if (_chunkBytesRemaining)
      std::move(
        _bypassPtr, _chunkBase + kChunkSize,
        _chunkBase + chunk_num_ae_bytes + 1);
  }

  // write out the length of the aec data
  _chunkBase[0] = uint8_t(chunk_num_ae_bytes);
}

//-----------------------------------------------------------------------------

inline void
ChunkStreamBuilder::startNextChunk()
{
  // start a new chunk
  if (_outputSizeRemaining < kChunkSize)
    throw std::runtime_error("Chunk buffer overflow");

  // NB: reserves one byte for the aec length
  _chunkBytesRemaining = kChunkSize - 1;
  _chunkBase += kChunkSize;
  _aecPtr = _chunkBase + 1;
  _bypassPtr = _chunkBase + kChunkSize - 1;
  _bypassBitIdx = 8;
  _bypassByteAllocCounter = -3;

  _outputSizeRemaining -= kChunkSize;
  _outputLength += kChunkSize;

  _bypassBitCounterHack = 0;
}

//=============================================================================

class ChunkStreamReader {
public:
  ChunkStreamReader(const ChunkStreamReader&) = delete;
  ChunkStreamReader(ChunkStreamReader&&) = delete;
  ChunkStreamReader& operator=(const ChunkStreamReader&) = delete;
  ChunkStreamReader& operator=(ChunkStreamReader&&) = delete;

  ChunkStreamReader()
    : _end(nullptr)
    , _aecBytesRemaining(0)
    , _aecNextChunk(nullptr)
    , _bypassNextChunk(nullptr)
    , _bypassAccumBitsRemaining(0)
    , _bypassBitsRemaining(0)
  {}

  ChunkStreamReader(uint8_t* buf, size_t size) { reset(buf, size); }

  void reset(const uint8_t* buf, size_t len);

  uint8_t readAecByte();
  bool readBypassBit();

private:
  static const int kChunkSize = 256;

  // the limit of the buffer.  Used for error checking
  const uint8_t* _end;

  // state for the aec substream
  int _aecBytesRemaining;
  const uint8_t* _aecByte;
  const uint8_t* _aecNextChunk;

  // state for the bypass substream
  const uint8_t* _bypassNextChunk;
  const uint8_t* _bypassByte;
  int _bypassAccumBitsRemaining;
  int _bypassBitsRemaining;
  uint8_t _bypassAccum;
};

//=============================================================================

inline void
ChunkStreamReader::reset(const uint8_t* buf, size_t size)
{
  _end = buf + size;
  _aecBytesRemaining = 0;
  _aecNextChunk = buf;
  _bypassNextChunk = buf;
  _bypassAccumBitsRemaining = 0;
  _bypassBitsRemaining = 0;
}

//-----------------------------------------------------------------------------

inline uint8_t
ChunkStreamReader::readAecByte()
{
  if (_aecBytesRemaining-- > 0)
    return *_aecByte++;

  const uint8_t* ptr = _aecNextChunk;
  int chunk_num_ae_bytes = 0;
  while (ptr < _end && !(chunk_num_ae_bytes = *ptr))
    ptr += kChunkSize;

  if (ptr + chunk_num_ae_bytes >= _end)
    return 0xff;
  //throw std::runtime_error("aec buffer exceeded");

  _aecNextChunk = ptr + kChunkSize;
  _aecByte = ptr + 1;
  _aecBytesRemaining = chunk_num_ae_bytes;
  _aecBytesRemaining--;

  return *_aecByte++;
}

//-----------------------------------------------------------------------------

inline bool
ChunkStreamReader::readBypassBit()
{
  // extract bit from accumulator
  if (_bypassAccumBitsRemaining-- > 0) {
    int bit = !!(_bypassAccum & 0x80);
    _bypassAccum <<= 1;
    return bit;
  }

  // try to refil accumulator
  _bypassBitsRemaining -= 8;
  if (_bypassBitsRemaining > 0) {
    _bypassAccum = *_bypassByte--;
    _bypassAccumBitsRemaining = std::min(_bypassBitsRemaining, 8);
    return readBypassBit();
  }

  // at end of current chunk, find next with bypass data
  const uint8_t* ptr = _bypassNextChunk;
  int chunk_num_ae_bytes = 0;
  while (ptr < _end && (chunk_num_ae_bytes = *ptr) == kChunkSize - 1)
    ptr += kChunkSize;

  // the last chunk may be truncated
  int chunkSize = kChunkSize;
  chunkSize = std::max(0, std::min(int(_end - ptr), chunkSize));

  if (ptr + chunkSize - 1 >= _end)
    throw std::runtime_error("bypass buffer exceeded");

  int chunk_bypass_num_flushed_bits = ptr[chunk_num_ae_bytes + 1] & 0x7;
  _bypassNextChunk = ptr + kChunkSize;
  _bypassByte = ptr + chunkSize - 1;
  _bypassAccum = *_bypassByte--;
  _bypassBitsRemaining =
    8 * (chunkSize - chunk_num_ae_bytes) - chunk_bypass_num_flushed_bits - 11;
  _bypassAccumBitsRemaining = std::min(_bypassBitsRemaining, 8);

  return readBypassBit();
}

//=============================================================================

}  // namespace pcc
