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

#ifndef PCCMath_h
#define PCCMath_h

#include <assert.h>
#include <cstddef>
#include <iostream>
#include <limits>
#include <math.h>
#include <string.h>

#include "PCCMisc.h"
#include "tables.h"

namespace pcc {
/// Vector dim 3
template<typename T>
class Vec3 {
public:
  T* begin() { return &data[0]; }
  const T* begin() const { return &data[0]; }

  T* end() { return &data[3]; }
  const T* end() const { return &data[3]; }

  T& operator[](size_t i)
  {
    assert(i < 3);
    return data[i];
  }
  const T& operator[](size_t i) const
  {
    assert(i < 3);
    return data[i];
  }
  size_t getElementCount() const { return 3; }
  T& r() { return data[0]; }
  T& g() { return data[1]; }
  T& b() { return data[2]; }
  const T& r() const { return data[0]; }
  const T& g() const { return data[1]; }
  const T& b() const { return data[2]; }
  T& x() { return data[0]; }
  T& y() { return data[1]; }
  T& z() { return data[2]; }
  const T& x() const { return data[0]; }
  const T& y() const { return data[1]; }
  const T& z() const { return data[2]; }

  T getNorm() const { return static_cast<T>(sqrt(getNorm2())); }
  T getNorm2() const { return (*this) * (*this); }
  T getNormInf() const
  {
    return std::max(data[2], std::max(abs(data[0]), abs(data[1])));
  }
  Vec3& operator=(const Vec3& rhs)
  {
    memcpy(data, rhs.data, sizeof(data));
    return *this;
  }
  Vec3& operator+=(const Vec3& rhs)
  {
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    return *this;
  }
  Vec3& operator-=(const Vec3& rhs)
  {
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    return *this;
  }
  Vec3& operator-=(const T a)
  {
    data[0] -= a;
    data[1] -= a;
    data[2] -= a;
    return *this;
  }
  Vec3& operator+=(const T a)
  {
    data[0] += a;
    data[1] += a;
    data[2] += a;
    return *this;
  }
  Vec3& operator<<=(int val)
  {
    data[0] <<= val;
    data[1] <<= val;
    data[2] <<= val;
    return *this;
  }
  Vec3& operator>>=(int val)
  {
    data[0] >>= val;
    data[1] >>= val;
    data[2] >>= val;
    return *this;
  }
  Vec3& operator/=(const T a)
  {
    assert(a != 0);
    data[0] /= a;
    data[1] /= a;
    data[2] /= a;
    return *this;
  }
  Vec3& operator*=(const T a)
  {
    data[0] *= a;
    data[1] *= a;
    data[2] *= a;
    return *this;
  }
  Vec3& operator=(const T a)
  {
    data[0] = a;
    data[1] = a;
    data[2] = a;
    return *this;
  }
  Vec3& operator=(const T* const rhs)
  {
    data[0] = rhs[0];
    data[1] = rhs[1];
    data[2] = rhs[2];
    return *this;
  }
  T operator*(const Vec3& rhs) const
  {
    return (
      data[0] * rhs.data[0] + data[1] * rhs.data[1] + data[2] * rhs.data[2]);
  }
  Vec3 operator-() const { return Vec3<T>(-data[0], -data[1], -data[2]); }
  friend Vec3 operator+(const Vec3& lhs, const Vec3& rhs)
  {
    return Vec3<T>(
      lhs.data[0] + rhs.data[0], lhs.data[1] + rhs.data[1],
      lhs.data[2] + rhs.data[2]);
  }
  friend Vec3 operator+(const T lhs, const Vec3& rhs)
  {
    return Vec3<T>(lhs + rhs.data[0], lhs + rhs.data[1], lhs + rhs.data[2]);
  }
  friend Vec3 operator+(const Vec3& lhs, const T rhs)
  {
    return Vec3<T>(lhs.data[0] + rhs, lhs.data[1] + rhs, lhs.data[2] + rhs);
  }
  friend Vec3 operator-(const Vec3& lhs, const Vec3& rhs)
  {
    return Vec3<T>(
      lhs.data[0] - rhs.data[0], lhs.data[1] - rhs.data[1],
      lhs.data[2] - rhs.data[2]);
  }
  friend Vec3 operator-(const T lhs, const Vec3& rhs)
  {
    return Vec3<T>(lhs - rhs.data[0], lhs - rhs.data[1], lhs - rhs.data[2]);
  }
  friend Vec3 operator-(const Vec3& lhs, const T rhs)
  {
    return Vec3<T>(lhs.data[0] - rhs, lhs.data[1] - rhs, lhs.data[2] - rhs);
  }
  friend Vec3 operator*(const T lhs, const Vec3& rhs)
  {
    return Vec3<T>(lhs * rhs.data[0], lhs * rhs.data[1], lhs * rhs.data[2]);
  }
  friend Vec3 operator*(const Vec3& lhs, const T rhs)
  {
    return Vec3<T>(lhs.data[0] * rhs, lhs.data[1] * rhs, lhs.data[2] * rhs);
  }
  friend Vec3 operator/(const Vec3& lhs, const T rhs)
  {
    assert(rhs != 0);
    return Vec3<T>(lhs.data[0] / rhs, lhs.data[1] / rhs, lhs.data[2] / rhs);
  }
  friend Vec3 operator<<(const Vec3& lhs, int val)
  {
    return Vec3<T>(lhs.data[0] << val, lhs.data[1] << val, lhs.data[2] << val);
  }
  friend Vec3 operator>>(const Vec3& lhs, int val)
  {
    return Vec3<T>(lhs.data[0] >> val, lhs.data[1] >> val, lhs.data[2] >> val);
  }
  bool operator<(const Vec3& rhs) const
  {
    if (data[0] == rhs.data[0]) {
      if (data[1] == rhs.data[1]) {
        return (data[2] < rhs.data[2]);
      }
      return (data[1] < rhs.data[1]);
    }
    return (data[0] < rhs.data[0]);
  }
  bool operator>(const Vec3& rhs) const
  {
    if (data[0] == rhs.data[0]) {
      if (data[1] == rhs.data[1]) {
        return (data[2] > rhs.data[2]);
      }
      return (data[1] > rhs.data[1]);
    }
    return (data[0] > rhs.data[0]);
  }
  bool operator==(const Vec3& rhs) const
  {
    return (
      data[0] == rhs.data[0] && data[1] == rhs.data[1]
      && data[2] == rhs.data[2]);
  }
  bool operator!=(const Vec3& rhs) const
  {
    return (
      data[0] != rhs.data[0] || data[1] != rhs.data[1]
      || data[2] != rhs.data[2]);
  }
  friend std::ostream& operator<<(std::ostream& os, const Vec3& vec)
  {
    os << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
    return os;
  }
  friend std::istream& operator>>(std::istream& is, Vec3& vec)
  {
    is >> vec[0] >> vec[1] >> vec[2];
    return is;
  }
  Vec3(const T a) { data[0] = data[1] = data[2] = a; }
  Vec3(const T x, const T y, const T z)
  {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  Vec3(const Vec3& vec)
  {
    data[0] = vec.data[0];
    data[1] = vec.data[1];
    data[2] = vec.data[2];
  }
  Vec3() = default;
  ~Vec3(void) = default;

private:
  T data[3];
};

template<typename T>
struct Box3 {
  Vec3<T> min;
  Vec3<T> max;
  bool contains(const Vec3<T> point) const
  {
    return !(
      point.x() < min.x() || point.x() > max.x() || point.y() < min.y()
      || point.y() > max.y() || point.z() < min.z() || point.z() > max.z());
  }

  Box3 merge(const Box3& box)
  {
    min.x() = std::min(min.x(), box.min.x());
    min.y() = std::min(min.y(), box.min.y());
    min.z() = std::min(min.z(), box.min.z());
    max.x() = std::max(max.x(), box.max.x());
    max.y() = std::max(max.y(), box.max.y());
    max.z() = std::max(max.z(), box.max.z());
    return box;
  }

  bool intersects(const Box3& box) const
  {
    return max.x() >= box.min.x() && min.x() <= box.max.x()
      && max.y() >= box.min.y() && min.y() <= box.max.y()
      && max.z() >= box.min.z() && min.z() <= box.max.z();
  }

  T getDist2(const Vec3<T>& point) const
  {
    const T dx = std::max(std::max(min[0] - point[0], 0.0), point[0] - max[0]);
    const T dy = std::max(std::max(min[1] - point[1], 0.0), point[1] - max[1]);
    const T dz = std::max(std::max(min[2] - point[2], 0.0), point[2] - max[2]);
    return dx * dx + dy * dy + dz * dz;
  }

  friend std::ostream& operator<<(std::ostream& os, const Box3& box)
  {
    os << box.min[0] << " " << box.min[1] << " " << box.min[2] << " "
       << box.max[0] << " " << box.max[1] << " " << box.max[2] << std::endl;
    return os;
  }
  friend std::istream& operator>>(std::istream& is, Box3& box)
  {
    is >> box.min[0] >> box.min[1] >> box.min[2] >> box.max[0] >> box.max[1]
      >> box.max[2];
    return is;
  }
};

//---------------------------------------------------------------------------

typedef DEPRECATED_MSVC Vec3<double> PCCVector3D DEPRECATED;
typedef DEPRECATED_MSVC Vec3<double> PCCPoint3D DEPRECATED;
typedef DEPRECATED_MSVC Box3<double> PCCBox3D DEPRECATED;
typedef DEPRECATED_MSVC Vec3<uint8_t> PCCColor3B DEPRECATED;

template<typename T>
using PCCVector3 DEPRECATED = Vec3<T>;

//---------------------------------------------------------------------------

template<typename T>
T
PCCClip(const T& n, const T& lower, const T& upper)
{
  return std::max(lower, std::min(n, upper));
}
template<typename T>
bool
PCCApproximatelyEqual(
  T a, T b, T epsilon = std::numeric_limits<double>::epsilon())
{
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

//---------------------------------------------------------------------------

inline int64_t
mortonAddr(const int32_t x, const int32_t y, const int32_t z)
{
  assert(x >= 0 && y >= 0 && z >= 0);
  int64_t answer = kMortonCode256X[(x >> 16) & 0xFF]
    | kMortonCode256Y[(y >> 16) & 0xFF] | kMortonCode256Z[(z >> 16) & 0xFF];
  answer = answer << 24 | kMortonCode256X[(x >> 8) & 0xFF]
    | kMortonCode256Y[(y >> 8) & 0xFF] | kMortonCode256Z[(z >> 8) & 0xFF];
  answer = answer << 24 | kMortonCode256X[x & 0xFF] | kMortonCode256Y[y & 0xFF]
    | kMortonCode256Z[z & 0xFF];
  return answer;
}

//---------------------------------------------------------------------------
// Convert a vector position (divided by 2^depth) to morton order address.

template<typename T>
int64_t
mortonAddr(const Vec3<T>& vec, int depth)
{
  int x = int(vec.x()) >> depth;
  int y = int(vec.y()) >> depth;
  int z = int(vec.z()) >> depth;
  return mortonAddr(x, y, z);
}

//---------------------------------------------------------------------------

inline int64_t
PCCClip(const int64_t& n, const int64_t& lower, const int64_t& upper)
{
  return std::max(lower, std::min(n, upper));
}

//---------------------------------------------------------------------------
// Integer division of @scalar by 2^shift, rounding intermediate half values
// away from zero.

inline int64_t
divExp2RoundHalfInf(int64_t scalar, int shift)
{
  if (!shift)
    return scalar;

  int64_t s0 = 1 << (shift - 1);
  return scalar >= 0 ? (s0 + scalar) >> shift : -((s0 - scalar) >> shift);
}

//---------------------------------------------------------------------------
// Integer division of @scalar by 2^shift, rounding intermediate half values
// away from zero.

inline uint64_t
divExp2RoundHalfInf(uint64_t scalar, int shift)
{
  if (!shift)
    return scalar;

  return ((1 << (shift - 1)) + scalar) >> shift;
}

//---------------------------------------------------------------------------
// Component-wise integer division of @vec by 2^shift, rounding intermediate
// half values away from zero.

template<typename T>
inline Vec3<T>
divExp2RoundHalfInf(Vec3<T> vec, int shift)
{
  for (int k = 0; k < 3; ++k)
    vec[k] = divExp2RoundHalfInf(vec[k], shift);

  return vec;
}

//---------------------------------------------------------------------------

} /* namespace pcc */

#endif /* PCCMath_h */
