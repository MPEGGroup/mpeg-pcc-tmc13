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
#include <type_traits>

#include "PCCMisc.h"
#include "tables.h"

#include <algorithm>

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

  T& x() { return data[0]; }
  T& y() { return data[1]; }
  T& z() { return data[2]; }
  const T& x() const { return data[0]; }
  const T& y() const { return data[1]; }
  const T& z() const { return data[2]; }

  T& s() { return data[0]; }
  T& t() { return data[1]; }
  T& v() { return data[2]; }
  const T& s() const { return data[0]; }
  const T& t() const { return data[1]; }
  const T& v() const { return data[2]; }

  template<typename ResultT>
  ResultT getNorm2() const
  {
    return Vec3<ResultT>(*this) * Vec3<ResultT>(*this);
  }

  T getNorm1() const
  {
    return std::abs(data[0]) + std::abs(data[1]) + std::abs(data[2]);
  }

  T getNormInf() const
  {
    return std::max(data[2], std::max(abs(data[0]), abs(data[1])));
  }

  // The minimum element
  T min() const { return std::min({data[0], data[1], data[2]}); }

  // The maximum element
  T max() const { return std::max({data[0], data[1], data[2]}); }

  // Applies std::abs to each element
  Vec3 abs() const
  {
    return {std::abs(data[0]), std::abs(data[1]), std::abs(data[2])};
  }

  Vec3& operator=(const Vec3& rhs)
  {
    memcpy(data, rhs.data, sizeof(data));
    return *this;
  }

  template<typename U>
  Vec3& operator+=(const typename pcc::Vec3<U>& rhs)
  {
    data[0] += rhs[0];
    data[1] += rhs[1];
    data[2] += rhs[2];
    return *this;
  }

  template<typename U>
  Vec3& operator-=(const typename pcc::Vec3<U>& rhs)
  {
    data[0] -= rhs[0];
    data[1] -= rhs[1];
    data[2] -= rhs[2];
    return *this;
  }

  template<typename U>
  Vec3& operator-=(U a)
  {
    data[0] -= a;
    data[1] -= a;
    data[2] -= a;
    return *this;
  }

  template<typename U>
  Vec3& operator+=(U a)
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

  Vec3& operator<<=(Vec3<int> val)
  {
    data[0] <<= val[0];
    data[1] <<= val[1];
    data[2] <<= val[2];
    return *this;
  }

  Vec3& operator>>=(int val)
  {
    data[0] >>= val;
    data[1] >>= val;
    data[2] >>= val;
    return *this;
  }

  Vec3& operator>>=(Vec3<int> val)
  {
    data[0] >>= val[0];
    data[1] >>= val[1];
    data[2] >>= val[2];
    return *this;
  }

  template<typename U>
  Vec3& operator/=(U a)
  {
    assert(a != 0);
    data[0] /= a;
    data[1] /= a;
    data[2] /= a;
    return *this;
  }

  template<typename U>
  Vec3& operator*=(U a)
  {
    data[0] *= a;
    data[1] *= a;
    data[2] *= a;
    return *this;
  }

  template<typename U>
  Vec3& operator=(const U a)
  {
    data[0] = a;
    data[1] = a;
    data[2] = a;
    return *this;
  }

  template<typename U>
  Vec3& operator=(const U* rhs)
  {
    data[0] = rhs[0];
    data[1] = rhs[1];
    data[2] = rhs[2];
    return *this;
  }

  Vec3 operator-() const { return Vec3<T>(-data[0], -data[1], -data[2]); }

  template<typename U>
  friend Vec3<typename std::common_type<T, U>::type>
  operator+(const Vec3& lhs, const typename pcc::Vec3<U>& rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
  }

  template<typename U>
  friend Vec3<typename std::enable_if<
    std::is_arithmetic<U>::value,
    typename std::common_type<T, U>::type>::type>
  operator+(const U lhs, const Vec3& rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs + rhs[0], lhs + rhs[1], lhs + rhs[2]);
  }

  template<typename U>
  friend Vec3<typename std::enable_if<
    std::is_arithmetic<U>::value,
    typename std::common_type<T, U>::type>::type>
  operator+(const Vec3& lhs, const U rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs[0] + rhs, lhs[1] + rhs, lhs[2] + rhs);
  }

  template<typename U>
  friend Vec3<typename std::common_type<T, U>::type>
  operator-(const Vec3& lhs, const typename pcc::Vec3<U>& rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
  }

  template<typename U>
  friend Vec3<typename std::enable_if<
    std::is_arithmetic<U>::value,
    typename std::common_type<T, U>::type>::type>
  operator-(const U lhs, const Vec3& rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs - rhs[0], lhs - rhs[1], lhs - rhs[2]);
  }

  template<typename U>
  friend Vec3<typename std::enable_if<
    std::is_arithmetic<U>::value,
    typename std::common_type<T, U>::type>::type>
  operator-(const Vec3& lhs, const U rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs[0] - rhs, lhs[1] - rhs, lhs[2] - rhs);
  }

  template<typename U>
  friend Vec3<typename std::enable_if<
    std::is_arithmetic<U>::value,
    typename std::common_type<T, U>::type>::type>
  operator*(const U lhs, const Vec3& rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs * rhs[0], lhs * rhs[1], lhs * rhs[2]);
  }

  // todo(df): make this elementwise?
  template<typename U>
  friend typename std::common_type<T, U>::type
  operator*(pcc::Vec3<T> lhs, const pcc::Vec3<U>& rhs)
  {
    return (lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]);
  }

  template<typename U>
  friend Vec3<typename std::enable_if<
    std::is_arithmetic<U>::value,
    typename std::common_type<T, U>::type>::type>
  operator*(const Vec3& lhs, const U rhs)
  {
    return Vec3<typename std::common_type<T, U>::type>(
      lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs);
  }

  template<typename U>
  friend Vec3<typename std::enable_if<
    std::is_arithmetic<U>::value,
    typename std::common_type<T, U>::type>::type>
  operator/(const Vec3& lhs, const U rhs)
  {
    assert(rhs != 0);
    return Vec3<typename std::common_type<T, U>::type>(
      lhs[0] / rhs, lhs[1] / rhs, lhs[2] / rhs);
  }

  friend Vec3 operator<<(const Vec3& lhs, int val)
  {
    return Vec3<T>(lhs[0] << val, lhs[1] << val, lhs[2] << val);
  }

  friend Vec3 operator<<(const Vec3& lhs, const Vec3<int>& val)
  {
    return Vec3<T>(lhs[0] << val[0], lhs[1] << val[1], lhs[2] << val[2]);
  }

  friend Vec3 operator>>(const Vec3& lhs, int val)
  {
    return Vec3<T>(lhs[0] >> val, lhs[1] >> val, lhs[2] >> val);
  }

  friend Vec3 operator>>(const Vec3& lhs, const Vec3<int>& val)
  {
    return Vec3<T>(lhs[0] >> val[0], lhs[1] >> val[1], lhs[2] >> val[2]);
  }

  bool operator<(const Vec3& rhs) const
  {
    if (data[0] == rhs[0]) {
      if (data[1] == rhs[1]) {
        return (data[2] < rhs[2]);
      }
      return (data[1] < rhs[1]);
    }
    return (data[0] < rhs[0]);
  }

  bool operator>(const Vec3& rhs) const
  {
    if (data[0] == rhs[0]) {
      if (data[1] == rhs[1]) {
        return (data[2] > rhs[2]);
      }
      return (data[1] > rhs[1]);
    }
    return (data[0] > rhs[0]);
  }

  bool operator==(const Vec3& rhs) const
  {
    return (data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2]);
  }

  bool operator!=(const Vec3& rhs) const
  {
    return (data[0] != rhs[0] || data[1] != rhs[1] || data[2] != rhs[2]);
  }

  friend std::ostream& operator<<(std::ostream& os, const Vec3& vec)
  {
    os << vec[0] << " " << vec[1] << " " << vec[2];
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

  template<typename U>
  Vec3(const typename pcc::Vec3<U>& vec)
  {
    data[0] = vec.data[0];
    data[1] = vec.data[1];
    data[2] = vec.data[2];
  }

  Vec3() = default;
  ~Vec3(void) = default;

private:
  T data[3];

  template<typename U>
  friend class pcc::Vec3;
};

template<typename T>
struct Box3 {
  Vec3<T> min;
  Vec3<T> max;

  Box3() = default;

  Box3(T min, T max) : min(min), max(max) {}

  Box3(const Vec3<T>& min, const Vec3<T>& max) : min(min), max(max) {}

  template<typename ForwardIt>
  Box3(ForwardIt begin, ForwardIt end)
    : Box3(std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest())
  {
    for (auto it = begin; it != end; ++it) {
      auto& pt = *it;
      for (int k = 0; k < 3; ++k) {
        if (pt[k] > max[k])
          max[k] = pt[k];
        if (pt[k] < min[k])
          min[k] = pt[k];
      }
    }
  }

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

  template<typename SquaredT>
  SquaredT getDist2(const Vec3<T>& point) const
  {
    using U = SquaredT;
    U dx = U(std::max(std::max(min[0] - point[0], T()), point[0] - max[0]));
    U dy = U(std::max(std::max(min[1] - point[1], T()), point[1] - max[1]));
    U dz = U(std::max(std::max(min[2] - point[2], T()), point[2] - max[2]));
    return dx * dx + dy * dy + dz * dz;
  }

  T getDist1(const Vec3<T>& point) const
  {
    T dx = T(std::max(std::max(min[0] - point[0], T()), point[0] - max[0]));
    T dy = T(std::max(std::max(min[1] - point[1], T()), point[1] - max[1]));
    T dz = T(std::max(std::max(min[2] - point[2], T()), point[2] - max[2]));
    return dx + dy + dz;
  }

  void insert(const Vec3<T>& point)
  {
    min.x() = std::min(min.x(), point.x());
    min.y() = std::min(min.y(), point.y());
    min.z() = std::min(min.z(), point.z());
    max.x() = std::max(max.x(), point.x());
    max.y() = std::max(max.y(), point.y());
    max.z() = std::max(max.z(), point.z());
  }

  friend std::ostream& operator<<(std::ostream& os, const Box3& box)
  {
    os << box.min[0] << " " << box.min[1] << " " << box.min[2] << " "
       << box.max[0] << " " << box.max[1] << " " << box.max[2];
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
// element-wise multiplication of two vectors

template<typename T, typename U>
Vec3<typename std::common_type<T, U>::type>
times(Vec3<T> lhs, const Vec3<U>& rhs)
{
  return Vec3<typename std::common_type<T, U>::type>(
    lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2]);
}

//---------------------------------------------------------------------------

typedef DEPRECATED_MSVC Vec3<double> PCCVector3D DEPRECATED;
typedef DEPRECATED_MSVC Vec3<double> PCCPoint3D DEPRECATED;
typedef DEPRECATED_MSVC Box3<double> PCCBox3D DEPRECATED;
typedef DEPRECATED_MSVC Vec3<uint8_t> PCCColor3B DEPRECATED;

template<typename T>
using PCCVector3 DEPRECATED = Vec3<T>;

//===========================================================================

struct Rational {
  int numerator;
  int denominator;

  Rational() : Rational(0, 1){};

  Rational(int numerator) : Rational(numerator, 1){};

  Rational(int numerator, int denominator)
    : numerator(numerator), denominator(denominator)
  {}

  Rational(float val);
  Rational(double val);

  operator double() const { return double(numerator) / double(denominator); }

  operator float() const { return float(numerator) / float(denominator); }
};

//---------------------------------------------------------------------------

inline Rational
reciprocal(const Rational x)
{
  return Rational(x.denominator, x.numerator);
}

//===========================================================================

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
// Convert a vector position to morton order address.

template<typename T>
int64_t
mortonAddr(const Vec3<T>& vec)
{
  return mortonAddr(int(vec.x()), int(vec.y()), int(vec.z()));
}

//---------------------------------------------------------------------------

inline int64_t
PCCClip(const int64_t& n, const int64_t& lower, const int64_t& upper)
{
  return std::max(lower, std::min(n, upper));
}

//---------------------------------------------------------------------------
// Integer division of @x by 2^shift, truncating towards zero.

template<typename T>
inline T
divExp2(T x, int shift)
{
  return x >= 0 ? x >> shift : -(-x >> shift);
}

//---------------------------------------------------------------------------
// Integer division of @x by 2^shift, rounding intermediate half values
// to +Inf.

inline int64_t
divExp2RoundHalfUp(int64_t x, int shift)
{
  if (!shift)
    return x;

  int64_t half = 1ll << (shift - 1);
  return (x + half) >> shift;
}

//---------------------------------------------------------------------------
// Integer division of @scalar by 2^shift, rounding intermediate half values
// away from zero.

inline int64_t
divExp2RoundHalfInf(int64_t scalar, int shift)
{
  if (!shift)
    return scalar;

  int64_t s0 = 1ll << (shift - 1);
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

  return ((1ull << (shift - 1)) + scalar) >> shift;
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

extern const uint16_t kDivApproxDivisor[256];

//---------------------------------------------------------------------------

inline int64_t
divInvDivisorApprox(const uint64_t b, int32_t& log2InvScale)
{
  assert(b > 0);
  const int32_t lutSizeLog2 = 8;

  const auto n = std::max(0, ilog2(b) + 1 - lutSizeLog2);
  const auto index = (b + ((1ull << n) >> 1)) >> n;
  assert(unsigned(index) <= (1 << lutSizeLog2));
  log2InvScale = n + (lutSizeLog2 << 1);
  return kDivApproxDivisor[index - 1] + 1;
}

//---------------------------------------------------------------------------

inline int64_t
divApprox(const int64_t a, const uint64_t b, const int32_t log2Scale)
{
  assert(abs(a) < (1ull << 46));
  int32_t log2InvScale;
  const int64_t invB = divInvDivisorApprox(b, log2InvScale);
  return (invB * a) >> (log2InvScale - log2Scale);
}

//---------------------------------------------------------------------------

template<unsigned NIter = 1>
inline int64_t
recipApprox(int64_t b, int32_t& log2Scale)
{
  int log2ScaleOffset = 0;
  int32_t log2bPlusOne = ilog2(uint64_t(b)) + 1;

  if (log2bPlusOne > 31) {
    b >>= log2bPlusOne - 31;
    log2ScaleOffset -= log2bPlusOne - 31;
  }

  if (log2bPlusOne < 31) {
    b <<= 31 - log2bPlusOne;
    log2ScaleOffset += 31 - log2bPlusOne;
  }

  // Initial approximation: 48/17 - 32/17 * b with 28 bits decimal prec
  int64_t bRecip = ((0x2d2d2d2dLL << 31) - 0x1e1e1e1eLL * b) >> 28;
  for (unsigned i = 0; i < NIter; ++i)
    bRecip += bRecip * ((1LL << 31) - (b * bRecip >> 31)) >> 31;

  log2Scale = (31 << 1) - log2ScaleOffset;
  return bRecip;
}

//---------------------------------------------------------------------------

inline Vec3<int64_t>
divApprox(const Vec3<int64_t> a, const uint64_t b, const int32_t log2Scale)
{
  assert(abs(a[0]) < (1ull << 46));
  assert(abs(a[1]) < (1ull << 46));
  assert(abs(a[2]) < (1ull << 46));
  int32_t log2InvScale;
  const int64_t invB = divInvDivisorApprox(b, log2InvScale);
  const int32_t n = log2InvScale - log2Scale;
  Vec3<int64_t> res;
  res[0] = (invB * a[0]) >> n;
  res[1] = (invB * a[1]) >> n;
  res[2] = (invB * a[2]) >> n;
  return res;
}

//---------------------------------------------------------------------------

inline Vec3<int64_t>
divApproxRoundHalfInf(
  const Vec3<int64_t> a, const uint64_t b, const int32_t log2Scale)
{
  assert(abs(a[0]) < (1ull << 46));
  assert(abs(a[1]) < (1ull << 46));
  assert(abs(a[2]) < (1ull << 46));
  int32_t log2InvScale;
  const int64_t invB = divInvDivisorApprox(b, log2InvScale);
  const int32_t n = log2InvScale - log2Scale;
  Vec3<int64_t> res;
  res[0] = divExp2RoundHalfInf(invB * a[0], n);
  res[1] = divExp2RoundHalfInf(invB * a[1], n);
  res[2] = divExp2RoundHalfInf(invB * a[2], n);
  return res;
}

//---------------------------------------------------------------------------

inline int32_t
isin0(const int32_t x, const int32_t log2Scale)
{
  assert(log2Scale >= kLog2ISineAngleScale);
  assert(x >= 0);
  assert(x <= (1 << (log2Scale - 2)));
  const auto ds = log2Scale - kLog2ISineAngleScale;
  const auto b = (1 << ds);
  const auto i0 = (x >> ds);
  const auto x0 = i0 << ds;
  const auto d1 = x - x0;
  assert(i0 <= (1 << kLog2ISineAngleScale) >> 2);
  return kISine[i0] + ((d1 * (kISine[i0 + 1] - kISine[i0]) + (b >> 1)) >> ds);
}

//---------------------------------------------------------------------------

inline int32_t
icos0(const int32_t x, const int32_t log2Scale)
{
  assert(x >= 0);
  assert(x <= (1 << (log2Scale - 2)));
  return isin0((1 << (log2Scale - 2)) - x, log2Scale);
}

//---------------------------------------------------------------------------

inline int32_t
isin(int32_t x, const int32_t log2Scale)
{
  const auto L = 1 << (log2Scale - 1);
  x = std::min(std::max(x, -L), L);
  assert(abs(x) <= (1 << (log2Scale - 1)));
  const auto Q0 = 1 << (log2Scale - 2);
  if (x >= Q0) {
    return isin0((1 << (log2Scale - 1)) - x, log2Scale);
  } else if (x >= 0) {
    return isin0(x, log2Scale);
  } else if (x >= -Q0) {
    return -isin0(-x, log2Scale);
  } else {
    return -isin0((1 << (log2Scale - 1)) + x, log2Scale);
  }
}

//---------------------------------------------------------------------------

inline int32_t
icos(int32_t x, const int32_t log2Scale)
{
  const auto Q0 = 1 << (log2Scale - 2);
  const auto ax = std::min(abs(x), (1 << (log2Scale - 1)));
  return ax <= Q0 ? icos0(ax, log2Scale)
                  : -icos0((1 << (log2Scale - 1)) - ax, log2Scale);
}

//---------------------------------------------------------------------------

} /* namespace pcc */

#endif /* PCCMath_h */
