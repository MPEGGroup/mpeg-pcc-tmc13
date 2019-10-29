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

#ifndef PCCPointSet_h
#define PCCPointSet_h

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "PCCMath.h"
#include "PCCMisc.h"

namespace pcc {
class PCCPointSet3 {
public:
  typedef Vec3<double> PointType;

  //=========================================================================
  // proxy object for use with iterator, allowing handling of PCCPointSet3's
  // structure-of-arrays as a single array.

  class iterator;
  class Proxy {
    friend class iterator;

    PCCPointSet3* parent_;
    size_t idx_;

  public:
    //-----------------------------------------------------------------------

    Proxy() : parent_(nullptr), idx_() {}

    Proxy(PCCPointSet3* parent, size_t idx) : parent_(parent), idx_(idx) {}

    //-----------------------------------------------------------------------

    PointType operator*() const { return (*parent_)[idx_]; }

    PointType& operator*() { return (*parent_)[idx_]; }

    //-----------------------------------------------------------------------
    // Swap the position of the current proxied point (including attributes)
    // with that of @other in the same PointSet.

    void swap(const Proxy& other) const
    {
      assert(parent_ == other.parent_);
      parent_->swapPoints(idx_, other.idx_);
    }

    //-----------------------------------------------------------------------
  };

  //=========================================================================
  // iterator for use with stl algorithms

  class iterator {
  private:
    Proxy p_;

  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef const Proxy value_type;
    typedef std::ptrdiff_t difference_type;
    typedef const Proxy* pointer;
    typedef const Proxy& reference;

    //-----------------------------------------------------------------------

    iterator() = default;
    iterator(const iterator&) = default;

    //-----------------------------------------------------------------------

    explicit iterator(PCCPointSet3* parent) : p_{parent, 0} {}

    explicit iterator(PCCPointSet3* parent, size_t idx) : p_{parent, idx} {}

    //-----------------------------------------------------------------------
    // :: Iterator

    reference operator*() const { return p_; }

    //-----------------------------------------------------------------------

    iterator& operator++()
    {
      p_.idx_++;
      return *this;
    }

    //-----------------------------------------------------------------------
    // :: ForwardIterator

    iterator operator++(int)
    {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    //-----------------------------------------------------------------------

    pointer operator->() const { return &p_; }

    //-----------------------------------------------------------------------

    bool operator==(const iterator& other) const
    {
      return p_.idx_ == other.p_.idx_;
    }

    //-----------------------------------------------------------------------

    bool operator!=(const iterator& other) const { return !(*this == other); }

    //-----------------------------------------------------------------------
    // :: BidirectionalIterator

    iterator& operator--()
    {
      p_.idx_--;
      return *this;
    }

    //-----------------------------------------------------------------------

    iterator operator--(int)
    {
      iterator retval = *this;
      --(*this);
      return retval;
    }

    //-----------------------------------------------------------------------
    // :: RandomAccessIterator

    value_type operator[](difference_type n)
    {
      return Proxy{p_.parent_, p_.idx_ + n};
    }

    //-----------------------------------------------------------------------

    iterator& operator+=(difference_type n)
    {
      p_.idx_ += n;
      return *this;
    }

    //-----------------------------------------------------------------------

    iterator operator+(difference_type n) const
    {
      iterator it(*this);
      it += n;
      return it;
    }

    //-----------------------------------------------------------------------

    iterator& operator-=(difference_type n)
    {
      p_.idx_ -= n;
      return *this;
    }

    //-----------------------------------------------------------------------

    iterator operator-(difference_type n) const
    {
      iterator it(*this);
      it -= n;
      return it;
    }

    //-----------------------------------------------------------------------

    difference_type operator-(const iterator& other) const
    {
      return p_.idx_ - other.p_.idx_;
    }

    //-----------------------------------------------------------------------
  };

  //=========================================================================

  PCCPointSet3()
  {
    withColors = false;
    withReflectances = false;
    withFrameIndex = false;
  }
  PCCPointSet3(const PCCPointSet3&) = default;
  PCCPointSet3& operator=(const PCCPointSet3& rhs) = default;
  ~PCCPointSet3() = default;

  void swap(PCCPointSet3& other)
  {
    using std::swap;
    swap(positions, other.positions);
    swap(colors, other.colors);
    swap(reflectances, other.reflectances);
    swap(frameidx, other.frameidx);
    swap(withColors, other.withColors);
    swap(withReflectances, other.withReflectances);
    swap(withFrameIndex, other.withFrameIndex);
  }

  Vec3<double> operator[](const size_t index) const
  {
    assert(index < positions.size());
    return positions[index];
  }
  Vec3<double>& operator[](const size_t index)
  {
    assert(index < positions.size());
    return positions[index];
  }
  Vec3<uint8_t> getColor(const size_t index) const
  {
    assert(index < colors.size() && withColors);
    return colors[index];
  }
  Vec3<uint8_t>& getColor(const size_t index)
  {
    assert(index < colors.size() && withColors);
    return colors[index];
  }
  void setColor(const size_t index, const Vec3<uint8_t> color)
  {
    assert(index < colors.size() && withColors);
    colors[index] = color;
  }
  uint16_t getReflectance(const size_t index) const
  {
    assert(index < reflectances.size() && withReflectances);
    return reflectances[index];
  }
  uint16_t& getReflectance(const size_t index)
  {
    assert(index < reflectances.size() && withReflectances);
    return reflectances[index];
  }
  void setReflectance(const size_t index, const uint16_t reflectance)
  {
    assert(index < reflectances.size() && withReflectances);
    reflectances[index] = reflectance;
  }

  bool hasReflectances() const { return withReflectances; }
  void addReflectances()
  {
    withReflectances = true;
    resize(getPointCount());
  }
  void removeReflectances()
  {
    withReflectances = false;
    reflectances.resize(0);
  }

  uint8_t getFrameIndex(const size_t index) const
  {
    assert(index < frameidx.size() && withFrameIndex);
    return frameidx[index];
  }

  uint8_t& getFrameIndex(const size_t index)
  {
    assert(index < frameidx.size() && withFrameIndex);
    return frameidx[index];
  }

  void setFrameIndex(const size_t index, const uint8_t frameindex)
  {
    assert(index < frameidx.size() && withFrameIndex);
    frameidx[index] = frameindex;
  }

  bool hasFrameIndex() const { return withFrameIndex; }
  void addFrameIndex()
  {
    withFrameIndex = true;
    resize(getPointCount());
  }
  void removeFrameIndex()
  {
    withFrameIndex = false;
    frameidx.resize(0);
  }

  bool hasColors() const { return withColors; }
  void addColors()
  {
    withColors = true;
    resize(getPointCount());
  }
  void removeColors()
  {
    withColors = false;
    colors.resize(0);
  }

  void addRemoveAttributes(bool withColors, bool withReflectances)
  {
    if (withColors)
      addColors();
    else
      removeColors();

    if (withReflectances)
      addReflectances();
    else
      removeReflectances();
  }

  size_t getPointCount() const { return positions.size(); }
  void resize(const size_t size)
  {
    positions.resize(size);
    if (hasColors()) {
      colors.resize(size);
    }
    if (hasReflectances()) {
      reflectances.resize(size);
    }
    if (hasFrameIndex()) {
      frameidx.resize(size);
    }
  }

  void reserve(const size_t size)
  {
    positions.reserve(size);
    if (hasColors()) {
      colors.reserve(size);
    }
    if (hasReflectances()) {
      reflectances.reserve(size);
    }
    if (hasFrameIndex()) {
      frameidx.reserve(size);
    }
  }
  void clear()
  {
    positions.clear();
    colors.clear();
    reflectances.clear();
    frameidx.clear();
  }

  size_t removeDuplicatePointInQuantizedPoint(int minGeomNodeSizeLog2)
  {
    for (int i = 0; i < positions.size(); i++) {
      Vec3<double> newPoint = positions[i];
      if (minGeomNodeSizeLog2 > 0) {
        uint32_t mask = ((uint32_t)-1) << minGeomNodeSizeLog2;
        positions[i].x() = ((int32_t)(positions[i].x()) & mask);
        positions[i].y() = ((int32_t)(positions[i].y()) & mask);
        positions[i].z() = ((int32_t)(positions[i].z()) & mask);
      }
    }
    positions.erase(
      std::unique(positions.begin(), positions.end()), positions.end());

    return positions.size();
  }

  void append(const PCCPointSet3& src)
  {
    if (!getPointCount())
      addRemoveAttributes(src.hasColors(), src.hasReflectances());

    int dstEnd = positions.size();
    int srcSize = src.positions.size();
    resize(dstEnd + srcSize);

    std::copy(
      src.positions.begin(), src.positions.end(),
      std::next(positions.begin(), dstEnd));

    if (hasColors() && src.hasColors())
      std::copy(
        src.colors.begin(), src.colors.end(),
        std::next(colors.begin(), dstEnd));

    if (hasReflectances() && src.hasReflectances())
      std::copy(
        src.reflectances.begin(), src.reflectances.end(),
        std::next(reflectances.begin(), dstEnd));
  }

  void swapPoints(const size_t index1, const size_t index2)
  {
    assert(index1 < getPointCount());
    assert(index2 < getPointCount());
    std::swap((*this)[index1], (*this)[index2]);
    if (hasColors()) {
      std::swap(getColor(index1), getColor(index2));
    }
    if (hasReflectances()) {
      std::swap(getReflectance(index1), getReflectance(index2));
    }
  }

  Box3<double> computeBoundingBox() const
  {
    Box3<double> bbox = {std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::lowest()};
    const size_t pointCount = getPointCount();
    for (size_t i = 0; i < pointCount; ++i) {
      const Vec3<double>& pt = (*this)[i];
      for (int k = 0; k < 3; ++k) {
        if (pt[k] > bbox.max[k]) {
          bbox.max[k] = pt[k];
        }
        if (pt[k] < bbox.min[k]) {
          bbox.min[k] = pt[k];
        }
      }
    }
    return bbox;
  }

  //--------------------------------------------------------------------------
  // Determine the bounding box of the set of points given by the indicies
  // given by iterating over [begin, end)

  template<typename ForwardIt>
  Box3<double> computeBoundingBox(ForwardIt begin, ForwardIt end) const
  {
    Box3<double> bbox = {std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::lowest()};

    for (auto it = begin; it != end; ++it) {
      int i = *it;
      const Vec3<double>& pt = (*this)[i];
      for (int k = 0; k < 3; ++k) {
        if (pt[k] > bbox.max[k]) {
          bbox.max[k] = pt[k];
        }
        if (pt[k] < bbox.min[k]) {
          bbox.min[k] = pt[k];
        }
      }
    }
    return bbox;
  }

  void convertRGBToYUV()
  {  // BT709
    for (auto& color : colors) {
      const uint8_t r = color[2];
      const uint8_t g = color[0];
      const uint8_t b = color[1];
      const double y = PCCClip(
        std::round(0.212600 * r + 0.715200 * g + 0.072200 * b), 0., 255.);
      const double u = PCCClip(
        std::round(-0.114572 * r - 0.385428 * g + 0.5 * b + 128.0), 0., 255.);
      const double v = PCCClip(
        std::round(0.5 * r - 0.454153 * g - 0.045847 * b + 128.0), 0., 255.);
      color[0] = static_cast<uint8_t>(y);
      color[1] = static_cast<uint8_t>(u);
      color[2] = static_cast<uint8_t>(v);
    }
  }

  void convertYUVToRGB()
  {  // BT709
    for (auto& color : colors) {
      const double y1 = color[0];
      const double u1 = color[1] - 128.0;
      const double v1 = color[2] - 128.0;
      const double r =
        PCCClip(round(y1 /*- 0.00000 * u1*/ + 1.57480 * v1), 0.0, 255.0);
      const double g =
        PCCClip(round(y1 - 0.18733 * u1 - 0.46813 * v1), 0.0, 255.0);
      const double b =
        PCCClip(round(y1 + 1.85563 * u1 /*+ 0.00000 * v1*/), 0.0, 255.0);
      color[2] = static_cast<uint8_t>(r);
      color[0] = static_cast<uint8_t>(g);
      color[1] = static_cast<uint8_t>(b);
    }
  }

private:
  std::vector<Vec3<double>> positions;
  std::vector<Vec3<uint8_t>> colors;
  std::vector<uint16_t> reflectances;
  std::vector<uint8_t> frameidx;
  bool withColors;
  bool withReflectances;
  bool withFrameIndex;
};

//===========================================================================
// Swap the position of two points (including attributes) in the PointSet
// as referenced by the proxies a and b.

inline void
swap(const PCCPointSet3::Proxy& a, const PCCPointSet3::Proxy& b)
{
  a.swap(b);
}

//---------------------------------------------------------------------------
// Swap two point clouds

inline void
swap(PCCPointSet3& a, PCCPointSet3& b)
{
  a.swap(b);
}

//---------------------------------------------------------------------------

} /* namespace pcc */

#endif /* PCCPointSet_h */
