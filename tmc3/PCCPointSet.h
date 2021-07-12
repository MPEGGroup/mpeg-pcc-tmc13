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

//============================================================================
// The type used for internally representing attribute data
typedef uint16_t attr_t;

// The type used for internally representing positions
typedef Vec3<int32_t> point_t;

//============================================================================

class PCCPointSet3 {
public:
  typedef point_t PointType;

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
    withLaserAngles = false;
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
    swap(laserAngles, other.laserAngles);
    swap(withColors, other.withColors);
    swap(withReflectances, other.withReflectances);
    swap(withFrameIndex, other.withFrameIndex);
    swap(withLaserAngles, other.withLaserAngles);
  }

  PointType operator[](const size_t index) const
  {
    assert(index < positions.size());
    return positions[index];
  }
  PointType& operator[](const size_t index)
  {
    assert(index < positions.size());
    return positions[index];
  }

  void swapPoints(std::vector<PointType>& other) { positions.swap(other); }

  Vec3<attr_t> getColor(const size_t index) const
  {
    assert(index < colors.size() && withColors);
    return colors[index];
  }
  Vec3<attr_t>& getColor(const size_t index)
  {
    assert(index < colors.size() && withColors);
    return colors[index];
  }
  void setColor(const size_t index, const Vec3<attr_t> color)
  {
    assert(index < colors.size() && withColors);
    colors[index] = color;
  }
  attr_t getReflectance(const size_t index) const
  {
    assert(index < reflectances.size() && withReflectances);
    return reflectances[index];
  }
  attr_t& getReflectance(const size_t index)
  {
    assert(index < reflectances.size() && withReflectances);
    return reflectances[index];
  }
  void setReflectance(const size_t index, const attr_t reflectance)
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

  int getLaserAngle(const size_t index) const
  {
    assert(index < laserAngles.size() && withLaserAngles);
    return laserAngles[index];
  }

  int& getLaserAngle(const size_t index)
  {
    assert(index < laserAngles.size() && withLaserAngles);
    return laserAngles[index];
  }

  void setLaserAngle(const size_t index, const int angle)
  {
    assert(index < laserAngles.size() && withLaserAngles);
    laserAngles[index] = angle;
  }

  bool hasLaserAngles() const { return withLaserAngles; }
  void addLaserAngles()
  {
    withLaserAngles = true;
    resize(getPointCount());
  }
  void removeLaserAngles()
  {
    withLaserAngles = false;
    laserAngles.resize(0);
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

  void addRemoveAttributes(const PCCPointSet3& ref)
  {
    ref.hasColors() ? addColors() : removeColors();
    ref.hasReflectances() ? addReflectances() : removeReflectances();
    ref.hasLaserAngles() ? addLaserAngles() : removeLaserAngles();
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
    if (hasLaserAngles()) {
      laserAngles.resize(size);
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
    if (hasLaserAngles()) {
      laserAngles.reserve(size);
    }
  }
  void clear()
  {
    positions.clear();
    colors.clear();
    reflectances.clear();
    frameidx.clear();
    laserAngles.clear();
  }

  size_t removeDuplicatePointInQuantizedPoint(int minGeomNodeSizeLog2)
  {
    for (int i = 0; i < positions.size(); i++) {
      PointType newPoint = positions[i];
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
      addRemoveAttributes(src);

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

    if (hasLaserAngles())
      std::copy(
        src.laserAngles.begin(), src.laserAngles.end(),
        std::next(laserAngles.begin(), dstEnd));
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
    if (hasLaserAngles()) {
      std::swap(getLaserAngle(index1), getLaserAngle(index2));
    }
  }

  Box3<int32_t> computeBoundingBox() const
  {
    Box3<int32_t> bbox(
      std::numeric_limits<int32_t>::max(),
      std::numeric_limits<int32_t>::lowest());
    const size_t pointCount = getPointCount();
    for (size_t i = 0; i < pointCount; ++i) {
      const auto& pt = (*this)[i];
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
  Box3<int32_t> computeBoundingBox(ForwardIt begin, ForwardIt end) const
  {
    Box3<int32_t> bbox(
      std::numeric_limits<int32_t>::max(),
      std::numeric_limits<int32_t>::lowest());

    for (auto it = begin; it != end; ++it) {
      int i = *it;
      const auto& pt = (*this)[i];
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

private:
  std::vector<PointType> positions;
  std::vector<Vec3<attr_t>> colors;
  std::vector<attr_t> reflectances;
  std::vector<uint8_t> frameidx;
  bool withColors;
  bool withReflectances;
  bool withFrameIndex;
  std::vector<int> laserAngles;
  bool withLaserAngles;
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
