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

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <type_traits>

namespace pcc {

//===========================================================================
// A ringbuf::iterator with RandomAccessIterator semantics.

template<typename T>
class ring_iterator {
  typedef typename std::remove_const<T>::type T_nonconst;

  typedef ring_iterator<T> iterator;
  typedef ring_iterator<const T_nonconst> const_iterator;

  // the backing store
  T* base_;

  // length of backing store
  size_t max_;

  // iterator's index into backing store
  size_t idx_;

  // notional start of ring buffer, used to determine index ordering
  const iterator* start_;

public:
  typedef std::random_access_iterator_tag iterator_category;
  typedef T value_type;
  typedef std::ptrdiff_t difference_type;
  typedef T* pointer;
  typedef T& reference;

  //-------------------------------------------------------------------------

  explicit ring_iterator(T* base, size_t max, const iterator* start)
    : base_(base), max_(max), idx_(), start_(start)
  {}

  operator ring_iterator<const T_nonconst>() const
  {
    return *(const_iterator*)this;
  }

  //-------------------------------------------------------------------------
  // Iterator

  reference operator*() const { return base_[idx_]; }

  iterator& operator++()
  {
    idx_++;
    idx_ = idx_ == max_ ? 0 : idx_;
    return *this;
  }

  //-------------------------------------------------------------------------
  // InputIterator

  bool operator==(const iterator& other) const { return idx_ == other.idx_; }
  bool operator!=(const iterator& other) const { return !(*this == other); }

  pointer operator->() const { return &base_[idx_]; }

  iterator operator++(int)
  {
    iterator retval = *this;
    ++(*this);
    return retval;
  }

  //-------------------------------------------------------------------------
  // BidirectionalIterator

  iterator& operator--()
  {
    idx_ = idx_ == 0 ? max_ : idx_;
    idx_--;
    return *this;
  }

  iterator operator--(int)
  {
    iterator retval = *this;
    --(*this);
    return retval;
  }

  //-------------------------------------------------------------------------
  // RandomAccessIterator

  reference operator[](size_t n)
  {
    size_t idx = idx_ + n;
    if (idx_ >= max_)
      idx_ -= max_;
    return base_[idx];
  }

  // Precondition: 0 < it + n < end;
  // ie, the calculated position is not modulo size().
  iterator& operator+=(difference_type n)
  {
    idx_ += n;
    if (idx_ >= max_) {
      if (n > 0)
        idx_ -= max_;
      else
        idx_ += max_;
    }
    return *this;
  }

  // Precondition: 0 < it + n < end;
  // ie, the calculated position is not modulo size().
  iterator operator+(difference_type n) const
  {
    iterator it(*this);
    it += n;
    return it;
  }

  // Precondition: 0 < it + n < end;
  // ie, the calculated position is not modulo size().
  iterator& operator-=(difference_type n)
  {
    *this += -n;
    return *this;
  }

  // Precondition: 0 < it + n < end;
  // ie, the calculated position is not modulo size().
  iterator operator-(difference_type n) const
  {
    iterator it(*this);
    it -= n;
    return it;
  }

  difference_type operator-(const iterator& other) const
  {
    size_t lin_pos_this = linear_idx(*start_, idx_);
    size_t lin_pos_other = linear_idx(*start_, other.idx_);
    return lin_pos_this - lin_pos_other;
  }

  bool operator<(const iterator& other) const
  {
    difference_type dist = operator-(other);
    return dist < 0;
  }

  bool operator<=(const iterator& other) const
  {
    difference_type dist = operator-(other);
    return dist <= 0;
  }

  bool operator>(const iterator& other) const { return !operator<=(other); }

  bool operator>=(const iterator& other) const { return !operator<(other); }

  //=========================================================================
private:
  // calculate the linear position of idx_ with respect to origin
  static size_t linear_idx(const iterator& origin, size_t idx_)
  {
    difference_type dist = difference_type(idx_ - origin.idx_);
    if (dist < 0)
      dist += origin.max_;
    return size_t(dist);
  }
};

//===========================================================================
// An indexed sequence cotainer allowing fast insertion and deletion at both
// the head and tail.
//
// As opposed to std::deque, a pcc::ringbuf uses a fixed-sized backing store
// that is not resized so as to avoid memory management calls when inserting
// or erasing elements.  To permit fast manipulation of the head and tail,
// elements are not guaranteed to be stored contigiously.
//
// ## The complexity (efficiency) of common operations is as follows:
//
//  - Random access -- constant O(1).
//  - Insertion / removal at the head or tail -- constant O(1).
//  - Insertion / removal of elements -- linear O(n).
//
// ## Incomplete implementation:
//
// The following methods are not implemented:
//
//  - push_front, emplace_front, insert, emplace, at, swap.
//
// ## Iterator validity:
//
//  - push_back, emplace_back do not invalidate any references to elements,
//    all iterators remain valid.
//  - pop_front, pop_back do not invalidate any references to non-erased
//    elements, and all iterators remain valid.

template<typename T>
class ringbuf {
public:
  typedef T value_type;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef size_t size_type;
  typedef std::ptrdiff_t difference_type;

  typedef ring_iterator<T> iterator;
  typedef ring_iterator<const T> const_iterator;

  //--------------------------------------------------------------------------

  ringbuf()
    : buf_(nullptr)
    , capacity_(0)
    , rd_it_(iterator(buf_.get(), capacity_, &rd_it_))
    , wr_it_(iterator(buf_.get(), capacity_, &rd_it_))
  {}

  //--------------------------------------------------------------------------

  ringbuf(size_t size)
    : buf_(static_cast<T*>(operator new[](sizeof(T) * (size + 1))))
    , capacity_(size + 1)
    , rd_it_(iterator(buf_.get(), capacity_, &rd_it_))
    , wr_it_(iterator(buf_.get(), capacity_, &rd_it_))
  {}

  //--------------------------------------------------------------------------

  ringbuf(ringbuf&& other) noexcept { *this = std::move(other); }

  //--------------------------------------------------------------------------

  ringbuf& operator=(ringbuf&& rhs) noexcept
  {
    // unfortunately the iterators cannot be moved, since they hold a
    // pointer to the ringbuffer's internal iterator.
    // it_zero is used to extract the true iterator idx
    iterator it_zero = iterator(nullptr, 0, &it_zero);

    auto lhs_rd_idx = -(it_zero - rd_it_);
    auto lhs_wr_idx = -(it_zero - wr_it_);

    auto rhs_rd_idx = -(it_zero - rhs.rd_it_);
    auto rhs_wr_idx = -(it_zero - rhs.wr_it_);

    std::swap(buf_, rhs.buf_);
    std::swap(capacity_, rhs.capacity_);

    rd_it_ = iterator(buf_.get(), capacity_, &rd_it_);
    wr_it_ = iterator(buf_.get(), capacity_, &rd_it_);
    rd_it_ += rhs_rd_idx;
    wr_it_ += rhs_wr_idx;

    rhs.rd_it_ = iterator(rhs.buf_.get(), rhs.capacity_, &rhs.rd_it_);
    rhs.wr_it_ = iterator(rhs.buf_.get(), rhs.capacity_, &rhs.rd_it_);
    rhs.rd_it_ += lhs_rd_idx;
    rhs.wr_it_ += lhs_wr_idx;

    return *this;
  }

  //--------------------------------------------------------------------------

  ~ringbuf()
  {
    while (rd_it_ != wr_it_) {
      pop_front();
    }
  }

  //--------------------------------------------------------------------------

  iterator begin() { return rd_it_; }
  iterator end() { return wr_it_; }

  const_iterator begin() const { return const_iterator(rd_it_); }
  const_iterator end() const { return const_iterator(wr_it_); }

  //--------------------------------------------------------------------------

  bool empty() const { return begin() == end(); }

  //--------------------------------------------------------------------------

  void push_back(const value_type& val) { emplace_back(val); }

  //--------------------------------------------------------------------------

  void push_back(value_type&& val) { emplace_back(std::move(val)); }

  //--------------------------------------------------------------------------

  template<class... Args>
  void emplace_back(Args&&... args)
  {
    new (&*wr_it_) T(args...);
    ++wr_it_;
  }

  //--------------------------------------------------------------------------

  void pop_back()
  {
    --wr_it_;
    wr_it_->~T();
  }

  //--------------------------------------------------------------------------

  void pop_front()
  {
    rd_it_->~T();
    ++rd_it_;
  }

  //--------------------------------------------------------------------------

  reference front() { return *rd_it_; }

  //--------------------------------------------------------------------------

  reference back() { return *std::prev(wr_it_); }

  //--------------------------------------------------------------------------

  reference operator[](size_type idx) { return *std::next(rd_it_, idx); }

  //--------------------------------------------------------------------------

  const_reference operator[](size_type idx) const
  {
    return *std::next(rd_it_, idx);
  }

  //--------------------------------------------------------------------------

  void clear()
  {
    while (!empty()) {
      pop_front();
    }
  }

  //--------------------------------------------------------------------------

  size_t size() const { return size_t(end() - begin()); }

  //--------------------------------------------------------------------------

  size_t capacity() const { return capacity_ - 1; }

  //--------------------------------------------------------------------------
private:
  struct operator_delete_arr {
    void operator()(T* ptr) { ::operator delete[](ptr); }
  };

  std::unique_ptr<T[], operator_delete_arr> buf_;
  size_t capacity_;
  iterator rd_it_;
  iterator wr_it_;
};

//===========================================================================

} /* namespace pcc */
