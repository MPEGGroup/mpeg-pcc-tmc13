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

#include <chrono>

//===========================================================================

namespace pcc {
namespace chrono {
  /**
 * Clock reporting elapsed user time of the current process and children.
 *
 * NB: under winapi, only child processes that have completed execution
 *     via pcc::system() are reported.
 */
  struct utime_inc_children_clock {
    typedef std::chrono::nanoseconds duration;
    typedef duration::rep rep;
    typedef duration::period period;
    typedef std::chrono::time_point<utime_inc_children_clock, duration>
      time_point;

    static constexpr bool is_steady = true;

    static time_point now() noexcept;
  };
}  // namespace chrono
}  // namespace pcc

//===========================================================================

namespace pcc {
namespace chrono {
  /**
 * Measurement of cumulative elapsed time intervals.
 *
 * This timer acts like a stopwatch and may be used to measure the
 * cumulative periods between successive calls to start() and stop().
 */
  template<typename Clock>
  class Stopwatch {
  public:
    typedef typename Clock::duration duration;

    /// Reset the accumulated interval count.
    void reset();

    /// Mark the beginning of a measurement period.
    void start();

    /// Mark the end of a measurement period.
    ///
    /// @return  the duration of the elapsed Clock time since start()
    duration stop();

    /// The sum of the previous elapsed time intervals.
    ///
    /// NB: this excludes any currently active period marked by start().
    ///
    /// @return  cumulative elapsed time
    constexpr duration count() const { return cumulative_time_; }

  private:
    typename Clock::time_point start_time_;
    duration cumulative_time_{duration::zero()};
  };
}  // namespace chrono
}  // namespace pcc

//---------------------------------------------------------------------------

template<typename Clock>
void
pcc::chrono::Stopwatch<Clock>::reset()
{
  cumulative_time_ = cumulative_time_.zero();
}

//---------------------------------------------------------------------------

template<typename Clock>
void
pcc::chrono::Stopwatch<Clock>::start()
{
  start_time_ = Clock::now();
}

//---------------------------------------------------------------------------

template<typename Clock>
typename pcc::chrono::Stopwatch<Clock>::duration
pcc::chrono::Stopwatch<Clock>::stop()
{
  const auto& delta = duration(Clock::now() - start_time_);
  cumulative_time_ += delta;
  return delta;
}

//===========================================================================
