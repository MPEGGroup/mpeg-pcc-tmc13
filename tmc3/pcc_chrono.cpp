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

#if _WIN32
#  define _UNICODE
#  include <windows.h>
#endif

#include "TMC3Config.h"
#include "pcc_chrono.h"

#if HAVE_GETRUSAGE
#  include <sys/time.h>
#  include <sys/resource.h>
#endif

//===========================================================================

#if _WIN32
namespace pcc {
namespace chrono {
  namespace detail {
    using hundredns = std::chrono::duration<int64_t, std::ratio<1, 10000000>>;

    // global state to emulate getrusage(RUSAGE_CHILDREN).
    hundredns g_cumulative_time_children{0};
  }  // namespace detail
}  // namespace chrono
}  // namespace pcc
#endif

//---------------------------------------------------------------------------

#if _WIN32
pcc::chrono::utime_inc_children_clock::time_point
pcc::chrono::utime_inc_children_clock::now() noexcept
{
  HANDLE hProcess = GetCurrentProcess();
  FILETIME dummy, userTime;

  GetProcessTimes(hProcess, &dummy, &dummy, &dummy, &userTime);

  ULARGE_INTEGER val;
  val.LowPart = userTime.dwLowDateTime;
  val.HighPart = userTime.dwHighDateTime;

  using namespace detail;
  using hundredns = std::chrono::duration<int64_t, std::ratio<1, 10000000>>;
  return time_point(hundredns(val.QuadPart) + g_cumulative_time_children);
}
#endif

//---------------------------------------------------------------------------

#if HAVE_GETRUSAGE
pcc::chrono::utime_inc_children_clock::time_point
pcc::chrono::utime_inc_children_clock::now() noexcept
{
  std::chrono::nanoseconds total;

  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  total = std::chrono::seconds(usage.ru_utime.tv_sec)
    + std::chrono::microseconds(usage.ru_utime.tv_usec);

  if (getrusage(RUSAGE_CHILDREN, &usage))
    return time_point(total);

  total += std::chrono::seconds(usage.ru_utime.tv_sec)
    + std::chrono::microseconds(usage.ru_utime.tv_usec);

  return time_point(total);
}
#endif

//===========================================================================
