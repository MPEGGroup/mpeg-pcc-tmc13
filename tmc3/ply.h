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
#include <string>

#include "PCCPointSet.h"

namespace pcc {
namespace ply {

  //============================================================================
  // This defines the the attribute names to be used when reading/writing ply
  // files.  The order of the names define the ordering used in the internal
  // point cloud representation.

  struct PropertyNameMap {
    // The names of the position attributes, typically {"x", "y", "z"}
    std::array<const char*, 3> position;
  };

  //============================================================================

  ///
  // Write @a pointCloud to a PLY file called @a fileName.
  // Each point position, pt, is converted prior to writing by:
  //  pt' = pt * positionScale + offset
  //
  // @param pointCloud  points to be written.
  // @param propertyNames  Describes ply property names of pointcloud attributes.
  // @param positionScale  scale factor for positions.
  // @param positionOffset  offset for positions (after scaling).
  // @param fileName  output filename.
  // @param asAscii  PLY writing format (true => ascii, false => binary).
  bool write(
    const PCCPointSet3& pointCloud,
    const PropertyNameMap& propertyNames,
    double positionScale,
    Vec3<double> positionOffset,
    const std::string& fileName,
    bool asAscii);

  ///
  // Read @a pointCloud to a PLY file called @a fileName.
  // Point positions are scaled by positionScale and converted to integers.
  //
  bool read(
    const std::string& fileName,
    const PropertyNameMap& propertyNames,
    double positionScale,
    PCCPointSet3& cloud);

  //============================================================================

}  // namespace ply
}  // namespace pcc