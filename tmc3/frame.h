/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2021, ISO/IEC
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

#include <vector>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "hls.h"

namespace pcc {

//============================================================================
// Represents a frame in the encoder or decoder.

struct CloudFrame {
  // The value of the decoder's reconstructed FrameCtr for this frame.
  int frameNum;

  // Defines the ordering of the position components (eg, xyz vs zyx)
  AxisOrder geometry_axis_order;

  // The length of the output cloud's unit vector.
  // The units of outputUnitLength is given by outputUnit.
  //
  // When outputUnit is ScaleUnit::kMetres, outputUnitLength is
  // measured in metres.
  //
  // When outputUnit is ScaleUnit::kDimensionless, outputUnitLength is
  // measured in units of an external coordinate system.
  double outputUnitLength;

  // The unit of the output cloud's unit vector
  ScaleUnit outputUnit;

  // The origin of the output cloud
  // NB: this respects geometry_axis_order.
  Vec3<int> outputOrigin;

  // Number of fractional bits in representaiton of cloud.positions.
  int outputFpBits;

  // Descriptions of each attribute.  Attribute parameters in the description
  // are only applicable to this frame -- they may change in a subsequent frame
  std::vector<AttributeDescription> attrDesc;

  // The output point cloud.  The coordinate system is defined by the
  // other parameters in this structure.
  // NB: Point positions respect geometry_axis_order.
  PCCPointSet3 cloud;

  // Determines parameters according to the sps.
  void setParametersFrom(const SequenceParameterSet& sps, int fixedPointBits);
};

//============================================================================
// Scale point cloud geometry by global scale factor

void scaleGeometry(
  PCCPointSet3& cloud,
  const SequenceParameterSet::GlobalScale& globalScale,
  int fixedPointFracBits);

//============================================================================

}  // namespace pcc
