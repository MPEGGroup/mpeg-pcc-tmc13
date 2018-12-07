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

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "PayloadBuffer.h"
#include "PCCMath.h"
#include "PCCPointSet.h"
#include "hls.h"

namespace pcc {

//============================================================================

enum class PartitionMethod
{
  // Don't partition input
  kNone = 0,

  // Partition according to uniform geometry
  kUniformGeom = 2,

  // Partition according to the depth of octree
  kOctreeUniform = 3,
};

//============================================================================

struct EncoderParams {
  SequenceParameterSet sps;
  GeometryParameterSet gps;

  // NB: information about attributes is split between the SPS and the APS.
  //  => The SPS enumerates the attributes, the APS controls coding params.
  std::vector<AttributeParameterSet> aps;

  // todo(df): this should go away
  std::map<std::string, int> attributeIdxMap;

  // Filename for saving recoloured point cloud.
  std::string postRecolorPath;

  // Method for partitioning the input cloud
  PartitionMethod partitionMethod;

  // Number of slices used by PartitionMethod::kUniformGeom
  int partitionNumUniformGeom;

  // Depth of octree used by PartitionMethod::kOctreeUniform
  int partitionOctreeDepth;
};

//============================================================================

class PCCTMC3Encoder3 {
public:
  PCCTMC3Encoder3() = default;
  PCCTMC3Encoder3(const PCCTMC3Encoder3&) = default;
  PCCTMC3Encoder3& operator=(const PCCTMC3Encoder3& rhs) = default;
  ~PCCTMC3Encoder3() = default;

  int compress(
    const PCCPointSet3& inputPointCloud,
    EncoderParams* params,
    std::function<void(const PayloadBuffer&)> outputFn,
    PCCPointSet3* reconstructedCloud = nullptr);

  void compressPartition(
    const PCCPointSet3& inputPointCloud,
    EncoderParams* params,
    std::function<void(const PayloadBuffer&)> outputFn,
    PCCPointSet3* reconstructedCloud = nullptr);

  static void fixupParameterSets(EncoderParams* params);

private:
  void appendReconstructedPoints(PCCPointSet3* reconstructedCloud);

  void encodeGeometryBrick(PayloadBuffer* buf);

  void quantization(const PCCPointSet3& inputPointCloud);

private:
  PCCPointSet3 pointCloud;

  // Position of the slice in the translated+scaled co-ordinate system.
  PCCVector3<int> _sliceOrigin;

  // Size of the current slice
  PCCVector3<int> _sliceBoxWhd;

  // The active parameter sets
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;
  std::vector<const AttributeParameterSet*> _aps;

  // Current identifier of payloads with the same geometry
  int _sliceId;

  // Identifies the current tile
  int _tileId;
};

//============================================================================

}  // namespace pcc
