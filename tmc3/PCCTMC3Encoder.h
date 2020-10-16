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

#include "Attribute.h"
#include "PayloadBuffer.h"
#include "PCCMath.h"
#include "PCCPointSet.h"
#include "geometry_params.h"
#include "pointset_processing.h"
#include "hls.h"
#include "partitioning.h"
#include "geometry.h"

namespace pcc {

//============================================================================

struct EncoderAttributeParams {
  // NB: this only makes sense for setting configurable parameters
  AttributeBrickHeader abh;

  // Threshold for choosing dist2 out of the population of nearest neighbour
  // distances.
  float dist2PercentileEstimate;
};

//----------------------------------------------------------------------------

struct EncoderParams {
  SequenceParameterSet sps;
  GeometryParameterSet gps;
  GeometryBrickHeader gbh;

  // NB: information about attributes is split between the SPS and the APS.
  //  => The SPS enumerates the attributes, the APS controls coding params.
  std::vector<AttributeParameterSet> aps;

  // Encoder specific parameters for attributes
  std::vector<EncoderAttributeParams> attr;

  // todo(df): this should go away
  std::map<std::string, int> attributeIdxMap;

  // Resolution of the source point cloud points per metres.
  // The value 0 indicates undefined.
  float srcResolution;

  // Quantisation scale factor for geometry pre-processing
  float geomPreScale;

  // Encoder specific parameters for geometry
  OctreeEncOpts geom;

  // Options for the predictive geometry coder
  PredGeomEncOpts predGeom;

  // Parameters that control partitioning
  PartitionParams partition;

  // attribute recolouring parameters
  RecolourParams recolour;

  // LiDAR head position
  Vec3<int> lidarHeadPosition;

  // number of expected lasers
  int numLasers;

  // floating Lasers' theta (have to be converted to fixed point in gps)
  std::vector<double> lasersTheta;

  // floating Lasers' H (have to be converted to fixed point in gps)
  std::vector<double> lasersZ;

  // per-slice trisoup node sizes
  std::vector<int> trisoupNodeSizesLog2;

  // Enable enforcement of level limits (encoder will abort if exceeded)
  bool enforceLevelLimits;

  // Qp used for IDCM quantisation (used to derive HLS values)
  int idcmQp;
};

//============================================================================

class PCCTMC3Encoder3 {
public:
  class Callbacks;

  PCCTMC3Encoder3();
  PCCTMC3Encoder3(const PCCTMC3Encoder3&) = delete;
  PCCTMC3Encoder3(PCCTMC3Encoder3&&) = default;
  PCCTMC3Encoder3& operator=(const PCCTMC3Encoder3& rhs) = delete;
  PCCTMC3Encoder3& operator=(PCCTMC3Encoder3&& rhs) = default;
  ~PCCTMC3Encoder3();

  int compress(
    const PCCPointSet3& inputPointCloud,
    EncoderParams* params,
    Callbacks*,
    PCCPointSet3* reconstructedCloud = nullptr);

  void compressPartition(
    const PCCPointSet3& inputPointCloud,
    const PCCPointSet3& originPartCloud,
    EncoderParams* params,
    Callbacks*,
    PCCPointSet3* reconstructedCloud = nullptr);

  static void deriveParameterSets(EncoderParams* params);
  static void fixupParameterSets(EncoderParams* params);

private:
  void appendReconstructedPoints(PCCPointSet3* reconstructedCloud);

  void encodeGeometryBrick(const EncoderParams*, PayloadBuffer* buf);

  SrcMappedPointSet quantization(const PCCPointSet3& src);

private:
  PCCPointSet3 pointCloud;

  // Point positions in spherical coordinates of the current slice
  std::vector<point_t> _posSph;

  // Scale factor used to quantise geometry during pre-processing
  float _geomPreScale;

  // Position of the slice in the translated+scaled co-ordinate system.
  Vec3<int> _sliceOrigin;

  // Size of the current slice
  Vec3<int> _sliceBoxWhd;

  // The active parameter sets
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;
  std::vector<const AttributeParameterSet*> _aps;

  // Cached copy of the curent _gbh (after encoding geometry)
  GeometryBrickHeader _gbh;

  // Indicates that this is the start of a new frame
  bool _firstSliceInFrame;

  // Current identifier of payloads with the same geometry
  int _sliceId;

  // Identifies the previous slice in bistream order
  int _prevSliceId;

  // Identifies the current tile
  int _tileId;

  // Current frame number.
  // NB: only the log2_max_frame_idx LSBs are sampled for frame_idx
  int _frameCounter;

  // Memorized context buffers
  std::unique_ptr<GeometryOctreeContexts> _ctxtMemOctreeGeom;
  std::unique_ptr<PredGeomContexts> _ctxtMemPredGeom;
  std::vector<AttributeContexts> _ctxtMemAttrs;
  std::vector<int> _ctxtMemAttrSliceIds;
};

//----------------------------------------------------------------------------

class PCCTMC3Encoder3::Callbacks {
public:
  virtual void onOutputBuffer(const PayloadBuffer&) = 0;
  virtual void onPostRecolour(const PCCPointSet3&) = 0;
};

//============================================================================

}  // namespace pcc
