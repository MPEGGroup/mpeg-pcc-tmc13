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
#include "frame.h"
#include "geometry.h"
#include "geometry_params.h"
#include "hls.h"
#include "partitioning.h"
#include "pointset_processing.h"

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

  // Determine the sequence bonuding box using the first input frame
  bool autoSeqBbox;

  // Length of the source point cloud unit vectors.
  double srcUnitLength;

  // Scale factor used to define the coordinate system used for coding.
  // This is the coordinate system where slicing is performed.
  //  P_cod = P_src * codedGeomScale
  double codedGeomScale;

  // Scale factor used to define the sequence coordinate system.
  //  P_seq = P_src * seqGeomScale
  double seqGeomScale;

  // Scale factor used to define the external coordinate system.
  //  P_ext = P_src * extGeomScale
  double extGeomScale;

  // Number of fractional bits used in output position representation.
  int outputFpBits;

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

  // precision expected for attributes after scaling with predgeom
  // and spherical coordinates
  int attrSphericalMaxLog2;
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
    CloudFrame* reconstructedCloud = nullptr);

  void compressPartition(
    const PCCPointSet3& inputPointCloud,
    const PCCPointSet3& originPartCloud,
    EncoderParams* params,
    Callbacks*,
    CloudFrame* reconstructedCloud = nullptr);

  static void deriveParameterSets(EncoderParams* params);
  static void fixupParameterSets(EncoderParams* params);

private:
  void appendSlice(PCCPointSet3& cloud);

  void encodeGeometryBrick(const EncoderParams*, PayloadBuffer* buf);

  SrcMappedPointSet quantization(const PCCPointSet3& src);

private:
  PCCPointSet3 pointCloud;

  // Point positions in spherical coordinates of the current slice
  std::vector<point_t> _posSph;

  // Scale factor used to decimate the input point cloud.
  // Decimation is performed as if the input were scaled by
  //   Round(P_src * inputDecimationScale)
  // and duplicate points removed.
  // todo: expose this parameter?
  double _inputDecimationScale;

  // Scale factor that defines coding coordinate system
  double _srcToCodingScale;

  // Sequence origin in terms of coding coordinate system
  Vec3<int> _originInCodingCoords;

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
  // NB: only the log2_max_frame_ctr LSBs are sampled for frame_ctr
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
