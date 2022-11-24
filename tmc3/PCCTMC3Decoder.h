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

#include "Attribute.h"
#include "PayloadBuffer.h"
#include "PCCMath.h"
#include "PCCPointSet.h"
#include "frame.h"
#include "framectr.h"
#include "geometry.h"
#include "hls.h"

namespace pcc {

//============================================================================

struct DecoderParams {
  // For partial decoding (aka, scalable bitstreams), the number of octree
  // layers to skip during the decode process (attribute coding must take
  // this into account)
  int minGeomNodeSizeLog2;

  // A maximum number of points to partially decode.
  int decodeMaxPoints;

  // Number of fractional bits used in output position representation.
  int outputFpBits;
};

//============================================================================

class PCCTMC3Decoder3 {
public:
  class Callbacks;

  PCCTMC3Decoder3(const DecoderParams& params);
  PCCTMC3Decoder3(const PCCTMC3Decoder3&) = delete;
  PCCTMC3Decoder3(PCCTMC3Decoder3&&) = default;
  PCCTMC3Decoder3& operator=(const PCCTMC3Decoder3& rhs) = delete;
  PCCTMC3Decoder3& operator=(PCCTMC3Decoder3&& rhs) = default;
  ~PCCTMC3Decoder3();

  void init();

  int decompress(const PayloadBuffer* buf, Callbacks* callback);

  //==========================================================================

  void storeSps(SequenceParameterSet&& sps);
  void storeGps(GeometryParameterSet&& gps);
  void storeAps(AttributeParameterSet&& aps);
  void storeTileInventory(TileInventory&& inventory);

  //==========================================================================
  void setMotionVectorFileName(std::string s) { motionVectorFileName = s; }

private:
  void activateParameterSets(const AttributeParamInventoryHdr& gbh);
  void activateParameterSets(const GeometryBrickHeader& gbh);
  int decodeGeometryBrick(const PayloadBuffer& buf);
  void decodeAttributeBrick(const PayloadBuffer& buf);
  void decodeConstantAttribute(const PayloadBuffer& buf);
  bool dectectFrameBoundary(const PayloadBuffer* buf);
  void outputCurrentCloud(Callbacks* callback);
  void storeCurrentCloudAsRef();

  void startFrame();
  void emplaceRefFrame(const SequenceParameterSet& sps);

  //==========================================================================

private:
  // Decoder specific parameters
  DecoderParams _params;

  // Indicates that pointcloud output should be suppressed at a frame boundary
  bool _suppressOutput;

  // Indicates that this is the start of a new frame.
  // NB: this is set to false quiet early in the decoding process
  bool _firstSliceInFrame;

  // Indicates whether the output has been initialised
  bool _outputInitialized;

  // Current identifier of payloads with the same geometry
  int _sliceId;

  // Identifies the previous slice in bistream order
  int _prevSliceId;

  // Cumulative frame counter
  FrameCtr _frameCtr;

  // Position of the slice in the translated+scaled co-ordinate system.
  Vec3<int> _sliceOrigin;

  // The point cloud currently being decoded
  PCCPointSet3 _currentPointCloud;

  // The accumulated decoded slices
  PCCPointSet3 _accumCloud;

  // The current output cloud
  CloudFrame _outCloud;

  // Point positions in spherical coordinates of the current slice
  std::vector<point_t> _posSph;

  // Received parameter sets, mapping parameter set id -> parameterset
  std::map<int, SequenceParameterSet> _spss;
  std::map<int, GeometryParameterSet> _gpss;
  std::map<int, AttributeParameterSet> _apss;

  // mapping sps id to reference Frame buffer
  std::map<int, CloudFrame> _refFrameSeq;

  // Metadata that allows slices/tiles to be indentified by their bounding box
  TileInventory _tileInventory;

  // The active SPS
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;

  // The active reference frame
  const CloudFrame* _refFrame;

  GeometryBrickHeader _gbh;

  // Memorized context buffers
  std::unique_ptr<GeometryOctreeContexts> _ctxtMemOctreeGeom;
  std::unique_ptr<PredGeomContexts> _ctxtMemPredGeom;
  std::vector<AttributeContexts> _ctxtMemAttrs;
  std::vector<int> _ctxtMemAttrSliceIds;

  // Attribute decoder for reuse between attributes of same slice
  std::unique_ptr<AttributeDecoderIntf> _attrDecoder;

  // Point positions in spherical coordinates of the reference frame
  PredGeomPredictor _refFrameSph;
  std::string motionVectorFileName;

  AttributeInterPredParams attrInterPredParams;

  pcc::point_t minPos_ref;
};

//----------------------------------------------------------------------------

class PCCTMC3Decoder3::Callbacks {
public:
  virtual void onOutputCloud(const CloudFrame&) = 0;
};

//============================================================================

}  // namespace pcc
