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

#include "PCCTMC3Decoder.h"

#include <algorithm>
#include <cassert>
#include <string>

#include "AttributeCommon.h"
#include "PayloadBuffer.h"
#include "PCCPointSet.h"
#include "coordinate_conversion.h"
#include "geometry.h"
#include "geometry_octree.h"
#include "geometry_predictive.h"
#include "hls.h"
#include "io_hls.h"
#include "io_tlv.h"
#include "pcc_chrono.h"
#include "osspecific.h"

namespace pcc {

//============================================================================

PCCTMC3Decoder3::PCCTMC3Decoder3(const DecoderParams& params) : _params(params)
{
  init();
}

//----------------------------------------------------------------------------

void
PCCTMC3Decoder3::init()
{
  _firstSliceInFrame = true;
  _outputInitialized = false;
  _suppressOutput = 1;
  _sps = nullptr;
  _gps = nullptr;
  _spss.clear();
  _gpss.clear();
  _apss.clear();

  _ctxtMemOctreeGeom.reset(new GeometryOctreeContexts);
  _ctxtMemPredGeom.reset(new PredGeomContexts);
}

//----------------------------------------------------------------------------

PCCTMC3Decoder3::~PCCTMC3Decoder3() = default;

//============================================================================

static bool
payloadStartsNewSlice(PayloadType type)
{
  return type == PayloadType::kGeometryBrick
    || type == PayloadType::kFrameBoundaryMarker;
}

//============================================================================

bool
PCCTMC3Decoder3::dectectFrameBoundary(const PayloadBuffer* buf)
{
  // This may be from either geometry brick or attr param inventory
  int frameCtrLsb;

  switch (buf->type) {
  case PayloadType::kFrameBoundaryMarker: {
    // the frame boundary data marker explcitly indicates a boundary
    // However, this implementation doesn't flush the output, rather
    // this happens naturally when the frame boundary is detected by
    // a change in frameCtr.
    auto fbm = parseFrameBoundaryMarker(*buf);
    frameCtrLsb = fbm.fbdu_frame_ctr_lsb;
    break;
  }

  case PayloadType::kGeometryBrick: {
    activateParameterSets(parseGbhIds(*buf));
    auto gbh = parseGbh(*_sps, *_gps, *buf, nullptr, nullptr);
    frameCtrLsb = gbh.frame_ctr_lsb;
    break;
  }

  case PayloadType::kGeneralizedAttrParamInventory: {
    auto apih = parseAttrParamInventoryHdr(*buf);
    activateParameterSets(apih);
    // todo(conf): check lsb_bits is same as sps
    frameCtrLsb = apih.attr_param_frame_ctr_lsb;
    break;
  }

  // other data units don't indicate a boundary
  default: return false;
  }

  auto bdry = _frameCtr.isDifferentFrame(frameCtrLsb, _sps->frame_ctr_bits);
  _frameCtr.update(frameCtrLsb, _sps->frame_ctr_bits);

  return bdry;
}

//============================================================================

void
PCCTMC3Decoder3::outputCurrentCloud(PCCTMC3Decoder3::Callbacks* callback)
{
  if (_suppressOutput)
    return;

  std::swap(_outCloud.cloud, _accumCloud);

  // Apply global scaling to output for integer conformance
  // todo: add other output scaling modes
  // NB: if accumCloud is reused for future inter-prediction, global scaling
  //     must be applied to a copy.
  scaleGeometry(_outCloud.cloud, _sps->globalScale, _outCloud.outputFpBits);

  callback->onOutputCloud(_outCloud);

  std::swap(_outCloud.cloud, _accumCloud);
  _accumCloud.clear();
}

//============================================================================

void
PCCTMC3Decoder3::startFrame()
{
  _outputInitialized = true;
  _firstSliceInFrame = true;
  _outCloud.frameNum = _frameCtr;

  // the following could be set once when the SPS is discovered
  _outCloud.setParametersFrom(*_sps, _params.outputFpBits);
}

//============================================================================

int
PCCTMC3Decoder3::decompress(
  const PayloadBuffer* buf, PCCTMC3Decoder3::Callbacks* callback)
{
  // Starting a new geometry brick/slice/tile, transfer any
  // finished points to the output accumulator
  if (!buf || payloadStartsNewSlice(buf->type)) {
    if (size_t numPoints = _currentPointCloud.getPointCount()) {
      for (size_t i = 0; i < numPoints; i++)
        for (int k = 0; k < 3; k++)
          _currentPointCloud[i][k] += _sliceOrigin[k];
      _accumCloud.append(_currentPointCloud);
    }
  }

  if (!buf) {
    // flush decoder, output pending cloud if any
    outputCurrentCloud(callback);
    return 0;
  }

  // process a frame boundary
  //  - this may update FrameCtr
  //  - this will activate the sps for GeometryBrick and AttrParamInventory
  //  - after outputing the current frame, the output must be reinitialized
  if (dectectFrameBoundary(buf)) {
    outputCurrentCloud(callback);
    _outputInitialized = false;
  }

  // process the buffer
  switch (buf->type) {
  case PayloadType::kSequenceParameterSet: {
    auto sps = parseSps(*buf);
    convertXyzToStv(&sps);
    storeSps(std::move(sps));
    return 0;
  }

  case PayloadType::kGeometryParameterSet: {
    auto gps = parseGps(*buf);
    // HACK: assume that an SPS has been received prior to the GPS.
    // This is not required, and parsing of the GPS is independent of the SPS.
    // todo(df): move GPS fixup to activation process
    _sps = &_spss.cbegin()->second;
    convertXyzToStv(*_sps, &gps);
    storeGps(std::move(gps));
    return 0;
  }

  case PayloadType::kAttributeParameterSet: {
    auto aps = parseAps(*buf);
    // HACK: assume that an SPS has been received prior to the APS.
    // This is not required, and parsing of the APS is independent of the SPS.
    // todo(df): move APS fixup to activation process
    _sps = &_spss.cbegin()->second;
    convertXyzToStv(*_sps, &aps);
    storeAps(std::move(aps));
    return 0;
  }

  case PayloadType::kFrameBoundaryMarker:
    if (!_outputInitialized)
      startFrame();
    return 0;

  case PayloadType::kGeometryBrick:
    if (!_outputInitialized)
      startFrame();

    // avoid accidents with stale attribute decoder on next slice
    _attrDecoder.reset();
    // Avoid dropping an actual frame
    _suppressOutput = false;
    return decodeGeometryBrick(*buf);

  case PayloadType::kAttributeBrick: decodeAttributeBrick(*buf); return 0;

  case PayloadType::kConstantAttribute:
    decodeConstantAttribute(*buf);
    return 0;

  case PayloadType::kTileInventory:
    // NB: the tile inventory is decoded in xyz order.  It may need
    //     conversion if it is used (it currently isn't).
    storeTileInventory(parseTileInventory(*buf));
    return 0;

  case PayloadType::kGeneralizedAttrParamInventory: {
    if (!_outputInitialized)
      startFrame();

    auto hdr = parseAttrParamInventoryHdr(*buf);
    assert(hdr.attr_param_sps_attr_idx < int(_sps->attributeSets.size()));
    auto& attrDesc = _outCloud.attrDesc[hdr.attr_param_sps_attr_idx];
    parseAttrParamInventory(attrDesc, *buf, attrDesc.params);
    return 0;
  }

  case PayloadType::kUserData: parseUserData(*buf); return 0;
  }

  // todo(df): error, unhandled payload type
  return 1;
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::storeSps(SequenceParameterSet&& sps)
{
  // todo(df): handle replacement semantics
  _spss.emplace(std::make_pair(sps.sps_seq_parameter_set_id, sps));
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::storeGps(GeometryParameterSet&& gps)
{
  // todo(df): handle replacement semantics
  _gpss.emplace(std::make_pair(gps.gps_geom_parameter_set_id, gps));
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::storeAps(AttributeParameterSet&& aps)
{
  // todo(df): handle replacement semantics
  _apss.emplace(std::make_pair(aps.aps_attr_parameter_set_id, aps));
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::storeTileInventory(TileInventory&& inventory)
{
  // todo(df): handle replacement semantics
  _tileInventory = inventory;
}

//==========================================================================

void
PCCTMC3Decoder3::activateParameterSets(const GeometryBrickHeader& gbh)
{
  // HACK: assume activation of the first SPS and GPS
  // todo(df): parse brick header here for propper sps & gps activation
  //  -- this is currently inconsistent between trisoup and octree
  assert(!_spss.empty());
  assert(!_gpss.empty());
  _sps = &_spss.cbegin()->second;
  _gps = &_gpss.cbegin()->second;
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::activateParameterSets(const AttributeParamInventoryHdr& hdr)
{
  // HACK: assume activation of the first SPS and GPS
  // todo(df): parse brick header here for propper sps & gps activation
  //  -- this is currently inconsistent between trisoup and octree
  assert(!_spss.empty());
  assert(!_gpss.empty());
  _sps = &_spss.cbegin()->second;
  _gps = &_gpss.cbegin()->second;
}

//==========================================================================
// Initialise the point cloud storage and decode a single geometry slice.

int
PCCTMC3Decoder3::decodeGeometryBrick(const PayloadBuffer& buf)
{
  assert(buf.type == PayloadType::kGeometryBrick);
  std::cout << "positions bitstream size " << buf.size() << " B\n";

  // todo(df): replace with attribute mapping
  bool hasColour = std::any_of(
    _sps->attributeSets.begin(), _sps->attributeSets.end(),
    [](const AttributeDescription& desc) {
      return desc.attributeLabel == KnownAttributeLabel::kColour;
    });

  bool hasReflectance = std::any_of(
    _sps->attributeSets.begin(), _sps->attributeSets.end(),
    [](const AttributeDescription& desc) {
      return desc.attributeLabel == KnownAttributeLabel::kReflectance;
    });

  _currentPointCloud.clear();
  _currentPointCloud.addRemoveAttributes(hasColour, hasReflectance);

  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
  clock_user.start();

  int gbhSize, gbfSize;
  _gbh = parseGbh(*_sps, *_gps, buf, &gbhSize, &gbfSize);
  _prevSliceId = _sliceId;
  _sliceId = _gbh.geom_slice_id;
  _sliceOrigin = _gbh.geomBoxOrigin;

  // sanity check for loss detection
  if (_gbh.entropy_continuation_flag) {
    assert(!_firstSliceInFrame);
    assert(_gbh.prev_slice_id == _prevSliceId);
  } else {
    // forget (reset) all saved context state at boundary
    _ctxtMemOctreeGeom->reset();
    _ctxtMemPredGeom->reset();
    for (auto& ctxtMem : _ctxtMemAttrs)
      ctxtMem.reset();
  }

  // set default attribute values (in case an attribute data unit is lost)
  // NB: it is a requirement that geom_num_points_minus1 is correct
  _currentPointCloud.resize(_gbh.footer.geom_num_points_minus1 + 1);
  if (hasColour) {
    auto it = std::find_if(
      _outCloud.attrDesc.cbegin(), _outCloud.attrDesc.cend(),
      [](const AttributeDescription& desc) {
        return desc.attributeLabel == KnownAttributeLabel::kColour;
      });

    Vec3<attr_t> defAttrVal = 1 << (it->bitdepth - 1);
    if (!it->params.attr_default_value.empty())
      for (int k = 0; k < 3; k++)
        defAttrVal[k] = it->params.attr_default_value[k];
    for (int i = 0; i < _currentPointCloud.getPointCount(); i++)
      _currentPointCloud.setColor(i, defAttrVal);
  }

  if (hasReflectance) {
    auto it = std::find_if(
      _outCloud.attrDesc.cbegin(), _outCloud.attrDesc.cend(),
      [](const AttributeDescription& desc) {
        return desc.attributeLabel == KnownAttributeLabel::kReflectance;
      });
    attr_t defAttrVal = 1 << (it->bitdepth - 1);
    if (!it->params.attr_default_value.empty())
      defAttrVal = it->params.attr_default_value[0];
    for (int i = 0; i < _currentPointCloud.getPointCount(); i++)
      _currentPointCloud.setReflectance(i, defAttrVal);
  }

  // Calculate a tree level at which to stop
  // It should result in at most max points being decoded
  if (_params.decodeMaxPoints && _gps->octree_point_count_list_present_flag) {
    if (_params.decodeMaxPoints > _gbh.footer.geom_num_points_minus1)
      _params.minGeomNodeSizeLog2 = 0;
    else {
      auto it = std::lower_bound(
        std::next(_gbh.footer.octree_lvl_num_points_minus1.begin()),
        _gbh.footer.octree_lvl_num_points_minus1.end(),
        _params.decodeMaxPoints);

      _params.minGeomNodeSizeLog2 =
        std::distance(it, _gbh.footer.octree_lvl_num_points_minus1.end()) + 1;
    }
  }

  EntropyDecoder aec;
  aec.setBuffer(buf.size() - gbhSize - gbfSize, buf.data() + gbhSize);
  aec.enableBypassStream(_sps->cabac_bypass_stream_enabled_flag);
  aec.start();

  if (_gps->predgeom_enabled_flag)
    decodePredictiveGeometry(
      *_gps, _gbh, _currentPointCloud, &_posSph, *_ctxtMemPredGeom, aec);
  else if (!_gps->trisoup_enabled_flag) {
    if (!_params.minGeomNodeSizeLog2) {
      decodeGeometryOctree(
        *_gps, _gbh, _currentPointCloud, *_ctxtMemOctreeGeom, aec);
    } else {
      decodeGeometryOctreeScalable(
        *_gps, _gbh, _params.minGeomNodeSizeLog2, _currentPointCloud,
        *_ctxtMemOctreeGeom, aec);
    }
  } else {
    decodeGeometryTrisoup(
      *_gps, _gbh, _currentPointCloud, *_ctxtMemOctreeGeom, aec);
  }

  // At least the first slice's geometry has been decoded
  _firstSliceInFrame = false;

  clock_user.stop();

  auto total_user =
    std::chrono::duration_cast<std::chrono::milliseconds>(clock_user.count());
  std::cout << "positions processing time (user): "
            << total_user.count() / 1000.0 << " s\n";
  std::cout << std::endl;

  return 0;
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::decodeAttributeBrick(const PayloadBuffer& buf)
{
  assert(buf.type == PayloadType::kAttributeBrick);
  // todo(df): replace assertions with error handling
  assert(_sps);
  assert(_gps);

  // verify that this corresponds to the correct geometry slice
  AttributeBrickHeader abh = parseAbhIds(buf);
  assert(abh.attr_geom_slice_id == _sliceId);

  // todo(df): validate that sps activation is not changed via the APS
  const auto it_attr_aps = _apss.find(abh.attr_attr_parameter_set_id);

  assert(it_attr_aps != _apss.cend());
  const auto& attr_aps = it_attr_aps->second;

  assert(abh.attr_sps_attr_idx < _sps->attributeSets.size());
  const auto& attr_sps = _sps->attributeSets[abh.attr_sps_attr_idx];
  const auto& label = attr_sps.attributeLabel;

  // sanity check for loss detection
  if (_gbh.entropy_continuation_flag)
    assert(_gbh.prev_slice_id == _ctxtMemAttrSliceIds[abh.attr_sps_attr_idx]);

  // Ensure context arrays are allocated context arrays
  // todo(df): move this to sps activation
  _ctxtMemAttrSliceIds.resize(_sps->attributeSets.size());
  _ctxtMemAttrs.resize(_sps->attributeSets.size());

  // In order to determinet hat the attribute decoder is reusable, the abh
  // must be inspected.
  int abhSize;
  abh = parseAbh(*_sps, attr_aps, buf, &abhSize);

  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;

  // replace the attribute decoder if not compatible
  if (!_attrDecoder || !_attrDecoder->isReusable(attr_aps, abh))
    _attrDecoder = makeAttributeDecoder();

  clock_user.start();

  // Convert cartesian positions to spherical for use in attribute coding.
  // NB: this retains the original cartesian positions to restore afterwards
  std::vector<pcc::point_t> altPositions;
  if (attr_aps.spherical_coord_flag) {
    // If predgeom was used, re-use the internal positions rather than
    // calculating afresh.
    Box3<int> bboxRpl;
    if (_gps->predgeom_enabled_flag) {
      altPositions = _posSph;
      bboxRpl = Box3<int>(altPositions.begin(), altPositions.end());
    } else {
      altPositions.resize(_currentPointCloud.getPointCount());

      auto laserOrigin = _gbh.geomAngularOrigin(*_gps);
      bboxRpl = convertXyzToRpl(
        laserOrigin, _gps->angularTheta.data(), _gps->angularTheta.size(),
        &_currentPointCloud[0],
        &_currentPointCloud[0] + _currentPointCloud.getPointCount(),
        altPositions.data());
    }

    offsetAndScale(
      bboxRpl.min, attr_aps.attr_coord_scale, altPositions.data(),
      altPositions.data() + altPositions.size());

    _currentPointCloud.swapPoints(altPositions);
  }

  auto& ctxtMemAttr = _ctxtMemAttrs.at(abh.attr_sps_attr_idx);
  _attrDecoder->decode(
    *_sps, attr_sps, attr_aps, abh, _gbh.footer.geom_num_points_minus1,
    _params.minGeomNodeSizeLog2, buf.data() + abhSize, buf.size() - abhSize,
    ctxtMemAttr, _currentPointCloud);

  if (attr_aps.spherical_coord_flag)
    _currentPointCloud.swapPoints(altPositions);

  // Note the current sliceID for loss detection
  _ctxtMemAttrSliceIds[abh.attr_sps_attr_idx] = _sliceId;

  clock_user.stop();

  std::cout << label << "s bitstream size " << buf.size() << " B\n";

  auto total_user =
    std::chrono::duration_cast<std::chrono::milliseconds>(clock_user.count());
  std::cout << label
            << "s processing time (user): " << total_user.count() / 1000.0
            << " s\n";
  std::cout << std::endl;
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::decodeConstantAttribute(const PayloadBuffer& buf)
{
  assert(buf.type == PayloadType::kConstantAttribute);
  // todo(df): replace assertions with error handling
  assert(_sps);
  assert(_gps);

  ConstantAttributeDataUnit cadu = parseConstantAttribute(*_sps, buf);

  // verify that this corresponds to the correct geometry slice
  assert(cadu.constattr_geom_slice_id == _sliceId);

  assert(cadu.constattr_sps_attr_idx < _sps->attributeSets.size());
  const auto& attrDesc = _sps->attributeSets[cadu.constattr_sps_attr_idx];
  const auto& label = attrDesc.attributeLabel;

  // todo(df): replace with proper attribute mapping
  if (label == KnownAttributeLabel::kColour) {
    Vec3<attr_t> defAttrVal;
    for (int k = 0; k < 3; k++)
      defAttrVal[k] = attrDesc.params.attr_default_value[k];
    for (int i = 0; i < _currentPointCloud.getPointCount(); i++)
      _currentPointCloud.setColor(i, defAttrVal);
  }

  if (label == KnownAttributeLabel::kReflectance) {
    attr_t defAttrVal = attrDesc.params.attr_default_value[0];
    for (int i = 0; i < _currentPointCloud.getPointCount(); i++)
      _currentPointCloud.setReflectance(i, defAttrVal);
  }
}

//============================================================================

}  // namespace pcc
