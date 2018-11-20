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

#include "PCCTMC3Encoder.h"

#include <cassert>
#include <set>

#include "AttributeEncoder.h"
#include "PCCPointSetProcessing.h"
#include "geometry.h"
#include "io_hls.h"
#include "osspecific.h"
#include "pcc_chrono.h"

namespace pcc {

//============================================================================

void
PCCTMC3Encoder3::init()
{
  minPositions = 0.0;
  boundingBox.min = uint32_t(0);
  boundingBox.max = uint32_t(0);
  pointCloud.clear();
}

//----------------------------------------------------------------------------

int
PCCTMC3Encoder3::compress(
  const PCCPointSet3& inputPointCloud,
  EncoderParams* params,
  std::function<void(const PayloadBuffer&)> outputFn,
  PCCPointSet3* reconstructedCloud)
{
  init();

  // fixup parameter set IDs
  params->sps.sps_seq_parameter_set_id = 0;
  params->gps.gps_seq_parameter_set_id = 0;
  params->gps.gps_geom_parameter_set_id = 0;
  for (int i = 0; i < params->aps.size(); i++) {
    params->aps[i].aps_seq_parameter_set_id = 0;
    params->aps[i].aps_attr_parameter_set_id = i;
  }

  // development level / header
  params->sps.profileCompatibility.profile_compatibility_flags = 0;
  params->sps.level = 0;

  // the encoder writes out the minPositions in the GBH:
  params->gps.geom_box_present_flag = true;

  // fixup attribute parameters
  for (auto it : params->attributeIdxMap) {
    auto& attr_sps = params->sps.attributeSets[it.second];
    attr_sps.attr_instance_id = it.second;

    attr_sps.cicp_colour_primaries_idx = 2;
    attr_sps.cicp_transfer_characteristics_idx = 2;
    attr_sps.cicp_matrix_coefficients_idx = 2;
    attr_sps.cicp_video_full_range_flag = true;

    if (it.first == "color") {
      attr_sps.attr_num_dimensions = 3;
      attr_sps.attributeLabel = KnownAttributeLabel::kColour;
    }

    if (it.first == "reflectance") {
      attr_sps.attr_num_dimensions = 1;
      attr_sps.attributeLabel = KnownAttributeLabel::kReflectance;
    }
  }

  // placeholder to "activate" the parameter sets
  _sps = &params->sps;
  _gps = &params->gps;
  _aps.clear();
  for (const auto& aps : params->aps) {
    _aps.push_back(&aps);
  }

  // write out all parameter sets prior to encoding
  outputFn(write(*_sps));
  outputFn(write(*_gps));
  for (const auto aps : _aps) {
    outputFn(write(*aps));
  }

  // geometry compression consists of the following stages:
  //  - prefilter/quantize geometry (non-normative)
  //  - encode geometry (single slice, id = 0)
  //  - recolour

  // The quantization process will determine the bounding box
  quantization(inputPointCloud);

  // geometry encoding
  _sliceId = 0;
  _tileId = 0;

  if (1) {
    PayloadBuffer payload(PayloadType::kGeometryBrick);

    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
    clock_user.start();

    encodeGeometryBrick(&payload);

    clock_user.stop();

    double bpp = double(8 * payload.size()) / inputPointCloud.getPointCount();
    std::cout << "positions bitstream size " << payload.size() << " B (" << bpp
              << " bpp)\n";

    auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << "positions processing time (user): "
              << total_user.count() / 1000.0 << " s" << std::endl;

    outputFn(payload);
  }

  // recolouring

  // NB: recolouring is required if points are added / removed
  bool recolourNeeded =
    _gps->geom_unique_points_flag || _gps->trisoup_node_size_log2 > 0;

  if (recolourNeeded) {
    recolour(
      inputPointCloud, _sps->seq_source_geom_scale_factor, minPositions,
      &pointCloud);
  }

  // dump recoloured point cloud
  if (!params->postRecolorPath.empty()) {
    PCCPointSet3 tempPointCloud(pointCloud);
    tempPointCloud.convertYUVToRGB();
    tempPointCloud.write(params->postRecolorPath);
  }

  // attributeCoding

  // for each attribute
  for (const auto& it : params->attributeIdxMap) {
    int attrIdx = it.second;
    const auto& attr_sps = _sps->attributeSets[attrIdx];
    const auto& attr_aps = *_aps[attrIdx];
    const auto& label = attr_sps.attributeLabel;

    PayloadBuffer payload(PayloadType::kAttributeBrick);

    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
    clock_user.start();

    // todo(df): move elsewhere?
    AttributeBrickHeader abh;
    abh.attr_attr_parameter_set_id = attr_aps.aps_attr_parameter_set_id;
    abh.attr_sps_attr_idx = attrIdx;
    abh.attr_geom_slice_id = _sliceId;
    write(abh, &payload);

    AttributeEncoder attrEncoder;
    attrEncoder.encode(attr_sps, attr_aps, pointCloud, &payload);
    clock_user.stop();

    int coded_size = int(payload.size());
    double bpp = double(8 * coded_size) / inputPointCloud.getPointCount();
    std::cout << label << "s bitstream size " << coded_size << " B (" << bpp
              << " bpp)\n";

    auto time_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << label
              << "s processing time (user): " << time_user.count() / 1000.0
              << " s" << std::endl;

    outputFn(payload);
  }

  // prevent re-use of this sliceId:  the next slice (geometry + attributes)
  // should be distinguishable from the current slice.
  _sliceId++;

  reconstructedPointCloud(reconstructedCloud);

  return 0;
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::encodeGeometryBrick(PayloadBuffer* buf)
{
  // todo(df): confirm minimum of 1 isn't needed
  uint32_t maxBB =
    std::max({1u, boundingBox.max[0], boundingBox.max[1], boundingBox.max[2]});

  // the current node dimension (log2) encompasing maxBB
  int nodeSizeLog2 = ceillog2(maxBB + 1);

  GeometryBrickHeader gbh;
  gbh.geom_geom_parameter_set_id = _gps->gps_geom_parameter_set_id;
  gbh.geomBoxOrigin.x() = int(minPositions.x());
  gbh.geomBoxOrigin.y() = int(minPositions.y());
  gbh.geomBoxOrigin.z() = int(minPositions.z());
  gbh.geom_slice_id = _sliceId;
  gbh.geom_tile_id = _tileId;
  gbh.geom_box_log2_scale = 0;
  gbh.geom_max_node_size_log2 = nodeSizeLog2;
  gbh.geom_num_points = int(pointCloud.getPointCount());
  write(*_gps, gbh, buf);

  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = int(pointCloud.getPointCount()) * 3 * 4 + 1024;
  EntropyEncoder arithmeticEncoder(maxAcBufLen, nullptr);
  arithmeticEncoder.start();

  if (_gps->trisoup_node_size_log2 == 0) {
    encodeGeometryOctree(*_gps, gbh, pointCloud, &arithmeticEncoder);
  } else {
    encodeGeometryTrisoup(*_gps, gbh, pointCloud, &arithmeticEncoder);
  }

  uint32_t dataLen = arithmeticEncoder.stop();
  std::copy_n(arithmeticEncoder.buffer(), dataLen, std::back_inserter(*buf));
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::reconstructedPointCloud(PCCPointSet3* reconstructedCloud)
{
  if (reconstructedCloud == nullptr) {
    return;
  }
  const size_t pointCount = pointCloud.getPointCount();

  reconstructedCloud->addRemoveAttributes(
    pointCloud.hasColors(), pointCloud.hasReflectances());
  reconstructedCloud->resize(pointCount);

  const double minPositionQuantizationScale = 0.0000000001;
  const double invScale =
    fabs(_sps->seq_source_geom_scale_factor) > minPositionQuantizationScale
    ? 1.0 / _sps->seq_source_geom_scale_factor
    : 1.0;
  for (size_t i = 0; i < pointCount; ++i) {
    const auto quantizedPoint = pointCloud[i];
    auto& point = (*reconstructedCloud)[i];
    for (size_t k = 0; k < 3; ++k) {
      point[k] = quantizedPoint[k] * invScale + minPositions[k];
    }
    if (pointCloud.hasColors()) {
      reconstructedCloud->setColor(i, pointCloud.getColor(i));
    }
    if (pointCloud.hasReflectances()) {
      reconstructedCloud->setReflectance(i, pointCloud.getReflectance(i));
    }
  }
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::computeMinPositions(const PCCPointSet3& inputPointCloud)
{
  const size_t inputPointCount = inputPointCloud.getPointCount();
  minPositions = inputPointCloud[0];
  for (size_t i = 1; i < inputPointCount; ++i) {
    const auto point = inputPointCloud[i];
    for (int k = 0; k < 3; ++k) {
      if (minPositions[k] > point[k]) {
        minPositions[k] = point[k];
      }
    }
  }
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::quantization(const PCCPointSet3& inputPointCloud)
{
  // todo(df): allow overriding of minposition from CLI
  // todo(df): remove trisoup hack
  minPositions = PCCVector3D{0.0};
  if (_gps->trisoup_node_size_log2 == 0)
    computeMinPositions(inputPointCloud);

  PCCBox3<int32_t> clampBox{{0, 0, 0}, {INT32_MAX, INT32_MAX, INT32_MAX}};
  if (_gps->geom_unique_points_flag) {
    quantizePositionsUniq(
      _sps->seq_source_geom_scale_factor, -minPositions, clampBox,
      inputPointCloud, &pointCloud);
  } else {
    quantizePositions(
      _sps->seq_source_geom_scale_factor, -minPositions, clampBox,
      inputPointCloud, &pointCloud);
  }

  const size_t pointCount = pointCloud.getPointCount();
  boundingBox.min = uint32_t(0);
  boundingBox.max = uint32_t(0);
  for (size_t i = 0; i < pointCount; ++i) {
    const PCCVector3D point = pointCloud[i];
    for (int k = 0; k < 3; ++k) {
      const uint32_t coord = uint32_t(point[k]);
      if (boundingBox.max[k] < coord) {
        boundingBox.max[k] = coord;
      }
    }
  }
}

//============================================================================

}  // namespace pcc
