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

#ifndef PCCTMC3Encoder_h
#define PCCTMC3Encoder_h

#include <assert.h>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include "ringbuf.h"

#include "AttributeEncoder.h"
#include "OctreeNeighMap.h"
#include "PayloadBuffer.h"
#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "PCCPointSetProcessing.h"
#include "PCCTMC3Common.h"
#include "hls.h"
#include "io_hls.h"
#include "osspecific.h"
#include "pcc_chrono.h"

#include "ArithmeticCodec.h"
#include "tables.h"

namespace pcc {

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
};

//============================================================================

class PCCTMC3Encoder3 {
public:
  PCCTMC3Encoder3() { init(); }
  PCCTMC3Encoder3(const PCCTMC3Encoder3&) = default;
  PCCTMC3Encoder3& operator=(const PCCTMC3Encoder3& rhs) = default;
  ~PCCTMC3Encoder3() = default;
  void init()
  {
    minPositions = 0.0;
    boundingBox.min = uint32_t(0);
    boundingBox.max = uint32_t(0);
    pointCloud.clear();
  }

  int compress(
    const PCCPointSet3& inputPointCloud,
    EncoderParams* params,
    std::function<void(const PayloadBuffer&)> outputFn,
    PCCPointSet3* reconstructedCloud = nullptr)
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
    //  - recolour
    //  - encode geometry
    if (_gps->geom_codec_type == GeometryCodecType::kOctree) {
      quantization(inputPointCloud);
    }
    if (_gps->geom_codec_type == GeometryCodecType::kTriSoup) {
      // todo(?): this ought to be replaced with quantization(...)
      if (!inputPointCloud.write("verticesInWorldCoords.ply", true)) {
        std::cerr << "Failed to write verticesInWorldCoords.ply\n";
        return -1;
      }
      std::string cmd =
        "TMC1_coordinateTransform"
        " -inworld verticesInWorldCoords.ply"
        " -outframe verticesInFrameCoords.ply"
        " -depth "
        + std::to_string(_gps->trisoup_depth) + " -scale "
        + std::to_string(_sps->donotuse_trisoup_int_to_orig_scale)
        + " -format binary_little_endian";
      if (int err = system(cmd.c_str())) {
        std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
        return -1;
      }
    }

    // geometry encoding

    if (1) {
      PayloadBuffer payload(PayloadType::kGeometryBrick);

      pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
      clock_user.start();

      if (_gps->geom_codec_type == GeometryCodecType::kOctree) {
        encodePositions(&payload);
      }
      if (_gps->geom_codec_type == GeometryCodecType::kTriSoup) {
        bool hasColor = inputPointCloud.hasColors();
        bool hasReflectance = inputPointCloud.hasReflectances();
        if (encodeTrisoup(hasColor, hasReflectance, &payload))
          return -1;
        if (recolorTrisoup(inputPointCloud))
          return -1;
      }

      clock_user.stop();

      double bpp =
        double(8 * payload.size()) / inputPointCloud.getPointCount();
      std::cout << "positions bitstream size " << payload.size() << " B ("
                << bpp << " bpp)\n";

      auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
        clock_user.count());
      std::cout << "positions processing time (user): "
                << total_user.count() / 1000.0 << " s" << std::endl;

      outputFn(payload);
    }

    // attributeCoding

    // dump recoloured point cloud
    if (!params->postRecolorPath.empty()) {
      PCCPointSet3 tempPointCloud(pointCloud);
      tempPointCloud.convertYUVToRGB();
      tempPointCloud.write(params->postRecolorPath);
    }

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
      abh.attr_geom_brick_id = 0;
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

    reconstructedPointCloud(reconstructedCloud);

    return 0;
  }

private:
  void reconstructedPointCloud(PCCPointSet3* reconstructedCloud)
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

  //--------------------------------------------------------------------------

  void encodePositions(PayloadBuffer* buf);

  //--------------------------------------------------------------------------

  int encodeTrisoup(bool hasColor, bool hasReflectance, PayloadBuffer* buf)
  {
    // todo(?): port TMC1 geometry encoder
    // todo(?): fix TMC1 to output where it is told not requiring an extra
    //          level of guesswork to find the output file.
    pcc::mkdir("outbin");
    std::string cmd =
      "TMC1_geometryEncode"
      " -inframe verticesInFrameCoords.ply"
      " -outbin outbin"
      " -outref refinedVertices.ply"
      " -depth "
      + std::to_string(_gps->trisoup_depth) + " -level "
      + std::to_string(_gps->trisoup_triangle_level)
      + " -format binary_little_endian";
    if (int err = system(cmd.c_str())) {
      std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
      return -1;
    }

    // todo(?): port interpolator
    cmd =
      "TMC1_voxelize"
      " -inref refinedVertices.ply"
      " -outvox voxelizedVertices.ply"
      " -depth "
      + std::to_string(_gps->trisoup_depth) + " -format binary_little_endian";
    if (int err = system(cmd.c_str())) {
      std::cerr << "Failed (" << err << ") to run :" << cmd << '\n';
      return -1;
    }

    // Read in quantizedPointCloud file.
    if (!pointCloud.read("voxelizedVertices.ply")) {
      std::cerr << "Failed to open voxelizedVertices.ply\n";
      return -1;
    }
    pointCloud.addRemoveAttributes(hasColor, hasReflectance);

    // encapsulate the TMC1 "bitstream"
    std::ifstream fin("outbin/outbin_0000.bin", std::ios::binary);
    if (!fin.is_open()) {
      std::cerr << "failed to open tmc1 bitstream\n";
      return -1;
    }

    std::copy(
      std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>(),
      std::back_inserter(*buf));
    return 0;
  }

  //--------------------------------------------------------------------------

  void computeMinPositions(const PCCPointSet3& inputPointCloud)
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

  int quantization(const PCCPointSet3& inputPointCloud)
  {
    computeMinPositions(inputPointCloud);

    bool hasColour = inputPointCloud.hasColors();
    bool hasReflectance = inputPointCloud.hasReflectances();

    pointCloud.resize(0);
    pointCloud.addRemoveAttributes(hasColour, hasReflectance);

    const size_t inputPointCount = inputPointCloud.getPointCount();
    if (_gps->geom_unique_points_flag) {  // retain unique quantized points
      // filter points based on quantized positions
      std::set<PCCVector3<int64_t>> retainedPoints;
      for (size_t i = 0; i < inputPointCount; ++i) {
        const PCCVector3D point = (inputPointCloud[i] - minPositions)
          * _sps->seq_source_geom_scale_factor;
        PCCVector3<int64_t> quantizedPoint;
        quantizedPoint[0] = int64_t(std::round(point[0]));
        quantizedPoint[1] = int64_t(std::round(point[1]));
        quantizedPoint[2] = int64_t(std::round(point[2]));
        const auto it = retainedPoints.find(quantizedPoint);
        if (it == retainedPoints.end()) {
          retainedPoints.insert(quantizedPoint);
        }
      }
      // compute reconstructed point cloud
      pointCloud.resize(retainedPoints.size());
      const double invScale = 1.0 / _sps->seq_source_geom_scale_factor;
      size_t pointCounter = 0;
      for (const auto& quantizedPoint : retainedPoints) {
        auto& point = pointCloud[pointCounter++];
        for (size_t k = 0; k < 3; ++k) {
          point[k] = double(quantizedPoint[k]) * invScale + minPositions[k];
        }
      }

      if (hasColour) {  // transfer colors
        if (!PCCTransfertColors(inputPointCloud, pointCloud)) {
          std::cout << "Error: can't transfer colors!" << std::endl;
          return -1;
        }
      }

      if (hasReflectance) {  // transfer reflectances
        if (!PCCTransfertReflectances(inputPointCloud, pointCloud)) {
          std::cout << "Error: can't transfer reflectances!" << std::endl;
          return -1;
        }
      }

      // compute quantized coordinates
      pointCounter = 0;
      for (const auto& quantizedPoint : retainedPoints) {
        auto& point = pointCloud[pointCounter++];
        for (size_t k = 0; k < 3; ++k) {
          point[k] = double(quantizedPoint[k]);
        }
      }
    } else {
      // quantized point cloud
      pointCloud.resize(inputPointCount);
      for (size_t i = 0; i < inputPointCount; ++i) {
        const PCCVector3D point = (inputPointCloud[i] - minPositions)
          * _sps->seq_source_geom_scale_factor;
        auto& quantizedPoint = pointCloud[i];
        quantizedPoint[0] = std::round(point[0]);
        quantizedPoint[1] = std::round(point[1]);
        quantizedPoint[2] = std::round(point[2]);
        if (hasColour) {
          pointCloud.setColor(i, inputPointCloud.getColor(i));
        }
        if (hasReflectance) {
          pointCloud.setReflectance(i, inputPointCloud.getReflectance(i));
        }
      }
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
    return 0;
  }

  //--------------------------------------------------------------------------

  int recolorTrisoup(const PCCPointSet3& inputPointCloud)
  {
    // Transform from integer to original coordinates.
    const size_t pointCount = pointCloud.getPointCount();
    for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
      pointCloud[pointIndex] *= _sps->donotuse_trisoup_int_to_orig_scale;
      // todo(pchou): should add intToOrigTranslation here.
    }

    // Recolour attributes.
    if (inputPointCloud.hasColors()) {
      if (!PCCTransfertColors(inputPointCloud, pointCloud)) {
        std::cout << "Error: can't transfer colors!" << std::endl;
        return -1;
      }
    }
    if (inputPointCloud.hasReflectances()) {
      if (!PCCTransfertReflectances(inputPointCloud, pointCloud)) {
        std::cout << "Error: can't transfer reflectance!" << std::endl;
        return -1;
      }
    }

    // Transform from original to integer coordinates.
    for (size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
      // todo(pchou): should subtract intToOrigTranslation here.
      // todo(?): don't use division here
      pointCloud[pointIndex] /= _sps->donotuse_trisoup_int_to_orig_scale;
    }

    return 0;
  }

private:
  PCCVector3D minPositions;
  PCCBox3<uint32_t> boundingBox;
  PCCPointSet3 pointCloud;

  // The active parameter sets
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;
  std::vector<const AttributeParameterSet*> _aps;
};
}  // namespace pcc

#endif /* PCCTMC3Encoder_h */
