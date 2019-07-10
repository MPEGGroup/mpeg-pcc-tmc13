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
#include "partitioning.h"
#include "pcc_chrono.h"

namespace pcc {

//============================================================================

int
PCCTMC3Encoder3::compress(
  const PCCPointSet3& inputPointCloud,
  EncoderParams* params,
  std::function<void(const PayloadBuffer&)> outputFn,
  PCCPointSet3* reconstructedCloud)
{
  fixupParameterSets(params);

  // Determine input bounding box (for SPS metadata) if not manually set
  if (params->sps.seq_bounding_box_whd == Vec3<int>{0}) {
    const auto& bbox = inputPointCloud.computeBoundingBox();
    for (int k = 0; k < 3; k++) {
      params->sps.seq_bounding_box_xyz0[k] = int(bbox.min[k]);

      // somehow determine the decoder's reconstructed points bounding box
      // and update sps accordingly.
      auto max_k = bbox.max[k] - bbox.min[k];
      max_k = std::round(max_k * params->sps.seq_source_geom_scale_factor);
      max_k = std::round(max_k / params->sps.seq_source_geom_scale_factor);

      // NB: plus one to convert to range
      params->sps.seq_bounding_box_whd[k] = int(max_k) + 1;
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

  // initial geometry IDs
  _tileId = 0;
  _sliceId = 0;
  _sliceOrigin = Vec3<int>{0};

  // Partition the input point cloud into tiles
  //  - quantize the input point cloud (without duplicate point removal)
  //  - inverse quantize the cloud above to get the initial-sized cloud
  //  - if tile partitioning is enabled,partitioning function produces
  //    vectors tileMaps which map tileIDs to point indexes.
  //    Compute the tile metadata for each partition.
  //  - if not, regard the whole input cloud as a single tile to facilitate
  //    slice partitioning subsequent
  //  todo(df):
  PartitionSet partitions;
  PCCPointSet3 quantizedInputCloud;
  quantizedInputCloud = quantization(inputPointCloud);

  std::vector<std::vector<int32_t>> tileMaps;
  if (params->partition.tileSize) {
    PCCPointSet3 inverseQuantizedCloud;
    Box3<int32_t> clampBox{{0, 0, 0}, {INT32_MAX, INT32_MAX, INT32_MAX}};
    quantizePositions(
      1.0 / _sps->seq_source_geom_scale_factor, 0, clampBox,
      quantizedInputCloud, &inverseQuantizedCloud);
    tileMaps = tilePartition(params->partition, inverseQuantizedCloud);

    // Get the bounding box of current tile and write it into tileInventory
    partitions.tileInventory.tiles.resize(tileMaps.size());
    for (int t = 0; t < tileMaps.size(); t++) {
      Box3<double> bbox = inverseQuantizedCloud.computeBoundingBox(
        tileMaps[t].begin(), tileMaps[t].end());

      auto& tileIvt = partitions.tileInventory.tiles[t];
      for (int k = 0; k < 3; k++) {
        tileIvt.tile_bounding_box_whd[k] = bbox.max[k] - bbox.min[k];
        tileIvt.tile_bounding_box_xyz0[k] =
          bbox.min[k] - _sps->seq_bounding_box_xyz0[k];
      }
    }
  } else {
    tileMaps.emplace_back();
    auto& tile = tileMaps.back();
    for (int i = 0; i < quantizedInputCloud.getPointCount(); i++)
      tile.push_back(i);
  }

  // If partitioning is not enabled, encode input as a single "partition"
  if (params->partition.method == PartitionMethod::kNone) {
    // todo(df): params->gps.geom_box_present_flag = false;
    _sliceOrigin = Vec3<int>{0};
    compressPartition(
      quantizedInputCloud, inputPointCloud, params, outputFn,
      reconstructedCloud);
    return 0;
  }

  // Partition the input point cloud
  //  - get the partitial cloud of each tile
  //  - partitioning function produces a list of point indexes, origin and
  //    optional tile metadata for each partition.
  //  - encode any tile metadata
  //  todo(df): consider requiring partitioning function to sort the input
  //            points and provide ranges rather than a set of indicies.
  do {
    Box3<int32_t> clampBox{{0, 0, 0}, {INT32_MAX, INT32_MAX, INT32_MAX}};

    for (int t = 0; t < tileMaps.size(); t++) {
      const auto& tile = tileMaps[t];

      // Get the point cloud of current tile and compute their bounding boxes
      PCCPointSet3 tileCloud;
      getSrcPartition(quantizedInputCloud, tileCloud, tile);
      Box3<double> bbox = tileCloud.computeBoundingBox();
      Vec3<int> tile_quantized_box_xyz0;
      for (int k = 0; k < 3; k++) {
        tile_quantized_box_xyz0[k] = int(bbox.min[k]);
      }

      // Move the tile cloud to coodinate origin
      // for the convenience of slice partitioning
      quantizePositions(
        1, tile_quantized_box_xyz0, clampBox, tileCloud, &tileCloud);

      //Slice partition of current tile
      std::vector<Partition> curSlices;
      switch (params->partition.method) {
      // NB: this method is handled earlier
      case PartitionMethod::kNone: return 1;

      case PartitionMethod::kUniformGeom:
        curSlices = partitionByUniformGeom(
          params->partition, tileCloud, t, _gps->trisoup_node_size_log2);
        break;

      case PartitionMethod::kOctreeUniform:
        curSlices = partitionByOctreeDepth(params->partition, tileCloud, t);
        break;
      }
      // Map slice indexes to tile indexes(the original indexes)
      for (int i = 0; i < curSlices.size(); i++) {
        for (int p = 0; p < curSlices[i].pointIndexes.size(); p++) {
          curSlices[i].pointIndexes[p] = tile[curSlices[i].pointIndexes[p]];
        }
      }
      // Adjust the point number of each slice
      // to the range between sliceMaxPoints and sliceMinPoints
      // by merge small slices and split large ones.
      refineSlices(params->partition, inputPointCloud, curSlices);

      partitions.slices.insert(
        partitions.slices.end(), curSlices.begin(), curSlices.end());
    }
    std::cout << "Slice number: " << partitions.slices.size() << std::endl;
  } while (0);

  if (partitions.tileInventory.tiles.size() > 1) {
    assert(partitions.tileInventory.tiles.size() == tileMaps.size());
    std::cout << "Tile number: " << tileMaps.size() << std::endl;
    outputFn(write(partitions.tileInventory));
  }

  // Encode each partition:
  //  - create a pointset comprising just the partitioned points
  //  - compress
  for (const auto& partition : partitions.slices) {
    // create partitioned point set
    PCCPointSet3 srcPartition;
    getSrcPartition(quantizedInputCloud, srcPartition, partition.pointIndexes);

    // Get the original partial cloud corresponding to each slice for recolor
    std::vector<int32_t> partitionOriginIdxes;
    for (int i = 0; i < partition.pointIndexes.size(); i++) {
      const auto& point = srcPartition[i];
      std::multimap<Vec3<double>, int32_t>::iterator pos;
      for (pos = quantizedToOrigin.lower_bound(point);
           pos != quantizedToOrigin.upper_bound(point); ++pos) {
        partitionOriginIdxes.push_back(pos->second);
      }
    }
    PCCPointSet3 partitionInOriginCloud;
    getSrcPartition(
      inputPointCloud, partitionInOriginCloud, partitionOriginIdxes);

    _sliceId = partition.sliceId;
    _tileId = partition.tileId;
    _sliceOrigin = partition.origin;
    compressPartition(
      srcPartition, partitionInOriginCloud, params, outputFn,
      reconstructedCloud);
  }

  return 0;
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::fixupParameterSets(EncoderParams* params)
{
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

  // the encoder writes out the slice origin in the GBH
  // NB: this may be disabled during encoding
  params->gps.geom_box_present_flag = true;
  params->gps.geom_box_log2_scale_present_flag = true;
  params->gps.gps_geom_box_log2_scale = 0;

  // fixup attribute parameters
  for (auto it : params->attributeIdxMap) {
    auto& attr_sps = params->sps.attributeSets[it.second];
    auto& attr_aps = params->aps[it.second];
    auto& attr_enc = params->attr[it.second];
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

    // the encoder options may not specify sufficient offsets for the number
    // of layers used by the sytax: extend with last value as appropriate
    if (!attr_enc.abh.attr_layer_qp_delta_luma.empty()) {
      int numLayers =
        attr_aps.attr_encoding == AttributeEncoding::kRAHTransform
        ? attr_aps.raht_depth + 1
        : attr_aps.num_detail_levels + 1;

      attr_enc.abh.attr_layer_qp_delta_luma.resize(
        numLayers, attr_enc.abh.attr_layer_qp_delta_luma.back());

      attr_enc.abh.attr_layer_qp_delta_chroma.resize(
        numLayers, attr_enc.abh.attr_layer_qp_delta_chroma.back());
    }
  }
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::compressPartition(
  const PCCPointSet3& inputPointCloud,
  const PCCPointSet3& originPartCloud,
  EncoderParams* params,
  std::function<void(const PayloadBuffer&)> outputFn,
  PCCPointSet3* reconstructedCloud)
{
  // geometry compression consists of the following stages:
  //  - prefilter/quantize geometry (non-normative)
  //  - encode geometry (single slice, id = 0)
  //  - recolour

  pointCloud.clear();
  pointCloud = inputPointCloud;

  // Offset the point cloud to account for (preset) _sliceOrigin.
  Vec3<double> sliceOriginD{double(_sliceOrigin[0]), double(_sliceOrigin[1]),
                            double(_sliceOrigin[2])};

  // The new maximum bounds of the offset cloud
  Vec3<int> maxBound{0};

  const size_t pointCount = pointCloud.getPointCount();
  for (size_t i = 0; i < pointCount; ++i) {
    const Vec3<double> point = (pointCloud[i] -= sliceOriginD);
    for (int k = 0; k < 3; ++k) {
      const int k_coord = int(point[k]);
      assert(k_coord >= 0);
      if (maxBound[k] < k_coord)
        maxBound[k] = k_coord;
    }
  }

  // todo(df): don't update maxBound if something is forcing the value?
  _sliceBoxWhd = maxBound;

  // geometry encoding
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
  if (_gps->geom_unique_points_flag || _gps->trisoup_node_size_log2 > 0) {
    recolour(
      params->recolour, originPartCloud, _sps->seq_source_geom_scale_factor,
      _sps->seq_bounding_box_xyz0, _sliceOrigin, &pointCloud);
  }

  // dump recoloured point cloud
  // todo(df): this needs to work with partitioned clouds
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
    const auto& attr_enc = params->attr[attrIdx];
    const auto& label = attr_sps.attributeLabel;

    PayloadBuffer payload(PayloadType::kAttributeBrick);

    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
    clock_user.start();

    // todo(df): move elsewhere?
    AttributeBrickHeader abh;
    abh.attr_attr_parameter_set_id = attr_aps.aps_attr_parameter_set_id;
    abh.attr_sps_attr_idx = attrIdx;
    abh.attr_geom_slice_id = _sliceId;
    abh.attr_qp_delta_luma = 0;
    abh.attr_qp_delta_chroma = 0;
    abh.attr_layer_qp_delta_luma = attr_enc.abh.attr_layer_qp_delta_luma;
    abh.attr_layer_qp_delta_chroma = attr_enc.abh.attr_layer_qp_delta_chroma;

    write(attr_aps, abh, &payload);

    AttributeEncoder attrEncoder;
    attrEncoder.encode(*_sps, attr_sps, attr_aps, abh, pointCloud, &payload);
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

  appendReconstructedPoints(reconstructedCloud);
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::encodeGeometryBrick(PayloadBuffer* buf)
{
  // todo(df): confirm minimum of 1 isn't needed
  int32_t maxBB =
    std::max({1, _sliceBoxWhd[0], _sliceBoxWhd[1], _sliceBoxWhd[2]});

  // the current node dimension (log2) encompasing maxBB
  int nodeSizeLog2 = ceillog2(maxBB + 1);

  GeometryBrickHeader gbh;
  gbh.geom_geom_parameter_set_id = _gps->gps_geom_parameter_set_id;
  gbh.geom_slice_id = _sliceId;
  gbh.geom_tile_id = std::max(0, _tileId);
  gbh.geomBoxOrigin = _sliceOrigin;
  gbh.geom_box_log2_scale = 0;
  gbh.geom_max_node_size_log2 = nodeSizeLog2;
  gbh.geom_num_points = int(pointCloud.getPointCount());
  write(*_gps, gbh, buf);

  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = int(pointCloud.getPointCount()) * 3 * 4 + 1024;
  EntropyEncoder arithmeticEncoder(maxAcBufLen, nullptr);
  arithmeticEncoder.enableBypassStream(_sps->cabac_bypass_stream_enabled_flag);
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
PCCTMC3Encoder3::appendReconstructedPoints(PCCPointSet3* reconstructedCloud)
{
  if (reconstructedCloud == nullptr) {
    return;
  }
  const size_t pointCount = pointCloud.getPointCount();
  size_t outIdx = reconstructedCloud->getPointCount();

  reconstructedCloud->addRemoveAttributes(
    pointCloud.hasColors(), pointCloud.hasReflectances());
  reconstructedCloud->resize(outIdx + pointCount);

  const double minPositionQuantizationScale = 0.0000000001;
  const double invScale =
    fabs(_sps->seq_source_geom_scale_factor) > minPositionQuantizationScale
    ? 1.0 / _sps->seq_source_geom_scale_factor
    : 1.0;
  for (size_t i = 0; i < pointCount; ++i, ++outIdx) {
    const auto quantizedPoint = pointCloud[i];
    auto& point = (*reconstructedCloud)[outIdx];
    for (size_t k = 0; k < 3; ++k) {
      point[k] = (quantizedPoint[k] + _sliceOrigin[k]) * invScale
        + _sps->seq_bounding_box_xyz0[k];
    }
    if (pointCloud.hasColors()) {
      reconstructedCloud->setColor(outIdx, pointCloud.getColor(i));
    }
    if (pointCloud.hasReflectances()) {
      reconstructedCloud->setReflectance(outIdx, pointCloud.getReflectance(i));
    }
  }
}

//----------------------------------------------------------------------------
// translates and scales inputPointCloud, storing the result in
// this->pointCloud for use by the encoding process.

PCCPointSet3
PCCTMC3Encoder3::quantization(const PCCPointSet3& inputPointCloud)
{
  PCCPointSet3 pointCloud0;
  pointCloud0.clear();
  quantizedToOrigin.clear();

  // Currently the sequence width/height/depth must be set
  assert(_sps->seq_bounding_box_whd != Vec3<int>{0});

  // Clamp all points to [clampBox.min, clampBox.max] after translation
  // and quantisation.
  Box3<int32_t> clampBox{{0, 0, 0}, {INT32_MAX, INT32_MAX, INT32_MAX}};
  if (_sps->seq_bounding_box_whd != Vec3<int>{0}) {
    // todo(df): this is icky (not to mention rounding issues)
    // NB: the sps seq_bounding_box_* uses unscaled co-ordinates => convert
    // NB: minus 1 to convert to max x/y/z position
    clampBox = Box3<int32_t>{{0, 0, 0}, _sps->seq_bounding_box_whd};
    for (int k = 0; k < 3; k++)
      clampBox.max[k] =
        int(ceil(clampBox.max[k] * _sps->seq_source_geom_scale_factor)) - 1;
  }

  if (_gps->geom_unique_points_flag) {
    quantizePositionsUniq(
      _sps->seq_source_geom_scale_factor, _sps->seq_bounding_box_xyz0,
      clampBox, inputPointCloud, &pointCloud0, quantizedToOrigin);
  } else {
    quantizePositions(
      _sps->seq_source_geom_scale_factor, _sps->seq_bounding_box_xyz0,
      clampBox, inputPointCloud, &pointCloud0);
  }

  // Offset the point cloud to account for (preset) _sliceOrigin.
  Vec3<double> sliceOriginD{double(_sliceOrigin[0]), double(_sliceOrigin[1]),
                            double(_sliceOrigin[2])};

  // The new maximum bounds of the offset cloud
  Vec3<int> maxBound{0};

  const size_t pointCount = pointCloud0.getPointCount();
  for (size_t i = 0; i < pointCount; ++i) {
    const Vec3<double> point = (pointCloud0[i] -= sliceOriginD);
    for (int k = 0; k < 3; ++k) {
      const int k_coord = int(point[k]);
      assert(k_coord >= 0);
      if (maxBound[k] < k_coord)
        maxBound[k] = k_coord;
    }
  }

  // todo(df): don't update maxBound if something is forcing the value?
  _sliceBoxWhd = maxBound;
  return pointCloud0;
}

//----------------------------------------------------------------------------
// get the partial point cloud according to required point indexes
void
PCCTMC3Encoder3::getSrcPartition(
  const PCCPointSet3& inputPointCloud,
  PCCPointSet3& srcPartition,
  std::vector<int32_t> Indexes)
{
  //PCCPointSet3 srcPartition;
  srcPartition.addRemoveAttributes(
    inputPointCloud.hasColors(), inputPointCloud.hasReflectances());

  int partitionSize = Indexes.size();
  srcPartition.resize(partitionSize);

  for (int i = 0; i < partitionSize; i++) {
    int inputIdx = Indexes[i];
    srcPartition[i] = inputPointCloud[inputIdx];

    if (inputPointCloud.hasColors())
      srcPartition.setColor(i, inputPointCloud.getColor(inputIdx));

    if (inputPointCloud.hasReflectances())
      srcPartition.setReflectance(i, inputPointCloud.getReflectance(inputIdx));
  }

  //return srcPartition;
}

//============================================================================

}  // namespace pcc
