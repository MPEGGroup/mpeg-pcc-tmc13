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
#include <stdexcept>

#include "Attribute.h"
#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"
#include "io_hls.h"
#include "osspecific.h"
#include "partitioning.h"
#include "pcc_chrono.h"
#include "ply.h"

namespace pcc {

//============================================================================

PCCTMC3Encoder3::PCCTMC3Encoder3() : _frameCounter(-1)
{}

//============================================================================

int
PCCTMC3Encoder3::compress(
  const PCCPointSet3& inputPointCloud,
  EncoderParams* params,
  PCCTMC3Encoder3::Callbacks* callback,
  PCCPointSet3* reconstructedCloud)
{
  // start of frame
  _frameCounter++;

  if (_frameCounter == 0) {
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

    // Determine the ladar head position relative to the sequence bounding box
    // NB: currently the sps box offset is unscaled
    if (params->gps.geom_angular_mode_enabled_flag) {
      auto origin = params->sps.seq_bounding_box_xyz0;
      auto scale = params->sps.seq_source_geom_scale_factor;
      params->gps.geom_angular_lidar_head_position -= origin;
      params->gps.geom_angular_lidar_head_position *= scale;
    }
  }

  // placeholder to "activate" the parameter sets
  _sps = &params->sps;
  _gps = &params->gps;
  _aps.clear();
  for (const auto& aps : params->aps) {
    _aps.push_back(&aps);
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

  // determine the dist2 parameters based upon the quantized point cloud
  // todo(df): do this only if no dist2 value is set but should be
  bool calcDist2 = false;
  for (auto& aps : params->aps)
    calcDist2 |= aps.num_detail_levels > 0;

  if (calcDist2) {
    int maxNodeSizeLog2 =
      ceillog2(std::max({_sliceBoxWhd[0], _sliceBoxWhd[1], _sliceBoxWhd[2]}));

    // workout an intrinsic dist2
    int baseDist2 = estimateDist2(quantizedInputCloud, maxNodeSizeLog2);

    // generate dist2 series for each aps
    for (auto& aps : params->aps) {
      if (aps.num_detail_levels == 0)
        continue;

      aps.dist2.resize(aps.num_detail_levels);

      int64_t d2 = baseDist2;
      for (int i = 0; i < aps.num_detail_levels; ++i) {
        aps.dist2[i] = d2;
        d2 = 4 * d2;
      }
    }
  }

  // write out all parameter sets prior to encoding
  callback->onOutputBuffer(write(*_sps));
  callback->onOutputBuffer(write(*_gps));
  for (const auto aps : _aps) {
    callback->onOutputBuffer(write(*aps));
  }

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
      Box3<int32_t> bbox = inverseQuantizedCloud.computeBoundingBox(
        tileMaps[t].begin(), tileMaps[t].end());

      auto& tileIvt = partitions.tileInventory.tiles[t];
      for (int k = 0; k < 3; k++) {
        tileIvt.tile_bounding_box_whd[k] = bbox.max[k] - bbox.min[k];
        tileIvt.tile_bounding_box_xyz0[k] =
          bbox.min[k] - _sps->seq_bounding_box_xyz0[k];
      }
    }
  } else {
    _sliceOrigin = quantizedInputCloud.computeBoundingBox().min;
    tileMaps.emplace_back();
    auto& tile = tileMaps.back();
    for (int i = 0; i < quantizedInputCloud.getPointCount(); i++)
      tile.push_back(i);
  }

  // don't partition if partitioning would result in a single slice.
  auto partitionMethod = params->partition.method;
  if (quantizedInputCloud.getPointCount() < params->partition.sliceMaxPoints)
    partitionMethod = PartitionMethod::kNone;

  // If partitioning is not enabled, encode input as a single "partition"
  if (partitionMethod == PartitionMethod::kNone) {
    compressPartition(
      quantizedInputCloud, inputPointCloud, params, callback,
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
      Box3<int32_t> bbox = tileCloud.computeBoundingBox();
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
      switch (partitionMethod) {
      // NB: this method is handled earlier
      case PartitionMethod::kNone: return 1;

      case PartitionMethod::kUniformGeom:
        curSlices = partitionByUniformGeom(
          params->partition, tileCloud, t, _gps->trisoup_node_size_log2);
        break;

      case PartitionMethod::kUniformSquare:
        curSlices = partitionByUniformSquare(
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
      if (partitionMethod == PartitionMethod::kUniformSquare)
        refineSlicesByAdjacentInfo(
          params->partition, quantizedInputCloud, curSlices);
      else
        refineSlices(params->partition, quantizedInputCloud, curSlices);

      partitions.slices.insert(
        partitions.slices.end(), curSlices.begin(), curSlices.end());
    }
    std::cout << "Slice number: " << partitions.slices.size() << std::endl;
  } while (0);

  if (partitions.tileInventory.tiles.size() > 1) {
    assert(partitions.tileInventory.tiles.size() == tileMaps.size());
    std::cout << "Tile number: " << tileMaps.size() << std::endl;
    callback->onOutputBuffer(write(partitions.tileInventory));
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
      std::multimap<point_t, int32_t>::iterator pos;
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
      srcPartition, partitionInOriginCloud, params, callback,
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

  // use one bit to indicate frame boundaries
  params->sps.log2_max_frame_idx = 1;

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

    // the encoder options may not specify sufficient offsets for the number
    // of layers used by the sytax: extend with last value as appropriate
    int numLayers = std::max(
      attr_enc.abh.attr_layer_qp_delta_luma.size(),
      attr_enc.abh.attr_layer_qp_delta_chroma.size());

    int lastDeltaLuma = 0;
    if (!attr_enc.abh.attr_layer_qp_delta_luma.empty())
      lastDeltaLuma = attr_enc.abh.attr_layer_qp_delta_luma.back();

    int lastDeltaChroma = 0;
    if (!attr_enc.abh.attr_layer_qp_delta_chroma.empty())
      lastDeltaChroma = attr_enc.abh.attr_layer_qp_delta_chroma.back();

    attr_enc.abh.attr_layer_qp_delta_luma.resize(numLayers, lastDeltaLuma);
    attr_enc.abh.attr_layer_qp_delta_chroma.resize(numLayers, lastDeltaChroma);
  }
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::compressPartition(
  const PCCPointSet3& inputPointCloud,
  const PCCPointSet3& originPartCloud,
  EncoderParams* params,
  PCCTMC3Encoder3::Callbacks* callback,
  PCCPointSet3* reconstructedCloud)
{
  // geometry compression consists of the following stages:
  //  - prefilter/quantize geometry (non-normative)
  //  - encode geometry (single slice, id = 0)
  //  - recolour

  pointCloud.clear();
  pointCloud = inputPointCloud;

  // Offset the point cloud to account for (preset) _sliceOrigin.
  // The new maximum bounds of the offset cloud
  Vec3<int> maxBound{0};

  const size_t pointCount = pointCloud.getPointCount();
  for (size_t i = 0; i < pointCount; ++i) {
    const point_t point = (pointCloud[i] -= _sliceOrigin);
    for (int k = 0; k < 3; ++k) {
      const int k_coord = int(point[k]);
      assert(k_coord >= 0);
      if (maxBound[k] < k_coord)
        maxBound[k] = k_coord;
    }
  }

  // todo(df): don't update maxBound if something is forcing the value?
  // NB: size is max - min + 1
  _sliceBoxWhd = maxBound + 1;

  // geometry encoding
  if (1) {
    PayloadBuffer payload(PayloadType::kGeometryBrick);

    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
    clock_user.start();

    encodeGeometryBrick(params, &payload);

    clock_user.stop();

    double bpp = double(8 * payload.size()) / inputPointCloud.getPointCount();
    std::cout << "positions bitstream size " << payload.size() << " B (" << bpp
              << " bpp)\n";

    auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << "positions processing time (user): "
              << total_user.count() / 1000.0 << " s" << std::endl;

    callback->onOutputBuffer(payload);
  }

  // verify that the per-level slice constraint has been met
  // todo(df): avoid hard coded value here (should be level dependent)
  if (params->enforceLevelLimits)
    if (pointCloud.getPointCount() > 1100000)
      throw std::runtime_error(
        std::string("level slice point count limit (1100000) exceeded: ")
        + std::to_string(pointCloud.getPointCount()));

  // recolouring

  // NB: recolouring is required if points are added / removed
  if (_gps->geom_unique_points_flag || _gps->trisoup_node_size_log2 > 0) {
    for (const auto& attr_sps : _sps->attributeSets) {
      recolour(
        attr_sps, params->recolour, originPartCloud,
        _sps->seq_source_geom_scale_factor, _sps->seq_bounding_box_xyz0,
        _sliceOrigin, &pointCloud);
    }
  }

  // dump recoloured point cloud
  // todo(df): this needs to work with partitioned clouds
  callback->onPostRecolour(pointCloud);

  // attributeCoding
  auto attrEncoder = makeAttributeEncoder();

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

    abh.attr_region_qp_present_flag = false;
    abh.attr_region_qp_origin = Vec3<int>{0};
    abh.attr_region_qp_whd = Vec3<int>{0};
    abh.attr_region_qp_delta = 0;

    write(attr_aps, abh, &payload);

    // replace the attribute encoder if not compatible
    if (!attrEncoder->isReusable(attr_aps))
      attrEncoder = makeAttributeEncoder();

    attrEncoder->encode(*_sps, attr_sps, attr_aps, abh, pointCloud, &payload);
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

    callback->onOutputBuffer(payload);
  }

  // prevent re-use of this sliceId:  the next slice (geometry + attributes)
  // should be distinguishable from the current slice.
  _sliceId++;

  appendReconstructedPoints(reconstructedCloud);
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::encodeGeometryBrick(
  const EncoderParams* params, PayloadBuffer* buf)
{
  GeometryBrickHeader gbh;
  gbh.geom_geom_parameter_set_id = _gps->gps_geom_parameter_set_id;
  gbh.geom_slice_id = _sliceId;
  gbh.geom_tile_id = std::max(0, _tileId);
  gbh.frame_idx = _frameCounter & ((1 << _sps->log2_max_frame_idx) - 1);
  gbh.geomBoxOrigin = _sliceOrigin;
  gbh.geom_box_log2_scale = 0;

  if (!_gps->implicit_qtbt_enabled_flag) {
    // NB: A minimum whd of 2 means there is always at least 1 tree level
    int32_t maxWhd =
      std::max({2, _sliceBoxWhd[0], _sliceBoxWhd[1], _sliceBoxWhd[2]});
    gbh.geom_max_node_size_log2 = ceillog2(maxWhd);
  } else {
    // different node dimension for xyz, for the purpose of implicit qtbt
    for (int k = 0; k < 3; k++)
      gbh.geom_max_node_size_log2_xyz[k] =
        ceillog2(std::max(2, _sliceBoxWhd[k]));
  }

  gbh.geom_num_points_minus1 = int(pointCloud.getPointCount()) - 1;
  gbh.geom_slice_qp_offset = params->gbh.geom_slice_qp_offset;
  gbh.geom_octree_qp_offset_depth = params->gbh.geom_octree_qp_offset_depth;
  gbh.geom_octree_parallel_max_node_size_log2 =
    params->gbh.geom_octree_parallel_max_node_size_log2;

  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = int(pointCloud.getPointCount()) * 3 * 4 + 1024;

  // The number of entropy substreams is 1 + parallel_max_node_size_log2
  // NB: the two substream case is syntactically prohibited
  std::vector<std::unique_ptr<EntropyEncoder>> arithmeticEncoders;
  for (int i = 0; i < 1 + gbh.geom_octree_parallel_max_node_size_log2; i++) {
    arithmeticEncoders.emplace_back(new EntropyEncoder(maxAcBufLen, nullptr));
    auto& aec = arithmeticEncoders.back();
    aec->enableBypassStream(_sps->cabac_bypass_stream_enabled_flag);
    aec->start();
  }

  if (_gps->trisoup_node_size_log2 == 0) {
    encodeGeometryOctree(*_gps, gbh, pointCloud, arithmeticEncoders);
  } else {
    // limit the number of points to the slice limit
    // todo(df): this should be derived from the level
    gbh.geom_num_points_minus1 = params->partition.sliceMaxPoints - 1;
    encodeGeometryTrisoup(*_gps, gbh, pointCloud, arithmeticEncoders);
  }

  // update the header with the actual number of points coded
  gbh.geom_num_points_minus1 = pointCloud.getPointCount() - 1;

  // determine the length of each sub-stream
  for (auto& arithmeticEncoder : arithmeticEncoders) {
    auto dataLen = arithmeticEncoder->stop();
    gbh.geom_octree_parallel_bitstream_offsets.push_back(dataLen);
  }

  // determine the number of bits to use for the offset fields
  // NB: don't include the last offset since it isn't signalled
  if (gbh.geom_octree_parallel_max_node_size_log2 > 0) {
    size_t maxOffset = *std::max_element(
      gbh.geom_octree_parallel_bitstream_offsets.begin(),
      std::prev(gbh.geom_octree_parallel_bitstream_offsets.end()));
    gbh.geom_octree_parallel_max_offset_log2 = ceillog2(maxOffset + 1);
  }

  // assemble data unit
  write(*_sps, *_gps, gbh, buf);
  for (int i = 0; i < 1 + gbh.geom_octree_parallel_max_node_size_log2; i++) {
    auto& aec = arithmeticEncoders[i];
    auto dataLen = gbh.geom_octree_parallel_bitstream_offsets[i];
    std::copy_n(aec->buffer(), dataLen, std::back_inserter(*buf));
  }
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

  for (size_t i = 0; i < pointCount; ++i, ++outIdx) {
    (*reconstructedCloud)[outIdx] = pointCloud[i] + _sliceOrigin;

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
  // The new maximum bounds of the offset cloud
  Vec3<int> maxBound{0};

  const size_t pointCount = pointCloud0.getPointCount();
  for (size_t i = 0; i < pointCount; ++i) {
    const point_t point = (pointCloud0[i] -= _sliceOrigin);
    for (int k = 0; k < 3; ++k) {
      const int k_coord = int(point[k]);
      assert(k_coord >= 0);
      if (maxBound[k] < k_coord)
        maxBound[k] = k_coord;
    }
  }

  // todo(df): don't update maxBound if something is forcing the value?
  // NB: size is max - min + 1
  _sliceBoxWhd = maxBound + 1;
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
