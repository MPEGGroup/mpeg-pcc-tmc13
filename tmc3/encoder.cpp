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
#include <limits>
#include <set>
#include <stdexcept>

#include "Attribute.h"
#include "AttributeCommon.h"
#include "coordinate_conversion.h"
#include "geometry_params.h"
#include "hls.h"
#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"
#include "geometry_predictive.h"
#include "io_hls.h"
#include "osspecific.h"
#include "partitioning.h"
#include "pcc_chrono.h"
#include "ply.h"

namespace pcc {

//============================================================================

PCCPointSet3
getPartition(const PCCPointSet3& src, const std::vector<int32_t>& indexes);

PCCPointSet3 getPartition(
  const PCCPointSet3& src,
  const SrcMappedPointSet& map,
  const std::vector<int32_t>& indexes);

//============================================================================

PCCTMC3Encoder3::PCCTMC3Encoder3() : _frameCounter(-1)
{
  _ctxtMemOctreeGeom.reset(new GeometryOctreeContexts);
  _ctxtMemPredGeom.reset(new PredGeomContexts);
}

//----------------------------------------------------------------------------

PCCTMC3Encoder3::~PCCTMC3Encoder3() = default;

//============================================================================

int
PCCTMC3Encoder3::compress(
  const PCCPointSet3& inputPointCloud,
  EncoderParams* params,
  PCCTMC3Encoder3::Callbacks* callback,
  CloudFrame* reconCloud)
{
  // start of frame
  _frameCounter++;

  if (_frameCounter == 0) {
    // Angular predictive geometry coding needs to determine spherical
    // positions.  To avoid quantization of the input disturbing this:
    //  - sequence scaling is replaced by decimation of the input
    //  - any user-specified global scaling is honoured
    _inputDecimationScale = 1.;
    if (params->gps.predgeom_enabled_flag) {
      _inputDecimationScale = params->codedGeomScale;
      params->codedGeomScale /= params->seqGeomScale;
      params->seqGeomScale = 1.;
    }

    deriveParameterSets(params);
    fixupParameterSets(params);

    _srcToCodingScale = params->codedGeomScale;

    // Determine input bounding box (for SPS metadata) if not manually set
    Box3<int> bbox;
    if (params->autoSeqBbox)
      bbox = inputPointCloud.computeBoundingBox();
    else {
      bbox.min = params->sps.seqBoundingBoxOrigin;
      bbox.max = bbox.min + params->sps.seqBoundingBoxSize - 1;
    }

    // Note whether the bounding box size is defined
    // todo(df): set upper limit using level
    bool bboxSizeDefined = params->sps.seqBoundingBoxSize > 0;
    if (!bboxSizeDefined)
      params->sps.seqBoundingBoxSize = (1 << 21) - 1;

    // Then scale the bounding box to match the reconstructed output
    for (int k = 0; k < 3; k++) {
      auto min_k = bbox.min[k];
      auto max_k = bbox.max[k];

      // the sps bounding box is in terms of the conformance scale
      // not the source scale.
      // NB: plus one to convert to range
      min_k = std::round(min_k * params->seqGeomScale);
      max_k = std::round(max_k * params->seqGeomScale);
      params->sps.seqBoundingBoxOrigin[k] = min_k;
      params->sps.seqBoundingBoxSize[k] = max_k - min_k + 1;

      // Compensate the sequence origin such that source point (0,0,0) coded
      // as P_c is reconstructed as (0,0,0):
      //   0 = P_c * globalScale + seqOrigin
      auto gs = Rational(params->sps.globalScale);
      int rem = params->sps.seqBoundingBoxOrigin[k] % gs.numerator;
      rem += rem < 0 ? gs.numerator : 0;
      params->sps.seqBoundingBoxOrigin[k] -= rem;
      params->sps.seqBoundingBoxSize[k] += rem;

      // Convert the origin to coding coordinate system
      _originInCodingCoords[k] = params->sps.seqBoundingBoxOrigin[k];
      _originInCodingCoords[k] /= double(gs);
    }

    // Determine the number of bits to signal the bounding box
    params->sps.sps_bounding_box_offset_bits =
      numBits(params->sps.seqBoundingBoxOrigin.abs().max());

    params->sps.sps_bounding_box_size_bits = bboxSizeDefined
      ? numBits(params->sps.seqBoundingBoxSize.abs().max())
      : 0;

    // Determine the lidar head position in coding coordinate system
    params->gps.gpsAngularOrigin *= _srcToCodingScale;
    params->gps.gpsAngularOrigin -= _originInCodingCoords;

    // determine the scale factors based on a characteristic of the
    // acquisition system
    if (params->gps.geom_angular_mode_enabled_flag) {
      auto gs = Rational(params->sps.globalScale);
      int maxX = (params->sps.seqBoundingBoxSize[0] - 1) / double(gs);
      int maxY = (params->sps.seqBoundingBoxSize[1] - 1) / double(gs);
      auto& origin = params->gps.gpsAngularOrigin;
      int rx = std::max(std::abs(origin[0]), std::abs(maxX - origin[0]));
      int ry = std::max(std::abs(origin[1]), std::abs(maxY - origin[1]));
      int r = std::max(rx, ry);
      int twoPi = 25735;
      int maxLaserIdx = params->gps.numLasers() - 1;

      if (params->gps.predgeom_enabled_flag) {
        auto& gps = params->gps;
        twoPi = 1 << (gps.geom_angular_azimuth_scale_log2_minus11 + 12);
        r >>= params->gps.geom_angular_radius_inv_scale_log2;
      }

      // todo(df): handle the single laser case better
      Box3<int> sphBox{0, {r, twoPi, maxLaserIdx}};
      int refScale = params->gps.azimuth_scaling_enabled_flag
        ? params->attrSphericalMaxLog2
        : 0;
      auto attr_coord_scale = normalisedAxesWeights(sphBox, refScale);
      for (auto& aps : params->aps)
        if (aps.spherical_coord_flag)
          aps.attr_coord_scale = attr_coord_scale;
    }

    // Allocate storage for attribute contexts
    _ctxtMemAttrs.resize(params->sps.attributeSets.size());
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
  _firstSliceInFrame = true;

  // Configure output coud
  if (reconCloud) {
    reconCloud->setParametersFrom(*_sps, params->outputFpBits);
    reconCloud->frameNum = _frameCounter;
  }

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
  SrcMappedPointSet quantizedInput = quantization(inputPointCloud);

  // write out all parameter sets prior to encoding
  callback->onOutputBuffer(write(*_sps));
  callback->onOutputBuffer(write(*_sps, *_gps));
  for (const auto aps : _aps) {
    callback->onOutputBuffer(write(*_sps, *aps));
  }

  std::vector<std::vector<int32_t>> tileMaps;
  if (params->partition.tileSize) {
    tileMaps = tilePartition(params->partition, quantizedInput.cloud);

    // To tag the slice with the tile id there must be sufficient bits.
    // todo(df): determine sps parameter from the paritioning?
    assert(numBits(tileMaps.size() - 1) <= _sps->slice_tag_bits);

    // Default is to use implicit tile ids (ie list index)
    partitions.tileInventory.tile_id_bits = 0;

    // all tileOrigins are relative to the sequence bounding box
    partitions.tileInventory.origin = _sps->seqBoundingBoxOrigin;

    // Get the bounding box of current tile and write it into tileInventory
    partitions.tileInventory.tiles.resize(tileMaps.size());

    // Convert tile bounding boxes to sequence coordinate system.
    // A position in the box must remain in the box after conversion
    // irrispective of how the decoder outputs positions (fractional | integer)
    //   => truncate origin (eg, rounding 12.5 to 13 would not allow all
    //      decoders to find that point).
    //   => use next integer for upper coordinate.
    double gs = Rational(_sps->globalScale);

    for (int t = 0; t < tileMaps.size(); t++) {
      Box3<int32_t> bbox = quantizedInput.cloud.computeBoundingBox(
        tileMaps[t].begin(), tileMaps[t].end());

      auto& tileIvt = partitions.tileInventory.tiles[t];
      tileIvt.tile_id = t;
      for (int k = 0; k < 3; k++) {
        auto origin = std::trunc(bbox.min[k] * gs);
        auto size = std::ceil(bbox.max[k] * gs) - origin + 1;
        tileIvt.tileOrigin[k] = origin;
        tileIvt.tileSize[k] = size;
      }
    }
  } else {
    tileMaps.emplace_back();
    auto& tile = tileMaps.back();
    for (int i = 0; i < quantizedInput.cloud.getPointCount(); i++)
      tile.push_back(i);
  }

  if (partitions.tileInventory.tiles.size() > 1) {
    auto& inventory = partitions.tileInventory;
    assert(inventory.tiles.size() == tileMaps.size());
    std::cout << "Tile number: " << tileMaps.size() << std::endl;
    inventory.ti_seq_parameter_set_id = _sps->sps_seq_parameter_set_id;
    inventory.ti_origin_bits_minus1 =
      numBits(inventory.origin.abs().max()) - 1;

    // The inventory comes into force on the first frame
    inventory.ti_frame_ctr_bits = _sps->frame_ctr_bits;
    inventory.ti_frame_ctr = _frameCounter & ((1 << _sps->frame_ctr_bits) - 1);

    // Determine the number of bits for encoding tile sizes
    int maxValOrigin = 1;
    int maxValSize = 1;
    for (const auto& entry : inventory.tiles) {
      maxValOrigin = std::max(maxValOrigin, entry.tileOrigin.max());
      maxValSize = std::max(maxValSize, entry.tileSize.max() - 1);
    }
    inventory.tile_origin_bits_minus1 = numBits(maxValOrigin) - 1;
    inventory.tile_size_bits_minus1 = numBits(maxValSize) - 1;

    callback->onOutputBuffer(write(*_sps, partitions.tileInventory));
  }

  // Partition the input point cloud
  //  - get the partitial cloud of each tile
  //  - partitioning function produces a list of point indexes, origin and
  //    optional tile metadata for each partition.
  //  - encode any tile metadata
  //  NB: the partitioning method is required to ensure that the output
  //      slices conform to any codec limits.
  //  todo(df): consider requiring partitioning function to sort the input
  //            points and provide ranges rather than a set of indicies.
  do {
    for (int t = 0; t < tileMaps.size(); t++) {
      const auto& tile = tileMaps[t];
      auto tile_id = partitions.tileInventory.tiles.empty()
        ? 0
        : partitions.tileInventory.tiles[t].tile_id;

      // Get the point cloud of current tile and compute their bounding boxes
      PCCPointSet3 tileCloud = getPartition(quantizedInput.cloud, tile);
      Box3<int32_t> bbox = tileCloud.computeBoundingBox();

      // Move the tile cloud to coodinate origin
      // for the convenience of slice partitioning
      for (int i = 0; i < tileCloud.getPointCount(); i++)
        tileCloud[i] -= bbox.min;

      // don't partition if partitioning would result in a single slice.
      auto partitionMethod = params->partition.method;
      if (tileCloud.getPointCount() < params->partition.sliceMaxPoints)
        partitionMethod = PartitionMethod::kNone;

      // use the largest trisoup node size as a partitioning boundary for
      // consistency between slices with different trisoup node sizes.
      int partitionBoundaryLog2 = 0;
      if (!params->trisoupNodeSizesLog2.empty())
        partitionBoundaryLog2 = *std::max_element(
          params->trisoupNodeSizesLog2.begin(),
          params->trisoupNodeSizesLog2.end());

      //Slice partition of current tile
      std::vector<Partition> curSlices;
      switch (partitionMethod) {
      case PartitionMethod::kNone:
        curSlices = partitionNone(params->partition, tileCloud, tile_id);
        break;

      case PartitionMethod::kUniformGeom:
        curSlices = partitionByUniformGeom(
          params->partition, tileCloud, tile_id, partitionBoundaryLog2);
        break;

      case PartitionMethod::kUniformSquare:
        curSlices = partitionByUniformSquare(
          params->partition, tileCloud, tile_id, partitionBoundaryLog2);
        break;

      case PartitionMethod::kOctreeUniform:
        curSlices =
          partitionByOctreeDepth(params->partition, tileCloud, tile_id);
        break;

      case PartitionMethod::kNpoints:
        curSlices = partitionByNpts(params->partition, tileCloud, tile_id);
        break;
      }
      // Map slice indexes to tile indexes(the original indexes)
      for (int i = 0; i < curSlices.size(); i++) {
        for (int p = 0; p < curSlices[i].pointIndexes.size(); p++) {
          curSlices[i].pointIndexes[p] = tile[curSlices[i].pointIndexes[p]];
        }
      }

      partitions.slices.insert(
        partitions.slices.end(), curSlices.begin(), curSlices.end());
    }
    std::cout << "Slice number: " << partitions.slices.size() << std::endl;
  } while (0);

  // Encode each partition:
  //  - create a pointset comprising just the partitioned points
  //  - compress
  for (const auto& partition : partitions.slices) {
    // create partitioned point set
    PCCPointSet3 sliceCloud =
      getPartition(quantizedInput.cloud, partition.pointIndexes);

    PCCPointSet3 sliceSrcCloud =
      getPartition(inputPointCloud, quantizedInput, partition.pointIndexes);

    _sliceId = partition.sliceId;
    _tileId = partition.tileId;
    _sliceOrigin = sliceCloud.computeBoundingBox().min;
    compressPartition(sliceCloud, sliceSrcCloud, params, callback, reconCloud);
  }

  // Apply global scaling to reconstructed point cloud
  if (reconCloud)
    scaleGeometry(
      reconCloud->cloud, _sps->globalScale, reconCloud->outputFpBits);

  return 0;
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::deriveParameterSets(EncoderParams* params)
{
  // fixup extGeomScale in the case that we're coding metres
  if (params->sps.seq_geom_scale_unit_flag == ScaleUnit::kMetre)
    params->extGeomScale = 0;

  if (params->extGeomScale == 0.)
    params->extGeomScale = params->srcUnitLength;

  // Derive the sps scale factor:  The sequence scale is normalised to an
  // external geometry scale of 1.
  //  - Ie, if the user specifies that extGeomScale=2 (1 seq point is equal
  //    to 2 external points), seq_geom_scale is halved.
  //
  //  - In cases where the sequence scale is in metres, the external system
  //    is defined to have a unit length of 1 metre, and srcUnitLength must
  //    be used to define the sequence scale.
  //  - The user may define the relationship to the external coordinate system.
  //
  // NB: seq_geom_scale is the reciprocal of unit length
  params->sps.seqGeomScale = params->seqGeomScale / params->extGeomScale;

  // Global scaling converts from the coded scale to the sequence scale
  // NB: globalScale is constrained, eg 1.1 is not representable
  // todo: consider adjusting seqGeomScale to make a valid globalScale
  // todo: consider adjusting codedGeomScale to make a valid globalScale
  params->sps.globalScale =
    Rational(params->seqGeomScale / params->codedGeomScale);
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
  params->sps.profile.main_profile_compatibility_flag = 0;
  params->sps.profile.reserved_profile_compatibility_21bits = 0;
  params->sps.level = 0;

  // constraints
  params->sps.profile.unique_point_positions_constraint_flag = false;
  params->sps.profile.slice_reordering_constraint_flag =
    params->sps.entropy_continuation_enabled_flag;

  // use one bit to indicate frame boundaries
  params->sps.frame_ctr_bits = 1;

  // number of bits for slice tag (tileid) if tiles partitioning enabled
  // NB: the limit of 64 tiles is arbritrary
  params->sps.slice_tag_bits = params->partition.tileSize > 0 ? 6 : 0;

  // slice origin parameters used by this encoder implementation
  params->gps.geom_box_log2_scale_present_flag = true;
  params->gps.gps_geom_box_log2_scale = 0;

  // don't code per-slice angular origin
  params->gps.geom_slice_angular_origin_present_flag = false;

  // derive the idcm qp offset from cli
  params->gps.geom_idcm_qp_offset = params->idcmQp - params->gps.geom_base_qp;

  // Feature dependencies
  if (!params->gps.neighbour_avail_boundary_log2_minus1) {
    params->gps.adjacent_child_contextualization_enabled_flag = 0;
    params->gps.intra_pred_max_node_size_log2 = 0;
  }

  if (params->gps.predgeom_enabled_flag)
    params->gps.geom_planar_mode_enabled_flag = false;

  // fixup attribute parameters
  for (auto it : params->attributeIdxMap) {
    auto& attr_sps = params->sps.attributeSets[it.second];
    auto& attr_aps = params->aps[it.second];
    auto& attr_enc = params->attr[it.second];

    // this encoder does not (yet) support variable length attributes
    // todo(df): add variable length attribute support
    attr_aps.raw_attr_variable_len_flag = 0;

    // sanitise any intra prediction skipping
    if (attr_aps.attr_encoding != AttributeEncoding::kPredictingTransform)
      attr_aps.intra_lod_prediction_skip_layers = attr_aps.kSkipAllLayers;
    if (attr_aps.intra_lod_prediction_skip_layers < 0)
      attr_aps.intra_lod_prediction_skip_layers = attr_aps.kSkipAllLayers;

    // avoid signalling overly large values
    attr_aps.intra_lod_prediction_skip_layers = std::min(
      attr_aps.intra_lod_prediction_skip_layers,
      attr_aps.maxNumDetailLevels() + 1);

    // dist2 is refined in the slice header
    //  - the encoder always writes them unless syntatically prohibited:
    attr_aps.aps_slice_dist2_deltas_present_flag =
      attr_aps.lodParametersPresent()
      && !attr_aps.scalable_lifting_enabled_flag
      && attr_aps.num_detail_levels_minus1
      && attr_aps.lod_decimation_type != LodDecimationMethod::kPeriodic;

    // disable dist2 estimation when decimating with centroid sampler
    if (attr_aps.lod_decimation_type == LodDecimationMethod::kCentroid)
      attr_aps.aps_slice_dist2_deltas_present_flag = false;

    // If the lod search ranges are negative, use a full-range search
    // todo(df): lookup level limit
    if (attr_aps.inter_lod_search_range < 0)
      attr_aps.inter_lod_search_range = 1100000;

    if (attr_aps.intra_lod_search_range < 0)
      attr_aps.intra_lod_search_range = 1100000;

    // If all intra prediction layers are skipped, don't signal a search range
    if (
      attr_aps.intra_lod_prediction_skip_layers
      > attr_aps.maxNumDetailLevels())
      attr_aps.intra_lod_search_range = 0;

    // If there are no refinement layers, don't signal an inter search range
    if (attr_aps.maxNumDetailLevels() == 1)
      attr_aps.inter_lod_search_range = 0;

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
  CloudFrame* reconCloud)
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

  // apply a custom trisoup node size
  params->gbh.trisoup_node_size_log2_minus2 = 0;
  if (_gps->trisoup_enabled_flag) {
    int idx = std::min(_sliceId, int(params->trisoupNodeSizesLog2.size()) - 1);
    params->gbh.trisoup_node_size_log2_minus2 =
      params->trisoupNodeSizesLog2[idx] - 2;
  }

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
  if (_gps->geom_unique_points_flag || _gps->trisoup_enabled_flag) {
    for (const auto& attr_sps : _sps->attributeSets) {
      recolour(
        attr_sps, params->recolour, originPartCloud, _srcToCodingScale,
        _originInCodingCoords + _sliceOrigin, &pointCloud);
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

    // NB: regionQpOrigin/regionQpSize use the STV axes, not XYZ.
    if (false) {
      abh.qpRegions.emplace_back();
      auto& region = abh.qpRegions.back();
      region.regionOrigin = 0;
      region.regionSize = 0;
      region.attr_region_qp_offset = {0, 0};
      abh.attr_region_bits_minus1 = -1
        + numBits(
            std::max(region.regionOrigin.max(), region.regionSize.max()));
    }
    // Number of regions is constrained to at most 1.
    assert(abh.qpRegions.size() <= 1);

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
        altPositions.resize(pointCloud.getPointCount());

        auto laserOrigin = _gbh.geomAngularOrigin(*_gps);
        bboxRpl = convertXyzToRpl(
          laserOrigin, _gps->angularTheta.data(), _gps->angularTheta.size(),
          &pointCloud[0], &pointCloud[0] + pointCloud.getPointCount(),
          altPositions.data());
      }

      offsetAndScale(
        bboxRpl.min, attr_aps.attr_coord_scale, altPositions.data(),
        altPositions.data() + altPositions.size());

      pointCloud.swapPoints(altPositions);
    }

    // calculate dist2 for this slice
    abh.attr_dist2_delta = 0;
    if (attr_aps.aps_slice_dist2_deltas_present_flag) {
      // todo(df): this could be set in the sps and refined only if necessary
      auto dist2 =
        estimateDist2(pointCloud, 100, 128, attr_enc.dist2PercentileEstimate);
      abh.attr_dist2_delta = dist2 - attr_aps.dist2;
    }

    // replace the attribute encoder if not compatible
    if (!attrEncoder->isReusable(attr_aps, abh))
      attrEncoder = makeAttributeEncoder();

    auto& ctxtMemAttr = _ctxtMemAttrs.at(abh.attr_sps_attr_idx);
    attrEncoder->encode(
      *_sps, attr_sps, attr_aps, abh, ctxtMemAttr, pointCloud, &payload);

    if (attr_aps.spherical_coord_flag)
      pointCloud.swapPoints(altPositions);

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

  // Note the current slice id for loss detection with entropy continuation
  _prevSliceId = _sliceId;

  // prevent re-use of this sliceId:  the next slice (geometry + attributes)
  // should be distinguishable from the current slice.
  _sliceId++;
  _firstSliceInFrame = false;

  if (reconCloud)
    appendSlice(reconCloud->cloud);
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::encodeGeometryBrick(
  const EncoderParams* params, PayloadBuffer* buf)
{
  GeometryBrickHeader gbh;
  gbh.geom_geom_parameter_set_id = _gps->gps_geom_parameter_set_id;
  gbh.geom_slice_id = _sliceId;
  gbh.prev_slice_id = _prevSliceId;
  // NB: slice_tag could be set to some other (external) meaningful value
  gbh.slice_tag = std::max(0, _tileId);
  gbh.frame_ctr_lsb = _frameCounter & ((1 << _sps->frame_ctr_bits) - 1);
  gbh.geomBoxOrigin = _sliceOrigin;
  gbh.gbhAngularOrigin = _gps->gpsAngularOrigin - _sliceOrigin;
  gbh.geom_box_origin_bits_minus1 = numBits(gbh.geomBoxOrigin.max()) - 1;
  gbh.geom_box_log2_scale = 0;
  gbh.geom_slice_qp_offset = params->gbh.geom_slice_qp_offset;
  gbh.geom_stream_cnt_minus1 = params->gbh.geom_stream_cnt_minus1;
  gbh.trisoup_node_size_log2_minus2 =
    params->gbh.trisoup_node_size_log2_minus2;

  gbh.geom_qp_offset_intvl_log2_delta =
    params->gbh.geom_qp_offset_intvl_log2_delta;

  // Entropy continuation is not permitted in the first slice of a frame
  gbh.entropy_continuation_flag = false;
  if (_sps->entropy_continuation_enabled_flag)
    gbh.entropy_continuation_flag = !_firstSliceInFrame;

  // inform the geometry coder what the root node size is
  for (int k = 0; k < 3; k++) {
    // NB: A minimum whd of 2 means there is always at least 1 tree level
    gbh.rootNodeSizeLog2[k] = ceillog2(std::max(2, _sliceBoxWhd[k]));

    // The root node size cannot be smaller than the trisoup node size
    // since this is how the root node size is defined at the decoder.
    // NB: the following isn't strictly necessary, but avoids accidents
    // involving the qtbt derivation.
    gbh.rootNodeSizeLog2[k] =
      std::max(gbh.trisoupNodeSizeLog2(*_gps), gbh.rootNodeSizeLog2[k]);
  }
  gbh.maxRootNodeDimLog2 = gbh.rootNodeSizeLog2.max();

  // use a cubic node if qtbt is disabled
  if (!_gps->predgeom_enabled_flag && !_gps->qtbt_enabled_flag)
    gbh.rootNodeSizeLog2 = gbh.maxRootNodeDimLog2;

  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = int(pointCloud.getPointCount()) * 3 * 4 + 1024;

  // allocate entropy streams
  std::vector<std::unique_ptr<EntropyEncoder>> arithmeticEncoders;
  for (int i = 0; i < 1 + gbh.geom_stream_cnt_minus1; i++) {
    arithmeticEncoders.emplace_back(new EntropyEncoder(maxAcBufLen, nullptr));
    auto& aec = arithmeticEncoders.back();
    aec->enableBypassStream(_sps->cabac_bypass_stream_enabled_flag);
    aec->start();
  }

  // forget (reset) all saved context state at boundary
  if (!gbh.entropy_continuation_flag) {
    _ctxtMemOctreeGeom->reset();
    _ctxtMemPredGeom->reset();
    for (auto& ctxtMem : _ctxtMemAttrs)
      ctxtMem.reset();
  }

  if (_gps->predgeom_enabled_flag)
    encodePredictiveGeometry(
      params->predGeom, *_gps, gbh, pointCloud, &_posSph, *_ctxtMemPredGeom,
      arithmeticEncoders[0].get());
  else if (!_gps->trisoup_enabled_flag)
    encodeGeometryOctree(
      params->geom, *_gps, gbh, pointCloud, *_ctxtMemOctreeGeom,
      arithmeticEncoders);
  else {
    // limit the number of points to the slice limit
    // todo(df): this should be derived from the level
    gbh.footer.geom_num_points_minus1 = params->partition.sliceMaxPoints - 1;
    encodeGeometryTrisoup(
      params->geom, *_gps, gbh, pointCloud, *_ctxtMemOctreeGeom,
      arithmeticEncoders);
  }

  // signal the actual number of points coded
  gbh.footer.geom_num_points_minus1 = pointCloud.getPointCount() - 1;

  // assemble data unit
  //  - record the position of each aec buffer for chunk concatenation
  std::vector<std::pair<size_t, size_t>> aecStreams;
  write(*_sps, *_gps, gbh, buf);
  for (auto& arithmeticEncoder : arithmeticEncoders) {
    auto aecLen = arithmeticEncoder->stop();
    auto aecBuf = arithmeticEncoder->buffer();
    aecStreams.emplace_back(buf->size(), aecLen);
    buf->insert(buf->end(), aecBuf, aecBuf + aecLen);
  }

  // This process is performed here from the last chunk to the first.  It
  // is also possible to implement this in a forwards direction too.
  if (_sps->cabac_bypass_stream_enabled_flag) {
    aecStreams.pop_back();
    for (auto i = aecStreams.size() - 1; i + 1; i--) {
      auto& stream = aecStreams[i];
      auto* ptr = reinterpret_cast<uint8_t*>(buf->data());
      auto* chunkA = ptr + stream.first + (stream.second & ~0xff);
      auto* chunkB = ptr + stream.first + stream.second;
      auto* end = ptr + buf->size();
      ChunkStreamBuilder::spliceChunkStreams(chunkA, chunkB, end);
    }
  }

  // append the footer
  write(*_gps, gbh, gbh.footer, buf);

  // Cache gbh for later reference
  _gbh = gbh;
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::appendSlice(PCCPointSet3& accumCloud)
{
  // offset current point cloud to be in coding coordinate system
  size_t numPoints = pointCloud.getPointCount();
  for (size_t i = 0; i < numPoints; i++)
    for (int k = 0; k < 3; k++)
      pointCloud[i][k] += _sliceOrigin[k];

  accumCloud.append(pointCloud);
}

//----------------------------------------------------------------------------
// translates and scales inputPointCloud, storing the result in
// this->pointCloud for use by the encoding process.

SrcMappedPointSet
PCCTMC3Encoder3::quantization(const PCCPointSet3& src)
{
  // Currently the sequence bounding box size must be set
  assert(_sps->seqBoundingBoxSize != Vec3<int>{0});

  // Clamp all points to [clampBox.min, clampBox.max] after translation
  // and quantisation.
  Box3<int32_t> clampBox(0, std::numeric_limits<int32_t>::max());

  // When using predictive geometry, sub-sample the point cloud and let
  // the predictive geometry coder quantise internally.
  if (_inputDecimationScale != 1.)
    return samplePositionsUniq(
      _inputDecimationScale, _srcToCodingScale, _originInCodingCoords, src);

  if (_gps->geom_unique_points_flag)
    return quantizePositionsUniq(
      _srcToCodingScale, _originInCodingCoords, clampBox, src);

  SrcMappedPointSet dst;
  quantizePositions(
    _srcToCodingScale, _originInCodingCoords, clampBox, src, &dst.cloud);
  return dst;
}

//----------------------------------------------------------------------------
// get the partial point cloud according to required point indexes

PCCPointSet3
getPartition(const PCCPointSet3& src, const std::vector<int32_t>& indexes)
{
  PCCPointSet3 dst;
  dst.addRemoveAttributes(src);

  int partitionSize = indexes.size();
  dst.resize(partitionSize);

  for (int i = 0; i < partitionSize; i++) {
    int inputIdx = indexes[i];
    dst[i] = src[inputIdx];

    if (src.hasColors())
      dst.setColor(i, src.getColor(inputIdx));

    if (src.hasReflectances())
      dst.setReflectance(i, src.getReflectance(inputIdx));

    if (src.hasLaserAngles())
      dst.setLaserAngle(i, src.getLaserAngle(inputIdx));
  }

  return dst;
}

//----------------------------------------------------------------------------
// get the partial point cloud according to required point indexes

PCCPointSet3
getPartition(
  const PCCPointSet3& src,
  const SrcMappedPointSet& map,
  const std::vector<int32_t>& indexes)
{
  // Without the list, do nothing
  if (map.idxToSrcIdx.empty())
    return {};

  // work out the destination size.
  // loop over each linked list until an element points to itself
  int size = 0;
  for (int idx : indexes) {
    int prevIdx, srcIdx = map.idxToSrcIdx[idx];
    do {
      size++;
      prevIdx = srcIdx;
      srcIdx = map.srcIdxDupList[srcIdx];
    } while (srcIdx != prevIdx);
  }

  PCCPointSet3 dst;
  dst.addRemoveAttributes(src);
  dst.resize(size);

  int dstIdx = 0;
  for (int idx : indexes) {
    int prevIdx, srcIdx = map.idxToSrcIdx[idx];
    do {
      dst[dstIdx] = src[srcIdx];

      if (src.hasColors())
        dst.setColor(dstIdx, src.getColor(srcIdx));

      if (src.hasReflectances())
        dst.setReflectance(dstIdx, src.getReflectance(srcIdx));

      if (src.hasLaserAngles())
        dst.setLaserAngle(dstIdx, src.getLaserAngle(srcIdx));

      dstIdx++;
      prevIdx = srcIdx;
      srcIdx = map.srcIdxDupList[srcIdx];
    } while (srcIdx != prevIdx);
  }

  return dst;
}

//============================================================================

}  // namespace pcc
