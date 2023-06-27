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
#include <numeric>

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
#include "TMC3.h"
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
  if (params->gps.biPredictionEnabledFlag)
    _frameCounter = biPredEncodeParams.currentFrameIndex;
  else
    _frameCounter++;

  if (_frameCounter == 0) {
    // Angular predictive geometry coding needs to determine spherical
    // positions.  To avoid quantization of the input disturbing this:
    //  - sequence scaling is replaced by decimation of the input
    //  - any user-specified global scaling is honoured
    _inputDecimationScale = 1.;
    if (
      params->gps.predgeom_enabled_flag
      && params->gps.geom_angular_mode_enabled_flag) {
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

    // Determine point number alone one laser beam
    params->gps.one_point_alone_laser_beam_flag = false;
    if (params->gps.geom_angular_mode_enabled_flag&&params->gps.geom_inter_idcm_enabled_flag) {
      int maxPointsNumPerTurn = 0;
      for (int i = 0; i < params->gps.angularNumPhiPerTurn.size(); i++) {
        maxPointsNumPerTurn += params->gps.angularNumPhiPerTurn[i];
      }
      if (inputPointCloud.getPointCount() / double(maxPointsNumPerTurn) < 2) {
        params->gps.one_point_alone_laser_beam_flag = true;
      }
      else {
        params->gps.one_point_alone_laser_beam_flag = false;
      }
    }
    else {
      params->gps.one_point_alone_laser_beam_flag = false;
    }
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

    if (params->gps.globalMotionEnabled) {
      if (params->gps.predgeom_enabled_flag) {
        _refFrameSph.parseMotionParams(motionVectorFileName, 1.0);
        if (params->gps.biPredictionEnabledFlag)
          biPredEncodeParams._refFrameSph2.parseMotionParams(motionVectorFileName, 1.0);
        params->interGeom.motionParams.parseFile(
          motionVectorFileName, params->codedGeomScale);
        deriveMotionParams(params);
      } else if (!params->gps.trisoup_enabled_flag) {
        params->interGeom.motionParams.parseFile(
          motionVectorFileName, params->codedGeomScale);
        deriveMotionParams(params);
      }
    }
  }
  if (params->interGeom.deriveGMThreshold && params->gps.geom_angular_mode_enabled_flag) {

    const auto pointCount = inputPointCloud.getPointCount();
    const auto bbox = inputPointCloud.computeBoundingBox();
    const int histScale = params->interGeom.gmThresholdHistScale;

    const int minZ = std::max(bbox.min.z(), params->interGeom.gmThresholdMinZ);
    const int maxZ = std::min(bbox.max.z(), params->interGeom.gmThresholdMaxZ);

    const int histRange = (maxZ - minZ) / histScale;
    std::vector<int> histInp(histRange), midVals(histRange);
    std::vector<float> probs(histRange);

    for (auto i = 0; i < pointCount; i++) {
      auto scaledZ = (inputPointCloud[i].z() - minZ) / histScale;
      if (scaledZ >= 0 && scaledZ < histRange)
        histInp[scaledZ]++;
    }

    for (size_t i = 0; i < histRange; ++i)
      midVals[i] = minZ + int((i + 0.5) * histScale);

    const float sum = std::accumulate(histInp.begin(), histInp.end(), 0.);

    for (size_t i = 0; i < histRange; ++i)
      probs[i] = histInp[i] / sum;

    const float mean =
      std::inner_product(probs.begin(), probs.end(), midVals.begin(), 0.);

    float variance = 0.0;
    for (size_t i = 0; i < histRange; ++i)
      variance += probs[i] * pow((midVals[i] - mean), 2);

    const int topIdx =
      std::max_element(histInp.begin(), histInp.end()) - histInp.begin();
    const int top = midVals[topIdx];

    const int deltaLeft =
      params->interGeom.gmThresholdLeftScale * std::sqrt(variance);
    const int deltaRight =
      params->interGeom.gmThresholdRightScale * std::sqrt(variance);
    const int leftOffset =
      std::max(minZ, top - deltaLeft) * params->codedGeomScale;
    const int rightOffset =
      std::min(maxZ, top + deltaRight) * params->codedGeomScale;

    if (params->gps.globalMotionEnabled) {
      if (params->gps.predgeom_enabled_flag) {
        _refFrameSph.updateThresholds(_frameCounter, leftOffset, rightOffset);
        if (params->gps.biPredictionEnabledFlag) {
          biPredEncodeParams._refFrameSph2.updateThresholds(_frameCounter, leftOffset, rightOffset);
        }
      }
      else if (!params->gps.trisoup_enabled_flag)
        params->interGeom.motionParams.updateThresholds(
          _frameCounter, leftOffset, rightOffset);
    } else {
      params->predGeom.splitter = rightOffset - bbox.min.z(); 
      params->predGeom.splitter =
        int(params->predGeom.splitter * _inputDecimationScale);
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

  // use the largest trisoup node size as a partitioning boundary for
  // consistency between slices with different trisoup node sizes.
  int partitionBoundaryLog2 = 0;
  if (!params->trisoupNodeSizesLog2.empty())
    partitionBoundaryLog2 = *std::max_element(
      params->trisoupNodeSizesLog2.begin(),
      params->trisoupNodeSizesLog2.end());

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
      if (!params->trisoupNodeSizesLog2.empty()) {
        for (int i = 0; i < curSlices.size(); i++) {
          for (int p = 0; p < curSlices[i].pointIndexesPadding.size(); p++) {
            curSlices[i].pointIndexesPadding[p] =
              tile[curSlices[i].pointIndexesPadding[p]];
          }
        }

        for (int i = 0; i < curSlices.size(); i++) {
          for (int p = 0; p < curSlices[i].pointIndexesPadding2.size(); p++) {
            curSlices[i].pointIndexesPadding2[p] =
              tile[curSlices[i].pointIndexesPadding2[p]];
          }
        }
      }

      partitions.slices.insert(
        partitions.slices.end(), curSlices.begin(), curSlices.end());
    }
    std::cout << "Slice number: " << partitions.slices.size() << std::endl;
  } while (0);

  if (_frameCounter) {
    if (params->gps.predgeom_enabled_flag) {
      if (params->gps.biPredictionEnabledFlag){
        if (!biPredEncodeParams.codeCurrentFrameAsBFrame){
          _refFrameSph.clearRefFrameCur();
          _refFrameSph.insert(biPredEncodeParams._refPosSph2);
          _refPosSph = biPredEncodeParams._refPosSph2;
        }
        _refFrameSph.setFrameCtr(_frameCounter - 2);    
        _refFrameSph.setRefFrameCtr(biPredEncodeParams.refFrameIndex - 1);
        if(previousAsInter)
          _refFrameSph.updateCurMovingStatue(_frameCounter, biPredEncodeParams.refFrameIndex);

        if (biPredEncodeParams.codeCurrentFrameAsBFrame){
          biPredEncodeParams._refFrameSph2.setFrameCtr(_frameCounter - 2);
          biPredEncodeParams._refFrameSph2.setRefFrameCtr(biPredEncodeParams.refFrameIndex2 - 1);
          if(previousAsInter)
            biPredEncodeParams._refFrameSph2.updateCurMovingStatue(_frameCounter, biPredEncodeParams.refFrameIndex2);
        } 
      }
      _refFrameSph.updateFrame(*_gps);
      if (params->gps.biPredictionEnabledFlag && biPredEncodeParams.codeCurrentFrameAsBFrame)
        biPredEncodeParams._refFrameSph2.updateFrame(*_gps);
    } else if (!params->gps.trisoup_enabled_flag) {
      params->interGeom.motionParams.AdvanceFrame();
    }
    if (
      params->gps.biPredictionEnabledFlag
      && !biPredEncodeParams.codeCurrentFrameAsBFrame) {
      attrInterPredParams.referencePointCloud =
        biPredEncodeParams.attrInterPredParams2.referencePointCloud;
      _refFrame.cloud = biPredEncodeParams.predPointCloud2;
    }
  } else {
    _refFrameSph.setGlobalMotionEnabled(_gps->globalMotionEnabled);
    biPredEncodeParams._refFrameSph2.setGlobalMotionEnabled(_gps->globalMotionEnabled);
  }

  // Encode each partition:
  //  - create a pointset comprising just the partitioned points
  //  - compress
  pcc::CloudFrame reconCloudAlt;

  for (const auto& partition : partitions.slices) {
    // create partitioned point set
    PCCPointSet3 sliceCloud =
      getPartition(quantizedInput.cloud, partition.pointIndexes);

    PCCPointSet3 sliceCloudPadding;
    if (!params->trisoupNodeSizesLog2.empty()) {
      sliceCloudPadding =
        getPartition(quantizedInput.cloud, partition.pointIndexesPadding);

      PCCPointSet3 sliceCloudPadding2 =
        getPartition(quantizedInput.cloud, partition.pointIndexesPadding2);

      sliceCloudPadding.append(sliceCloudPadding2);
    }

    PCCPointSet3 sliceSrcCloud =
      getPartition(inputPointCloud, quantizedInput, partition.pointIndexes);

    _sliceId = partition.sliceId;
    _tileId = partition.tileId;
    _sliceOrigin = sliceCloud.computeBoundingBox().min;

    compressPartition(
      sliceCloud, sliceSrcCloud, sliceCloudPadding, params, callback,
      reconCloud, &reconCloudAlt);
  }

  const PCCPointSet3 accurrentPointCloud =(reconCloud) ? reconCloud->cloud : PCCPointSet3{};
  if (_sps->inter_frame_prediction_enabled_flag) {
    // buffer the current frame for potential use in predictionif (_gps->biPredictionEnabledFlag){
    if (_gps->biPredictionEnabledFlag) {
      if (!_gps->predgeom_enabled_flag) {
        if (!biPredEncodeParams.codeCurrentFrameAsBFrame) {
          biPredEncodeParams.predPointCloud2 = accurrentPointCloud;
        } else {
          _refFrame = *reconCloud;
          _refFrameAlt = reconCloudAlt;
        }
      }
    } else {
      _refFrame = *reconCloud;
      _refFrameAlt = reconCloudAlt;
    }
  }

  biPredEncodeParams.boundingBoxSize =
    accurrentPointCloud.computeBoundingBox().max
    - accurrentPointCloud.computeBoundingBox().min;

   if (_gps->geom_z_compensation_enabled_flag && reconCloud) {
    auto plyScale = reciprocal(_sps->seqGeomScale);
    plyScale.numerator *= 1000;
    auto laserOrigin = _gps->gpsAngularOrigin;
    compensateZCoordinate(
      reconCloud->cloud, _gps, plyScale, laserOrigin,
      reconCloud->outputUnitLength, reconCloud->outputOrigin);
  }

  // Apply global scaling to reconstructed point cloud
  if (reconCloud)
    scaleGeometry(
      reconCloud->cloud, _sps->globalScale, reconCloud->outputFpBits);

  return 0;
}

//----------------------------------------------------------------------------

int
PCCTMC3Encoder3::compressHGOF(
  const PCCPointSet3& inputPointCloud,
  EncoderParams* params,
  PCCTMC3Encoder3::Callbacks* callback,
  CloudFrame* reconCloud)
{
  auto& gof = hGOFEncodeParams.gof;
  auto& gofSpherical = hGOFEncodeParams.gof_spherical;
  auto& gofPosSph = hGOFEncodeParams.gof_posSph;
  auto& refPointCloud1 = attrInterPredParams.referencePointCloud;
  auto& refPointCloud2 = biPredEncodeParams.attrInterPredParams2.referencePointCloud;
  auto& predPointCloud1 = _refFrame.cloud;
  auto& predPointCloud2 = biPredEncodeParams.predPointCloud2;
  auto& refPosSph1 = _refPosSph;
  auto& refPosSph2 = biPredEncodeParams._refPosSph2;
  auto& refFrameSph1 = _refFrameSph;
  auto& refFrameSph2 = biPredEncodeParams._refFrameSph2;

  if (!biPredEncodeParams.codeCurrentFrameAsBFrame) {
    if (gof.size()) {
      predPointCloud2 = gof[gof.size() - 1];
      refPointCloud2 = gofSpherical[gof.size() - 1];
      refPosSph2 = gofPosSph[gof.size() - 1];
      hGOFEncodeParams.clearGofs();
    }
  } else {
    if (!gof.size()) {
      hGOFEncodeParams.initializeGof(
        biPredEncodeParams.refTimesList.size(), predPointCloud1,
        predPointCloud2, refPointCloud1, refPointCloud2, refPosSph1, refPosSph2);
    }

    const auto idx1 = hGOFEncodeParams.currFrameIndexInGOF
      + biPredEncodeParams.refFrameIndex
      - biPredEncodeParams.currentFrameIndex;
    const auto idx2 = hGOFEncodeParams.currFrameIndexInGOF
      + biPredEncodeParams.refFrameIndex2
      - biPredEncodeParams.currentFrameIndex;
    hGOFEncodeParams.updateReferenceFrames(
      predPointCloud1, predPointCloud2, refPointCloud1, refPointCloud2, refPosSph1, refPosSph2,
      refFrameSph1, refFrameSph2, idx1,idx2);

    if (!(--(biPredEncodeParams.refTimesList[idx1]))) {
      hGOFEncodeParams.clearFrame(idx1);
    }

    if (!(--(biPredEncodeParams.refTimesList[idx2]))) {
      hGOFEncodeParams.clearFrame(idx2);
    }
  }

  int ret = compress(inputPointCloud, params, callback, reconCloud);
  if (biPredEncodeParams.codeCurrentFrameAsBFrame) {
    hGOFEncodeParams.storeReferenceFrame(hGOFEncodeParams.currFrameIndexInGOF, predPointCloud1, refPointCloud1, refPosSph1);
  }

  return ret;
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
    params->sps.entropy_continuation_enabled_flag
    || params->sps.inter_entropy_continuation_enabled_flag;

  // use one bit to indicate frame boundaries
  params->sps.frame_ctr_bits = 1;
  if (params->gps.biPredictionEnabledFlag) {
    int bits = 1;
    int predictionPeriod = params->gps.biPredictionPeriod;
    while (predictionPeriod >> bits) {
      bits++;
    }
    params->sps.frame_ctr_bits = bits + 1;
  } 

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

    if (attr_aps.rahtPredParams.raht_prediction_search_range < 0)
      attr_aps.rahtPredParams.raht_prediction_search_range = 1100000;

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


//---------------------------------------------------------------------------
// motion parameter derivation
// Setup the block sizes b
void
PCCTMC3Encoder3::deriveMotionParams(EncoderParams* params)
{
  auto scaleFactor = params->codedGeomScale;
  params->interGeom.th_dist = int(1000 * scaleFactor);
  params->gbh.motion_block_size = {0, 0, 0};
  for (int i = 0; i < 3; i++) {
    if (params->interGeom.motion_block_size[i] > 0)
      params->gbh.motion_block_size[i] = std::max(
        64,
        int(std::round(params->interGeom.motion_block_size[i] * scaleFactor)));
    else
      params->gbh.motion_block_size[i] = 0;
  }
  params->interGeom.motion_window_size = std::max(
    2, int(std::round(params->interGeom.motion_window_size * scaleFactor)));
}
//---------------------------------------------------------------------------
// set the bi-prediction params of one frame
void
PCCTMC3Encoder3::setBiPredEncodeParams(
  const bool codeCurrFrameAsBFrame,
  const int currFrameIndex,
  const int refFrameIndex,
  const int refFrameIndex2,
  const int qpShift)
{
  biPredEncodeParams.codeCurrentFrameAsBFrame = codeCurrFrameAsBFrame;
  biPredEncodeParams.currentFrameIndex = currFrameIndex;
  biPredEncodeParams.refFrameIndex = refFrameIndex;
  biPredEncodeParams.refFrameIndex2 = refFrameIndex2;
  biPredEncodeParams.qpShiftTimes = qpShift;
}

//---------------------------------------------------------------------------
// initalization the bi-prediction params of one GOF when the hierarchical GOF structure is applied
void
PCCTMC3Encoder3::initBiPredEncodeParamsGOF(const int currPredPeriod)
{
  hGOFEncodeParams.clearLists();
  hGOFEncodeParams.refTimesList.resize(currPredPeriod + 1, 0);
  hGOFEncodeParams.refTimesList[currPredPeriod]++;
  hGOFEncodeParams.GenerateList(0, currPredPeriod, 0, 1, 1);
}

//----------------------------------------------------------------------------

//---------------------------------------------------------------------------
//Global motion constrain for B-frame structure
bool
PCCTMC3Encoder3::biPredictionEligibility(
  const int curPicIndex, const int prePicIndex, EncoderParams* params)
{
  if (params->gps.globalMotionEnabled) {
    std::vector<std::vector<double>> gm = {
      {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
    if (!params->gps.trisoup_enabled_flag)
      params->interGeom.motionParams.getMotionParamsMultiple(
        gm, curPicIndex, prePicIndex);

    const auto& bBoxSize = biPredEncodeParams.boundingBoxSize;
    const double thr1_tan_pre_gof = tan(M_PI * 0.05 / 180);
    const double thr1_sin_pre_gof = sin(M_PI * 0.05 / 180);
    const double thr2_pre_gof = 0.001;
    const double SxPercent = abs(gm[3][0]) / double(bBoxSize[0]);
    const double SyPercent = abs(gm[3][1]) / double(bBoxSize[1]);
    const double SzPercent = abs(gm[3][2]) / double(bBoxSize[2]);
    const double Rx = abs(gm[1][2] / gm[2][2]);
    const double Ry = abs(gm[0][2]);
    const double Rz = abs(gm[0][1] / gm[0][0]);

    return SxPercent < thr2_pre_gof && SyPercent < thr2_pre_gof
      && SzPercent < thr2_pre_gof && Rx < thr1_tan_pre_gof
      && Ry < thr1_sin_pre_gof && Rz < thr1_tan_pre_gof;
  } else {
    return true;
  }
}
//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::compressPartition(
  const PCCPointSet3& inputPointCloud,
  const PCCPointSet3& originPartCloud,
  const PCCPointSet3& paddingPointCloud,
  EncoderParams* params,
  PCCTMC3Encoder3::Callbacks* callback,
  CloudFrame* reconCloud,
  CloudFrame* reconCloudAlt)
{
  // geometry compression consists of the following stages:
  //  - prefilter/quantize geometry (non-normative)
  //  - encode geometry (single slice, id = 0)
  //  - recolour
  pointCloud.clear();
  pointCloud = inputPointCloud;
  if (paddingPointCloud.getPointCount() != 0)
    pointCloudPadding = paddingPointCloud;

  // apply a custom trisoup node size
  params->gbh.trisoup_node_size_log2_minus2 = 0;
  if (_gps->trisoup_enabled_flag) {
    int idx = std::min(_sliceId, int(params->trisoupNodeSizesLog2.size()) - 1);
    params->gbh.trisoup_node_size_log2_minus2 =
      params->trisoupNodeSizesLog2[idx] - 2;
  }

  // Offset the point cloud to account for (preset) _sliceOrigin.
  // The new maximum bounds of the offset cloud
  const size_t pointCount = pointCloud.getPointCount();
  for (size_t i = 0; i < pointCount; ++i) {
    pointCloud[i] -= _sliceOrigin;
  }
  if (pointCloudPadding.getPointCount() != 0) {
    const size_t pointCountPadding = pointCloudPadding.getPointCount();
    for (size_t i = 0; i < pointCountPadding; ++i) {
      const point_t point = (pointCloudPadding[i] -= _sliceOrigin);
    }
  }

  // the encoded coordinate system is non-negative.
  Box3<int32_t> bbox = pointCloud.computeBoundingBox();

  if ( _gps->trisoup_enabled_flag ){
    params->gbh.slice_bb_pos   = 0;
    params->gbh.slice_bb_width = 0;
    int mask = ( 1 << params->gbh.trisoupNodeSizeLog2(params->gps) ) - 1;
    params->gbh.slice_bb_pos_bits = 0;
    params->gbh.slice_bb_pos_log2_scale = 0;
    params->gbh.slice_bb_width_bits = 0;
    params->gbh.slice_bb_width_log2_scale = 0;
    if( params->gps.non_cubic_node_start_edge ) {
      params->gbh.slice_bb_pos   = bbox.min;
      if( ( bbox.min[0] & mask ) ||
          ( bbox.min[1] & mask ) ||
          ( bbox.min[2] & mask ) ) {
        params->gbh.slice_bb_pos_bits = numBits( params->gbh.slice_bb_pos.max() );
        params->gbh.slice_bb_pos_log2_scale = 0;
      }
    }
    if( params->gps.non_cubic_node_end_edge ) {
      params->gbh.slice_bb_width = bbox.max - params->gbh.slice_bb_pos;
      if( ( bbox.max[0] & mask ) ||
          ( bbox.max[1] & mask ) ||
          ( bbox.max[2] & mask ) ) {
        params->gbh.slice_bb_width_bits = numBits( params->gbh.slice_bb_width.max() );
        params->gbh.slice_bb_width_log2_scale = 0;
      }
    }
  }

  // todo(df): don't update maxBound if something is forcing the value?
  // NB: size is max - min + 1
  _sliceBoxWhd = bbox.max + 1;

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

  PCCPointSet3 reconSliceAltPositions;

  bool currFrameNotCodedAsB =
    (_gps->biPredictionEnabledFlag
      && !biPredEncodeParams.codeCurrentFrameAsBFrame);
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
    abh.attr_raht_ac_coeff_qp_delta_luma = {};
    abh.attr_raht_ac_coeff_qp_delta_chroma = {};

    if(_gbh.interPredictionEnabledFlag)
      abh.attr_qp_delta_luma = attr_aps.qpShiftStep;

    if (_gps->biPredictionEnabledFlag)
      abh.attr_qp_delta_luma *= biPredEncodeParams.qpShiftTimes;

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
    
    abh.enableAttrInterPred =
      (attr_aps.attr_encoding == AttributeEncoding::kRAHTransform)
      ? _codeCurrFrameAsInter
      : movingState;
    attrInterPredParams.enableAttrInterPred =
      attr_aps.attrInterPredictionEnabled && abh.enableAttrInterPred;
 

    abh.attrInterPredSearchRange = attr_aps.attrInterPredSearchRange;
    abh.disableAttrInterPredForRefFrame2 = !biPredEncodeParams.movingState2;
    biPredEncodeParams.attrInterPredParams2.enableAttrInterPred =
      _gps->biPredictionEnabledFlag && attr_aps.attrInterPredictionEnabled
      && !abh.disableAttrInterPredForRefFrame2;

    if (_gps->biPredictionEnabledFlag) {
      //FrameMerge for attribute inter prediction
      if (
        biPredEncodeParams.codeCurrentFrameAsBFrame
        && biPredEncodeParams.attrInterPredParams2.enableAttrInterPred) {
        if (attrInterPredParams.enableAttrInterPred) {
          attrInterPredParams.referencePointCloud.append(
            biPredEncodeParams.attrInterPredParams2.referencePointCloud);
          //The times of neighbor search process is halfed
          abh.attrInterPredSearchRange /= 2;
        } else {
          attrInterPredParams = biPredEncodeParams.attrInterPredParams2;
          attrInterPredParams.frameDistance = 1;
        }
      }
    }

    attrInterPredParams.paramsForInterRAHT.raht_inter_prediction_depth_minus1 =
      attr_aps.raht_inter_prediction_depth_minus1;
    attrInterPredParams.paramsForInterRAHT.raht_inter_prediction_enabled =
      attr_aps.attrInterPredictionEnabled;
    attrInterPredParams.paramsForInterRAHT.enableFilterEstimation = attr_aps.raht_send_inter_filters;
    attrInterPredParams.paramsForInterRAHT.skipInitLayersForFiltering = attr_aps.raht_inter_skip_layers;

    attrInterPredParams.paramsForInterRAHT.raht_enable_inter_intra_layer_RDO =
      attr_aps.raht_enable_code_layer;
    attrInterPredParams.attr_layer_code_mode.clear();
    if (attr_sps.attr_num_dimensions_minus1 == 0 && !_gps->geom_angular_mode_enabled_flag)
      attrInterPredParams.enableSkipCode = false;
    else
      attrInterPredParams.enableSkipCode = true;

    attrInterPredParams.attrInterIntraSliceRDO =
      attr_enc.attrInterIntraSliceRDO && attr_aps.attrInterPredictionEnabled;




    // Convert cartesian positions to spherical for use in attribute coding.
    // NB: this retains the original cartesian positions to restore afterwards
    std::vector<pcc::point_t> altPositions;
    if (attr_aps.spherical_coord_flag) {
      // If predgeom was used, re-use the internal positions rather than
      // calculating afresh.
      Box3<int> bboxRpl;

      pcc::point_t minPos = 0;

      if (_gps->predgeom_enabled_flag) {
        altPositions = _posSph;
        bboxRpl = Box3<int>(altPositions.begin(), altPositions.end());
        minPos = bboxRpl.min;
        if (
          attrInterPredParams.enableAttrInterPred
          || biPredEncodeParams.attrInterPredParams2.enableAttrInterPred) {
          for (auto i = 0; i < 3; i++)
            minPos[i] = minPos[i] < minPos_ref[i] ? minPos[i] : minPos_ref[i];
          auto minPos_shift = minPos_ref - minPos;

          if (minPos_shift[0] || minPos_shift[1] || minPos_shift[2])
            offsetAndScaleShift(
              minPos_shift, attr_aps.attr_coord_scale,
              &attrInterPredParams.referencePointCloud[0],
              &attrInterPredParams.referencePointCloud[0]
                + attrInterPredParams.getPointCount());
        }
        minPos_ref = minPos;
      } else {
        altPositions.resize(pointCloud.getPointCount());

        auto laserOrigin = _gbh.geomAngularOrigin(*_gps);
        bboxRpl = convertXyzToRpl(
          laserOrigin, _gps->angularTheta.data(), _gps->angularTheta.size(),
          &pointCloud[0], &pointCloud[0] + pointCloud.getPointCount(),
          altPositions.data());

        if(!attr_aps.attrInterPredictionEnabled){
          minPos = bboxRpl.min;
        }
      }

      offsetAndScale(
        minPos, attr_aps.attr_coord_scale, altPositions.data(),      
        //bboxRpl.min, attr_aps.attr_coord_scale, altPositions.data(),
        altPositions.data() + altPositions.size());

      pointCloud.swapPoints(altPositions);
    }

    // calculate dist2 for this slice
    abh.attr_dist2_delta = 0;
    if (attr_aps.aps_slice_dist2_deltas_present_flag || attrInterPredParams.enableAttrInterPred) {
      // todo(df): this could be set in the sps and refined only if necessary
      auto dist2 =
        estimateDist2(pointCloud, 100, 128, attr_enc.dist2PercentileEstimate);
      abh.attr_dist2_delta = dist2 - attr_aps.dist2;
    }

    // replace the attribute encoder if not compatible
    if (!attrEncoder->isReusable(attr_aps, abh))
      attrEncoder = makeAttributeEncoder();
    if (!attr_aps.spherical_coord_flag)
      for (auto i = 0; i < pointCloud.getPointCount(); i++)
        pointCloud[i] += _sliceOrigin; 

   if (attr_aps.attrInterPredictionEnabled) {
      if (attr_aps.attr_encoding != AttributeEncoding::kRAHTransform) {
        Box3<int> currentFrameBox = pointCloud.computeBoundingBox();
        attrInterPredParams.referencePointCloud = _refFrameAlt.cloud;
        int count = 0;
        for (int i = 0; i < attrInterPredParams.getPointCount(); i++) {
          point_t p = attrInterPredParams.referencePointCloud[i];
          if (currentFrameBox.contains(p)) {
            attrInterPredParams.referencePointCloud[count] = p;
            if (attrInterPredParams.referencePointCloud.hasReflectances())
              attrInterPredParams.referencePointCloud.setReflectance(
                count,
                attrInterPredParams.referencePointCloud.getReflectance(i));
            if (attrInterPredParams.referencePointCloud.hasColors())
              attrInterPredParams.referencePointCloud.setColor(
                count, attrInterPredParams.referencePointCloud.getColor(i));
            count++;
          }
        }
        attrInterPredParams.referencePointCloud.resize(count);
      }
    }

    auto& ctxtMemAttr = _ctxtMemAttrs.at(abh.attr_sps_attr_idx);
    attrEncoder->encode(
      *_sps, attr_sps, attr_aps, abh, ctxtMemAttr, pointCloud, &payload, attrInterPredParams);

    reconSliceAltPositions = pointCloud;
    bool currFrameNotCodedAsB =
      (_gps->biPredictionEnabledFlag
        && !biPredEncodeParams.codeCurrentFrameAsBFrame);
    auto& refCloud = currFrameNotCodedAsB
      ? biPredEncodeParams.attrInterPredParams2.referencePointCloud
      : attrInterPredParams.referencePointCloud;
    if (attr_aps.spherical_coord_flag) {
      refCloud = pointCloud;
      pointCloud.swapPoints(altPositions);
    } else {
      refCloud = pointCloud;
    }

    if (!attr_aps.spherical_coord_flag)
      for (auto i = 0; i < pointCloud.getPointCount(); i++)
        pointCloud[i] -= _sliceOrigin;

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
  auto& refPosSph = currFrameNotCodedAsB
    ? biPredEncodeParams._refPosSph2
    : _refPosSph;
  auto& refFrameSph = currFrameNotCodedAsB
    ? biPredEncodeParams._refFrameSph2
    : _refFrameSph;
  if (_gps->interPredictionEnabledFlag) {
    if (_gps->predgeom_enabled_flag){
      if (currFrameNotCodedAsB){
        biPredEncodeParams._refPosSph2 = _posSph;
        biPredEncodeParams._refFrameSph2.insert(_posSph);
      }
      else{
        _refPosSph = _posSph;
        _refFrameSph.insert(_posSph);
      }
    }
  }

  // Note the current slice id for loss detection with entropy continuation
  _prevSliceId = _sliceId;

  // prevent re-use of this sliceId:  the next slice (geometry + attributes)
  // should be distinguishable from the current slice.
  _sliceId++;
  _firstSliceInFrame = false;

  if (reconCloud)
    appendSlice(reconCloud->cloud);

  if (reconCloudAlt)
    reconCloudAlt->cloud.append(reconSliceAltPositions);
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
  gbh.trisoup_vertex_quantization_bits =
    params->gbh.trisoup_vertex_quantization_bits;
  gbh.trisoup_centroid_vertex_residual_flag =
    params->gbh.trisoup_centroid_vertex_residual_flag;
  gbh.trisoup_face_vertex_flag = params->gbh.trisoup_face_vertex_flag;

  gbh.trisoup_halo_flag =
    params->gbh.trisoup_halo_flag;
  gbh.trisoup_adaptive_halo_flag =
    params->gbh.trisoup_adaptive_halo_flag;
  gbh.trisoup_fine_ray_tracing_flag =
    params->gbh.trisoup_fine_ray_tracing_flag;
  gbh.slice_bb_pos   = params->gbh.slice_bb_pos;
  gbh.slice_bb_width = params->gbh.slice_bb_width;
  gbh.slice_bb_pos_bits = params->gbh.slice_bb_pos_bits;
  gbh.slice_bb_pos_log2_scale = params->gbh.slice_bb_pos_log2_scale;
  gbh.slice_bb_width_bits = params->gbh.slice_bb_width_bits;
  gbh.slice_bb_width_log2_scale = params->gbh.slice_bb_width_log2_scale;
  gbh.interPredictionEnabledFlag = _codeCurrFrameAsInter;
  gbh.geom_qp_offset_intvl_log2_delta =
    params->gbh.geom_qp_offset_intvl_log2_delta;
  gbh.interFrameRefGmcFlag = false;
  gbh.interFrameRefGmcFlag2 = false;
  gbh.gm_matrix = {65536, 0, 0, 0, 65536, 0, 0, 0, 65536};
  gbh.gm_trans = 0;
  gbh.gm_thresh = {0, 0};
  gbh.biPredictionEnabledFlag = biPredEncodeParams.codeCurrentFrameAsBFrame;
  gbh.gm_matrix2 = {65536, 0, 0, 0, 65536, 0, 0, 0, 65536};
  gbh.gm_trans2 = 0;
  gbh.gm_thresh2 = {0, 0};
  if (gbh.interPredictionEnabledFlag && !_gps->predgeom_enabled_flag) {
    gbh.lpu_type = params->interGeom.lpuType;
    gbh.motion_block_size = params->gbh.motion_block_size;
  }

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
    aec->setBypassBinCodingWithoutProbUpdate(_sps->bypass_bin_coding_without_prob_update);
    aec->start();
  }

  // forget (reset) all saved context state at boundary
  if (!gbh.entropy_continuation_flag) {
    if (
      !_sps->inter_entropy_continuation_enabled_flag
      || !gbh.interPredictionEnabledFlag) {
      _ctxtMemOctreeGeom->reset();
      _ctxtMemPredGeom->reset();
      for (auto& ctxtMem : _ctxtMemAttrs)
        ctxtMem.reset();
    }
  }

  if (_gps->predgeom_enabled_flag) {
    _refFrameSph.setInterEnabled(gbh.interPredictionEnabledFlag);
    biPredEncodeParams._refFrameSph2.setInterEnabled(gbh.interPredictionEnabledFlag && gbh.biPredictionEnabledFlag);
    if (!gbh.interPredictionEnabledFlag) {    
      if (_gps->globalMotionEnabled){
        _refFrameSph.updateNextMovingStatus(false);
      }
      _refFrameSph.clearRefFrame(); 
      if (_gps->biPredictionEnabledFlag){
        if (_gps->globalMotionEnabled)
          biPredEncodeParams._refFrameSph2.updateNextMovingStatus(false);
        biPredEncodeParams._refFrameSph2.clearRefFrame();         
      }
      previousAsInter = false;
    }
    else{
      previousAsInter = true;
    }
    encodePredictiveGeometry(
      params->predGeom, *_gps, gbh, pointCloud, &_posSph, _refFrameSph, biPredEncodeParams._refFrameSph2,
      *_ctxtMemPredGeom, arithmeticEncoders[0].get());
  } else if (!_gps->trisoup_enabled_flag) {
    encodeGeometryOctree(
      params->geom, *_gps, gbh, pointCloud, *_ctxtMemOctreeGeom,
      arithmeticEncoders, _refFrame, *_sps, params->interGeom,
      biPredEncodeParams);
  }
  else
  {
    // limit the number of points to the slice limit
    // todo(df): this should be derived from the level
    gbh.footer.geom_num_points_minus1 = params->partition.sliceMaxPointsTrisoup - 1;
    encodeGeometryTrisoup(
      params->trisoup, params->geom, *_gps, gbh, pointCloud, pointCloudPadding,
      *_ctxtMemOctreeGeom, arithmeticEncoders, _refFrame,
      *_sps, params->interGeom);
  }
  // signal the actual number of points coded
  gbh.footer.geom_num_points_minus1 = pointCloud.getPointCount() - 1;

  if ( _gps->predgeom_enabled_flag && gbh.interPredictionEnabledFlag && _gps->globalMotionEnabled) {
    gbh.interFrameRefGmcFlag = _refFrameSph.getFrameMovingState();
    if (gbh.biPredictionEnabledFlag)
      gbh.interFrameRefGmcFlag2 = biPredEncodeParams._refFrameSph2.getFrameMovingState();
    if (_gps->biPredictionEnabledFlag){
      _refFrameSph.getMotionParamsMultiple(gbh.gm_thresh, gbh.gm_matrix, gbh.gm_trans,
        biPredEncodeParams.currentFrameIndex, biPredEncodeParams.refFrameIndex);
      if (gbh.biPredictionEnabledFlag)
        biPredEncodeParams._refFrameSph2.getMotionParamsMultiple(gbh.gm_thresh2, gbh.gm_matrix2, gbh.gm_trans2,
        biPredEncodeParams.currentFrameIndex, biPredEncodeParams.refFrameIndex2);
    }
    else
      _refFrameSph.getMotionParams(gbh.gm_thresh, gbh.gm_matrix, gbh.gm_trans);
  }

  attrInterPredParams.frameDistance = 1;
  movingState = false;
  biPredEncodeParams.movingState2 = false;

  if (gbh.interPredictionEnabledFlag && params->attributeIdxMap.size()
    && _aps[params->attributeIdxMap.begin()->second]->attr_encoding
    != AttributeEncoding::kRAHTransform) {
    auto checkMovingState =
      [&](const VecInt& mat, const pcc::point_t& tran) -> bool {
      const double scale = 65536.;
      const double thr1_pre = 0.1 / attrInterPredParams.frameDistance;
      const double thr2_pre = params->attrInterPredTranslationThreshold;

      const double thr1_tan_pre = tan(M_PI * thr1_pre / 180);
      const double thr1_sin_pre = sin(M_PI * thr1_pre / 180);
      const double Rx = std::abs((mat[5] / scale) / (1. + mat[8] / scale));
      const double Ry = std::abs(mat[2] / scale);
      const double Rz = std::abs((mat[1] / scale) / (1. + mat[0] / scale));
      const double Sx = std::abs(tran[0]);
      const double Sy = std::abs(tran[1]);
      const double Sz = std::abs(tran[2]);

      return (
        Rx < thr1_tan_pre&& Ry < thr1_sin_pre&& Rz < thr1_tan_pre
        && Sx < thr2_pre&& Sy < thr2_pre&& Sz < thr2_pre);
    };

    movingState = checkMovingState(gbh.gm_matrix, gbh.gm_trans);

    if (gbh.biPredictionEnabledFlag) {
      biPredEncodeParams.movingState2 =
        checkMovingState(gbh.gm_matrix2, gbh.gm_trans2);
    }
  }

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
