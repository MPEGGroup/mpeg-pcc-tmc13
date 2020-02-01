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

#include "BitReader.h"
#include "BitWriter.h"
#include "hls.h"
#include "io_hls.h"

#include <iomanip>
#include <iterator>
#include <algorithm>

namespace pcc {

//============================================================================

std::ostream&
operator<<(std::ostream& os, const AttributeLabel& label)
{
  switch (KnownAttributeLabel(label.attribute_label_four_bytes)) {
  case KnownAttributeLabel::kColour: os << "color"; break;
  case KnownAttributeLabel::kReflectance: os << "reflectance"; break;
  default:
    auto iosFlags = os.flags(std::ios::hex);
    os << std::setw(8) << label.attribute_label_four_bytes;
    os.flags(iosFlags);
  }

  return os;
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps)
{
  PayloadBuffer buf(PayloadType::kSequenceParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUn(24, sps.profileCompatibility.profile_compatibility_flags);
  bs.writeUn(8, sps.level);
  bs.writeUe(sps.sps_seq_parameter_set_id);

  bool seq_bounding_box_present_flag = true;
  bs.write(seq_bounding_box_present_flag);
  if (seq_bounding_box_present_flag) {
    bs.writeSe(sps.seq_bounding_box_xyz0.x());
    bs.writeSe(sps.seq_bounding_box_xyz0.y());
    bs.writeSe(sps.seq_bounding_box_xyz0.z());
    int seq_bounding_box_scale = 1;
    bs.writeUe(seq_bounding_box_scale);
    bs.writeUe(sps.seq_bounding_box_whd.x());
    bs.writeUe(sps.seq_bounding_box_whd.y());
    bs.writeUe(sps.seq_bounding_box_whd.z());
  }
  // todo(df): determine encoding of scale factor
  bs.writeF(sps.seq_source_geom_scale_factor);

  int num_attribute_sets = int(sps.attributeSets.size());
  bs.writeUe(num_attribute_sets);
  for (const auto& attr : sps.attributeSets) {
    // todo(df): should be attr_num_dimensions_minus1
    bs.writeUe(attr.attr_num_dimensions);
    bs.writeUe(attr.attr_instance_id);
    bs.writeUe(attr.attr_bitdepth);
    if (attr.attr_num_dimensions > 1)
      bs.writeUe(attr.attr_bitdepth_secondary);

    bs.writeUe(attr.cicp_colour_primaries_idx);
    bs.writeUe(attr.cicp_transfer_characteristics_idx);
    bs.writeUe(attr.cicp_matrix_coefficients_idx);
    bs.write(attr.cicp_video_full_range_flag);

    const auto& label = attr.attributeLabel;

    bool known_attribute_label_flag = label.known_attribute_label_flag();
    bs.write(known_attribute_label_flag);
    if (known_attribute_label_flag) {
      bs.writeUe(int(label.known_attribute_label()));
    } else {
      bs.writeUn(32, label.attribute_label_four_bytes);
    }
  }

  bs.writeUn(5, sps.log2_max_frame_idx);
  bs.writeUn(3, sps.geometry_axis_order);
  bs.write(sps.cabac_bypass_stream_enabled_flag);

  bool sps_extension_flag = false;
  bs.write(sps_extension_flag);
  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

SequenceParameterSet
parseSps(const PayloadBuffer& buf)
{
  SequenceParameterSet sps;
  assert(buf.type == PayloadType::kSequenceParameterSet);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(24, &sps.profileCompatibility.profile_compatibility_flags);
  bs.readUn(8, &sps.level);
  bs.readUe(&sps.sps_seq_parameter_set_id);

  bool seq_bounding_box_present_flag = bs.read();
  if (seq_bounding_box_present_flag) {
    bs.readSe(&sps.seq_bounding_box_xyz0.x());
    bs.readSe(&sps.seq_bounding_box_xyz0.y());
    bs.readSe(&sps.seq_bounding_box_xyz0.z());
    int seq_bounding_box_scale;
    bs.readUe(&seq_bounding_box_scale);
    bs.readUe(&sps.seq_bounding_box_whd.x());
    bs.readUe(&sps.seq_bounding_box_whd.y());
    bs.readUe(&sps.seq_bounding_box_whd.z());
  }
  bs.readF(&sps.seq_source_geom_scale_factor);

  int num_attribute_sets = int(bs.readUe());
  for (int i = 0; i < num_attribute_sets; i++) {
    sps.attributeSets.emplace_back();
    auto& attr = sps.attributeSets.back();
    bs.readUe(&attr.attr_num_dimensions);
    bs.readUe(&attr.attr_instance_id);
    bs.readUe(&attr.attr_bitdepth);
    if (attr.attr_num_dimensions > 1)
      bs.readUe(&attr.attr_bitdepth_secondary);

    bs.readUe(&attr.cicp_colour_primaries_idx);
    bs.readUe(&attr.cicp_transfer_characteristics_idx);
    bs.readUe(&attr.cicp_matrix_coefficients_idx);
    bs.read(&attr.cicp_video_full_range_flag);

    auto& label = attr.attributeLabel;

    bool known_attribute_label_flag = bs.read();
    if (known_attribute_label_flag) {
      KnownAttributeLabel known_attribute_label;
      bs.readUe(&known_attribute_label);
      label = known_attribute_label;
    } else {
      bs.readUn(32, &label.attribute_label_four_bytes);
    }
  }

  bs.readUn(5, &sps.log2_max_frame_idx);
  bs.readUn(3, &sps.geometry_axis_order);
  bs.read(&sps.cabac_bypass_stream_enabled_flag);

  bool sps_extension_flag = bs.read();
  if (sps_extension_flag) {
    // todo(df): sps_extension_data;
    assert(!sps_extension_flag);
  }
  bs.byteAlign();

  return sps;
}

//============================================================================

PayloadBuffer
write(const GeometryParameterSet& gps)
{
  PayloadBuffer buf(PayloadType::kGeometryParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUe(gps.gps_geom_parameter_set_id);
  bs.writeUe(gps.gps_seq_parameter_set_id);
  bs.write(gps.geom_box_present_flag);
  if (gps.geom_box_present_flag) {
    bs.write(gps.geom_box_log2_scale_present_flag);
    if (!gps.geom_box_log2_scale_present_flag)
      bs.writeUe(gps.gps_geom_box_log2_scale);
  }
  bs.write(gps.geom_unique_points_flag);
  bs.write(gps.implicit_qtbt_enabled_flag);
  bs.write(gps.neighbour_context_restriction_flag);
  bs.write(gps.inferred_direct_coding_mode_enabled_flag);
  bs.write(gps.bitwise_occupancy_coding_flag);
  bs.write(gps.adjacent_child_contextualization_enabled_flag);

  bs.write(gps.geom_planar_mode_enabled_flag);
  if (gps.geom_planar_mode_enabled_flag) {
    bs.writeUe(gps.geom_planar_threshold0);
    bs.writeUe(gps.geom_planar_threshold1);
    bs.writeUe(gps.geom_planar_threshold2);
    bs.writeUe(gps.geom_planar_idcm_threshold);
    bs.write(gps.geom_angular_mode_enabled_flag);
  }

  if (gps.geom_angular_mode_enabled_flag) {
    bs.writeUe(gps.geom_angular_lidar_head_position[0]);
    bs.writeUe(gps.geom_angular_lidar_head_position[1]);
    bs.writeUe(gps.geom_angular_lidar_head_position[2]);
    bs.writeUe(gps.geom_angular_num_lidar_lasers());
    for (int i = 0; i < gps.geom_angular_num_lidar_lasers(); i++) {
      bs.writeUe(gps.geom_angular_theta_laser[i] + 1048576);
      bs.writeUe(gps.geom_angular_z_laser[i] + 1048576);
    }
    bs.write(gps.planar_buffer_disabled_flag);
  }

  if (gps.geom_angular_mode_enabled_flag && gps.implicit_qtbt_enabled_flag) {
    bs.writeUe(gps.implicit_qtbt_angular_max_node_min_dim_log2_to_split_z);
    bs.writeUe(gps.implicit_qtbt_angular_max_diff_to_split_z);
  }

  bs.writeUe(gps.geom_occupancy_ctx_reduction_factor);
  bs.writeUe(gps.neighbour_avail_boundary_log2);
  bs.writeUe(gps.intra_pred_max_node_size_log2);
  bs.writeUe(gps.trisoup_node_size_log2);
  bs.write(gps.geom_scaling_enabled_flag);
  if (gps.geom_scaling_enabled_flag)
    bs.writeUe(gps.geom_base_qp_minus4);

  if (gps.implicit_qtbt_enabled_flag)
    if (!gps.trisoup_node_size_log2) {
      bs.writeUe(gps.max_num_implicit_qtbt_before_ot);
      bs.writeUe(gps.min_implicit_qtbt_size_log2);
    }

  bool gps_extension_flag = false;
  bs.write(gps_extension_flag);
  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

GeometryParameterSet
parseGps(const PayloadBuffer& buf)
{
  GeometryParameterSet gps;
  assert(buf.type == PayloadType::kGeometryParameterSet);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUe(&gps.gps_geom_parameter_set_id);
  bs.readUe(&gps.gps_seq_parameter_set_id);
  bs.read(&gps.geom_box_present_flag);
  if (gps.geom_box_present_flag) {
    bs.read(&gps.geom_box_log2_scale_present_flag);
    if (!gps.geom_box_log2_scale_present_flag)
      bs.readUe(&gps.gps_geom_box_log2_scale);
  }
  bs.read(&gps.geom_unique_points_flag);
  bs.read(&gps.implicit_qtbt_enabled_flag);
  bs.read(&gps.neighbour_context_restriction_flag);
  bs.read(&gps.inferred_direct_coding_mode_enabled_flag);
  bs.read(&gps.bitwise_occupancy_coding_flag);
  bs.read(&gps.adjacent_child_contextualization_enabled_flag);

  gps.geom_angular_mode_enabled_flag = false;
  bs.read(&gps.geom_planar_mode_enabled_flag);
  if (gps.geom_planar_mode_enabled_flag) {
    bs.readUe(&gps.geom_planar_threshold0);
    bs.readUe(&gps.geom_planar_threshold1);
    bs.readUe(&gps.geom_planar_threshold2);
    bs.readUe(&gps.geom_planar_idcm_threshold);
    bs.read(&gps.geom_angular_mode_enabled_flag);
  }

  gps.planar_buffer_disabled_flag = false;
  if (gps.geom_angular_mode_enabled_flag) {
    bs.readUe(&gps.geom_angular_lidar_head_position[0]);
    bs.readUe(&gps.geom_angular_lidar_head_position[1]);
    bs.readUe(&gps.geom_angular_lidar_head_position[2]);

    int geom_angular_num_lidar_lasers;
    bs.readUe(&geom_angular_num_lidar_lasers);
    gps.geom_angular_theta_laser.resize(geom_angular_num_lidar_lasers);
    gps.geom_angular_z_laser.resize(geom_angular_num_lidar_lasers);
    for (int i = 0; i < geom_angular_num_lidar_lasers; i++) {
      bs.readUe(&gps.geom_angular_theta_laser[i]);
      bs.readUe(&gps.geom_angular_z_laser[i]);
      gps.geom_angular_theta_laser[i] -= 1048576;
      gps.geom_angular_z_laser[i] -= 1048576;
    }
    bs.read(&gps.planar_buffer_disabled_flag);
  }

  gps.implicit_qtbt_angular_max_node_min_dim_log2_to_split_z = 0;
  gps.implicit_qtbt_angular_max_diff_to_split_z = 0;
  if (gps.geom_angular_mode_enabled_flag && gps.implicit_qtbt_enabled_flag) {
    bs.readUe(&gps.implicit_qtbt_angular_max_node_min_dim_log2_to_split_z);
    bs.readUe(&gps.implicit_qtbt_angular_max_diff_to_split_z);
  }

  bs.readUe(&gps.geom_occupancy_ctx_reduction_factor);
  bs.readUe(&gps.neighbour_avail_boundary_log2);
  bs.readUe(&gps.intra_pred_max_node_size_log2);
  bs.readUe(&gps.trisoup_node_size_log2);

  bs.read(&gps.geom_scaling_enabled_flag);
  if (gps.geom_scaling_enabled_flag)
    bs.readUe(&gps.geom_base_qp_minus4);

  gps.max_num_implicit_qtbt_before_ot = 0;
  gps.min_implicit_qtbt_size_log2 = 0;
  if (gps.implicit_qtbt_enabled_flag)
    if (gps.trisoup_node_size_log2 == 0) {
      bs.readUe(&gps.max_num_implicit_qtbt_before_ot);
      bs.readUe(&gps.min_implicit_qtbt_size_log2);
    }

  bool gps_extension_flag = bs.read();
  if (gps_extension_flag) {
    // todo(df): gps_extension_data;
    assert(!gps_extension_flag);
  }
  bs.byteAlign();

  return gps;
}

//============================================================================

PayloadBuffer
write(const AttributeParameterSet& aps)
{
  PayloadBuffer buf(PayloadType::kAttributeParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUe(aps.aps_attr_parameter_set_id);
  bs.writeUe(aps.aps_seq_parameter_set_id);
  bs.writeUe(aps.attr_encoding);

  bs.writeUe(aps.init_qp);
  bs.writeSe(aps.aps_chroma_qp_offset);
  bs.write(aps.aps_slice_qp_deltas_present_flag);

  if (aps.lodParametersPresent()) {
    bs.writeUe(aps.num_pred_nearest_neighbours);
    bs.writeUe(aps.search_range);
    bs.writeUe(aps.lod_neigh_bias.x());
    bs.writeUe(aps.lod_neigh_bias.y());
    bs.writeUe(aps.lod_neigh_bias.z());

    if (aps.attr_encoding == AttributeEncoding::kLiftingTransform)
      bs.write(aps.scalable_lifting_enabled_flag);

    if (!aps.scalable_lifting_enabled_flag) {
      bs.writeUe(aps.num_detail_levels);
      if (aps.num_detail_levels > 0)
        bs.write(aps.lod_decimation_enabled_flag);

      for (int idx = 0; idx < aps.num_detail_levels; idx++) {
        // todo(??): is this an appropriate encoding?
        bs.writeUe64(aps.dist2[idx]);
      }
    }
  }

  if (aps.attr_encoding == AttributeEncoding::kPredictingTransform) {
    bs.writeUe(aps.adaptive_prediction_threshold);
    bs.write(aps.intra_lod_prediction_enabled_flag);
    bs.write(aps.inter_component_prediction_enabled_flag);
    bs.writeUe(aps.max_num_direct_predictors);
  }

  if (aps.attr_encoding == AttributeEncoding::kRAHTransform) {
    bs.write(aps.raht_prediction_enabled_flag);
    if (aps.raht_prediction_enabled_flag) {
      bs.writeUe(aps.raht_prediction_threshold0);
      bs.writeUe(aps.raht_prediction_threshold1);
    }
  }

  bool aps_extension_flag = false;
  bs.write(aps_extension_flag);
  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

AttributeParameterSet
parseAps(const PayloadBuffer& buf)
{
  AttributeParameterSet aps;
  assert(buf.type == PayloadType::kAttributeParameterSet);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUe(&aps.aps_attr_parameter_set_id);
  bs.readUe(&aps.aps_seq_parameter_set_id);
  bs.readUe(&aps.attr_encoding);

  bs.readUe(&aps.init_qp);
  bs.readSe(&aps.aps_chroma_qp_offset);
  bs.read(&aps.aps_slice_qp_deltas_present_flag);

  if (aps.lodParametersPresent()) {
    bs.readUe(&aps.num_pred_nearest_neighbours);
    bs.readUe(&aps.search_range);
    bs.readUe(&aps.lod_neigh_bias.x());
    bs.readUe(&aps.lod_neigh_bias.y());
    bs.readUe(&aps.lod_neigh_bias.z());

    aps.scalable_lifting_enabled_flag = false;
    if (aps.attr_encoding == AttributeEncoding::kLiftingTransform)
      bs.read(&aps.scalable_lifting_enabled_flag);

    if (!aps.scalable_lifting_enabled_flag) {
      bs.readUe(&aps.num_detail_levels);
      if (aps.num_detail_levels > 0)
        bs.read(&aps.lod_decimation_enabled_flag);

      aps.dist2.resize(aps.num_detail_levels);
      for (int idx = 0; idx < aps.num_detail_levels; idx++) {
        bs.readUe(&aps.dist2[idx]);
      }
    }
  }

  aps.intra_lod_prediction_enabled_flag = false;
  if (aps.attr_encoding == AttributeEncoding::kPredictingTransform) {
    bs.readUe(&aps.adaptive_prediction_threshold);
    bs.read(&aps.intra_lod_prediction_enabled_flag);
    bs.read(&aps.inter_component_prediction_enabled_flag);
    bs.readUe(&aps.max_num_direct_predictors);
  }

  if (aps.attr_encoding == AttributeEncoding::kRAHTransform) {
    bs.read(&aps.raht_prediction_enabled_flag);
    if (aps.raht_prediction_enabled_flag) {
      bs.readUe(&aps.raht_prediction_threshold0);
      bs.readUe(&aps.raht_prediction_threshold1);
    }
  }

  bool aps_extension_flag = bs.read();
  if (aps_extension_flag) {
    // todo(df): aps_extension_data;
    assert(!aps_extension_flag);
  }
  bs.byteAlign();

  return aps;
}

//============================================================================

void
write(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kGeometryBrick);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  bs.writeUn(24, gbh.geom_num_points);
  bs.writeUe(gbh.geom_geom_parameter_set_id);
  bs.writeUe(gbh.geom_tile_id);
  bs.writeUe(gbh.geom_slice_id);
  bs.writeUn(sps.log2_max_frame_idx, gbh.frame_idx);

  if (gps.geom_box_present_flag) {
    int geomBoxLog2Scale = gbh.geomBoxLog2Scale(gps);
    int geom_box_origin_x = gbh.geomBoxOrigin.x() >> geomBoxLog2Scale;
    int geom_box_origin_y = gbh.geomBoxOrigin.y() >> geomBoxLog2Scale;
    int geom_box_origin_z = gbh.geomBoxOrigin.z() >> geomBoxLog2Scale;

    if (gps.geom_box_log2_scale_present_flag)
      bs.writeUe(gbh.geom_box_log2_scale);
    bs.writeUe(geom_box_origin_x);
    bs.writeUe(geom_box_origin_y);
    bs.writeUe(geom_box_origin_z);
  }

  if (!gps.implicit_qtbt_enabled_flag) {
    bs.writeUe(gbh.geom_max_node_size_log2);
  } else {
    auto& geom_max_node_size_log2_x = gbh.geom_max_node_size_log2_xyz[0];

    int geom_max_node_size_log2_delta_y =
      gbh.geom_max_node_size_log2_xyz[1] - gbh.geom_max_node_size_log2_xyz[0];

    int geom_max_node_size_log2_delta_z =
      gbh.geom_max_node_size_log2_xyz[2] - gbh.geom_max_node_size_log2_xyz[1];

    bs.writeUe(geom_max_node_size_log2_x);
    bs.writeSe(geom_max_node_size_log2_delta_y);
    bs.writeSe(geom_max_node_size_log2_delta_z);
  }

  if (gbh.geom_octree_parallel_max_node_size_log2 == 0)
    bs.writeUe(0);
  else {
    bs.writeUe(gbh.geom_octree_parallel_max_node_size_log2 - 1);
    bs.writeUn(6, gbh.geom_octree_parallel_max_offset_log2);

    // NB: the length of the last substream is not signalled
    for (int i = 0; i < gbh.geom_octree_parallel_max_node_size_log2; i++) {
      int width = gbh.geom_octree_parallel_max_offset_log2;
      bs.writeUn(width, gbh.geom_octree_parallel_bitstream_offsets[i]);
    }
  }

  if (gps.geom_scaling_enabled_flag) {
    bs.writeSe(gbh.geom_slice_qp_offset);
    bs.writeUe(gbh.geom_octree_qp_offset_depth);
  }

  bs.byteAlign();
}

//----------------------------------------------------------------------------

GeometryBrickHeader
parseGbh(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const PayloadBuffer& buf,
  int* bytesRead)
{
  GeometryBrickHeader gbh;
  assert(buf.type == PayloadType::kGeometryBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(24, &gbh.geom_num_points);
  bs.readUe(&gbh.geom_geom_parameter_set_id);
  bs.readUe(&gbh.geom_tile_id);
  bs.readUe(&gbh.geom_slice_id);
  bs.readUn(sps.log2_max_frame_idx, &gbh.frame_idx);

  if (gps.geom_box_present_flag) {
    if (gps.geom_box_log2_scale_present_flag)
      bs.readUe(&gbh.geom_box_log2_scale);

    int geom_box_origin_x;
    int geom_box_origin_y;
    int geom_box_origin_z;
    bs.readUe(&geom_box_origin_x);
    bs.readUe(&geom_box_origin_y);
    bs.readUe(&geom_box_origin_z);
    int geomBoxLog2Scale = gbh.geomBoxLog2Scale(gps);
    gbh.geomBoxOrigin.x() = geom_box_origin_x << geomBoxLog2Scale;
    gbh.geomBoxOrigin.y() = geom_box_origin_y << geomBoxLog2Scale;
    gbh.geomBoxOrigin.z() = geom_box_origin_z << geomBoxLog2Scale;
  }

  if (!gps.implicit_qtbt_enabled_flag) {
    bs.readUe(&gbh.geom_max_node_size_log2);
  } else {
    int geom_max_node_size_log2_x;
    int geom_max_node_size_log2_delta_y;
    int geom_max_node_size_log2_delta_z;

    bs.readUe(&geom_max_node_size_log2_x);
    bs.readSe(&geom_max_node_size_log2_delta_y);
    bs.readSe(&geom_max_node_size_log2_delta_z);

    gbh.geom_max_node_size_log2_xyz[0] = geom_max_node_size_log2_x;

    gbh.geom_max_node_size_log2_xyz[1] =
      geom_max_node_size_log2_delta_y + gbh.geom_max_node_size_log2_xyz[0];

    gbh.geom_max_node_size_log2_xyz[2] =
      geom_max_node_size_log2_delta_z + gbh.geom_max_node_size_log2_xyz[1];
  }

  bs.readUe(&gbh.geom_octree_parallel_max_node_size_log2);
  if (gbh.geom_octree_parallel_max_node_size_log2 > 0) {
    gbh.geom_octree_parallel_max_node_size_log2++;

    gbh.geom_octree_parallel_bitstream_offsets.resize(
      gbh.geom_octree_parallel_max_node_size_log2);

    bs.readUn(6, &gbh.geom_octree_parallel_max_offset_log2);
    for (int i = 0; i < gbh.geom_octree_parallel_max_node_size_log2; i++) {
      int width = gbh.geom_octree_parallel_max_offset_log2;
      bs.readUn(width, &gbh.geom_octree_parallel_bitstream_offsets[i]);
    }
  }

  if (gps.geom_scaling_enabled_flag) {
    bs.readSe(&gbh.geom_slice_qp_offset);
    bs.readUe(&gbh.geom_octree_qp_offset_depth);
  }

  bs.byteAlign();

  if (bytesRead)
    *bytesRead = int(std::distance(buf.begin(), bs.pos()));

  return gbh;
}

//----------------------------------------------------------------------------

GeometryBrickHeader
parseGbhIds(const PayloadBuffer& buf)
{
  GeometryBrickHeader gbh;
  assert(buf.type == PayloadType::kGeometryBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUn(24, &gbh.geom_num_points);
  bs.readUe(&gbh.geom_geom_parameter_set_id);
  bs.readUe(&gbh.geom_tile_id);
  bs.readUe(&gbh.geom_slice_id);

  /* NB: this function only decodes ids at the start of the header. */
  /* NB: do not attempt to parse any further */

  return gbh;
}

//============================================================================

void
write(
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kAttributeBrick);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  bs.writeUe(abh.attr_attr_parameter_set_id);
  bs.writeUe(abh.attr_sps_attr_idx);
  bs.writeUe(abh.attr_geom_slice_id);

  if (aps.aps_slice_qp_deltas_present_flag) {
    bs.writeSe(abh.attr_qp_delta_luma);
    bs.writeSe(abh.attr_qp_delta_chroma);
  }

  bool attr_layer_qp_present_flag = !abh.attr_layer_qp_delta_luma.empty();
  bs.write(attr_layer_qp_present_flag);
  if (attr_layer_qp_present_flag) {
    int attr_num_qp_layers_minus1 = abh.attr_num_qp_layers_minus1();
    bs.writeUe(attr_num_qp_layers_minus1);
    for (int i = 0; i <= attr_num_qp_layers_minus1; i++) {
      bs.writeSe(abh.attr_layer_qp_delta_luma[i]);
      bs.writeSe(abh.attr_layer_qp_delta_chroma[i]);
    }
  }

  bs.write(abh.attr_region_qp_present_flag);
  if (abh.attr_region_qp_present_flag) {
    bs.writeUe(abh.attr_region_qp_origin.x());
    bs.writeUe(abh.attr_region_qp_origin.y());
    bs.writeUe(abh.attr_region_qp_origin.z());
    bs.writeUe(abh.attr_region_qp_whd.x());
    bs.writeUe(abh.attr_region_qp_whd.y());
    bs.writeUe(abh.attr_region_qp_whd.z());
    bs.writeSe(abh.attr_region_qp_delta);
  }
  bs.byteAlign();
}

//----------------------------------------------------------------------------

AttributeBrickHeader
parseAbhIds(const PayloadBuffer& buf)
{
  AttributeBrickHeader abh;
  assert(buf.type == PayloadType::kAttributeBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUe(&abh.attr_attr_parameter_set_id);
  bs.readUe(&abh.attr_sps_attr_idx);
  bs.readUe(&abh.attr_geom_slice_id);

  /* NB: this function only decodes ids at the start of the header. */
  /* NB: do not attempt to parse any further */

  return abh;
}

//----------------------------------------------------------------------------

AttributeBrickHeader
parseAbh(
  const AttributeParameterSet& aps, const PayloadBuffer& buf, int* bytesRead)
{
  AttributeBrickHeader abh;
  assert(buf.type == PayloadType::kAttributeBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUe(&abh.attr_attr_parameter_set_id);
  bs.readUe(&abh.attr_sps_attr_idx);
  bs.readUe(&abh.attr_geom_slice_id);

  if (aps.aps_slice_qp_deltas_present_flag) {
    bs.readSe(&abh.attr_qp_delta_luma);
    bs.readSe(&abh.attr_qp_delta_chroma);
  }

  bool attr_layer_qp_present_flag;
  bs.read(&attr_layer_qp_present_flag);
  if (attr_layer_qp_present_flag) {
    int attr_num_qp_layers_minus1;
    bs.readUe(&attr_num_qp_layers_minus1);
    abh.attr_layer_qp_delta_luma.resize(attr_num_qp_layers_minus1 + 1);
    abh.attr_layer_qp_delta_chroma.resize(attr_num_qp_layers_minus1 + 1);
    for (int i = 0; i <= attr_num_qp_layers_minus1; i++) {
      bs.readSe(&abh.attr_layer_qp_delta_luma[i]);
      bs.readSe(&abh.attr_layer_qp_delta_chroma[i]);
    }
  }

  bs.read(&abh.attr_region_qp_present_flag);
  if (abh.attr_region_qp_present_flag) {
    bs.readUe(&abh.attr_region_qp_origin.x());
    bs.readUe(&abh.attr_region_qp_origin.y());
    bs.readUe(&abh.attr_region_qp_origin.z());
    bs.readUe(&abh.attr_region_qp_whd.x());
    bs.readUe(&abh.attr_region_qp_whd.y());
    bs.readUe(&abh.attr_region_qp_whd.z());
    bs.readSe(&abh.attr_region_qp_delta);
  }

  bs.byteAlign();

  if (bytesRead)
    *bytesRead = int(std::distance(buf.begin(), bs.pos()));

  return abh;
}

//============================================================================

PayloadBuffer
write(const TileInventory& inventory)
{
  PayloadBuffer buf(PayloadType::kTileInventory);
  auto bs = makeBitWriter(std::back_inserter(buf));

  int num_tiles = inventory.tiles.size();
  bs.writeUe(num_tiles);
  for (const auto& entry : inventory.tiles) {
    bs.writeSe(entry.tile_bounding_box_xyz0.x());
    bs.writeSe(entry.tile_bounding_box_xyz0.y());
    bs.writeSe(entry.tile_bounding_box_xyz0.z());
    bs.writeUe(entry.tile_bounding_box_whd.x());
    bs.writeUe(entry.tile_bounding_box_whd.y());
    bs.writeUe(entry.tile_bounding_box_whd.z());
  }

  bs.byteAlign();

  return buf;
}

//----------------------------------------------------------------------------

TileInventory
parseTileInventory(const PayloadBuffer& buf)
{
  TileInventory inventory;
  assert(buf.type == PayloadType::kTileInventory);
  auto bs = makeBitReader(buf.begin(), buf.end());

  int num_tiles;
  bs.readUe(&num_tiles);
  for (int i = 0; i < num_tiles; i++) {
    TileInventory::Entry entry;
    bs.readSe(&entry.tile_bounding_box_xyz0.x());
    bs.readSe(&entry.tile_bounding_box_xyz0.y());
    bs.readSe(&entry.tile_bounding_box_xyz0.z());
    bs.readUe(&entry.tile_bounding_box_whd.x());
    bs.readUe(&entry.tile_bounding_box_whd.y());
    bs.readUe(&entry.tile_bounding_box_whd.z());
    inventory.tiles.push_back(entry);
  }

  bs.byteAlign();

  return inventory;
}

//============================================================================

}  // namespace pcc
