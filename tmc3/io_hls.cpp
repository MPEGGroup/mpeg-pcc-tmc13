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
#include "PCCMisc.h"
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
    auto sps_bounding_box_offset_xyz =
      toXyz(sps.geometry_axis_order, sps.seqBoundingBoxOrigin);

    bs.writeSe(sps_bounding_box_offset_xyz.x());
    bs.writeSe(sps_bounding_box_offset_xyz.y());
    bs.writeSe(sps_bounding_box_offset_xyz.z());

    int seq_bounding_box_offset_log2_scale = 0;
    bs.writeUe(seq_bounding_box_offset_log2_scale);

    auto seq_bounding_box_whd =
      toXyz(sps.geometry_axis_order, sps.seqBoundingBoxSize);

    bs.writeUe(seq_bounding_box_whd.x());
    bs.writeUe(seq_bounding_box_whd.y());
    bs.writeUe(seq_bounding_box_whd.z());
  }
  // todo(df): determine encoding of scale factor
  bs.writeF(sps.seq_source_geom_scale_factor);

  int num_attribute_sets = int(sps.attributeSets.size());
  bs.writeUe(num_attribute_sets);
  for (const auto& attr : sps.attributeSets) {
    bs.writeUe(attr.attr_num_dimensions_minus1);
    bs.writeUe(attr.attr_instance_id);

    int attr_bitdepth_minus1 = attr.bitdepth - 1;
    bs.writeUe(attr_bitdepth_minus1);

    if (attr.attr_num_dimensions_minus1) {
      int attr_bitdepth_secondary_minus1 = attr.bitdepthSecondary - 1;
      bs.writeUe(attr_bitdepth_secondary_minus1);
    }

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
    Vec3<int> seq_bounding_box_offset;
    bs.readSe(&seq_bounding_box_offset.x());
    bs.readSe(&seq_bounding_box_offset.y());
    bs.readSe(&seq_bounding_box_offset.z());

    int seq_bounding_box_offset_log2_scale;
    bs.readUe(&seq_bounding_box_offset_log2_scale);
    seq_bounding_box_offset *= 1 << seq_bounding_box_offset_log2_scale;

    Vec3<int> seq_bounding_box_whd;
    bs.readUe(&seq_bounding_box_whd.x());
    bs.readUe(&seq_bounding_box_whd.y());
    bs.readUe(&seq_bounding_box_whd.z());

    // NB: these are in XYZ axis order until the SPS is converted to STV
    sps.seqBoundingBoxOrigin = seq_bounding_box_offset;
    sps.seqBoundingBoxSize = seq_bounding_box_whd;
  }
  bs.readF(&sps.seq_source_geom_scale_factor);

  int num_attribute_sets = int(bs.readUe());
  for (int i = 0; i < num_attribute_sets; i++) {
    sps.attributeSets.emplace_back();
    auto& attr = sps.attributeSets.back();
    bs.readUe(&attr.attr_num_dimensions_minus1);
    bs.readUe(&attr.attr_instance_id);

    int attr_bitdepth_minus1;
    bs.readUe(&attr_bitdepth_minus1);
    attr.bitdepth = attr_bitdepth_minus1 + 1;

    if (attr.attr_num_dimensions_minus1) {
      int attr_bitdepth_secondary_minus1 = 0;
      bs.readUe(&attr_bitdepth_secondary_minus1);
      attr.bitdepthSecondary = attr_bitdepth_secondary_minus1 + 1;
    }

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

//----------------------------------------------------------------------------

void
convertXyzToStv(SequenceParameterSet* sps)
{
  // permute the bounding box from xyz to internal stv order
  sps->seqBoundingBoxOrigin =
    fromXyz(sps->geometry_axis_order, sps->seqBoundingBoxOrigin);

  sps->seqBoundingBoxSize =
    fromXyz(sps->geometry_axis_order, sps->seqBoundingBoxSize);
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps, const GeometryParameterSet& gps)
{
  PayloadBuffer buf(PayloadType::kGeometryParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUe(gps.gps_geom_parameter_set_id);
  bs.writeUe(gps.gps_seq_parameter_set_id);
  bs.write(gps.geom_box_log2_scale_present_flag);
  if (!gps.geom_box_log2_scale_present_flag)
    bs.writeUe(gps.gps_geom_box_log2_scale);
  bs.write(gps.predgeom_enabled_flag);
  bs.write(gps.geom_unique_points_flag);

  if (!gps.predgeom_enabled_flag) {
    bs.write(gps.qtbt_enabled_flag);
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
    }

    bs.write(gps.geom_angular_mode_enabled_flag);
    if (gps.geom_angular_mode_enabled_flag) {
      auto geom_angular_origin =
        toXyz(sps.geometry_axis_order, gps.geomAngularOrigin);
      bs.writeUe(geom_angular_origin.x());
      bs.writeUe(geom_angular_origin.y());
      bs.writeUe(geom_angular_origin.z());
      bs.writeUe(gps.geom_angular_num_lidar_lasers());

      if (gps.geom_angular_num_lidar_lasers()) {
        bs.writeSe(gps.geom_angular_theta_laser[0]);
        bs.writeSe(gps.geom_angular_z_laser[0]);
        bs.writeUe(gps.geom_angular_num_phi_per_turn[0]);
      }

      for (int i = 1; i < gps.geom_angular_num_lidar_lasers(); i++) {
        int geom_angular_theta_laser_diff = gps.geom_angular_theta_laser[i]
          - gps.geom_angular_theta_laser[i - 1];

        int geom_angular_z_laser_diff =
          gps.geom_angular_z_laser[i] - gps.geom_angular_z_laser[i - 1];

        // NB: angles must be in increasing monotonic order
        assert(geom_angular_theta_laser_diff >= 0);
        bs.writeUe(geom_angular_theta_laser_diff);
        bs.writeSe(geom_angular_z_laser_diff);
        bs.writeUe(gps.geom_angular_num_phi_per_turn[i]);
      }
      bs.write(gps.planar_buffer_disabled_flag);
    }

    bs.writeUe(gps.neighbour_avail_boundary_log2);
    bs.writeUe(gps.intra_pred_max_node_size_log2);
    bs.writeUe(gps.trisoup_node_size_log2);
    bs.write(gps.geom_scaling_enabled_flag);
    if (gps.geom_scaling_enabled_flag) {
      bs.writeUe(gps.geom_base_qp);
      bs.writeSe(gps.geom_idcm_qp_offset);
    }
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
  bs.read(&gps.geom_box_log2_scale_present_flag);
  if (!gps.geom_box_log2_scale_present_flag)
    bs.readUe(&gps.gps_geom_box_log2_scale);
  bs.read(&gps.predgeom_enabled_flag);
  bs.read(&gps.geom_unique_points_flag);

  if (!gps.predgeom_enabled_flag) {
    bs.read(&gps.qtbt_enabled_flag);
    bs.read(&gps.neighbour_context_restriction_flag);
    bs.read(&gps.inferred_direct_coding_mode_enabled_flag);
    bs.read(&gps.bitwise_occupancy_coding_flag);
    bs.read(&gps.adjacent_child_contextualization_enabled_flag);

    bs.read(&gps.geom_planar_mode_enabled_flag);
    if (gps.geom_planar_mode_enabled_flag) {
      bs.readUe(&gps.geom_planar_threshold0);
      bs.readUe(&gps.geom_planar_threshold1);
      bs.readUe(&gps.geom_planar_threshold2);
      bs.readUe(&gps.geom_planar_idcm_threshold);
    }

    gps.planar_buffer_disabled_flag = false;
    bs.read(&gps.geom_angular_mode_enabled_flag);
    if (gps.geom_angular_mode_enabled_flag) {
      Vec3<int> geom_angular_origin;
      bs.readUe(&geom_angular_origin.x());
      bs.readUe(&geom_angular_origin.y());
      bs.readUe(&geom_angular_origin.z());

      // NB: this is in XYZ axis order until the GPS is converted to STV
      gps.geomAngularOrigin = geom_angular_origin;

      int geom_angular_num_lidar_lasers;
      bs.readUe(&geom_angular_num_lidar_lasers);
      gps.geom_angular_theta_laser.resize(geom_angular_num_lidar_lasers);
      gps.geom_angular_z_laser.resize(geom_angular_num_lidar_lasers);
      gps.geom_angular_num_phi_per_turn.resize(geom_angular_num_lidar_lasers);

      if (geom_angular_num_lidar_lasers) {
        bs.readSe(&gps.geom_angular_theta_laser[0]);
        bs.readSe(&gps.geom_angular_z_laser[0]);
        bs.readUe(&gps.geom_angular_num_phi_per_turn[0]);
      }

      for (int i = 1; i < geom_angular_num_lidar_lasers; i++) {
        int geom_angular_theta_laser_diff;
        int geom_angular_z_laser_diff;
        bs.readUe(&geom_angular_theta_laser_diff);
        bs.readSe(&geom_angular_z_laser_diff);
        bs.readUe(&gps.geom_angular_num_phi_per_turn[i]);

        gps.geom_angular_theta_laser[i] =
          gps.geom_angular_theta_laser[i - 1] + geom_angular_theta_laser_diff;

        gps.geom_angular_z_laser[i] =
          gps.geom_angular_z_laser[i - 1] + geom_angular_z_laser_diff;
      }
      bs.read(&gps.planar_buffer_disabled_flag);
    }

    bs.readUe(&gps.neighbour_avail_boundary_log2);
    bs.readUe(&gps.intra_pred_max_node_size_log2);
    bs.readUe(&gps.trisoup_node_size_log2);

    gps.geom_base_qp = 0;
    gps.geom_idcm_qp_offset = 0;
    bs.read(&gps.geom_scaling_enabled_flag);
    if (gps.geom_scaling_enabled_flag) {
      bs.readUe(&gps.geom_base_qp);
      bs.readSe(&gps.geom_idcm_qp_offset);
    }
  }

  bool gps_extension_flag = bs.read();
  if (gps_extension_flag) {
    // todo(df): gps_extension_data;
    assert(!gps_extension_flag);
  }
  bs.byteAlign();

  return gps;
}

//----------------------------------------------------------------------------

void
convertXyzToStv(const SequenceParameterSet& sps, GeometryParameterSet* gps)
{
  gps->geomAngularOrigin =
    fromXyz(sps.geometry_axis_order, gps->geomAngularOrigin);
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps, const AttributeParameterSet& aps)
{
  PayloadBuffer buf(PayloadType::kAttributeParameterSet);
  auto bs = makeBitWriter(std::back_inserter(buf));

  bs.writeUe(aps.aps_attr_parameter_set_id);
  bs.writeUe(aps.aps_seq_parameter_set_id);
  bs.writeUe(aps.attr_encoding);

  bs.writeUe(aps.init_qp_minus4);
  bs.writeSe(aps.aps_chroma_qp_offset);
  bs.write(aps.aps_slice_qp_deltas_present_flag);

  if (aps.lodParametersPresent()) {
    bs.writeUe(aps.num_pred_nearest_neighbours_minus1);
    bs.writeUe(aps.search_range);

    auto lod_neigh_bias = toXyz(sps.geometry_axis_order, aps.lodNeighBias);
    bs.writeUe(lod_neigh_bias.x());
    bs.writeUe(lod_neigh_bias.y());
    bs.writeUe(lod_neigh_bias.z());

    if (aps.attr_encoding == AttributeEncoding::kLiftingTransform) {
      bs.write(aps.scalable_lifting_enabled_flag);
      if (aps.scalable_lifting_enabled_flag)
        bs.writeUe(aps.max_neigh_range);
    }

    if (!aps.scalable_lifting_enabled_flag) {
      bs.writeUe(aps.num_detail_levels);
      if (!aps.num_detail_levels)
        bs.write(aps.canonical_point_order_flag);
      else {
        bs.write(aps.lod_decimation_enabled_flag);

        if (aps.lod_decimation_enabled_flag) {
          for (int idx = 0; idx < aps.num_detail_levels; idx++) {
            auto lod_sampling_period_minus2 = aps.lodSamplingPeriod[idx] - 2;
            bs.writeUe(lod_sampling_period_minus2);
          }
        } else {
          for (int idx = 0; idx < aps.num_detail_levels; idx++) {
            auto numerator = aps.dist2[idx];
            auto denominator = idx > 0 ? aps.dist2[idx - 1] : 1;
            int lod_sampling_scale_minus1 = (numerator / denominator) - 1;
            int lod_sampling_offset = numerator % denominator;
            bs.writeUe(lod_sampling_scale_minus1);
            if (idx > 0)
              bs.writeUe(lod_sampling_offset);
          }
        }
      }
    }
  }

  if (aps.attr_encoding == AttributeEncoding::kPredictingTransform) {
    bs.writeUe(aps.max_num_direct_predictors);
    if (aps.max_num_direct_predictors)
      bs.writeUe(aps.adaptive_prediction_threshold);
    bs.write(aps.intra_lod_prediction_enabled_flag);
    bs.write(aps.inter_component_prediction_enabled_flag);
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

  bs.readUe(&aps.init_qp_minus4);
  bs.readSe(&aps.aps_chroma_qp_offset);
  bs.read(&aps.aps_slice_qp_deltas_present_flag);

  if (aps.lodParametersPresent()) {
    bs.readUe(&aps.num_pred_nearest_neighbours_minus1);
    bs.readUe(&aps.search_range);

    Vec3<int> lod_neigh_bias;
    bs.readUe(&lod_neigh_bias.x());
    bs.readUe(&lod_neigh_bias.y());
    bs.readUe(&lod_neigh_bias.z());
    // NB: this is in XYZ axis order until the GPS is converted to STV
    aps.lodNeighBias = lod_neigh_bias;

    aps.scalable_lifting_enabled_flag = false;
    if (aps.attr_encoding == AttributeEncoding::kLiftingTransform) {
      bs.read(&aps.scalable_lifting_enabled_flag);
      if (aps.scalable_lifting_enabled_flag)
        bs.readUe(&aps.max_neigh_range);
    }

    aps.canonical_point_order_flag = false;
    if (!aps.scalable_lifting_enabled_flag) {
      bs.readUe(&aps.num_detail_levels);
      if (!aps.num_detail_levels)
        bs.read(&aps.canonical_point_order_flag);
      else {
        bs.read(&aps.lod_decimation_enabled_flag);

        if (aps.lod_decimation_enabled_flag) {
          aps.lodSamplingPeriod.resize(aps.num_detail_levels);
          for (int idx = 0; idx < aps.num_detail_levels; idx++) {
            int lod_sampling_period_minus2;
            bs.readUe(&lod_sampling_period_minus2);
            aps.lodSamplingPeriod[idx] = lod_sampling_period_minus2 + 2;
          }
        } else {
          aps.dist2.resize(aps.num_detail_levels);
          for (int idx = 0; idx < aps.num_detail_levels; idx++) {
            int lod_sampling_scale_minus1;
            int lod_sampling_offset = 0;
            bs.readUe(&lod_sampling_scale_minus1);
            if (idx == 0)
              aps.dist2[idx] = lod_sampling_scale_minus1 + 1;
            else {
              int lod_sampling_offset;
              bs.readUe(&lod_sampling_offset);
              aps.dist2[idx] =
                aps.dist2[idx - 1] * (lod_sampling_scale_minus1 + 1)
                + lod_sampling_offset;
            }
          }
        }
      }
    }
  }

  aps.intra_lod_prediction_enabled_flag = false;
  if (aps.attr_encoding == AttributeEncoding::kPredictingTransform) {
    bs.readUe(&aps.max_num_direct_predictors);
    aps.adaptive_prediction_threshold = 0;
    if (aps.max_num_direct_predictors)
      bs.readUe(&aps.adaptive_prediction_threshold);
    bs.read(&aps.intra_lod_prediction_enabled_flag);
    bs.read(&aps.inter_component_prediction_enabled_flag);
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

//----------------------------------------------------------------------------

void
convertXyzToStv(const SequenceParameterSet& sps, AttributeParameterSet* aps)
{
  aps->lodNeighBias = fromXyz(sps.geometry_axis_order, aps->lodNeighBias);
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

  bs.writeUn(24, gbh.geom_num_points_minus1);
  bs.writeUe(gbh.geom_geom_parameter_set_id);
  bs.writeUe(gbh.geom_tile_id);
  bs.writeUe(gbh.geom_slice_id);
  bs.writeUn(sps.log2_max_frame_idx, gbh.frame_idx);

  int geomBoxLog2Scale = gbh.geomBoxLog2Scale(gps);
  auto geom_box_origin = toXyz(sps.geometry_axis_order, gbh.geomBoxOrigin);
  geom_box_origin.x() >>= geomBoxLog2Scale;
  geom_box_origin.y() >>= geomBoxLog2Scale;
  geom_box_origin.z() >>= geomBoxLog2Scale;

  if (gps.geom_box_log2_scale_present_flag)
    bs.writeUe(gbh.geom_box_log2_scale);
  bs.writeUe(geom_box_origin.x());
  bs.writeUe(geom_box_origin.y());
  bs.writeUe(geom_box_origin.z());

  if (!gps.predgeom_enabled_flag) {
    int tree_depth_minus1 = gbh.tree_lvl_coded_axis_list.size() - 1;
    bs.writeUe(tree_depth_minus1);
    if (gps.qtbt_enabled_flag)
      for (int i = 0; i <= tree_depth_minus1; i++)
        bs.writeUn(3, gbh.tree_lvl_coded_axis_list[i]);

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

  bs.readUn(24, &gbh.geom_num_points_minus1);
  bs.readUe(&gbh.geom_geom_parameter_set_id);
  bs.readUe(&gbh.geom_tile_id);
  bs.readUe(&gbh.geom_slice_id);
  bs.readUn(sps.log2_max_frame_idx, &gbh.frame_idx);

  if (gps.geom_box_log2_scale_present_flag)
    bs.readUe(&gbh.geom_box_log2_scale);

  Vec3<int> geom_box_origin;
  bs.readUe(&geom_box_origin.x());
  bs.readUe(&geom_box_origin.y());
  bs.readUe(&geom_box_origin.z());
  gbh.geomBoxOrigin = fromXyz(sps.geometry_axis_order, geom_box_origin);
  gbh.geomBoxOrigin *= 1 << gbh.geomBoxLog2Scale(gps);

  if (!gps.predgeom_enabled_flag) {
    int tree_depth_minus1;
    bs.readUe(&tree_depth_minus1);

    gbh.tree_lvl_coded_axis_list.resize(tree_depth_minus1 + 1, 7);
    if (gps.qtbt_enabled_flag)
      for (int i = 0; i <= tree_depth_minus1; i++)
        bs.readUn(3, &gbh.tree_lvl_coded_axis_list[i]);

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

  bs.readUn(24, &gbh.geom_num_points_minus1);
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
  const SequenceParameterSet& sps,
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
    auto attr_region_qp_origin =
      toXyz(sps.geometry_axis_order, abh.regionQpOrigin);

    auto attr_region_qp_whd = toXyz(sps.geometry_axis_order, abh.regionQpSize);

    bs.writeUe(attr_region_qp_origin.x());
    bs.writeUe(attr_region_qp_origin.y());
    bs.writeUe(attr_region_qp_origin.z());
    bs.writeUe(attr_region_qp_whd.x());
    bs.writeUe(attr_region_qp_whd.y());
    bs.writeUe(attr_region_qp_whd.z());
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
  const SequenceParameterSet& sps,
  const AttributeParameterSet& aps,
  const PayloadBuffer& buf,
  int* bytesRead)
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
    Vec3<int> attr_region_qp_origin;
    bs.readUe(&attr_region_qp_origin.x());
    bs.readUe(&attr_region_qp_origin.y());
    bs.readUe(&attr_region_qp_origin.z());
    abh.regionQpOrigin =
      fromXyz(sps.geometry_axis_order, attr_region_qp_origin);

    Vec3<int> attr_region_qp_whd;
    bs.readUe(&attr_region_qp_whd.x());
    bs.readUe(&attr_region_qp_whd.y());
    bs.readUe(&attr_region_qp_whd.z());
    abh.regionQpSize = fromXyz(sps.geometry_axis_order, attr_region_qp_whd);

    bs.readSe(&abh.attr_region_qp_delta);
  }

  bs.byteAlign();

  if (bytesRead)
    *bytesRead = int(std::distance(buf.begin(), bs.pos()));

  return abh;
}

//============================================================================

PayloadBuffer
write(const SequenceParameterSet& sps, const TileInventory& inventory)
{
  PayloadBuffer buf(PayloadType::kTileInventory);
  auto bs = makeBitWriter(std::back_inserter(buf));

  int num_tiles = inventory.tiles.size();
  bs.writeUn(16, num_tiles);

  // calculate the maximum size of any values
  int maxVal = 1;
  for (const auto& entry : inventory.tiles) {
    for (int k = 0; k < 3; k++) {
      maxVal = std::max(maxVal, entry.tileOrigin[k]);
      maxVal = std::max(maxVal, entry.tileSize[k]);
    }
  }

  int tile_bounding_box_bits = ceillog2(uint32_t(maxVal));
  bs.writeUn(8, tile_bounding_box_bits);

  for (const auto& entry : inventory.tiles) {
    auto tile_origin = toXyz(sps.geometry_axis_order, entry.tileOrigin);
    bs.writeSn(tile_bounding_box_bits, tile_origin.x());
    bs.writeSn(tile_bounding_box_bits, tile_origin.y());
    bs.writeSn(tile_bounding_box_bits, tile_origin.z());

    auto tile_size = toXyz(sps.geometry_axis_order, entry.tileSize);
    bs.writeUn(tile_bounding_box_bits, tile_size.x());
    bs.writeUn(tile_bounding_box_bits, tile_size.y());
    bs.writeUn(tile_bounding_box_bits, tile_size.z());
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
  bs.readUn(16, &num_tiles);

  int tile_bounding_box_bits;
  bs.readUn(8, &tile_bounding_box_bits);

  for (int i = 0; i < num_tiles; i++) {
    Vec3<int> tile_origin;
    bs.readSn(tile_bounding_box_bits, &tile_origin.x());
    bs.readSn(tile_bounding_box_bits, &tile_origin.y());
    bs.readSn(tile_bounding_box_bits, &tile_origin.z());
    Vec3<int> tile_size;
    bs.readUn(tile_bounding_box_bits, &tile_size.x());
    bs.readUn(tile_bounding_box_bits, &tile_size.y());
    bs.readUn(tile_bounding_box_bits, &tile_size.z());

    // NB: this is in XYZ axis order until the inventory is converted to STV
    TileInventory::Entry entry;
    entry.tileOrigin = tile_origin;
    entry.tileSize = tile_size;
    inventory.tiles.push_back(entry);
  }

  bs.byteAlign();

  return inventory;
}

//----------------------------------------------------------------------------

void
convertXyzToStv(const SequenceParameterSet& sps, TileInventory* inventory)
{
  for (auto& tile : inventory->tiles) {
    tile.tileOrigin = fromXyz(sps.geometry_axis_order, tile.tileOrigin);
    tile.tileSize = fromXyz(sps.geometry_axis_order, tile.tileSize);
  }
}

//============================================================================

}  // namespace pcc
