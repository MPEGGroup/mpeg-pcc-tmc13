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

  bs.writeUe(sps.sps_seq_parameter_set_id);

  int num_attribute_sets = int(sps.attributeSets.size());
  bs.writeUe(num_attribute_sets);
  for (const auto& attr : sps.attributeSets) {
    // todo(df): should be attr_num_dimensions_minus1
    bs.writeUe(attr.attr_num_dimensions);
    bs.writeUe(attr.attr_instance_id);
    bs.writeUe(attr.attr_bitdepth);
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

  bs.readUe(&sps.sps_seq_parameter_set_id);

  int num_attribute_sets = int(bs.readUe());
  for (int i = 0; i < num_attribute_sets; i++) {
    sps.attributeSets.emplace_back();
    auto& attr = sps.attributeSets.back();
    bs.readUe(&attr.attr_num_dimensions);
    bs.readUe(&attr.attr_instance_id);
    bs.readUe(&attr.attr_bitdepth);
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
  bs.write(gps.geom_unique_points_flag);
  bs.write(gps.neighbour_context_restriction_flag);
  bs.write(gps.inferred_direct_coding_mode_enabled_flag);
  bs.write(gps.bitwise_occupancy_coding_flag);
  bs.writeUe(gps.geom_occupancy_ctx_reduction_factor);
  bs.writeUe(gps.neighbour_avail_boundary_log2);
  bs.writeUe(gps.intra_pred_max_node_size_log2);
  bs.writeUe(gps.trisoup_node_size_log2);

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
  bs.read(&gps.geom_unique_points_flag);
  bs.read(&gps.neighbour_context_restriction_flag);
  bs.read(&gps.inferred_direct_coding_mode_enabled_flag);
  bs.read(&gps.bitwise_occupancy_coding_flag);
  bs.readUe(&gps.geom_occupancy_ctx_reduction_factor);
  bs.readUe(&gps.neighbour_avail_boundary_log2);
  bs.readUe(&gps.intra_pred_max_node_size_log2);
  bs.readUe(&gps.trisoup_node_size_log2);

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

  bool isLifting = aps.attr_encoding == AttributeEncoding::kLiftingTransform
    || aps.attr_encoding == AttributeEncoding::kPredictingTransform;
  if (isLifting) {
    bs.writeUe(aps.num_pred_nearest_neighbours);
    bs.writeUe(aps.max_num_direct_predictors);
    bs.writeUe(aps.search_range);
    bs.writeUe(aps.quant_step_size_luma);
    bs.writeUe(aps.quant_step_size_chroma);
    bs.write(aps.lod_binary_tree_enabled_flag);

    bs.writeUe(aps.num_detail_levels);
    for (int idx = 0; idx < aps.num_detail_levels; idx++) {
      // todo(??): is this an appropriate encoding?
      bs.writeUe64(aps.dist2[idx]);
    }
  }

  if (aps.attr_encoding == AttributeEncoding::kPredictingTransform) {
    bs.writeUe(aps.adaptive_prediction_threshold);
  }

  if (aps.attr_encoding == AttributeEncoding::kRAHTransform) {
    bs.writeUe(aps.raht_depth);
    bs.writeUe(aps.raht_binary_level_threshold);
    bs.writeUe(aps.quant_step_size_luma);
    // todo(?): raht chroma quant_step_size?
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

  bool isLifting = aps.attr_encoding == AttributeEncoding::kLiftingTransform
    || aps.attr_encoding == AttributeEncoding::kPredictingTransform;
  if (isLifting) {
    bs.readUe(&aps.num_pred_nearest_neighbours);
    bs.readUe(&aps.max_num_direct_predictors);
    bs.readUe(&aps.search_range);
    bs.readUe(&aps.quant_step_size_luma);
    bs.readUe(&aps.quant_step_size_chroma);
    bs.read(&aps.lod_binary_tree_enabled_flag);

    aps.num_detail_levels = int(bs.readUe());
    aps.dist2.resize(aps.num_detail_levels);
    for (int idx = 0; idx < aps.num_detail_levels; idx++) {
      bs.readUe(&aps.dist2[idx]);
    }
  }

  if (aps.attr_encoding == AttributeEncoding::kPredictingTransform) {
    bs.readUe(&aps.adaptive_prediction_threshold);
  }

  if (aps.attr_encoding == AttributeEncoding::kRAHTransform) {
    bs.readUe(&aps.raht_depth);
    bs.readUe(&aps.raht_binary_level_threshold);
    bs.readUe(&aps.quant_step_size_luma);
    // todo(?): raht chroma quant_step_size
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
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kGeometryBrick);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  bs.writeUe(gbh.geom_geom_parameter_set_id);
  bs.writeUe(gbh.geom_tile_id);
  bs.writeUe(gbh.geom_slice_id);

  if (gps.geom_box_present_flag) {
    int geom_box_origin_x = gbh.geomBoxOrigin.x() >> gbh.geom_box_log2_scale;
    int geom_box_origin_y = gbh.geomBoxOrigin.y() >> gbh.geom_box_log2_scale;
    int geom_box_origin_z = gbh.geomBoxOrigin.z() >> gbh.geom_box_log2_scale;

    bs.writeUe(gbh.geom_box_log2_scale);
    bs.writeUe(geom_box_origin_x);
    bs.writeUe(geom_box_origin_y);
    bs.writeUe(geom_box_origin_z);
  }

  bs.writeUe(gbh.geom_max_node_size_log2);
  bs.writeUe(gbh.geom_num_points);
  bs.byteAlign();
}

//----------------------------------------------------------------------------

GeometryBrickHeader
parseGbh(
  const GeometryParameterSet& gps, const PayloadBuffer& buf, int* bytesRead)
{
  GeometryBrickHeader gbh;
  assert(buf.type == PayloadType::kGeometryBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUe(&gbh.geom_geom_parameter_set_id);
  bs.readUe(&gbh.geom_tile_id);
  bs.readUe(&gbh.geom_slice_id);

  if (gps.geom_box_present_flag) {
    bs.readUe(&gbh.geom_box_log2_scale);

    int geom_box_origin_x;
    int geom_box_origin_y;
    int geom_box_origin_z;
    bs.readUe(&geom_box_origin_x);
    bs.readUe(&geom_box_origin_y);
    bs.readUe(&geom_box_origin_z);
    gbh.geomBoxOrigin.x() = geom_box_origin_x << gbh.geom_box_log2_scale;
    gbh.geomBoxOrigin.y() = geom_box_origin_y << gbh.geom_box_log2_scale;
    gbh.geomBoxOrigin.z() = geom_box_origin_z << gbh.geom_box_log2_scale;
  }

  bs.readUe(&gbh.geom_max_node_size_log2);
  bs.readUe(&gbh.geom_num_points);
  bs.byteAlign();

  if (bytesRead)
    *bytesRead = int(std::distance(buf.begin(), bs.pos()));

  return gbh;
}

//============================================================================

void
write(const AttributeBrickHeader& abh, PayloadBuffer* buf)
{
  assert(buf->type == PayloadType::kAttributeBrick);
  auto bs = makeBitWriter(std::back_inserter(*buf));

  bs.writeUe(abh.attr_attr_parameter_set_id);
  bs.writeUe(abh.attr_sps_attr_idx);
  bs.writeUe(abh.attr_geom_slice_id);

  bs.byteAlign();
}

//----------------------------------------------------------------------------

AttributeBrickHeader
parseAbh(const PayloadBuffer& buf, int* bytesRead)
{
  AttributeBrickHeader abh;
  assert(buf.type == PayloadType::kAttributeBrick);
  auto bs = makeBitReader(buf.begin(), buf.end());

  bs.readUe(&abh.attr_attr_parameter_set_id);
  bs.readUe(&abh.attr_sps_attr_idx);
  bs.readUe(&abh.attr_geom_slice_id);

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
