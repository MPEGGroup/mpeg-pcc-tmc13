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

#pragma once

#include "PCCMath.h"

#include <cstdint>
#include <ostream>
#include <vector>

namespace pcc {

//============================================================================

enum class PayloadType
{
  kSequenceParameterSet = 0,
  kGeometryParameterSet = 1,
  kGeometryBrick = 2,
  kAttributeParameterSet = 3,
  kAttributeBrick = 4,
  kTileInventory = 5,
};

//============================================================================

enum class KnownAttributeLabel : uint32_t
{
  kColour = 0,
  kReflectance = 1,
};

//============================================================================

struct AttributeLabel {
  uint32_t attribute_label_four_bytes = 0xffffffffu;

  //--------------------------------------------------------------------------

  AttributeLabel() = default;

  AttributeLabel(KnownAttributeLabel known_attribute_label)
  {
    attribute_label_four_bytes = int(known_attribute_label);
  }

  //--------------------------------------------------------------------------

  friend bool
  operator==(const AttributeLabel& lhs, const KnownAttributeLabel& rhs)
  {
    return uint32_t(rhs) == lhs.attribute_label_four_bytes;
  }

  //--------------------------------------------------------------------------

  bool known_attribute_label_flag() const
  {
    switch (KnownAttributeLabel(attribute_label_four_bytes)) {
    case KnownAttributeLabel::kColour:
    case KnownAttributeLabel::kReflectance: return true;
    }

    return false;
  }

  //--------------------------------------------------------------------------

  KnownAttributeLabel known_attribute_label() const
  {
    return KnownAttributeLabel(attribute_label_four_bytes);
  }
};

//============================================================================

std::ostream& operator<<(std::ostream& os, const AttributeLabel& label);

//============================================================================

enum class AttributeEncoding
{
  kPredictingTransform = 0,
  kRAHTransform = 1,
  kLiftingTransform = 2,
};

//============================================================================

// invariant properties
struct AttributeDescription {
  int attr_num_dimensions;
  int attr_instance_id;
  int attr_bitdepth;
  int cicp_colour_primaries_idx;
  int cicp_transfer_characteristics_idx;
  int cicp_matrix_coefficients_idx;
  bool cicp_video_full_range_flag;

  AttributeLabel attributeLabel;
};

//============================================================================

struct ProfileCompatibility {
  int profile_compatibility_flags;
};

//============================================================================

struct SequenceParameterSet {
  int sps_seq_parameter_set_id;

  ProfileCompatibility profileCompatibility;
  int level;

  // todo(df): encode the following
  PCCVector3<int> seq_bounding_box_xyz0;
  PCCVector3<int> seq_bounding_box_whd;
  //int seq_bounding_box_scale_log2;

  // A value describing the scaling of the source positions prior to encoding.
  float seq_source_geom_scale_factor;

  // NB: attributeSets.size() = num_attribute_sets
  std::vector<AttributeDescription> attributeSets;
};

//============================================================================

struct GeometryParameterSet {
  int gps_geom_parameter_set_id;
  int gps_seq_parameter_set_id;

  // Indicates that the GeometryBrickHeader contains a valid
  // geom_box_origin_xyz.
  int geom_box_present_flag;

  // Controls the ability to represent multiple points (with associated
  // attributes) at the same spatial position.
  bool geom_unique_points_flag;

  // Controls the use of neighbour based contextualisation of octree
  // occupancy during geometry coding.  When true, only neighbours that
  // are direct siblings are available.
  bool neighbour_context_restriction_flag;

  // Defines the size of the neighbour availiability volume (aka
  // look-ahead cube size) for occupancy searches.  A value of 0
  // indicates that the feature is disabled.
  int neighbour_avail_boundary_log2;

  // Controls the use of early termination of the geometry tree
  // by directly coding the position of isolated points.
  bool inferred_direct_coding_mode_enabled_flag;

  // Selects between bitwise and bytewise occupancy coding
  bool bitwise_occupancy_coding_flag;

  // Experimental knob to control the number of contexts used
  // for occupancy coding.
  int geom_occupancy_ctx_reduction_factor;

  // Maximum node size where intra prediction is enabled
  int intra_pred_max_node_size_log2;

  // size of triangle nodes (reconstructed surface) in trisoup geometry.
  // a value of zero disables the feature
  int trisoup_node_size_log2;
};

//============================================================================

struct GeometryBrickHeader {
  int geom_geom_parameter_set_id;
  int geom_tile_id;
  int geom_slice_id;

  // derived from geom_box_origin_{x,y,z} * (1 << geom_box_log2_scale)
  PCCVector3<int> geomBoxOrigin;
  int geom_box_log2_scale;

  // todo(df): minus1?
  int geom_max_node_size_log2;
  int geom_num_points;
};

//============================================================================

struct AttributeParameterSet {
  int aps_attr_parameter_set_id;
  int aps_seq_parameter_set_id;
  AttributeEncoding attr_encoding;

  //--- lifting/predicting transform parameters

  bool lod_binary_tree_enabled_flag;
  int num_pred_nearest_neighbours;
  int max_num_direct_predictors;
  int adaptive_prediction_threshold;
  int search_range;

  // NB: derived from num_detail_levels_minus1
  int num_detail_levels;
  std::vector<int64_t> dist2;

  // NB: these parameters are shared by raht and lift
  int quant_step_size_luma;
  int quant_step_size_chroma;

  //--- raht parameters

  int raht_depth;
  int raht_binary_level_threshold;
};

//============================================================================

struct AttributeBrickHeader {
  int attr_sps_attr_idx;
  int attr_attr_parameter_set_id;
  int attr_geom_slice_id;
};

//============================================================================

struct TileInventory {
  struct Entry {
    PCCVector3<int> tile_bounding_box_xyz0;
    PCCVector3<int> tile_bounding_box_whd;
  };
  std::vector<Entry> tiles;
};

//============================================================================

}  // namespace pcc
