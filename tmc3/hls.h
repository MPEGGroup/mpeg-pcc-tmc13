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
  kFrameBoundaryMarker = 6,
  kConstantAttribute = 7,
  kGeneralizedAttrParamInventory = 8,
  kUserData = 9,
};

//============================================================================

enum class KnownAttributeLabel : uint32_t
{
  kColour = 0,
  kReflectance = 1,
  kOpacity = 2,
  kFrameIndex = 3,
  kFrameNumber = 4,
  kMaterialId = 5,
  kNormal = 6,

  // Indicates that the attrabute label is described by an Oid
  kOid = 0xffffffff,
};

//============================================================================

struct Oid {
  // A sequence of encoded subidentifiers according to Rec. ITU-T X.690 |
  // ISO/IEC 8825-1.  NB: this does not include any identifier octets, length
  // octets or end-of-content octets of the basic encoding rules.
  std::vector<uint8_t> contents;

  Oid() = default;
  Oid(const std::string& str);

  // Convert the oid to a string representation
  operator std::string() const;

  friend bool operator==(const Oid& lhs, const Oid& rhs);
};

//============================================================================

struct AttributeLabel {
  KnownAttributeLabel known_attribute_label;
  Oid oid;

  //--------------------------------------------------------------------------

  AttributeLabel() = default;

  AttributeLabel(KnownAttributeLabel known_attribute_label)
    : known_attribute_label(known_attribute_label)
  {}

  //--------------------------------------------------------------------------

  friend bool
  operator==(const AttributeLabel& lhs, const KnownAttributeLabel& rhs)
  {
    return lhs.known_attribute_label == rhs;
  }

  //--------------------------------------------------------------------------

  bool known_attribute_label_flag() const
  {
    return known_attribute_label != KnownAttributeLabel::kOid;
  }
};

//============================================================================

std::ostream& operator<<(std::ostream& os, const AttributeLabel& label);

//============================================================================

enum class AttributeEncoding
{
  kRAHTransform = 0,
  kPredictingTransform = 1,
  kLiftingTransform = 2,
  kRaw = 3,
};

//============================================================================

enum class LodDecimationMethod
{
  kNone = 0,
  kPeriodic = 1,
  kCentroid = 2,
};

//============================================================================

enum class AxisOrder
{
  kZYX = 0,
  kXYZ = 1,
  kXZY = 2,
  kYZX = 3,
  kZYX_4 = 4,
  kZXY = 5,
  kYXZ = 6,
  kXYZ_7 = 7,
};

// Permute the internal STV axes to XYZ order.
template<typename T>
Vec3<T>
toXyz(AxisOrder order, const Vec3<T>& stv)
{
  switch (order) {
  case AxisOrder::kZYX: return {stv.v(), stv.t(), stv.s()};
  case AxisOrder::kXYZ: return {stv.s(), stv.t(), stv.v()};
  case AxisOrder::kXZY: return {stv.s(), stv.v(), stv.t()};
  case AxisOrder::kYZX: return {stv.v(), stv.s(), stv.t()};
  case AxisOrder::kZYX_4: return {stv.v(), stv.t(), stv.s()};
  case AxisOrder::kZXY: return {stv.t(), stv.v(), stv.s()};
  case AxisOrder::kYXZ: return {stv.t(), stv.s(), stv.v()};
  case AxisOrder::kXYZ_7: return {stv.s(), stv.t(), stv.v()};
  }
}

// Permute the an XYZ axis order to the internal STV order.
template<typename T>
Vec3<T>
fromXyz(AxisOrder order, const Vec3<T>& xyz)
{
  switch (order) {
  case AxisOrder::kZYX: return {xyz.z(), xyz.y(), xyz.x()};
  case AxisOrder::kXYZ: return {xyz.x(), xyz.y(), xyz.z()};
  case AxisOrder::kXZY: return {xyz.x(), xyz.z(), xyz.y()};
  case AxisOrder::kYZX: return {xyz.y(), xyz.z(), xyz.x()};
  case AxisOrder::kZYX_4: return {xyz.z(), xyz.y(), xyz.x()};
  case AxisOrder::kZXY: return {xyz.z(), xyz.x(), xyz.y()};
  case AxisOrder::kYXZ: return {xyz.y(), xyz.x(), xyz.z()};
  case AxisOrder::kXYZ_7: return {xyz.x(), xyz.y(), xyz.z()};
  }
}

//============================================================================
// ISO/IEC 23001-8 codec independent code points
enum class ColourMatrix : uint8_t
{
  kIdentity = 0,
  kBt709 = 1,
  kUnspecified = 2,
  kReserved_3 = 3,
  kUsa47Cfr73dot682a20 = 4,
  kBt601 = 5,
  kSmpte170M = 6,
  kSmpte240M = 7,
  kYCgCo = 8,
  kBt2020Ncl = 9,
  kBt2020Cl = 10,
  kSmpte2085 = 11,
};

//============================================================================

enum class AttributeParameterType : uint8_t
{
  kItuT35 = 0,
  kOid = 1,
  kCicp = 2,
  kScaling = 3,
  kDefaultValue = 4,
  /* [3, 127] are reserved for future use */
  /* [128, 255] are specified according to the attribute label */
};

//============================================================================

struct OpaqueAttributeParameter {
  // the type of the data
  AttributeParameterType attr_param_type;

  // identifies the type of attr_param_byte data when attr_param_type = 0.
  int attr_param_itu_t_t35_country_code;
  int attr_param_itu_t_t35_country_code_extension;

  // identifies the type of attr_param_byte data when attr_param_type = 1.
  Oid attr_param_oid;

  // the attribute data excluding type0/type1 identification bytes */
  std::vector<uint8_t> attr_param_byte;
};

//============================================================================

struct AttributeParameters {
  // indicates if the cicp attribute parameter is valid
  bool cicpParametersPresent;
  int cicp_colour_primaries_idx;
  int cicp_transfer_characteristics_idx;
  ColourMatrix cicp_matrix_coefficients_idx;
  bool cicp_video_full_range_flag;

  // attribute scaling
  bool scalingParametersPresent;
  int attr_scale_minus1;
  int attr_frac_bits;
  int attr_offset;

  // soft default attribute values.
  // If empty, use hard default values.
  std::vector<int> attr_default_value;

  // Unknown attribute parameters
  std::vector<OpaqueAttributeParameter> opaqueParameters;

  int numParams() const
  {
    int count = opaqueParameters.size() + cicpParametersPresent
      + scalingParametersPresent + !attr_default_value.empty();
    return count;
  }

  void clear()
  {
    cicpParametersPresent = false;
    scalingParametersPresent = false;
    attr_default_value.clear();
    opaqueParameters.clear();
  }
};

//============================================================================

struct AttributeDescription {
  int attr_num_dimensions_minus1;

  // NB: the instance id is not the attribute id / attrId used in the decoding
  // process.  The instance id is used to distinguish between, in the decoded
  // output, multiple attributes with the same label.  Eg, rgb0 and rgb1.
  int attr_instance_id;

  int bitdepth;

  AttributeLabel attributeLabel;

  AttributeParameters params;
};

//============================================================================

struct AttributeParamInventoryHdr {
  int attr_param_seq_parameter_set_id;

  // LSBs of FrameCtr used to identify the frame the parameters apply to
  int attr_param_frame_ctr_lsb;

  // The attribute index in the SPS.
  int attr_param_sps_attr_idx;
};

//---------------------------------------------------------------------------

struct AttributeParamInventory : public AttributeParamInventoryHdr {
  AttributeParameters params;
};

//============================================================================

struct ProfileCompatibility {
  // indicates conformance with the "main" profile
  bool main_profile_compatibility_flag;

  // reserved for future use
  int reserved_profile_compatibility_21bits;

  // indicates that the bistream may break if slices are reordered
  bool slice_reordering_constraint_flag;

  // indicates that there are no duplicate points in the reconstructed frames
  bool unique_point_positions_constraint_flag;

  // during development, no profile bits are set
  bool isDraftProfile() const
  {
    return main_profile_compatibility_flag == 0
      && reserved_profile_compatibility_21bits == 0;
  }
};

//============================================================================

enum class ScaleUnit : bool
{
  kDimensionless = 0,
  kMetre = 1,
};

//============================================================================

struct SequenceParameterSet {
  int sps_seq_parameter_set_id;

  ProfileCompatibility profile;
  int level;

  // Number of bits used to code seqBoundingBoxOrigin
  int sps_bounding_box_offset_bits;

  // the bounding box origin (in stv axis order).
  Vec3<int> seqBoundingBoxOrigin;

  // Number of bits used to code seqBoundingBoxSize
  int sps_bounding_box_size_bits;

  // the size of the bounding box (in stv axis order).
  Vec3<int> seqBoundingBoxSize;

  // A value describing the scaling of the source positions prior to encoding.
  Rational seqGeomScale;

  // Indicates that units used to interpret seqGeomScale.
  ScaleUnit seq_geom_scale_unit_flag;

  // Represents a coded factorisation of a rational:
  //  (2^denominatorLog2 + numeratorModDenominator) * 2^numeratorMulLog2
  //  ------------------------------------------------------------------
  //                       2^denominatorLog2
  struct GlobalScale {
    int numeratorMulLog2 = 0;
    int numeratorModDenominator = 0;
    int denominatorLog2 = 0;

    // Initialized to 1/1
    GlobalScale() = default;

    // Convert rational to global scale representation.
    // NB: may throw an exception if not possible
    GlobalScale(Rational x);

    // Convert to simplified Rational form
    operator Rational() const;
  };

  // Scale factor applied to the coded geometry to get the output geometry
  GlobalScale globalScale;

  int& global_scale_mul_log2() { return globalScale.numeratorMulLog2; }
  int& global_scale_rem() { return globalScale.numeratorModDenominator; }
  int& global_scale_fp_bits() { return globalScale.denominatorLog2; }

  int global_scale_mul_log2() const { return globalScale.numeratorMulLog2; }
  int global_scale_rem() const { return globalScale.numeratorModDenominator; }
  int global_scale_fp_bits() const { return globalScale.denominatorLog2; }

  // NB: attributeSets.size() = num_attribute_sets
  std::vector<AttributeDescription> attributeSets;

  // The number of bits to use for frame_ctr
  int frame_ctr_bits;

  // The number of bits to use for slice_tag
  int slice_tag_bits;

  // Defines the ordering of the position components (eg, xyz vs zyx)
  AxisOrder geometry_axis_order;

  // Controls whether bypass bins are written to a seperate sub-stream, or
  // encoded as ep bins via CABAC.
  bool cabac_bypass_stream_enabled_flag;

  // Indicates that context state may be propagated between slices.
  bool entropy_continuation_enabled_flag;
};

//============================================================================

struct GeometryParameterSet {
  int gps_geom_parameter_set_id;
  int gps_seq_parameter_set_id;

  // Indicates the presence of gps_geom_box_log2_scale and
  // geom_box_log2_scale.
  bool geom_box_log2_scale_present_flag;

  // Default scaling factor for per-slice geometry box origin
  int gps_geom_box_log2_scale;

  // Selects between predictive and octree geometry coding methods.
  bool predgeom_enabled_flag;

  // Controls the ability to represent multiple points (with associated
  // attributes) at the same spatial position.
  bool geom_unique_points_flag;

  // Defines the size of the neighbour availiability volume (aka
  // look-ahead cube size) for occupancy searches.  A value of 0
  // indicates that only neighbours that are direct siblings are available.
  int neighbour_avail_boundary_log2_minus1;

  // Controls the use of early termination of the geometry tree
  // by directly coding the position of isolated points.
  int inferred_direct_coding_mode;

  // Permits coding the common prefix of two idcm points
  bool joint_2pt_idcm_enabled_flag;

  // Selects between bitwise and bytewise occupancy coding
  bool bitwise_occupancy_coding_flag;

  // Controlls contextualization of occupancy bits and refinement of
  // the neighbour pattern according to the occupancy of adjacent
  // children in neighbouring nodes.
  bool adjacent_child_contextualization_enabled_flag;

  // Maximum node size where intra prediction is enabled
  int intra_pred_max_node_size_log2;

  // Enables trisoup
  bool trisoup_enabled_flag;

  // sampling value of trisoup decoding process
  // a value of zero set the automatic sampling value setting to avoid over point of slice MAX points(sliceMaxPoints)
  int trisoup_sampling_value;

  // controls the ability to perform in-loop geometry scaling
  bool geom_scaling_enabled_flag;

  // factor by which to shift geometry QPs before use
  int geom_qp_multiplier_log2;

  // intial qp for geometry scaling, scaled by the qp multiplier
  int geom_base_qp;

  // initial qp (offset) for idcm nodes, scaled by the qp multiplier
  int geom_idcm_qp_offset;

  // Enables/disables non-cubic geometry nodes
  bool qtbt_enabled_flag;

  // Controls the use of planar mode
  bool geom_planar_mode_enabled_flag;
  int geom_planar_threshold0;
  int geom_planar_threshold1;
  int geom_planar_threshold2;
  int geom_idcm_rate_minus1;

  // Enables angular coding in octree/predgeom
  bool geom_angular_mode_enabled_flag;

  // Indicates whether the angualed origin is signalled in the gps or slice
  bool geom_slice_angular_origin_present_flag;

  // Origin for angular mode computations relative to sequence bounding box
  // (in stv axis order).
  Vec3<int> gpsAngularOrigin;

  int numLasers() const { return angularTheta.size(); }

  std::vector<int> angularTheta;
  std::vector<int> angularZ;
  std::vector<int> angularNumPhiPerTurn;

  int geomAngularThetaPred(int i) const
  {
    if (!--i)
      return angularTheta[i];
    return 2 * angularTheta[i] - angularTheta[i - 1];
  }

  // disable the use of planar buffer when angular mode is enabled
  bool planar_buffer_disabled_flag;

  // block size (i.e. number of points per block) in predictive geometry coding
  int geom_qp_offset_intvl_log2;

  // scale factor for azimuth in coding predictive geometry coding
  int geom_angular_azimuth_scale_log2_minus11;
  int geom_angular_azimuth_speed_minus1;

  // inverse scale factor for radius coding in predictive geometry coding
  int geom_angular_radius_inv_scale_log2;

  // Indicates that the geometry footer contains a count of point
  // in each octree level.
  bool octree_point_count_list_present_flag;

  // Quantize azimuth residuals in algular predictive geometry coder.
  bool azimuth_scaling_enabled_flag;
};

//============================================================================

struct GeometryBrickFooter {
  // The actual number of points present in the slice
  int geom_num_points_minus1;

  // The number of points that can be decoded at a particular octree level
  std::vector<int> octree_lvl_num_points_minus1;
};

//============================================================================

struct GeometryBrickHeader {
  int geom_geom_parameter_set_id;
  int slice_tag;
  int geom_slice_id;

  // The LSBs of the system frame counter
  int frame_ctr_lsb;

  // Origin of the reconstructed geometry, relative to sequence bounding box
  // (in stv axis order).
  Vec3<int> geomBoxOrigin;
  int geom_box_log2_scale;

  // Number of bits to represent geomBoxOrigin >> geom_box_log2_scale
  int geom_box_origin_bits_minus1;

  // Per-slice version of gpsAngularOrigin, but relative to slice bounding box
  // (in stv axis order).
  Vec3<int> gbhAngularOrigin;

  // Origin for angular mode computations relative to slice bounding box
  Vec3<int> geomAngularOrigin(const GeometryParameterSet& gps) const
  {
    if (gps.geom_slice_angular_origin_present_flag)
      return gbhAngularOrigin;
    return gps.gpsAngularOrigin - geomBoxOrigin;
  }

  // the size of the root geometry node
  // NB: this is only needed for the initial node size determination at
  //     the encoder
  Vec3<int> rootNodeSizeLog2;

  Vec3<int> pgeom_resid_abs_log2_bits;

  // the largest dimension of the root geometry node
  mutable int maxRootNodeDimLog2;

  std::vector<int8_t> tree_lvl_coded_axis_list;

  int tree_depth_minus1() const { return tree_lvl_coded_axis_list.size() - 1; }

  // qp offset for geometry scaling (if enabled)
  int geom_slice_qp_offset;

  int sliceQp(const GeometryParameterSet& gps) const
  {
    return (gps.geom_base_qp + geom_slice_qp_offset)
      << gps.geom_qp_multiplier_log2;
  }

  // block size offset for predictive geometry coding (if enabled)
  int geom_qp_offset_intvl_log2_delta;

  // number of entropy streams used to encode the octree
  int geom_stream_cnt_minus1;

  int geomBoxLog2Scale(const GeometryParameterSet& gps) const
  {
    if (!gps.geom_box_log2_scale_present_flag)
      return gps.gps_geom_box_log2_scale;
    return geom_box_log2_scale;
  }

  // size of triangle nodes (reconstructed surface) in trisoup geometry.
  int trisoup_node_size_log2_minus2;

  int trisoupNodeSizeLog2(const GeometryParameterSet& gps) const
  {
    return gps.trisoup_enabled_flag ? trisoup_node_size_log2_minus2 + 2 : 0;
  }

  // downsampling rate used in tringle voxelisation
  int trisoup_sampling_value_minus1;

  int num_unique_segments_minus1;

  // Number of bits to represent num_unique_segments_minus1
  int num_unique_segments_bits_minus1;

  // 'Header' information that appears at the end of the data unit
  GeometryBrickFooter footer;

  // Indicates the current slice reuses contexts from the prevous slice
  bool entropy_continuation_flag;

  // The id of the previous slice in bitsream order
  int prev_slice_id;

  // minimum radius for predictive geometry coding with angular mode
  int pgeom_min_radius;
};

//============================================================================
// NB: when updating this, remember to update AttributeLods::isReusable(...)

struct AttributeParameterSet {
  int aps_attr_parameter_set_id;
  int aps_seq_parameter_set_id;
  AttributeEncoding attr_encoding;

  bool lodParametersPresent() const
  {
    return attr_encoding == AttributeEncoding::kLiftingTransform
      || attr_encoding == AttributeEncoding::kPredictingTransform;
  }

  //--- lifting/predicting transform parameters

  LodDecimationMethod lod_decimation_type;
  bool canonical_point_order_flag;
  int num_pred_nearest_neighbours_minus1;
  int max_num_direct_predictors;
  bool direct_avg_predictor_disabled_flag;
  int adaptive_prediction_threshold;
  int intra_lod_search_range;
  int inter_lod_search_range;

  // Neighbour contribution weights used to calculate quantization weights
  Vec3<uint32_t> quant_neigh_weight;

  int adaptivePredictionThreshold(const AttributeDescription& desc) const
  {
    return adaptive_prediction_threshold << std::max(0, desc.bitdepth - 8);
  }

  // NB: in stv order
  Vec3<int32_t> lodNeighBias;

  // The number of detail levels, starting from the finest, to skip intra
  // prediction.
  int intra_lod_prediction_skip_layers;

  // A large value for intra_lod_prediction_skip_layers that causes all layers
  // to be skipped.
  static const int kSkipAllLayers = 0x7fffffff;

  bool inter_component_prediction_enabled_flag;
  bool last_component_prediction_enabled_flag;

  // Whether prediction weights are blended according to spatial positions
  bool pred_weight_blending_enabled_flag;

  // The number of refinement layers
  int num_detail_levels_minus1;

  // The configured number of detail levels (there is always at least one).
  // NB: the actual number of detail levels generated may be lower.
  int maxNumDetailLevels() const
  {
    return scalable_lifting_enabled_flag ? 21 /* MaxRootNodeDimLog2 */
                                         : num_detail_levels_minus1 + 1;
  }

  std::vector<int> lodSamplingPeriod;

  int dist2;
  bool aps_slice_dist2_deltas_present_flag;

  // NB: these parameters are shared by all transform implementations
  int init_qp_minus4;
  int aps_chroma_qp_offset;
  bool aps_slice_qp_deltas_present_flag;

  //--- raht parameters
  bool raht_prediction_enabled_flag;
  int raht_prediction_threshold0;
  int raht_prediction_threshold1;

  //--- lifting parameters
  bool scalable_lifting_enabled_flag;
  int max_neigh_range_minus1;

  // indicates that attribute coding should be performed in
  // pseudo-spherical domain
  bool spherical_coord_flag;

  // (r, phi, laserid) scale factors for domain conversion
  Vec3<int> attr_coord_scale;

  // Whether raw attribute are coded as fixed width or variable length.
  bool raw_attr_variable_len_flag;
};

//============================================================================

struct AttributeBrickHeader {
  int attr_sps_attr_idx;
  int attr_attr_parameter_set_id;
  int attr_geom_slice_id;

  // Last component prediction coefficients.  Only present for lifting
  // transform with three components.
  std::vector<int8_t> attrLcpCoeffs;

  // indicates whether last component prediction coefficients are present.
  bool lcpPresent(
    const AttributeDescription& desc, const AttributeParameterSet& aps) const
  {
    if (aps.attr_encoding != AttributeEncoding::kLiftingTransform)
      return false;
    if (!aps.last_component_prediction_enabled_flag)
      return false;
    if (desc.attr_num_dimensions_minus1 != 2)
      return false;
    return true;
  }

  // Inter-component prediction coefficients.
  std::vector<Vec3<int8_t>> icpCoeffs;

  // indicates whether inter component prediction coefficients are present.
  bool icpPresent(
    const AttributeDescription& desc, const AttributeParameterSet& aps) const
  {
    if (aps.attr_encoding != AttributeEncoding::kPredictingTransform)
      return false;
    if (!aps.inter_component_prediction_enabled_flag)
      return false;
    if (desc.attr_num_dimensions_minus1 == 0)
      return false;
    return true;
  }

  int attr_qp_delta_luma;
  int attr_qp_delta_chroma;

  std::vector<int> attr_layer_qp_delta_luma;
  std::vector<int> attr_layer_qp_delta_chroma;

  bool attr_layer_qp_present_flag() const
  {
    return !attr_layer_qp_delta_luma.empty();
  }

  int attr_num_qp_layers_minus1() const
  {
    return attr_layer_qp_delta_luma.size() - 1;
  }

  struct QpRegion {
    // NB: in stv order
    Vec3<int> regionOrigin;

    // NB: in stv order
    Vec3<int> regionSize;

    std::array<int, 2> attr_region_qp_offset;
  };

  std::vector<QpRegion> qpRegions;

  // Number of bits to represent regionOrigin and regionSize
  int attr_region_bits_minus1;

  int32_t attr_dist2_delta;
};

//============================================================================

struct ConstantAttributeDataUnit {
  int constattr_sps_attr_idx;
  int constattr_attr_parameter_set_id;
  int constattr_geom_slice_id;

  std::vector<int> constattr_default_value;
};

//============================================================================

struct FrameBoundaryMarker {
  // Identifies the frame to be terminated
  int fbdu_frame_ctr_lsb;
};

//============================================================================

struct TileInventory {
  struct Entry {
    // The tile id (either manually specified, or the implicit value).
    int tile_id;

    // NB: in stv order
    Vec3<int> tileOrigin;

    // NB: in stv order
    Vec3<int> tileSize;
  };

  // id of an applicable sequence parameter set
  int ti_seq_parameter_set_id;

  // Number of bits for frame idx
  int ti_frame_ctr_bits;

  // Frame idx when tile inventory comes into force
  int ti_frame_ctr;

  // Number of bits, if any, used to signal tile_id
  int tile_id_bits;

  // the origin of the tiles (in stv axis order).  Likely the sps origin
  Vec3<int> origin;

  // Number of bits to represent the inventory origin
  int ti_origin_bits_minus1;

  std::vector<Entry> tiles;

  // Number of bits to represent each tile's origin
  int tile_origin_bits_minus1;

  // Number of bits to represent each tile's size
  int tile_size_bits_minus1;
};

//============================================================================

struct UserData {
  // Identifies the type of user data
  Oid user_data_oid;
};

//============================================================================

void convertXyzToStv(SequenceParameterSet*);
void convertXyzToStv(const SequenceParameterSet&, GeometryParameterSet*);
void convertXyzToStv(const SequenceParameterSet&, AttributeParameterSet*);
void convertXyzToStv(const SequenceParameterSet&, TileInventory*);

//============================================================================

}  // namespace pcc
