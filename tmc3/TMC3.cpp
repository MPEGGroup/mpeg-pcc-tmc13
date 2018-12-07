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

#include "TMC3.h"

#include <memory>

#include "PCCTMC3Encoder.h"
#include "PCCTMC3Decoder.h"
#include "constants.h"
#include "program_options_lite.h"
#include "io_tlv.h"
#include "version.h"

using namespace std;
using namespace pcc;

//============================================================================

struct Parameters {
  bool isDecoder;

  // command line parsing should adjust dist2 values according to PQS
  bool positionQuantizationScaleAdjustsDist2;

  // output mode for ply writing (binary or ascii)
  bool outputBinaryPly;

  // when true, configure the encoder as if no attributes are specified
  bool disableAttributeCoding;

  std::string uncompressedDataPath;
  std::string compressedStreamPath;
  std::string reconstructedDataPath;

  // Filename for saving pre inverse scaled point cloud.
  std::string preInvScalePath;

  pcc::EncoderParams encoder;
  pcc::DecoderParams decoder;

  // todo(df): this should be per-attribute
  ColorTransform colorTransform;

  // todo(df): this should be per-attribute
  int reflectanceScale;
};

//============================================================================

int
main(int argc, char* argv[])
{
  cout << "MPEG PCC tmc3 version " << ::pcc::version << endl;

  Parameters params;
  if (!ParseParameters(argc, argv, params)) {
    return -1;
  }

  // Timers to count elapsed wall/user time
  pcc::chrono::Stopwatch<std::chrono::steady_clock> clock_wall;
  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;

  clock_wall.start();

  int ret = 0;
  if (params.isDecoder) {
    ret = Decompress(params, clock_user);
  } else {
    ret = Compress(params, clock_user);
  }

  clock_wall.stop();

  using namespace std::chrono;
  auto total_wall = duration_cast<milliseconds>(clock_wall.count()).count();
  auto total_user = duration_cast<milliseconds>(clock_user.count()).count();
  std::cout << "Processing time (wall): " << total_wall / 1000.0 << " s\n";
  std::cout << "Processing time (user): " << total_user / 1000.0 << " s\n";

  return ret;
}

//---------------------------------------------------------------------------
// :: Command line / config parsing helpers

template<typename T>
static std::istream&
readUInt(std::istream& in, T& val)
{
  unsigned int tmp;
  in >> tmp;
  val = T(tmp);
  return in;
}

static std::istream&
operator>>(std::istream& in, ColorTransform& val)
{
  return readUInt(in, val);
}

namespace pcc {
static std::istream&
operator>>(std::istream& in, AttributeEncoding& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, PartitionMethod& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const AttributeEncoding& val)
{
  switch (val) {
  case AttributeEncoding::kPredictingTransform: out << "0 (Pred)"; break;
  case AttributeEncoding::kRAHTransform: out << "1 (RAHT)"; break;
  case AttributeEncoding::kLiftingTransform: out << "2 (Lift)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const PartitionMethod& val)
{
  switch (val) {
  case PartitionMethod::kNone: out << "0 (None)"; break;
  case PartitionMethod::kUniformGeom: out << "0 (UniformGeom)"; break;
  case PartitionMethod::kOctreeUniform: out << "0 (UniformOctree)"; break;
  default: out << int(val) << " (Unknown)"; break;
  }
  return out;
}
}  // namespace pcc

namespace df {
namespace program_options_lite {
  template<typename T>
  struct option_detail<pcc::PCCVector3<T>> {
    static constexpr bool is_container = true;
    static constexpr bool is_fixed_size = true;
    typedef T* output_iterator;

    static void clear(pcc::PCCVector3<T>& container){};
    static output_iterator make_output_iterator(pcc::PCCVector3<T>& container)
    {
      return &container[0];
    }
  };
}  // namespace program_options_lite
}  // namespace df

//---------------------------------------------------------------------------
// :: Command line / config parsing

bool
ParseParameters(int argc, char* argv[], Parameters& params)
{
  namespace po = df::program_options_lite;

  struct {
    AttributeDescription desc;
    AttributeParameterSet aps;
  } params_attr;

  bool print_help = false;

  // a helper to set the attribute
  std::function<po::OptionFunc::Func> attribute_setter =
    [&](po::Options&, const std::string& name, po::ErrorReporter) {
      // copy the current state of parsed attribute parameters
      //
      // NB: this does not cause the default values of attr to be restored
      // for the next attribute block.  A side-effect of this is that the
      // following is allowed leading to attribute foo having both X=1 and
      // Y=2:
      //   "--attr.X=1 --attribute foo --attr.Y=2 --attribute foo"
      //

      // NB: insert returns any existing element
      const auto& it = params.encoder.attributeIdxMap.insert(
        {name, int(params.encoder.attributeIdxMap.size())});

      if (it.second) {
        params.encoder.sps.attributeSets.push_back(params_attr.desc);
        params.encoder.aps.push_back(params_attr.aps);
        return;
      }

      // update existing entry
      params.encoder.sps.attributeSets[it.first->second] = params_attr.desc;
      params.encoder.aps[it.first->second] = params_attr.aps;
    };

  /* clang-format off */
  // The definition of the program/config options, along with default values.
  //
  // NB: when updating the following tables:
  //      (a) please keep to 80-columns for easier reading at a glance,
  //      (b) do not vertically align values -- it breaks quickly
  //
  po::Options opts;
  opts.addOptions()
  ("help", print_help, false, "this help text")
  ("config,c", po::parseConfigFile, "configuration file name")

  (po::Section("General"))

  ("mode", params.isDecoder, false,
    "The encoding/decoding mode:\n"
    "  0: encode\n"
    "  1: decode")

  // i/o parameters
  ("reconstructedDataPath",
    params.reconstructedDataPath, {},
    "The ouput reconstructed pointcloud file path (decoder only)")

  ("uncompressedDataPath",
    params.uncompressedDataPath, {},
    "The input pointcloud file path")

  ("compressedStreamPath",
    params.compressedStreamPath, {},
    "The compressed bitstream path (encoder=output, decoder=input)")

  ("postRecolorPath",
    params.encoder.postRecolorPath, {},
    "Recolored pointcloud file path (encoder only)")

  ("preInvScalePath",
    params.preInvScalePath, {},
    "Pre inverse scaled pointcloud file path (decoder only)")

  ("outputBinaryPly",
    params.outputBinaryPly, false,
    "Output ply files using binary (or otherwise ascii) format")

  // general
  // todo(df): this should be per-attribute
  ("colorTransform",
    params.colorTransform, COLOR_TRANSFORM_RGB_TO_YCBCR,
    "The colour transform to be applied:\n"
    "  0: none\n"
    "  1: RGB to YCbCr (Rec.709)")

  // todo(df): this should be per-attribute
  ("hack.reflectanceScale",
    params.reflectanceScale, 1,
    "scale factor to be applied to reflectance "
    "pre encoding / post reconstruction")

  // NB: if adding decoder options, uncomment the Decoder section marker
  // (po::Section("Decoder"))

  (po::Section("Encoder"))

  ("seq_bounding_box_xyz0",
    params.encoder.sps.seq_bounding_box_xyz0, {0},
    "seq_bounding_box_xyz0.  NB: seq_bounding_box_whd must be set for this "
    "parameter to have an effect")

  ("seq_bounding_box_whd",
    params.encoder.sps.seq_bounding_box_whd, {0},
    "seq_bounding_box_whd")

  ("positionQuantizationScale",
    params.encoder.sps.seq_source_geom_scale_factor, 1.f,
    "Scale factor to be applied to point positions during quantization process")

  ("positionQuantizationScaleAdjustsDist2",
    params.positionQuantizationScaleAdjustsDist2, false,
    "Scale dist2 values by squared positionQuantizationScale")

  ("mergeDuplicatedPoints",
    params.encoder.gps.geom_unique_points_flag, true,
    "Enables removal of duplicated points")

  ("partitionMethod",
    params.encoder.partitionMethod, PartitionMethod::kNone,
    "Method used to partition input point cloud into slices/tiles:\n"
    "  0: none\n"
    "  1: none (deprecated)\n"
    "  2: n Uniform-Geometry partition bins along the longest edge\n"
    "  3: Uniform Geometry partition at n octree depth")

  ("partitionNumUniformGeom",
    params.encoder.partitionNumUniformGeom, 0,
    "Number of bins for partitionMethod=2:\n"
    "  0: slice partition with adaptive-defined bins\n"
    "  >=1: slice partition with user-defined bins\n")

  ("partitionOctreeDepth",
    params.encoder.partitionOctreeDepth, 2,
    "Depth of octree partition for partitionMethod=3")

  ("disableAttributeCoding",
    params.disableAttributeCoding, false,
    "Ignore attribute coding configuration")

  (po::Section("Geometry"))

  // tools
  ("bitwiseOccupancyCoding",
    params.encoder.gps.bitwise_occupancy_coding_flag, true,
    "Selects between bitwise and bytewise occupancy coding:\n"
    "  0: bytewise\n"
    "  1: bitwise")

  ("neighbourContextRestriction",
    params.encoder.gps.neighbour_context_restriction_flag, false,
    "Limit geometry octree occupancy contextualisation to sibling nodes")

  ("neighbourAvailBoundaryLog2",
    params.encoder.gps.neighbour_avail_boundary_log2, 0,
    "Defines the avaliability volume for neighbour occupancy lookups."
    " 0: unconstrained")

  ("inferredDirectCodingMode",
    params.encoder.gps.inferred_direct_coding_mode_enabled_flag, true,
    "Permits early termination of the geometry octree for isolated points")

  ("intra_pred_max_node_size_log2",
    params.encoder.gps.intra_pred_max_node_size_log2, 0,
    "octree nodesizes eligible for occupancy intra prediction")

  ("ctxOccupancyReductionFactor",
     params.encoder.gps.geom_occupancy_ctx_reduction_factor, 3,
     "Adjusts the number of contexts used in occupancy coding")

  ("trisoup_node_size_log2",
    params.encoder.gps.trisoup_node_size_log2, 0,
    "Size of nodes for surface triangulation.\n"
    "  0: disabled\n")

  (po::Section("Attributes"))

  // attribute processing
  //   NB: Attribute options are special in the way they are applied (see above)
  ("attribute",
    attribute_setter,
    "Encode the given attribute (NB, must appear after the"
    "following attribute parameters)")

  ("bitdepth",
    params_attr.desc.attr_bitdepth, 8,
    "Attribute bitdepth")

  ("transformType",
    params_attr.aps.attr_encoding, AttributeEncoding::kPredictingTransform,
    "Coding method to use for attribute:\n"
    "  0: Hierarchical neighbourhood prediction\n"
    "  1: Region Adaptive Hierarchical Transform (RAHT)\n"
    "  2: Hierarichical neighbourhood prediction as lifting transform")

  ("rahtLeafDecimationDepth",
    params_attr.aps.raht_binary_level_threshold, 3,
    "Sets coefficients to zero in the bottom n levels of RAHT tree. "
    "Used for chroma-subsampling in attribute=color only.")

  ("rahtQuantizationStep",
    params_attr.aps.quant_step_size_luma, 0,
    "deprecated -- use quantizationStepsLuma")

  ("rahtDepth",
    params_attr.aps.raht_depth, 21,
    "Number of bits for morton representation of an RAHT co-ordinate"
    "component")

  ("numberOfNearestNeighborsInPrediction",
    params_attr.aps.num_pred_nearest_neighbours, 3,
    "Attribute's maximum number of nearest neighbors to be used for prediction")

  ("adaptivePredictionThreshold",
    params_attr.aps.adaptive_prediction_threshold, -1,
    "Neighbouring attribute value difference that enables choice of "
    "single|multi predictors. Applies to transformType=2 only.\n"
    "  -1: auto = 2**(bitdepth-2)")

  ("attributeSearchRange",
    params_attr.aps.search_range, 128,
    "Range for nearest neighbor search")

  ("lodBinaryTree",
    params_attr.aps.lod_binary_tree_enabled_flag, false,
    "Controls LoD generation method:\n"
    " 0: distance based subsampling\n"
    " 1: binary tree")

  ("max_num_direct_predictors",
    params_attr.aps.max_num_direct_predictors, 3,
    "Maximum number of nearest neighbour candidates used in direct"
    "attribute prediction")

  ("levelOfDetailCount",
    params_attr.aps.num_detail_levels, 1,
    "Attribute's number of levels of detail")

  ("quantizationStepLuma",
    params_attr.aps.quant_step_size_luma, 0,
    "Attribute's luma quantization step size")

  ("quantizationStepChroma",
    params_attr.aps.quant_step_size_chroma, 0,
    "Attribute's chroma quantization step size")

  ("dist2",
    params_attr.aps.dist2, {},
    "Attribute's list of squared distances, or initial value for automatic"
    "derivation")
  ;
  /* clang-format on */

  po::setDefaults(opts);
  po::ErrorReporter err;
  const list<const char*>& argv_unhandled =
    po::scanArgv(opts, argc, (const char**)argv, err);

  for (const auto arg : argv_unhandled) {
    err.warn() << "Unhandled argument ignored: " << arg << "\n";
  }

  if (argc == 1 || print_help) {
    po::doHelp(std::cout, opts, 78);
    return false;
  }

  // Certain coding modes are not available when trisoup is enabled.
  // Disable them, and warn if set (they may be set as defaults).
  if (params.encoder.gps.trisoup_node_size_log2 > 0) {
    if (!params.encoder.gps.geom_unique_points_flag)
      err.warn() << "TriSoup geometry does not preserve duplicated points\n";

    if (params.encoder.gps.inferred_direct_coding_mode_enabled_flag)
      err.warn() << "TriSoup geometry is incompatable with IDCM\n";

    params.encoder.gps.geom_unique_points_flag = true;
    params.encoder.gps.inferred_direct_coding_mode_enabled_flag = false;
  }

  // support disabling attribute coding (simplifies configuration)
  if (params.disableAttributeCoding) {
    params.encoder.attributeIdxMap.clear();
    params.encoder.sps.attributeSets.clear();
    params.encoder.aps.clear();
  }

  // fixup any per-attribute settings
  for (const auto& it : params.encoder.attributeIdxMap) {
    auto& attr_sps = params.encoder.sps.attributeSets[it.second];
    auto& attr_aps = params.encoder.aps[it.second];

    // Avoid wasting bits signalling chroma quant step size for reflectance
    if (it.first == "reflectance") {
      attr_aps.quant_step_size_chroma = 0;
    }

    bool isLifting =
      attr_aps.attr_encoding == AttributeEncoding::kPredictingTransform
      || attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform;

    // derive the dist2 values based on an initial value
    if (isLifting) {
      if (attr_aps.dist2.size() > attr_aps.num_detail_levels) {
        attr_aps.dist2.resize(attr_aps.num_detail_levels);
      } else if (
        attr_aps.dist2.size() < attr_aps.num_detail_levels
        && !attr_aps.dist2.empty()) {
        if (attr_aps.dist2.size() < attr_aps.num_detail_levels) {
          attr_aps.dist2.resize(attr_aps.num_detail_levels);
          const double distRatio = 4.0;
          uint64_t d2 = attr_aps.dist2[0];
          for (int i = 0; i < attr_aps.num_detail_levels; ++i) {
            attr_aps.dist2[i] = d2;
            d2 = uint64_t(std::round(distRatio * d2));
          }
        }
      }
    }
    // In order to simplify specification of dist2 values, which are
    // depending on the scale of the coded point cloud, the following
    // adjust the dist2 values according to PQS.  The user need only
    // specify the unquantised PQS value.
    if (params.positionQuantizationScaleAdjustsDist2) {
      double pqs = params.encoder.sps.seq_source_geom_scale_factor;
      double pqs2 = pqs * pqs;
      for (auto& dist2 : attr_aps.dist2)
        dist2 = int64_t(std::round(pqs2 * dist2));
    }

    // Set default threshold based on bitdepth
    if (attr_aps.adaptive_prediction_threshold == -1) {
      attr_aps.adaptive_prediction_threshold = 1
        << (attr_sps.attr_bitdepth - 2);
    }

    if (attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform) {
      attr_aps.adaptive_prediction_threshold = 0;
    }

    // For RAHT, ensure that the unused lod count = 0 (prevents mishaps)
    if (attr_aps.attr_encoding == AttributeEncoding::kRAHTransform) {
      attr_aps.num_detail_levels = 0;
      attr_aps.adaptive_prediction_threshold = 0;

      // todo(df): suggest chroma quant_step_size for raht
      attr_aps.quant_step_size_chroma = 0;
    }
  }

  // sanity checks

  if (params.encoder.gps.intra_pred_max_node_size_log2)
    if (!params.encoder.gps.neighbour_avail_boundary_log2)
      err.error() << "Geometry intra prediction requires finite"
                     "neighbour_avail_boundary_log2\n";

  for (const auto& it : params.encoder.attributeIdxMap) {
    const auto& attr_sps = params.encoder.sps.attributeSets[it.second];
    const auto& attr_aps = params.encoder.aps[it.second];

    bool isLifting =
      attr_aps.attr_encoding == AttributeEncoding::kPredictingTransform
      || attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform;

    if (it.first == "color") {
      // todo(??): permit relaxing of the following constraint
      if (attr_sps.attr_bitdepth > 8)
        err.error() << it.first << ".bitdepth must be less than 9\n";
    }

    if (it.first == "reflectance") {
      if (attr_sps.attr_bitdepth > 16)
        err.error() << it.first << ".bitdepth must be less than 17\n";
    }

    if (isLifting) {
      int lod = attr_aps.num_detail_levels;
      if (lod > 255 || lod < 0) {
        err.error() << it.first
                    << ".levelOfDetailCount must be in the range [0,255]\n";
      }
      if (attr_aps.dist2.size() != lod) {
        err.error() << it.first << ".dist2 does not have " << lod
                    << " entries\n";
      }

      if (attr_aps.adaptive_prediction_threshold < 0) {
        err.error() << it.first
                    << ".adaptivePredictionThreshold must be positive\n";
      }

      if (
        attr_aps.num_pred_nearest_neighbours
        > kAttributePredictionMaxNeighbourCount) {
        err.error() << it.first
                    << ".numberOfNearestNeighborsInPrediction must be <= "
                    << kAttributePredictionMaxNeighbourCount << "\n";
      }
    }
  }

  // check required arguments are specified

  if (!params.isDecoder && params.uncompressedDataPath.empty())
    err.error() << "uncompressedDataPath not set\n";

  if (params.isDecoder && params.reconstructedDataPath.empty())
    err.error() << "reconstructedDataPath not set\n";

  if (params.compressedStreamPath.empty())
    err.error() << "compressedStreamPath not set\n";

  // report the current configuration (only in the absence of errors so
  // that errors/warnings are more obvious and in the same place).
  if (err.is_errored)
    return false;

  // Dump the complete derived configuration
  cout << "+ Effective configuration parameters\n";

  po::dumpCfg(cout, opts, "General", 4);
  if (params.isDecoder) {
    po::dumpCfg(cout, opts, "Decoder", 4);
  } else {
    po::dumpCfg(cout, opts, "Encoder", 4);
    po::dumpCfg(cout, opts, "Geometry", 4);

    for (const auto& it : params.encoder.attributeIdxMap) {
      // NB: when dumping the config, opts references params_attr
      params_attr.desc = params.encoder.sps.attributeSets[it.second];
      params_attr.aps = params.encoder.aps[it.second];
      cout << "    " << it.first << "\n";
      po::dumpCfg(cout, opts, "Attributes", 8);
    }
  }

  cout << endl;

  return true;
}

int
Compress(Parameters& params, Stopwatch& clock)
{
  PCCPointSet3 pointCloud;
  if (
    !pointCloud.read(params.uncompressedDataPath)
    || pointCloud.getPointCount() == 0) {
    cout << "Error: can't open input file!" << endl;
    return -1;
  }

  // Sanitise the input point cloud
  // todo(df): remove the following with generic handling of properties
  bool codeColour = params.encoder.attributeIdxMap.count("color");
  if (!codeColour)
    pointCloud.removeColors();
  assert(codeColour == pointCloud.hasColors());

  bool codeReflectance = params.encoder.attributeIdxMap.count("reflectance");
  if (!codeReflectance)
    pointCloud.removeReflectances();
  assert(codeReflectance == pointCloud.hasReflectances());

  ofstream fout(params.compressedStreamPath, ios::binary);
  if (!fout.is_open()) {
    return -1;
  }

  clock.start();

  if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
    pointCloud.convertRGBToYUV();
  }

  if (params.reflectanceScale > 1 && pointCloud.hasReflectances()) {
    const auto pointCount = pointCloud.getPointCount();
    for (size_t i = 0; i < pointCount; ++i) {
      int val = pointCloud.getReflectance(i) / params.reflectanceScale;
      pointCloud.setReflectance(i, val);
    }
  }

  PCCTMC3Encoder3 encoder;

  // The reconstructed point cloud
  std::unique_ptr<PCCPointSet3> reconPointCloud;
  if (!params.reconstructedDataPath.empty()) {
    reconPointCloud.reset(new PCCPointSet3);
  }

  int ret = encoder.compress(
    pointCloud, &params.encoder,
    [&](const PayloadBuffer& buf) { writeTlv(buf, fout); },
    reconPointCloud.get());
  if (ret) {
    cout << "Error: can't compress point cloud!" << endl;
    return -1;
  }

  std::cout << "Total bitstream size " << fout.tellp() << " B" << std::endl;
  fout.close();

  clock.stop();

  if (!params.reconstructedDataPath.empty()) {
    if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
      reconPointCloud->convertYUVToRGB();
    }

    if (params.reflectanceScale > 1 && reconPointCloud->hasReflectances()) {
      const auto pointCount = reconPointCloud->getPointCount();
      for (size_t i = 0; i < pointCount; ++i) {
        int val = reconPointCloud->getReflectance(i) * params.reflectanceScale;
        reconPointCloud->setReflectance(i, val);
      }
    }

    reconPointCloud->write(
      params.reconstructedDataPath, !params.outputBinaryPly);
  }

  return 0;
}
int
Decompress(Parameters& params, Stopwatch& clock)
{
  ifstream fin(params.compressedStreamPath, ios::binary);
  if (!fin.is_open()) {
    return -1;
  }

  clock.start();

  PayloadBuffer buf;
  PCCTMC3Decoder3 decoder;

  while (true) {
    PayloadBuffer* buf_ptr = &buf;
    readTlv(fin, &buf);

    // at end of file (or other error), flush decoder
    if (!fin)
      buf_ptr = nullptr;

    int ret = decoder.decompress(
      params.decoder, buf_ptr, [&](const PCCPointSet3& decodedPointCloud) {
        PCCPointSet3 pointCloud(decodedPointCloud);

        if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
          pointCloud.convertYUVToRGB();
        }

        if (params.reflectanceScale > 1 && pointCloud.hasReflectances()) {
          const auto pointCount = pointCloud.getPointCount();
          for (size_t i = 0; i < pointCount; ++i) {
            int val = pointCloud.getReflectance(i) * params.reflectanceScale;
            pointCloud.setReflectance(i, val);
          }
        }

        // Dump the decoded colour using the pre inverse scaled geometry
        if (!params.preInvScalePath.empty()) {
          pointCloud.write(params.preInvScalePath, !params.outputBinaryPly);
        }

        decoder.inverseQuantization(pointCloud);

        clock.stop();

        if (!pointCloud.write(
              params.reconstructedDataPath, !params.outputBinaryPly)) {
          cout << "Error: can't open output file!" << endl;
        }

        clock.start();
      });

    if (ret) {
      cout << "Error: can't decompress point cloud!" << endl;
      return -1;
    }

    if (!buf_ptr)
      break;
  }

  fin.clear();
  fin.seekg(0, ios_base::end);
  std::cout << "Total bitstream size " << fin.tellg() << " B" << std::endl;

  clock.stop();

  return 0;
}
