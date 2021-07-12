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
#include "ply.h"
#include "pointset_processing.h"
#include "program_options_lite.h"
#include "io_tlv.h"
#include "version.h"

using namespace std;
using namespace pcc;

//============================================================================

enum class OutputSystem
{
  // Output after global scaling, don't convert to external system
  kConformance = 0,

  // Scale output to external coordinate system
  kExternal = 1,
};

//----------------------------------------------------------------------------

struct Parameters {
  bool isDecoder;

  // command line parsing should adjust dist2 values according to PQS
  bool positionQuantizationScaleAdjustsDist2;

  // Scale factor to apply when loading the ply before integer conversion.
  // Eg, If source point positions are in fractional metres converting to
  // millimetres will allow some fidelity to be preserved.
  double inputScale;

  // Length of the output point clouds unit vectors.
  double outputUnitLength;

  // output mode for ply writing (binary or ascii)
  bool outputBinaryPly;

  // Fractional fixed-point bits retained in conformance output
  int outputFpBits;

  // Output coordinate system to use
  OutputSystem outputSystem;

  // when true, configure the encoder as if no attributes are specified
  bool disableAttributeCoding;

  // Frame number of first file in input sequence.
  int firstFrameNum;

  // Number of frames to process.
  int frameCount;

  std::string uncompressedDataPath;
  std::string compressedStreamPath;
  std::string reconstructedDataPath;

  // Filename for saving recoloured point cloud (encoder).
  std::string postRecolorPath;

  // Filename for saving pre inverse scaled point cloud (decoder).
  std::string preInvScalePath;

  pcc::EncoderParams encoder;
  pcc::DecoderParams decoder;

  // perform attribute colourspace conversion on ply input/output.
  bool convertColourspace;

  // resort the input points by azimuth angle
  bool sortInputByAzimuth;
};

//----------------------------------------------------------------------------

class SequenceCodec {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceCodec(Parameters* params) : params(params) {}

  // Perform conversions and write output point cloud
  //  \params cloud  a mutable copy of reconFrame.cloud
  void writeOutputFrame(
    const std::string& postInvScalePath,
    const std::string& preInvScalePath,
    const CloudFrame& reconFrame,
    PCCPointSet3& cloud);

  // determine the output ply scale factor
  double outputScale(const CloudFrame& cloud) const;

  // the output ply origin, scaled according to output coordinate system
  Vec3<double> outputOrigin(const CloudFrame& cloud) const;

  void scaleAttributesForInput(
    const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

  void scaleAttributesForOutput(
    const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

protected:
  Parameters* params;
};

//----------------------------------------------------------------------------

class SequenceEncoder
  : public SequenceCodec
  , PCCTMC3Encoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceEncoder(Parameters* params);

  int compress(Stopwatch* clock);

protected:
  int compressOneFrame(Stopwatch* clock);

  void onOutputBuffer(const PayloadBuffer& buf) override;
  void onPostRecolour(const PCCPointSet3& cloud) override;

private:
  ply::PropertyNameMap _plyAttrNames;

  // The raw origin used for input sorting
  Vec3<int> _angularOrigin;

  PCCTMC3Encoder3 encoder;

  std::ofstream bytestreamFile;

  int frameNum;
};

//----------------------------------------------------------------------------

class SequenceDecoder
  : public SequenceCodec
  , PCCTMC3Decoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceDecoder(Parameters* params);

  int decompress(Stopwatch* clock);

protected:
  void onOutputCloud(const CloudFrame& cloud) override;

private:
  PCCTMC3Decoder3 decoder;

  std::ofstream bytestreamFile;

  Stopwatch* clock;
};

//============================================================================

void convertToGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

void convertFromGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

//============================================================================

int
main(int argc, char* argv[])
{
  cout << "MPEG PCC tmc3 version " << ::pcc::version << endl;

  Parameters params;

  try {
    if (!ParseParameters(argc, argv, params))
      return 1;
  }
  catch (df::program_options_lite::ParseFailure& e) {
    std::cerr << "Error parsing option \"" << e.arg << "\" with argument \""
              << e.val << "\"." << std::endl;
    return 1;
  }

  // Timers to count elapsed wall/user time
  pcc::chrono::Stopwatch<std::chrono::steady_clock> clock_wall;
  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;

  clock_wall.start();

  int ret = 0;
  if (params.isDecoder) {
    ret = SequenceDecoder(&params).decompress(&clock_user);
  } else {
    ret = SequenceEncoder(&params).compress(&clock_user);
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

std::array<const char*, 3>
axisOrderToPropertyNames(AxisOrder order)
{
  static const std::array<const char*, 3> kAxisOrderToPropertyNames[] = {
    {"z", "y", "x"}, {"x", "y", "z"}, {"x", "z", "y"}, {"y", "z", "x"},
    {"z", "y", "x"}, {"z", "x", "y"}, {"y", "x", "z"}, {"x", "y", "z"},
  };

  return kAxisOrderToPropertyNames[int(order)];
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

namespace pcc {
static std::istream&
operator>>(std::istream& in, ScaleUnit& val)
{
  try {
    readUInt(in, val);
  }
  catch (...) {
    in.clear();
    std::string str;
    in >> str;

    val = ScaleUnit::kDimensionless;
    if (str == "metre")
      val = ScaleUnit::kMetre;
    else if (!str.empty())
      throw std::runtime_error("Cannot parse unit");
  }
  return in;
}
}  // namespace pcc

static std::istream&
operator>>(std::istream& in, OutputSystem& val)
{
  return readUInt(in, val);
}

namespace pcc {
static std::istream&
operator>>(std::istream& in, ColourMatrix& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, AxisOrder& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, AttributeEncoding& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, LodDecimationMethod& val)
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
static std::istream&
operator>>(std::istream& in, PredGeomEncOpts::SortMode& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, OctreeEncOpts::QpMethod& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

static std::ostream&
operator<<(std::ostream& out, const OutputSystem& val)
{
  switch (val) {
  case OutputSystem::kConformance: out << "0 (Conformance)"; break;
  case OutputSystem::kExternal: out << "1 (External)"; break;
  }
  return out;
}

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const ScaleUnit& val)
{
  switch (val) {
  case ScaleUnit::kDimensionless: out << "0 (Dimensionless)"; break;
  case ScaleUnit::kMetre: out << "1 (Metre)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const ColourMatrix& val)
{
  switch (val) {
  case ColourMatrix::kIdentity: out << "0 (Identity)"; break;
  case ColourMatrix::kBt709: out << "1 (Bt709)"; break;
  case ColourMatrix::kUnspecified: out << "2 (Unspecified)"; break;
  case ColourMatrix::kReserved_3: out << "3 (Reserved)"; break;
  case ColourMatrix::kUsa47Cfr73dot682a20:
    out << "4 (Usa47Cfr73dot682a20)";
    break;
  case ColourMatrix::kBt601: out << "5 (Bt601)"; break;
  case ColourMatrix::kSmpte170M: out << "6 (Smpte170M)"; break;
  case ColourMatrix::kSmpte240M: out << "7 (Smpte240M)"; break;
  case ColourMatrix::kYCgCo: out << "8 (kYCgCo)"; break;
  case ColourMatrix::kBt2020Ncl: out << "9 (Bt2020Ncl)"; break;
  case ColourMatrix::kBt2020Cl: out << "10 (Bt2020Cl)"; break;
  case ColourMatrix::kSmpte2085: out << "11 (Smpte2085)"; break;
  default: out << "Unknown"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const AxisOrder& val)
{
  switch (val) {
  case AxisOrder::kZYX: out << "0 (zyx)"; break;
  case AxisOrder::kXYZ: out << "1 (xyz)"; break;
  case AxisOrder::kXZY: out << "2 (xzy)"; break;
  case AxisOrder::kYZX: out << "3 (yzx)"; break;
  case AxisOrder::kZYX_4: out << "4 (zyx)"; break;
  case AxisOrder::kZXY: out << "5 (zxy)"; break;
  case AxisOrder::kYXZ: out << "6 (yxz)"; break;
  case AxisOrder::kXYZ_7: out << "7 (xyz)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const AttributeEncoding& val)
{
  switch (val) {
  case AttributeEncoding::kRAHTransform: out << "0 (RAHT)"; break;
  case AttributeEncoding::kPredictingTransform: out << "1 (Pred)"; break;
  case AttributeEncoding::kLiftingTransform: out << "2 (Lift)"; break;
  case AttributeEncoding::kRaw: out << "3 (Raw)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const LodDecimationMethod& val)
{
  switch (val) {
  case LodDecimationMethod::kNone: out << "0 (None)"; break;
  case LodDecimationMethod::kPeriodic: out << "1 (Periodic)"; break;
  case LodDecimationMethod::kCentroid: out << "2 (Centroid)"; break;
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
  case PartitionMethod::kUniformGeom: out << "2 (UniformGeom)"; break;
  case PartitionMethod::kOctreeUniform: out << "3 (UniformOctree)"; break;
  case PartitionMethod::kUniformSquare: out << "4 (UniformSquare)"; break;
  case PartitionMethod::kNpoints: out << "5 (NPointSpans)"; break;
  default: out << int(val) << " (Unknown)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const PredGeomEncOpts::SortMode& val)
{
  switch (val) {
    using SortMode = PredGeomEncOpts::SortMode;
  case SortMode::kNoSort: out << int(val) << " (None)"; break;
  case SortMode::kSortMorton: out << int(val) << " (Morton)"; break;
  case SortMode::kSortAzimuth: out << int(val) << " (Azimuth)"; break;
  case SortMode::kSortRadius: out << int(val) << " (Radius)"; break;
  case SortMode::kSortLaserAngle: out << int(val) << " (LaserAngle)"; break;
  default: out << int(val) << " (Unknown)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const OctreeEncOpts::QpMethod& val)
{
  switch (val) {
    using Method = OctreeEncOpts::QpMethod;
  case Method::kUniform: out << int(val) << " (Uniform)"; break;
  case Method::kRandom: out << int(val) << " (Random)"; break;
  case Method::kByDensity: out << int(val) << " (ByDensity)"; break;
  default: out << int(val) << " (Unknown)"; break;
  }
  return out;
}
}  // namespace pcc

namespace df {
namespace program_options_lite {
  template<typename T>
  struct option_detail<pcc::Vec3<T>> {
    static constexpr bool is_container = true;
    static constexpr bool is_fixed_size = true;
    typedef T* output_iterator;

    static void clear(pcc::Vec3<T>& container){};
    static output_iterator make_output_iterator(pcc::Vec3<T>& container)
    {
      return &container[0];
    }
  };
}  // namespace program_options_lite
}  // namespace df

//---------------------------------------------------------------------------
// :: Command line / config parsing

void sanitizeEncoderOpts(
  Parameters& params, df::program_options_lite::ErrorReporter& err);

//---------------------------------------------------------------------------

bool
ParseParameters(int argc, char* argv[], Parameters& params)
{
  namespace po = df::program_options_lite;

  struct {
    AttributeDescription desc;
    AttributeParameterSet aps;
    EncoderAttributeParams encoder;
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
        params.encoder.attr.push_back(params_attr.encoder);
        return;
      }

      // update existing entry
      params.encoder.sps.attributeSets[it.first->second] = params_attr.desc;
      params.encoder.aps[it.first->second] = params_attr.aps;
      params.encoder.attr[it.first->second] = params_attr.encoder;
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
  ("firstFrameNum",
     params.firstFrameNum, 0,
     "Frame number for use with interpolating %d format specifiers "
     "in input/output filenames")

  ("frameCount",
     params.frameCount, 1,
     "Number of frames to encode")

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
    params.postRecolorPath, {},
    "Recolored pointcloud file path (encoder only)")

  ("preInvScalePath",
    params.preInvScalePath, {},
    "Pre inverse scaled pointcloud file path (decoder only)")

  ("convertPlyColourspace",
    params.convertColourspace, true,
    "Convert ply colourspace according to attribute colourMatrix")

  ("outputBinaryPly",
    params.outputBinaryPly, true,
    "Output ply files using binary (or ascii) format")

  ("outputUnitLength",
    params.outputUnitLength, 0.,
    "Length of reconstructed point cloud x,y,z unit vectors\n"
    " 0: use srcUnitLength")

  ("outputScaling",
    params.outputSystem, OutputSystem::kExternal,
    "Output coordnate system scaling\n"
    " 0: Conformance\n"
    " 1: External")

  ("outputPrecisionBits",
    params.outputFpBits, -1,
    "Fractional bits in conformance output (prior to external scaling)\n"
    " 0: integer,  -1: automatic (full)")

  // This section controls all general geometry scaling parameters
  (po::Section("Coordinate system scaling"))

  ("srcUnitLength",
    params.encoder.srcUnitLength, 1.,
    "Length of source point cloud x,y,z unit vectors in srcUnits")

  ("srcUnit",
    params.encoder.sps.seq_geom_scale_unit_flag, ScaleUnit::kDimensionless,
    " 0: dimensionless\n 1: metres")

  ("inputScale",
    params.inputScale, 1.,
    "Scale input while reading src ply. "
    "Eg, 1000 converts metres to integer millimetres")

  ("codingScale",
    params.encoder.codedGeomScale, 1.,
    "Scale used to represent coded geometry. Relative to inputScale")

  ("sequenceScale",
    params.encoder.seqGeomScale, 1.,
    "Scale used to obtain sequence coordinate system. "
    "Relative to inputScale")

  // Alias for compatibility with old name.
  ("positionQuantizationScale", params.encoder.seqGeomScale, 1.,
   "(deprecated)")

  ("externalScale",
    params.encoder.extGeomScale, 1.,
    "Scale used to define external coordinate system.\n"
    "Meaningless when srcUnit = metres\n"
    "  0: Use srcUnitLength\n"
    " >0: Relative to inputScale")

  (po::Section("Decoder"))

  ("skipOctreeLayers",
    params.decoder.minGeomNodeSizeLog2, 0,
    "Partial decoding of octree and attributes\n"
    " 0   : Full decode\n"
    " N>0 : Skip the bottom N layers in decoding process")

  ("decodeMaxPoints",
    params.decoder.decodeMaxPoints, 0,
    "Partially decode up to N points")

  (po::Section("Encoder"))

  ("sortInputByAzimuth",
    params.sortInputByAzimuth, false,
    "Sort input points by azimuth angle")

  ("geometry_axis_order",
    params.encoder.sps.geometry_axis_order, AxisOrder::kXYZ,
    "Sets the geometry axis coding order:\n"
    "  0: (zyx)\n  1: (xyz)\n  2: (xzy)\n"
    "  3: (yzx)\n  4: (zyx)\n  5: (zxy)\n"
    "  6: (yxz)\n  7: (xyz)")

  ("autoSeqBbox",
    params.encoder.autoSeqBbox, true,
    "Calculate seqOrigin and seqSizeWhd automatically.")

  // NB: the underlying variable is in STV order.
  //     Conversion happens during argument sanitization.
  ("seqOrigin",
    params.encoder.sps.seqBoundingBoxOrigin, {0},
    "Origin (x,y,z) of the sequence bounding box "
    "(in input coordinate system). "
    "Requires autoSeqBbox=0")

  // NB: the underlying variable is in STV order.
  //     Conversion happens during argument sanitization.
  ("seqSizeWhd",
    params.encoder.sps.seqBoundingBoxSize, {0},
    "Size of the sequence bounding box "
    "(in input coordinate system). "
    "Requires autoSeqBbox=0")

  ("mergeDuplicatedPoints",
    params.encoder.gps.geom_unique_points_flag, true,
    "Enables removal of duplicated points")

  ("partitionMethod",
    params.encoder.partition.method, PartitionMethod::kUniformSquare,
    "Method used to partition input point cloud into slices/tiles:\n"
    "  0: none\n"
    "  2: n Uniform-geometry partition bins along the longest edge\n"
    "  3: Uniform geometry partition at n octree depth\n"
    "  4: Uniform square partition\n"
    "  5: n-point spans of input")

  ("partitionOctreeDepth",
    params.encoder.partition.octreeDepth, 1,
    "Depth of octree partition for partitionMethod=4")

  ("sliceMaxPoints",
    params.encoder.partition.sliceMaxPoints, 1100000,
    "Maximum number of points per slice")

  ("sliceMinPoints",
    params.encoder.partition.sliceMinPoints, 550000,
    "Minimum number of points per slice (soft limit)")

  ("tileSize",
    params.encoder.partition.tileSize, 0,
    "Partition input into cubic tiles of given size")

  ("cabac_bypass_stream_enabled_flag",
    params.encoder.sps.cabac_bypass_stream_enabled_flag, false,
    "Controls coding method for ep(bypass) bins")

  ("entropyContinuationEnabled",
    params.encoder.sps.entropy_continuation_enabled_flag, false,
    "Propagate context state between slices")

  ("disableAttributeCoding",
    params.disableAttributeCoding, false,
    "Ignore attribute coding configuration")

  ("enforceLevelLimits",
    params.encoder.enforceLevelLimits, true,
    "Abort if level limits exceeded")

  (po::Section("Geometry"))

  ("geomTreeType",
    params.encoder.gps.predgeom_enabled_flag, false,
    "Selects the tree coding method:\n"
    "  0: octree\n"
    "  1: predictive")

  ("qtbtEnabled",
    params.encoder.gps.qtbt_enabled_flag, true,
    "Enables non-cubic geometry bounding box")

  ("maxNumQtBtBeforeOt",
    params.encoder.geom.qtbt.maxNumQtBtBeforeOt, 4,
    "Max number of qtbt partitions before ot")

  ("minQtbtSizeLog2",
    params.encoder.geom.qtbt.minQtbtSizeLog2, 0,
    "Minimum size of qtbt partitions")

  ("numOctreeEntropyStreams",
    // NB: this is adjusted by minus 1 after the arguments are parsed
    params.encoder.gbh.geom_stream_cnt_minus1, 1,
    "Number of entropy streams for octree coding")

  ("bitwiseOccupancyCoding",
    params.encoder.gps.bitwise_occupancy_coding_flag, true,
    "Selects between bitwise and bytewise occupancy coding:\n"
    "  0: bytewise\n"
    "  1: bitwise")

  ("neighbourAvailBoundaryLog2",
    // NB: this is adjusted by minus 1 after the arguments are parsed
    params.encoder.gps.neighbour_avail_boundary_log2_minus1, 0,
    "Defines the avaliability volume for neighbour occupancy lookups:\n"
    "<2: Limited to sibling nodes only")

  ("inferredDirectCodingMode",
    params.encoder.gps.inferred_direct_coding_mode, 1,
    "Early termination of the geometry octree for isolated points:"
    " 0: disabled\n"
    " 1: fully constrained\n"
    " 2: partially constrained\n"
    " 3: unconstrained (fastest)")

  ("jointTwoPointIdcm",
    params.encoder.gps.joint_2pt_idcm_enabled_flag, true,
    "Jointly code common prefix of two IDCM points")

  ("adjacentChildContextualization",
    params.encoder.gps.adjacent_child_contextualization_enabled_flag, true,
    "Occupancy contextualization using neighbouring adjacent children")

  ("intra_pred_max_node_size_log2",
    params.encoder.gps.intra_pred_max_node_size_log2, 0,
    "octree nodesizes eligible for occupancy intra prediction")

  ("planarEnabled",
    params.encoder.gps.geom_planar_mode_enabled_flag, true,
    "Use planar mode for geometry coding")

  ("planarModeThreshold0",
    params.encoder.gps.geom_planar_threshold0, 77,
    "Activation threshold (0-127) of first planar mode. "
    "Lower values imply more use of the first planar mode")

  ("planarModeThreshold1",
    params.encoder.gps.geom_planar_threshold1, 99,
    "Activation threshold (0-127) of second planar mode. "
    "Lower values imply more use of the first planar mode")

  ("planarModeThreshold2",
    params.encoder.gps.geom_planar_threshold2, 113,
    "Activation threshold (0-127) of third planar mode. "
    "Lower values imply more use of the third planar mode")

   ("planarModeIdcmUse",
    // NB: this is adjusted by minus1 after thearguments are parsed
    params.encoder.gps.geom_idcm_rate_minus1, 0,
    "Degree (1/32%) of IDCM activation when planar mode is enabled\n"
    "  0 => never, 32 => always")

  ("trisoupNodeSizeLog2",
    params.encoder.trisoupNodeSizesLog2, {0},
    "Node size for surface triangulation\n"
    " <2: disabled")

  ("trisoup_sampling_value",
    params.encoder.gps.trisoup_sampling_value, 0,
    "Trisoup voxelisation sampling rate\n"
    "  0: automatic")

  ("positionQuantisationEnabled",
    params.encoder.gps.geom_scaling_enabled_flag, false,
    "Enable in-loop quantisation of positions")

  ("positionQuantisationMethod",
    params.encoder.geom.qpMethod, OctreeEncOpts::QpMethod::kUniform,
    "Method used to determine per-node QP:\n"
    "  0: uniform\n"
    "  1: random\n"
    "  2: by node point density")

  ("positionQpMultiplierLog2",
    params.encoder.gps.geom_qp_multiplier_log2, 0,
    "Granularity of QP to step size mapping:\n"
    "  n: 2^n QPs per doubling interval, n in 0..3")

  ("positionBaseQp",
    params.encoder.gps.geom_base_qp, 0,
    "Base QP used in position quantisation (0 = lossless)")

  ("positionIdcmQp",
    params.encoder.idcmQp, 0,
    "QP used in position quantisation of IDCM nodes")

  ("positionSliceQpOffset",
    params.encoder.gbh.geom_slice_qp_offset, 0,
    "Per-slice QP offset used in position quantisation")

  ("positionQuantisationOctreeSizeLog2",
    params.encoder.geom.qpOffsetNodeSizeLog2, -1,
    "Octree node size used for signalling position QP offsets "
    "(-1 => disabled)")

  ("positionQuantisationOctreeDepth",
    params.encoder.geom.qpOffsetDepth, -1,
    "Octree depth used for signalling position QP offsets (-1 => disabled)")

  ("positionBaseQpFreqLog2",
    params.encoder.gps.geom_qp_offset_intvl_log2, 8,
    "Frequency of sending QP offsets in predictive geometry coding")

  // NB: this will be corrected to be relative to base value later
  ("positionSliceQpFreqLog2",
    params.encoder.gbh.geom_qp_offset_intvl_log2_delta, 0,
    "Frequency of sending QP offsets in predictive geometry coding")

  ("angularEnabled",
    params.encoder.gps.geom_angular_mode_enabled_flag, false,
    "Controls angular contextualisation of occupancy")

  // NB: the underlying variable is in STV order.
  //     Conversion happens during argument sanitization.
  ("lidarHeadPosition",
    params.encoder.gps.gpsAngularOrigin, {0, 0, 0},
    "laser head position (x,y,z) in angular mode")

  ("numLasers",
    params.encoder.numLasers, 0,
    "Number of lasers in angular mode")

  ("lasersTheta",
    params.encoder.lasersTheta, {},
    "Vertical laser angle in angular mode")

  ("lasersZ",
    params.encoder.lasersZ, {},
    "Vertical laser offset in angular mode")

  ("lasersNumPhiPerTurn",
    params.encoder.gps.angularNumPhiPerTurn, {},
    "Number of sampling poisitions in a complete laser turn in angular mode")

  ("planarBufferDisabled",
    params.encoder.gps.planar_buffer_disabled_flag, false,
    "Disable planar buffer (when angular mode is enabled)")

  ("predGeomAzimuthQuantization",
    params.encoder.gps.azimuth_scaling_enabled_flag, true,
    "Quantize azimuth according to radius in predictive geometry coding")

  ("positionAzimuthScaleLog2",
    params.encoder.gps.geom_angular_azimuth_scale_log2_minus11, 5,
    "Additional bits to represent azimuth angle in predictive geometry coding")

  // NB: this will be corrected to be minus 1 later
  ("positionAzimuthSpeed",
    params.encoder.gps.geom_angular_azimuth_speed_minus1, 363,
    "Scale factor applied to azimuth angle in predictive geometry coding")

  ("positionRadiusInvScaleLog2",
    params.encoder.gps.geom_angular_radius_inv_scale_log2, 0,
    "Inverse scale factor applied to radius in predictive geometry coding")

  ("predGeomSort",
    params.encoder.predGeom.sortMode, PredGeomEncOpts::kSortMorton,
    "Predictive geometry tree construction order")

  ("predGeomAzimuthSortPrecision",
    params.encoder.predGeom.azimuthSortRecipBinWidth, 0.0f,
    "Reciprocal precision used in azimuthal sorting for tree construction")

  ("predGeomTreePtsMax",
    params.encoder.predGeom.maxPtsPerTree, 1100000,
    "Maximum number of points per predictive geometry tree")

  ("pointCountMetadata",
    params.encoder.gps.octree_point_count_list_present_flag, false,
    "Add octree layer point count metadata")

  (po::Section("Attributes"))

  // attribute processing
  //   NB: Attribute options are special in the way they are applied (see above)
  ("attribute",
    attribute_setter,
    "Encode the given attribute (NB, must appear after the"
    "following attribute parameters)")

  // NB: the cli option sets +1, the minus1 will be applied later
  ("attrScale",
    params_attr.desc.params.attr_scale_minus1, 1,
    "Scale factor used to interpret coded attribute values")

  ("attrOffset",
    params_attr.desc.params.attr_offset, 0,
    "Offset used to interpret coded attribute values")

  ("bitdepth",
    params_attr.desc.bitdepth, 8,
    "Attribute bitdepth")

  ("defaultValue",
    params_attr.desc.params.attr_default_value, {},
    "Default attribute component value(s) in case of data omission")

  // todo(df): this should be per-attribute
  ("colourMatrix",
    params_attr.desc.params.cicp_matrix_coefficients_idx, ColourMatrix::kBt709,
    "Matrix used in colourspace conversion\n"
    "  0: none (identity)\n"
    "  1: ITU-T BT.709\n"
    "  8: YCgCo")

  ("transformType",
    params_attr.aps.attr_encoding, AttributeEncoding::kPredictingTransform,
    "Coding method to use for attribute:\n"
    "  0: Region Adaptive Hierarchical Transform (RAHT)\n"
    "  1: Hierarchical neighbourhood prediction\n"
    "  2: Hierarchical neighbourhood prediction as lifting transform")

  ("rahtPredictionEnabled",
    params_attr.aps.raht_prediction_enabled_flag, true,
    "Controls the use of transform-domain prediction")

  ("rahtPredictionThreshold0",
    params_attr.aps.raht_prediction_threshold0, 2,
    "Grandparent threshold for early transform-domain prediction termination")

  ("rahtPredictionThreshold1",
    params_attr.aps.raht_prediction_threshold1, 6,
    "Parent threshold for early transform-domain prediction termination")

  // NB: the cli option sets +1, the minus1 will be applied later
  ("numberOfNearestNeighborsInPrediction",
    params_attr.aps.num_pred_nearest_neighbours_minus1, 3,
    "Attribute's maximum number of nearest neighbors to be used for prediction")

  ("adaptivePredictionThreshold",
    params_attr.aps.adaptive_prediction_threshold, 1 << 6,
    "Neighbouring attribute value difference that enables direct "
    "prediction. 8-bit value scaled to attribute bitdeph. "
    "Applies to transformType=0 only")

  ("intraLodSearchRange",
    params_attr.aps.intra_lod_search_range, -1,
    "Intra LoD nearest neighbor search range\n"
    " -1: Full-range")

  ("interLodSearchRange",
    params_attr.aps.inter_lod_search_range, -1,
    "Inter LoD nearest neighbor search range\n"
    " -1: Full-range")

  // NB: the underlying variable is in STV order.
  //     Conversion happens during argument sanitization.
  ("lod_neigh_bias",
    params_attr.aps.lodNeighBias, {1, 1, 1},
    "Attribute's (x,y,z) component intra prediction weights")

  ("lodDecimator",
    params_attr.aps.lod_decimation_type, LodDecimationMethod::kNone,
    "LoD decimation method:\n"
    " 0: none\n"
    " 1: periodic subsampling using lodSamplingPeriod\n"
    " 2: centroid subsampling using lodSamplingPeriod")

  ("max_num_direct_predictors",
    params_attr.aps.max_num_direct_predictors, 3,
    "Maximum number of nearest neighbour candidates used in direct"
    "attribute prediction")

  ("direct_avg_predictor_disabled_flag",
    params_attr.aps.direct_avg_predictor_disabled_flag, false,
    "Disable average predictor")

  ("predWeightBlending",
    params_attr.aps.pred_weight_blending_enabled_flag, false,
    "Blend prediction weights according to neigbour distances. "
    "Applies to transformType=0 only")

  // NB: this parameter actually represents the number of refinement layers
  ("levelOfDetailCount",
    params_attr.aps.num_detail_levels_minus1, 1,
    "Attribute's number of levels of detail")

  ("dist2",
    params_attr.aps.dist2, 0,
    "Initial squared distance used in LoD generation")

  ("dist2PercentileEstimate",
    params_attr.encoder.dist2PercentileEstimate, 0.85f,
    "Percentile for dist2 estimation during nearest neighbour search")

  ("positionQuantizationScaleAdjustsDist2",
    params.positionQuantizationScaleAdjustsDist2, false,
    "Scale dist2 values by squared positionQuantizationScale")

  ("lodSamplingPeriod",
    params_attr.aps.lodSamplingPeriod, {4},
    "List of per LoD sampling periods used in LoD generation")

  ("intraLodPredictionSkipLayers",
    params_attr.aps.intra_lod_prediction_skip_layers, -1,
    "Number of finest detail levels that skip intra prediction\n"
    " -1: skip all (disables intra pred)")

  ("interComponentPredictionEnabled",
    params_attr.aps.inter_component_prediction_enabled_flag, false,
    "Use primary attribute component to predict values of subsequent "
    "components")

  ("lastComponentPredictionEnabled",
    params_attr.aps.last_component_prediction_enabled_flag, true,
    "Use second attribute component to predict value of the final component")

  ("canonical_point_order_flag",
    params_attr.aps.canonical_point_order_flag, false,
    "Enable skipping morton sort in case of number of LoD equal to 1")

  ("spherical_coord_flag",
     params_attr.aps.spherical_coord_flag, false,
     "Code attributes in spherical domain")

  ("attrSphericalMaxLog2",
    params.encoder.attrSphericalMaxLog2, 0,
    "Override spherical coordinate normalisation factor")

  ("aps_scalable_enable_flag",
    params_attr.aps.scalable_lifting_enabled_flag, false,
    "Enable scalable attritube coding")

  ("max_neigh_range",
    // NB: this is adjusted by minus 1 after the arguments are parsed
    params_attr.aps.max_neigh_range_minus1, 5,
    "maximum nearest neighbour range for scalable lifting")

  ("qp",
    // NB: this is adjusted with minus 4 after the arguments are parsed
    params_attr.aps.init_qp_minus4, 4,
    "Attribute's luma quantisation parameter")

  ("qpChromaOffset",
    params_attr.aps.aps_chroma_qp_offset, 0,
    "Attribute's chroma quantisation parameter offset (relative to luma)")

  ("aps_slice_qp_deltas_present_flag",
    params_attr.aps.aps_slice_qp_deltas_present_flag, false,
    "Enable signalling of per-slice QP values")

  ("qpLayerOffsetsLuma",
    params_attr.encoder.abh.attr_layer_qp_delta_luma, {},
      "Attribute's per layer luma QP offsets")

  ("qpLayerOffsetsChroma",
      params_attr.encoder.abh.attr_layer_qp_delta_chroma, {},
      "Attribute's per layer chroma QP offsets")

  ("quantNeighWeight",
    params_attr.aps.quant_neigh_weight, {16, 8, 4},
    "Factors used to derive quantization weights (transformType=1)")

  // This section is just dedicated to attribute recolouring (encoder only).
  // parameters are common to all attributes.
  (po::Section("Recolouring"))

  ("recolourSearchRange",
    params.encoder.recolour.searchRange, 1,
    "")

  ("recolourNumNeighboursFwd",
    params.encoder.recolour.numNeighboursFwd, 8,
    "")

  ("recolourNumNeighboursBwd",
    params.encoder.recolour.numNeighboursBwd, 1,
    "")

  ("recolourUseDistWeightedAvgFwd",
    params.encoder.recolour.useDistWeightedAvgFwd, true,
    "")

  ("recolourUseDistWeightedAvgBwd",
    params.encoder.recolour.useDistWeightedAvgBwd, true,
    "")

  ("recolourSkipAvgIfIdenticalSourcePointPresentFwd",
    params.encoder.recolour.skipAvgIfIdenticalSourcePointPresentFwd, true,
    "")

  ("recolourSkipAvgIfIdenticalSourcePointPresentBwd",
    params.encoder.recolour.skipAvgIfIdenticalSourcePointPresentBwd, false,
    "")

  ("recolourDistOffsetFwd",
    params.encoder.recolour.distOffsetFwd, 4.,
    "")

  ("recolourDistOffsetBwd",
    params.encoder.recolour.distOffsetBwd, 4.,
    "")

  ("recolourMaxGeometryDist2Fwd",
    params.encoder.recolour.maxGeometryDist2Fwd, 1000.,
    "")

  ("recolourMaxGeometryDist2Bwd",
    params.encoder.recolour.maxGeometryDist2Bwd, 1000.,
    "")

  ("recolourMaxAttributeDist2Fwd",
    params.encoder.recolour.maxAttributeDist2Fwd, 1000.,
    "")

  ("recolourMaxAttributeDist2Bwd",
    params.encoder.recolour.maxAttributeDist2Bwd, 1000.,
    "")

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

  // set default output units (this works for the decoder too)
  if (params.outputUnitLength <= 0.)
    params.outputUnitLength = params.encoder.srcUnitLength;
  params.encoder.outputFpBits = params.outputFpBits;
  params.decoder.outputFpBits = params.outputFpBits;

  if (!params.isDecoder)
    sanitizeEncoderOpts(params, err);

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
    po::dumpCfg(cout, opts, "Coordinate system scaling", 4);
    po::dumpCfg(cout, opts, "Encoder", 4);
    po::dumpCfg(cout, opts, "Geometry", 4);
    po::dumpCfg(cout, opts, "Recolouring", 4);

    for (const auto& it : params.encoder.attributeIdxMap) {
      // NB: when dumping the config, opts references params_attr
      params_attr.desc = params.encoder.sps.attributeSets[it.second];
      params_attr.aps = params.encoder.aps[it.second];
      params_attr.encoder = params.encoder.attr[it.second];
      cout << "    " << it.first << "\n";
      po::dumpCfg(cout, opts, "Attributes", 8);
    }
  }

  cout << endl;

  return true;
}

//----------------------------------------------------------------------------

void
sanitizeEncoderOpts(
  Parameters& params, df::program_options_lite::ErrorReporter& err)
{
  // Input scaling affects the definition of the source unit length.
  // eg, if the unit length of the source is 1m, scaling by 1000 generates
  // a cloud with unit length 1mm.
  params.encoder.srcUnitLength /= params.inputScale;

  // global scale factor must be positive
  if (params.encoder.codedGeomScale > params.encoder.seqGeomScale) {
    err.warn() << "codingScale must be <= sequenceScale, adjusting\n";
    params.encoder.codedGeomScale = params.encoder.seqGeomScale;
  }

  // fix the representation of various options
  params.encoder.gbh.geom_stream_cnt_minus1--;
  params.encoder.gps.geom_idcm_rate_minus1--;
  params.encoder.gps.geom_angular_azimuth_speed_minus1--;
  params.encoder.gps.neighbour_avail_boundary_log2_minus1 =
    std::max(0, params.encoder.gps.neighbour_avail_boundary_log2_minus1 - 1);
  for (auto& attr_sps : params.encoder.sps.attributeSets) {
    attr_sps.params.attr_scale_minus1--;
  }
  for (auto& attr_aps : params.encoder.aps) {
    attr_aps.init_qp_minus4 -= 4;
    attr_aps.num_pred_nearest_neighbours_minus1--;
    attr_aps.max_neigh_range_minus1--;
  }

  // Config options are absolute, but signalling is relative
  params.encoder.gbh.geom_qp_offset_intvl_log2_delta -=
    params.encoder.gps.geom_qp_offset_intvl_log2;

  // If idcm rate is configured as 0, disable idcm
  // NB: if user has requested less contrained idcm, warn
  if (params.encoder.gps.geom_idcm_rate_minus1 < 0) {
    if (params.encoder.gps.inferred_direct_coding_mode == 1)
      params.encoder.gps.inferred_direct_coding_mode = 0;
  }

  if (params.encoder.gps.geom_idcm_rate_minus1 < 31) {
    if (params.encoder.gps.inferred_direct_coding_mode > 1) {
      params.encoder.gps.geom_idcm_rate_minus1 = 31;
      err.warn() << "ignoring planarModeIdcmUse < 32: "
                    "contradicts inferredDirectCodingMode > 1\n";
    }
  }

  // convert coordinate systems if the coding order is different from xyz
  convertXyzToStv(&params.encoder.sps);
  convertXyzToStv(params.encoder.sps, &params.encoder.gps);
  for (auto& aps : params.encoder.aps)
    convertXyzToStv(params.encoder.sps, &aps);

  // Trisoup is enabled when a node size is specified
  // sanity: don't enable if only node size is 0.
  // todo(df): this needs to take into account slices where it is disabled
  if (params.encoder.trisoupNodeSizesLog2.size() == 1)
    if (params.encoder.trisoupNodeSizesLog2[0] < 2)
      params.encoder.trisoupNodeSizesLog2.clear();

  for (auto trisoupNodeSizeLog2 : params.encoder.trisoupNodeSizesLog2)
    if (trisoupNodeSizeLog2 < 2)
      err.error() << "Trisoup node size must be greater than 1\n";

  params.encoder.gps.trisoup_enabled_flag =
    !params.encoder.trisoupNodeSizesLog2.empty();

  // Certain coding modes are not available when trisoup is enabled.
  // Disable them, and warn if set (they may be set as defaults).
  if (params.encoder.gps.trisoup_enabled_flag) {
    if (!params.encoder.gps.geom_unique_points_flag)
      err.warn() << "TriSoup geometry does not preserve duplicated points\n";

    if (params.encoder.gps.inferred_direct_coding_mode)
      err.warn() << "TriSoup geometry is incompatable with IDCM\n";

    params.encoder.gps.geom_unique_points_flag = true;
    params.encoder.gps.inferred_direct_coding_mode = 0;
  }

  // tweak qtbt generation when trisoup is /isn't enabled
  params.encoder.geom.qtbt.trisoupEnabled =
    params.encoder.gps.trisoup_enabled_flag;

  // Planar coding mode is not available for bytewise coding
  if (!params.encoder.gps.bitwise_occupancy_coding_flag) {
    if (params.encoder.gps.geom_planar_mode_enabled_flag)
      err.warn() << "Bytewise geometry coding does not support planar mode\n";
    params.encoder.gps.geom_planar_mode_enabled_flag = false;
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
    auto& attr_enc = params.encoder.attr[it.second];

    // default values for attribute
    attr_sps.attr_instance_id = 0;
    auto& attrMeta = attr_sps.params;
    attrMeta.cicp_colour_primaries_idx = 2;
    attrMeta.cicp_transfer_characteristics_idx = 2;
    attrMeta.cicp_video_full_range_flag = true;
    attrMeta.cicpParametersPresent = false;
    attrMeta.attr_frac_bits = 0;
    attrMeta.scalingParametersPresent = false;

    // Enable scaling if a paramter has been set
    //  - pre/post scaling is only currently supported for reflectance
    attrMeta.scalingParametersPresent = attrMeta.attr_offset
      || attrMeta.attr_scale_minus1 || attrMeta.attr_frac_bits;

    // todo(df): remove this hack when scaling is generalised
    if (it.first != "reflectance" && attrMeta.scalingParametersPresent) {
      err.warn() << it.first << ": scaling not supported, disabling\n";
      attrMeta.scalingParametersPresent = 0;
    }

    if (it.first == "reflectance") {
      // Avoid wasting bits signalling chroma quant step size for reflectance
      attr_aps.aps_chroma_qp_offset = 0;
      attr_enc.abh.attr_layer_qp_delta_chroma.clear();

      // There is no matrix for reflectace
      attrMeta.cicp_matrix_coefficients_idx = ColourMatrix::kUnspecified;
      attr_sps.attr_num_dimensions_minus1 = 0;
      attr_sps.attributeLabel = KnownAttributeLabel::kReflectance;
    }

    if (it.first == "color") {
      attr_sps.attr_num_dimensions_minus1 = 2;
      attr_sps.attributeLabel = KnownAttributeLabel::kColour;
      attrMeta.cicpParametersPresent = true;
    }

    // Assume that YCgCo is actually YCgCoR for now
    // This requires an extra bit to represent chroma (luma will have a
    // reduced range)
    if (attrMeta.cicp_matrix_coefficients_idx == ColourMatrix::kYCgCo)
      attr_sps.bitdepth++;

    // Extend the default attribute value to the correct width if present
    if (!attrMeta.attr_default_value.empty())
      attrMeta.attr_default_value.resize(
        attr_sps.attr_num_dimensions_minus1 + 1,
        attrMeta.attr_default_value.back());

    // In order to simplify specification of dist2 values, which are
    // depending on the scale of the coded point cloud, the following
    // adjust the dist2 values according to PQS.  The user need only
    // specify the unquantised PQS value.
    if (params.positionQuantizationScaleAdjustsDist2) {
      auto delta = log2(params.encoder.codedGeomScale);
      attr_aps.dist2 =
        std::max(0, int32_t(std::round(attr_aps.dist2 + delta)));
    }

    // derive samplingPeriod values based on initial value
    if (
      !attr_aps.lodParametersPresent()
      || (attr_aps.lod_decimation_type == LodDecimationMethod::kNone)) {
      attr_aps.lodSamplingPeriod.clear();
    } else if (!attr_aps.lodSamplingPeriod.empty()) {
      auto i = attr_aps.lodSamplingPeriod.size();
      attr_aps.lodSamplingPeriod.resize(attr_aps.num_detail_levels_minus1);
      // add any extra values as required
      for (; i < attr_aps.num_detail_levels_minus1; i++)
        attr_aps.lodSamplingPeriod[i] = attr_aps.lodSamplingPeriod[i - 1];
    }

    if (attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform) {
      attr_aps.adaptive_prediction_threshold = 0;
      attr_aps.intra_lod_prediction_skip_layers = -1;
    }

    // For RAHT, ensure that the unused lod count = 0 (prevents mishaps)
    if (attr_aps.attr_encoding == AttributeEncoding::kRAHTransform) {
      attr_aps.num_detail_levels_minus1 = 0;
      attr_aps.adaptive_prediction_threshold = 0;
    }

    if (!params.encoder.gps.geom_angular_mode_enabled_flag) {
      if (attr_aps.spherical_coord_flag)
        err.warn() << it.first
                   << ".spherical_coord_flag=1 requires angularEnabled=1, "
                      "disabling\n";
      attr_aps.spherical_coord_flag = false;
    }
  }

  // convert floating point values of Lasers' Theta and H to fixed point
  if (params.encoder.gps.geom_angular_mode_enabled_flag) {
    if (params.encoder.numLasers == 0)
      err.error() << "numLasers must be at least 1\n";

    for (auto val : params.encoder.lasersTheta) {
      int one = 1 << 18;
      params.encoder.gps.angularTheta.push_back(round(val * one));
    }

    for (auto val : params.encoder.lasersZ) {
      int one = 1 << 3;
      auto scale = params.encoder.codedGeomScale;
      if (params.encoder.gps.predgeom_enabled_flag)
        scale = params.encoder.codedGeomScale / params.encoder.seqGeomScale;

      params.encoder.gps.angularZ.push_back(round(val * scale * one));
    }

    if (params.encoder.gps.angularTheta.size() != params.encoder.numLasers)
      err.error() << "lasersZ.size() != numLasers\n";

    if (params.encoder.gps.angularZ.size() != params.encoder.numLasers)
      err.error() << "lasersTheta.size() != numLasers\n";

    if (
      params.encoder.gps.angularNumPhiPerTurn.size()
      != params.encoder.numLasers)
      err.error() << "lasersNumPhiPerTurn.size() != numLasers\n";

    if (params.encoder.gps.qtbt_enabled_flag) {
      params.encoder.geom.qtbt.angularMaxNodeMinDimLog2ToSplitV =
        std::max<int>(0, 8 + log2(params.encoder.codedGeomScale));
      params.encoder.geom.qtbt.angularMaxDiffToSplitZ =
        std::max<int>(0, 1 + log2(params.encoder.codedGeomScale));
    }

    if (params.encoder.gps.predgeom_enabled_flag) {
      auto& gps = params.encoder.gps;
      int maxSpeed = 1 << (gps.geom_angular_azimuth_scale_log2_minus11 + 12);
      if (params.encoder.gps.geom_angular_azimuth_speed_minus1 + 1 > maxSpeed)
        err.error() << "positionAzimuthSpeed > max (" << maxSpeed << ")\n";
    }
  }

  // tweak qtbt when angular is / isn't enabled
  params.encoder.geom.qtbt.angularTweakEnabled =
    params.encoder.gps.geom_angular_mode_enabled_flag;

  if (!params.encoder.geom.qtbt.angularTweakEnabled) {
    // NB: these aren't used in this condition
    params.encoder.geom.qtbt.angularMaxNodeMinDimLog2ToSplitV = 0;
    params.encoder.geom.qtbt.angularMaxDiffToSplitZ = 0;
  }

  // sanity checks

  if (params.encoder.gps.geom_qp_multiplier_log2 & ~3)
    err.error() << "positionQpMultiplierLog2 must be in the range 0..3\n";

  if (!params.encoder.gps.geom_angular_mode_enabled_flag) {
    if (params.encoder.gps.planar_buffer_disabled_flag) {
      params.encoder.gps.planar_buffer_disabled_flag = 0;
      err.warn() << "ignoring planarBufferDisabled without angularEnabled\n";
    }
  }

  // The following featues depend upon the occupancy atlas
  if (!params.encoder.gps.neighbour_avail_boundary_log2_minus1) {
    if (params.encoder.gps.adjacent_child_contextualization_enabled_flag)
      err.warn() << "ignoring adjacentChildContextualization when"
                    " neighbourAvailBoundaryLog2=0\n";

    if (params.encoder.gps.intra_pred_max_node_size_log2)
      err.warn() << "ignoring intra_pred_max_node_size_log2 when"
                    " neighbourAvailBoundaryLog2=0\n";

    params.encoder.gps.adjacent_child_contextualization_enabled_flag = 0;
    params.encoder.gps.intra_pred_max_node_size_log2 = 0;
  }

  if (
    params.encoder.partition.sliceMaxPoints
    < params.encoder.partition.sliceMinPoints)
    err.error()
      << "sliceMaxPoints must be greater than or equal to sliceMinPoints\n";

  for (const auto& it : params.encoder.attributeIdxMap) {
    const auto& attr_sps = params.encoder.sps.attributeSets[it.second];
    const auto& attr_aps = params.encoder.aps[it.second];
    auto& attr_enc = params.encoder.attr[it.second];

    if (it.first == "color") {
      if (
        attr_enc.abh.attr_layer_qp_delta_luma.size()
        != attr_enc.abh.attr_layer_qp_delta_chroma.size()) {
        err.error() << it.first
                    << ".qpLayerOffsetsLuma length != .qpLayerOffsetsChroma\n";
      }
    }

    if (attr_sps.bitdepth > 16)
      err.error() << it.first << ".bitdepth must be less than 17\n";

    if (attr_aps.lodParametersPresent()) {
      int lod = attr_aps.num_detail_levels_minus1;
      if (lod > 255 || lod < 0) {
        err.error() << it.first
                    << ".levelOfDetailCount must be in the range [0,255]\n";
      }

      // if zero, values are derived automatically
      if (attr_aps.dist2 < 0 || attr_aps.dist2 > 20) {
        err.error() << it.first << ".dist2 must be in the range [0,20]\n";
      }

      if (lod > 0 && attr_aps.canonical_point_order_flag) {
        err.error() << it.first
                    << "when levelOfDetailCount > 0, "
                       "canonicalPointOrder must be 0\n";
      }

      if (
        attr_aps.attr_encoding == AttributeEncoding::kPredictingTransform
        && lod == 0 && attr_aps.intra_lod_prediction_skip_layers != 0) {
        err.error()
          << "when transformType == 0 (Pred) and levelOfDetailCount == 0, "
             "intraLodPredictionSkipLayers must be 0\n";
      }

      if (
        (attr_aps.lod_decimation_type != LodDecimationMethod::kNone)
        && attr_aps.lodSamplingPeriod.empty()) {
        err.error() << it.first
                    << ".lodSamplingPeriod must contain at least one entry\n";
      }

      for (auto samplingPeriod : attr_aps.lodSamplingPeriod) {
        if (samplingPeriod < 2)
          err.error() << it.first << ".lodSamplingPeriod values must be > 1\n";
      }

      if (attr_aps.adaptive_prediction_threshold < 0) {
        err.error() << it.first
                    << ".adaptivePredictionThreshold must be positive\n";
      }

      if (
        attr_aps.num_pred_nearest_neighbours_minus1
        >= kAttributePredictionMaxNeighbourCount) {
        err.error() << it.first
                    << ".numberOfNearestNeighborsInPrediction must be <= "
                    << kAttributePredictionMaxNeighbourCount << "\n";
      }

      if (attr_aps.scalable_lifting_enabled_flag) {
        if (attr_aps.lod_decimation_type != LodDecimationMethod::kNone) {
          err.error() << it.first << ".lod_decimation_type must be 0\n";
        }

        if (params.encoder.gps.trisoup_enabled_flag) {
          err.error() << it.first
                      << " trisoup_enabled_flag must be disabled\n";
        }

        if (params.encoder.gps.geom_qp_multiplier_log2 != 3)
          err.error() << it.first << " positionQpMultiplierLog2 must be 3\n";
      }
    }

    if (attr_aps.init_qp_minus4 < 0 || attr_aps.init_qp_minus4 + 4 > 51)
      err.error() << it.first << ".qp must be in the range [4,51]\n";

    if (std::abs(attr_aps.aps_chroma_qp_offset) > 51 - 4) {
      err.error() << it.first
                  << ".qpChromaOffset must be in the range [-47,47]\n";
    }
  }
}

//============================================================================

SequenceEncoder::SequenceEncoder(Parameters* params) : SequenceCodec(params)
{
  // determine the naming (ordering) of ply properties
  _plyAttrNames.position =
    axisOrderToPropertyNames(params->encoder.sps.geometry_axis_order);

  // NB: this is the raw origin before the encoder tweaks it
  _angularOrigin = params->encoder.gps.gpsAngularOrigin;
}

//----------------------------------------------------------------------------

int
SequenceEncoder::compress(Stopwatch* clock)
{
  bytestreamFile.open(params->compressedStreamPath, ios::binary);
  if (!bytestreamFile.is_open()) {
    return -1;
  }

  const int lastFrameNum = params->firstFrameNum + params->frameCount;
  for (frameNum = params->firstFrameNum; frameNum < lastFrameNum; frameNum++) {
    if (compressOneFrame(clock))
      return -1;
  }

  std::cout << "Total bitstream size " << bytestreamFile.tellp() << " B\n";
  bytestreamFile.close();

  return 0;
}

//----------------------------------------------------------------------------

int
SequenceEncoder::compressOneFrame(Stopwatch* clock)
{
  std::string srcName{expandNum(params->uncompressedDataPath, frameNum)};
  PCCPointSet3 pointCloud;
  if (
    !ply::read(srcName, _plyAttrNames, params->inputScale, pointCloud)
    || pointCloud.getPointCount() == 0) {
    cout << "Error: can't open input file!" << endl;
    return -1;
  }

  // Some evaluations wish to scan the points in azimuth order to simulate
  // real-time acquisition (since the input has lost its original order).
  // NB: because this is trying to emulate the input order, binning is disabled
  if (params->sortInputByAzimuth)
    sortByAzimuth(
      pointCloud, 0, pointCloud.getPointCount(), 0., _angularOrigin);

  // Sanitise the input point cloud
  // todo(df): remove the following with generic handling of properties
  bool codeColour = params->encoder.attributeIdxMap.count("color");
  if (!codeColour)
    pointCloud.removeColors();
  assert(codeColour == pointCloud.hasColors());

  bool codeReflectance = params->encoder.attributeIdxMap.count("reflectance");
  if (!codeReflectance)
    pointCloud.removeReflectances();
  assert(codeReflectance == pointCloud.hasReflectances());

  clock->start();

  if (params->convertColourspace)
    convertFromGbr(params->encoder.sps.attributeSets, pointCloud);

  scaleAttributesForInput(params->encoder.sps.attributeSets, pointCloud);

  // The reconstructed point cloud
  CloudFrame recon;
  auto* reconPtr = params->reconstructedDataPath.empty() ? nullptr : &recon;

  auto bytestreamLenFrameStart = bytestreamFile.tellp();

  int ret = encoder.compress(pointCloud, &params->encoder, this, reconPtr);
  if (ret) {
    cout << "Error: can't compress point cloud!" << endl;
    return -1;
  }

  auto bytestreamLenFrameEnd = bytestreamFile.tellp();
  int frameLen = bytestreamLenFrameEnd - bytestreamLenFrameStart;
  std::cout << "Total frame size " << frameLen << " B" << std::endl;

  clock->stop();

  if (reconPtr)
    writeOutputFrame(params->reconstructedDataPath, {}, recon, recon.cloud);

  return 0;
}

//----------------------------------------------------------------------------

void
SequenceEncoder::onOutputBuffer(const PayloadBuffer& buf)
{
  writeTlv(buf, bytestreamFile);
}

//----------------------------------------------------------------------------

void
SequenceEncoder::onPostRecolour(const PCCPointSet3& cloud)
{
  if (params->postRecolorPath.empty()) {
    return;
  }

  // todo(df): don't allocate if conversion is not required
  PCCPointSet3 tmpCloud(cloud);
  CloudFrame frame;
  frame.setParametersFrom(params->encoder.sps, params->encoder.outputFpBits);
  frame.cloud = cloud;
  frame.frameNum = frameNum - params->firstFrameNum;

  writeOutputFrame(params->postRecolorPath, {}, frame, tmpCloud);
}

//============================================================================

SequenceDecoder::SequenceDecoder(Parameters* params)
  : SequenceCodec(params), decoder(params->decoder)
{}

//----------------------------------------------------------------------------

int
SequenceDecoder::decompress(Stopwatch* clock)
{
  ifstream fin(params->compressedStreamPath, ios::binary);
  if (!fin.is_open()) {
    return -1;
  }

  this->clock = clock;
  clock->start();

  PayloadBuffer buf;
  while (true) {
    PayloadBuffer* buf_ptr = &buf;
    readTlv(fin, &buf);

    // at end of file (or other error), flush decoder
    if (!fin)
      buf_ptr = nullptr;

    if (decoder.decompress(buf_ptr, this)) {
      cout << "Error: can't decompress point cloud!" << endl;
      return -1;
    }

    if (!buf_ptr)
      break;
  }

  fin.clear();
  fin.seekg(0, ios_base::end);
  std::cout << "Total bitstream size " << fin.tellg() << " B" << std::endl;

  clock->stop();

  return 0;
}

//----------------------------------------------------------------------------

void
SequenceDecoder::onOutputCloud(const CloudFrame& frame)
{
  clock->stop();

  // copy the point cloud in order to modify it according to the output options
  PCCPointSet3 pointCloud(frame.cloud);
  writeOutputFrame(
    params->reconstructedDataPath, params->preInvScalePath, frame, pointCloud);

  clock->start();
}

//============================================================================

double
SequenceCodec::outputScale(const CloudFrame& frame) const
{
  switch (params->outputSystem) {
  case OutputSystem::kConformance: return 1.;

  case OutputSystem::kExternal:
    // The scaling converts from the frame's unit length to configured output.
    // In terms of specification this is the external coordinate system.
    return frame.outputUnitLength / params->outputUnitLength;
  }
}

//----------------------------------------------------------------------------

Vec3<double>
SequenceCodec::outputOrigin(const CloudFrame& frame) const
{
  switch (params->outputSystem) {
  case OutputSystem::kConformance: return 0.;

  case OutputSystem::kExternal: return frame.outputOrigin * outputScale(frame);
  }
}

//----------------------------------------------------------------------------

void
SequenceCodec::writeOutputFrame(
  const std::string& postInvScalePath,
  const std::string& preInvScalePath,
  const CloudFrame& frame,
  PCCPointSet3& cloud)
{
  if (postInvScalePath.empty() && preInvScalePath.empty())
    return;

  scaleAttributesForOutput(frame.attrDesc, cloud);

  if (params->convertColourspace)
    convertToGbr(frame.attrDesc, cloud);

  // the order of the property names must be determined from the sps
  ply::PropertyNameMap attrNames;
  attrNames.position = axisOrderToPropertyNames(frame.geometry_axis_order);

  // offset frame number
  int frameNum = frame.frameNum + params->firstFrameNum;

  // Dump the decoded colour using the pre inverse scaled geometry
  if (!preInvScalePath.empty()) {
    std::string filename{expandNum(preInvScalePath, frameNum)};
    ply::write(cloud, attrNames, 1.0, 0.0, filename, !params->outputBinaryPly);
  }

  auto plyScale = outputScale(frame) / (1 << frame.outputFpBits);
  auto plyOrigin = outputOrigin(frame);
  std::string decName{expandNum(postInvScalePath, frameNum)};
  if (!ply::write(
        cloud, attrNames, plyScale, plyOrigin, decName,
        !params->outputBinaryPly)) {
    cout << "Error: can't open output file!" << endl;
  }
}

//============================================================================

const AttributeDescription*
findColourAttrDesc(const std::vector<AttributeDescription>& attrDescs)
{
  // todo(df): don't assume that there is only one colour attribute in the sps
  for (const auto& desc : attrDescs) {
    if (desc.attributeLabel == KnownAttributeLabel::kColour)
      return &desc;
  }
  return nullptr;
}

//----------------------------------------------------------------------------

void
convertToGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  const AttributeDescription* attrDesc = findColourAttrDesc(attrDescs);
  if (!attrDesc)
    return;

  switch (attrDesc->params.cicp_matrix_coefficients_idx) {
  case ColourMatrix::kBt709: convertYCbCrBt709ToGbr(cloud); break;

  case ColourMatrix::kYCgCo:
    // todo(df): select YCgCoR vs YCgCo
    // NB: bitdepth is the transformed bitdepth, not the source
    convertYCgCoRToGbr(attrDesc->bitdepth - 1, cloud);
    break;

  default: break;
  }
}

//----------------------------------------------------------------------------

void
convertFromGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  const AttributeDescription* attrDesc = findColourAttrDesc(attrDescs);
  if (!attrDesc)
    return;

  switch (attrDesc->params.cicp_matrix_coefficients_idx) {
  case ColourMatrix::kBt709: convertGbrToYCbCrBt709(cloud); break;

  case ColourMatrix::kYCgCo:
    // todo(df): select YCgCoR vs YCgCo
    // NB: bitdepth is the transformed bitdepth, not the source
    convertGbrToYCgCoR(attrDesc->bitdepth - 1, cloud);
    break;

  default: break;
  }
}

//============================================================================

const AttributeDescription*
findReflAttrDesc(const std::vector<AttributeDescription>& attrDescs)
{
  // todo(df): don't assume that there is only one in the sps
  for (const auto& desc : attrDescs) {
    if (desc.attributeLabel == KnownAttributeLabel::kReflectance)
      return &desc;
  }
  return nullptr;
}

//----------------------------------------------------------------------------

struct AttrFwdScaler {
  template<typename T>
  T operator()(const AttributeParameters& params, T val) const
  {
    int scale = params.attr_scale_minus1 + 1;
    return ((val - params.attr_offset) << params.attr_frac_bits) / scale;
  }
};

//----------------------------------------------------------------------------

struct AttrInvScaler {
  template<typename T>
  T operator()(const AttributeParameters& params, T val) const
  {
    int scale = params.attr_scale_minus1 + 1;
    return ((val * scale) >> params.attr_frac_bits) + params.attr_offset;
  }
};

//----------------------------------------------------------------------------

template<typename Op>
void
scaleAttributes(
  const std::vector<AttributeDescription>& attrDescs,
  PCCPointSet3& cloud,
  Op scaler)
{
  // todo(df): extend this to other attributes
  const AttributeDescription* attrDesc = findReflAttrDesc(attrDescs);
  if (!attrDesc || !attrDesc->params.scalingParametersPresent)
    return;

  auto& params = attrDesc->params;

  // Parameters present, but nothing to do
  bool unityScale = !params.attr_scale_minus1 && !params.attr_frac_bits;
  if (unityScale && !params.attr_offset)
    return;

  const auto pointCount = cloud.getPointCount();
  for (size_t i = 0; i < pointCount; ++i) {
    auto& val = cloud.getReflectance(i);
    val = scaler(params, val);
  }
}

//----------------------------------------------------------------------------

void
SequenceCodec::scaleAttributesForInput(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  scaleAttributes(attrDescs, cloud, AttrFwdScaler());
}

//----------------------------------------------------------------------------

void
SequenceCodec::scaleAttributesForOutput(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  scaleAttributes(attrDescs, cloud, AttrInvScaler());
}

//============================================================================
