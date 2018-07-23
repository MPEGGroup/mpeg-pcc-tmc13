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
#include "program_options_lite.h"
#include "version.h"

using namespace std;
using namespace pcc;

int main(int argc, char *argv[]) {
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
  if (params.mode == CODEC_MODE_ENCODE || params.mode == CODEC_MODE_ENCODE_LOSSLESS_GEOMETRY ||
      params.mode == CODEC_MODE_ENCODE_TRISOUP_GEOMETRY) {
    ret = Compress(params, clock_user);
  } else {
    ret = Decompress(params, clock_user);
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

template <typename T>
static std::istream& readUInt(std::istream &in, T &val) {
  unsigned int tmp;
  in >> tmp;
  val = T(tmp);
  return in;
}

static std::istream& operator>>(std::istream &in, CodecMode &val) {
  return readUInt(in, val);
}

static std::istream& operator>>(std::istream &in, ColorTransform &val) {
  return readUInt(in, val);
}

namespace pcc {
static std::istream& operator>>(std::istream &in, TransformType &val) {
  return readUInt(in, val);
}}

namespace pcc {
static std::istream& operator>>(std::istream &in, GeometryCodecType &val) {
  return readUInt(in, val);
}}

namespace pcc {
static std::ostream& operator<<(std::ostream &out, const TransformType &val) {
  switch (val) {
  case TransformType::kIntegerLift: out << "0 (IntegerLifting)"; break;
  case TransformType::kRAHT:        out << "1 (RAHT)"; break;
  }
  return out;
}}

namespace pcc {
static std::ostream& operator<<(std::ostream &out, const GeometryCodecType &val) {
  switch (val) {
  case GeometryCodecType::kBypass:  out << "0 (Bypass)"; break;
  case GeometryCodecType::kOctree:  out << "1 (TMC1 Octree)"; break;
  case GeometryCodecType::kTriSoup: out << "2 (TMC3 TriSoup)"; break;
  }
  return out;
}}

//---------------------------------------------------------------------------
// :: Command line / config parsing

bool ParseParameters(int argc, char *argv[], Parameters &params) {

  namespace po = df::program_options_lite;

  PCCAttributeEncodeParamaters params_attr;
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
      params.encodeParameters.attributeEncodeParameters[name] = params_attr;
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

  ("mode", params.mode, CODEC_MODE_ENCODE,
    "The encoding/decoding mode:\n"
    "  0: encode\n"
    "  1: decode\n"
    // NB: the following forms are deprecated
    "  2: encode with lossless geometry\n"
    "  3: decode with lossless geometry")

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
    params.encodeParameters.postRecolorPath, {},
    "Recolored pointcloud file path (encoder only)")

  ("preInvScalePath",
    params.decodeParameters.preInvScalePath, {},
    "Pre inverse scaled pointcloud file path (decoder only)")

  // general
  ("colorTransform",
    params.colorTransform, COLOR_TRANSFORM_RGB_TO_YCBCR,
    "The colour transform to be applied:\n"
    "  0: none\n"
    "  1: RGB to YCbCr (Rec.709)")

  ("positionQuantizationScale",
    params.encodeParameters.positionQuantizationScale, 1.,
    "Scale factor to be applied to point positions during quantization process")

  ("mergeDuplicatedPoints",
    params.encodeParameters.mergeDuplicatedPoints, true,
    "Enables removal of duplicated points")

  ("roundOutputPositions",
    params.decodeParameters.roundOutputPositions, false,
    "todo(kmammou)")

  // tools
  ("geometryCodec",
    params.encodeParameters.geometryCodec, GeometryCodecType::kOctree,
    "Controls the method used to encode geometry:"
    "  0: bypass (a priori)\n"
    "  1: octree (TMC3)\n"
    "  2: trisoup (TMC1)")

  ("neighbourContextualisation",
    params.encodeParameters.neighbourContextsEnabled, true,
    "Contextualise geometry octree occupancy based on neighbour patterns")

  ("inferredDirectCodingMode",
    params.encodeParameters.inferredDirectCodingModeEnabled, true,
    "Permits early termination of the geometry octree for isolated points")

  // (trisoup) geometry parameters
  ("triSoupDepth",  // log2(maxBB+1), where maxBB+1 is analogous to image width
    params.encodeParameters.triSoup.depth, 10,
    "Depth of voxels (reconstructed points) in trisoup geometry")

  ("triSoupLevel",
    params.encodeParameters.triSoup.level, 7,
    "Level of triangles (reconstructed surface) in trisoup geometry")

  ("triSoupIntToOrigScale",  // reciprocal of positionQuantizationScale
    params.encodeParameters.triSoup.intToOrigScale, 1.,
    "orig_coords = integer_coords * intToOrigScale + intToOrigTranslation")

  ("triSoupIntToOrigTranslation",
    params.encodeParameters.triSoup.intToOrigTranslation, {0., 0., 0.},
    "orig_coords = integer_coords * intToOrigScale + intToOrigTranslation")

  // attribute processing
  //   NB: Attribute options are special in the way they are applied (see above)
  ("attribute",
    attribute_setter,
    "Encode the given attribute (NB, must appear after the"
    "following attribute parameters)")

  ("searchRange",
    params_attr.searchRange, size_t(0),
    "Attribute's todo(kmammou)")

  ("transformType",
    params_attr.transformType, TransformType::kIntegerLift,
    "Coding method to use for attribute:\n"
    "  0: Nearest neighbour prediction with integer lifting transform\n"
    "  1: Region Adaptive Hierarchical Transform (RAHT)")

  ("rahtLeafDecimationDepth",
    params_attr.binaryLevelThresholdRaht, 3,
    "Sets coefficients to zero in the bottom n levels of RAHT tree. "
    "Used for chroma-subsampling in attribute=color only.")

  ("rahtQuantizationStep",
    params_attr.quantizationStepRaht, 1,
    "Quantization step size used in RAHT")

  ("rahtDepth",
    params_attr.depthRaht, 21,
    "Number of bits for morton representation of an RAHT co-ordinate"
    "component")

  ("numberOfNearestNeighborsInPrediction",
    params_attr.numberOfNearestNeighborsInPrediction, size_t(4),
    "Attribute's maximum number of nearest neighbors to be used for prediction")

  ("levelOfDetailCount",
    params_attr.levelOfDetailCount, size_t(6),
    "Attribute's number of levels of detail")

  ("quantizationSteps",
    params_attr.quantizationSteps, {},
    "Attribute's list of quantization step sizes (one for each LoD)")

  ("dist2", params_attr.dist2, {},
    "Attribute's list of squared distances (one for each LoD)")
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

  // Set GeometryCodecType according to codec mode
  // NB: for bypass, the decoder must load the a priori geometry
  if (params.mode == 2 || params.mode == 3) {
    params.encodeParameters.geometryCodec = GeometryCodecType::kBypass;
  }
  if (params.mode == 4) {
    params.encodeParameters.geometryCodec = GeometryCodecType::kTriSoup;
  }

  // Restore params.mode to be encode vs decode
  params.mode = CodecMode(params.mode & 1);

  // For trisoup, ensure that positionQuantizationScale is the exact inverse of intToOrigScale.
  if (params.encodeParameters.geometryCodec == GeometryCodecType::kTriSoup) {
    params.encodeParameters.positionQuantizationScale =
        1.0 / params.encodeParameters.triSoup.intToOrigScale;
  }

  // For RAHT, ensure that the unused lod count = 0 (prevents mishaps)
  for (auto &attr : params.encodeParameters.attributeEncodeParameters) {
    if (attr.second.transformType != TransformType::kIntegerLift) {
      attr.second.levelOfDetailCount = 0;
    }
  }

  // sanity checks
  //  - validate that quantizationSteps, dist2
  //    of each attribute contain levelOfDetailCount elements.
  for (const auto &attr : params.encodeParameters.attributeEncodeParameters) {
    if (attr.second.transformType == TransformType::kIntegerLift) {
      int lod = attr.second.levelOfDetailCount;

      if (lod > 255) {
        err.error() << attr.first
          << ".levelOfDetailCount must be less than 256\n";
      }
      if (attr.second.dist2.size() != lod) {
        err.error() << attr.first << ".dist2 does not have " << lod << " entries\n";
      }
      if (attr.second.quantizationSteps.size() != lod) {
        err.error() << attr.first << ".quantizationSteps does not have " << lod << " entries\n";
      }
      if (attr.second.numberOfNearestNeighborsInPrediction >
          PCCTMC3MaxPredictionNearestNeighborCount)
      {
        err.error() << attr.first
          << ".numberOfNearestNeighborsInPrediction must be less than "
          << PCCTMC3MaxPredictionNearestNeighborCount << "\n";
      }
    }
  }

  // check required arguments are specified

  const bool encode = params.mode == CODEC_MODE_ENCODE;

  if (encode && params.uncompressedDataPath.empty())
    err.error() << "uncompressedDataPath not set\n";

  if (!encode && params.reconstructedDataPath.empty())
    err.error() << "reconstructedDataPath not set\n";

  if (params.compressedStreamPath.empty())
    err.error() << "compressedStreamPath not set\n";

  // currently the attributes with lossless geometry require the source data
  // todo(?): remove this dependency by improving reporting
  if (params.encodeParameters.geometryCodec == GeometryCodecType::kBypass
   && params.uncompressedDataPath.empty())
    err.error() << "uncompressedDataPath not set\n";

  // report the current configuration (only in the absence of errors so
  // that errors/warnings are more obvious and in the same place).
  if (err.is_errored)
    return false;

  cout << "+ Parameters" << endl;
  if (params.mode == CODEC_MODE_ENCODE) {
    cout << "\t mode                        encode\n";
  } else if (params.mode == CODEC_MODE_DECODE) {
    cout << "\t mode                        decode\n";
  }
  cout << "\t uncompressedDataPath        " << params.uncompressedDataPath << endl;
  cout << "\t compressedStreamPath        " << params.compressedStreamPath << endl;
  cout << "\t reconstructedDataPath       " << params.reconstructedDataPath << endl;
  cout << "\t colorTransform              " << params.colorTransform << endl;
  if (encode) {
    cout << "\t geometryCodec               "
         << params.encodeParameters.geometryCodec << endl;

    cout << "\t positionQuantizationScale   " << params.encodeParameters.positionQuantizationScale
         << endl;

    if (params.encodeParameters.geometryCodec == GeometryCodecType::kOctree) {
      cout << "\t mergeDuplicatedPoints       "
           << params.encodeParameters.mergeDuplicatedPoints << '\n';
      cout << "\t neighbourContextualisation  "
           << params.encodeParameters.neighbourContextsEnabled << '\n';
      cout << "\t inferredDirectCodingMode    "
           << params.encodeParameters.inferredDirectCodingModeEnabled << '\n';
    }

    if (params.encodeParameters.geometryCodec == GeometryCodecType::kTriSoup) {
      cout << "\t triSoupDepth                " << params.encodeParameters.triSoup.depth << endl;
      cout << "\t triSoupLevel                " << params.encodeParameters.triSoup.level << endl;
      cout << "\t triSoupIntToOrigScale       " << params.encodeParameters.triSoup.intToOrigScale << endl;
      cout << "\t triSoupIntToOrigTranslation ";
      for (const auto tr : params.encodeParameters.triSoup.intToOrigTranslation) {
        cout << tr << " ";
      }
      cout << endl;
    }
    for (const auto & attributeEncodeParameters : params.encodeParameters.attributeEncodeParameters) {
      cout << "\t " << attributeEncodeParameters.first << endl;
      cout << "\t\t transformType                          "
           << attributeEncodeParameters.second.transformType << endl;

      if (attributeEncodeParameters.second.transformType == TransformType::kRAHT) {
        cout << "\t\t rahtQuantizationStep                   "
             << attributeEncodeParameters.second.quantizationStepRaht << endl;
        cout << "\t\t rahtDepth                              "
             << attributeEncodeParameters.second.depthRaht << endl;
      }

      if (attributeEncodeParameters.second.transformType == TransformType::kIntegerLift) {
        cout << "\t\t numberOfNearestNeighborsInPrediction   "
             << attributeEncodeParameters.second.numberOfNearestNeighborsInPrediction << endl;
        cout << "\t\t searchRange                            "
             << attributeEncodeParameters.second.searchRange << endl;
        cout << "\t\t levelOfDetailCount                     "
             << attributeEncodeParameters.second.levelOfDetailCount << endl;
        cout << "\t\t dist2                                  ";
        for (const auto qs : attributeEncodeParameters.second.dist2) {
          cout << qs << " ";
        }
        cout << endl;
        cout << "\t\t quantizationSteps                      ";
        for (const auto qs : attributeEncodeParameters.second.quantizationSteps) {
          cout << qs << " ";
        }
        cout << endl;
      }
    }
    cout << endl;
  } else {
      cout << "\t roundOutputPositions        "
           << params.decodeParameters.roundOutputPositions << '\n';
  }

  return true;
}

int Compress(const Parameters &params, Stopwatch& clock) {
  PCCPointSet3 pointCloud;
  if (!pointCloud.read(params.uncompressedDataPath) || pointCloud.getPointCount() == 0) {
    cout << "Error: can't open input file!" << endl;
    return -1;
  }

  clock.start();

  if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
    pointCloud.convertRGBToYUV();
  }
  PCCTMC3Encoder3 encoder;
  PCCBitstream bitstream = {};
  const size_t predictedBitstreamSize =
      encoder.estimateBitstreamSize(pointCloud, params.encodeParameters);
  std::unique_ptr<uint8_t[]> buffer(new uint8_t[predictedBitstreamSize]);
  bitstream.buffer = buffer.get();
  bitstream.capacity = predictedBitstreamSize;
  bitstream.size = 0;

  std::unique_ptr<PCCPointSet3> reconstructedPointCloud;
  if (!params.reconstructedDataPath.empty()) {
    reconstructedPointCloud.reset(new PCCPointSet3);
  }

  int ret = encoder.compress(
      pointCloud, params.encodeParameters, bitstream,
      reconstructedPointCloud.get());
  if (ret) {
    cout << "Error: can't compress point cloud!" << endl;
    return -1;
  }

  clock.stop();

  assert(bitstream.size <= bitstream.capacity);
  std::cout << "Total bitstream size " << bitstream.size << " B" << std::endl;
  ofstream fout(params.compressedStreamPath, ios::binary);
  if (!fout.is_open()) {
    return -1;
  }
  fout.write(reinterpret_cast<const char *>(bitstream.buffer), bitstream.size);
  fout.close();

  if (!params.reconstructedDataPath.empty()) {
    if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
      reconstructedPointCloud->convertYUVToRGB();
    }
    reconstructedPointCloud->write(params.reconstructedDataPath, true);
  }

  return 0;
}
int Decompress(const Parameters &params, Stopwatch &clock) {
  PCCBitstream bitstream = {};
  ifstream fin(params.compressedStreamPath, ios::binary);
  if (!fin.is_open()) {
    return -1;
  }
  fin.seekg(0, std::ios::end);
  uint64_t bitStreamSize = fin.tellg();
  fin.seekg(0, std::ios::beg);
  unique_ptr<uint8_t[]> buffer(new uint8_t[bitStreamSize]);
  bitstream.buffer = buffer.get();
  bitstream.capacity = bitStreamSize;
  bitstream.size = 0;
  fin.read(reinterpret_cast<char *>(bitstream.buffer), bitStreamSize);
  if (!fin) {
    return -1;
  }
  fin.close();

  clock.start();

  PCCTMC3Decoder3 decoder;
  PCCPointSet3 pointCloud;

  // read a priori geometry from input file for bypass case
  if (params.encodeParameters.geometryCodec == GeometryCodecType::kBypass) {
    if (!pointCloud.read(params.uncompressedDataPath) || pointCloud.getPointCount() == 0) {
      cout << "Error: can't open input file!" << endl;
      return -1;
    }
    pointCloud.removeReflectances();
    pointCloud.removeColors();
  }

  int ret = decoder.decompress(params.decodeParameters, bitstream, pointCloud);
  if (ret) {
    cout << "Error: can't decompress point cloud!" << endl;
    return -1;
  }
  assert(bitstream.size <= bitstream.capacity);
  std::cout << "Total bitstream size " << bitstream.size << " B" << std::endl;

  if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
    pointCloud.convertYUVToRGB();
  }

  clock.stop();

  if (!pointCloud.write(params.reconstructedDataPath, true)) {
    cout << "Error: can't open output file!" << endl;
    return -1;
  }
  return 0;
}
