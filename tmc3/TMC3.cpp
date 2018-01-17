/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * <OWNER> = Apple Inc.
 * <ORGANIZATION> = Apple Inc.
 * <YEAR> = 2017
 *
 * Copyright (c) 2017, Apple Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "TMC3.h"

using namespace std;
using namespace pcc;

int main(int argc, char *argv[]) {
  Parameters params;
  Usage();
  if (!ParseParameters(argc, argv, params)) {
    std::cout << "Error: missing parameters!" << std::endl;
    return -1;
  }
  const auto start = std::chrono::high_resolution_clock::now();
  int ret = 0;
  if (params.mode == CODEC_MODE_ENCODE || params.mode == CODEC_MODE_ENCODE_LOSSLESS_GEOMETRY) {
    ret = Compress(params);
  } else {
    ret = Decompress(params);
  }

  const auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Processing time: "
            << std::chrono::duration<double, std::milli>(end - start).count() / 1000.0 << " s"
            << std::endl;
  return ret;
}
void Usage() {
  std::cout << "tmc3 v" << TMC3_VERSION_MAJOR << "." << TMC3_VERSION_MAJOR << std::endl
            << std::endl;

  std::cout << "+ Usage" << std::endl;
  std::cout << "\t Encode example: \n tmc3 --mode 0 --mergeDuplicatedPoints 1 "
               "--uncompressedDataPath Ford_01-0100.ply --compressedStreamPath compressed.bin "
               "--colorTransform 1 --numberOfNearestNeighborsInPrediction 8 "
               "--positionQuantizationScale 500.0 --levelOfDetailCount 6 --dist2 16777216 4194304 "
               "1048576 262144 65536 0 --quantizationSteps 1 1 1 2 2 2 --quantizationDeadZoneSizes "
               "1 1 1 2 2 2 --searchRange 2 --attribute color --attribute reflectance"
            << std::endl
            << std::endl;
  std::cout << "\t Decode example: \n tmc3 --mode 1 --compressedStreamPath compressed.bin "
               "--reconstructedDataPath reconstructed_dec.ply --colorTransform 1"
            << std::endl
            << std::endl;
  std::cout << std::endl;
}

bool ParseParameters(int argc, char *argv[], Parameters &params) {
  PCCAttributeEncodeParamaters attributeEncodeParams;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--mode")) {
      if (++i < argc) params.mode = atoi(argv[i]);
    } else if (!strcmp(argv[i], "--reconstructedDataPath")) {
      if (++i < argc) params.reconstructedDataPath = argv[i];
    } else if (!strcmp(argv[i], "--uncompressedDataPath")) {
      if (++i < argc) params.uncompressedDataPath = argv[i];
    } else if (!strcmp(argv[i], "--compressedStreamPath")) {
      if (++i < argc) params.compressedStreamPath = argv[i];
    } else if (!strcmp(argv[i], "--attribute")) {
      if (++i < argc)
        params.encodeParameters.attributeEncodeParameters[argv[i]] = attributeEncodeParams;
    } else if (!strcmp(argv[i], "--numberOfNearestNeighborsInPrediction")) {
      if (++i < argc) attributeEncodeParams.numberOfNearestNeighborsInPrediction = atoi(argv[i]);
    } else if (!strcmp(argv[i], "--levelOfDetailCount")) {
      if (++i < argc) {
        attributeEncodeParams.levelOfDetailCount = atoi(argv[i]);
        attributeEncodeParams.quantizationSteps.resize(attributeEncodeParams.levelOfDetailCount);
        attributeEncodeParams.quantizationDeadZoneSizes.resize(
            attributeEncodeParams.levelOfDetailCount);
        attributeEncodeParams.dist2.resize(attributeEncodeParams.levelOfDetailCount);
      }
    } else if (!strcmp(argv[i], "--quantizationSteps")) {
      for (size_t k = 0; k < attributeEncodeParams.levelOfDetailCount; ++k) {
        if (++i < argc) attributeEncodeParams.quantizationSteps[k] = atoi(argv[i]);
      }
    } else if (!strcmp(argv[i], "--quantizationDeadZoneSizes")) {
      for (size_t k = 0; k < attributeEncodeParams.levelOfDetailCount; ++k) {
        if (++i < argc) attributeEncodeParams.quantizationDeadZoneSizes[k] = atoi(argv[i]);
      }
    } else if (!strcmp(argv[i], "--dist2")) {
      for (size_t k = 0; k < attributeEncodeParams.levelOfDetailCount; ++k) {
        if (++i < argc) attributeEncodeParams.dist2[k] = atoi(argv[i]);
      }
    } else if (!strcmp(argv[i], "--searchRange")) {
      if (++i < argc) attributeEncodeParams.searchRange = atoi(argv[i]);
    } else if (!strcmp(argv[i], "--colorTransform")) {
      if (++i < argc) params.colorTransform = static_cast<ColorTransform>(atoi(argv[i]));
    } else if (!strcmp(argv[i], "--positionQuantizationScale")) {
      if (++i < argc) params.encodeParameters.positionQuantizationScale = atof(argv[i]);
    } else if (!strcmp(argv[i], "--mergeDuplicatedPoints")) {
      if (++i < argc) params.encodeParameters.mergeDuplicatedPoints = atoi(argv[i]) != 0;
    } else if (!strcmp(argv[i], "--roundOutputPositions")) {
        if (++i < argc) params.roundOutputPositions = atoi(argv[i]) != 0;
    }
  }

  const bool encode =
      (params.mode == CODEC_MODE_ENCODE || params.mode == CODEC_MODE_ENCODE_LOSSLESS_GEOMETRY);
  cout << "+ Parameters" << endl;
  cout << "\t mode                        ";
  if (params.mode == CODEC_MODE_ENCODE) {
    cout << "encode" << endl;
  } else if (params.mode == CODEC_MODE_ENCODE_LOSSLESS_GEOMETRY) {
    cout << "encode with lossless geometry" << endl;
  } else if (params.mode == CODEC_MODE_DECODE) {
    cout << "decode" << endl;
  } else if (params.mode == CODEC_MODE_DECODE_LOSSLESS_GEOMETRY) {
    cout << "decode with lossless geometry" << endl;
  }
  cout << "\t uncompressedDataPath        " << params.uncompressedDataPath << endl;
  cout << "\t compressedStreamPath        " << params.compressedStreamPath << endl;
  cout << "\t reconstructedDataPath       " << params.reconstructedDataPath << endl;
  cout << "\t colorTransform              " << params.colorTransform << endl;
  if (encode) {
    cout << "\t mergeDuplicatedPoints       " << params.encodeParameters.mergeDuplicatedPoints
         << endl;
    cout << "\t positionQuantizationScale   " << params.encodeParameters.positionQuantizationScale
         << endl;
    for (const auto & attributeEncodeParameters : params.encodeParameters.attributeEncodeParameters) {
      cout << "\t " << attributeEncodeParameters.first << endl;
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
      cout << "\t\t quantizationDeadZoneSizes              ";
      for (const auto dz : attributeEncodeParameters.second.quantizationDeadZoneSizes) {
        cout << dz << " ";
      }
      cout << endl;
    }
    cout << endl;
  } else {
      cout << "\t roundOutputPositions        " << params.roundOutputPositions << endl;
  }

  const bool test1 =
      encode && (params.uncompressedDataPath.empty() || params.compressedStreamPath.empty());
  const bool test2 =
      !encode && (params.reconstructedDataPath.empty() || params.compressedStreamPath.empty());
  const bool test3 =
      params.mode == CODEC_MODE_DECODE_LOSSLESS_GEOMETRY && params.uncompressedDataPath.empty();
  if (test1 || test2 || test3) {
    return false;
  }
  return true;
}
int Compress(const Parameters &params) {
  PCCPointSet3 pointCloud;
  if (!pointCloud.read(params.uncompressedDataPath) || pointCloud.getPointCount() == 0) {
    cout << "Error: can't open input file!" << endl;
    return -1;
  }

  if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
    pointCloud.convertRGBToYUV();
  }
  PCCTMC3Encoder3 encoder;
  PCCBitstream bitstream = {};
  const size_t predictedBitstreamSize =
      encoder.estimateBitstreamSize(pointCloud, params.encodeParameters);
  std::unique_ptr<uint8_t> buffer(new uint8_t[predictedBitstreamSize]);
  bitstream.buffer = buffer.get();
  bitstream.capacity = predictedBitstreamSize;
  bitstream.size = 0;

  std::unique_ptr<PCCPointSet3> reconstructedPointCloud;
  if (!params.reconstructedDataPath.empty()) {
    reconstructedPointCloud.reset(new PCCPointSet3);
  }

  if ((params.mode == CODEC_MODE_ENCODE &&
       encoder.compress(pointCloud, params.encodeParameters, bitstream,
                        reconstructedPointCloud.get())) ||
      (params.mode == CODEC_MODE_ENCODE_LOSSLESS_GEOMETRY &&
       encoder.compressWithLosslessGeometry(pointCloud, params.encodeParameters, bitstream,
                                            reconstructedPointCloud.get()))) {
    cout << "Error: can't compress point cloud!" << endl;
    return -1;
  }

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
int Decompress(const Parameters &params) {
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

  PCCTMC3Decoder3 decoder;
  PCCPointSet3 pointCloud;
  if (params.mode == CODEC_MODE_DECODE_LOSSLESS_GEOMETRY) {  // read geometry from input file
    if (!pointCloud.read(params.uncompressedDataPath) || pointCloud.getPointCount() == 0) {
      cout << "Error: can't open input file!" << endl;
      return -1;
    }
    pointCloud.removeReflectances();
    pointCloud.removeColors();
  }

  if ((params.mode == CODEC_MODE_DECODE && decoder.decompress(bitstream, pointCloud, params.roundOutputPositions)) ||
      (params.mode == CODEC_MODE_DECODE_LOSSLESS_GEOMETRY &&
       decoder.decompressWithLosslessGeometry(bitstream, pointCloud))) {
    cout << "Error: can't decompress point cloud!" << endl;
    return -1;
  }
  assert(bitstream.size <= bitstream.capacity);
  std::cout << "Total bitstream size " << bitstream.size << " B" << std::endl;

  if (params.colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
    pointCloud.convertYUVToRGB();
  }

  if (!pointCloud.write(params.reconstructedDataPath, true)) {
    cout << "Error: can't open output file!" << endl;
    return -1;
  }
  return 0;
}
