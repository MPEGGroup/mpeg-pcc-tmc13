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

#include "ply.h"

#include "PCCMisc.h"
#include "PCCPointSet.h"

#include <fstream>
#include <string>
#include <vector>

namespace pcc {

//============================================================================

static bool
compareSeparators(char aChar, const char* const sep)
{
  int i = 0;
  while (sep[i] != '\0') {
    if (aChar == sep[i])
      return false;
    i++;
  }
  return true;
}

//============================================================================

static bool
getTokens(
  const char* str, const char* const sep, std::vector<std::string>& tokens)
{
  if (!tokens.empty())
    tokens.clear();
  std::string buf = "";
  size_t i = 0;
  size_t length = ::strlen(str);
  while (i < length) {
    if (compareSeparators(str[i], sep)) {
      buf += str[i];
    } else if (buf.length() > 0) {
      tokens.push_back(buf);
      buf = "";
    }
    i++;
  }
  if (!buf.empty())
    tokens.push_back(buf);
  return !tokens.empty();
}

//============================================================================

bool
ply::write(
  const PCCPointSet3& cloud,
  const PropertyNameMap& attributeNames,
  double positionScale,
  Vec3<double> positionOffset,
  const std::string& fileName,
  bool asAscii)
{
  std::ofstream fout(fileName, std::ofstream::out);
  if (!fout.is_open()) {
    return false;
  }

  const size_t pointCount = cloud.getPointCount();
  fout << "ply" << std::endl;

  if (asAscii) {
    fout << "format ascii 1.0" << std::endl;
  } else {
    PCCEndianness endianess = PCCSystemEndianness();
    if (endianess == PCC_BIG_ENDIAN) {
      fout << "format binary_big_endian 1.0" << std::endl;
    } else {
      fout << "format binary_little_endian 1.0" << std::endl;
    }
  }
  fout << "element vertex " << pointCount << std::endl;
  if (asAscii) {
    fout << "property float " << attributeNames.position[0] << std::endl;
    fout << "property float " << attributeNames.position[1] << std::endl;
    fout << "property float " << attributeNames.position[2] << std::endl;
  } else {
    fout << "property float64 " << attributeNames.position[0] << std::endl;
    fout << "property float64 " << attributeNames.position[1] << std::endl;
    fout << "property float64 " << attributeNames.position[2] << std::endl;
  }

  if (cloud.hasColors()) {
    fout << "property uchar green" << std::endl;
    fout << "property uchar blue" << std::endl;
    fout << "property uchar red" << std::endl;
  }
  if (cloud.hasReflectances()) {
    fout << "property uint16 refc" << std::endl;
  }
  if (cloud.hasFrameIndex()) {
    fout << "property uint8 frameindex" << std::endl;
  }
  fout << "element face 0" << std::endl;
  fout << "property list uint8 int32 vertex_index" << std::endl;
  fout << "end_header" << std::endl;
  if (asAscii) {
    //      fout << std::setprecision(std::numeric_limits<double>::max_digits10);
    fout << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < pointCount; ++i) {
      Vec3<double> position = cloud[i] * positionScale + positionOffset;
      fout << position.x() << " " << position.y() << " " << position.z();
      if (cloud.hasColors()) {
        const Vec3<attr_t>& color = cloud.getColor(i);
        fout << " " << static_cast<int>(color[0]) << " "
             << static_cast<int>(color[1]) << " "
             << static_cast<int>(color[2]);
      }
      if (cloud.hasReflectances()) {
        fout << " " << static_cast<int>(cloud.getReflectance(i));
      }
      if (cloud.hasFrameIndex()) {
        fout << " " << static_cast<int>(cloud.getFrameIndex(i));
      }
      fout << std::endl;
    }
  } else {
    fout.clear();
    fout.close();
    fout.open(fileName, std::ofstream::binary | std::ofstream::app);
    for (size_t i = 0; i < pointCount; ++i) {
      Vec3<double> position = cloud[i] * positionScale + positionOffset;
      fout.write(
        reinterpret_cast<const char* const>(&position), sizeof(double) * 3);
      if (cloud.hasColors()) {
        const Vec3<attr_t>& c = cloud.getColor(i);
        Vec3<uint8_t> val8b{uint8_t(c[0]), uint8_t(c[1]), uint8_t(c[2])};
        fout.write(reinterpret_cast<const char*>(&val8b), sizeof(uint8_t) * 3);
      }
      if (cloud.hasReflectances()) {
        const attr_t& reflectance = cloud.getReflectance(i);
        fout.write(
          reinterpret_cast<const char*>(&reflectance), sizeof(uint16_t));
      }
      if (cloud.hasFrameIndex()) {
        const uint16_t& findex = cloud.getFrameIndex(i);
        fout.write(reinterpret_cast<const char*>(&findex), sizeof(uint16_t));
      }
    }
  }
  fout.close();
  return true;
}

//============================================================================

bool
ply::read(
  const std::string& fileName,
  const PropertyNameMap& attributeNames,
  double positionScale,
  PCCPointSet3& cloud)
{
  std::ifstream ifs(fileName, std::ifstream::in | std::ifstream::binary);
  if (!ifs.is_open()) {
    return false;
  }
  enum AttributeType
  {
    ATTRIBUTE_TYPE_FLOAT64 = 0,
    ATTRIBUTE_TYPE_FLOAT32 = 1,
    ATTRIBUTE_TYPE_UINT64 = 2,
    ATTRIBUTE_TYPE_UINT32 = 3,
    ATTRIBUTE_TYPE_UINT16 = 4,
    ATTRIBUTE_TYPE_UINT8 = 5,
    ATTRIBUTE_TYPE_INT64 = 6,
    ATTRIBUTE_TYPE_INT32 = 7,
    ATTRIBUTE_TYPE_INT16 = 8,
    ATTRIBUTE_TYPE_INT8 = 9,
  };
  struct AttributeInfo {
    std::string name;
    AttributeType type;
    size_t byteCount;
  };

  std::vector<AttributeInfo> attributesInfo;
  attributesInfo.reserve(16);
  const size_t MAX_BUFFER_SIZE = 4096;
  char tmp[MAX_BUFFER_SIZE];
  const char* sep = " \t\r";
  std::vector<std::string> tokens;

  ifs.getline(tmp, MAX_BUFFER_SIZE);
  getTokens(tmp, sep, tokens);
  if (tokens.empty() || tokens[0] != "ply") {
    std::cout << "Error: corrupted file!" << std::endl;
    return false;
  }
  bool isAscii = false;
  double version = 1.0;
  size_t pointCount = 0;
  bool isVertexProperty = true;
  while (1) {
    if (ifs.eof()) {
      std::cout << "Error: corrupted header!" << std::endl;
      return false;
    }
    ifs.getline(tmp, MAX_BUFFER_SIZE);
    getTokens(tmp, sep, tokens);
    if (tokens.empty() || tokens[0] == "comment") {
      continue;
    }
    if (tokens[0] == "format") {
      if (tokens.size() != 3) {
        std::cout << "Error: corrupted format info!" << std::endl;
        return false;
      }
      isAscii = tokens[1] == "ascii";
      version = atof(tokens[2].c_str());
    } else if (tokens[0] == "element") {
      if (tokens.size() != 3) {
        std::cout << "Error: corrupted element info!" << std::endl;
        return false;
      }
      if (tokens[1] == "vertex") {
        pointCount = atoi(tokens[2].c_str());
      } else {
        isVertexProperty = false;
      }
    } else if (tokens[0] == "property" && isVertexProperty) {
      if (tokens.size() != 3) {
        std::cout << "Error: corrupted property info!" << std::endl;
        return false;
      }
      const std::string& propertyType = tokens[1];
      const std::string& propertyName = tokens[2];
      const size_t attributeIndex = attributesInfo.size();
      attributesInfo.resize(attributeIndex + 1);
      AttributeInfo& attributeInfo = attributesInfo[attributeIndex];
      attributeInfo.name = propertyName;
      if (propertyType == "float64") {
        attributeInfo.type = ATTRIBUTE_TYPE_FLOAT64;
        attributeInfo.byteCount = 8;
      } else if (propertyType == "float" || propertyType == "float32") {
        attributeInfo.type = ATTRIBUTE_TYPE_FLOAT32;
        attributeInfo.byteCount = 4;
      } else if (propertyType == "uint64") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT64;
        attributeInfo.byteCount = 8;
      } else if (propertyType == "uint32") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT32;
        attributeInfo.byteCount = 4;
      } else if (propertyType == "uint16") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT16;
        attributeInfo.byteCount = 2;
      } else if (propertyType == "uchar" || propertyType == "uint8") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT8;
        attributeInfo.byteCount = 1;
      } else if (propertyType == "int64") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT64;
        attributeInfo.byteCount = 8;
      } else if (propertyType == "int32") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT32;
        attributeInfo.byteCount = 4;
      } else if (propertyType == "int16") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT16;
        attributeInfo.byteCount = 2;
      } else if (propertyType == "char" || propertyType == "int8") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT8;
        attributeInfo.byteCount = 1;
      }
    } else if (tokens[0] == "end_header") {
      break;
    }
  }
  if (version != 1.0) {
    std::cout << "Error: non-supported version!" << std::endl;
    return false;
  }

  size_t indexX = PCC_UNDEFINED_INDEX;
  size_t indexY = PCC_UNDEFINED_INDEX;
  size_t indexZ = PCC_UNDEFINED_INDEX;
  size_t indexR = PCC_UNDEFINED_INDEX;
  size_t indexG = PCC_UNDEFINED_INDEX;
  size_t indexB = PCC_UNDEFINED_INDEX;
  size_t indexReflectance = PCC_UNDEFINED_INDEX;
  size_t indexFrame = PCC_UNDEFINED_INDEX;
  size_t indexNX = PCC_UNDEFINED_INDEX;
  size_t indexNY = PCC_UNDEFINED_INDEX;
  size_t indexNZ = PCC_UNDEFINED_INDEX;
  size_t indexLaserAngle = PCC_UNDEFINED_INDEX;
  const size_t attributeCount = attributesInfo.size();
  for (size_t a = 0; a < attributeCount; ++a) {
    const auto& attributeInfo = attributesInfo[a];
    if (
      attributeInfo.name == attributeNames.position[0]
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexX = a;
    } else if (
      attributeInfo.name == attributeNames.position[1]
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexY = a;
    } else if (
      attributeInfo.name == attributeNames.position[2]
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexZ = a;
    } else if (attributeInfo.name == "red" && attributeInfo.byteCount == 1) {
      indexR = a;
    } else if (attributeInfo.name == "green" && attributeInfo.byteCount == 1) {
      indexG = a;
    } else if (attributeInfo.name == "blue" && attributeInfo.byteCount == 1) {
      indexB = a;
    } else if (
      (attributeInfo.name == "reflectance" || attributeInfo.name == "refc")
      && attributeInfo.byteCount <= 2) {
      indexReflectance = a;
    } else if (
      attributeInfo.name == "frameindex" && attributeInfo.byteCount <= 2) {
      indexFrame = a;
    } else if (
      attributeInfo.name == "nx"
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexNX = a;
    } else if (
      attributeInfo.name == "ny"
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexNY = a;
    } else if (
      attributeInfo.name == "nz"
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexNZ = a;
    } else if (attributeInfo.name == "laserangle") {
      indexLaserAngle = a;
    }
  }
  if (
    indexX == PCC_UNDEFINED_INDEX || indexY == PCC_UNDEFINED_INDEX
    || indexZ == PCC_UNDEFINED_INDEX) {
    std::cout << "Error: missing coordinates!" << std::endl;
    return false;
  }
  bool withColors = indexR != PCC_UNDEFINED_INDEX
    && indexG != PCC_UNDEFINED_INDEX && indexB != PCC_UNDEFINED_INDEX;
  bool withReflectances = indexReflectance != PCC_UNDEFINED_INDEX;
  bool withFrameIndex = indexFrame != PCC_UNDEFINED_INDEX;
  bool withLaserAngles = indexLaserAngle != PCC_UNDEFINED_INDEX;

  cloud.addRemoveAttributes(withColors, withReflectances);
  if (withFrameIndex)
    cloud.addFrameIndex();
  else
    cloud.removeFrameIndex();

  if (withLaserAngles)
    cloud.addLaserAngles();
  else
    cloud.removeLaserAngles();

  cloud.resize(pointCount);
  if (isAscii) {
    size_t pointCounter = 0;
    while (!ifs.eof() && pointCounter < pointCount) {
      ifs.getline(tmp, MAX_BUFFER_SIZE);
      getTokens(tmp, sep, tokens);
      if (tokens.empty()) {
        continue;
      }
      if (tokens.size() < attributeCount) {
        return false;
      }
      auto& position = cloud[pointCounter];
      position[0] = atof(tokens[indexX].c_str()) * positionScale;
      position[1] = atof(tokens[indexY].c_str()) * positionScale;
      position[2] = atof(tokens[indexZ].c_str()) * positionScale;
      if (cloud.hasColors()) {
        auto& color = cloud.getColor(pointCounter);
        color[0] = atoi(tokens[indexG].c_str());
        color[1] = atoi(tokens[indexB].c_str());
        color[2] = atoi(tokens[indexR].c_str());
      }
      if (cloud.hasReflectances()) {
        cloud.getReflectance(pointCounter) =
          uint16_t(atoi(tokens[indexReflectance].c_str()));
      }
      if (cloud.hasFrameIndex()) {
        cloud.getFrameIndex(pointCounter) =
          uint8_t(atoi(tokens[indexFrame].c_str()));
      }
      if (cloud.hasLaserAngles()) {
        cloud.getLaserAngle(pointCounter) =
          std::round(atof(tokens[indexLaserAngle].c_str()));
      }
      ++pointCounter;
    }
  } else {
    for (size_t pointCounter = 0; pointCounter < pointCount && !ifs.eof();
         ++pointCounter) {
      auto& position = cloud[pointCounter];
      for (size_t a = 0; a < attributeCount && !ifs.eof(); ++a) {
        const auto& attributeInfo = attributesInfo[a];
        if (a == indexX) {
          if (attributeInfo.byteCount == 4) {
            float x;
            ifs.read(reinterpret_cast<char*>(&x), sizeof(float));
            position[0] = x * positionScale;
          } else {
            double x;
            ifs.read(reinterpret_cast<char*>(&x), sizeof(double));
            position[0] = x * positionScale;
          }
        } else if (a == indexY) {
          if (attributeInfo.byteCount == 4) {
            float y;
            ifs.read(reinterpret_cast<char*>(&y), sizeof(float));
            position[1] = y * positionScale;
          } else {
            double y;
            ifs.read(reinterpret_cast<char*>(&y), sizeof(double));
            position[1] = y * positionScale;
          }
        } else if (a == indexZ) {
          if (attributeInfo.byteCount == 4) {
            float z;
            ifs.read(reinterpret_cast<char*>(&z), sizeof(float));
            position[2] = z * positionScale;
          } else {
            double z;
            ifs.read(reinterpret_cast<char*>(&z), sizeof(double));
            position[2] = z * positionScale;
          }
        } else if (a == indexR && attributeInfo.byteCount == 1) {
          uint8_t val8b;
          ifs.read(reinterpret_cast<char*>(&val8b), sizeof(uint8_t));
          cloud.getColor(pointCounter)[2] = val8b;
        } else if (a == indexG && attributeInfo.byteCount == 1) {
          uint8_t val8b;
          ifs.read(reinterpret_cast<char*>(&val8b), sizeof(uint8_t));
          cloud.getColor(pointCounter)[0] = val8b;
        } else if (a == indexB && attributeInfo.byteCount == 1) {
          uint8_t val8b;
          ifs.read(reinterpret_cast<char*>(&val8b), sizeof(uint8_t));
          cloud.getColor(pointCounter)[1] = val8b;
        } else if (a == indexReflectance && attributeInfo.byteCount <= 2) {
          if (attributeInfo.byteCount == 1) {
            uint8_t reflectance;
            ifs.read(reinterpret_cast<char*>(&reflectance), sizeof(uint8_t));
            cloud.getReflectance(pointCounter) = reflectance;
          } else {
            auto& reflectance = cloud.getReflectance(pointCounter);
            ifs.read(reinterpret_cast<char*>(&reflectance), sizeof(uint16_t));
          }
        } else if (a == indexFrame && attributeInfo.byteCount <= 2) {
          if (attributeInfo.byteCount == 1) {
            auto& findex = cloud.getFrameIndex(pointCounter);
            ifs.read(reinterpret_cast<char*>(&findex), sizeof(uint8_t));
          } else {
            uint16_t findex;
            ifs.read(reinterpret_cast<char*>(&findex), sizeof(uint16_t));
            cloud.getFrameIndex(pointCounter) = uint8_t(findex);
          }
        } else {
          char buffer[128];
          ifs.read(buffer, attributeInfo.byteCount);
        }
      }
    }
  }
  return true;
}

//============================================================================

}  // namespace pcc
