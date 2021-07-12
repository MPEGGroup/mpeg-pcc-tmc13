/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2019, ISO/IEC
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

#include <climits>
#include <exception>

#include "PCCMisc.h"
#include "PCCPointSet.h"
#include "ply.h"
#include "program_options_lite.h"
#include "version.h"

using namespace std;
using namespace pcc;

//============================================================================

struct Options;

bool parseParameters(int argc, char* argv[], Options& opts);
void runMerge(const Options& opts);
void runSplit(const Options& opts);

//============================================================================

struct Options {
  enum class Mode
  {
    Merge,
    Split
  } mode;

  // output mode for ply writing (binary or ascii)
  bool outputBinaryPly;

  // path (using %d to indicate frame number) of input files
  std::string srcPath;

  // path (using %d to indicate frame number) of output files
  std::string outPath;

  // number of src frames to process;
  int frameCount;

  // first frame number in input sequence
  int firstFrameNum;

  // first frame number in output sequence
  int firstOutputFrameNum;

  // number of src frames per merged output frame
  int groupSize;
};

//============================================================================

static std::istream&
operator>>(std::istream& in, Options::Mode& val)
{
  std::string word;
  in >> word;
  if (word == "merge")
    val = Options::Mode::Merge;
  else if (word == "split")
    val = Options::Mode::Split;
  else
    throw std::exception();
  return in;
}

//----------------------------------------------------------------------------

static std::ostream&
operator<<(std::ostream& out, const Options::Mode& val)
{
  switch (val) {
  case Options::Mode::Merge: out << "merge"; break;
  case Options::Mode::Split: out << "split"; break;
  }
  return out;
}

//============================================================================

int
main(int argc, char* argv[])
{
  cout << "MPEG PCC ply merge/split tool from Test Model C13" << endl;

  Options opts;
  if (!parseParameters(argc, argv, opts)) {
    return 1;
  }

  try {
    switch (opts.mode) {
    case Options::Mode::Merge: runMerge(opts); break;
    case Options::Mode::Split: runSplit(opts); break;
    }
  }
  catch (const exception& e) {
    cerr << "Error:" << e.what() << endl;
  }

  return 0;
}

//---------------------------------------------------------------------------
// :: Command line / config parsing

bool
parseParameters(int argc, char* argv[], Options& params)
{
  namespace po = df::program_options_lite;
  bool print_help = false;

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

  ("mode", params.mode, Options::Mode::Merge,
    "The combine/split mode:\n"
    "  split: extract n ply files from one\n"
    "  merge: combine n ply files into one")

  ("srcPath",
    params.srcPath, {},
    "Input pointcloud file path")

  ("outPath",
    params.outPath, {},
    "Output pointcloud file path")

  ("outputBinaryPly",
    params.outputBinaryPly, false,
    "Output ply files using binary (or otherwise ascii) format")

  ("firstFrameNum",
    params.firstFrameNum, 0,
    "Number of first frame of input sequence (used in %d interpolation)")

  ("firstOutputFrameNum",
    params.firstOutputFrameNum, 0,
    "Number of first frame of output sequence (used in %d interpolation)")

  ("frameCount",
    params.frameCount, 8,
    "Number of source ply files to process")

  ("groupSize",
    params.groupSize, 8,
    "Number of source ply files per combined output (merge mode only)")
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

  po::dumpCfg(cout, opts, 4);

  return !err.is_errored;
}

//---------------------------------------------------------------------------
// Read a series of frames, write a single merged frame.
//

void
runMerge(const Options& opts)
{
  ply::PropertyNameMap propNames;
  propNames.position = {"x", "y", "z"};

  // iterate over all input frames in groups
  int outFrameNum = opts.firstOutputFrameNum;
  for (int i = 0; i < opts.frameCount; outFrameNum++) {
    vector<PCCPointSet3> srcClouds;
    srcClouds.reserve(opts.groupSize);

    // read source frames
    for (int j = 0; j < opts.groupSize && i < opts.frameCount; j++, i++) {
      int frameNum = opts.firstFrameNum + i;
      string srcName{expandNum(opts.srcPath, frameNum)};

      srcClouds.emplace_back();
      auto& srcCloud = srcClouds.back();
      if (
        !ply::read(srcName, propNames, 1.0, srcCloud)
        || srcCloud.getPointCount() == 0) {
        throw runtime_error("failed to read input file: " + srcName);
      }
    }

    int totalPoints = 0;
    for (const auto& cloud : srcClouds)
      totalPoints += cloud.getPointCount();

    // merge sources, setting the frameIndex of each point in the merged
    // cloud to the group index of the corresponding source frame.
    PCCPointSet3 outCloud;
    outCloud.addFrameIndex();
    if (srcClouds[0].hasColors())
      outCloud.addColors();
    if (srcClouds[0].hasReflectances())
      outCloud.addReflectances();
    outCloud.resize(totalPoints);

    int outPtIdx = 0;
    for (int j = 0; j < srcClouds.size(); j++) {
      const auto& cloud = srcClouds[j];

      for (int ptIdx = 0; ptIdx < cloud.getPointCount(); ptIdx++, outPtIdx++) {
        outCloud[outPtIdx] = cloud[ptIdx];
        outCloud.setFrameIndex(outPtIdx, j);

        if (cloud.hasColors())
          outCloud.setColor(outPtIdx, cloud.getColor(ptIdx));

        if (cloud.hasReflectances())
          outCloud.setReflectance(outPtIdx, cloud.getReflectance(ptIdx));
      }
    }

    string outName{expandNum(opts.outPath, outFrameNum)};
    if (!ply::write(outCloud, propNames, 1, 0, outName, !opts.outputBinaryPly))
      throw runtime_error("failed to write output file: " + outName);
    cout << outName << endl;
  }
}

//---------------------------------------------------------------------------
// Read a merged frame, write out each component frame.
//

void
runSplit(const Options& opts)
{
  ply::PropertyNameMap propNames;
  propNames.position = {"x", "y", "z"};

  int outFrameNum = opts.firstOutputFrameNum;
  int srcFrameNum = opts.firstFrameNum;

  for (int i = 0; i < opts.frameCount; i++, srcFrameNum++, outFrameNum++) {
    string srcName{expandNum(opts.srcPath, srcFrameNum)};

    PCCPointSet3 srcCloud;
    if (
      !ply::read(srcName, propNames, 1.0, srcCloud)
      || srcCloud.getPointCount() == 0) {
      throw runtime_error("failed to read input file: " + srcName);
    }

    int numSrcPoints = srcCloud.getPointCount();

    if (!srcCloud.hasFrameIndex())
      throw runtime_error("missing frameindex property: " + srcName);

    // Extract each frame based on the frame index (assumed to start at 0)
    int frameIdx = 0;
    do {
      PCCPointSet3 outCloud;
      if (srcCloud.hasColors())
        outCloud.addColors();
      if (srcCloud.hasReflectances())
        outCloud.addReflectances();
      outCloud.resize(srcCloud.getPointCount());

      int nextFrameIdx = INT_MAX;
      int outPtIdx = 0;
      for (int ptIdx = 0; ptIdx < numSrcPoints; ptIdx++) {
        int ptFrameIdx = srcCloud.getFrameIndex(ptIdx);
        if (ptFrameIdx > frameIdx)
          nextFrameIdx = min(ptFrameIdx, nextFrameIdx);

        if (ptFrameIdx != frameIdx)
          continue;

        outCloud[outPtIdx] = srcCloud[ptIdx];

        if (srcCloud.hasColors())
          outCloud.setColor(outPtIdx, srcCloud.getColor(ptIdx));

        if (srcCloud.hasReflectances())
          outCloud.setReflectance(outPtIdx, srcCloud.getReflectance(ptIdx));

        outPtIdx++;
      }
      outCloud.resize(outPtIdx);

      string outName{expandNum(opts.outPath, outFrameNum)};
      if (outCloud.getPointCount() > 0)
        if (!ply::write(
              outCloud, propNames, 1, 0, outName, !opts.outputBinaryPly))
          throw runtime_error("failed to write output file: " + outName);
      cout << outName << endl;

      if (nextFrameIdx != INT_MAX)
        outFrameNum += nextFrameIdx - frameIdx;
      frameIdx = nextFrameIdx;
    } while (frameIdx != INT_MAX);
  }
}
