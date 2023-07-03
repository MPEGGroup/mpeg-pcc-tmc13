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

#include "RAHT.h"

#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>
#include <stdio.h>

#include "PCCTMC3Common.h"
#include "PCCMisc.h"

namespace pcc {

//============================================================================

void
PCCRAHTACCoefficientEntropyEstimate::init()
{
  for (int k = 0; k < 3; k++)
    probResGt0[k] = probResGt1[k] = (scaleRes >> 1);
  sumCostBits = 0.;
}

//---------------------------------------------------------------------------

void
PCCRAHTACCoefficientEntropyEstimate::updateCostBits(int32_t value, int k)
{
  int log2scaleRes = ilog2(uint32_t(scaleRes));
  double bits = 0;
  bits += value ? log2scaleRes - log2(probResGt0[k])
                : log2scaleRes - log2(scaleRes - probResGt0[k]);  //Gt0
  int mag = abs(value);
  if (mag) {
    bits += mag > 1 ? log2scaleRes - log2(probResGt1[k])
                    : log2scaleRes - log2(scaleRes - probResGt1[k]);  //Gt1
    bits += 1;  //sign bit.
    if (mag > 1)
      bits += 2.0 * log2(mag - 1.0) + 1.0;  //EG0 approximation.
  }
  sumCostBits += bits;
}

//----------------------------------------------------------------------------

void
PCCRAHTACCoefficientEntropyEstimate::resStatUpdate(int32_t value, int k)
{
  probResGt0[k] += value ? (scaleRes - probResGt0[k]) >> windowLog2
                         : -((probResGt0[k]) >> windowLog2);
  if (value)
    probResGt1[k] += abs(value) > 1 ? (scaleRes - probResGt1[k]) >> windowLog2
                                    : -((probResGt1[k]) >> windowLog2);
}

//============================================================================

struct UrahtNode {
  int64_t pos;
  int weight;
  Qps qp;

  uint8_t occupancy;
  std::vector<UrahtNode>::iterator firstChild;
  std::vector<UrahtNode>::iterator lastChild;
};

//============================================================================
// remove any non-unique leaves from a level in the uraht tree

int
reduceUnique(
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,
  std::vector<UrahtNode>* weightsOut,
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut,
  bool integer_haar_enable_flag)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInWrIt = weightsIn->begin();
  auto weightsInRdIt = weightsIn->cbegin();
  auto attrsInWrIt = attrsIn->begin();
  auto attrsInRdIt = attrsIn->begin();
  for (int i = 0; i < numNodes; i++) {
    const auto& node = *weightsInRdIt++;

    // copy across unique nodes
    if (node.pos != posPrev) {
      posPrev = node.pos;
      *weightsInWrIt++ = node;
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;
      continue;
    }

    // duplicate node
    (weightsInWrIt - 1)->weight += node.weight;
    weightsOut->push_back(node);
    for (int k = 0; k < numAttrs; k++) {
      if (integer_haar_enable_flag) {
        attrsOut->push_back(*attrsInRdIt++ - *(attrsInWrIt - numAttrs + k));
        *(attrsInWrIt - numAttrs + k) += attrsOut->back() >> 1;
      } else {
        *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
        attrsOut->push_back(*attrsInRdIt++);
      }
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Split a level of values into sum and difference pairs.

int
reduceLevel(
  int level,
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,
  std::vector<UrahtNode>* weightsOut,
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut,
  bool integer_haar_enable_flag
  )
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInWrIt = weightsIn->begin();
  auto weightsInRdIt = weightsIn->cbegin();
  auto attrsInWrIt = attrsIn->begin();
  auto attrsInRdIt = attrsIn->begin();

  for (int i = 0; i < numNodes; i++) {
    auto& node = *weightsInRdIt++;
    bool newPair = (posPrev ^ node.pos) >> level != 0;
    posPrev = node.pos;
    if (newPair) {
      *weightsInWrIt++ = node;
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;
    } else {
      auto& left = *(weightsInWrIt - 1);
      left.weight += node.weight;
      left.qp[0] = (left.qp[0] + node.qp[0]) >> 1;
      left.qp[1] = (left.qp[1] + node.qp[1]) >> 1;
      weightsOut->push_back(node);

      for (int k = 0; k < numAttrs; k++) {
        if (integer_haar_enable_flag) {
          attrsOut->push_back(*attrsInRdIt++ - *(attrsInWrIt - numAttrs + k));
          *(attrsInWrIt - numAttrs + k) += attrsOut->back() >> 1;
        } else {
          *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
          attrsOut->push_back(*attrsInRdIt++);
        }
      }
    }
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================
// Merge sum and difference values to form a tree.

void
expandLevel(
  int level,
  int numNodes,
  int numAttrs,
  std::vector<UrahtNode>* weightsIn,  // expand by numNodes before expand
  std::vector<UrahtNode>* weightsOut,  // shrink after expand
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut,
  bool integer_haar_enable_flag)
{
  if (numNodes == 0)
    return;

  // process a single level of the tree
  auto weightsInWrIt = weightsIn->rbegin();
  auto weightsInRdIt = std::next(weightsIn->crbegin(), numNodes);
  auto weightsOutRdIt = weightsOut->crbegin();
  auto attrsInWrIt = attrsIn->rbegin();
  auto attrsInRdIt = std::next(attrsIn->crbegin(), numNodes * numAttrs);
  auto attrsOutRdIt = attrsOut->crbegin();

  for (int i = 0; i < numNodes;) {
    bool isPair = (weightsOutRdIt->pos ^ weightsInRdIt->pos) >> level == 0;
    if (!isPair) {
      *weightsInWrIt++ = *weightsInRdIt++;
      for (int k = 0; k < numAttrs; k++)
        *attrsInWrIt++ = *attrsInRdIt++;
      
      continue;
    }

    // going to process a pair
    i++;

    // Out node is inserted before In node.
    const auto& nodeDelta = *weightsInWrIt++ = *weightsOutRdIt++;
    auto curAttrIt = attrsInWrIt;
    for (int k = 0; k < numAttrs; k++)
      *attrsInWrIt++ = *attrsOutRdIt++;

    // move In node to correct position, subtracting delta
    *weightsInWrIt = *weightsInRdIt++;
    (weightsInWrIt++)->weight -= nodeDelta.weight;
    for (int k = 0; k < numAttrs; k++) {
      *attrsInWrIt = *attrsInRdIt++;
      if (integer_haar_enable_flag) {
        *attrsInWrIt -= *curAttrIt >> 1;
        *curAttrIt++ += *attrsInWrIt++;
      } else {
        *attrsInWrIt++ -= *curAttrIt++;
      }
    }
  }
}

//============================================================================
// Search for neighbour with @value in the ordered list [first, last).
//
// If distance is positive, search [from, from+distance].
// If distance is negative, search [from-distance, from].

template<typename It, typename T, typename T2, typename Cmp>
It
findNeighbour(It first, It last, It from, T value, T2 distance, Cmp compare)
{
  It start = first;
  It end = last;

  if (distance >= 0) {
    start = from;
    if ((distance + 1) < std::distance(from, last))
      end = std::next(from, distance + 1);
  } else {
    end = from;
    if ((-distance) < std::distance(first, from))
      start = std::prev(from, -distance);
  }

  auto found = std::lower_bound(start, end, value, compare);
  if (found == end)
    return last;
  return found;
}

//============================================================================
// Find the neighbours of the node indicated by @t between @first and @last.
// The position weight of each found neighbour is stored in two arrays.

template<typename It>
void
findNeighbours(
  It first,
  It last,
  It it,
  It firstChild,
  It lastChild,
  int level,
  uint8_t occupancy,
  int parentNeighIdx[19],
  int childNeighIdx[12][8],
  const bool rahtSubnodePredictionEnabled,
  const int& raht_prediction_search_range)
{
  static const uint8_t neighMasks[19] = {255, 240, 204, 170, 192, 160, 136,
                                         3,   5,   15,  17,  51,  85,  10,
                                         34,  12,  68,  48,  80};

  // current position (discard extra precision)
  int64_t cur_pos = it->pos >> level;

  // the position of the parent, offset by (-1,-1,-1)
  int64_t base_pos = morton3dAdd(cur_pos, -1ll);

  // these neighbour offsets are relative to base_pos
  static const uint8_t neighOffset[19] = {0, 35, 21, 14, 49, 42, 28, 1,  2, 3,
                                          4, 5,  6,  10, 12, 17, 20, 33, 34};

  // special case for the direct parent (no need to search);
  parentNeighIdx[0] = std::distance(first, it);

  for (int i = 1; i < 19; i++) {
    // Only look for neighbours that have an effect
    if (!(occupancy & neighMasks[i])) {
      parentNeighIdx[i] = -1;
      continue;
    }

    // compute neighbour address to look for
    // the delta between it and the current position is
    int64_t neigh_pos = morton3dAdd(base_pos, neighOffset[i]);

    int64_t delta = neigh_pos - cur_pos;
    ///< in there will limit the prediction nearset neighbors searchRange 
    if (delta >= 0)
      delta =
      delta >= raht_prediction_search_range ? raht_prediction_search_range : delta;
    else
      delta =
      (-delta) >= raht_prediction_search_range ? -raht_prediction_search_range : delta;
    // find neighbour
    auto found = findNeighbour(
      first, last, it, neigh_pos, delta,
      [=](decltype(*it)& candidate, int64_t neigh_pos) {
        return (candidate.pos >> level) < neigh_pos;
      });

    if (found == last) {
      parentNeighIdx[i] = -1;
      continue;
    }

    if ((found->pos >> level) != neigh_pos) {
      parentNeighIdx[i] = -1;
      continue;
    }

    parentNeighIdx[i] = std::distance(first, found);
  }

  if (rahtSubnodePredictionEnabled) {
    //initialize the childNeighIdx
    for (int *p = (int*)childNeighIdx, i = 0; i < 96; p++, i++)
      *p = -1;

    static const uint8_t occuMasks[12] = {3,  5,  15, 17, 51, 85,
                                          10, 34, 12, 68, 48, 80};
    static const uint8_t occuShift[12] = {6, 5, 4, 3, 2, 1, 3, 1, 2, 1, 2, 3};

    int curLevel = level - 3;
    for (int i = 0; i < 9; i++) {
      if (parentNeighIdx[7 + i] == -1)
        continue;

      auto neiIt = first + parentNeighIdx[7 + i];
      uint8_t mask =
        (neiIt->occupancy >> occuShift[i]) & occupancy & occuMasks[i];
      if (!mask)
        continue;

      for (auto it = neiIt->firstChild; it != neiIt->lastChild; it++) {
        int nodeIdx = ((it->pos >> curLevel) & 0x7) - occuShift[i];
        if ((nodeIdx >= 0) && ((mask >> nodeIdx) & 1)) {
          childNeighIdx[i][nodeIdx] = std::distance(firstChild, it);
        }
      }
    }

    for (int i = 9; i < 12; i++) {
      if (parentNeighIdx[7 + i] == -1)
        continue;

      auto neiIt = first + parentNeighIdx[7 + i];
      uint8_t mask =
        (neiIt->occupancy << occuShift[i]) & occupancy & occuMasks[i];
      if (!mask)
        continue;

      for (auto it = neiIt->firstChild; it != neiIt->lastChild; it++) {
        int nodeIdx = ((it->pos >> curLevel) & 0x7) + occuShift[i];
        if ((nodeIdx < 8) && ((mask >> nodeIdx) & 1)) {
          childNeighIdx[i][nodeIdx] = std::distance(firstChild, it);
        }
      }
    }
  }
}

//============================================================================
// Generate the spatial prediction of a block.

template<bool isEncoder, bool rahtExtension, typename It>
void
intraDcPred(
  int numAttrs,
  const int parentNeighIdx[19],
  const int childNeighIdx[12][8],
  int occupancy,
  It first,
  It firstChild,
  It intraFirstChild,
  FixedPoint predBuf[][8],
  FixedPoint intraPredBuf[][8],
  const RahtPredictionParams &rahtPredParams, 
  int64_t& limitLow,
  int64_t& limitHigh,
  const bool& enableACInterPred)
{
  static const uint8_t predMasks[19] = {255, 240, 204, 170, 192, 160, 136,
                                        3,   5,   15,  17,  51,  85,  10,
                                        34,  12,  68,  48,  80};

  const auto& predWeightParent = rahtPredParams.predWeightParent;
  const auto& predWeightChild = rahtPredParams.predWeightChild;

  static const int kDivisors[64] = {
    32768, 16384, 10923, 8192, 6554, 5461, 4681, 4096, 3641, 3277, 2979,
    2731,  2521,  2341,  2185, 2048, 1928, 1820, 1725, 1638, 1560, 1489,
    1425,  1365,  1311,  1260, 1214, 1170, 1130, 1092, 1057, 1024, 993,
    964,   936,   910,   886,  862,  840,  819,  799,  780,  762,  745,
    728,   712,   697,   683,  669,  655,  643,  630,  618,  607,  596,
    585,   575,   565,   555,  546,  537,  529,  520,  512};

  int weightSum[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

  std::fill_n(&predBuf[0][0], 8 * numAttrs, FixedPoint(0));
  if (isEncoder && enableACInterPred)
    std::fill_n(&intraPredBuf[0][0], 8 * numAttrs, FixedPoint(0));

  int64_t neighValue[3];
  int64_t childNeighValue[3];
  int64_t intraChildNeighValue[3];

  const auto parentOnlyCheckMaxIdx =
    rahtPredParams.raht_subnode_prediction_enabled_flag ? 7 : 19;
  for (int i = 0; i < parentOnlyCheckMaxIdx; i++) {
    if (parentNeighIdx[i] == -1)
      continue;

    auto neighValueIt = std::next(first, numAttrs * parentNeighIdx[i]);
    for (int k = 0; k < numAttrs; k++)
      neighValue[k] = *neighValueIt++;

    // skip neighbours that are outside of threshold limits
    if (i) {
      if (10 * neighValue[0] <= limitLow || 10 * neighValue[0] >= limitHigh)
        continue;
    } else {
      constexpr int ratioThreshold1 = 2;
      constexpr int ratioThreshold2 = 25;
      limitLow = ratioThreshold1 * neighValue[0];
      limitHigh = ratioThreshold2 * neighValue[0];
    }

    // apply weighted neighbour value to masked positions
    for (int k = 0; k < numAttrs; k++)
      if (rahtExtension)
        neighValue[k] *= predWeightParent[i];
      else
        neighValue[k] *= predWeightParent[i] << pcc::FixedPoint::kFracBits;

    int mask = predMasks[i] & occupancy;
    for (int j = 0; mask; j++, mask >>= 1) {
      if (mask & 1) {
        weightSum[j] += predWeightParent[i];
        for (int k = 0; k < numAttrs; k++) {
          predBuf[k][j].val += neighValue[k];
          if (isEncoder && enableACInterPred)
            intraPredBuf[k][j].val += neighValue[k];
        }
      }
    }
  }
  if (rahtPredParams.raht_subnode_prediction_enabled_flag) {
    for (int i = 0; i < 12; i++) {
      if (parentNeighIdx[7 + i] == -1)
        continue;

      auto neighValueIt = std::next(first, numAttrs * parentNeighIdx[7 + i]);
      for (int k = 0; k < numAttrs; k++)
        neighValue[k] = *neighValueIt++;

      // skip neighbours that are outside of threshold limits
      if (10 * neighValue[0] <= limitLow || 10 * neighValue[0] >= limitHigh)
        continue;

      // apply weighted neighbour value to masked positions
      for (int k = 0; k < numAttrs; k++)
        if (rahtExtension)
          neighValue[k] *= predWeightParent[7 + i];
        else
          neighValue[k] *= predWeightParent[7 + i] << pcc::FixedPoint::kFracBits;

      int mask = predMasks[7 + i] & occupancy;
      for (int j = 0; mask; j++, mask >>= 1) {
        if (mask & 1) {
          if (childNeighIdx[i][j] != -1) {
            weightSum[j] += predWeightChild[i];
            auto childNeighValueIt =
              std::next(firstChild, numAttrs * childNeighIdx[i][j]);
            for (int k = 0; k < numAttrs; k++)
              if (rahtExtension)
                childNeighValue[k] = (*childNeighValueIt++)
                  * predWeightChild[i];
              else
                childNeighValue[k] = (*childNeighValueIt++)
                  * (predWeightChild[i] << pcc::FixedPoint::kFracBits);

            for (int k = 0; k < numAttrs; k++)
              predBuf[k][j].val += childNeighValue[k];

            if (isEncoder && enableACInterPred) {
              auto intraChildNeighValueIt =
                std::next(intraFirstChild, numAttrs * childNeighIdx[i][j]);
              for (int k = 0; k < numAttrs; k++)
                if (rahtExtension)
                  intraChildNeighValue[k] =
                    (*intraChildNeighValueIt++) * predWeightChild[i];
                else
                  intraChildNeighValue[k] = (*intraChildNeighValueIt++)
                    * (predWeightChild[i] << pcc::FixedPoint::kFracBits);
              for (int k = 0; k < numAttrs; k++)
                intraPredBuf[k][j].val += intraChildNeighValue[k];
            }
          } else {
            weightSum[j] += predWeightParent[7 + i];
            for (int k = 0; k < numAttrs; k++) {
              predBuf[k][j].val += neighValue[k];
              if (isEncoder && enableACInterPred)
                intraPredBuf[k][j].val += neighValue[k];
            }
          }
        }
      }
    }
  }

  // normalise
  FixedPoint div;
  for (int i = 0; i < 8; i++, occupancy >>= 1) {
    if (occupancy & 1) {
      div.val = kDivisors[weightSum[i]];
      for (int k = 0; k < numAttrs; k++) {
        predBuf[k][i] *= div;
        if (isEncoder && enableACInterPred)
          intraPredBuf[k][i] *= div;
      }
      if (rahtPredParams.integer_haar_enable_flag) {
        for (int k = 0; k < numAttrs; k++) {
          predBuf[k][i].val = (predBuf[k][i].val >> predBuf[k][i].kFracBits)
            << predBuf[k][i].kFracBits;
          if (isEncoder && enableACInterPred)
            intraPredBuf[k][i].val =
              (intraPredBuf[k][i].val >> intraPredBuf[k][i].kFracBits)
              << intraPredBuf[k][i].kFracBits;
        }
      }
    }
  }
}

//============================================================================
// Encapsulation of a RAHT transform stage.

class RahtKernel {
public:
  RahtKernel(int weightLeft, int weightRight)
  {
    uint64_t w = weightLeft + weightRight;
    uint64_t isqrtW = irsqrt(w);
    _a.val =
      (isqrt(uint64_t(weightLeft) << (2 * _a.kFracBits)) * isqrtW) >> 40;
    _b.val =
      (isqrt(uint64_t(weightRight) << (2 * _b.kFracBits)) * isqrtW) >> 40;
  }

  void fwdTransform(
    FixedPoint left, FixedPoint right, FixedPoint* lf, FixedPoint* hf)
  {
    FixedPoint a = _a, b = _b;
    // lf = left * a + right * b
    // hf = right * a - left * b

    *lf = right;
    *lf *= b;
    *hf = right;
    *hf *= a;

    a *= left;
    b *= left;

    *lf += a;
    *hf -= b;
  }

  void invTransform(
    FixedPoint lf, FixedPoint hf, FixedPoint* left, FixedPoint* right)
  {
    FixedPoint a = _a, b = _b;

    *left = lf;
    *left *= a;
    *right = lf;
    *right *= b;

    b *= hf;
    a *= hf;

    *left -= b;
    *right += a;
  }

private:
  FixedPoint _a, _b;
};

//============================================================================
// Encapsulation of an Integer Haar transform stage.

class HaarKernel {
public:
  HaarKernel(int weightLeft, int weightRight) {}

  void fwdTransform(
    FixedPoint left, FixedPoint right, FixedPoint* lf, FixedPoint* hf)
  {
    hf->val = right.val - left.val;
    lf->val = left.val + ((hf->val >> (1 + hf->kFracBits)) << hf->kFracBits);
  }

  void invTransform(
    FixedPoint lf, FixedPoint hf, FixedPoint* left, FixedPoint* right)
  {
    left->val = lf.val - (((hf.val >> (1 + hf.kFracBits)) << hf.kFracBits));
    right->val = hf.val + left->val;
  }
};

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

template<class Kernel>
void
fwdTransformBlock222(
  int numBufs, FixedPoint buf[][8], int weights[8 + 8 + 8 + 8])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 0, iw = 0; i < 12; i++, iw += 2) {
    int i0 = a[i];
    int i1 = b[i];

    if (weights[iw] + weights[iw + 1] == 0)
      continue;

    // only one occupied, propagate to next level
    if (!weights[iw] || !weights[iw + 1]) {
      if (!weights[iw]) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

    // actual transform
    Kernel kernel(weights[iw], weights[iw + 1]);
    for (int k = 0; k < numBufs; k++) {
      auto& bufk = buf[k];
      kernel.fwdTransform(bufk[i0], bufk[i1], &bufk[i0], &bufk[i1]);
    }
  }
}

//============================================================================
// In-place inverse transform a set of sparse 2x2x2 blocks each using the
// same weights

template<class Kernel>
void
invTransformBlock222(
  int numBufs, FixedPoint buf[][8], int weights[8 + 8 + 8 + 8])
{
  static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
  static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
  for (int i = 11, iw = 22; i >= 0; i--, iw -= 2) {
    int i0 = a[i];
    int i1 = b[i];

    if (weights[iw] + weights[iw + 1] == 0)
      continue;

    // only one occupied, propagate to next level
    if (!weights[iw] || !weights[iw + 1]) {
      if (!weights[iw]) {
        for (int k = 0; k < numBufs; k++)
          std::swap(buf[k][i0], buf[k][i1]);
      }
      continue;
    }

    // actual transform
    Kernel kernel(weights[iw], weights[iw + 1]);
    for (int k = 0; k < numBufs; k++) {
      auto& bufk = buf[k];
      kernel.invTransform(bufk[i0], bufk[i1], &bufk[i0], &bufk[i1]);
    }
  }
}

//============================================================================
// expand a set of eight weights into three levels

void
mkWeightTree(int weights[8 + 8 + 8 + 8])
{
  int* in = &weights[0];
  int* out = &weights[8];

  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
  out += 4;
  for (int i = 0; i < 4; i++) {
    out[0] = out[4] = in[0] + in[1];
    if (!in[0] || !in[1])
      out[4] = 0;  // single node, no high frequencies
    in += 2;
    out++;
  }
}

//============================================================================
// Invoke mapFn(coefIdx) for each present coefficient in the transform

template<class T>
void
scanBlock(int weights[8 + 8 + 8 + 8], T mapFn)
{
  static const int8_t kRahtScanOrder[] = {0, 4, 2, 1, 6, 5, 3, 7};

  // there is always the DC coefficient (empty blocks are not transformed)
  mapFn(0);

  for (int i = 1; i < 8; i++) {
    if (!weights[24 + kRahtScanOrder[i]])
      continue;

    mapFn(kRahtScanOrder[i]);
  }
}

//============================================================================
// Tests if two positions are siblings at the given tree level

static bool
isSibling(int64_t pos0, int64_t pos1, int level)
{
  return ((pos0 ^ pos1) >> level) == 0;
}

//============================================================================
// estimate filter tap by a binary search

int getFilterTap (int64_t autocorr, int64_t crosscorr)
{
  //binary search to replace divison. ( returns 128*crosscorr/autocorr)
  if (crosscorr == 0)
    return 0;
  
  bool isneg = crosscorr < 0;
  crosscorr = abs(crosscorr);
  
  if (crosscorr == autocorr)
    return (isneg ? -128 : 128);
  
  int tapint = 0, tapfrac = 0;
  
  //determine integer part by repeated subtraction
  while (crosscorr >= autocorr) {
    crosscorr -= autocorr;
    tapint += 128;
  }
  if (crosscorr == 0) {
    return (isneg ? -tapint : tapint);
  }
  
  int min = 0, max = 128;
  
  while (min < (max - 1)) {
    int mid = (min+max) >> 1;
    int64_t midval = (mid*autocorr)>>7;
    if (crosscorr == midval) {
      tapfrac = mid; 
      return (isneg ? -(tapint + tapfrac) : (tapint + tapfrac));
    }
    else if (crosscorr < midval) {
      max = mid;
    }
    else {
      min = mid;
    }
  }
  tapfrac = min;
  return  (isneg ? -(tapint+tapfrac) : (tapint+tapfrac));
  
}

int
estimate_layer_filter(
  const std::vector<UrahtNode>& weightsLf,
  const std::vector<UrahtNode>& weightsLf_ref,
  const std::vector<int>& attrsLf,
  const std::vector<int>& attrsLf_ref,
  int level, int level_ref, int numAttrs, bool inheritDc,  bool rahtExtension
){
  int64_t autocorr = 0, crosscorr = 0;
  int layerFilter = 128;
  for (int i = 0, j = 0, iLast, jLast, iEnd = weightsLf.size(), jEnd = weightsLf_ref.size(); i < iEnd; i = iLast) {
    FixedPoint transformBuf[6][8] = {};
    FixedPoint transformInterPredBuf[3][8] = {};
    int weights[8 + 8 + 8 + 8] = {};
    Qps nodeQp[8] = {};
    uint8_t occupancy = 0;
    int nodeCnt = 0;
    
    int weights_ref[8 + 8 + 8 + 8] = {};
    bool interNode = false;
    
    const auto cur_pos = weightsLf[i].pos >> (level + 3);
    auto ref_pos = j < jEnd-1 ? weightsLf_ref[j].pos >> (level_ref + 3) : 0x7FFFFFFFFFFFFFFFLL;
    while ((j < jEnd-1) && (cur_pos > ref_pos)) {
      j++;
      ref_pos = weightsLf_ref[j].pos >> (level_ref + 3);
    }
    if (cur_pos == ref_pos) {
      interNode = true;
    }
    
    if (interNode) {
      for (jLast = j; jLast < jEnd; jLast++) {
        int nextNode = jLast > j
        && !isSibling(weightsLf_ref[jLast].pos, weightsLf_ref[j].pos, level_ref + 3);
        if (nextNode)
          break;
        int nodeIdx = (weightsLf_ref[jLast].pos >> level_ref) & 0x7;
        weights_ref[nodeIdx] = weightsLf_ref[jLast].weight;
        for (int k = 0; k < numAttrs; k++)
          transformInterPredBuf[k][nodeIdx] = attrsLf_ref[jLast * numAttrs + k];
      }
    }
    
    for (iLast = i; iLast < iEnd; iLast++) {
      int nextNode = iLast > i
      && !isSibling(weightsLf[iLast].pos, weightsLf[i].pos, level + 3);
      if (nextNode)
        break;
      
      int nodeIdx = (weightsLf[iLast].pos >> level) & 0x7;
      weights[nodeIdx] = weightsLf[iLast].weight;
      
      occupancy |= 1 << nodeIdx;
      
      if (rahtExtension)
        nodeCnt++;
  
      for (int k = 0; k < numAttrs; k++)
        transformBuf[k][nodeIdx] = attrsLf[iLast * numAttrs + k];
    }
    
    mkWeightTree(weights);
    mkWeightTree(weights_ref);
    
    if (rahtExtension && nodeCnt == 1) {
      interNode = false;
    }
    
    if (interNode) {
      for (int childIdx = 0; childIdx < 8; childIdx++) {
        if (weights_ref[childIdx] <= 1)
          continue;
        if(weights_ref[childIdx]) {
          FixedPoint rsqrtWeight;
          uint64_t w = weights_ref[childIdx];
          int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformInterPredBuf[k][childIdx].val >>= shift;
            transformInterPredBuf[k][childIdx] *= rsqrtWeight;
          }
        }
      }
    }
    
    for (int childIdx = 0; childIdx < 8; childIdx++) {
      if (weights[childIdx] <= 1)
        continue;
      
      // Summed attribute values
      if (true) {
        FixedPoint rsqrtWeight;
        uint64_t w = weights[childIdx];
        int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
        rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
        for (int k = 0; k < numAttrs; k++) {
          transformBuf[k][childIdx].val >>= shift;
          transformBuf[k][childIdx] *= rsqrtWeight;
        }
      }
    }
    
    if (interNode) {
      fwdTransformBlock222<RahtKernel>(numAttrs, transformBuf, weights);
      fwdTransformBlock222<RahtKernel>(numAttrs, transformInterPredBuf, weights_ref);
      scanBlock(weights, [&](int idx) {
        if (inheritDc && !idx)
          return;
        
        int shiftbits = transformBuf[0][idx].kFracBits;
        int64_t refVal = transformInterPredBuf[0][idx].val;
        if (refVal) {
          autocorr += (refVal*refVal) >> shiftbits;
          crosscorr += (refVal*transformBuf[0][idx].val) >> shiftbits;
        }
      });
    }
  }
  if (autocorr) {
    layerFilter = getFilterTap(autocorr, crosscorr);
  }
  return layerFilter;
}

//============================================================================
// Core transform process (for encoder/decoder)

template<bool isEncoder, bool rahtExtension>
void
uraht_process(
  const RahtPredictionParams& rahtPredParams,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int numPoints,
  int numAttrs,
  int64_t* positions,
  int* attributes,
  int32_t* coeffBufIt,
  AttributeInterPredParams& attrInterPredParams)
{
  // coefficients are stored in three planar arrays.  coeffBufItK is a set
  // of iterators to each array.
  int32_t* coeffBufItK[3] = {
    coeffBufIt,
    coeffBufIt + numPoints,
    coeffBufIt + numPoints * 2,
  };

  if (numPoints == 1) {
    auto quantizers = qpset.quantizers(0, pointQpOffsets[0]);
    for (int k = 0; k < numAttrs; k++) {
      auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

      if (isEncoder) {
        auto coeff = attributes[k];
        assert(coeff <= INT_MAX && coeff >= INT_MIN);
        *coeffBufItK[k]++ = coeff =
          q.quantize(coeff << kFixedPointAttributeShift);
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
      } else {
        int64_t coeff = *coeffBufItK[k]++;
        attributes[k] =
          divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
      }
    }
    return;
  }

  std::vector<UrahtNode> weightsLf, weightsHf;
  std::vector<int> attrsLf, attrsHf;

  std::vector<UrahtNode> weightsLf_ref, weightsHf_ref;
  std::vector<int> attrsLf_ref, attrsHf_ref;

  bool enableACInterPred = attrInterPredParams.enableAttrInterPred;
  const bool enableDCInterPred = attrInterPredParams.enableAttrInterPred;
  bool enableFilterEstimation = attrInterPredParams.paramsForInterRAHT.enableFilterEstimation;
  bool enableACRDOInterPred =
    attrInterPredParams.paramsForInterRAHT.raht_enable_inter_intra_layer_RDO
    && enableACInterPred;

  int refIdx = 0;
  int treeDepth = 0;
  int treeDepthLimit =
    attrInterPredParams.paramsForInterRAHT.raht_inter_prediction_depth_minus1
    + 1;

  std::vector<int64_t> fixedFilterTaps = {128, 128, 128, 127, 125, 121, 115};
  int skipInitLayersForFiltering = attrInterPredParams.paramsForInterRAHT.skipInitLayersForFiltering;

  weightsLf.reserve(numPoints);
  attrsLf.reserve(numPoints * numAttrs);


  int regionQpShift = 4;

  const int maxAcCoeffQpOffsetLayers = qpset.rahtAcCoeffQps.size() - 1;

  // copy positions into internal form
  // todo(df): lift to api
  for (int i = 0; i < numPoints; i++) {
    weightsLf.emplace_back(UrahtNode{
      positions[i],
      1,
      {pointQpOffsets[i][0] << regionQpShift,
       pointQpOffsets[i][1] << regionQpShift}});
    for (int k = 0; k < numAttrs; k++) {
      attrsLf.push_back(attributes[i * numAttrs + k]);
    }
  }

  weightsHf.reserve(numPoints);
  attrsHf.reserve(numPoints * numAttrs);

  if (enableACInterPred) {
    weightsLf_ref.reserve(attrInterPredParams.paramsForInterRAHT.voxelCount);
    attrsLf_ref.reserve(attrInterPredParams.paramsForInterRAHT.voxelCount * numAttrs);

    for (int i = 0; i < attrInterPredParams.paramsForInterRAHT.voxelCount; i++) {
      weightsLf_ref.emplace_back(UrahtNode{attrInterPredParams.paramsForInterRAHT.mortonCode[i],
         1, {0, 0}});
      for (int k = 0; k < numAttrs; k++) {
        attrsLf_ref.push_back(attrInterPredParams.paramsForInterRAHT.attributes[i * numAttrs + k]);
      }
    }

    weightsHf_ref.reserve(attrInterPredParams.paramsForInterRAHT.voxelCount);
    attrsHf_ref.reserve(attrInterPredParams.paramsForInterRAHT.voxelCount * numAttrs);
  }
  
  // ascend tree
  std::vector<int> levelHfPos;
  std::vector<int> levelHfPos_ref;

  int numDupNodes = numPoints;
  for (int level = 0, numNodes = weightsLf.size(); numNodes > 1; level++) {
    levelHfPos.push_back(weightsHf.size());
    if (level == 0) {
      // process any duplicate points
      numNodes = reduceUnique(
        numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf,
        rahtPredParams.integer_haar_enable_flag);
      numDupNodes -= numNodes;
    } else {
      // normal level reduction
      numNodes = reduceLevel(
        level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf,
        rahtPredParams.integer_haar_enable_flag);
    }
  }
  

  if (enableACInterPred){
    for (int level = 0, numNodes = weightsLf_ref.size(); numNodes > 1; level++) {
      levelHfPos_ref.push_back(weightsHf_ref.size());
      if (level == 0) {
        // process any duplicate points
        numNodes = reduceUnique(
          numNodes, numAttrs, &weightsLf_ref, &weightsHf_ref,
          &attrsLf_ref, &attrsHf_ref, rahtPredParams.integer_haar_enable_flag);
      } else {
        // normal level reduction
        numNodes = reduceLevel(
          level, numNodes, numAttrs, &weightsLf_ref, &weightsHf_ref,
          &attrsLf_ref, &attrsHf_ref, rahtPredParams.integer_haar_enable_flag);
      }
    }    
  }

  assert(weightsLf[0].weight == numPoints);

  // reconstruction buffers
  std::vector<int64_t> attrRec, attrRecParent, intraAttrRec;
  attrRec.resize(numPoints * numAttrs);
  attrRecParent.resize(numPoints * numAttrs);
  if (isEncoder && enableACRDOInterPred)
    intraAttrRec.resize(numPoints * numAttrs);

  std::vector<int64_t> attrRecUs, attrRecParentUs, intraAttrRecUs;
  attrRecUs.resize(numPoints * numAttrs);
  attrRecParentUs.resize(numPoints * numAttrs);
  if (isEncoder && enableACRDOInterPred)
    intraAttrRecUs.resize(numPoints * numAttrs);

  std::vector<UrahtNode> weightsParent;
  weightsParent.reserve(numPoints);

  std::vector<int> numParentNeigh, numGrandParentNeigh;
  numParentNeigh.resize(numPoints);
  numGrandParentNeigh.resize(numPoints);

  // quant layer selection
  auto qpLayer = 0;
  // AC coeff QP offset laer
  auto acCoeffQpLayer = -1;

  // descend tree
  weightsLf.resize(1);
  attrsLf.resize(numAttrs);

  weightsLf_ref.resize(1);
  attrsLf_ref.resize(numAttrs);

  int trainZeros = 0;
  int sumNodes = 0;
  int intraTrainZeros = 0;
  PCCRAHTACCoefficientEntropyEstimate intraEstimate;
  PCCRAHTACCoefficientEntropyEstimate curEstimate;
  int depth = 0;
  std::vector<int> intraACCoeffcients;
  if (isEncoder && enableACRDOInterPred) {
    intraACCoeffcients.resize(numPoints * numAttrs);
  }
  std::vector<NodeInfoRAHT> intraNodes;
  for (int level = levelHfPos.size() - 1, level_ref = levelHfPos_ref.size() - 1, isFirst = 1; level > 0; /*nop*/) {
    int numNodes = weightsHf.size() - levelHfPos[level];
    sumNodes += numNodes;
    weightsLf.resize(weightsLf.size() + numNodes);
    attrsLf.resize(attrsLf.size() + numNodes * numAttrs);
    expandLevel(
      level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf,
      rahtPredParams.integer_haar_enable_flag);

    weightsHf.resize(levelHfPos[level]);
    attrsHf.resize(levelHfPos[level] * numAttrs);


    if (level_ref <= 0) {
      enableACInterPred = false;
    }

    if (treeDepth >= treeDepthLimit)
      enableACInterPred = false;

    if (enableACInterPred) {
      int numNodes_ref = weightsHf_ref.size() - levelHfPos_ref[level_ref];
      weightsLf_ref.resize(weightsLf_ref.size() + numNodes_ref);
      attrsLf_ref.resize(attrsLf_ref.size() + numNodes_ref * numAttrs);
      expandLevel(
        level_ref, numNodes_ref, numAttrs, &weightsLf_ref, &weightsHf_ref,
        &attrsLf_ref, &attrsHf_ref, rahtPredParams.integer_haar_enable_flag);
      weightsHf_ref.resize(levelHfPos_ref[level_ref]);
      attrsHf_ref.resize(levelHfPos_ref[level_ref] * numAttrs);
    }

    enableACRDOInterPred =
      attrInterPredParams.paramsForInterRAHT.raht_enable_inter_intra_layer_RDO
      && enableACInterPred;

    // expansion of level is complete, processing is now on the next level
    level--;
    level_ref--;

    // every three levels, perform transform
    if (level % 3)
      continue;
    ///< current level nodes number is equal to previous nodes level
    if (sumNodes == 0)
      continue;

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;
    bool enablePredictionInLvl =
      inheritDc && rahtPredParams.raht_prediction_enabled_flag;
    isFirst = 0;

    bool predictionInter = false;
    int64_t position = 0;

    if (enablePredictionInLvl) {
      for (auto& ele : weightsParent)
        ele.occupancy = 0;

      const int parentCount = weightsParent.size();
      auto it = weightsLf.begin();
      for (auto i = 0; i < parentCount; i++) {
        weightsParent[i].firstChild = it++;

        while (it != weightsLf.end()
               && !((it->pos ^ weightsParent[i].pos) >> (level + 3)))
          it++;
        weightsParent[i].lastChild = it;
      }
    }

    int32_t* intraCoeffBufItK[3] = {
      intraACCoeffcients.data(),
      intraACCoeffcients.data() + sumNodes,
      intraACCoeffcients.data() + sumNodes * 2,
    };
    int32_t* intraCoeffBufItBeginK[3] = {
      intraCoeffBufItK[0],
      intraCoeffBufItK[1],
      intraCoeffBufItK[2],
    };

    int32_t* coeffBufItBeginK[3] = {
      coeffBufItK[0],
      coeffBufItK[1],
      coeffBufItK[2],
    };
    bool curLevelEnableACInterPred = false;
    if (isEncoder)
      curLevelEnableACInterPred =
        enablePredictionInLvl && enableACRDOInterPred;
    else {
      curLevelEnableACInterPred = enablePredictionInLvl && enableACRDOInterPred
        && attrInterPredParams.attr_layer_code_mode[depth];
    }

    // select quantiser according to transform layer
    qpLayer = std::min(qpLayer + 1, int(qpset.layers.size()) - 1);
    acCoeffQpLayer++;
    
    int64_t interFilterTap = 128;
    if ((!enableFilterEstimation) && (enableACInterPred) && (treeDepth < treeDepthLimit)) {
      int filtexidx = treeDepth < fixedFilterTaps.size() ? treeDepth : (fixedFilterTaps.size()-1);
      interFilterTap = fixedFilterTaps[filtexidx];
    }

    // prepare reconstruction buffers
    //  previous reconstruction -> attrRecParent
    std::swap(attrRec, attrRecParent);
    std::swap(attrRecUs, attrRecParentUs);
    std::swap(numParentNeigh, numGrandParentNeigh);
    auto attrRecParentUsIt = attrRecParentUs.cbegin();
    auto attrRecParentIt = attrRecParent.cbegin();
    auto weightsParentIt = weightsParent.begin();
    auto numGrandParentNeighIt = numGrandParentNeigh.cbegin();
    
    bool enableEstimateLayer = enableFilterEstimation && enableACInterPred && (treeDepth<treeDepthLimit) && (treeDepth >= skipInitLayersForFiltering ) ;
    //begin filter estimation at encoder
    if (isEncoder && enableEstimateLayer ) {
      
      int origFilterTap =  estimate_layer_filter (weightsLf, weightsLf_ref,  attrsLf,attrsLf_ref, level, level_ref, numAttrs, inheritDc,  rahtExtension  );
      int residueFilterTap = 128 - origFilterTap;
      auto quantizers = qpset.quantizers(qpLayer, {0,0});
      auto& q = quantizers[0];
      int64_t quantizedResFilterTap = q.quantize(residueFilterTap << kFixedPointAttributeShift);
      int64_t recResidueFilterTap = divExp2RoundHalfUp(q.scale(quantizedResFilterTap), kFixedPointAttributeShift);
      attrInterPredParams.paramsForInterRAHT.FilterTaps.push_back(quantizedResFilterTap);
      interFilterTap = 128-recResidueFilterTap;
    } //end filter estimation
    
    //get filter tap at the decoder 
    bool enabledecoderparsing = !isEncoder && (enableFilterEstimation) && (treeDepth < (attrInterPredParams.paramsForInterRAHT.FilterTaps.size() + skipInitLayersForFiltering))  && (treeDepth >= skipInitLayersForFiltering);
    if (enabledecoderparsing) {
      auto quantizers = qpset.quantizers(qpLayer, {0,0});
      auto& q = quantizers[0];
      int idxtoread= treeDepth - skipInitLayersForFiltering;
      int64_t recResidueFilterTap = divExp2RoundHalfUp(q.scale(attrInterPredParams.paramsForInterRAHT.FilterTaps[idxtoread]), kFixedPointAttributeShift);
      interFilterTap = 128 - recResidueFilterTap;
    }
    for (int i = 0, j = 0, iLast, jLast, iEnd = weightsLf.size(), jEnd = weightsLf_ref.size(); i < iEnd; i = iLast) {
      // todo(df): hoist and dynamically allocate
      FixedPoint transformBuf[6][8] = {};
      FixedPoint(*transformPredBuf)[8] = &transformBuf[numAttrs];
      FixedPoint transformInterPredBuf[3][8] = {};
      FixedPoint transformIntraBuf[3][8] = {};
      FixedPoint transformIntraPredBuf[3][8] = {};
      int weights[8 + 8 + 8 + 8] = {};
      Qps nodeQp[8] = {};
      uint8_t occupancy = 0;

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
      int nodeCnt = 0;

      int weights_ref[8 + 8 + 8 + 8] = {};
      bool interNode = false;
      if (curLevelEnableACInterPred
          || (enableACInterPred && !enablePredictionInLvl)) {
        const auto cur_pos = weightsLf[i].pos >> (level + 3);
        auto ref_pos = weightsLf_ref[j].pos >> (level_ref + 3);
        while ((j < weightsLf_ref.size() - 1) && (cur_pos > ref_pos)) {
          j++;
          ref_pos = weightsLf_ref[j].pos >> (level_ref + 3);
        }
        if (cur_pos == ref_pos) {
          interNode = true;
        }
      }

      if (interNode) {
        for (jLast = j; jLast < jEnd; jLast++) {
          int nextNode = jLast > j
            && !isSibling(weightsLf_ref[jLast].pos, weightsLf_ref[j].pos, level_ref + 3);
          if (nextNode)
            break;
          int nodeIdx = (weightsLf_ref[jLast].pos >> level_ref) & 0x7;
          weights_ref[nodeIdx] = weightsLf_ref[jLast].weight;
          for (int k = 0; k < numAttrs; k++)
            transformInterPredBuf[k][nodeIdx] = attrsLf_ref[jLast * numAttrs + k];
        }
      }

      for (iLast = i; iLast < iEnd; iLast++) {
        int nextNode = iLast > i
          && !isSibling(weightsLf[iLast].pos, weightsLf[i].pos, level + 3);
        if (nextNode)
          break;

        int nodeIdx = (weightsLf[iLast].pos >> level) & 0x7;
        weights[nodeIdx] = weightsLf[iLast].weight;
        nodeQp[nodeIdx][0] = weightsLf[iLast].qp[0] >> regionQpShift;
        nodeQp[nodeIdx][1] = weightsLf[iLast].qp[1] >> regionQpShift;

        occupancy |= 1 << nodeIdx;

        if (rahtExtension)
          nodeCnt++;

        if (isEncoder) {
          for (int k = 0; k < numAttrs; k++)
            transformBuf[k][nodeIdx] = attrsLf[iLast * numAttrs + k];
        }
      }

      mkWeightTree(weights);

      mkWeightTree(weights_ref);   

      if (!inheritDc) {
        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNeigh[j++] = 19;
        }
      }
      if (rahtExtension && nodeCnt == 1){
        interNode = false;
      }

      // Inter-level prediction:
      //  - Find the parent neighbours of the current node
      //  - Generate prediction for all attributes into transformPredBuf
      //  - Subtract transformed coefficients from forward transform
      //  - The transformPredBuf is then used for reconstruction
      bool enablePrediction = enablePredictionInLvl;
      if (enablePredictionInLvl) {
        weightsParentIt->occupancy = occupancy;
        // indexes of the neighbouring parents
        int parentNeighIdx[19];
        int childNeighIdx[12][8];

        int parentNeighCount = 0;
        if (rahtExtension && nodeCnt == 1) {
          enablePrediction = false;
          parentNeighCount = 19;
        } else if (
          *numGrandParentNeighIt < rahtPredParams.raht_prediction_threshold0) {
          enablePrediction = false;
        } else {
          findNeighbours(
            weightsParent.begin(), weightsParent.end(), weightsParentIt,
            weightsLf.begin(), weightsLf.begin() + i, level + 3, occupancy,
            parentNeighIdx, childNeighIdx,
            rahtPredParams.raht_subnode_prediction_enabled_flag, rahtPredParams.raht_prediction_search_range);
          for (int i = 0; i < 19; i++) {
            parentNeighCount += (parentNeighIdx[i] != -1);
          }
          if (parentNeighCount < rahtPredParams.raht_prediction_threshold1) {
            enablePrediction = false;
          } else {
            int64_t limitLow = 0;
            int64_t limitHigh = 0;
            intraDcPred<isEncoder, rahtExtension>(
              numAttrs, parentNeighIdx, childNeighIdx, occupancy,
              attrRecParent.begin(), attrRec.begin(), intraAttrRec.begin(),
              transformPredBuf, transformIntraPredBuf, rahtPredParams,
              limitLow, limitHigh, curLevelEnableACInterPred);
          }
        }

        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNeigh[j++] = parentNeighCount;
        }
      }

      int parentWeight = 0;
      if (inheritDc) {
        parentWeight = weightsParentIt->weight;
        weightsParentIt++;
        numGrandParentNeighIt++;
      }

      bool enableIntraPrediction =
        curLevelEnableACInterPred && enablePrediction;
      bool enableInterPrediction = curLevelEnableACInterPred;

      if (!rahtPredParams.integer_haar_enable_flag) {
      // normalise coefficients
      if (interNode){
        for (int childIdx = 0; childIdx < 8; childIdx++) {
          if (weights_ref[childIdx] <= 1)
            continue;
          if(weights_ref[childIdx]){
            FixedPoint rsqrtWeight;
            uint64_t w = weights_ref[childIdx];
            int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
            rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits); 
            for (int k = 0; k < numAttrs; k++) {
              transformInterPredBuf[k][childIdx].val >>= shift;
              transformInterPredBuf[k][childIdx] *= rsqrtWeight;
            }
          }
        }
        if (!isEncoder)
          enablePrediction = false;
      }


      for (int childIdx = 0; childIdx < 8; childIdx++) {

        if (weights[childIdx] <= 1)
          continue;

        // Summed attribute values
        if (isEncoder) {
          FixedPoint rsqrtWeight;
          uint64_t w = weights[childIdx];
          int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][childIdx].val >>= shift;
            transformBuf[k][childIdx] *= rsqrtWeight;
          }
        }

        // Predicted attribute values
        FixedPoint sqrtWeight;
        if (enablePrediction) {
          sqrtWeight.val =
            isqrt(uint64_t(weights[childIdx]) << (2 * FixedPoint::kFracBits));
          for (int k = 0; k < numAttrs; k++) {
            transformPredBuf[k][childIdx] *= sqrtWeight;
          }
        }
        if (isEncoder && enableIntraPrediction) {
          for (int k = 0; k < numAttrs; k++) {
            transformIntraPredBuf[k][childIdx] *= sqrtWeight;
          }
        }
      }
      }

      // forward transform:
      //  - encoder: transform both attribute sums and prediction
      //  - decoder: just transform prediction
      if (rahtPredParams.integer_haar_enable_flag) {
        if (isEncoder && enablePrediction)
          fwdTransformBlock222<HaarKernel>(2 * numAttrs, transformBuf, weights);
        else if (isEncoder)
          fwdTransformBlock222<HaarKernel>(numAttrs, transformBuf, weights);
        else if (enablePrediction)
          fwdTransformBlock222<HaarKernel>(numAttrs, transformPredBuf, weights);

        if (interNode) {
          fwdTransformBlock222<HaarKernel>(numAttrs, transformInterPredBuf, weights_ref);
          for (int childIdx = 0; childIdx < 8; childIdx++) {
            for (int k = 0; k < numAttrs; k++){
              int64_t refVal = transformInterPredBuf[k][childIdx].val;
              // NB: integer haar is not compatible with inter filter.
              transformPredBuf[k][childIdx].val = refVal;
            }
          }
          enablePrediction = true;
        }

        if (isEncoder && enableIntraPrediction)
          fwdTransformBlock222<HaarKernel>(numAttrs, transformIntraPredBuf, weights);
      }
      else {
        if (isEncoder && enablePrediction)
          fwdTransformBlock222<RahtKernel>(2 * numAttrs, transformBuf, weights);
        else if (isEncoder)
          fwdTransformBlock222<RahtKernel>(numAttrs, transformBuf, weights);
        else if (enablePrediction)
          fwdTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);

        if (interNode) {
          fwdTransformBlock222<RahtKernel>(numAttrs, transformInterPredBuf, weights_ref);
          for (int childIdx = 0; childIdx < 8; childIdx++) {
            for (int k = 0; k < numAttrs; k++){
              int64_t refVal = transformInterPredBuf[k][childIdx].val;
              int64_t filteredVal = (treeDepth<skipInitLayersForFiltering) ? refVal : (refVal*interFilterTap)>>7;
              transformPredBuf[k][childIdx].val = filteredVal;
            }
          }
          enablePrediction = true;
        }

        if (isEncoder && enableIntraPrediction)
          fwdTransformBlock222<RahtKernel>(numAttrs, transformIntraPredBuf, weights);
      }

      if (isEncoder && curLevelEnableACInterPred)
        std::copy_n(&transformBuf[0][0], 8 * numAttrs, &transformIntraBuf[0][0]);

      // per-coefficient operations:
      //  - subtract transform domain prediction (encoder)
      //  - write out/read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      scanBlock(weights, [&](int idx) {
        // skip the DC coefficient unless at the root of the tree
        if (inheritDc && !idx)
          return;

        if (isEncoder && enablePrediction) {
          // subtract transformed prediction (skipping DC)
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][idx] -= transformPredBuf[k][idx];
          }
        }

        if (isEncoder && enableIntraPrediction) {
          for (int k = 0; k < numAttrs; k++) {
            transformIntraBuf[k][idx] -= transformIntraPredBuf[k][idx];
          }
        }

        // decision for RDOQ
        int64_t sumCoeff = 0;
        static const int LUTlog[16] = {0,   256, 406, 512, 594, 662, 719,  768,
                                812, 850, 886, 918, 947, 975, 1000, 1024};
        static const int LUTbins[11] = {1, 2, 3, 5, 5, 7, 7, 9, 9, 11, 11};
        bool flagRDOQ = false;
        int64_t intraSumCoeff = 0;
        bool intraFlagRDOQ = false;
        if (isEncoder && !rahtPredParams.integer_haar_enable_flag) {
          int64_t Dist2 = 0;
          int Ratecoeff = 0;
          int64_t lambda0;

          int64_t intraDist2 = 0;
          int intraRatecoeff = 0;

          for (int k = 0; k < numAttrs; k++) {
            //auto q = Quantizer(qpLayer[std::min(k, int(quantizers.size()) - 1)] + nodeQp[idx]);
            auto quantizers = qpset.quantizers(qpLayer, nodeQp[idx]);
            auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
            auto coeff = transformBuf[k][idx].round();
            Dist2 += coeff * coeff;
            auto Qcoeff = q.quantize(coeff << kFixedPointAttributeShift);
            sumCoeff += std::abs(Qcoeff);
            //Ratecoeff += !!Qcoeff; // sign
            Ratecoeff +=
              std::abs(Qcoeff) < 15 ? LUTlog[std::abs(Qcoeff)] : LUTlog[15];
            if (!k)
              lambda0 = q.scale(1);
            if (curLevelEnableACInterPred) {
              auto intraCoeff = transformIntraBuf[k][idx].round();
              intraDist2 += intraCoeff * intraCoeff;
              auto intraQcoeff =
                q.quantize(intraCoeff << kFixedPointAttributeShift);
              intraSumCoeff += std::abs(intraQcoeff);
              //Ratecoeff += !!Qcoeff; // sign
              intraRatecoeff += std::abs(intraQcoeff) < 15
                ? LUTlog[std::abs(intraQcoeff)]
                : LUTlog[15];
            }
          }
          const int64_t lambda = lambda0 * lambda0 * (numAttrs == 1 ? 25 : 35);
          if (sumCoeff < 3) {
            int Rate = LUTbins[trainZeros > 10 ? 10 : trainZeros];
            if (trainZeros > 10) {
              int temp = trainZeros - 11;
              // prefix k =2
              temp += 1;
              int a = 0;
              while (temp) {
                a++;
                temp >>= 1;
              }
              Rate += 2 * a - 1;
              // suffix  k=2
              Rate += 2;
            }
            //Rate = Rate / std::max(1, trainZeros);
            Rate += (Ratecoeff + 128) >> 8;

            flagRDOQ = (Dist2 << 26) < lambda * Rate;
          }
          if (curLevelEnableACInterPred && intraSumCoeff < 3) {
            int intraRate =
              LUTbins[intraTrainZeros > 10 ? 10 : intraTrainZeros];
            if (intraTrainZeros > 10) {
              int temp = intraTrainZeros - 11;
              // prefix k =2
              temp += 1;
              int a = 0;
              while (temp) {
                a++;
                temp >>= 1;
              }
              intraRate += 2 * a - 1;
              // suffix  k=2
              intraRate += 2;
            }
            intraRate += (intraRatecoeff + 128) >> 8;
            intraFlagRDOQ = (intraDist2 << 26) < lambda * intraRate;
          }

          // Track RL for RDOQ
          if (flagRDOQ || sumCoeff == 0)
            trainZeros++;
          else
            trainZeros = 0;

          if (curLevelEnableACInterPred) {
            if (intraFlagRDOQ || intraSumCoeff == 0)
              intraTrainZeros++;
            else
              intraTrainZeros = 0;
          }
        }

        // Check if QP offset is to be applied to AC coeffiecients
        Qps coeffQPOffset = (acCoeffQpLayer <= maxAcCoeffQpOffsetLayers && idx)
          ? qpset.rahtAcCoeffQps[acCoeffQpLayer][idx - 1]
          : Qps({0, 0});

        Qps nodeQPOffset = {
          nodeQp[idx][0] + coeffQPOffset[0],
          nodeQp[idx][1] + coeffQPOffset[1]};

        // The RAHT transform
        auto quantizers = qpset.quantizers(qpLayer, nodeQPOffset);
        for (int k = 0; k < numAttrs; k++) {

          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          if (isEncoder) {
            if (flagRDOQ)  // apply RDOQ
              transformBuf[k][idx].val = 0;

            if (intraFlagRDOQ)  // apply RDOQ
              transformIntraBuf[k][idx].val = 0;
            auto coeff = transformBuf[k][idx].round();
            assert(coeff <= INT_MAX && coeff >= INT_MIN);
            coeff = q.quantize(coeff << kFixedPointAttributeShift);

            //DC inter prediction at encoder
            auto coeff_tmp = coeff;
            if (curLevelEnableACInterPred)
              curEstimate.updateCostBits(coeff, k);
            *coeffBufItK[k]++ = coeff;
            transformPredBuf[k][idx] += divExp2RoundHalfUp(
              q.scale(coeff_tmp), kFixedPointAttributeShift);
            if (curLevelEnableACInterPred)
              curEstimate.resStatUpdate(coeff, k);
            if (curLevelEnableACInterPred) {  //< estimate
              auto intraCoeff = transformIntraBuf[k][idx].round();
              assert(intraCoeff <= INT_MAX && intraCoeff >= INT_MIN);
              intraCoeff = q.quantize(intraCoeff << kFixedPointAttributeShift);

              intraEstimate.updateCostBits(intraCoeff, k);
              *intraCoeffBufItK[k]++ = intraCoeff;
              transformIntraPredBuf[k][idx] += divExp2RoundHalfUp(
                q.scale(intraCoeff), kFixedPointAttributeShift);
              intraEstimate.resStatUpdate(intraCoeff, k);
            }
          } else {
            int64_t coeff = *coeffBufItK[k]++;

            transformPredBuf[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
          }
        }
      });

      // replace DC coefficient with parent if inheritable
      if (inheritDc) {
        for (int k = 0; k < numAttrs; k++) {
          attrRecParentIt++;
          int64_t val = *attrRecParentUsIt++;
          if (rahtExtension)
            transformPredBuf[k][0].val = val;
          else if (val > 0)
            transformPredBuf[k][0].val = val << (15 - 2);
          else
            transformPredBuf[k][0].val = -((-val) << (15 - 2));
          if (isEncoder && curLevelEnableACInterPred) {
            ///< inherit the parent DC coefficients
            transformIntraPredBuf[k][0].val = transformPredBuf[k][0].val;
          }
        }
      }

      if (rahtPredParams.integer_haar_enable_flag) {
        invTransformBlock222<HaarKernel>(numAttrs, transformPredBuf, weights);
        if (isEncoder && curLevelEnableACInterPred)
          invTransformBlock222<HaarKernel>(numAttrs, transformIntraPredBuf, weights);
      } else {
        invTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);
        if (isEncoder && curLevelEnableACInterPred)
          invTransformBlock222<RahtKernel>(numAttrs, transformIntraPredBuf, weights);
      }

      for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
        if (!weights[nodeIdx])
          continue;

        for (int k = 0; k < numAttrs; k++)
          if (rahtExtension) {
            attrRecUs[j * numAttrs + k] = transformPredBuf[k][nodeIdx].val;
            if (isEncoder && curLevelEnableACInterPred) {
              intraAttrRecUs[j * numAttrs + k] =
                transformIntraPredBuf[k][nodeIdx].val;
            }
          }
          else {
            FixedPoint temp = transformPredBuf[k][nodeIdx];
            temp.val <<= 2;
            attrRecUs[j * numAttrs + k] = temp.round();
            if (isEncoder && curLevelEnableACInterPred) {
              temp = transformIntraPredBuf[k][nodeIdx];
              temp.val <<= 2;
              intraAttrRecUs[j * numAttrs + k] = temp.round();
            }
          }

        // scale values for next level
        if (!rahtPredParams.integer_haar_enable_flag) {
        if (weights[nodeIdx] > 1) {
          FixedPoint rsqrtWeight;
          uint64_t w = weights[nodeIdx];
          int shift = w > 1024 ? ilog2(w - 1) >> 1 : 0;
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);
          for (int k = 0; k < numAttrs; k++) {
            transformPredBuf[k][nodeIdx].val >>= shift;
            transformPredBuf[k][nodeIdx] *= rsqrtWeight;
            if (isEncoder && curLevelEnableACInterPred) {
              transformIntraPredBuf[k][nodeIdx].val >>= shift;
              transformIntraPredBuf[k][nodeIdx] *= rsqrtWeight;
            }
          }
        }
        }

        for (int k = 0; k < numAttrs; k++) {
          attrRec[j * numAttrs + k] = rahtExtension
            ? transformPredBuf[k][nodeIdx].val
            : transformPredBuf[k][nodeIdx].round();
          if (isEncoder && curLevelEnableACInterPred) {
            intraAttrRec[j * numAttrs + k] = rahtExtension
              ? transformIntraPredBuf[k][nodeIdx].val
              : transformIntraPredBuf[k][nodeIdx].round();
          }
        }
        j++;
      }
      // increment reference buffer index if inter prediction is enabled
    }


    if (isEncoder && curLevelEnableACInterPred) {
      double curCost = curEstimate.costBits();
      double intraCost = intraEstimate.costBits();
      if (intraCost < curCost) {
        for (int k = 0; k < numAttrs; ++k)
          std::copy_n(intraCoeffBufItBeginK[k], sumNodes, coeffBufItBeginK[k]);
        std::swap(intraAttrRec, attrRec);
        std::swap(intraAttrRecUs, attrRecUs);
        curEstimate = intraEstimate;
        attrInterPredParams.attr_layer_code_mode.push_back(
          0);
        trainZeros = intraTrainZeros;
      } else {
        intraEstimate = curEstimate;
        attrInterPredParams.attr_layer_code_mode.push_back(1);
        intraTrainZeros = trainZeros;
      }
      curEstimate.resetCostBits();
      intraEstimate.resetCostBits();
    }
    if (enablePredictionInLvl && enableACRDOInterPred)
      ++depth;
    sumNodes = 0;
    // preserve current weights/positions for later search
    weightsParent = weightsLf;
    // increment tree depth
    treeDepth++;
  }
  // process duplicate points at level 0
  if (numDupNodes) {
    std::swap(attrRec, attrRecParent);
    auto attrRecParentIt = attrRecParent.cbegin();
    auto attrsHfIt = attrsHf.cbegin();

    for (int i = 0, out = 0, iEnd = weightsLf.size(); i < iEnd; i++) {
      int weight = weightsLf[i].weight;
      // unique points have weight = 1
      if (weight == 1) {
        for (int k = 0; k < numAttrs; k++)
          attrRec[out++] = *attrRecParentIt++;
        continue;
      }
      Qps nodeQp = {
        weightsLf[i].qp[0] >> regionQpShift,
        weightsLf[i].qp[1] >> regionQpShift};

      // duplicates
      FixedPoint attrSum[3];
      FixedPoint attrRecDc[3];
      FixedPoint sqrtWeight;
      sqrtWeight.val = isqrt(uint64_t(weight) << (2 * FixedPoint::kFracBits));

      int64_t sumCoeff = 0;
      for (int k = 0; k < numAttrs; k++) {
        if (isEncoder)
          attrSum[k] = attrsLf[i * numAttrs + k];
        if (rahtExtension)
          attrRecDc[k].val = *attrRecParentIt++;
        else
          attrRecDc[k] = *attrRecParentIt++;
        if (!rahtPredParams.integer_haar_enable_flag) {
          attrRecDc[k] *= sqrtWeight;
        }
      }

      FixedPoint rsqrtWeight;
      for (int w = weight - 1; w > 0; w--) {
        RahtKernel kernel(w, 1);
        HaarKernel haarkernel(w, 1);
        int shift = w > 1024 ? ilog2(uint32_t(w - 1)) >> 1 : 0;
        if (isEncoder)
          rsqrtWeight.val = irsqrt(w) >> (40 - shift - FixedPoint::kFracBits);

        auto quantizers = qpset.quantizers(qpLayer, nodeQp);
        for (int k = 0; k < numAttrs; k++) {
          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          FixedPoint transformBuf[2];
          if (isEncoder) {
            // invert the initial reduction (sum)
            // NB: read from (w-1) since left side came from attrsLf.
            transformBuf[1] = attrsHfIt[(w - 1) * numAttrs + k];
            if (rahtPredParams.integer_haar_enable_flag) {
              attrSum[k].val -= transformBuf[1].val >> 1;
              transformBuf[1].val += attrSum[k].val;
              transformBuf[0] = attrSum[k];
            } else {
              attrSum[k] -= transformBuf[1];
              transformBuf[0] = attrSum[k];

              // NB: weight of transformBuf[1] is by construction 1.
              transformBuf[0].val >>= shift;
              transformBuf[0] *= rsqrtWeight;
            }

            if (rahtPredParams.integer_haar_enable_flag) {
              haarkernel.fwdTransform(
                transformBuf[0], transformBuf[1], &transformBuf[0],
                &transformBuf[1]);
            } else {
              kernel.fwdTransform(
                transformBuf[0], transformBuf[1], &transformBuf[0],
                &transformBuf[1]);
            }

            auto coeff = transformBuf[1].round();
            assert(coeff <= INT_MAX && coeff >= INT_MIN);
            *coeffBufItK[k]++ = coeff =
              q.quantize(coeff << kFixedPointAttributeShift);
            transformBuf[1] =
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);

            sumCoeff +=
              std::abs(q.quantize(coeff << kFixedPointAttributeShift));
          } else {
            int64_t coeff = *coeffBufItK[k]++;
            transformBuf[1] =
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
          }

          // inherit the DC value
          transformBuf[0] = attrRecDc[k];

          if (rahtPredParams.integer_haar_enable_flag) {
            haarkernel.invTransform(
              transformBuf[0], transformBuf[1], &transformBuf[0],
              &transformBuf[1]);
          } else {
            kernel.invTransform(
              transformBuf[0], transformBuf[1], &transformBuf[0],
              &transformBuf[1]);
          }

          attrRecDc[k] = transformBuf[0];
          attrRec[out + w * numAttrs + k] =
            rahtExtension ? transformBuf[1].val : transformBuf[1].round();
          if (w == 1)
            attrRec[out + k] =
              rahtExtension ? transformBuf[0].val : transformBuf[0].round();
        }

        // Track RL for RDOQ
        if (isEncoder) {
          if (sumCoeff == 0)
            trainZeros++;
          else
            trainZeros = 0;
        }
      }

      attrsHfIt += (weight - 1) * numAttrs;
      out += weight * numAttrs;
    }
  }
  

  // write-back reconstructed attributes
  assert(attrRec.size() == numAttrs * numPoints);
  if(rahtExtension)
    for (auto& attr : attrRec) {
      attr += FixedPoint::kOneHalf;
      *(attributes++) = attr >> FixedPoint::kFracBits;
    }
  else
    std::copy(attrRec.begin(), attrRec.end(), attributes);
}

//============================================================================
/*
 * RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 *
 * Outputs:
 * weights = list of 'voxelCount' weights associated with each transform coefficient
 * coefficients = quantized transformed attributes array, in column-major order
 * binaryLayer = binary layer where each coefficient was generated
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalTransform(
  const RahtPredictionParams &rahtPredParams,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients,
  const bool rahtExtension,
  AttributeInterPredParams& attrInterPredParams)
{
  if (rahtExtension)
    uraht_process<true, true>(
      rahtPredParams, qpset, pointQpOffsets, voxelCount, attribCount, mortonCode,
      attributes, coefficients, attrInterPredParams);
  else
    uraht_process<true, false>(
      rahtPredParams, qpset, pointQpOffsets, voxelCount, attribCount, mortonCode,
      attributes, coefficients, attrInterPredParams);
}

//============================================================================
/*
 * inverse RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 * coefficients = quantized transformed attributes array, in column-major order
 *
 * Outputs:
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalInverseTransform(
  const RahtPredictionParams &rahtPredParams,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients,
  const bool rahtExtension,
  AttributeInterPredParams& attrInterPredParams)
{
  if (rahtExtension)
    uraht_process<false, true>(
      rahtPredParams, qpset, pointQpOffsets, voxelCount, attribCount,
      mortonCode, attributes, coefficients, attrInterPredParams);
  else
    uraht_process<false, false>(
      rahtPredParams, qpset, pointQpOffsets, voxelCount, attribCount,
      mortonCode, attributes, coefficients, attrInterPredParams);
}

//============================================================================

}  // namespace pcc
