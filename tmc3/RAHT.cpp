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

struct UrahtNode {
  int64_t pos;
  int weight;
  int qp;
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
  std::vector<int>* attrsOut)
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
      *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
      attrsOut->push_back(*attrsInRdIt++);
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
  std::vector<int>* attrsOut)
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
      (weightsInWrIt - 1)->weight += node.weight;
      (weightsInWrIt - 1)->qp = ((weightsInWrIt - 1)->qp + node.qp) >> 1;
      weightsOut->push_back(node);

      for (int k = 0; k < numAttrs; k++) {
        *(attrsInWrIt - numAttrs + k) += *attrsInRdIt;
        attrsOut->push_back(*attrsInRdIt++);
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
  std::vector<UrahtNode>* weightsIn,   // expand by numNodes before expand
  std::vector<UrahtNode>* weightsOut,  // shrink after expand
  std::vector<int>* attrsIn,
  std::vector<int>* attrsOut)
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
      *attrsInWrIt++ -= *curAttrIt++;
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
  int level,
  uint8_t occupancy,
  int parentNeighIdx[19],
  int parentNeighWeights[19])
{
  static const uint8_t neighMasks[19] = {255, 15, 240, 51, 204, 85,  170,
                                         3,   12, 5,   10, 48,  192, 80,
                                         160, 17, 34,  68, 136};

  // current position (discard extra precision)
  int64_t cur_pos = it->pos >> level;

  // the position of the parent, offset by (-1,-1,-1)
  int64_t base_pos = morton3dAdd(cur_pos, -1ll);

  // these neighbour offsets are relative to base_pos
  static const uint8_t neighOffset[19] = {0,  3,  35, 5,  21, 6, 14, 1,  17, 2,
                                          10, 33, 49, 34, 42, 4, 12, 20, 28};

  // special case for the direct parent (no need to search);
  parentNeighIdx[0] = std::distance(first, it);
  parentNeighWeights[0] = it->weight;

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
    parentNeighWeights[i] = found->weight;
  }
}

//============================================================================
// Generate the spatial prediction of a block.

template<typename It>
void
intraDcPred(
  int numAttrs,
  const int neighIdx[19],
  const int neighWeights[19],
  int occupancy,
  It first,
  FixedPoint predBuf[][8])
{
  static const uint8_t predMasks[19] = {255, 15, 240, 51, 204, 85,  170,
                                        3,   12, 5,   10, 48,  192, 80,
                                        160, 17, 34,  68, 136};

  static const FixedPoint predWeight[19] = {3.8, 2,   2,   2,   2,   2,   2,
                                            0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,
                                            0.7, 0.7, 0.7, 0.7, 0.7};

  FixedPoint weightSum[8] = {};
  std::fill_n(&predBuf[0][0], 8 * numAttrs, FixedPoint(0));

  for (int i = 0; i < 19; i++) {
    if (neighIdx[i] == -1)
      continue;

    // apply weighted neighbour value to masked positions
    auto neighValueIt = std::next(first, numAttrs * neighIdx[i]);
    FixedPoint neighValue[3];
    for (int k = 0; k < numAttrs; k++) {
      neighValue[k] = *neighValueIt++;
      neighValue[k] *= predWeight[i];
    }

    uint8_t mask = predMasks[i] & occupancy;
    for (int j = 0; j < 8; j++) {
      if (!(mask & (1 << j)))
        continue;

      weightSum[j] += predWeight[i];
      for (int k = 0; k < numAttrs; k++)
        predBuf[k][j] += neighValue[k];
    }
  }

  // normalise
  for (int i = 0; i < 8; i++) {
    if (!(occupancy & (1 << i)))
      continue;

    for (int k = 0; k < numAttrs; k++)
      predBuf[k][i] /= weightSum[i];
  }
}

//============================================================================
// Encapsulation of a RAHT transform stage.

class RahtKernel {
public:
  RahtKernel(int weightLeft, int weightRight)
  {
    int w = weightLeft + weightRight;
    _a.val = isqrt((int64_t(weightLeft) << (2 * _a.kFracBits)) / w);
    _b.val = isqrt((int64_t(weightRight) << (2 * _b.kFracBits)) / w);
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
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

template<class Kernel>
void
fwdTransformBlock222(int numBufs, FixedPoint buf[][8], int weights[8 + 4 + 2])
{
  for (int iw = 0, stride = 1; stride < 8; stride <<= 1) {
    for (int i0 = 0; i0 < 8; i0 += 2 * stride, iw += 2) {
      int i1 = i0 + stride;

      // the following tests can be simplified using an occupancy mask
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
}

//============================================================================
// In-place inverse transform a set of sparse 2x2x2 blocks each using the
// same weights

template<class Kernel>
void
invTransformBlock222(int numBufs, FixedPoint buf[][8], int weights[8 + 4 + 2])
{
  for (int iw = 12, stride = 4; stride > 0; stride >>= 1) {
    for (int i1 = 8 - stride; i1 > 0; i1 -= 2 * stride, iw -= 2) {
      int i0 = i1 - stride;

      // the following tests can be simplified using an occupancy mask
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
}

//============================================================================
// expand a set of eight weights into three levels

void
mkWeightTree(int weights[8 + 4 + 2])
{
  int* in = &weights[0];
  int* out = &weights[8];
  for (int i = 0; i < 6; i++) {
    *out++ = in[0] + in[1];
    in += 2;
  }
}

//============================================================================
// Invoke mapFn(coefIdx) for each present coefficient in the transform

template<class T>
void
scanBlock(int weights[8 + 4 + 2], T mapFn)
{
  static const int8_t kRahtScanOrder[] = {0, 4, 6, 2, 7, 5, 3, 1};

  // there is always the DC coefficient (empty blocks are not transformed)
  mapFn(0);

  for (int i = 1, iw = 12; iw >= 0; i++, iw -= 2) {
    if (!weights[iw] || !weights[iw + 1])
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
// Core transform process (for encoder/decoder)

template<bool isEncoder>
void
uraht_process(
  bool raht_prediction_enabled_flag,
  const int predictionThreshold[2],
  const std::vector<Qps>& qpLayers,
  int numPoints,
  int numAttrs,
  int64_t* positions,
  int* attributes,
  int32_t* coeffBufIt,
  int* regionQpOffset)
{
  // coefficients are stored in three planar arrays.  coeffBufItK is a set
  // of iterators to each array.
  int32_t* coeffBufItK[3] = {
    coeffBufIt,
    coeffBufIt + numPoints,
    coeffBufIt + numPoints * 2,
  };

  std::vector<UrahtNode> weightsLf, weightsHf;
  std::vector<int> attrsLf, attrsHf;

  weightsLf.reserve(numPoints);
  attrsLf.reserve(numPoints * numAttrs);

  int regionQpShift = 4;

  // copy positions into internal form
  // todo(df): lift to api
  for (int i = 0; i < numPoints; i++) {
    weightsLf.emplace_back(
      UrahtNode{positions[i], 1, regionQpOffset[i] << regionQpShift});
    for (int k = 0; k < numAttrs; k++) {
      attrsLf.push_back(attributes[i * numAttrs + k]);
    }
  }

  weightsHf.reserve(numPoints);
  attrsHf.reserve(numPoints * numAttrs);

  // ascend tree
  std::vector<int> levelHfPos;

  for (int level = 0, numNodes = weightsLf.size(); numNodes > 1; level++) {
    levelHfPos.push_back(weightsHf.size());
    if (level == 0) {
      // process any duplicate points
      numNodes = reduceUnique(
        numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    } else {
      // normal level reduction
      numNodes = reduceLevel(
        level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    }
  }

  assert(weightsLf[0].weight == numPoints);

  // reconstruction buffers
  std::vector<int> attrRec, attrRecParent;
  attrRec.resize(numPoints * numAttrs);
  attrRecParent.resize(numPoints * numAttrs);

  std::vector<int> attrRecUs, attrRecParentUs;
  attrRecUs.resize(numPoints * numAttrs);
  attrRecParentUs.resize(numPoints * numAttrs);

  std::vector<UrahtNode> weightsParent;
  weightsParent.reserve(numPoints);

  std::vector<int> numParentNiegh, numGrandParentNeigh;
  numParentNiegh.reserve(numPoints);
  numGrandParentNeigh.reserve(numPoints);

  // quant layer selection
  auto qpLayerIt = qpLayers.begin();

  // descend tree
  weightsLf.resize(1);
  attrsLf.resize(numAttrs);
  for (int level = levelHfPos.size() - 1, isFirst = 1; level > 0; /*nop*/) {
    int numNodes = weightsHf.size() - levelHfPos[level];
    weightsLf.resize(weightsLf.size() + numNodes);
    attrsLf.resize(attrsLf.size() + numNodes * numAttrs);
    expandLevel(
      level, numNodes, numAttrs, &weightsLf, &weightsHf, &attrsLf, &attrsHf);
    weightsHf.resize(levelHfPos[level]);
    attrsHf.resize(levelHfPos[level] * numAttrs);

    // expansion of level is complete, processing is now on the next level
    level--;

    // every three levels, perform transform
    if (level % 3)
      continue;

    // initial scan position of the coefficient buffer
    //  -> first level = all coeffs
    //  -> otherwise = ac coeffs only
    bool inheritDc = !isFirst;
    bool enablePredictionInLvl = inheritDc && raht_prediction_enabled_flag;
    isFirst = 0;

    // select quantiser according to transform layer
    const auto& qp = *qpLayerIt;
    if (std::next(qpLayerIt) != qpLayers.end())
      qpLayerIt++;

    // prepare reconstruction buffers
    //  previous reconstruction -> attrRecParent
    std::swap(attrRec, attrRecParent);
    std::swap(attrRecUs, attrRecParentUs);
    std::swap(numParentNiegh, numGrandParentNeigh);
    auto attrRecParentUsIt = attrRecParentUs.cbegin();
    auto attrRecParentIt = attrRecParent.cbegin();
    auto weightsParentIt = weightsParent.cbegin();
    auto numGrandParentNeighIt = numGrandParentNeigh.cbegin();

    for (int i = 0, iLast, iEnd = weightsLf.size(); i < iEnd; i = iLast) {
      // todo(df): hoist and dynamically allocate
      FixedPoint transformBuf[6][8] = {};
      FixedPoint(*transformPredBuf)[8] = &transformBuf[numAttrs];
      int weights[8 + 4 + 2] = {};
      int nodeQp[8 + 4 + 2] = {};
      uint8_t occupancy = 0;

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
      for (iLast = i; iLast < iEnd; iLast++) {
        int nextNode = iLast > i
          && !isSibling(weightsLf[iLast].pos, weightsLf[i].pos, level + 3);
        if (nextNode)
          break;

        int nodeIdx = (weightsLf[iLast].pos >> level) & 0x7;
        weights[nodeIdx] = weightsLf[iLast].weight;
        nodeQp[nodeIdx] = weightsLf[iLast].qp >> regionQpShift;

        occupancy |= 1 << nodeIdx;

        if (isEncoder) {
          for (int k = 0; k < numAttrs; k++)
            transformBuf[k][nodeIdx] = attrsLf[iLast * numAttrs + k];
        }
      }

      mkWeightTree(weights);

      if (!inheritDc) {
        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNiegh[j++] = 19;
        }
      }

      // Inter-level prediction:
      //  - Find the parent neighbours of the current node
      //  - Generate prediction for all attributes into transformPredBuf
      //  - Subtract transformed coefficients from forward transform
      //  - The transformPredBuf is then used for reconstruction
      bool enablePrediction = enablePredictionInLvl;
      if (enablePredictionInLvl) {
        // indexes of the neighbouring parents
        int parentNeighIdx[19];
        int parentNeighWeights[19];

        int parentNeighCount = 0;
        if (*numGrandParentNeighIt < predictionThreshold[0]) {
          enablePrediction = false;
        } else {
          findNeighbours(
            weightsParent.cbegin(), weightsParent.cend(), weightsParentIt,
            level + 3, occupancy, parentNeighIdx, parentNeighWeights);
          for (int i = 0; i < 19; i++) {
            parentNeighCount += (parentNeighIdx[i] != -1);
          }
          if (parentNeighCount < predictionThreshold[1]) {
            enablePrediction = false;
          } else
            intraDcPred(
              numAttrs, parentNeighIdx, parentNeighWeights, occupancy,
              attrRecParent.begin(), transformPredBuf);
        }

        for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
          if (!weights[nodeIdx])
            continue;
          numParentNiegh[j++] = parentNeighCount;
        }
      }

      int parentWeight = 0;
      if (inheritDc) {
        parentWeight = weightsParentIt->weight;
        weightsParentIt++;
        numGrandParentNeighIt++;
      }

      // normalise coefficients
      for (int childIdx = 0; childIdx < 8; childIdx++) {
        if (weights[childIdx] <= 1)
          continue;

        FixedPoint sqrtWeight;
        sqrtWeight.val =
          isqrt(uint64_t(weights[childIdx]) << (2 * FixedPoint::kFracBits));

        // Summed attribute values
        if (isEncoder) {
          for (int k = 0; k < numAttrs; k++)
            transformBuf[k][childIdx] /= sqrtWeight;
        }

        // Predicted attribute values
        if (enablePrediction) {
          for (int k = 0; k < numAttrs; k++)
            transformPredBuf[k][childIdx] *= sqrtWeight;
        }
      }

      // forward transform:
      //  - encoder: transform both attribute sums and prediction
      //  - decoder: just transform prediction
      if (isEncoder && enablePrediction)
        fwdTransformBlock222<RahtKernel>(2 * numAttrs, transformBuf, weights);
      else if (isEncoder)
        fwdTransformBlock222<RahtKernel>(numAttrs, transformBuf, weights);
      else if (enablePrediction)
        fwdTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);

      // per-coefficient operations:
      //  - subtract transform domain prediction (encoder)
      //  - write out/read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      scanBlock(weights, [&](int idx) {
        // skip the DC coefficient unless at the root of the tree
        if (inheritDc && !idx)
          return;

        // subtract transformed prediction (skipping DC)
        if (isEncoder && enablePrediction) {
          for (int k = 0; k < numAttrs; k++) {
            transformBuf[k][idx] -= transformPredBuf[k][idx];
          }
        }

        // The RAHT transform
        for (int k = 0; k < numAttrs; k++) {
          // todo: hoist to preallocated array
          auto q =
            Quantizer(qp[std::min(k, int(qp.size()) - 1)] + nodeQp[idx]);

          if (isEncoder) {
            auto coeff = transformBuf[k][idx].round();
            assert(coeff <= INT_MAX && coeff >= INT_MIN);
            *coeffBufItK[k]++ = coeff =
              q.quantize(coeff << kFixedPointAttributeShift);
            transformPredBuf[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
          } else {
            int64_t coeff = *coeffBufItK[k]++;
            transformPredBuf[k][idx] +=
              divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
          }
        }
      });

      // replace DC coefficient with parent if inheritable
      if (inheritDc) {
        FixedPoint sqrtWeight;
        sqrtWeight.val =
          isqrt(uint64_t(parentWeight) << (2 * FixedPoint::kFracBits));

        for (int k = 0; k < numAttrs; k++) {
          attrRecParentIt++;
          transformPredBuf[k][0] = *attrRecParentUsIt++;
        }
      }

      invTransformBlock222<RahtKernel>(numAttrs, transformPredBuf, weights);

      for (int j = i, nodeIdx = 0; nodeIdx < 8; nodeIdx++) {
        if (!weights[nodeIdx])
          continue;

        for (int k = 0; k < numAttrs; k++)
          attrRecUs[j * numAttrs + k] = transformPredBuf[k][nodeIdx].round();

        // scale values for next level
        if (weights[nodeIdx] > 1) {
          FixedPoint sqrtWeight;
          sqrtWeight.val =
            isqrt(uint64_t(weights[nodeIdx]) << (2 * FixedPoint::kFracBits));
          for (int k = 0; k < numAttrs; k++)
            transformPredBuf[k][nodeIdx] /= sqrtWeight;
        }

        for (int k = 0; k < numAttrs; k++)
          attrRec[j * numAttrs + k] = transformPredBuf[k][nodeIdx].round();
        j++;
      }
    }

    // preserve current weights/positions for later search
    weightsParent = weightsLf;
  }

  // process duplicate points at level 0
  std::swap(attrRec, attrRecParent);
  auto attrRecParentIt = attrRecParent.cbegin();
  auto attrsHfIt = attrsHf.cbegin();
  const auto& qp = qpLayers.back();

  for (int i = 0, out = 0, iEnd = weightsLf.size(); i < iEnd; i++) {
    int weight = weightsLf[i].weight;
    int nodeQp = weightsLf[i].qp >> regionQpShift;

    // unique points have weight = 1
    if (weight == 1) {
      for (int k = 0; k < numAttrs; k++)
        attrRec[out++] = *attrRecParentIt++;
      continue;
    }

    // duplicates
    FixedPoint attrSum[3];
    FixedPoint attrRecDc[3];
    FixedPoint sqrtWeight;
    sqrtWeight.val = isqrt(uint64_t(weight) << (2 * FixedPoint::kFracBits));
    for (int k = 0; k < numAttrs; k++) {
      if (isEncoder)
        attrSum[k] = attrsLf[i * numAttrs + k];
      attrRecDc[k] = *attrRecParentIt++;
      attrRecDc[k] *= sqrtWeight;
    }

    for (int w = weight - 1; w > 0; w--) {
      RahtKernel kernel(w, 1);
      sqrtWeight.val = isqrt(uint64_t(w) << (2 * FixedPoint::kFracBits));

      for (int k = 0; k < numAttrs; k++) {
        auto q = Quantizer(qp[std::min(k, int(qp.size()) - 1)] + nodeQp);

        FixedPoint transformBuf[2];
        if (isEncoder) {
          // invert the initial reduction (sum)
          // NB: read from (w-1) since left side came from attrsLf.
          transformBuf[1] = attrsHfIt[(w - 1) * numAttrs + k];
          attrSum[k] -= transformBuf[1];
          transformBuf[0] = attrSum[k];

          // NB: weight of transformBuf[1] is by construction 1.
          transformBuf[0] /= sqrtWeight;

          kernel.fwdTransform(
            transformBuf[0], transformBuf[1], &transformBuf[0],
            &transformBuf[1]);

          auto coeff = transformBuf[1].round();
          assert(coeff <= INT_MAX && coeff >= INT_MIN);
          *coeffBufItK[k]++ = coeff =
            q.quantize(coeff << kFixedPointAttributeShift);
          transformBuf[1] =
            divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
        } else {
          int64_t coeff = *coeffBufItK[k]++;
          transformBuf[1] =
            divExp2RoundHalfUp(q.scale(coeff), kFixedPointAttributeShift);
        }

        // inherit the DC value
        transformBuf[0] = attrRecDc[k];

        kernel.invTransform(
          transformBuf[0], transformBuf[1], &transformBuf[0],
          &transformBuf[1]);

        attrRecDc[k] = transformBuf[0];
        attrRec[out + w * numAttrs + k] = transformBuf[1].round();
        if (w == 1)
          attrRec[out + k] = transformBuf[0].round();
      }
    }

    attrsHfIt += (weight - 1) * numAttrs;
    out += weight * numAttrs;
  }

  // write-back reconstructed attributes
  assert(attrRec.size() == numAttrs * numPoints);
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
  bool raht_prediction_enabled_flag,
  const int predictionThreshold[2],
  const std::vector<Qps>& qpLayers,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients,
  int* regionQPOffset)
{
  uraht_process<true>(
    raht_prediction_enabled_flag, predictionThreshold, qpLayers, voxelCount,
    attribCount, mortonCode, attributes, coefficients, regionQPOffset);
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
  bool raht_prediction_enabled_flag,
  const int predictionThreshold[2],
  const std::vector<Qps>& qpLayers,
  int64_t* mortonCode,
  int* attributes,
  const int attribCount,
  const int voxelCount,
  int* coefficients,
  int* regionQPOffset)
{
  uraht_process<false>(
    raht_prediction_enabled_flag, predictionThreshold, qpLayers, voxelCount,
    attribCount, mortonCode, attributes, coefficients, regionQPOffset);
}

//============================================================================

}  // namespace pcc
