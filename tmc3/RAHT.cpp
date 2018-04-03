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

#include <vector>

#include "PCCMath.h"
#include "PCCMisc.h"

namespace pcc {

/*
 * RAHT
 *
 * Inputs:
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 * depth = number of bits of precision (e.g., 10) for each coordinate in mortonCode
 *
 * Outputs:
 * attributes = 'voxelCount' x 'attribCount' array of transformed attributes, in row-major order
 * weights = list of 'voxelCount' weights associated with each transform coefficient
 *
 * Note transform is done 'in-place' (i.e., array 'attributes' contains
 * attributes on input, and transformed attributes on output).
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void regionAdaptiveHierarchicalTransform(
  long long *mortonCode, float *attributes, float *weight, int *binaryLayer,
  int attribCount, int voxelCount, int depth)
{
  // Prologue, common to both RAHT and IRAHT:

  //% Process one level at a time, from leaves(b = 1) to root(b = 3 * depth).
  // for b = 1:3 * depth
  //    if b == 1 % initialize indices of coeffs at level b
  //        I{ b } = uint32(1:N)'; % vector of indices from 1 to N
  //    else % define indices of coeffs at level b
  //        I{ b } = I{ b - 1 }(~[0; F{ b - 1 }]);
  //    end
  //    Mb = M1(I{ b }); % Morton codes at level b
  //    W{ b } = double([I{ b }(2:end); N + 1] - I{ b }); % weights
  //    D = bitxor(Mb(1:end - 1), Mb(2:end)); % path diffs
  //    F{ b } = (bitand(D, (bitshift(c1, 3 * depth) - bitshift(c1, b)))) == 0; % is left sibling
  // end

  // Allocate I, W, F.
  // I[b] = indices of Morton codes of lead voxels at binary level b
  // W[b] = weights of cells at binary level b
  // F[b] = boolean flags to indicate left siblings at binary level b
  int **I = new int *[3 * depth];
  int **W = new int *[3 * depth];
  bool **F = new bool *[3 * depth];
  int *N = new int[3 * depth];                // length of above arrays at each level
  long long *Mb = new long long[voxelCount];  // Morton codes of lead voxels at binary level b

  // Process one level at a time, from leaves (b=0) to root (b=3*depth-1).
  int Nparent = voxelCount;
  for (int b = 0; b < 3 * depth; b++) {
    // Get indices for each node at this level.
    int Nb = Nparent;
    N[b] = Nb;
    int *Ib = new int[Nb];
    I[b] = Ib;
    if (b == 0) {
      // Initialize indices of coefficients at level b.
      for (int n = 0; n < Nb; n++) Ib[n] = n;
    } else {
      // Define indices of coefficients at level b.
      int Nchild = N[b - 1];    // number of nodes at child level
      int *Ichild = I[b - 1];   // indices of nodes at child level
      bool *Fchild = F[b - 1];  // left-sibling flags at child level
      int i = 0;
      Ib[i++] = Ichild[0];  // always promote index of child[0] to parent
      for (int n = 1; n < Nchild; n++) {
        // If node to left of child[n] is not a left child,
        // then child[n] is not a right child, i.e.,
        // child[n] is a left child or an orphan,
        // in which case promote its index to its parent
        if (Fchild[n - 1] == false) {
          Ib[i++] = Ichild[n];
        }
      }
      // assert i==Nb
    }

    // Collect Morton codes for each node at this level.
    for (int n = 0; n < Nb; n++) Mb[n] = mortonCode[Ib[n]];

    // Compute Weights for each node at this level.
    int *Wb = new int[Nb];
    W[b] = Wb;
    for (int n = 0; n < Nb - 1; n++) Wb[n] = Ib[n + 1] - Ib[n];
    Wb[Nb - 1] = voxelCount - Ib[Nb - 1];

    // Flag whether this is a left sibling for each node at this level.
    bool *Fb = new bool[Nb - 1];
    F[b] = Fb;
    long long mask =
        ((long long)1 << (3 * depth)) - ((long long)1 << (b + 1));  // turn on all bits above b
    Nparent = 1;  // number of nodes at parent level
    for (int n = 0; n < Nb - 1; n++) {
      long long D = Mb[n + 1] ^ Mb[n];  // path difference
      Fb[n] = (D & mask) == 0;          // is left sibling if prefix is shared
      if (Fb[n] == false) Nparent++;    // count non-left-siblings
    }
  }

  // RAHT-specific part:

  //% Initialize weights to all ones.
  // weights = ones(size(transformedAttributes,1),1);
  //
  //% Process one level at a time, from leaves (b=1) to root (b=3*depth).
  // for b = 1:3*depth
  //    i0 = I{b}([F{b};0] == 1); % left sibling indices
  //    if isempty(i0)
  //        continue;
  //    end
  //    i1 = I{b}([0;F{b}] == 1); % right sibling indices
  //    w0 = W{b}([F{b};0] == 1); % left sibling weights
  //    w1 = W{b}([0;F{b}] == 1); % right sibling weights
  //    x0 = transformedAttributes(i0,:); % left sibling coefficients
  //    x1 = transformedAttributes(i1,:); % right sibling coefficients
  //    alpha = repmat(single(sqrt(w0 ./ (w0+w1))),1,size(transformedAttributes,2));
  //    beta = repmat(single(sqrt(w1 ./ (w0+w1))),1,size(transformedAttributes,2));
  //    transformedAttributes(i0,:) = alpha .* x0 + beta .* x1;
  //    transformedAttributes(i1,:) = -beta .* x0 + alpha .* x1;
  //    weights(i0) = weights(i0) + weights(i1);
  //    weights(i1) = weights(i0); % same as i0
  // end

  // Allocate local workspace for attributes.
  float *x0 = new float[attribCount];
  float *x1 = new float[attribCount];

  // Initialize weights.
  for (int n = 0; n < voxelCount; n++) {
    weight[n] = 1.0f;
    binaryLayer[n] = 0;
  }

  // Process one level at a time, from leaves (b=0) to root (b=3*depth-1).
  for (int b = 0; b < 3 * depth; b++) {
    int Nb = N[b];
    int *Ib = I[b];
    int *Wb = W[b];
    bool *Fb = F[b];

    // Process all (left,right) siblings.
    for (int n = 0; n < Nb - 1; n++) {
      if (Fb[n] == true) {
        // Found a (left,right) sibling.
        // Copy attributes into local workspace.
        int i0 = Ib[n];
        int i1 = Ib[n + 1];
        for (int k = 0; k < attribCount; k++) {
          x0[k] = attributes[attribCount * i0 + k];
          x1[k] = attributes[attribCount * i1 + k];
        }
        // Perform Givens rotation on attributes of (left,right) siblings.
        int w0 = Wb[n];
        int w1 = Wb[n + 1];
        float frac = (float)w0 / (float)(w0 + w1);
        float alpha = sqrt(frac);
        float beta = sqrt(1.0f - frac);
        for (int k = 0; k < attribCount; k++) {
          attributes[attribCount * i0 + k] = alpha * x0[k] + beta * x1[k];
          attributes[attribCount * i1 + k] = -beta * x0[k] + alpha * x1[k];
        }
        // Set weights.
        weight[i0] = (w0 + w1);
        weight[i1] = (w0 + w1);

        binaryLayer[i0] = b + 1;

        n++;  // skip over right sibling
      }
    }
  }

  // Free arrays.
  delete[] x0;
  delete[] x1;
  for (int b = 0; b < 3 * depth; b++) {
    delete[] I[b];
    delete[] W[b];
    delete[] F[b];
  }
  delete[] I;
  delete[] W;
  delete[] F;
  delete[] N;
  delete[] Mb;
}

/*
 * inverse RAHT

 * Inputs:
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attributes = 'voxelCount' x 'attribCount' array of transformed attributes, in row-major order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 * depth = number of bits of precision (e.g., 10) for each coordinate in mortonCode
 *
 * Outputs:
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 *
 * Note transform is done 'in-place' (i.e., array 'attributes' contains
 * transformed attributes on input, and attributes on output).
 */
void regionAdaptiveHierarchicalInverseTransform(
  long long *mortonCode, float *attributes,
  int attribCount, int voxelCount, int depth)
{
  // Prologue, common to both RAHT and IRAHT:

  //% Process one level at a time, from leaves(b = 1) to root(b = 3 * depth).
  // for b = 1:3 * depth
  //    if b == 1 % initialize indices of coeffs at level b
  //        I{ b } = uint32(1:N)'; % vector of indices from 1 to N
  //    else % define indices of coeffs at level b
  //        I{ b } = I{ b - 1 }(~[0; F{ b - 1 }]);
  //    end
  //    Mb = M1(I{ b }); % Morton codes at level b
  //    W{ b } = double([I{ b }(2:end); N + 1] - I{ b }); % weights
  //    D = bitxor(Mb(1:end - 1), Mb(2:end)); % path diffs
  //    F{ b } = (bitand(D, (bitshift(c1, 3 * depth) - bitshift(c1, b)))) == 0; % is left sibling
  // end

  // Allocate I, W, F.
  // I[b] = indices of Morton codes of lead voxels at binary level b
  // W[b] = weights of cells at binary level b
  // F[b] = boolean flags to indicate left siblings at binary level b
  int **I = new int *[3 * depth];
  int **W = new int *[3 * depth];
  bool **F = new bool *[3 * depth];
  int *N = new int[3 * depth];                // length of above arrays at each level
  long long *Mb = new long long[voxelCount];  // Morton codes of lead voxels at binary level b

  // Process one level at a time, from leaves (b=0) to root (b=3*depth-1).
  int Nparent = voxelCount;
  for (int b = 0; b < 3 * depth; b++) {
    // Get indices for each node at this level.
    int Nb = Nparent;
    N[b] = Nb;
    int *Ib = new int[Nb];
    I[b] = Ib;
    if (b == 0) {
      // Initialize indices of coefficients at level b.
      for (int n = 0; n < Nb; n++) Ib[n] = n;
    } else {
      // Define indices of coefficients at level b.
      int Nchild = N[b - 1];    // number of nodes at child level
      int *Ichild = I[b - 1];   // indices of nodes at child level
      bool *Fchild = F[b - 1];  // left-sibling flags at child level
      int i = 0;
      Ib[i++] = Ichild[0];  // always promote index of child[0] to parent
      for (int n = 1; n < Nchild; n++) {
        // If node to left of child[n] is not a left child,
        // then child[n] is not a right child, i.e.,
        // child[n] is a left child or an orphan,
        // in which case promote its index to its parent
        if (Fchild[n - 1] == false) {
          Ib[i++] = Ichild[n];
        }
      }
      // assert i==Nb
    }

    // Collect Morton codes for each node at this level.
    for (int n = 0; n < Nb; n++) Mb[n] = mortonCode[Ib[n]];

    // Compute Weights for each node at this level.
    int *Wb = new int[Nb];
    W[b] = Wb;
    for (int n = 0; n < Nb - 1; n++) Wb[n] = Ib[n + 1] - Ib[n];
    Wb[Nb - 1] = voxelCount - Ib[Nb - 1];

    // Flag whether this is a left sibling for each node at this level.
    bool *Fb = new bool[Nb - 1];
    F[b] = Fb;
    long long mask =
        ((long long)1 << (3 * depth)) - ((long long)1 << (b + 1));  // turn on all bits above b
    Nparent = 1;  // number of nodes at parent level
    for (int n = 0; n < Nb - 1; n++) {
      long long D = Mb[n + 1] ^ Mb[n];  // path difference
      Fb[n] = (D & mask) == 0;          // is left sibling if prefix is shared
      if (Fb[n] == false) Nparent++;    // count non-left-siblings
    }
  }

  // IRAHT-specific part:

  //% Process one level at a time, from root(b = 3 * depth) to leaves(b = 1).
  // for b = 3 * depth:-1 : (3 * (depth - level) + 1)
  //    i0 = I{ b }([F{ b }; 0] == 1); % left sibling indices
  //    if isempty(i0)
  //        continue;
  //    end
  //    i1 = I{ b }([0; F{ b }] == 1); % right sibling indices
  //    w0 = W{ b }([F{ b }; 0] == 1); % left sibling weights
  //    w1 = W{ b }([0; F{ b }] == 1); % right sibling weights
  //    x0 = attributes(i0, :); % left sibling coefficients
  //    x1 = attributes(i1, :); % right sibling coefficients
  //    alpha = repmat(single(sqrt(w0./(w0 + w1))), 1, size(attributes, 2));
  //    beta = repmat(single(sqrt(w1./(w0 + w1))), 1, size(attributes, 2));
  //    attributes(i0, :) = alpha.*x0 - beta.*x1;
  //    attributes(i1, :) = beta.*x0 + alpha.*x1;
  // end

  // Allocate local workspace for attributes.
  float *x0 = new float[attribCount];
  float *x1 = new float[attribCount];

  // Process one level at a time, root (b=3*depth-1) to from leaves (b=0).
  for (int b = 3 * depth - 1; b >= 0; b--) {
    int Nb = N[b];
    int *Ib = I[b];
    int *Wb = W[b];
    bool *Fb = F[b];

    // Process all (left,right) siblings.
    for (int n = 0; n < Nb - 1; n++) {
      if (Fb[n] == true) {
        // Found a (left,right) sibling.
        // Copy attributes into local workspace.
        int i0 = Ib[n];
        int i1 = Ib[n + 1];
        for (int k = 0; k < attribCount; k++) {
          x0[k] = attributes[attribCount * i0 + k];
          x1[k] = attributes[attribCount * i1 + k];
        }
        // Perform Givens rotation on attributes of (left,right) siblings.
        int w0 = Wb[n];
        int w1 = Wb[n + 1];
        float frac = (float)w0 / (float)(w0 + w1);
        float alpha = sqrt(frac);
        float beta = sqrt(1.0f - frac);
        for (int k = 0; k < attribCount; k++) {
          attributes[attribCount * i0 + k] = alpha * x0[k] - beta * x1[k];
          attributes[attribCount * i1 + k] = beta * x0[k] + alpha * x1[k];
        }
        n++;  // skip over right sibling
      }
    }
  }

  // Free arrays.
  delete[] x0;
  delete[] x1;
  for (int b = 0; b < 3 * depth; b++) {
    delete[] I[b];
    delete[] W[b];
    delete[] F[b];
  }
  delete[] I;
  delete[] W;
  delete[] F;
  delete[] N;
  delete[] Mb;
}

} /* namespace pcc */
