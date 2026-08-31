# Charge Solving Algorithm

## Overview

The charge solving pipeline determines blob charges from wire measurements using
LASSO (Least Absolute Shrinkage and Selection Operator) regression. The core idea:
wire measurements are linear combinations of blob charges, so we solve an
over-determined system `A * x = m` where `A` is the blob-measurement association
matrix, `x` is blob charges, and `m` is wire measurements.

## 1. ChargeSolving (`src/ChargeSolving.cxx`)

### Pipeline Orchestration

The `operator()` method runs:

```
Input Cluster -> Unpack -> [Weight -> Solve -> Prune] x N_strategies -> Repack -> Output
```

1. **Unpack** (`CS::unpack`): Decompose cluster graph into per-slice bipartite
   blob-measurement subgraphs. Each subgraph is further split into connected
   components so each solve is independent.

2. **Weight + Solve + Prune**: For each configured weighting strategy (default: "uboone"):
   ```cpp
   std::transform(sgs.begin(), sgs.end(), sgs.begin(),
       [&](graph_t& sg) {
           blob_weighter(in_graph, sg);
           auto tmp_csg = solve(sg, sparams);
           return prune(tmp_csg, blob_threshold[ind]);
       });
   ```

3. **Repack** (`CS::repack`): Replace blob values in original cluster graph with
   solved charges, removing pruned blobs.

### Blob Weighting Strategies

Weight is stored in `vtx.value.uncertainty()` (overloading the uncertainty field).

**`uniform`** (line 53-61):
- All blobs get weight 9.0
- Simplest strategy, no topology awareness

**`simple`** (line 64-85):
- Weight = number of unique slice idents connected via b-b edges + 1 (current slice)
- More connected blobs get higher weight (less regularization)
- The "+1" comes from inserting the current slice ident first

**`uboone`** (line 90-138):
- Base weight = 9.0
- If blob connects to a blob in next slice with charge > 300: divide by 3
- If blob connects to a blob in prev slice with charge > 300: divide by 3
- Result: isolated blobs get weight 9, one-sided connected get 3, both-sided get 1
- This strongly penalizes isolated blobs during LASSO solving

| Connectivity | Weight |
|-------------|--------|
| No connections | 9.0 |
| One direction (prev or next) | 3.0 |
| Both directions | 1.0 |

---

## 2. CSGraph (`src/CSGraph.cxx`, `inc/WireCellImg/CSGraph.h`)

### CS::unpack (lines 291-385)

Converts the cluster graph into per-slice bipartite subgraphs:

1. Find all blob nodes (`'b'`), for each find its slice (`'s'`) and measures (`'m'`)
2. Group by slice, creating a graph per slice with blob and measure nodes
3. Apply measurement threshold: skip measures with value < threshold or uncertainty > threshold
4. Critical check (line 351): `if (!(msum.uncertainty()>0))` -- skip measures with
   non-positive uncertainty (logs warning instead of throwing)
5. Split each per-slice graph into connected components via `connected_subgraphs()`

### CS::solve (lines 31-211)

The core LASSO solver:

1. **Setup**: Extract blob source values, weights, measurement values, and covariance
   - `source(bind)` = initial blob charge estimate
   - `weight(bind)` = blob weight from weighting strategy
   - `measure(mind)` = wire measurement value
   - `mcov(mind, mind)` = measurement uncertainty squared (diagonal covariance)

2. **Special case** (lines 93-104): If `config == simple` and only 1 blob, average
   all measurement values instead of running LASSO:
   ```cpp
   csg_out[nbdesc].value = val / nmeas;
   ```

3. **Build A matrix** (lines 106-132): Binary association matrix from blob-measure edges.
   `A(mind, bind) = 1` if blob `bind` is connected to measure `mind`.

4. **LASSO parameters** (lines 142-147, uboone config):
   ```cpp
   lambda = 3.0 / total_wire_charge / 2.0 * params.scale
   tolerance = total_wire_charge / 3.0 / params.scale / R_mat.cols() * 0.005
   ```
   - Lambda controls L1 regularization strength
   - Tolerance controls convergence criterion
   - Both scale inversely with total wire charge

5. **Cholesky whitening** (lines 157-170, if `params.whiten == true`):
   ```cpp
   Eigen::LLT<double_matrix_t> llt(mcov.inverse());
   double_matrix_t U = llt.matrixL().transpose();
   m_vec = U * measure;          // whitened measurements
   R_mat = params.scale * U * A; // whitened response
   ```
   This transforms the problem into a whitened basis where measurement errors are
   unit-variance, improving LASSO conditioning.

6. **Solve** (line 182): `Ress::solve(R_mat, m_vec, rparams, source, weight)`
   - Uses external Ress library for LASSO optimization
   - `source` = initial guess, `weight` = regularization weights

7. **Update results** (lines 197-204):
   ```cpp
   bvalue.value = solution[ind] * params.scale;  // note: drops weight
   ```
   The comment "drops weight" means uncertainty is not updated from the solution.

### CS::prune (lines 213-252)

Simple threshold-based blob removal:
- Copy all non-blob nodes and blobs above threshold to output graph
- Copy edges where both endpoints survive

### CS::repack (lines 400-491)

Reconstruct original cluster graph with solved values:
1. Index surviving blobs and measures from solved subgraphs
2. Copy non-{b,m} nodes as-is, replace live blobs with solved values
3. Transfer edges where both endpoints exist in output

---

## 3. BlobSolving (`src/BlobSolving.cxx`)

A simpler, standalone LASSO solver variant (not using the CSGraph pipeline):
- Same weight scheme: base=9, divide by 3 per connection, max 2 slices checked
- Directly calls `Ress::solve` per slice
- Uncertainty always set to 0.0 (FIXME noted in code)

---

## 4. ChargeErrorFrameEstimator (`src/ChargeErrorFrameEstimator.cxx`)

Estimates per-channel charge uncertainty:
- Uses pre-computed error waveforms indexed by ROI (Region of Interest) length
- Applies plane-specific fudge factors
- Three modes based on time position:
  - Before `time_limits.first`: one error model
  - Between limits: standard model
  - After `time_limits.second`: another error model
- Error scales as `sqrt(waveform_value * fudge_factor)`

---

## Mathematical Summary

The charge solving problem per time slice:

Given:
- `m` measurements (wire signals): vector `m_vec` of size `nmeas`
- `n` blobs (unknowns): vector `x` of size `nblob`
- Association matrix `A` (nmeas x nblob): `A[i,j] = 1` if blob j covers measure i
- Measurement covariance `Sigma` (diagonal, `sigma_i^2`)
- Regularization weights `w` per blob

Solve (LASSO):
```
minimize ||U*A*x - U*m||_2^2 + lambda * sum(w_i * |x_i|)
```
where `U` is the Cholesky factor of `Sigma^{-1}` (whitening transform).

The regularization drives small/uncertain blob charges to zero, effectively
performing model selection -- blobs with insufficient evidence are removed.
