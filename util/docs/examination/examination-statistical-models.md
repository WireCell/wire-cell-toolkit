# Statistical / Regression Models - Code Examination

## Overview

This group of code implements linear regression models with L1 (LASSO) and L2 (Ridge) regularization, their combination (Elastic Net), and the Kolmogorov-Smirnov (KS) statistical test for comparing distributions. The regression models use coordinate descent optimization, while the KS test code provides both one-sample and two-sample variants with asymptotic p-value computation.

## Files Examined

- `util/inc/WireCellUtil/LinearModel.h`
- `util/src/LinearModel.cxx`
- `util/inc/WireCellUtil/ElasticNetModel.h`
- `util/src/ElasticNetModel.cxx`
- `util/inc/WireCellUtil/LassoModel.h`
- `util/src/LassoModel.cxx`
- `util/inc/WireCellUtil/KSTest.h`
- `util/src/KSTest.cxx`

## Algorithm Description

### LinearModel

Base class for linear regression solving `y = X * beta`. Provides:
- `Predict()`: returns `X * beta`.
- `chi2_base()`: returns the squared L2 norm of the residual `||y - X*beta||_2^2`.
- `MeanResidual()`: returns `||y - X*beta||_1_norm / N` (see BUG-5 below on naming).

### ElasticNetModel

Minimizes: `(1/2) * ||y - X*beta||_2^2 + N*lambda*(alpha*||beta||_1 + 0.5*(1-alpha)*||beta||_2^2)`

Uses coordinate descent with active-set management:

1. **Column norm precomputation** (lines 50-58): computes `||X_j||_2^2` for each feature. Columns with near-zero norm are flagged and their norm is set to 1 to avoid division by zero.

2. **Coordinate update** (lines 64-85): for each active feature `j`, the full residual excluding feature `j` is computed as `r_j = y - X * beta_tmp` where `beta_tmp(j)=0`. The update is:
   ```
   beta(j) = S(X_j' * r_j / ||X_j||^2, lambda*alpha*w_j) / (1 + lambda*(1-alpha))
   ```
   where `S()` is the soft-thresholding operator.

3. **Active-set strategy** (lines 79-81, 86-102): coefficients that fall below `1e-6` are deactivated. When convergence is detected, all coefficients are reactivated for one final pass (`double_check` mechanism) to confirm the solution.

4. **Soft-thresholding** (lines 110-128): implements `S(x, lambda) = sign(x)*max(|x|-lambda, 0)`, with an option to enforce non-negativity by only returning the positive branch.

**Key algorithmic note on the Ridge denominator**: The standard Elastic Net coordinate update divides by `(||X_j||^2 + lambda*(1-alpha))`. Here, the code divides `delta_j` by `norm(j)` first (inside `_soft_thresholding`), then divides the thresholded result by `(1 + lambda*(1-alpha))`. This is mathematically equivalent only if `norm(j) == 1` (i.e., columns are pre-normalized). See BUG-1 for details.

### LassoModel

Specializes ElasticNet with `alpha = 1` (pure L1 penalty). Overrides `Fit()` with a more efficient approach:

1. **Precomputed Gram matrix** (lines 64-105): computes `X'X` as a sparse matrix and `X'y` as a dense vector before iteration, avoiding repeated matrix-vector products.

2. **Sparse coordinate update** (lines 116-136): for each feature `j`:
   ```
   beta(j) = S( (X'y_j - sum_{k!=j} X'X_{j,k} * beta(k)) / ||X_j||^2, lambda*w_j )
   ```
   This iterates only over nonzero entries of `X'X` column `j`, which is efficient when `X` is sparse.

3. **Same active-set and convergence strategy** as ElasticNet.

4. **`chi2_l1()`** (line 162): returns `2 * lambda * ||beta||_1 * N`. See BUG-2 regarding the factor of 2.

### KS Test (1-sample)

`KSTest1Sample` computes the KS D-statistic by comparing the empirical CDF (ECDF) of a sorted sample set against a theoretical CDF:

1. For each sample point `x_i` (line 27-48), the ECDF steps from `i/N` to `(i+1)/N`. The D-statistic checks both `|ECDF(x_i) - F(x_i)|` and `|ECDF(x_i^-) - F(x_i)|`, taking the supremum over all points.

2. Accepts a `std::set<double>` which guarantees uniqueness and sorted order.

### KS Test (2-sample)

`KSTest2Sample` computes the D-statistic between two empirical distributions:

1. Merges both sorted sample vectors and iterates over all unique values (lines 82-106).
2. At each value, computes both ECDFs using `std::upper_bound` binary search.
3. The D-statistic is the maximum absolute difference between the two ECDFs.

### KS P-value (Asymptotic Kolmogorov Distribution)

`ks_pvalue()` (lines 120-172) computes the survival function of the Kolmogorov distribution:

1. Computes effective sample size: `n_eff = n1` (one-sample) or `n_eff = n1*n2/(n1+n2)` (two-sample).
2. Scales: `z = d * sqrt(n_eff)`.
3. Evaluates the alternating series: `P = 2 * sum_{k=1}^{inf} (-1)^{k-1} * exp(-2*k^2*z^2)`.
4. Clamps result to `[0, 1]`.

### kslike_compare

A simpler KS-like comparison (lines 215-237) for histogram/PDF vectors: normalizes both vectors to unit sum, then computes the max absolute difference between the running cumulative sums. Mimics ROOT's `TH1::KolmogorovTest` with option "M".

### discrete_cdf

Numerically integrates a PDF into a discrete CDF using trapezoidal rule (lines 175-213), then normalizes so the final value is 1.0.

## Potential Bugs

### [BUG-1] Non-standard Ridge denominator in ElasticNet coordinate update

- **File**: `util/src/ElasticNetModel.cxx`:75-76
- **Severity**: High
- **Description**: The standard Elastic Net coordinate descent update for coefficient `j` is:
  ```
  beta_j = S(X_j' * r_j, lambda*alpha*w_j) / (||X_j||^2 + lambda*(1-alpha))
  ```
  The code instead computes:
  ```
  beta_j = S(X_j' * r_j / ||X_j||^2, lambda*alpha*w_j) / (1 + lambda*(1-alpha))
  ```
  Dividing `delta_j` by `norm(j)` before soft-thresholding and then dividing by `(1 + lambda*(1-alpha))` instead of `(norm(j) + lambda*(1-alpha))` is only correct when `norm(j) == 1` (i.e., when columns of `X` are pre-normalized to unit L2 norm). If columns are not normalized, the effective L2 penalty strength becomes `lambda*(1-alpha)/norm(j)` rather than the intended `lambda*(1-alpha)`, causing features with larger norms to receive relatively weaker Ridge regularization.
- **Impact**: Regression results will be incorrect when `X` columns have varying norms and `alpha < 1`. If columns are always pre-normalized before calling `Fit()`, this is benign. However, there is no enforcement or documentation of this requirement.

### [BUG-2] Suspicious factor of 2 in chi2_l1

- **File**: `util/src/LassoModel.cxx`:162
- **Severity**: Medium
- **Description**: `chi2_l1()` returns `2 * lambda * ||beta||_1 * N`. The LASSO objective function as documented in the comment (line 11) is `(1/2)*||y-X*beta||^2 + N*lambda*||beta||_1`. The L1 penalty term is `N * lambda * ||beta||_1`, without a factor of 2. The factor of 2 may be intentional (e.g., to match the derivative or to combine with `chi2_base()` which computes the un-halved squared norm), but it is not documented and could lead to confusion. If the intent is to compute the full objective `||y-X*beta||^2 + 2*N*lambda*||beta||_1` (i.e., multiplying the objective by 2 to cancel the 1/2), the naming `chi2_l1` suggests it should return only the L1 penalty portion, making the factor of 2 unexpected.
- **Impact**: Users combining `chi2_base() + chi2_l1()` would get `||r||^2 + 2*N*lambda*||beta||_1`, which does not equal `2 * objective` (that would be `||r||^2 + 2*N*lambda*||beta||_1`, which happens to match). So the factor of 2 is consistent with `2 * objective = chi2_base() + chi2_l1()`, but this convention is nowhere documented.

### [BUG-3] Division by zero in kslike_compare for zero-sum vectors

- **File**: `util/src/KSTest.cxx`:227-228
- **Severity**: Medium
- **Description**: `kslike_compare()` computes `norm1 = 1.0/sum1` and `norm2 = 1.0/sum2` without checking whether `sum1` or `sum2` is zero. If either input vector sums to zero (e.g., all-zero entries, or positive and negative entries that cancel), this produces `inf` or `nan`, which propagates through subsequent calculations.
- **Impact**: Undefined behavior (inf/nan propagation) when called with zero-sum histogram vectors. The function documentation says inputs should be PDFs, but there is no validation.

### [BUG-4] Redundant sort in KSTest2Sample constructor and measure()

- **File**: `util/src/KSTest.cxx`:65, 75
- **Severity**: Low
- **Description**: In the `KSTest2Sample` constructor (line 65), `ref_samples_` is sorted after being copied from a `std::set<double>`, which is already sorted by definition. Similarly, in `measure()` (line 75), `test_samples_vec` is sorted after being copied from a `std::set<double>`. Both sorts are redundant.
- **Impact**: No correctness issue, but unnecessary O(N log N) work on already-sorted data. The code comment on line 65 acknowledges this: "std::set already provides this, but vector copy is explicit."

### [BUG-5] MeanResidual naming inconsistency

- **File**: `util/src/LinearModel.cxx`:14
- **Severity**: Low
- **Description**: `MeanResidual()` computes `(y - Predict()).norm() / y.size()`, which is the L2 norm (Euclidean norm) of the residual divided by N. The name "MeanResidual" suggests the mean of the residuals (i.e., `sum(|r_i|) / N` or `sum(r_i) / N`). The actual computation is `sqrt(sum(r_i^2)) / N`, which is neither the mean absolute residual nor the root-mean-square error (RMSE would be `sqrt(sum(r_i^2)/N)`). It is `||r||_2 / N`, an unusual metric.
- **Impact**: Potential confusion for users expecting standard residual statistics. The value returned is smaller than the mean absolute error by a factor related to the residual distribution shape, and smaller than RMSE by a factor of `1/sqrt(N)`.

## Efficiency Considerations

### [EFF-1] Full residual recomputation in ElasticNet coordinate descent

- **File**: `util/src/ElasticNetModel.cxx`:69-72
- **Severity**: High (performance)
- **Description**: Inside the inner loop over features `j` (line 64), the code computes the full residual `r_j = y - X * beta_tmp` at every coordinate step. This requires a full matrix-vector multiplication costing O(N*p) per feature, leading to O(N*p^2) per iteration. The standard approach is to maintain a running residual vector and perform rank-1 updates: when `beta(j)` changes by `delta`, update `residual -= delta * X.col(j)`, which costs only O(N) per feature and O(N*p) per iteration.
- **Suggestion**: Maintain a residual vector `r = y - X*beta` and update it incrementally:
  ```cpp
  double old_beta_j = beta(j);
  // ... compute new beta(j) ...
  r += (old_beta_j - beta(j)) * X.col(j);
  ```
  This is exactly what `LassoModel::Fit()` achieves via the precomputed `X'X` approach. The same optimization should be applied to `ElasticNetModel::Fit()`.

### [EFF-2] X'X sparsity estimation in LassoModel overestimates storage

- **File**: `util/src/LassoModel.cxx`:87-88
- **Severity**: Low
- **Description**: The estimated number of nonzeros is `nbeta * (nbeta / 2)`, which is O(p^2/2) -- essentially half of the dense matrix size. For truly sparse design matrices, the actual number of nonzeros in `X'X` can be much smaller. The `reserve` call on line 88 also reserves `nbeta/2` entries per column, which may over-allocate for very sparse problems or under-allocate for dense ones.

    Additionally, the triplet-based construction (lines 85-105) still computes every `X.col(i).dot(X.col(j))` for all `i,j` pairs (O(p^2 * N) work), even when many dot products will be zero. For truly sparse X, exploiting the sparsity structure of X directly (e.g., using Eigen's sparse matrix multiplication `X.transpose() * X`) would be more efficient.
- **Suggestion**: Use Eigen's built-in sparse or dense matrix multiplication for the Gram matrix computation: `XdX = X.transpose() * X`. This leverages optimized BLAS routines and avoids the manual double loop. If X is stored as a sparse matrix, the result automatically captures the sparsity pattern.

### [EFF-3] Convergence tolerance in ks_pvalue is excessively tight

- **File**: `util/src/KSTest.cxx`:146
- **Severity**: Low
- **Description**: The convergence tolerance is set to `1e-15`, which is near the limit of double-precision floating-point arithmetic (~1e-16 epsilon). For a p-value computation, this level of precision is rarely meaningful -- a tolerance of `1e-10` or even `1e-8` would be more than sufficient for statistical testing. Each unnecessary term requires an `exp()` evaluation.
- **Suggestion**: Relax the tolerance to `1e-10`. In practice, the exponential decay of terms means this has minimal impact on iteration count, but it documents the intended precision level more accurately. The `exp(-700)` guard on line 150 already provides the practical convergence bound for large z.

### [EFF-4] Data copies in ElasticNetModel::Fit() and setter methods

- **File**: `util/src/ElasticNetModel.cxx`:44-45; `util/inc/WireCellUtil/LinearModel.h`:19-26
- **Severity**: Medium
- **Description**: In `Fit()`, lines 44-45 create full copies of `_y` and `_X` via the getter methods, which return by reference but are then assigned to local values. For large datasets, this doubles memory usage unnecessarily. Similarly, the setter methods `SetData()`, `Sety()`, `SetX()`, and `Setbeta()` all take their Eigen arguments by value rather than by const reference, causing unnecessary copies on every call.
- **Suggestion**: Change local aliases to references: `Eigen::VectorXd& y = _y; Eigen::MatrixXd& X = _X;`. Change setter signatures to take `const Eigen::MatrixXd&` and `const Eigen::VectorXd&` parameters instead of by-value.

## Summary

| ID | Title | Category | Severity | File |
|----|-------|----------|----------|------|
| BUG-1 | Non-standard Ridge denominator in ElasticNet | Bug | High | ElasticNetModel.cxx:75-76 |
| BUG-2 | Suspicious factor of 2 in chi2_l1 | Bug | Medium | LassoModel.cxx:162 |
| BUG-3 | Division by zero in kslike_compare | Bug | Medium | KSTest.cxx:227-228 |
| BUG-4 | Redundant sort in KSTest2Sample | Bug | Low | KSTest.cxx:65,75 |
| BUG-5 | MeanResidual naming inconsistency | Bug | Low | LinearModel.cxx:14 |
| EFF-1 | Full residual recomputation in ElasticNet | Efficiency | High | ElasticNetModel.cxx:69-72 |
| EFF-2 | X'X sparsity estimation in Lasso | Efficiency | Low | LassoModel.cxx:87-88 |
| EFF-3 | Excessively tight convergence tolerance in ks_pvalue | Efficiency | Low | KSTest.cxx:146 |
| EFF-4 | Data copies in setters and Fit() | Efficiency | Medium | ElasticNetModel.cxx:44-45, LinearModel.h:19-26 |
