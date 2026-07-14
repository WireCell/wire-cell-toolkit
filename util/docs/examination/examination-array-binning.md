# Array & Binning - Code Examination

## Overview

This group of code provides foundational array/data structure utilities for the Wire-Cell Toolkit. It includes 2D array operations (downsample, upsample, masking, baseline subtraction, PCA), a uniform binning abstraction for discretizing linear spaces, a templated 2D object array container, and a streaming buffered 2D histogram that can dynamically expand and be sliced along the X axis.

## Files Examined

- `util/inc/WireCellUtil/Array.h`
- `util/src/Array.cxx`
- `util/inc/WireCellUtil/Binning.h`
- `util/inc/WireCellUtil/ObjectArray2d.h`
- `util/inc/WireCellUtil/BufferedHistogram2D.h`
- `util/src/BufferedHistogram2D.cxx`

## Algorithm Description

### Downsample (Array.cxx:11-28)

Reduces a 2D float array along a specified axis by a factor of `k`. A zero-initialized output array of size `floor(rows/k) x cols` (dim=0) or `rows x floor(cols/k)` (dim=1) is created. Every input row/column `i` is accumulated into output row/column `i/k` (integer division). The accumulated result is divided by `k` to produce an average. Excess rows/columns beyond the last full group of `k` are included in the accumulation but the output array is only sized for complete groups, causing out-of-bounds access for the remainder.

### Upsample (Array.cxx:30-47)

Expands a 2D float array along a specified axis by a factor of `k`. A zero-initialized output of size `rows*k x cols` (dim=0) or `rows x cols*k` (dim=1) is created. Each output row/column `i` copies the content of input row/column `i/k`, effectively repeating each input element `k` times.

### Threshold Mask (Array.cxx:49-57)

Applies an element-wise mask to an input array. Both arrays must have identical shape. Elements where the mask value exceeds the threshold `th` are kept; all others are zeroed. Uses Eigen's `.select()` for efficient conditional assignment.

### Linear Baseline Subtraction (Array.cxx:59-81)

Operates column-by-column. Scans each column for contiguous regions of interest (ROI) delimited by zero-valued samples. For each ROI bounded by `[sta, end]`, a linear interpolation between the values at `sta` and `end` is subtracted from each sample within the ROI. The algorithm uses exact floating-point equality (`== 0`) to detect ROI boundaries.

### PCA (Array.h:92-97)

Performs Principal Component Analysis on a matrix where rows are observations and columns are features. The data is mean-centered per column, a covariance matrix is computed with Bessel's correction `1/(n-1)`, and Eigen's `SelfAdjointEigenSolver` extracts eigenvalues (ascending order) and corresponding eigenvectors (as columns).

### Uniform Binning (Binning.h:15-103)

Provides a uniform 1D binning abstraction: `nbins` equally-spaced bins over `[minval, maxval)`. Key operations include bin lookup (`bin()`), center/edge calculation, range membership testing (`inside()`), sub-range extraction (`subset()`), and Gaussian bin-integrated distribution calculation (`gaussian()`).

### ObjectArray2d (ObjectArray2d.h:14-76)

A simple templated 2D container backed by a flat `std::vector<Thing>`. Elements are stored in row-major order and accessed via `operator()(irow, icol)` with bounds checking (`at()`). Provides copy construction/assignment and range-based iteration.

### BufferedHistogram2D (BufferedHistogram2D.h/cxx)

A dynamically-expandable 2D histogram backed by a `deque<vector<double>>`. The X axis grows by appending entries to the deque; each X-bin's Y axis grows independently via `vector::resize()`. The `fill()` method auto-expands both axes. `popx()` removes and returns the lowest X-bin, advancing `xmin`. The `ysize()` method scans all X-bins to find the maximum Y-bin count.

## Potential Bugs

### [BUG-1] Out-of-bounds write in downsample() when input size is not divisible by k

- **File**: `util/src/Array.cxx:15-16` (dim=0), `util/src/Array.cxx:22-23` (dim=1)
- **Severity**: Critical
- **Description**: The output array is sized to `in.rows() / k` (or `in.cols() / k`), which is `floor(N/k)`. However, the loop iterates over all input rows `i` from `0` to `in.rows()-1`. For the remainder rows where `i >= floor(N/k) * k`, the expression `i / k` yields `floor(N/k)`, which equals the output array's row count and is therefore out of bounds. For example, with 10 rows and k=3, the output has 3 rows (indices 0-2), but when `i=9`, `i/k = 3`, accessing `out.row(3)` which does not exist.
- **Impact**: Undefined behavior -- memory corruption or crash. Any input whose row/column count is not an exact multiple of `k` will trigger this bug.

### [BUG-2] Exact float equality in baseline_subtraction() ROI boundary detection

- **File**: `util/src/Array.cxx:66`
- **Severity**: Medium
- **Description**: The condition `in(it, ich) == 0` uses exact floating-point equality to detect zero-valued samples that delimit regions of interest. Due to floating-point arithmetic, values that are conceptually zero (e.g., `1e-16` from prior operations) will not match, causing ROI boundaries to be missed. Conversely, legitimate near-zero signal values may be incorrectly treated as boundaries.
- **Impact**: Incorrect ROI detection leading to wrong baseline subtraction results. The severity depends on upstream data precision; if the input is known to contain exact zeros (e.g., from integer conversion), this is acceptable. Otherwise, an epsilon-based comparison would be more robust.

### [BUG-3] Division by zero in baseline_subtraction() when sta == end

- **File**: `util/src/Array.cxx:70`
- **Severity**: Medium
- **Description**: When an ROI consists of a single non-zero sample (i.e., a lone non-zero value between two zeros), `sta` and `end` will both point to the same index. The guard `sta < end` (line 67) prevents the subtraction loop from executing in this case, so the single-sample ROI is silently dropped (its output remains zero). While this avoids an actual division by zero, the single non-zero sample is lost rather than preserved or baseline-corrected, which may not be the intended behavior.
- **Impact**: Single-sample ROIs are zeroed out in the output. If the intent is to preserve them, this is a data loss bug.

### [BUG-4] ysize() iterates by value causing unnecessary copies and potentially incorrect results

- **File**: `util/src/BufferedHistogram2D.cxx:19`
- **Severity**: Low
- **Description**: The loop `for (auto v : m_xbindeque)` copies each `XBin` (a `std::vector<double>`) by value. While this does not produce incorrect results per se, it is wasteful. The `ysize()` method itself returns the maximum Y-bin vector length across all X-bins. If `m_xbindeque` is empty, it returns 0, which is correct. However, `ymax()` (line 15) depends on `ysize()`, and when no data has been filled, `ymax()` returns `ymin()`, which may be unexpected for callers expecting a meaningful range.
- **Impact**: Semantic issue: `ymax() == ymin()` when the histogram is empty may confuse callers. Performance is addressed in EFF-1.

### [BUG-5] Binning::bin() returns negative values for inputs below minval

- **File**: `util/inc/WireCellUtil/Binning.h:75`
- **Severity**: Low
- **Description**: The `bin()` method computes `int((val - m_minval) / m_binsize)` with no range checking. For `val < m_minval`, this returns a negative bin index. The documentation notes "no range checking is performed," so this is by design, but callers that use `bin()` without checking `inside()` first could produce out-of-bounds array access. The `sample_bin_range()` method (line 101) does clamp against 0 and `m_nbins`, providing a safe alternative.
- **Impact**: Callers must be aware that `bin()` can return out-of-range values. Not a bug in the class itself but a potential source of bugs in calling code.

### [BUG-6] ObjectArray2d copy constructor does not handle self-assignment

- **File**: `util/inc/WireCellUtil/ObjectArray2d.h:20-23`
- **Severity**: Low
- **Description**: The copy constructor delegates to `operator=`, and `operator=` (line 25) clears `m_things` before copying from `other`. If somehow self-assignment occurs (unlikely for the copy constructor but possible for `operator=`), the data would be destroyed before being read. There is no self-assignment guard (`if (this == &other) return *this;`).
- **Impact**: Self-assignment of an `ObjectArray2d` instance would result in data loss. This is an unusual usage pattern but a standard C++ correctness issue.

### [BUG-7] popx() returns by value instead of by move

- **File**: `util/src/BufferedHistogram2D.cxx:48-59`
- **Severity**: Low
- **Description**: `popx()` copies the front vector into `ret` on line 56 (`ret = m_xbindeque.front()`), then pops the front. Since `front()` returns a reference to an element about to be destroyed, a move would be more appropriate: `ret = std::move(m_xbindeque.front())`. Most modern compilers may optimize this, but the copy is semantically unnecessary.
- **Impact**: Performance overhead for large Y-bin vectors. Addressed further in EFF-2.

## Efficiency Considerations

### [EFF-1] ysize() is O(n) and copies every vector

- **File**: `util/src/BufferedHistogram2D.cxx:16-23`
- **Description**: `ysize()` iterates over all X-bins to find the maximum Y-bin count. The loop variable `auto v` copies each `std::vector<double>` by value. Even with the correct `const auto&`, this remains O(n) in the number of X-bins. Since `ysize()` is called by `ymax()`, and these may be called frequently, this linear scan is inefficient.
- **Suggestion**: Use `const auto& v` to eliminate copies. Consider caching the maximum Y size as a member variable, updated in `fill()` and `popx()`, to make `ysize()` O(1).

### [EFF-2] popx() copies vector by value instead of moving

- **File**: `util/src/BufferedHistogram2D.cxx:56`
- **Description**: The assignment `ret = m_xbindeque.front()` copies the front vector. Since the front is immediately popped and destroyed on the next line, the copy is wasted work.
- **Suggestion**: Use `ret = std::move(m_xbindeque.front())` to transfer ownership without copying. This avoids allocation and element-wise copy for potentially large vectors.

### [EFF-3] ObjectArray2d lacks move constructor and move assignment

- **File**: `util/inc/WireCellUtil/ObjectArray2d.h:14-76`
- **Description**: The class defines a copy constructor and copy assignment operator but does not provide move equivalents. For types `Thing` that are expensive to copy (e.g., `std::vector`, `std::string`), moving an `ObjectArray2d` will fall back to copying, which is unnecessarily expensive.
- **Suggestion**: Add a move constructor and move assignment operator, or simplify by following the Rule of Zero (remove the custom copy operations and let the compiler generate all five special member functions, since the class only holds `size_t` and `std::vector`).

### [EFF-4] ObjectArray2d::reset() uses emplace_back in a loop instead of resize

- **File**: `util/inc/WireCellUtil/ObjectArray2d.h:50-57`
- **Description**: The `reset()` method clears the vector, reserves capacity, and then calls `emplace_back()` in a loop to default-construct `m_nrows * m_ncols` elements. This could be simplified to a single `m_things.resize(m_nrows * m_ncols)` call, which default-constructs all elements in one operation and may be more efficient.
- **Suggestion**: Replace the clear/reserve/emplace_back loop with `m_things.clear(); m_things.resize(m_nrows * m_ncols);` or simply `m_things.assign(m_nrows * m_ncols, Thing{});`.

### [EFF-5] downsample() iterates all input rows/cols instead of only full groups

- **File**: `util/src/Array.cxx:15-16`
- **Description**: Even after fixing BUG-1, the loop iterates over all input rows/columns. The documentation says "extra rows/cols are ignored," so the loop should only iterate up to `(in.rows() / k) * k` to skip the remainder elements, avoiding unnecessary additions to the last output row that would then be divided out.
- **Suggestion**: Change the loop bound from `in.rows()` to `(in.rows() / k) * k` (and similarly for dim=1), which both fixes BUG-1 and avoids wasted computation on remainder elements.

## Summary

| ID | Type | Severity | File | Title |
|----|------|----------|------|-------|
| BUG-1 | Bug | Critical | Array.cxx:15-16,22-23 | Out-of-bounds write in downsample() when size not divisible by k |
| BUG-2 | Bug | Medium | Array.cxx:66 | Exact float equality in baseline_subtraction() ROI detection |
| BUG-3 | Bug | Medium | Array.cxx:70 | Single-sample ROIs silently dropped in baseline_subtraction() |
| BUG-4 | Bug | Low | BufferedHistogram2D.cxx:19 | ysize() copies vectors by value; empty histogram edge case |
| BUG-5 | Bug | Low | Binning.h:75 | bin() returns negative indices without range checking |
| BUG-6 | Bug | Low | ObjectArray2d.h:25 | No self-assignment guard in operator= |
| BUG-7 | Bug | Low | BufferedHistogram2D.cxx:56 | popx() copies front vector instead of moving |
| EFF-1 | Efficiency | Medium | BufferedHistogram2D.cxx:16-23 | ysize() is O(n) with per-element vector copies |
| EFF-2 | Efficiency | Low | BufferedHistogram2D.cxx:56 | popx() copies instead of moves |
| EFF-3 | Efficiency | Low | ObjectArray2d.h:14-76 | Missing move constructor and move assignment |
| EFF-4 | Efficiency | Low | ObjectArray2d.h:50-57 | reset() loop could be simplified to resize() |
| EFF-5 | Efficiency | Low | Array.cxx:15-16 | downsample() iterates remainder elements unnecessarily |
