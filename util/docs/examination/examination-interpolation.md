# Interpolation - Code Examination

## Overview

This file provides linear interpolation utilities for the Wire-Cell Toolkit. It includes two main interpolator classes (`irrterp` for irregularly-spaced samples and `linterp` for regularly-spaced samples) plus a convenience function (`linterpolate`) for one-shot resampling. Both classes support single-point evaluation and batch sequence output. The irregular interpolator uses a `std::map` for O(log n) point lookup with clamped boundary behavior; the regular interpolator uses a `std::vector` with O(1) arithmetic index computation.

## Files Examined

- `util/inc/WireCellUtil/Interpolate.h`

## Algorithm Description

### irrterp - Irregular Linear Interpolation

Stores sample points in a `std::map<X, Y>` which keeps them sorted by X coordinate. Single-point evaluation (line 115-132) uses `upper_bound` for O(log n) lookup to find the bracketing interval, then performs linear interpolation. Values outside the domain are clamped to the nearest boundary value (lines 120-125).

Batch sequence output (lines 137-200) iterates through `num` output points at regular spacing (`xstart`, `xstep`). It maintains an iterator `ub` that advances forward through the map as `xcur` increases, avoiding repeated O(log n) lookups. This gives amortized O(1) per output point when stepping forward. Three regions are handled: before the domain (clamp to first value), after the domain (clamp to last value), and within the domain (linear interpolation via `std::lerp`).

### linterp - Regular Linear Interpolation

Stores sample values in a `std::vector<Y>` with a starting X value (`m_le`), ending X value (`m_re`), and step size (`m_step`). Single-point evaluation (lines 313-323) computes the integer bin index via arithmetic: `ind = int((x - m_le) / m_step)`, then interpolates between `m_dat[ind]` and `m_dat[ind+1]` using `std::lerp`. Values outside the domain are clamped (lines 315-316). This is O(1) per point.

Batch sequence output (lines 325-345) simply loops and calls the single-point operator for each output point. No iterator optimization is performed.

### linterpolate - One-Shot Resampling

A free function (lines 350-361) that resamples an input range to an output range by mapping both to the normalized domain [0, 1]. It constructs a temporary `linterp` object and calls its batch operator.

### std::lerp Polyfill

Lines 53-60 provide a naive `std::lerp` implementation for pre-C++20 compilers. It uses the formula `a + t * (b - a)`.

## Potential Bugs

### [BUG-1] Compilation error in map copy constructor

- **File**: `util/inc/WireCellUtil/Interpolate.h`:104
- **Severity**: Critical
- **Description**: Line 104 reads `points(pts.begin(), pts.end)` -- missing parentheses on `pts.end`. This should be `pts.end()`. This is a function call, not a member pointer. This code will fail to compile if this constructor overload is ever instantiated.
- **Impact**: Any code attempting to construct an `irrterp` from a `const std::map<X,Y>&` via this constructor will produce a compilation error. The other constructors (iterator-range, map-move) still work, so this may be latent.

### [BUG-2] Type truncation when Y is an integer type in linterp

- **File**: `util/inc/WireCellUtil/Interpolate.h`:321
- **Severity**: High
- **Description**: The `delta` variable is declared as type `Y` (line 321): `Y delta = (x-x0)/m_step;`. The expression `(x-x0)/m_step` produces a fractional value of type `X`, but if `Y` is an integer type (e.g., `int`), the fractional part is truncated. The subsequent `std::lerp(m_dat[ind], m_dat[ind+1], delta)` would then receive an integer `delta` that is always 0 (since `0 <= (x-x0)/m_step < 1` within a bin), meaning interpolation effectively becomes a floor/step function.
- **Impact**: When `linterp<double, int>` or similar integer-Y instantiations are used, interpolation returns only the left sample value with no blending. The same bug exists in `irrterp` at line 188.

### [BUG-3] Undefined behavior when m_dat is empty in linterp constructor

- **File**: `util/inc/WireCellUtil/Interpolate.h`:290
- **Severity**: High
- **Description**: The constructor at line 285-291 computes `m_re = m_le + m_step * (m_dat.size() - 1)`. If the iterator range `[beg, end)` is empty, `m_dat.size()` is 0 and `m_dat.size() - 1` wraps around to `SIZE_MAX` (unsigned underflow), producing a nonsensical `m_re`. Furthermore, any subsequent call to `operator()` will call `m_dat.front()` or `m_dat.back()` on an empty vector, which is undefined behavior.
- **Impact**: Constructing an `linterp` from an empty range causes silent corruption of `m_re` and eventual undefined behavior on evaluation. Unlike `irrterp`, which checks for emptiness and throws, `linterp` has no such guard.

### [BUG-4] Potential out-of-bounds access in linterp near right edge

- **File**: `util/inc/WireCellUtil/Interpolate.h`:318-322
- **Severity**: Medium
- **Description**: The clamping guard at line 316 is `if (x >= m_re) return m_dat.back();`. For values just below `m_re`, the index computation `int ind = int((x - m_le) / m_step)` could yield `ind == m_dat.size() - 1` due to floating-point rounding. The subsequent access `m_dat[ind+1]` at line 322 would then read one past the last element. This can occur when `x` is very close to but slightly less than `m_re`, and floating-point arithmetic rounds `(x - m_le) / m_step` up to exactly `m_dat.size() - 1`.
- **Impact**: Out-of-bounds read on the data vector. In practice this is rare due to the `>=` guard, but floating-point edge cases near `m_re` can trigger it.

### [BUG-5] irrterp batch mode assumes forward-stepping only

- **File**: `util/inc/WireCellUtil/Interpolate.h`:159-198
- **Severity**: Medium
- **Description**: The batch operator advances the `ub` iterator only forward (lines 167-173, 192-197) using `while (xcur > ub->first) { ++ub; }`. If `xstep` is negative (reverse traversal), `xcur` decreases each iteration but `ub` never moves backward. After the first point, the iterator will remain stuck at its current position, and subsequent points will be interpolated using the wrong interval or will incorrectly appear to be "before the domain."
- **Impact**: Batch interpolation with negative step produces incorrect results. Single-point `operator()` works correctly for any x value since it uses `upper_bound` each time.

### [BUG-6] Off-by-one in linterpolate resampling grid

- **File**: `util/inc/WireCellUtil/Interpolate.h`:356-360
- **Severity**: Low
- **Description**: The input step is computed as `olddx = (xmax-xmin)/ilen` (line 356) and the output step as `newdx = (xmax-xmin)/olen` (line 359). The `linterp` constructor interprets its data as `ilen` points covering `xmin` to `xmin + olddx * (ilen-1)`. So the effective domain right edge is `xmin + (xmax-xmin) * (ilen-1)/ilen`, which is slightly less than `xmax`. Meanwhile, the output may request points up to `xmin + newdx * (olen-1)`, which for large `olen` approaches `xmax`. Points beyond the effective domain are clamped to the last sample value rather than being interpolated, causing a small plateau at the end of the resampled output.
- **Impact**: The last few output samples may be clamped rather than interpolated, producing a subtle inaccuracy in resampled waveforms. Using `(xmax-xmin)/(ilen-1)` and `(xmax-xmin)/(olen-1)` as step sizes would map endpoints correctly.

## Efficiency Considerations

### [EFF-1] linterp batch mode does not optimize iterator advancement

- **File**: `util/inc/WireCellUtil/Interpolate.h`:339-343
- **Description**: The `linterp` batch operator (lines 325-345) simply calls the single-point `operator()` in a loop. While each call is O(1), it redundantly recomputes the clamping checks and index each time. By contrast, `irrterp`'s batch mode (lines 137-200) maintains an advancing iterator to avoid repeated O(log n) lookups.
- **Suggestion**: For `linterp`, the overhead per call is small (arithmetic, two comparisons), so this is a minor concern. However, for hot paths a specialized loop that tracks the current bin index and advances it could eliminate the repeated division and branching.

### [EFF-2] Redundant copy/move constructors and assignment operators

- **File**: `util/inc/WireCellUtil/Interpolate.h`:83-101, 237-282
- **Description**: Both `irrterp` and `linterp` explicitly define copy constructors, move constructors, copy assignment, and move assignment operators that do exactly what the compiler-generated defaults would do. These add code bulk with no behavioral benefit.
- **Suggestion**: Remove the explicit definitions and rely on the compiler-generated defaults (Rule of Zero), or use `= default` declarations. The move constructor for `linterp` (lines 244-255) also unnecessarily zeroes out the moved-from object's primitives, which is not required by the standard.

### [EFF-3] irrterp batch mode copies the map key/value on each interpolation

- **File**: `util/inc/WireCellUtil/Interpolate.h`:185-189
- **Description**: In the inner loop of the batch operator, `lb` and `ub` iterators are dereferenced to access `first` and `second` members. For simple numeric types this is fine, but the pattern of `auto lb = ub; --lb;` creates a new iterator copy each iteration.
- **Suggestion**: This is a minor concern for numeric types. The current implementation is appropriate for the typical use case.

## Summary

| ID | Title | Severity | Type |
|----|-------|----------|------|
| BUG-1 | Compilation error in map copy constructor (`pts.end` missing parens) | Critical | Bug |
| BUG-2 | Type truncation when Y is integer type (delta cast to Y) | High | Bug |
| BUG-3 | Undefined behavior when m_dat is empty in linterp constructor | High | Bug |
| BUG-4 | Potential out-of-bounds in linterp near right edge | Medium | Bug |
| BUG-5 | irrterp batch mode assumes forward-stepping only | Medium | Bug |
| BUG-6 | Off-by-one in linterpolate resampling grid | Low | Bug |
| EFF-1 | linterp batch mode does not optimize iterator advancement | Low | Efficiency |
| EFF-2 | Redundant copy/move constructors and assignment operators | Low | Efficiency |
| EFF-3 | irrterp batch mode iterator copy per iteration | Low | Efficiency |
