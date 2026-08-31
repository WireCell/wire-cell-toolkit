# Response Functions - Code Examination

## Overview

This code provides electronics response functions and field response data structures for the Wire-Cell Toolkit signal processing pipeline. It models the shaping behavior of cold electronics (DUNE), parameterized cold electronics (with tunable undershoot/overshoot), warm electronics (ICARUS), and provides utilities for loading field response data from JSON, averaging over wire regions, and applying frequency-domain filters (high-pass, low-pass). These response functions are convolved with field responses to simulate or deconvolve detector signals.

## Files Examined

- `util/inc/WireCellUtil/Response.h`
- `util/src/Response.cxx`

## Algorithm Description

### Cold Electronics Response (`coldelec`)

The `coldelec` function (line 311-335, Response.cxx) computes the time-domain impulse response of the DUNE cold electronics front-end ASIC. It is derived from the inverse Laplace transform of the electronics transfer function (provided by Hucheng). The function:

1. Returns 0 for `time <= 0` or `time >= 10 us` (validity range).
2. Normalizes time by the shaping time: `reltime = time / shaping`.
3. Applies a gain scaling factor of `10 * 1.012` (noted as slightly shaping-dependent).
4. Evaluates a sum of 12 terms, each a product of an exponential decay and trigonometric functions with hardcoded coefficients (from Mathematica symbolic inversion). The expression has the general form: `sum_i A_i * exp(-alpha_i * reltime) * [cos/sin terms]`.

The five poles of the transfer function produce three exponential decay rates (2.94809, 2.82833, 2.40318) and two oscillation frequencies (1.19361, 2.5928) in normalized units.

### Parameterized Cold Electronics Response (`paramscoldelec`)

The `paramscoldelec` function (lines 341-432, Response.cxx) is a generalization of `coldelec` with four additional parameters (`k3`, `k4`, `k5`, `k6`) controlling undershoot and overshoot behavior. The algorithm:

1. Returns 0 for `time <= 0` or `time >= 50 us`.
2. Converts time to microseconds and computes derived pole/residue parameters from the shaping time via calibration constants (`CT = 1/1.996`, `p0 = 1.477/tp/CT`, etc.).
3. Evaluates a massive closed-form expression (~70 lines of dense algebra) representing the partial-fraction decomposition of a transfer function with poles at `p0`, `pr1 +/- i*pi1`, `pr2 +/- i*pi2`, `k4`, and `k6`. The expression has three pure-exponential terms plus two damped-sinusoidal terms (cos and sin at each of pi1 and pi2 frequencies).

The formula is the result of symbolic computation (likely Mathematica) and represents the exact inverse Laplace transform of a rational transfer function with the given poles and two additional real poles/zeros at k4 and k6.

### Warm Electronics Response (`warmelec`)

The `warmelec` function (lines 494-506, Response.cxx) models the ICARUS warm electronics using a Bessel-function-inspired approximation from arXiv:1805.03931. The algorithm:

1. Returns 0 for `time <= 0` or `time >= 10 us`.
2. Computes `reltime = time / shaping`.
3. Returns `gain * (1 - exp(-0.5 * (0.9*reltime)^2)) * exp(-0.5 * (0.5*reltime)^2)`.

This produces a pulse shape with a Gaussian-like rise (controlled by the left factor) and Gaussian decay (right factor). The factors 0.9 and 0.5 are empirical fits to match SPICE simulation and CERN test data.

### High-Frequency Filter (`hf_filter`)

Lines 435-442: A Gaussian-type filter `exp(-0.5 * (freq/sigma)^power)` with optional zero-frequency removal (`flag`). When `flag` is true and `freq == 0`, returns 0.

### Low-Frequency Filter (`lf_filter`)

Line 444: A complementary filter `1 - exp(-(freq/tau)^2)` that suppresses low frequencies and passes high frequencies.

### Field Response Loading (`Schema::load`)

Lines 31-71: Loads field response data from a JSON file via the `Persist::load` utility. Parses a nested JSON structure (`FieldResponse -> planes[] -> PlaneResponse -> paths[] -> PathResponse -> current/pitchpos/wirepos`) into C++ schema objects. Returns a default-constructed `FieldResponse` on failure.

### Wire Region Averaging (`wire_region_average`)

Lines 78-195: Averages path responses within each wire region. The algorithm:

1. For each plane, maps each path to an integer key `eff_num = pitchpos / (0.01 * pitch)` and mirrors it to negative positions (exploiting symmetry).
2. When multiple paths map to the same key, averages them pairwise.
3. Computes pitch-position ranges for each response by taking midpoints between adjacent positions.
4. Determines wire regions by rounding pitch positions to nearest wire number.
5. For each wire region, integrates each response contribution by multiplying by the overlap fraction `(high_limit - low_limit) / pitch`.

### 1D Averaging (`average_1D`)

Lines 197-224: First calls `wire_region_average`, then sums all wire-region responses into a single aggregate response per plane. This produces a single 1D response function per plane.

### SimpleRC Response

Lines 517-536: Models a simple RC circuit response: a delta function at the offset time plus an exponential decay `-tick/width * exp(-(time - offset)/width)`. The tick factor normalizes the integral per time bin.

### SysResp (Systematic Response)

Lines 540-561: A Gaussian smearing function for field response systematics studies. When `smear > 0`, returns a normalized Gaussian `tick * exp(-0.5*((time-offset)/smear)^2) / smear * 0.3989422804` (where 0.3989422804 = 1/sqrt(2*pi)). Otherwise returns a rectangular pulse of width `tick`.

### `as_array`

Lines 226-262: Converts a `PlaneResponse` into an Eigen 2D array where each row is a path's current waveform. Two overloads: one auto-sized, one with explicit dimensions (zero-padded).

## Potential Bugs

### [BUG-1] Division by zero in `paramscoldelec` when poles coincide
- **File**: `util/src/Response.cxx`:359-361
- **Severity**: High
- **Description**: The denominator of the first three terms of the `paramscoldelec` expression contains factors like `(k4 - k6)`, `(k4 - p0)`, `(k6 - p0)`, and quadratic terms `(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))`. If `k4 == k6`, `k4 == p0`, or `k6 == p0`, a division by zero occurs. Similarly, the sinusoidal terms (lines 392-393, 429-430) divide by `pi1` and `pi2`; if the imaginary parts of the poles are zero, this produces division by zero.
- **Impact**: Undefined behavior (infinity or NaN) propagating through the signal processing chain. With default parameters (`k3=0.1, k4=0.1, k5=0.0, k6=0.0`), `k6=0.0` could coincide with `p0` if `tp*CT` yields `p0=0`, though with typical shaping times this is unlikely. More dangerous if users pass arbitrary parameter values.

### [BUG-2] Floating-point equality comparison in `hf_filter`
- **File**: `util/src/Response.cxx`:438
- **Severity**: Low
- **Description**: The condition `if (freq == 0)` uses exact floating-point equality. If `freq` is a very small but nonzero number (e.g., due to floating-point arithmetic), the zero-removal will not trigger, and the filter will return a value near 1.0 instead of 0.0.
- **Impact**: The DC component may not be properly removed. In practice, the frequency is often computed from an FFT bin index where bin 0 is exactly 0.0, so this may work correctly in typical usage.

### [BUG-3] `nsamples` may be zero in `wire_region_average`
- **File**: `util/src/Response.cxx`:97-108
- **Severity**: Medium
- **Description**: The variable `nsamples` is initialized to 0 (line 97) and only set inside the `else` branch (line 105) when a duplicate `eff_num` is found. If all paths map to unique `eff_num` values (no duplicates), `nsamples` remains 0, and the averaging loop at line 177 (`for (int k = 0; k != nsamples; k++)`) will do nothing, producing all-zero averaged responses. Also, `avgs[wire_no] = realseq_t(nsamples)` at line 160 would create zero-length vectors.
- **Impact**: Silent production of empty/zero response arrays if the input has unique pitch positions (no mirroring collisions).

### [BUG-4] Incorrect pairwise averaging for duplicate paths
- **File**: `util/src/Response.cxx`:106-108
- **Severity**: Medium
- **Description**: When a duplicate `eff_num` is found, the code averages the existing value with the new one: `(existing + new) / 2`. If there are three or more paths mapping to the same `eff_num`, each subsequent path is averaged with the running result, giving unequal weights. For example, with values A, B, C: result = ((A+B)/2 + C)/2 = A/4 + B/4 + C/2, instead of the correct (A+B+C)/3.
- **Impact**: Incorrect field response weighting when more than two paths share the same discretized pitch position.

### [BUG-5] Variable shadowing of `nsamples`
- **File**: `util/src/Response.cxx`:114
- **Severity**: Low
- **Description**: Inside the negative-mirror `else` block (line 114), `int nsamples = path.current.size()` declares a new local variable that shadows the outer `nsamples` declared at line 97. The outer `nsamples` is only updated in the positive-key else branch (line 105). If the first duplicate is found only in the negative-mirror branch, the outer `nsamples` remains 0.
- **Impact**: Same as BUG-3 -- the outer `nsamples` may not be set, leading to zero-length averaging.

### [BUG-6] `as_array` crashes on empty `paths` vector
- **File**: `util/src/Response.cxx`:229, 252
- **Severity**: Medium
- **Description**: Both `as_array` overloads access `pr.paths[0].current.size()` without checking if `pr.paths` is empty. If a `PlaneResponse` has no paths, this is undefined behavior (out-of-bounds access).
- **Impact**: Segmentation fault if called with an empty plane response.

### [BUG-7] Redundant bounds check in `as_array` with explicit dimensions
- **File**: `util/src/Response.cxx`:237-238
- **Severity**: Low
- **Description**: The check `if (irow < set_nrows)` at line 238 is always true because the loop runs `irow < nrows` and we already verified `set_nrows >= nrows` at line 232. Similarly, `if (icol < set_ncols)` at line 241 is always true. Not a bug per se, but dead code indicating possible confusion about the intended logic.
- **Impact**: No runtime impact; code clarity issue only.

### [BUG-8] `dump` function is a no-op
- **File**: `util/src/Response.cxx`:73
- **Severity**: Low
- **Description**: The `dump` function has an empty body `{}`. It is declared in the header (line 149) but never implemented, silently doing nothing when called.
- **Impact**: Users expecting to serialize a `FieldResponse` to a file will get no output and no error.

### [BUG-9] SimpleRC response not properly bounded for negative times
- **File**: `util/src/Response.cxx`:525-536
- **Severity**: Low
- **Description**: The `SimpleRC::operator()` does not check for `time < offset`. For `time < offset`, the exponential term `exp(-(time - offset) / _width)` becomes `exp(positive)`, which grows exponentially. Combined with the negative sign, this produces large negative values for times well before the offset. The delta-function term (line 532-534) also triggers for `time < _offset + _tick`, which includes all negative times when `_offset = 0`.
- **Impact**: Unphysical exponentially growing response for times before the impulse, if the caller does not restrict the time range. The `Generator::generate` functions typically start from `domain.first` or `tbins.center(0)` which may include negative times.

### [BUG-10] `coldelec` gain scaling comment says "slightly dependent on shaping time"
- **File**: `util/src/Response.cxx`:320-321
- **Severity**: Low
- **Description**: The comment on line 319-320 acknowledges that the scaling factor `10 * 1.012` is slightly dependent on shaping time but uses a fixed constant regardless. For non-default shaping times, the peak gain will not exactly match the requested gain.
- **Impact**: Small gain calibration error (likely < 1%) for non-default shaping times.

## Efficiency Considerations

### [EFF-1] Triple-nested loop in `wire_region_average`
- **File**: `util/src/Response.cxx`:157-182
- **Severity**: Medium
- **Description**: The averaging loop iterates over `wire_regions x fresp_map x nsamples`. For each wire region, it scans all response entries to find overlapping ones. With W wire regions, R response entries, and N samples, this is O(W * R * N).
- **Suggestion**: Pre-compute which response entries overlap each wire region to avoid iterating all entries for every wire. Alternatively, since the pitch positions are sorted (from `std::map`), use binary search to find the overlapping range for each wire region.

### [EFF-2] Repeated map lookups in `wire_region_average`
- **File**: `util/src/Response.cxx`:100-118
- **Severity**: Low
- **Description**: The code calls `fresp_map.find(eff_num)` and `fresp_map.find(-eff_num)` separately, then accesses `fresp_map[eff_num]` again inside the else branch (via `.at(k)`). Each `find` + `operator[]` is two separate tree traversals. Lines 165-166 also perform separate lookups for `pitch_pos_range_map[resp_num]`.
- **Suggestion**: Use the iterator returned by `find()` instead of re-looking up the key. Use `auto [it, inserted] = fresp_map.emplace(...)` pattern to combine lookup and insertion.

### [EFF-3] Element-by-element array copy in `as_array`
- **File**: `util/src/Response.cxx`:240-242, 257-258
- **Severity**: Medium
- **Description**: The `as_array` functions copy current waveforms into an Eigen array one element at a time. The code even has a comment asking "maybe there is a fast way to do this copy?" (line 242, 258).
- **Suggestion**: Use `Eigen::Map` to wrap the `std::vector` data and assign entire rows at once: `ret.row(irow).head(ncols) = Eigen::Map<const Eigen::VectorXf>(path.current.data(), ncols);` (assuming compatible scalar types).

### [EFF-4] Copying large structures by value in range-for loops
- **File**: `util/src/Response.cxx`:47, 84, 99, 206, 213
- **Severity**: Medium
- **Description**: Several range-for loops copy `PlaneResponse` and `PathResponse` objects by value: `for (auto plane : fr.planes)` (line 84), `for (auto path : plane.paths)` (lines 99, 213), `for (auto plane : fr["planes"])` (line 47). Each `PathResponse` contains a `realseq_t` (vector) that gets deep-copied.
- **Suggestion**: Use `const auto&` references: `for (const auto& plane : fr.planes)`, `for (const auto& path : plane.paths)`.

### [EFF-5] `push_back` in JSON loading loop without `reserve`
- **File**: `util/src/Response.cxx`:54-56
- **Severity**: Low
- **Description**: The `current` vector is built by `push_back` in a loop over JSON elements (line 55) without a prior `reserve()`. The JSON array shape information (`par["current"]["array"]["shape"]`) is available and could be used to pre-allocate.
- **Suggestion**: Read the shape and call `current.reserve(expected_size)` before the loop, or use `std::transform` with a pre-sized vector.

### [EFF-6] Repeated `pow()` calls in `paramscoldelec`
- **File**: `util/src/Response.cxx`:359-430
- **Severity**: Low
- **Description**: The massive expression in `paramscoldelec` repeatedly computes `pow(pi1,2)`, `pow(pi2,2)`, `pow(pr1,2)`, `pow(pr2,2)`, `pow(k4,2)`, `pow(k6,2)`, `pow(p0,2)`, etc. -- hundreds of times. While the compiler may optimize some of these via CSE (common subexpression elimination), it is not guaranteed, especially at lower optimization levels.
- **Suggestion**: Pre-compute squared values (`double pi1_2 = pi1*pi1;`, etc.) and intermediate subexpressions before the main formula. This would also improve readability.

### [EFF-7] `average_1D` creates temporary `as_array` just to get column count
- **File**: `util/src/Response.cxx`:209
- **Severity**: Low
- **Description**: `int nsamples = Response::as_array(plane).cols()` constructs a full 2D Eigen array just to determine the number of columns (samples). This allocates and fills an entire matrix only to read one dimension.
- **Suggestion**: Use `plane.paths[0].current.size()` directly (with an empty-check), avoiding the full array construction.

## Summary

| ID | Type | Severity | Title |
|----|------|----------|-------|
| BUG-1 | Bug | High | Division by zero in `paramscoldelec` when poles coincide |
| BUG-2 | Bug | Low | Floating-point equality comparison in `hf_filter` |
| BUG-3 | Bug | Medium | `nsamples` may be zero in `wire_region_average` |
| BUG-4 | Bug | Medium | Incorrect pairwise averaging for duplicate paths |
| BUG-5 | Bug | Low | Variable shadowing of `nsamples` |
| BUG-6 | Bug | Medium | `as_array` crashes on empty `paths` vector |
| BUG-7 | Bug | Low | Redundant bounds check in `as_array` |
| BUG-8 | Bug | Low | `dump` function is a no-op |
| BUG-9 | Bug | Low | SimpleRC response not bounded for negative times |
| BUG-10 | Bug | Low | Fixed gain scaling factor is shaping-time dependent |
| EFF-1 | Efficiency | Medium | Triple-nested loop in `wire_region_average` |
| EFF-2 | Efficiency | Low | Repeated map lookups in `wire_region_average` |
| EFF-3 | Efficiency | Medium | Element-by-element array copy in `as_array` |
| EFF-4 | Efficiency | Medium | Copying large structures by value in range-for loops |
| EFF-5 | Efficiency | Low | `push_back` without `reserve` in JSON loading |
| EFF-6 | Efficiency | Low | Repeated `pow()` calls in `paramscoldelec` |
| EFF-7 | Efficiency | Low | Full array construction just to get column count |
