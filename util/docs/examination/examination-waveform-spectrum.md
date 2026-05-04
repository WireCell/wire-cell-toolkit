# Waveform & Spectrum Processing - Code Examination

## Overview

This group of code provides primitive operations for time-domain waveforms and frequency-domain spectra in the Wire-Cell Toolkit. `Waveform.h/.cxx` supplies types, statistics (mean, RMS, median, percentile), bin-range merging, channel-mask merging, time-domain resampling, and element-wise arithmetic on sequences. `Spectrum.h` supplies frequency-domain operations including interpolation, extrapolation, aliasing, resampling, and Hermitian symmetry enforcement. `EigenFFT.h` is a thin wrapper that pulls in the Eigen unsupported FFT module.

## Files Examined

- `util/inc/WireCellUtil/Waveform.h`
- `util/src/Waveform.cxx`
- `util/inc/WireCellUtil/Spectrum.h`
- `util/inc/WireCellUtil/EigenFFT.h`
- `util/inc/WireCellUtil/Interpolate.h` (dependency of Spectrum.h)

## Algorithm Description

### Waveform Resampling (Waveform.h:91-115)

Time-domain resampling via linear interpolation. For each output sample, it computes the corresponding fractional index in the input waveform and performs a weighted average of the two bracketing input samples. Boundary clamping is applied: values outside the original domain are set to the first or last sample. The interpolation weights `d1` and `d2` are computed from the fractional position.

### Statistics: mean_rms (Waveform.cxx:19-39)

Computes population mean and RMS from a real-valued waveform. To avoid float precision issues, the input float vector is first uplifted to double. Uses the identity `rms = sqrt((sum2 - sum^2/n) / n)` where `sum2` is sum of squares and `sum` is the sum.

### Statistics: percentile / percentile_binned (Waveform.cxx:89-134)

`percentile()` uses `std::nth_element` for O(n) selection on a copy of the waveform. `percentile_binned()` builds a histogram with `nbins = wave.size()` bins and walks the cumulative histogram to find the bin at the target percentile.

### Channel Mask Merging (Waveform.cxx:155-216)

`merge(BinRangeList)` sorts ranges and merges overlapping ones in a single pass. `merge(BinRangeList, BinRangeList)` concatenates two lists and delegates to the single-list merge. `merge(ChannelMasks)` and `merge(ChannelMaskMap)` apply the same logic per-channel and per-label.

### Spectrum Interpolation (Spectrum.h:162-172)

Delegates to `linterpolate()` from `Interpolate.h`, which maps both input and output onto the domain [0,1] and uses `linterp` for linear interpolation. The result is scaled by `sqrt(newsize/oldsize)` to preserve energy (RMS normalization).

### Spectrum Extrapolation (Spectrum.h:193-230)

Extends a spectrum to a larger size by copying the lower and upper frequency halves to the appropriate positions in the output and filling the gap (new high-frequency bins) with a constant value. Handles both even and odd-length spectra differently with respect to the Nyquist bin.

### Spectrum Aliasing (Spectrum.h:250-280)

Downsamples a spectrum by folding higher-frequency bins back onto lower-frequency bins (spectral aliasing). Computes only the positive half-spectrum, then applies `hermitian_mirror()` to fill the negative frequencies. Normalizes by `sqrt(newsize/oldsize)`.

### Hermitian Mirror Symmetry (Spectrum.h:56-89)

Enforces Hermitian symmetry on a complex spectrum so that IFFT produces a real-valued signal. Sets DC bin to its real part, and for even-length spectra also sets the Nyquist bin to its real part. Copies positive-frequency bins to negative-frequency positions in reverse order, then conjugates the negative-frequency half.

### Spectrum Resample (Spectrum.h:293-309)

Two-step operation: first interpolates the spectrum to an intermediate size `ceil(olen*relperiod)`, then either aliases (if relperiod > 1, i.e., period increases / downsampling in time) or extrapolates (if relperiod < 1, period decreases / upsampling in time).

### most_frequent (Waveform.cxx:218-227)

Computes the most frequent value in a vector of `short` by building a histogram of size 2^16 indexed by treating values as `unsigned short`. Returns the smallest index at the maximum count.

## Potential Bugs

### [BUG-1] Incorrect interpolation weights in Waveform::resample

- **File**: `util/inc/WireCellUtil/Waveform.h:109-111`
- **Severity**: High
- **Description**: The interpolation weight `d1` is computed as `oldfracsteps - oldstep * oldind`. The variable `oldfracsteps` is `(cursor - domain.first) / oldstep`, which is already in units of sample indices (dimensionless). But `oldstep * oldind` has units of time (since `oldstep` has units of time and `oldind` is an integer index). These are dimensionally inconsistent. The correct fractional part should be `oldfracsteps - oldind` (both dimensionless). Furthermore, the interpolation formula `(wave[oldind] * d1 + wave[oldind+1] * d2) / oldstep` divides by `oldstep`, but `d1 + d2 = oldstep` only under the current (incorrect) dimensional mixing. If the intent is linear interpolation, the formula should be `wave[oldind] * (1 - frac) + wave[oldind+1] * frac` where `frac = oldfracsteps - oldind`. The current code happens to produce correct results because the errors cancel: `d1 = oldfracsteps - oldstep*oldind` (wrong dimensionally but numerically equals `(cursor - domain.first)/oldstep - oldstep*oldind`), then dividing by `oldstep` yields ... the math does reduce but the weights are swapped. Let `f = oldfracsteps - oldind` be the true fractional index. Then `d1 = f*oldstep - oldstep*oldind + oldind*oldstep`... Actually, expanding carefully: `d1 = oldfracsteps - oldstep*oldind`. Since `oldfracsteps = (cursor - domain.first)/oldstep`, we get `d1 = (cursor - domain.first)/oldstep - oldstep*oldind`. This is not dimensionally consistent. If `oldstep` is, say, `0.5 us` and `oldind = 3`, `oldfracsteps = 3.7`, then `d1 = 3.7 - 0.5*3 = 2.2`, `d2 = 0.5 - 2.2 = -1.7`. This produces negative weights and incorrect interpolation. The formula is only correct when `oldstep = 1.0`.
- **Impact**: Produces incorrect resampled waveforms whenever the sample width is not 1.0, which is typical for real physics data (e.g., 0.5 us tick size). This can silently corrupt waveform data.

### [BUG-2] merge(BinRangeList) crashes on empty input

- **File**: `util/src/Waveform.cxx:160`
- **Severity**: Medium
- **Description**: `merge(const BinRangeList& brl)` accesses `tmp[0]` without checking if `tmp` is empty. If an empty `BinRangeList` is passed, this results in undefined behavior (out-of-bounds access).
- **Impact**: Crash or undefined behavior when merging an empty bin range list.

### [BUG-3] merge(BinRangeList) does not take maximum of endpoints for overlapping ranges

- **File**: `util/src/Waveform.cxx:165-166`
- **Severity**: Medium
- **Description**: When merging overlapping ranges, line 166 unconditionally sets `out.back().second = this_br.second`. If the current range is entirely contained within the previous range (e.g., [0,10) and [2,5)), the merged range incorrectly shrinks to [0,5) instead of remaining [0,10). The fix is `out.back().second = std::max(out.back().second, this_br.second)`.
- **Impact**: Merged bin ranges can shrink rather than grow, leading to loss of masked regions and incorrect downstream masking.

### [BUG-4] percentile_binned division by zero when all values are identical

- **File**: `util/src/Waveform.cxx:114`
- **Severity**: Medium
- **Description**: When `vmin == vmax`, `binsize = 0`, leading to division by zero at line 117 `(val - vmin) / binsize`. The function does not handle the degenerate case of a constant waveform.
- **Impact**: Floating-point exception or NaN/Inf propagation for constant-valued waveforms.

### [BUG-5] Spectrum::alias uses incorrect fold stride

- **File**: `util/inc/WireCellUtil/Spectrum.h:265-269`
- **Severity**: Medium
- **Description**: The aliasing loop computes `L = ceil(ilen/olen)` as the number of folds and uses stride `M = olen/2` via `oldind = m + l*M`. However, the correct stride for spectral aliasing should be `olen` (the decimation factor), not `olen/2`. Using `M = olen/2` means each output bin incorrectly accumulates from input bins spaced by `olen/2` instead of by `olen`, effectively double-folding. Furthermore, the loop condition `oldind > half` uses the input half-spectrum size as the bound, which is correct for zero-extending, but the stride error means bins are accumulated from wrong positions.
- **Impact**: Incorrect spectral aliasing results when downsampling. The aliased spectrum will not correctly represent the folded higher-frequency content.

### [BUG-6] Spectrum::extrap does not normalize output

- **File**: `util/inc/WireCellUtil/Spectrum.h:192-230`
- **Severity**: Low
- **Description**: The docstring at line 188-190 states "The new spectrum is normalized by sqrt(newsize/oldsize)" but the implementation does not apply any normalization factor. The `interp()` function applies normalization, but `extrap()` does not, which is inconsistent with the documented behavior.
- **Impact**: Energy is not conserved as documented when using `extrap()` directly. Users relying on the documented normalization will get incorrect amplitudes. However, when called via `resample()`, the preceding `interp()` step applies its own normalization, which may partially compensate.

### [BUG-7] Spectrum::resample applies double normalization via interp

- **File**: `util/inc/WireCellUtil/Spectrum.h:294-309`
- **Severity**: Low
- **Description**: `resample()` calls `interp()` which normalizes by `sqrt(tmpsiz/ilen)`, then calls either `alias()` which normalizes by `sqrt(olen/tmpsiz)` or `extrap()` which does not normalize. When using `alias()`, the combined normalization is `sqrt(tmpsiz/ilen) * sqrt(olen/tmpsiz) = sqrt(olen/ilen)`, which is correct. But when using `extrap()` (relperiod <= 1), only the interp normalization `sqrt(tmpsiz/ilen)` is applied, where `tmpsiz = ceil(olen*relperiod)`. This may not equal the desired `sqrt(olen/ilen)`.
- **Impact**: Potentially incorrect normalization when resampling with `relperiod <= 1`.

### [BUG-8] Spectrum::hermitian_mirror second form may read uninitialized output data

- **File**: `util/inc/WireCellUtil/Spectrum.h:80-89`
- **Severity**: Low
- **Description**: The two-iterator form copies `ibeg..ihalf` to `obeg`, then calls the in-place `hermitian_mirror(obeg, obeg+len)`. The in-place version reads from `beg+1..mid` to fill `mid+1..end`. But only the first `len/2+1` elements of the output were initialized by the copy. The in-place call then reads these initialized elements and writes to the rest, so this is actually correct -- the uninitialized elements are written before being read. This is safe but relies on the specific read/write pattern of `reverse_copy`.
- **Impact**: No actual bug, but the design is fragile and could break if the in-place implementation changes.

### [BUG-9] most_frequent reinterprets signed short as unsigned

- **File**: `util/src/Waveform.cxx:218-227`
- **Severity**: Low
- **Description**: The function takes `const std::vector<short>& vals` but the loop iterates with `unsigned short val`, causing implicit conversion. Negative short values (e.g., -1 = 0xFFFF as unsigned short) map to high histogram indices, which is actually the intended behavior for hashing into the 2^16 histogram. However, the return value `it - hist.begin()` is cast back to `short`, so the returned "most frequent value" correctly round-trips. The implicit signed-to-unsigned conversion may trigger compiler warnings.
- **Impact**: Functional but may produce compiler warnings. The semantics are correct due to the two's complement round-trip.

## Efficiency Considerations

### [EFF-1] Waveform::resample uses push_back without reserve

- **File**: `util/inc/WireCellUtil/Waveform.h:96-113`
- **Description**: The `ret` vector grows via `push_back` in a loop of `nsamples` iterations without calling `reserve(nsamples)` first, causing multiple reallocations.
- **Suggestion**: Add `ret.reserve(nsamples)` before the loop at line 97.

### [EFF-2] percentile makes a full copy of the waveform

- **File**: `util/src/Waveform.cxx:103`
- **Description**: `percentile()` copies the entire waveform to use `nth_element`. This is necessary to avoid modifying the input, but for large waveforms called repeatedly, this is costly.
- **Suggestion**: Consider accepting a mutable reference overload or using a move-based API for cases where the caller does not need the original.

### [EFF-3] percentile_binned does two passes over data

- **File**: `util/src/Waveform.cxx:110-134`
- **Description**: First calls `std::minmax_element` (one pass), then iterates again to build the histogram (second pass). This is inherent to the algorithm but could be combined into a single pass that tracks min/max while building an initial estimate, then refines.
- **Suggestion**: For very large waveforms, consider a single-pass approach or streaming percentile estimation.

### [EFF-4] merge(ChannelMasks) copies the entire first map

- **File**: `util/src/Waveform.cxx:189`
- **Description**: `ChannelMasks out = one` copies the entire first map. If the first map is large and only a few channels overlap, this is wasteful.
- **Suggestion**: Minor concern; acceptable for typical use cases. Could accept by const-ref and build output incrementally if profiling shows this is a bottleneck.

### [EFF-5] most_frequent allocates 256KB histogram on stack

- **File**: `util/src/Waveform.cxx:221`
- **Description**: Allocates a `std::vector<unsigned int>` of size 65536 (256 KB) on every call. For repeated calls this allocation and zero-initialization is non-trivial.
- **Suggestion**: Consider using a static thread-local vector or passing in a pre-allocated buffer if called in a hot loop.

### [EFF-6] Spectrum::interp creates intermediate linterp object

- **File**: `util/inc/WireCellUtil/Interpolate.h:351-361` (called from `Spectrum.h:165`)
- **Description**: `linterpolate()` creates a `linterp` object that copies the input range into an internal `std::vector`. When the input is already a vector, this is a redundant copy.
- **Suggestion**: Consider providing an overload that takes a data pointer/span to avoid the copy.

### [EFF-7] Spectrum::alias computes only half spectrum but normalizes full spectrum

- **File**: `util/inc/WireCellUtil/Spectrum.h:277-279`
- **Description**: After `hermitian_mirror` fills the full spectrum, the normalization transform iterates over all `olen` elements. The normalization could be applied to just the half-spectrum before mirroring, saving roughly half the multiplications.
- **Suggestion**: Apply normalization before `hermitian_mirror()`.

### [EFF-8] mean_rms copies entire waveform to double vector

- **File**: `util/src/Waveform.cxx:31`
- **Description**: The entire float waveform is copied to a `std::vector<double>` to avoid precision issues. For large waveforms (e.g., 10k+ samples), this doubles memory usage.
- **Suggestion**: Use Kahan/compensated summation on the float data directly, or accumulate into double accumulators without copying the entire vector. A simple loop with `double sum = 0; for (auto v : wf) sum += v;` avoids the copy.

## Summary

| ID | Category | Severity | Title |
|----|----------|----------|-------|
| BUG-1 | Bug | High | Incorrect interpolation weights in Waveform::resample |
| BUG-2 | Bug | Medium | merge(BinRangeList) crashes on empty input |
| BUG-3 | Bug | Medium | merge(BinRangeList) does not take max of endpoints for overlapping ranges |
| BUG-4 | Bug | Medium | percentile_binned division by zero when all values identical |
| BUG-5 | Bug | Medium | Spectrum::alias uses incorrect fold stride |
| BUG-6 | Bug | Low | Spectrum::extrap does not normalize output despite documentation |
| BUG-7 | Bug | Low | Spectrum::resample applies inconsistent normalization via interp+extrap path |
| BUG-8 | Bug | Low | Spectrum::hermitian_mirror second form has fragile read/write ordering |
| BUG-9 | Bug | Low | most_frequent signed-to-unsigned implicit conversion |
| EFF-1 | Efficiency | Medium | resample uses push_back without reserve |
| EFF-2 | Efficiency | Low | percentile makes full copy of waveform |
| EFF-3 | Efficiency | Low | percentile_binned does two passes over data |
| EFF-4 | Efficiency | Low | merge(ChannelMasks) copies entire first map |
| EFF-5 | Efficiency | Low | most_frequent allocates 256KB histogram per call |
| EFF-6 | Efficiency | Low | Spectrum::interp creates redundant copy via linterp |
| EFF-7 | Efficiency | Low | Spectrum::alias normalizes full spectrum after mirror |
| EFF-8 | Efficiency | Medium | mean_rms copies entire waveform to double vector |
