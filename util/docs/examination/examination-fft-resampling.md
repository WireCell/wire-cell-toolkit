# FFT & Resampling - Code Examination

## Overview

This group of code provides two core signal-processing utilities for the Wire-Cell Toolkit:

1. **FFT optimal length selection** -- given an arbitrary sample count, return the nearest highly-composite number that yields efficient FFT execution.
2. **LMN rational resampling** -- implement frequency-domain resampling by computing a floating-point GCD to find rational relationships between sampling periods, resizing real-valued arrays, and resampling complex spectra by splitting positive/negative frequency halves.

## Files Examined

| File | Description |
|------|-------------|
| `util/inc/WireCellUtil/FFTBestLength.h` | Header declaring `fft_best_length()` |
| `util/src/FFTBestLength.cxx` | Implementation with hard-coded highly-composite number tables |
| `util/inc/WireCellUtil/LMN.h` | Header declaring LMN resampling functions (`gcd`, `rational`, `nhalf`, `nbigger`, `resize`, `resample`) |
| `util/src/LMN.cxx` | Implementation of LMN resampling functions |

## Algorithm Description

### FFT Best Length (`fft_best_length`)

Selects an FFT-friendly size greater than or equal to the requested `window_length`. Three separate lookup tables of highly-composite numbers are maintained: one general-purpose table, one for even lengths, and one for odd lengths. The `keep_odd_even` flag controls whether the returned length must preserve the odd/even parity of the input. A linear scan finds the first table entry >= the requested size; if no entry qualifies the original length is returned unchanged.

### Floating-Point GCD (`LMN::gcd`)

Implements a recursive Euclidean algorithm using `fmod()` instead of the integer modulus operator, allowing GCD computation on floating-point values. Recursion terminates when the first argument falls below `eps`.

### Rational Resampling Size (`LMN::rational`)

Given source period `Ts` and target period `Tr`, computes the minimum number of source samples that allows exact rational resampling. Uses the floating-point GCD of `Tr` and `|Ts - Tr|` to derive the ratio, with two error checks to ensure the result is sufficiently close to an integer.

### Spectrum Half-Size (`LMN::nhalf`)

Returns the number of positive-frequency bins in a DFT of length N, excluding the DC bin and any Nyquist bin. For odd N: `(N-1)/2`; for even N: `(N-2)/2`.

### Padded Size (`LMN::nbigger`)

Returns the smallest multiple of `Nrat` that is >= N.

### Real Array Resize (`LMN::resize`)

Copies the first `min(Nr, Ns)` samples from the source into a zero-initialized output of size `Nr`. Provides both `std::vector<float>` and Eigen `array_xxf` overloads.

### Complex Spectrum Resample (`LMN::resample`)

Resamples a frequency-domain spectrum by copying positive-frequency and negative-frequency halves from the source into a zero-initialized output. For upsampling, zeros are inserted around the Nyquist frequency; for downsampling, high-frequency content near Nyquist is discarded. Provides both `std::vector<complex<float>>` and Eigen `array_xxc` overloads.

## Potential Bugs

### [BUG-1] `gcd()` does not propagate `eps` in recursive call

- **File**: `util/src/LMN.cxx:19`
- **Severity**: High
- **Description**: The recursive call `LMN::gcd(fmod(b, a), a)` omits the third argument `eps`. Although the header declares a default value of `1e-6`, if a caller invokes `gcd()` with a non-default `eps` (e.g., `1e-3`), the recursive call silently reverts to the default `1e-6`. This is almost certainly unintentional.
- **Impact**: When `eps` differs from the default, the recursion termination threshold changes mid-computation, potentially causing either too many recursion levels (if the caller's eps was larger) or incorrect early termination (if the caller's eps was smaller). In the worst case with extremely small eps values this could lead to deep recursion.

### [BUG-2] Operator precedence error in `nbigger()` -- `!` binds tighter than `%`

- **File**: `util/inc/WireCellUtil/LMN.h:36`
- **Severity**: Critical
- **Description**: The expression `if (! N%Nrat)` is parsed as `if ((!N) % Nrat)` due to C++ operator precedence. The unary `!` operator binds tighter than `%`. The intended expression is `if (!(N % Nrat))`. Because `!N` evaluates to `0` for any nonzero N (the normal case), `0 % Nrat` is always `0` (falsy), so the early-return is never taken. The function always computes `Nrat * (N/Nrat + 1)`, which overshoots by `Nrat` when `N` is already a multiple of `Nrat`.
- **Impact**: Any call to `nbigger()` where `N` is already divisible by `Nrat` returns a value that is `Nrat` too large. This would cause unnecessary over-padding of arrays.

### [BUG-3] Duplicate assignment lines in `resize()` for Eigen arrays

- **File**: `util/src/LMN.cxx:68-69` and `util/src/LMN.cxx:72-73`
- **Severity**: Low
- **Description**: Both the `axis == 0` and `axis == 1` branches contain an identical assignment statement repeated twice:
  ```
  rs(safe, all) = in(safe, all);
  rs(safe, all) = in(safe, all);   // duplicate
  ```
  and
  ```
  rs(all, safe) = in(all, safe);
  rs(all, safe) = in(all, safe);   // duplicate
  ```
  The second assignment in each pair is redundant -- it overwrites with the same data.
- **Impact**: No functional impact (the result is correct), but it performs an unnecessary memory copy on every call. This appears to be a copy-paste artifact.

### [BUG-4] FIXME: Nyquist bin not handled in `resample()`

- **File**: `util/src/LMN.cxx:152` and `util/src/LMN.cxx:176`
- **Severity**: Medium
- **Description**: Both the Eigen and `std::vector` overloads of `resample()` contain the comment `// FIXME: deal with Nyquist bin.` For even-length DFTs, the Nyquist bin (index `N/2`) has special symmetry requirements. Currently the code copies `N_half + 1` positive-frequency bins (indices `0..N_half`) and `N_half` negative-frequency bins. When the source length is even, the Nyquist bin at index `Ns/2` is neither explicitly copied nor split. The `nhalf()` function already excludes the Nyquist bin from its count (returning `(N-2)/2` for even N), so the Nyquist bin's energy may be silently zeroed or double-counted depending on the direction of resampling.
- **Impact**: When resampling spectra of even length, energy at the Nyquist frequency may be lost (upsampling) or incorrectly included (downsampling). This could introduce subtle amplitude errors or aliasing artifacts in resampled waveforms.

### [BUG-5] Typo in `nhalf()` comment -- "event" should be "even"

- **File**: `util/inc/WireCellUtil/LMN.h:32`
- **Severity**: Low (cosmetic)
- **Description**: The comment reads `// event` but should read `// even`.
- **Impact**: No functional impact; documentation-only issue.

### [BUG-6] Wrong variable in error message in `rational()`

- **File**: `util/src/LMN.cxx:39`
- **Severity**: Low
- **Description**: The error message for the second check uses `err1` instead of `err2`: `raise<ValueError>("gcd error two too big %f > %f", err1, eps)`. The conditional correctly checks `err2 > eps`, but the formatted message reports `err1`'s value.
- **Impact**: Misleading diagnostic output when this error path is triggered. The exception is still raised correctly, but the reported error magnitude is wrong.

### [BUG-7] Potential unbounded recursion in `gcd()` for adversarial inputs

- **File**: `util/src/LMN.cxx:14-20`
- **Severity**: Medium
- **Description**: The recursive `gcd()` has no explicit recursion depth limit. While the Euclidean algorithm converges, floating-point `fmod()` can exhibit slow convergence for certain input ratios (e.g., values very close to irrational ratios like the golden ratio). Combined with BUG-1 (eps not propagated), if the effective eps is very small, recursion depth could become large, risking stack overflow.
- **Impact**: Potential stack overflow for pathological inputs. In typical usage with well-behaved sampling periods this is unlikely, but the lack of a depth guard is a latent risk.

## Efficiency Considerations

### [EFF-1] Linear scan in `fft_best_length()` instead of binary search

- **File**: `util/src/FFTBestLength.cxx:53-55`
- **Description**: The lookup tables are pre-sorted, but the code performs a linear scan (`for` loop with sequential comparison) to find the first entry >= `window_length`.
- **Suggestion**: Use `std::lower_bound()` for O(log N) lookup instead of O(N). With table sizes of ~70-80 entries the absolute cost is small, but the fix is trivial and semantically clearer.

### [EFF-2] Lookup table reconstructed on every call

- **File**: `util/src/FFTBestLength.cxx:14-50`
- **Description**: The `edges` vector is constructed and populated via `boost::assign` on every invocation of `fft_best_length()`. This involves heap allocation and element-by-element insertion on each call.
- **Suggestion**: Make the tables `static const` (either as `static const std::vector` initialized once, or preferably as `static constexpr std::array`). This eliminates repeated allocation and would also allow `std::lower_bound()` to work on a compile-time constant array.

### [EFF-3] Recursive `gcd()` could be iterative

- **File**: `util/src/LMN.cxx:14-20`
- **Description**: The floating-point GCD is implemented via tail recursion. While many compilers optimize tail recursion, it is not guaranteed by the C++ standard.
- **Suggestion**: Convert to an iterative `while` loop:
  ```cpp
  while (a >= eps) {
      double t = fmod(b, a);
      b = a;
      a = t;
  }
  return b;
  ```
  This also eliminates BUG-7 (unbounded recursion risk) and naturally fixes BUG-1 (eps propagation) since the parameter stays in scope.

### [EFF-4] `fill_constant()` reimplements `std::fill`

- **File**: `util/src/LMN.cxx:78-86`
- **Description**: `fill_constant()` is a manual while-loop that duplicates `std::fill(begin, end, value)`.
- **Suggestion**: Replace with `std::fill(begin, end, value)` which may benefit from compiler vectorization and specialization.

### [EFF-5] `fill_linear()` recomputes step divisor each iteration

- **File**: `util/src/LMN.cxx:88-96`
- **Description**: The expression `first + ind*((last-first)/N)` recomputes `(last-first)/N` on every iteration, though the compiler may hoist this.
- **Suggestion**: Pre-compute `const float step = (last - first) / N;` before the loop for clarity and to guarantee the division is performed once.

## Summary

| ID | Category | Severity | Title |
|----|----------|----------|-------|
| BUG-1 | Bug | High | `gcd()` does not propagate `eps` in recursive call |
| BUG-2 | Bug | Critical | Operator precedence error in `nbigger()` -- `!` vs `%` |
| BUG-3 | Bug | Low | Duplicate assignment lines in Eigen `resize()` |
| BUG-4 | Bug | Medium | Nyquist bin not handled in `resample()` (FIXME) |
| BUG-5 | Bug | Low | Typo: "event" should be "even" in `nhalf()` comment |
| BUG-6 | Bug | Low | Wrong variable (`err1` vs `err2`) in error message |
| BUG-7 | Bug | Medium | Potential unbounded recursion in `gcd()` |
| EFF-1 | Efficiency | Low | Linear scan instead of binary search in `fft_best_length()` |
| EFF-2 | Efficiency | Medium | Lookup table rebuilt on every call to `fft_best_length()` |
| EFF-3 | Efficiency | Low | Recursive `gcd()` could be iterative (also fixes BUG-1 and BUG-7) |
| EFF-4 | Efficiency | Low | `fill_constant()` reimplements `std::fill` |
| EFF-5 | Efficiency | Low | `fill_linear()` recomputes step divisor each iteration |
