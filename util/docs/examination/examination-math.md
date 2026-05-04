# Math Utilities - Code Examination

## Overview

This group of code provides basic mathematical and statistical utility
functions for the Wire-Cell Toolkit: integer GCD computation via the
Euclidean algorithm, a nearest-coprime search, a running statistics
accumulator (count, mean, RMS), and simple histogram generation from
accumulated samples.

## Files Examined

- `util/inc/WireCellUtil/Math.h`
- `util/src/Math.cxx`
- `util/inc/WireCellUtil/Stats.h`

## Algorithm Description

### GCD (Greatest Common Divisor)

Implemented in `Math.cxx:3-15` using the Euclidean algorithm in an
iterative two-step form.  Each iteration computes `a = a % b` and then
`b = b % a`, returning whichever operand first reaches zero.  This is
equivalent to the standard Euclidean algorithm but performs two modulo
operations per loop iteration rather than one.

### nearest_coprime

Implemented in `Math.cxx:17-44`.  Given a `number` and a `target`,
searches outward from `target` in both directions (below and above) to
find the nearest integer that is coprime to `number` (GCD == 1).  The
search is bounded to the half of `[0, number)` that contains `target`:
either `[0, number/2)` or `[number/2, number)`.  The search alternates
between checking `target - step` and `target + step` for increasing
`step`, returning the first coprime found.  Returns 0 if no coprime is
found.

### Stats Accumulator

Implemented in `Stats.h:13-83`.  Maintains running sums S0 (count), S1
(sum of values), and S2 (sum of squares).  Provides:

- `mean()` (line 19): Returns `S1/S0`, with a special-case return of
  0.0 when `S0 == 0`.
- `rms()` (line 24): Computes sample standard deviation using
  `sqrt((S2 - S1*S1/S0) / (S0-1))`, the Bessel-corrected formula.
  Returns -1.0 when `S0 <= 1`.  Clamps negative discriminants (from
  floating-point round-off) to zero.
- `operator()` (line 32): Accumulates a value into S0, S1, S2.
  Optionally appends to a samples vector for later histogramming.

### Histogram Generation

Implemented in `Stats.h:44-81`.  Builds a histogram from stored samples.
If the binning has min == max == 0, it auto-detects the range from the
sample min/max.  Optionally adds two extra bins for underflow/overflow.
Uses `Binning::bin()` to map each sample to a bin index.

## Potential Bugs

### BUG-1: GCD crashes on zero input

- **File**: `util/src/Math.cxx:6`
- **Severity**: High
- **Description**: `GCD(a, b)` performs `a = a % b` on line 6 as its
  first operation.  If `b == 0`, this is undefined behavior (division by
  zero).  Similarly, if `a` is initially 0, line 7 returns `b`, which is
  correct, but the function provides no documented contract about zero
  inputs.  If both `a == 0` and `b == 0`, the modulo on line 6 is
  undefined behavior.
- **Impact**: Calling `GCD(x, 0)` or `GCD(0, 0)` causes a crash or
  undefined behavior.  This can happen if `nearest_coprime` is called
  with `number == 0`.

### BUG-2: nearest_coprime returns 0 as error sentinel but 0 is also a valid search value

- **File**: `util/src/Math.cxx:43`
- **Severity**: Medium
- **Description**: The function returns 0 to indicate "no coprime
  found."  However, when `target` is 0 and `step` is 0, it tests
  `below = target - step = 0` and checks `GCD(number, 0)`.  Since
  `GCD(n, 0) = n` (not 1 for n > 1), the value 0 is never coprime, so
  it never returns 0 as a "found" result.  But callers cannot distinguish
  "searched and found nothing" from a hypothetical coprime value of 0.
  The sentinel value is confusing and potentially error-prone for callers
  that do not carefully check.
- **Impact**: Callers that use the return value without checking for 0
  may use an invalid value.

### BUG-3: nearest_coprime unsigned underflow when target is 0

- **File**: `util/src/Math.cxx:27`
- **Severity**: High
- **Description**: `below = target - step` is computed with `size_t`
  (unsigned) arithmetic.  When `step > target`, this wraps around to a
  very large value.  The guard `below >= lo` (line 28) would usually
  catch this since the wrapped value is enormous, but if `lo` is 0 the
  guard `below >= 0` is always true for unsigned types.  On the first
  iteration `step == 0` so `below == target`, but on subsequent
  iterations when `target < step`, the wrapped value passes the
  `below >= lo` check (since `lo == 0` and any `size_t` is >= 0) and
  then `GCD(number, huge_value)` is called with an unintended argument.
- **Impact**: Incorrect coprime search results or unexpected behavior
  when `target` is small relative to the search range.  In practice this
  means calling with small `target` values in the lower half
  (`lo == 0`) can produce wrong results.

### BUG-4: Integer division in target > number/2 creates boundary ambiguity

- **File**: `util/src/Math.cxx:20`
- **Severity**: Low
- **Description**: The comparison `target > number/2` uses integer
  division.  For odd `number` (e.g., `number = 7`), `number/2 = 3`, so
  `target = 3` goes to the lower half `[0, 3)` and `target = 4` goes to
  the upper half `[3, 7)`.  The boundary value `target == number/2` is
  always placed in the lower half with `hi = number/2`, but the upper
  half starts at `lo = number/2`.  This means `target == number/2` is
  searched in `[0, number/2)` where `hi = number/2` is exclusive (the
  loop checks `above < hi`), and the target itself is at the boundary.
  The `below` path handles it (since `target >= lo` with `lo = 0`), but
  the `above` path at `step = 0` checks `target < number/2` which is
  false, so the target value itself is only checked via the `below`
  path.  This is correct but fragile and non-obvious.
- **Impact**: Minor; the logic happens to be correct but the boundary
  handling is subtle and could easily break under modification.

### BUG-5: Stats::rms() returns -1.0 as error sentinel

- **File**: `util/inc/WireCellUtil/Stats.h:25`
- **Severity**: Medium
- **Description**: When `S0 <= 1`, `rms()` returns -1.0.  This is an
  in-band error sentinel that callers must know to check.  A negative
  RMS is mathematically impossible, but callers that use the return
  value in arithmetic (e.g., computing error bars, tolerances) without
  checking for -1.0 will get incorrect results.  Returning NaN or
  throwing an exception would be safer.
- **Impact**: Silent propagation of the sentinel value -1.0 into
  downstream calculations.  The `operator<<` on line 85 will print
  "-1" for RMS without any warning.

### BUG-6: Stats::mean() returns 0.0 for empty accumulator

- **File**: `util/inc/WireCellUtil/Stats.h:20`
- **Severity**: Low
- **Description**: When no values have been accumulated (`S0 == 0`),
  `mean()` returns 0.0.  This is indistinguishable from a genuine mean
  of 0.0 (e.g., accumulating a single value of 0).  Callers have no way
  to tell whether the mean is truly zero or the accumulator was empty,
  short of separately checking `S0`.
- **Impact**: Callers that do not independently check `S0` may
  misinterpret an empty accumulator as having a mean of zero.

### BUG-7: operator<< defined in header without inline

- **File**: `util/inc/WireCellUtil/Stats.h:85`
- **Severity**: High
- **Description**: The free function `operator<<(std::ostream&, const
  Stats&)` is defined directly in the header file without the `inline`
  specifier and is not a template.  If this header is included by more
  than one translation unit, the linker will see multiple definitions of
  this function, violating the One Definition Rule (ODR) and causing a
  linker error (multiple definition of
  `WireCell::operator<<(std::ostream&, WireCell::Stats const&)`).
- **Impact**: Link failure if `Stats.h` is included in more than one
  `.cxx` file in the same binary.

## Efficiency Considerations

### EFF-1: Stats accumulator uses naive sum-of-squares formula

- **File**: `util/inc/WireCellUtil/Stats.h:24-29`
- **Description**: The RMS computation uses `S2 - S1*S1/S0`, the naive
  two-pass-equivalent formula.  For large datasets with values far from
  zero, catastrophic cancellation can make `S2 - S1*S1/S0` lose
  significant precision, potentially going negative (handled by the
  clamp on line 28).  Welford's online algorithm would be numerically
  superior.
- **Suggestion**: Consider Welford's algorithm for numerically stable
  online variance computation, especially if the accumulator is used
  with large counts or values with large magnitude.

### EFF-2: Samples vector grows unbounded

- **File**: `util/inc/WireCellUtil/Stats.h:37`
- **Description**: When `sample == true`, every accumulated value is
  appended to `samples` with no size limit.  For long-running
  accumulations this can consume unbounded memory.
- **Suggestion**: Consider adding a maximum sample count, reservoir
  sampling, or documenting the expected usage pattern.  Alternatively,
  allow callers to enable sampling only when histogramming is needed and
  disable it afterward.

### EFF-3: Histogram auto-range recomputes min/max

- **File**: `util/inc/WireCellUtil/Stats.h:51`
- **Description**: When auto-detecting the range, `std::minmax_element`
  performs a linear scan over all samples.  This is fine for occasional
  use but could be avoided by tracking min/max incrementally in
  `operator()`.
- **Suggestion**: Track running min and max in the accumulator to avoid
  the extra pass.

## Summary

| ID    | Type | Severity | Summary                                           |
|-------|------|----------|---------------------------------------------------|
| BUG-1 | Bug  | High     | GCD crashes on zero input (division by zero)       |
| BUG-2 | Bug  | Medium   | nearest_coprime uses 0 as ambiguous error sentinel |
| BUG-3 | Bug  | High     | Unsigned underflow in nearest_coprime subtraction  |
| BUG-4 | Bug  | Low      | Integer division creates subtle boundary behavior  |
| BUG-5 | Bug  | Medium   | rms() returns -1.0 as in-band error sentinel       |
| BUG-6 | Bug  | Low      | mean() returns 0.0 for empty, indistinguishable    |
| BUG-7 | Bug  | High     | operator<< in header without inline causes ODR     |
| EFF-1 | Eff  | Medium   | Naive sum-of-squares formula, cancellation risk    |
| EFF-2 | Eff  | Low      | Unbounded samples vector growth                    |
| EFF-3 | Eff  | Low      | Auto-range histogram rescans for min/max           |
