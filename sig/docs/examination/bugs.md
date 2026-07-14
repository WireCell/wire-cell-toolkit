# Potential Bugs in sig/

Code examined: `sig/src/Decon2DResponse.cxx`, `sig/src/Decon2DFilter.cxx`, `sig/src/Util.cxx`, `sig/test/check_tensorset_filter.cxx`

---

## BUG-1: Reversed linear interpolation in `init_overall_response()` (Decon2DResponse.cxx:195-196)

**Severity: High**

The redigitization interpolation has the weights swapped:

```cpp
wfs.at(i) = ((ctime - ftbins.at(fcount - 1)) / fravg.period * arr(irow, fcount - 1) +
             (ftbins.at(fcount) - ctime) / fravg.period * arr(irow, fcount));
```

Standard linear interpolation between points `A` (at `t_lo`) and `B` (at `t_hi`) for query time `t` is:

```
result = (t_hi - t) / (t_hi - t_lo) * A  +  (t - t_lo) / (t_hi - t_lo) * B
```

The code uses `(ctime - ftbins[fcount-1])` as the weight for `arr(irow, fcount-1)` (the left point), but it should be the weight for `arr(irow, fcount)` (the right point). The two weights are swapped.

When `ctime` is close to `ftbins[fcount-1]`, the weight `(ctime - ftbins[fcount-1])` is small, but the code multiplies it with `arr(irow, fcount-1)` (which should have large weight). This means the interpolation gives more weight to the *farther* point rather than the *closer* point.

**Note:** If the response function is smooth and the fine/coarse bin sizes are similar, the numerical impact may be small but the logic is incorrect. Additionally, the denominator uses `fravg.period` instead of `(ftbins[fcount] - ftbins[fcount-1])` -- these are equal by construction so that part is correct.

---

## BUG-2: Division by complex zero without protection in Decon2DResponse.cxx:331

**Severity: Medium**

```cpp
c_data = c_data / c_resp;
```

This performs element-wise division of the signal by the response in the 2D frequency domain. If any element of `c_resp` is zero (or near-zero), this produces `inf` or `NaN`. The NaN/inf cleanup happens later (lines 341-345), but division by near-zero values can produce extremely large floating-point values that are not `inf` but are still numerically harmful. The NaN/inf check only catches the extreme cases, not large-but-finite garbage values.

By contrast, the per-channel response correction (lines 300-308) properly checks `std::abs(four) != 0` before dividing.

**Recommendation:** Apply a threshold-based protection (similar to the per-channel case) or use a Wiener-style regularized division: `c_data * conj(c_resp) / (|c_resp|^2 + epsilon)`.

---

## BUG-3: CMM data cast to `double*` without type verification (Decon2DFilter.cxx:130-131)

**Severity: Medium**

```cpp
Eigen::Map<Eigen::ArrayXXd> ranges_arr((double *) cmm_range->data(), nranges, 2);
Eigen::Map<Eigen::ArrayXd> channels_arr((double *) cmm_channel->data(), nranges);
```

The code casts the raw tensor data to `double*` without verifying the tensor's actual element type. If the upstream producer stores CMM data as `float` or `int`, this reinterprets the memory incorrectly, producing garbage values. The waveform tensor is correctly handled using `itensor_to_eigen_array<>` which is type-aware, but the CMM tensors bypass this mechanism.

---

## BUG-4: CMM channel used as direct row index without channel-to-row mapping (Decon2DFilter.cxx:132-134)

**Severity: Medium**

```cpp
auto ch = channels_arr(ind);
if (!(ch < r_data.rows())) {
    continue;
}
```

The `ch` value from the CMM is treated as a direct row index into `r_data`. However, the tensor rows correspond to the channels listed in `ch_arr` (the channel tensor), which may not start at 0 or be contiguous. If the CMM stores actual channel IDs (e.g., channel 2400), but the tensor only has rows 0..N for a subset of channels, `ch` would be out of range and silently skipped (or index the wrong row).

The code should map from channel ID to row index using `ch_arr`.

---

## BUG-5: `m_pad_nwires` always zero, making wire-padding dead code (Decon2DFilter.cxx:87-88,122)

**Severity: Low**

```cpp
int m_pad_nwires = 0;
// ...
Array::array_xxf r_data = tm_r_data.block(m_pad_nwires, 0, m_nwires, m_nticks);
```

`m_pad_nwires` is hardcoded to 0 and `m_nwires = nchans`. The `.block()` call is thus equivalent to just using the full array. The FIXME on line 89 suggests this is known incomplete. If wire-padding is ever needed, this will silently produce wrong results because the padding logic was never wired to configuration.

---

## BUG-6: Negative time_shift not handled in Decon2DResponse.cxx:372-379

**Severity: Medium**

```cpp
int time_shift = (m_coarse_time_offset + m_intrinsic_time_offset) / m_period;
if (time_shift > 0) {
    // ... circular shift logic
}
```

The code only handles positive `time_shift`. With default `m_coarse_time_offset = -8.0 us` and typical `m_intrinsic_time_offset` values, `time_shift` could be negative. A negative shift is silently ignored, meaning the time alignment is not applied. The same issue exists for `fine_time_shift` in `init_overall_response()` (line 174).

---

## BUG-7: Unused variable `filter_names` in both files

**Severity: Low (code quality)**

`Decon2DFilter.cxx:84`:
```cpp
const std::vector<std::string> filter_names{"Wire_ind", "Wire_ind", "Wire_col"};
```

This vector is declared but never used in `Decon2DFilter`. It *is* used in `Decon2DResponse.cxx:246,336` via `filter_names[iplane]`, but in `Decon2DFilter` it is dead code.

---

## BUG-8: `std::abs()` on complex vs `float` ambiguity (Decon2DResponse.cxx:341)

**Severity: Low**

```cpp
float val = abs(c_data(irow, icol));
```

`abs()` (unqualified) on a `std::complex<float>` may call the integer `abs()` on some compilers/platforms rather than `std::abs()` for complex types, which returns the magnitude. Should use `std::abs()` explicitly to ensure correct complex magnitude computation.

---

## BUG-9: Test does not configure the node before use (check_tensorset_filter.cxx:46-51)

**Severity: Low (test-only)**

```cpp
auto node = std::make_shared<Sig::Decon2DResponse>();
auto cfg = node->default_configuration();
// cfg is never passed to node->configure(cfg)
(*node)(iitens, oitens);
```

`default_configuration()` is called but `configure()` is never called with it. The node is used unconfigured. Since `operator()` accesses `m_anode`, `m_fresp`, `m_dft` (all nullptr), this test would crash if actually run. This suggests it is a skeleton/placeholder test that has never been executed.

---

## BUG-10: Raw `new` for ITensor::vector without RAII (Decon2DFilter.cxx:151, Decon2DResponse.cxx:391)

**Severity: Low**

```cpp
ITensor::vector *itv = new ITensor::vector;
```

If any operation between this `new` and the wrapping in `ITensor::shared_vector(itv)` throws, the memory leaks. This is a minor concern since the intervening code is unlikely to throw, but it is a pattern that should use `std::make_shared` or direct construction.
