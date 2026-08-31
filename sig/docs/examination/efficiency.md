# Efficiency Observations in sig/

Code examined: `sig/src/Decon2DResponse.cxx`, `sig/src/Decon2DFilter.cxx`, `sig/src/Util.cxx`

---

## EFF-1: Element-wise filter application uses scalar loop instead of Eigen vectorization (Decon2DFilter.cxx:114-118)

**Impact: Medium**

```cpp
for (int irow = 0; irow < c_data.rows(); ++irow) {
    for (int icol = 0; icol < c_data.cols(); ++icol) {
        c_data_afterfilter(irow, icol) = c_data(irow, icol) * filter.at(icol);
    }
}
```

This nested loop applies a 1D filter (same across all rows) to each column. Eigen supports broadcasting:

```cpp
Eigen::Map<const Eigen::ArrayXf> filter_arr(filter.data(), filter.size());
c_data_afterfilter = c_data.rowwise() * filter_arr.transpose().cast<std::complex<float>>();
```

This enables SIMD vectorization. The same pattern applies to Decon2DResponse.cxx:164-168 (response * electronics multiplication) and Decon2DResponse.cxx:338-349 (wire filter application).

---

## EFF-2: `c_data_afterfilter` is an unnecessary full copy (Decon2DFilter.cxx:113)

**Impact: Medium**

```cpp
Array::array_xxc c_data_afterfilter(c_data.rows(), c_data.cols());
```

A new array the same size as `c_data` is allocated, filled element by element, then `c_data` is never used again. The filter could be applied in-place:

```cpp
for (int irow = 0; irow < c_data.rows(); ++irow) {
    for (int icol = 0; icol < c_data.cols(); ++icol) {
        c_data(irow, icol) *= filter.at(icol);
    }
}
```

This halves peak memory usage for this array (which can be large: nchans x nticks complex floats).

---

## EFF-3: `init_overall_response()` is recomputed on every call (Decon2DResponse.cxx:267)

**Impact: High**

```cpp
auto m_overall_resp = init_overall_response(in);
```

`operator()` calls `init_overall_response()` every time a new TensorSet arrives. This function:
- Reads the field response and averages it
- Generates cold electronics response
- Does forward FFT, element-wise multiplication, inverse FFT
- Performs redigitization with linear interpolation

The response depends only on configuration parameters (gain, shaping time, field response, etc.), not on the input data (it only reads `nticks` and `tick` from the input metadata). If these are constant across calls (which they typically are), the response can be computed once in `configure()` or lazily cached after the first call.

---

## EFF-4: Scalar redigitization loop with redundant search (Decon2DResponse.cxx:184-203)

**Impact: Medium**

```cpp
for (int irow = 0; irow < fine_nwires; ++irow) {
    size_t fcount = 1;
    for (int i = 0; i != m_fft_nticks; i++) {
        double ctime = ctbins.at(i);
        if (fcount < fine_nticks)
            while (ctime > ftbins.at(fcount)) { ... }
```

The inner `while` loop re-searches for the correct bin from where it last left off. Because both `ctbins` and `ftbins` are uniformly spaced, the bin index can be computed directly:

```cpp
size_t fcount = static_cast<size_t>(ctime / fravg.period);
```

This eliminates the linear search entirely. Additionally, `fcount` is reset to 1 for each wire row, but since `ctbins` is identical for all wires, the index mapping is the same -- it could be precomputed once outside the row loop.

---

## EFF-5: Temporary array copies for circular shift (Decon2DResponse.cxx:362-368, 372-379)

**Impact: Medium**

Wire shift:
```cpp
Array::array_xxf arr1(m_wire_shift, ncols);
arr1 = r_data.block(nrows - m_wire_shift, 0, m_wire_shift, ncols);
Array::array_xxf arr2(nrows - m_wire_shift, ncols);
arr2 = r_data.block(0, 0, nrows - m_wire_shift, ncols);
r_data.block(0, 0, m_wire_shift, ncols) = arr1;
r_data.block(m_wire_shift, 0, nrows - m_wire_shift, ncols) = arr2;
```

Two large temporary arrays are allocated for a circular shift. Alternatives:
1. Use `Eigen::PermutationMatrix` for row permutation
2. Perform the circular shift in the frequency domain (multiply by a phase ramp before the inverse FFT) -- zero extra memory
3. At minimum, only one temporary is needed (copy the smaller segment, shift the larger in-place, paste back)

The same pattern appears for time shift (lines 372-379) and in `init_overall_response()` (lines 175-181).

---

## EFF-6: Response array `r_resp` zero-initialized then partially filled (Decon2DResponse.cxx:317-322)

**Impact: Low**

```cpp
Array::array_xxf r_resp = Array::array_xxf::Zero(r_data.rows(), m_fft_nticks);
for (size_t i = 0; i != m_overall_resp.size(); i++) {
    for (int j = 0; j != m_fft_nticks; j++) {
        r_resp(i, j) = m_overall_resp.at(i).at(j);
    }
}
```

The array has `r_data.rows()` rows but only `m_overall_resp.size()` rows are filled (typically a small number like 10-20 paths). The remaining rows stay zero. This is correct for the FFT but allocates much more memory than needed. The FFT on axis 0 pads implicitly if the response is placed in a smaller array. However, the zero-fill is needed for the 2D FFT so this is more of a minor note.

The inner copy loop could use `Eigen::Map` to avoid the scalar loop:
```cpp
Eigen::Map<const Eigen::ArrayXf> row(m_overall_resp[i].data(), m_fft_nticks);
r_resp.row(i) = row;
```

---

## EFF-7: `restore_baseline()` iterates over columns 3 times (Util.cxx:5-34)

**Impact: Low-Medium**

For each row:
1. First pass (lines 10-15): collect non-zero values into `signal`
2. Compute median of `signal`
3. Second pass (lines 21-26): collect values within 500 of median into `temp_signal`
4. Compute median of `temp_signal`
5. Third pass (lines 31-33): subtract baseline from non-zero entries

Each `median()` call internally sorts the data (O(n log n)). Using `std::nth_element` (O(n)) would be faster for median computation, though this depends on the `Waveform::median` implementation.

The initial allocation of `signal` and `temp_signal` to full column width then resize is fine since `resize` doesn't reallocate when shrinking, but the pattern is repeated per row. Pre-allocating a single buffer outside the row loop and reusing it would avoid repeated allocation.

---

## EFF-8: `filter.at()` bounds checking in hot loop (Decon2DFilter.cxx:107-109, 116)

**Impact: Low**

```cpp
filter.at(i) *= wave.at(i);
// ...
c_data_afterfilter(irow, icol) = c_data(irow, icol) * filter.at(icol);
```

`.at()` performs bounds checking on every access. In tight inner loops over arrays whose bounds are already validated, `operator[]` avoids this overhead. The same applies to `m_overall_resp.at(i).at(j)` (Decon2DResponse.cxx:320), `elec.at(icol)` (line 166), `ch_elec.at(icol)` (line 301), etc.

---

## EFF-9: Per-channel FFT in electronics response correction (Decon2DResponse.cxx:293-309)

**Impact: Medium-High**

```cpp
for (int irow = 0; irow != c_data.rows(); irow++) {
    Waveform::realseq_t tch_resp = m_cresp->channel_response(ch_arr[irow]);
    tch_resp.resize(m_fft_nticks, 0);
    const WireCell::Waveform::compseq_t ch_elec = fwd_r2c(m_dft, tch_resp);
    // ... element-wise division
}
```

For every channel (potentially thousands), an FFT is computed independently. If many channels share the same response (which is common for channels on the same ASIC/FEMB), the FFT results could be cached by response identity. Even without caching, batching the FFTs into a single 2D FFT would be more efficient.

---

## EFF-10: Final forward FFT seems architecturally wasteful (Decon2DResponse.cxx:380)

**Impact: Medium (architectural)**

```cpp
c_data = fwd_r2c(m_dft, r_data, 1);
```

After the full inverse FFT pipeline (producing `r_data` in the time domain), the code immediately does a forward FFT again before outputting. This means the downstream `Decon2DFilter` receives frequency-domain data and must inverse-FFT it again. If the two components are always used together, the intermediate inverse+forward FFT pair is redundant -- the wire-domain frequency-domain data could be passed directly from Response to Filter, saving two full 2D FFTs.
