# Response & Convolution Group -- Code Examination

## Summary of Findings

This examination covers 8 source/header pairs implementing detector response
modeling: field response interpolation across impact positions, electronics
response (cold, warm, RC, JSON-based, systematic), and the 2D FFT-based
convolution engine that combines charge depositions with response functions.

**Bugs found: 4** (2 confirmed, 2 likely)

| ID | File | Severity | Description |
|----|------|----------|-------------|
| B1 | ResponseSys.cxx:12 | High | Constructor stores `tick` into `"nticks"` config key instead of `nticks` |
| B2 | WarmElecResponse.h:34 / .cxx:28 | Medium | Raw `new` without `delete`; memory leak on every `configure()` call |
| B3 | ResponseSys.h:34 / .cxx:28 | Medium | Raw `new` without `delete`; memory leak on every `configure()` call |
| B4 | PlaneImpactResponse.cxx:417 | Low | Off-by-one: `irind > m_ir.size()` should be `>=` |

**Efficiency concerns: 6**

| ID | File | Description |
|----|------|-------------|
| E1 | ImpactTransform.cxx | Massive code duplication between paired and central convolution paths |
| E2 | ImpactTransform.cxx:158-196 | Repeated inv-FFT + resize + fwd-FFT of response spectra per group; could be precomputed |
| E3 | PlaneImpactResponse.cxx:291-305 | FR waveform copy + resize + FFT per path (npaths can be 210+) |
| E4 | ImpactTransform.cxx:66 | Copies spectrum from ImpactResponse but never uses it (dead code) |
| E5 | PlaneImpactResponse.cxx:80 | `m_overall_short_padding` read from config before being initialized by constructor default |
| E6 | ImpactData.cxx:51 | No bounds check on `absbin` before indexing `m_waveform` |

---

## 1. PlaneImpactResponse

**Files:** `gen/inc/WireCellGen/PlaneImpactResponse.h`, `gen/src/PlaneImpactResponse.cxx`

### Algorithm overview

PlaneImpactResponse builds a collection of per-impact-position response
functions for one wire plane. It loads a field response (FR) from an
`IFieldResponse` component, resamples it from the FR's native time binning to
the digitization tick via frequency-domain interpolation (`spectrum_resize`),
convolves with "short" electronics responses (e.g., cold electronics), and
optionally prepares a separate "long" response (e.g., RC filter) for later
convolution. The responses are indexed by a 2D map (`m_bywire`) that maps
(wire, impact) pairs to `ImpactResponse` objects, exploiting the mirror symmetry
of the field response about each wire.

### Potential bugs

**B4 -- Off-by-one in bounds check (line 417)**

```cpp
if (irind < 0 || irind > (int) m_ir.size()) {
```

This should be `>=` not `>`. When `irind == m_ir.size()`, the subsequent
`m_ir[irind]` is an out-of-bounds access. In practice this is unlikely to be
triggered because index values come from `m_bywire` which is constructed from
valid path indices, but it is a latent bug.

**E5 -- Uninitialized member used as default (line 80)**

```cpp
m_overall_short_padding = get(cfg, "overall_short_padding", m_overall_short_padding);
```

`m_overall_short_padding` is not initialized in the constructor (line 27-34),
so the default value passed to `get()` is indeterminate. The same applies to
`m_long_padding` on line 81. In practice, `default_configuration()` sets these
values, and `configure()` is typically called with a config that includes them,
but if a config omits the key, behavior is undefined.

### Efficiency concerns

**E3 -- Per-path FFT pipeline (lines 291-349)**

For each of the ~210 paths, the code:
1. Copies the current waveform
2. Resizes it
3. Performs a forward FFT
4. Calls `spectrum_resize` (which does another FFT implicitly via `hermitian_mirror`)
5. Multiplies with the short response spectrum
6. Performs an inverse FFT
7. Resizes
8. Performs another forward FFT

This amounts to roughly 3 FFTs per path, or ~630 FFTs total. The short response
convolution (step 5) could be done once in a combined spectrum, but the
resampling and convolution are interleaved per-path, so this is inherent to the
algorithm. However, the final resize + re-FFT on lines 346-349 (going from
`n_short_length` to `m_nbins` and back to frequency domain) could be avoided if
downstream consumers accepted the shorter length.

**`spectrum_resize` intermediate allocation (lines 96-140)**

Each call allocates a new `compseq_t` of size `newsize` and copies data. For
210 paths, this is 210 temporary allocations. A reusable buffer would reduce
allocator pressure.

### Key algorithmic details

- The resampling from FR tick to digitization tick is done in the frequency
  domain. The `spectrum_resize` function truncates or zero-pads the Hermitian
  spectrum and applies a normalization factor. The `norm = rawresp_tick` choice
  (line 300) converts instantaneous current samples to integrated charge by
  multiplying by the FR sampling period.

- The `oopsilon` nudge (line 259) prevents floating-point `ceil()` from
  misclassifying wire numbers when pitch positions are near exact multiples of
  the pitch.

- Mirror symmetry construction (lines 360-371): for each relative wire, the
  impact indices are the "direct" set from `wire_to_ind[irelwire]` concatenated
  with the reverse of `wire_to_ind[-irelwire]` (skipping the first element to
  avoid double-counting the on-wire impact position).

---

## 2. ImpactTransform

**Files:** `gen/inc/WireCellGen/ImpactTransform.h`, `gen/src/ImpactTransform.cxx`

### Algorithm overview

ImpactTransform performs the core 2D convolution that turns diffused charge
distributions (from `BinnedDiffusion_transform`) into wire waveforms. It
exploits the fact that the field response is the same for symmetric impact
positions by packing two charge groups into the real and imaginary parts of a
complex array, performing a single 2D FFT, multiplying with the response in the
frequency domain, and inverse-transforming. The central (on-wire) impact group
is handled separately as a real-only convolution. After the 2D inverse FFT, the
real and imaginary parts are separated and summed with appropriate column
reversal to recover both contributions.

The `waveform()` method extracts per-wire results and optionally convolves with
the "long" auxiliary response (e.g., RC filter) in 1D.

### Potential bugs

**No confirmed bugs, but several robustness concerns:**

- **Line 130-131, 141, 227-228**: Array indexing into `c_data` and `data_t_w`
  uses charge tuple values directly offset by `npad_wire - start_ch` and
  `- m_start_tick`. There is no bounds checking. If a charge tuple has an
  out-of-range channel or tick, this will produce an out-of-bounds Eigen matrix
  access (undefined behavior). The correctness depends entirely on
  `BinnedDiffusion_transform` producing valid ranges.

- **Line 52-63**: When `m_pir->closest()` throws `ValueError`, the code catches
  it, prints to `cerr`, and does `continue`. This means `map_resp` will be
  missing the entry for that wire offset. Later, when the response is accessed
  via `m_vec_map_resp.at(i)[key]`, a missing key in `std::map` will
  default-construct a null `shared_ptr`, leading to a null dereference on
  line 158 or 176 when `->spectrum()` is called.

- **Line 66**: After the catch block, the code accesses
  `map_resp[j - m_num_pad_wire]->spectrum()` to copy the spectrum, but this
  result is never used. This is dead code (E4) but also means the successful
  path does an unnecessary copy.

### Efficiency concerns

**E1 -- Massive code duplication (lines 123-212 vs 217-301)**

The "paired groups" loop (lines 123-212) and the "central group" block
(lines 217-301) contain nearly identical response-spectrum construction logic
(the inv-FFT, resize, fwd-FFT pipeline for each wire's response). This code
is duplicated verbatim. Extracting a helper function would reduce ~90 lines
of duplication and make maintenance less error-prone.

**E2 -- Repeated response FFT reprocessing**

For each of the `num_double` groups plus the central group, the response
spectra are:
1. Inverse-FFT'd from stored spectrum to time domain
2. Truncated/zero-padded to the working size
3. Forward-FFT'd back to frequency domain

This is done for `m_num_pad_wire * 2 + 1` responses per group, and there are
`num_double + 1` groups. The responses at a given wire offset are often the
same across groups (they come from `m_vec_map_resp`), so this work could be
cached. Even within a single group, the FFT size is constant, so the resized
spectra could be precomputed once.

**Memory: large temporary Eigen arrays**

Each iteration of the paired-groups loop allocates two `array_xxc` matrices
of size `(nwires_padded x nticks_padded)`. For typical parameters (e.g.,
200 wires x 10000 ticks), each complex float matrix is ~16 MB. The `resp_f_w`
matrix is re-allocated inside the loop body. Moving it outside and using
`.setZero()` would avoid repeated allocation.

### Key algorithmic details

- **Real+imaginary packing trick (lines 130-147)**: Two charge groups that are
  symmetric about the wire center are packed as real and imaginary parts of a
  single complex array. After the 2D FFT and element-wise multiplication with
  the response, the inverse FFT produces results where the real part corresponds
  to one group and the imaginary part (with column reversal, line 305) to the
  other. This halves the number of 2D FFTs needed.

- **Line 305**: `acc_data_f_w.imag().colwise().reverse()` -- the column reversal
  is needed because the imaginary-part charge was stored with reversed wire
  indexing (line 141: `end_ch + npad_wire - 1 - get<0>(...)`).

---

## 3. ImpactData

**Files:** `gen/inc/WireCellGen/ImpactData.h`, `gen/src/ImpactData.cxx`

### Algorithm overview

ImpactData accumulates the charge waveform at a single impact position by
summing contributions from overlapping `GaussianDiffusion` patches. Each
diffusion has a 2D patch (pitch x time) and the ImpactData extracts the row
corresponding to its pitch bin. It also computes a "weight" waveform used for
interpolation between impact positions.

### Potential bugs

**E6 -- No bounds check on absbin (line 50-51)**

```cpp
const int absbin = tbin + toffset_bin;
m_waveform[absbin] += patch(pbin, tbin);
```

If `absbin < 0` or `absbin >= nticks`, this is an out-of-bounds write into
`m_waveform`. The `std::vector::operator[]` does not check bounds. The
correctness depends on `GaussianDiffusion` producing valid time offsets relative
to the allocated `nticks`. A defensive `at()` or bounds check would prevent
potential memory corruption.

### Efficiency concerns

- **Line 33**: `for (auto diff : m_diffusions)` copies each `shared_ptr`. Using
  `const auto&` would avoid the atomic reference count increment/decrement per
  diffusion.

- The `calculate()` method uses an early-return idempotency check (line 27-29)
  based on waveform size. This is fragile -- if `calculate()` is called with
  different `nticks` values, the second call silently returns stale data.

### Key algorithmic details

- The weight computation `qweight[pbin] * patch(pbin, tbin)` gives, for each
  time bin, a weighted measure of where charge is concentrated along the pitch
  direction. This enables linear interpolation between neighboring impact
  positions in the ImpactTransform.

---

## 4. ColdElecResponse

**Files:** `gen/inc/WireCellGen/ColdElecResponse.h`, `gen/src/ColdElecResponse.cxx`

### Algorithm overview

Implements the ICARUS cold electronics response function as an `IWaveform`.
Wraps `Response::ColdElec` (from WireCellUtil) which implements a 4th-order
semi-Gaussian shaping amplifier transfer function. Configured with gain,
shaping time, and an optional post-gain scaling factor.

### Potential bugs

None found. The implementation is clean and straightforward.

### Efficiency concerns

None significant. The waveform is generated once at configuration time and
cached. The `waveform_samples(Binning)` overload regenerates on each call,
which is appropriate for rebinning.

### Key algorithmic details

- Uses `std::unique_ptr<Response::ColdElec>` for proper RAII memory management.
  This is the correct pattern (contrast with WarmElecResponse and ResponseSys).

---

## 5. WarmElecResponse

**Files:** `gen/inc/WireCellGen/WarmElecResponse.h`, `gen/src/WarmElecResponse.cxx`

### Algorithm overview

Identical structure to ColdElecResponse but for ICARUS warm electronics. Wraps
`Response::WarmElec`.

### Potential bugs

**B2 -- Memory leak (line 28)**

```cpp
m_warmresp = new Response::WarmElec(m_cfg["gain"].asDouble(), m_cfg["shaping"].asDouble());
```

`m_warmresp` is a raw pointer (`Response::WarmElec*`, declared at line 34 of
the header). There is no `delete` in the destructor (no destructor is defined),
and no `delete` before reassignment if `configure()` is called multiple times.
Each call to `configure()` leaks the previous `Response::WarmElec` object.

**Fix:** Change `m_warmresp` to `std::unique_ptr<Response::WarmElec>` as done
in ColdElecResponse and RCResponse.

### Efficiency concerns

None beyond the memory leak.

---

## 6. RCResponse

**Files:** `gen/inc/WireCellGen/RCResponse.h`, `gen/src/RCResponse.cxx`

### Algorithm overview

Implements a simple RC (resistor-capacitor) filter response as an `IWaveform`.
Wraps `Response::SimpleRC`. Used as a "long" response in PlaneImpactResponse
for modeling long-timescale effects.

### Potential bugs

None found. Uses `std::unique_ptr` correctly.

### Efficiency concerns

None significant.

### Key algorithmic details

- The comment on line 22-23 ("fixme: why give SimpleRC tick twice?") suggests
  the API of `Response::SimpleRC` is slightly awkward -- the tick is passed both
  at construction and implicitly through the `Binning` in `generate()`. This is
  a design smell but not a bug.

---

## 7. JsonElecResponse

**Files:** `gen/inc/WireCellGen/JsonElecResponse.h`, `gen/src/JsonElecResponse.cxx`

### Algorithm overview

Loads an electronics response function from a JSON file containing `"times"`
and `"amplitudes"` arrays, and generates the waveform by linear interpolation
at the requested time bins.

### Potential bugs

**Null-check ordering (lines 59-66)**

```cpp
if (jtimes.size() != jamps.size()) {      // line 59
    THROW(...)
}
// ...
if (jtimes.isNull()) {                     // line 66
    THROW(...)
}
```

The null check on line 66 comes *after* `jtimes.size()` is already called on
line 59. If `jtimes` is null, calling `.size()` on a null Json::Value returns 0,
so the size check passes, and then the null check fires. This is not a crash
bug (Json::Value handles null gracefully), but the logic is inverted -- the null
check should come first for clarity and correctness.

**Unused variable (line 69)**

```cpp
const double tick = waveform_period();
```

This `tick` variable is declared but not used until line 77 where it is used
in the `Binning` constructor. Actually it is used, so this is fine. No bug.

**`lower_bound` edge case (line 27-29)**

```cpp
auto it = std::lower_bound(times.begin(), times.end(), time);
auto k = it - times.begin();
if (*it == time) ret.at(ind) = amps.at(k);
```

If `time` equals `times.back()`, `lower_bound` returns an iterator to the last
element, and the exact-match branch fires correctly. If `time` is slightly less
than `times.back()` but `lower_bound` returns `times.end()` (impossible since
we already checked `time <= times.back()`), there would be a problem. The
bounds check on line 24 prevents this, so this is safe.

However, there is an implicit assumption that `times` is sorted. If the JSON
file provides unsorted times, `lower_bound` produces incorrect results silently.
No validation of sort order is performed.

### Efficiency concerns

- The `generate()` function is a free function (not a method) in an anonymous
  namespace equivalent (it is at file scope). It performs a `lower_bound` search
  for each sample bin. For large waveforms, this is O(n log m) where n is the
  number of output bins and m is the number of input samples. Since the output
  bins are uniformly spaced and the input is sorted, a linear scan would be
  O(n + m).

---

## 8. ResponseSys

**Files:** `gen/inc/WireCellGen/ResponseSys.h`, `gen/src/ResponseSys.cxx`

### Algorithm overview

Generates a Gaussian-shaped systematic response function for studying field
response uncertainties. The Gaussian is parameterized by magnitude, time smear
(sigma), and time offset. Used as a convolution kernel to model systematic
distortions.

### Potential bugs

**B1 -- Constructor stores `tick` into `"nticks"` (line 12)**

```cpp
Gen::ResponseSys::ResponseSys(int nticks, double start, double tick, ...) {
    m_cfg["nticks"] = tick;   // BUG: should be nticks
    m_cfg["start"] = start;
    m_cfg["tick"] = tick;
```

The first config assignment uses `tick` (a `double`, the sampling period)
instead of `nticks` (an `int`, the number of ticks). This means
`default_configuration()` will return a config where `"nticks"` has the value
of `tick` (e.g., 250 ns = 2.5e-7 in system units), which when read as `int`
on line 30 would yield 0, producing an empty waveform.

In practice, if `configure()` is always called with a user-supplied config that
includes the correct `"nticks"`, this bug is masked. But if anyone relies on
the default configuration, the response will be broken.

**B3 -- Memory leak (line 28)**

```cpp
m_sysresp = new Response::SysResp(tick, m_cfg["magnitude"].asDouble(), sigma, offset);
```

`m_sysresp` is a raw pointer (`Response::SysResp*`, line 34 of header). No
destructor is defined, and no `delete` before reassignment. Same issue as
WarmElecResponse (B2).

### Efficiency concerns

None beyond the memory leak.

### Key algorithmic details

- The systematic response is a Gaussian convolution kernel:
  `magnitude * exp(-0.5 * ((t - offset) / sigma)^2)`. When convolved with the
  field response, it models timing uncertainties and gain variations.

---

## Cross-cutting Observations

### Inconsistent memory management patterns

Three different patterns are used across the response classes:
- `std::unique_ptr` (ColdElecResponse, RCResponse) -- correct
- Raw `new` without `delete` (WarmElecResponse, ResponseSys) -- memory leak
- No heap allocation needed (JsonElecResponse stores data directly) -- correct

The WarmElecResponse and ResponseSys should be updated to use `unique_ptr`
consistent with the other classes.

### Configuration pattern

All response classes follow the same configure pattern: store defaults in the
constructor via `m_cfg`, return `m_cfg` from `default_configuration()`, and
replace `m_cfg` entirely in `configure()`. The ResponseSys constructor bug (B1)
shows a weakness of this pattern -- typos in constructor initialization are hard
to catch since `configure()` typically overwrites everything.

### ImpactTransform complexity

The ImpactTransform constructor is approximately 300 lines of dense numerical
code with no helper functions. The real+imaginary packing optimization is
clever but makes the code fragile. The massive duplication between the paired
and central paths (E1) is a significant maintenance burden. Any fix to the
response-spectrum reprocessing logic must be applied in two places.
