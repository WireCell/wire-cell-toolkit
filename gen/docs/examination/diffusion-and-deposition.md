# Diffusion & Deposition Group -- Code Examination

Examination of 10 source/header pairs in `gen/` related to diffusion modeling
and deposition-to-frame conversion.

---

## Summary of findings

| Severity | Count | Description |
|----------|-------|-------------|
| Bug | 5 | Dimensional error in Diffuser sigmaT, division-by-zero paths, uninitialized offset bins, index overflow in get_charge_vec, DenseAccumulator tbin not updated |
| Logic | 4 | Floating-point equality in comparator, magic numbers in DepoSplat, early-return before bounds check, `break` vs `continue` |
| Efficiency | 6 | Unnecessary copies, redundant map lookups, repeated tbins construction, per-depo face search |
| Robustness | 3 | Raw `new` without RAII, cerr instead of logging, no null-depo check in DepoFramer |

---

## 1. GaussianDiffusion (.h + .cxx)

### Algorithm overview

`GausDesc` describes a 1D Gaussian (center, sigma). `GaussianDiffusion` forms a
2D charge patch by taking the outer product of two independent 1D Gaussian
distributions (time and pitch), then optionally applying binomial fluctuations
that preserve total charge.

The 1D distributions are computed via `binint()` which integrates erf() across
bin edges -- this is the correct approach for bin-integrated Gaussian sampling.

### Potential bugs

**B1. Floating-point equality comparison in `GausDiffTimeCompare` (line 249)**

```cpp
if (lhs->depo_time() != rhs->depo_time()) {
```

Comparing floating-point times with `!=` can produce non-deterministic ordering
when two depos have times that differ by floating-point rounding. This comparator
is used in `std::set`, so inconsistent ordering can violate the strict weak
ordering contract and cause undefined behavior. A tolerance-based comparison or
bit-exact comparison (via `std::memcmp` or integer reinterpretation) would be
safer.

**B2. `set_sampling` returns early before checking `npss` if `ntss == 0` (line 132)**

The check `if (!ntss)` returns early, but this happens *after* `m_toffset_bin`
has already been set (line 128). If later code queries `toffset_bin()` without
checking `patch().size()`, it will get a stale bin offset with no patch data.
The same pattern applies to the `npss` check at line 147.

**B3. Uninitialized `m_toffset_bin` / `m_poffset_bin` if `set_sampling` is never called**

They are initialized to `-1` in the constructor (line 111-112), which is a
sentinel value but is not documented as such. Calling `toffset_bin()` or
`poffset_bin()` before `set_sampling()` returns `-1`, which could silently
produce wrong indices downstream.

### Efficiency concerns

**E1. `weights()` returns by value (line 242)**

```cpp
const std::vector<double> Gen::GaussianDiffusion::weights() const { return m_qweights; }
```

This copies the vector every time. Should return `const std::vector<double>&`.

**E2. `GausDesc::weight()` takes `pvec` by value (line 61)**

```cpp
std::vector<double> weight(double start, double step, int nbins, std::vector<double> pvec) const;
```

Should be `const std::vector<double>&`.

### Key algorithmic details

- `binint()` uses the erf integral approach: for each bin `[a,b]`, the integrated
  Gaussian is `0.5 * erf(b') - 0.5 * erf(a')` where `x' = (x - center) / (sqrt(2) * sigma)`.
- `weight()` computes linear interpolation weights for redistributing charge
  to the two nearest impact positions. The formula involves the Gaussian PDF
  evaluated at bin edges and the integrated charge per bin.
- Binomial fluctuation (lines 186-213) approximates a multinomial distribution
  by applying independent binomial draws then renormalizing. The comment notes
  this introduces ~1% error at 10000 electrons.
- When `sigma == 0` (point source), all charge goes to a single bin.
  `GausDesc::distance()` returns raw distance (not in sigma) when sigma is zero,
  which is correct for the boundary-check logic in BinnedDiffusion.

---

## 2. BinnedDiffusion (.h + .cxx)

### Algorithm overview

Manages a collection of `GaussianDiffusion` objects, associating each with
impact position bins. Each depo is added as a GaussianDiffusion covering the
range of impact bins within `nsigma`. When `impact_data()` is called, it
lazily triggers `set_sampling()` on all diffusions for that bin and then
`calculate()` on the `ImpactData`.

### Potential bugs

**B4. Dangling reference to `m_dft` (line 101)**

```cpp
const IDFT::pointer& m_dft;
```

The member stores a *reference* to the `IDFT::pointer` passed to the
constructor. If the original shared_ptr goes out of scope or is reassigned,
this reference dangles. Should store by value (`IDFT::pointer m_dft`).

**B5. Dangling references to `m_pimpos` and `m_tbins` (lines 100, 102)**

```cpp
const Pimpos& m_pimpos;
const Binning& m_tbins;
```

Same issue: these store references to externally-owned objects. If the caller's
objects are destroyed before the BinnedDiffusion, these dangle. This is a design
risk, not necessarily a current bug if lifetimes are managed correctly by callers.

### Efficiency concerns

**E3. `pitch_range()` and `time_range()` rebuild a `vector<GausDesc>` every call (lines 157-163)**

Each call copies all GausDesc values from all diffusions into a temporary vector.
This is O(N) allocation + copy for N depos, and these functions may be called
multiple times.

**E4. `gausdesc_range()` takes vector by value (line 136)**

```cpp
static std::pair<double, double> gausdesc_range(const std::vector<Gen::GausDesc> gds, double nsigma)
```

Missing `&` -- copies the entire vector. Should be `const std::vector<Gen::GausDesc>& gds`.

### Key algorithmic details

- Uses a `std::set` with `GausDiffTimeCompare` ordering, which (per the
  floating-point comparison issue in B1) may cause subtle ordering differences
  across platforms.
- `impact_data()` triggers lazy evaluation: diffusions are sampled on first
  access. This means the `m_fluctuate` random number generator is consumed in
  access order, which may affect reproducibility if access order varies.

---

## 3. BinnedDiffusion_transform (.h + .cxx)

### Algorithm overview

A variant of BinnedDiffusion that accumulates all diffusions into a vector
(not per-impact-bin), then distributes charge into sparse matrices or charge
vectors via `get_charge_matrix()` / `get_charge_vec()`. This is the "transform"
path used by `DepoTransform` + `ImpactTransform`.

### Potential bugs

**B6. Hard-coded index hash overflow in `get_charge_vec()` (line 320)**

```cpp
long int index1 = channel * 100000 + abs_tbin;
```

If `abs_tbin >= 100000` (which is 100k ticks, possible at 0.5us tick over 50ms),
the hash collides with the next channel. This is a latent bug for long readouts.

**B7. `map_imp_redimp[abs_pbin] + 1` may be out of `map_redimp_vec` range (line 181, 296)**

In `get_charge_matrix()`:
```cpp
vec_spmatrix.at(map_redimp_vec[map_imp_redimp[abs_pbin] + 1])
```

If `abs_pbin` is the last impact of a wire, `map_imp_redimp[abs_pbin] + 1` may
not exist in `map_redimp_vec`, causing `map_redimp_vec[]` to default-construct
a zero entry. Then `vec_spmatrix.at(0)` may or may not be valid. Whether this is
exercised depends on the Pimpos layout.

**B8. Periodic clearing of lookup maps every 5000 iterations (line 363-377)**

```cpp
if (counter % 5000 == 0) {
    for (auto it = vec_map_pair_pos.begin(); it != vec_map_pair_pos.end(); it++) {
        it->clear();
    }
}
```

This clears the de-duplication maps, meaning that charges from depo N and
depo N+5001 that land in the same (channel, tick) cell will create *duplicate*
entries in `vec_vec_charge` instead of being accumulated. This is a correctness
issue: the output vector will have duplicate (channel, tick) tuples with partial
charges rather than a single accumulated value.

### Efficiency concerns

**E5. Multiple `std::map` lookups per (pbin, tbin) pair**

The nested maps `map_imp_ch`, `map_imp_redimp`, `map_redimp_vec` are all looked
up repeatedly in the inner loop. These could be precomputed into direct arrays
indexed by impact bin number.

**E6. `diff->clear_sampling()` frees patch memory inside the loop (line 208, 379)**

This is actually a *good* pattern for memory management -- it releases the patch
immediately after use. No issue here.

### Dangling references

Same issue as BinnedDiffusion: `m_pimpos` and `m_tbins` are stored as `const` references.

---

## 4. Diffuser (.h + .cxx)

### Algorithm overview

Older code that models longitudinal and transverse diffusion of drifted charge.
Takes IDepo objects, computes diffusion sigmas from drift distance and diffusion
coefficients (DL, DT), then creates an IDiffusion patch via erf-integrated 1D
Gaussians. Uses an internal buffer (`m_input`, an `IDiffusionSet`) to maintain
time ordering.

### Potential bugs

**B9. Dimensional error in `sigmaT` calculation (line 124)**

```cpp
const double sigmaT = sqrt(2 * m_DT * drift_time / units::centimeter2) * units::centimeter2;
```

The `sqrt()` produces units of `[cm]` (since `DT * time / cm^2` is dimensionless,
then `sqrt` gives dimensionless, then `* cm` gives cm). But the code multiplies
by `units::centimeter2` (cm^2) instead of `units::centimeter`, giving dimensions
of `[cm^2]` for what should be a length. Compare with `sigmaL`:

```cpp
const double sigmaL = sqrt(tmpcm2) * units::centimeter / m_drift_velocity;
```

Here `tmpcm2` is dimensionless and `* cm / velocity` gives time units, which is
correct for the longitudinal (time) direction.

For `sigmaT`, the result should be a length (pitch direction), so it should be:
```cpp
const double sigmaT = sqrt(2 * m_DT * drift_time / units::centimeter2) * units::centimeter;
```

**B10. `operator()` unconditionally prints to cerr (line 129)**

```cpp
cerr << "Diffuser: "
     << " drift distance=" << drift_distance << ...
```

This is not guarded by any debug flag and will produce output for every single
depo processed.

**B11. `oned()` produces zero-length vector if sigma is zero**

If `sigma == 0`, `oned()` will evaluate `0.5 * (x - mean) / 0` which is NaN/Inf.
The caller `diffuse()` checks for empty vectors but not NaN content.

### Efficiency concerns

**E7. `diffuse()` uses raw `new Diffusion(...)` (line 207)**

```cpp
Diffusion* smear = new Diffusion(depo, l_bins.size(), t_bins.size(), ...);
```

This is immediately wrapped in a shared_ptr (line 217), so no leak, but
`std::make_shared<Diffusion>(...)` would be better practice (single allocation).

### Key algorithmic details

- The buffer drain logic (lines 137-150) releases diffusions from the front of
  the ordered set only when the center of the newest diffusion is far enough
  ahead that no future diffusion can leapfrog the current leader. The threshold
  is `m_max_sigma_l * m_nsigma`.
- Diffusion sigmas follow the standard formula: `sigma = sqrt(2 * D * t)` where
  D is the diffusion coefficient and t is the drift time.

---

## 5. Diffusion (.h + .cxx)

### Algorithm overview

Simple 2D array wrapper over `boost::multi_array<double, 2>`. Stores a
diffusion patch with accessor methods for getting/setting values and computing
positions from indices.

### Potential bugs

**B12. `operator=` does not copy `lbin` and `tbin` (lines 26-37)**

```cpp
Diffusion& Diffusion::operator=(const Diffusion& other)
{
    m_depo = other.m_depo;
    array.resize(boost::extents[other.lsize()][other.tsize()]);
    array = other.array;
    lmin = other.lmin;
    tmin = other.tmin;
    lmax = other.lmax;
    tmax = other.tmax;
    return *this;
}
```

`lbin` and `tbin` are computed in the constructor from
`(lmax - lmin) / nlong` but are *not* copied in `operator=`. After assignment,
`lpos()` and `tpos()` will use the *old* `lbin`/`tbin` values of the target
object since they are recomputed only in the constructor. The copy constructor
(lines 16-24) also fails to copy `lbin` and `tbin`.

**Fix:** Add `lbin = other.lbin; tbin = other.tbin;` to both copy constructor
and `operator=`, or recompute them from the newly assigned lmin/lmax/nlong.

### Efficiency concerns

None significant -- this is a simple data container.

### Key algorithmic details

- `lpos(ind, 0.5)` and `tpos(ind, 0.5)` give bin centers due to the
  `offset` parameter.

---

## 6. DepoTransform (.h + .cxx)

### Algorithm overview

Converts an IDepoSet to an IFrame using `BinnedDiffusion_transform` and
`ImpactTransform`. For each face/plane combination, it bins all depositions,
convolves with the plane impact response, and extracts per-wire waveforms.

### Potential bugs

**B13. `m_pirs` indexed by `iplane` (local loop counter), not by plane index (line 213)**

```cpp
auto pir = m_pirs.at(iplane);
```

`iplane` is incremented for every plane in the face, starting from -1 then
pre-incremented. If `process_planes` skips a plane, `iplane` still increments
(it counts all planes, not just processed ones). So if planes [0,2] are
processed, `iplane` will be 0 for plane 0 and 2 for plane 2, and
`m_pirs.at(2)` is used for the collection plane. This works if `m_pirs` has
3 entries, but it breaks if `m_pirs` only has entries for the processed planes.

**B14. Redundant `tbins` declaration inside the plane loop (line 203)**

```cpp
Binning tbins(m_readout_time / m_tick, m_start_time, m_start_time + m_readout_time);
```

This shadows the identically-constructed `tbins` at line 183. Wasteful but
not a bug.

### Efficiency concerns

**E8. `depo->extent_long() / m_drift_speed` computed without caching**

For each depo, the sigma_time is `extent_long / drift_speed`. This division is
cheap but the pattern is duplicated across DepoTransform, DepoSplat, etc.

---

## 7. DepoSplat (.h + .cxx)

### Algorithm overview

"Splats" depositions directly into a frame, approximating combined
simulation + signal processing. Each depo is Gaussian-diffused and directly
added to channel waveforms without field/electronics response convolution.

### Potential bugs

**B15. Hard-coded magic numbers for additional smearing (lines 281-297)**

```cpp
const int time_offset = 2;  // # of ticks -- "why????"
```

```cpp
double add_sigma_L = 1.428249 * time_slice_width / nrebin / (m_tick / units::us);
```

```cpp
if (iplane == 0) add_sigma_T *= (0.402993 * 0.3);
else if (iplane == 1) add_sigma_T *= (0.402993 * 0.5);
else if (iplane == 2) add_sigma_T *= (0.188060 * 0.2);
```

These are unexplained detector-specific constants (likely MicroBooNE). They make
this component non-portable to other detectors. The code itself acknowledges this
with the comment "What are these magic numbers???" and "why????".

**B16. `time_offset = 2` shifts patch indexing (line 363)**

```cpp
auto icol = it - gd->toffset_bin() + time_offset;
```

This 2-tick offset is applied when reading from the patch but not when computing
the patch via `set_sampling()`. This means the patch is read with a systematic
2-tick shift, which is presumably compensating for some response delay, but it
could cause the last 2 columns of the patch to be read out-of-bounds (though
the `icol >= patch.cols()` check at line 364 guards against this).

**B17. `tend` clamped to `tbins.nbins() - 1` instead of `tbins.nbins()` (line 314)**

```cpp
const int tend = std::min(tbins.bin(tcen + twid) + 1, tbins.nbins() - 1);
```

This truncates the last valid bin. If a depo's charge falls in the last time
bin, it will be excluded. Should be `tbins.nbins()`.

**B18. `charge.resize(tend, 0.0)` should be `tend + 1` or use `tend` as exclusive bound (line 359)**

```cpp
if ((int) charge.size() < tend) {
    charge.resize(tend, 0.0);
}
```

Since the loop is `for (int it = tbeg; it < tend; ++it)` and `charge[it]` is
accessed, `charge` needs size `>= tend`. The resize to `tend` makes the valid
range `[0, tend)`, so `charge[tend-1]` is valid. This is correct.

**B19. Charge sign: `std::abs()` forces positive (line 365)**

```cpp
charge[it] += std::abs(patch(irow, icol) * charge_scaler);
```

The patch already has the depo charge baked in (which is typically negative for
ionization electrons). Taking `abs()` flips the sign convention. This is
intentional (matching signal processing output convention) but should be
documented.

### Efficiency concerns

**E9. Per-depo GaussianDiffusion allocation (line 334)**

Each depo creates a new `shared_ptr<GaussianDiffusion>`. For a splat operation,
a stack-allocated object (like DepoFluxSplat does) would avoid heap allocation.

---

## 8. DepoFluxSplat (.h + .cxx)

### Algorithm overview

A cleaner, more general version of DepoSplat. Converts IDepoSet to IFrame with
configurable smearing. Supports both sparse and dense output frames. No
detector-specific magic numbers -- smearing is configured via `smear_long` and
`smear_tran` parameters (in units of tick and pitch, respectively).

### Potential bugs

**B20. `break` instead of `continue` in pitch range check (line 352)**

```cpp
if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) {
    break;  // <--- breaks out of plane loop, skipping remaining planes
}
```

This `break` exits the plane loop for this depo, meaning if a depo is outside
the pitch range of plane 1, plane 2 will be skipped entirely. Should be
`continue` to skip only the current plane.

**B21. DenseAccumulator does not update `tbin` when extending trace (line 234-243)**

```cpp
if (need.first < have.first || need.second > have.second) {
    std::vector<float> newcharge(need.second-need.first, 0);
    std::copy(oldcharge.begin(), oldcharge.end(),
              newcharge.begin() + have.first-need.first);
    ...
    tp->charge() = newcharge;
}
```

When the range is extended to an earlier `tbin`, the trace's `tbin()` is not
updated to `need.first`. The `SimpleTrace` still reports the old starting tbin,
but the charge vector now has data starting at `need.first`. This causes
misalignment between `tbin()` and the charge array.

Additionally, when the ranges overlap but the new data does not extend the
existing range (the `else` branch is missing), the charge from the new patch
is silently dropped.

**B22. `eff_nsigma` uses `depo->extent_long()` instead of `sigma_L` (line 325)**

```cpp
double eff_nsigma = depo->extent_long() > 0 ? m_nsigma : 0;
```

After computing `sigma_L` (which may include extra smearing), the range check
still uses the raw `depo->extent_long()` to decide the effective nsigma. If
`extent_long() == 0` but `smear_long > 0`, then `eff_nsigma` is 0 and the check
becomes `nmin_sigma > 0 || nmax_sigma < 0`, which could reject depos that
actually have a valid smeared extent.

Same issue for the transverse direction at line 349:
```cpp
double eff_nsigma = depo->extent_tran() > 0 ? m_nsigma : 0;
```

### Efficiency concerns

**E10. `find_face()` iterates all faces for every depo (line 159-166)**

For N depos and F faces, this is O(N*F). Could sort depos by face first.

**E11. Stack-allocated `GaussianDiffusion` (line 357)**

```cpp
Gen::GaussianDiffusion gd(depo, time_desc, pitch_desc);
```

Good: avoids heap allocation, unlike DepoSplat.

### Key algorithmic details

- The `DenseAccumulator` pattern merges overlapping traces by computing the
  union of their tick ranges and using `std::transform` with `std::plus<float>`
  for accumulation.
- Nominal time is computed by shifting depo time to the collection plane:
  `nominal_depo_time = depo->time() + m_origin / m_speed`.

---

## 9. DepoPlaneX (.h + .cxx)

### Algorithm overview

Drifts depositions to a given X plane, maintaining them in causal time order.
Uses a working queue (sorted set) and a frozen queue (deque). Depos are "frozen
out" when new arrivals guarantee no future depo can arrive earlier.

### Potential bugs

None significant. The logic is straightforward.

### Efficiency concerns

**E12. `drain()` and `freezeout()` copy doomed depos then erase (lines 29-41, 43-53)**

```cpp
for (auto depo : m_queue) {
    if (depo->time() < time) {
        m_frozen.push_back(depo);
        doomed.push_back(depo);
    }
}
for (auto depo : doomed) {
    m_queue.erase(depo);
}
```

This iterates the entire queue twice and builds a temporary vector. Since
`m_queue` is a sorted set, a `lower_bound`-based range erase would be more
efficient: find the first element with `time >= time` and erase everything
before it in one sweep.

### Key algorithmic details

- `drain()` uses `depo->time() < time` (strict less-than), so depos exactly at
  the drain time are *not* frozen. This is consistent with the causal argument
  that depos at exactly the current time could still be reordered.
- `freezeout_time()` returns `-1.0 * units::second` when the frozen queue is
  empty, which is a sentinel value.

---

## 10. DepoFramer (.h + .cxx)

### Algorithm overview

Composes a drifter and a ductor: takes an IDepoSet, drifts all depos, then
feeds them through a ductor to produce frames, and sums the partial frames.

### Potential bugs

**B23. No null-check on `in` (line 33)**

```cpp
bool Gen::DepoFramer::operator()(const input_pointer& in, output_pointer& out)
{
    const int ident = in->ident();
```

If `in` is nullptr (EOS), this dereferences a null pointer. Unlike DepoTransform
and DepoFluxSplat, which check for null input, DepoFramer does not.

**B24. `partial_frames` may contain nullptr (line 63)**

The ductor may push a nullptr (EOS marker) onto the output queue. These nullptrs
are collected into `partial_frames` and passed to `Aux::sum()`. If `sum()` does
not handle nullptrs, this could crash.

**B25. Warning uses cerr instead of logging (line 52)**

```cpp
std::cerr << "Gen::DepoFramer: warning: failed to get null on last drifted depo\n";
```

Should use the logging framework.

### Efficiency concerns

**E13. Sorts depos by descending time then pushes nullptr (lines 39-40)**

```cpp
std::sort(depos.begin(), depos.end(), descending_time);
depos.push_back(nullptr);
```

This is fine for the drifter protocol (depos arrive in descending time order,
terminated by nullptr). No issue.

---

## Cross-cutting concerns

### Dangling references
BinnedDiffusion and BinnedDiffusion_transform both store references to
externally-owned `Pimpos` and `Binning` objects. This is safe only if the
owning code guarantees lifetime. In DepoTransform (line 205), the
`BinnedDiffusion_transform` is a local variable that uses `*pimpos` (a pointer
from the plane), which is safe as long as the plane outlives the local variable
(it does, since the plane comes from the anode face).

### Code duplication
`BinnedDiffusion` and `BinnedDiffusion_transform` share significant duplicated
code (the `gausdesc_range` static function is defined identically in both .cxx
files, as are `pitch_range`, `time_range`, `impact_bin_range`, `time_bin_range`).

### Sign convention
DepoSplat and DepoFluxSplat both force charge positive via `std::abs()`.
DepoTransform does not, relying on the downstream electronics response
convolution to handle the sign.
