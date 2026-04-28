# Detector-Specific Implementations Examination

Files examined:
- `src/Microboone.cxx` (1,352 lines)
- `src/Protodune.cxx` (1,020 lines)
- `src/ProtoduneHD.cxx` (1,043 lines)
- `src/ProtoduneVD.cxx` (1,354 lines)
- `src/DuneCrp.cxx` (369 lines)
- `src/Icarus.cxx` (110 lines)
- Corresponding headers

---

## Algorithm Overview

### Common Per-Channel Noise Filtering Pattern

All detectors follow the same general per-channel pipeline:

1. Nominal baseline subtraction
2. Gain correction
3. Detector-specific preprocessing (chirp, sticky codes, etc.)
4. FFT to frequency domain
5. RC undershoot correction (spectral division)
6. Spectral noise mask application
7. DC removal, IFFT
8. Robust baseline correction (sigma-clipping + median)
9. Adaptive baseline for partial RC channels
10. Noisy channel tagging (RMS-based)

### Coherent Noise Subtraction Pattern

Shared across MicroBooNE, PDHD, PDVD:

1. Compute median waveform across channel group
2. Signal protection: detect real signal regions (via deconvolution or amplitude
   thresholds) and exclude them from the median computation
3. Scaled subtraction: compute per-channel correlation coefficient with median,
   scale the median, and subtract

### Flag Encoding System

ADC values are overloaded to carry status flags:
- MicroBooNE (12-bit): signal flag = 4096, dead = 10000, protected = +20000
- PDHD/PDVD/DuneCrp (14-bit): signal flag = 16384, dead = 100000, protected = +200000

### Cross-Detector Comparison

| Feature | MicroBooNE | Protodune-SP | PDHD | PDVD | DuneCrp | ICARUS |
|---------|-----------|-------------|------|------|---------|--------|
| ADC bits (flag threshold) | 12 (4096) | 12 (via uB) | 14 (16384) | 14 (16384) | 14 (16384) | N/A |
| Sticky code mitigation | No | Yes | No | No | No | No |
| Chirp detection | Yes | No | No | No | No | No |
| Harmonic noise filter | No | Yes (coll.) | No | No | No | No |
| Per-channel freqmask | Yes | Yes (harm.) | Yes | Yes | No | No |
| RC undershoot correction | Yes | Yes | Commented out | Disabled | Yes | Yes |
| Adaptive baseline window | 20 | 20 (via uB) | 512 | 512 | 512 | No |
| Coherent noise sub | Yes | No | Yes | Yes | No | No |
| FEMB noise detection | No | No | Yes | Yes | No | No |
| Shield coupling sub | No | No | No | Yes | No | No |
| ADC bit shift correction | Yes | No | No | No | No | No |
| Low-freq noise detection | Yes | No | No | No | No | No |
| Relative gain calibration | No | Yes | No | No | No | No |
| Noisy channel tagging | Yes | Yes | Commented out | Yes | Commented out | No |

---

## Potential Bugs

### BUG-DET-1: Swapped min/max variable names in `LinearInterpSticky` (HIGH)
**File:** Protodune.cxx:343-346

```cpp
auto min = std::max_element(digits.begin(), digits.end());
auto max = std::min_element(digits.begin(), digits.end());
double max_value = *max;  // actually the minimum!
double min_value = *min;  // actually the maximum!
```
The iterator named `min` points to the max element, and vice versa. The subsequent
signal-like detection logic at lines 356-363 tests whether the **minimum** exceeds
the threshold (instead of the maximum), inverting the intended behavior. This causes
the algorithm to misclassify which sticky regions are signal-like vs. interpolatable.

### BUG-DET-2: Division by zero in `RawAdapativeBaselineAlg` (MEDIUM)
**Files:** Microboone.cxx:722, ProtoduneHD.cxx:620, ProtoduneVD.cxx:621, DuneCrp.cxx:184

```cpp
baselineVec[j] = ((j - downIndex) * baselineVec[downIndex] + (upIndex - j) * baselineVec[upIndex])
                 / ((double) upIndex - downIndex);
```
If `upIndex == downIndex` (both searches converge to the same index), divides by zero.

### BUG-DET-3: Hardcoded tick count 6000 in `LedgeIdentify` (MEDIUM)
**File:** Protodune.cxx:321,268,278

```cpp
if (LedgeEnd == 0) LedgeEnd = 6000;
if (LedgeStart < 5750) ...
if (height > 30 && LedgeStart < 5900) ...
```
Hardcoded for 6000-tick readout. The improved `LedgeIdentify1` (lines 37-199)
correctly uses `signal.size()`, indicating this is a known legacy issue.

### BUG-DET-4: ProtoDUNE-SP uses MicroBooNE's 4096 ADC flag threshold (MEDIUM)
**File:** Protodune.cxx:890-892

`Protodune::OneChannelNoise::apply()` calls `Microboone::SignalFilter`,
`Microboone::RawAdapativeBaselineAlg`, and `Microboone::RemoveFilterFlags` which
use 4096/10000/20000 constants. If ProtoDUNE-SP has 14-bit ADCs, real ADC values
above 4096 would be incorrectly treated as flagged.

### BUG-DET-5: Missing bounds check in harmonic noise removal (LOW)
**File:** Protodune.cxx:845-848

```cpp
spectrum.at(mag.size() + 1 - j).real(0);
```
When `j == 1`, accesses `mag.size()` which is out of bounds. The conjugate mirror
index should be `spectrum.size() - j`.

### BUG-DET-6: PDHD/PDVD `Subtract_WScaling` hardcodes threshold=4 (LOW)
**Files:** ProtoduneHD.cxx:87, ProtoduneVD.cxx:93

The Microboone version uses the configurable `correlation_threshold` parameter.
The PDHD/PDVD variants dropped it during copy-paste, hardcoding to 4.

### BUG-DET-7: PDVD `Subtract_WScaling` checks `ave_coef > 0` instead of `!= 0` (LOW)
**File:** ProtoduneVD.cxx:125

MicroBooNE and PDHD use `ave_coef != 0`. If all channels are anti-correlated with
the median, PDVD sets `scaling = 0` instead of computing the ratio.

### BUG-DET-8: Static warning flag shared across all instances (LOW)
**File:** Protodune.cxx:775

```cpp
static bool warned = false;
```
Only the first mismatched channel across the entire program lifetime gets a warning.
Masks configuration errors in multi-APA setups.

### BUG-DET-9: `RelGainCalib::apply` unchecked vector access (LOW)
**File:** Protodune.cxx:991

`m_rel_gain.at(ch)` uses channel number as direct index. If the JSON array is
smaller than the channel number, throws with no clear error message.

---

## Efficiency Issues

### EFF-DET-1: Per-channel FFT in all detectors (MEDIUM)
Every channel performs its own forward/inverse FFT. The DFT plan setup overhead
could be amortized across channels (plans are likely cached by the DFT
implementation, but the function call overhead remains).

### EFF-DET-2: Redundant waveform copy for baseline calculation (MEDIUM)
All detectors copy the entire waveform to compute a robust baseline:
```cpp
auto temp_signal = signal;  // full copy
```
A histogram-based or in-place approach would avoid the copy.

### EFF-DET-3: `CalcRMSWithFlags` incremental push_back (LOW)
**Files:** All detectors

```cpp
WireCell::Waveform::realseq_t temp;
for (...) if (sig.at(i) < 4096) temp.push_back(sig.at(i));
```
Growing a vector without `reserve()`. Pre-reserving would avoid reallocations.

### EFF-DET-4: ROIs iterated by value (LOW)
All detectors: `for (auto roi : rois)` copies each `vector<int>` ROI.
Should be `for (const auto& roi : rois)`.

### EFF-DET-5: PDVD per-bin temporary vectors in shield coupling (LOW)
**File:** ProtoduneVD.cxx:1159-1185

For each time bin, a new `temp` vector is created, filled, and passed to
`median_binned`. Thousands of short-lived vectors. Should pre-allocate.

### EFF-DET-6: PDVD `ShieldCouplingSub` makes 3 full passes (LOW)
**File:** ProtoduneVD.cxx:1255-1338

Pass 1: scale down by strip_length. Pass 2: subtract medians. Pass 3: scale back up.
Could be combined into fewer passes.

---

## Detector-Specific Algorithm Notes

### MicroBooNE Unique Features
- **ADC bit shift detection/correction**: Identifies lower bits stuck or shifted
- **Chirp noise detection** (`Diagnostics::Chirp`): Window-based RMS analysis
- **Partial RC undershoot detection** (`Diagnostics::Partial`): Spectral shape
- **Low-frequency noise identification** (`OneChannelStatus::ID_lf_noisy`)
- **freqmask schema (re-verified)**: `cfg/pgrapher/experiment/uboone/chndb-base.jsonnet`
  uses the legacy `{value, lobin, hibin}` bin-baked-at-jsonnet-time form.  This is
  safe because `Microboone::OneChannelNoise` uses a half-complex `fwd_r2c` FFT —
  only bins 0..nticks/2 (~4797) are consumed, so the positive-frequency notches at
  bins 169-173 and 513-516 (~35 kHz, ~107 kHz) hit correctly without explicit
  conjugate-mirror.  The intentional `daq.nticks=9595` vs `nf.nsamples=9592` offset
  was re-verified safe under `OmniChannelNoiseDB::set_nsamples()`: the auto-rebuild
  widens the spectrum to 9595, but the notch bin indices land at the same physical
  frequencies (±0.03%) and the consumer reads only the first 4798 bins regardless.
  New notches should use `wc.freqmasks_phys([freqs], delta)` from `wirecell.jsonnet`.

### SBND Unique Features
- **freqmask status**: `cfg/pgrapher/experiment/sbnd/chndb-base.jsonnet` has a
  `freqbinner.freqmasks(harmonic_freqs, ...)` call, but `harmonic_freqs == []`
  (all candidate frequencies commented out), so the effective mask is a pass-through.
  SBND's single production frame size (3400 ticks) matches `params.nf.nsamples`, so
  `OmniChannelNoiseDB::set_nsamples()` rebuilds to the same size — a no-op.  When
  harmonic notches are re-enabled after analysis, use `wc.freqmasks_phys([freqs], delta)`
  instead of `freqbinner.freqmasks`, which does not auto-mirror conjugate-frequency bins.

### Other experiments (freqmask audit summary)

The following experiments were audited against the runtime auto-rebuild and the
new `{value, flo, fhi}` freqmasks schema.  All instantiate `OmniChannelNoiseDB`,
so the `OmnibusNoiseFilter::set_nsamples()` dynamic_cast hits in every flow.
None require code or config-data migration; in-line comment-only updates were
applied where useful for future maintainers.

| Experiment | nticks/nsamples | freqmask state | Consumer | Notes |
|---|---|---|---|---|
| `dune10kt-1x2x6` | 6000/6000 (Reframer) | Legacy `{lobin,hibin}` U/V notches @ 169-173, 513-516 | `Protodune.cxx` (W-only) | U/V dead config (consumer gated on iplane==2); rebuild is no-op |
| `dune-vd` | tracks `daq.nticks` | inherits dune10kt configs; `single` filter commented out in `nf.jsonnet` | `Microboone.cxx` (inactive) | No active consumer |
| `dune-vd-coldbox` | 6000/6000 | empty | `Protodune.cxx` | No freqmasks; rebuild is no-op |
| `dunevd-crp2` | 6000/6000 | empty | `DuneCrp.cxx` (does not call `noise(ch)`) | No freqmask consumer at all |
| `icarus` | 4096/4096 | empty | `Icarus.cxx` (does not call `noise(ch)`); `Microboone::CoherentNoiseSub` reads `response()` only | No freqmask consumer |
| `iceberg` | 8256/8256 (splusn variants override `RUN_NTICKS`) | active empty; commented legacy blocks | `Microboone.cxx` (active, all-channel `noise(ch)` consumer) | Most relevant to the auto-rebuild — splusn variants align both nticks and nsamples to `RUN_NTICKS`; rebuild safety-net protects future drift |
| `pcbro-50liter` | 6000/6000 | Legacy `{lobin,hibin}` U/V notches | `Protodune.cxx` (W-only) | Same dead-U/V pattern as pdsp/dune10kt |

In-line comments documenting these invariants were added to each experiment's
`chndb-base.jsonnet`, recommending `wc.freqmasks_phys([freqs], delta)` if/when
notches are reactivated.

### ProtoDUNE-SP Unique Features
- **Sticky code mitigation**: Two-stage repair:
  1. Linear interpolation for non-signal-like sticky regions
  2. FFT-based interpolation using even/odd subsample prediction
- **Ledge artifact removal**: Step-function + exponential recovery detection
- **50 kHz harmonic removal**: Iterative (5 passes) statistical outlier detection
  in frequency domain, using adaptive baseline on magnitude spectrum
- **FEMB clock correction** (`FftScaling`): Frequency-domain resampling
- **Relative gain calibration**: Pulse-area-based gain correction from JSON
- **freqmask schema (re-verified)**: `cfg/pgrapher/experiment/pdsp/chndb-base.jsonnet`
  defines U/V plane notches at FFT bins 169-173 and 513-516 in the legacy
  `{value, lobin, hibin}` form.  pdsp has a single production frame size
  (`daq.nticks = nf.nsamples = 6000`, no Resampler upstream of NF), so the
  `OmniChannelNoiseDB::set_nsamples()` runtime rebuild is an idempotent
  no-op and bin indices stay physically correct.  Note that the current
  `pdOneChannelNoise` consumer (`Protodune.cxx`, gated on `iplane==2`) reads
  `m_noisedb->noise(ch)` only for W-plane channels, so the U/V freqmasks are
  effectively dead in today's code path — the rebuild does not regress
  anything because there is no live consumer.  When notches are reactivated,
  prefer `wc.freqmasks_phys([freqs], delta)` from `wirecell.jsonnet`.

### ProtoDUNE-HD Unique Features
- **FEMB negative pulse detection** (`FEMBNoiseSub`): Projects all channels in
  a group to 1D, finds wide ROIs below -3.5 sigma. Tags affected time ranges.
- **Wide adaptive baseline**: 512-tick window (vs 20 for MicroBooNE)
- **Per-channel frequency mask** (`ProtoduneHD.cxx`, `PDHD::OneChannelNoise::apply()`):
  Same mechanism as PDVD (see ProtoDUNE-VD section below).  Pre-existing
  U/V plane notches at FFT bins 169-173 (~57 kHz) and 513-516 (~171 kHz) in
  `cfg/pgrapher/experiment/pdhd/chndb-base.jsonnet` were never actually
  applied (no consumer); they have been cleared to `freqmasks: []` pending
  re-analysis (see in-line comments).  Toggle infrastructure is identical to
  PDVD: `use_freqmask` TLA on `pdhd/wct-nf-sp.jsonnet` (default `true`).
  Use `wc.freqmasks_phys([freqs], delta)` when populating new entries (see
  ProtoDUNE-VD section for the schema).  PDHD also emits a
  `PDHD OneChannelNoise ch=...: freqmask size N != FFT size M, mask SKIPPED`
  warning if the chndb spectrum size and runtime FFT size still mismatch
  after the auto-rebuild (defensive check; should never fire in practice).

### ProtoDUNE-VD Unique Features
- **Per-channel frequency mask** (`ProtoduneVD.cxx`, `PDVD::OneChannelNoise::apply()`):
  After the forward FFT, the waveform spectrum is element-wise multiplied by the
  per-channel `noise` spectrum read from the chndb (`m_noisedb->noise(ch)`).  The
  spectrum is populated from the `freqmasks` field in `chndb-base.jsonnet` via
  `OmniChannelNoiseDB::parse_freqmasks`.  An empty spectrum (the default) is a
  no-op.  Use `wc.freqmasks_phys([freqs], delta)` in `channel_info[]` entries
  to specify notch frequencies; the helper emits
  `{value: 0.0, flo: f-delta, fhi: f+delta}` records, and `parse_freqmasks`
  resolves them to bin indices at runtime from the live `m_tick`/`m_nsamples`
  *and* automatically zeros the conjugate-mirror bins.  This makes the same
  chndb work correctly for any frame size (eg PDVD has both 6400- and
  8000-tick raw frames in production).  Currently populated for the W-plane
  on anode 0 (channels 2188-2195 and 2480-2485, harmonics of 23.5 kHz from
  47 kHz to 282 kHz, half-width 1 kHz) — diagnosed from
  `magnify-run040475-evt0-anode0.root.rms.root`.  Controlled at run time with
  the `use_freqmask` TLA in `wct-nf-sp.jsonnet` (default `true`); pass
  `--tla-code use_freqmask=false` to disable.

  The legacy `{value, lobin, hibin}` form (used by sbnd, pdsp, uboone,
  iceberg, dune10kt-1x2x6, pcbro-50liter) is still accepted but does *not*
  auto-mirror — it trusts the JSON to enumerate both halves explicitly, and
  the bin indices are interpreted against the runtime `m_nsamples` as-is.

  Two related runtime safety nets work alongside the freqmask:
  1. **Auto-rebuild on size mismatch** — `OmnibusNoiseFilter` discovers the
     real frame size from the first input frame and pushes it to the chndb
     via `OmniChannelNoiseDB::set_nsamples`.  If it differs from the
     jsonnet-configured `nsamples`, all spectra (rcrc, config, noise/freqmasks,
     response) are rebuilt at the correct size and a single info-level log
     `OmniChannelNoiseDB: nsamples X -> Y, rebuilding spectra` is emitted.
  2. **Loud-skip warning** — if `noise(ch).size() != spectrum.size()` despite
     the rebuild (eg a non-`OmniChannelNoiseDB` chndb implementation that
     ignores the size handshake), the consumer logs
     `PDVD OneChannelNoise ch=...: freqmask size N != FFT size M, mask SKIPPED`
     instead of silently dropping the mask.
- **Shield coupling subtraction**: Novel noise removal for capacitive coupling
  between TDE U-plane strips and shield/grid:
  1. Scale each channel's signal by inverse strip length (capacitance weighting)
  2. Positive-only signal protection with wide padding (70 bins)
  3. Compute median excluding flagged/outlier bins
  4. Subtract median (fixed scaling=1)
  5. Scale back by strip length
