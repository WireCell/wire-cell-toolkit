# L1SPFilterPD — PDHD/PDVD Unipolar-Induction L1SP Component

This document covers the implementation of `L1SPFilterPD`, a PDHD/PDVD-specific
fork of `L1SPFilter` that handles **unipolar induction-plane signals** using a
per-ROI polarity-detecting LASSO (Strategy B).

For physical motivation, the three strategies (A/B/C), and the full uBooNE L1SP
algorithm background, see [`README.md`](./README.md).

---

## Relationship to `L1SPFilter`

`L1SPFilterPD` is a **full fork** — `L1SPFilter.{h,cxx}` are byte-identical
after this component was added.  No shared helpers were extracted from the
production file; any code that looks similar to `L1SPFilter.cxx` was
deliberately duplicated rather than refactored.

Key structural differences:

| | `L1SPFilter` | `L1SPFilterPD` |
|---|---|---|
| Physical problem | uBooNE shorted U/Y wires | PDHD/PDVD unipolar induction |
| Response basis | {collection W, induction V} | {bipolar, ±unipolar} polarity-selected |
| Polarity | always flag=1 (net-positive) | flag=+1 / −1 / 0 (Strategy B, per-ROI ratio) |
| Layer 4 cross-channel cleaning | ✓ (ch±1, ±3 ticks) | removed (shorted-wire only) |
| Propagation polarity tracking | `bool flag_shorted` | removed (flag=2 not ported) |
| Response pointers | raw `linterp<double>*` | `std::map<int, unique_ptr<linterp<double>>>` (per-plane) |
| Response init | builds from `IFieldResponse` on every `operator()` call | lazy (once per configure cycle); loads from pre-built JSON+bz2 via `Persist::load` |
| Response source | `IFieldResponse` + in-C++ FR⊗ER convolution | pre-built kernel file (`kernels_file`); generated offline by `wirecell-sigproc gen-l1sp-kernels` |

---

## Source files

```
sigproc/inc/WireCellSigProc/L1SPFilterPD.h
sigproc/src/L1SPFilterPD.cxx
```

Registered WCT component type: `L1SPFilterPD`
(implements `IFrameFilter` + `IConfigurable`; cxx:28–30).

Build: `sigproc/wscript_build` globs `src/*.cxx`, so no build-system edit was
needed.

---

## Algorithm structure

### Anonymous helpers (cxx, `namespace {}`)

**`build_G(nrow, ncol, row_offset, t_lo, t_hi, overall_offset, basis1_offset, scaling, resp_scale, basis0, basis1)`**

Builds an `nrow × 2·ncol` response matrix:
- Columns `[0, ncol)` — `basis0(dt + overall_offset)` at each `(meas-tick, signal-tick)` pair.
- Columns `[ncol, 2·ncol)` — `basis1(dt + overall_offset − basis1_offset)`.

`row_offset = (W_first_tick − β_first_tick)` decouples the W tick range
from the β tick range so the W vector can be **padded** by pad_L=30 ticks
before / pad_R=20 (positive) or 30 (negative) ticks after the β-coverage
span (matching the `(−15 µs, +10/+15 µs)` response window — see `t_hi`
below). Without that padding, the boundary β coefficients see only one
half of the kernel and can be inflated arbitrarily by the LASSO to soak
up unexplained ADC residuals.

`overall_offset` is the global LASSO frame origin (kernel-file metadata
`frame_origin_us`, in WCT time units): the kernel native time at which
"source signal at t = 0" sits in the LASSO frame.  By convention this is
the reference plane's bipolar zero-crossing (V plane for PDHD/PDVD), so
β at LASSO tick j ↔ "charge passed V plane at tick j".  A single global
value is used for both U and V channel fits; the U/V geometric arrival
difference is already encoded in each plane's kernel shape.

`basis1_offset` is the per-plane W shift loaded from the kernel file
(`unipolar_time_offset_us = zero_crossing − W_peak`).  Subtracting it puts
the W kernel peak at LASSO `dt = zero_crossing` (= bipolar zero-crossing),
i.e. inside the response window.  Negative-polarity case has
`basis1_offset = 0` (no shift): the `negative.unipolar` basis is
neg-half(bipolar), whose trough sits at native time +12 µs (PDVD V) /
+10 µs (PDHD V).  The caller widens `t_hi` to +15 µs for negative ROIs
so the trough is well inside the window; β then fires a few ticks
before the trough, recovering the physical "source signal arrives
before the ADC dip" timing.

The caller passes whichever `{basis0, basis1}` pair is appropriate for the
detected polarity (see `l1_fit` below).

**`lasso_solve(G, W, lambda, niter, eps)`**

Thin wrapper around `WireCell::LassoModel`: sets data, fits, returns `Getbeta()`.

### `init_resp()` (cxx:667–)

Lazy, idempotent.  On the first `operator()` call after construction or
reconfigure, loads `m_kernels_file` (resolved via `WIRECELL_PATH`) using
`Persist::load`.  For each plane entry in `"planes"`:
- `m_lin_bipolar[plane]` — bipolar kernel (positive case; shared with negative).
- `m_lin_pos_unipolar[plane]` — W-plane unipolar kernel (collection-on-induction).
- `m_lin_neg_unipolar[plane]` — neg-half(bipolar) (anode-induction).
- `m_unipolar_toff_pos[plane]` — W shift in WCT time units (positive case).
  Negative case has no shift; handled inline in `l1_fit()`.

After loading, every kernel sample is multiplied by `m_kernels_scale`
(= `kernels_scale` config key, default 1.0) before being handed to the
linterp.  The kernels are stored in ADC/electron at the reference 14 mV/fC
FE gain; `kernels_scale = params.elec.gain / (14 mV/fC)` corrects the
amplitudes to the runtime gain, keeping the LASSO bases consistent with
the actual raw-ADC data.

The kernel file is generated offline by `wirecell-sigproc gen-l1sp-kernels`
(see `wire-cell-python/wirecell/sigproc/l1sp.py`).  The PDHD file is
`pdhd_l1sp_kernels.json.bz2` (in `wire-cell-data`); it contains per-plane
(U=0, V=1) kernels for both positive and negative cases.

### `compute_asym()` (anonymous helper, cxx)

Single helper that returns an `AsymRecord` with the per-ROI features used
by the trigger, the adjacency-expansion pass, and the calibration dump.
A `bool fill_dump_fields` parameter gates the dump-only computations so
production (non-dump) runs only pay for what `decide_trigger()` and the
adjacency pass actually consume.

Always computed:

- `gmax` — drives `decide_trigger()` and the adjacency loose precondition.
- `raw_asym_wide` — wide-window raw asymmetry over
  `[roi_start − raw_asym_pad, roi_end + raw_asym_pad]`, gated by
  `raw_asym_eps`.  Used by the adjacency sign-aligned precondition.
- `sub_windows` (`std::vector<SubInfo>`) — every `|gauss| > core_g_thr`
  contiguous run inside the ROI with its `(run_len, abs_sum, fill, fwhm,
  aw)`.  Walked once by `enumerate_subwindows()` and consumed by both
  `decide_trigger()` (first-firing) and the best-by-score selection
  below.  `ef` (energy fraction) is intentionally **not** in `SubInfo`:
  it requires the wide ±`energy_pad` pad and is only consulted for runs
  that survive `min_length`, so `decide_trigger` computes it on demand.
- `core_length`, `core_raw_asym_wide` — selected by walking `sub_windows`
  and picking `argmax(run_len · |aw|)`; used by the adjacency loose
  precondition.

Dump-only (computed only when `fill_dump_fields == true`):

- In-ROI accumulators: `temp_sum`, `temp1_sum`, `temp2_sum`, `max_val`,
  `min_val`, `temp_sum_pos/neg`, `n_above_pos/neg`, `argmax/argmin_tick`,
  `sig_peak`, `sig_integral`, `gauss_abs_sum_roi`.
- Per-ROI gauss shape: `gauss_fill`, `gauss_fwhm_frac`.
- Wide-window decon energy fraction `roi_energy_frac` over
  `[roi_start − energy_pad, roi_end + energy_pad]`.
- Remaining best-sub-window fields: `core_lo`, `core_hi`, `core_fill`,
  `core_fwhm_frac`.

Definitions match `pdhd/nf_plot/find_long_decon_artifacts.py` (iter-7) so
Python ↔ C++ values can be compared bit-for-bit during validation.

### `decide_trigger()` (anonymous helper, cxx)

Walks the precomputed `AsymRecord::sub_windows` (no second pass over the
ROI ticks).  For each sub-window the multi-arm gate fires if all of:

- `gmax >= l1_gmax_min`
- `run_len >= l1_min_length`
- per-sub-window `ef >= l1_energy_frac_thr`
- ANY of the five arms:
  - `|aw| >= l1_asym_strong`, OR
  - `run_len >= l1_len_long_mod   AND |aw| >= l1_asym_mod`, OR
  - `run_len >= l1_len_long_loose AND |aw| >= l1_asym_loose`, OR
  - `run_len >= l1_len_fill_shape AND fill <= l1_fill_shape_fill_thr
                                  AND fwhm <= l1_fill_shape_fwhm_thr
                                  AND |aw| >= l1_asym_mod`, OR
  - `run_len >= l1_len_very_long  AND |aw| >= l1_asym_very_long`
    (very-long arm — OFF by C++ default; PDHD enables at `(140, 0.35)` via
    `cfg/pgrapher/experiment/pdhd/sp.jsonnet`).

A fired sub-window is then handed to `pdvd_track_veto_hits()` (default
no-op — see "PDVD-only opt-in track veto" below).  If the veto rejects
it, the loop continues to the next sub-window.

Polarity = `sign(aw)` of the first firing (and non-vetoed) sub-window.
Returns `{−1, 0, +1}`.

### PDVD-only opt-in track veto (`pdvd_track_veto_hits()`)

PDVD bottom anode 0 turns up four prolonged-track FPs that the five
trigger arms cannot cleanly suppress with threshold tuning alone:
they fire `asym_strong` at `|aw|` in 0.69–0.85 (matching real
high-asym TPs in the same range), or fire the length-only arms with
multi-channel spread.  See `pdvd/sp_plot/handscan_039324_anode0.csv`
(run 39324 events 0-5, anode 0) for the four reference cases.

The veto is a per-sub-window **post-trigger** rejection that adds
shape evidence to the existing arms.  It is OFF by C++ default (so
PDHD `decide_trigger` is bit-identical to before this change) and
enabled via the PDVD `data:` block (`l1_pdvd_track_veto_enable:
true`).  The check is:

```
|aw| < l1_pdvd_track_high_asym  AND  (
  run_len >= l1_pdvd_track_long_cl
  OR (run_len >= l1_pdvd_track_med_cl
      AND fill   >= l1_pdvd_track_med_fill
      AND fwhm   >= l1_pdvd_track_med_fwhm)
)
```

Defaults (see "Configuration parameters" below): `high_asym=0.85,
long_cl=170, med_cl=100, med_fill=0.40, med_fwhm=0.40`.  The first
arm catches long sub-windows whose `|aw|` is just above the trigger
floor; the second arm catches medium-length sub-windows that look
multi-peak (high fill+fwhm) at moderate asym.  Real L1SP unipolar
lobes either land above `0.85 |aw|` (escape) or stay below 100 ticks
(too short to trigger the veto).

Validated against `handscan_039324_anode0.csv`: 7 prolonged-track
fires suppressed (4 standalone clusters + 3 adjacency-coupled
sub-fires); 0 GT-positive hits lost (each cluster always has another
firing ROI that escapes the veto).  Aggregate per-ROI metrics on the
PDVD bottom anode 0 NPZ dumps moved from `F1=0.862` (4 jsonnet
threshold tweaks alone) to `F1=0.923` (with veto enabled).

### `l1_fit()` (cxx)

Called per ROI as the LASSO-applier; the trigger decision is precomputed
in `operator()` pass 2 (so each ROI's features are computed exactly once)
and passed in via the `polarity` argument.  Sequence:

1. **Basis selection** — `polarity = +1` → `{m_lin_bipolar, m_lin_pos_unipolar}`;
   `polarity = -1` → `{m_lin_bipolar, m_lin_neg_unipolar}`; `polarity = 0`
   → pass-through (no LASSO).
2. **Pass-through guard** — if the matching unipolar basis is null
   (config keys empty), early-return without modifying `newtrace`.
3. **Segmented LASSO solve** — segments of `l1_seg_length` ticks; for each
   segment call `build_G` then `lasso_solve`.  The `build_G` window is
   `(−15 µs, +10 µs)` for positive polarity and `(−15 µs, +15 µs)` for
   negative polarity (the wider upper edge admits the neg-half-bipolar
   trough at native +12 µs).  The W vector for each segment is loaded
   with pad_L=30 / pad_R=20 (positive) or 30 (negative) ticks of raw-ADC
   context around the segment's β span (clipped to the trace bounds) so
   boundary β coefficients have full kernel support; the padded raw ADC
   is fit context only and never written back. Combine the two basis
   coefficients and rescale to electron units in one step:
   `l1_signal[t] = (β₀ * basis0_scale + β₁ * basis1_scale) * scaling_factor`.
4. **Post-processing** — apply Gaussian smearing on `l1_signal` (which is
   already in electron units; the smearing kernel is sum-normalised to 1
   so the integral is preserved), then floor at `l1_decon_limit` (also in
   electrons), remove peaks below `peak_threshold`/`mean_threshold`, write
   back into `newtrace`.

The propagation/hint-polarity path (uBooNE Layer 3) and the `flag_l1 = 2`
zero-out branch are deliberately not ported — see Design decisions below.

### `operator()` (cxx)

1. Retrieve `adctag` and `sigtag` traces.  Sizes must match.
2. Build per-channel tick sets from decon signal (positive samples) and raw ADC
   (samples above the noise threshold, padded ±`raw_pad` ticks).
3. Merge and pad into ROI pairs per channel (±`roi_pad`).
4. **Pass 2 — decide.** For every in-scope ROI compute the `AsymRecord`
   features (with `fill_dump_fields = m_dump_mode`) and run
   `decide_trigger` once over the precomputed `sub_windows`.  Cache the
   per-ROI `(polarity, AsymRecord)` keyed by `(channel, roi_index)`.
5. **Pass 3 — adjacency expansion** (gated by `l1_adj_enable`, default
   ON; see "Cross-channel adjacency expansion" below).
6. **Pass 4 — apply.**  For each ROI with a non-zero post-adjacency
   polarity, call `l1_fit` with that polarity so the LASSO writeback
   honours the cached / promoted decision instead of recomputing it.
   Negative decon samples inside every ROI are zeroed in `newtrace`
   regardless of trigger result (matching `L1SPFilter` behaviour).
7. In `dump_mode`, write a per-frame NPZ with all ROI records, the
   pre-adjacency `flag_l1`, the post-adjacency `flag_l1_adj`, and
   `adj_donor_ch` (channel of the donor ROI, or −1 if none).
8. Emit `IFrame` with traces tagged `outtag`.

Layer 3 propagation and Layer 4 cross-channel cleaning (uBooNE shorted-wire
heuristic) are intentionally omitted for PDHD/PDVD — see Design decisions.
The adjacency-expansion pass is a *different* mechanism: it admits more
ROIs into LASSO based on a triggered neighbour, rather than zeroing one
out.

---

## Per-ROI trigger gate (Strategy B, retuned)

The legacy uBooNE single-ratio gate
(`temp_sum / (temp1_sum * rescaling / nbin_fit) > adc_ratio_threshold`) is
preserved in the dump record as `flag` / `ratio` for diagnostics, but no
longer drives the live `flag_l1`.  The retuned trigger uses six per-ROI
shape features computed per `|gauss| > l1_core_g_thr` sub-window, then
fires polarised on `sign(raw_asym_wide)`.

### Why per-sub-window and not per-ROI

The C++ ROI window (`gauss > 0` plus raw-ADC noise hits, padded
±`raw_pad` ticks) is wider than the iter-7 reference window
(`|gauss| > 50`).  A single C++ ROI can wrap several iter-7 candidates.
Computing features over the whole padded ROI dilutes asymmetry and length;
computing them per `|gauss| > l1_core_g_thr` run mirrors iter-7 and
recovers per-candidate granularity.

### Trigger logic

Fixed ROI preconditions: `gmax >= l1_gmax_min`.  Per sub-window
preconditions: `run_len >= l1_min_length`, `ef >= l1_energy_frac_thr`.
Then ANY of:

- **strong-asym arm**: `|raw_asym_wide| >= l1_asym_strong`
- **long-moderate arm**: `run_len >= l1_len_long_mod`
  AND `|raw_asym_wide| >= l1_asym_mod`
- **very-long-loose arm**: `run_len >= l1_len_long_loose`
  AND `|raw_asym_wide| >= l1_asym_loose`
- **fill-shape arm**: `run_len >= l1_len_fill_shape`
  AND `gauss_fill <= l1_fill_shape_fill_thr`
  AND `gauss_fwhm_frac <= l1_fill_shape_fwhm_thr`
  AND `|raw_asym_wide| >= l1_asym_mod`
- **very-long arm** (default OFF; PDHD enables at `(140, 0.35)`):
  `run_len >= l1_len_very_long` AND `|raw_asym_wide| >= l1_asym_very_long`.
  Targets long sub-windows whose per-sub-window asym sits below
  `l1_asym_mod` (0.40) but whose great length still indicates a unipolar
  artifact.  Calibrated 2026-05-02 against APA3 V ch 8753 in 027409:0
  (core_length=144, |craw|=0.36 — fails mod arm by 0.04).  Multi-event
  scan across 027409 events 0–7 measured ~1 new promotion/event globally
  on a baseline of ~74 — surgical enough to ship default-on for PDHD
  without rerunning the iter-7 calibration.  C++ default keeps the arm
  OFF (`l1_len_very_long = INT_MAX`, `l1_asym_very_long = 1.0`) so other
  experiments inherit prior behaviour bit-identically.

The four core arms map 1:1 onto the four iter-7 `cluster_pass()` rules in
`pdhd/nf_plot/find_long_decon_artifacts.py`.  Polarity = sign of the
firing sub-window's `raw_asym_wide` (more stable than the in-ROI ratio:
the wide window has a well-defined denominator gate via `l1_raw_asym_eps`
and looks at the surrounding context).

### Gain-scale convention

Raw-ADC trigger knobs (`l1_raw_asym_eps`, `raw_ROI_th_adclimit`,
`adc_sum_threshold`) are tuned at the 14 mV/fC reference FE gain.  The
PDHD jsonnet (`cfg/pgrapher/experiment/pdhd/sp.jsonnet`) multiplies them
by `gain_scale = params.elec.gain / (14.0 * wc.mV / wc.fC)` at configure
time, mirroring the same convention used for `adc_limit`, `min_rms_cut`,
and `max_rms_cut` in `chndb-base.jsonnet`.

Deconvolved-domain knobs (`l1_gmax_min`, `l1_core_g_thr`, all asym
ratios, lengths, energy fraction) operate on gain-normalised signals
and are gain-invariant — they are NOT scaled by `gain_scale`.

### Validation status

Tuned and verified against the iter-7 offline detector
(`pdhd/nf_plot/find_long_decon_artifacts.py`) on R=27409 evts 0–7,12,
U-plane APA 0–3:

| Metric | Value | Target |
|---|---|---|
| Recall vs iter-7 (clusters hit) | **90.0 %** | ≥ 90 % |
| Extras / cpp_fired (over-triggers) | **7.7 %** | ≤ 10 % |
| iter-7 clusters (reference)       | 230   | — |
| C++ fired ROIs                    | 432   | — |

Validators live next to the iter-7 detector:

- `pdhd/nf_plot/eval_l1sp_trigger.py` — compares C++ `flag_l1` vs a
  hand-scan CSV of (ch_lo,ch_hi,t_lo,t_hi) ground-truth boxes
  (`pdhd/nf_plot/handscan_27409.csv`).
- `pdhd/nf_plot/compare_trigger_vs_iter7.py` — compares C++ `flag_l1` vs
  iter-7 cluster CSVs over multi-event/multi-APA aggregates with
  `--show-misses` / `--show-extras` for spot-checks.

### Calibration inputs

**Kernel file**: `pdhd_l1sp_kernels.json.bz2` in `wire-cell-data` (PDHD).
Contains per-plane (U=0, V=1) bipolar + positive/negative unipolar kernels
with per-plane W shifts auto-derived from the zero-crossing calculation,
plus a top-level `meta.frame_origin_us` (= reference plane's bipolar
zero-crossing, V plane by default) used as the global LASSO frame origin.
Regenerate when the PDHD field response or electronics parameters change.
Use the `-d/--detector` preset (recommended; reads detector-specific
defaults from `wirecell/sigproc/track_response_defaults.jsonnet`):

```
wirecell-sigproc gen-l1sp-kernels -d pdhd        pdhd_l1sp_kernels.json.bz2
wirecell-sigproc gen-l1sp-kernels -d pdvd-bottom pdvd_bottom_l1sp_kernels.json.bz2
wirecell-sigproc gen-l1sp-kernels -d pdvd-top    pdvd_top_l1sp_kernels.json.bz2
```

`pdvd-top` uses a JsonElecResponse JSON.bz2 (cold-box characterised
electronics); `pdhd` and `pdvd-bottom` use the parametric cold model
with `gain` / `shaping`.  All three add an output-window pad on PDVD
(160 µs) so the bipolar induction tail does not wrap circularly.

The flat-flag form is still supported when overriding individual
parameters; the FR file moves to a `--fr` option:

```
wirecell-sigproc gen-l1sp-kernels \
    --fr dune-garfield-1d565.json.bz2 \
    --gain '14*mV/fC' --shaping '2.2*us' \
    --postgain 1.0 --adc-per-mv 11.7028571 \
    --coarse-time-offset '-8*us' \
    pdhd_l1sp_kernels.json.bz2
```

Validate the output with:
```
python3 pdhd/nf_plot/track_response_l1sp_kernels.py \
    --from-file pdhd_l1sp_kernels.json.bz2
```

---

## Configuration parameters

### Input / output tags

| Key | Default | Meaning |
|-----|---------|---------|
| `adctag` | `"raw"` | Input tag for raw ADC traces (post-NF) |
| `sigtag` | `"gauss"` | Input tag for decon signal traces (post-OmnibusSigProc) |
| `outtag` | `"l1sp"` | Output tag for corrected signal traces |

### Kernel file

| Key | Default | Meaning |
|-----|---------|---------|
| `kernels_file` | `""` | Path (resolved via `WIRECELL_PATH`) to the pre-built JSON+bz2 kernel file |
| `kernels_scale` | `1.0` | Amplitude multiplier applied to every loaded kernel sample. Set to `params.elec.gain / (14 mV/fC)` when the detector runs at a gain other than the reference 14 mV/fC. |

The kernel file is **required** — `init_resp()` throws if the key is empty.
Generate it offline with:
```
wirecell-sigproc gen-l1sp-kernels -d <detector>  <out>_l1sp_kernels.json.bz2
```
where `<detector>` is one of `pdhd`, `pdvd-bottom`, `pdvd-top`,
`uboone`, `sbnd` (presets in
`wirecell/sigproc/track_response_defaults.jsonnet`).
See `wire-cell-python/wirecell/sigproc/l1sp.py` for the schema.
The detector kernels live in `wire-cell-data/`:
`pdhd_l1sp_kernels.json.bz2`, `pdvd_top_l1sp_kernels.json.bz2`,
`pdvd_bottom_l1sp_kernels.json.bz2`.

### ROI building

| Key | Default | Meaning |
|-----|---------|---------|
| `roi_pad` | `3` | Tick padding added to each side of a decon ROI |
| `raw_pad` | `15` | Tick padding added around raw-ADC above-threshold samples |
| `raw_ROI_th_nsigma` | `4` | Raw-ADC threshold in units of estimated σ |
| `raw_ROI_th_adclimit` | `10` | Minimum absolute raw-ADC threshold (ADC counts) |

### Trigger thresholds — Strategy B retuned

Defaults seeded from the iter-7 offline detector
(`pdhd/nf_plot/find_long_decon_artifacts.py`) and validated on the C++
dump corpus (R=27409 evts 0–7,12 U-plane APA 0–3).  See "Per-ROI trigger
gate" above for how each knob is combined.

| Key | Default | Meaning |
|-----|---------|---------|
| `l1_min_length`              | `30`    | Min sub-window run length (ticks) |
| `l1_gmax_min`                | `1500`  | Min ROI peak `\|gauss\|` (electron units) |
| `l1_energy_frac_thr`         | `0.66`  | Min sub-window decon energy fraction (isolated-lobe gate) |
| `l1_energy_pad_ticks`        | `500`   | Wide-window pad for `roi_energy_frac` |
| `l1_raw_asym_pad_ticks`      | `20`    | Wide-window pad for `raw_asym_wide` |
| `l1_raw_asym_eps`            | `20.0`  | Per-tick raw-ADC gate (sign-routed) |
| `l1_core_g_thr`              | `50.0`  | Per-tick `\|gauss\|` gate defining sub-windows |
| `l1_asym_strong`             | `0.65`  | Strong-asym arm threshold |
| `l1_asym_mod`                | `0.40`  | Moderate-asym arm threshold |
| `l1_asym_loose`              | `0.30`  | Loose-asym arm threshold |
| `l1_len_long_mod`            | `100`   | Length needed to enable moderate-asym arm |
| `l1_len_long_loose`          | `200`   | Length needed to enable loose-asym arm |
| `l1_len_fill_shape`          | `50`    | Length needed to enable fill-shape arm |
| `l1_fill_shape_fill_thr`     | `0.38`  | `gauss_fill` ceiling for fill-shape arm |
| `l1_fill_shape_fwhm_thr`     | `0.30`  | `gauss_fwhm_frac` ceiling for fill-shape arm |
| `l1_len_very_long`           | `INT_MAX` | Length needed to enable very-long arm (OFF default) — PDHD overrides to `140` |
| `l1_asym_very_long`          | `1.0`   | Very-long arm asym threshold (OFF default) — PDHD overrides to `0.35` |
| `l1_pdvd_track_veto_enable`  | `false` | Master switch for the PDVD-only post-trigger track veto (see "PDVD-only opt-in track veto" above).  PDVD bottom anodes set this `true` in `cfg/pgrapher/experiment/protodunevd/sp.jsonnet`; PDHD never sets it, so PDHD's `decide_trigger` is bit-identical to the pre-veto behaviour |
| `l1_pdvd_track_high_asym`    | `0.85`  | `\|aw\|` escape: sub-windows above this never veto |
| `l1_pdvd_track_long_cl`      | `170`   | Long-arm: any sub-window with `run_len ≥ this` (and `\|aw\|` below escape) is rejected |
| `l1_pdvd_track_med_cl`       | `100`   | Shape-arm length floor |
| `l1_pdvd_track_med_fill`     | `0.40`  | Shape-arm `fill` floor (high fill suggests multi-peak track) |
| `l1_pdvd_track_med_fwhm`     | `0.40`  | Shape-arm `fwhm_frac` floor |

Legacy uBooNE knobs retained for diagnostics only (drive the `flag` /
`ratio` fields in the calibration dump but no longer affect `flag_l1`):

| Key | Default | Meaning |
|-----|---------|---------|
| `adc_l1_threshold`           | `6`     | Min `\|ADC\|` for in-ROI accumulators |
| `adc_sum_threshold`          | `160`   | Legacy Σ`\|ADC\|` floor (gain-scaled in pdhd jsonnet) |
| `adc_sum_rescaling`          | `90`    | Legacy ratio denominator |
| `adc_ratio_threshold`        | `0.2`   | Legacy asymmetry-ratio cut |

### LASSO solve

| Key | Default | Meaning |
|-----|---------|---------|
| `l1_seg_length` | `120` | Segment length (ticks) for the segmented solve |
| `l1_scaling_factor` | `500` | Numerical conditioning on G; cancels in linear algebra |
| `l1_lambda` | `10` | LASSO L1 regularization weight; per-coefficient sparsity threshold in electrons = `l1_lambda × l1_scaling_factor` (default 5000 e) |
| `l1_epsilon` | `0.05` | Convergence tolerance |
| `l1_niteration` | `100000` | Maximum LASSO iterations |
| `l1_resp_scale` | `1.0` | Kernel amplitude scale; must be 1.0 for ADC/electron kernels |
| `overall_time_offset` | `0` | Additive override on top of the kernel-file `frame_origin_us` (default 0; tuning only) |

Per-plane unipolar time offsets (positive case) are read from the kernel file
(`unipolar_time_offset_us` per plane); the negative case has no shift.

The global LASSO frame origin is loaded from kernel-file
`meta.frame_origin_us` (= reference plane's bipolar zero-crossing) and used
for *all* induction planes.  `overall_time_offset` is added on top as an
additive override (default 0); typical operation should leave it at 0.

### Output reconstruction

| Key | Default | Meaning |
|-----|---------|---------|
| `l1_decon_limit` | `100` | Floor (electrons) applied per-tick after smearing |
| `l1_basis0_scale` | `1.0` | Weight for the bipolar (basis0) component (β₀ already in electrons) |
| `l1_basis1_scale` | `1.0` | Weight for the unipolar (basis1) component (β₁ already in electrons) |
| `peak_threshold` | `1000` | Drop an output ROI if its peak < this |
| `mean_threshold` | `500` | Drop an output ROI if its mean < this |
| `filter` | `[]` | Explicit smearing kernel taps (overrides auto-derivation if non-empty) |

#### Smearing kernel — auto-derived from `Gaus_wide`

When `filter` is empty (the PDHD/PDVD default), the time-domain smearing
kernel is built once at `configure()` time from the SP `Gaus_wide` HfFilter:

1. Fetch `IFilterWaveform` named by `gauss_filter` (`"HfFilter:Gaus_wide"`).
2. Build the frequency-domain spectrum: `H[k] = exp(−½(|f_k|/σ)²)`, DC zeroed.
3. IFFT (1/N normalisation) → `h[n]`, peak at index 0.
4. Scan outward until `|h[k]| < kernel_threshold · h[0]` to find `n_half`.
5. Extract `kernel[i] = h[(i+N) % N]` for `i ∈ [−n_half, n_half]` and sum-normalise.

The IFFT bin spacing is `1/(2·max_freq)`. With `HfFilter`'s default
`max_freq = 1 MHz` this gives exactly 500 ns per bin, matching the SP tick on
both uBooNE and PDHD (post-resampler from 512 ns native).

The uBooNE explicit 21-tap `filter` array is numerically equivalent to the
IFFT-derived kernel (σ = 0.111408 MHz, 500 ns tick) to within max |Δ| ≈ 5×10⁻⁶.

**PDVD per-side filter**: the PDVD jsonnet registers `Gaus_wide_b`
(bottom CRP, `anode.data.ident < 4`) and `Gaus_wide_t` (top CRP) as
separate instances, and `make_sigproc` sets `gauss_filter` to
`'HfFilter:Gaus_wide_b'` or `'HfFilter:Gaus_wide_t'` per anode.  PDHD has
a single `Gaus_wide` shared across all four APAs.

| Key | Default | Meaning |
|-----|---------|---------|
| `gauss_filter` | `"HfFilter:Gaus_wide"` | Type-name of the frequency-domain filter to IFFT |
| `kernel_threshold` | `1e-3` | Truncation threshold as a fraction of the peak |
| `kernel_max_half` | `64` | Maximum half-width of the derived kernel (ticks) |
| `kernel_nticks` | `4096` | FFT size used for the IFFT derivation |

Set `gauss_filter` to `""` (and leave `filter` empty) to disable smearing entirely.

### DFT

| Key | Default | Meaning |
|-----|---------|---------|
| `dft` | `"FftwDFT"` | IDFT component type-name (used for smearing kernel derivation) |

Electronics response parameters (gain, shaping, postgain, ADC_mV,
coarse/fine time offsets) are no longer C++ config keys — they are baked
into the kernel file at build time by `gen-l1sp-kernels`.

---

## Run-to-run determinism

`L1SPFilterPD` is bit-deterministic for a fixed input frame (verified
2026-05-02 on PDHD event 027409:0, APA0 and APA1, two consecutive
runs each with `setarch -R`).  All `frame_gauss{N}` / `frame_wiener{N}`
arrays compared bit-identical.  The component has no RNG, no
threading, no hash-keyed containers (only `std::map<int,…>` /
`std::set<int>` keyed by channel/tick), constructs a fresh
`LassoModel` per ROI (no carry-over), and rides on the
`FFTW_ESTIMATE | FFTW_UNALIGNED` plans established by
`aux/src/FftwDFT.cxx` for upstream NF/SP determinism.  Eigen runs
single-threaded (no OpenMP/TBB/MKL linked).

---

## Performance profile

Measured 2026-05-02 with gperftools (`libprofiler` LD_PRELOAD,
`CPUPROFILE_FREQUENCY=400`) on a single-event NF+SP run of PDHD event
027409:0 APA1 (U+V planes both processed by L1SP).  Total 4953
samples ≈ 12 s of CPU.  Reported as percent of total wall-clock to
calibrate where future optimisation effort is worth spending.

| Component                          | %total | Notes                                                                          |
|------------------------------------|--------|--------------------------------------------------------------------------------|
| `OmnibusSigProc::operator()`       | 40.7%  | dominated by FFTs (`FftwDFT::inv1b` 12.9%, `fwd1b` 7.8%) inside `decon_2D_*`   |
| `FrameFileSink` output             | 22.7%  | `BZ2_compressBlock` 17.1% — pure I/O cost, codec-bound                         |
| `Main::initialize`                 | 16.7%  | one-time JSON / FieldResponse load (`Persist::load` 10.9%)                     |
| `OmnibusNoiseFilter::operator()`   | 9.4%   | NF stage                                                                       |
| `FrameFileSource` input            | 5.5%   | bzip2 decompression of input frames                                            |
| **`L1SPFilterPD::operator()`**     | **1.6%** | breakdown below                                                              |

Within `L1SPFilterPD::operator()` (78 samples = 1.6% of total):

| Sub-cost                           | % of L1SP | % of total |
|------------------------------------|-----------|------------|
| `std::nth_element` (noise percentile in operator() setup) | 33%  | 0.52% |
| `l1_fit` total                     | 22%       | 0.34%      |
|  └ `LassoModel::Fit`               | 19%       | 0.30%      |
| `std::set::insert` (`init_map` per-tick fills) | 21%  | 0.32%      |
| `compute_asym` + `enumerate_subwindows`        | 2.6% | 0.04%      |

### What this tells us

L1SPFilterPD is a small slice of total runtime.  The shipped
optimisations (single-copy noise percentile, walk-once `SubInfo`
refactor, dump-only field gating, plane cache, single sigtraces
sweep) extracted the bulk of the addressable budget here.

Items that were considered but **deliberately dropped** after the
profile because the addressable improvement is too small to justify:

- `init_map` `std::set<int>` → bitmap.  Targets the 0.32% set-insert
  bucket; saves ≤ 0.25% of total runtime even in the best case.
- Consolidating `compute_asym` pass 1 with the sub-window walk.
  Targets a 0.04% bucket; effectively zero return.
- Replacing `.at(idx)` with `[idx]` in the per-ROI inner loops.
  Same 0.04% ceiling; trades the bounds check for a sub-millisecond
  win.

If a future pass wants real wall-clock improvements, the targets
are outside this component:

- **`OmnibusSigProc` (40.7%)** — FFTW-bound; FFT-plan caching,
  batched DFT calls, or replacing `inv_c2r`/`fwd_r2c` with a
  faster path are the levers.  See follow-up below.
- **`FrameFileSink` bzip2 compression (17.1%)** — swapping bzip2
  for zstd or lz4 would likely save 10-15% of total wall-clock
  outright.  Reader-side coordination required.
- **Output tag duplication** — L1SP feeds the corrected gauss
  trace under both the `gauss%d` and `wiener%d` output tags
  (`cfg/pgrapher/experiment/pdhd/sp.jsonnet` `final_merger`).
  Dropping the wiener copy if downstream no longer consumes it
  would cut both output size and `FrameFileSink` time.

Profile artifacts (raw `.prof` + `pprof --pdf` callgraph) live
under `/home/xqian/tmp/l1sp_efficiency/` on the dev host.

### Follow-up: `OmnibusSigProc` Tier 1 + 2 (2026-05-02)

Two FFT-side cleanups landed against the 40.7% OmnibusSigProc
budget identified above.  Neither touched signal-processing
algorithm; both are pure plumbing.

**Tier 1 — eliminate redundant 2D-array copies in `DftTools`
(commit `4e37cd83`).**  `DftTools::fwd`/`inv` (axis form) at
`aux/src/DftTools.cxx:67,76` allocated a fresh complex array and
copied the input into it before running the FFT in-place.  Added
`fwd_inplace`/`inv_inplace` overloads that skip the copy when the
caller is about to overwrite the input anyway, and rewired
`fwd_r2c(2D)`/`inv_c2r(2D)` to use them internally (drops a second
copy that the cast-then-`fwd` chain was incurring).  Three
overwrite call sites in `OmnibusSigProc::decon_2D_init` (lines
1070/1083/1107) switched to the in-place overload.  Bit-identical
(same FFTW plan signature → same cached plan → same butterflies);
verified 41/41 dump-NPZ arrays equal via `np.array_equal` on
027409:0 APA0+APA1.

**Tier 2 — batch the per-channel electronics FFT loop (commit
`42705188`).**  `decon_2D_init` ran ~960 separate 1D r2c FFTs per
plane via the per-channel correction loop at `cxx:1050-1066`.
Replaced with a single batched `fwd_r2c(2D)` on a stacked response
matrix, followed by a per-row vectorized correction.  The FFTW plan
signature changes (howmany goes from 1 to nchans), so this is *not*
guaranteed bit-identical in general — but on the verification run
across all four PDHD APAs of 027409:0 the dump-mode L1SP NPZ
outputs were 41/41 bit-identical and `flag_l1` / `flag_l1_adj` sums
matched exactly: FFTW_ESTIMATE picks the same butterfly path for
both decompositions at these sizes.

**Wall-clock impact** (APA1 NF+SP, 3-run averages on the dev host):

| Build                       | Wall (s)   | Δ vs pre-Tier-1 |
|-----------------------------|------------|------------------|
| Pre-Tier-1 (`bc7817f7`)     | 19.73 ± 0.10 | —              |
| Tier 1 alone (`4e37cd83`)   | 18.96 ± 0.05 | −0.77 s (−3.9%)|
| Tier 1 + Tier 2 (`42705188`) | 18.83 ± 0.08 | −0.90 s (−4.6%)|

Tier 1 captured most of the gain; Tier 2 added another ~0.7%.

**Item considered and dropped — field-response frequency-domain
cache.**  The plan called for caching `c_resp` per plane to avoid
the 2D field-response FFT inside every `decon_2D_init`.
Invariance check fails: `overall_resp[plane]` is mutated in-place
at `OmnibusSigProc.cxx:1653` by the filter-response multiplication
*before* the first `decon_2D_init` call (pass 1), and the second
pass (APA0 only, line 1818) reloads it fresh.  Within-event caching
therefore needs state tracking, and cross-event caching depends on
`m_fft_nticks` derived from input frame trace size (cxx:802).
Combined complexity outweighs the ~2-3% potential gain.

**Items NOT pursued (out of scope for "low-hanging fruit"):**

- **r2c at the `IDFT` interface level.**  Using a true r2c plan
  instead of zero-imag c2c saves ~half the butterfly work on real
  inputs and would be the largest single FFT-side win, but
  requires extending `IDFT` (`iface/inc/WireCellIface/IDFT.h`)
  with new method families and updating `FftwDFT` plus any other
  backends.  Future architectural work.
- **Wire-axis transpose-then-FFT.**  FFTW handles strided data
  internally; not a measured bottleneck.
- **Output-tag duplication (FrameFileSink budget).**  Confirmed
  not an FFT lever — all OmnibusSigProc output tags share the
  same time-domain `m_r_data[plane]`; tags are labels on the same
  trace, no extra IFFTs.  The remaining FrameFileSink savings
  come from codec/IO work (zstd swap), which is a separate
  effort.

---

## Cross-channel adjacency expansion

Long unipolar artifacts have spatial coherence: when a neighbouring
channel clearly fires L1SP and a sub-threshold ROI on the current
channel time-overlaps it with similar extent, that ROI is almost
certainly the same artifact and should be processed.  The Strategy B
multi-arm gate is computed per-sub-window inside `decide_trigger` and
can miss such ROIs when the per-sub-window asymmetry is diluted, the
core sub-window is shorter than `l1_min_length`, or the ROI's `gmax`
is just below `l1_gmax_min`.

To recover those candidates, `operator()` runs a pass (default ON
since 2026-05-02) between the trigger decision and the LASSO apply
that *promotes* a ROI's polarity to that of an adjacent neighbour
ROI when the criteria below are met.  The expansion is **iterative**
(BFS layer by layer, default since 2026-05-03): originally-triggered
ROIs are hop 0, direct ±1 neighbours hop 1, neighbours-of-neighbours
hop 2, and so on, capped at `l1_adj_max_hops` (default 3 ⇔ ±3
channels from any original donor).  The loose preconditions are
re-applied independently at every hop, so a transitive promotion
only happens when each link in the chain looks unipolar in its own
right.  Within one hop, all promotions are computed from the donor
state at the *end of the previous hop* (no in-layer chaining), which
makes the result independent of map iteration order.

Set `l1_adj_enable=false` to disable the pass entirely (pre-2026-05-02
behaviour); set `l1_adj_max_hops=1` to recover the pre-2026-05-03
single-hop behaviour (donors must be originally-triggered).

### Promotion criteria

For candidate ROI `R_c` on channel `c` (with `hop == −1`) and donor
ROI `R_d` on channel `c±1` (with `0 ≤ donor_hop < current_hop`):

| Knob | Default | Check |
|------|---------|-------|
| `l1_adj_overlap_pad` | 3 ticks | `(R_c.start − pad) ≤ R_d.end + pad` AND `(R_c.end + pad) ≥ R_d.start − pad` |
| `l1_adj_gap_max` | 100 ticks | `|R_c.start − R_d.start| ≤ gap_max` (sanity check) |
| `l1_adj_max_hops` | 3 | maximum BFS depth from any original trigger |
| `l1_adj_len_ratio` | 0.40 | `min(len_c, len_d) / max(len_c, len_d) ≥ ratio` |
| `l1_adj_loose_gmax` | 300 (gain-normalised) | `R_c.gmax ≥ this` |
| `l1_adj_loose_core_len` | 2 ticks | `R_c.core_length ≥ this` |
| `l1_adj_loose_asym_abs` | 0.30 | `max(aw_aligned, craw_aligned) ≥ this`, where `aw_aligned = sign(polarity_d) · raw_asym_wide` and likewise for `craw_aligned`. The candidate must have at least one wide- or core-window asymmetry pointing in the donor's polarity direction with magnitude ≥ this. |

The loose preconditions on the candidate itself prevent noise ROIs
from being swept up — they apply at *every* hop, so a chain can only
extend through ROIs that each independently look unipolar.  Defaults
were first tuned on event 027409:0 APA0 (the original screenshot
reproducer); the iterative cap was added after event 027409:0 APA1
showed a 4-channel chain (ch 4071 ⇆ 4072 ⇆ 4073) where the second
hop was needed to catch ch 4071.  The asymmetry precondition was
loosened from `|core_raw_asym_wide|` to a sign-aligned
`max(aw_aligned, craw_aligned)` after event 027409:0 APA3 U-plane
showed ch 8354 sandwiched between two originally-triggered chs
(8353, 8355) with `aw=+0.45` (matching donor sign) but
`craw=+0.28` (just below the legacy 0.30 cut).  The sign-alignment
also tightens the gate against wrong-sign promotions: in the same
revision, three APA0 ev0 channels (223, 588, 590) that the legacy
sign-agnostic check accepted are now rejected because their
strongest asymmetry indicator points opposite to the donor's
polarity (their LASSO fits would have been wrong-signed).

### Honest framing on prior art

`L1SPFilter`'s `cxx:432-476` neighbour-channel pass is *not* the
inspiration: that path is a shorted-wire suppression heuristic that
zeroes a flag-2 seed ROI when an adjacent flag-1 ROI overlaps in
time, which does not apply to PDHD/PDVD geometry.  The borrowing
from MicroBooNE here is limited to the ±3-tick overlap convention.

### Validation

**Single-hop (2026-05-02), event 027409:0 APA0.**  With adjacency
disabled, the long unipolar tail on channel 324 (ticks ~5830-5945)
was passed through unchanged with peak ~4750 ADC.  With adjacency
enabled, the tail beyond the leading peak is zeroed by the LASSO
fit; only the genuine ~5830-5840 collection-induction lobe survives.
The number of triggered ROIs across APA0 increases from 8 → 17.

**Iterative expansion (2026-05-03), event 027409:0 APA1 V-plane.**
Channels 4073 and 4074 trigger originally; ±1 hop catches 4072
(donor 4073) and 4075 (donor 4074); the second hop catches 4071
(donor 4072), which has `|core_raw_asym_wide|=0.75`, length 154 vs
donor 180 (ratio 0.86), and meets every loose precondition.  Across
the whole event APA1: 24 originals, 20 hop-1 adjacency promotions
under the old behaviour (`l1_adj_max_hops=1`), 35 promotions under
the new default — i.e. iterative adds 15 hop-2/3 promotions on top
of the 20 hop-1 ones, all gated by the per-candidate loose
preconditions.

Setting `l1_adj_enable=false` reverts to the pre-2026-05-02 output,
bit-identical at the inner-`.npy` level.  Setting `l1_adj_max_hops=1`
reverts to the pre-2026-05-03 single-hop output.

### Toggling

The pass is on by default with `l1_adj_max_hops=3`.  To disable it
entirely (pre-2026-05-02 output):

```bash
wire-cell --tla-code l1sp_pd_adj_enable=false ... -c wct-nf-sp.jsonnet
```

To restrict to a single hop (pre-2026-05-03 output):

```bash
wire-cell --tla-code l1sp_pd_adj_max_hops=1 ... -c wct-nf-sp.jsonnet
```

Or in `pgrapher/experiment/pdhd/sp.jsonnet`:

```jsonnet
sp.make_sigproc(anode, l1sp_pd_adj_enable=false)
sp.make_sigproc(anode, l1sp_pd_adj_max_hops=1)
```

The threshold knobs above are exposed via `L1SPFilterPD`'s direct
config (`l1_adj_overlap_pad`, etc.) and can be overridden through the
component's `data:` block in jsonnet.

---

## Calibration dump schema

When `dump_mode=true` a single NPZ file per frame is written to `dump_path`.
All per-ROI arrays are parallel (same length = total ROI count across all
in-scope channels).

### Frame scalars

| Key | Type | Meaning |
|-----|------|---------|
| `frame_ident` | int32 | Frame identifier |
| `frame_time` | float64 | Frame timestamp |
| `call_count` | int32 | `operator()` invocation index |
| `n_rois` | int32 | Total number of ROIs in this frame |

### Per-ROI locator

| Key | Type | Meaning |
|-----|------|---------|
| `channel` | int32\[\] | Channel ID |
| `roi_start` | int32\[\] | First tick of the ROI (inclusive) |
| `roi_end` | int32\[\] | Last tick of the ROI (inclusive) |
| `nbin_fit` | int32\[\] | ROI width in ticks (`roi_end - roi_start + 1`) |

### Per-ROI asymmetry scalars (threshold-gated, `|adc| > adc_l1_threshold`)

| Key | Type | Meaning |
|-----|------|---------|
| `temp_sum` | float64\[\] | Σ ADC (signed) |
| `temp1_sum` | float64\[\] | Σ \|ADC\| |
| `temp2_sum` | float64\[\] | Σ \|gauss\| (sig absolute sum, gated) |
| `max_val` | float64\[\] | Max ADC in ROI (ungated) |
| `min_val` | float64\[\] | Min ADC in ROI (ungated) |

### Per-ROI same-channel adjacency

| Key | Type | Meaning |
|-----|------|---------|
| `prev_roi_end` | int32\[\] | `roi_end` of preceding ROI on same channel (−1 if none) |
| `next_roi_start` | int32\[\] | `roi_start` of following ROI on same channel (−1 if none) |
| `prev_gap` | int32\[\] | Gap in ticks to preceding ROI (−1 if none) |
| `next_gap` | int32\[\] | Gap in ticks to following ROI (−1 if none) |

### Per-ROI polarity classification (Tier 1)

| Key | Type | Meaning |
|-----|------|---------|
| `flag` | int32\[\] | `{0, +1, -1}` under current threshold config |
| `ratio` | float64\[\] | `temp_sum / (temp1_sum * adc_sum_rescaling / nbin_fit)`; 0 when `temp1_sum==0` |
| `temp_sum_pos` | float64\[\] | Σ ADC for positive thresholded samples only |
| `temp_sum_neg` | float64\[\] | Σ ADC for negative thresholded samples only (≤ 0) |
| `n_above_pos` | int32\[\] | Count of samples with `adc > +adc_l1_threshold` |
| `n_above_neg` | int32\[\] | Count of samples with `adc < −adc_l1_threshold` |

A pure bipolar signal has `temp_sum_pos ≈ −temp_sum_neg`; a pure unipolar signal
has one of them ≈ 0.  `n_above_pos / n_above_neg` carries similar information as
a sample-count ratio.

### Per-ROI peak location and decon scalars (Tier 2)

| Key | Type | Meaning |
|-----|------|---------|
| `argmax_tick` | int32\[\] | Absolute tick of `max_val` within the ROI |
| `argmin_tick` | int32\[\] | Absolute tick of `min_val` within the ROI |
| `sig_peak` | float64\[\] | Max of gauss (decon) signal within the ROI, ungated |
| `sig_integral` | float64\[\] | Σ gauss within the ROI, ungated (signed sum) |

`sig_peak` and `sig_integral` are independent of the ADC threshold gate — they
capture the full decon content without the raw-ADC bias.  `argmax_tick` /
`argmin_tick` are useful for `unipolar_time_offset` calibration once unipolar
field-response files are available.

### Per-ROI Strategy-B features (Tier 3)

Computed by `compute_asym()` and matching the iter-7 detector's feature
definitions bit-for-bit.  Used by `decide_trigger()` and as offline
analysis inputs.

| Key | Type | Meaning |
|-----|------|---------|
| `gmax`                | float64\[\] | `max(\|gauss[t]\|)` for `t ∈ [roi_start, roi_end]` |
| `gauss_fill`          | float64\[\] | `Σ\|gauss\| / (gmax · nbin_fit)` (full ROI) |
| `gauss_fwhm_frac`     | float64\[\] | `count(\|gauss\| > 0.5·gmax) / nbin_fit` (full ROI) |
| `roi_energy_frac`     | float64\[\] | `Σ\|gauss\|_ROI / Σ\|gauss\|_(ROI±l1_energy_pad_ticks)` |
| `raw_asym_wide`       | float64\[\] | `(pos+neg)/(pos−neg)` over raw ADC in ROI±`l1_raw_asym_pad_ticks`, gated by `±l1_raw_asym_eps` |
| `core_lo`             | int32\[\]   | First tick of longest `\|gauss\|>l1_core_g_thr` sub-window (−1 if none) |
| `core_hi`             | int32\[\]   | Last tick of that sub-window (−1 if none) |
| `core_length`         | int32\[\]   | `core_hi − core_lo + 1` |
| `core_fill`           | float64\[\] | `gauss_fill` recomputed on the core sub-window |
| `core_fwhm_frac`      | float64\[\] | `gauss_fwhm_frac` recomputed on the core sub-window |
| `core_raw_asym_wide`  | float64\[\] | `raw_asym_wide` recomputed around the core sub-window |
| `flag_l1`             | int32\[\]   | Pre-adjacency `decide_trigger()` result `{−1, 0, +1}` under current config |
| `flag_l1_adj`         | int32\[\]   | Post-adjacency polarity actually used to drive the LASSO branch.  With `l1_adj_enable=true` (the default) this can differ from `flag_l1`; equals `flag_l1` when `l1_adj_enable=false` |
| `adj_donor_ch`        | int32\[\]   | Channel of the adjacent donor ROI when an originally `flag_l1==0` ROI was promoted, else −1 |

`flag_l1_adj` is the trigger that drove the LASSO branch for this ROI.
Compare against `flag_l1` to see ROIs the cross-channel adjacency pass
recovered, and against legacy `flag` to see how the new gate diverges
from the uBooNE single-ratio decision.

---

## Wiring into a PDHD/PDVD graph

The wiring pattern mirrors uBooNE (`cfg/pgrapher/experiment/uboone/sp.jsonnet:65–163`).
The main differences are:

- `ChannelSelector` is optional; Strategy B does not require a static channel list.
  If the artifact is geometrically localized, a selector may still help reduce CPU.
- The `FrameMerger` feeding L1SP must supply both `adctag` and `sigtag` on the
  same frame.
- The output merger writes `outtag` back over `gauss` (and optionally `wiener`).

Skeleton (jsonnet):

```jsonnet
local l1spfilterpd = g.pnode({
    type: "L1SPFilterPD",
    data: {
        dft: wc.tn(tools.dft),
        kernels_file:  "pdhd_l1sp_kernels.json.bz2",  // resolved via WIRECELL_PATH
        kernels_scale: params.elec.gain / (14.0 * wc.mV / wc.fC),
        adctag: "raw%d" % n,    // post-NF raw ADC
        sigtag: "gauss%d" % n,  // post-OmnibusSigProc decon
        outtag: "gauss%d" % n,
        process_planes: [0, 1],  // 0=U, 1=V; skip W
    }
}, nin=1, nout=1, uses=[tools.dft, anode]),
```

See `cfg/pgrapher/experiment/pdhd/sp.jsonnet` for the live PDHD wiring
(the `l1sp_pd_mode != ''` branch) and the per-APA `process_planes` defaults
(APA0 → `[0]` U only; APA1-3 → `[0, 1]` U+V).

### PDVD wiring

`cfg/pgrapher/experiment/protodunevd/sp.jsonnet` mirrors the PDHD wiring with
two structural differences driven by PDVD's dual electronics:

| Knob | Bottom (anodes 0–3) | Top (anodes 4–7) |
|------|---------------------|------------------|
| `kernels_file` | `pdvd_bottom_l1sp_kernels.json.bz2` | `pdvd_top_l1sp_kernels.json.bz2` |
| `kernels_scale` | `params.elec.gain / (7.8 mV/fC)` (runtime scalar gain knob) | `1.0` (fixed `JsonElecResponse`, no runtime gain knob) |
| `gauss_filter` | `HfFilter:Gaus_wide_b` | `HfFilter:Gaus_wide_t` |
| ADC threshold reference | 7.8 mV/fC bottom electronics | top JsonElecResponse + 1.52 postgain |

Both kernel files live in `wire-cell-data/` and resolve via `WIRECELL_PATH`.
Generated offline with `wirecell-sigproc gen-l1sp-kernels -d pdvd-bottom <out>.json.bz2`
(and `-d pdvd-top` for the top region).

The dispatch primitive is `local sfx = if anode.data.ident < 4 then '_b' else '_t';` —
the same suffix used elsewhere in the PDVD jsonnet (e.g. `nf.jsonnet:24`,
`sp.jsonnet:24`) for region-specific filter / channel-noise-DB / electronics-response
selection.

**Default mode is `'dump'`**, not `'process'` as on PDHD: the ROI tagger runs
and writes per-event NPZ records, but `l1_fit()` is bypassed and the LASSO
output is not written back into `gauss`/`wiener`.  This lets the user
validate the tagger and the per-region kernel files independently before
turning on full replacement.  To support this, `operator()` guards
`init_resp()` with `if (!m_dump_mode)` so the kernel JSON is not loaded
until the user opts into process mode (added 2026-05-03).  Switching to
`'process'` is a runtime decision (`-w <wf_dir>` on
`pdvd/run_nf_sp_evt.sh`); no jsonnet change is needed since the
`kernels_file` paths are already wired.

PDHD's `l1_len_very_long=140` / `l1_asym_very_long=0.35` overrides are **not**
applied on PDVD; the very-long arm is left at the C++ default (OFF) until
PDVD-side calibration justifies enabling it.

PDVD bottom anodes (`anode.data.ident < 4`) additionally override four
trigger thresholds and enable the PDVD-only track veto, tuned against
`pdvd/sp_plot/handscan_039324_anode0.csv`:

| Knob | C++ default | PDVD-bottom override |
|---|---|---|
| `l1_len_long_mod` | 100 | 180 |
| `l1_len_fill_shape` | 50 | 90 |
| `l1_fill_shape_fill_thr` | 0.38 | 0.30 |
| `l1_fill_shape_fwhm_thr` | 0.30 | 0.25 |
| `l1_pdvd_track_veto_enable` | false | **true** |
| `l1_pdvd_track_{high_asym,long_cl,med_cl,med_fill,med_fwhm}` | (defaults) | left at defaults (0.85, 170, 100, 0.40, 0.40) |

PDVD top anodes (`anode.data.ident >= 4`) inherit all C++ defaults
pending their own hand-scan validation.

The runtime entry point (`pdvd/run_nf_sp_evt.sh`) provides flags that mirror
PDHD's: `-c <calib_dir>` overrides the auto-generated calibration dump
location, `-w <wf_dir>` switches to process mode + per-ROI waveform dump
(requires kernels), `-x` disables L1SP entirely.

---

## Test coverage

Added 2026-05-02 alongside the PDHD NF+SP chain.  All sit under
`sigproc/test/` and are auto-discovered by `./wcb build --tests`.

### Atomic doctest tests (no external fixtures)

- `doctest_dfttools_inplace.cxx` — pins the `Aux::DftTools::fwd_inplace` /
  `inv_inplace` overloads added in commit 4e37cd83 (Tier 1 FFT path).
  Bit-equal vs the allocating overload, both axes, both directions.
- `doctest_framemerger.cxx` — `replace` and `include` rules, per-trace
  summary propagation, channel-mask union (commit 7d05a7cb).
- `doctest_framesplitter.cxx` — broadcast + EOS pass-through.
- `doctest_femb_noise.cxx` — `PDHD::Is_FEMB_noise` width gate, padding,
  `nsigma` threshold, polarity (negative-only), boundary clamping
  (commit 139d782f).
- `doctest_l1sp_structural.cxx` — config schema, lazy `init_resp()`
  (configure does not throw without a kernel), explicit `ValueError`
  when `kernels_file` is empty.

### Kernel-gated doctest (auto-skips if `WIRECELL_PATH` cannot resolve
the production kernel)

- `doctest_l1sp_kernel.cxx`:
  - Loads `pdhd_l1sp_kernels.json.bz2`, runs `operator()` against a
    minimal synthetic `raw`+`gauss` frame, asserts no exception and
    that the output carries the L1SP output tag.
  - Pins the `dump_mode` NPZ schema: every documented per-frame
    scalar (`frame_ident`, `frame_time`, `call_count`, `n_rois`)
    and per-ROI Tier 1/2/3 column (Tier 3 includes `gmax`,
    `gauss_fill`, `gauss_fwhm_frac`, `roi_energy_frac`,
    `raw_asym_wide`, `core_*`, `flag_l1`, `flag_l1_adj`,
    `adj_donor_ch`).  Schema-pin only; no value assertions.

### End-to-end APA1 regression (BATS)

`check_pdhd_apa1_nf_sp.bats` runs the production
`wct-nf-sp.jsonnet` graph on the APA1 raw frame
(`run027409` / `evt 0`) and compares the SP output (`gauss1`,
`wiener1` tags) against the in-tree fixture
`sigproc/test/data/protodunehd-sp-frames-anode1.tar.bz2`.
Tolerances: ≤1% per-channel RMS on ≥99% of channels, ≤0.1% per-tag
absolute integral, soft `np.allclose(atol=1.0_ADC, rtol=1e-3)`.
Observed run-to-run spread is bit-zero on this branch (commit
47d16673 fixed the FFTW plan-cache nondeterminism), so the tolerances
are conservative against future drift.

The reference fixture was generated from current HEAD on 2026-05-02
after the L1SP→`wiener1` routing change in commits 616cfcb2 and
eab4e6ac.  Refresh it (and update `WCT_PDHD_REF` if you don't want
to overwrite the in-tree copy) only when a deliberate algorithmic
change is being baselined.

Skip semantics: the BATS test skips cleanly when any of
`WCT_PDHD_DATA` (raw input), `WCT_PDHD_DEPLOY` (jsonnet entry
point), or the in-tree reference fixture is absent.

### What's deliberately not tested at unit level

- **`OmnibusSigProc` Tier 1 + 2 batched FFT path** — already
  paired-run verified bit-identical on all 4 APAs (commits 4e37cd83,
  42705188).  Covered empirically by the BATS test.
- **L1SP Strategy-B trigger branches** (synthetic polarity /
  adjacency) — modelling the trigger gates well enough to fabricate
  inputs that hit each branch is brittle and high-cost.  The BATS
  test on real data is strictly stronger evidence and already
  caught the wiener-routing change on first run.
- **`CoherentNoiseDump` round-trip** — diagnostic-only emitter, not
  in the data flow; format is not a stable contract.

## Pending work

1. **LASSO body tuning** — calibrate `l1_basis0_scale`, `l1_basis1_scale`
   and verify the auto-derived smearing kernel against post-fit reconstruction.
   The trigger gate (trigger thresholds, dump validation) is independent and
   does not need to be re-run.  Kernel regeneration (if PDHD electronics
   parameters change) uses `wirecell-sigproc gen-l1sp-kernels`; the new file
   must be placed in `wire-cell-data` and its name updated in
   `cfg/.../pdhd/sp.jsonnet:kernels_file`.

2. **PDVD jsonnet port** — the trigger logic is detector-symmetric; wire
   `L1SPFilterPD` into `cfg/pgrapher/experiment/protodunevd/sp.jsonnet`
   with a PDVD-specific `kernels_file` (generate with the PDVD field response).
   The same `gain_scale` block applies.

## Design decisions

- **`flag=2` (zero-out) not ported** — uBooNE's flag=2 branch zeros ROIs with
  very low ADC content or "phantom decon, flat ADC" signatures.  Both conditions
  are specific to the uBooNE shorted-wire pathology and have no analog in PDHD/PDVD.
  PDHD/PDVD uses only `{0, +1, -1}` flags.

- **Propagation layers not ported** — uBooNE's forward/reverse same-channel sweep
  (rescues flag=2 ROIs near a confident fit) and `ch±1` cross-channel cleanup
  (suppresses phantom signals on shorted neighbors) were removed.  Both exist
  exclusively because of shorted U/Y wires.  The Layer 4 comment in `operator()`
  documents the cross-channel decision.  BUG-L1-1 (reverse-scan index mismatch in
  the propagation loops) is moot since those loops no longer exist.
