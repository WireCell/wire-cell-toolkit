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
| Response pointers | raw `linterp<double>*` | `unique_ptr<linterp<double>>` |
| Response init | `init_resp()` at every `operator()` call | lazy (once per configure cycle) |

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

**`build_G(nbin, t_lo, t_hi, overall_offset, basis1_offset, scaling, resp_scale, basis0, basis1)`**

Builds the N×2N response matrix for one segment:
- Columns `[0, N)` — `basis0` evaluated at each `(meas-tick, signal-tick)` pair.
- Columns `[N, 2N)` — `basis1` evaluated at the same pairs, shifted by `basis1_offset`.

The caller passes whichever `{basis0, basis1}` pair is appropriate for the
detected polarity (see `l1_fit` below).

**`lasso_solve(G, W, lambda, niter, eps)`**

Thin wrapper around `WireCell::LassoModel`: sets data, fits, returns `Getbeta()`.

### `build_response(fields_tn, plane_index)` (cxx:225–251)

Loads the named `IFieldResponse`, averages it across all wires in `plane_index`,
convolves with `Response::ColdElec` at the configured gain/shaping, and returns
a `unique_ptr<linterp<double>>`.  Called by `init_resp()`.

### `init_resp()` (cxx:254–265)

Lazy, idempotent.  On the first `operator()` call after construction or
reconfigure:
- Always builds `m_lin_bipolar` from `fields` / `bipolar_plane`.
- Builds `m_lin_pos_unipolar` from `fields_pos_unipolar` / `unipolar_plane` only
  if that config key is non-empty.
- Builds `m_lin_neg_unipolar` from `fields_neg_unipolar` / `unipolar_plane` only
  if that config key is non-empty.

If either unipolar config key is empty, the corresponding pointer stays null and
the trigger stub (which always returns 0) prevents the missing-pointer path from
being exercised.

### `compute_asym()` (anonymous helper, cxx)

Single helper that runs five passes over the ROI region and returns an
`AsymRecord` with all per-ROI features used by both the trigger and the
calibration dump:

1. In-ROI accumulators (`temp_sum`, `temp1_sum`, `temp2_sum`, `max_val`,
   `min_val`, threshold-gated by `adc_l1_threshold`).
2. Per-ROI gauss `gmax`, `gauss_fill`, `gauss_fwhm_frac`.
3. Wide-window decon energy fraction `roi_energy_frac` over
   `[roi_start − energy_pad, roi_end + energy_pad]`.
4. Wide-window raw asymmetry `raw_asym_wide` over
   `[roi_start − raw_asym_pad, roi_end + raw_asym_pad]`, with positive /
   negative samples gated by `raw_asym_eps`.
5. Core sub-window scan: identifies `|gauss| > core_g_thr` runs and reports
   the longest one's `(core_lo, core_hi, core_length, core_fill,
   core_fwhm_frac, core_raw_asym_wide)` for the dump record.

Definitions match `pdhd/nf_plot/find_long_decon_artifacts.py` (iter-7) so
Python ↔ C++ values can be compared bit-for-bit during validation.

### `decide_trigger()` (anonymous helper, cxx)

Independent of `compute_asym()`'s record so the dump and the live trigger
can disagree without coupling.  Walks every `|gauss| > core_g_thr`
sub-window in the ROI; for each sub-window computes its own
`(run_len, fill, fwhm, ef, aw)` and fires if all of:

- `gmax >= l1_gmax_min`
- `run_len >= l1_min_length`
- per-sub-window `ef >= l1_energy_frac_thr`
- ANY of the four arms:
  - `|aw| >= l1_asym_strong`, OR
  - `run_len >= l1_len_long_mod   AND |aw| >= l1_asym_mod`, OR
  - `run_len >= l1_len_long_loose AND |aw| >= l1_asym_loose`, OR
  - `run_len >= l1_len_fill_shape AND fill <= l1_fill_shape_fill_thr
                                  AND fwhm <= l1_fill_shape_fwhm_thr
                                  AND |aw| >= l1_asym_mod`.

Polarity = `sign(aw)` of the first firing sub-window.  Returns
`{−1, 0, +1}`.

### `l1_fit()` (cxx)

Called per ROI.  Sequence:

1. **Feature extraction** — `compute_asym()` once, populating the dump
   record (used in dump mode and never mutated afterwards).
2. **Trigger** — `decide_trigger()` returns `flag_l1 ∈ {0, ±1}`.
3. **Basis selection** — `flag_l1 = +1` → `{m_lin_bipolar, m_lin_pos_unipolar}`;
   `flag_l1 = -1` → `{m_lin_bipolar, m_lin_neg_unipolar}`.
4. **Pass-through guard** — if the matching unipolar basis is null
   (config keys empty), early-return without modifying `newtrace`.
5. **Segmented LASSO solve** — segments of `l1_seg_length` ticks; for each
   segment call `build_G` then `lasso_solve`.  Scatter beta into
   `l1_signal[t] += m_l1_basis0_scale * beta[j] + m_l1_basis1_scale * beta[N+j]`.
6. **Post-processing** — Gaussian smearing (if `filter` key is non-empty),
   zero samples below `l1_decon_limit`, remove peaks below `peak_threshold` /
   `mean_threshold`, write back into `newtrace`.

The propagation/hint-polarity path (uBooNE Layer 3) and the `flag_l1 = 2`
zero-out branch are deliberately not ported — see Design decisions below.

### `operator()` (cxx)

1. Retrieve `adctag` and `sigtag` traces.  Sizes must match.
2. Build per-channel tick sets from decon signal (positive samples) and raw ADC
   (samples above the noise threshold, padded ±`raw_pad` ticks).
3. Merge and pad into ROI pairs per channel (±`roi_pad`).
4. Call `l1_fit` on each in-scope ROI (channel filtered by `process_planes`
   and the optional `eligible_channels` whitelist).  Negative decon samples
   inside the ROI are zeroed in `newtrace` regardless of trigger result
   (matching `L1SPFilter` behaviour).
5. In `dump_mode`, write a per-frame NPZ with all ROI records plus the
   live `flag_l1` decision.
6. Emit `IFrame` with traces tagged `outtag`.

Layer 3 propagation and Layer 4 cross-channel cleaning are intentionally
omitted for PDHD/PDVD — see Design decisions.

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

The four arms map 1:1 onto the four iter-7 `cluster_pass()` rules in
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

### Calibration inputs still needed

**Unipolar field-response files**: the `fields_pos_unipolar` and
`fields_neg_unipolar` config keys accept any `IFieldResponse` type-name.
Until these are provided, the trigger fires correctly but `l1_fit()`
falls through to pass-through (the unipolar-nullptr early-exit at
the top of the LASSO block).

---

## Configuration parameters

### Input / output tags

| Key | Default | Meaning |
|-----|---------|---------|
| `adctag` | `"raw"` | Input tag for raw ADC traces (post-NF) |
| `sigtag` | `"gauss"` | Input tag for decon signal traces (post-OmnibusSigProc) |
| `outtag` | `"l1sp"` | Output tag for corrected signal traces |

### Field response

| Key | Default | Meaning |
|-----|---------|---------|
| `fields` | `"FieldResponse"` | IFieldResponse type-name for the bipolar basis |
| `bipolar_plane` | `1` | Plane index within `fields` (0=U, 1=V, 2=W) |
| `fields_pos_unipolar` | `""` | IFieldResponse for arriving-only (+unipolar) basis |
| `fields_neg_unipolar` | `""` | IFieldResponse for leaving-only (−unipolar) basis |
| `unipolar_plane` | `1` | Plane index within the unipolar field responses |

Leave both `fields_*_unipolar` keys empty to keep the component as a pass-through.

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
| `l1_scaling_factor` | `500` | Global scaling applied to the G matrix |
| `l1_lambda` | `5` | LASSO L1 regularization weight |
| `l1_epsilon` | `0.05` | Convergence tolerance |
| `l1_niteration` | `100000` | Maximum LASSO iterations |
| `l1_resp_scale` | `0.5` | Additional response scaling inside `build_G` |
| `unipolar_time_offset` | `3.0 µs` | Time offset of unipolar basis relative to bipolar |
| `overall_time_offset` | `0` | Global time offset applied to both bases |

### Output reconstruction

| Key | Default | Meaning |
|-----|---------|---------|
| `l1_decon_limit` | `100` | Zero samples below this value after solve |
| `l1_basis0_scale` | `1.15` | Weight for the bipolar (basis0) component |
| `l1_basis1_scale` | `0.5` | Weight for the unipolar (basis1) component |
| `peak_threshold` | `1000` | Drop an output ROI if its peak < this |
| `mean_threshold` | `500` | Drop an output ROI if its mean < this |
| `filter` | `[]` | Gaussian smearing kernel (empty = no smearing) |

### Electronics response

| Key | Default | Meaning |
|-----|---------|---------|
| `gain` | `14 mV/fC` | Preamp gain |
| `shaping` | `2.2 µs` | Shaping time |
| `postgain` | `1.2` | Post-amp gain |
| `ADC_mV` | `4096/2000` | ADC-to-mV ratio |
| `fine_time_offset` | `0` | Fine time offset |
| `coarse_time_offset` | `−8.0 µs` | Coarse time offset |
| `dft` | `"FftwDFT"` | IDFT component type-name |

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
| `flag_l1`             | int32\[\]   | Live `decide_trigger()` result `{−1, 0, +1}` under current config |

`flag_l1` is the trigger that drove the LASSO branch for this ROI.
Compare against legacy `flag` to see how the new gate diverges from the
uBooNE single-ratio decision.

---

## Wiring into a PDHD/PDVD graph

The wiring pattern mirrors uBooNE (`cfg/pgrapher/experiment/uboone/sp.jsonnet:65–163`).
The main differences are:

- `ChannelSelector` is optional; Strategy B does not require a static channel list.
  If the artifact is geometrically localized, a selector may still help reduce CPU.
- The `FrameMerger` feeding L1SP must supply both `adctag` and `sigtag` on the
  same frame.
- The output merger writes `outtag` back over `gauss` (and optionally `wiener`).

Skeleton (jsonnet, to be finalized after trigger calibration):

```jsonnet
local l1spfilterpd = g.pnode({
    type: "L1SPFilterPD",
    data: {
        dft: wc.tn(tools.dft),
        fields: wc.tn(tools.field),
        // fields_pos_unipolar: wc.tn(tools.field_pos_unipolar),  // uncomment when ready
        // fields_neg_unipolar: wc.tn(tools.field_neg_unipolar),
        adctag: "orig",    // post-NF raw ADC
        sigtag: "gauss",   // post-OmnibusSigProc decon
        outtag: "l1sp",
    }
}, nin=1, nout=1, uses=[tools.dft, tools.field]),
```

Gate the feature with a boolean in the params file (default `false`) so existing
production configs remain bit-identical with the feature off.

---

## Pending work

1. **Unipolar field-response inputs** — collect or derive the truncated-drift
   field responses for the PDHD/PDVD induction plane; supply them via
   `fields_pos_unipolar` / `fields_neg_unipolar` config keys.  Until these are
   provided the component is a graceful pass-through (trigger fires but the
   LASSO early-exits; see `l1_fit()` unipolar-nullptr check).

2. **LASSO body tuning** — once unipolar bases are populated, calibrate
   `unipolar_time_offset`, `l1_basis0_scale`, `l1_basis1_scale`, and the
   smearing `filter` kernel.  Trigger gate is independent of these and does
   not need to be re-run.

3. **PDVD jsonnet port** — the trigger logic itself is detector-symmetric;
   the only change needed in `cfg/pgrapher/experiment/protodunevd/sp.jsonnet`
   is the same `gain_scale` block applied here for PDHD.

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
