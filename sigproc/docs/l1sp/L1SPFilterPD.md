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
| Polarity | always flag=1 (net-positive) | flag=+1 / −1 / 0 (stub, see below) |
| Layer 4 cross-channel cleaning | ✓ (ch±1, ±3 ticks) | removed |
| Propagation polarity tracking | `bool flag_shorted` | `int active_polarity` (±1 / 0) |
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

### `l1_fit()` (cxx:267–438)

Called per ROI.  Sequence:

1. **ADC accumulation** — compute `temp_sum` (signed), `temp1_sum` (absolute),
   `temp2_sum` (absolute decon), `max_val`, `min_val` over samples above
   `adc_l1_threshold`.

2. **Trigger (STUB)** — currently always returns `flag_l1 = 0` (pass-through).
   See [stub](#trigger-stub-and-strategy-b-implementation) below.

3. **Polarity override** — if `hint_polarity ≠ 0` (propagation pass), that
   polarity is used instead of the stub result.

4. **Basis selection** — `flag_l1 = 1` → `{m_lin_bipolar, m_lin_pos_unipolar}`;
   `flag_l1 = -1` → `{m_lin_bipolar, m_lin_neg_unipolar}`.

5. **Segmented LASSO solve** — segments of `l1_seg_length` ticks; for each
   segment call `build_G` then `lasso_solve`.  Scatter beta into
   `l1_signal[t] += m_l1_basis0_scale * beta[j] + m_l1_basis1_scale * beta[N+j]`.

6. **Post-processing** — Gaussian smearing (if `filter` key is non-empty),
   zero samples below `l1_decon_limit`, remove peaks below `peak_threshold` /
   `mean_threshold`, write back into `newtrace`.

7. **Zero-out path** — `flag_l1 = 2`: zero the ROI in `newtrace`.

### `operator()` (cxx:440–631)

1. Retrieve `adctag` and `sigtag` traces.  Sizes must match.
2. Pre-compute `ntot_ticks` (max trace length) — avoids repeated per-trace passes.
3. Build per-channel tick sets from decon signal (positive samples) and raw ADC
   (samples above noise threshold, padded ±`raw_pad` ticks).
4. Merge and pad into ROI pairs per channel (±`roi_pad`).
5. **First pass**: call `l1_fit` on each ROI; record `flag_rois`.  Also zero any
   negative decon values inside the ROI (same behaviour as `L1SPFilter`).
6. **Forward propagation** (Layer 3): iterate ROIs in time order; maintain
   `active_polarity` (±1 / 0).  A flag-2 ROI within 20 ticks of a flag-±1 ROI
   is re-fitted with that polarity hint.
7. **Reverse propagation**: same, iterating in reverse time order.
8. **Layer 4 omitted** — cross-channel cleaning (ch±1, ±3 ticks) is not applicable
   to PDHD/PDVD geometry.
9. Emit `IFrame` with traces tagged `outtag`.

---

## Trigger stub and Strategy B implementation

### Current state

The trigger in `l1_fit` (cxx:298–341) always returns `flag_l1 = 0` (pass-through).
The ADC-sum variables are computed but suppressed with `(void)` casts, and the
comment enumerates the planned cuts.

### When implementing Strategy B

Replace the stub with:

```cpp
// Guard: only evaluate if there is enough signal.
if (temp1_sum > m_adc_sum_threshold && temp1_sum > 0) {
    double ratio = m_adc_sum_rescaling / nbin_fit;  // ~rescaling factor
    double asym  = temp_sum / (temp1_sum * ratio);

    if (asym > m_adc_ratio_threshold)        flag_l1 =  1;  // +unipolar
    else if (-asym > m_adc_ratio_threshold)  flag_l1 = -1;  // −unipolar
    else if (temp1_sum * ratio < m_adc_sum_rescaling_limit) flag_l1 = 2; // artifact
}
```

The sign of `asym` selects the polarity; `m_adc_ratio_threshold` is the asymmetry
cut to tune.  All threshold parameters are already in `default_configuration()`.

### Calibration inputs needed

Before filling in the trigger:

1. **Event-scan + ROI tagging**: scan PDHD/PDVD data, tag ROIs where the raw ADC
   is clearly unipolar (by eye or with a prototype classifier), and histogram
   `asym` for the tagged vs untagged populations.  This gives the operating point
   for `adc_ratio_threshold`.

2. **Unipolar field-response files**: the `fields_pos_unipolar` and
   `fields_neg_unipolar` config keys accept any `IFieldResponse` type-name.  Until
   these are provided (and the trigger is enabled), `L1SPFilterPD` acts as a
   pass-through on every frame.

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

### Trigger thresholds

| Key | Default | Meaning |
|-----|---------|---------|
| `adc_l1_threshold` | `6` | Min |ADC| to include a sample in the sum |
| `adc_sum_threshold` | `160` | Min Σ|ADC| required to evaluate the trigger |
| `adc_sum_rescaling` | `90` | Rescaling denominator for the asymmetry ratio |
| `adc_sum_rescaling_limit` | `50` | Rescaling guard for artifact detection |
| `adc_ratio_threshold` | `0.2` | Asymmetry ratio cut for flag=±1 |

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

1. **Trigger cuts** — run the event-scan + ROI tagging analysis on PDHD/PDVD data;
   fill in the Strategy B condition in `l1_fit()` (cxx:298–341, see stub comment).

2. **Unipolar field-response inputs** — collect or derive the truncated-drift
   field responses for the PDHD/PDVD induction plane; supply them via
   `fields_pos_unipolar` / `fields_neg_unipolar` config keys.

3. **Threshold calibration** — retune `adc_ratio_threshold`, `unipolar_time_offset`,
   `l1_basis0_scale`, `l1_basis1_scale`, and the smearing `filter` kernel from
   PDHD/PDVD field responses.

4. **Jsonnet wiring** — add `L1SPFilterPD` to the PDHD/PDVD graph in
   `cfg/pgrapher/experiment/pdhd/sp.jsonnet` (and the pdvd equivalent) behind a
   boolean gate (default off).

5. **Review BUG-L1-1** from `sigproc/docs/examination/03-l1sp-filter.md` before
   relying on the propagation layers in production — the reverse-scan index
   mismatch in `L1SPFilter.cxx:378–397` was reproduced in the propagation logic
   here; it should be verified or fixed before shipping.
