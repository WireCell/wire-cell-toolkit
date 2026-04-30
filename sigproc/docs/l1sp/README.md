# L1SP Filter

A compressed-sensing / LASSO-based sparse signal-processing stage for LArTPC
wire channels in **shorted-wire regions**, where induction- and collection-plane
wires are electrically connected and the readout observes a mixture of unipolar
and bipolar detector responses simultaneously.

**Canonical reference**: MicroBooNE Collaboration, _Ionization Electron Signal
Processing in Single Phase LArTPCs II: Data/Simulation Comparison and Performance
in MicroBooNE_, arXiv:1802.08709.

See also:

- [`index.org`](./index.org) — original brief intro and test-harness instructions.
- [`../examination/03-l1sp-filter.md`](../examination/03-l1sp-filter.md) — detailed
  code audit including six known bugs (BUG-L1-1..6) and six efficiency issues
  (EFF-L1-1..6). Any new integration should review that list before shipping.

---

## Why L1, not standard deconvolution

Standard 2D deconvolution inside `OmnibusSigProc` assumes each channel carries
**one** detector response. In a shorted region both collection (unipolar) and
induction (bipolar) responses are superimposed on the same wire. Deconvolving
with either response alone leaves large artifacts.

L1SP solves this by simultaneously fitting **both** response functions. For each
ROI segment of length N ticks it solves:

```
minimize  ‖ G · β − W ‖²  +  λ · ‖ β ‖₁
```

where:

- **W** ∈ ℝᴺ — observed ADC waveform in the segment.
- **β** ∈ ℝ²ᴺ — unknown signal coefficients: first half for collection (W-plane),
  second half for induction (V-plane).
- **G** ∈ ℝᴺˣ²ᴺ — forward model: each column is the detector response (field
  response convolved with cold-electronics shaping) for a unit impulse at that
  time, computed separately for W- and V-planes.
- **λ** (`l1_lambda`) — L1 regularization weight enforcing sparsity; physically
  appropriate because ionization from tracks occupies only a few ticks.

The solver is `WireCell::LassoModel` (`util/inc/WireCellUtil/LassoModel.h`).
The response interpolators `lin_W` / `lin_V` are built once in `init_resp()`
from `IFieldResponse`, averaged across the relevant plane wires, and convolved
with `Response::ColdElec` at the configured gain/shaping.

---

## Algorithm structure

**Source files:**

- `sigproc/inc/WireCellSigProc/L1SPFilter.h`
- `sigproc/src/L1SPFilter.cxx`

Registered as WCT component `L1SPFilter` (implements `IFrameFilter` +
`IConfigurable`; cxx:27-29).

### `init_resp()` (cxx:59-107)

Called at the start of each `operator()` invocation (idempotent guard inside).
Fetches the `IFieldResponse` named by `fields`, averages V- and W-plane
response vectors, convolves with `Response::ColdElec` at the configured
electronics parameters, and stores the result in `lin_V` / `lin_W`
linear-interpolation objects.

### `operator()` (cxx:220-494)

Main pipeline entry point. Receives one input `IFrame`, emits one output `IFrame`.

1. **Trace selection** (cxx:242-243): retrieve traces tagged `adctag` ("raw",
   post-NF ADC) and `sigtag` ("gauss", post-OmnibusSigProc decon signal).
2. **Decon ROI detection** (cxx:259-273): mark ticks where the gauss signal is
   positive charge.
3. **Raw ADC ROI detection** (cxx:280-306): mark ticks where raw ADC exceeds
   `raw_ROI_th_nsigma × noise_sigma` (percentile-estimated per channel), padded
   ± `raw_pad` ticks.
4. **Merge ROIs** (cxx:311-356): union decon and raw tick sets; merge gaps
   ≤ `roi_pad` ticks into contiguous ROI intervals.
5. **First L1 pass** (cxx:368-381): call `L1_fit` on each ROI. Returns flag:
   - **0** — no action (leave gauss signal unchanged).
   - **1** — shorted: L1 solve applied, output replaces gauss.
   - **2** — artifact: zero the output.
6. **Forward propagation** (cxx:383-399): scan ROIs in time order; if a
   flag-2 ROI is within 20 ticks of a flag-1 ROI, re-fit as shorted.
7. **Reverse propagation** (cxx:400-426): same scan in reverse.
   ⚠ Known bug BUG-L1-1: index logic is scrambled — see audit doc.
8. **Cross-channel cleaning** (cxx:432-476): zero any flag-2 ROI if either
   neighbour channel (ch ± 1) has a flag-1 ROI overlapping within 3 ticks.
9. **Emit output** (cxx:480-494): build new `IFrame` with traces tagged `outtag`
   ("l1sp") carrying the processed signal.

### `L1_fit()` (cxx:496-720)

Called per ROI. Signature: `int L1_fit(newtrace, adctrace, start_tick, end_tick, flag_shorted)`.

**Classification** (cxx:525-536):

- `flag = 1` (do L1 solve) if `sum/rescaled_abs_sum > adc_ratio_threshold`
  AND `abs_sum > adc_sum_threshold` — significant net-positive charge.
- `flag = 2` (zero) if rescaled sum is below limit, or large decon charge with
  small raw ADC.
- `flag = 0` otherwise (pass through).

**L1 solve** (cxx:538-591):

- Divide ROI into sections of ≤ `l1_seg_length` ticks.
- For each section: build `G (N×2N)` from `lin_W` (collection columns) and
  `lin_V` (induction columns), windowed to [−15 µs, +10 µs] around each tick.
- Run LASSO: `LassoModel(l1_lambda, l1_niteration, l1_epsilon)`.

**Post-processing** (cxx:593-680):

- Reconstruct signal: `β_col × l1_col_scale + β_ind × l1_ind_scale`,
  scaled by `l1_scaling_factor`.
- Convolve with `filter` smearing kernel.
- Zero bins below `l1_decon_limit / l1_scaling_factor`.
- Remove ROI if peak < `peak_threshold` or mean < `mean_threshold`.

---

## Configuration parameters

Defaults from `default_configuration()` (cxx:109-175). All are overridable in
the jsonnet config block.

### I/O & infrastructure

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `fields` | `"FieldResponse"` | WCT component name providing `IFieldResponse` |
| `filter` | `[]` | Smearing kernel (array of doubles). **Must be set by user.** |
| `adctag` | `"raw"` | Trace tag for NF-output raw ADC |
| `sigtag` | `"gauss"` | Trace tag for OmnibusSigProc gauss signal |
| `outtag` | `"l1sp"` | Trace tag on output waveforms |
| `dft` | `"FftwDFT"` | WCT DFT component name |

### ROI detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `raw_ROI_th_nsigma` | 4 | N-σ threshold for raw ADC ROI detection |
| `raw_ROI_th_adclimit` | 10 | Upper ADC limit (same threshold) |
| `roi_pad` | 3 | Gap tolerance (ticks) when merging ROI intervals |
| `raw_pad` | 15 | Padding added on each side of raw-ADC ROI hits |
| `overall_time_offset` | 0 | Global time offset (µs); ⚠ known unit bug — see audit |
| `collect_time_offset` | 3.0 | Collection-plane time offset relative to induction (µs) |

### ROI flag classification

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `adc_l1_threshold` | 6 | Per-tick ADC threshold entering sum |
| `adc_sum_threshold` | 160 | Minimum absolute ADC sum for flag-1 |
| `adc_sum_rescaling` | 90 | Rescaling divisor for ratio test |
| `adc_sum_rescaling_limit` | 50 | Lower bound of rescaled sum for flag-2 |
| `adc_ratio_threshold` | 0.2 | Minimum `sum/rescaled_abs_sum` for flag-1 |

### LASSO solver

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `l1_seg_length` | 120 | Max ticks per LASSO segment |
| `l1_scaling_factor` | 500 | Response matrix scale (numerical conditioning) |
| `l1_lambda` | 5 | L1 regularization strength |
| `l1_epsilon` | 0.05 | LASSO convergence tolerance |
| `l1_niteration` | 100000 | Max LASSO iterations |
| `l1_decon_limit` | 100 | Min output signal (electrons) after solve |

### Reconstruction & cleanup

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `l1_resp_scale` | 0.5 | Response scaling for matrix construction |
| `l1_col_scale` | 1.15 | Collection-component weight in output reconstruction |
| `l1_ind_scale` | 0.5 | Induction-component weight in output reconstruction |
| `peak_threshold` | 1000 | Drop ROI if peak < this (electrons) |
| `mean_threshold` | 500 | Drop ROI if mean < this (electrons) |

### Electronics response (used in `init_resp`)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `gain` | 14 mV/fC | Preamp gain |
| `shaping` | 2.2 µs | Shaping time |
| `postgain` | 1.2 | Post-amplifier gain |
| `ADC_mV` | 4096/2000 | ADC-to-mV conversion |
| `fine_time_offset` | 0 | Fine time offset (µs) |
| `coarse_time_offset` | −8.0 µs | Coarse time offset (µs) |

---

## How MicroBooNE wires L1SP in

**Config files:**

- `cfg/pgrapher/experiment/uboone/sp.jsonnet` (primary; lines 65-163)
- `cfg/layers/mids/uboone/api/sp.jsonnet` (layered-API mirror)

### Physical motivation

In the MicroBooNE detector a subset of U-plane (induction) and Y-plane
(collection) wires were electrically shorted, creating a ~740-channel region
(channels 3566–4305) where each readout wire carries a superposition of
collection and induction responses. Standard deconvolution leaves large
correlated residuals in this region; L1SP was developed to recover usable
signal there.

### Pipeline graph

```
rawsplit ──────────────────────────────────────────────────► rawsigmerge (port 1)
rawsplit ──► sigproc (OmnibusSigProc) ──► sigsplit ──► rawsigmerge (port 0)
                                                   └──► l1merge (port 1)   [wiener/gauss passthrough]
rawsigmerge ──► chsel ──► l1spfilter ──► l1merge (port 0)
                                                l1merge ──► [downstream]
```

`rawsigmerge` (FrameMerger, rule=`replace`) reassembles the "raw" tag from the
NF-output split with the "gauss" tag from OmnibusSigProc, giving L1SP both
inputs on a single frame.

`l1merge` (FrameMerger, rule=`replace`) writes the L1SP "l1sp" output back over
the "gauss" and "wiener" tags on the affected channels:

```jsonnet
mergemap: [
    ["raw",  "raw",    "raw"],
    ["l1sp", "gauss",  "gauss"],
    ["l1sp", "wiener", "wiener"],
]
```

Channels NOT in the `chsel` range keep their original OmnibusSigProc gauss/wiener
traces unchanged.

### Channel gating

```jsonnet
// cfg/pgrapher/experiment/uboone/sp.jsonnet:65-74
local chsel = g.pnode({
    type: "ChannelSelector",
    data: {
        channels: std.range(3566, 4305),  // shorted U/Y region
        tags: ["raw", "gauss"]
    }
}, nin=1, nout=1),
```

**There is no runtime predicate** — L1SP runs every event, but only on this
static channel list. Within the component, `L1_fit` applies the per-ROI
{0,1,2} flag logic that gates whether each individual ROI receives an L1 solve.

### MicroBooNE production overrides

Compared to `default_configuration()`, uBooNE tightens two ROI thresholds and
adds an explicit 21-tap Gaussian smearing kernel:

```jsonnet
raw_ROI_th_nsigma:    4.2,   // default: 4
raw_ROI_th_adclimit:  9,     // default: 10
filter: [0.000305453, 0.000978027, 0.00277049, 0.00694322, 0.0153945,
         0.0301973,   0.0524048,   0.0804588,  0.109289,   0.131334,
         0.139629,    0.131334,    0.109289,   0.0804588,  0.0524048,
         0.0301973,   0.0153945,   0.00694322, 0.00277049, 0.000978027,
         0.000305453],
```

All other parameters match the defaults.

---

## Applicability to pdhd / pdvd

### Current state

Both `cfg/pgrapher/experiment/pdhd/sp.jsonnet` and
`cfg/pgrapher/experiment/protodunevd/sp.jsonnet` instantiate **only
`OmnibusSigProc`** — no `FrameSplitter`, no `FrameMerger`, no `L1SPFilter`, no
channel-restricted branch. The same is true for `pdsp`.

`OmnibusSigProc` includes its own ROI refinement (`use_roi_refinement: true` in
pdhd configs, with tight/tighter ROI filters and multi-plane protection). This is
a 2D Wiener-based ROI refinement, **not** the L1-norm sparse joint-response fit
that `L1SPFilter` performs. The two are complementary, not duplicates.

### Questions to answer before deciding to integrate

1. **Does a shorted-wire (or equivalent mixed-response) region exist in pdhd or
   pdvd?** The MicroBooNE U/Y short is a hardware defect; it does not exist by
   design in ProtoDUNE. If no channel sees a physical superposition of two
   distinct response functions, L1SP solves a problem that does not exist.

2. **3-view geometry in pdvd**: ProtoDUNE-VD uses CRP anodes with two induction
   views (U, V) and one collection view (X). L1SP's response model is
   hard-coded to two components (W-plane collection, V-plane induction). It is
   not directly applicable to a 3-response mix without code modification.

3. **Response rescaling**: `l1_col_scale = 1.15` and `l1_ind_scale = 0.5` were
   derived from MicroBooNE field responses. If L1SP is used on pdhd/pdvd, these
   must be re-derived from the appropriate `IFieldResponse` objects. The same
   applies to `gain`, `shaping`, `coarse_time_offset`, and `collect_time_offset`.

### Prerequisites for integration (if motivated by a real detector effect)

- Identify the channel range(s) that require the joint-response treatment, and
  confirm the physical cause.
- Mirror the uBooNE `rawsplit / sigsplit / rawsigmerge / chsel / l1spfilter /
  l1merge` subgraph in the relevant cfg jsonnet (model on
  `cfg/pgrapher/experiment/uboone/sp.jsonnet:65-163`).
- Re-derive `l1_col_scale`, `l1_ind_scale`, `collect_time_offset`, and the
  smearing `filter` kernel from pdhd/pdvd detector responses.
- Gate the feature with a jsonnet boolean (default `false`) so existing
  production configs remain bit-identical with the feature off.
- Review (and ideally fix) **BUG-L1-1** (reverse-scan index mismatch) in
  `L1SPFilter.cxx:378-397` before shipping on new data.

---

## Test harness

From `index.org` — requires Magnify-format ROOT input:

```
wire-cell -V input=mag.root -V output=l1sp.root \
          -c l1sp/mag-l1sp-mag.jsonnet
```

where `mag-l1sp-mag.jsonnet` lives in the `wire-cell-cfg` repository (not this
toolkit). It sets up the `ChannelSelector → L1SPFilter → FrameMerger` sub-graph
driven by real MicroBooNE magnify data.
