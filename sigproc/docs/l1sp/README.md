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
  Because the writeback line multiplies β by `l1_scaling_factor`, the
  per-coefficient sparsity threshold/shrinkage in electron units is
  `λ_in_e = l1_lambda × l1_scaling_factor` (default 10 × 500 = 5000 e).

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
- For each section: build `G (M×2N)` from `lin_W` (collection columns) and
  `lin_V` (induction columns), windowed to [−15 µs, +10 µs] around each tick.
  M ≥ N: the W vector is padded by 30 ticks before / 20 ticks after the
  segment's β span (clipped to the trace bounds) so boundary β
  coefficients see the full kernel response and cannot grow to fit
  imaginary out-of-window signal. The padded raw ADC is fit context
  only — β positions and writeback range stay strictly inside the ROI.
- Run LASSO: `LassoModel(l1_lambda, l1_niteration, l1_epsilon)`.

**Post-processing** (cxx:850-905):

- Combine basis coefficients and rescale to electrons in one step:
  `l1_signal[t] = (β₀ × l1_basis0_scale + β₁ × l1_basis1_scale) × l1_scaling_factor`.
- Convolve with `filter` smearing kernel (sum-normalised to 1, so the
  integral is preserved; the result stays in electron units).
- Zero bins below `l1_decon_limit` (in electrons).
- Remove ROI if peak < `peak_threshold` or mean < `mean_threshold`.

---

## How L1SP is triggered: the flag logic in detail

The decision "should we run an L1 solve on this waveform?" is taken in **four
layers**, applied in sequence. None of them is a generic two-sided asymmetry
test — MicroBooNE's trigger is **one-sided, biased toward positive net
charge**, because the artifact L1SP was designed to fix is collection-plane
charge leaking onto a shorted induction wire.

### Layer 1 — Static channel selection (graph level)

Before any ROI logic runs, the WCT graph already restricts L1SP's input to a
fixed channel range via `ChannelSelector(channels = std.range(3566, 4305))`.
Channels outside that range never reach `L1SPFilter::operator()`. There is no
runtime predicate at this layer — the channel list is part of the cfg.

### Layer 2 — Per-ROI flag classification in `L1_fit()`

For each ROI on each selected channel, `L1_fit` (cxx:530-564) computes three
ADC sums over samples whose `|ADC| > adc_l1_threshold` (default 6):

```
temp_sum   = Σ ADC_i        (signed sum)
temp1_sum  = Σ |ADC_i|      (unsigned sum)
temp2_sum  = Σ |decon_i|    (unsigned sum of the gauss decon signal)
```

It then assigns one of three flags:

```cpp
// flag = 1  (run L1 solve)
if (temp1_sum > adc_sum_threshold &&
    temp_sum / (temp1_sum * adc_sum_rescaling / nbin_fit) > adc_ratio_threshold)

// flag = 2  (zero output -- artifact)
else if (temp1_sum * adc_sum_rescaling / nbin_fit < adc_sum_rescaling_limit)
else if (temp2_sum > 30*nbin_fit && temp1_sum < 2*nbin_fit && (max_W - min_W) < 22)

// flag = 0  (pass through, leave gauss decon unchanged)
else
```

Reading the flag-1 inequality: it is essentially `temp_sum / temp1_sum >
adc_ratio_threshold * adc_sum_rescaling / nbin_fit`. With defaults
(0.2 × 90 / nbin_fit ≈ 18/nbin_fit) it asks **"does the signed ADC sum
exceed a fraction of the unsigned sum, scaled by ROI length?"**

This **is** an asymmetry-style test, but **only one-sided**: `temp_sum` is
signed and the comparison is `> threshold`. A purely *negative* unipolar ROI
gives `temp_sum ≈ −temp1_sum` and fails the test; the ratio is large but
negative. So flag-1 fires on **net-positive** ROIs only. A balanced bipolar
ROI gives `temp_sum ≈ 0` and also fails (correctly — there is nothing to
fix). A fully positive collection-leaked ROI gives `temp_sum/temp1_sum ≈ 1`
and passes.

The flag-2 conditions are pure-noise cleanup: either the ROI has too little
unsigned ADC activity, or the decon claims a large signal where the raw ADC
is essentially flat (a deconvolution artifact).

### Layer 3 — Same-channel time propagation (forward + reverse)

`operator()` then walks the channel's ROI list **twice** — once in time
order, once in reverse (cxx:383-402 and cxx:404-426). Each pass carries a
`flag_shorted` boolean state:

- A flag-1 ROI **sets** `flag_shorted = true`.
- A flag-2 ROI **within 20 ticks** of the previous ROI, *while*
  `flag_shorted` is true, gets **rescued**: `L1_fit` is called again on it
  and its flag is rewritten to 0 (handled). The rationale: a real shorted
  region produces a run of ROIs along the wire; if one of them looked like
  an artifact (flag-2) but its neighbour was clearly shorted (flag-1), it
  was probably also shorted but failed the flag-1 cut by chance.
- A flag-0 ROI **clears** `flag_shorted` (the run has ended).
- A gap of more than 20 ticks between adjacent ROIs **clears**
  `flag_shorted` (different runs).

The reverse pass repeats the same logic scanning ROIs from latest to
earliest, so a flag-2 ROI can be rescued by a flag-1 neighbour on either
side. (Note: BUG-L1-1 in the audit doc reports the reverse pass has an
index-mismatch bug; the *intent* is what is described here.)

So nearby ROIs **are** considered — temporally, on the same channel,
within ±20 ticks.

### Layer 4 — Cross-channel cleaning

After all per-channel propagation finishes, a final loop (cxx:432-476)
looks at neighbouring channels (`ch ± 1`):

- Take each ROI on channel `ch` whose flag is still 2.
- If channel `ch + 1` *or* channel `ch − 1` has a flag-1 ROI whose time
  range overlaps within ±3 ticks → **zero this ROI's output**.
- Otherwise leave it alone.

This is the spatial counterpart to Layer 3: a flag-2 ROI sandwiched between
truly-shorted neighbours is treated as confirmed-artifact and silenced. So
nearby ROIs are also considered **across channels**, within ±1 channel and
±3 ticks.

### Summary of the trigger

| Question | Answer |
|----------|--------|
| Is the trigger asymmetry-based? | Yes, but **one-sided** (net-positive only). |
| Does it look at neighbouring ROIs in time? | Yes — ±20 ticks, same channel, both directions. |
| Does it look at neighbouring channels? | Yes — ±1 channel, ±3 tick overlap, for artifact cleanup. |
| Is it event-by-event or static? | Static channel range; per-ROI flag within that range. |

### Implication for PDHD/PDVD adaptation

`L1SPFilterPD` (see `L1SPFilterPD.md`) is the PD-specific implementation. It
keeps Layers 1 and 2 with modifications; Layers 3 and 4 are intentionally
dropped because their sole purpose was handling the uBooNE shorted-wire pathology:

- **Layer 1 (channel selection):** kept — configurable `process_planes` (default
  U+V) plus optional `eligible_channels` whitelist for a known positive-half list.
- **Layer 2 (per-ROI flag):** modified to two-sided `{0, +1, -1}` — positive
  asymmetry → `+1` (arriving electrons); negative asymmetry → `-1` (leaving
  electrons). No `flag=2` zero-out: PD has no shorted-wire phantom-signal pathology.
- **Layer 3 (time propagation):** dropped — it only rescued `flag=2` ROIs from
  neighbouring confident hits. With `flag=2` gone, Layer 3 is a no-op and has
  been removed.
- **Layer 4 (cross-channel cleaning):** dropped — it suppressed phantom signals on
  `ch±1` neighbours of confirmed shorted-region hits. PD channels are not shorted.

---

## Configuration parameters

Defaults from `default_configuration()` (cxx:109-175). All are overridable in
the jsonnet config block.

### I/O & infrastructure

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `fields` | `"FieldResponse"` | WCT component name providing `IFieldResponse` |
| `filter` | `[]` | Smearing kernel (array of doubles). PDHD/PDVD: leave empty to auto-derive from `Gaus_wide`; uBooNE keeps explicit 21-tap array. |
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
| `overall_time_offset` | 0 | Additive override (µs) on top of kernel-file `frame_origin_us`. The global LASSO frame origin is loaded from `meta.frame_origin_us` (= reference plane's bipolar zero crossing); this knob is for tuning only and should normally stay 0. |
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
| `l1_scaling_factor` | 500 | Numerical conditioning on G (cancels in linear algebra) |
| `l1_lambda` | 10 | L1 regularization strength; per-coefficient threshold in electrons = `l1_lambda × l1_scaling_factor` (default 5000 e) |
| `l1_epsilon` | 0.05 | LASSO convergence tolerance |
| `l1_niteration` | 100000 | Max LASSO iterations |
| `l1_decon_limit` | 100 | Per-tick floor (electrons) applied after smearing |

### Reconstruction & cleanup

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `l1_resp_scale` | 1.0 | Kernel amplitude scale; must be 1.0 for ADC/electron kernels |
| `l1_basis0_scale` | 1.0 | Bipolar (basis0) component weight (β₀ already in electrons) |
| `l1_basis1_scale` | 1.0 | Unipolar (basis1) component weight (β₁ already in electrons) |
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

The explicit 21-tap array is numerically equivalent to the IFFT-derived kernel
(σ = 0.111408 MHz, 500 ns tick) to within max |Δ| ≈ 5×10⁻⁶; see
`pdhd/nf_plot/plot_l1sp_smearing_kernel.py` (in `wcp-porting-validation`) for validation.

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

## Adapting L1SP for pdhd / pdvd unipolar-induction artifacts

### The physical problem

ProtoDUNE-HD and ProtoDUNE-VD induction planes see two distinct failure modes
that standard 2D deconvolution cannot handle, because it assumes the detector
response is **bipolar with zero net integral** (∫response dt ≈ 0).

**Case 1 — Anode-induction (electrons only leaving).** When ionization
originates at or very near the induction plane, the electrons only drift
*away* from the wire — they never approach it. Only the **negative** lobe of
the induced-current waveform is observed. The raw ADC on the induction wire
carries a purely negative unipolar signal.

**Case 2 — Imperfect-geometry collection on induction (electrons only
arriving).** In regions of the detector where mechanical geometry deviates
from design, some electrons that should pass through the induction plane
instead are collected on it. Only the **positive** (approaching) lobe of the
response is induced before the electrons stop. The raw ADC carries a purely
positive unipolar signal.

**Why standard deconvolution fails.** The inverse filter is the Fourier-space
reciprocal of the bipolar response. Applied to a unipolar input, the filter
tries to "close" the missing lobe, producing a long monotone tail that
stretches hundreds of ticks past the real hit location. The artifact is
visible in two reference examples from real data:

- V-plane channel 4340 (`pdvd/pics/Picture1.png`): strong positive/negative
  asymmetry in the raw waveform maps to a large oscillating artifact in the
  deconvolved signal.
- U-plane channel 48 (`pdvd/pics/Picture2.png`): a brief unipolar raw pulse
  deconvolves to a tall positive peak followed by a plateau extending
  hundreds of ticks — the 2D image shows it as long vertical colour streaks
  stretching well beyond the true track region.

### What L1SP contributes

L1SP's core machinery — a LASSO solver over a configurable G-matrix response
basis, per-ROI gating, and the existing WireCell graph plumbing (ChannelSelector,
FrameMerger) — is exactly the framework needed. Only two things require
adaptation:

1. **The response basis** in the G matrix. Instead of {bipolar-induction,
   collection-unipolar} as in MicroBooNE, we need one or more
   *truncated-unipolar* basis functions derived from the PDHD/PDVD field
   response.
2. **The ROI gating logic** inside `L1_fit()`. The current flag-1 test
   (`temp_sum / rescaled_abs_sum > adc_ratio_threshold`) detects net-positive
   charge as the proxy for "shorted". The PDHD/PDVD problem requires a
   **polarity-asymmetry** test that fires on either strongly positive *or*
   strongly negative ROIs.

### Strategy A — Static channel list (uBooNE-style port)

Identify the geometric channel ranges where the collection-on-induction case
is concentrated (e.g. channels near the CRP/APA edge where the field
geometry distorts). Construct a "arriving-only" response (positive lobe of
the field response, zeroed past the zero-crossing) and wire L1SP exactly as
MicroBooNE does:

```
G ∈ ℝᴺˣ²ᴺ  with  col[0..N] = standard bipolar induction response
                   col[N..2N] = arriving-only (positive-half) response
```

**Pros:** smallest code change; the entire uBooNE graph fragment
(`rawsplit / sigsplit / rawsigmerge / chsel / l1spfilter / l1merge`) ports
without modification.

**Cons:** the **anode-induction case is event-by-event** (any track that
starts near the anode anywhere in the active volume produces it). A static
channel list cannot capture it.

### Strategy B — Per-ROI polarity-asymmetry detection *(implemented in `L1SPFilterPD`)*

Run `OmnibusSigProc` on all channels as today. Then apply a downstream L1SP
pass only on ROIs where the **raw ADC** polarity asymmetry exceeds a
threshold.

**Implemented** in `L1SPFilterPD::l1_fit()`:

```cpp
double ratio = (temp1_sum > 0)
             ? temp_sum / (temp1_sum * adc_sum_rescaling / nbin_fit)
             : 0.0;
int flag_l1 = 0;
if (temp1_sum > adc_sum_threshold) {
    if      (ratio >  adc_ratio_threshold) flag_l1 = +1;  // positive unipolar
    else if (ratio < -adc_ratio_threshold) flag_l1 = -1;  // negative unipolar
}
```

When flagged, the G matrix is assembled with the polarity-appropriate basis:

- `flag_l1 = +1` → col[N..2N] = arriving-only (`fields_pos_unipolar`) response.
- `flag_l1 = -1` → col[N..2N] = leaving-only (`fields_neg_unipolar`) response.

col[0..N] retains the standard bipolar induction response in both cases.

**Threshold values** (`adc_sum_threshold=160`, `adc_sum_rescaling=90`,
`adc_ratio_threshold=0.2`) are inherited from uBooNE and need re-tuning from
PDHD/PDVD dump-mode hand-scan data (see `L1SPFilterPD.md`, Pending Work).

**Pros:** no static channel list; handles event-by-event unipolar occurrence;
detection test is cheap (ADC sum over ROI window).

**Cons:** threshold calibration required; unipolar field-response files for
PDHD/PDVD needed before the fit actually runs (component degrades gracefully
to pass-through without them).

### Strategy C — 3-basis LASSO (most flexible, highest cost)

Extend G to three response blocks:

```
G ∈ ℝᴺˣ³ᴺ  with  col[0..N]   = standard bipolar induction
                   col[N..2N]  = positive-unipolar (arriving-only)
                   col[2N..3N] = negative-unipolar (leaving-only)
```

The LASSO solver discovers the mixture per tick without requiring an upstream
detector. Normal hits get solved by the bipolar block; unipolar hits get
solved by whichever unipolar block matches.

**Pros:** handles both unipolar cases plus normal bipolar uniformly without
any classifier.

**Cons:** the 30–50% larger matrix increases solve time; the positive-unipolar
and bipolar responses share the positive lobe, making the G columns
correlated and the LASSO convergence slower or less stable. Needs careful
λ retuning.

### Building the unipolar response basis

Two options, in increasing accuracy:

1. **Approximate (fast start):** take the existing `IFieldResponse` for the
   induction plane, load the average response into a temporary vector, and
   zero all samples after the zero-crossing. This gives the arriving-only
   (positive) half. Negate and zero the other half for the leaving-only
   (negative) version. No new Garfield/COMSOL run needed.

2. **Physical (higher accuracy):** re-run the field-response calculation
   starting the drift path at the induction plane (no approaching phase) or
   terminating it on the induction plane (no leaving phase). This correctly
   captures the field shaping near the boundary; the approximate approach
   may mis-model the onset ramp.

Recommendation: start with (1) to validate the L1SP path end-to-end; switch
to (2) if the reconstruction residuals are larger than acceptable.

### PDVD 3-view geometry

PDVD has **two induction views** (U, V) and one collection view (X). Both U
and V can suffer both unipolar cases. L1SP's response storage currently
hard-codes two `linterp<double>` objects (`lin_V`, `lin_W` in `L1SPFilter.h`)
and averages `fravg.planes[1]` and `fravg.planes[2]` in `init_resp()` (cxx:80-81).

Two practical options for supporting both induction views:

- **Two L1SP instances:** instantiate `L1SPFilter` twice in the WCT graph,
  once for U channels (with U-plane field response) and once for V channels
  (with V-plane field response). Each uses the same code; only the
  `field_response`, `collect_time_offset`, and channel-selection inputs
  differ. This is the lowest-risk option.
- **Generalize `init_resp()`:** add a `plane_index` config key that selects
  which plane from `IFieldResponse` to average into `lin_V`. This removes the
  duplication but requires a small code change.

The unipolar-basis construction above applies identically to both views.

### Recommended first step: diagnostic-only asymmetry scan

Before touching any signal-processing path, implement a standalone diagnostic
that reads existing PDHD/PDVD processed frames and, for each induction-plane
ROI, computes and logs the polarity asymmetry. This is safe (read-only), fast
to write (~50 lines), and answers the key scoping questions:

- What fraction of induction ROIs are significantly asymmetric?
- Are asymmetric ROIs concentrated in predictable geometric regions (→
  supports Strategy A) or distributed across all channels and events (→
  motivates Strategy B)?
- Is the negative-unipolar (anode-induction) case, the positive-unipolar
  (collection-on-induction) case, or both common enough to matter?

The diagnostic drives the choice between Strategies A, B, and C without any
risk of altering production outputs.

### Open questions before implementation

1. **Unipolar-response basis files** — do the PDHD/PDVD field-response sets
   include a truncated-drift variant, or should we start with the
   zero-the-second-lobe approximation described above?
2. **Geometric localization** — is the collection-on-induction case
   concentrated in specific known channel ranges (edge channels, dead-wire
   neighbours), or is it scattered? If localized, Strategy A may be
   sufficient for that case.
3. **Validation data** — can existing data or simulation provide
   hand-labelled unipolar ROIs to tune and validate the asymmetry threshold?
4. **Threshold calibration** — the polarity-asymmetry cut and the modified
   flag-1 condition both need dedicated tuning on pdhd/pdvd field responses,
   since the MicroBooNE values (net-positive proxy for shorted wires) are
   physically different from the new detection target.

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
