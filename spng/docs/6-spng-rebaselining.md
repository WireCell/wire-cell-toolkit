# #6 — ROI rebaselining in SPNG: the problem and the `SPNGMaskRebaseliner` fix

Written 2026-10-18. Companion to [#3](3-spng-vs-osp-study.md)
§"Charge normalisation" (where the problem was found) and the study
[README](../studies/spcmp/README.md). Code is on `wire-cell-toolkit` branch
`modernize_osp`, **uncommitted** at the time of writing (file list at the end).

## 1. What "rebaselining" does in SPNG

After deconvolution and the gauss time filter, SPNG multiplies the waveform
by the Boolean ROI mask from the DNN (`SPNGThreshold` on the network
output, upsampled ×4 with `repeat_interleave`, so ROI edges lie on a 4-tick
grid). The product is then "rebaselined": inside every ROI a straight line
is subtracted so that the ROI's charge is measured relative to a local
baseline rather than relative to zero. This is meant to remove the slowly
varying residual the deconvolution leaves near signals — mainly on the
induction planes, where the low-frequency remnant of the bipolar response
survives the filters (measured ±5–10 % of pulse peak near ROI edges on the
smoke event; ≈ 0 on the collection plane). It is **not** a noise filter:
electronics noise is zero-mean and sums out over an ROI.

OSP does the same thing in `ROI_refinement::apply_roi`
(`sigproc/src/ROI_refinement.cxx:119-218`): a line between the two end
samples of the ROI — but of the *extended/padded* ROI, so the anchors sit
in the padding outside the tight signal region — and it does **not**
rebaseline the collection plane at all (l.179-217).

## 2. The bug

`SPNGRebaseliner` (`spng/src/Rebaseliner.cxx` → `rebaseline_zero()` in
`spng/src/Rebaseline.cxx`) only ever receives `gauss * mask`. It finds ROIs
as runs of nonzero samples (delimited by ≥ 2 exact zeros) and subtracts the
line through the **first and last samples inside the ROI**
(`apply_rebaseline_to_roi_like_original`, l.34-64). With SPNG's tight,
pure ROIs those two samples are still on the pulse — at ~8 % of peak on the
smoke event, ~13 % for a ±2σ Gaussian — so the subtracted "baseline" is a
positive trapezoid under the whole pulse and real charge ∝ (edge amplitude
× ROI width) is removed:

| W plane, smoke event (muon θxz=45 θyz=30) | charge / truth |
|---|---|
| gauss decon, no ROI, no rebaseline (in true pixels) | 1.01 |
| decon inside SPNG's ROI mask (mask loss only) | 0.996 |
| after `SPNGRebaseliner` | **0.94** |

The loss grows with ROI tightness (threshold 0.7 > 0.5), is larger on V/W
than U, and is −13…−18 % per channel for electron showers (wide, low,
fragmented signals). An emulation with anchors 4 ticks outside the ROI
recovered 1.00–1.01 (#3 §"Direct test of the mechanism").

The old CUDA path had a second problem: two `.item<float>()` host syncs plus
an `arange` and `index_put_` **per ROI** (l.57-58), so it took 350 ms per
APA on an RTX 4090 against 8 ms for the raw-pointer CPU loop.

## 3. The fix: `SPNGMaskRebaseliner`

A new two-input node — port 0 the unmasked waveform, port 1 the Boolean
mask — replaces the pair `SPNGReduce(mul)` + `SPNGRebaseliner` in
`applyroi_view`. Because it sees the unmasked waveform, the anchors can be
placed *outside* the ROI. The old node and the free function
`rebaseline_zero()` are untouched and remain selectable.

Files:
- `spng/inc/WireCellSpng/Rebaseline.h`, `spng/src/Rebaseline.cxx`:
  `MaskRebaselineOptions`, `rebaseline_masked()`, `rebaseline_masked_rows()`.
- `spng/inc/WireCellSpng/MaskRebaseliner.h`, `spng/src/MaskRebaseliner.cxx`:
  the node (`FaninBase<ITorchTensor>`, Hana config, factory name
  `SPNGMaskRebaseliner`).
- `spng/test/doctest_rebaseline_masked.cxx` (10 cases, 169 assertions with
  the legacy suites), `spng/test/test_rebaseline_masked_timing.cxx`.
- `spng/cfg/spng/subgraphs.jsonnet`, `adc-to-spng.jsonnet`,
  `detector.jsonnet`, `detconfs/pdhd.jsonnet`, `detconfs/pdhd_badAPA1.jsonnet`.

### 3.1 Algorithm

Inputs: waveform `W` (float) and mask `M` (bool; nonzero for other dtypes),
same shape, ROIs along `dim` (ticks).

1. **Dilate** the mask by `anchor_pad` samples on each side (max-pool with
   kernel `2·pad+1`, stride 1). ROIs closer than `2·pad` merge; that is
   intended, an anchor in such a gap would be on signal anyway.
2. **ROI geometry on the dilated mask**: for each dilated ROI `[s, e]` the
   anchors are its first and last samples, i.e. `pad` samples outside the
   tight ROI. With `anchor_width > 1` each anchor is the mean of `width`
   samples extending away from the ROI, clipped so the window never reaches
   into a neighbouring dilated ROI or past the tensor edge (an ROI touching
   the tensor edge falls back to its own edge sample there).
3. **Baseline** = the line through `(s, a)` and `(e, b)`; ROIs with
   `e − s ≤ min_roi_size − 1` are left untouched.
4. **Output** = `(W − baseline)` inside the output mask (the tight mask, or
   the dilated one if `keep_pad`), 0 elsewhere; optional `remove_negative`.

`anchor_pad = 0, anchor_width = 1` reproduces `rebaseline_zero(W·M)`
exactly (tested on the legacy fixtures and randomised inputs; max |diff|
1.2e-4 on a 2560×6000 workflow tensor).

### 3.2 Two implementations, one semantics

- **Vectorised torch path** (`rebaseline_masked_core`): mask shifts, running
  `cummax`/`cummin` scans for per-sample ROI start/end and neighbour
  bounds, a float64 prefix sum + `gather` for the anchor-window means, one
  `where`. No per-ROI loop, no host/device synchronisation. Used on CUDA
  (and on CPU with `force_vectorized`, a testing aid). Peak transient memory
  ≈ 1 GB for 2560×6000; the node processes row ranges separately, which
  bounds it per view.
- **Direct CPU path** (`rebaseline_masked_core_cpu`): one raw-pointer scan
  per row, rows in parallel with `at::parallel_for`. Used for float32 CPU
  tensors. Cross-checked against the vectorised path over 90 parameter
  combinations (`pad ∈ {0,1,3,4,7}`, `width ∈ {1,2,5}`, `keep_pad`,
  `min_roi_size`) at 1e-4.

Timing, 2560×6000, mean of 3 (`test_rebaseline_masked_timing`):

| | legacy `rebaseline_zero` | new, pad 0 | new, pad 4 | new, torch path |
|---|---|---|---|---|
| CPU (64 cores) | 8.3 ms | 11.4 ms | 11.8 ms | 263 ms |
| CUDA (RTX 4090) | 347 ms | 11.5 ms | 11.8 ms | 11.8 ms |

So the new node is ~1.4× the legacy cost on CPU (with padding and anchor
averaging) and ~30× faster than legacy on GPU. The ~12 ms is small against
the ~3 s of a full `adc-to-spng` job.

### 3.3 Node configuration

```
SPNGMaskRebaseliner data:
  dim: -1               # tick dimension
  anchor_pad: 4         # samples outside the tight ROI for the anchors (0 = legacy)
  anchor_width: 1       # samples averaged per anchor
  min_roi_size: 1
  keep_pad: false       # output the dilated ROI instead of the tight one
  remove_negative: false
  rows: [[b,e], ...]    # channel-row ranges to rebaseline; empty = all
  enabled: true         # false -> output = waveform*mask
  tag: "signal"
  datapath_format: ""
```
The output metadata derives from the **waveform** input only (the mask's
`tag:"roi"`/datapath must not leak into the signal tensor). The node raises
on shape mismatch, bad row ranges or a range beyond the channel count.

### 3.4 Jsonnet plumbing

- `subgraphs.jsonnet`
  - `rebaseline_default` = `{anchor_pad:4, anchor_width:1, min_roi_size:1, keep_pad:false, remove_negative:false, views:[0,1,2]}`.
  - `view_rows(sel_views, cat_views)` turns a list of views into channel-row
    ranges of the concatenated tensor using the new `tpc.view_nchans`
    (`[800,800,960]` for pdhd; set in `detconfs/pdhd.jsonnet` and
    `pdhd_badAPA1.jsonnet`, declared in `detector.jsonnet tpc()`).
  - `applyroi_view(view, extra_name, rebaseline=null, rows=[], enabled=null)`:
    `rebaseline == null` → the legacy `mul`+`SPNGRebaseliner` pair, exactly
    as before; otherwise one `SPNGMaskRebaseliner`. `applyroi_views`
    forwards. `dnnroi_inference_cat(..., rebaseline=$.rebaseline_default)`
    passes `rows = view_rows(rebaseline.views)` so the per-view selection
    works on its single all-planes node. The per-view graphs
    (`dnnroi_inference`, `_simple`, `_pdhd_apa1_regres`) and `roi.jsonnet`
    still build the legacy pair — switching them is a one-argument change.
- `adc-to-spng.jsonnet` TLAs: `rebaseline_pad=4`, `rebaseline_width=1`,
  `rebaseline_keep_pad=false`, `rebaseline_views=[0,1,2]`,
  `rebaseline_legacy=false` (the last three via `--tla-code`).
  `dump=rebaseline` finds the new node too.

```
wire-cell -A input=adc.npz -A model_file=m.ts -A threshold=0.7 -A device=gpu \
          -A rebaseline_pad=4 --tla-code "rebaseline_views=[0,1]" spng/adc-to-spng.jsonnet
```

## 4. Validation

Smoke event (muon θxz=45 θyz=30, SPNG threshold 0.7, truth threshold 400 e⁻):

| variant | U q_ratio / purity | V | W |
|---|---|---|---|
| legacy (old node) | 0.983 / 0.992 | 0.928 / 0.993 | 0.938 / 0.998 |
| **default** (pad 4, all views) | **1.011 / 0.992** | **0.991 / 0.993** | **1.013 / 0.998** |
| noW (views [0,1]) | 1.011 / 0.992 | 0.991 / 0.993 | 1.002 / 0.998 |
| keeppad | 1.011 / 0.859 | 0.995 / 0.758 | 1.018 / 0.740 |
| OSP+DNNROI (reference) | 1.022 / 0.867 | 1.015 / 0.781 | 1.016 / 0.728 |

The charge deficit is gone with the ROI (purity) untouched; `keep_pad`
shows what including the tails costs in purity, which is what DNNROI's
wider ROIs also pay. GPU and CPU runs of the full job agree to 4e-5
relative. The remaining ~1 % excess on W (1.013) is the same +1–2 % seen in
the un-rebaselined decon and in DNNROI, i.e. not from this node.

Full sample (294 events, `studies/spcmp` config `spng_rebaseline` +
`variants`; details in #3 §"Rebaseliner fix: full-sample results"):
muon Q/truth per plane 0.96/0.93/0.98 → **0.995/0.990/1.011**, electron
0.92/0.90/0.93 → **0.97/0.96/0.98**, per-channel bias −5…−7 % → ≤ 1 %
(muons) and −13…−18 % → −5 % (electrons), all ROI metrics unchanged, plane
spread 0.067 → 0.034. `noW` is marginally better on plane consistency;
`keeppad` raises the 3-view trio fraction 0.94 → 0.98 at DNNROI-like
purity.

Job-level cost, per-node wire-cell timers in real `adc-to-spng` jobs
(`rebaseline_legacy=true` selects the old node; one pdhd APA). The old
GPU node's cost is per ROI (two host syncs each), so it only shows on
ROI-dense events:

| event | ROIs | old `mul`+`Rebaseliner` CPU / GPU | new `MaskRebaseliner` CPU / GPU | whole job CPU / GPU (old → new) |
|---|---|---|---|---|
| single muon (1 frame) | few hundred | 24 ms / 30 ms | 23 ms / 62 ms¹ | 5.0 → 5.0 s / 1.4 → 1.7 s |
| cosmics (10 frames, ~3.7k ROIs and 1.8k active channels per frame) | ~37k | 226 ms / **2.16 s** | 181 ms / **93 ms** | 37.6 → 36.8 s / 13.0 → 10.8 s |

¹ one call per job, so first-call CUDA kernel setup is included.

On the cosmics event the new node is 23× faster on GPU (216 → 9 ms per
frame) and 20 % faster on CPU, worth 17 % of the whole GPU job. Files:
`/nfs/data/1/calcuttj/data/spcmp-angled/cosmics-timing/spng-{legacy,fixed}-{gpu,cpu}.log`.

## 5. Open points

- `anchor_pad = 4` = one DNN output pixel. It was chosen from the emulation
  (anchors at 2 ticks already recover 0.98, 4 ticks 1.00–1.01); a scan of
  pad and width on the full sample would pin it down.
- Whether to rebaseline W at all (OSP does not): the `noW` variant tests
  it; on this simulation W has no residual to remove, real data may differ.
- The legacy pair is still what the per-view graphs use; once the new node
  is accepted they should switch, and `SPNGRebaseliner` could be retired
  or get the vectorised CUDA path.
- Not touched: the ~1-tick SPNG/OSP timing offset (#5).

## 6. Toolkit files changed (to commit)

```
spng/inc/WireCellSpng/Rebaseline.h            + MaskRebaselineOptions, rebaseline_masked(), rebaseline_masked_rows()
spng/src/Rebaseline.cxx                       + implementations (vectorised + CPU direct)
spng/inc/WireCellSpng/MaskRebaseliner.h       new node
spng/src/MaskRebaseliner.cxx                  new node
spng/test/doctest_rebaseline_masked.cxx       new
spng/test/test_rebaseline_masked_timing.cxx   new
spng/cfg/spng/subgraphs.jsonnet               rebaseline_default, view_rows, applyroi_view(s), dnnroi_inference_cat (+ threshold param from earlier)
spng/cfg/spng/adc-to-spng.jsonnet             rebaseline_* and threshold TLAs, dump finder
spng/cfg/spng/detector.jsonnet                tpc(view_nchans)
spng/cfg/spng/detconfs/pdhd.jsonnet           view_nchans=[800,800,960]
spng/cfg/spng/detconfs/pdhd_badAPA1.jsonnet   view_nchans=[800,800,960]
```
