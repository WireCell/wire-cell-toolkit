# Stage 2 — WCT-native PDHD light reconstruction (DUNE-method port)

Status: complete; chain runs end-to-end on the example data.  Quantitative
comparison against the LArSoft reference products is in `validation.md`
(stage 3).

The chain ports the DUNE LArSoft optical reconstruction (the modules the PDHD
master job stages but does not yet run, see
`pdhd/docs/photon-detector-chain.md` §3) into three `flash/` components:

```
PDHDOpWaveformSource ("raw" snippets, OpChannel basis)
  -> Flash::OpDecon       [duneopdet Deconvolution, Wiener + SPE template]
  -> Flash::OpHitFinder   [duneopdet OpHitFinderDeco / larana AlgoSlidingWindow]
  -> Flash::OpFlashFinder [larana OpFlashAlg, protodune_opflash]
  -> opflash tensor set (design.md §3.4) -> opflash_pdhd-wct.tar.gz
```

## fcl ↔ config correspondence

### `Flash::OpDecon` (`IFrameFilter`, "raw" → "decon")

Port of `duneopdet/OpticalDetector/Deconvolution/Deconvolution_module.cc`
(duneopdet v10_20_09d00) with the `protodunehd_deconvolution` data values:

| fcl | OpDecon config | default |
|---|---|---|
| `Samples` | `samples` | 1024 |
| `PreTrigger` | `pre_trigger` | 50 |
| `PedestalBuffer` | `pedestal_buffer` | 30 |
| `LineNoiseRMS` | `line_noise_rms` | 4.5 (flat-N² fallback only, see noise templates) |
| `InputPolarity` | `input_polarity` | -1 |
| `WfmFilter` = Wiener, `AutoScale: true` | `auto_scale` | true |
| `ApplyPostfilter`, `WfmPostfilter.Cutoff: 1.5` | `apply_postfilter`, `postfilter_cutoff` | true, 1.5 MHz |
| `ApplyPostBLCorrection` | `apply_post_blcorr` | true |
| `SPETemplateFiles` + `SPETemplateMap*` | `spe_file` (JSON) | `pdhd-spe-templates.json` |
| `NoiseTemplateFiles` + `NoiseTemplateMap*` | `noise_file` (JSON, empty = flat N²) | `pdhd-noise-templates.json` (set by `flash.jsonnet`) |

Algorithm (per snippet, verified bit-equal against an independent numpy
replication, max |Δ| ~6e-7):

1. pedestal = mean of the first `pre_trigger − pedestal_buffer` (= 20) samples;
2. `xv = polarity·(ADC − pedestal)`;
3. Wiener filter `G = conj(H)·S²/(|H|²·S² + N²)` with `H` = FFT of the
   channel's SPE template, `S² = (max(xv)/SPE_amplitude)²` (the module's
   δ-function input guess), `N²` per frequency bin from the channel's noise
   power spectrum (`noise_file`; flat `LineNoiseRMS²·N` when unmapped/absent);
4. auto-scale from the filtered SPE response (integral of the positive region
   around its peak, clamped to ≥1);
5. baseline of the deconvolved waveform from the first 20 samples (pre
   post-filter), subtracted after;
6. Gauss post-filter (cutoff 1.5 MHz at 62.5 MHz sampling).

Differences from the module: a shorter-than-`samples` input is zero-padded
(LArSoft pads with random noise — never triggered for these 1024-tick
snippets); `kScint` input shape is not ported (unused in the PDHD
configuration).

### SPE + noise templates — default is the 2024 averages (NOT v1)

Two template sets exist; **the default `pdhd-spe-templates.json` is the 2024
NP04 FBK/HPK *average* templates** (`extract_pdhd_spe_templates.py`,
`SPE_NP04_{FBK,HPK}_2024_without_pretrigger.dat`) with **flat Wiener N²**
(`noise_file: ''`).  This is a deliberate choice over the LArSoft production
calibration, made after the per-channel study in `pdhd/pics/pd/`:

- the **v1 per-channel** templates (preserved in `pdhd-spe-templates-v1.json`,
  `extract_pdhd_spe_templates_v1.py`; the 113 run28368
  `..._Jun2025.txt` from DUNE StashCache
  `/cvmfs/dune.osgstorage.org/.../ProtoDUNE/HD/opdetresponse/v1/`) reproduce
  the in-file LArSoft `deconv` **exactly** (correlation 1.0000, peak 1.000,
  shift 0) — but the LArSoft v1 deconvolution **over-subtracts the slow tail
  below zero on ~every channel** (their net template area is large-positive,
  not DC-balanced);
- the **2024 averages** are DC-balanced (total area ≈ 0) and deconvolve to a
  **flat tail** (real C++ output: median late-tail ≈ 0 % of peak across the
  example events), at the cost of a +1–2 tick shift vs the v1/LArSoft peak.

So `flash/` deliberately diverges from the LArSoft production deconvolution in
favour of the physically flat 2024 result.  The optional `noise_file`
(run27950 per-channel spectra, `pdhd-noise-templates.json`,
`FFT_Noise_Template_run27950_..._{fbk,hpk}_Dec2024.txt`) is a second-order
effect and slightly less flat with these templates, so it is off by default.
A data-driven SPE template extracted from these snippets does **not** beat the
2024 average (the single-PE pulses sit at the noise floor; see
`pdhd/pics/pd/README.md`).  Switching sets is a one-line `flash.jsonnet`
change (`spe_file` / `noise_file`).

### `Flash::OpHitFinder` (`IFrameTensorSet`, "decon" → ophits tensor)

Port of `OpHitFinderDeco` (`dune_ophit_finder_deco`) essentials:

- deconvolved snippet × `ScalingFactor` (=100) cast to `short`
  (`RunHitFinder_deco`), so all thresholds/areas are in those units;
- `PedAlgoEdges` head method (`standard_algo_pedestal_edges`): pedestal
  mean/σ from the first 3 samples;
- `AlgoSlidingWindow` (direct port, public static `sliding_window()` for unit
  tests): ADCThreshold 3.0, NSigma 1, Tail/End thresholds 1/1, MinPulseWidth 1,
  NumPre/PostSample 2/2, positive polarity;
- `ConstructHit`: keep pulses with peak ≥ `HitThreshold` (3.0);
  `PE = area/SPEArea` (=100); peak/start times in trigger-relative WCT ns via
  the frame time/tick; `flash_id = −1`.

**Overlapping-pulse splitting (extension, not in larana).** A second light pulse
riding on the first's slow scintillation tail never dips below the (very low)
tail threshold, so `sliding_window` merges the two into one OpHit and the second
pulse's PE is absorbed.  `split_pulse()` (public static, post-processes each
found pulse without touching the state machine) splits a pulse at prominent
valleys between its sub-peaks — the per-channel analogue of the prototype
`ToyLightReco` flash-level KS/re-trigger separation.  A valley splits two peaks
only when it is deep both **relatively** (`split_min_prominence` 0.4 of the
smaller flanking peak) and **absolutely** (`split_min_prominence_abs`, scaled
units; PDHD 100 = 1 PE/tick), so a smooth tail is not fragmented; each sub-pulse
recomputes its own area/peak/times.  Gated by `split_enable` — **component
default off / bit-identical** to larana (the disabled path returns the pulse
verbatim; covered by the `split disabled` doctest), **on in PDHD's
`flash.jsonnet`**.  On run 27305 evt 150: 827 → 852 OpHits, 55 → 60 OpFlashes,
PE conserved to 0.01 %.  See `pdhd/docs/pdhd-light-raw-data.md` §3.4.

### `Flash::OpFlashFinder` (`ITensorSetFilter`, ophits → opflash tensor set)

Port of `larana/OpticalDetector/OpFlashAlg.cxx` `RunFlashFinder` with
`protodune_opflash` (= `dunefd_opflash` = `standard_opflash` +
`FlashThreshold: 3.5`) values:

- double-offset 1 µs accumulators (`bin_width`), flash candidates ≥
  `flash_threshold` (3.5 PE), largest-first hit claiming;
- `RefineHitsInFlash` with `width_tolerance` 0.5;
- `ConstructFlash`: PE-weighted flash time, per-OpDet PE vector (OpChannel ==
  OpDet), y/z centroid/width from `pdhd-opdet-geom.json` (mm);
- `RemoveLateLight` (configurable `remove_late_light`, default true; Argon
  constant 1.6 µs).

Not ported: wire-plane centers/widths, Frame/InBeamFrame/OnBeamTime
bookkeeping (LArSoft frame conventions don't apply).

## Running

```
./run_light_evt.sh -m reco  <run> <evt>    # reconstruction only
./run_light_evt.sh -m both  <run> <evt>    # stage-1 convert + reco
```

Outputs in `pdhd/work/<RUN6>_<EVT>/`:

- `opflash_pdhd-wct.tar.gz` — opflash tensor set (`producer: "wct-flash"`),
  same schema/archive convention as the stage-1 converted LArSoft flashes,
  directly consumable by `FlashTensorToOpticalPCs{nchan:160}` /
  `QLMatching{nchan:160}`;
- `light-frames-wct.tar.bz2` — "raw" + "decon" frames for waveform-level
  validation.

Unit tests: `flash/test/doctest_ophitfinder.cxx` (sliding-window pulse
finding); `build/flash/wcdoctest-flash`.
