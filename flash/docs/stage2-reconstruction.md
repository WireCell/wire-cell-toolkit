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
| `LineNoiseRMS` | `line_noise_rms` | 4.5 |
| `InputPolarity` | `input_polarity` | -1 |
| `WfmFilter` = Wiener, `AutoScale: true` | `auto_scale` | true |
| `ApplyPostfilter`, `WfmPostfilter.Cutoff: 1.5` | `apply_postfilter`, `postfilter_cutoff` | true, 1.5 MHz |
| `ApplyPostBLCorrection` | `apply_post_blcorr` | true |
| `SPETemplateFiles` + `SPETemplateMap*` | `spe_file` (JSON) | `pdhd-spe-templates.json` |

Algorithm (per snippet, verified bit-equal against an independent numpy
replication, max |Δ| ~6e-7):

1. pedestal = mean of the first `pre_trigger − pedestal_buffer` (= 20) samples;
2. `xv = polarity·(ADC − pedestal)`;
3. Wiener filter `G = conj(H)·S²/(|H|²·S² + N²)` with `H` = FFT of the
   channel's SPE template, `S² = (max(xv)/SPE_amplitude)²` (the module's
   δ-function input guess), `N² = LineNoiseRMS²·N`;
4. auto-scale from the filtered SPE response (integral of the positive region
   around its peak, clamped to ≥1);
5. baseline of the deconvolved waveform from the first 20 samples (pre
   post-filter), subtracted after;
6. Gauss post-filter (cutoff 1.5 MHz at 62.5 MHz sampling).

Differences from the module: a shorter-than-`samples` input is zero-padded
(LArSoft pads with random noise — never triggered for these 1024-tick
snippets); `kScint` input shape and noise templates are not ported (unused in
the PDHD configuration).

### SPE templates — known fidelity limit

`cfg/pgrapher/experiment/pdhd/pdhd-spe-templates.json`
(`flash/test/extract_pdhd_spe_templates.py`) carries the two official 2024
templates `SPE_NP04_{FBK,HPK}_2024_without_pretrigger.dat` (duneopdet
`config_data/`) with the FBK/HPK channel split of
`protodunehd_pds_channels_mc` (`dune_opdet_channels.fcl`).

The in-file reference `deconv` waveforms were however produced with the
**per-channel run28368 v1 templates** (`protodunehd_template_list_v1`,
`ProtoDUNE/HD/opdetresponse/v1/run28368_APA4_CH*_Jun2025.txt`), which are not
available on the accessible cvmfs products.  Consequence (run 27305 evt 150):
our deconvolution agrees with the reference at correlation 0.99 but with a
~2-tick peak shift, ~10% amplitude difference and a residual ~0.1 PE/tick
baseline, which inflates downstream hit areas and counts.  Swapping in the v1
templates when obtainable only requires regenerating the JSON (the
channel→template map supports per-channel templates as-is).

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
