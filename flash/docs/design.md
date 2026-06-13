# PDHD Light (Photon Detector) Reconstruction — Design

Status: stages 1–3 complete — ROOT → WCT conversion
(stage1-root-conversion.md), WCT-native DUNE-method reconstruction
(stage2-reconstruction.md), validation vs the LArSoft reference
(validation.md).  This document is the living design record for
the toolkit-native PDHD light-signal-processing subsystem.  It will be updated at each
implementation stage; per-stage details land in sibling documents:

- `stage1-root-conversion.md` — ROOT → WCT format conversion (I/O components, dump job)
- `stage2-reconstruction.md` — WCT-native optical reconstruction (decon → ophit → opflash)
- `validation.md` — comparison against the LArSoft reference products

## 1. Goal and scope

Reconstruct PDHD photon-detector (PDS) light from DAPHNE waveforms inside the
wire-cell toolkit, ending in **flash** objects that the existing SBND-style
charge-light (Q-L) matching machinery (`FlashTensorToOpticalPCs` → `QLMatching`)
consumes without modification.  The method template is the DUNE LArSoft optical chain
(see `pdhd/docs/photon-detector-chain.md`):

```
raw OpDetWaveform → deconvolution (Wiener + SPE template)   [duneopdet Deconvolution]
                 → OpHit finding                            [dune_ophit_finder_deco]
                 → OpFlash assembly                         [protodune_opflash]
```

Package layout:

- `root/` — ROOT-based file I/O **only**.  The current ROOT files are a *temporary*
  exchange format; everything ROOT-specific is isolated in two source components so a
  future format change only replaces the sources, nothing downstream.
- `flash/` (this package) — detector-file-format-agnostic light reconstruction
  components operating on `IFrame` / `ITensorSet`.
- `cfg/pgrapher/experiment/pdhd/` — configuration (jsonnet library + data JSONs).
- `pdhd/` (separate `wcp-porting-img` repo) — job entry points, run scripts, work area.

## 2. Input data (temporary ROOT format)

Example files: `pdhd/input_data_7p8_new_coh_grouping/run<RUN6>/np04hd_raw_run<RUN6>_*.root`
(runs 27305, 27980, 28084, 29107 with 24–31 events each).  Each file holds four
TDirectories:

| Directory | Content |
|---|---|
| `opflashana/` | LArSoft `opflashana` analyzer trees: `PerOpHitTree` (OpHits), `PerFlashTree` (flash summaries), `FlashBreakdownTree` (per-flash per-OpChannel PE), `PerEventFlashTree` (per-event vectors incl. flat `[nflash×256]` PE matrix), `FlashHitMatchTree` (flash↔hit assoc.) |
| `flashopdet/` | `opdet_geo` (160 OpDets: ident + x,y,z in **cm**), `flash_opdet` (per-flash **per-OpDet** PE with positions), `opch_map` (OpChannel→OpDet; **only present in `onevent_run27305_final.root`**, 224 of 256 channels, 4 ganged OpChannels per OpDet window) |
| `trigoff/` | `trigger_offset`: per event `rd_timestamp` (readout start, DTS 16 ns ticks), `tc_time_candidate` (trigger candidate), `tc_type` (14/29/31), `offset_us = (tc − rd)·16 ns ≈ 249.7–250.1 µs`.  May contain multiple candidates per event. |
| `decoana/` | `run_<R>_evt_<E>/ch<N>/{raw,deconv}/waveform_<i>`: TH1D snippets, 1024 bins @ 16 ns, x-axis in µs on the absolute DTS clock.  `raw` = ADC (pedestal ≈ 8180, **negative pulse polarity**); `deconv` = SPE-normalised output of the LArSoft `Deconvolution` module — the golden reference for our stage-2 port.  Only self-triggered channels appear (~59 per event, OpChannel 0–82 in the examples). |

Key numerology:

- 160 physical OpDets (X-ARAPUCA acceptance windows: 4 APAs × 10 bars × 4 windows),
  positions known for all 160 from `opdet_geo`.
- 256 electronics **OpChannel** slots (224 mapped in the available map); 4 OpChannels
  gang onto one OpDet window position.
- Optical clock: 62.5 MHz → 16 ns ticks (DTS timestamps count these ticks).
- TPC readout window opens ≈250 µs *before* the trigger (`offset_us`); this is the
  same −250 µs anchor as `protodunehd_detectorclocks.TriggerOffsetTPC` and the
  commented `time_offset` in `cfg/pgrapher/experiment/pdhd/clus.jsonnet`.

## 3. Conventions

### 3.1 Time

All times in WCT products are **ns (WCT units) relative to the event trigger**:

- Trigger anchor: `tc_time_candidate` from `trigoff/trigger_offset`.  When an event has
  several candidates, the one whose `offset_us` is closest to a configurable
  `nominal_offset_us` (default 250.0) is used; the chosen `tc_type`, `offset_us` and
  the DTS tick values are recorded in output metadata so the choice is auditable.
- Flash time: `t = (FlashTime_dts − tc_dts) · 16 ns`.
- Waveform frame: `IFrame::time() = (rd_timestamp − tc_dts) · 16 ns ≈ −250 µs`
  (readout start relative to trigger); each trace `tbin` = snippet start relative to
  `rd_timestamp` in 16 ns ticks; `IFrame::tick() = 16 ns`.
- Relation to TPC: charge imaging runs with `time_offset = 0` (the raw `x` scope stays
  offset-free; see `clus.jsonnet`).  Rather than baking the −250 µs anchor into `x_raw`
  at imaging time, the per-event `offset_us` is applied **downstream** at Q/L matching:
  `QLMatching` folds it into the per-flash drift `x` correction (`trigger_offset` config)
  and `T0Correction` adds it to the `x_t0cor` scope (`trigger_offset` DV-metadata key).
  Both default 0 ⇒ bit-identical when unused.  See `match/docs/joint-qlmatching-design.md`
  and `clus/docs/clustering_with_t0.md`.

### 3.2 Channels

- **Flash-level products use the 160 OpDet (geometry) basis** (`nchan = 160`):
  Q-L matching needs a per-channel position (semi-analytical light model), which only
  exists per OpDet.  PE from ganged OpChannels is summed onto their OpDet.
- **Waveform-level products keep the OpChannel (electronics) basis** — faithful to the
  data; the OpChannel→OpDet reduction happens once, in flash assembly.
- Maps are committed as JSON config data (extracted once from the ROOT files):
  - `cfg/pgrapher/experiment/pdhd/pdhd-opdet-geom.json` — 160 OpDets, positions in mm;
  - `cfg/pgrapher/experiment/pdhd/pdhd-opch-map.json` — OpChannel→OpDet (with a
    coverage caveat: 224/256 channels, from the one file carrying the map).

### 3.3 Data model

No new `iface` IData classes.  Light data uses existing toolkit currency:

- waveforms: `IFrame` with tagged trace sets (`"raw"`, `"deconv"`/`"decon"`),
  trace channel = OpChannel, multiple snippets per channel allowed;
- hits and flashes: `ITensorSet` with the documented schema below;
- at Q-L matching time the opflash tensor is turned into the canonical
  `flash`/`light`/`flashlight` point clouds by the existing
  `FlashTensorToOpticalPCs` (`aux/`), configured with `nchan: 160`.

### 3.4 Opflash tensor-set schema

One schema shared by stage 1 (converted LArSoft flashes) and stage 2 (WCT-native
reconstruction) so every downstream consumer is identical:

- tensor 0 (`opflash`): f8, shape `[nflash, 1+160]`, row-major;
  column 0 = flash time (ns, trigger-relative), columns 1..160 = PE per OpDet.
  This is exactly the `FlashTensorToOpticalPCs` input contract.  Flashes may be
  built per cathode side (`group_by_side`, on for PDHD; see
  `stage2-reconstruction.md`) — the schema is unchanged, the side is implicit in
  which OpDet columns carry PE.
- tensor 1 (`flash_summary`): f8, shape `[nflash, 8]`:
  `flash_id, total_pe, y_center, z_center, y_width, z_width, abs_time_dts, nhits`
  (lengths in mm; `abs_time_dts` in 16 ns DTS ticks as double; `nhits = -1` if unknown).
- tensor 2 (`ophits`, optional): f8, shape `[nhit, 9]`:
  `channel, peak_time_ns, width_ns, area, amplitude, pe, start_time_ns, flash_id, fast_to_total`
  (`channel` = OpChannel; `flash_id = -1` when unassociated).  With OpHit splitting
  enabled (PDHD; see `stage2-reconstruction.md`) a row may be one sub-pulse of a
  merged pulse — the schema is unchanged, there are just more rows.
- tensor-set metadata: `run, subrun, event, tc_type, tc_time_dts (string),
  rd_timestamp_dts (string), offset_us, producer` (`"opflashana"` for converted
  LArSoft products, `"wct-flash"` for our reconstruction); tensor-set `ident` = event.
  The stage-1 converter writes the full set; the stage-2 `"wct-flash"` producer writes
  `event, nchan, producer` plus `offset_us` (the WCT-native reco can't recover the
  readout-vs-trigger offset from the self-triggered snippets, so `OpFlashFinder` takes
  it as a config value — supplied by `run_light_evt.sh` from the ROOT `trigoff` tree —
  and stamps it verbatim).  Downstream Q/L matching reads this `offset_us`.
- On disk: `TensorFileSink` archive (`opflash_*.tar.gz`, prefix `opflash_`), the same
  convention as SBND's `opflash_apa<n>.tar.gz`, so
  `TensorFileSource → FlashTensorToOpticalPCs{nchan:160} → QLMatching{nchan:160}`
  consume it unmodified.

## 4. Components

### Stage 1 — ROOT → WCT conversion (`root/`)

- `Root::PDHDOpFlashSource` (`ITensorSetSource`): `flash_opdet` + `PerFlashTree` +
  `PerOpHitTree` + `trigoff` → one opflash tensor set per configured event.
- `Root::PDHDOpWaveformSource` (`IFrameSource`): `decoana` TH1D snippets → `IFrame`
  with `"raw"` and `"deconv"` tagged traces (raw ADC values kept faithful).

A dump job (`pdhd/wct-light-convert.jsonnet` + `pdhd/run_light_evt.sh`) writes
`opflash_pdhd.tar.gz` (TensorFileSink) and `light-frames.tar.bz2` (FrameFileSink)
into `pdhd/work/<RUN>_<EVT>/`.

### Stage 2 — reconstruction (`flash/`)

Port of the DUNE chain with parameters from
`/cvmfs/dune.opensciencegrid.org/products/dune/duneopdet/v10_20_09d00/`:

- `Flash::OpDecon` (`IFrameFilter`) — `Deconvolution_module.cc` math: pedestal
  subtraction, polarity flip, FFT, Wiener filter built from the per-channel
  run28368 v1 SPE templates and the run27950 noise power spectra (the
  `protodunehd_pds_channels_data_v1` production set; flat `LineNoiseRMS²·N`
  as fallback when no noise file/channel is mapped), auto-scale,
  Gauss post-filter (cutoff 1.5 MHz), post baseline correction.
  PDHD fcl values: `Samples 1024, Pedestal 8180, PreTrigger 50, PedestalBuffer 30,
  InputPolarity −1`.
- `Flash::OpHitFinder` (`IFrameTensorSet`) — `dune_ophit_finder_deco` essentials:
  SlidingWindow hit finding on deconvolved snippets, `HitThreshold 3.0` (Wiener),
  PE = area / `ScalingFactor 100`.
- `Flash::OpFlashFinder` (`ITensorSetFilter`) — `protodune_opflash` essentials:
  time clustering of hits into flashes, per-OpDet PE summation via the channel map,
  y/z centroid/width from OpDet positions; emits the §3.4 schema.

### Stage 3 — validation

Compare, on common channels, our decon vs the in-file `deconv` snippets (same
algorithm and parameters ⇒ near-identical expected), our hits vs `PerOpHitTree`, our
flashes vs `PerFlashTree`/`flash_opdet`.  Results in `validation.md`.

## 5. Risks / open items

1. **decoana x-axis semantics** — RESOLVED in stage 1: the axis low edge encodes
   the snippet start in 16 ns ticks relative to the event's earliest saved snippet
   (`t_first`, not recorded in the file; recovered per event from OpHit↔deconv-peak
   coincidence voting — see stage1-root-conversion.md).  All absolute timestamps in
   the dump are doubles with a 16-tick ULP, so converted times carry a ±8-tick
   (±128 ns) quantization.
2. **Multiple trigger candidates** — nearest-to-250 µs selection is a heuristic;
   metadata records the choice.  Physics-trigger `tc_type` may need expert input
   (run 29107 mixes tc_type 14/29/31; the selection picks 14 and validates).
3. **opch_map coverage** — present in one file only, 224/256 channels (~54 OpDets).
   Stage-2 flashes will under-populate the 160-OpDet matrix relative to `opflashana`
   on events with unmapped active channels; unmapped channels warn-and-drop.
4. **SPE-template↔channel mapping** — RESOLVED in stage 3: the production
   `protodunehd_pds_channels_data_v1` set (run28368 v1 per-channel SPE
   templates + run27950 noise spectra, retrieved from the DUNE StashCache
   cvmfs) is carried in `pdhd-spe-templates.json` / `pdhd-noise-templates.json`
   (`extract_pdhd_spe_templates_v1.py`); deconvolution matches the reference
   exactly.
5. **`protodune_opflash` parameter fidelity** — exact `dunefd_opflash` values are
   copied at implementation time; tuned only if flash-count comparisons demand.
6. **Future data format** — when the temporary ROOT format is replaced, only the two
   `root/` sources are rewritten; `flash/` components and all schemas are unchanged.
