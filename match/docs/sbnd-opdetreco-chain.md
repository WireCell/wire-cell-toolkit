# SBND optical (PMT) reconstruction chain: raw waveform → OpHit NPE → OpFlash

How `sbndcode/OpDetReco/` (plus the legacy raw-path finder in `OpDetSim/`) turns a
digitized `raw::OpDetWaveform` into reconstructed photo-electrons (`recob::OpHit::PE()`)
and time-clustered `recob::OpFlash`es.

> Source tree: `sbndcode/sbndcode/OpDetReco/` and
> `sbndcode/sbndcode/OpDetSim/opHitFinderSBND_module.cc` (linked checkout at `./sbndcode`).
> Companion to [`sbnd-opdetsim-chain.md`](sbnd-opdetsim-chain.md) (the simulation side).

## The chain at a glance

```
raw::OpDetWaveform (ADC)
   │  SBNDOpDeconvolution  (Wiener)        ← OpDetReco/OpDeconvolution
   ▼
deconvolved OpDetWaveform (PE-density)
   │  SBNDOpHitFinder      (SlidingWindow) ← OpDetReco/OpHit
   ▼
recob::OpHit  (PE = Area / SPEArea)
   │  SBNDFlashFinder      (SimpleFlashAlgo) ← OpDetReco/OpFlash
   ▼
recob::OpFlash (total PE, time, Y/Z barycenter)
```

The default production fcl `OpDeconvolution/job/run_decohitfinder.fcl` runs the producer
path **`opdecopmt → ophitpmt → opflashtpc0 / opflashtpc1`** (the X-ARAPUCA siblings exist
but are commented out by default).

There is **also a legacy direct path** that skips deconvolution and finds hits on raw
waveforms (`OpDetSim/opHitFinderSBND_module.cc`); both compute a PE, with different
calibration constants. Both are documented below.

## Stage 1 — Deconvolution (raw ADC → PE-density waveform)

`OpDeconvolution/SBNDOpDeconvolution_module.cc` reads raw `OpDetWaveform`s (label
`opdaq`), filters by PD type + electronics, and calls a Wiener-filter tool:

- `OpDeconvolutionAlgWiener_tool.cc` (MC) / `OpDeconvolutionAlgWienerData_tool.cc` (data,
  per-channel SPE templates + filter parameters from the calibration DB).
- Loads the **single-PE response (SPR)** template from `OpDetSim/digi_pmt_sbnd_v2int0.root`.
- Builds a light-signal hypothesis (2-exponential scintillation, fast/slow constants
  `HypoSignalTimeConsts: [25, 1500]` ns, weights `[0.3, 0.7]`).
- Per waveform: estimate baseline mean + RMS from the front samples, build a Wiener kernel
  (SPR ⊗ signal-hypothesis vs noise PSD), FFT-convolve, then rescale by
  `DecoWaveformPrecision`.

Output: a deconvolved `OpDetWaveform` whose amplitude is ~photo-electron density rather
than ADC. This is the input to the deco-path hit finder.

## Stage 2 — OpHit finding and the PE (NPE) formula

### Deco path (default) — `OpDetReco/OpHit/SBNDOpHitFinder_module.cc`

Wraps the standard larana `RunHitFinder` with a `SlidingWindow` pulse algorithm on the
deconvolved waveform. The PE conversion is a **standard photon calibrator**:

```cpp
// SBNDOpHitFinder_module.cc:149-161
bool  areaToPE = pset.get<bool> ("AreaToPE");
float SPEArea  = pset.get<float>("SPEArea");
...
fCalib = new calib::PhotonCalibratorStandard(SPEArea, SPEShift, areaToPE);
```

With `AreaToPE = true`, larana computes **`PE = OpHit.Area / SPEArea`**, where `Area` is
the integral of the (deconvolved) waveform over the hit window. Config
(`OpDeconvolution/job/sbnd_ophitfinder_deco.fcl`):

| param | PMT | X-ARAPUCA |
|-------|-----|-----------|
| `SPEArea` | 200 | 500 |
| `InputModule` | `opdecopmt` | `opdecoxarapuca` |
| `HitAlgoPset.Name` | `SlidingWindow` | `SlidingWindow` |
| `HitAlgoPset.ADCThreshold` (start / end) | 12 / 8 | 20 / 20 |
| `PedAlgoPset` | Edges (200/200) | Edges (50/50) |

(The data variant `…PMT_data` adds a per-channel `ADCThresholdVector` of 312 entries.)

### Raw path (legacy) — `OpDetSim/opHitFinderSBND_module.cc`

Operates directly on raw `OpDetWaveform`s (no deconvolution):

1. `subtractBaseline()` — average the first `BaselineSample` (95) ticks, subtract, and
   apply polarity (`PulsePolarityPMT: -1`, `PulsePolarityArapuca: +1`).
2. `findAndSuppressPeak()` — find the max sample; if above threshold, walk out to the
   threshold-crossing edges and integrate; convert ADC·ticks → **ADC·ns** by dividing by
   the sampling rate (`/ fSampling`, or `/ fSampling_Daphne` for DAPHNE electronics).
3. PE from a fixed SPE area (`opHitFinderSBND_module.cc:201-204`):

```cpp
if (pmt_coated || pmt_uncoated) phelec = Area / fArea1pePMT;   // 132.66 ADC*ns
else /* xarapuca */             phelec = Area / fArea1peSiPM;  // 8541.  ADC*ns
```

Config (`OpDetSim/ophitfinder_sbnd.fcl`): `Area1pePMT: 132.66`, `Area1peSiPM: 8541`,
`ThresholdPMT: 8`, `ThresholdArapuca: 20`, `BaselineSample: 95`.

### The two paths side by side

| aspect | deco path (default) | raw path (legacy) |
|--------|--------------------|--------------------|
| input | deconvolved waveform | raw ADC waveform |
| algorithm | larana `SlidingWindow` | local peak-find + suppress |
| `Area` units | deconvolved (~PE) | ADC·ns |
| PE = | `Area / SPEArea` (200 / 500) | `Area / Area1pe{PMT,SiPM}` (132.66 / 8541) |
| baseline | larana `Edges` ped algo | front-95-tick mean |

In both, **NPE is an integral-over-a-calibration-area**: integrate the pulse, divide by
the per-detector-class single-PE area.

## Stage 3 — OpFlash (cluster OpHits in time)

`OpFlash/SBNDFlashFinder_module.cc` collects OpHits, re-calibrates each to PE, and clusters
them in time:

- **Per-channel recalibration** (`OpFlash/FlashFinder/PECalib.cxx:48`):
  ```cpp
  double area_pe = area / _spe_area_gain_v[opdet] * _relative_qe_v[opdet];
  ```
  i.e. `Area / SPEAreaGain × relativeQE`, with `SPEAreaGain` = 200 (PMT) / 500 (Arapuca)
  and a per-channel relative-QE vector (defaults to 1.0).
- **Time clustering**: `SimpleFlashAlgo` groups hits inside an integration window subject
  to min-PE / coincidence cuts.
- **Geometry**: `FlashTools/FlashGeoBarycenter_tool.cc` — PE-weighted Y/Z barycenter and
  RMS widths over the channels' positions.
- **Drift / propagation**: `FlashTools/DriftEstimatorPMTRatio_tool.cc` estimates drift from
  the PMT/Arapuca PE ratio and corrects the light-propagation time.
- **Output**: `recob::OpFlash` (total PE, time + width, Y/Z center + width) plus
  OpHit↔OpFlash associations; produced per TPC (`opflashtpc0`, `opflashtpc1`).

## Standard job ordering

`OpDeconvolution/job/run_decohitfinder.fcl` (MC) includes `opdeconvolution_sbnd.fcl`,
`sbnd_ophitfinder_deco.fcl`, `sbnd_flashfinder_deco.fcl` and runs:

```
reco: [ opdecopmt, ophitpmt, opflashtpc0, opflashtpc1 ]
```

(`run_decohitfinder_data.fcl` is the data analogue using the `_data` tool variants.)

## Where this lives in code

- Deconvolution: `OpDetReco/OpDeconvolution/SBNDOpDeconvolution_module.cc`,
  `Alg/OpDeconvolutionAlgWiener_tool.cc`, `Alg/OpDeconvolutionAlgWienerData_tool.cc`.
- OpHit (deco): `OpDetReco/OpHit/SBNDOpHitFinder_module.cc`,
  `OpDeconvolution/job/sbnd_ophitfinder_deco.fcl`.
- OpHit (raw): `OpDetSim/opHitFinderSBND_module.cc`, `OpDetSim/ophitfinder_sbnd.fcl`.
- OpFlash: `OpDetReco/OpFlash/SBNDFlashFinder_module.cc`,
  `FlashFinder/SimpleFlashAlgo.*`, `FlashFinder/PECalib.cxx`, `FlashTools/*`.
- Jobs: `OpDetReco/OpDeconvolution/job/run_decohitfinder.fcl`.

## Known failure mode: closely-spaced flashes (≲8 µs) and cross-TPC inconsistency

When two genuine flashes arrive within ~8 µs of each other, `SimpleFlashAlgo` cannot keep
both **in the same TPC**, and — because the finder runs **independently per TPC with no
cross-TPC reconciliation** — the two TPCs can end up keeping *different* members of the
pair. A cross-TPC object that physically shares one flash time then matches two
non-coincident flashes and fails cross-TPC (xTPC) matching downstream.

### Symptom (data event 59415)

Two cluster halves of a cathode-crossing cosmic — at `(x,y,z) = (31.8, -57.0, 451.1)` and
`(-2.5, -81.2, 450.1)`, i.e. nearly identical y/z, opposite sides of the cathode — did **not**
get xTPC matched. Each half matched a different flash, and the two flashes are **not
coincident in time**, even though the geometry says both halves belong to the same flash
group.

### Root cause

The finder is `SBNDFlashFinder_module.cc` driving `SimpleFlashAlgo::RecoFlash`
(`OpFlash/FlashFinder/SimpleFlashAlgo.cxx`). Four structural facts combine:

1. **Per-TPC, independent.** `run_flashfinder.fcl:60-61` instantiates two producers,
   `opflashtpc0` (`SimpleFlashTPC0.TPC: 0`) and `opflashtpc1` (`TPC: 1`); each algorithm
   only sees its own TPC's PMTs (`Configure()` → `ListOpChannelsByTPC(_tpc)`,
   SimpleFlashAlgo.cxx:63-72). The two `recob::OpFlash` collections are never reconciled in
   time anywhere in sbndcode.

2. **Hard veto, brightest-first.** Candidate time-bins are processed in **descending PE
   order** (`pesum_idx_map` keyed by `1./pe`, SimpleFlashAlgo.cxx:190-197). The brightest
   claims a flash; every later candidate whose start lies within
   `veto_ctr = VetoSize / TimeResolution` of an already-accepted flash is **skipped entirely**
   (SimpleFlashAlgo.cxx:211-250, the `skip=true` branches at lines 224-231). Defaults
   (`job/sbnd_flashalgo.fcl`): `VetoSize = 8 µs`, `IntegralTime = 8 µs`,
   `TimeResolution = 0.01 µs`, `PEThreshold = 20`, `MinPECoinc = 6`, `MinMultCoinc = 3`.
   So **within one TPC, two flashes < 8 µs apart cannot both exist** — only the single
   brightest survives.

3. **The bug.** For two real flashes within 8 µs, each TPC keeps only whichever flash is
   brighter **in that TPC**. If flash A is brighter in TPC0 and flash B is brighter in TPC1,
   TPC0 keeps A and vetoes B's coincident partner; TPC1 keeps B and vetoes A's. The two
   surviving flash times then differ by the A–B separation, so the cross-TPC pair splits onto
   two non-coincident flashes. (Secondary effect: the surviving flash's `IntegralTime = 8 µs`
   window — `pesum` loop at SimpleFlashAlgo.cxx:253-256 and PE fill at 278-284 — also absorbs
   the vetoed partner's PE, corrupting total PE and the per-channel pattern used for the Y/Z
   barycenter and the downstream Q/L predicted light.) Note also the flash time is the single
   **brightest 10 ns bin**, not a PE-weighted centroid (`flash_time_v` = `idx`,
   SimpleFlashAlgo.cxx:266, 310-314).

4. **No spatial discrimination.** Time-bins are cut on PE-sum + hit-multiplicity only
   (SimpleFlashAlgo.cxx:191-197). There is **no analog of MicroBooNE's per-PMT
   Kolmogorov–Smirnov spatial-pattern test** — the very information (different PMT
   illumination of A vs B) that would distinguish two close flashes is discarded.

### Evidence — reconstructed flashes for event 59415

From `opflash_apa{0,1}.tar.gz` (col 0 = time [ns], cols 1..312 = per-channel PE), the region
near −1.26 ms:

| flash | TPC0 (apa0) | TPC1 (apa1) |
|-------|-------------|-------------|
| near −1.2688 ms | **−1268808 ns, 22712 PE** | absent (nearest −1259707, 212 PE) |
| near −1.2640 ms | absent | **−1264017 ns, 12394 PE** |

The two bright flashes are **4.79 µs apart**, each bright in exactly one TPC and **missing in
the other**. TPC0 has no flash at −1264017 (it is 4.79 µs < 8 µs from its brighter −1268808 →
vetoed); TPC1 has no flash at −1268808 (same reason). This veto signature is the direct cause
of the 59415 xTPC failure.

### Contrast with the MicroBooNE reference (mechanism, not threshold)

`prototype_base/2dtoy/src/ToyLightReco.cxx` (`Process_beam_wfs`) bins 1500 samples → 250 bins
(6 samples/bin ≈ 93.75 ns/bin at 64 MHz). Its 78-bin auto-new-flash gap / default window is
**≈ 7.3 µs** — *comparable* to SBND's 8 µs, so the real difference is the **mechanism**:

- MicroBooNE separates close-in-time flashes with **gap + brightening + a per-PMT KS
  spatial-pattern test** (the 4–15 bin / ~0.4–1.4 µs regime), truncates each flash window at
  the next flash's start (no PE double-counting), and runs on **one coherent PMT system**, so
  a flash has a single consistent time.
- SBND uses a **flat PE+multiplicity threshold + a hard veto, independently per TPC**, so it
  (a) cannot split two close flashes inside a TPC and (b) produces **inconsistent cross-TPC
  flash times** for one physical event.

### Possible remedies (not implemented — diagnostic note only)

- **Reconcile TPC0 + TPC1 flashes into cross-TPC flash groups** by coincident time before
  matching (cleanest fix; mirrors how the toolkit already groups via `flash_t0_window`). This
  fixes the xTPC symptom without touching the per-TPC time clustering.
- **Reduce `VetoSize` / `IntegralTime`** so two flashes a few µs apart can both be
  reconstructed within a TPC.
- **Add a per-PMT spatial (KS-style) test** to split close-in-time flashes that differ in
  spatial pattern, à la MicroBooNE.
- Trade-off: a smaller veto risks splitting a single flash's late (slow-scintillation) light
  into two spurious flashes, so any veto reduction needs validation against single-flash late
  light.
