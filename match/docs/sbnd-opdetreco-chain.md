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
