# Signal Processing Pipeline Overview

## What This Module Does

The `sigproc` module implements the full signal processing chain for Liquid Argon Time
Projection Chamber (LArTPC) wire readout data. It converts raw ADC waveforms from
thousands of sense wires into calibrated charge signals suitable for 3D reconstruction.

Reference: arXiv:1802.08709 (MicroBooNE collaboration).

## Pipeline Architecture

A typical signal processing job runs two major stages in sequence:

```
Raw ADC Frame
    |
    v
[Noise Filtering]  -- OmnibusNoiseFilter / OmnibusPMTNoiseFilter
    |
    v
Clean ADC Frame
    |
    v
[Signal Processing] -- OmnibusSigProc
    |
    v
Deconvolved Charge Frame (with ROI masks)
```

### Stage 1: Noise Filtering (OmnibusNoiseFilter)

Three-pass architecture:

1. **Per-channel filters**: Applied independently to each wire. Typical operations:
   - Stuck-bit / sticky-code mitigation (ProtoDUNE)
   - Chirp noise detection (MicroBooNE)
   - RC undershoot correction (spectral division by RC response)
   - Spectral noise mask application (notch filters for known noise frequencies)
   - Robust baseline subtraction
   - Noisy channel identification (RMS-based)

2. **Grouped filters**: Applied to groups of channels sharing electronics
   (e.g., channels on the same ASIC/motherboard). Primary operation:
   - Coherent noise subtraction: compute median waveform across group,
     subtract it from each channel (with signal protection)

3. **Per-channel status filters**: Final channel quality determination.

An optional **PMT noise filter** (OmnibusPMTNoiseFilter) removes cross-talk from
photomultiplier tubes into TPC wires by identifying large negative excursions on
collection wires and correlated signals on induction wires.

### Stage 2: Signal Processing (OmnibusSigProc)

Multi-step pipeline per wire plane (U induction, V induction, W collection):

1. **Response initialization**: Compute average field response convolved with
   electronics response, redigitized to the frame's tick period.

2. **Load and rebase**: Copy trace data into padded Eigen matrices, optional
   linear baseline subtraction.

3. **2D deconvolution** (`decon_2D_init`): The core signal recovery step.
   - FFT in time (per wire) -> per-channel response correction
   - FFT in wire direction -> divide by 2D field response
   - Apply wire-direction filter -> IFFT(wire) -> IFFT(time)
   - Time and wire shift corrections

4. **ROI formation**: Identify signal regions in the deconvolved waveforms.
   - **Tight ROIs**: High-threshold regions (clear signals)
   - **Loose ROIs**: Low-threshold regions (potential signals, wider coverage)
   - Multi-plane geometric consistency checks

5. **ROI refinement**: Iterative multi-stage cleanup:
   - `CleanUpROIs`: Remove isolated ROIs not connected to tight ROIs (BFS flood-fill)
   - `generate_merge_ROIs`: Promote uncovered tight ROIs
   - `MP3ROI`/`MP2ROI`: Multi-plane geometric protection (3-plane and 2-plane consistency)
   - `BreakROIs`: Peak finding to split wide ROIs (TSpectrum-like algorithm)
   - `ShrinkROIs`: Tighten ROI boundaries using neighbor evidence
   - `CleanUpCollectionROIs`/`CleanUpInductionROIs`: Remove fake signals
   - `ExtendROIs`: Extend to cover neighbor ROI time ranges

6. **Final filtering and output**: Apply Wiener and Gaussian filters within ROIs,
   produce output traces with charge estimates.

### Optional: L1SP Filter

For **shorted wire regions** where induction and collection plane wires are
electrically connected, standard deconvolution fails. The L1SP filter uses
L1-norm regularized sparse deconvolution (LASSO/compressed sensing) to jointly
decompose the mixed signal into collection and induction components:

  minimize ||G * beta - W||_2^2 + lambda * ||beta||_1

where G encodes both response functions and the L1 penalty enforces sparsity
(most time bins have zero signal).

## Key Data Structures

| Structure | Role |
|-----------|------|
| `IFrame` / `SimpleFrame` | Container of traces (waveforms) with metadata |
| `ITrace` / `SimpleTrace` | Single channel waveform: channel ID + tick offset + samples |
| `Eigen::Array2D` (m_r_data) | 2D matrix [wire x tick] for per-plane processing |
| `SignalROI` | Region of interest: channel, time range, contained charge |
| `SignalROIList` | Per-wire list of SignalROIs |
| `front_rois` / `back_rois` | Maps tracking wire-to-wire ROI adjacency |
| `contained_rois` | Map of tight ROIs contained within loose ROIs |
| `IChannelNoiseDatabase` | Per-channel noise parameters and filter spectra |

## Detector-Specific Implementations

Each detector has its own noise filtering functions registered as WCT components:

| Detector | Key Noise Patterns | Unique Features |
|----------|--------------------|-----------------|
| MicroBooNE | Chirp, ADC bit shift, partial RC, coherent | Foundational implementation |
| ProtoDUNE-SP | Sticky codes, ledge artifacts, 50kHz harmonic, FEMB clock | FFT-based sticky repair, relative gain calibration |
| ProtoDUNE-HD | Baseline drift, FEMB negative pulses | Wide adaptive baseline (512 ticks) |
| ProtoDUNE-VD | Shield coupling, FEMB negative pulses | Novel shield coupling subtraction |
| DUNE CRP | RC undershoot, partial RC | Minimal/clean implementation |
| ICARUS | RC undershoot only | Simplest; uses IWaveform for RC response |

## Configuration

The pipeline is configured via Jsonnet (a JSON-generating language). Key
configuration entry points:

- `adc-noise-sig.jsonnet`: Full pipeline test (ADC voltage + noise + NF + SP)
- `check_pdsp_sim_sp.jsonnet`: ProtoDUNE simulation with signal processing

Configuration specifies:
- Detector geometry (anode, field response files)
- Noise database parameters (per-channel gains, filters, masks)
- Signal processing parameters (thresholds, ROI parameters, filter types)
- Pipeline graph topology (which components connect to which)
