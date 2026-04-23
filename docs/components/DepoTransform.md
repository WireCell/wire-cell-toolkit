# DepoTransform

Converts a set of drifted depositions into a simulated readout frame by convolving diffused charge with per-plane field and electronics response functions using FFT-based transforms.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoTransform` |
| Concrete class | `WireCell::Gen::DepoTransform` |
| Node category | function |
| Primary interface | `IDepoFramer` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type-name of the IAnodePlane component |
| `nsigma` | number of Gaussian sigma at which to truncate diffusion (default: 3.0) |
| `fluctuate` | if true apply charge fluctuation (default: false) |
| `rng` | type-name of the IRandom component used when fluctuate is true |
| `dft` | type-name of the IDFT component for FFT transforms (default: "FftwDFT") |
| `start_time` | readout gate start time in WCT time units (default: 0 ns) |
| `readout_time` | readout gate duration in WCT time units (default: 5 ms) |
| `tick` | waveform sample period in WCT time units (default: 0.5 us) |
| `drift_speed` | nominal electron drift speed (default: 1.0 mm/us) |
| `first_frame_number` | starting frame counter value (default: 0) |
| `pirs` | JSON array of type-names of IPlaneImpactResponse components, one per plane |
| `process_planes` | optional array of plane indices to process; defaults to all three |
