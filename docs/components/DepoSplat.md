# DepoSplat

A simplified ductor that converts drifted depositions into wire-plane signal frames using Gaussian diffusion splatting without applying field or electronics response functions.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoSplat` |
| Concrete class | `WireCell::Gen::DepoSplat` |
| Node category | queuedout |
| Primary interface | `IDuctor` |
| Input type(s) | `IDepo` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type-name of the IAnodePlane component (default: "AnodePlane") |
| `nsigma` | number of Gaussian sigma at which to truncate the diffusion (default: 3.0) |
| `start_time` | initial readout start time in WCT time units (default: 0 ns) |
| `readout_time` | duration of each readout frame in WCT time units (default: 5 ms) |
| `tick` | waveform sample period in WCT time units (default: 0.5 us) |
| `drift_speed` | nominal electron drift speed (default: 1.0 mm/us) |
| `continuous` | if true use continuous readout mode (default: true) |
| `fixed` | if true use fixed time-window readout mode (default: false) |
| `first_frame_number` | starting frame counter value (default: 0) |
| `frame_tag` | tag string applied to output frames (default: "") |
| `fluctuate` | if true apply Poisson charge fluctuation (default: false) |
| `rng` | type-name of the IRandom component used when fluctuate is true (default: "Random") |
