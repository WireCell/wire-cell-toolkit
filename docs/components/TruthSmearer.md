# TruthSmearer

Produces signal-truth frames by diffusing depositions in time and wire dimensions without full field-response convolution; applies configurable Gaussian time smearing and a discrete wire smearing filter.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TruthSmearer` |
| Concrete class | `WireCell::Gen::TruthSmearer` |
| Node category | queuedout |
| Primary interface | `IDuctor` |
| Input type(s) | `IDepo` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component (default: "AnodePlane") |
| `rng` | type:name of the IRandom component (default: "Random") |
| `dft` | type:name of the IDFT component (default: "FftwDFT") |
| `start_time` | initial readout start time (default: 0.0 ns) |
| `readout_time` | duration of each readout frame (default: 5.0 ms) |
| `tick` | ADC sample period (default: 0.5 us) |
| `drift_speed` | nominal electron drift speed (default: 1.0 mm/us) |
| `continuous` | if true produce frames continuously (default: true) |
| `nsigma` | number of Gaussian sigma to retain for diffusion truncation (default: 3.0) |
| `fluctuate` | if true apply Poisson fluctuations (default: true) |
| `time_smear` | additional Gaussian longitudinal smearing sigma in time (default: 1.4 us) |
| `wire_smear_ind` | discrete wire smearing weight for induction planes (default: 0.75) |
| `wire_smear_col` | discrete wire smearing weight for collection plane (default: 0.95) |
| `truth_gain` | multiplicative gain applied to output charge (default: -1.0) |
| `first_frame_number` | starting frame counter (default: 0) |
