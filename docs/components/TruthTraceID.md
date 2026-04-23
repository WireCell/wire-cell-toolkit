# TruthTraceID

Produces truth-level frames by diffusing depositions and optionally applying configurable high-frequency wire and time filter responses; supports "Bare" (no filter) or filtered truth types for signal-processing studies.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TruthTraceID` |
| Concrete class | `WireCell::Gen::TruthTraceID` |
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
| `truth_type` | "Bare" for no filter, or any other string to enable HF filter convolution (default: "Bare") |
| `truth_gain` | multiplicative gain on output charge (default: -1.0) |
| `fluctuate` | if true apply fluctuations (default: true) |
| `first_frame_number` | starting frame counter (default: 0) |
