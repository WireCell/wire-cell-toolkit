# Drifter

Drifts depositions to a response plane along the X axis, applying longitudinal and transverse diffusion, electron lifetime absorption, optional charge fluctuation, and a time offset.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Drifter` |
| Concrete class | `WireCell::Gen::Drifter` |
| Node category | queuedout |
| Primary interface | `IDrifter` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `xregions` | JSON array of X-region objects with "anode", "response", and "cathode" sub-objects (required) |
| `DL` | longitudinal diffusion coefficient in cm²/s (default: 7.2) |
| `DT` | transverse diffusion coefficient in cm²/s (default: 12.0) |
| `lifetime` | electron absorption lifetime in WCT time units (default: 8 ms) |
| `fluctuate` | if true fluctuate the number of absorbed electrons binomially (default: true) |
| `drift_speed` | electron drift speed in WCT velocity units (default: 1.6 mm/us) |
| `time_offset` | constant time offset added to drifted deposition times (default: 0) |
| `charge_scale` | multiplicative scale factor applied to the drifted charge (default: 1.0) |
| `rng` | type-name of the IRandom component used for fluctuation (default: "Random") |
