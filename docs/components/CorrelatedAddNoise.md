# CorrelatedAddNoise

Adds frequency-banded inter-channel correlated electronics noise to a frame; draws correlated noise per frequency bin across wires via a coloring matrix loaded from a JSON model file.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `CorrelatedAddNoise` |
| Concrete class | `WireCell::Gen::CorrelatedAddNoise` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `model_file` | path to the JSON (or .bz2-compressed JSON) correlated noise model file |
| `rng` | type:name of the IRandom service (default: "Random:default") |
| `nsamples` | number of time samples per waveform (default: 2128) |
| `dt` | sample period in WCT time units; values less than 1.0 are interpreted as seconds (default: 0.5 us) |
| `ifft_scale` | multiplicative scale applied after the IFFT normalization (default: 1.0) |
