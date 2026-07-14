# UncorrelatedAddNoise

Adds per-channel uncorrelated electronics noise to each trace in a frame; noise is synthesized in frequency space using independent complex Gaussian coefficients scaled to match a per-wire mean-magnitude spectrum from a JSON model file.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `UncorrelatedAddNoise` |
| Concrete class | `WireCell::Gen::UncorrelatedAddNoise` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `rng` | type:name of the IRandom component (default: "Random:default") |
| `model_file` | path to the JSON or JSON.bz2 noise model file (default: "uncorrelated_noise_model.json.bz2") |
| `nsamples` | number of time samples expected per trace (default: 2128) |
| `dt` | ADC sample period; values less than 1.0 are interpreted as seconds (default: 0.5 us) |
| `ifft_scale` | additional scale factor applied after the IFFT normalization (default: 1.0) |
