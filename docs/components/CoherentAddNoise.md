# CoherentAddNoise

Adds coherent (per-group shared) noise to a frame using IGroupSpectrum noise models; all channels in the same group receive an identical noise waveform drawn once per group per frame.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `CoherentAddNoise` |
| Concrete class | `WireCell::Gen::CoherentAddNoise` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `model` | type:name of a single IGroupSpectrum noise model, or an array of such names |
| `models` | array of IGroupSpectrum type:name strings (alternative to "model") |
| `rng` | type:name of the IRandom service (default: "Random") |
| `dft` | type:name of the IDFT service (default: "FftwDFT") |
| `nsamples` | number of time samples expected per noise waveform (default: 9600) |
| `replacement_percentage` | fraction of the recycled random buffer refreshed each call (default: 0.02) |
