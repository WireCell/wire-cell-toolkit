# AddNoise

Historical alias for IncoherentAddNoise; adds independent per-channel noise to every trace in a frame using IChannelSpectrum noise models.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `AddNoise` |
| Concrete class | `WireCell::Gen::IncoherentAddNoise` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `model` | type:name of a single IChannelSpectrum noise model, or an array of such names |
| `models` | array of IChannelSpectrum type:name strings |
| `rng` | type:name of the IRandom service (default: "Random") |
| `dft` | type:name of the IDFT service (default: "FftwDFT") |
| `nsamples` | number of time samples expected per noise waveform (default: 9600) |
| `replacement_percentage` | fraction of the recycled random buffer refreshed each call (default: 0.02) |
