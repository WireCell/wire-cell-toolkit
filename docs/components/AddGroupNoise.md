# AddGroupNoise

Adds grouped correlated noise to a frame by assigning channels to noise groups and generating one shared noise waveform per group via inverse FFT.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `AddGroupNoise` |
| Concrete class | `WireCell::Gen::AddGroupNoise` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `spectra_file` | path to JSON file with per-group noise amplitude spectra |
| `map_file` | path to JSON file mapping channels to groups |
| `rng` | type:name of the IRandom service (default: "Random") |
| `dft` | type:name of the IDFT service (default: "FftwDFT") |
| `nsamples` | number of time samples in the output noise waveform (default: 4096) |
| `spec_scale` | multiplicative scale factor applied to spectral amplitudes (default: 1.0) |
