# Misconfigure

Simulates misconfigured front-end electronics by replacing the existing electronics response on each trace with a different one; convolves each waveform with the target cold-electronics response while deconvolving the original.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Misconfigure` |
| Concrete class | `WireCell::Gen::Misconfigure` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `from.gain` | gain of the original amplifier in mV/fC (default: 14.0) |
| `from.shaping` | shaping time of the original amplifier in us (default: 2.2) |
| `to.gain` | gain of the target (misconfigured) amplifier in mV/fC (default: 4.7) |
| `to.shaping` | shaping time of the target amplifier in us (default: 1.1) |
| `nsamples` | number of samples used to represent the response functions (default: 50) |
| `tick` | sample period for the response functions in us (default: 0.5) |
| `truncate` | if true, clip extra samples added by convolution so waveform length is preserved (default: true) |
| `dft` | type:name of the IDFT component (default: "FftwDFT") |
