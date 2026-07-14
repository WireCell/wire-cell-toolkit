# SilentNoise

Generates frames of zero-valued (silent) noise traces, useful as a placeholder or for testing; produces a configurable number of output frames each containing a fixed number of all-zero traces.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `SilentNoise` |
| Concrete class | `WireCell::Gen::SilentNoise` |
| Node category | source |
| Primary interface | `IFrameSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `noutputs` | number of frames to produce; 0 means run forever (default: 0) |
| `nchannels` | number of zero-valued traces per frame (default: 0) |
| `traces_tag` | tag string applied to the trace set in each output frame (default: "") |
