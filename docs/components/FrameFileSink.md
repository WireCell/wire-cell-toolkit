# FrameFileSink

Saves frame waveform data to a streaming archive file as NumPy arrays, with configurable tag selection, optional digitization to 16-bit integers, baseline/scale/offset transforms, and dense-frame padding.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameFileSink` |
| Concrete class | `WireCell::Sio::FrameFileSink` |
| Node category | sink |
| Primary interface | `IFrameSink` |
| Input type(s) | `IFrame` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `outname` | output archive file path; extension determines compression (required) |
| `tags` | array of trace tags to save; use ["*"] to save all traces (default all) |
| `digitize` | if true, cast waveform samples to 16-bit integers (default false) |
| `baseline` | value to initialize the waveform array before filling (default 0.0) |
| `scale` | multiplicative factor applied to each sample (default 1.0) |
| `offset` | additive offset applied after scaling (default 0.0) |
| `masks` | if true, also save channel mask maps as NumPy arrays (default true) |
| `dense` | optional object with chbeg, chend, tbbeg, tbend to force a dense channel/time-bin range |
