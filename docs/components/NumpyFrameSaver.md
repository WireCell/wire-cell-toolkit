# NumpyFrameSaver

Saves frame waveform data as 2D NumPy arrays appended to a NumPy zip file for each configured frame tag, with optional digitization and waveform scaling, while passing the frame through unchanged.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `NumpyFrameSaver` |
| Concrete class | `WireCell::Sio::NumpyFrameSaver` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path to the output NumPy zip file (default "wct-frame.npz") |
| `frame_tags` | array of frame/trace tags to save; if empty, all traces are saved |
| `digitize` | if true, cast samples to 16-bit integers before saving (default false) |
| `baseline` | value subtracted before filling (default 0.0) |
| `scale` | multiplicative scale applied to each sample (default 1.0) |
| `offset` | additive offset applied after scaling (default 0.0) |
