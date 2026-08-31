# HDF5FrameTap

Passes an IFrame through unchanged while writing the selected tagged traces to an HDF5 file, storing per-tag 2-D sample arrays, channel arrays, and tick-info arrays under a group named by the frame ident.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `HDF5FrameTap` |
| Concrete class | `WireCell::Hio::HDF5FrameTap` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path of the HDF5 output file (default "wct-frame.hdf5") |
| `trace_tags` | array of trace tags to save; if empty, all traces are saved |
| `digitize` | if true, save samples as 16-bit integers; otherwise as 32-bit floats (default false) |
| `baseline` | value added to the sample array before scaling (default 0.0) |
| `scale` | multiplier applied to each sample after baseline addition (default 1.0) |
| `offset` | additive offset applied after scaling (default 0.0) |
| `gzip` | gzip compression level 0-9; 0 disables compression (default 0) |
| `shuffle` | if true, apply HDF5 byte-shuffle filter before compression (default false) |
| `chunk` | integer or two-element array setting the HDF5 chunk dimensions |
| `tick0` | override for the first tick index |
| `nticks` | override for the number of ticks per channel |
