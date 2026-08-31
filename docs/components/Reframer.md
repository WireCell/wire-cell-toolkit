# Reframer

Re-samples input frame traces onto a fixed rectangular grid of channels (from an anode plane) and a fixed time window, accumulating overlapping trace samples; supports selection by trace tags and optional frame-level output tag.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Reframer` |
| Concrete class | `WireCell::Gen::Reframer` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane providing the channel list |
| `tags` | array of trace tag strings to select from the input frame; if empty all traces are used |
| `frame_tag` | tag string to apply to the output frame (default: "") |
| `tbin` | starting time bin index of the output window (default: 0) |
| `nticks` | number of ticks in the output window |
| `toffset` | additional time offset added to the output frame's time (default: 0.0) |
| `fill` | fill value for output samples with no input coverage (default: 0.0) |
