# SumSlices

Slices a frame into individual ISlice objects by summing non-zero charge samples from tagged traces within each tick-span window, emitting each slice individually; optionally pads with an empty slice when no activity is found.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `SumSlices` |
| Concrete class | `WireCell::Img::SumSlices` |
| Node category | function |
| Primary interface | `IFrameSlices` |
| Input type(s) | `IFrame` |
| Output type(s) | `ISliceSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component (user must set) |
| `tick_span` | number of ticks summed into each output slice (default: 4) |
| `tag` | trace tag selecting which traces to use; empty string uses all traces |
| `slice_eos` | if true, append a null (EOS) slice after each frame's slices (default: false) |
| `pad_empty` | if true, emit a single empty slice when no nominal slices are produced from a frame (default: true) |
