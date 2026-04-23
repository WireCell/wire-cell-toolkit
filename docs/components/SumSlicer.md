# SumSlicer

Slices a frame into time slices by summing all non-zero charge samples from tagged (or all) traces within each tick-span window; outputs slices for one frame as a single ISliceFrame.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `SumSlicer` |
| Concrete class | `WireCell::Img::SumSlicer` |
| Node category | queuedout |
| Primary interface | `IFrameSlicer` |
| Input type(s) | `IFrame` |
| Output type(s) | `ISlice` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component (user must set) |
| `tick_span` | number of ticks summed into each output slice (default: 4) |
| `tag` | trace tag selecting which traces to use; empty string uses all traces |
| `slice_eos` | if true, append a null (EOS) slice at the end of slices from each input frame (default: false) |
