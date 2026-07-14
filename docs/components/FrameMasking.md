# FrameMasking

Applies a channel mask to frame traces by zeroing out charge samples in masked time ranges, producing a new frame with the same structure but masked regions set to zero.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameMasking` |
| Concrete class | `WireCell::Img::FrameMasking` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component (default: "AnodePlane") |
| `cm_tag` | name of the channel mask map entry to use for masking |
| `trace_tags` | array of trace tag strings; only traces with these tags are processed and masked |
