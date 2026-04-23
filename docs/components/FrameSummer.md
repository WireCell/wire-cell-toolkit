# FrameSummer

Joins two input frames by summing their traces into a single output frame; supports optional time alignment and a constant time offset applied to the second frame.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameSummer` |
| Concrete class | `WireCell::Gen::FrameSummer` |
| Node category | join |
| Primary interface | `IFrameJoiner` |
| Input type(s) | `(IFrame, IFrame)` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `align` | if true, the second frame's time is replaced with the first frame's time before summing (default: false) |
| `offset` | constant time offset added to the second frame's time (default: 0.0) |
