# FrameSummerYZ

Sums a fixed number of input frames (e.g. from Y and Z drift volumes) into a single output frame by accumulating all traces; multiplicity must be set to the expected number of inputs.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameSummerYZ` |
| Concrete class | `WireCell::Gen::FrameSummerYZ` |
| Node category | fanin |
| Primary interface | `IFrameFanin` |
| Input type(s) | `vector<IFrame>` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of input frames to sum (must be positive) |
