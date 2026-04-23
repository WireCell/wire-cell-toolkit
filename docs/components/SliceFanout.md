# SliceFanout

Fans out a single input slice to multiple output ports by broadcasting the same slice pointer to each port; also propagates end-of-stream to all outputs.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `SliceFanout` |
| Concrete class | `WireCell::Img::SliceFanout` |
| Node category | fanout |
| Primary interface | `ISliceFanout` |
| Input type(s) | `ISlice` |
| Output type(s) | `vector<ISlice>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of output ports to fan out to (must be positive) |
