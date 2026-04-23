# TaggedFrameTensorSet

Converts an IFrame to an ITensorSet by extracting per-tag waveform, channel-index, and optional summary tensors for each specified trace tag.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TaggedFrameTensorSet` |
| Concrete class | `WireCell::Aux::TaggedFrameTensorSet` |
| Node category | function |
| Primary interface | `IFrameTensorSet` |
| Input type(s) | `IFrame` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `tensors` | Array of objects each with a "tag" key and optional "pad" key, specifying which tagged trace groups to convert to tensors. |
