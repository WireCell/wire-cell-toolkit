# TaggedTensorSetFrame

Converts an ITensorSet to an IFrame by reconstructing tagged traces from per-tag waveform, channel, and optional summary tensors for each specified tag.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TaggedTensorSetFrame` |
| Concrete class | `WireCell::Aux::TaggedTensorSetFrame` |
| Node category | function |
| Primary interface | `ITensorSetFrame` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `tensors` | Array of objects each with a "tag" key identifying which tagged tensor groups to reconstruct as traces in the output frame. |
