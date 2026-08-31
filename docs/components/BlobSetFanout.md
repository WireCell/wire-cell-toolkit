# BlobSetFanout

Fans out a single IBlobSet to N identical output ports, where N is set by the multiplicity configuration parameter.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobSetFanout` |
| Concrete class | `WireCell::Img::BlobSetFanout` |
| Node category | fanout |
| Primary interface | `IBlobSetFanout` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `vector<IBlobSet>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of output ports to fan out to (must be positive; required) |
