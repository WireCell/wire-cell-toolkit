# BlobSetSync

Collects N IBlobSets from parallel input ports and combines their blobs into a single output IBlobSet, using the earliest slice start time, without merging underlying slice activity.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobSetSync` |
| Concrete class | `WireCell::Img::BlobSetSync` |
| Node category | fanin |
| Primary interface | `IBlobSetFanin` |
| Input type(s) | `vector<IBlobSet>` |
| Output type(s) | `IBlobSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of input ports to synchronise (must be positive; required) |
