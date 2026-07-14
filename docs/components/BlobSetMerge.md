# BlobSetMerge

Merges N synchronised IBlobSets (one per input port) into a single IBlobSet by combining their blobs and merging the underlying slice activity; input slices must be time-synchronised.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobSetMerge` |
| Concrete class | `WireCell::Img::BlobSetMerge` |
| Node category | fanin |
| Primary interface | `IBlobSetFanin` |
| Input type(s) | `vector<IBlobSet>` |
| Output type(s) | `IBlobSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of input ports to merge (must be positive; required) |
