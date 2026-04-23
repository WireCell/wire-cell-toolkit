# DeadLiveMerging

Merges multiple cluster graphs (e.g. dead-region and live-region clusters) into a single unified cluster graph by combining all nodes and edges from each input.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DeadLiveMerging` |
| Concrete class | `WireCell::Img::DeadLiveMerging` |
| Node category | fanin |
| Primary interface | `IClusterFanin` |
| Input type(s) | `vector<ICluster>` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of input cluster ports (default: 2, must be > 0) |
| `tags` | array of string tags, one per input port, labeling each input cluster's traces |
