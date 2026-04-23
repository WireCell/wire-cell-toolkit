# ClusterFanout

Fans out a single ICluster to N identical output ports, where N is set by the multiplicity configuration parameter.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ClusterFanout` |
| Concrete class | `WireCell::Img::ClusterFanout` |
| Node category | fanout |
| Primary interface | `IClusterFanout` |
| Input type(s) | `ICluster` |
| Output type(s) | `vector<ICluster>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of output ports to fan the cluster out to (must be positive; required) |
