# ClusterScopeFilter

Filters a cluster graph by removing all blob nodes that do not belong to the configured anode face index, passing through all other node types (wires, channels, slices, measures) unchanged.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ClusterScopeFilter` |
| Concrete class | `WireCell::Img::ClusterScopeFilter` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `face_index` | index of the anode face whose blobs are retained; blobs on other faces are removed (default: -1, meaning no filtering) |
