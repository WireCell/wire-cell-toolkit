# GlobalGeomClustering

Re-clusters blobs globally across the entire detector volume by removing existing blob-blob edges and creating new ones using a geometric clustering policy applied without regard to existing cluster groupings.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `GlobalGeomClustering` |
| Concrete class | `WireCell::Img::GlobalGeomClustering` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `clustering_policy` | string name of the geometric clustering policy to apply when forming new blob-blob edges |
