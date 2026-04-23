# LocalGeomClustering

Re-clusters blobs locally by preserving existing cluster groupings, removing old blob-blob edges, and adding new ones using geometric clustering constrained within each original cluster group.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `LocalGeomClustering` |
| Concrete class | `WireCell::Img::LocalGeomClustering` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `dryrun` | if true, pass input cluster through unchanged (default: false) |
| `clustering_policy` | string name of the geometric clustering policy to apply within each cluster group |
