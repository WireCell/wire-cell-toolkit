# ShadowGhosting

Computes blob and cluster shadow graphs from the input cluster as a diagnostic/development component; currently passes the input cluster through unchanged without modifying any blobs.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ShadowGhosting` |
| Concrete class | `WireCell::Img::ShadowGhosting` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `shadow_type` | single-character string specifying the shadow projection type, either "w" (wire) or "c" (channel) |
