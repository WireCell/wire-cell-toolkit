# ClusterFileSource

Reads ICluster graphs from a streaming archive file in JSON cluster graph schema or NumPy cluster array schema format, reconstructing full cluster graphs using the provided anode planes.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ClusterFileSource` |
| Concrete class | `WireCell::Sio::ClusterFileSource` |
| Node category | source |
| Primary interface | `IClusterSource` |
| Input type(s) | `(none)` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `inname` | input archive file path (required) |
| `prefix` | expected string prefix of cluster entry names within the archive |
| `anodes` | array of type:name strings for IAnodePlane components (required for numpy format) |
