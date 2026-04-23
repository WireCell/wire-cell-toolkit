# ClusterFileSink

Writes ICluster graphs to a streaming archive file (tar, zip, etc.) in one of several formats: JSON cluster graph schema, GraphViz DOT format, or NumPy array schema.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ClusterFileSink` |
| Concrete class | `WireCell::Sio::ClusterFileSink` |
| Node category | sink |
| Primary interface | `IClusterSink` |
| Input type(s) | `ICluster` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `outname` | output archive file path; extension controls compression (required) |
| `prefix` | string prefix prepended to each entry name within the archive |
| `format` | serialization format, "json", "dot", "numpy", or "dummy" |
