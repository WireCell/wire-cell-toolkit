# BlobGrouping

Takes a cluster with blob-wire-channel topology and adds measure nodes by grouping channels into connected components per slice per plane, accumulating per-plane signal into each measure.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobGrouping` |
| Concrete class | `WireCell::Img::BlobGrouping` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

This component is configurable but has no documented parameters.
