# BlobSolving

Solves for blob charge values within a cluster using a LASSO regression per slice, weighting blobs by their connectivity to adjacent slices, and returns a new cluster with updated blob charge assignments.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobSolving` |
| Concrete class | `WireCell::Img::BlobSolving` |
| Node category | function |
| Primary interface | `IClusterFilter` |
| Input type(s) | `ICluster` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

This component is configurable but has no documented parameters.
