# BlobSampling

Samples blobs from a blob set using named IBlobSampler instances and outputs the resulting point cloud tree as a tensor set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `BlobSampling` |
| Concrete class | `WireCell::Clus::BlobSampling` |
| Node category | function |
| Primary interface | `IBlobSampling` |
| Input type(s) | `IBlobSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `datapath` | tensor datapath template for the output point cloud tree; "%d" is interpolated with the blob-set ident (default "pointtrees/%d") |
| `samplers` | object mapping point-cloud names to IBlobSampler type/name strings; at least one entry is required |
