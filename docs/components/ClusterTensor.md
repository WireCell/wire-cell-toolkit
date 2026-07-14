# ClusterTensor

Converts an ICluster graph to an ITensorSet using the tensor data model representation of clusters.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ClusterTensor` |
| Concrete class | `WireCell::Aux::ClusterTensor` |
| Node category | function |
| Primary interface | `IClusterTensorSet` |
| Input type(s) | `ICluster` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `datapath` | Tensor data model path (with optional %d format for cluster ident) at which cluster tensors are stored; default "". |
