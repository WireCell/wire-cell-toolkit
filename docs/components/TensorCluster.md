# TensorCluster

Converts an ITensorSet to an ICluster graph by locating the cluster tensor at a configurable datapath regex and reconstructing the cluster using provided anode planes.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TensorCluster` |
| Concrete class | `WireCell::Aux::TensorCluster` |
| Node category | function |
| Primary interface | `ITensorSetCluster` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `ICluster` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `dpre` | Regular expression matched against tensor datapaths to find the cluster root tensor; default "clusters/[[:digit:]]+$". |
| `anodes` | Array of IAnodePlane type-name strings required to reconstruct the cluster geometry; default []. |
