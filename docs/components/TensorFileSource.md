# TensorFileSource

Reads ITensorSet objects from a streaming archive file by parsing paired metadata JSON files and NumPy array files per tensor, grouped by tensor set ident, and reconstructing them into ITensorSet objects.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TensorFileSource` |
| Concrete class | `WireCell::Sio::TensorFileSource` |
| Node category | source |
| Primary interface | `ITensorSetSource` |
| Input type(s) | `(none)` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `inname` | input archive file path (required) |
| `prefix` | expected string prefix of tensor entry names within the archive |
