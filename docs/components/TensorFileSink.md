# TensorFileSink

Writes each ITensorSet to a streaming archive file as a collection of NumPy array files (tensor data) and JSON files (tensor set and individual tensor metadata), with an optional dump mode that discards data for debugging.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TensorFileSink` |
| Concrete class | `WireCell::Sio::TensorFileSink` |
| Node category | sink |
| Primary interface | `ITensorSetSink` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `outname` | output archive file path; extension determines compression (required) |
| `prefix` | string prefix prepended to all entry names within the archive |
| `dump_mode` | if true, discard all data without writing (default false) |
