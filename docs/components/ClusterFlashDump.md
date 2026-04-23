# ClusterFlashDump

Logs a summary of a tensor-data-model PC tree holding cluster and flash information, serving as a debugging/inspection sink.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ClusterFlashDump` |
| Concrete class | `WireCell::Aux::ClusterFlashDump` |
| Node category | sink |
| Primary interface | `ITensorSetSink` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `datapath` | Tensor data model path (with optional %d format for ident) locating the point-cloud tree in the ITensorSet; required, no default. |
