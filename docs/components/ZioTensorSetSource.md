# ZioTensorSetSource

Receives ITensorSet objects from a remote producer over a ZIO/ZeroMQ flow-inject port and injects them into the data-flow graph, treating an empty payload or EOT as end-of-stream.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ZioTensorSetSource` |
| Concrete class | `WireCell::Zio::TensorSetSource` |
| Node category | source |
| Primary interface | `ITensorSetSource` |
| Input type(s) | `(none)` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `nodename` | ZIO node nickname used for Zyre discovery |
| `portname` | name of the ZIO flow port (default "flow") |
| `timeout` | ZeroMQ send/receive timeout in milliseconds (default 1000) |
| `credit` | ZIO flow credit (window size) for back-pressure (default 10) |
| `binds` | array of addresses the port should bind to |
| `connects` | array of addresses the port should connect to |
| `headers` | object of additional Zyre header key/value pairs to advertise |
