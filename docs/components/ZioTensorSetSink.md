# ZioTensorSetSink

Receives ITensorSet objects and forwards them over a ZIO/ZeroMQ flow-extract port to a remote server, handling end-of-stream and timeout signalling.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ZioTensorSetSink` |
| Concrete class | `WireCell::Zio::TensorSetSink` |
| Node category | sink |
| Primary interface | `ITensorSetSink` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `(none)` |
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
| `attributes` | object of key/value pairs placed in the BOT flow label |
