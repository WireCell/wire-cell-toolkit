# ZioTorchScript

Forwards an input tensor set to a remote TorchScript inference service over ZIO/ZeroMQ using a Majordomo client pattern and returns the resulting tensor set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ZioTorchScript` |
| Concrete class | `WireCell::Pytorch::ZioTorchScript` |
| Node category | function |
| Primary interface | `ITensorSetFilter` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `address` | ZeroMQ endpoint address of the broker or inference server (default "tcp://localhost:5555") |
| `service` | Majordomo service name to request (default "torch:dnnroi") |
| `model` | model identifier passed in configuration (default "model.ts") |
| `gpu` | hint to the remote service to run on GPU (default true) |
| `wait_time` | milliseconds to wait between retries on failure (default 500) |
