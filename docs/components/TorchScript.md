# TorchScript

Loads a TorchScript model file and runs it directly in-process on an input tensor set to produce an output tensor set, with optional retry on transient runtime failures.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TorchScript` |
| Concrete class | `WireCell::Pytorch::TorchScript` |
| Node category | function |
| Primary interface | `ITensorSetFilter` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `model` | filesystem path to the TorchScript (.ts) model file (default "model.ts") |
| `gpu` | if true, load and run the model on a CUDA GPU; otherwise use CPU (default false) |
| `wait_time` | time in milliseconds to wait before retrying after a model execution failure |
