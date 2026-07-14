# TensorPacker

Packs a fixed number of individual ITensor inputs from parallel ports into a single ITensorSet output.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TensorPacker` |
| Concrete class | `WireCell::Aux::TensorPacker` |
| Node category | fanin |
| Primary interface | `ITensorPacker` |
| Input type(s) | `vector<ITensor>` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | Number of input ports (individual ITensors) to pack into one ITensorSet; default 1. |
