# TensorSetFanin

Merges multiple ITensorSet inputs into a single ITensorSet by concatenating tensors from selected input ports in a configurable order.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TensorSetFanin` |
| Concrete class | `WireCell::Aux::TensorSetFanin` |
| Node category | fanin |
| Primary interface | `ITensorSetFanin` |
| Input type(s) | `vector<ITensorSet>` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | Number of input ITensorSet ports; required. |
| `ident_port` | Index of the input port whose ident is used for the output; default 0. |
| `tensor_order` | Array of port indices specifying concatenation order; defaults to all ports in order. |
