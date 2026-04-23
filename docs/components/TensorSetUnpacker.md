# TensorSetUnpacker

Unpacks individual ITensors from an ITensorSet onto separate output ports, selecting each tensor by a (tag, type) pair.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TensorSetUnpacker` |
| Concrete class | `WireCell::Aux::TensorSetUnpacker` |
| Node category | fanout |
| Primary interface | `ITensorSetUnpacker` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `vector<ITensor>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `tags` | Array of tag strings (one per output port) used to look up each tensor; length sets output multiplicity. |
| `types` | Array of type strings (one per output port) paired with tags to identify tensors. |
| `multiplicity` | Number of output ports; set implicitly by the length of the tags array. |
