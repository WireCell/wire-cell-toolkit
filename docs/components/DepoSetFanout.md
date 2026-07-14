# DepoSetFanout

Fans out a single IDepoSet to a configurable number of identical output ports, passing EOS through to all outputs.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoSetFanout` |
| Concrete class | `WireCell::Gen::DepoSetFanout` |
| Node category | fanout |
| Primary interface | `IDepoSetFanout` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `vector<IDepoSet>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of output ports to fan out to (must be positive) |
