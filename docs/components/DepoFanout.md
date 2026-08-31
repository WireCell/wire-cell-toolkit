# DepoFanout

Fans out each input IDepo (or EOS marker) to N identical output ports; used to send the same deposition stream to multiple downstream consumers.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoFanout` |
| Concrete class | `WireCell::Gen::DepoFanout` |
| Node category | fanout |
| Primary interface | `IDepoFanout` |
| Input type(s) | `IDepo` |
| Output type(s) | `vector<IDepo>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | integer number of output ports to fan out to (must be positive) |
