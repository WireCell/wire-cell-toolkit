# DepoFramer

Converts an IDepoSet into a voltage-level IFrame by delegating to an IDrifter to drift depositions and an IDuctor to convolve charge with field and electronics response.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoFramer` |
| Concrete class | `WireCell::Gen::DepoFramer` |
| Node category | function |
| Primary interface | `IDepoFramer` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `Drifter` | type:name of the IDrifter component to use (default: "Drifter") |
| `Ductor` | type:name of the IDuctor component to use (default: "Ductor") |
