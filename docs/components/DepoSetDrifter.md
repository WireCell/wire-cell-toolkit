# DepoSetDrifter

Adapts an IDrifter component to operate on a full set of depositions at once, drifting each depo through the detector medium and returning the surviving drifted set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoSetDrifter` |
| Concrete class | `WireCell::Gen::DepoSetDrifter` |
| Node category | function |
| Primary interface | `IDepoSetFilter` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `drifter` | type-name of the IDrifter component to use (default: "Drifter") |
