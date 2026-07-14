# DepoSetScaler

Adapts an IDrifter-based scaler component to operate on a full set of depositions at once, applying per-depo charge or position scaling to every deposition in the set.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoSetScaler` |
| Concrete class | `WireCell::Gen::DepoSetScaler` |
| Node category | function |
| Primary interface | `IDepoSetFilter` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `scaler` | type-name of the IDrifter-based scaler component to use (default: "Scaler") |
