# DepoSetFilter

Filters depositions in a set to retain only those whose position falls within the sensitive volume of any face of a given anode plane.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoSetFilter` |
| Concrete class | `WireCell::Gen::DepoSetFilter` |
| Node category | function |
| Primary interface | `IDepoSetFilter` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type-name of the IAnodePlane component used to define the sensitive volume (default: "AnodePlane") |
