# DepoSetFilterYZ

Filters depositions in a set by checking their Y-Z position against a 2D acceptance map loaded from a file, retaining only depositions that match a specified response bin value.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoSetFilterYZ` |
| Concrete class | `WireCell::Gen::DepoSetFilterYZ` |
| Node category | function |
| Primary interface | `IDepoSetFilter` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `yzmap_filename` | path to the JSON file containing the Y-Z acceptance map |
| `bin_width` | bin width in the Z direction (in WCT length units) |
| `tpc_width` | width of the TPC in the X direction |
| `bin_height` | bin height in the Y direction |
| `yoffset` | Y offset applied when computing the Y bin index |
| `zoffset` | Z offset applied when computing the Z bin index |
| `nbinsy` | number of bins along Y |
| `nbinsz` | number of bins along Z |
| `resp` | response bin value to select |
| `anode` | type-name of the IAnodePlane component (default: "AnodePlane") |
| `plane` | wire plane index used to look up the map entry |
