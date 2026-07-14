# Scaler

Scales the charge of each deposition by a spatially-varying factor read from a 2D YZ-binned map file, filtering depositions outside the anode sensitive volume; intended for applying space-charge or recombination corrections.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Scaler` |
| Concrete class | `WireCell::Gen::Scaler` |
| Node category | queuedout |
| Primary interface | `IDrifter` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `yzmap_scale_filename` | path to a JSON file containing the 2D YZ charge-scale map indexed by anode name and plane |
| `anode` | type:name of the IAnodePlane component |
| `plane` | wire plane index used to look up the scale map (default: 0) |
| `bin_width` | z-bin width for the YZ map |
| `bin_height` | y-bin height for the YZ map |
| `n_ybin` | number of bins along Y (default: 31) |
| `n_zbin` | number of bins along Z (default: 180) |
| `yoffset` | offset added to deposition Y coordinate before binning |
| `zoffset` | offset added to deposition Z coordinate before binning |
