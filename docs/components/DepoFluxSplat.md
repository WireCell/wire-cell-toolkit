# DepoFluxSplat

Converts an IDepoSet into a true-signal IFrame by applying Gaussian smearing in longitudinal and transverse directions then binning the resulting charge flux onto wire channels.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoFluxSplat` |
| Concrete class | `WireCell::Gen::DepoFluxSplat` |
| Node category | function |
| Primary interface | `IDepoFramer` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component |
| `field_response` | type:name of an IFieldResponse used to obtain drift_speed and response_plane |
| `drift_speed` | explicit drift velocity in WCT units (overrides field_response) |
| `response_plane` | explicit response plane distance in WCT length units |
| `tick` | sample period for time binning (default: 0.5 us) |
| `window_start` | start time of the acceptance window (default: 0) |
| `window_duration` | duration of the acceptance window |
| `sparse` | bool; if true produce a sparse frame (default: false) |
| `nsigma` | number of sigma at which to truncate the depo Gaussian (default: 3.0) |
| `smear_long` | extra longitudinal smearing in ticks (default: 0.0) |
| `smear_tran` | extra transverse smearing in wire pitches (default: 0.0) |
| `process_planes` | array of plane indices to process; defaults to all three |
