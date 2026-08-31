# Diffuser

Models longitudinal and transverse diffusion of a drifted deposition, projecting charge onto a wire pitch grid and producing time-ordered IDiffusion patch objects with Gaussian spreading in both drift-time and pitch directions.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Diffuser` |
| Concrete class | `WireCell::Diffuser` |
| Node category | queuedout |
| Primary interface | `IDiffuser` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDiffusion` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `pitch_origin` | 3D point giving the origin of the wire pitch coordinate system |
| `pitch_direction` | unit vector giving the wire pitch direction |
| `pitch_distance` | wire pitch spacing in WCT length units (default: 5 mm) |
| `timeslice` | longitudinal bin size (time slice width) in WCT time units (default: 2 us) |
| `timeoffset` | time offset added to each deposition time (default: 0) |
| `starttime` | longitudinal grid origin in WCT time units (default: 0) |
| `DL` | longitudinal diffusion coefficient in cm²/s (default: 5.3) |
| `DT` | transverse diffusion coefficient in cm²/s (default: 12.8) |
| `drift_velocity` | electron drift velocity (default: 1.6 mm/us) |
| `max_sigma_l` | maximum longitudinal sigma for output ordering buffer (default: 5 us) |
| `nsigma` | number of Gaussian sigma at which to truncate the diffusion patch (default: 3.0) |
