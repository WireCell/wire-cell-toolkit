# WireSource

Convenience source that combines WireParams and WireGenerator to produce a wire geometry directly from geometry parameters without requiring a separate source of IWireParameters.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `WireSource` |
| Concrete class | `WireCell::WireSource` |
| Node category | source |
| Primary interface | `IWireSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IWire::vector` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `center_mm` | object with x, y, z keys giving the detector center in mm |
| `size_mm` | object with x, y, z keys giving the full detector size in mm |
| `pitch_mm` | object with u, v, w keys giving the wire pitch for each plane in mm |
| `angle_deg` | object with u, v, w keys giving the wire angle w.r.t. the +Y axis in degrees |
| `offset_mm` | object with u, v, w keys giving the pitch-direction offset from center to first wire in mm |
| `plane_mm` | object with u, v, w keys giving the X position of each wire plane in mm |
