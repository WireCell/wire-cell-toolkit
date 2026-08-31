# GridTiling

Tiles a time slice into blobs by projecting channel activity onto a ray-grid coordinate system for a single anode face, producing a set of geometrically consistent 3D blobs.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `GridTiling` |
| Concrete class | `WireCell::Img::GridTiling` |
| Node category | function |
| Primary interface | `ITiling` |
| Input type(s) | `ISlice` |
| Output type(s) | `IBlobSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component (user must set) |
| `face` | integer index of the anode face to tile (default: 0) |
| `threshold` | minimum activity value required for a channel to be considered active |
| `nudge` | small nudge value applied during blob boundary construction |
