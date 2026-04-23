# DepoSetRotate

Applies an axis permutation and per-axis scaling to the positions of all depositions in a set, enabling coordinate rotations or reflections of the deposition cloud.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DepoSetRotate` |
| Concrete class | `WireCell::Gen::DepoSetRotate` |
| Node category | function |
| Primary interface | `IDepoSetFilter` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IDepoSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `rotate` | boolean flag to enable the rotation/transpose transform (default: false) |
| `transpose` | 3-element integer array specifying the axis permutation [ix, iy, iz] |
| `scale` | 3-element float array specifying per-axis scale factors applied after the transpose |
