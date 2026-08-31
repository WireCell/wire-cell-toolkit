# TrackDepos

Generates point depositions along one or more straight tracks specified as time-stamped rays; deposits are spaced by a fixed step size at a configurable fraction of the speed of light.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `TrackDepos` |
| Concrete class | `WireCell::Gen::TrackDepos` |
| Node category | source |
| Primary interface | `IDepoSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `step_size` | spacing between successive depositions along a track (default: 1.0 mm) |
| `clight` | track speed as a fraction of the speed of light (default: 1.0) |
| `tracks` | array of track objects each with "time", "ray" (start/end points), and "charge" keys |
| `group_time` | if positive, chunk the deposition stream into sub-streams each spanning at most this time duration (default: -1, disabled) |
