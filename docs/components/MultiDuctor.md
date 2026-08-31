# MultiDuctor

Routes each deposition to one of several sub-ductors selected by configurable rules, merges their output frames, and produces readout-time-windowed output frames; supports continuous or data-driven readout timing.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `MultiDuctor` |
| Concrete class | `WireCell::Gen::MultiDuctor` |
| Node category | queuedout |
| Primary interface | `IDuctor` |
| Input type(s) | `IDepo` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component |
| `chain` | array of chain objects; each chain is an array of rule objects with "ductor" (type:name), "rule" ("wirebounds" or "bool"), and "args" keys |
| `tick` | sample period in us; must match sub-ductors (default: 0.5) |
| `start_time` | initial readout start time in ns (default: 0.0) |
| `readout_time` | duration of each readout window in ms (default: 5.0) |
| `continuous` | if true produce frames continuously even with no depositions (default: false) |
| `first_frame_number` | starting frame counter value (default: 0) |
