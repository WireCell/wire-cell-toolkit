# WireBoundedDepos

Filters depositions by whether they project onto specified wire-number ranges in given planes; can operate in accept or reject mode, supporting multiple wire regions.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `WireBoundedDepos` |
| Concrete class | `WireCell::Gen::WireBoundedDepos` |
| Node category | queuedout |
| Primary interface | `IDrifter` |
| Input type(s) | `IDepo` |
| Output type(s) | `IDepo` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component |
| `regions` | array of region objects; each region is an array of constraint objects with "plane", "min" (minimum wire number), and "max" (maximum wire number) keys |
| `mode` | "accept" to pass depositions matching a region, any other value to reject them (default: "accept") |
