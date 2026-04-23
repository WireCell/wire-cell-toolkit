# PlaneSelector

Filters an IFrame to retain only traces belonging to channels of a specified wire plane on a given anode, preserving and remapping trace tags via configurable tag rules.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `PlaneSelector` |
| Concrete class | `WireCell::Aux::PlaneSelector` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | Type-name of the IAnodePlane instance used to determine channel membership; default "AnodePlane". |
| `plane` | Wire plane index (0, 1, or 2) selecting which plane's channels to keep; default 0. |
| `tags` | Array of input trace tags to consider; default []. |
| `tag_rules` | Array of rules mapping input frame/trace tags to output tags; default []. |
