# CelltreeSource

Reads waveform data from a ROOT celltree file (Event/Sim TTree), loading one event entry and producing a tagged IFrame with configurable branch name to trace-tag mappings, optional per-channel thresholds, and bad-channel masks.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `CelltreeSource` |
| Concrete class | `WireCell::Root::CelltreeSource` |
| Node category | source |
| Primary interface | `IFrameSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `filename` | path to the input ROOT celltree file (required) |
| `entry` | which TTree entry (event) to load (default 0) |
| `in_branch_base_names` | array of ROOT branch base names to read |
| `out_trace_tags` | array of output trace tags corresponding to each branch base name |
| `in_branch_thresholds` | array of branch names for per-channel thresholds (empty string to skip) |
| `time_scale` | integer factor to expand each histogram time bin into ticks (default 4) |
