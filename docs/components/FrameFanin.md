# FrameFanin

Merges multiple input frames into a single output frame, concatenating all traces; supports per-port trace tagging and tag transformation rules; multiplicity may be fixed or dynamic.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameFanin` |
| Concrete class | `WireCell::Gen::FrameFanin` |
| Node category | fanin |
| Primary interface | `IFrameFanin` |
| Input type(s) | `vector<IFrame>` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of input ports; 0 enables dynamic multiplicity (default: 0) |
| `tags` | array of tag strings, one per input port; applies that tag to all traces arriving on that port in the output frame |
| `tag_rules` | array of per-port rule objects with "frame" and/or "trace" keys mapping regex patterns to output tags |
