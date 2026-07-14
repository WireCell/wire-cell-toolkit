# FrameFanout

Fans a single input frame out to multiple output ports; in trivial mode the same frame is sent to all ports; in rule-driven mode per-port tag transformation rules determine the frame and trace tags on each output.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameFanout` |
| Concrete class | `WireCell::Gen::FrameFanout` |
| Node category | fanout |
| Primary interface | `IFrameFanout` |
| Input type(s) | `IFrame` |
| Output type(s) | `vector<IFrame>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `multiplicity` | number of output ports (must be positive) |
| `tag_rules` | array of per-output-port rule objects with "frame" and/or "trace" keys mapping regex patterns to output tags; empty selects trivial fanout |
