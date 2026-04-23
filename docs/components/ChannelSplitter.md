# ChannelSplitter

Splits an input frame into multiple output frames, one per configured anode plane, by routing each trace to the output port whose anode owns that channel; supports per-port tag transformation rules.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ChannelSplitter` |
| Concrete class | `WireCell::SigProc::ChannelSplitter` |
| Node category | fanout |
| Primary interface | `IFrameFanout` |
| Input type(s) | `IFrame` |
| Output type(s) | `vector<IFrame>` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anodes` | array of type:name strings for IAnodePlane components, one per output port (required) |
| `tag_rules` | array of per-port tag-rule objects mapping input frame/trace tags to output tags |
