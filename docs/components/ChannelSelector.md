# ChannelSelector

Filters frame traces to retain only those whose channel ID appears in a configured set and/or whose tag matches a configured list of tags; supports trace summary propagation and optional tag transformation rules.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `ChannelSelector` |
| Concrete class | `WireCell::SigProc::ChannelSelector` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `channels` | array of channel IDs to retain; only traces on these channels pass through |
| `tags` | array of trace tags to consider; if empty, untagged traces are used |
| `tag_rules` | array of tag-transformation rule objects; if present, output tags are derived from input tags via regex rules |
