# Retagger

Transforms frame and trace tags according to configurable regex-based rules, optionally merging tagged trace sets under new combined tags and remapping channel mask map entries by name.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Retagger` |
| Concrete class | `WireCell::Gen::Retagger` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `tag_rules` | array with a single element (index 0) containing an object with "frame", "trace", and/or "merge" keys; each maps regex patterns to output tag strings or arrays |
| `maskmap` | object mapping input channel mask map names to output names (default: {"bad": "bad"}) |
