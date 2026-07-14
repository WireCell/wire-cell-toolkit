# OmnibusNoiseFilter

Applies a comprehensive noise filtering pipeline to frame traces: runs per-channel filters, coherent noise subtraction via grouped multi-channel filters, and per-channel status filters, propagating and updating channel mask maps throughout.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `OmnibusNoiseFilter` |
| Concrete class | `WireCell::SigProc::OmnibusNoiseFilter` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `nticks` | number of time samples to use per waveform; 0 means infer from data (default 0) |
| `maskmap` | object mapping CMM key names to output CMM key names (default {"chirp":"bad","noisy":"bad"}) |
| `channel_filters` | array of type:name strings for per-channel IChannelFilter components |
| `channel_status_filters` | array of type:name strings for per-channel status IChannelFilter components |
| `grouped_filters` | array of type:name strings for coherent-group IChannelFilter components |
| `noisedb` | type:name of an IChannelNoiseDatabase component (required) |
| `intraces` | input trace tag (default "orig") |
| `outtraces` | output trace tag (default "raw") |
