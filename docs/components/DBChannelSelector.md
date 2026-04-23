# DBChannelSelector

Selects traces from an input frame by channel using a channel list obtained from an IChannelNoiseDatabase; can select either "bad" channels or "misconfigured" channels depending on the configured type.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DBChannelSelector` |
| Concrete class | `WireCell::SigProc::DBChannelSelector` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `type` | which channel list to use: "bad" or "misconfigured" (default "misconfigured") |
| `channelDB` | type:name of an IChannelNoiseDatabase component (required) |
