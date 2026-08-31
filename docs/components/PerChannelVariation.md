# PerChannelVariation

Applies per-channel electronics response variation by replacing the nominal cold-electronics response on each trace with an individually calibrated channel response from an IChannelResponse service.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `PerChannelVariation` |
| Concrete class | `WireCell::Gen::PerChannelVariation` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `gain` | nominal amplifier gain in mV/fC for the response being replaced (default: 14.0) |
| `shaping` | nominal shaping time in us for the response being replaced (default: 2.2) |
| `nsamples` | number of response samples (default: 310) |
| `tick` | sample period in us (default: 0.5) |
| `truncate` | if true, clip waveforms back to original length after convolution (default: true) |
| `per_chan_resp` | type:name of an IChannelResponse component providing per-channel responses; if empty the frame is passed through unchanged |
| `dft` | type:name of the IDFT component (default: "FftwDFT") |
