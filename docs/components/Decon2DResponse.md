# Decon2DResponse

Performs 2D deconvolution (wire and time) of waveform data in a tensor set using field and electronics responses, with optional per-channel electronics response correction, wire and time shift corrections, outputting a complex FFT tensor ready for downstream filtering.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Decon2DResponse` |
| Concrete class | `WireCell::Sig::Decon2DResponse` |
| Node category | function |
| Primary interface | `ITensorSetFilter` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `dft` | type:name of the IDFT component (default "FftwDFT") |
| `anode` | type:name of the IAnodePlane component (required) |
| `field_response` | type:name of the IFieldResponse component (required) |
| `per_chan_resp` | type:name of an IChannelResponse for per-channel electronics correction |
| `tag` | trace tag identifying the waveform tensor (default "trace_tag") |
| `iplane` | plane index (default 0) |
| `gain` | electronics gain in mV/fC (default 14.0) |
| `shaping_time` | electronics shaping time (default 2.2 us) |
| `inter_gain` | intermediate gain factor (default 1.2) |
| `ADC_mV` | ADC-to-mV conversion factor |
| `ftoffset` | fine time offset |
| `ctoffset` | coarse time offset (default -8 us) |
