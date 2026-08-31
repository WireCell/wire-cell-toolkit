# OmnibusSigProc

Performs comprehensive 2D deconvolution signal processing on all wire planes: applies field and electronics response deconvolution, ROI formation and refinement with multiple filter stages, multi-plane protection, and outputs tagged deconvolved charge and signal waveforms.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `OmnibusSigProc` |
| Concrete class | `WireCell::SigProc::OmnibusSigProc` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component (required) |
| `dft` | type:name of the IDFT component (default "FftwDFT") |
| `field_response` | type:name of IFieldResponse component |
| `per_chan_resp` | type:name of IChannelResponse for per-channel corrections |
| `gain` | electronics gain |
| `shaping` | electronics shaping time |
| `sparse` | if true use sparse frame output (default false) |
| `process_planes` | array of plane indices to process (default all) |
| `wiener_tag` | output tag for Wiener-filtered traces |
| `gauss_tag` | output tag for Gaussian-filtered traces |
| `decon_charge_tag` | output tag for deconvolved charge traces |
| `frame_tag` | output frame-level tag |
