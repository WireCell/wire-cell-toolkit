# FrameTensor

Converts an IFrame to an ITensorSet using the tensor data model, with configurable encoding mode and optional ADC digitization.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `FrameTensor` |
| Concrete class | `WireCell::Aux::FrameTensor` |
| Node category | function |
| Primary interface | `IFrameTensorSet` |
| Input type(s) | `IFrame` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `datapath` | Tensor data model path for the frame tensor, supports %d format for frame ident; default "frames/%d". |
| `mode` | Mapping mode for frame-to-tensor conversion; one of "unified", "tagged", or "sparse"; default "tagged". |
| `digitize` | If true, store trace waveforms as short integers instead of floats; default false. |
| `baseline` | Baseline value added to samples before encoding; default 0. |
| `scale` | Scale factor applied after baseline shift; default 1. |
| `offset` | Offset added after scaling; default 0. |
