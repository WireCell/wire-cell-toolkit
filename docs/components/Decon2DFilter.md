# Decon2DFilter

Applies a configurable set of software time-domain filter waveforms to a 2D (wire x time) complex FFT array in a tensor set, then inverse-FFTs in the wire dimension and applies channel mask map zeroing to produce the filtered real waveform tensor.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Decon2DFilter` |
| Concrete class | `WireCell::Sig::Decon2DFilter` |
| Node category | function |
| Primary interface | `ITensorSetFilter` |
| Input type(s) | `ITensorSet` |
| Output type(s) | `ITensorSet` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `dft` | type:name of the IDFT component to use for FFT operations (default "FftwDFT") |
| `tag` | trace tag identifying the waveform tensor within the tensor set (default "trace_tag") |
| `iplane` | plane index (0=U, 1=V, 2=W) used to select a wire filter name (default 0) |
| `filters` | array of type:name strings for IFilterWaveform components to multiply together into the time filter |
