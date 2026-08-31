# Resampler

Resamples all waveforms in an IFrame from the input sampling period to a target period using FFT-based rational resampling via the LMN library.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Resampler` |
| Concrete class | `WireCell::Aux::Resampler` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `period` | Target output sampling period (in WCT time units) to resample waveforms to; required. |
| `dft` | Type-name of the IDFT service to use; default "FftwDFT". |
| `time_pad` | Padding strategy before resampling; one of "zero", "first", "last", or "linear"; default "zero". |
| `time_sizing` | Output waveform length strategy; "duration" or "padded"; default "padded". |
