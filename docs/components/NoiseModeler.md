# NoiseModeler

Accumulates noise spectra from voltage-level frames, classifies traces using an ITraceRanker, collects FFT amplitudes per channel group, and at termination writes mean half-spectra to a JSON output file for use as a noise model.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `NoiseModeler` |
| Concrete class | `WireCell::SigProc::NoiseModeler` |
| Node category | sink |
| Primary interface | `IFrameSink` |
| Input type(s) | `IFrame` |
| Output type(s) | `(none)` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `dft` | type:name of the IDFT component (default "FftwDFT") |
| `isnoise` | type:name of the ITraceRanker component that ranks traces for noise likelihood (default "NoiseRanker") |
| `threshold` | minimum noise rank to consider a trace as noise (default 0.9) |
| `groups` | array of {groupID, channels} objects, or a path to a JSON file; null for a single group |
| `nfft` | number of FFT samples; if null, set from first trace size |
| `outname` | output file path for the noise model JSON (required) |
