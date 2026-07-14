# NoiseSource

Generates uncorrelated electronics noise frames by producing random waveforms for every channel of an anode plane using a per-channel frequency-domain noise spectrum model.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `NoiseSource` |
| Concrete class | `WireCell::Gen::NoiseSource` |
| Node category | source |
| Primary interface | `IFrameSource` |
| Input type(s) | `(none)` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component |
| `model` | type:name of the IChannelSpectrum component providing per-channel noise spectra |
| `rng` | type:name of the IRandom component (default: "Random") |
| `dft` | type:name of the IDFT component (default: "FftwDFT") |
| `start_time` | simulation start time (default: 0.0 ns) |
| `stop_time` | simulation stop time (default: 1.0 ms) |
| `readout_time` | duration of each output frame (default: 5.0 ms) |
| `sample_period` | ADC sample period (default: 0.5 us) |
| `first_frame_number` | starting frame counter value (default: 0) |
| `nsamples` | number of time samples per trace (default: 9600) |
| `replacement_percentage` | fraction of random numbers replaced each frame (default: 0.02) |
