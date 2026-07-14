# Digitizer

Converts a voltage-level waveform frame into an ADC-count frame by applying a gain, per-plane baselines, and a linear ADC digitization with configurable resolution and full-scale range.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Digitizer` |
| Concrete class | `WireCell::Gen::Digitizer` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type-name of the IAnodePlane component (default: "AnodePlane") |
| `resolution` | ADC bit resolution, giving 2^resolution possible values |
| `gain` | voltage-to-ADC gain factor applied before digitization |
| `fullscale` | 2-element array [min_voltage, max_voltage] defining the ADC full-scale range |
| `baselines` | 3-element array of per-plane baseline voltages indexed by wire plane |
| `round` | if true round to nearest ADC count; if false truncate toward zero |
| `frame_tag` | tag string applied to the output frame (default: "") |
