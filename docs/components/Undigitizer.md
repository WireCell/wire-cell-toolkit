# Undigitizer

Inverts the ADC digitization transform to convert integer ADC sample values back to analog voltage values, reversing the gain, baseline subtraction, and full-scale mapping applied by the Digitizer component.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `Undigitizer` |
| Concrete class | `WireCell::SigProc::Undigitizer` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type:name of the IAnodePlane component (required) |
| `resolution` | ADC bit resolution (required) |
| `gain` | linear gain factor of the electronics (required) |
| `fullscale` | two-element array [Vmin, Vmax] giving the full-scale voltage range (required) |
| `baselines` | array of per-plane baseline voltages (required) |
