# DNNROIFinding

Applies a deep-neural-network ROI finder to deconvolved wire-plane traces by running them through a TorchScript model and produces a new frame with signal-processed charge traces.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DNNROIFinding` |
| Concrete class | `WireCell::Pytorch::DNNROIFinding` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane that defines the channel range (default "AnodePlane") |
| `plane` | wire-plane index (0=U, 1=V, 2=W) selecting which channels to process (default 0) |
| `forward` | type/name of the ITensorForward service running the model (default "TorchService") |
| `intags` | array of input trace tags fed as network input channels (required) |
| `decon_charge_tag` | trace tag for the deconvolved charge used as the signal source after ROI masking (required) |
| `outtag` | output trace tag applied to the resulting signal traces |
| `input_scale` | scale factor applied to input charge values before the model (default 1/4000) |
| `output_scale` | scale factor applied to output charge after ROI masking (default 1.0) |
| `tick0` | first tick index of the data window (default 0) |
| `nticks` | number of ticks in the data window (default 6000) |
| `mask_thresh` | probability threshold for converting the network ROI mask to binary (default 0.5) |
