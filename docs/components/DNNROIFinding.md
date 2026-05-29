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
| `input_scale` | scale factor applied uniformly to all input charge values before the model (default 1/4000); set to 1.0 when the TorchScript model bakes per-channel normalization internally (e.g. 6-channel models) |
| `input_offset` | constant offset added to input charge values after scaling (default 0.0) |
| `output_scale` | scale factor applied to output charge after ROI masking (default 1.0) |
| `tick0` | first tick index of the data window (default 0) |
| `nticks` | number of ticks in the data window (default 6000) |
| `mask_thresh` | probability threshold for converting the network ROI mask to binary (default 0.5) |
| `tick_per_slice` | model time-axis rebin factor; raw ticks are downsampled by this stride before inference and upsampled after (default 10; current PDHD/PDVD configs use 4) |
| `tick_pad_multiple` | round model input up to a multiple of this value before inference (default 0 → fall back to `tick_per_slice`). Set larger (e.g. 128) for models whose post-rebin width must be divisible by 2^N for some N stride-2 down/up cycles |

## Tick padding

The node rounds the input tick count up to a multiple of `pad_mult` before
inference, then crops the output back to the original `input_ticks` length:

```
pad_mult     = tick_pad_multiple > 0 ? tick_pad_multiple : tick_per_slice
model_ticks  = ceil(input_ticks / pad_mult) * pad_mult
output       = inference(...).leftCols(input_ticks)   // crop back
```

So padding is invisible to downstream nodes — it only zero-fills the tail of
the model's input to satisfy the network's alignment constraint.

How to pick `tick_pad_multiple` for a given `.ts`:

| Model topology (in the tick axis, post-rebin) | Required `tick_pad_multiple` |
|---|---|
| Plain UNet with no internal stride-2 cascade (PDHD MobileNetV3-large UNet) | unset / `tick_per_slice` |
| UNet with N stride-2 down+up levels post-rebin | `tick_per_slice · 2^N` |
| PDVD MobileNetV3-UNet (5 stride-2 levels, `tick_per_slice=4`) | **128** (= 4·32) |

Mismatch surfaces as a runtime tensor-shape error inside the model, not a
toolkit-side check, so set this knob from the model card.

## Input Normalization

Before being passed to the model, all input channels are transformed uniformly:

```
model_input = raw_charge * input_scale + input_offset
```

The location of normalization depends on the model:

- **3-channel models** — normalization is external: set `input_scale = 1/4000` (the default) so C++ scales the charge before inference.
- **6-channel models** — per-channel z-score normalization is baked into the TorchScript graph itself, so `input_scale` should be `1.0` (pass-through) and the model handles the rest internally.

