# L1SPFilter

Applies L1-sparse-penalty signal processing to collection-plane waveforms using both raw ADC and deconvolved signal traces as inputs, identifying and refitting regions of interest to suppress shorted-wire artifacts and improve signal recovery.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `L1SPFilter` |
| Concrete class | `WireCell::SigProc::L1SPFilter` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `fields` | type:name of the IFieldResponse component (default "FieldResponse") |
| `filter` | array of doubles giving the smearing filter waveform applied to L1-fitted signals |
| `adctag` | trace tag for raw ADC waveforms (default "raw") |
| `sigtag` | trace tag for deconvolved signal waveforms (default "gauss") |
| `outtag` | trace tag for the output L1SP-processed waveforms (default "l1sp") |
| `gain` | electronics gain |
| `shaping` | electronics shaping time |
| `postgain` | post-gain factor |
| `ADC_mV` | ADC-to-mV conversion |
| `fine_time_offset` | fine time offset (default 0) |
| `coarse_time_offset` | coarse time offset |
| `dft` | type:name of IDFT component (default "FftwDFT") |
| `roi_pad` | ROI padding in ticks (default 3) |
| `raw_pad` | raw ROI padding in ticks (default 15) |
| `adc_l1_threshold` | ADC threshold for L1 trigger (default 6) |
| `adc_sum_threshold` | ADC sum threshold (default 160) |
| `l1_seg_length` | L1 segment length in ticks (default 120) |
| `l1_lambda` | L1 regularization parameter (default 5) |
| `l1_epsilon` | L1 convergence epsilon (default 0.05) |
| `l1_niteration` | maximum L1 iterations (default 100000) |
| `l1_decon_limit` | L1 deconvolution electron threshold (default 100) |
