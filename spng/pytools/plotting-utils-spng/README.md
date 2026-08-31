osp-spng-frame-visualizers
==========================

YAML-configured frame visualization for comparing OSP/SPNG-style NPZ and
HDF5 frame data.

The project installs a `frame-visualizers` command that loads one YAML file,
opens up to three active frame datasets at a time, and provides an interactive
Matplotlib UI for image views plus projected 1D traces.

Architecture
------------

- `frame_visualizers.py`: command-line entry point.
- `frame_config.py`: YAML loading, option dumping, and validation.
- `frame_models.py`: shared dataclasses and view constants.
- `frame_data_loaders.py`: NPZ/HDF5 loading and 1D-to-2D reshape handling.
- `frame_collection.py`: active-frame selection and view slicing.
- `frame_visualizer_ui.py`: Matplotlib UI, image panels, controls, selection,
  zoom, and projected 1D trace plotting.

The code is intentionally split so configuration parsing, data loading, data
selection, and UI rendering can change independently.

Installation
------------

From this repository:

```bash
python3 -m pip install -e .
```

This installs the package in editable mode and exposes:

```bash
frame-visualizers
```

Python 3.10 or newer is required. Runtime dependencies are declared in
`pyproject.toml`:

- `h5py`
- `matplotlib`
- `numpy`
- `PyYAML`

Usage
-----

Run the visualizer with a YAML configuration:

```bash
frame-visualizers 3unets-frame-visualizers.yaml
```

Select a specific frame ordinal when a dataset pattern or repeated HDF5 basename
matches multiple frames:

```bash
frame-visualizers 3unets-frame-visualizers.yaml --frame-index 0
```

Change the half-window around the selected tick in the projected 1D trace:

```bash
frame-visualizers 3unets-frame-visualizers.yaml --trace-half-window 200
```

Print all YAML options supported by the current installed code:

```bash
frame-visualizers --dump-yaml-options
```

YAML Configuration
------------------

The YAML file must contain a top-level `files` list. Each entry describes one
input frame source.

```yaml
files:
  - label: spng
    name: 3UNets/spng-sigproc-cosmics_0_wct-depos.npz
    dset: "frame_*_1"
    format: npz
    xaxis_label: tick
    yaxis_label: channel within view
    trace_yaxis_label: ADC values
    channel_mapping:
      U: [0, 800]
      V: [800, 1600]
      W: [1600, 2560]
      tick_info: 6000
```

Supported file fields:

- `label`: display name. Defaults to the input file stem.
- `name`: input file path. Relative paths resolve from the YAML file location.
- `format`: `npz`, `h5`, or `hdf5`.
- `dset`: NPZ key, HDF5 dataset path, or HDF5 dataset basename. Required for
  HDF5. For NPZ, omitting it selects sorted `frame_*` keys by `--frame-index`.
- `xaxis_label`: x-axis label for the 2D image panel.
- `yaxis_label`: y-axis label for the 2D image panel.
- `trace_yaxis_label`: left y-axis label for normal projected 1D traces.
- `channel_mapping`: view-to-channel mapping for `U`, `V`, and `W`.
- `channel_mapping.tick_info`: tick count used to reshape 1D arrays into
  channel-by-tick frames.
- `separate_axis`: optional mapping for drawing that file's projected 1D trace
  on its own right/top axes.

Separate Axis Traces
--------------------

Use `separate_axis` when one trace has a different scale and needs its own axis.

```yaml
files:
  - label: drifted
    name: 3UNets/drifted-cosmics_0_wct-depos.npz
    dset: "frame_*_1"
    format: npz
    trace_yaxis_label: collected charge
    channel_mapping:
      U: [0, 800]
      V: [800, 1600]
      W: [1600, 2560]
      tick_info: 6000
    separate_axis:
      enabled: true
      xaxis_label: tick number
      yaxis_label: ADC values
```

`separate_axis: true` is still accepted. The mapping form is preferred because it
can name the separate top x-axis and right y-axis.

Interactive UI
--------------

The UI shows:

- Three image panels for the active files.
- U/V/W view buttons.
- A reset button for zoom.
- An absolute-channel text box.
- Active-file checkboxes with a maximum of three active files.
- A projected 1D trace panel for the selected channel around the selected tick.

Click a bin in an image panel to update the projected trace. Drag on an image
panel to zoom into a tick/channel region. Use the channel box to jump to an
absolute channel.
