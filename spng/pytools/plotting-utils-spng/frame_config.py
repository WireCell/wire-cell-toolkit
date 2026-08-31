"""YAML configuration parsing for frame visualizer inputs."""

import copy
import re
from pathlib import Path

import yaml

from frame_models import DEFAULT_VIEW_RANGES, FileSpec, SeparateAxisConfig, VIEW_ORDER


YAML_OPTIONS = {
    "files": [
        {
            "label": "optional display name; defaults to the input file stem",
            "name": "required input file path; relative paths resolve from this YAML file",
            "format": "required input format: npz, h5, or hdf5",
            "dset": (
                "optional NPZ key, HDF5 dataset path, or HDF5 dataset basename; "
                "required for HDF5"
            ),
            "xaxis_label": "optional x-axis label; default: tick; alias: xaxis-label",
            "yaxis_label": (
                "optional y-axis label; default: channel within view; "
                "alias: yaxis-label"
            ),
            "trace_yaxis_label": (
                "optional left y-axis label for normal projected 1D traces; "
                "alias: trace-yaxis-label"
            ),
            "channel_mapping": {
                "U": "optional [start, stop] channel range; default: [0, 800]",
                "V": "optional [start, stop] channel range; default: [800, 1600]",
                "W": "optional [start, stop] channel range; default: [1600, 2560]",
                "tick_info": (
                    "optional positive tick count for reshaping 1D inputs; "
                    "alias: tick-info"
                ),
            },
            "separate_axis": {
                "enabled": "optional boolean; defaults to true when separate_axis is a mapping",
                "xaxis_label": (
                    "optional label for this file's separate top x-axis; "
                    "alias: xaxis-label"
                ),
                "yaxis_label": (
                    "optional label for this file's separate right y-axis; "
                    "alias: yaxis-label"
                ),
            },
        }
    ]
}


def available_yaml_options():
    """Return a description of supported YAML configuration fields."""
    return copy.deepcopy(YAML_OPTIONS)


def dump_yaml_options(stream=None):
    """Dump supported YAML configuration fields as YAML.

    If *stream* is omitted, return the YAML text. Otherwise write to *stream*.
    """
    options_yaml = yaml.safe_dump(available_yaml_options(), sort_keys=False)
    if stream is None:
        return options_yaml
    stream.write(options_yaml)
    if not options_yaml.endswith("\n"):
        stream.write("\n")


class Configuration:
    """Load and validate the YAML configuration file."""

    def __init__(self, config_path):
        self.path = Path(config_path)
        self.base_dir = self.path.resolve().parent
        self.file_specs = self._load_file_specs()

    def _load_file_specs(self):
        with open(self.path, "r", encoding="utf-8") as config_file:
            config = yaml.safe_load(config_file)

        if not isinstance(config, dict):
            raise ValueError("YAML config must contain a top-level mapping")
        file_entries = config.get("files")
        if not isinstance(file_entries, list) or not file_entries:
            raise ValueError("YAML config must contain a non-empty files list")

        specs = []
        seen_labels = set()
        for index, entry in enumerate(file_entries, start=1):
            spec = self._parse_file_entry(index, entry)
            if spec.label in seen_labels:
                raise ValueError(f"Duplicate file label {spec.label!r}")
            seen_labels.add(spec.label)
            specs.append(spec)
        return specs

    def _parse_file_entry(self, index, entry):
        if not isinstance(entry, dict):
            raise ValueError(f"files[{index}] must be a mapping")

        raw_name = entry.get("name")
        if not raw_name:
            raise ValueError(f"files[{index}] is missing required field 'name'")

        raw_format = entry.get("format")
        if not raw_format:
            raise ValueError(f"files[{index}] is missing required field 'format'")
        file_format = str(raw_format).lower()
        if file_format not in ("npz", "h5", "hdf5"):
            raise ValueError(
                f"files[{index}].format must be npz, h5, or hdf5, got {raw_format!r}"
            )

        path = Path(raw_name)
        if not path.is_absolute():
            path = self.base_dir / path

        channel_mapping, tick_info = self._parse_channel_mapping(
            self._get_config_value(entry, "channel_mapping", "channel-mapping")
        )
        return FileSpec(
            label=str(entry.get("label") or path.stem),
            path=path,
            file_format=file_format,
            dataset=entry.get("dset"),
            xaxis_label=self._get_config_value(
                entry,
                "xaxis_label",
                "xaxis-label",
                default="tick",
            ),
            yaxis_label=self._get_config_value(
                entry,
                "yaxis_label",
                "yaxis-label",
                default="channel within view",
            ),
            channel_mapping=channel_mapping,
            tick_info=tick_info,
            trace_yaxis_label=self._optional_string(
                self._get_config_value(
                    entry,
                    "trace_yaxis_label",
                    "trace-yaxis-label",
                )
            ),
            separate_axis=self._parse_separate_axis(entry.get("separate_axis")),
        )

    @staticmethod
    def _get_config_value(mapping, *keys, default=None):
        for key in keys:
            if key in mapping:
                return mapping[key]
        return default

    @classmethod
    def _parse_separate_axis(cls, value):
        if value is None:
            return SeparateAxisConfig()
        if isinstance(value, dict):
            enabled = cls._parse_bool(value.get("enabled"), default=True)
            return SeparateAxisConfig(
                enabled=enabled,
                xaxis_label=cls._optional_string(
                    cls._get_config_value(value, "xaxis_label", "xaxis-label")
                ),
                yaxis_label=cls._optional_string(
                    cls._get_config_value(value, "yaxis_label", "yaxis-label")
                ),
            )
        return SeparateAxisConfig(enabled=cls._parse_bool(value, default=False))

    @staticmethod
    def _optional_string(value):
        if value is None:
            return None
        value = str(value).strip()
        return value or None

    @classmethod
    def _parse_channel_mapping(cls, raw_mapping):
        if raw_mapping is None:
            return dict(DEFAULT_VIEW_RANGES), None
        if not isinstance(raw_mapping, dict):
            raise ValueError("channel_mapping must be a mapping")

        channel_mapping = {}
        for view_name in VIEW_ORDER:
            raw_range = raw_mapping.get(view_name)
            if raw_range is None:
                raw_range = raw_mapping.get(view_name.lower())
            if raw_range is None:
                channel_mapping[view_name] = DEFAULT_VIEW_RANGES[view_name]
            else:
                channel_mapping[view_name] = cls._parse_range(
                    raw_range,
                    f"channel_mapping.{view_name}",
                )

        tick_info = cls._get_config_value(raw_mapping, "tick_info", "tick-info")
        if tick_info is not None:
            tick_info = int(tick_info)
            if tick_info <= 0:
                raise ValueError("channel_mapping.tick_info must be positive")
        return channel_mapping, tick_info

    @staticmethod
    def _parse_range(value, label):
        if isinstance(value, (list, tuple)) and len(value) == 2:
            start, stop = value
        elif isinstance(value, str):
            match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", value)
            if not match:
                raise ValueError(
                    f"{label} must be [start, stop] or 'start - stop', got {value!r}"
                )
            start, stop = match.groups()
        else:
            raise ValueError(
                f"{label} must be [start, stop] or 'start - stop', got {value!r}"
            )

        start = int(start)
        stop = int(stop)
        if start < 0 or stop <= start:
            raise ValueError(f"{label} has invalid channel range {start}:{stop}")
        return start, stop

    @staticmethod
    def _parse_bool(value, default=False):
        if value is None:
            return default
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            lowered = value.strip().lower()
            if lowered in ("true", "yes", "1", "on"):
                return True
            if lowered in ("false", "no", "0", "off"):
                return False
        raise ValueError(f"Expected boolean value, got {value!r}")
