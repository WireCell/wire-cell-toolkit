"""Shared models and utility functions for frame visualization."""

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


DEFAULT_VIEW_RANGES = {
    "U": (0, 800),
    "V": (800, 1600),
    "W": (1600, 2560),
}
VIEW_ORDER = ("U", "V", "W")
MAX_ACTIVE_FILES = 3


@dataclass(frozen=True)
class SeparateAxisConfig:
    enabled: bool = False
    xaxis_label: str | None = None
    yaxis_label: str | None = None


@dataclass(frozen=True)
class FileSpec:
    label: str
    path: Path
    file_format: str
    dataset: str | None
    xaxis_label: str
    yaxis_label: str
    channel_mapping: dict[str, tuple[int, int]]
    tick_info: int | None
    trace_yaxis_label: str | None
    separate_axis: SeparateAxisConfig


@dataclass
class LoadedFrame:
    spec: FileSpec
    source_dataset: str
    frame: np.ndarray


@dataclass
class UISelection:
    active_view_name: str = VIEW_ORDER[0]
    selected_channel: int = 0
    selected_tick: int = 0
    zoom_limits: tuple[int, int, int, int] | None = None


def frame_sort_key(name):
    matches = re.findall(r"\d+", name)
    return tuple(int(match) for match in matches) if matches else (-1,)


def finite_limits(arrays):
    finite_arrays = [array[np.isfinite(array)] for array in arrays]
    finite_arrays = [array for array in finite_arrays if array.size]
    if not finite_arrays:
        return 0.0, 1.0

    combined = np.concatenate(finite_arrays)
    vmin = float(np.percentile(combined, 1.0))
    vmax = float(np.percentile(combined, 99.5))
    if vmin == vmax:
        vmax = vmin + 1.0
    return vmin, vmax
