"""Command-line interface for the YAML-configured frame visualizer."""

import argparse
import sys

from frame_config import Configuration, dump_yaml_options
from frame_data_loaders import DataLoaderFactory
from frame_visualizer_ui import FrameVisualizerUI


def load_config(config_path):
    """Backward-compatible helper for callers that used the old function."""
    return Configuration(config_path).file_specs


def load_frame(spec, frame_index):
    """Backward-compatible helper for callers that used the old function."""
    return DataLoaderFactory.create(spec, frame_index).load()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize and compare NPZ/HDF5 frame datasets from YAML."
    )
    parser.add_argument(
        "config",
        nargs="?",
        help="YAML config file with a top-level files list",
    )
    parser.add_argument(
        "--dump-yaml-options",
        action="store_true",
        help="print all supported YAML configuration options and exit",
    )
    parser.add_argument(
        "--frame-index",
        type=int,
        default=-1,
        help="frame ordinal for sorted frame_* NPZ keys or repeated HDF5 basename matches",
    )
    parser.add_argument(
        "--trace-half-window",
        type=int,
        default=100,
        help="number of ticks shown on each side of the selected tick",
    )
    args = parser.parse_args()
    if not args.dump_yaml_options and not args.config:
        parser.error("config is required unless --dump-yaml-options is used")
    return args


def main():
    args = parse_args()
    if args.dump_yaml_options:
        dump_yaml_options(sys.stdout)
        return

    config = Configuration(args.config)
    loaded_frames = [
        DataLoaderFactory.create(spec, args.frame_index).load()
        for spec in config.file_specs
    ]
    for loaded in loaded_frames:
        print(
            f"Loaded {loaded.spec.label} from "
            f"{loaded.spec.path}:{loaded.source_dataset} shape={loaded.frame.shape}"
        )
    FrameVisualizerUI(
        loaded_frames,
        trace_half_window=args.trace_half_window,
    ).show()


if __name__ == "__main__":
    main()
