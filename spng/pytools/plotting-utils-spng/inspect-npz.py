"""Inspect the contents of an NPZ file.

Prints each stored entry with shape, dtype, size, and optional numeric summary
statistics. Metadata entries such as JSON bytes are shown without forcing them
into arrays.
"""

import argparse
from pathlib import Path

import numpy as np


def numeric_summary(array):
    if not np.issubdtype(array.dtype, np.number):
        return ""
    if array.size == 0:
        return " empty"

    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return " finite=0"

    return (
        f" min={float(np.min(finite)):.6g}"
        f" max={float(np.max(finite)):.6g}"
        f" mean={float(np.mean(finite)):.6g}"
        f" finite={finite.size}/{array.size}"
    )


def describe_value(key, value, include_summary):
    if isinstance(value, bytes):
        print(f"{key}: type=bytes length={len(value)}")
        return

    array = np.asarray(value)
    summary = numeric_summary(array) if include_summary else ""
    print(
        f"{key}: shape={array.shape} dtype={array.dtype} "
        f"ndim={array.ndim} size={array.size}{summary}"
    )


def inspect_npz(filename, include_summary=False, allow_pickle=False):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"NPZ file does not exist: {path}")

    with np.load(path, allow_pickle=allow_pickle) as npz_file:
        print(f"file: {path}")
        print(f"entries: {len(npz_file.files)}")
        for key in npz_file.files:
            describe_value(key, npz_file[key], include_summary)


def parse_args():
    parser = argparse.ArgumentParser(description="Inspect keys and arrays in an NPZ file.")
    parser.add_argument("filename", help="NPZ file to inspect")
    parser.add_argument(
        "--summary",
        action="store_true",
        help="print min/max/mean/finite counts for numeric arrays",
    )
    parser.add_argument(
        "--allow-pickle",
        action="store_true",
        help="allow object arrays when loading the NPZ file",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    inspect_npz(
        args.filename,
        include_summary=args.summary,
        allow_pickle=args.allow_pickle,
    )


if __name__ == "__main__":
    main()
