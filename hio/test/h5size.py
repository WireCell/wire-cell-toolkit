#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "click",
#     "h5py",
#     "numpy",
#     "hdf5plugin",
#     "matplotlib",
# ]
# ///

"""
HDF5 file analysis and repacking tool.

Provides utilities to analyze HDF5 file storage characteristics and
repack files with different chunking, compression, and transformation options.
"""

import os
import sys
import math
import time
import resource
import tempfile
import itertools
import json
from pathlib import Path
import click
import h5py
import numpy as np

try:
    import hdf5plugin
except ImportError:
    hdf5plugin = None


def format_bytes(num_bytes):
    """Format bytes into human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if abs(num_bytes) < 1024.0:
            return f"{num_bytes:.2f} {unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:.2f} PB"


def calculate_chunk_shape(dataset_shape, chunk_size_bytes, element_size, method='major'):
    """
    Calculate chunk shape based on dataset shape, target chunk size, and method.

    Args:
        dataset_shape: Tuple of dataset dimensions
        chunk_size_bytes: Target chunk size in bytes
        element_size: Size of each element in bytes
        method: 'major' or 'square'

    Returns:
        Tuple of chunk dimensions
    """
    ndim = len(dataset_shape)
    if ndim == 0:
        return None

    max_elements = chunk_size_bytes // element_size
    if max_elements < 1:
        max_elements = 1

    if method == 'major':
        # Last dimension gets min(chunk_size, dim_size)
        # Second-to-last gets sized to fit chunk_size
        # Earlier dimensions get size 1
        chunk_shape = [1] * ndim

        # Set last dimension
        chunk_shape[-1] = min(dataset_shape[-1], max_elements)

        if ndim > 1:
            # Calculate remaining elements budget
            remaining_elements = max_elements // chunk_shape[-1]
            if remaining_elements > 0:
                chunk_shape[-2] = min(dataset_shape[-2], remaining_elements)

        return tuple(chunk_shape)

    elif method == 'square':
        # All dimensions get equal chunk size (as close as possible)
        target_dim_size = int(max_elements ** (1.0 / ndim))
        if target_dim_size < 1:
            target_dim_size = 1

        chunk_shape = []
        for dim_size in dataset_shape:
            chunk_shape.append(min(dim_size, target_dim_size))

        # Verify we don't exceed chunk_size_bytes
        total_elements = np.prod(chunk_shape)
        if total_elements * element_size > chunk_size_bytes:
            # Adjust down
            scale = (chunk_size_bytes / (total_elements * element_size)) ** (1.0 / ndim)
            chunk_shape = [max(1, int(c * scale)) for c in chunk_shape]

        return tuple(chunk_shape)

    else:
        raise ValueError(f"Unknown chunk method: {method}")


def convert_to_voxels(array):
    """
    Convert array to WarpConvNet Voxels sparse representation.

    Returns a dict with:
        - features: 1D array of non-zero values
        - coords: 2D array of coordinates (N x ndim)
        - offsets: indices for batch boundaries (if applicable)
    """
    # Find non-zero elements
    nonzero_coords = np.array(np.nonzero(array)).T  # Shape: (N, ndim)
    features = array[tuple(np.nonzero(array))]  # 1D array of values

    # For single batch, offsets would be [0, len(features)]
    offsets = np.array([0, len(features)], dtype=np.int64)

    return {
        'features': features,
        'coords': nonzero_coords.astype(np.int32),
        'offsets': offsets
    }


def visit_all_datasets(h5obj, func, path=''):
    """
    Recursively visit all datasets in HDF5 file/group.

    Args:
        h5obj: h5py.File or h5py.Group
        func: Function to call on each dataset, receives (path, dataset)
        path: Current path (used for recursion)
    """
    if isinstance(h5obj, h5py.Dataset):
        func(path, h5obj)
    elif isinstance(h5obj, h5py.Group):
        for key in h5obj.keys():
            new_path = f"{path}/{key}" if path else key
            visit_all_datasets(h5obj[key], func, new_path)


@click.group()
def cli():
    """HDF5 file analysis and repacking tool."""
    pass


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
def stats(filename):
    """
    Read HDF5 file and print storage statistics.

    Shows logical storage size, size on disk, and compression ratio
    for each dataset and the overall file.
    """
    click.echo(f"Analyzing: {filename}")
    click.echo("=" * 70)

    total_logical = 0
    total_storage = 0

    def print_dataset_stats(path, dataset):
        nonlocal total_logical, total_storage

        # Logical size (total elements * element size in bytes)
        logical_size = dataset.size * dataset.dtype.itemsize
        # Storage size (actual space allocated on disk)
        storage_size = dataset.id.get_storage_size()

        total_logical += logical_size
        total_storage += storage_size

        ratio = logical_size / storage_size if storage_size > 0 else float('inf')

        click.echo(f"\nDataset: {path}")
        click.echo(f"  Shape: {dataset.shape}, dtype: {dataset.dtype}")
        click.echo(f"  Logical Size: {format_bytes(logical_size)} ({logical_size} bytes)")
        click.echo(f"  Storage Size: {format_bytes(storage_size)} ({storage_size} bytes)")
        click.echo(f"  Ratio: {ratio:.2f}x")

        if dataset.chunks:
            click.echo(f"  Chunks: {dataset.chunks}")
        if dataset.compression:
            click.echo(f"  Compression: {dataset.compression}")
            if dataset.compression_opts:
                click.echo(f"  Compression opts: {dataset.compression_opts}")
        if dataset.shuffle:
            click.echo(f"  Shuffle: enabled")
        if dataset.fillvalue is not None:
            click.echo(f"  Fill value: {dataset.fillvalue}")

    with h5py.File(filename, 'r') as f:
        visit_all_datasets(f, print_dataset_stats)

    # Overall file statistics
    physical_size = os.path.getsize(filename)

    click.echo("\n" + "=" * 70)
    click.echo("SUMMARY")
    click.echo(f"  Total Logical Size: {format_bytes(total_logical)} ({total_logical} bytes)")
    click.echo(f"  Total Storage Size: {format_bytes(total_storage)} ({total_storage} bytes)")
    click.echo(f"  Physical File Size: {format_bytes(physical_size)} ({physical_size} bytes)")

    if total_storage > 0:
        logical_ratio = total_logical / total_storage
        click.echo(f"  Logical/Storage Ratio: {logical_ratio:.2f}x")

    if physical_size > 0:
        physical_ratio = total_logical / physical_size
        click.echo(f"  Logical/Physical Ratio: {physical_ratio:.2f}x")


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
@click.option('--chunk-size', default=1024, type=int,
              help='Target chunk size in bytes (default: 1024)')
@click.option('--chunk-method', type=click.Choice(['major', 'square']), default='major',
              help='Chunking strategy: major (last dims prioritized) or square (equal dims)')
@click.option('--fill-value', default=None, type=float,
              help='HDF5 fill value for sparse datasets')
@click.option('--shuffle/--no-shuffle', default=False,
              help='Enable HDF5 shuffle filter')
@click.option('--compression-method', type=click.Choice(['gzip', 'zstd']), default=None,
              help='Compression method (default: None)')
@click.option('--compression-level', default=None, type=int,
              help='Compression level (gzip: 0-9, zstd: 1-22)')
@click.option('--transform', type=click.Choice(['voxels']), default=None,
              help='Data transformation: voxels (WarpConvNet sparse representation)')
def repack(input_file, output_file, chunk_size, chunk_method, fill_value,
           shuffle, compression_method, compression_level, transform):
    """
    Repack HDF5 file with different chunking, compression, and transformation options.

    Reads INPUT_FILE and writes to OUTPUT_FILE with specified options.
    Can optionally transform data (e.g., to sparse voxel representation).
    """
    # Validate compression options
    if compression_method == 'zstd' and hdf5plugin is None:
        click.echo("Error: zstd compression requires 'hdf5plugin' package", err=True)
        click.echo("Install with: pip install hdf5plugin", err=True)
        sys.exit(1)

    # Set compression options
    compression = compression_method
    compression_opts = compression_level

    if compression == 'gzip' and compression_opts is None:
        compression_opts = 4  # Default gzip level
    elif compression == 'zstd':
        if compression_opts is None:
            compression_opts = 3  # Default zstd level
        compression = hdf5plugin.Zstd(clevel=compression_opts)
        compression_opts = None  # Plugin handles opts

    click.echo(f"Repacking: {input_file} -> {output_file}")
    click.echo(f"  Chunk size: {chunk_size} bytes")
    click.echo(f"  Chunk method: {chunk_method}")
    click.echo(f"  Fill value: {fill_value}")
    click.echo(f"  Shuffle: {shuffle}")
    click.echo(f"  Compression: {compression_method}")
    if compression_level:
        click.echo(f"  Compression level: {compression_level}")
    click.echo(f"  Transform: {transform}")
    click.echo()

    with h5py.File(input_file, 'r') as fin, h5py.File(output_file, 'w') as fout:

        def copy_dataset(path, dataset):
            click.echo(f"Processing: {path}")

            # Read data
            data = dataset[...]

            # Apply transformation if requested
            if transform == 'voxels':
                voxels = convert_to_voxels(data)

                # Create group with original dataset name
                grp = fout.create_group(path)

                # Write voxel components
                for name, array in voxels.items():
                    element_size = array.dtype.itemsize
                    chunks = calculate_chunk_shape(array.shape, chunk_size,
                                                   element_size, chunk_method)

                    grp.create_dataset(
                        name,
                        data=array,
                        chunks=chunks,
                        compression=compression,
                        compression_opts=compression_opts,
                        shuffle=shuffle,
                        fillvalue=fill_value
                    )
                    click.echo(f"  -> {path}/{name}: shape={array.shape}, chunks={chunks}")
            else:
                # Calculate chunking
                element_size = dataset.dtype.itemsize
                chunks = calculate_chunk_shape(dataset.shape, chunk_size,
                                              element_size, chunk_method)

                # Create dataset with new options
                fout.create_dataset(
                    path,
                    data=data,
                    chunks=chunks,
                    compression=compression,
                    compression_opts=compression_opts,
                    shuffle=shuffle,
                    fillvalue=fill_value,
                    dtype=dataset.dtype
                )

                click.echo(f"  Shape: {dataset.shape} -> chunks={chunks}")

        # Copy all datasets
        visit_all_datasets(fin, copy_dataset)

        # Copy attributes at root level
        for key, value in fin.attrs.items():
            fout.attrs[key] = value

    # Report results
    input_size = os.path.getsize(input_file)
    output_size = os.path.getsize(output_file)
    ratio = input_size / output_size if output_size > 0 else float('inf')

    click.echo()
    click.echo("=" * 70)
    click.echo("REPACK COMPLETE")
    click.echo(f"  Input size:  {format_bytes(input_size)} ({input_size} bytes)")
    click.echo(f"  Output size: {format_bytes(output_size)} ({output_size} bytes)")
    click.echo(f"  Ratio: {ratio:.2f}x")
    if output_size < input_size:
        savings = (1 - output_size / input_size) * 100
        click.echo(f"  Savings: {savings:.1f}%")
    elif output_size > input_size:
        overhead = (output_size / input_size - 1) * 100
        click.echo(f"  Overhead: {overhead:.1f}%")


def create_baseline(input_file, output_file):
    """
    Create a baseline HDF5 file with no compression, no chunking, no transforms.

    Args:
        input_file: Source HDF5 file
        output_file: Destination baseline file

    Returns:
        File size in bytes
    """
    with h5py.File(input_file, 'r') as fin, h5py.File(output_file, 'w') as fout:
        def copy_uncompressed(path, dataset):
            data = dataset[...]
            fout.create_dataset(
                path,
                data=data,
                chunks=None,  # No chunking
                compression=None,
                shuffle=False,
                dtype=dataset.dtype
            )

        visit_all_datasets(fin, copy_uncompressed)

        # Copy attributes
        for key, value in fin.attrs.items():
            fout.attrs[key] = value

    return os.path.getsize(output_file)


def run_repack_timed(input_file, output_file, chunk_size, chunk_method,
                     fill_value, shuffle, compression_method, compression_level,
                     transform):
    """
    Run repack with timing and return statistics.

    Returns:
        Dict with 'wall_time' (seconds), 'core_time' (seconds), 'size' (bytes), 'ratio' (vs input)
    """
    # Start timing - both wall clock and CPU time
    start_wall = time.time()
    start_usage = resource.getrusage(resource.RUSAGE_SELF)

    # Set up compression
    compression = compression_method
    compression_opts = compression_level

    if compression == 'gzip' and compression_opts is None:
        compression_opts = 4
    elif compression == 'zstd':
        if compression_opts is None:
            compression_opts = 3
        if hdf5plugin:
            compression = hdf5plugin.Zstd(clevel=compression_opts)
            compression_opts = None

    # Perform repack
    with h5py.File(input_file, 'r') as fin, h5py.File(output_file, 'w') as fout:
        def copy_dataset(path, dataset):
            data = dataset[...]

            if transform == 'voxels':
                voxels = convert_to_voxels(data)
                grp = fout.create_group(path)

                for name, array in voxels.items():
                    element_size = array.dtype.itemsize
                    chunks = calculate_chunk_shape(array.shape, chunk_size,
                                                   element_size, chunk_method)
                    grp.create_dataset(
                        name,
                        data=array,
                        chunks=chunks,
                        compression=compression,
                        compression_opts=compression_opts,
                        shuffle=shuffle,
                        fillvalue=fill_value
                    )
            else:
                element_size = dataset.dtype.itemsize
                chunks = calculate_chunk_shape(dataset.shape, chunk_size,
                                              element_size, chunk_method)
                fout.create_dataset(
                    path,
                    data=data,
                    chunks=chunks,
                    compression=compression,
                    compression_opts=compression_opts,
                    shuffle=shuffle,
                    fillvalue=fill_value,
                    dtype=dataset.dtype
                )

        visit_all_datasets(fin, copy_dataset)

        for key, value in fin.attrs.items():
            fout.attrs[key] = value

    # Calculate elapsed times
    wall_time = time.time() - start_wall
    end_usage = resource.getrusage(resource.RUSAGE_SELF)

    # Core time = user time + system time (across all cores/threads)
    user_time = end_usage.ru_utime - start_usage.ru_utime
    sys_time = end_usage.ru_stime - start_usage.ru_stime
    core_time = user_time + sys_time

    output_size = os.path.getsize(output_file)
    input_size = os.path.getsize(input_file)
    ratio = input_size / output_size if output_size > 0 else float('inf')

    return {
        'wall_time': wall_time,
        'core_time': core_time,
        'user_time': user_time,
        'sys_time': sys_time,
        'size': output_size,
        'ratio': ratio
    }


def get_file_stats(filename):
    """
    Get comprehensive stats for an HDF5 file.

    Returns:
        Dict with 'logical_size', 'storage_size', 'physical_size'
    """
    total_logical = 0
    total_storage = 0

    with h5py.File(filename, 'r') as f:
        def accumulate_stats(path, dataset):
            nonlocal total_logical, total_storage
            total_logical += dataset.size * dataset.dtype.itemsize
            total_storage += dataset.id.get_storage_size()

        visit_all_datasets(f, accumulate_stats)

    physical_size = os.path.getsize(filename)

    return {
        'logical_size': total_logical,
        'storage_size': total_storage,
        'physical_size': physical_size
    }


def read_file_timed(filename):
    """
    Read all datasets from HDF5 file into memory with timing.

    Returns:
        Dict with 'read_wall_time', 'read_core_time', 'read_user_time', 'read_sys_time'
    """
    # Start timing - both wall clock and CPU time
    start_wall = time.time()
    start_usage = resource.getrusage(resource.RUSAGE_SELF)

    # Read all datasets into memory
    with h5py.File(filename, 'r') as f:
        def read_dataset(path, dataset):
            # Force full read into memory
            _ = dataset[...]

        visit_all_datasets(f, read_dataset)

    # Calculate elapsed times
    read_wall_time = time.time() - start_wall
    end_usage = resource.getrusage(resource.RUSAGE_SELF)

    # Core time = user time + system time (across all cores/threads)
    read_user_time = end_usage.ru_utime - start_usage.ru_utime
    read_sys_time = end_usage.ru_stime - start_usage.ru_stime
    read_core_time = read_user_time + read_sys_time

    return {
        'read_wall_time': read_wall_time,
        'read_core_time': read_core_time,
        'read_user_time': read_user_time,
        'read_sys_time': read_sys_time
    }


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output-dir', type=click.Path(), default=None,
              help='Directory for temporary files (default: system temp)')
@click.option('--keep-files/--no-keep-files', default=False,
              help='Keep intermediate repacked files')
@click.option('--results', type=click.Path(), default=None,
              help='Output JSON file for full results')
def benchmark(input_file, output_dir, keep_files, results):
    """
    Benchmark different repack options and report compression ratios and timing.

    Creates a baseline uncompressed file, then tests the full outer product of
    chunking, compression, and transformation options. Reports a table showing
    the compression ratio and time for each configuration.
    """
    click.echo(f"Benchmarking: {input_file}")
    click.echo("=" * 100)

    # Create working directory
    if output_dir:
        work_dir = Path(output_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        temp_context = None
    else:
        temp_context = tempfile.TemporaryDirectory()
        work_dir = Path(temp_context.name)

    try:
        # Create baseline
        click.echo("\nCreating baseline (uncompressed, no chunking)...")
        baseline_file = work_dir / "baseline.h5"
        baseline_size = create_baseline(input_file, str(baseline_file))
        baseline_stats = get_file_stats(str(baseline_file))
        click.echo(f"  Baseline size: {format_bytes(baseline_size)} ({baseline_size} bytes)")
        click.echo(f"  Logical size:  {format_bytes(baseline_stats['logical_size'])}")
        click.echo(f"  Storage size:  {format_bytes(baseline_stats['storage_size'])}")

        # Define parameter space dimensions.
        chunk_sizes = [1024, 65536, 1048576]
        chunk_methods = ['major', 'square']
        #fill_values = [None, 0]
        fill_values = [None]
        shuffle_opts = [False, True]
        compression_configs = [
            (None, None),
            ('gzip', 1),
            ('gzip', 2),
            ('gzip', 3),

        ]
        if hdf5plugin:
            compression_configs.extend([
                ('zstd', 1),
                ('zstd', 5),
                ('zstd', 9),
            ])
        #transforms = [None, 'voxels']
        transforms = [None]

        # Calculate total combinations
        total_combinations = (len(chunk_sizes) * len(chunk_methods) * len(fill_values) *
                            len(shuffle_opts) * len(compression_configs) * len(transforms))

        click.echo(f"\nTesting {total_combinations} combinations...")
        click.echo(f"  Chunk sizes: {chunk_sizes}")
        click.echo(f"  Chunk methods: {chunk_methods}")
        click.echo(f"  Fill values: {fill_values}")
        click.echo(f"  Shuffle: {shuffle_opts}")
        click.echo(f"  Compressions: {[f'{m}:{l}' if m else 'None' for m, l in compression_configs]}")
        click.echo(f"  Transforms: {transforms}")
        click.echo()

        # Generate all combinations using itertools.product
        benchmark_results = []

        for i, (chunk_size, chunk_method, fill_value, shuffle,
                (comp_method, comp_level), transform) in enumerate(
            itertools.product(chunk_sizes, chunk_methods, fill_values,
                            shuffle_opts, compression_configs, transforms), 1):

            output_file = work_dir / f"test_{i:04d}.h5"

            try:
                result = run_repack_timed(
                    str(baseline_file),
                    str(output_file),
                    chunk_size=chunk_size,
                    chunk_method=chunk_method,
                    fill_value=fill_value,
                    shuffle=shuffle,
                    compression_method=comp_method,
                    compression_level=comp_level,
                    transform=transform
                )

                # Measure read performance
                read_timing = read_file_timed(str(output_file))
                result.update(read_timing)

                # Store configuration and results
                config = {
                    'chunk_size': chunk_size,
                    'chunk_method': chunk_method,
                    'fill_value': fill_value,
                    'shuffle': shuffle,
                    'compression_method': comp_method,
                    'compression_level': comp_level,
                    'transform': transform
                }

                result['config'] = config
                result['baseline_ratio'] = baseline_size / result['size'] if result['size'] > 0 else float('inf')
                benchmark_results.append(result)

                # One-line progress summary
                comp_str = f"{comp_method}:{comp_level}" if comp_method else "None"
                fill_str = str(fill_value) if fill_value is not None else "None"
                trans_str = transform if transform else "None"

                # Calculate parallelism factors
                write_parallelism = result['core_time'] / result['wall_time'] if result['wall_time'] > 0 else 0
                read_parallelism = result['read_core_time'] / result['read_wall_time'] if result['read_wall_time'] > 0 else 0

                click.echo(f"[{i:4d}/{total_combinations:4d}] "
                          f"chunk={chunk_size:6d}/{chunk_method:6s} "
                          f"fill={fill_str:4s} shuffle={str(shuffle):5s} "
                          f"comp={comp_str:9s} trans={trans_str:6s} | "
                          f"ratio={result['baseline_ratio']:5.2f}x "
                          f"W:{result['wall_time']:5.3f}s/{write_parallelism:4.1f}x "
                          f"R:{result['read_wall_time']:5.3f}s/{read_parallelism:4.1f}x "
                          f"size={format_bytes(result['size']):>10s}")

                # Clean up unless keeping files
                if not keep_files:
                    output_file.unlink()

            except Exception as e:
                click.echo(f"[{i:4d}/{total_combinations:4d}] Error: {e}", err=True)

        # Sort results by compression ratio (best first)
        benchmark_results.sort(key=lambda x: x['baseline_ratio'], reverse=True)

        # Write JSON results if requested
        if results:
            click.echo(f"\nWriting results to {results}...")
            json_data = {
                'input_file': input_file,
                'baseline': {
                    'size': baseline_size,
                    'logical_size': baseline_stats['logical_size'],
                    'storage_size': baseline_stats['storage_size']
                },
                'total_combinations': total_combinations,
                'results': benchmark_results
            }
            with open(results, 'w') as f:
                json.dump(json_data, f, indent=2)
            click.echo(f"Results written to {results}")

        # Print summary table
        click.echo("\n" + "=" * 145)
        click.echo("BENCHMARK SUMMARY (sorted by compression ratio)")
        click.echo("=" * 145)
        click.echo(f"{'Rank':>4s} {'Chunk':>6s} {'Method':>8s} {'Comp':>6s} {'Lvl':>4s} "
                  f"{'Shuf':>5s} {'Fill':>5s} {'Trans':>6s} {'Ratio':>7s} "
                  f"{'W-Wall':>7s} {'W-Core':>7s} {'W-Para':>6s} "
                  f"{'R-Wall':>7s} {'R-Core':>7s} {'R-Para':>6s} {'Size':>12s}")
        click.echo("-" * 145)

        # Show top 20 results
        # display_count = min(20, len(benchmark_results))
        display_count = len(benchmark_results)
        for rank, result in enumerate(benchmark_results[:display_count], 1):
            cfg = result['config']
            trans_str = cfg['transform'] if cfg['transform'] else "None"
            write_parallelism = result['core_time'] / result['wall_time'] if result['wall_time'] > 0 else 0
            read_parallelism = result['read_core_time'] / result['read_wall_time'] if result['read_wall_time'] > 0 else 0

            click.echo(f"{rank:4d} "
                      f"{cfg['chunk_size']:6d} "
                      f"{cfg['chunk_method']:>8s} "
                      f"{str(cfg['compression_method']):>6s} "
                      f"{str(cfg['compression_level']):>4s} "
                      f"{str(cfg['shuffle']):>5s} "
                      f"{str(cfg['fill_value']):>5s} "
                      f"{trans_str:>6s} "
                      f"{result['baseline_ratio']:7.2f}x "
                      f"{result['wall_time']:6.3f}s "
                      f"{result['core_time']:6.3f}s "
                      f"{write_parallelism:5.2f}x "
                      f"{result['read_wall_time']:6.3f}s "
                      f"{result['read_core_time']:6.3f}s "
                      f"{read_parallelism:5.2f}x "
                      f"{format_bytes(result['size']):>12s}")

        if len(benchmark_results) > display_count:
            click.echo(f"... ({len(benchmark_results) - display_count} more results)")
            click.echo(f"Use --results option to save full results to JSON")

        # Best configuration
        if benchmark_results:
            best = benchmark_results[0]
            best_write_parallelism = best['core_time'] / best['wall_time'] if best['wall_time'] > 0 else 0
            best_read_parallelism = best['read_core_time'] / best['read_wall_time'] if best['read_wall_time'] > 0 else 0

            click.echo("\n" + "=" * 100)
            click.echo("BEST CONFIGURATION")
            click.echo(f"  Compression ratio: {best['baseline_ratio']:.2f}x")
            click.echo(f"  Size: {format_bytes(best['size'])}")
            click.echo(f"  Write timing:")
            click.echo(f"    Wall time: {best['wall_time']:.3f}s")
            click.echo(f"    Core time: {best['core_time']:.3f}s (user: {best['user_time']:.3f}s, sys: {best['sys_time']:.3f}s)")
            click.echo(f"    Parallelism: {best_write_parallelism:.2f}x")
            click.echo(f"  Read timing:")
            click.echo(f"    Wall time: {best['read_wall_time']:.3f}s")
            click.echo(f"    Core time: {best['read_core_time']:.3f}s (user: {best['read_user_time']:.3f}s, sys: {best['read_sys_time']:.3f}s)")
            click.echo(f"    Parallelism: {best_read_parallelism:.2f}x")
            click.echo(f"  Options:")
            for key, value in best['config'].items():
                click.echo(f"    --{key.replace('_', '-')}: {value}")

    finally:
        # Clean up temporary directory
        if temp_context:
            temp_context.cleanup()
        elif not keep_files and work_dir.exists():
            # Clean up baseline if not keeping files
            if baseline_file.exists():
                baseline_file.unlink()


def parse_param_filter(param_specs):
    """
    Parse parameter filter specifications.

    Args:
        param_specs: List of strings like "chunk_size=1024,4096" or "compression_method=gzip"

    Returns:
        Dict mapping parameter names to lists of allowed values
    """
    filters = {}
    for spec in param_specs:
        if '=' not in spec:
            click.echo(f"Warning: Invalid param spec '{spec}', should be 'name=value1,value2,...'", err=True)
            continue

        name, values_str = spec.split('=', 1)
        name = name.strip()

        # Parse values - handle different types
        values = []
        for v in values_str.split(','):
            v = v.strip()
            # Try to convert to appropriate type
            if v.lower() == 'none':
                values.append(None)
            elif v.lower() in ('true', 'false'):
                values.append(v.lower() == 'true')
            elif v.isdigit():
                values.append(int(v))
            else:
                values.append(v)

        filters[name] = values

    return filters


def filter_results(results, param_filters):
    """
    Filter benchmark results based on parameter constraints.

    Args:
        results: List of result dicts with 'config' field
        param_filters: Dict mapping parameter names to allowed values

    Returns:
        Filtered list of results
    """
    if not param_filters:
        return results

    # Check for invalid parameter names once
    warned_params = set()
    if results:
        sample_config = results[0]['config']
        for param_name in param_filters.keys():
            if param_name not in sample_config:
                click.echo(f"Warning: Parameter '{param_name}' not found in config, ignoring filter", err=True)
                warned_params.add(param_name)

    filtered = []
    for result in results:
        config = result['config']
        matches = True

        for param_name, allowed_values in param_filters.items():
            if param_name in warned_params:
                continue  # Skip invalid parameters

            if config[param_name] not in allowed_values:
                matches = False
                break

        if matches:
            filtered.append(result)

    return filtered


def get_compression_marker(compression_method):
    """Get matplotlib marker style for compression method."""
    markers = {
        None: 'o',      # circle
        'gzip': 's',    # square
        'zstd': 'D',    # diamond
    }
    return markers.get(compression_method, 'o')


def format_config_label(config, exclude_keys=None):
    """
    Format configuration as a compact label.

    Args:
        config: Configuration dict
        exclude_keys: Set of keys to exclude from label

    Returns:
        Compact string representation
    """
    if exclude_keys is None:
        exclude_keys = set()

    parts = []

    # Add chunk info if not excluded
    if 'chunk_size' not in exclude_keys:
        cs = config['chunk_size']//1024
        parts.append(f"c={cs}kB")
    if 'chunk_method' not in exclude_keys:
        parts.append(f"{config['chunk_method'][:3]}")

    # Add compression info if not excluded
    if 'compression_method' not in exclude_keys:
        comp = config['compression_method']
        level = config['compression_level']
        if comp:
            comp_str = f"{comp}"
            if 'compression_level' not in exclude_keys and level is not None:
                comp_str += f":{level}"
            parts.append(comp_str)
        else:
            parts.append("nocomp")

    # Add other params if not excluded
    if 'shuffle' not in exclude_keys and config.get('shuffle'):
        parts.append("shuf")
    if 'fill_value' not in exclude_keys and config.get('fill_value') is not None:
        parts.append(f"fill={config['fill_value']}")
    if 'transform' not in exclude_keys and config.get('transform'):
        parts.append(f"t={config['transform']}")

    return ' '.join(parts)


def plot_bar_chart(ax, results, value_key, title, xlabel, ylabel='Configuration'):
    """
    Create a horizontal bar chart sorted by value.

    Args:
        ax: Matplotlib axes
        results: List of result dicts
        value_key: Key to extract value from results
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
    """
    # Sort by value (descending)
    sorted_results = sorted(results, key=lambda x: x[value_key], reverse=True)

    # Extract values and labels
    values = [r[value_key] for r in sorted_results]
    labels = [format_config_label(r['config']) for r in sorted_results]

    # Create bar chart
    y_pos = np.arange(len(labels))
    ax.barh(y_pos, values)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis='x', alpha=0.3)


def plot_scatter(ax, results, x_key, y_key, color_key, x_label, y_label,
                color_label, title, cmap='gray'):
    """
    Create a scatter plot with colored points and different markers.

    Args:
        ax: Matplotlib axes
        results: List of result dicts
        x_key: Key for x-axis values
        y_key: Key for y-axis values
        color_key: Key for color values
        x_label: X-axis label
        y_label: Y-axis label
        color_label: Color bar label
        title: Plot title
        cmap: Colormap name
    """
    import matplotlib.pyplot as plt

    # Group by compression method for different markers
    comp_groups = {}
    for result in results:
        comp_method = result['config']['compression_method']
        if comp_method not in comp_groups:
            comp_groups[comp_method] = []
        comp_groups[comp_method].append(result)

    # Plot each group with its marker
    scatter_objs = []
    for comp_method, group_results in comp_groups.items():
        x_vals = [r[x_key] for r in group_results]
        y_vals = [r[y_key] for r in group_results]
        c_vals = [r[color_key] for r in group_results]
        marker = get_compression_marker(comp_method)

        # Determine what to exclude from labels
        exclude = {x_key.replace('_', ' '), y_key.replace('_', ' '),
                  color_key.replace('_', ' '), 'compression_method'}
        labels = [format_config_label(r['config'], exclude) for r in group_results]

        scatter = ax.scatter(x_vals, y_vals, c=c_vals, marker=marker,
                           s=100, alpha=0.7, cmap=cmap,
                           label=f"{comp_method or 'none'}")
        scatter_objs.append(scatter)

        # Add labels for points
        for x, y, label in zip(x_vals, y_vals, labels):
            if label:  # Only add label if there's info not shown elsewhere
                ax.annotate(label, (x, y), fontsize=6, alpha=0.6,
                          xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(title='Compression', loc='best')

    # Add colorbar using the last scatter object
    if scatter_objs:
        cbar = plt.colorbar(scatter_objs[-1], ax=ax)
        cbar.set_label(color_label)


@cli.command()
@click.argument('json_file', type=click.Path(exists=True))
@click.option('-o', '--output', required=True, type=click.Path(),
              help='Output PDF filename')
@click.option('-p', '--param', 'params', multiple=True,
              help='Filter parameter: name=value1,value2,... (can be used multiple times)')
@click.option('--cmap', default='gray',
              help='Matplotlib colormap name (default: gray)')
def plot(json_file, output, params, cmap):
    """
    Create visualization plots from benchmark results.

    Reads a JSON file produced by the benchmark command and generates
    a multi-page PDF with various plots analyzing compression performance.

    Use --param to filter results, e.g., --param compression_method=gzip,zstd
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    click.echo(f"Reading benchmark results from {json_file}...")

    # Load JSON data
    with open(json_file, 'r') as f:
        data = json.load(f)

    results = data['results']
    baseline_size = data['baseline']['size']

    click.echo(f"Loaded {len(results)} benchmark results")

    # Parse and apply filters
    param_filters = parse_param_filter(params)
    if param_filters:
        click.echo(f"Applying filters: {param_filters}")
        results = filter_results(results, param_filters)
        click.echo(f"Filtered to {len(results)} results")

    if not results:
        click.echo("Error: No results after filtering", err=True)
        return

    # Calculate derived metrics
    for result in results:
        # Write speed (bytes/sec)
        result['write_speed'] = baseline_size / result['wall_time'] if result['wall_time'] > 0 else 0
        # Read speed (bytes/sec)
        result['read_speed'] = baseline_size / result['read_wall_time'] if result['read_wall_time'] > 0 else 0

    click.echo(f"Creating plots in {output}...")

    # Create multi-page PDF
    with PdfPages(output) as pdf:

        # Page 1: Bar chart - Compression Ratio
        fig, ax = plt.subplots(figsize=(11, 8.5))
        plot_bar_chart(ax, results, 'baseline_ratio',
                      'Compression Ratio (Baseline / Compressed)',
                      'Compression Ratio', 'Configuration')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2: Bar chart - Write Speed
        fig, ax = plt.subplots(figsize=(11, 8.5))
        plot_bar_chart(ax, results, 'write_speed',
                      'Write Speed (Compression)',
                      'Write Speed (bytes/sec)', 'Configuration')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 3: Bar chart - Read Speed
        fig, ax = plt.subplots(figsize=(11, 8.5))
        plot_bar_chart(ax, results, 'read_speed',
                      'Read Speed (Decompression)',
                      'Read Speed (bytes/sec)', 'Configuration')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 4: Scatter - Write vs Read Time, colored by compression ratio
        fig, ax = plt.subplots(figsize=(11, 8.5))
        plot_scatter(ax, results, 'wall_time', 'read_wall_time', 'baseline_ratio',
                    'Write Time (seconds)', 'Read Time (seconds)',
                    'Compression Ratio',
                    'Read vs Write Performance (colored by compression ratio)',
                    cmap=cmap)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 5: Scatter - Write vs Read Time, colored by chunk size
        fig, ax = plt.subplots(figsize=(11, 8.5))
        # Extract chunk sizes for color mapping
        for result in results:
            result['chunk_size_val'] = result['config']['chunk_size']

        plot_scatter(ax, results, 'wall_time', 'read_wall_time', 'chunk_size_val',
                    'Write Time (seconds)', 'Read Time (seconds)',
                    'Chunk Size (bytes)',
                    'Read vs Write Performance (colored by chunk size)',
                    cmap=cmap)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Add metadata
        d = pdf.infodict()
        d['Title'] = 'HDF5 Compression Benchmark Results'
        d['Author'] = 'h5size.py'
        d['Subject'] = f'Analysis of {json_file}'
        d['Keywords'] = 'HDF5, compression, benchmark'

    click.echo(f"PDF created: {output}")


if __name__ == '__main__':
    cli()
