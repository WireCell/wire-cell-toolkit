#!/usr/bin/env -S uv run --script
# -*- python -*-
# /// script
# requires-python = ">=3.12"
# dependencies = ["pandas", "matplotlib"]
# ///
#!/usr/bin/env python3
"""
Make PDF plots from csv file produced by check_dft_measure.cxx
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
import math
import numpy as np # For handling potential NaN values
from itertools import product

cpu_linestyle='None'
cpu_marker='.'
gpu_linestyle='None'
gpu_marker='x'

#fft_func_names = ('fft', 'rfft', 'ifft', 'irfft', 'fft2', 'rfft2', 'ifft2', 'irfft2')
fft_func_names = ('fft', 'ifft')
fft_device_names = ('cpu', 'gpu')
default_columns = tuple([f'{f}_{d}_ms' for d,f in product(fft_device_names, fft_func_names)])

def get_faster_absolute_size_data(df, function_names, cuda_available):
    """
    Analyzes the DataFrame to find, for each original size, the 'faster' size
    (i.e., a size >= original_size with the smallest time) for each FFT function and device.
    This is the original implementation of 'get_faster_size_data'.
    """
    faster_size_map = {}
    print("\nAnalyzing data to determine 'faster absolute' sizes...")
    for func_name in function_names:
        for device_type in ['cpu', 'gpu']:
            if device_type == 'gpu' and not cuda_available:
                continue

            col_name = f'{func_name}_{device_type}_ms'
            if col_name not in df.columns:
                continue

            current_faster_data = []
            for i in range(len(df)):
                original_size = df['size'].iloc[i]
                
                # Slice the DataFrame from the current index onwards to find the min time
                slice_df = df.iloc[i:]
                if not slice_df.empty:
                    # idxmin() returns the original index from the DataFrame
                    # .dropna() to ensure min is not affected by NaN, if any
                    min_time_idx_global = slice_df[col_name].dropna().idxmin()
                    faster_size = df['size'].loc[min_time_idx_global]
                else:
                    faster_size = original_size 
                
                current_faster_data.append((original_size, faster_size))
            faster_size_map[(func_name, device_type)] = current_faster_data
    print("Finished analyzing 'faster absolute' sizes.")
    return faster_size_map

def get_faster_nearest_size_data(df, function_names, cuda_available):
    """
    Analyzes the DataFrame to find, for each original size N, the 'first faster' size.
    The 'first faster' size is the size at the first valley in the raw timing measurements
    at a size greater or equal to the original size.
    """
    first_faster_size_map = {}
    print("\nAnalyzing data to determine 'first faster (nearest valley)' sizes...")

    for func_name in function_names:
        for device_type in ['cpu', 'gpu']:
            if device_type == 'gpu' and not cuda_available:
                continue

            col_name = f'{func_name}_{device_type}_ms'
            if col_name not in df.columns:
                continue

            current_first_faster_data = []
            for i in range(len(df)):
                original_size = df['size'].iloc[i]
                original_time = df[col_name].iloc[i]

                candidate_size = original_size
                candidate_time = original_time

                # Start checking from the next size
                for j in range(i + 1, len(df)):
                    next_size = df['size'].iloc[j]
                    next_time = df[col_name].iloc[j]

                    # Handle NaN values gracefully
                    if pd.isna(candidate_time):
                        # If current candidate time is NaN, move to the next valid point if available
                        candidate_size = next_size
                        candidate_time = next_time
                        continue
                    if pd.isna(next_time):
                        # If next_time is NaN, it's not smaller, so consider candidate_time as the valley
                        break 
                    
                    if next_time < candidate_time:
                        candidate_size = next_size
                        candidate_time = next_time
                    else: # next_time >= candidate_time, we found a valley (or plateau)
                        break # Stop at the first increase or plateau

                current_first_faster_data.append((original_size, candidate_size))
            first_faster_size_map[(func_name, device_type)] = current_first_faster_data
    print("Finished analyzing 'first faster (nearest valley)' sizes.")
    return first_faster_size_map

def largest_prime_factor(n):
    """
    Finds the largest prime factor of a given integer n.
    """
    max_prime = -1

    # Handle the case of factor 2
    while n % 2 == 0:
        max_prime = 2
        n //= 2  # Equivalent to n = n / 2

    # Handle odd factors starting from 3
    # We only need to check up to the square root of n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            max_prime = i
            n //= i

    # If n is still greater than 2 after all divisions,
    # it means the remaining n is a prime factor itself
    if n > 2:
        max_prime = n

    return max_prime


def max_prime_colors(sizes):
    #         0       1    2        3      4     5       6     7
    colors = ["black",None,"red", "green", None, "blue", None, "orange"]
    primes = [largest_prime_factor(s) for s in sizes]
    return [colors[0 if p > 7 else p] for p in primes]

def plot_raw_timings(df, function_names, cuda_available, pdf):
    """
    Generates and saves 'Raw timing plots' to the PDF.
    """
    print("\nGenerating 'Raw timing plots':")
    for func_name in function_names:
        fig, ax = plt.subplots(figsize=(10, 6)) # Create a new figure for each plot
        
        # Plot CPU data
        cpu_col_name = f'{func_name}_cpu_ms'
        if cpu_col_name in df.columns:
            colors = max_prime_colors(df['size'])
            #ax.plot(df['size'], df[cpu_col_name], marker=cpu_marker, linestyle=cpu_linestyle, label='CPU')
            ax.scatter(df['size'], df[cpu_col_name], c=colors,
                    marker=cpu_marker, linestyle=cpu_linestyle, label='CPU')

        
        # Plot GPU data if available
        gpu_col_name = f'{func_name}_gpu_ms'
        if cuda_available and gpu_col_name in df.columns:
            ax.plot(df['size'], df[gpu_col_name], marker=gpu_marker, linestyle=gpu_linestyle, label='GPU')
        
        ax.set_title(f'Raw Timing for {func_name.upper()} Function')
        ax.set_xlabel('Tensor Dimension Size')
        ax.set_ylabel('Execution Time (ms)')
        ax.legend()
        ax.grid(True, which="both", ls="--", c='0.7')
        
        # Use log scale for Y-axis if the data range is very large
        y_data = []
        if cpu_col_name in df.columns: y_data.extend(df[cpu_col_name].dropna().tolist())
        if cuda_available and gpu_col_name in df.columns: y_data.extend(df[gpu_col_name].dropna().tolist())
        
        if y_data:
            min_val = min(val for val in y_data if val > 0)
            max_val = max(y_data)
            ax.set_yscale('log')
                
        pdf.savefig(fig) # Save the current figure to the PDF
        plt.close(fig)   # Close the figure to free memory
        print(f" - Plot created: Raw Timing for {func_name.upper()}")

def plot_best_timings_speedup(df, function_names, cuda_available, faster_size_map, pdf):
    """
    Generates and saves 'Best timing plots' (Speedup Ratio using absolute fastest size) to the PDF.
    """
    print("\nGenerating 'Best timing plots' (Speedup Ratio - Absolute Fastest):")
    for func_name in function_names:
        fig, ax = plt.subplots(figsize=(10, 6)) # Create a new figure for each plot
        
        # Calculate and plot CPU speedup ratio
        cpu_col_name = f'{func_name}_cpu_ms'
        if cpu_col_name in df.columns and (func_name, 'cpu') in faster_size_map:
            cpu_speedup_ratios = []
            for i, (original_size, faster_size) in enumerate(faster_size_map[(func_name, 'cpu')]):
                t_original = df[cpu_col_name].iloc[i]
                # Find the actual time for the faster_size
                t_faster_idx = df[df['size'] == faster_size].index[0]
                t_faster = df[cpu_col_name].loc[t_faster_idx]
                
                # Handle division by zero or negative times for robustness
                if t_original > 0 and not pd.isna(t_original) and not pd.isna(t_faster): 
                    ratio = t_original / t_faster
                else: # t_original is non-positive or NaN, or t_faster is NaN
                    ratio = 0.0 
                cpu_speedup_ratios.append(ratio)
            ax.plot(df['size'], cpu_speedup_ratios, marker=cpu_marker, linestyle=cpu_linestyle, label='CPU Speedup Ratio (Absolute Fastest)')
        
        # Calculate and plot GPU speedup ratio if available
        gpu_col_name = f'{func_name}_gpu_ms'
        if cuda_available and gpu_col_name in df.columns and (func_name, 'gpu') in faster_size_map:
            gpu_speedup_ratios = []
            for i, (original_size, faster_size) in enumerate(faster_size_map[(func_name, 'gpu')]):
                t_original = df[gpu_col_name].iloc[i]
                t_faster_idx = df[df['size'] == faster_size].index[0]
                t_faster = df[gpu_col_name].loc[t_faster_idx]
                
                if t_original > 0 and not pd.isna(t_original) and not pd.isna(t_faster):
                    ratio = t_original / t_faster
                else:
                    ratio = 0.0
                gpu_speedup_ratios.append(ratio)
            ax.plot(df['size'], gpu_speedup_ratios, marker=gpu_marker, linestyle=gpu_linestyle, label='GPU Speedup Ratio (Absolute Fastest)')

        ax.set_title(f'Speedup Ratio for {func_name.upper()} Function (Absolute Fastest)')
        ax.set_xlabel('Tensor Dimension Size')
        ax.set_ylabel('Speedup Ratio: (Original Time - Faster Time) / Original Time')
        ax.legend()
        ax.grid(True, which="both", ls="--", c='0.7')
        ax.axhline(0, color='red', linestyle=':', linewidth=0.8, label='No Speedup/Slowdown') # Baseline for no speedup
        ax.set_ylim(1, ax.get_ylim()[1])
        
        pdf.savefig(fig) # Save the current figure to the PDF
        plt.close(fig)   # Close the figure to free memory
        print(f" - Plot created: Speedup Ratio (Absolute Fastest) for {func_name.upper()}")

def plot_memory_ratios(df, function_names, cuda_available, faster_size_map, pdf):
    """
    Generates and saves 'Memory ratio plots' (using absolute fastest size) to the PDF.
    """
    print("\nGenerating 'Memory ratio plots' (Relative to Absolute Fastest Size):")
    for func_name in function_names:
        fig, ax = plt.subplots(figsize=(10, 6)) # Create a new figure for each plot

        # Plot CPU memory ratio
        if (func_name, 'cpu') in faster_size_map:
            cpu_faster_data = faster_size_map[(func_name, 'cpu')]
            cpu_original_sizes = [pair[0] for pair in cpu_faster_data]
            cpu_faster_sizes = [pair[1] for pair in cpu_faster_data]
            
            cpu_memory_ratios = []
            for orig_size, fast_size in zip(cpu_original_sizes, cpu_faster_sizes):
                if orig_size != 0:
                    #ratio = (orig_size - fast_size) / orig_size
                    ratio = fast_size / orig_size
                else:
                    ratio = 0.0 # Handle case where original_size is 0
                cpu_memory_ratios.append(ratio)
            ax.plot(cpu_original_sizes, cpu_memory_ratios, marker=cpu_marker, linestyle=cpu_linestyle, label='CPU Memory Ratio (Absolute Fastest)')

        # Plot GPU memory ratio if available
        if cuda_available and (func_name, 'gpu') in faster_size_map:
            gpu_faster_data = faster_size_map[(func_name, 'gpu')]
            gpu_original_sizes = [pair[0] for pair in gpu_faster_data]
            gpu_faster_sizes = [pair[1] for pair in gpu_faster_data]
            
            gpu_memory_ratios = []
            for orig_size, fast_size in zip(gpu_original_sizes, gpu_faster_sizes):
                if orig_size != 0:
                    #ratio = (orig_size - fast_size) / orig_size
                    ratio = fast_size / orig_size
                else:
                    ratio = 0.0
                gpu_memory_ratios.append(ratio)
            ax.plot(gpu_original_sizes, gpu_memory_ratios, marker=gpu_marker, linestyle=gpu_linestyle, label='GPU Memory Ratio (Absolute Fastest)')

        ax.set_title(f'Memory Ratio for {func_name.upper()} Function (Absolute Fastest)')
        ax.set_xlabel('Original Tensor Dimension Size')
        ax.set_ylabel('Memory Ratio: Faster Size / Original Size')
        ax.legend()
        ax.grid(True, which="both", ls="--", c='0.7')
        ax.axhline(0, color='red', linestyle=':', linewidth=0.8, label='No Size Change')

        #ax.set_ylim(0.001, ax.get_ylim()[1])
        ax.set_yscale('log')
        
        pdf.savefig(fig)
        plt.close(fig)
        print(f" - Plot created: Memory Ratio (Absolute Fastest) for {func_name.upper()}")

def plot_speedup_nearest(df, function_names, cuda_available, faster_nearest_size_map, pdf):
    """
    Generates and saves 'Speedup nearest' plots (Speedup Ratio using 'first faster' size) to the PDF.
    """
    print("\nGenerating 'Speedup nearest' plots (Speedup Ratio - First Valley):")
    for func_name in function_names:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Calculate and plot CPU speedup ratio
        cpu_col_name = f'{func_name}_cpu_ms'
        if cpu_col_name in df.columns and (func_name, 'cpu') in faster_nearest_size_map:
            cpu_speedup_ratios = []
            for i, (original_size, first_faster_size) in enumerate(faster_nearest_size_map[(func_name, 'cpu')]):
                t_original = df[cpu_col_name].iloc[i]
                # Find the actual time for the first_faster_size
                t_faster_idx = df[df['size'] == first_faster_size].index[0]
                t_faster = df[cpu_col_name].loc[t_faster_idx]
                
                if t_original > 0 and not pd.isna(t_original) and not pd.isna(t_faster): 
                    ratio = t_original / t_faster
                else: 
                    ratio = 0.0 
                cpu_speedup_ratios.append(ratio)
            ax.plot(df['size'], cpu_speedup_ratios, marker=cpu_marker, linestyle=cpu_linestyle, label='CPU Speedup Ratio (First Valley)')
        
        # Calculate and plot GPU speedup ratio if available
        gpu_col_name = f'{func_name}_gpu_ms'
        if cuda_available and gpu_col_name in df.columns and (func_name, 'gpu') in faster_nearest_size_map:
            gpu_speedup_ratios = []
            for i, (original_size, first_faster_size) in enumerate(faster_nearest_size_map[(func_name, 'gpu')]):
                t_original = df[gpu_col_name].iloc[i]
                t_faster_idx = df[df['size'] == first_faster_size].index[0]
                t_faster = df[gpu_col_name].loc[t_faster_idx]
                
                if t_original > 0 and not pd.isna(t_original) and not pd.isna(t_faster):
                    ratio = t_original /  t_faster
                else:
                    ratio = 0.0
                gpu_speedup_ratios.append(ratio)
            ax.plot(df['size'], gpu_speedup_ratios, marker=gpu_marker, linestyle=gpu_linestyle, label='GPU Speedup Ratio (First Valley)')

        ax.set_title(f'Speedup Ratio for {func_name.upper()} Function (First Valley)')
        ax.set_xlabel('Tensor Dimension Size')
        ax.set_ylabel('Speedup Ratio: (Original Time - First Valley Time) / Original Time')
        ax.legend()
        ax.grid(True, which="both", ls="--", c='0.7')
        ax.axhline(0, color='red', linestyle=':', linewidth=0.8, label='No Speedup/Slowdown')
        # ax.set_ylim(-1.1, 1.1) 
        ax.set_ylim(1.0, ax.get_ylim()[1])

        pdf.savefig(fig)
        plt.close(fig)
        print(f" - Plot created: Speedup Ratio (First Valley) for {func_name.upper()}")

def plot_memory_nearest(df, function_names, cuda_available, faster_nearest_size_map, pdf):
    """
    Generates and saves 'Memory nearest' plots (using 'first faster' size) to the PDF.
    """
    print("\nGenerating 'Memory nearest' plots (First Valley):")
    for func_name in function_names:
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot CPU memory ratio
        if (func_name, 'cpu') in faster_nearest_size_map:
            cpu_faster_data = faster_nearest_size_map[(func_name, 'cpu')]
            cpu_original_sizes = [pair[0] for pair in cpu_faster_data]
            cpu_first_faster_sizes = [pair[1] for pair in cpu_faster_data]
            
            cpu_memory_ratios = []
            for orig_size, fast_size in zip(cpu_original_sizes, cpu_first_faster_sizes):
                if orig_size != 0:
                    ratio = fast_size / orig_size
                else:
                    ratio = 0.0
                cpu_memory_ratios.append(ratio)
            ax.plot(cpu_original_sizes, cpu_memory_ratios, marker=cpu_marker, linestyle=cpu_linestyle, label='CPU Memory Ratio (First Valley)')

        # Plot GPU memory ratio if available
        if cuda_available and (func_name, 'gpu') in faster_nearest_size_map:
            gpu_faster_data = faster_nearest_size_map[(func_name, 'gpu')]
            gpu_original_sizes = [pair[0] for pair in gpu_faster_data]
            gpu_first_faster_sizes = [pair[1] for pair in gpu_faster_data]
            
            gpu_memory_ratios = []
            for orig_size, fast_size in zip(gpu_original_sizes, gpu_first_faster_sizes):
                if orig_size != 0:
                    ratio = fast_size / orig_size
                else:
                    ratio = 0.0
                gpu_memory_ratios.append(ratio)
            ax.plot(gpu_original_sizes, gpu_memory_ratios, marker=gpu_marker, linestyle=gpu_linestyle, label='GPU Memory Ratio (First Valley)')

        ax.set_title(f'Memory Ratio for {func_name.upper()} Function (First Valley)')
        ax.set_xlabel('Original Tensor Dimension Size')
        ax.set_ylabel('Memory Ratio: (Original Size - First Valley Size) / Original Size')
        ax.legend()
        ax.grid(True, which="both", ls="--", c='0.7')
        ax.axhline(0, color='red', linestyle=':', linewidth=0.8, label='No Size Change')
        #ax.set_ylim(min(-0.1, ax.get_ylim()[0]), max(1.1, ax.get_ylim()[1]))
        ax.set_ylim(1.0, ax.get_ylim()[1])

        pdf.savefig(fig)
        plt.close(fig)
        print(f" - Plot created: Memory Ratio (First Valley) for {func_name.upper()}")


def create_plots(csv_filename, pdf_filename, columns=None):
    """
    Reads a CSV file containing FFT benchmark data and generates five families of plots
    using matplotlib, saving them to a single PDF document.
    """
    if not columns:
        columns = default_columns

    if not os.path.exists(csv_filename):
        print(f"Error: CSV file '{csv_filename}' not found.")
        sys.exit(1)

    try:
        df = pd.read_csv(csv_filename)
    except Exception as e:
        print(f"Error reading CSV file '{csv_filename}': {e}")
        sys.exit(1)

    if df.empty:
        print("Error: CSV file is empty or could not be parsed.")
        sys.exit(1)
    
    # Check for 'size' column and ensure it's sorted
    if 'size' not in df.columns:
        print("Error: CSV file must contain a 'size' column.")
        sys.exit(1)
    df = df.sort_values(by='size').reset_index(drop=True)


    # Extract base FFT function names from column headers (e.g., 'fft', 'rfft', 'fft2')
    function_names = sorted(list(set([col.split('_')[0] for col in df.columns if col in columns])))
    if not function_names:
        print("Error: No FFT function data columns found in the CSV (e.g., 'fft_cpu_ms').")
        sys.exit(1)

    # Determine if GPU data is present in the CSV
    cuda_available = any("_gpu_ms" in col for col in df.columns)
    if not cuda_available:
        print("Warning: No GPU data found in CSV. GPU plots will be skipped.")

    # Initialize PDF document
    with PdfPages(pdf_filename) as pdf:
        # Family 1: Raw Timing plots
        plot_raw_timings(df, function_names, cuda_available, pdf)

        # Pre-calculate faster sizes for existing plots (Absolute Fastest)
        faster_absolute_size_map = get_faster_absolute_size_data(df, function_names, cuda_available)
        
        # Family 2: Best Timing (Speedup Ratio - Absolute Fastest) plots
        plot_best_timings_speedup(df, function_names, cuda_available, faster_absolute_size_map, pdf)

        # Family 3: Memory Ratio (Relative to Absolute Fastest Size) plots
        plot_memory_ratios(df, function_names, cuda_available, faster_absolute_size_map, pdf)

        # Pre-calculate faster sizes for new plots (First Valley)
        faster_nearest_size_map = get_faster_nearest_size_data(df, function_names, cuda_available)

        # Family 4: Speedup nearest (Speedup Ratio - First Valley) plots
        plot_speedup_nearest(df, function_names, cuda_available, faster_nearest_size_map, pdf)

        # Family 5: Memory nearest (Relative to First Valley Size) plots
        plot_memory_nearest(df, function_names, cuda_available, faster_nearest_size_map, pdf)

    print(f"\nAll plots successfully saved to '{pdf_filename}'")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plot_fft_benchmarks.py <input_csv_file> <output_pdf_file>")
        sys.exit(1)

    input_csv_file = sys.argv[1]
    output_pdf_file = sys.argv[2]
    create_plots(input_csv_file, output_pdf_file)
