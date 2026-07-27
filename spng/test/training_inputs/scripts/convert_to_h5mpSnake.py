"""
This script converts DNNROI training tar files to HDF5 format.
Works as a Snakemake script.

Input: tar file containing:
  - dnnroi-training-fodder.npz (dense, mp2, mp3 for U and V views)
  - dnnroi-training-truth.npz (truth for U, V, W views)

Output: HDF5 files following naming convention:
  - <prefix>-g4-rec-0.h5 (fodder data for U view)
  - <prefix>-g4-rec-1.h5 (fodder data for V view)
  - <prefix>-g4-tru-0.h5 (truth data for U view)
  - <prefix>-g4-tru-1.h5 (truth data for V view)
"""

import tarfile
import numpy as np
import h5py
import tempfile
from pathlib import Path


def extract_tar(tar_path, extract_dir):
    """Extract tar file to temporary directory."""
    with tarfile.open(tar_path, 'r') as tar:
        tar.extractall(extract_dir)
    return extract_dir


def load_npz_arrays(npz_path):
    """Load arrays from npz file, organizing by event and index."""
    data = np.load(npz_path)
    
    print(f"DEBUG: Files in {npz_path}:")
    for key in sorted(data.files):
        if '_array' in key:
            print(f"  {key}: shape={data[key].shape if hasattr(data[key], 'shape') else 'N/A'}")
    
    # Group arrays by event number (tensor_<event>_<index>_array)
    events = {}
    for key in data.files:
        # Look for keys with pattern: tensor_X_Y_array (may or may not have .npy extension)
        if 'tensor_' in key and '_array' in key:
            # Extract the numeric parts
            try:
                # Remove any .npy extension
                clean_key = key.replace('.npy', '')
                # Split on underscore
                parts = clean_key.split('_')
                
                # Find tensor, event, index positions
                if 'tensor' in parts:
                    tensor_idx = parts.index('tensor')
                    if len(parts) > tensor_idx + 2:
                        event_num = int(parts[tensor_idx + 1])
                        array_idx = int(parts[tensor_idx + 2])
                        
                        if event_num not in events:
                            events[event_num] = {}
                        events[event_num][array_idx] = data[key]
                        print(f"DEBUG: Loaded event={event_num}, idx={array_idx}, shape={data[key].shape}")
            except (ValueError, IndexError) as e:
                print(f"DEBUG: Could not parse key '{key}': {e}")
                continue
    
    print(f"DEBUG: Total events loaded: {len(events)}")
    for event_num in sorted(events.keys()):
        print(f"  Event {event_num}: indices {sorted(events[event_num].keys())}")
    
    return events


def process_fodder(npz_path, output_rec0=None, output_rec1=None, output_rec2=None):
    """
    Process fodder npz file and create HDF5 files for U and V views.
    
    Fodder array indices:
    0 --> Dense array (U) -> g4-rec-0.h5
    1 --> MP2 Information (U)
    2 --> MP3 Information (U)
    3 --> Dense array (V) -> g4-rec-1.h5
    4 --> MP2 Information (V)
    5 --> MP3 Information (V)
    """
    events = load_npz_arrays(npz_path)
    
    if not events:
        print("WARNING: No events found in fodder file!")
        return
    
    # Create HDF5 files for U (0) and V (1) views with new naming convention
    if output_rec0 is not None:
        h5_u = h5py.File(output_rec0, 'w')
    if output_rec1 is not None:
        h5_v = h5py.File(output_rec1, 'w')
    if output_rec2 is not None:
        h5_w = h5py.File(output_rec2, 'w')

    try:
        for event_num in sorted(events.keys()):
            event_data = events[event_num]
            print(f"Processing fodder event {event_num}")
            
            # Process U view (indices 0, 1, 2)
            if output_rec0 is not None:
                if 0 in event_data:
                    print(f"  Creating rec-0/{event_num}/frame_loose_lf0, shape={event_data[0].shape}")
                    h5_u.create_dataset(f"{event_num}/frame_loose_lf0", 
                                    data=event_data[0], 
                                    compression='gzip')
                if 1 in event_data:
                    print(f"  Creating rec-0/{event_num}/frame_mp2_roi0, shape={event_data[1].shape}")
                    h5_u.create_dataset(f"{event_num}/frame_mp2_roi0", 
                                    data=event_data[1], 
                                    compression='gzip')
                if 2 in event_data:
                    print(f"  Creating rec-0/{event_num}/frame_mp3_roi0, shape={event_data[2].shape}")
                    h5_u.create_dataset(f"{event_num}/frame_mp3_roi0", 
                                    data=event_data[2], 
                                    compression='gzip')
            
            # Process V view (indices 3, 4, 5)
            if output_rec1 is not None:
                if 3 in event_data:
                    print(f"  Creating rec-1/{event_num}/frame_loose_lf0, shape={event_data[3].shape}")
                    h5_v.create_dataset(f"{event_num}/frame_loose_lf0", 
                                    data=event_data[3], 
                                    compression='gzip')
                if 4 in event_data:
                    print(f"  Creating rec-1/{event_num}/frame_mp2_roi0, shape={event_data[4].shape}")
                    h5_v.create_dataset(f"{event_num}/frame_mp2_roi0", 
                                    data=event_data[4], 
                                    compression='gzip')
                if 5 in event_data:
                    print(f"  Creating rec-1/{event_num}/frame_mp3_roi0, shape={event_data[5].shape}")
                    h5_v.create_dataset(f"{event_num}/frame_mp3_roi0", 
                                    data=event_data[5], 
                                    compression='gzip')
    
            # Process V view (indices 3, 4, 5)
            if output_rec2 is not None:
                if 6 in event_data:
                    print(f"  Creating rec-2/{event_num}/frame_loose_lf0, shape={event_data[6].shape}")
                    h5_w.create_dataset(f"{event_num}/frame_loose_lf0", 
                                    data=event_data[6], 
                                    compression='gzip')
                if 7 in event_data:
                    print(f"  Creating rec-2/{event_num}/frame_mp2_roi0, shape={event_data[7].shape}")
                    h5_w.create_dataset(f"{event_num}/frame_mp2_roi0", 
                                    data=event_data[7], 
                                    compression='gzip')
                if 8 in event_data:
                    print(f"  Creating rec-2/{event_num}/frame_mp3_roi0, shape={event_data[8].shape}")
                    h5_w.create_dataset(f"{event_num}/frame_mp3_roi0", 
                                    data=event_data[8], 
                                    compression='gzip')
    finally:
        if output_rec0 is not None:
            h5_u.close()
        if output_rec1 is not None:
            h5_v.close()
        if output_rec2 is not None:
            h5_w.close()


def process_truth(npz_path, output_tru0=None, output_tru1=None, output_tru2=None):
    """
    Process truth npz file and create HDF5 files for U and V views.
    
    Truth array indices:
    0 --> Dense array (U) -> g4-tru-0.h5
    1 --> Dense array (V) -> g4-tru-1.h5
    2 --> Dense array (W1)
    3 --> Dense array (W2)
    """
    events = load_npz_arrays(npz_path)
    
    if not events:
        print("WARNING: No events found in truth file!")
        return
    
    # Create HDF5 files for U (0) and V (1) views with new naming convention
    if output_tru0 is not None:
        h5_u = h5py.File(output_tru0, 'w')
    if output_tru1 is not None:
        h5_v = h5py.File(output_tru1, 'w')
    if output_tru2 is not None:
        h5_w = h5py.File(output_tru2, 'w')
    
    try:
        for event_num in sorted(events.keys()):
            event_data = events[event_num]
            print(f"Processing truth event {event_num}")
            
            # Process U view (index 0) - use frame_ductor0 as dataset name
            if 0 in event_data and output_tru0 is not None:
                print(f"  Creating tru-0/{event_num}/frame_ductor0, shape={event_data[0].shape}")
                h5_u.create_dataset(f"{event_num}/frame_ductor0", 
                                   data=event_data[0], 
                                   compression='gzip')
            
            # Process V view (index 1) - use frame_ductor0 as dataset name
            if 1 in event_data and output_tru1 is not None:
                print(f"  Creating tru-1/{event_num}/frame_ductor0, shape={event_data[1].shape}")
                h5_v.create_dataset(f"{event_num}/frame_ductor0", 
                                   data=event_data[1], 
                                   compression='gzip')
            
            if 2 in event_data and 3 in event_data and output_tru2 is not None:
                w = np.concat([event_data[2], event_data[3]])
                print(f"  Creating tru-2/{event_num}/frame_ductor0, shape={w.shape}")
                h5_w.create_dataset(f"{event_num}/frame_ductor0", 
                                    data=w, 
                                    compression='gzip')
    finally:
        if output_tru0 is not None: h5_u.close()
        if output_tru1 is not None: h5_v.close()
        if output_tru2 is not None: h5_w.close()


def get_output(name, output):
    return (
        None if name not in output.keys()
        else output[name]
    )

# Main execution for Snakemake
if 'snakemake' in globals():
    tar_path = snakemake.input[0]

    output_rec0 = get_output('rec0', snakemake.output)
    output_rec1 = get_output('rec1', snakemake.output)
    output_rec2 = get_output('rec2', snakemake.output)
    output_tru0 = get_output('tru0', snakemake.output)
    output_tru1 = get_output('tru1', snakemake.output)
    output_tru2 = get_output('tru2', snakemake.output)

    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Extracting {tar_path} to {temp_dir}...")
        extract_tar(tar_path, temp_dir)
        
        # Find npz files - look for the exact filenames
        temp_path = Path(temp_dir)
        fodder_file = temp_path / "dnnroi-training-fodder.npz"
        truth_file = temp_path / "dnnroi-training-truth.npz"
        
        print(f"Looking for:")
        print(f"  fodder: {fodder_file} (exists: {fodder_file.exists()})")
        print(f"  truth: {truth_file} (exists: {truth_file.exists()})")
        
        # Fallback: search for files with these patterns
        if not fodder_file.exists():
            fodder_files = list(temp_path.glob("*fodder*.npz"))
            if fodder_files:
                fodder_file = fodder_files[0]
                print(f"  Found alternative fodder file: {fodder_file}")
        
        if not truth_file.exists():
            truth_files = list(temp_path.glob("*truth*.npz"))
            if truth_files:
                truth_file = truth_files[0]
                print(f"  Found alternative truth file: {truth_file}")
        
        if not fodder_file.exists() or not truth_file.exists():
            print("Error: Could not find both fodder and truth npz files in tar")
            print("Contents of tar:")
            for f in temp_path.glob("*"):
                print(f"  {f.name}")
            raise FileNotFoundError("Missing fodder or truth npz files")
        
        print(f"\nProcessing fodder file: {fodder_file}")
        process_fodder(str(fodder_file), output_rec0, output_rec1, output_rec2)
        print(f"Created {output_rec0} and {output_rec1}")
        
        print(f"\nProcessing truth file: {truth_file}")
        process_truth(str(truth_file), output_tru0, output_tru1, output_tru2)
        print(f"Created {output_tru0} and {output_tru1}")
    
    print("\nConversion complete!")