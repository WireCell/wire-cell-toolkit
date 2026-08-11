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
def load_remap(filename):
    with open(filename, 'r') as f:
        remap = {i:l.strip('\n') for i, l in enumerate(f.readlines()) if '#' not in l}
    return remap

def process_fodder_dense(
        npz_path,
        output_rec0=None,
        output_rec1=None,
        output_rec2=None,
        remap=None):
    """
    Process fodder npz file and create HDF5 files for all views concatted (NO MP INFO).
    
    Fodder array indices:
    0 --> Dense array (U) -> g4-rec-0.h5
    1 --> Dense array (V)
    2 --> Dense array (W)
    """
    events = load_npz_arrays(npz_path)
    
    if not events:
        print("WARNING: No events found in fodder file!")
        return

    if remap is not None:
        remap_map = load_remap(remap)
    
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
            event_num_out = event_num
            if remap is not None:
                event_num_out = remap_map[event_num]
                print(f'\tRemapping to {event_num_out}')

            # # If any are missing that's bad
            # if not all(i in event_data for i in [0, 1, 2]):
            #     raise RuntimeError(f'Could not find all keys in input file {npz_path}, event {event_num}')
            # print(event_data[0].shape, event_data[1].shape, event_data[2].shape)
            # # Risky assumption about the shape/order of chans/ticks
            # full_data = np.concat([event_data[i] for i in range(3)])
            # print(f"  Creating rec-0/{event_num}/frame_loose_lf0, shape={full_data.shape}")
            # h5_file.create_dataset(f"{event_num}/frame_loose_lf0", 
            #                     data=full_data, 
            #                     compression='gzip')
            # Process U view (indices 0, 1, 2)
            if output_rec0 is not None:
                if 0 in event_data:
                    print(f"  Creating rec-0/{event_num_out}/frame_loose_lf0, shape={event_data[0].shape}")
                    h5_u.create_dataset(f"{event_num_out}/frame_loose_lf0", 
                                    data=event_data[0], 
                                    compression='gzip')
                if 1 in event_data:
                    print(f"  Creating rec-1/{event_num_out}/frame_loose_lf0, shape={event_data[1].shape}")
                    h5_v.create_dataset(f"{event_num_out}/frame_loose_lf0", 
                                    data=event_data[1], 
                                    compression='gzip')
                if 2 in event_data:
                    print(f"  Creating rec-2/{event_num_out}/frame_loose_lf0, shape={event_data[2].shape}")
                    h5_w.create_dataset(f"{event_num_out}/frame_loose_lf0", 
                                    data=event_data[2], 
                                    compression='gzip')
    finally:
        if output_rec0 is not None:
            h5_u.close()
        if output_rec1 is not None:
            h5_v.close()
        if output_rec2 is not None:
            h5_w.close()


if __name__ == '__main__':
    from argparse import ArgumentParser as ap
    parser = ap()
    parser.add_argument('-i', type=str, help='Input', required=True)
    parser.add_argument('-o', type=str, help='Output', required=True)
    parser.add_argument('--remap', type=str, help='File to remap index to known event nums')
    args = parser.parse_args()

    fodder_file = args.i
    output_rec0 = args.o.replace('placeholder', '0')
    output_rec1 = args.o.replace('placeholder', '1')
    output_rec2 = args.o.replace('placeholder', '2')

    print(f"\nProcessing fodder file: {fodder_file}")
    process_fodder_dense(fodder_file, output_rec0, output_rec1, output_rec2, remap=args.remap)
    print(f"Created {output_rec0}, {output_rec1}, {output_rec2}")