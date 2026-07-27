"""
Plot frames from DNNROI training HDF5 files.

This script reads fodder and truth HDF5 files and creates comparison plots
showing the dense (loose_lf), mp2, mp3 ROIs from fodder alongside the truth mask.

Usage:
    python plot_hdf5.py <prefix> <view> [event_num]
    
    prefix: The output prefix (e.g., "cosmics_10" for files like cosmics_10-g4-rec-0.h5)
    view: Either "0" or "1" (0=U view, 1=V view)
    event_num: Optional event number to plot (default: plots all events)

Example:
    python plot_hdf5.py cosmics_10 0 0
    python plot_hdf5.py cosmics_10 1 0
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
import sys
from pathlib import Path


def load_event_data(fodder_file, truth_file, event_num):
    """Load data for a specific event from both fodder and truth files."""
    data = {}
    
    with h5py.File(fodder_file, 'r') as f:
        event_str = str(event_num)
        if event_str not in f:
            print(f"Warning: Event {event_num} not found in {fodder_file}")
            return None
        
        event_group = f[event_str]
        
        # Load fodder data
        if 'frame_loose_lf0' in event_group:
            data['loose_lf'] = event_group['frame_loose_lf0'][:]
        if 'frame_mp2_roi0' in event_group:
            data['mp2'] = event_group['frame_mp2_roi0'][:]
        if 'frame_mp3_roi0' in event_group:
            data['mp3'] = event_group['frame_mp3_roi0'][:]
    
    with h5py.File(truth_file, 'r') as f:
        event_str = str(event_num)
        if event_str not in f:
            print(f"Warning: Event {event_num} not found in {truth_file}")
            return None
        
        event_group = f[event_str]
        
        # Load truth data
        if 'frame_ductor0' in event_group:
            data['truth'] = event_group['frame_ductor0'][:]
    
    return data


def plot_event(data, event_num, view, output_file=None):
    """Create a comparison plot of fodder and truth data."""
    
    # Map view number to label
    view_label = 'U' if view == '0' else 'V'
    
    # Create figure with 2 rows and 4 columns
    fig = plt.figure(figsize=(20, 10))
    gs = GridSpec(2, 4, figure=fig, hspace=0.3, wspace=0.3)
    
    fig.suptitle(f'DNNROI Training Data - Event {event_num} - View {view_label} (plane {view})', 
                 fontsize=16, fontweight='bold')
    
    # Determine vmin/vmax for consistent color scaling
    if 'loose_lf' in data:
        lf_min, lf_max = np.percentile(data['loose_lf'], [1, 99])
    else:
        lf_min, lf_max = 0, 1
    
    # Row 1: Fodder data
    # Plot 1: Dense (loose_lf)
    ax1 = fig.add_subplot(gs[0, 0])
    if 'loose_lf' in data:
        im1 = ax1.imshow(data['loose_lf'], aspect='auto', cmap='viridis',
                        vmin=lf_min, vmax=lf_max)
        plt.colorbar(im1, ax=ax1, label='ADC')
    ax1.set_title('Fodder: Dense (loose_lf)')
    ax1.set_xlabel('Time Tick')
    ax1.set_ylabel('Channel')
    
    # Plot 2: MP2 ROI
    ax2 = fig.add_subplot(gs[0, 1])
    if 'mp2' in data:
        im2 = ax2.imshow(data['mp2'], aspect='auto', cmap='Blues', vmin=0, vmax=1)
        plt.colorbar(im2, ax=ax2, label='MP2 ROI Mask')
    ax2.set_title('Fodder: MP2 ROI')
    ax2.set_xlabel('Time Tick')
    ax2.set_ylabel('Channel')
    
    # Plot 3: MP3 ROI
    ax3 = fig.add_subplot(gs[0, 2])
    if 'mp3' in data:
        im3 = ax3.imshow(data['mp3'], aspect='auto', cmap='Oranges', vmin=0, vmax=1)
        plt.colorbar(im3, ax=ax3, label='MP3 ROI Mask')
    ax3.set_title('Fodder: MP3 ROI')
    ax3.set_xlabel('Time Tick')
    ax3.set_ylabel('Channel')
    
    # Plot 4: Combined MP2 & MP3 with different colors
    ax4 = fig.add_subplot(gs[0, 3])
    if 'mp2' in data and 'mp3' in data:
        # Create a combined image: 0=background, 1=MP2 only, 2=MP3 only, 3=both
        combined = np.zeros_like(data['mp2'])
        mp2_mask = data['mp2'] > 0
        mp3_mask = data['mp3'] > 0
        
        combined[mp2_mask & ~mp3_mask] = 1  # MP2 only
        combined[~mp2_mask & mp3_mask] = 2  # MP3 only
        combined[mp2_mask & mp3_mask] = 3   # Both
        
        # Create custom colormap: white, blue, orange, red
        colors = ['white', 'blue', 'orange', 'red']
        n_bins = 4
        cmap = ListedColormap(colors)
        
        im4 = ax4.imshow(combined, aspect='auto', cmap=cmap, vmin=0, vmax=3)
        cbar = plt.colorbar(im4, ax=ax4, ticks=[0.375, 1.125, 1.875, 2.625])
        cbar.ax.set_yticklabels(['None', 'MP2', 'MP3', 'Both'])
    ax4.set_title('Fodder: MP2 (blue) | MP3 (orange)')
    ax4.set_xlabel('Time Tick')
    ax4.set_ylabel('Channel')
    
    # Row 2: Truth and comparisons
    # Plot 5: Truth - with same style as dense
    ax5 = fig.add_subplot(gs[1, 0])
    if 'truth' in data:
        im5 = ax5.imshow(data['truth'], aspect='auto', cmap='viridis',
                        vmin=0, vmax=data['truth'].max())
        plt.colorbar(im5, ax=ax5, label='Truth Values')
    ax5.set_title('Truth')
    ax5.set_xlabel('Time Tick')
    ax5.set_ylabel('Channel')
    
    # Plot 6: Overlay - Dense with Truth contours
    ax6 = fig.add_subplot(gs[1, 1])
    if 'loose_lf' in data and 'truth' in data:
        im6 = ax6.imshow(data['loose_lf'], aspect='auto', cmap='viridis',
                        vmin=lf_min, vmax=lf_max)
        # Add truth contours with same style
        if data['truth'].max() > 0:
            ax6.contour(data['truth'], levels=np.linspace(data['truth'].min(), 
                        data['truth'].max(), 5), colors='red', linewidths=2, alpha=0.7)
        plt.colorbar(im6, ax=ax6, label='ADC')
    ax6.set_title('Dense + Truth Contours')
    ax6.set_xlabel('Time Tick')
    ax6.set_ylabel('Channel')
    
    # Plot 7: Overlay - MP2/MP3 with Truth
    ax7 = fig.add_subplot(gs[1, 2])
    if 'mp2' in data and 'mp3' in data and 'truth' in data:
        # Create RGB image
        rgb = np.zeros((*data['truth'].shape, 3))
        rgb[:, :, 0] = data['truth'] / data['truth'].max() if data['truth'].max() > 0 else 0  # Red channel for truth
        rgb[:, :, 1] = data['mp2']    # Green channel for MP2
        rgb[:, :, 2] = data['mp3']    # Blue channel for MP3
        ax7.imshow(rgb, aspect='auto')
    ax7.set_title('Overlay: Truth(R) MP2(G) MP3(B)')
    ax7.set_xlabel('Time Tick')
    ax7.set_ylabel('Channel')
    
    # Plot 8: Statistics
    ax8 = fig.add_subplot(gs[1, 3])
    ax8.axis('off')
    
    # Calculate statistics
    stats_text = "Statistics:\n\n"
    
    if 'loose_lf' in data:
        stats_text += f"Dense (loose_lf):\n"
        loose_lf_coverage = (data['loose_lf']>0).sum() / data['loose_lf'].size * 100
        stats_text += f"  Coverage: {loose_lf_coverage:.2f}%\n"
        stats_text += f"  Shape: {data['loose_lf'].shape}\n"
        stats_text += f"  Min: {data['loose_lf'].min():.2f}\n"
        stats_text += f"  Max: {data['loose_lf'].max():.2f}\n"
        stats_text += f"  Mean: {data['loose_lf'].mean():.2f}\n\n"
    
    if 'mp2' in data:
        mp2_coverage = (data['mp2'] > 0).sum() / data['mp2'].size * 100
        stats_text += f"MP2 ROI:\n"
        stats_text += f"  Coverage: {mp2_coverage:.2f}%\n\n"
    
    if 'mp3' in data:
        mp3_coverage = (data['mp3'] > 0).sum() / data['mp3'].size * 100
        stats_text += f"MP3 ROI:\n"
        stats_text += f"  Coverage: {mp3_coverage:.2f}%\n\n"
    
    if 'truth' in data:
        truth_coverage = (data['truth'] > 0).sum() / data['truth'].size * 100
        stats_text += f"Truth:\n"
        stats_text += f"  Shape: {data['truth'].shape}\n"
        stats_text += f"  Coverage: {truth_coverage:.2f}%\n\n"
    
    ax8.text(0.1, 0.9, stats_text, transform=ax8.transAxes,
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {output_file}")
    else:
        plt.show()
    
    plt.close()


def main(prefix, view, event_num=None):
    """Main function to create plots."""
    
    # Construct file paths with new naming convention
    fodder_file = f"{prefix}-g4-rec-{view}.h5"
    truth_file = f"{prefix}-g4-tru-{view}.h5"
    
    # Check if files exist
    if not Path(fodder_file).exists():
        print(f"Error: Fodder file not found: {fodder_file}")
        return 1
    
    if not Path(truth_file).exists():
        print(f"Error: Truth file not found: {truth_file}")
        return 1
    
    print(f"Reading from:")
    print(f"  Fodder: {fodder_file}")
    print(f"  Truth: {truth_file}")
    
    # Get list of events
    with h5py.File(fodder_file, 'r') as f:
        event_nums = sorted([int(k) for k in f.keys()])
    
    print(f"Found {len(event_nums)} events: {event_nums}")
    
    # Plot specific event or all events
    if event_num is not None:
        if event_num not in event_nums:
            print(f"Error: Event {event_num} not found in files")
            return 1
        
        events_to_plot = [event_num]
    else:
        events_to_plot = event_nums
    
    # Create plots
    for evt in events_to_plot:
        print(f"\nPlotting event {evt}...")
        data = load_event_data(fodder_file, truth_file, evt)
        
        if data is None:
            print(f"Skipping event {evt} due to missing data")
            continue
        
        output_file = f"{prefix}_view{view}_event{evt}.png"
        plot_event(data, evt, view, output_file)
    
    print(f"\nPlotting complete!")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_hdf5.py <prefix> <view> [event_num]")
        print("  prefix: Output file prefix (e.g., 'cosmics_10')")
        print("  view: '0' or '1' (0=U view, 1=V view)")
        print("  event_num: Optional specific event to plot (default: all events)")
        print("\nExamples:")
        print("  python plot_hdf5_spngmp.py <path-to-h5-files>/cosmics_10 0      # Plot all events for U view (plane 0)")
        print("  python plot_hdf5_spngmp.py <path-to-h5-files>/cosmics_10 1 0    # Plot event 0 for V view (plane 1)")
        sys.exit(1)
    
    prefix = sys.argv[1]
    view = sys.argv[2]
    
    if view not in ['0', '1']:
        print("Error: view must be '0' (U) or '1' (V)")
        sys.exit(1)
    
    event_num = int(sys.argv[3]) if len(sys.argv) > 3 else None
    
    sys.exit(main(prefix, view, event_num))