#!/usr/bin/env python3
"""
Interactive viewer for .pt files with u/v/w plane separation.
Displays all files for a given plane simultaneously with interactive crop selection.
"""
import torch
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons, RectangleSelector, Button
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path
import re
import argparse
import sys
import traceback


class PTFileViewer:
    def __init__(self, directory=".", urange=None, vrange=None, wrange=None):
        self.all_files = sorted(glob.glob(f"{directory}/*.pt"))
        if not self.all_files:
            raise ValueError(f"No .pt files found in {directory}")

        # Load first file to get total channel count
        sample_data = torch.load(self.all_files[0])
        total_channels = sample_data['y'].shape[1]

        # Channel regions for each plane (default to full range)
        self.regions = {
            'u': urange if urange else (0, total_channels),
            'v': vrange if vrange else (0, total_channels),
            'w': wrange if wrange else (0, total_channels)
        }

        # Group files by plane
        self.files_by_plane = {'u': [], 'v': [], 'w': []}
        for filepath in self.all_files:
            plane = self.get_plane_from_filename(filepath)
            self.files_by_plane[plane].append(filepath)

        # Sort v-plane by t2 (second angle)
        self.files_by_plane['v'] = sorted(
            self.files_by_plane['v'],
            key=lambda f: self.extract_t_values(f)[1]
        )

        self.current_plane = 'u'
        self.current_item = 'y'
        self.n_ticks = sample_data['y'].shape[2]

        # Crop ranges (None means use full range)
        self.crop_rect = None  # (tick_start, tick_end, ch_start, ch_end)

        # Setup the plot
        self.setup_plot()

    def get_plane_from_filename(self, filepath):
        """Extract plane (u/v/w) from filename"""
        filename = Path(filepath).name.lower()
        for plane in ['uplane', 'vplane', 'wplane']:
            if plane in filename:
                return plane[0]
        return 'u'

    def extract_t_values(self, filepath):
        """Extract t1 and t2 values from filename"""
        filename = Path(filepath).name
        match = re.search(r't1-(\d+)-t2-(\d+)', filename)
        if match:
            return int(match.group(1)), int(match.group(2))
        return 0, 0

    def setup_plot(self):
        """Setup the matplotlib figure with controls"""
        self.fig = plt.figure(figsize=(18, 10))

        # Control panel on the left
        control_left = 0.02
        control_width = 0.10

        # Radio buttons for selecting item
        ax_item_radio = self.fig.add_axes([control_left, 0.75, control_width, 0.15])
        ax_item_radio.set_title('Item', fontsize=10, pad=5)
        self.item_radio = RadioButtons(ax_item_radio, ('y', 'labels', 'feat'), active=0)
        self.item_radio.on_clicked(self.on_item_change)

        # Radio buttons for selecting plane
        ax_plane_radio = self.fig.add_axes([control_left, 0.55, control_width, 0.15])
        ax_plane_radio.set_title('Plane', fontsize=10, pad=5)
        self.plane_radio = RadioButtons(ax_plane_radio, ('u', 'v', 'w'), active=0)
        self.plane_radio.on_clicked(self.on_plane_change)

        # Reset crop button
        ax_reset_button = self.fig.add_axes([control_left, 0.40, control_width, 0.05])
        self.reset_button = Button(ax_reset_button, 'Reset Crop')
        self.reset_button.on_clicked(self.reset_crop)

        # Instructions text
        self.fig.text(control_left, 0.30,
                     'Draw a box on\n the top left image to\ncrop all plots',
                     fontsize=9, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Main plot area
        self.plot_left = 0.14
        self.plot_right = 0.95
        self.plot_bottom = 0.08
        self.plot_top = 0.92

        self.selector = None
        self.title_text = None

        # Initial plot
        self.update_plot()

    def on_select(self, eclick, erelease):
        """Handle rectangle selection"""
        x1, x2 = sorted([eclick.xdata, erelease.xdata])
        y1, y2 = sorted([eclick.ydata, erelease.ydata])

        # Convert to integer indices
        tick_start = max(0, int(x1))
        tick_end = min(self.n_ticks, int(x2))

        # Get channel size for current plane
        start_ch, end_ch = self.regions[self.current_plane]
        channel_size = end_ch - start_ch

        ch_start = max(0, int(y1))
        ch_end = min(channel_size, int(y2))

        # Store crop rectangle
        self.crop_rect = (tick_start, tick_end, ch_start, ch_end)

        # Update all plots
        self.update_plot()

    def reset_crop(self, event):
        """Reset crop to full view"""
        self.crop_rect = None
        self.update_plot()

    def update_plot(self):
        """Update the plot with current selections"""
        # Remove all subplot axes but keep control widgets
        control_axes = [self.item_radio.ax, self.plane_radio.ax, self.reset_button.ax]

        for ax in self.fig.axes[:]:
            if ax not in control_axes:
                self.fig.delaxes(ax)

        files = self.files_by_plane[self.current_plane]
        if not files:
            return

        # Calculate grid layout
        n_files = len(files)
        n_cols = min(4, n_files)
        n_rows = (n_files + n_cols - 1) // n_cols

        # Get channel range for this plane
        start_ch, end_ch = self.regions[self.current_plane]
        channel_size = end_ch - start_ch

        # Apply crop
        if self.crop_rect:
            tick_start, tick_end, ch_start, ch_end = self.crop_rect
        else:
            tick_start, tick_end = 0, self.n_ticks
            ch_start, ch_end = 0, channel_size

        # Adjust to absolute channel indices
        abs_ch_start = start_ch + ch_start
        abs_ch_end = start_ch + ch_end

        # Load all data
        all_data = []
        for filepath in files:
            data = torch.load(filepath)
            item_data = data[self.current_item][0, abs_ch_start:abs_ch_end, tick_start:tick_end].cpu().numpy()
            all_data.append(item_data)

        # Set colormap range based on item
        if self.current_item in ['y', 'labels']:
            vmin, vmax = 0, 1
        else:  # feat - use automatic range
            vmin = min(d.min() for d in all_data)
            vmax = max(d.max() for d in all_data)

        # Calculate subplot positions
        plot_width = self.plot_right - self.plot_left
        plot_height = self.plot_top - self.plot_bottom

        subplot_width = plot_width / n_cols * 0.85
        subplot_height = plot_height / n_rows * 0.85

        h_spacing = plot_width / n_cols * 0.15
        v_spacing = plot_height / n_rows * 0.15

        # Create subplots
        axes_list = []
        for idx, (filepath, plot_data) in enumerate(zip(files, all_data)):
            row = idx // n_cols
            col = idx % n_cols

            # Calculate position
            left = self.plot_left + col * (subplot_width + h_spacing)
            bottom = self.plot_top - (row + 1) * (subplot_height + v_spacing)

            ax = self.fig.add_axes([left, bottom, subplot_width, subplot_height])
            axes_list.append(ax)

            # Plot heatmap
            im = ax.imshow(plot_data, aspect='auto', cmap='viridis',
                          interpolation='nearest', origin='lower',
                          vmin=vmin, vmax=vmax,
                          extent=[tick_start, tick_end, ch_start, ch_end])

            # Extract t1, t2 from filename
            t1, t2 = self.extract_t_values(filepath)
            ax.set_title(f't1={t1}, t2={t2}', fontsize=9)

            ax.set_xlabel('Tick', fontsize=8)
            ax.set_ylabel(f'Ch offset', fontsize=8)
            ax.tick_params(labelsize=7)

            # Add rectangle selector to first plot
            if idx == 0:
                self.selector = RectangleSelector(
                    ax, self.on_select,
                    useblit=True,
                    button=[1],  # Left mouse button
                    minspanx=5, minspany=5,
                    spancoords='data',
                    interactive=False,
                    props=dict(facecolor='red', edgecolor='red', alpha=0.3, fill=True)
                )

        # Add a single colorbar with fixed normalization
        if all_data:
            cbar_ax = self.fig.add_axes([self.plot_right + 0.01, self.plot_bottom, 0.02, self.plot_top - self.plot_bottom])
            # Create a ScalarMappable with the fixed normalization
            norm = Normalize(vmin=0, vmax=1)
            sm = ScalarMappable(cmap='viridis', norm=norm)
            sm.set_array([])
            cbar = self.fig.colorbar(sm, cax=cbar_ax, label='Value')
            cbar.ax.tick_params(labelsize=8)

        # Remove old title if it exists
        if self.title_text is not None:
            self.title_text.remove()

        # Add title
        crop_info = f"| Ticks: {tick_start}-{tick_end} | Channels: {abs_ch_start}-{abs_ch_end}" if self.crop_rect else "| Full view"
        self.title_text = self.fig.text(
            0.55, 0.96,
            f"Plane: {self.current_plane.upper()} | Item: {self.current_item} {crop_info}",
            ha='center', fontsize=10, fontweight='bold'
        )

        self.fig.canvas.draw_idle()

    def on_item_change(self, label):
        """Handle item selection change"""
        self.current_item = label
        self.update_plot()

    def on_plane_change(self, label):
        """Handle plane selection change"""
        self.current_plane = label
        # Reset crop when changing planes
        self.crop_rect = None
        self.update_plot()

    def show(self):
        """Display the interactive plot"""
        plt.show()


def parse_range(range_str):
    """Parse a range string like '0-800' or '0:800' into a tuple (start, end)"""
    if not range_str:
        return None

    # Support both '-' and ':' as separators
    if '-' in range_str:
        parts = range_str.split('-')
    elif ':' in range_str:
        parts = range_str.split(':')
    else:
        raise ValueError(f"Invalid range format: {range_str}. Use 'start-end' or 'start:end'")

    if len(parts) != 2:
        raise ValueError(f"Invalid range format: {range_str}. Expected 2 values.")

    try:
        start = int(parts[0])
        end = int(parts[1])
        return (start, end)
    except ValueError:
        raise ValueError(f"Invalid range format: {range_str}. Values must be integers.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Interactive PT file viewer')
    parser.add_argument('directory', nargs='?', default='.',
                       help='Directory containing .pt files (default: current directory)')
    parser.add_argument('--urange', type=str, default=None,
                       help='Channel range for u-plane (e.g., 0-800 or 0:800)')
    parser.add_argument('--vrange', type=str, default=None,
                       help='Channel range for v-plane (e.g., 800-1600 or 800:1600)')
    parser.add_argument('--wrange', type=str, default=None,
                       help='Channel range for w-plane (e.g., 1600-2560 or 1600:2560)')

    args = parser.parse_args()

    try:
        # Parse ranges
        urange = parse_range(args.urange)
        vrange = parse_range(args.vrange)
        wrange = parse_range(args.wrange)

        viewer = PTFileViewer(args.directory, urange=urange, vrange=vrange, wrange=wrange)
        print(f"\nLoaded {len(viewer.all_files)} files")
        print(f"  u-plane: {len(viewer.files_by_plane['u'])} files")
        print(f"  v-plane: {len(viewer.files_by_plane['v'])} files")
        print(f"  w-plane: {len(viewer.files_by_plane['w'])} files")
        print("\nControls:")
        print("  - Item radio buttons: Select item (y, labels, feat)")
        print("  - Plane radio buttons: Select plane (u, v, w)")
        print("  - Draw rectangle: Click and drag on the FIRST image to select crop area")
        print("  - Reset Crop button: Return to full view")
        print("\nPlane channel ranges:")
        print(f"  u: {viewer.regions['u'][0]}-{viewer.regions['u'][1]}")
        print(f"  v: {viewer.regions['v'][0]}-{viewer.regions['v'][1]}")
        print(f"  w: {viewer.regions['w'][0]}-{viewer.regions['w'][1]}")
        viewer.show()
    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()
        sys.exit(1)
