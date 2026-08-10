"""Matplotlib UI for frame visualization."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, CheckButtons, RectangleSelector, TextBox

from frame_collection import FrameCollection
from frame_models import UISelection, VIEW_ORDER, finite_limits


class ProjectedTracePlot:
    """Draw projected 1D channel traces and their optional separate axes."""

    def __init__(self, figure, axis, frames, trace_half_window):
        self.figure = figure
        self.axis = axis
        self.frames = frames
        self.trace_half_window = trace_half_window
        self.extra_axes = []
        self.colors = plt.rcParams["axes.prop_cycle"].by_key().get(
            "color",
            ["C0", "C1", "C2", "C3", "C4", "C5"],
        )

    def trace_tick_limits(self, view_name, selected_tick):
        _, tick_count = self.frames.require_compatible_view_shape(view_name)
        selected_tick = int(np.clip(selected_tick, 0, tick_count - 1))
        tick_min = max(0, selected_tick - self.trace_half_window)
        tick_max = min(tick_count - 1, selected_tick + self.trace_half_window)
        return tick_min, tick_max

    def clear(self):
        for axis in self.extra_axes:
            axis.remove()
        self.extra_axes.clear()
        self.axis.clear()

    def plot(self, view_name, channel, selected_tick):
        self.clear()
        tick_min, tick_max = self.trace_tick_limits(view_name, selected_tick)
        tick_values = np.arange(tick_min, tick_max + 1)
        line_handles, line_labels = self._draw_lines(
            view_name,
            channel,
            selected_tick,
            tick_min,
            tick_max,
            tick_values,
        )
        self._format_axis(
            view_name,
            channel,
            selected_tick,
            tick_min,
            tick_max,
            line_handles,
            line_labels,
        )

    def _draw_lines(
        self,
        view_name,
        channel,
        selected_tick,
        tick_min,
        tick_max,
        tick_values,
    ):
        line_handles = []
        line_labels = []
        separate_axis_count = 0

        for color_index, loaded in enumerate(self.frames.active_frames()):
            array = self.frames.view_array(loaded, view_name)
            if not (0 <= channel < array.shape[0]):
                continue

            color = self.colors[color_index % len(self.colors)]
            trace_axis = self.axis
            if loaded.spec.separate_axis.enabled:
                trace_axis = self._new_separate_axis(separate_axis_count, color)
                separate_axis_count += 1

            line = trace_axis.plot(
                tick_values,
                array[channel, tick_min : tick_max + 1],
                label=loaded.spec.label,
                linewidth=1.0,
                color=color,
            )[0]
            trace_axis.plot(
                selected_tick,
                array[channel, selected_tick],
                marker="o",
                markersize=4,
                color=color,
            )
            if loaded.spec.separate_axis.enabled:
                self._style_separate_axis(
                    trace_axis,
                    loaded.spec,
                    color,
                    separate_axis_count - 1,
                )
            line_handles.append(line)
            line_labels.append(loaded.spec.label)
        return line_handles, line_labels

    def _new_separate_axis(self, axis_index, color):
        trace_axis = self.axis.twinx()
        trace_axis.spines["right"].set_position(("axes", 1.0 + 0.035 * axis_index))
        trace_axis.spines["top"].set_position(("axes", 1.0 + 0.15 * axis_index))
        trace_axis.spines["right"].set_color(color)
        trace_axis.spines["top"].set_color(color)
        self.extra_axes.append(trace_axis)
        return trace_axis

    @staticmethod
    def _style_separate_axis(axis, spec, color, axis_index):
        axis.set_ylabel(spec.separate_axis.yaxis_label or f"{spec.label} value")
        axis.yaxis.set_label_position("right")
        axis.yaxis.tick_right()
        axis.tick_params(axis="y", colors=color, pad=2)
        axis.yaxis.label.set_color(color)

        axis.set_xlabel(spec.separate_axis.xaxis_label or "tick number")
        axis.xaxis.set_label_position("top")
        axis.xaxis.tick_top()
        axis.tick_params(
            axis="x",
            colors=color,
            bottom=False,
            labelbottom=False,
            top=True,
            labeltop=False,
            pad=2,
        )
        axis.xaxis.label.set_color(color)
        axis.xaxis.labelpad = 6 + 8 * axis_index

    def _default_trace_yaxis_label(self):
        for loaded in self.frames.active_frames():
            if loaded.spec.separate_axis.enabled:
                continue
            if loaded.spec.trace_yaxis_label:
                return loaded.spec.trace_yaxis_label
        first_spec = self.frames.first_active_frame().spec
        return first_spec.trace_yaxis_label or "value"

    def _format_axis(
        self,
        view_name,
        channel,
        selected_tick,
        tick_min,
        tick_max,
        line_handles,
        line_labels,
    ):
        self.axis.axvline(
            selected_tick,
            color="black",
            linestyle="--",
            linewidth=1.0,
            alpha=0.65,
        )
        for axis in [self.axis, *self.extra_axes]:
            axis.set_xlim(tick_min, tick_max)
            axis.grid(axis="x", alpha=0.25)

        absolute_channel = self.frames.view_offset(
            self.frames.active_labels[0],
            view_name,
        ) + channel
        self.axis.set_xlabel("tick number")
        self.axis.set_ylabel(self._default_trace_yaxis_label())
        self.axis.yaxis.set_label_position("left")
        self.axis.yaxis.tick_left()
        self.axis.set_title(
            f"{view_name} value by tick for channel {channel} "
            f"(absolute channel={absolute_channel})",
            pad=22 if self.extra_axes else None,
        )
        if line_handles:
            self.axis.legend(line_handles, line_labels, loc="upper right")
        self.axis.grid(True, alpha=0.25)


class FrameVisualizerUI:
    """Matplotlib UI for browsing and comparing loaded frame data."""

    def __init__(self, loaded_frames, trace_half_window=100):
        self.frames = FrameCollection(loaded_frames)
        self.selection = UISelection()
        self.trace_half_window = trace_half_window
        self.click_pixel_tolerance = 5
        self.mouse_press_event = None
        self.updating_channel_box = False
        self.updating_checks = False
        self.axes_to_label = {}
        self.image_artists = []
        self.selectors = []
        self.view_buttons = []
        self.trace_plot = None
        self.colorbar = None

        self._create_figure()
        self.trace_plot = ProjectedTracePlot(
            self.figure,
            self.trace_axis,
            self.frames,
            self.trace_half_window,
        )
        self._build_controls()
        self._build_selectors()
        self._connect_events()
        self._initial_draw()

    def _create_figure(self):
        self.figure = plt.figure(figsize=(18, 10))
        self.grid = self.figure.add_gridspec(
            3,
            5,
            height_ratios=(0.12, 1.0, 0.75),
            width_ratios=(1.0, 1.0, 1.0, 0.16, 1.05),
        )
        self.button_axes = [
            self.figure.add_subplot(self.grid[0, index]) for index in range(3)
        ]
        self.reset_button_axis = self.figure.add_subplot(self.grid[0, 3])
        self.channel_box_axis = self.figure.add_subplot(self.grid[0, 4])
        self.image_axes = [
            self.figure.add_subplot(self.grid[1, index]) for index in range(3)
        ]
        self.colorbar_axis = self.figure.add_subplot(self.grid[1, 3])
        self.check_axis = self.figure.add_subplot(self.grid[1:, 4])
        self.trace_axis = self.figure.add_subplot(self.grid[2, :4])
        self.status_text = self.figure.text(
            0.01,
            0.01,
            "Select 1-3 files, switch U/V/W, click a bin for a trace, "
            "drag to zoom, or enter an absolute channel.",
        )
        self.figure.subplots_adjust(
            left=0.055,
            right=0.98,
            bottom=0.085,
            top=0.94,
            wspace=0.32,
            hspace=0.62,
        )

    def _build_controls(self):
        self._build_view_buttons()
        self._build_reset_button()
        self._build_channel_box()
        self._build_file_checks()
        self._update_view_button_style()

    def _build_view_buttons(self):
        for button_axis, view_name in zip(self.button_axes, VIEW_ORDER):
            button = Button(button_axis, view_name)
            button.on_clicked(
                lambda _event, selected_view=view_name: self._set_active_view(
                    selected_view
                )
            )
            self.view_buttons.append(button)

    def _build_reset_button(self):
        self.reset_button = Button(self.reset_button_axis, "reset")
        self.reset_button.on_clicked(lambda _event: self._reset_zoom())

    def _build_channel_box(self):
        self.channel_box = TextBox(self.channel_box_axis, "channel", initial="0")
        self.channel_box.on_submit(lambda _text: self._plot_from_channel_box())

    def _build_file_checks(self):
        check_states = [
            label in self.frames.active_labels for label in self.frames.labels
        ]
        self.file_checks = CheckButtons(self.check_axis, self.frames.labels, check_states)
        self.file_checks.on_clicked(self._on_file_check)
        self.check_axis.set_title("active files")

    def _build_selectors(self):
        for axis in self.image_axes:
            selector = RectangleSelector(
                axis,
                self._on_select,
                button=[1],
                minspanx=self.click_pixel_tolerance,
                minspany=self.click_pixel_tolerance,
                spancoords="pixels",
                interactive=False,
                useblit=True,
            )
            self.selectors.append(selector)

    def _connect_events(self):
        self.figure.canvas.mpl_connect("button_press_event", self._on_mouse_press)
        self.figure.canvas.mpl_connect("button_release_event", self._on_mouse_release)

    def _initial_draw(self):
        self._redraw_images()
        self._plot_channel(
            self.selection.active_view_name,
            channel=0,
            selected_tick=0,
        )

    def _active_view_shape(self, view_name=None):
        return self.frames.require_compatible_view_shape(
            view_name or self.selection.active_view_name
        )

    def _active_view_arrays(self):
        return self.frames.active_view_arrays(self.selection.active_view_name)

    def _display_limits(self):
        return finite_limits(list(self._active_view_arrays().values()))

    def _update_view_button_style(self):
        for button, view_name in zip(self.view_buttons, VIEW_ORDER):
            color = (
                "#d7ecff"
                if view_name == self.selection.active_view_name
                else "#f0f0f0"
            )
            button.ax.set_facecolor(color)

    def _on_file_check(self, clicked_label):
        if self.updating_checks:
            return
        statuses = self.file_checks.get_status()
        selected = [
            label for label, status in zip(self.frames.labels, statuses) if status
        ]

        try:
            self.frames.set_active_labels(selected, self.selection.active_view_name)
        except ValueError as error:
            self._revert_check(clicked_label)
            self.status_text.set_text(str(error))
            self.figure.canvas.draw_idle()
            return

        self._clamp_selected_bin()
        self._set_channel_box_to_current_channel()
        self._refresh_all_plots()
        self.status_text.set_text(
            f"Showing {len(self.frames.active_labels)} active file(s): "
            f"{', '.join(self.frames.active_labels)}."
        )
        self.figure.canvas.draw_idle()

    def _revert_check(self, label):
        try:
            index = self.frames.labels.index(label)
        except ValueError:
            return
        self.updating_checks = True
        try:
            self.file_checks.set_active(index)
        finally:
            self.updating_checks = False

    def _set_active_view(self, view_name):
        old_view = self.selection.active_view_name
        self.selection.active_view_name = view_name
        try:
            self._active_view_shape(view_name)
        except ValueError as error:
            self.selection.active_view_name = old_view
            self.status_text.set_text(str(error))
            self.figure.canvas.draw_idle()
            return

        self._clamp_selected_bin()
        self._set_channel_box_to_current_channel()
        self._refresh_all_plots()
        self._update_view_button_style()
        self.status_text.set_text(f"Showing {view_name} for active files.")
        self.figure.canvas.draw_idle()

    def _clamp_selected_bin(self):
        channel_count, tick_count = self._active_view_shape()
        self.selection.selected_channel = min(
            self.selection.selected_channel,
            channel_count - 1,
        )
        self.selection.selected_tick = min(self.selection.selected_tick, tick_count - 1)

    def _set_channel_box_to_current_channel(self):
        absolute_channel = (
            self.frames.view_offset(
                self.frames.active_labels[0],
                self.selection.active_view_name,
            )
            + self.selection.selected_channel
        )
        self._set_channel_box_value(absolute_channel)

    def _set_channel_box_value(self, absolute_channel):
        self.updating_channel_box = True
        try:
            self.channel_box.set_val(str(int(absolute_channel)))
        finally:
            self.updating_channel_box = False

    def _refresh_all_plots(self):
        self._redraw_images()
        self._plot_channel(
            self.selection.active_view_name,
            self.selection.selected_channel,
            selected_tick=self.selection.selected_tick,
        )

    def _redraw_images(self):
        self.axes_to_label.clear()
        self.image_artists.clear()
        vmin, vmax = self._display_limits()

        for axis in self.image_axes:
            axis.clear()
            axis.set_visible(False)
        self.colorbar_axis.clear()

        image = None
        for axis, loaded in zip(self.image_axes, self.frames.active_frames()):
            image = self._draw_image_axis(axis, loaded, vmin, vmax)

        if image is not None:
            self.colorbar = self.figure.colorbar(
                image,
                cax=self.colorbar_axis,
                label="value",
            )
        self._apply_zoom_limits()

    def _draw_image_axis(self, axis, loaded, vmin, vmax):
        array = self.frames.view_array(loaded, self.selection.active_view_name)
        image = axis.imshow(
            array,
            aspect="auto",
            origin="lower",
            cmap="plasma",
            vmin=vmin,
            vmax=vmax,
        )
        self.image_artists.append(image)
        axis.set_visible(True)
        axis.set_title(f"{loaded.spec.label} {self.selection.active_view_name}")
        axis.set_xlabel(loaded.spec.xaxis_label)
        axis.set_ylabel(loaded.spec.yaxis_label)
        self.axes_to_label[axis] = loaded.spec.label
        return image

    def _selection_from_events(self, press_event, release_event):
        coordinates = (
            press_event.xdata,
            press_event.ydata,
            release_event.xdata,
            release_event.ydata,
        )
        if any(coordinate is None for coordinate in coordinates):
            return None

        label = self.axes_to_label.get(press_event.inaxes)
        if label is None:
            return None
        array = self.frames.view_array(
            self.frames.frames_by_label[label],
            self.selection.active_view_name,
        )
        tick_min = int(np.floor(min(press_event.xdata, release_event.xdata)))
        tick_max = int(np.ceil(max(press_event.xdata, release_event.xdata)))
        channel_min = int(np.floor(min(press_event.ydata, release_event.ydata)))
        channel_max = int(np.ceil(max(press_event.ydata, release_event.ydata)))

        tick_min = int(np.clip(tick_min, 0, array.shape[1] - 1))
        tick_max = int(np.clip(tick_max, 0, array.shape[1] - 1))
        channel_min = int(np.clip(channel_min, 0, array.shape[0] - 1))
        channel_max = int(np.clip(channel_max, 0, array.shape[0] - 1))
        if tick_min == tick_max or channel_min == channel_max:
            return None
        return tick_min, tick_max, channel_min, channel_max

    def _is_click_sized_motion(self, press_event, release_event):
        dx = release_event.x - press_event.x
        dy = release_event.y - press_event.y
        return (dx * dx + dy * dy) <= self.click_pixel_tolerance**2

    def _on_select(self, press_event, release_event):
        if self._is_click_sized_motion(press_event, release_event):
            return
        selection = self._selection_from_events(press_event, release_event)
        if selection is None:
            return

        self.selection.zoom_limits = selection
        self._apply_zoom_limits()
        tick_min, tick_max, channel_min, channel_max = selection
        self.status_text.set_text(
            f"Selected ticks {tick_min}-{tick_max}, "
            f"channels {channel_min}-{channel_max} within "
            f"{self.selection.active_view_name}."
        )
        self.figure.canvas.draw_idle()

    def _limits_for_active_view(self):
        channel_count, tick_count = self._active_view_shape()
        if self.selection.zoom_limits is None:
            return 0, tick_count - 1, 0, channel_count - 1

        tick_min, tick_max, channel_min, channel_max = self.selection.zoom_limits
        return (
            int(np.clip(tick_min, 0, tick_count - 1)),
            int(np.clip(tick_max, 0, tick_count - 1)),
            int(np.clip(channel_min, 0, channel_count - 1)),
            int(np.clip(channel_max, 0, channel_count - 1)),
        )

    def _apply_zoom_limits(self):
        tick_min, tick_max, channel_min, channel_max = self._limits_for_active_view()
        for axis in self.image_axes:
            if axis.get_visible():
                axis.set_xlim(tick_min - 0.5, tick_max + 0.5)
                axis.set_ylim(channel_min - 0.5, channel_max + 0.5)

    def _reset_zoom(self):
        self.selection.zoom_limits = None
        self._apply_zoom_limits()
        self.status_text.set_text(
            f"Selection reset; showing the full "
            f"{self.selection.active_view_name} range."
        )
        self.figure.canvas.draw_idle()

    def _event_bin(self, event):
        label = self.axes_to_label.get(event.inaxes)
        if label is None or event.xdata is None or event.ydata is None:
            return None

        array = self.frames.view_array(
            self.frames.frames_by_label[label],
            self.selection.active_view_name,
        )
        tick = int(round(event.xdata))
        channel = int(round(event.ydata))
        if not (0 <= tick < array.shape[1] and 0 <= channel < array.shape[0]):
            return None
        return self.selection.active_view_name, label, channel, tick, array[channel, tick]

    def _on_mouse_press(self, event):
        if event.button != 1 or event.inaxes not in self.image_axes:
            self.mouse_press_event = None
            return
        self.mouse_press_event = event

    def _on_mouse_release(self, event):
        if event.button != 1 or self.mouse_press_event is None:
            self.mouse_press_event = None
            return

        press_event = self.mouse_press_event
        self.mouse_press_event = None
        if event.inaxes != press_event.inaxes or event.inaxes not in self.image_axes:
            return
        if not self._is_click_sized_motion(press_event, event):
            return

        selected = self._event_bin(event)
        if selected is None:
            return

        view_name, label, channel, tick, pixel_value = selected
        self.selection.selected_channel = channel
        self.selection.selected_tick = tick
        self._set_channel_box_value(self.frames.view_offset(label, view_name) + channel)
        self._plot_selected_bin(view_name, label, channel, tick, pixel_value)

    def _plot_selected_bin(self, view_name, label, channel, tick, pixel_value):
        tick_min, tick_max = self._trace_tick_limits(view_name, tick)
        absolute_channel = self.frames.view_offset(label, view_name) + channel
        self.status_text.set_text(
            f"Selected absolute channel {absolute_channel} "
            f"({view_name} channel {channel}); showing ticks {tick_min}-{tick_max} "
            f"around tick {tick}; {label} value={pixel_value:.6g}."
        )
        self._plot_channel(view_name, channel, selected_tick=tick)
        self.figure.canvas.draw_idle()

    def _plot_from_channel_box(self):
        if self.updating_channel_box:
            return
        try:
            absolute_channel = int(self.channel_box.text.strip())
        except ValueError:
            self.status_text.set_text("Channel entry must be an integer.")
            self.figure.canvas.draw_idle()
            return

        view_lookup = self.frames.view_for_absolute_channel(absolute_channel)
        if view_lookup is None:
            self._show_invalid_channel_status(absolute_channel)
            return

        view_name, channel = view_lookup
        if view_name != self.selection.active_view_name:
            self._set_active_view(view_name)

        _, _, channel_min, channel_max = self._limits_for_active_view()
        if not channel_min <= channel <= channel_max:
            self.status_text.set_text(
                f"Channel {absolute_channel} is outside the current selection."
            )
            self.figure.canvas.draw_idle()
            return

        self.selection.selected_channel = channel
        self._plot_channel(view_name, channel, selected_tick=self.selection.selected_tick)
        tick_min, tick_max = self._trace_tick_limits(
            view_name,
            self.selection.selected_tick,
        )
        self.status_text.set_text(
            f"Selected absolute channel {absolute_channel} "
            f"({view_name} channel {channel}); showing ticks {tick_min}-{tick_max} "
            f"around tick {self.selection.selected_tick}."
        )
        self.figure.canvas.draw_idle()

    def _show_invalid_channel_status(self, absolute_channel):
        loaded = self.frames.first_active_frame()
        valid_ranges = ", ".join(
            f"{name}={start}:{stop}"
            for name, (start, stop) in loaded.spec.channel_mapping.items()
        )
        self.status_text.set_text(
            f"Channel {absolute_channel} is outside active ranges: {valid_ranges}."
        )
        self.figure.canvas.draw_idle()

    def _trace_tick_limits(self, view_name, selected_tick):
        return self.trace_plot.trace_tick_limits(view_name, selected_tick)

    def _plot_channel(self, view_name, channel, selected_tick=None):
        self.selection.selected_channel = channel
        if selected_tick is None:
            selected_tick = self.selection.selected_tick
        self.selection.selected_tick = selected_tick
        self.trace_plot.plot(view_name, channel, selected_tick)

    def show(self):
        plt.show()
