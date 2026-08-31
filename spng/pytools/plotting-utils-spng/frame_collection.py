"""Loaded-frame collection and view-slicing logic."""

from frame_models import MAX_ACTIVE_FILES, VIEW_ORDER


class FrameCollection:
    """Loaded-frame data model used by the UI."""

    def __init__(self, loaded_frames):
        if not loaded_frames:
            raise ValueError("At least one loaded frame is required")
        self.loaded_frames = loaded_frames
        self.labels = [loaded.spec.label for loaded in loaded_frames]
        self.frames_by_label = {loaded.spec.label: loaded for loaded in loaded_frames}
        self.active_labels = self.labels[:MAX_ACTIVE_FILES]

    def active_frames(self):
        return [self.frames_by_label[label] for label in self.active_labels]

    def set_active_labels(self, labels, view_name):
        if not labels:
            raise ValueError("At least one file must remain selected.")
        if len(labels) > MAX_ACTIVE_FILES:
            raise ValueError("Select at most 3 files at once.")
        old_labels = self.active_labels
        self.active_labels = list(labels)
        try:
            self.require_compatible_view_shape(view_name)
        except Exception:
            self.active_labels = old_labels
            raise

    def view_slice(self, loaded, view_name):
        start, stop = loaded.spec.channel_mapping[view_name]
        if stop > loaded.frame.shape[0]:
            raise ValueError(
                f"{loaded.spec.label} view {view_name} uses channels {start}:{stop}, "
                f"but frame has only {loaded.frame.shape[0]} channels"
            )
        return slice(start, stop), start

    def view_array(self, loaded, view_name):
        channel_slice, _ = self.view_slice(loaded, view_name)
        return loaded.frame[channel_slice, :]

    def active_view_arrays(self, view_name):
        return {
            loaded.spec.label: self.view_array(loaded, view_name)
            for loaded in self.active_frames()
        }

    def require_compatible_view_shape(self, view_name):
        arrays = [self.view_array(loaded, view_name) for loaded in self.active_frames()]
        shapes = {array.shape for array in arrays}
        if len(shapes) != 1:
            readable = ", ".join(
                f"{loaded.spec.label}: {self.view_array(loaded, view_name).shape}"
                for loaded in self.active_frames()
            )
            raise ValueError(
                f"Active files have incompatible {view_name} shapes: {readable}"
            )
        return arrays[0].shape

    def view_offset(self, label, view_name):
        loaded = self.frames_by_label[label]
        _, offset = self.view_slice(loaded, view_name)
        return offset

    def first_active_frame(self):
        return self.frames_by_label[self.active_labels[0]]

    def view_for_absolute_channel(self, absolute_channel):
        loaded = self.first_active_frame()
        for view_name in VIEW_ORDER:
            start, stop = loaded.spec.channel_mapping[view_name]
            if start <= absolute_channel < stop:
                return view_name, absolute_channel - start
        return None
