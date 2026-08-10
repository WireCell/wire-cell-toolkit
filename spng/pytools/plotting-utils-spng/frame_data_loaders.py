"""Data loader classes for NPZ and HDF5 frame inputs."""

from abc import ABC, abstractmethod

import numpy as np

from frame_models import LoadedFrame, frame_sort_key

try:
    import h5py
except ImportError:
    h5py = None


class DataLoader(ABC):
    """Base class for loading one configured dataset into a 2D frame."""

    def __init__(self, spec, frame_index):
        self.spec = spec
        self.frame_index = frame_index

    def load(self):
        if not self.spec.path.exists():
            raise FileNotFoundError(
                f"{self.spec.label} file does not exist: {self.spec.path}"
            )
        source_dataset, array = self.load_array()
        return LoadedFrame(
            self.spec,
            source_dataset,
            self.as_2d_frame(array),
        )

    @abstractmethod
    def load_array(self):
        """Return (source_dataset_name, array)."""

    def as_2d_frame(self, array):
        array = np.asarray(array)
        array = np.squeeze(array)
        if array.ndim == 2:
            return array
        if array.ndim != 1:
            raise ValueError(
                f"{self.spec.label} must be a 1D or 2D dataset after squeeze, "
                f"got {array.shape}"
            )

        max_channel = max(stop for _, stop in self.spec.channel_mapping.values())
        tick_count = self.spec.tick_info
        if tick_count is None:
            if array.size % max_channel != 0:
                raise ValueError(
                    f"{self.spec.label} is 1D with {array.size} values. "
                    "Provide tick_info so it can be reshaped into channel x tick."
                )
            tick_count = array.size // max_channel

        expected_size = max_channel * tick_count
        if array.size != expected_size:
            raise ValueError(
                f"{self.spec.label} cannot reshape {array.size} values into "
                f"({max_channel}, {tick_count})"
            )
        return array.reshape(max_channel, tick_count)


class NpzDataLoader(DataLoader):
    """Load frames from NPZ files."""

    def load_array(self):
        with np.load(self.spec.path) as npz_file:
            if self.spec.dataset:
                return self._load_named_dataset(npz_file)
            return self._load_frame_index_dataset(npz_file)

    def _load_named_dataset(self, npz_file):
        if self.spec.dataset not in npz_file.files:
            available = ", ".join(npz_file.files)
            raise ValueError(
                f"Dataset {self.spec.dataset!r} not found in {self.spec.path}; "
                f"available keys: {available}"
            )
        return self.spec.dataset, np.asarray(npz_file[self.spec.dataset])

    def _load_frame_index_dataset(self, npz_file):
        frame_names = sorted(
            (name for name in npz_file.files if name.startswith("frame_")),
            key=frame_sort_key,
        )
        if not frame_names:
            raise ValueError(f"No frame_* datasets found in {self.spec.path}")
        try:
            frame_name = frame_names[self.frame_index]
        except IndexError as error:
            raise ValueError(
                f"Frame index {self.frame_index} is unavailable in {self.spec.path}; "
                f"the file contains {len(frame_names)} frame(s)"
            ) from error
        return frame_name, np.asarray(npz_file[frame_name])


class Hdf5DataLoader(DataLoader):
    """Load frames from HDF5 files."""

    def load_array(self):
        if h5py is None:
            raise ImportError(
                "HDF5 input requires h5py. Install it with "
                "`python3 -m pip install h5py` or use only NPZ inputs."
            )
        if not self.spec.dataset:
            raise ValueError(f"{self.spec.label} is HDF5 and requires a dset value")

        with h5py.File(self.spec.path, "r") as h5_file:
            if self._is_exact_dataset_path(h5_file):
                return self.spec.dataset, np.asarray(h5_file[self.spec.dataset])
            return self._load_by_basename(h5_file)

    def _is_exact_dataset_path(self, h5_file):
        return (
            self.spec.dataset in h5_file
            and isinstance(h5_file[self.spec.dataset], h5py.Dataset)
        )

    def _load_by_basename(self, h5_file):
        dataset_paths = []

        def collect_dataset(name, obj):
            if isinstance(obj, h5py.Dataset) and name.split("/")[-1] == self.spec.dataset:
                dataset_paths.append(name)

        h5_file.visititems(collect_dataset)
        dataset_paths.sort(key=frame_sort_key)
        if not dataset_paths:
            raise ValueError(f"No {self.spec.dataset} dataset found in {self.spec.path}")

        try:
            dataset_path = dataset_paths[self.frame_index]
        except IndexError as error:
            raise ValueError(
                f"Frame index {self.frame_index} is unavailable in {self.spec.path}; "
                f"the file contains {len(dataset_paths)} {self.spec.dataset} dataset(s)"
            ) from error
        return dataset_path, np.asarray(h5_file[dataset_path])


class DataLoaderFactory:
    """Create the correct loader for a FileSpec."""

    @staticmethod
    def create(spec, frame_index):
        if spec.file_format == "npz":
            return NpzDataLoader(spec, frame_index)
        if spec.file_format in ("h5", "hdf5"):
            return Hdf5DataLoader(spec, frame_index)
        raise ValueError(f"Unsupported file format: {spec.file_format}")
