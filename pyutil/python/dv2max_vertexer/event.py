"""Event container (vendored from sbndvtx/event.py, inference subset): x-fixed reco inputs, Wire-Cell reco rows.

All reco coordinates in an Event are float32, reco frame.  Truth is never filled at inference (event.truth is
None); `strip_truth()` is kept because the vendored feature / pool code calls it.  Dropped from the original:
load_event, meta_dict / from_meta (cache and truth I/O) and split_segments.
"""
import dataclasses
from dataclasses import dataclass, field

import numpy as np

def _empty3():
    return np.zeros((0, 3), np.float32)


@dataclass
class Event:
    key: str                                   # "{sample}:{run}:{subrun}:{event}"
    sample: str                                # 'nue' | 'bnb'
    index: int                                 # row in the r3 file
    rse: np.ndarray                            # (3,) int64
    fold: int
    # clustering_pr: whole-event 3D imaging points (cosmics included)
    points_xyz: np.ndarray = field(default_factory=_empty3)
    points_q: np.ndarray = field(default_factory=lambda: np.zeros(0, np.float32))
    points_cluster_id: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    # track_fit: ordered trajectory-fit points
    fit_xyz: np.ndarray = field(default_factory=_empty3)
    # q is Bee-encoded: q = 0.1*dQ - 1000 clipped at 0 -> dQ = 10*(q + 1000) e; use features.fitgraph.decode_dq
    fit_q: np.ndarray = field(default_factory=lambda: np.zeros(0, np.float32))
    fit_cluster_id: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    fit_segment_id: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))  # PR segment, -1 = vertex pt
    fit_index: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))   # index along segment (-1)
    # shower_track: q is a tag (15000 on shower-tagged points, 0 on track-tagged points)
    shower_xyz: np.ndarray = field(default_factory=_empty3)
    shower_q: np.ndarray = field(default_factory=lambda: np.zeros(0, np.float32))
    shower_cluster_id: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    # Wire-Cell PR vertex candidates and the production vertex
    cand_xyz: np.ndarray = field(default_factory=_empty3)
    cand_q: np.ndarray = field(default_factory=lambda: np.zeros(0, np.float32))
    tag_nu: np.ndarray = field(default_factory=lambda: np.full(3, np.nan, np.float32))
    wc_reco: list = field(default_factory=list)        # list[WCRecoNu], x-fixed (INPUT)
    wc_status: str = "missing"                          # 'ok' | 'none' | 'missing'
    truth: object = None                                # LABELS ONLY (always None at inference)
    features: dict = field(default_factory=dict)        # name -> array, see sbndvtx.features

    # -- convenience ----------------------------------------------------------
    @property
    def has_tag_nu(self):
        return bool(np.all(np.isfinite(self.tag_nu)))

    @property
    def has_target(self):
        return self.truth is not None and self.truth.has_target

    @property
    def target_xyz(self):
        return self.truth.target_xyz if self.truth is not None else np.full(3, np.nan, np.float32)

    @property
    def n_points(self):
        return len(self.points_xyz)

    @property
    def n_fit(self):
        return len(self.fit_xyz)

    @property
    def n_shower(self):
        return len(self.shower_xyz)

    @property
    def n_cand(self):
        return len(self.cand_xyz)

    def level_size(self, level):
        return {"point": self.n_points, "fit": self.n_fit, "shower": self.n_shower,
                "cand": self.n_cand, "event": 1}[level]

    def strip_truth(self):
        """Shallow copy without MC truth (what feature code and models receive)."""
        return dataclasses.replace(self, truth=None, features=dict(self.features))
