"""Detector geometry in the x-fixed reco frame (measured extents of clustering_pr on the dev cache).

SBND: two TPCs drifting along x, cathode at x=0, anodes at |x|=X_ANODE.
Measured on 800 dev events: x in [-201.18, 201.18], y in [-199.96, 199.96], z in [0.15, 501.0].
"""
import numpy as np

X_ANODE = 201.2
Y_HALF = 200.0
Z_LO, Z_HI = 0.0, 501.0
BEAM_DIR = np.array([0.0, 0.0, 1.0])


def wall_distances(xyz):
    """(n,3) -> (n,4) float32 distances (cm, >=0): anode plane, cathode, nearest y face, nearest z face."""
    xyz = np.asarray(xyz, np.float64).reshape(-1, 3)
    ax = np.abs(xyz[:, 0])
    d = np.stack([X_ANODE - ax, ax, Y_HALF - np.abs(xyz[:, 1]),
                  np.minimum(xyz[:, 2] - Z_LO, Z_HI - xyz[:, 2])], 1)
    return np.maximum(d, 0.0).astype(np.float32)


def boundary_distance(xyz):
    """(n,3) -> (n,) distance to the nearest TPC boundary a track can cross (anode, y, z faces; not the cathode)."""
    d = wall_distances(xyz)
    return np.min(d[:, [0, 2, 3]], axis=1)
