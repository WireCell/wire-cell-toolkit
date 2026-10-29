#!/usr/bin/env python3
"""pr_ml doc 33 perf round 1: the vectorized steiner_truth.py paths against the legacy forms they replace, on synthetic tracks --
OwnerIndex._eval (packed-key unique vs np.unique(..., axis=0)), true_vertices (one batched query per track + bounding-box cut vs the
per-pair loop), the bounded 2 cm scope query and batched k-d queries vs per-track queries.  Exact equality required; prints PASS/FAIL."""
import os, sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import steiner_truth as S

def legacy_eval(self, Q, qi, oi, D, Sa, T, O):
    pr = np.unique(np.column_stack([oi, qi]), axis=0)
    st = np.r_[0, np.flatnonzero(np.diff(pr[:, 0])) + 1]
    en = np.r_[st[1:], len(pr)]
    for s0, e0 in zip(st, en):
        o, qs = int(pr[s0, 0]), pr[s0:e0, 1]
        d, s, t = self.polys[o].nearest(Q[qs])
        better = d < D[qs]
        q = qs[better]
        D[q], Sa[q], T[q], O[q] = d[better], s[better], t[better], o

def legacy_true_vertices(trajs, lmin, dv):
    keys = [k for k in sorted(trajs) if trajs[k].length >= lmin]
    out = []
    for k in keys:
        c = trajs[k].start[None]
        for k2 in keys:
            if k2 != k and trajs[k2].nearest(c)[0][0] <= dv:
                out.append((k, trajs[k].start)); break
    return out

def tracks(rng, n):
    T = {}
    for i in range(n):
        p0 = rng.uniform(-300, 300, 3)
        if i and i % 3 == 0:
            p0 = T[i].units[0].V[-1].copy()                          # shares the previous track's end
        elif i and i % 5 == 0:                                       # a near miss: 7 mm off the previous track's end (dv = 10)
            u = rng.normal(0, 1, 3); p0 = T[i].units[0].V[-1] + 7.0 * u / np.linalg.norm(u)
        steps = rng.normal(0, 4, (int(rng.integers(1, 40)), 3)) + rng.normal(0, 3, 3)
        V = np.vstack([p0, p0 + np.cumsum(steps, axis=0)])
        T[i + 1] = S.KeyTraj(i + 1, [S.Polyline(V)], [i + 1], [0.0])
    return T

def main():
    rng = np.random.default_rng(20261028); ok = True
    T = tracks(rng, 120)
    a, b = S.true_vertices(T, 10.0, 10.0), legacy_true_vertices(T, 10.0, 10.0)
    same = len(a) == len(b) and all(x[0] == y[0] and np.array_equal(x[1], y[1]) for x, y in zip(a, b))
    print("true_vertices", len(a), "PASS" if same and len(a) > 0 else "FAIL"); ok &= same and len(a) > 0
    oi = S.OwnerIndex([t.units[0] for t in T.values()])
    Q = rng.uniform(-320, 320, (5000, 3))
    new = oi.nearest(Q, cap=80.0)
    S.OwnerIndex._eval, keep = legacy_eval, S.OwnerIndex._eval
    old = oi.nearest(Q, cap=80.0)
    S.OwnerIndex._eval = keep
    same = all(np.array_equal(x, y) for x, y in zip(new, old)) and np.isfinite(new[0]).sum() > 100
    print("OwnerIndex.nearest", int(np.isfinite(new[0]).sum()), "PASS" if same else "FAIL"); ok &= same
    P = rng.uniform(-300, 300, (4000, 3)); kd = cKDTree(P); X = rng.uniform(-310, 310, (20000, 3))
    full = kd.query(X)[0]; bnd = kd.query(X, distance_upper_bound=20.0 * (1 + 1e-9))[0]
    same = np.array_equal(full <= 20.0, bnd <= 20.0) and np.array_equal(full[full <= 20.0], bnd[bnd <= 20.0])
    print("bounded scope query", int((full <= 20.0).sum()), "PASS" if same else "FAIL"); ok &= same
    parts = np.array_split(X, 7); same = np.array_equal(np.concatenate([kd.query(p)[0] for p in parts]), full)
    print("batched k-d query", "PASS" if same else "FAIL"); ok &= same
    return ok

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
