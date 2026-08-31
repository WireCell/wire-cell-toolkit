# ClusteringRegular — Underlying Function Audit

Functions called by `clustering_regular` in `clus/src/clustering_regular.cxx`.

---

## Underlying Functions

| Function | Source | Called where |
|---|---|---|
| `Find_Closest_Points` | `clustering_extend.cxx` | line 75 (in `Clustering_2nd_round`) |
| `get_strategic_points` | `clustering_extend.cxx` | inside `Find_Closest_Points` |
| `compute_wireplane_params` | `Clustering_Util.cxx` | wire-angle setup block |
| `sort_clusters` | `Facade_Cluster.cxx` | line 473 |
| `merge_clusters` | `ClusteringFuncs.cxx` | line 504 |

---

## Randomness

### `Find_Closest_Points` / `get_strategic_points`

The randomness fix applied to `get_strategic_points` (replacing pointer-ordered
`std::set<std::pair<geo_point_t, const Blob*>>` with a sorted+deduped `std::vector`)
directly benefits `clustering_regular` since it calls `Find_Closest_Points`.
See `extend-underlying.md` for details.

### `compute_wireplane_params`

Iterates over `wpids` which is `std::set<WirePlaneId>` — keyed by `WirePlaneId` struct
(APA/face/layer integers), deterministic hardware-topology ordering. No pointer
dependence. No issues.

### `sort_clusters` / `merge_clusters`

See `extend-underlying.md`. Both deterministic. No issues.

---

## Efficiency

No additional improvements found.

- `compute_wireplane_params` fills `wpid_params` and `apas` which are declared in
  `clustering_regular` but never read after the call (dead within this file). These
  cannot be removed without changing the shared helper's signature — other callers may
  read them. Noted but not changed (same as `clustering_parallel_prolong`).

---

## Notes

- `wpid_U_dir` / `wpid_V_dir` / `wpid_W_dir` maps remain `std::map<WirePlaneId, ...>` —
  keyed by `WirePlaneId` struct, deterministic. No change needed.
- `ilive2desc` was already `std::unordered_map<int, int>`.
