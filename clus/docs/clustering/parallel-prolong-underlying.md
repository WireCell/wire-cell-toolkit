# ClusteringParallelProlong — Underlying Function Audit

Functions called by `clustering_parallel_prolong` in
`clus/src/clustering_parallel_prolong.cxx`.

---

## Underlying Functions

| Function | Source | Called where |
|---|---|---|
| `Find_Closest_Points` | `clustering_extend.cxx` | line 73 (in `Clustering_2nd_round`) |
| `get_strategic_points` | `clustering_extend.cxx` | inside `Find_Closest_Points` |
| `compute_wireplane_params` | `Clustering_Util.cxx` | line 305 |
| `sort_clusters` | `Facade_Cluster.cxx` | line 316 |
| `merge_clusters` | `ClusteringFuncs.cxx` | line 347 |

---

## Randomness

### `Find_Closest_Points` / `get_strategic_points`

The randomness fix applied to `get_strategic_points` (see `extend-underlying.md`)
directly benefits `clustering_parallel_prolong`.

### `compute_wireplane_params` / wire-direction maps

`wpid_U_dir` / `wpid_V_dir` / `wpid_W_dir` are `std::map<WirePlaneId, ...>` — keyed by
`WirePlaneId` struct (APA/face/layer integers). Ordering is deterministic
hardware-topology order. No pointer dependence. No issues.

### `sort_clusters` / `merge_clusters`

See `extend-underlying.md`. Both deterministic. No issues.

---

## Efficiency

No additional improvements found.

- `wpid_params` and `apas` (output params of `compute_wireplane_params`) are filled but
  never read in this file — they are dead locals. Cannot be removed without changing the
  shared helper signature. Noted but not changed (same situation as
  `clustering_regular`).

---

## Notes

- `ilive2desc` was already `std::unordered_map<int, int>`.
