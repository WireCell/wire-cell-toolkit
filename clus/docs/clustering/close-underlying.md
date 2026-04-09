# ClusteringClose — Underlying Function Audit

Functions called by `clustering_close` / `Clustering_3rd_round` in
`clus/src/clustering_close.cxx`.

---

## Underlying Functions

| Function | Source | Called where |
|---|---|---|
| `Find_Closest_Points` | `clustering_extend.cxx` | line 64 (in `Clustering_3rd_round`) |
| `get_strategic_points` | `clustering_extend.cxx` | inside `Find_Closest_Points` |
| `sort_clusters` | `Facade_Cluster.cxx` | line 180 |
| `merge_clusters` | `ClusteringFuncs.cxx` | line 225 |

---

## Randomness

### `Find_Closest_Points` / `get_strategic_points`

The randomness fix applied to `get_strategic_points` (see `extend-underlying.md`)
directly benefits `clustering_close`.

### `sort_clusters` / `merge_clusters`

See `extend-underlying.md`. Both deterministic. No issues.

---

## Efficiency

No additional improvements found.

- `Clustering_3rd_round` also calls `cluster1.vhough_transform`, `cluster1.ndipole`,
  `cluster1.nnearby`, `cluster1.calc_ave_pos` — these are Cluster facade methods backed
  by KD-tree and Hough transforms. No randomness or efficiency issues found at this
  level; they operate on the cluster's internal point cloud directly.

---

## Notes

- `ilive2desc` was already `std::unordered_map<int, int>`.
- `used_clusters` and `map_cluster_index` were converted to `unordered_set` /
  `unordered_map` in the prior session (see `close-analysis.md`).
