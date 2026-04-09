# ClusteringLiveDead — Underlying Function Audit

Functions called by `clustering_live_dead` in `clus/src/clustering_live_dead.cxx`.

---

## Underlying Functions

| Function | Source | Called where |
|---|---|---|
| `sort_clusters` | `Facade_Cluster.cxx` | lines 100, 103 (live + dead cluster ordering) |
| `merge_clusters` | `ClusteringFuncs.cxx` | line 345 |

`Find_Closest_Points` is **not** called by `clustering_live_dead`.

---

## Efficiency

No issues found.

- `sort_clusters` uses `std::sort` with `cluster_less` — O(N log N), already optimal.
- `merge_clusters` uses `boost::connected_components` + `std::map<int, std::set<int>>` —
  integer keys, already O(N log N).

---

## Randomness

No issues found.

- `sort_clusters` is deterministic (see `cluster_less` analysis in
  `extend-underlying.md`).
- `merge_clusters` is deterministic given deterministically-ordered vertex insertion
  (guaranteed by `sort_clusters` on `live_clusters` before graph construction).

---

## Notes

- Dead and live cluster vectors are each sorted with `sort_clusters` before use — already
  applied in prior session. Ordering of both is therefore deterministic.
