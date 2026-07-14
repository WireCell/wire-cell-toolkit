# ClusteringExtendLoop Analysis

Comparison of `ClusteringExtendLoop::visit` in `clus/src/clustering_extend.cxx`
against the caller loop in `prototype_base/2dtoy/src/ToyClustering.cxx`
(lines 293–323).

---

## Functional Equivalence

The toolkit `ClusteringExtendLoop::visit` is identical to the prototype loop:

| Step | Prototype | Toolkit |
|------|-----------|---------|
| Busy-event guard | `if (size > 1100) num_try = 1` | `if (nchildren() > 1100) num_try = 1` |
| flag=1 (prolong) | `Clustering_extend(..., 1, 150 cm, 0)` | `clustering_extend(..., 1, 150 cm, 0)` |
| flag=2 (parallel) | `Clustering_extend(..., 2, 30 cm, 0)` | `clustering_extend(..., 2, 30 cm, 0)` |
| flag=3 (regular) | `Clustering_extend(..., 3, 15 cm, 0)` | `clustering_extend(..., 3, 15 cm, 0)` |
| flag=4 i=0 | `Clustering_extend(..., 4, 60 cm, 0)` | `clustering_extend(..., 4, 60 cm, 0)` |
| flag=4 i>0 | `Clustering_extend(..., 4, 35 cm, i)` | `clustering_extend(..., 4, 35 cm, i)` |

The `num_try` default (0 in toolkit, 3 hardcoded in prototype) is a configuration
difference: the pipeline JSON sets `num_try: 3` explicitly, so runtime behavior
is identical.

---

## Bugs Found

None.

---

## Efficiency and Randomness

No changes needed. `ClusteringExtendLoop::visit` is a thin wrapper that only
calls `clustering_extend` in a loop — no containers of its own. All
efficiency/randomness improvements to `clustering_extend` were already applied
in a prior session (sorted `live_clusters`, `unordered_map` for
`map_cluster_index`, `unordered_set` for `used_clusters`).

---

## Notes

- The prototype also calls `Clustering_extend(..., 4, 60 cm, 0, 15 cm, 1)` once
  *before* the main loop (line 264). This pre-loop call is handled separately in
  the toolkit pipeline as a standalone `ClusteringExtend` visitor with `flag=4`,
  `length_cut=60 cm`, `num_dead_try=1` — not part of `ClusteringExtendLoop`.
