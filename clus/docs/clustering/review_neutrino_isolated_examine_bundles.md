# Code Review: `Clustering_neutrino` / `Clustering_isolated` / `ExamineBundles`

**Date:** 2026-04-10  
**Reviewer:** Claude Code  
**Scope:** Three toolkit clustering stages versus their prototype counterparts, focusing
on (1) functional parity, (2) bugs, (3) algorithmic efficiency, (4) nondeterminism from
pointer-keyed containers, (5) callees, and (6) multi-APA/face correctness.

---

## 1. Overview and Methodology

### Toolkit source files reviewed
| Function | File | Lines |
|---|---|---|
| `Clustering_neutrino` | `clus/src/clustering_neutrino.cxx` | 56-993 |
| `Clustering_isolated` | `clus/src/clustering_isolated.cxx` | 59-540 |
| `ExamineBundles` | `clus/src/clustering_examine_bundles.cxx` | 71-183 |

### Prototype originals
| Function | File | Lines |
|---|---|---|
| `Clustering_neutrino` | `prototype_base/2dtoy/src/ToyClustering_neutrino.h` | 46-1006 |
| `Clustering_isolated` | `prototype_base/2dtoy/src/ToyClustering_isolated.h` | 3-326 |
| `ExamineBundles` / `ExamineBundle` | `prototype_base/2dtoy/src/ExamineBundles.cxx` | 5-123 |

### Terminology
- **wpid** — `WirePlaneId`: identifies a specific wire plane (APA, face, plane-layer). `WirePlaneId(0)` is the synthetic "all planes" sentinel used to query global detector metadata.
- **wpid_all** — `WirePlaneId wpid_all(0)`, the global-detector-envelope sentinel.
- **cluster_less_functor** — content-based `Cluster*` comparator (`Facade_Cluster.cxx:2498`) that orders by (length, nchildren, npoints, per-wpid range), providing deterministic pointer-keyed containers.

---

## 2. `Clustering_neutrino`

### 2.1 Function map (prototype → toolkit)

| Prototype block | Prototype lines | Toolkit lines | Notes |
|---|---|---|---|
| Geometry/plane params | 47-54 | 63-128 | Prototype: single `TPCParams` singleton. Toolkit: per-wpid loop. |
| Sort live clusters by length | 60-70 | 95-99 | Both sort descending. |
| Containment box / FV tests | 73-189 | 158-276 | See §2.3 for a bug at toolkit line 178. |
| Closest-cluster map | 191-217 | 279-310 | Logic identical. |
| `to_be_merged_pairs`, `used_clusters`, cloud map | 222-229 | 314-321 | Identical. |
| Pre-PCA warmup loop | 231-234 | 323-326 (commented) | Removed; toolkit uses lazy PCA. |
| Candidate × contained main loop | 239-925 | 331-953 | Heart of function; logic matches line-for-line (see §2.3 for one intentional divergence). |
| Separate_2 + extended cloud (cluster1) | 304-469 | 383-529 | API change: blob-id-map instead of clone/delete. |
| Separate_2 + extended cloud (cluster2) | 510-669 | 566-711 | Same. |
| merge_type 1/2/3 decisions | 676-875 | 718-916 | Identical. |
| `used_clusters` bookkeeping | 881-917 | 918-945 | One intentional divergence; see §2.3. |
| Final merge — graph | 929-1003 | 959-993 | Prototype: hand-written union-find. Toolkit: Boost `connected_components` via `merge_clusters()`. |

### 2.2 Logic parity with prototype

The port is structurally faithful. All merge-type decision cascades, threshold values, and
loop structure map 1:1. The `Separate_2` API change (returning a blob-to-group-id `vector<int>`
instead of cloning new clusters) is transparent to callers in `clustering_neutrino.cxx` because
the extended-point-cloud fill loop is identical — it just receives `b2id` instead of raw pointers.

The prototype's manual union-find at lines 929-1003 is replaced by building a
`cluster_connectivity_graph_t` and calling `merge_clusters(g, live_grouping)` at line 975.
The factored-out helper in `ClusteringFuncs.cxx:48-119` uses `std::map<int, std::set<int>>`
internally for connected-component bookkeeping, which is deterministic. However, it inherits
the vertex-index ordering from the caller's `map_cluster_index` (pointer-keyed; see §2.5).

### 2.3 Bugs and divergences

#### BUG — FV-x gate uses per-wpid map instead of global detector envelope (`clustering_neutrino.cxx:178`)

```cpp
// CURRENT (wrong for multi-TPC):
if (el_wcps.first.x()  < map_FV_xmin.begin()->second - map_FV_xmin_margin.begin()->second ||
    el_wcps.second.x() > map_FV_xmax.begin()->second + map_FV_xmax_margin.begin()->second ||
    cluster->get_length() < 6.0 * units::cm)
    continue;
```

`map_FV_xmin.begin()->second` picks an arbitrary per-wpid entry. The file already reads the
global-detector-envelope bounds at lines 119-125:

```cpp
WirePlaneId wpid_all(0);
double det_FV_xmin = dv->metadata(wpid_all)["FV_xmin"].asDouble();
double det_FV_xmax = dv->metadata(wpid_all)["FV_xmax"].asDouble();
// (margins are currently commented-out at lines 126-127)
```

This is exactly the pattern used in `clustering_separate.cxx:360-380`:

```cpp
WirePlaneId wpid_all(0);
double det_FV_xmin        = dv->metadata(wpid_all)["FV_xmin"].asDouble();
double det_FV_xmax        = dv->metadata(wpid_all)["FV_xmax"].asDouble();
double det_FV_xmin_margin = dv->metadata(wpid_all)["FV_xmin_margin"].asDouble();
double det_FV_xmax_margin = dv->metadata(wpid_all)["FV_xmax_margin"].asDouble();
```

**Fix (two steps):**

1. Un-comment lines 126-127 in `clustering_neutrino.cxx`:
   ```cpp
   double det_FV_xmin_margin = dv->metadata(wpid_all)["FV_xmin_margin"].asDouble();
   double det_FV_xmax_margin = dv->metadata(wpid_all)["FV_xmax_margin"].asDouble();
   ```
2. Replace line 178:
   ```cpp
   if (el_wcps.first.x()  < det_FV_xmin - det_FV_xmin_margin ||
       el_wcps.second.x() > det_FV_xmax + det_FV_xmax_margin ||
       cluster->get_length() < 6.0 * units::cm)
       continue;
   ```

This makes the FV-x check a single global detector envelope covering all TPCs, identical to the
`clustering_separate` convention. The per-wpid `map_FV_xmin` / `map_FV_xmax` maps (built at
lines 108-118) are then unused — they can be removed if no other use is found.

#### INTENTIONAL DIVERGENCE — extra `used_clusters.insert(cluster2)` (`clustering_neutrino.cxx:936-938`)

```cpp
// Toolkit only:
if (cluster2->get_length() < 15*units::cm &&
    cluster1->get_length() > 2 * cluster2->get_length())
    used_clusters.insert(cluster2);
```

The prototype (`ToyClustering_neutrino.h:902-903`) inserts only `cluster1`. The toolkit also
locks out `cluster2` for further merges when it is a very short fragment relative to `cluster1`.
This will produce modestly different merge graphs on events with many short fragments.
**Confirmed intentional by the developer — no change required.** Documented here for traceability.

#### OBSERVATION — `drift_dir_abs(1,0,0)` hard-coded (`clustering_neutrino.cxx:135`)

The comment correctly states "we do not care about the actual direction of drift_dir, so just
picking up the first instance." All dot-product uses of `drift_dir_abs` downstream are
sign-insensitive (looking for the axis, not the orientation). This is safe for single-face and
dual-face (opposite-drift) geometries as long as no branch takes a signed x-comparison against
`drift_dir_abs`. Spot-checks through the main loop confirm this is the case.

### 2.4 Efficiency

- **Dropped pre-PCA warmup** (`clustering_neutrino.cxx:323-326`, commented out): correct.
  The facade caches PCA lazily on first access, making the prototype's explicit `Calc_PCA`
  warmup loop redundant.
- **`Separate_2` blob-id-map API**: avoids heap allocation and `delete` of per-sub-cluster
  copies. Significant memory savings for large clusters.
- **`merge_clusters` factored out**: eliminates ~75 lines of hand-written union-find (prototype
  lines 929-1003). Shared code with the other clustering stages.
- **No remaining efficiency concern** in the main loop — geometric calculations are identical in
  complexity to the prototype.

### 2.5 Determinism — pointer-keyed containers

All of the following containers are keyed on raw `Cluster*` and are iterated in address order,
which varies across runs.

| Line | Container | Role |
|---|---|---|
| 279 | `std::map<Cluster*, std::pair<Cluster*, double>>` | closest-cluster map |
| 314 | `std::set<std::pair<Cluster*, Cluster*>>` | merge pairs |
| 316 | `std::set<Cluster*>` | `used_clusters` |
| 318 | `std::map<Cluster*, shared_ptr<Simple3DPointCloud>>` | extended cloud map |
| 320-321 | `std::map<Cluster*, geo_point_t>` (×2) | dir maps |
| 963 | `std::map<const Cluster*, int>` | Boost vertex index assignment |

**Impact:** The most consequential is line 963 (`map_cluster_index`). The vertex indices it
assigns to `boost::add_vertex` determine which component gets which id in `connected_components`,
and therefore which existing cluster "becomes" the merged cluster in `merge_clusters()`. On
events with more than two-cluster merges this can shift which cluster inherits the output ident.

**Recommended fix pattern (matching the established project convention):**
- Replace raw-pointer maps/sets with `cluster_less_functor`-ordered variants:
  ```cpp
  std::map<Cluster*, ..., cluster_less_functor>
  std::set<Cluster*, cluster_less_functor>
  ```
- For `map_cluster_index` (line 963), build a `std::vector` from
  `live_grouping.children()` (already deterministically ordered) and use array index
  directly instead of a map lookup.

### 2.6 Multi-APA/face handling

| Item | Verdict |
|---|---|
| `wpid_params` build loop (lines 63-92) | Correct: iterates all wpids from `live_grouping.wpids()`. |
| Per-wpid time-slice-width map (lines 107-118) | Correct. |
| FV-x gate (line 178) | **Bug** — see §2.3 above. |
| Y/Z FV bounds (lines 181-276) | Use `det_FV_ymax / ymin / zmax / zmin` from `wpid_all(0)` — correct. |
| `drift_dir_abs(1,0,0)` (line 135) | Intentional; sign-insensitive; safe. |
| `Separate_2` (clustering_separate.cxx:1600) | Loops all `(apa, face)` pairs internally — correct. |
| `judge_vertex(dv)` (Facade_Cluster.cxx:2393) | Accepts `dv`; no singleton — correct. |
| `merge_clusters(g, live_grouping)` | Geometry-agnostic Boost CC; correct. |

### 2.7 Callees reviewed

| Callee | File:line | Multi-TPC | Deterministic | Notes |
|---|---|---|---|---|
| `Separate_2` | `clustering_separate.cxx:1600` | Yes — loops `(apa,face)` | Yes (int-keyed ids) | Returns blob-to-group-id `vector<int>` |
| `JudgeSeparateDec_1` | `clustering_separate.cxx:318` | Yes — uses 3D x-extent, not per-TPC tick count | Yes | Replaces prototype's `get_num_time_slices()*time_slice_width` |
| `Cluster::judge_vertex(dv)` | `Facade_Cluster.cxx:2393` | Yes — `dv` replaces singleton | n/a | Compares against `dv` metadata per-wpid |
| `merge_clusters(g, live_grouping)` | `ClusteringFuncs.cxx:48` | Yes | Internally yes (`std::map<int,std::set<int>>`) — but inherits caller's vertex-index ordering | Boost CC |
| `Simple3DPointCloud` | used at lines 360, 539 | n/a | n/a | Replaces prototype's `ToyPointCloud` |

### 2.8 Suggested fixes

Priority order:

1. **[Bug, multi-TPC]** Line 178: replace `map_FV_xmin.begin()->second` with `det_FV_xmin`,
   uncomment margin variables at lines 126-127. (See §2.3 for full diff.)
2. **[Determinism]** Lines 279, 314, 316, 318, 320-321: add `cluster_less_functor` to
   pointer-keyed containers.
3. **[Determinism]** Line 963: replace `map_cluster_index` with a `live_grouping.children()`
   index vector (same pattern as `clustering_isolated.cxx:474`).

---

## 3. `Clustering_isolated`

### 3.1 Function map (prototype → toolkit)

| Prototype block | Prototype lines | Toolkit lines | Notes |
|---|---|---|---|
| Parameter setup | 5-27 | 66-131 | Prototype: single `TPCParams`. Toolkit: per-wpid loop. |
| Pre-sort (toolkit only) | — | 97-101 | Prototype does not sort; toolkit sorts descending by length. |
| Scope pre-pass (toolkit only) | — | 103-108 | No prototype analogue. |
| Classify big vs small | 28-83 | 132-211 | Identical thresholds; see §3.3 for `Separate_2` semantics. |
| Small ↔ big pairing | 85-112 | 217-260 | Identical. |
| Small ↔ small (first pass, 5 cm) | 113-138 | 262-282 | Identical. |
| Small ↔ small (second pass, 50 cm) | 140-162 | 283-300 | Identical. |
| Big ↔ big pairing | 164-224 | 303-382 | **Bug at toolkit line 331**; see §3.3. |
| Union-find merge + results map | 226-313 | 385-467 | Same logic; `results` uses `cluster_less_functor` in toolkit. |
| Final merge | 314-326 | 470-489 | Prototype returns map; toolkit merges in-place via `merge_clusters`. |

### 3.2 Logic parity with prototype

All threshold values match exactly:

| Threshold | Toolkit | Prototype |
|---|---|---|
| `range_cut` | 150 (:130) | 150 (:15) |
| `length_cut` | 20 cm (:131) | 20 cm (:16) |
| `Separate_2 dis_cut` | 2.5 cm (:162) | 2.5 cm (:45) |
| split length threshold | 60 cm (:159) | 60 cm (:43) |
| `small_big_dis_cut` | 80 cm (:218) | 80 cm (:85) |
| first `small_small_dis_cut` | 5 cm (:250) | 5 cm (:113) |
| second `small_small_dis_cut` | 50 cm (:282) | 50 cm (:145) |
| `big_dis_cut` | 3 cm (:304) | 3 cm (:164) |
| `big_dis_range_cut` | 16 cm (:305) | 16 cm (:165) |
| outside-fraction cut | `0.125*N \|\| >400` (:362) | same (:214) |
| pca_ratio override | `>60 cm && <0.0015` (:365) | same (:218) |

### 3.3 Bugs and divergences

#### BUG — `cluster1` mutated across inner-loop iterations (`clustering_isolated.cxx:307-331`)

```cpp
for (size_t i = 0; i != big_clusters.size(); i++) {
    Cluster *cluster1 = big_clusters.at(i);           // line 307: bound once per OUTER iteration
    for (size_t j = i + 1; j != big_clusters.size(); j++) {
        Cluster *cluster2 = big_clusters.at(j);       // line 310: fresh each INNER iteration
        ...
        if (!(cluster1->get_length() > cluster2->get_length()))
            std::swap(cluster1, cluster2);             // line 331: MUTATES cluster1!
        if (used_big_clusters.find(cluster2) != ...) continue;
        ...
    }
}
```

After the first inner iteration that executes the `swap`, `cluster1` on the next inner
iteration is whatever the previous pair's shorter cluster was — not `big_clusters.at(i)`. This
breaks the invariant that `cluster1` is always the i-th big cluster.

The prototype (`ToyClustering_isolated.h:174-188`) avoided this by handling both orderings with
explicit `if/else` branches that never mutated the outer variable.

**Fix:** scope the sorted pair to fresh locals inside the inner loop:

```cpp
for (size_t i = 0; i != big_clusters.size(); i++) {
    for (size_t j = i + 1; j != big_clusters.size(); j++) {
        Cluster *c1 = big_clusters.at(i);
        Cluster *c2 = big_clusters.at(j);
        if (c2->get_length() > c1->get_length()) std::swap(c1, c2); // c1 is always longer
        if (used_big_clusters.find(c2) != used_big_clusters.end()) continue;
        double pca_ratio = c2->get_pca().values.at(1) / c2->get_pca().values.at(0);
        double small_cluster_length = c2->get_length();
        ...
    }
}
```

Then replace all subsequent `cluster1`/`cluster2` references in the inner body with `c1`/`c2`.

#### MINOR — uninitialized `max_cluster` (`clustering_isolated.cxx:443`)

```cpp
Cluster *max_cluster;   // uninitialized — UB if cluster_set is empty
```

`cluster_set` is always non-empty by construction, so this is not reachable in practice.
Still worth `= nullptr` for defensive coding, matching common project style.

#### OBSERVATION — scope gating (`clustering_isolated.cxx:137`)

```cpp
if (!live_clusters.at(i)->get_scope_filter(scope)) continue;
```

This silently skips clusters the prototype would have processed. Ensure that integration tests
set scope filters correctly (all clusters pass the filter for the relevant scope) or the
toolkit will report different big/small counts from the prototype. No code change needed but
document in parity test setup.

### 3.4 Efficiency

- **Lazy point-cloud construction** replaces prototype's explicit `Create_point_cloud()` calls
  at prototype lines 89, 95, 168, 171, 306, 313. No functional change; the facade constructs
  the cloud on first access via `get_closest_points` / `get_closest_dis`.
- **`Separate_2` blob-id-map API** avoids heap allocation and `delete` cycles for all
  intermediate sub-clusters.
- **Sorted input vector** (`clustering_isolated.cxx:97-101`): makes the classification loop
  deterministic and matches the spirit of the prototype (which assumed length-sorted input from
  the caller).

### 3.5 Determinism — pointer-keyed containers

| Line | Container | Role |
|---|---|---|
| 215 | `std::set<std::pair<Cluster*, Cluster*>>` | merge pairs |
| 219 | `std::set<Cluster*>` | `used_small_clusters` |
| 303 | `std::set<Cluster*>` | `used_big_clusters` |
| 385 | `std::vector<std::set<Cluster*>>` | union-find sets |
| 386 | `std::set<Cluster*>` | `used_clusters` |
| 432 | `std::set<Cluster*>` | `temp_clusters` |
| 473 | `std::map<const Cluster*, int>` | Boost vertex index |

The `results` map at line 439 correctly uses `cluster_less_functor` — this is the established
project pattern. Apply it consistently to the other six containers above.

**Most consequential:** the element sets inside `merged_clusters` (line 385) are
`std::set<Cluster*>` and are iterated at lines 444 and 456 to build the `results` entry. Since
`cluster_less_functor` is used for `results` itself, the `max_cluster` tie-break (longest by
`get_length`) is deterministic. However the `(temp_cluster, dis)` pairs appended at line 465
are in pointer-iteration order — edge insertion order into the Boost graph, which can shift
component numbering. Use `std::set<Cluster*, cluster_less_functor>` for the inner sets.

**Sort at line 99** is length-only and uses an unstable comparator. Ties are broken by
`std::sort`'s internal ordering, which is nondeterministic. Use `cluster_less_functor` or
wrap with `std::stable_sort`.

### 3.6 Multi-APA/face handling

| Item | Verdict |
|---|---|
| Per-wpid param build (lines 66-95) | Correct — iterates all wpids. |
| `get_uvwt_range()` per-wpid sum (lines 138-146) | Correct — sums across all wpids, dividing time component by per-wpid `nticks_live_slice`. |
| `Facade::get_uvwt_range(cluster, b2id, id)` (lines 172-180) | Correct — per-wpid partitioned sum, see `Facade_Cluster.cxx:1788`. |
| `JudgeSeparateDec_1` | Correct — uses `get_earliest_latest_points()` (3D x-extent) instead of prototype's per-TPC tick count. |
| `Separate_2` | Correct — loops all `(apa, face)` pairs internally. |
| `Cluster::get_length()` (`Facade_Cluster.cxx:1759`) | Correct — iterates over every `(apa, face)` from `get_uvwt_range()`. |
| `drift_dir_abs(1,0,0)` (line 127) | Intentional — `JudgeSeparateDec_1` folds to `[0, π/2]`; sign-insensitive. |
| No hard-coded `face=0`/`apa=0` | Confirmed — no explicit zero-index in `clustering_isolated.cxx`. |

The multi-APA/face generalization of `clustering_isolated` is **correct**.

### 3.7 Callees reviewed

| Callee | File:line | Multi-TPC | Deterministic |
|---|---|---|---|
| `Separate_2` | `clustering_separate.cxx:1600` | Yes | Yes |
| `JudgeSeparateDec_1` | `clustering_separate.cxx:318` | Yes | Yes |
| `Cluster::get_length()` | `Facade_Cluster.cxx:1759` | Yes — loops wpid | n/a |
| `Facade::get_uvwt_range(b2id,id)` | `Facade_Cluster.cxx:1788` | Yes | Yes |
| `Facade::get_length(b2id,id)` | `Facade_Cluster.cxx:1826` | Yes | Yes |
| `cluster_less_functor` | `Facade_Cluster.cxx:2498` | n/a | Yes (content-based) |
| `merge_clusters(g, grouping, "isolated")` | `ClusteringFuncs.cxx:48` | Yes | Internally yes; inherits vertex-index order from caller |

### 3.8 Suggested fixes

Priority order:

1. **[Bug]** Lines 307-331: fix big-big swap scoping (see §3.3 for full diff).
2. **[Determinism]** Lines 215, 219, 303, 385 (inner sets), 386, 432: add `cluster_less_functor`.
3. **[Determinism]** Line 473: replace `map_cluster_index` with vector index.
4. **[Determinism]** Line 99: use `cluster_less_functor` or `std::stable_sort`.
5. **[Minor]** Line 443: `Cluster *max_cluster = nullptr;`.

---

## 4. `ExamineBundles`

### 4.1 Function map (prototype → toolkit)

The toolkit is a **functional rewrite**, not a line-by-line port.

| Prototype concept | Toolkit equivalent |
|---|---|
| Collect all mcells from main + other clusters → new `PR3DCluster` | Already in `Cluster` with blobs; no aggregation needed. |
| `Examine_graph(ct_point_cloud)` — rebuild graph, run CC | `Cluster::connected_blobs(dv, pcts, graph_name)` — same semantics, configurable graph flavor. |
| Main cluster: max overlap of mcells with old `main_cluster` | Blobs previously tagged `−1` in `("isolated","perblob")` pcarray; max-overlap of those tags with new CC component ids. |
| Length fallback | None in prototype. Toolkit adds length-based fallback when previous tag array is missing or stale. |
| Construct new `PR3DCluster` + `FlashTPCBundle` | Mutate per-blob tag in place via `put_pcarray`. No allocation. |
| Delete old bundles and `other_clusters` | Not needed; toolkit owns memory via the tree. |

### 4.2 Logic parity with prototype

This function is intentionally redesigned; it is not a line-for-line port. The semantic
equivalent is:

- **Prototype overlap criterion:** `SlimMergeGeomCell*` set intersection between each new CC
  component and the old main cluster's mcell set.
- **Toolkit overlap criterion:** count how many blobs in each new CC component id were previously
  tagged `−1` in the `isolated/perblob` pcarray. Semantically equivalent, provided the
  `perblob` array is in sync with the current blob-vector ordering.

The toolkit's new **length-based fallback** (lines 139-172, triggered when `flag_largest=true`)
has no prototype analogue. It kicks in when:
- The previous `isolated` pcarray size differs from the current blob count (stale/uninitialized), or
- No blob was previously tagged `−1` (first call on this cluster).

This is a deliberate behavioral improvement, not a divergence from the prototype's intent.

### 4.3 Bugs

None found in `clustering_examine_bundles.cxx`. This file is clean.

### 4.4 Efficiency

- **No heap allocation:** the prototype allocates/deletes `PR3DCluster` and `FlashTPCBundle`
  objects per bundle. The toolkit tags per-blob integers in place — O(1) memory vs O(N blobs).
- **`unique_ids` + `cluster_lengths` double pass** (lines 142-153): currently iterates
  `b2groupid` twice (once to build `std::set<int> unique_ids`, once implied via the map).
  Can be collapsed into a single sweep:
  ```cpp
  std::map<int, double> cluster_lengths;
  for (const auto& id : b2groupid) {
      if (id >= 0 && cluster_lengths.find(id) == cluster_lengths.end())
          cluster_lengths[id] = get_length(live_clusters.at(i), b2groupid, id);
  }
  ```
  Minor; not a bottleneck, but cleaner.

### 4.5 Determinism

All containers in the driver are integer-keyed:
- `std::map<int,int> overlap_counts` (line 108)
- `std::map<int,double> cluster_lengths` (line 141)
- `std::set<int> unique_ids` (line 142)

Iteration order is deterministic. Tie-breaking uses strict `>`, so the smallest component id
wins ties. This is deterministic and reproducible across runs.

The prototype used `std::set<SlimMergeGeomCell*>` for membership testing (not iteration) —
no nondeterminism there either.

**`ExamineBundles` is the deterministic reference implementation among the three reviewed functions.**

### 4.6 Multi-APA/face handling

The driver itself contains **no** `face`, `apa`, or `wpid` references. All geometry is
delegated to helpers:

| Helper | Multi-TPC status |
|---|---|
| `Cluster::connected_blobs(dv, pcts, graph_name)` | Passes `dv` and `pcts` — geometry-aware. The driver is clean. |
| `graph_algorithms(flavor, dv, pcts)` → `make_graph_relaxed` / `make_graph_ctpc` | **Audit target** — see §4.7 |
| `Facade::get_length(cluster, b2id, id)` (`Facade_Cluster.cxx:1826`) | Correct — iterates `wpid` from `get_uvwt_range()`, uses per-face tick/drift/pitch. |

**No hard-coded `face=0` / `apa=0` in `clustering_examine_bundles.cxx`.** Confirmed.

### 4.7 Follow-up audit target: `make_graph_relaxed` / `make_graph_ctpc`

`clustering_examine_bundles.cxx` delegates all graph construction to `graph_algorithms(flavor,
dv, pcts)` at `Facade_Cluster.cxx:2827`. This factory calls `make_graph_relaxed` or
`make_graph_ctpc` depending on the `graph_name_` config key. These factories are the most
likely location for hidden multi-APA/face assumptions (e.g. `wire_angles(0,0)`, `get_xyz(0)`,
`get_tpc_params(0,0)`).

**Recommended follow-up:** search `make_graph_relaxed` and `make_graph_ctpc` for:
```
grep -n "face(0)\|apa=0\|face=0\|wire_angles(0\|get_tpc_params(0" clus/src/Facade_Cluster.cxx
```
Any hits should be replaced with per-blob wpid dispatch.

### 4.8 Suggested fixes

1. **[Minor, efficiency]** Lines 142-153: collapse `unique_ids` + `cluster_lengths` into a
   single sweep (see §4.4).
2. **[Follow-up audit]** `make_graph_relaxed` / `make_graph_ctpc` — check for hard-coded
   `(0,0)` APA/face (outside scope of this file, tracked here for completeness).

---

## 5. Cross-cutting recommendations

### 5.1 Adopt `cluster_less_functor` as the project standard for pointer-keyed containers

`cluster_less_functor` (`Facade_Cluster.cxx:2498`) is already used in `clustering_isolated.cxx:439`
and in several other clustering stages. It provides content-based ordering (length, nchildren,
npoints, per-wpid range tuples) that is deterministic across runs. Every `std::set<Cluster*>`
or `std::map<Cluster*, ...>` in the three reviewed files that is **iterated** should use it.

Containers that are only accessed via `find()` (membership testing) do not affect algorithmic
output and are lower priority.

### 5.2 Use `WirePlaneId(0)` global envelope for detector-wide FV checks

Functions that test whether a cluster is inside the fiducial volume of the **whole detector**
(not a specific TPC) should use:
```cpp
WirePlaneId wpid_all(0);
double det_FV_xmin = dv->metadata(wpid_all)["FV_xmin"].asDouble();
```
not `map_FV_xmin.begin()->second`. The established pattern is in `clustering_separate.cxx:360-380`.
`clustering_neutrino.cxx:178` needs this fix (see §2.3).

### 5.3 Document the `scope_filter` contract

Both `clustering_neutrino` and `clustering_isolated` gate on `cluster->get_scope_filter(scope)`.
The prototype has no such concept. For parity tests, all clusters in the test fixture must pass
the relevant scope filter or the toolkit will silently skip them, producing different results
from the prototype without any error. A comment or assert in the test harness is recommended.

### 5.4 Replace `map_cluster_index` pattern with vector index

Both `clustering_neutrino.cxx:963` and `clustering_isolated.cxx:473` build
`std::map<const Cluster*, int>` just to assign consecutive vertex indices for the Boost graph.
Since `live_grouping.children()` returns a deterministically-ordered vector, the index can be
looked up from the vector directly:
```cpp
auto live_clusters_all = live_grouping.children();  // already deterministic
for (size_t i = 0; i < live_clusters_all.size(); i++) {
    ilive2desc[i] = boost::add_vertex(i, g);
}
// lookup: std::find(live_clusters_all.begin(), live_clusters_all.end(), ptr) - live_clusters_all.begin()
```
Or build a `std::unordered_map<Cluster*, int>` keyed on pointer but populated in
`children()` order (so values are deterministic even if the map itself is unordered).

---

## 6. Suggested fix backlog

Grouped by priority. Each item lists the exact file and line(s) to change.

### P0 — Correctness bugs

| # | File:line | Description |
|---|---|---|
| 1 | `clustering_neutrino.cxx:126-127, 178` | Un-comment margin vars; replace `map_FV_xmin.begin()->second` with `det_FV_xmin - det_FV_xmin_margin` (global `wpid_all` envelope). |
| 2 | `clustering_isolated.cxx:307-334` | Fix big-big swap scoping: use fresh inner-loop locals `c1, c2` instead of mutating outer `cluster1`. |

### P1 — Determinism (affects reproducibility across runs)

| # | File:line | Description |
|---|---|---|
| 3 | `clustering_neutrino.cxx:279, 314, 316, 318, 320-321` | Add `cluster_less_functor` to all pointer-keyed containers. |
| 4 | `clustering_neutrino.cxx:963` | Replace `map<const Cluster*,int>` with vector-index approach. |
| 5 | `clustering_isolated.cxx:215, 219, 303, 385-inner, 386, 432` | Add `cluster_less_functor` to all pointer-keyed containers. |
| 6 | `clustering_isolated.cxx:473` | Replace `map<const Cluster*,int>` with vector-index approach. |
| 7 | `clustering_isolated.cxx:99` | Use `cluster_less_functor` comparator or `std::stable_sort` for initial sort. |

### P2 — Minor / code quality

| # | File:line | Description |
|---|---|---|
| 8 | `clustering_isolated.cxx:443` | Initialize `max_cluster = nullptr`. |
| 9 | `clustering_examine_bundles.cxx:142-153` | Collapse `unique_ids` + `cluster_lengths` into one loop. |
| 10 | `make_graph_relaxed` / `make_graph_ctpc` in `Facade_Cluster.cxx` | Audit for hard-coded `face=0`/`apa=0`. (Follow-up task.) |
