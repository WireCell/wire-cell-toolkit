# Clustering 6-Function Review

**Review date:** 2026-04-10  
**Scope:** `ClusteringLiveDead`, `ClusteringExtend`, `ClusteringRegular`,
`ClusteringParallelProlong`, `ClusteringClose`, `ClusteringExtendLoop`  
**Criteria:** (1) logic fidelity vs. prototype, (2) bugs, (3) efficiency,
(4) determinism/randomness, (5) multi-APA/face correctness

**Prototype:** `prototype_base/2dtoy/src/ToyClustering_*.h`  
**Toolkit:** `clus/src/clustering_*.cxx`

---

## Overall Status

All 6 functions are functionally correct and have had all major prior fixes
applied. One minor code fix was made in this session; one stale doc was updated.

| Function | Logic | Bugs | Efficiency | Determinism | Multi-APA | Net |
|---|---|---|---|---|---|---|
| `Clustering_live_dead` | ✓ | ✓ | ✓ | ✓ | ✓ | Pass |
| `Clustering_extend` | ✓ | Fixed* | ✓ | ✓ | ✓ | Pass |
| `Clustering_regular` | ✓ | ✓ | ✓ | ✓ | ✓ | Pass |
| `Clustering_parallel_prolong` | ✓ | ✓ | ✓ | ✓ | ✓ | Pass |
| `Clustering_close` | ✓ | ✓ | ✓ | ✓ | N/A | Pass |
| `ClusteringExtendLoop` | ✓ | ✓ | ✓ | ✓ | ✓ | Pass |

\* One minor fix applied in this session (see §2.2 below).

---

## 1. ClusteringLiveDead

**Files:** `clus/src/clustering_live_dead.cxx` ↔ `prototype_base/2dtoy/src/ToyClustering_dead_live.h`

### 1.1 Logic Fidelity — Pass

Core algorithm is identical to prototype:
1. Build `dead → [live clusters]` mapping via `is_connected` (wire+time overlap)
2. For each dead cluster bridging ≥2 live clusters, test each pair
3. Convergence loop: alternating nearest-neighbour search for closest points
4. Compute Hough-transform directions, distance, angle metrics
5. Four merge-criteria branches by cluster length (both-short / one-short / both-long)
6. Transitive merge via Boost connected components

**Toolkit enhancement (intentional):** Direction flip protection (D1). After
`dir1 = vhough_transform(mcell1_center, 30cm)`, the toolkit also computes
`dir5 = vhough_transform(p1, 30cm)` and flips `dir1 *= -1` if they disagree
by more than 120°. Not in prototype. Prevents ambiguous Hough results from
producing a backwards direction vector.

**Prototype bug fixed:** In the `length_2 > 12cm && length_1 ≤ 12cm` branch the
prototype has a copy-paste error (`angle_diff2 <= 45 || angle_diff2 <= 45` — same
variable twice). The toolkit correctly uses `angle_diff1 <= 45 || angle_diff2 <= 45`.

### 1.2 Bugs — None

No bugs found in current code.

### 1.3 Efficiency

Prior fixes applied and verified:
- `tested_pairs`: replaced `std::set<pair<Cluster*,Cluster*>>` with index-based
  `unordered_set<size_t>` keyed on `min_idx * nlive + max_idx`. O(1) lookup,
  half the entries, no pointer-address dependence.
- `map_cluster_index`: `unordered_map` — O(1) lookup.
- `dead_live_cluster_mapping`: `emplace` return avoids double lookup.

### 1.4 Determinism

- `sort_clusters(live_clusters)` and `sort_clusters(dead_clusters)` — consistent
  full tie-breaking (length → nblobs → npoints → wire ranges → PCA center →
  pointer last resort).
- `map_cluster_index` uses integer indices from sorted order → symmetric pair key
  is content-based, not address-based.
- `wpid_params` / `wpid_U/V/W_dir` maps keyed by `WirePlaneId` (struct with
  apa/face/layer) — hardware-topology deterministic.

### 1.5 Multi-APA/Face — Pass

Correct per-wpid pattern (fixed in prior session, confirmed here):
- Lines 71–78: calls `compute_wireplane_params` over all grouping wpids → builds
  four per-wpid maps.
- Lines 218–221: looks up geometry per-cluster:
  ```cpp
  auto wpid_1 = cluster_1->wpid(mcell1_center);
  auto wpid_2 = cluster_2->wpid(mcell2_center);
  const auto& [drift_dir_1, angle_u_1, angle_v_1, angle_w_1] = wpid_params.at(wpid_1);
  const auto& [drift_dir_2, angle_u_2, angle_v_2, angle_w_2] = wpid_params.at(wpid_2);
  ```
- Each `is_angle_consistent` call uses the correct APA's angles for that cluster.
- `is_connected` uses integer wire ranges — already APA-agnostic.

---

## 2. ClusteringExtend (and helpers Clustering_4th_prol/para/reg/dead)

**Files:** `clus/src/clustering_extend.cxx` ↔ `prototype_base/2dtoy/src/ToyClustering_extend.h`

### 2.1 Logic Fidelity — Pass

All four flag modes (1=prolong, 2=parallel, 3=regular, 4=dead) follow prototype
flow identically. Decision trees in each `Clustering_4th_*` helper are functionally
identical to prototype. Toolkit enhancement: per-wpid geometry maps instead of
MicroBooNE hardcoded angles (correct multi-detector generalization).

### 2.2 Bug Fixed (this session) — flag=2 lowest_p loop missing `used_clusters` guard

**Location:** `clustering_extend.cxx`, flag=2 parallel branch, lowest_p loop.

The flag=2 branch tests both the `highest_p` and `lowest_p` endpoints of
`cluster_1`. The `highest_p` loop (lines 727–743) correctly:
- Checks `if (used_clusters.find(cluster_2) != used_clusters.end()) continue;`
- Inserts short merged clusters: `if (cluster_2->get_length() < 15cm) used_clusters.insert(cluster_2);`

The `lowest_p` loop previously had neither check. The same omission exists in the
prototype (inherited behaviour), so this is not a regression. The boost graph
deduplicates via connected components so duplicate edges caused no data corruption.
However, the asymmetry was inconsistent with the `highest_p` and `flag=1`/`flag=3`
patterns. **Fix applied:** added the `find` guard and conditional insert to the
`lowest_p` loop, bringing it in line with all other endpoint loops.

### 2.3 Efficiency

Prior fixes verified:
- `live_clusters` sorted with `sort_clusters()` before the outer loop.
- `used_clusters`: `unordered_set` — O(1) find/insert.
- `map_cluster_index`: `unordered_map` — O(1) lookup.
- `Find_Closest_Points` / `get_strategic_points`: sorted by `geo_point_t` content,
  deduplicated before use; no pointer-address ordering.

No further algorithmic improvements found. The O(N²) outer loop is inherent to the
algorithm (same in prototype). KD-tree caching is already used by the Facade layer.

### 2.4 Determinism

- Sorted `live_clusters` → outer loop order is content-based.
- `used_clusters` / `map_cluster_index` are insertion-order independent.
- `wpid_U/V/W_dir` maps keyed by `WirePlaneId` — content-based deterministic.

### 2.5 Multi-APA/Face — Pass

- Lines 570–595: builds `wpid_U/V/W_dir` maps for all grouping wpids.
- Each helper looks up geometry per-cluster via `cluster.wpid(point)`.
- `drift_dir_abs(1,0,0)` hardcoding in `Clustering_4th_reg`/`_dead` is correct:
  used only in `fabs(angle - π/2)` checks which are sign-insensitive w.r.t. drift
  direction, so ±X drift faces both work.
- Cross-APA pairs resolved via `get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv)`
  (two-point overload, line 732 of `Facade_Util.cxx`) — returns the APA whose
  bounding box has the longer ray intersection. Both wpids are always in the maps
  (built from `grouping.wpids()`), so no out-of-range risk.

---

## 3. ClusteringRegular

**Files:** `clus/src/clustering_regular.cxx` ↔ `prototype_base/2dtoy/src/ToyClustering_reg.h`

### 3.1 Logic Fidelity — Pass

Decision tree in `Clustering_1st_round` is functionally identical to prototype.
The prototype's no-op `for (int kk=0;kk!=1;kk++)` outer wrapper is correctly omitted.
Unused `angle5`/`angle5_1` variables from prototype are correctly dropped.

**Toolkit enhancement:** Cross-APA merge special case (lines 286–288): for tracks
> 100 cm spanning different APAs/faces (PDHD geometry), the toolkit applies a
relaxed merge criterion. Not in prototype. Intentional.

### 3.2 Bugs — None

### 3.3 Efficiency

`Find_Closest_Points`, `calc_ave_pos`, `vhough_transform` all use Facade's cached
KD-tree. No redundant rebuilds. O(N²) pair loop is inherent.

### 3.4 Determinism

- `sort_clusters(live_clusters)` at line 473.
- `unordered_map<const Cluster*, int> map_cluster_index`.
- `wpid_U/V/W_dir` maps keyed by `WirePlaneId`.

### 3.5 Multi-APA/Face — Pass

Per-wpid maps built at lines 432–457. Cluster-pair geometry resolved via the
safe two-point `get_wireplaneid` overload. Cross-APA wpid choice uses the
longer-bounding-box criterion, documented in `multi-apa-analysis.md`.

---

## 4. ClusteringParallelProlong

**Files:** `clus/src/clustering_parallel_prolong.cxx` ↔ `prototype_base/2dtoy/src/ToyClustering_para_prol.h`

### 4.1 Logic Fidelity — Pass

`Clustering_2nd_round` decision tree is functionally identical to prototype.
Unused `flag_para_U`/`flag_para_V` variables from prototype are correctly dropped.

### 4.2 Bugs — None

### 4.3 Efficiency

Same cached KD-tree usage as regular. No issues.

### 4.4 Determinism

- `sort_clusters(live_clusters)`.
- `unordered_map<const Cluster*, int> map_cluster_index`.
- `compute_wireplane_params` builds deterministic `WirePlaneId`-keyed maps.

### 4.5 Multi-APA/Face — Pass

Per-wpid geometry via `compute_wireplane_params`. Single dominant `wpid_ps`
(longer bounding-box APA) used for both direction checks — safe because
cross-APA parallel-prolong candidates are geometrically unlikely (they would
require the same Y-Z position on opposite sides of a cathode).

Dead output parameters `wpid_params` and `apas` from `compute_wireplane_params`
are unused after the call (same in `clustering_regular`). Cannot be removed
without changing the shared helper's signature; other callers may use them.
Acceptable.

---

## 5. ClusteringClose

**Files:** `clus/src/clustering_close.cxx` ↔ `prototype_base/2dtoy/src/ToyClustering_close.h`

### 5.1 Logic Fidelity — Pass

`Clustering_3rd_round` merge logic is identical. Thresholds (0.5cm, 1.0cm, 2.0cm,
length_cut), Hough/dipole/nearby checks, and short-cluster skip all match prototype.

**Prototype UB fixed:** `num_p1`, `num_p2`, `num_tp1`, `num_tp2` left uninitialized
in prototype when `dis >= 2cm`. Toolkit uses brace-init `{0}`.

### 5.2 Bugs — None

### 5.3 Efficiency

One pass through N(N-1)/2 pairs; no inner rescans. No KD-tree caching issues.

### 5.4 Determinism

- `sort_clusters(live_clusters)`.
- `unordered_set<const Cluster*> used_clusters`.
- `unordered_map<const Cluster*, int> map_cluster_index`.

### 5.5 Multi-APA/Face — N/A

`Clustering_3rd_round` decisions use only distance, cluster lengths, Hough
directions, and local point counts — no wire-plane geometry. Multi-APA extension
is automatic and requires no per-wpid maps.

---

## 6. ClusteringExtendLoop

**Files:** `clus/src/clustering_extend.cxx` (`ClusteringExtendLoop::visit`)
↔ `prototype_base/2dtoy/src/ToyClustering.cxx` lines 293–323

### 6.1 Logic Fidelity — Pass

| Step | Prototype | Toolkit |
|------|-----------|---------|
| Busy-event guard | `if (size > 1100) num_try = 1` | `if (nchildren() > 1100) num_try = 1` |
| flag=1 (prolong) | `Clustering_extend(..., 1, 150cm, 0)` | `clustering_extend(..., 1, 150cm, 0)` |
| flag=2 (parallel) | `Clustering_extend(..., 2, 30cm, 0)` | `clustering_extend(..., 2, 30cm, 0)` |
| flag=3 (regular) | `Clustering_extend(..., 3, 15cm, 0)` | `clustering_extend(..., 3, 15cm, 0)` |
| flag=4, i=0 | `Clustering_extend(..., 4, 60cm, 0)` | `clustering_extend(..., 4, 60cm, 0)` |
| flag=4, i>0 | `Clustering_extend(..., 4, 35cm, i)` | `clustering_extend(..., 4, 35cm, i)` |

`num_try` default: 0 in toolkit config; pipeline JSON sets `num_try: 3` explicitly
— runtime identical to prototype's hardcoded 3.

The prototype also calls `Clustering_extend(..., 4, 60cm, 0, 15cm, 1)` once before
the main loop. This pre-loop call is handled in the toolkit pipeline as a standalone
`ClusteringExtend` visitor (`flag=4`, `length_cut=60cm`, `num_dead_try=1`) — not
part of `ClusteringExtendLoop`.

### 6.2 Bugs — None

### 6.3 Efficiency

Each call to `clustering_extend` rebuilds from scratch: geometry maps, sorted
cluster list, Boost graph. This is O(N²) × 5 calls × num_try iterations. The
same pattern is in the prototype; no incremental update was ever implemented.
Acceptable for current event sizes.

### 6.4 Determinism

Inherits all determinism fixes from `clustering_extend`. No containers of its own.

### 6.5 Multi-APA/Face — Pass

Inherits from `clustering_extend`.

---

## 7. Changes Made in This Session

### Code change

**`clus/src/clustering_extend.cxx`** — flag=2 lowest_p loop: added
`used_clusters.find` guard and conditional insert after successful merge,
matching the `highest_p` loop and all other flag patterns.

### Documentation update

**`clus/docs/clustering/multi-apa-analysis.md`** — The "Bug: ClusteringLiveDead
— Hard Exception for Multi-APA" section was rewritten to describe the fix that
was applied. Summary table updated from `BUG` to `Fixed`.

---

## 8. Functions / Utilities Verified Correct

| Name | Location | Notes |
|---|---|---|
| `sort_clusters` | `Facade_Cluster.cxx` | Full tie-breaking; pointer is last resort |
| `cluster_less` | `Facade_Cluster.cxx` | PCA-center self-comparison bug fixed (prior session) |
| `merge_clusters` | `ClusteringFuncs.cxx` | Boost connected components; deterministic |
| `compute_wireplane_params` | `Clustering_Util.cxx` | Builds per-wpid maps; deterministic |
| `get_wireplaneid` (2-pt) | `Facade_Util.cxx:732` | Ray–bbox intersection; safe for cross-APA |
| `Find_Closest_Points` | `clustering_extend.cxx` | Sorted strategic points; early exit |
| `get_strategic_points` | `clustering_extend.cxx` | Sorted by `geo_point_t`; no pointer order |
| `is_angle_consistent` | `Facade_Util.cxx` | 8-param version with `num_cut`; correct |

**Note on `get_wireplaneid` single-point overload** (`Facade_Util.cxx:724`):
This overload calls `dv->contained_by(point)` and could in principle return a wpid
not in the per-APA maps if the point falls in an unmapped region. However, none of
the 6 functions reviewed here use this overload — they all call the two-point version.
Callers in `connect_graph_relaxed.cxx`, `clustering_protect_overclustering.cxx`, etc.
are outside the scope of this review.
