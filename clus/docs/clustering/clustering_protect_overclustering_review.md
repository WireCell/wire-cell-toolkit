# Review: `ClusteringProtectOverclustering` / `Clustering_protect_overclustering`

Review date: 2026-04-10

## 1. Purpose and algorithm

`Clustering_protect_overclustering` detects clusters that are physically disjoint but were incorrectly merged by earlier clustering stages, and splits them. It works by building a local boost graph over the cluster's points (one vertex per point, edges for spatially nearby points sharing wire-plane indices), computing connected components, and calling `grouping->separate()` if more than one component survives after path-quality checks and MST validation.

The prototype equivalent is `WCP2dToy::Examine_overclustering` + `WCP2dToy::Clustering_protect_overclustering` in `prototype_base/2dtoy/src/ToyClustering_protect_overclustering.h` (lines 5–911).

---

## 2. Port status

| Item | Status |
|---|---|
| Class wiring (`IConfigurable`, `IEnsembleVisitor`) | Ported |
| Config mixins (`NeedDV`, `NeedPCTS`, `NeedScope`) | Ported |
| Pipeline registration (`WIRECELL_FACTORY`) | Ported |
| `clus.jsonnet` helper `protect_overclustering(name)` | Present (`clus.jsonnet:343`) |
| Active in uboone flow | **Commented out** (`root/test/uboone.jsonnet:268`) |
| Active in sbnd flow | **Commented out** (`cfg/pgrapher/experiment/sbnd/clus.jsonnet:149`) |
| Active in dune-vd flow | **Commented out** (three sites in `dune-vd/clus.jsonnet:281,427,583`) |

The component is built but not scheduled in any experiment config. Enable with `protect_overclustering("...name...")` in the relevant jsonnet.

---

## 3. Side-by-side walkthrough

### 3.1 Driver function

Prototype `Clustering_protect_overclustering` (lines 5–45):
- Iterates `live_clusters`, calls `Examine_overclustering`.
- If sub-clusters returned, schedules original for deletion and sub-clusters for addition.
- Recomputes `cluster_length_map` for new clusters using single-TPC `TPCParams` pitches/angles.

Toolkit `clustering_protect_overclustering` (lines 809–843 original; now ~797–831):
- Iterates a copy of `live_grouping.children()` — safety: `grouping->separate(..., true)` destroys the original cluster, but the copy was taken before the loop so subsequent pointers are unaffected.
- Calls `Separate_overclustering(cluster, dv, pcts, scope)` which performs the separation in-place via `grouping->separate`. No length-map bookkeeping needed (toolkit manages lengths elsewhere).
- Sets `cluster->set_default_scope(scope)` before the call, meaning the cluster is re-interpreted in the requested scope for this pass. This differs from the prototype (no scope concept). Intent: correct.

### 3.2 Workhorse function (`Separate_overclustering` vs `Examine_overclustering`)

Both functions follow the same five-phase structure:

**Phase 1 — Intra-blob edges.** For every blob, find point pairs sharing nearby wire indices (within `max_wire_interval` / `min_wire_interval`) and add edges. Faithful port: identical logic, same `set_intersection` approach. ✓

**Phase 2 — Inter-blob edges.** Identify overlapping blob pairs within the same time slice (and ±1 or ±2 slices), then connect their points similarly. Toolkit adds multi-APA/face handling here — see §5.

**Phase 3 — Connected components.** `boost::connected_components`. If `num <= 1`, return early (no split needed).

**Phase 4 — Per-component point clouds, pairwise distance + Hough direction + path quality.** For each component pair (j,k):
  - Compute closest-point distance `index_index_dis[j][k]`.
  - If both clouds are large enough (Hough condition), compute directional closest points `index_index_dis_dir1/2[j][k]`.
  - Walk a 1 cm-step path between the closest points; count `num_bad` steps that fail `is_good_point`. Invalidate the pair if too many bad steps.
  - Build MSTs on the surviving pairs.

**Phase 5 — Re-component and separate.** Add MST edges back to the main graph, run `connected_components` again. If `num1 > 1`, call `grouping->separate`.

---

## 4. Logic fidelity findings

### F1. `closest_index` semantics: intentional toolkit extension (confirmed by user)

Prototype (`ToyClustering_protect_overclustering.h:275, 368–376, 508–517`):
- `std::map<std::pair<int,int>, std::pair<int,double>> closest_index` — one nearest-neighbor per `(point_index1, time2)` key.
- Adds exactly **one edge** per key (the single nearest).

Toolkit (original lines 270–461):
- `std::map<std::pair<int,int>, std::set<std::pair<double,int>>> closest_index` with `max_num_nodes = 5` — keeps the five nearest neighbors per key.
- In the final loop: adds all five edges unless the nearest is already `> 0.25 cm`, in which case stops after the first.
- Net effect: for tight pairs (nearest ≤ 0.25 cm), up to 5 edges are inserted instead of 1, adding extra inter-component connectivity. This is a deliberate toolkit improvement to reduce sensitivity to a single nearest-point measurement.

The `if (it5 == begin && dis > 0.25*cm) break` idiom means: "when the closest point is already well-connected (< 0.25 cm gap), keep all 5; when it is already far (> 0.25 cm), one edge is enough." **This is intentional and preserved.**

### F2. Discarded return value of `Separate_overclustering`

`clustering_protect_overclustering` (line 825 original) calls `Separate_overclustering(cluster, ...)` and discards the returned `map<int,Cluster*>`. The split already happened via `grouping->separate(..., true)` inside the function. Functionally equivalent to the prototype pattern. No issue.

### F3. Intra-blob edge insertion differences

Prototype (line 201): `auto edge = add_edge(index1, index2, *graph)` — unweighted edge; weight assigned separately via `(*graph)[edge.first].dist = ...`.  
Toolkit (line 177 original): `add_edge(index1, index2, dis, *graph)` — weighted edge directly.  
Same logical effect; different graph representation. ✓

---

## 5. Multi-APA / face handling

### M1. `af_time_slices` — correct

Toolkit lines 198–209 build `map<apa, map<face, vector<time_slices>>>` from `cluster->time_blob_map()`. The prototype uses a flat `map<int, SMGCSet> time_cells_set_map` (single TPC). The toolkit correctly decomposes by `(apa, face)` before the adjacency loop.

### M2. `nticks_per_slice` per (apa, face) — correct

Adjacency tests (±1 / ±2 time slices) use `grouping->get_nticks_per_slice().at(apa).at(face)` as the multiplier. The prototype hard-codes `1` and `2`. This correctly adapts to variable tick widths in multi-face geometries.

### M3. Intra-face-only policy — deliberate

`connected_mcells` is populated only for blob pairs within the same `(apa, face)`. Cross-face blob connectivity is not considered. This is the chosen policy: blobs from different faces are never connected as candidates in this function. For PDHD and dune-vd the invariant is: a cluster that spans faces will already be fragmented at the face boundary by earlier clustering; protect_overclustering only reconnects within each face. **Documented; no change.**

### M4. Wire-plane-type per blob — safe

`get_max_wire_type()` and `get_min_wire_type()` are properties of each blob (not shared per detector). Since inter-blob connectivity is restricted to same-face pairs, no cross-face wire-type mismatch can occur.

### M5. Path-check `get_wireplaneid` + `transform->backward` — correct

The path-check at lines 590/625/660 (original) uses `Facade::get_wireplaneid(test_p, wpid_p1, wpid_p2, dv)` which resolves which APA/face `test_p` belongs to. The `backward` transform signature is `(pos, t0, face, apa)` matching `IPCTransform.h:24`. Argument order verified correct. ✓

---

## 6. Bugs / correctness

### B1. Dead variable `scope_name` — **FIXED**

Original line 799: `auto scope_name = cluster->get_scope_transform(scope);` — result never used. Removed.

### B2. Unused `num_edges` counter — **FIXED**

`int num_edges = 0` was declared, incremented in the intra-blob loop and the final edge loop, and then cast away `(void)num_edges`. Removed all references.

### B3. Early return after unnecessary sort — **FIXED**

Original: `std::sort(ordered_components...)` was called before `if (num <= 1) return {}`. For the common case of a single connected component the sort was wasted work. Moved the guard before the sort.

### B4. Iteration-invalidation safety — safe, documented

`live_clusters` is a copy taken before the loop. `grouping->separate(..., true)` destroys the original cluster but does not modify the copy vector. Subsequent iterations access pointers to *other* clusters that are not touched by the current separation. Safe. Added comment in code.

---

## 7. Determinism

### D1. `map<const Blob*, ...> map_mcell_wind_wcps[3]` — **FIXED**

Original declaration used `const Blob*` as map key, making iteration order address-dependent. Although the outer map was never iterated (only indexed directly), the pointer key made correctness fragile and non-obvious.

**Fix:** Changed to `std::vector<std::map<int,std::set<int>>>` sized by `mcells.size()`, indexed by blob position in `mcells`. Added `std::map<const Blob*, int> blob_to_idx` built from the same vector. All six lookup sites updated to `blob_to_idx.at(mcell)`.

`children()` order (from `NaryTree::nursery_`) is insertion-order and therefore deterministic. `mcells` vector is deterministic. `blob_to_idx` is content-keyed (pointer values used as map key but the indices are the stable values). ✓

### D2. `BlobSet` iteration — safe (no change needed)

`BlobSet = std::set<const Blob*, BlobLess>` where `BlobLess` delegates to `blob_less(a, b)` (defined in `Facade_Blob.cxx:248`) comparing by `wpid`, `npoints`, `charge`, etc. — **content-based, deterministic**. Iteration in the time-slice adjacency loop is deterministic.

### D3. `connected_mcells` vector and `ordered_components` sort — deterministic

`connected_mcells` is insertion-ordered, with iteration driven by the deterministic `BlobSet` (D2). `ordered_components` is sorted by `min_vertex` (stable integer key). Both are deterministic. ✓

### D4. `closest_index` — deterministic

Key type `pair<int,int>` is deterministic. Values are `set<pair<double,int>>` — sorted by distance first (stable for unequal distances; tie-broken by `index2` integer, deterministic). ✓

---

## 8. Efficiency

### E1. Scope-transform hoist — **FIXED**

The three path-check blocks each re-evaluated `cluster->get_default_scope().hash() != cluster->get_raw_scope().hash()` and re-called `pcts->pc_transform(...)` and `cluster->get_cluster_t0()` on every step of the 1 cm walk. These are constant for the entire call. Hoisted once before the `for (j)` loop:

```cpp
const bool needs_scope_transform = cluster->get_default_scope().hash() != cluster->get_raw_scope().hash();
const auto scope_transform = needs_scope_transform
    ? pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()))
    : std::shared_ptr<IPCTransform>{};
const double cluster_t0 = needs_scope_transform ? cluster->get_cluster_t0() : 0.0;
```

### E2. `closest_index` double lookup — **FIXED**

Both update sites used `find(key) == end` + `operator[key]` (2–3 lookups) + a manual `it++` loop to advance the iterator for trimming. Replaced with:
```cpp
auto& ci_entry = closest_index[key];
ci_entry.insert(std::make_pair(dis, index2));
if (ci_entry.size() > max_num_nodes) {
    ci_entry.erase(std::next(ci_entry.begin(), max_num_nodes), ci_entry.end());
}
```
`operator[]` on a missing key inserts a default-constructed empty set — correct. Single lookup, no manual loop.

### E3. Operator-precedence ambiguity — **FIXED**

The Hough-condition gate:
```cpp
if (num < 100 && ... > 400 || ... > 500 && ...) {
```
(`&&` binds tighter than `||`) was previously papered over with `#pragma GCC diagnostic ignored "-Wparentheses"`. Explicit parentheses added; pragma removed.

### E4. Deferred efficiency items (for follow-up PR)

- **E.d1** — Two ~150-line inter-blob blocks ("test 2 against 1" / "test 1 against 2") are near-identical. Extract into one helper `accumulate_closest(mcell_src, mcell_dst, wcps_src, winds, blob_to_idx, map_mcell_wind_wcps, closest_index, max_num_nodes)`.
- **E.d2** — Three path-check blocks (for `index_index_dis`, `dir1`, `dir2`) share the same body modulo which matrix entry is updated. Extract into `check_path(p1, p2, wpid1, wpid2, dv, cluster, scope_transform, needs_scope_transform, cluster_t0) -> int num_bad`.
- **E.d3** — `cluster->point3d(index)` is called twice for the same `index1`/`index2` pair inside the intra-blob and inter-blob loops. If `point3d` involves non-trivial lookup, cache as `std::vector<geo_point_t>` once.
- **E.d4** — `af_time_slices` rebuilds `vector<int> time_slices_vec` per (apa,face) pair when the source map already has the keys in sorted order. Iterate `time_blob_map()` directly.

---

## 9. Summary of applied fixes

| ID | Category | Change |
|---|---|---|
| B1 | Bug | Delete dead `scope_name` variable |
| B2 | Bug | Delete `num_edges` counter and `(void)num_edges` |
| B3 | Bug | Move `if (num <= 1) return {}` before `std::sort` |
| D1 | Determinism | Replace `map<Blob*, ...>` with `vector<...>` indexed by blob position |
| E1 | Efficiency | Hoist scope-transform check and `pc_transform` call out of per-step inner loop |
| E2 | Efficiency | Replace double-lookup `find+operator[]` with single `operator[]+std::next` |
| E3 | Clarity | Add explicit parentheses to Hough condition; remove `#pragma GCC diagnostic` |

---

## 10. Verification

1. **Build:** `./wcb --target=WireCellClus` from the toolkit root.
2. **Enable the visitor:** in `cfg/pgrapher/experiment/uboone/clus.jsonnet:268` (or equivalent) uncomment the `protect_overclustering(...)` call. Run on a small input:
   ```
   wire-cell -c cfg/pgrapher/experiment/uboone/main.jsonnet <input.npz>
   ```
   Confirm no crash and inspect the live grouping cluster count before/after.
3. **Determinism:** run twice on the same input, compare output cluster counts byte-by-byte.
4. **Regression check for top-5 edge behavior:** compare cluster count with the old single-nearest behavior (revert `closest_index` temporarily) on a sample event to confirm the top-5 extension does not over-connect unrelated tracks.
