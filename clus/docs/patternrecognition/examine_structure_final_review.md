# Review: `examine_structure_final` / `examine_structure_final_{1,1p,2,3}` / `search_for_vertex_activities`

**Reviewed by:** Claude Code  
**Date:** 2026-04-11  
**Branch:** `apply-pointcloud`  
**Prototype ref:** `prototype_base/pid/src/NeutrinoID_final_structure.h` (697 lines), `prototype_base/pid/src/NeutrinoID_improve_vertex.h` (lines 1039–1211)  
**Function map:** `clus/docs/porting/neutrino_id_function_map.md` lines 97–103, 143

---

## Implementation status

| Item | Description | Status | File |
|---|---|---|---|
| B.1 | ESF2 / ESF3: `flag_start` determined from `sg` but applied to walk `sg1` | ✅ Applied | `NeutrinoStructureExaminer.cxx:3066, 3253` |
| B.2 | `search_for_vertex_activities` Round 2: raw `boost::edges` instead of `ordered_edges` | ✅ Applied | `NeutrinoVertexFinder.cxx:149` |

---

## 1. Scope

This review audits:

| Prototype | Toolkit | Toolkit file:line |
|---|---|---|
| `examine_structure_final(cluster)` | `examine_structure_final(graph, main_vertex, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:3492` |
| `examine_structure_final_1(cluster)` | `examine_structure_final_1(graph, main_vertex, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:2627` |
| `examine_structure_final_1p(cluster)` | `examine_structure_final_1p(graph, main_vertex, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:2759` |
| `examine_structure_final_2(cluster)` | `examine_structure_final_2(graph, main_vertex, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:3008` |
| `examine_structure_final_3(cluster)` | `examine_structure_final_3(graph, main_vertex, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:3197` |
| `search_for_vertex_activities(vtx, sg_set, cluster, range)` | `search_for_vertex_activities(graph, vtx, segments_set, cluster, track_fitter, dv, range)` | `NeutrinoVertexFinder.cxx:22` |

Five review goals:
1. Functional equivalence with the prototype.
2. Bugs in the toolkit.
3. Algorithmic efficiency.
4. Determinism (pointer-keyed / unordered containers).
5. Multi-APA / multi-face correctness.

---

## 2. Per-function walk-through

### `examine_structure_final` (ESF)

**Prototype** (`NeutrinoID_final_structure.h:2–13`): calls `examine_structure_final_1`, then `_1p`, then `_2`, then `_3` in sequence. Always returns `true`. Does **not** call `do_multi_tracking` at this level; each sub-function does so internally when it makes changes.

**Toolkit** (`NeutrinoStructureExaminer.cxx:3492–3498`): identical control flow. Returns `true`. Equivalent.

---

### `examine_structure_final_1` (ESF1)

**Purpose:** For each degree-2 vertex (excluding `main_vertex`): if its two segments are duplicate (same endpoints within 0.1 cm), remove one; otherwise, if a straight line between the far endpoints passes a good-point scan (≤1 bad point at 0.6 cm steps, 0.2 cm tolerance), merge the two segments into one new straight segment.

**Iteration:** Toolkit uses `ordered_nodes(graph)` (line 2646) — deterministic. Prototype iterates `map_vertex_segments` (a `std::map<ProtoVertex*, ...>`) — pointer-address order, non-deterministic. **Improvement.**

**Duplicate check:** Toolkit compares wcpt distances (`< 0.1 cm`) in all four endpoint pairing combinations (lines 2675–2688). Prototype compares `.index` equality. The distance-based check is more robust.

**Good-point check:** Toolkit resolves `dv->contained_by(test_p)` → `transform->backward()` → `grouping->is_good_point(apa, face, 0.2 cm)`. Points with `face==-1` (outside all TPCs) are silently skipped and not counted as bad. Prototype uses single-TPC `ct_point_cloud->is_good_point(test_p, 0.2 cm, 0, 0)`, which returns `false` for out-of-volume points (thus counting them as bad). **The toolkit is more permissive for inter-TPC gaps — reasonable for multi-TPC.**

**Merge:** Toolkit delegates to `merge_two_segments_into_one(graph, sg1, vtx, sg2, dv)` helper. Prototype builds the new segment inline. Functionally equivalent.

**`do_multi_tracking`:** Both call it after any update with `(true, true, ...)`. Equivalent.

---

### `examine_structure_final_1p` (ESF1p)

**Purpose:** When `main_vertex` has exactly 2 segments that are nearly collinear (angle > 175°), and one is short (< 6 cm and shorter than the other), merge the short segment into the longer one, relocating `main_vertex` to the short segment's far endpoint.

**Degree check:** Both prototype and toolkit require exactly 2 connections at main_vertex. Equivalent.

**Angle computation:** Prototype uses ROOT `TVector3::Angle()`. Toolkit uses `std::acos(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()))`. Equivalent.

**Endpoint orientation:** Prototype uses `.index` equality (`vec_wcps.front().index == main_vertex->get_wcpt().index`). Toolkit uses distance threshold `< 0.01 cm` (lines 2815, 2818). Reasonable adaptation — toolkit WCPoints are not guaranteed to share index semantics.

**Four orientation cases** (ff / fb / bf / bb): Faithfully reproduced. Same deduplication guard on push (`> 0.01 cm`).

**Post-merge updates:**
- Toolkit calls `create_segment_point_cloud(sg, merged_pts, dv, "main")` after updating `sg->wcpts()`. Prototype does not rebuild the point cloud (relies on `do_multi_tracking` to refit). **Improvement.**
- Toolkit sets `main_vertex->fit(vtx->fit())` when valid (lines 2867–2869). Prototype only sets `wcpt`. **Improvement.**

**Reconnection:** Toolkit collects vtx's other segments into a vector before re-adding to avoid invalidating iterators during graph mutation (lines 2872–2890). Correct.

**`do_multi_tracking`:** Called after update with `(true, true, false, false, false, &cluster)`. Prototype calls with `(true, true, true)` (includes `flag_exclusion=true`). The `flag_exclusion=false` difference is codebase-wide (see `examine_structure_review.md §B.1`).

---

### `examine_structure_final_2` (ESF2)

**Purpose:** For each segment from `main_vertex` to a nearby vertex `vtx1` (distance < 2.0 cm): check that all of vtx1's other segments have good connectivity to `main_vertex`. If so, merge `vtx1` into `main_vertex`.

**Connectivity check logic:**
1. Find point on `sg1` at ~3 cm from `vtx1`.
2. Walk from that point toward `vtx1` via `is_good_point` (0.2 cm tolerance).
3. If connected, test straight-line path from `min_point` to `main_vtx_point` (0.3 cm steps, 0.2 cm tolerance). Any bad point vetoes the merge.
4. If vtx1 has exactly 2 connections and no merge so far: check if `sg` itself has any bad points or midpoints. If yes, force-merge (the connecting segment is questionable).

Faithfully reproduced in toolkit. Good-point checks use multi-TPC path (`dv->contained_by` + `transform->backward` + `grouping->is_good_point`). Outside-TPC points skipped.

**Merge:** `merge_vertex_into_another(graph, vtx1, main_vertex, dv)` helper. Prototype does the WCPt splicing inline. Functionally equivalent.

**`do_multi_tracking`:** Equivalent.

**§B.1 — `flag_start` from wrong segment (fixed):** See §B.1 for full analysis. The toolkit now correctly uses `sg1_wcpts.front().point` vs `vtx1_wcpt_point` (line 3066).

---

### `examine_structure_final_3` (ESF3)

**Purpose:** Reverse of ESF2. Merges `main_vertex` into a nearby vertex `vtx1` (distance < 2.5 cm). All of `main_vertex`'s other segments must have good connectivity to `vtx1`. On merge, segment wcpts near `main_vertex` are rebuilt to start from `vtx1`'s position.

**Connectivity check:** same structure as ESF2 but thresholds differ — 2.5 cm distance limit (vs 2.0 cm), 0.3 cm good-point tolerance on the straight-line path (vs 0.2 cm). Both match the prototype.

**§B.1 — `flag_start` from wrong segment (fixed):** Same bug as ESF2. Fixed similarly: `sg1_wcpts.front().point` vs `main_wcpt_point` (line 3253).

**Wcpt rebuilding:** After the merge decision, for each non-`sg` segment of `main_vertex`, the toolkit rebuilds the near-vtx1 portion of wcpts by interpolating through the steiner point cloud at 1 cm steps (lines 3373–3398). Then calls `create_segment_point_cloud()` to update the point cloud. Prototype does the same inline (lines 126–168). **Toolkit also calls `create_segment_point_cloud()` — prototype does not. Improvement.**

**Removal order:** Prototype removes `vtx1` first then `sg` (lines 185–186). Toolkit removes `sg` first then `vtx1` (lines 3472–3473) with an explanatory comment: removing `vtx1` first via `boost::remove_vertex` would silently free the `sg` edge descriptor, causing a double-free in the subsequent `remove_segment(sg)`. **Correct fix over the prototype.**

**`main_vertex` update:** Toolkit sets both `wcpt` and `fit` from `vtx1` when valid (lines 3439–3442). Prototype sets only `wcpt`. **Improvement.**

**`do_multi_tracking`:** Equivalent.

---

### `search_for_vertex_activities` (SfVA)

**Purpose:** Search for unreconstructed track activity near a given vertex. Finds steiner terminal points within `search_range` that are far from all existing segments. If found, create a new short segment from the vertex to the best candidate.

**Two-round search:**

**Round 1** (strict): candidates with `min_dis > 0.6 cm` AND `sum_2d > 1.2 cm`. Scored by `sum_angle * (sum_charge + 1e-9)`.

**Round 2** (relaxed, only if Round 1 yields nothing): candidates with `min_dis > 0.36 cm` AND `sum_2d > 0.8 cm` AND `sum_charge > 20000`. Scored by `min_dis + sum_2d / sqrt(3)`.

All thresholds faithfully ported.

**Segment iteration scope:**
- Prototype iterates `map_segment_vertices` (all segments in the whole neutrino event, not just current cluster). This is a prototype quirk — it uses the global segment map for distance minimization, meaning activity near other clusters' segments is also suppressed. The toolkit filters by `sg->cluster() != &cluster` to match only the current cluster. This is a **correct improvement** for the multi-cluster toolkit architecture.

**Round 1 iteration:** `ordered_edges(graph)` (line 84). Deterministic. Correct.

**Round 2 iteration (B.2 — fixed):** Was `boost::edges(graph).first / .second` (line 149) — non-deterministic. Fixed to `ordered_edges(graph)` matching Round 1. The min-distance computation is floating-point; with different edge orders, two nearly-equal distances could resolve differently across runs. Since Round 2 already represents a degraded fallback, non-determinism here is particularly undesirable.

**Identity check:** Prototype skips candidates whose `.index == vtx->get_wcpt().index`. Toolkit skips candidates within `0.01 cm` of `vertex->wcpt().point` (line 70). Equivalent.

**Charge query:**
- Prototype: `ct_point_cloud->get_closest_dead_chs(test_p, 0)` — checks exact channel (range = 0 wires).
- Toolkit: `grouping->get_closest_dead_chs(test_p_raw, 1, apa, face, plane)` — checks ±1 wire range.
The wider range means more points are classified as near dead channels. This may slightly reduce `sum_charge` by excluding a plane from the average (larger `ncount` denominator). This is a deliberate multi-TPC adaptation (the raw-coordinate `test_p_raw` is used; a ±1 wire buffer guards against boundary effects from the coordinate transform).

**New segment creation:** Toolkit uses `make_vertex`, `create_segment_for_cluster`, `add_segment`. The 3-point path (vertex → midpoint via steiner kNN → candidate) matches the prototype. Functionally equivalent.

---

## 3. Bugs

### B.1 — `flag_start` from wrong segment (ESF2 line 3066, ESF3 line 3253)

**Status:** ✅ Fixed.

The prototype (and original toolkit) set `flag_start` by checking whether vtx1 is at the front of `sg` (the connecting segment between main_vertex and vtx1), then used it to control the walk direction through `sg1` (a completely different segment). Since `sg` and `sg1` are independent, their orientations relative to the shared vertex are unrelated.

**Impact analysis:**

`flag_start` determines the walk direction from `min_index` (the point on `sg1` closest to ~3 cm from the pivot vertex) toward the pivot vertex. If wrong, the walk goes in the opposite direction.

The `flag_connect / flag_update` interaction matters:
- `flag_connect=false` → straight-path check **skipped** → `flag_update` left unchanged (stays `true`)
- `flag_connect=true` + bad straight path → `flag_update=false` (merge vetoed)

This gives two real failure modes:

**False positive merge** (wrong direction finds good far-end, misses gapped near-vtx1 portion):
1. vtx1 at front of sg1, but flag_start=false → walk toward end (away from vtx1)
2. Far end is clean → `flag_connect=true`
3. Straight path from min_point to main_vertex is clean → `flag_update` stays true → **merge proceeds**
4. Correct walk (toward front=vtx1) would have found a gap → `flag_connect=false` → straight path skipped → `flag_update` unchanged — same outcome? No: correct case gives `flag_connect=false` so veto skipped too. But consider a segment that is gapped near vtx1 AND the straight-path to main_vertex is blocked:
   - Wrong code: far end good → `flag_connect=true` → straight path checked → bad → **merge vetoed** ✓ (happens to be correct)
   - OR: far end good → flag_connect=true → straight path clear → merge proceeds (false positive) if the straight path happens to be clear.

**False negative veto** (wrong direction finds gapped far-end, skips checking near-vtx1 which is good):
1. vtx1 at front of sg1, flag_start=false → walk toward end (away from vtx1)  
2. Far end has a gap → `flag_connect=false` → straight-path check **skipped** → flag_update unchanged (true)
3. Correct walk (toward front=vtx1) would find clean near-vtx1 → `flag_connect=true` → straight path checked → if bad → `flag_update=false` (merge vetoed)
4. Result: correct code vetoes, wrong code allows — **false positive merge**.

**Fix:** use `sg1_wcpts.front().point` (not `sg_wcpts.front().point`) for the orientation check. ESF2: line 3066; ESF3: line 3253. The old `const auto& sg_wcpts = sg->wcpts()` binding and its empty-check guard were also removed since they are no longer needed for this purpose.

---

### B.2 — `search_for_vertex_activities` Round 2 non-determinism (line 149)

**Status:** ✅ Fixed.

Round 2 was iterating edges with `boost::edges(graph).first / .second`, which does not guarantee a consistent order across runs. Changed to `ordered_edges(graph)`, consistent with Round 1 (line 84). The minimum-distance updates (`min_dis`, `min_dis_u/v/w`) are not affected by iteration order in exact arithmetic, but floating-point tie-breaking at segment boundaries depends on evaluation order.

---

## 4. Efficiency

- **ESF1:** `merge_two_segments_into_one()` encapsulates steiner kNN lookups in a single helper call. No redundant computation compared to prototype.
- **ESF2 / ESF3:** The `while(flag_continue)` loop recomputes `out_edges` on every restart. This is unavoidable since the graph mutates on each merge.
- **SfVA:** Both rounds iterate all graph edges per candidate — O(candidates × edges). Pre-computing segment distance data into a cache would help but is also true of the prototype. Not a regression.
- **ESF1p:** The two branches (sg1-short vs sg2-short) are symmetric copies. Could be refactored into a helper, but this is cosmetic and not a correctness issue.

---

## 5. Determinism

| Function | Container / iteration | Status |
|---|---|---|
| ESF1 outer | `ordered_nodes(graph)` | ✅ OK |
| ESF1 inner edge fetch | `boost::out_edges` (always exactly 2 edges, order-stable) | ✅ OK |
| ESF1p edge fetch | `boost::out_edges` (exactly 2 edges; both branches symmetric) | OK (low risk) |
| ESF2 outer | `boost::out_edges(main_vd, graph)` — breaks after first match, restarts | Low risk (documented) |
| ESF2 inner | `boost::out_edges(vtx1_vd, graph)` — examines all, no tie-breaking needed | OK |
| ESF3 outer | `boost::out_edges(main_vd, graph)` — same as ESF2 | Low risk (documented) |
| ESF3 inner | `boost::out_edges(main_vd, graph)` — examines all | OK |
| SfVA Round 1 | `ordered_edges(graph)` | ✅ OK |
| SfVA Round 2 | ~~`boost::edges(graph)`~~ → `ordered_edges(graph)` | ✅ Fixed (B.2) |

**Note on ESF2/ESF3 outer loop:** The `break`-and-restart pattern means only one merge per outer-loop pass. Non-deterministic edge order means a different eligible segment may be selected first on each run. However, the `while(flag_continue)` loop guarantees all eligible merges eventually happen (unless they are mutually exclusive — which is a general correctness limitation of the sequential greedy approach, not a determinism bug per se). Risk is low and matches the prototype's behavior exactly.

---

## 6. Multi-APA / multi-face correctness

All six functions correctly replace the prototype's single-TPC `ct_point_cloud->is_good_point(p, tol, 0, 0)` with the three-step multi-TPC pattern:
1. `dv->contained_by(test_p)` → get `(face, apa)` for the point's detector volume.
2. `transform->backward(test_p, cluster_t0, face, apa)` → raw detector coordinates.
3. `grouping->is_good_point(raw_p, apa, face, tol, 0, 0)` → signal presence check.

Points outside all TPCs (`face==-1` or `apa==-1`) are silently skipped (not counted as bad) across all six functions. This is the correct behavior for inter-APA gaps in multi-TPC geometry.

`search_for_vertex_activities` additionally:
- Skips terminal steiner point candidates entirely if they are outside all TPCs (line 82, 147).
- Resolves the candidate point's APA/face once before querying charge per plane (line 111, 166).
- Passes `(apa, face)` to `segment_get_closest_2d_distances` so 2D distances are computed in the correct wire plane coordinate system for the candidate point's TPC.

---

## 7. Improvements over prototype

| # | Description | Functions |
|---|---|---|
| I.1 | `ordered_nodes(graph)` replaces pointer-order `map_vertex_segments` iteration | ESF1 |
| I.2 | `merge_two_segments_into_one()` helper replaces inline steiner-cloud wcpt building | ESF1 |
| I.3 | Duplicate segment check (distance-based, 4 pairing combinations) | ESF1 |
| I.4 | `create_segment_point_cloud()` called after wcpt mutation | ESF1p, ESF3 |
| I.5 | `main_vertex->fit()` updated from vtx when valid after relocate | ESF1p, ESF3 |
| I.6 | `remove_segment` before `remove_vertex` to avoid double-free on boost graph | ESF3 |
| I.7 | `ordered_edges(graph)` in Round 1 of SfVA (also now Round 2 via B.2 fix) | SfVA |
| I.8 | Segment distance queries filtered to current cluster only (`sg->cluster() != &cluster`) | SfVA |
| I.9 | Multi-TPC `dv->contained_by` + `transform->backward` + `grouping->is_good_point` throughout | All |
