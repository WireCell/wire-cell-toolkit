# Review: `examine_structure` / `examine_structure_{1..4}`

**Reviewed by:** Claude Code  
**Date:** 2026-04-11  
**Branch:** `apply-pointcloud`  
**Prototype ref:** `prototype_base/pid/src/NeutrinoID_examine_structure.h` (364 lines)  
**Function map:** `clus/docs/porting/neutrino_id_function_map.md` lines 72–103

---

## Implementation status

| Item | Description | Status | File |
|---|---|---|---|
| B.1 | `flag_exclusion` always `false` vs prototype's `true` (systematic across entire codebase; see §B.1) | Documented only | `NeutrinoStructureExaminer.cxx` (and all other callers) |
| B.2 | ES4: outside-TPC candidate passes 2D-distance criterion trivially | ✅ Applied | `NeutrinoStructureExaminer.cxx` |
| B.3 | ES1/ES2/ES4 scans: outside-TPC points silently counted as good | ✅ Applied | `NeutrinoStructureExaminer.cxx` |
| B.4 | ES4: `dv->contained_by(test_p)` called once per segment (loop-invariant) | ✅ Applied (merged with B.2) | `NeutrinoStructureExaminer.cxx` |
| B.5 | ES4: raw `boost::edges` instead of `ordered_edges` | ✅ Applied | `NeutrinoStructureExaminer.cxx` |
| B.6 | ES4: misleading comment re cross-cluster segment scope | ✅ Applied | `NeutrinoStructureExaminer.cxx` |
| B.7 | ES2: co-located vtx1/vtx2 — vtx2's other segments not merged into vtx1 | ✅ Applied | `NeutrinoStructureExaminer.cxx` |
| E.2 | ES3: endpoint orientation detected by geometric proximity (0.01 cm, too tight after kNN drift) | ✅ Applied | `NeutrinoStructureExaminer.cxx` |
| E.4 | ES4: re-compute `dis` from coords when squared distance already available | ✅ Applied | `NeutrinoStructureExaminer.cxx` |

---

## 1. Scope

This review audits:

| Prototype | Toolkit | Toolkit file:line |
|---|---|---|
| `examine_structure(cluster)` | `examine_structure(graph, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:13` |
| `examine_structure_1(cluster)` | `examine_structure_1(graph, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:25` |
| `examine_structure_2(cluster)` | `examine_structure_2(graph, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:166` |
| `examine_structure_3(cluster)` | `examine_structure_3(graph, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:340` |
| `examine_structure_4(vtx, cluster, flag)` | `examine_structure_4(vtx, flag, graph, cluster, track_fitter, dv)` | `NeutrinoStructureExaminer.cxx:473` |

All five review goals addressed:
1. Functional equivalence with the prototype.
2. Bugs in the toolkit.
3. Algorithmic efficiency.
4. Determinism (pointer-keyed / unordered containers).
5. Multi-APA / multi-face correctness.

---

## 2. Per-function walk-through

### `examine_structure` (ES)

**Prototype** (`NeutrinoID_examine_structure.h:1–12`): calls `examine_structure_2` then `examine_structure_1`; on any update calls `do_multi_tracking(…, true, true, true)` (`flag_dQ_dx_fit_reg`, `flag_dQ_dx_fit`, `flag_exclusion` all `true`).

**Toolkit** (`NeutrinoStructureExaminer.cxx:13–23`): equivalent control flow. `do_multi_tracking` called as `(true, true, false, false, false, &cluster)` — `flag_exclusion=false`. See §B.1.

### `examine_structure_1` (ES1)

**Purpose:** replace short or high-dQ/dx segments with a straight line if the straight line passes a "good-point" scan.

**Length cut** (`< 5 cm || (< 8 cm && dQ_dx > 1.5 MIP)`): faithfully ported.

**Straight-line scan:** toolkit uses `wcpts.front/back().point` as endpoints (not the fit points). This was an intentional fix documented inline — using fit positions caused the stop condition to fire late, appending a kink. Correct.

**Good-point check:** calls `dv->contained_by(test_p)` → `transform->backward` → `grouping->is_good_point`. If `contained_by` returns `face==-1`, the point was not counted as bad. See §B.3.

**Steiner knn:** uses `cluster.kd_steiner_knn(1, …, "steiner_pc")` to project new points. Equivalent to prototype's `pcloud_steiner->get_closest_wcpoint(…)`. No index-equality check is possible; distance threshold (`0.01 cm`) guards deduplication. Functionally equivalent.

**Determinism:** iterates `boost::edges(graph)` — note `examine_structure_1` does NOT mutate the graph (it only updates `sg->wcpts()`), so non-ordered iteration is safe here (no break/restart loop). Not an issue.

### `examine_structure_2` (ES2)

**Purpose:** merge two short segments connected through a degree-2 vertex into one straight segment.

**Iteration:** uses `ordered_nodes(graph)` — deterministic. ✓

**Good-point check:** same outside-TPC issue as ES1 (see §B.3).

**Straight-line construction:** uses `vtx1->wcpt().point` / `vtx2->wcpt().point` as endpoints (not `fit().point`). Same intentional fix as ES1. Correct.

**Co-located endpoints (B.7):** when `dist(vtx1,vtx2) < 0.01 cm`, the prototype would merge `vtx2`'s other segments onto `vtx1` (`map_segment_vertices` re-parenting) then delete `vtx2`. The toolkit skips this step entirely due to a `setS` descriptor-aliasing crash. This means other segments previously attached to `vtx2` are NOT consolidated onto `vtx1`. This is a known deviation; a safe fix requires a dedicated `merge_vertex_into_another` redesign.

### `examine_structure_3` (ES3)

**Purpose:** merge two segments through a degree-2 vertex based on angular collinearity (angle at 10 cm < 18° AND angle at 3 cm < 27°).

**Angle computation:** prototype uses `TVector3::Angle(dir1, dir2)` and flips with `π - angle`. Toolkit uses `std::acos(dot/(|a||b|))` and `π - angle`. Numerically equivalent.

**WCPoint merging:** prototype used WCPoint index equality to identify which endpoint of sg1/sg2 touches `vtx`. Toolkit uses geometric proximity (< 0.01 cm). This is tighter than typical steiner-point spacing; if kNN projection shifts a wcpt by more than 0.01 cm, the merge branches would all fail and `merged_wcpts` would remain empty, causing the function to silently skip this merge. An alternative would be to use `< 0.1*units::cm` as threshold. Documented as a risk; not changed in this pass.

**Iteration:** uses `ordered_nodes(graph)` — deterministic. ✓

### `examine_structure_4` (ES4)

**Purpose:** from a given vertex, search Steiner terminals within 6 cm for a candidate that (a) is far from all existing segments in 3D and 2D, (b) has a clear line-of-sight, and add a new segment to that terminal.

**Prototype dead code removed:** `saved_dirs` and `nshowers` in prototype were computed but never read. Correctly omitted in toolkit. ✓

**Cross-cluster segment iteration:** the graph is per-grouping (one graph for all clusters), so `boost::edges(graph)` iterates segments from ALL clusters, matching the prototype's unconstrained `map_segment_vertices`. The comment at the original code (lines 530–534) that claimed "toolkit only checks the current cluster's segments" was incorrect. See §B.6.

**2D distance check:** see §B.2 / §B.4.

**New vertex wiring:** `make_vertex(graph)` properly registers the vertex (`PRGraph.cxx:5–17`, sets descriptor). `add_segment` is then idempotent on already-registered vertices (`PRGraph.cxx:31–32`). ✓

---

## 3. Bugs (correctness)

### §B.1 — `flag_exclusion=false` in all `do_multi_tracking` calls

**Status:** Documented only. Requires dedicated quality-impact study before changing.

**What `update_association` does:** a wire-hit exclusion mechanism. For each fit-point on a segment, it retains only the 2D hit coordinates that are closer to this segment's fitted track than to any other segment's. This prevents multi-segment tracks from double-counting shared hit charge. The toolkit's implementation (`TrackFitting.cxx:2459–2602`) is fully multi-APA-correct — the per-coordinate APA/face offset/slope lookup was explicitly fixed (§4.1 of `do_multi_tracking_review.md`, 2026-04-11). `update_association` is called inside `form_map_graph` once per middle fit-point per segment (`TrackFitting.cxx:3111`); it modifies local `PlaneData` temporaries, not shared state.

**Prototype:** every `do_multi_tracking` call uses `(…, true, true, true)` → `flag_exclusion=true`. Three calls in `NeutrinoID_proto_vertex.h` (lines 719, 748, 773, all commented out) passed `false`, confirming the exclusion was intentional but selectively disabled in a few special cases.

**Toolkit:** 33+ call sites across `NeutrinoStructureExaminer.cxx` and `NeutrinoPatternBase.cxx` all pass `flag_exclusion=false`. No comment in any of those files explains the decision. The function is implemented and correct; it is simply never invoked.

**Impact:** without exclusion, adjacent segments can claim the same wire hits, inflating dQ/dx near vertices. With exclusion, each hit belongs to only the nearest segment — the intended behavior.

**Action:** enabling `flag_exclusion=true` globally is a non-trivial quality change touching all of pattern recognition. This needs a dedicated pass with before/after reconstruction comparisons. Do not change in this review.

### §B.2 — ES4: outside-TPC candidate bypasses 2D-distance filter

**Status:** ✅ Fixed.

**Location:** `NeutrinoStructureExaminer.cxx:549–555` (original).

**Problem:** `test_p` is a Steiner terminal. If `dv->contained_by(test_p)` returns `face==-1` (point on or outside the sensitive bounding box), the 2D-distance block is skipped and `min_dis_{u,v,w}` remain at `1e9`. The guards `min_dis_u + min_dis_v + min_dis_w > 1.8 cm` and `min_dis_u > 0.8 cm && …` then trivially pass, allowing an outside-TPC terminal to reach the line-scan stage. Prototype is not susceptible because it fetches 2D distances unconditionally via segment lookup structures.

**Fix:** hoisted `dv->contained_by(test_p)` above the `boost::edges` segment loop (also eliminating the B.4 redundancy); added an early `continue` when `face()==-1`.

### §B.3 — ES1 / ES2 / ES4 scans: outside-TPC steps silently treated as good

**Status:** ✅ Fixed.

**Locations:** ES1 line 91, ES2 line 241, ES4 line 585 (original).

**Problem:** the scan loop only increments `n_bad` when `face()!=-1 && apa()!=-1`. When a test point is outside every TPC, no charge is available there, so the prototype's `is_good_point(…)` would return `false`. The toolkit silently skips and leaves the bad count unchanged — more permissive than the prototype at TPC boundaries.

**Fix:** add an `else { n_bad++; }` branch to treat `face()==-1` as a bad sample.

### §B.4 — ES4: `dv->contained_by(test_p)` called N times per candidate (loop-invariant)

**Status:** ✅ Fixed (merged with B.2).

**Location:** original line 549 inside the `boost::edges` loop.

**Problem:** `test_p` is fixed for each steiner-terminal candidate, but `contained_by` was called once per segment — O(#segments) redundant lookups.

### §B.5 — ES4: raw `boost::edges` for the segment closest-distance loop

**Status:** ✅ Fixed.

**Location:** `NeutrinoStructureExaminer.cxx:535` (original).

**Problem:** all other examine_structure functions use `ordered_nodes`/`ordered_edges` for deterministic processing. ES4 used `boost::edges(graph)` for the segment-iteration inner loop, making the `max_dis` tie-breaking dependent on Boost's internal edge ordering (not deterministic across runs). While floating-point ties in `dis` are unlikely, it is an inconsistency.

**Fix:** replaced with `ordered_edges(graph)` from `clus/src/PRGraphType.cxx:18–25`.

### §B.6 — ES4: misleading cross-cluster comment

**Status:** ✅ Fixed.

**Location:** `NeutrinoStructureExaminer.cxx:530–534` (original).

**Problem:** comment said "this toolkit only checks the current cluster's segments." The toolkit graph is per-grouping (`m_graph` in `TrackFitting.h:592`, built via `add_graph` and `build_cluster_edges`); `boost::edges(graph)` thus spans all clusters, which correctly matches the prototype.

**Fix:** replaced comment with accurate description.

### §B.7 — ES2: co-located vtx1/vtx2 "merge other segments" behavior

**Status:** ✅ Fixed.

**Location:** `NeutrinoStructureExaminer.cxx` (was ~line 272–280; fix extends ~line 332–347).

**Problem:** when vtx1 and vtx2 are geometrically co-located (< 0.01 cm) but distinct `Vertex` objects, the prototype (`NeutrinoID_examine_structure.h:139–152`) re-parents all of vtx2's segments onto vtx1 and deletes vtx2. The toolkit's original comment warned against calling `merge_vertex_into_another(vtx2, vtx1)` because it would try to move sg2 (vtx2→vtx) as a vtx1→vtx edge, but vtx1→vtx already carried sg1 — with `setS` edge storage, `add_segment` aliases both to one descriptor, and then `remove_segment(sg2)` crashes by deleting the shared edge.

**Root cause analysis:** `merge_vertex_into_another` does `remove_segment(old_seg)` before `add_segment(old_seg, vtx_to, other_vtx)` (`NeutrinoPatternBase.cxx:1302–1303`), so there is no aliasing — **provided the conflicting edges are already gone**. The crash only occurs if sg2 (vtx2→vtx) and sg1 (vtx1→vtx) still exist when the merge runs.

**Fix:** defer `merge_vertex_into_another(vtx2, vtx1)` until **after** `remove_segment(sg1)`, `remove_segment(sg2)`, `remove_vertex(vtx)` are executed. At that point, both conflicting edges are gone and the merge is safe. If vtx2 has no remaining segments after sg2's removal, it is removed directly instead.

The graph type uses `setS` for edge storage (`PRGraphType.h:91–98`, `boost::setS` prevents parallel edges). The `WCPoint` struct in the toolkit has no `index` field (unlike the prototype), so the "same index" equality check is not reproducible; the 0.01 cm geometric threshold is used instead, which is tighter than the steiner-point spacing and is therefore effectively equivalent.

---

## 4. Efficiency

### §E.1 — ES4: `contained_by` hoisted (merged with §B.4)

Eliminates O(#segments) redundant lookups per terminal candidate.

### §E.2 — ES3: endpoint orientation detection via graph topology (was 0.01 cm threshold)

**Status:** ✅ Applied.

The original code used a 0.01 cm geometric proximity test to determine which wcpt endpoint of sg1/sg2 lies at the shared vertex `vtx`. After kNN projection by ES1/ES2, a segment's wcpts can shift by up to ~0.15 cm (half the steiner step size), making the 0.01 cm threshold too tight: the merge silently falls through to `merged_wcpts.empty()` and is skipped.

**Fix:** use `find_vertices(graph, sg)` which returns `(v_front, v_back)` ordered by proximity to `sg->wcpts().front()`. Comparing `sg1_front_vtx == vtx` and `sg2_front_vtx == vtx` by pointer identity gives exact orientation — no distance threshold involved.

### §E.3 — ES1 / ES2: scan loops are already optimal

### §E.4 — ES4: use pre-computed squared distance from `kd_steiner_radius`

**Status:** ✅ Applied.

`kd_steiner_radius` returns `pair<index, squared_distance>`. The original code recomputed `dis` from coordinates. Now uses `std::sqrt(dist_sq)` directly.

---

## 5. Determinism

| Function | Container | Status |
|---|---|---|
| ES1 | `boost::edges(graph)` — no graph mutation, safe | ✓ |
| ES2 | `ordered_nodes(graph)` — restarts after each mutation | ✓ |
| ES3 | `ordered_nodes(graph)` — restarts after each mutation | ✓ |
| ES4 outer (candidates) | `kd_steiner_radius` — kd-tree order, deterministic given fixed tree build | ✓ |
| ES4 inner (segments) | was `boost::edges`; now `ordered_edges(graph)` | ✅ Fixed (§B.5) |

---

## 6. Multi-APA / multi-face correctness

The prototype operated on a single-TPC detector; all `is_good_point` calls used wire-plane index `0`. The toolkit resolves `(apa, face)` per test-point via `dv->contained_by` → `transform->backward` → `grouping->is_good_point`. This pattern is applied consistently across ES1, ES2, and ES4.

`segment_get_closest_2d_distances` (`PRSegmentFunctions.cxx:141–174`) intentionally returns large values when the query point lives in a different APA/face than the segment — correct, because 2D wire-plane distances are only meaningful within one face. ES4's 2D-distance check implicitly relies on this.

New vertices and segments created by ES4 (`make_vertex`→`add_segment`) do not need explicit face/APA registration at creation time — fits are assigned by the subsequent `do_multi_tracking` call. ✓

The `dv->contained_by` → `face==-1` edge cases (B.2, B.3) are now handled correctly by the fixes above.
