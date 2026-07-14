# do_multi_tracking Port Review

**Reviewed:** 2026-04-11  
**Fixes applied (2026-04-11):** §2.2 (FIXED — pixel-frame projection in Gaussian weighting), §2.4 (FIXED — re-query wpid after fit_point for vertex projection), §4.1 (FIXED — pass pre-built segment list into update_association), §4.2 (FIXED — `FaceGeom` cache in `multi_trajectory_fit` avoids per-vertex double map lookup). §2.1 deferred (missing port). §2.3 NOT A BUG. §4.3 NOT A BUG (third `form_map_graph` call genuinely needed after `reset_fit_prop` + `organize_segments_path_3rd`). §5.1 NOT APPLICABLE (U/V slopes are non-zero in all real detector geometries; W already uses correct slope).  
**Entry point:** `TrackFitting::do_multi_tracking` — `clus/src/TrackFitting.cxx:7423`  
**Prototype entry:** `WCPPID::PR3DCluster::do_multi_tracking` — `prototype-dev/wire-cell/pid/src/PR3DCluster_multi_track_fitting.h:2`  
**Scope:** Functional equivalence, bugs, efficiency, determinism, multi-APA/face correctness.  
**Method:** Direct line-by-line comparison of all multi-tracking-specific functions against the prototype, plus cross-referencing the single-tracking review for shared helpers.  
All five review goals addressed:
1. Functional identity with prototype.
2. Bugs in toolkit.
3. Algorithmic efficiency.
4. Determinism (pointer-keyed containers).
5. Multi-APA / multi-face correctness.

**Key reference notes:**
- Toolkit blob wire ranges are half-open `[min, max)`. The `<` vs prototype's `<=` is intentional — see `do_single_tracking_review.md`.
- Prototype multi-tracking is `PR3DCluster_multi_track_fitting.h` (trajectory) + `PR3DCluster_multi_dQ_dx_fit.h` (dQ/dx). The toolkit merges both into `TrackFitting`.
- The prototype's `form_map_multi_segments` became `form_map_graph` in the toolkit (shared between single and multi).
- Prototype `map_vertex_segments` / `map_segment_vertices` are pointer-keyed maps; the toolkit equivalent is the Boost graph `m_graph` with `m_ordered_nodes_vec` for stable vertex iteration.

**Cross-reference:** The following shared helpers were thoroughly reviewed in `do_single_tracking_review.md` and are only summarised here with new multi-specific concerns:

| Shared helper | Single review sections | Multi-specific note |
|---|---|---|
| `prepare_data` | §2.8, §2.9 | No new issues |
| `form_point_association` | §2.4–2.6, §5.1 | Multi-face fix already applied |
| `examine_point_association` | §2.17, §5.5 | No new issues |
| `skip_trajectory_point` | §2.1, §2.2 | §2.1 index bug affects `examine_segment_trajectory` — see §2.3 below |
| `calculate_compact_matrix_multi` | §3.D1, §4.5 | Stable int key fix applied |
| `dQ_dx_multi_fit` | §2.10, §2.11, §5.2 | No new multi-specific issues beyond what is documented there |

---

## 1. Function Map (prototype → toolkit)

| Prototype function | Prototype file:line | Toolkit function | Toolkit file:line | Notes |
|---|---|---|---|---|
| `do_multi_tracking` | `PR3DCluster_multi_track_fitting.h:2` | `do_multi_tracking` | `TrackFitting.cxx:7423` | Missing `collect_charge_multi_trajectory` call — see §2.1 |
| `collect_charge_multi_trajectory` | `PR3DCluster_multi_track_fitting.h:1721` | **MISSING** | — | Not ported — see §2.1 |
| `organize_segments_path` | `PR3DCluster_multi_track_fitting.h:~600` | `organize_segments_path` | `TrackFitting.cxx:1551` | Ported cleanly; multi-APA uses `generate_fits_with_projections` |
| `organize_segments_path_2nd` | `PR3DCluster_multi_track_fitting.h:1324` | `organize_segments_path_2nd` | `TrackFitting.cxx:1382` | Toolkit adds degenerate-collapse guard — improvement |
| `organize_segments_path_3rd` | `PR3DCluster_multi_track_fitting.h:1100` | `organize_segments_path_3rd` | `TrackFitting.cxx:1230` | Toolkit adds `check_and_reset_close_vertices()` — improvement |
| `form_map_multi_segments` | `PR3DCluster_multi_track_fitting.h:715` | `form_map_graph` | `TrackFitting.cxx:2945` | Shared with single; stable int key fix applied |
| `multi_trajectory_fit` | `PR3DCluster_multi_track_fitting.h:207` | `multi_trajectory_fit` | `TrackFitting.cxx:3516` | Frame-mixing bug in charge_div_method==2 — see §2.2 |
| `fit_point` | `PR3DCluster_multi_track_fitting.h:548` | `fit_point` | `TrackFitting.cxx:3291` | Per-point wpid lookup correct; see §2.4 for vertex projection issue |
| `examine_trajectory` | `PR3DCluster_multi_track_fitting.h:422` | `examine_segment_trajectory` | `TrackFitting.cxx:4472` | Inherits skip_trajectory_point §2.1 index bug |
| `skip_trajectory_point` | (shared) | `skip_trajectory_point` | `TrackFitting.cxx:4632` | Cross-ref single review §2.1, §2.2 |
| `update_association` | `PR3DCluster_multi_track_fitting.h:970` | `update_association` | `TrackFitting.cxx:2426` | Multi-APA geometry correct — see §5.1 |
| `get_pos_multi` | `PR3DCluster_multi_track_fitting.h:194` | (inlined) | — | Prototype helper inlined into map_index_pss lookup |

---

## 2. Bugs (S1 — Correctness)

### 2.1 `collect_charge_multi_trajectory` not ported — DEFERRED

**Toolkit:** not present  
**Prototype:** `PR3DCluster_multi_track_fitting.h:17` calls it at the very start of `do_multi_tracking`, before `prepare_data`.

`collect_charge_multi_trajectory` (prototype `:1721`) collects charge from nearby 2D clusters that are *outside* the cluster's own blobs, storing them in `collected_charge_map`. This map is then available during `prepare_data` to include additional charge hits not already in the cluster geometry. Without this call, multi-segment tracks near APA boundaries or Bragg peaks may miss charge deposited in adjacent blobs.

The prototype also has hardcoded uBooNE volume bounds at `:1773`:

```cpp
if (pts.at(i).y < -120*units::cm || pts.at(i).y > 120*units::cm ||
    pts.at(i).z < -5*units::cm   || pts.at(i).z > 1070*units::cm) continue;
```

**Porting notes:** If this function is ported, these bounds must be replaced with a detector-agnostic active volume query (DUNE FD has different y/z extents per APA). The toolkit's `m_dv->contained_by(p)` infrastructure provides this.

**Fix:** Deferred. Same status as single review §2.16.

---

### 2.2 `multi_trajectory_fit`: frame-mixing bug in Gaussian charge division (charge_div_method==2) — FIXED

**File:** `TrackFitting.cxx:3618–3631`  
**Prototype:** `PR3DCluster_multi_track_fitting.h:269–283` (single-TPC, no apa/face ambiguity)

In the Gaussian weighting path (charge_div_method == 2), the outer loop over `coord_idx_fac` is keyed by `(apa_face, coord_idx)`, where `apa_face` is the **pixel**'s (apa, face). The code fetches the pixel's geometry:

```cpp
// outer loop: apa_face = pixel's (apa, face)
auto offset_t = std::get<0>(offset_it->second);  // pixel's offset
auto slope_x  = std::get<0>(slope_it->second);   // pixel's slope
```

Then transforms the **3D point** using the point's own wpid:

```cpp
auto test_wpid = m_dv->contained_by(pss_vec_it->second.first);
auto p_raw = transform->backward(pss_vec_it->second.first, cluster_t0,
                                  test_wpid.face(), test_wpid.apa());   // point's own face
double central_t = slope_x * p_raw.x() + offset_t;  // pixel's slope × point's raw x
```

If the 3D point lies in a face different from the pixel's face (e.g. a track crossing an APA boundary), `p_raw.x()` is in the point's local frame while `slope_x` and `offset_t` are in the pixel's frame. This produces an incorrect `central_t` and `central_ch`, biasing the Gaussian weights.

This is the same root cause as single review §2.3 (`trajectory_fit`). The fix is the same: re-project `p_raw` into the **pixel's** (apa, face) using `transform->backward(p, cluster_t0, pixel_face, pixel_apa)`, then apply the pixel's slopes and offsets.

**Fix:** Inside the `for (auto& [coord_idx, fac] : coord_idx_fac)` loop, replace:
```cpp
auto test_wpid = m_dv->contained_by(pss_vec_it->second.first);
auto p_raw = transform->backward(pss_vec_it->second.first, cluster_t0,
                                  test_wpid.face(), test_wpid.apa());
```
with:
```cpp
// Project the 3D point into the PIXEL's face for this charge-division factor
int pixel_apa = apa_face.first, pixel_face = apa_face.second;
auto p_raw = transform->backward(pss_vec_it->second.first, cluster_t0,
                                  pixel_face, pixel_apa);
```
Then use the already-loaded `slope_x`, `offset_t`, etc. (which are the pixel's) consistently.

---

### 2.3 `examine_segment_trajectory` inherits index-semantics bug from `skip_trajectory_point` — NOT A BUG (comment already added)

**File:** `TrackFitting.cxx:4500`  
**Prototype:** `PR3DCluster_multi_track_fitting.h:429` passes `init_indices.at(i)` (the form_map count) as the `index` argument of `skip_trajectory_point`, separate from `i` (the position index). The toolkit `examine_segment_trajectory` passes `i` as both:

```cpp
bool flag_skip = skip_trajectory_point(p, apa_face, i, pss_vec, fine_tracking_path);
```

As documented in single review §2.1, the toolkit's `skip_trajectory_point` uses `i` for `m_3d_to_2d.find(i)` dead-plane lookup, which works correctly only if `i` equals the form_map count index. After `multi_trajectory_fit`, the `fits` vector stored per-segment has `fit.index` set to the form_map count — but `examine_segment_trajectory` loops over `final_ps_vec`/`init_ps_vec` by position, so `i` here IS a position index, not a form_map count.

However, the `m_3d_to_2d` used by `skip_trajectory_point` is the same member map that was populated by the most-recent `form_map_graph` call, whose count values correspond to segment fit indices. In practice both `i` and `fit.index` usually agree for continuous segments, so the silent mismatch rarely fires. The issue is identical to single review §2.1 — a clarifying comment covers it.

---

### 2.4 `multi_trajectory_fit`: vertex projection uses pre-fit wpid — FIXED

**File:** `TrackFitting.cxx:3690–3730`  
**Prototype:** `PR3DCluster_multi_track_fitting.h:342` projects `init_p` (which by that point IS the fitted position after `fit_point`) using global single-TPC geometry — no wpid ambiguity.

In the toolkit vertex loop, `test_wpid` is obtained from `init_p` (the **pre-fit** position) at `:3690`, then reused to backward-transform the **fitted** `fitted_p` at `:3724`:

```cpp
auto test_wpid = m_dv->contained_by(init_p);    // ← from PRE-fit position
// ...
auto fitted_p_raw = transform->backward(fitted_p, cluster_t0,
                                         test_wpid.face(), test_wpid.apa());  // ← applied to fitted_p
vertex_fit.pu = offset_u + (slope_yu * fitted_p_raw.y() + slope_zu * fitted_p_raw.z());
```

If `fit_point` moves the vertex across an APA boundary, `test_wpid` is wrong for `fitted_p`, producing incorrect `pu/pv/pw/pt` projections stored in `vertex_fit`. The geometry scalars `offset_u`, `slope_yu`, etc. at `:3699–3710` are also looked up using the pre-fit wpid — these would also be wrong.

**Fix:** Re-query `m_dv->contained_by(fitted_p)` after `fit_point` returns, and re-look up `wpid_offsets`/`wpid_slopes` from the new wpid before computing `fitted_p_raw` and the projections.

---

## 3. Determinism Issues (S2)

### 3.D1 Vertex iteration uses `m_ordered_nodes_vec` — NOT A BUG

**File:** `TrackFitting.cxx:3661` (multi_trajectory_fit vertex loop), `TrackFitting.cxx:7436` (do_multi_tracking vertex reset)

Both loops iterate `m_ordered_nodes_vec`, which is built from `build_cluster_edges()` using a topological BFS/DFS order that is deterministic across runs for a given graph topology. This correctly replaces the prototype's `map_vertex_segments` (pointer-keyed, non-deterministic iteration order) with a stable ordering.

---

### 3.D2 `get_segment_edges()` iteration order — NOT A BUG (verify if needed)

**File:** `TrackFitting.cxx:3521` (multi_trajectory_fit segment loop), `TrackFitting.cxx:2984` (form_map_graph segment loop)

`get_segment_edges()` returns edges from `build_cluster_edges()`, which iterates `boost::edges(*m_graph)`. With Boost `adjacency_list<listS, vecS>`, edges are stored in per-vertex adjacency lists and iterated in insertion order. Edge insertion happens during graph construction (fixed per run), so iteration order is deterministic within a run. No action needed.

---

### 3.D3 `m_2d_to_3d` key ordering — NOT A BUG

**File:** `TrackFitting.cxx:3079–3086` (form_map_graph fills m_2d_to_3d)

`m_2d_to_3d` is `std::map<Coord2D, std::set<int>>`. `Coord2D` is compared by its integer fields (apa, face, time, wire, channel, plane) in lexicographic order — fully deterministic. The `std::set<int>` values are also ordered. No pointer-keyed ordering issues.

---

### 3.D4 `map_index_pss` in `multi_trajectory_fit` uses stable int keys — NOT A BUG

**File:** `TrackFitting.cxx:3520–3530`

```cpp
std::map<int, std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> map_index_pss;
for (const auto& ed : get_segment_edges()) {
    for (const auto& fit : fits) {
        int idx = fit.index;          // stable int assigned by form_map_graph
        map_index_pss[idx] = ...;
    }
}
```

Keys are `fit.index` integers from `form_map_graph`'s sequential `count`. Iteration is in deterministic int order. The equivalent in the prototype was `map_3D_tuple` (also int-keyed). No issues.

---

## 4. Efficiency Wins (S3)

### 4.1 `update_association` rebuilds segment list per call

**File:** `TrackFitting.cxx:2437–2443`

```cpp
std::vector<std::shared_ptr<PR::Segment>> all_segments;
for (const auto& ed : get_segment_edges()) {
    auto& edge_bundle = (*m_graph)[ed];
    if (edge_bundle.segment && edge_bundle.segment != segment)
        all_segments.push_back(edge_bundle.segment);
}
```

`update_association` is called once per middle fit point per segment inside `form_map_graph`. With N segments × M points each, this rebuilds an O(N) list O(N×M) times. The `get_segment_edges()` call itself is cheap (cached vector), but the vector allocation and fill is repeated.

**Fix:** Pass the pre-built `all_segments` vector (from `form_map_graph`'s own local `segments` vector at `:2972`) into `update_association` as a `const` reference parameter, avoiding the per-call rebuild. The signature would become:
```cpp
void update_association(std::shared_ptr<PR::Segment> segment,
                        const std::vector<std::shared_ptr<PR::Segment>>& all_segments,
                        PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt);
```

---

### 4.2 `multi_trajectory_fit` vertex loop re-looks up wpid/offsets/slopes per vertex

**File:** `TrackFitting.cxx:3690–3710`

Each vertex in the loop calls `m_dv->contained_by(init_p)` and then two `std::map::find` calls into `wpid_offsets` and `wpid_slopes`. For events with many vertices (N-prong interactions), this is O(N log F) where F = number of faces. Consider caching the per-face geometry in a small flat array indexed by `(apa, face)` at the top of `multi_trajectory_fit` (same pattern as already done in `trajectory_fit`'s outer loop).

---

### 4.3 `form_map_graph` called redundantly before `dQ_dx_multi_fit` — NOT A BUG

**File:** `TrackFitting.cxx:7919` (third `form_map_graph` call in `do_multi_tracking`)

In `do_multi_tracking`, `form_map_graph` is called three times: before `multi_trajectory_fit` (1st pass), before `multi_trajectory_fit` (2nd pass, after vertex/segment reset), and before `dQ_dx_multi_fit`. The third call is genuinely necessary for two reasons:

1. **`organize_segments_path_3rd` moves fit points** (redistributes them at finer 0.6 cm spacing), so the `m_3d_to_2d` / `m_2d_to_3d` maps from round 2 are stale.
2. **`reset_fit_prop()` is called on all vertices and segments** immediately before the third `form_map_graph`, clearing the index assignments from round 2. The third call re-assigns fresh sequential `fit.index` values needed by `dQ_dx_multi_fit`.

No action needed.

---

## 5. Multi-APA / Multi-Face Correctness (S5)

### 5.1 `update_association`: per-Coord2D geometry lookup — CORRECT ✓

**File:** `TrackFitting.cxx:2447–2566`

For each 2D coordinate, the toolkit reads:
```cpp
int apa = coord.apa;
int face = coord.face;
// ...
auto offset_it = wpid_offsets.find(WirePlaneId(kAllLayers, face, apa));
auto slope_it  = wpid_slopes.find(WirePlaneId(kAllLayers, face, apa));
```
and inverts the 2D coordinate to a local 3D test point using the coordinate's own (apa, face) geometry. The prototype used a single global offset/slope for all coordinates (uBooNE single-face). The toolkit is strictly better here.

One note: the reconstruction from wire index to (x, y or z) is:

- U plane: `raw_y = (coord.wire - offset_u) / slope_yu`  
- V plane: `raw_y = (coord.wire - offset_v) / slope_yv`  
- W plane: `raw_z = (coord.wire - offset_w) / slope_zw`

For U and V, z is set to 0; for W, y is set to 0. The resulting `test_point` has one coordinate absent. Distance is then computed by `segment_get_closest_2d_distances` which projects the segment points into the same reduced 2D space. This matches the prototype and is correct for this nearest-neighbour threshold decision.

**Caution:** `slope_yv` and `slope_yu` could be zero for certain detector geometries (rotated wires), causing division by zero. This is a latent issue inherited from the prototype. Consider adding a `std::abs(slope) > epsilon` guard.

---

### 5.2 `organize_segments_path` / `_2nd` / `_3rd`: endpoint extension is geometry-agnostic — CORRECT ✓

**File:** `TrackFitting.cxx:1551, 1382, 1230`

All three `organize_segments_path` variants compute endpoint extensions using pure 3D Euclidean arithmetic (no projection to wire/time). Multi-APA segments with endpoints near APA boundaries are handled correctly because no frame-specific offset is needed for the spacing computation. The toolkit calls `generate_fits_with_projections` at the end of each pass to compute per-point 2D projections using per-point `contained_by` lookups — this is the correct multi-APA generalisation.

---

### 5.3 `multi_trajectory_fit` vertex projection uses pre-fit wpid — FIXED (see §2.4)

Cross-reference §2.4. Fixed: re-query `contained_by(fitted_p)` after `fit_point` and reload geometry from the new wpid before computing `fitted_p_raw` and the projections.

---

### 5.4 `multi_trajectory_fit` Gaussian weighting mixes coordinate frames — FIXED (see §2.2)

Cross-reference §2.2.

---

### 5.5 `form_map_graph` vertex range: per-point `(apa, face)` for form_point_association — CORRECT ✓

**File:** `TrackFitting.cxx:3058`

```cpp
form_point_association(segment, fits[i].point, temp_2dut, temp_2dvt, temp_2dwt, dis_cut, nlevel, time_tick_cut);
```

`form_point_association` (single review §2.4) already uses the multi-face logic: it groups nearby blobs by `(apa, face)` and projects the 3D point into each face's frame. No extra work needed here.

---

## 6. Parameter Parity Checklist (S4)

| Parameter | Prototype value | Toolkit value | Match |
|---|---|---|---|
| `low_dis_limit` (1st pass) | `1.2*units::cm` | `m_params.low_dis_limit` = `1.2*units::cm` | ✓ |
| `end_point_limit` (1st pass) | `0.6*units::cm` | `m_params.end_point_limit` = `0.6*units::cm` | ✓ |
| `low_dis_limit` (2nd pass) | `0.6*units::cm` | `m_params.low_dis_limit / 2.` = `0.6*units::cm` | ✓ |
| `end_point_limit` (2nd pass) | `0.3*units::cm` | `m_params.end_point_limit / 2.` = `0.3*units::cm` | ✓ |
| `step_size` (3rd pass) | `0.6*units::cm` | `0.6*units::cm` (hardcoded) | ✓ |
| `charge_div_method` | 1 (equal div) | 1 (equal div) | ✓ |
| `div_sigma` | `m_params.div_sigma` | `m_params.div_sigma` | ✓ |
| `flag_exclusion` | passed in | passed in | ✓ |
| `+0.5` wire offset in pu/pv/pw | present (uBooNE convention) | absent | intentional — toolkit uses local wire indices per face |
| `+2400/+4800` plane bias in pv/pw | present (uBooNE channel numbering) | absent | intentional — uBooNE-specific channel convention |

---

## 7. Execution Order (severity priority)

| Phase | Work | Severity | Status |
|---|---|---|---|
| Now | Fix §2.2: frame-mixing in `multi_trajectory_fit` Gaussian weighting | S1 correctness + S5 multi-APA | **FIXED** |
| Now | Fix §2.4: vertex projection uses pre-fit wpid | S1/S5 | **FIXED** |
| When convenient | §4.1: pass `all_segments` to `update_association` to avoid per-call rebuild | S3 efficiency | **FIXED** |
| When convenient | §4.2: `FaceGeom` cache in `multi_trajectory_fit` to avoid per-vertex double map lookup | S3 efficiency | **FIXED** |
| Deferred | §2.1: port `collect_charge_multi_trajectory` with generalized volume bounds | S1 feature gap | DEFERRED |
| Closed | §4.3: third `form_map_graph` before `dQ_dx_multi_fit` | — | NOT A BUG |
| Closed | §5.1: U/V divide-by-zero in `update_association` | — | NOT APPLICABLE |
| Already addressed | §2.3: `skip_trajectory_point` index semantics | — | NOT A BUG |
| Already fixed | S2 determinism: stable int keys in `form_map_graph` / `dQ_dx_multi_fit` | — | FIXED (single review §3.D1) |

---

## 8. Verification

1. **Single-APA regression**: Run the same single-face test event used for the single-tracking review. Confirm dQ/dx results are unchanged after fixing §2.2 and §2.4 (these should only change results for boundary-crossing tracks).

2. **Multi-APA / boundary track**: Find or construct an event where at least one reconstructed segment has fit points that span two faces (confirmed by `fit.paf` values differing along the segment). Verify that `vertex_fit.paf` and `vertex_fit.pu/pv/pw/pt` after §2.4 fix are consistent with the fitted point's actual face.

3. **Gaussian weighting check (§2.2)**: Call `do_multi_tracking` with `charge_div_method=2`. For a boundary-crossing segment, compare `central_t` values before and after the fix to confirm cross-face pixels now use the correct projection.

4. **Determinism sweep**: Run the same multi-vertex event twice with ASLR disabled (or fixed seed). Confirm segment `fit.index` values and `vertex_fit.pu/pv/pw/pt` are bit-identical across runs. This validates the S2 items.

5. **`update_association` slope check (§5.1, NOT APPLICABLE)**: The U-plane division `raw_y = (coord.wire - offset_u) / slope_yu` is safe for all real detector geometries. U and V planes are defined to have angled wires (non-zero `slope_yu`, `slope_yv`); only W is horizontal (`slope_yw = 0`), and W already divides by `slope_zw` (see code comment at `:2541`). No epsilon guard needed.

6. **Efficiency benchmark (§4.1)**: Profile `form_map_graph` on a high-multiplicity event (many segments meeting at one vertex) and confirm that total time in `update_association` is reduced by the pre-built segment vector change.
