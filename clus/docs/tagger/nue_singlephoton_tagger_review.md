# Code Review: `singlephoton_tagger`

**Date:** 2026-04-09  
**Reviewer:** Claude Code  
**Scope:** Logic fidelity vs prototype, bug hunting, efficiency, determinism.

---

## Files Examined

| Role | File |
|---|---|
| Toolkit implementation | `clus/src/NeutrinoTaggerSinglePhoton.cxx` (2485 lines) |
| Prototype implementation | `prototype_pid/src/NeutrinoID_singlephoton_tagger.h` (4275 lines) |

---

## Functions Reviewed

| Function | Toolkit Location | Prototype Location |
|---|---|---|
| `singlephoton_tagger` (entry point) | line 2166 | line 2 |
| `bad_reconstruction_sp` | line 135 | line 3766 |
| `bad_reconstruction_1_sp` | line 481 | line 4170 |
| `bad_reconstruction_2_sp` | line 611 | line 3461 |
| `bad_reconstruction_3_sp` | line 904 | line 3223 |
| `mip_identification_sp` | line 1113 | line 2102 |
| `high_energy_overlapping_sp` | line 1578 | line 2596 |
| `low_energy_overlapping_sp` | line 1765 | line 2797 |
| `pi0_identification_sp` | line 1956 | line 2955 |
| `low_energy_michel_sp` | line 2093 | line 545 |

Supporting file-local helpers reviewed: `vtx_fit_pt`, `vtx_degree`, `shower_energy`, `seg_is_shower`, `seg_endpoint_near` (lines 61–94); `SpContext` struct (lines 102–116).

---

## Summary of Findings

Two logic bugs found in the toolkit entry point. Both fixed in this review.

- **Bug 1** — First-pass shower loop fills the real `TaggerInfo` (with vector fields) for every electron shower instead of a throw-away object.  Prototype passed `flag_fill=false` to all bad-reconstruction helpers in the first pass.  Result: vector fields (`shw_sp_br3_3_v_*`, `shw_sp_br3_5_v_*`, `shw_sp_br3_6_v_*`) over-populated; scalar br1/br2/br3/br4 fields overwritten by last-processed shower rather than max\_shower.
- **Bug 2** — First-pass call to `bad_reconstruction_1_sp` hardcoded `num_valid_tracks=0` instead of computing the per-shower value.  Prototype computed `num_valid_tracks` inside the loop, skipping each shower's own start segment.  Result: more aggressive first-pass br2 rejection when valid tracks exist at the main vertex, potentially selecting the wrong `max_shower`.

No bugs found in any helper function (bad_reconstruction_*, mip_identification_sp, high_energy_overlapping_sp, low_energy_overlapping_sp, pi0_identification_sp, low_energy_michel_sp).

---

## Structure Overview

The prototype is a single large function (~542 lines for the entry point) that calls ten helper member functions scattered throughout the header.  The toolkit restructures this as:

1. **File-local helpers** — `bad_reconstruction_sp`, `bad_reconstruction_1_sp`, `bad_reconstruction_2_sp`, `bad_reconstruction_3_sp`, `mip_identification_sp`, `high_energy_overlapping_sp`, `low_energy_overlapping_sp`, `pi0_identification_sp`, `low_energy_michel_sp`.
2. **Public entry point** — `PatternAlgorithms::singlephoton_tagger()` (line 2166).
3. **SpContext** — bundle of shared state, replacing prototype member variables.

The entry point follows the same two-pass structure as the prototype:
- **First pass**: classify all electron showers at main\_vertex (badreco1–4), accumulate counters, select `max_shower`.
- **Second pass**: run the full suite of checks on `max_shower` only, fill TaggerInfo, update `flag_sp`.

---

## Per-Function Findings

### ✅ `bad_reconstruction_sp` (br1_1, br1_2, br1_3)

- **br1_1** (stem/topology check, prototype lines 3789–3808): Energy, type, vertex degree, segment count, topology/trajectory flags, and length comparisons all match exactly. ✅
- **br1_2** (muon-track extension via `find_cont_muon_segment_nue`, prototype lines 3813–3974):  Length thresholds, the 6 cm offset for topology/off-cluster segments, `n_connected`/`n_connected1` accumulation, and the full energy-binned cut table all match exactly.  The `check_len2` lambda replaces the prototype's verbose if-else chain — same logic. ✅
- **br1_3** (long straight track at far end of stem, prototype lines 3976–4161): The three-way dis1/dis2/else case structure, the inner sg2 loop (angle1 > 170, 0.75·min(dis1,dis2) threshold), the energy-binned cut table (`check_len3`), and the "only flag if max\_length > sg length" guard all match. ✅

### ✅ `bad_reconstruction_1_sp` (br2)

- PCA axis computation, `seg_endpoint_near` for vertex-side endpoint, `dir_shower` / `dir1` / `dir2` / `dir3` angle computation, the `flag_single_shower || num_valid_tracks == 0` guard, energy-binned cuts (500/1000 MeV), and the two additional unconditional cuts (angle>40; trajectory+length cut) all match. ✅
- `max_angle` loop over far-end vertex edges: prototype iterates `map_vertex_segments[other_vertex]`; toolkit uses `boost::out_edges(other_vertex, graph)` — accepted graph-API translation. ✅

### ✅ `bad_reconstruction_2_sp` (br3_1 … br3_8)

- **br3_1** (four straight-shower conditions): The third condition in the prototype (line 3497) has a subtle operator-precedence issue:
  ```cpp
  // prototype:
  if (... && (direct_length/length > 0.95 || fabs(...)<5 && sg->get_flag_shower_trajectory()) && length > 25*units::cm)
  // reads as:  >0.95  ||  (<5 && trajectory)  because && binds tighter than ||
  ```
  The toolkit (line 645–650) explicitly parenthesizes this as `(>0.95 || (<5 && trajectory))` — correct and unambiguous. ✅
- **br3_2** (segment composition + fiducial): n\_ele/n\_other counting, the two bad2 conditions, and `other_fid` via `FiducialUtilsPtr` all match. ✅
- **br3_3 / br3_4** (backward segments): 🐛 **Bug fixed 2026-04-12.** The `angle > 105` check and all TaggerInfo vector pushes (`shw_sp_br3_3_v_*`) were incorrectly placed inside the `if (dir1.magnitude() > 10*units::cm)` guard. In the prototype (line 3590), `angle > 105 && len > 15*cm` is NOT guarded by `dir1.Mag() > 10*cm` — only the `angle > 90` accumulation and `angle > 150` check are. Also, the TaggerInfo vectors are pushed for EVERY segment (prototype lines 3592-3597), not only those with `dir1.Mag() > 10*cm`. Fix: moved `angle > 105` check and vector pushes outside the guard, computing `angle` unconditionally (defaults to 0 when `dir1.magnitude() <= 10*cm`). The rest of the loop (main cluster filter, `acc_length`/`total_main_len2` accumulation) matches.
- **br3_5** (average non-stem point): `side_total_length` replaces prototype's `total_length` (reset at line 3623); `dir1 = ave_p - other_point`; the `flag_bad5 = false` override for multi-cluster showers uses `shower->get_total_length(sg->cluster())`. ✅
- **br3_6** (segments at far-end vertex): Prototype iterates `map_vertex_segments[other_vertex]` (graph-wide); toolkit iterates `boost::out_edges(other_vertex, graph)` — correct translation.  `n_other_vtx_segs = boost::out_degree(...)` matches `map_vertex_segments[other_vertex].size()`. ✅
- **br3_7** (`shower_main_len = vertex->cluster() ? shower->get_total_length(vertex->cluster()) : 0`; prototype uses `vertex->get_cluster_id()`): equivalent. ✅
- **br3_8** (sliding dQ/dx window): `sg1->fits().size()-5` inner loop matches `int(sg1->get_point_vec().size())-5`; `segment_median_dQ_dx(sg1, i, i+5)` is the per-range overload matching `sg1->get_medium_dQ_dx(i, i+5)`. ✅

### ✅ `bad_reconstruction_3_sp` (br4_1, br4_2)

- **br4_1** (farthest main-cluster vertex, closest off-cluster segment): `segment_get_closest_point` replaces `sg->get_closest_point`; six threshold conditions, two false-override conditions, and the `num_main_segs>=4` guard all match exactly. ✅
- **br4_2** (angular distribution of shower fit points): The prototype's nested loop counts interior fit points then iterates `it->second` for vertex positions.  The toolkit uses `find_vertices(ctx.graph, sg1)` to get the two endpoints and adds them separately.  For a well-formed graph each segment has exactly two endpoint vertices, so the count matches.  The full `flag_bad2` expression (two compound conditions) matches exactly. ✅
- `n_vtx_segs` for the `size()==1` guard in br4_1: prototype uses `map_vertex_segments[shower->get_start_vertex().first].size()`; toolkit uses `boost::out_degree(start_vtx, graph)` — correct. ✅

### ✅ `mip_identification_sp`

This is the largest helper (~450 lines).  All of the following match the prototype exactly:

- `dQ_dx_cut` energy-binned table (1.45 / 1.6 / 1.85 / 1.3). ✅
- `n_end_reduction`, `n_first_mip`, `n_first_non_mip`, `n_first_non_mip_1`, `n_first_non_mip_2` scan loops. ✅
- Primary MIP classification condition (9 sub-clauses). ✅
- Refinements for `mip_id==-1` cases (7013, 6640 events). ✅
- `flag_strong_check` branch. ✅
- `n_good_tracks` loop (boost edges replacing map\_vertex\_segments). ✅
- Energy-dependent `vtx_deg==1 / >1` corrections.
  - Prototype uses `map_vertex_segments[vertex].size()`; toolkit uses `vtx_degree(vertex, ctx.graph)` = `boost::out_degree(vertex, graph)`. ✅
- Angular cuts on low-energy showers (ang\_beam > 40 / >30 blocks). ✅
- `min_dQ_dx_5` single-shower condition; shape cuts (`lowest_dQ_dx`, `highest_dQ_dx`, `iso_angle`). ✅
- Other-shower loop (`ctx.showers`; prototype iterates `showers`). ✅
- `min_dis` computation: prototype iterates `map_vtx_segs` for main-cluster vertices; toolkit builds `main_cl_vtxs` from `shower->fill_sets` and does the same pairwise distance. ✅
- `n_other_vertex`: `vtx_degree(other_vertex, graph)` replacing `map_vertex_segments[other_vertex].size()`. ✅
- Median/mean dQ/dx → dedx conversion (Bethe-Bloch, first 7 values, same formula). ✅
- All 20 `shw_sp_vec_dQ_dx_*` fills. ✅

### ✅ `high_energy_overlapping_sp` (hol_1, hol_2) | **Minor fix 2026-04-12**

- **hol_1**: `n_valid_tracks`, `min_angle`, `flag_all_showers`, `num_showers` computation; the three `flag_overlap1` conditions; dQ/dx boost conditions. 🐛 **Minor bug fixed 2026-04-12:** When `dir2.magnitude() == 0` for an electron/weak-muon segment, the toolkit incorrectly set `flag_all_showers = false` before `continue`. The prototype (line 2641) simply calls `continue` without modifying `flag_all_showers`. Fix: removed the `flag_all_showers = false` assignment.
- **hol_2**: `min_ang2` / `min_sg` scan; `ncount` loop (consecutive fit points within 0.6 cm); `medium_dQ_dx` from `min_sg` using near-endpoint index — prototype checks `wcpt` front/back equality; toolkit uses geometric distance. Both yield the same near endpoint under the accepted convention. ✅

### ✅ `low_energy_overlapping_sp` (lol_1, lol_2, lol_3)

- **lol_1**: Prototype iterates `map_vtx_segs` (shower's internal vertex map), selecting vertices with `it->second.size()==2` — the toolkit collects shower segments at each shower vertex via `boost::out_edges`, checks `vtx_ss.size()==2`.  The cut fires only on the 2-segment case; `front()`/`back()` are deterministic. ✅
- **lol_2**: Prototype checks `(*it)->get_particle_type() == 13`; toolkit checks `sg1->particle_info()->pdg() == 13` — equivalent. ✅
- **lol_3**: `n_out/n_sum` counts use shower fit points (interior points) and shower vertices — matches prototype's `map_seg_vtxs` + `map_vtx_segs` loops. ✅

### ✅ `pi0_identification_sp` (pio_1, pio_2)

- **pio_1**: mass window conditions, asymmetric-pair veto (7049), `min(Eshower_1,Eshower_2) > max(10 MeV, threshold)` — all match. ✅
- **pio_2**: Prototype iterates `map_vertex_segments` (all TPC vertices); toolkit iterates `graph_nodes(ctx.graph)`.  Both see all graph vertices.  `acc_length` per-cluster precomputed from `boost::edges(ctx.graph)` vs prototype's `map_segment_vertices` scan — same set. ✅
- Return value: `pi0_flag_pio` (is shower in pio map, not flag\_pi0) — correctly returned in both. ✅

### ✅ `low_energy_michel_sp`

- `n_3seg` counting: prototype uses `it1->second.size()>=3` from `map_vtx_segs` (shower internal vertex-to-segment map); toolkit counts segments in `sh_segs` that are connected to each shower vertex via boost edges — equivalent. ✅
- Two `flag_bad` conditions match exactly. ✅
- Prototype also computes `E_range` (and recalculates if zero) but never uses it in the decision logic; toolkit omits it cleanly. ✅

---

## 🐛 Bug 1 — First-pass TaggerInfo fill pollutes vector fields

**Location:** `singlephoton_tagger` entry point, first-pass shower loop (toolkit lines 2305–2318 before fix).

**Prototype behaviour (lines 177–213):** All four `bad_reconstruction_*_sp` helpers are called with `flag_fill=false` in the per-shower classification loop.  TaggerInfo is only filled (with `flag_fill=true`) in the second pass for `max_shower` alone.

**Toolkit bug:** The toolkit helpers fill TaggerInfo unconditionally (no `flag_fill` parameter).  In the first pass they were called with the real `ti` reference, causing:
1. **Scalar fields** (`shw_sp_br1_*`, `shw_sp_br2_*`, `shw_sp_br3_1/2/4/7/8_*`, `shw_sp_br4_*`) overwritten once per electron shower; the final value reflects whichever shower was processed last — not `max_shower`.
2. **Vector fields** (`shw_sp_br3_3_v_*`, `shw_sp_br3_5_v_*`, `shw_sp_br3_6_v_*`) populated once per shower in the first pass, then again for `max_shower` in the second pass → vectors contain entries from all showers plus a duplicate of max\_shower.

**Impact:** Medium–high.  Scalar fields get corrected by the second-pass overwrite; vector fields (which are BDT inputs) are wrong when more than one electron shower is present at the main vertex.

**Fix (applied):** Introduce a local `TaggerInfo tmp_ti{}` inside the per-shower loop and direct the first-pass calls to it:

```cpp
// Before fix (line 2306):
bool badreco1 = !bad_reconstruction_sp(ctx, shower, ti);
// ...
bool badreco3 = !bad_reconstruction_2_sp(ctx, shw_vtx_main, shower, ti);
bool badreco4 = !bad_reconstruction_3_sp(ctx, shw_vtx_main, shower, ti);

// After fix:
TaggerInfo tmp_ti{};
bool badreco1 = !bad_reconstruction_sp(ctx, shower, tmp_ti);
// ...
bool badreco3 = !bad_reconstruction_2_sp(ctx, shw_vtx_main, shower, tmp_ti);
bool badreco4 = !bad_reconstruction_3_sp(ctx, shw_vtx_main, shower, tmp_ti);
```

---

## 🐛 Bug 2 — First-pass `bad_reconstruction_1_sp` uses hardcoded `num_valid_tracks=0`

**Location:** `singlephoton_tagger` entry point, first-pass shower loop (toolkit line 2310 before fix).

**Prototype behaviour (lines 183–194):** For each electron shower in the classification loop, the prototype computes `num_valid_tracks` by iterating `map_vertex_segments[main_vertex]` and counting non-shower segments that are long enough (> 8 cm, or > 5 cm if not dir-weak), excluding the current shower's start segment.  This value is then passed to `bad_reconstruction_1_sp`.

**Toolkit bug:** The first-pass call hardcoded `num_valid_tracks=0`.  Inside `bad_reconstruction_1_sp`, the cut block at line 556:
```cpp
if (flag_single_shower || num_valid_tracks == 0) { ... }
```
With `num_valid_tracks=0`, this block always fires.  With the correct value, it is skipped when the vertex has valid tracks (and the event is not single-shower).  As a result, for multi-shower, multi-track events, the first-pass `badreco2` is more conservative than the prototype: more showers are incorrectly classified as having bad stem direction.  This can change which shower is selected as `max_shower`.

**Impact:** Medium.  Affects multi-track events with > 1 electron shower at the main vertex, which are likely to be background events anyway — but the wrong shower could still be selected.

**Fix (applied):** Compute `first_pass_valid_tracks` inside the loop body before calling `bad_reconstruction_1_sp`:

```cpp
// Compute num_valid_tracks for this shower (mirrors prototype lines 183-190)
int first_pass_valid_tracks = 0;
{
    auto vd = main_vertex->get_descriptor();
    for (auto [eit, eend] = boost::out_edges(vd, graph); eit != eend; ++eit) {
        SegmentPtr sg1 = graph[*eit].segment;
        if (!sg1 || sg1 == sg) continue;
        double len1 = segment_track_length(sg1);
        if (!seg_is_shower(sg1) &&
            (len1 > 8*units::cm || (!sg1->dir_weak() && len1 > 5*units::cm)))
            ++first_pass_valid_tracks;
    }
}
// ...
bool badreco2 = sg_at_main
    ? !bad_reconstruction_1_sp(ctx, shower, flag_single_shower,
                                first_pass_valid_tracks, tmp_ti)
    : true;
```

---

## Summary Table

| Function | Status | Notes |
|---|---|---|
| `bad_reconstruction_sp` (br1) | ✅ Correct | All three sub-checks match |
| `bad_reconstruction_1_sp` (br2) | ✅ Correct | PCA logic, angle cuts match |
| `bad_reconstruction_2_sp` (br3) | 🐛 **br3_3 scoping bug fixed** | br3_3 `angle>105` and vector fills were inside `dir1.mag>10cm` guard |
| `bad_reconstruction_3_sp` (br4) | ✅ Correct | br4_1 + br4_2 match |
| `mip_identification_sp` | ✅ Correct | All scan loops, classification, cuts match |
| `high_energy_overlapping_sp` (hol) | 🐛 **Minor fix** | hol_1: removed incorrect `flag_all_showers=false` on zero-mag dir |
| `low_energy_overlapping_sp` (lol) | ✅ Correct | lol_1/2/3 match |
| `pi0_identification_sp` (pio) | ✅ Correct | pio_1 + pio_2 match |
| `low_energy_michel_sp` | ✅ Correct | Matches; dead `E_range` code omitted cleanly |
| `singlephoton_tagger` (entry) | 🐛🐛 **2 bugs fixed** | See above |

---

## Changes Made

Both bugs are in `clus/src/NeutrinoTaggerSinglePhoton.cxx`, first-pass shower classification loop (approximately lines 2305–2318 before the fix, lines 2305–2336 after):

1. Replaced `TaggerInfo& ti` with a local `TaggerInfo tmp_ti{}` as the destination for all four `bad_reconstruction_*_sp` calls in the first-pass loop.
2. Added a `first_pass_valid_tracks` computation (14 lines) matching prototype lines 183–190, and passed it instead of the hardcoded `0` to `bad_reconstruction_1_sp`.

No changes were made to any helper function in the initial review.

**2026-04-12:** Two additional helper bugs found and fixed:
3. `bad_reconstruction_2_sp` br3_3: Moved `angle > 105` check and TaggerInfo vector pushes outside the `dir1.magnitude() > 10*units::cm` guard to match prototype scope.
4. `high_energy_overlapping_sp` hol_1: Removed incorrect `flag_all_showers = false` when `dir2.magnitude() == 0` (prototype does not set this flag in this case).

---

## Prototype Bugs Fixed by Toolkit

1. **PROTO-BUG-1 (Low) — `num_muons`/`num_pions` parenthesis bug**: Prototype lines 127-128 have `abs(sg->get_particle_type()==13)` which evaluates the comparison first (bool), then abs. Anti-particles (pdg=-13, -211) count in `num_mip_tracks` but not in `num_muons`/`num_pions`. Toolkit correctly uses `std::abs(pdg) == 13` / `== 211`.

---

## Efficiency / Structure Improvements

1. **SpContext bundle**: Shared state packaged into a file-local struct, replacing scattered member-variable access.
2. **Static helpers**: All 9 helpers are file-local `static` functions.
3. **flag_fill elimination**: Throw-away `TaggerInfo tmp_ti{}` replaces the prototype's `flag_fill=false` parameter threading.
4. **pi0 acc_length precomputation**: O(E) edge scan + O(1) lookup replaces prototype's O(V*S) per-vertex recomputation.
5. **Dead code removal**: `E_range`, unused `n_segs`, `flag_overlap_4/5`, `max_dQ_dx` in br1_2 — all removed.
6. **Debug output removal**: All `std::cout`/`flag_print` output removed.
7. **Overall**: 4275 lines → 2511 lines (41% reduction).

---

## Determinism

| Container | Prototype | Toolkit | Issue? |
|-----------|-----------|---------|--------|
| Shower loop at main vertex | `map_vertex_to_shower[v]` (pointer set) | `VertexShowerSetMap` = `IndexedShowerSet` | **Fixed** — index-ordered |
| Segment/vertex sets in helpers | `map_seg_vtxs`/`map_vtx_segs` (pointer-keyed) | `IndexedSegmentSet`, `IndexedVertexSet` | **Fixed** — index-ordered |
| `cluster_acc_length` (pi0) | `map_segment_vertices` (pointer-keyed) | `std::map<Cluster*, double>` | **OK** — lookup only, never iterated |
| pio_2 vector population | `map_vertex_segments` (pointer-ordered) | `graph_nodes()` (pointer-ordered) | **Inherited** — non-deterministic vector order in both; BDT treats as unordered |

---

## Multi-APA / Multi-Face

The single-photon tagger operates on graph topology and segment-level quantities. Drift direction is hardcoded as `Vector(1,0,0)` (same as prototype), but all drift-angle computations use `fabs(angle - 90°)` (symmetric w.r.t. drift sign). APA/face correctly derived from main vertex position via `dv->contained_by()`. No wire-plane geometry references. **No multi-APA issues found.**

---

## Angle Safety

All angle computations use `D3Vector::angle()` which internally clamps the cosine to [-1, 1] before `std::acos`. No NaN risk.
