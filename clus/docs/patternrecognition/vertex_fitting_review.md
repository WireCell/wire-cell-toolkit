# Vertex Fitting & Overall Main Vertex — Porting Review

Review of 11 functions ported from WCP prototype to WCT toolkit, covering vertex
fitting (MyFCN, fit_vertex, improve_vertex), overall main vertex determination
(determine_overall_main_vertex, determine_overall_main_vertex_DL), cluster
management (examine_main_vertices, check_switch_main_cluster, swap_main_cluster,
compare_main_vertices_global, calc_dir_cluster), and 2D charge collection
(collect_2D_charge).

## 1. Function Map

| # | Prototype Name | Toolkit Name | Prototype File:Line | Toolkit File:Line |
|---|---|---|---|---|
| 1 | `improve_vertex` | `PatternAlgorithms::improve_vertex` | `NeutrinoID_improve_vertex.h:44` | `NeutrinoVertexFinder.cxx:1911` |
| 2 | `MyFCN` (class) | `MyFCN` (class) | `NeutrinoID_improve_vertex.h:504` | `MyFCN.cxx:14` / `MyFCN.h:17` |
| 3 | `fit_vertex` | `PatternAlgorithms::fit_vertex` | `NeutrinoID_improve_vertex.h:11` | `NeutrinoVertexFinder.cxx:1844` |
| 4 | `determine_overall_main_vertex_DL` | same | `NeutrinoID_DL.h:4` | `NeutrinoVertexFinder.cxx:3166` |
| 5 | `determine_overall_main_vertex` | same | `NeutrinoID.cxx:333` | `NeutrinoVertexFinder.cxx:3380` |
| 6 | `collect_2D_charges` | `TrackFitting::collect_2D_charge` | `NeutrinoID_energy_reco.h:1` | `TrackFitting.cxx:897` |
| 7 | `examine_main_vertices()` | same | `NeutrinoID.cxx:409` | `NeutrinoPatternBase.cxx:1922` |
| 8 | `check_switch_main_cluster()` | same + `_2` variant | `NeutrinoID.cxx:736,953` | `NeutrinoVertexFinder.cxx:3033,3127` |
| 9 | `swap_main_cluster` | same | `NeutrinoID.cxx:728` | `NeutrinoPatternBase.cxx:1899` |
| 10 | `calc_dir_cluster` | same | `NeutrinoID.cxx:684` | `NeutrinoPatternBase.cxx:1839` |
| 11 | `compare_main_vertices_global` | same | `NeutrinoID.cxx:794` | `NeutrinoVertexFinder.cxx:2871` |

## 2. Bugs Found and Fixed

### Bug 1 (CRITICAL) — `MyFCN::AddSegment` eigenvalue ordering

**File:** `MyFCN.cxx:128-143`

Eigen's `SelfAdjointEigenSolver` returns eigenvalues in ascending order (smallest
first), while ROOT's `TMatrixDEigen` used in the prototype returns them in
descending order (largest first).  The rest of the toolkit codebase already handles
this — `Facade_Cluster.cxx:2307` explicitly reverses with `[2-i]`,
`NeutrinoPatternBase.cxx:1826` uses column 2 for the largest eigenvalue, and
`ClusteringFuncs.cxx:171` does the same.  Only `MyFCN.cxx` was missing the
reversal.

`FitVertex` zeros out row 0 of the R matrix (no constraint along that direction)
and weights the remaining rows by `sqrt(lambda_0/lambda_k)`.  With the prototype's
descending order, row 0 corresponds to the track direction (largest eigenvalue) and
the perpendicular directions get amplified weights >1.  With ascending order,
row 0 is the smallest-eigenvalue direction, the wrong direction is zeroed out, and
perpendicular weights become <1 — a ~30x difference in constraint strength for
typical track geometry.

**Fix:** Reversed the eigenvalue/eigenvector assignment order so that index 0 holds
the largest eigenvalue, matching the prototype convention.

### Bug 2 (minor) — `improve_vertex` non-deterministic `fitted_vertices`

**File:** `NeutrinoVertexFinder.cxx:1914`

`std::set<VertexPtr> fitted_vertices` is iterated at line 2204 (else branch) to
apply direction determination.  Since `VertexPtr = shared_ptr<Vertex>`, the set
orders by pointer address — non-deterministic across runs.

**Fix:** Changed to `IndexedVertexSet` (ordered by stable `get_graph_index()`).

### Bug 3 (cosmetic) — redundant `ParticleInfo` creation

**File:** `NeutrinoVertexFinder.cxx:3477-3481`

In `determine_overall_main_vertex`, the proton-tagging block created a
`ParticleInfo`, set pdg=2212 and mass, then immediately replaced it with a new
`ParticleInfo` containing the 4-momentum.  The first creation was dead code.
Compare with `determine_overall_main_vertex_DL:3332-3336` which does it correctly
in one step.

**Fix:** Removed the redundant creation; only the final `ParticleInfo` with
4-momentum is now constructed.

## 3. Determinism Improvements (toolkit over prototype)

| Function | Improvement |
|---|---|
| `compare_main_vertices_global` | Sorts vertex_candidates by `cluster_id` before scoring; final max-selection iterates sorted vector |
| `check_switch_main_cluster` | Sorts candidates by `cluster_id` before comparison |
| `determine_overall_main_vertex` | Sorts `all_candidate_clusters` by `cluster_id`, deduplicates before length scan |
| `examine_main_vertices` | Sorts phase-2 candidates by `cluster_id` for deterministic evaluation order |
| `improve_vertex` | Uses `ordered_nodes(graph)` for vertex iteration instead of `map_vertex_segments` |
| `determine_overall_main_vertex_DL` | Uses `ordered_edges(graph)` for segment point collection |
| `improve_vertex` (fixed) | `fitted_vertices` changed from `std::set<VertexPtr>` to `IndexedVertexSet` |

## 4. Multi-APA / Multi-Face Correctness

| Function | Multi-APA Handling |
|---|---|
| `MyFCN::UpdateInfo` | Uses `dv->contained_by(fit_pos)` for APA/face, looks up per-APA geometry via `WirePlaneId(kAllLayers, face, apa)`, uses `PC transform` for coordinate conversion |
| `fit_vertex` | Validates charge at new vertex position via APA/face-aware `grouping->average_3d_charge()` |
| `collect_2D_charge` | Fundamentally redesigned: iterates per-APA `m_charge_data`, handles wrapped wires via `get_wires_for_channel(apa, channel)`, builds `(apa, channel) -> (face, plane, wire)` geometry map |
| Others (`calc_dir_cluster`, `compare_main_vertices_global`, etc.) | Operate on 3D coordinates; no explicit APA/face needed |

## 5. Per-Function Equivalence Notes

### `MyFCN` (class)

- **Constructor, destructor, `update_fit_range`, `get_fittable_tracks`, `get_seg_info`:** Identical logic.
- **`AddSegment`:** Same point filtering (vertex_protect_dis / fit_dis), same PCA computation, same regularization (+0.15cm^2). Toolkit uses Eigen's `SelfAdjointEigenSolver` instead of ROOT's `TMatrixDEigen`; eigenvalue ordering now fixed (see Bug 1). Uses `sg->wcpts()` (steiner points) consistent with prototype's `get_point_vec()`.
- **`FitVertex`:** Mathematically identical — same angle check (15 deg), same BiCGSTAB solver, same vertex constraint logic (`1/vtx_constraint_range * sqrt(npoints)`). Toolkit uses `Eigen::Vector3d` dot product instead of `TVector3::Angle()`.
- **`UpdateInfo`:** Enhanced for multi-APA. Prototype uses singleton `TPCParams` for wire geometry; toolkit uses `IDetectorVolumes` + `TrackFitting` geometry maps. Prototype determines segment front/back by wcpt index matching; toolkit uses distance comparison (more robust when indices change). Steiner path interpolation logic identical (2cm step, 0.3cm rejection threshold).

### `fit_vertex`

Identical logic: creates `MyFCN` with `flag_vtx_constraint=true`, `vtx_constraint_range=0.43cm`, adds segments, calls `FitVertex`. Toolkit adds charge validation: queries `grouping->average_3d_charge()` at old and new positions, rejects if new_charge < 5000 && < 0.4*old or new_charge < 8000 && < 0.6*old. Uses `dv->contained_by()` for APA/face resolution.

### `improve_vertex`

Four-phase structure matches prototype:
1. Setup: counts tracks/showers, optionally calls `examine_structure_4()`
2. Vertex fitting: iterates vertices via `ordered_nodes(graph)`, applies `fit_vertex`, refits if moved >0.5cm
3. Activity search: `search_for_vertex_activities()`, refit 3-segment vertices, `eliminate_short_vertex_activities()`
4. Main vertex segment handling: direction determination for shower/track classification

Key difference: prototype iterates `map_vertex_segments` (pointer-ordered); toolkit iterates `ordered_nodes(graph)` (deterministic). The shower classification, topology checks, and proton tagging thresholds are all identical.

### `determine_overall_main_vertex`

Identical logic flow: find max-length cluster, call `examine_main_vertices()`, check cluster switch (dev_chain → `check_switch_main_cluster`, frozen → `check_switch_main_cluster_2`), proton-tag short high-dQ/dx stubs, clean up long muon sets.

Toolkit sorts `all_candidate_clusters` by `cluster_id` and deduplicates for deterministic max-length search. Proton tagging uses same thresholds (daughter_count==1, length<1.5cm, median_dQdx>1.6).

### `determine_overall_main_vertex_DL`

Same algorithm: build (x,y,z,q) point cloud from vertices + segment interior points, run DL inference, find nearest vertex, validate (direction sanity check, distance cut), optionally switch main cluster and re-examine directions.

Toolkit improvements: `try/catch` around DL inference, parameterized weights path and DL vertex cut (prototype hardcoded), `ordered_edges(graph)` for deterministic segment iteration. Prototype finds cluster to swap by iterating `other_clusters` with cluster_id comparison; toolkit uses `min_vertex->cluster()` pointer comparison.

### `collect_2D_charge`

Fundamentally redesigned for multi-APA support. Prototype: collects time/channel ranges from cluster point clouds, calls `ct_point_cloud->get_overlap_good_ch_charge()` per plane, then `fill_2d_charge_dead_chs()`. Toolkit: iterates `m_charge_data` (pre-populated unordered_map of `CoordReadout -> ChargeMeasurement`), classifies by plane using `get_wires_for_channel()`, builds geometry map `(apa, channel) -> vector<(face, plane, wire)>`. Two-pass design: first pass splits charges + collects unique channels, second pass builds geometry map (avoids redundant lookups).

### `examine_main_vertices`

Two-pass structure matches prototype:
1. Remove short clusters without strong tracks (same thresholds: `cluster_length_cut = min(main_length*0.6, 6cm)`, track-like = not_shower && dirsign!=0 && (!dir_weak || dQdx>2)); for survived-but-short (<5cm) clusters far from main (>100cm), remove anyway
2. If main vertex is all-showers: check for cluster swap (angle>160, length>10cm, length>0.5*main) or removal (angle<10, small+close)

Toolkit uses KD-tree `kd_steiner_knn(1, ...)` for closest-distance queries (O(log N)) vs prototype's `pcloud->get_closest_dis()` (O(N)). Phase-2 candidates sorted by `cluster_id`.

### `check_switch_main_cluster` / `check_switch_main_cluster_2`

Two overloads map to the prototype's two overloads. First (no-arg equivalent): checks if main vertex has all showers, collects all vertices from map, calls `compare_main_vertices_global`, swaps if better vertex found. Second (2-arg equivalent): simple shower check + swap to max_length_cluster.

Toolkit sorts candidates by `cluster_id` before comparison. Logic identical.

### `swap_main_cluster`

Identical logic: unset main_cluster flag on old, push old to other_clusters, set flag on new, erase new from other_clusters. Toolkit takes references instead of pointers and manages `Facade::Flags::main_cluster` flags.

### `calc_dir_cluster`

Identical algorithm: accumulate fitted points and vertex positions within `dis_cut` of `orig_p`, compute mean direction, normalize. Toolkit uses squared-distance comparison (avoids sqrt per point) — minor efficiency gain. Skips first/last fitted points same as prototype.

### `compare_main_vertices_global`

Same multi-component scoring:
- Z-position penalty: `-(z - min_z) / 200cm`
- Segment scoring: shower +0.125, track +0.25, clear proton +0.25, track with direction +0.125
- Main cluster bonus: +0.25
- Fiducial volume: +0.5
- Direction pointing: angle<15 → +0.25, angle<30 → +0.125
- Isolation penalty: delta==0 && not_main && total_length<6cm → -0.25*num_tracks

Toolkit sorts candidates by `cluster_id` before scoring; iterates sorted vector for max selection (deterministic). Prototype uses `map_vertex_num` pointer-keyed map iterated by `vertex_candidates` pointer vector.

## 6. Efficiency Notes

| Function | Improvement |
|---|---|
| `calc_dir_cluster` | Squared-distance comparison avoids `sqrt()` per point |
| `examine_main_vertices` | KD-tree nearest-neighbor (O(log N)) vs linear scan (O(N)) |
| `determine_overall_main_vertex` | Deduplicates candidate clusters to avoid redundant `get_length()` calls |
| `collect_2D_charge` | Two-pass design avoids redundant `get_wires_for_channel()` calls |

## 7. Summary

| Function | Status | Notes |
|---|---|---|
| `MyFCN::AddSegment` | **Bug fixed** | Eigenvalue ordering reversed to match prototype |
| `MyFCN::FitVertex` | Equivalent | Same BiCGSTAB solver and constraint logic |
| `MyFCN::UpdateInfo` | Better | Multi-APA support, distance-based front/back |
| `fit_vertex` | Better | Charge validation at new position |
| `improve_vertex` | **Bug fixed** + Better | `IndexedVertexSet` for determinism; `ordered_nodes` for iteration |
| `determine_overall_main_vertex` | **Cleaned up** + Better | Removed dead code; deterministic cluster sorting |
| `determine_overall_main_vertex_DL` | Better | try/catch, parameterized weights, ordered_edges |
| `collect_2D_charge` | Redesigned | Multi-APA native; wrapped wire handling |
| `examine_main_vertices` | Better | KD-tree, sorted candidates |
| `check_switch_main_cluster` | Better | Sorted candidates |
| `swap_main_cluster` | Equivalent | Flag management added |
| `calc_dir_cluster` | Better | Squared-distance optimization |
| `compare_main_vertices_global` | Better | Sorted candidates, deterministic max selection |
