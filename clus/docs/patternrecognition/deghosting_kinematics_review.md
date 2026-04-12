# Review: `deghosting` / `deghost_clusters` / `deghost_segments` / `calculate_shower_kinematics` / `calculate_kinematics_long_muon` / `cal_kine_charge` / `calculate_kinematics`

**Reviewed by:** Claude Code  
**Date:** 2026-04-12  
**Branch:** `apply-pointcloud`  
**Prototype refs:**
- `prototype_base/pid/src/NeutrinoID_deghost.h` (486 lines)
- `prototype_base/pid/src/NeutrinoID_energy_reco.h` (455 lines)
- `prototype_base/pid/src/NeutrinoID_shower_clustering.h` (lines 1407–1427)
- `prototype_base/pid/src/WCShower.cxx` (lines 288–527)

---

## Implementation status

| Item | Description | Status | File |
|---|---|---|---|
| B.1 | `deghosting()`: `map_cluster_main_vertices` passed by value — cleanup never propagates back | ✅ Fixed | `NeutrinoDeghoster.cxx:589`, `NeutrinoPatternBase.h:214` |
| B.2 | `deghost_segments` vertex cleanup: `boost::vertices(graph)` (pointer-order) instead of `ordered_nodes` | ✅ Fixed | `NeutrinoDeghoster.cxx:575` |
| L.1 | `cal_corr_factor()` returns 1.0 — U-plane gain correction not yet implemented | Documented only | `NeutrinoEnergyReco.cxx:14` |

---

## 1. Scope

This review audits:

| Prototype | Toolkit | Toolkit file:line |
|---|---|---|
| `NeutrinoID::deghosting()` | `PatternAlgorithms::deghosting(graph, map_cluster_main_vertices, all_clusters, track_fitter, dv)` | `NeutrinoDeghoster.cxx:589` |
| `NeutrinoID::deghost_clusters()` | `PatternAlgorithms::deghost_clusters(graph, all_clusters, track_fitter, dv)` | `NeutrinoDeghoster.cxx:74` |
| `NeutrinoID::deghost_segments()` | `PatternAlgorithms::deghost_segments(graph, map_cluster_main_vertices, all_clusters, track_fitter, dv)` | `NeutrinoDeghoster.cxx:353` |
| `NeutrinoID::calculate_shower_kinematics()` | `PatternAlgorithms::calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon, graph, track_fitter, dv, particle_data, recomb_model)` | `NeutrinoEnergyReco.cxx:269` |
| `WCShower::calculate_kinematics_long_muon(segments_in_muons)` | `Shower::calculate_kinematics_long_muon(segments_in_muons, particle_data, recomb_model)` | `PRShower.cxx:1126` |
| `WCShower::calculate_kinematics()` | `Shower::calculate_kinematics(particle_data, recomb_model)` | `PRShower.cxx:870` |
| `NeutrinoID::cal_kine_charge(WCShower*)` | `PatternAlgorithms::cal_kine_charge(ShowerPtr, ...)` | `NeutrinoEnergyReco.cxx:182` |
| `NeutrinoID::cal_kine_charge(ProtoSegment*)` | `PatternAlgorithms::cal_kine_charge(SegmentPtr, ...)` | `NeutrinoEnergyReco.cxx:235` |

Review goals: (1) functional equivalence, (2) bug detection, (3) efficiency, (4) determinism, (5) multi-APA/face correctness.

Helper functions also reviewed: `order_clusters`, `order_segments`, `kine_charge_from_maps`, `cal_corr_factor`, `fill_kine_tree`.

---

## 2. Per-function walk-through

### `deghosting`

**Prototype** (`NeutrinoID_deghost.h:13–43`): calls `deghost_clusters()`, then `deghost_segments()`, then cleans up `map_cluster_main_vertices` by erasing clusters whose main vertex is no longer in `map_vertex_segments` (i.e., was deleted during deghosting). Since `map_cluster_main_vertices` is a member variable, the erasure persists.

**Toolkit** (`NeutrinoDeghoster.cxx:589`): same sequence, but the cleanup tries to erase from `map_cluster_main_vertices`. **Bug B.1**: the parameter was originally declared as `ClusterVertexMap map_cluster_main_vertices` (by value), so erasure only affected the local copy. Fixed by changing to `ClusterVertexMap& map_cluster_main_vertices`. The cleanup logic itself (checking `boost::out_degree(vit, graph) > 0`) is correct.

The timing instrumentation (`DG_Clock`, `DG_ms`) present in the prototype is omitted in the toolkit; replaced by structured debug logging. This is a toolkit-wide convention.

---

### `order_clusters`

**Prototype** (`NeutrinoID_deghost.h:455–485`): iterates `map_segment_vertices` (a `std::map<ProtoSegment*, ...>`) — pointer-address order, non-deterministic. Cluster length is a running total.

**Toolkit** (`NeutrinoDeghoster.cxx:27–72`): iterates `ordered_edges(graph)` — deterministic. **Improvement.** The sort comparator `sortbysec` adds a tiebreaker on `cluster->ident()` (line 17), making the output fully deterministic when two clusters have equal total length. **Improvement.**

---

### `deghost_clusters`

**Purpose:** Detect clusters that are 2D-overlap duplicates of already-accepted clusters and remove them.

**Algorithm (identical in prototype and toolkit):**
1. Sort clusters by total segment length, longest first.
2. Seed three global 2D point clouds from all clusters that have no segments (noise-only).
3. For the longest cluster (i=0): unconditionally add its point cloud, steiner cloud, and skeleton cloud.
4. For each subsequent cluster: count per-plane dead-channel and unique (non-overlapping) points against the three global clouds. Apply three ghosting criteria (detailed below). If ghosted, collect for removal. If not ghosted, add to global clouds.

**Ghosting criteria (lines 265–276, matching prototype exactly):**
- Condition 1: high dead fraction (≥80%) AND low unique fraction
- Condition 2: very low unique fraction across all planes
- Condition 3: medium dead fraction (70–80%) AND low unique fraction
- Extra check: one plane entirely dead, two planes have zero unique points, AND `max_unique_percent < 0.75`

**Distance thresholds (matching prototype):**
- Global point cloud: `dis_cut/2` (0.6 cm)
- Global steiner cloud: `dis_cut*2/3` (0.8 cm)
- Global skeleton cloud: `dis_cut*6/4` (1.8 cm)

**Key differences from prototype (all improvements):**

1. **Division-by-zero protection** (line 250): toolkit guards `if (num_total_points > 0)` before computing percentages. Prototype would compute `NaN` for empty clusters (undefined behavior).

2. **Skeleton cloud population** (lines 130–139): toolkit uses `seg->wcpts()` for the skeleton cloud. Prototype uses `sg->get_point_vec()` (all blob 3D points). The choice of wcpts avoids t0-corrected x coordinates that could land outside the detector volume; for 2D wire-overlap checking the wire projections of wcpts are equivalent.

3. **Multi-TPC 2D distances**: toolkit passes `(plane, face, apa)` to `global_point_cloud->get_closest_2d_point_info(test_point, plane, face, apa)`, correctly selecting the wire geometry for the test point's APA/face. Prototype used a single global wire geometry.

4. **Points outside detector**: toolkit checks `apa == -1 || face == -1` (line 161) and skips such points, also not counting them in `num_total_points`. This is correct — points between APAs have ambiguous wire projections. Prototype counted all points.

5. **Orphan vertex cleanup** (line 320): uses `ordered_nodes(graph)` — deterministic.

---

### `deghost_segments`

**Purpose:** Within the surviving clusters (after `deghost_clusters`), identify and remove individual segments that are 2D-overlap duplicates of the global skeleton built from longer segments.

**Algorithm (identical in prototype and toolkit):**
1. Seed `global_point_cloud` from clusters with no segments.
2. Early return if `global_point_cloud` is empty (no reference to compare against).
3. For each cluster (longest first), for each segment (longest first):
   - Only analyze if segment is terminal (`start_n==1 || end_n==1`), low dQ/dx (<1.1 MIP), and long (>3.6 cm).
   - Count per-plane dead-channel and unique (non-overlapping) points against the three global clouds.
   - If all planes have zero unique points: mark for removal.
   - Protection: if segment is the only connection to the cluster's main vertex, keep it regardless.
4. After processing all segments in a cluster: add the cluster's full point cloud and steiner cloud to the globals.

**Distance thresholds (matching prototype):**
- Global point cloud: `dis_cut*2/3` (0.8 cm)
- Global steiner cloud: `dis_cut*2/3` (0.8 cm)
- Global skeleton cloud: `dis_cut*3/4` (0.9 cm)

**Key differences from prototype (all improvements):**

1. **Deterministic segment ordering** (line 346): `std::stable_sort` preserves the relative order of equal-length segments (which come from `ordered_edges` — deterministic). Prototype's `std::sort` is not stable.

2. **Multi-TPC 2D distances**: same as `deghost_clusters` — passes `(plane, face, apa)` to each cloud query.

3. **Main vertex protection** (lines 525–538): toolkit checks `(v1 == main_vtx && out_degree(source)==1) || (v2 == main_vtx && out_degree(target)==1)`, equivalent to prototype's `map_vertex_segments[tmp_vtx].find(sg)!=end && size()==1`.

4. **Orphan vertex cleanup** (lines 575–581): **Bug B.2** (fixed): was using raw `boost::vertices(graph)` (pointer-order iteration), now changed to `ordered_nodes(graph)` for consistency with `deghost_clusters`.

5. **Unused variable suppression** (line 512): `(void) num_total_points;` — `num_total_points` is accumulated but not used in the removal decision for `deghost_segments` (only the `num_unique` sums matter). This is correct behavior, matching prototype.

---

### `calculate_shower_kinematics`

**Prototype** (`NeutrinoID_shower_clustering.h:1407–1427`): iterates `showers` vector; for each not-yet-done shower: dispatches to `shower->calculate_kinematics()` or `calculate_kinematics_long_muon(segments_in_long_muon)` based on particle type, then calls `cal_kine_charge(shower)` and sets the result.

**Toolkit** (`NeutrinoEnergyReco.cxx:269–336`): functionally identical control flow, with improvements:

1. **Dependency injection**: passes `particle_data` and `recomb_model` to kinematics functions instead of singleton access. **Improvement.**

2. **Pre-collected charge maps**: uses `m_charge_2d_u/v/w` (populated once by `shower_clustering_with_nv` via `collect_charge_maps()`), falling back to on-demand collection if empty. This avoids re-reading 2D charge maps for every shower. **Efficiency improvement.**

3. **`vertices_in_long_muon` parameter**: declared but unused (`(void)vertices_in_long_muon` at line 270). The prototype does not use this parameter at this stage either — it's used inside `calculate_kinematics_long_muon` only via `segments_in_long_muon`. The extra parameter exists for future extension.

4. **`get_flag_shower()` return path**: toolkit calls `shower->get_flag_kinematics()` and skips already-computed showers. Prototype does the same (`if (!shower->get_flag_kinematics())`). Equivalent.

---

### `Shower::calculate_kinematics`

**Prototype** (`WCShower.cxx:339–527`): three cases — single segment, multi-segment single track (all segments connected), multi-segment multiple tracks (disconnected pieces).

**Toolkit** (`PRShower.cxx:870–1124`): identical three-case structure. Differences:

1. **Deterministic end_point**: toolkit uses `ordered_nodes(*this, m_full_graph)` (line 954, 1061) — deterministic index-based ordering. Prototype iterates `map_vtx_segs` (pointer-keyed map, non-deterministic on ties). **Improvement.**

2. **Start-point fallback**: multi-segment case (line 1026): prototype tests `start_point == (0,0,0)` as a sentinel; toolkit uses `sgcp_dist >= 0` returned by `shower_get_closest_point` as a validity signal. The `(0,0,0)` test would falsely trigger fallback for a valid hit at the origin; the `sgcp_dist` test is correct. **Bug fix** (labeled B16.1 in code comment).

3. **Zero init_dir fallback**: toolkit explicitly checks `if (data.init_dir.magnitude() == 0)` (lines 977, 1047) and falls back to `shower_cal_dir_3vector`. Prototype does not have this guard; it would silently produce a zero direction. **Improvement.**

4. **`flag_shower` definition**: toolkit mirrors prototype's three-part definition:
   `kShowerTrajectory || kShowerTopology || (|pdg|==11)`.

---

### `Shower::calculate_kinematics_long_muon`

**Prototype** (`WCShower.cxx:288–336`): iterates `map_seg_vtxs` and `map_vtx_segs` (pointer-keyed, non-deterministic). Finds length from muon segments, collects all dQ/dx, finds farthest muon-connected vertex.

**Toolkit** (`PRShower.cxx:1126–1190`): equivalent logic with improvements:

1. **Deterministic edge iteration**: uses `ordered_edges(*this, m_full_graph)` (line 1149). **Improvement.**

2. **Deterministic vertex deduplication**: uses `std::map<size_t, VertexPtr> muon_vertices_by_index` (line 1147) keyed by graph index — stable across runs. Prototype used pointer-keyed containers. **Improvement.**

3. **Single pass**: prototype used two separate loops (one to accumulate length, one to find farthest vertex). Toolkit combines accumulation and vertex collection into one loop. **Efficiency improvement.**

---

### `cal_kine_charge` (shower and segment overloads)

**Prototype** (`NeutrinoID_energy_reco.h:44–455`): large function for each overload, iterating `charge_2d_u/v/w` maps. For each 2D charge hit: converts time/channel to 2D coordinates, finds closest 3D point in the shower/segment's associated and fitted point clouds (within 0.6 cm), applies `cal_corr_factor`, accumulates per-plane charge sums.

**Toolkit** (`NeutrinoEnergyReco.cxx:182–266`): both overloads delegate to `kine_charge_from_maps()` helper (anonymous namespace, lines 47–178). The core algorithm is preserved; key differences:

1. **2D coordinate lookup**: toolkit uses `grouping->convert_time_wire_2Dpoint(time_slice, local_wire, apa, face, plane_id)` (line 84), passing the plane-local wire index. Prototype passed the global channel number directly. **Bug fix**: V channels start at ~2400 and W at ~4800; using the global channel as a wire index would produce enormous wrong distances. The toolkit correctly uses the local wire index extracted from `map_apa_ch_plane_wires`.

2. **Direct 2D query**: uses `get_closest_2d_point_info_direct(p2d.first, p2d.second, plane_id, face, apa)` (line 94) — the `_direct` variant avoids applying the angle projection again (since `convert_time_wire_2Dpoint` already provides coordinates in wire-perpendicular space). **Correctness fix.**

3. **Charge weighting / asymmetry logic**: identical to prototype — finds min/med/max plane by charge, computes max asymmetry, uses weighted average unless one plane is significantly higher (max_asy > 0.04). Numerics match.

4. **Fudge/recombination factors**: identical constants (shower: 0.8/0.5; proton: 0.95/0.35; default: 0.95/0.7).

5. **`cal_corr_factor` is a stub** (limitation L.1): returns 1.0. Prototype applies two corrections: (a) hard-coded U-plane gain correction for specific wire ranges (296–327, 336–337, etc.) using `factor /= 0.7`, and (b) a configurable correction factor from `TPCParams::get_corr_factor()`. These detector-specific corrections are not yet ported. The impact is a systematic ~43% energy over-estimation in U-plane regions where prototype applies the factor, and a further gain correction mismatch if the configurable correction was non-trivial.

---

### `fill_kine_tree` (helper)

**Toolkit** (`NeutrinoKinematics.cxx:43`): fills the `KineInfo` output structure from the graph traversal. This has no direct 1:1 prototype function — it is a consolidation of inline code from the prototype's main processing loop. The BFS traversal and energy accumulation logic match prototype behavior, with the following toolkit improvements:

1. **Deterministic traversal**: uses `ordered_out_edges` via `boost::out_edges` (prototype used pointer-keyed `map_vertex_segments`). Since graph edge insertion order is deterministic (set up by `add_segment`), this provides deterministic results.

2. **Structured output**: wraps all output into `KineInfo` struct instead of scattered member variables.

3. **SCE correction**: uses `geom_helper->get_corrected_point(nu_vtx, IClusGeomHelper::SCE, apa, face)` for the neutrino vertex position. Prototype used a different correction mechanism.

---

## 3. Bugs found and fixed

### B.1 — `deghosting`: `map_cluster_main_vertices` by-value prevents cleanup from propagating

**Files:** `NeutrinoDeghoster.cxx:589`, `NeutrinoPatternBase.h:214`

**Problem:** The prototype's `deghosting()` erases stale cluster entries from `map_cluster_main_vertices` after deghost runs. In the toolkit, the parameter was `ClusterVertexMap map_cluster_main_vertices` (by value). The cleanup loop at lines 598–620 erased from the local copy only; the caller's map (`TaggerCheckNeutrino.cxx:233`) retained entries for clusters whose main vertices had been removed. Downstream calls to `determine_overall_main_vertex` at line 250–253 use `map_cluster_main_vertices` and could encounter stale vertex pointers.

**Fix:** Changed parameter to `ClusterVertexMap& map_cluster_main_vertices` in both the header and implementation. The cleanup logic itself was already correct.

### B.2 — `deghost_segments`: vertex cleanup uses `boost::vertices` (pointer-order)

**File:** `NeutrinoDeghoster.cxx:575`

**Problem:** The orphan vertex cleanup in `deghost_segments` used `auto [vbegin, vend] = boost::vertices(graph)` — vertex descriptor iteration order in BGL depends on insertion order, which can depend on pointer addresses. `deghost_clusters` already uses `ordered_nodes(graph)` for the same cleanup.

**Fix:** Changed to `for (auto v : ordered_nodes(graph))` — consistent with `deghost_clusters` and deterministic.

---

## 4. Known limitations (not fixed)

### L.1 — `cal_corr_factor` is a stub

**File:** `NeutrinoEnergyReco.cxx:14`

`cal_corr_factor()` returns 1.0. The prototype applies:
1. A per-wire-range U-plane gain correction (`factor /= 0.7` for U wires 296–327, 336–337, 343–351, 376–400, 410–484, 501–524, 536–671)
2. A configurable overall gain correction from `TPCParams`

The wire-range list is MicroBooNE-specific. The toolkit's multi-detector design requires this to be implemented via the detector-geometry API rather than hard-coded wire ranges. Until implemented, `cal_kine_charge` results will be systematically higher than prototype in the affected wire regions.

---

## 5. Efficiency notes

1. **`calculate_shower_kinematics`**: pre-collected 2D charge maps (`m_charge_2d_u/v/w`) are computed once per event via `collect_charge_maps()`. Each `cal_kine_charge` call reuses these maps instead of re-querying the 2D data store. **O(n_showers) improvement** over prototype which re-read the maps per shower.

2. **`calculate_kinematics_long_muon`**: single-pass over edges accumulates both length and dQ/dx. Prototype used two separate loops.

3. **`kine_charge_from_maps`**: shared implementation eliminates code duplication between shower and segment overloads of `cal_kine_charge` and `calculate_shower_kinematics`.

---

## 6. Determinism improvements (toolkit vs. prototype)

| Location | Prototype (non-deterministic) | Toolkit (deterministic) |
|---|---|---|
| `order_clusters` | iterates `map_segment_vertices` (pointer-keyed map) | `ordered_edges(graph)` |
| `order_clusters` sort tiebreaker | none (unstable for equal lengths) | `cluster->ident()` tiebreaker |
| `order_segments` sort | `std::sort` (not stable) | `std::stable_sort` |
| `deghost_clusters` vertex cleanup | N/A (member-variable-based) | `ordered_nodes(graph)` |
| `deghost_segments` vertex cleanup | N/A | `ordered_nodes(graph)` (after B.2 fix) |
| `calculate_kinematics` end_point | iterates `map_vtx_segs` (pointer order) | `ordered_nodes(*this, m_full_graph)` |
| `calculate_kinematics_long_muon` edge iteration | `map_seg_vtxs` (pointer order) | `ordered_edges(*this, m_full_graph)` |
| `calculate_kinematics_long_muon` vertex dedup | pointer-keyed set | `map<size_t, VertexPtr>` by graph index |

---

## 7. Multi-APA / multi-face correctness

All three deghosting functions use per-APA/face 2D geometry:

- **Point cloud queries**: `global_*_cloud->get_closest_2d_point_info(test_point, plane, face, apa)` — the face/apa arguments select the correct wire geometry (wire spacing, angle, offset) for the test point's detector region.

- **Dead channel queries**: `cluster->grouping()->get_closest_dead_chs(test_point, 1, apa, face, plane)` — selects the correct dead-wire map.

- **Out-of-detector guard**: `if (apa == -1 || face == -1) continue` — points between APAs are excluded from both the numerator and denominator of the unique/dead-fraction calculation. This prevents spurious "ghosting" for clusters that span the inter-APA region.

- **Skeleton cloud population**: `dv->contained_by(wcpt.point)` per wcpt, used as the WirePlaneId for `make_points_direct`, so each skeleton point is registered with its correct APA/face.

The kinematics functions (`cal_kine_charge`, `calculate_shower_kinematics`) use `grouping->convert_time_wire_2Dpoint(time_slice, local_wire, apa, face, plane_id)` and `get_closest_2d_point_info_direct` with explicit APA/face. This is correct for multi-TPC.

---

## 8. Summary of improvements over prototype

1. Division-by-zero protection in `deghost_clusters` when `num_total_points == 0`
2. `order_clusters`: deterministic via `ordered_edges` + ident tiebreaker
3. `order_segments`: deterministic via `stable_sort`
4. `deghost_clusters` / `deghost_segments`: full multi-TPC 2D overlap checking with APA/face routing
5. `deghost_segments`: `ordered_nodes` for vertex cleanup (B.2 fix)
6. `deghosting`: by-reference `map_cluster_main_vertices` so cleanup propagates (B.1 fix)
7. `calculate_kinematics`: `ordered_nodes` for deterministic end_point; `sgcp_dist >= 0` sentinel (not zero-point test) for start_point fallback; zero init_dir guard
8. `calculate_kinematics_long_muon`: single-pass accumulation; `ordered_edges` + index-keyed vertex map
9. `cal_kine_charge`: local-wire vs. global-channel for 2D coordinate lookup; `_direct` query avoids double projection
10. `calculate_shower_kinematics`: pre-collected charge maps avoid per-shower re-read; dependency injection for particle_data and recomb_model
