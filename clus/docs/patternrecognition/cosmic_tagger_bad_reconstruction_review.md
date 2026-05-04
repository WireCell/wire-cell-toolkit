# Code Review: cosmic_tagger, bad_reconstruction, calculate_num_daughter_tracks, find_cont_muon_segment_nue

**Toolkit files reviewed:**
- `clus/src/NeutrinoTaggerCosmic.cxx` (lines 97-1321): `bad_reconstruction`, `cosmic_tagger`
- `clus/src/NeutrinoTrackShowerSep.cxx` (lines 221-325): `calculate_num_daughter_tracks`, `find_cont_muon_segment_nue`

**Prototype files compared:**
- `NeutrinoID_cosmic_tagger.h` (lines 1-865): `cosmic_tagger`
- `NeutrinoID_nue_tagger.h` (lines 3450-3852): `bad_reconstruction`
- `NeutrinoID_track_shower.h` (lines 724-763): `calculate_num_daughter_tracks`
- `NeutrinoID_track_shower.h` (lines 2372-2441): `find_cont_muon_segment_nue`

**Porting dictionary reference:** `clus/docs/porting/porting_dictionary.md`

---

## 1. calculate_num_daughter_tracks

### Prototype (NeutrinoID_track_shower.h:724-763)
```
std::pair<int,double> WCPPID::NeutrinoID::calculate_num_daughter_tracks(
    ProtoVertex *vtx, ProtoSegment *sg, bool flag_count_shower, double length_cut)
```

### Toolkit (NeutrinoTrackShowerSep.cxx:221-269)
```
std::pair<int,double> PatternAlgorithms::calculate_num_daughter_tracks(
    Graph& graph, VertexPtr vtx, SegmentPtr sg,
    bool flag_count_shower, double length_cut)
```

### 1.1 Functional Equivalence

**Verdict: Functionally equivalent, with correct improvements.**

| Aspect | Prototype | Toolkit | Match? |
|---|---|---|---|
| BFS traversal | `map_vertex_segments[vtx]` | `boost::out_edges(vd, graph)` | Yes (equivalent) |
| Shower test | `sg->get_flag_shower()` | `flags_any(kShowerTrajectory) \|\| flags_any(kShowerTopology) \|\| abs(pdg)==11` | Yes |
| Counting logic | Count if `!shower \|\| flag_count_shower`, length > cut | Same | Yes |
| Segment insertion | After counting | Before counting (line 238) | **Different - see below** |

**Segment insertion order difference:** In the toolkit, `used_segments.insert(current_sg)` is done atomically via `.insert().second` at line 238, _before_ the counting block. In the prototype, `used_segments.insert(current_sg)` happens at line 750, _after_ counting. This is **functionally identical** because the insertion guard only affects the _next_ time the segment is encountered, and the counting happens on the same iteration regardless. The toolkit's approach is slightly cleaner.

### 1.2 Bugs
None found.

### 1.3 Efficiency
The BFS is efficient. The `std::set` lookups are O(log n). No improvements needed.

### 1.4 Determinism
**Issue:** The toolkit uses `std::set<VertexPtr>` and `std::set<SegmentPtr>` (lines 228-229) for internal tracking (`used_vertices`, `used_segments`). These are ordered by raw pointer value, which is non-deterministic across runs.

**Impact: Low.** These sets are only used for membership testing (`count`/`insert`), never iterated for output. The BFS traversal order is determined by the graph edge order from `boost::out_edges`, not by these sets. So the non-deterministic pointer ordering in these sets does **not** affect the result.

### 1.5 Multi-APA/Face
No APA/face-specific logic. This function operates purely on the graph topology and is APA-agnostic. **OK.**

---

## 2. find_cont_muon_segment_nue

### Prototype (NeutrinoID_track_shower.h:2372-2441)
### Toolkit (NeutrinoTrackShowerSep.cxx:275-325)

### 2.1 Functional Equivalence

**Verdict: Functionally equivalent. Toolkit is a clean simplification.**

| Aspect | Prototype | Toolkit | Match? |
|---|---|---|---|
| Direction at 15cm | `sg->cal_dir_3vector(vtx->get_fit_pt(), 15cm)` | `segment_cal_dir_3vector(sg, vtx_pt, 15cm)` | Yes |
| Direction at 30cm | Same threshold (>30cm for either) | Same | Yes |
| Angle calculation | `(pi - dir1.Angle(dir2))/pi*180` | `(M_PI - dir1.angle(dir2))/M_PI*180` | Yes |
| Angle thresholds | <12.5, or <15 when sg<6cm | Same | Yes |
| dQ/dx threshold | ratio < 1.3 | Same | Yes |
| Selection criterion | `length * cos(angle)` maximized | Same | Yes |

**Simplification in toolkit:** The prototype tracks `max_ratio`, `max_ratio1`, `max_ratio1_length`, and `flag_cont` — but `flag_cont` and `max_ratio1*` are effectively unused (the final check is just `if(flag_cont) return (sg1,vtx1) else return (0,0)`, and `flag_cont` is set iff `sg1` is set). The toolkit correctly removes these dead variables and returns `{sg1, vtx1}` directly. This is a **correct improvement**.

### 2.2 Bugs
None found.

### 2.3 Efficiency
Iterates all edges at a vertex — minimal, no improvement possible.

### 2.4 Determinism
The function selects the segment with maximum `proj = length * cos(angle)`. If two segments have exactly the same projection, the result depends on iteration order from `boost::out_edges`. This is the same situation as the prototype (which iterates `map_vertex_segments[vtx]`). **Acceptable** — exact ties are extremely unlikely with floating-point geometry.

### 2.5 Multi-APA/Face
No APA/face-specific logic. Operates purely on 3D geometry. **OK.**

---

## 3. bad_reconstruction

### Prototype (NeutrinoID_nue_tagger.h:3450-3852)
### Toolkit (NeutrinoTaggerCosmic.cxx:97-449)

### 3.1 Functional Equivalence

**Verdict: Functionally equivalent. All thresholds match.**

#### Sub-check 1 (flag_bad_shower_1)

| Aspect | Prototype | Toolkit | Match? |
|---|---|---|---|
| start_type==1, degree==1, E<120MeV, nsegs<=3 | `pair_result.second==1 && map_vertex_segments[vtx].size()==1` | `start_type==1 && vtx_degree(vtx,graph)==1` | Yes |
| stem >80cm → bad | Yes | Yes | Yes |
| Topology/trajectory check | `!get_flag_shower_topology() && !get_flag_shower_trajectory()` | `!flags_any(kShowerTopology) && !flags_any(kShowerTrajectory)` | Yes |
| Shower energy | `get_kine_best() ?: get_kine_charge()` | `shower->get_kine_best()` | **See note** |

**Energy note:** The prototype explicitly falls back from `get_kine_best()` to `get_kine_charge()` when best is 0. The toolkit calls `shower->get_kine_best()` directly. Per the comment at line 82-84, `PRShower::get_kine_best()` already implements this fallback internally. **Functionally equivalent.**

#### Sub-check 2 (flag_bad_shower_2)

All energy-dependent length thresholds match exactly:

| Energy Range | Prototype thresholds | Toolkit thresholds | Match? |
|---|---|---|---|
| <200 MeV | 38, 42, 46, 50 | 38, 42, 46, 50 | Yes |
| <400 MeV | 42, 49, 52, 55 | 42, 49, 52, 55 | Yes |
| <600 MeV | 45, 48, 54, 62 | 45, 48, 54, 62 | Yes |
| <800 MeV | 51, 52, 56, 62 | 51, 52, 56, 62 | Yes |
| <1500 MeV | 55, 60, 65, 75 | 55, 60, 65, 75 | Yes |
| >=1500 MeV | 55, 65, 70, 75 | 55, 65, 70, 75 | Yes |

Special conditions also match:
- 400MeV: `n_connected + n_connected1 > 4 && max_length <= 72cm` override
- 800MeV: `vtx degree==1 && max_length < 68cm || n_connected >= 6 && max_length < 76cm` override
- 800MeV: `num_segments >= 15 && max_length < 60cm` override
- >1000MeV: `max_length_ratio < 0.95` override
- `max_length > 0.75 * total_length && max_length > 35cm` final override

All match.

#### Sub-check 3 (flag_bad_shower_3)

All thresholds match. The conditional logic for the final guard (lines 421-431 toolkit vs 3819-3828 prototype) also matches exactly.

### 3.2 Bugs
None found.

### 3.3 Efficiency

**Minor improvement opportunity in sub-check 3:** The inner nested loop (lines 353-382) iterates `shower_segs` for each `sg1` in `shower_segs`, giving O(n^2) where n is the number of shower segments. In practice, n is small (typically <20 segments per shower), so this is not a performance concern.

### 3.4 Determinism

The `IndexedSegmentSet shower_segs` at line 119 uses `SegmentIndexCmp` (ordered by `get_graph_index()`), so iteration is deterministic. **Good.**

The `max_length` selection could depend on iteration order in case of ties, but floating-point ties are extremely unlikely. **Acceptable.**

### 3.5 Multi-APA/Face

**Line 210-212:** The cluster-ID comparison `sg1->cluster()->get_cluster_id() != start_cl` is used to add a 6cm length offset for cross-cluster segments. This is the same logic as the prototype's `sg1->get_cluster_id() != shower->get_start_segment()->get_cluster_id()`. In multi-APA, clusters from different APAs will have different IDs, which correctly triggers the offset. **OK.**

---

## 4. cosmic_tagger

### Prototype (NeutrinoID_cosmic_tagger.h:1-865)
### Toolkit (NeutrinoTaggerCosmic.cxx:471-1321)

### 4.1 Functional Equivalence

**Verdict: Functionally equivalent with one noted difference in flag_cosmic_1.**

#### Flag 1: Main vertex outside FV

| Aspect | Prototype | Toolkit |
|---|---|---|
| FV check | `fid->inside_fiducial_volume(test_p, offset_x, &stm_tol_vec)` with stm_tol_vec = -1.5cm all around | `fiducial_utils->inside_fiducial_volume(test_p)` |

**Difference:** The prototype uses a tighter fiducial volume (shrunk by 1.5cm in each direction via `stm_tol_vec`). The toolkit uses the default FV without the STM tolerance shrinkage. There is a TODO comment at line 538 acknowledging this. **This is a known, documented deviation.** The toolkit's FV boundary is slightly _more permissive_ (a point 1cm inside the boundary will pass in the toolkit but fail in the prototype). The impact is that the toolkit will tag slightly fewer events as cosmic via flag_1 — a small signal-efficiency gain at the cost of slightly more cosmic background leaking through this particular cut. The downstream BDT should absorb most of this difference, but it is worth tracking.

**Recommendation:** Extend `FiducialUtils::inside_fiducial_volume` to accept an optional tolerance vector, then pass `-1.5cm` here to match the prototype exactly. Alternatively, if the BDT is retrained on toolkit output, this difference is acceptable.

#### Flags 2-3: Single muon / long muon direction

All logic matches. The muon-finding criteria, `dQ_dx_cut` formula, `valid_tracks` counting, `connected_showers` counting, and the cosmic decision thresholds all match.

**One subtle difference in muon-finding (line 54 prototype):** The prototype has operator precedence issue:
```cpp
if (sg->get_particle_type() == 13 || (!sg->get_flag_shower()) && medium_dQ_dx < dQ_dx_cut * ... && sg->get_particle_type()!=211)
```
Due to C++ precedence, `||` is lower than `&&`, so this is parsed as:
```
(pdg==13) || ((!shower) && dQ_dx_ok && pdg!=211)
```
The toolkit writes this explicitly with proper grouping at line 564:
```cpp
bool muon_like = (pdg == 13) || (!is_shower && med_dqdx < dQ_dx_cut * 1.05 * 43e3/units::cm && pdg != 211);
```
**Functionally identical** due to C++ operator precedence rules, but the toolkit is clearer.

**Energy check in valid_tracks counting (line 872 toolkit):** `Emi < 25 * units::cm` — this compares an energy (in MeV units) against a length (25cm). The prototype (line 398) has the same expression: `Eshower < 25*units::cm`. This appears to be a **bug inherited from the prototype**: the threshold should likely be `25 * units::MeV`, not `25 * units::cm`. With standard WCT units where `units::cm = 1.0` and `units::MeV = 1.0`, this evaluates to the same numerical value (25.0), so **there is no functional difference at runtime**. However, if unit systems were ever changed, this would break. Worth fixing for correctness.

#### Flags 4-5: Upstream-pointing muon

All logic matches.

#### Flags 6-8: Stopped muon + Michel electron

All logic matches. The muon-finding, michel-finding, `muon_2nd` selection, valid_tracks counting, angular adjustments, `flag_sec` determination, and cosmic decision criteria all match.

#### Flag 9: PCA-based cluster direction analysis

| Aspect | Prototype | Toolkit | Match? |
|---|---|---|---|
| Small piece counting | `pts.front().y > 50cm` | `pca.center.y() > 50cm` | **Different** |
| Cluster point collection | `sg->get_point_vec()` (raw 3D points) | `sg->fits()` (fit points) | **Different but better** |
| PCA computation | `main_cluster->Calc_PCA(pts)` then `get_PCA_axis(0)` | `calc_PCA_main_axis(pts).second` | Equivalent |
| Centroid direction | Sum of `(pt - main_vertex_pt)` | Same | Yes |
| Angle thresholds | All match | All match | Yes |
| Shower filter | `shower->get_start_segment()->get_particle_type()!=11` | Same check | Yes |

**Difference: Small-piece Y check.** The prototype checks `pts.front().y` (first point of the first blob in the cluster). The toolkit checks `pca.center.y()` (PCA centroid of the cluster). For small clusters (<3cm), the PCA center is a more robust representative point than the first point. **Toolkit is an improvement.**

**Difference: Fit points vs raw points.** The toolkit uses `sg->fits()` (trajectory-fitted points) rather than raw 3D cloud points. For the purpose of PCA direction analysis, fitted points are smoother and more representative. **Toolkit is an improvement.**

#### Flag 10: Front-end vertex check

All logic matches. Iterates all vertices in the main cluster, checks z<15cm and outside FV, then checks track direction/length.

### 4.2 Bugs

1. **Units mismatch (line 872):** `Emi < 25 * units::cm` should be `Emi < 25 * units::MeV`. As noted above, numerically identical with standard units but semantically wrong. Inherited from prototype.

### 4.3 Efficiency

**Repeated computations in flags 2-3 and 4-5:** The muon/long_muon finding, `valid_tracks` counting, and `connected_showers` counting at lines 551-621 are computed once and reused for flags 2-5. This is already efficient.

**Repeated `find_other_vertex` in flags 2 and 4-5 (lines 625 and 728):** For the muon case, `find_other_vertex(graph, muon, main_vertex)` is called once in the flag 2 block and again in the flag 4 block. These blocks are sequential in the same scope and `muon` doesn't change. The duplicate call could be hoisted, but the cost is negligible (single graph edge lookup).

**Flag 6-8 duplicates muon-finding logic from flags 2-5 (lines 761-808).** This is a faithful port of the prototype, which also duplicates this logic. The second block has a slightly different criterion (requiring `pdg==13` for muon, and looking for Michel electrons). Merging these would complicate the logic for minimal gain.

**Flag 9 point collection (lines 1083-1115):** Iterates all edges and all vertices in the graph to collect per-cluster point clouds. This is O(E+V) which is optimal.

### 4.4 Determinism

**Flag 9: `std::map<int,...>` for cluster maps (lines 1055-1059):** Uses `int` keys (cluster IDs), so iteration order is deterministic. **Good.**

**Flag 10: `boost::vertices(graph)` iteration (line 1270):** The vertex iteration order from BGL depends on the graph's internal storage. Since the graph uses `vecS` for vertex storage, iteration is deterministic (by vertex descriptor index). **Good.**

**Flag 2-8: `segs_at_vtx` lambda:** Uses `boost::out_edges` which iterates edges in insertion order for `vecS` edge storage. Deterministic. **Good.**

### 4.5 Multi-APA/Face

**Flag 1: FV check.** The fiducial volume check via `FiducialUtils` should already handle multi-APA geometry, as it uses `IDetectorVolumes` which defines the FV boundaries per TPC. **Assuming FiducialUtils is multi-APA aware, this is OK.**

**Flag 9: Cluster-level analysis.** The `all_clusters` vector is iterated, which includes clusters from all APAs. The PCA direction analysis is done per cluster. If a neutrino interaction spans multiple APAs (e.g., a muon track crossing an APA boundary split into two clusters), these would be analyzed as separate clusters. This is the correct behavior — the PCA analysis is meant to characterize each cluster independently.

**Flag 10: Front-end vertex check.** The z<15cm check is hard-coded for a single-APA assumption where the detector front face is at z~0. In a multi-APA layout (e.g., DUNE), the front face of each APA module is at a different z location. **This could be a multi-APA issue.** If the detector has APAs at different z positions, the check should use the APA-local z coordinate or the APA front face position rather than a global z<15cm. However, this depends on whether the detector geometry is already transformed to per-APA coordinates. If the FV check (`inside_fv(vpt)`) already excludes points outside the current APA, then the z<15cm check is redundant and harmless. **Recommend verifying that this works correctly for multi-APA.**

---

## 5. calculate_num_daughter_showers

This function (NeutrinoTrackShowerSep.cxx:164-214) is called by `cosmic_tagger` and was also reviewed for completeness.

### 5.1 Functional Equivalence
Matches the prototype `calculate_num_daughter_showers` pattern. The counting logic counts shower segments (or all segments if `flag_count_shower` is false).

### 5.2 Determinism
Same as `calculate_num_daughter_tracks` — internal sets use pointer ordering but only for membership checks. **OK.**

---

## Summary of Findings

### Functional Differences from Prototype

| ID | Location | Description | Status |
|---|---|---|---|
| D1 | `cosmic_tagger` flag 1 (line 539) | Missing STM tolerance shrinkage (-1.5cm) in FV check | **Fixed** — implemented tolerance support in `FiducialUtils::inside_fiducial_volume` and passed `-1.5cm` tolerance vector |
| D2 | `cosmic_tagger` flag 9 (line 1075) | Small-piece Y check uses PCA center instead of first-point Y | **Intentional improvement** |
| D3 | `cosmic_tagger` flag 9 (lines 1089-1095) | Uses fit points instead of raw points for PCA | **Intentional improvement** |

### Bugs

| ID | Location | Description | Status |
|---|---|---|---|
| B1 | `cosmic_tagger` line 872 | `Emi < 25 * units::cm` should be `25 * units::MeV` (units mismatch, inherited from prototype; numerically equivalent with standard units) | **Fixed** |

### Efficiency Improvements

No significant efficiency issues found. All algorithms are linear or near-linear in the number of segments/vertices.

### Determinism

| ID | Location | Description | Severity |
|---|---|---|---|
| R1 | `calculate_num_daughter_tracks` lines 228-229 | Internal `std::set<VertexPtr>` / `std::set<SegmentPtr>` use pointer ordering | **None** — used only for membership checks, not iterated for output |
| R2 | `calculate_num_daughter_showers` lines 168-169 | Same as R1 | **None** |

**All externally-visible sets/maps use indexed comparators (`SegmentIndexCmp`, `VertexIndexCmp`, `ShowerIndexCmp`), providing deterministic iteration.** The toolkit has correctly addressed the prototype's use of pointer-ordered `std::set`s in the function interfaces.

### Multi-APA/Face Concerns

| ID | Location | Description | Status |
|---|---|---|---|
| M1 | `cosmic_tagger` flag 10 | Hard-coded `vpt.z() < 15cm` assumes detector front face at z~0 | **Fixed** — now computes `z_front` from `IDetectorVolumes::inner_bounds()` and checks `vpt.z() < z_front + 15cm` |
| M2 | `cosmic_tagger` flag 9 (line 1075) | `pca.center.y() > 50cm` and `highest_y > 80/100/102cm` thresholds assume single-APA vertical geometry | **Low** — these are conservative cuts that should still work for multi-APA (cosmic rejection remains valid for vertical tracks) |

---

## Fixes Applied

1. **B1 (NeutrinoTaggerCosmic.cxx):** Changed `Emi < 25 * units::cm` to `Emi < 25 * units::MeV` — corrects a units mismatch inherited from the prototype.

2. **D1 (FiducialUtils.cxx + NeutrinoTaggerCosmic.cxx):** Implemented tolerance support in `FiducialUtils::inside_fiducial_volume()`. The method now accepts tolerance vectors of size 1 (uniform), 3 (per-axis), or 6 (per-face). For negative tolerances (shrinking), it checks that the point offset outward by |tolerance| in each axis direction is still contained. The `cosmic_tagger` flag 1 now passes `stm_tol_vec = {-1.5cm}` (6 elements, all -1.5cm) to match the prototype's STM boundary shrinkage.

3. **M1 (NeutrinoTaggerCosmic.cxx):** Replaced hard-coded `vpt.z() < 15cm` with `vpt.z() < z_front + 15cm`, where `z_front` is computed from the minimum z of all sensitive volumes via `IDetectorVolumes::inner_bounds()`. Falls back to `z_front=0` if no detector volumes are available, preserving prototype behavior for single-APA.
