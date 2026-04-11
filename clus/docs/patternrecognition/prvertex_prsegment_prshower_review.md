# PRVertex / PRSegment / PRShower Port Review

**Reviewed:** 2026-04-11  
**S1 fixes applied (2026-04-11):** F2 (ClusterIdCmp in clustering_points_segments), F3 (SegmentIndexCmp for seg_dpc/pts caches), F5 (paf populated in clear_fit), F6 (n2 clamp in segment_track_direct_length), F7 (empty-wcpts guard in find_vertices), F8 (HasCluster<Vertex> CRTP fix), F9 (0,0,0 sentinel fallback in two missing calculate_kinematics branches), F10 (IndexedSegmentSet for complete_structure and get_last_segment_vertex_long_muon), F11 (null guards in set_start_vertex and calculate_kinematics_long_muon). F1 left as-is: fill_maps() returning *this is correct for all existing callers (update_shower_maps handles global map population).  
**S2 fixes applied (2026-04-11):** F15 (fuse twin calculate_kinematics branches into one shared block, init_dir zero-fallback now also covers single-track case), F16 (lazy cache for get_total_length() and get_num_main_segments() — invalidated in add_segment/add_shower/set_start_segment), F17 (global per-(plane,apa,face) NFKDVec::Tree<double,IndexDynamic>(2) built after seg_pts2d pre-build; inner O(S) segment loop in ghost-removal replaced with three O(log N) knn(1) queries — PRSegmentFunctions.cxx:1768–1783, 1930–1952), F18 (pre-count reserve in fill_point_vector), F19 (set_fit_associate_vec takes std::vector by value + std::move at call site in TrackFitting.cxx).  
**S3 fixes applied (2026-04-11):** F4 (get_angles_cached replaced with eager pre-population from seg_pts2d keys — eliminates silent (0,0,0) caching when the first segment in cache order lacked the queried (apa,face), and removes iteration-order dependence. Lookup is now O(1) via ang_cache.emplace). Items 4.1/4.2 (F2/F3/F10) were already applied in S1; items 4.3/4.4 are documentation-only.  
**Multi-APA / Minor / Doc fixes applied (2026-04-11):** F14 (DynamicPointCloud::merge_wpid_params() + get_wpid_params() added; set_start_segment, add_segment, add_shower now call merge_wpid_params before add_points — fixes crash when shower DPC is queried for a face not present in the first segment; add_shower fixed to use TrajectoryView::add_segment to eliminate double-add of DPC points). F20 (create_segment_fit_point_cloud: removed origin-exclusion guard — all fit points included unconditionally). F21 (find_vertices: removed redundant boost::edge re-lookup after valid stored descriptor). F22 (porting_dictionary.md: filled in WCT equivalents for ProtoSegment, ProtoVertex, WCShower). F12 documented as non-bug for ±x-drift detectors (|angle-90°| is drift-sign-symmetric); F13 documented as intentionally single-face (query point is a specific-face 2D measurement).  
**Toolkit entry files:**
- `clus/inc/WireCellClus/PRVertex.h` (107 lines)
- `clus/inc/WireCellClus/PRSegment.h` + `clus/src/PRSegment.cxx` (176 + 120)
- `clus/inc/WireCellClus/PRSegmentFunctions.h` + `clus/src/PRSegmentFunctions.cxx` (123 + 2544)
- `clus/inc/WireCellClus/PRShower.h` + `clus/src/PRShower.cxx` (247 + 1266)
- `clus/inc/WireCellClus/PRShowerFunctions.h` + `clus/src/PRShowerFunctions.cxx` (26 + 178)
- `clus/inc/WireCellClus/PRGraph.h` + `clus/src/PRGraph.cxx` (141 + 142)
- `clus/inc/WireCellClus/PRCommon.h` (261) — `WCPoint`, `Fit`, `Flagged`, `Graphed`, `HasCluster`, `HasDPCs`
- `clus/inc/WireCellClus/PRTrajectoryView.h` + `clus/src/PRTrajectoryView.cxx`

**Prototype counterparts:**
- `prototype_base/pid/inc/WCPPID/ProtoVertex.h` + `pid/src/ProtoVertex.cxx`
- `prototype_base/pid/inc/WCPPID/ProtoSegment.h` + `pid/src/ProtoSegment.cxx`
- `prototype_base/pid/inc/WCPPID/WCShower.h` + `pid/src/WCShower.cxx`

**Scope:** Functional equivalence, bugs, algorithmic efficiency, determinism (pointer-ordered containers), and multi-APA / multi-face correctness.

**Key reference notes:**
- Wire ranges are half-open `[min, max)` in the toolkit; the `<` vs prototype's `<=` is intentional — see `do_single_tracking_review.md`.
- The prototype `WCShower` is the analogue of toolkit `PR::Shower`. There is no type called `ProtoShower` in the prototype.
- Toolkit model is graph-centric: a `PR::Graph` owns `Vertex` nodes and `Segment` edges; `Shower` is a `TrajectoryView` over a subgraph. Prototype used pointer maps (`map_vtx_segs`, `map_seg_vtxs`).
- `PRSegmentFunctions.cxx` is the toolkit equivalent of `ProtoSegment.cxx`'s algorithmic body (moved to free functions taking `SegmentPtr`).
- The prototype assumed a single APA/face. Multi-face support is the primary new dimension in the toolkit port; that is examined separately in §5.

---

## §1 Function Map (prototype → toolkit)

### 1.1 ProtoVertex → PR::Vertex

| Prototype method | Toolkit equivalent | Toolkit file:line | Notes |
|---|---|---|---|
| `get_wcpt() / set_wcpt()` | `wcpt() / wcpt(v)` | `PRVertex.h:55-65` | Chainable setter added — improvement |
| `get_fit_pt() / set_fit_pt()` | `fit().point` | `PRVertex.h:59-66` | Consolidated into `Fit` struct |
| `set_fit(pt,dQ,dx,pu,pv,pw,pt,chi2)` | `fit()` + caller sets fields | `PRCommon.h:119-148` | No aggregate setter; callers set `Fit` fields directly. Not a bug but a style change. |
| `get/set_fit_flag()` | `SegmentFlags::kFit` (on Segment) | — | Vertex-level fit flag not ported; exists only on Segment. Acceptable. |
| `get/set_fit_index()` | `fit_index() / fit_index(int)` | `PRVertex.h:69-70` | Equivalent |
| `get/set_flag_fit_fix()` | `flag_fix() / flag_fix(bool)` | `PRVertex.h:73-74` | Renamed, equivalent |
| `get/set_fit_range()` | `fit_range() / fit_range(double)` | `PRVertex.h:71-72` | Equivalent |
| `reset_fit_prop()` | `reset_fit_prop()` | `PRVertex.h:77-81` | Equivalent; clears `Fit.index/range/flag_fix` |
| `get/set_dQ/dx/pu/pv/pw/pt/reduced_chi2` | `fit().dQ` etc. | `PRCommon.h:121` | Consolidated into `Fit` struct — improvement |
| `get_dQ_dx()` | `fit().dQ / (fit().dx + 1e-9)` inline | — | No dedicated helper; callers compute inline |
| `get_id()` | `get_graph_index()` | `PRVertex.h:92` | Semantically different: `id` was an auto-incremented int; `graph_index` is a stable index within the graph. Equivalent in uniqueness. |
| `get_cluster_id()` | `cluster()->get_cluster_id()` via `HasCluster<>` | `PRCommon.h:30` | See §2.1 bug — `HasCluster<Segment>` CRTP typo |
| `get_dis(Point)` | `get_dis(Point)` | `PRVertex.h:88-90` | Equivalent |
| `get_fit_init_dis()` | `fit_distance()` | `PRVertex.h:83-86` | Equivalent |
| `set/get_neutrino_vertex()` | `VertexFlags::kNeutrinoVertex` via `Flagged<>` | `PRVertex.h:23` | Equivalent; cleaner flag model |
| `flag_fit_fix` (separate bool member) | merged into `Fit::flag_fix` | `PRCommon.h:126` | Improvement — single source of truth |

### 1.2 ProtoSegment → PR::Segment + PRSegmentFunctions

Prototype allocated 8 parallel `std::vector` members for fit data; toolkit uses `std::vector<Fit>` (struct-of-value). All algorithmic methods moved to free functions in `PRSegmentFunctions.cxx`.

| Prototype method | Toolkit equivalent | Toolkit file:line | Notes |
|---|---|---|---|
| `get_wcpt_vec()` | `wcpts()` | `PRSegment.h:92` | Equivalent |
| `get_point_vec() / set_point_vec()` | `fits()` (the `point` field) | `PRSegment.h:85-88` | Consolidated |
| `get_dQ_vec() / get_dx_vec()` etc. | `fits()` (individual fields) | `PRCommon.h:121` | Consolidated into `Fit` struct |
| `build_pcloud_fit()` | `create_segment_fit_point_cloud(seg, dv, name)` | `PRSegmentFunctions.cxx:25` | Multi-face aware — improvement |
| `set_fit_associate_vec(pts,idx,skip)` | `set_fit_associate_vec(fits, dv, name)` | `PRSegment.cxx:37` | Takes pre-built `Fit` vector; DI for `dv` |
| `reset_fit_prop()` | `reset_fit_prop()` | `PRSegment.cxx:6` | Equivalent |
| `clear_fit()` | `clear_fit(dv, cloud_name)` | `PRSegment.cxx:80` | Needs `dv` for DPC rebuild; see §2.3 bug |
| `get_fit_index_vec()` | `fit_index(i)` / `fit_index(i,idx)` | `PRSegment.cxx:12-22` | Equivalent (indexed access) |
| `get_fit_flag_skip()` | `fit_flag_skip(i)` / `fit_flag_skip(i,flag)` | `PRSegment.cxx:24-34` | Equivalent |
| `get_length()` / `get_length(n1,n2)` | `segment_track_length(seg,n1,n2)` | `PRSegmentFunctions.cxx:~650` | Equivalent |
| `get_direct_length(n1,n2)` | `segment_track_direct_length(seg,n1,n2)` | `PRSegmentFunctions.cxx:~720` | See §2.4 inconsistent-clamp bug |
| `get_max_deviation(n1,n2)` | `segment_track_max_deviation(seg,n1,n2)` | `PRSegmentFunctions.cxx:759` | Equivalent |
| `get_medium_dQ_dx()` | `segment_median_dQ_dx(seg)` | `PRSegmentFunctions.cxx:~823` | Equivalent |
| `get_rms_dQ_dx()` | `segment_rms_dQ_dx(seg)` | `PRSegmentFunctions.cxx:~867` | Equivalent |
| `get_closest_wcpt(Point)` | **MISSING** | — | No toolkit equivalent; callers scan `seg->wcpts()` manually |
| `get_closest_point(Point)` | `segment_get_closest_point(seg, pt, name)` | `PRSegmentFunctions.cxx:~100` | Equivalent |
| `get_closest_2d_dis(Point)` | `segment_get_closest_2d_distances(seg,pt,apa,face,name)` | `PRSegmentFunctions.cxx:137` | Signature extended for multi-face; see §5 |
| `get_closest_2d_dis(x,y,plane)` | `segment_get_closest_2d_distance(seg,pt,apa,face,plane,name)` | `PRSegmentFunctions.cxx:172` | Equivalent |
| `break_segment_at_point(p)` | `break_segment(seg, pt, graph, name)` | `PRSegmentFunctions.h:24` | Equivalent; adds cloud rebuild |
| `search_kink(start_p)` | `segment_search_kink(seg, start_p, name)` | `PRSegmentFunctions.cxx:190` | See §5.3 hard-coded drift direction |
| `reset_associate_points()` | `seg->dpcloud(name, nullptr)` per cloud | — | Callers clear by cloud name |
| `add_associate_point(...)` | `clustering_points_segments(...)` | `PRSegmentFunctions.cxx:1712` | Redesigned for multi-face — improvement |
| `is_shower_trajectory()` | `segment_is_shower_trajectory(seg)` | `PRSegmentFunctions.cxx:~979` | Equivalent |
| `is_shower_topology()` | `segment_is_shower_topology(seg)` | `PRSegmentFunctions.cxx:~2276` | Equivalent |
| `determine_shower_direction()` | `segment_determine_shower_direction(seg)` | `PRSegmentFunctions.cxx` | Equivalent |
| `determine_dir_track()` | `segment_determine_dir_track(seg,...)` | `PRSegmentFunctions.cxx` | Equivalent |
| `do_track_pid()` | `segment_do_track_pid(seg,...)` | `PRSegmentFunctions.cxx` | Equivalent; explicit `ParticleDataSet` DI |
| `cal_kine_range / cal_kine_dQdx` | `segment_cal_kine_range/dQdx` | `PRSegmentFunctions.cxx` | Explicit `IRecombinationModel` DI; improvement over `TPCParams` singleton |
| `cal_4mom()` | `segment_cal_4mom(seg,...)` | `PRSegmentFunctions.cxx` | Equivalent with DI |
| `get_flag_shower()` returning PDG-11 for shower | **Behavior MISSING** | — | Prototype coerced particle type to 11 (electron) in `get_flag_shower()`; toolkit stores PDG in `ParticleInfo` but callers must re-implement the coercion if needed |
| `print_dis()` | **MISSING** | — | Debug helper; not a functionality gap |
| `pcloud_associated_steiner` / `add_associate_point_steiner` | **Intentionally removed** | `PRSegment.h:171` | Comment confirms it was never used — improvement |

### 1.3 WCShower → PR::Shower

| Prototype method | Toolkit equivalent | Toolkit file:line | Notes |
|---|---|---|---|
| `set_start_vertex(vtx, type)` | `set_start_vertex(vtx, type)` | `PRShower.cxx:65` | Equivalent; see §2.7 null-check |
| `set_start_segment(seg)` | `set_start_segment(seg, flag_inc, fit_name, assoc_name)` | `PRShower.cxx:82` | Extended signature for cloud names |
| `set_start_segment(seg, map_seg_vtxs)` | same (graph provides topology) | `PRShower.cxx:82` | Graph replaces explicit maps |
| `complete_structure_with_start_segment(map_vtx_segs, map_seg_vtxs, used_segs)` | `complete_structure_with_start_segment(used_segs, fit_name, assoc_name)` | `PRShower.cxx:265` | Maps eliminated (graph provides topology); see §4.2 pointer-ordered parameter |
| `add_segment(seg, map_seg_vtxs)` | `add_segment(seg, flag_inc_vtxs, fit_name, assoc_name)` | `PRShower.cxx:139` | Equivalent |
| `add_shower(other)` | `add_shower(other, fit_name, assoc_name)` | `PRShower.cxx:187` | DPC merge batched via `add_points` — improvement over `rebuild_point_clouds` |
| `fill_maps(map_vtx_shower, map_seg_shower)` | `fill_maps()` (no-op, no out-params) | `PRShower.cxx:415` | **BUG — see §2.1** |
| `fill_sets(used_vtxs, used_segs)` | `fill_sets(used_vtxs, used_segs, excl_start)` | `PRShower.cxx` | Equivalent |
| `fill_point_vec(pts, flag_main)` | `fill_point_vector(pts, flag_main)` | `PRShower.cxx:371` | Equivalent; see §3.3 missing reserve |
| `get_connected_pieces(seg)` | `get_connected_pieces(seg)` | `PRShower.cxx:456` | Equivalent; returns `IndexedSegmentSet` (deterministic) |
| `calculate_kinematics()` | `calculate_kinematics(particle_data, recomb)` | `PRShower.cxx:863` | DI added; equivalent logic |
| `calculate_kinematics_long_muon(segs)` | `calculate_kinematics_long_muon(IndexedSegmentSet, pd, rm)` | `PRShower.cxx:1185` | Single-pass improvement; DI added |
| `get_last_segment_vertex_long_muon(segs)` | same signature | `PRShower.cxx:504` | See §4.2 pointer-ordered parameter |
| `get_stem_dQ_dx(vtx, seg, limit)` | same | `PRShower.cxx:709` | Equivalent; slight efficiency improvement |
| `get_total_length()` / `get_total_length(cluster_id)` | same | `PRShower.cxx` | Equivalent; see §3.2 cache recommendation |
| `get_total_track_length()` | same | `PRShower.cxx` | Equivalent |
| `get_num_segments()` | same | `PRShower.cxx` | Equivalent |
| `get_num_main_segments()` | same | `PRShower.cxx` | Equivalent |
| `update_particle_type()` | `update_particle_type(particle_data, recomb)` | `PRShower.cxx:655` | DI added — improvement |
| `cal_dir_3vector(p, dis_cut)` | `shower_cal_dir_3vector(shower, p, dis_cut)` | `PRShowerFunctions.cxx:122` | Free function; zero-vector guard added — fixes prototype divide-by-zero (WCShower.cxx:591) |
| `rebuild_point_clouds()` / `build_point_clouds()` | `add_segment` / `add_shower` batch `add_points` | `PRShower.cxx:117, 172` | Superseded — improvement |
| `get_fit_pcloud()` / `get_associated_pcloud()` | `get_pcloud(name)` | `PRShower.h:178-183` | Generalized to named DPCs |
| `get_map_seg_vtxs()` / `get_map_vtx_segs()` | graph traversal via `TrajectoryView` | — | Maps eliminated in favor of graph; `fill_maps` API broken (§2.1) |

---

## §2 Bugs (S1 — Correctness)

### 2.1 `PR::Shower::fill_maps` is a no-op  ⚠ CRITICAL

**Prototype:** `WCShower::fill_maps(map<ProtoVertex*, WCShower*>&, map<ProtoSegment*, WCShower*>&)` populates two caller-supplied global lookup tables that map every vertex/segment in the shower back to the shower owning it.

**Toolkit:** `PRShower.cxx:415-417`
```cpp
TrajectoryView& Shower::fill_maps() {
    return *this;
}
```
The method takes no out-parameters and does nothing. `ShowerVertexMap` and `ShowerSegmentMap` (declared in `PRShower.h:233-234`) are defined as `std::map<VertexPtr, ShowerPtr, VertexIndexCmp>` and `std::map<SegmentPtr, ShowerPtr, SegmentIndexCmp>` but the method to populate them is a stub.

**Impact:** Any caller that relied on `fill_maps` to build a vertex→shower or segment→shower lookup (e.g. to find which shower a given segment belongs to) silently gets empty maps. Likely affects downstream PID and neutrino tagger logic that navigates shower ownership.

**Fix:** Implement `fill_maps` with two out-parameter maps, or provide a separate `fill_shower_maps(ShowerVertexMap&, ShowerSegmentMap&)` that iterates `this->edges()` and `this->nodes()` and inserts `{ptr → this->shared_from_this()}` pairs. Alternatively, remove the method and require callers to build these maps themselves using `TrajectoryView` iteration — but that requires auditing all call sites.

---

### 2.2 `PR::Vertex` inherits `HasCluster<Segment>` instead of `HasCluster<Vertex>`  ⚠ BUG

**File:** `PRVertex.h:48`
```cpp
class Vertex
    : public Flagged<VertexFlags>
    , public Graphed<node_descriptor>
    , public HasCluster<Segment>      // ← should be HasCluster<Vertex>
```
`HasCluster<Subclass>::cluster(Facade::Cluster* cptr)` at `PRCommon.h:35` does:
```cpp
return *dynamic_cast<Subclass*>(this);
```
With `Subclass=Segment`, `dynamic_cast<Segment*>(vertex_ptr)` returns `nullptr` (different class hierarchy), so the chainable setter `vertex->cluster(ptr)` is broken. The non-chainable getter `vertex->cluster()` works fine since it only reads `m_cluster`.

**Impact:** Only manifests when a caller chains `.cluster(ptr)` on a `Vertex`; currently harmless. Becomes a silent null-deref if caller stores the returned `Segment*` reference.

**Fix:** `class Vertex : public HasCluster<Vertex>`

---

### 2.3 `find_vertices` dereferences `wcpts().front()` without an empty guard  ⚠ BUG

**File:** `PRGraph.cxx:87`
```cpp
auto ept = seg->wcpts().front().point;
```
If a segment was constructed with fits but no wcpts (or if wcpts were cleared), this is undefined behavior.

**Fix:** Add a guard:
```cpp
if (seg->wcpts().empty()) {
    // fall back to the fit point
    ept = seg->fits().empty() ? Point{} : seg->fits().front().point;
}
```

---

### 2.4 `segment_track_direct_length` clamps `n2` inconsistently  ⚠ BUG

**File:** `PRSegmentFunctions.cxx:739` vs `PRSegmentFunctions.cxx:775`

`segment_track_direct_length` clamps a standalone negative `n2` to 0:
```cpp
if (n2 < 0) n2 = 0;   // line 739
```
`segment_track_max_deviation` clamps to `fits.size()-1`:
```cpp
if (n2 < 0) n2 = static_cast<int>(fits.size()) - 1;   // line 775
```
The pair `(n1 < 0 && n2 < 0)` is handled identically (both default to `[0, fits.size()-1]`), but calling `segment_track_direct_length(seg, 5, -1)` silently yields the distance from index 5 to index 0 (backwards), while `segment_track_max_deviation(seg, 5, -1)` yields the expected full-range deviation.

**Impact:** Callers that pass `n2=-1` as a sentinel for "end of segment" to `segment_track_direct_length` get wrong results. Check prototype `ProtoSegment::get_direct_length(n1,n2)` for intended default.

**Fix:** Align to `segment_track_max_deviation`'s convention (line 739 → `n2 = static_cast<int>(fits.size()) - 1`).

---

### 2.5 `Segment::clear_fit` drops per-fit `paf` (APA/face)  ⚠ BUG (multi-APA)

**File:** `PRSegment.cxx:87-99`

`clear_fit` rebuilds `m_fits` from `m_wcpts` but does not populate `fit.paf`:
```cpp
for (size_t i = 0; i != m_wcpts.size(); i++) {
    m_fits.at(i).point = m_wcpts.at(i).point;
    m_fits.at(i).dQ = -1;
    // ... other fields ...
    // m_fits.at(i).paf is never set → stays {-1,-1}
}
```
After `clear_fit`, every `fit.paf = {-1,-1}`, breaking any downstream code that uses `fit.paf` to select wire-plane projections (e.g. `TrackFitting.cxx:6201-6208`, which compares `paf` between adjacent fits to detect face crossings).

**Fix:** Populate `paf` during `clear_fit` using `dv->contained_by(m_wcpts[i].point)`:
```cpp
auto wpid = dv->contained_by(m_wcpts.at(i).point);
m_fits.at(i).paf = {wpid.apa(), wpid.face()};
```

---

### 2.6 `get_angles_cached` first-hit short-circuit may cache `(0,0,0)`  ⚠ BUG (multi-APA)

**File:** `PRSegmentFunctions.cxx:1755-1768`

The lambda iterates `seg_dpc_cache` (a pointer-ordered `std::map`) and breaks on the first entry that returns a non-zero angle triple for the requested `(apa, face)`:
```cpp
for (auto& [s, dpc] : seg_dpc_cache) {
    a = dpc->get_angles(face, apa);
    if (a[0] != 0.0 || a[1] != 0.0 || a[2] != 0.0) break;
}
ang_cache[key] = a;
```
If the first segment in address order does not have the requested `(apa, face)` in its `wpid_params`, `get_angles` returns `(0,0,0)` and the loop moves on. But if *no* segment has it, `(0,0,0)` is stored and all subsequent 2D-distance calculations for that `(apa, face)` silently return garbage (projection angles of zero degenerate 2D coordinates).

This is also a determinism issue since the iteration order is pointer-order (see §4.1).

**Fix:** Change the loop to iterate until a genuinely populated `(apa, face)` entry is found, or collect all `wpid_params` into a global map once at the start of `clustering_points_segments`.

---

### 2.7 `calculate_kinematics` `(0,0,0)` sentinel not handled in two branches  ⚠ BUG

**File:** `PRShower.cxx`

`shower_get_closest_point` returns `{-1, (0,0,0)}` when the shower has no "fit" DPC. The nseg==1/conn_type!=1 branch handles this at `PRShower.cxx:937-942`:
```cpp
if (data.start_point.x() == 0 && data.start_point.y() == 0 && data.start_point.z() == 0) {
    // fallback to first fit point
}
```
But the `nseg>1 / nseg==nconnected` branch at PRShower.cxx:~1009 and the `else` (multi-track) branch at ~1088 both call `shower_get_closest_point` without checking for the `(0,0,0)` sentinel. If those branches are reached with an empty "fit" cloud, `data.start_point` is left at the origin `(0,0,0)`.

**Fix:** Apply the same fallback pattern to the two unchecked branches.

---

### 2.8 `set_start_vertex` and `calculate_kinematics_long_muon` missing null checks  ⚠ MINOR

- `PRShower.cxx:67`: `vtx->descriptor_valid()` — if `vtx` is `nullptr`, this is a null dereference. Add `if (!vtx)` guard before accessing any member.
- `PRShower.cxx:1187`: `m_start_segment->particle_info()->pdg()` — both `m_start_segment` and `particle_info()` can be null. This is inherited from prototype `WCShower.cxx`. Add null guards.

---

### 2.9 `create_segment_fit_point_cloud` silently drops a fit at the origin  ⚠ MINOR

**File:** `PRSegment.cxx:84` (approximate — inside `create_segment_fit_point_cloud`):

A fit is included only if `fit.valid() || (fit.point.x() != 0 || ...)`. A legitimate `(0,0,0)` fit point at the detector origin is silently excluded. Extremely unlikely in real data but a latent edge-case.

**Fix:** Remove the `(0,0,0)` check; rely on `fit.valid()` (index ≥ 0) alone.

---

## §3 Efficiency (S2)

### 3.1 `calculate_kinematics` double-pass over shower edges  ⚠ EFFICIENCY

`PRShower.cxx:994-1071` and `1073-1158` are nearly identical twin branches (`nseg==nconnected` vs else). Each iterates `ordered_edges(*this, m_full_graph)` twice — once to compute `total_length` and once to fill `vec_dQ/vec_dx`. These should be fused:

```cpp
double total_length = 0;
std::vector<double> vec_dQ, vec_dx;
for (auto edesc : ordered_edges(*this, m_full_graph)) {
    // ... accumulate total_length, vec_dQ, vec_dx in one pass ...
}
```

### 3.2 `get_total_length`, `get_num_segments`, `get_num_main_segments` recomputed on every call  ⚠ EFFICIENCY

These traverse all shower edges every call. In `calculate_kinematics` they are called multiple times. Cache the result in `ShowerData` after structure finalization (analogous to prototype's in-memory cached values). Mark cache dirty in `add_segment`, `add_shower`, `complete_structure_with_start_segment`.

### 3.3 `fill_point_vector` missing `reserve`  ⚠ MINOR

`PRShower.cxx:~371`: `points.push_back(...)` in a loop with no prior `reserve`. Calls to `fill_point_vector` in kinematics allocate repeatedly. Add:
```cpp
points.reserve(this->edges_count() * 50); // rough estimate
```

### 3.4 `clustering_points_segments` ghost-removal inner loop is O(N·S·P)  ⚠ EFFICIENCY

`PRSegmentFunctions.cxx:~1896-1905`: for each query point, the ghost-removal loop scans all segments globally. Building a single KD-tree over all fit points once at the start of `clustering_points_segments` would reduce this to O(N log N).

### 3.5 `Segment::set_fit_associate_vec` copies the fit vector  ⚠ MINOR

`PRSegment.cxx:39`: `m_fits = tmp_fit_vec;` is a full copy. Change the parameter to `std::vector<PR::Fit>&& tmp_fit_vec` and use `m_fits = std::move(tmp_fit_vec)` to avoid the copy. Update callers accordingly.

### 3.6 Redundant `boost::edge` check in `find_vertices`  ⚠ MINOR

`PRGraph.cxx:81-82`: after obtaining `ed = seg->get_descriptor()` and extracting source/target, the code re-calls `boost::edge(vd1, vd2, graph)` to verify the edge is in the graph. Since `ed` is already a valid descriptor (guarded by `seg->descriptor_valid()` at line 74), this re-check is redundant.

### 3.7 `add_shower` already improved over prototype  ✓ IMPROVEMENT

`PRShower.cxx:117, 172` batch DPC points via `add_points` rather than the prototype's full `rebuild_point_clouds` on every segment addition (`WCShower.cxx:681-692`). Good.

### 3.8 `shower_cal_dir_3vector` zero-vector early return  ✓ IMPROVEMENT

`PRShowerFunctions.cxx:153` returns early with a zero vector when `ncount == 0`. Prototype `WCShower.cxx:591` would divide by zero in the same case.

---

## §4 Determinism (S3 — pointer-keyed containers)

### 4.1 `clustering_points_segments` pointer-keyed maps  ⚠ DETERMINISM BUG

**File:** `PRSegmentFunctions.cxx:1717, 1730-1732`

```cpp
std::map<Facade::Cluster*, std::vector<SegmentPtr>> map_cluster_segs;  // line 1717 — address-ordered
std::map<SegmentPtr, ...>  seg_dpc_cache;                               // line 1730 — address-ordered
std::map<SegmentPtr, ...>  seg_pts3d;                                   // line 1731 — address-ordered
std::map<SegmentPtr, ...>  seg_pts2d;                                   // line 1732 — address-ordered
```

`map_cluster_segs` is iterated at line ~1803 in a loop that calls `create_segment_point_cloud` with side effects (modifying each segment's DPC). The iteration order depends on `Cluster*` addresses, which vary between runs. This is a real determinism bug.

The three `SegmentPtr`-keyed maps are used primarily for lookup (`.find`), which is address-order-independent. However, `seg_dpc_cache` is iterated in `get_angles_cached` (§2.6), coupling address order into the cached angle result.

**Fix:**
```cpp
std::map<Facade::Cluster*, std::vector<SegmentPtr>, ClusterPtrCmp> map_cluster_segs; // ClusterPtrCmp is already defined in PRShower.h:238
```
For `seg_dpc_cache`/`seg_pts3d`/`seg_pts2d`, use `SegmentIndexCmp` as the comparator or key by `seg->get_graph_index()` (integer). Both comparators are already defined in `PRShower.h:29-40`.

### 4.2 `complete_structure_with_start_segment` and `get_last_segment_vertex_long_muon` take pointer-ordered `std::set<SegmentPtr>`  ⚠ DETERMINISM

**File:** `PRShower.h:188, 199`

```cpp
void complete_structure_with_start_segment(std::set<SegmentPtr>& used_segments, ...);           // line 188 — pointer-ordered
std::pair<SegmentPtr, VertexPtr> get_last_segment_vertex_long_muon(std::set<SegmentPtr>& segs); // line 199 — pointer-ordered
```

These parameters accept the default `std::set<SegmentPtr>` (pointer comparator), while everywhere else in `PRShower.h` the corresponding type is `IndexedSegmentSet` (with `SegmentIndexCmp`). Callers building these sets via pointer order will see non-deterministic iteration in external loops.

The methods themselves only use `find`/`insert` on these sets, so the *internal* behavior is deterministic, but having an inconsistent API invites callers to pass sets with the wrong comparator.

**Fix:** Change both parameters to `IndexedSegmentSet&`.

Note: `calculate_kinematics_long_muon` already uses `IndexedSegmentSet&` (`PRShower.h:214`) — the inconsistency is clear.

### 4.3 Internal `unordered_set` usage — decision-neutral  ✓ OK

`PRShower.cxx:422-423, 459-460, 512`: `std::unordered_set<SegmentPtr>` and `std::unordered_set<VertexPtr>` are used inside `count_connected_segments`, `get_connected_pieces`, and `get_last_segment_vertex_long_muon` for membership tests only. The returned `IndexedSegmentSet`/`IndexedVertexSet` are always index-ordered, so the API is deterministic. The internal unordered sets are benign.

### 4.4 Correctly deterministic patterns  ✓ GOOD

- `VertexIndexCmp` / `SegmentIndexCmp` / `ShowerIndexCmp` defined at `PRShower.h:29-40` and used in all public-facing containers.
- `ordered_nodes` / `ordered_edges` / `sorted_out_edges` used in all farthest-vertex decision loops (`PRShower.cxx:947, 1030, 1128`, `PRShower.cxx:296, 356, 364, 384, 401`).
- `std::map<size_t, VertexPtr>` (keyed on `graph_index`) in `calculate_kinematics_long_muon` at `PRShower.cxx:1200`.
- Angle cache in `clustering_points_segments` keyed on `std::pair<int,int>` at `PRSegmentFunctions.cxx:1754`.

---

## §5 Multi-APA / Multi-Face Audit

### 5.1 Face resolution API  ✓ CONSISTENT

All face/APA assignment goes through `IDetectorVolumes::contained_by(Point) → WirePlaneId` (`aux/src/DetectorVolumes.cxx:336`). This is called in `create_segment_point_cloud` (`PRSegmentFunctions.cxx:32`) per path point and in `create_segment_fit_point_cloud` similarly. The resulting `WirePlaneId` is stored per DPC point and per `Fit::paf`.

### 5.2 Correctly multi-face aware spots  ✓ GOOD

- `Fit::paf{apa, face}` stored per-point (`PRCommon.h:122`). Allows `TrackFitting.cxx` to detect face crossings.
- `create_segment_point_cloud` builds per-`(apa,face)` kd-trees and projection angles from `cluster.grouping()->wpids()` (`PRSegmentFunctions.cxx:37-43`). Each segment correctly supports multiple faces.
- Fast 2D distance path in `clustering_points_segments` (`PRSegmentFunctions.cxx:1744-1768, 1873-1874`): caches per-`(apa,face)` angles and per-segment 2D projections keyed by `(plane, apa, face)`, returns `{1e9, 1e9, 1e9}` for cross-APA queries. This is the most thoroughly multi-face-correct code in the file.
- `PR::Shower` contains zero face/anode references — coordinate-agnostic. All face logic is delegated to the segment helpers. Correct by design.

### 5.3 Suspicious / potentially single-face  ⚠ AUDIT

**5.3.1 `segment_get_closest_2d_distances` uses a single `(apa,face)` caller-supplied**  
`PRSegmentFunctions.cxx:137-188`: the function accepts one `apa` and one `face` and queries kd2d for all three planes at that face only. If a segment spans two faces (which the toolkit now permits), the portion of the segment on the other face is invisible.

This is a direct port of the prototype's single-face `get_closest_2d_dis(Point)`. The prototype had no multi-face issue because only one face existed.

**Fix:** Either iterate all `(apa,face)` entries in the DPC's `wpid_params` and return the minimum, or update the caller to pass the correct face for each query point.

**5.3.2 `segment_search_kink` hard-codes the drift direction**  
`PRSegmentFunctions.cxx:194`:
```cpp
WireCell::Vector drift_dir_abs(1,0,0);
```
Prototype `ProtoSegment.cxx` did the same (single detector). For a multi-APA detector where opposing TPC volumes have anti-parallel drift directions (e.g. SBND, DUNE FD), this gives wrong kink detection for one of the two drift orientations.

**Fix:** Resolve the drift direction from `IDetectorVolumes::contained_by(seg start/end point)` and use the APA's actual drift axis.

**5.3.3 `Segment::clear_fit` drops `paf`**  
Already flagged as §2.5. This is also a multi-face correctness issue: downstream code that reads `fit.paf` to select wire-plane geometry after `clear_fit`+refit will read `{-1,-1}` and use incorrect projections.

**5.3.4 `Shower::set_start_segment` and `add_segment` share the first segment's DPC pointer**  
`PRShower.cxx:111-113, 166-168` (explicit TODO comments in the code):
```cpp
// For now, just share the pointer - we'll need to modify this if independent DPCs are needed
this->dpcloud(cloud_name_fit, seg_dpc_fit);
```
When a shower is composed of segments from different APAs/faces, the shower-level DPC inherits only the first segment's `wpid_params`. Any subsequent call to `shower_dpc->get_angles(face, apa)` for a face not present in the first segment returns `(0,0,0)`.

**Fix:** Expose a `wpid_params()` accessor on `DynamicPointCloud` and construct a genuine merged shower DPC whose `wpid_params` is the union of all constituent segment DPCs.

**5.3.5 dQ/dx aggregators are face-blind**  
`segment_median_dQ_dx`, `segment_rms_dQ_dx`, `segment_cal_kine_dQdx`, `segment_is_shower_trajectory`, `Shower::get_stem_dQ_dx`: all accumulate `fit.dQ/fit.dx` across all fits without checking `fit.paf`. If different APAs have different gain calibrations, a segment crossing a face boundary will silently mix calibration scales in its dQ/dx estimate.

In current MicroBooNE usage (one face), this is not an issue. For SBND/DUNE FD, this should be addressed before physics analysis.

**5.3.6 `DetectorVolumes::contained_by` tie-break is insertion-order dependent**  
`aux/src/DetectorVolumes.cxx:336-407`: when a 3D point lies near a face boundary, the grid cell may contain BBs from two faces. The first matching BB in `m_grid[i][j][k]` wins, and that order is insertion-order. For PR code this means a fit point at a geometric face boundary may be assigned to either face in a run-to-run non-deterministic fashion.

**Fix:** Document (or enforce) a deterministic tie-break rule, e.g. prefer the face with lower APA index, or prefer the face whose center is closer to the point.

### 5.4 Single-face sub-case preserved  ✓ CORRECT

For a detector with exactly one `(apa, face)`, all the new multi-face machinery degenerates gracefully to the prototype:
- `create_segment_point_cloud` builds one entry in `wpid_params`, producing behavior identical to the prototype's single `pcloud_fit`.
- `clustering_points_segments`'s fast-2D path finds only one `(apa, face)` per query and returns the same distance.
- `seg_dpc_cache` and `map_cluster_segs` have the same iteration-order as prototype since there is effectively one face bucket.
- Shower DPC merge sharing the first (and only) segment's `wpid_params` is a no-op difference.
- `drift_dir_abs(1,0,0)` is correct for MicroBooNE/SBND TPC1.

---

## §6 Called Helpers — Correctness Notes

| Helper | Toolkit file | Assessment |
|---|---|---|
| `create_segment_fit_point_cloud` | `PRSegmentFunctions.cxx:25` | Correct. Multi-face aware. See §2.9 origin-drop edge case. |
| `create_segment_point_cloud` | `PRSegmentFunctions.cxx:~1600` | Correct. Populates per-`(apa,face)` kd-trees. |
| `clustering_points_segments` | `PRSegmentFunctions.cxx:1712` | Good algorithm; has determinism bug (§4.1) and `get_angles_cached` bug (§2.6). |
| `segment_track_length` | `PRSegmentFunctions.cxx:~650` | Correct. Face-independent (geometric). |
| `segment_track_direct_length` | `PRSegmentFunctions.cxx:~720` | n2 clamp inconsistency (§2.4). |
| `segment_track_max_deviation` | `PRSegmentFunctions.cxx:759` | Correct. |
| `segment_median_dQ_dx / _rms` | `PRSegmentFunctions.cxx:~823,867` | Correct for single-face; face-blind (§5.3.5). |
| `segment_search_kink` | `PRSegmentFunctions.cxx:190` | Hard-coded drift direction (§5.3.2). |
| `segment_do_track_pid / _comp` | `PRSegmentFunctions.cxx` | Port of prototype logic; KS-ratio thresholds inherited. |
| `segment_cal_kine_range / _dQdx / _4mom` | `PRSegmentFunctions.cxx` | Correct. Explicit DI is improvement. |
| `segment_determine_dir_track / _shower` | `PRSegmentFunctions.cxx` | Correct port; 1:1 with prototype thresholds. |
| `shower_cal_dir_3vector` | `PRShowerFunctions.cxx:122` | Correct. Zero-vector guard added vs prototype. |
| `shower_get_closest_point` | `PRShowerFunctions.cxx:9` | Correct for single-face; see §2.7 for sentinel handling. |
| `find_vertices` | `PRGraph.cxx:72` | See §2.3 empty-wcpts bug. |
| `find_other_vertex` | `PRGraph.cxx:98` | Correct; null-safe. |
| `add_vertex / add_segment` | `PRGraph.cxx:5,28` | Correct; stable graph indices assigned. Redundant `boost::edge` check (§3.6). |

---

## §7 Priority-Ordered Fix List

| # | Severity | Location | Fix |
|---|---|---|---|
| F1 | CRITICAL | `PRShower.cxx:415-417` | Implement `fill_maps` to populate `ShowerVertexMap` / `ShowerSegmentMap` out-params (or refactor API) |
| F2 | BUG | `PRSegmentFunctions.cxx:1717` | Replace `std::map<Cluster*, ...>` with `ClusterPtrCmp`-keyed map |
| F3 | BUG | `PRSegmentFunctions.cxx:1730-1732` | Replace `std::map<SegmentPtr, ...>` with `SegmentIndexCmp`-keyed map; also fixes `get_angles_cached` non-determinism |
| F4 | BUG | `PRSegmentFunctions.cxx:1755-1768` | Fix `get_angles_cached` to not rely on iteration order for `(0,0,0)` detection |
| F5 | BUG | `PRSegment.cxx:87-99` | Populate `fit.paf` in `clear_fit` via `dv->contained_by` |
| F6 | BUG | `PRSegmentFunctions.cxx:739` | Align `segment_track_direct_length` n2 default clamp to `fits.size()-1` |
| F7 | BUG | `PRGraph.cxx:87` | Guard `wcpts().front()` against empty wcpts |
| F8 | BUG | `PRVertex.h:48` | Change `HasCluster<Segment>` to `HasCluster<Vertex>` |
| F9 | BUG | `PRShower.cxx:~1009, ~1088` | Add `(0,0,0)` sentinel check for `shower_get_closest_point` in two missing branches |
| F10 | BUG | `PRShower.h:188, 199` | Change `std::set<SegmentPtr>&` parameters to `IndexedSegmentSet&` |
| F11 | BUG | `PRShower.cxx:67, 1187` | Add null guards for `vtx` and `m_start_segment->particle_info()` |
| F12 | MULTI-APA | `PRSegmentFunctions.cxx:194` | Resolve drift direction from `IDetectorVolumes` per-segment |
| F13 | MULTI-APA | `PRSegmentFunctions.cxx:137-188` | `segment_get_closest_2d_distances` — iterate all `(apa,face)` in DPC |
| F14 | MULTI-APA | `PRShower.cxx:111-113, 166-168` | Implement proper DPC merge in `set_start_segment` / `add_segment` (resolve the TODO) |
| F15 | EFFICIENCY | `PRShower.cxx:994-1158` | Fuse twin edge passes in `calculate_kinematics` |
| F16 | EFFICIENCY | `PRShower.cxx` | Cache `get_total_length` / `get_num_segments` / `get_num_main_segments` |
| F17 | EFFICIENCY | `PRSegmentFunctions.cxx:1896-1905` | Global KD-tree for ghost-removal inner loop |
| F18 | MINOR | `PRShower.cxx:~371` | Add `reserve` in `fill_point_vector` |
| F19 | MINOR | `PRSegment.cxx:39` | Move-semantics for `set_fit_associate_vec` |
| F20 | MINOR | `PRSegment.cxx:84` | Remove `(0,0,0)` fit-point rejection; use `fit.valid()` only |
| F21 | MINOR | `PRGraph.cxx:81` | Remove redundant `boost::edge` re-check in `find_vertices` |
| F22 | DOC | `clus/docs/porting/porting_dictionary.md:211-217` | Fill in WCT equivalents for ProtoVertex / ProtoSegment / WCShower stubs |

---

## §8 Verification Plan

1. **Unit tests** (`clus/test/doctest_prvertex.cxx`, `clus/test/doctest_prsegment.cxx`):
   - Add a test that calls `vertex->cluster(ptr)` and verifies the returned reference is a `Vertex&`, not garbage. This will catch the `HasCluster<Segment>` CRTP fix.
   - Add a test constructing a segment with fits but no wcpts and calling `find_vertices` — verifies the empty-guard fix.
   - Add a test calling `clear_fit` on a multi-face segment and checking that `fit.paf != {-1,-1}` afterwards.
   - Add a test for `segment_track_direct_length(seg, 5, -1)` vs `segment_track_max_deviation(seg, 5, -1)` — they should both interpret `n2=-1` as the last index.

2. **Integration** (single-face regression):
   - Run the full `do_multi_tracking` → `PRShower::calculate_kinematics` pipeline on a MicroBooNE sample before and after applying fixes F1–F21.
   - Compare `ShowerData` output (kenergy_range, kenergy_dQdx, start_point, end_point) and verify no regression in reconstructed kinematic distributions.

3. **Determinism test**:
   - Run the pipeline twice on the same event with different heap allocation patterns (or `LD_PRELOAD` a shuffling allocator).
   - Diff `ShowerData` and `IndexedSegmentSet` membership between the two runs.
   - Fixes F2–F4 and F10 must eliminate all diffs.

4. **Multi-face smoke test** (after F5, F13, F14):
   - Construct a synthetic segment spanning two simulated APA faces, run `create_segment_point_cloud`, and verify both `(apa0,face0)` and `(apa1,face1)` kd-trees are populated.
   - Call `segment_get_closest_2d_distances` with a query point on each face and verify sensible distances are returned for both.
