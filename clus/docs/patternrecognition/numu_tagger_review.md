# NuMu Tagger, BDT Scorer & match\_isFC Review

Review of the following ported functions:

| Prototype (WCP) | Toolkit (WCT) | Source file |
|---|---|---|
| `NeutrinoID::numu_tagger()` | `PatternAlgorithms::numu_tagger()` | `clus/src/NeutrinoTaggerNuMu.cxx:161` |
| `NeutrinoID::count_daughters(ProtoSegment*)` | `PatternAlgorithms::count_daughters(Graph&, SegmentPtr, VertexPtr)` | `clus/src/NeutrinoTaggerNuMu.cxx:82` |
| `NeutrinoID::count_daughters(WCShower*)` | `PatternAlgorithms::count_daughters(Graph&, ShowerPtr, VertexPtr, IndexedSegmentSet&)` | `clus/src/NeutrinoTaggerNuMu.cxx:118` |
| `NeutrinoID::find_cont_muon_segment()` | `PatternAlgorithms::find_cont_muon_segment()` | `clus/src/NeutrinoVertexFinder.cxx:920` |
| `NeutrinoID::cal_numu_bdts_xgboost()` | `UbooneNumuBDTScorer::cal_numu_bdts_xgboost()` | `root/src/UbooneNumuBDTScorer.cxx:275` |
| `match_isFC` (via `fid->check_fully_contained`) | `Facade::cluster_fc_check()` | `clus/src/Clustering_Util.cxx:56` |

Prototype sources: `NeutrinoID_numu_tagger.h`, `NeutrinoID_numu_bdts.h`,
`NeutrinoID_track_shower.h:2250`, `prod-nusel-eval.cxx:847`.

---

## 1. Functional Equivalence

### 1.1 `numu_tagger()`

**Verdict: Equivalent.**

The function identifies numu CC events via three checks (flag 1: direct muon at
main vertex, flag 2: long-muon shower, flag 3: indirect muon in main cluster)
and returns `{flag_long_muon, max_muon_length}`.

| Aspect | Prototype | Toolkit | Assessment |
|---|---|---|---|
| Flag 1 loop | `map_vertex_segments[main_vertex]` (pointer-keyed map) | `boost::out_edges(vd, graph)` | Equivalent; toolkit is deterministic |
| Flag 2 loop | `showers` iterator | Range-for over `IndexedShowerSet` | Equivalent; toolkit iteration order is index-based |
| Flag 3 loop | `map_segment_vertices` iterator | `boost::edges(graph)` filtered by cluster id | Equivalent |
| dQ/dx cut formula | `0.8866 + 0.9533 * pow(18cm/length, 0.4234)` | Same | Identical |
| Daughter cut | `!(n_daughter_tracks > 1 \|\| n_daughter_all - n_daughter_tracks > 2)` | Same | Identical |
| Flag 3 condition (line 138 prototype) | `sg->get_flag_shower() && !sg->get_flag_shower_topology() \|\| !sg->get_flag_shower() \|\| length > 50cm` | `!is_shower_topo \|\| length > 50cm` | Toolkit correctly simplifies (since `kShowerTopology` implies `is_shower`) |
| Flag 3 threshold (line 181 prototype) | `max_length > 25cm && acc_track_length > 0 \|\| max_length > 30cm && acc_track_length == 0` | Same with explicit parentheses | Equivalent; toolkit clarifies precedence |
| `neutrino_type` bits | Sets `1UL<<2` (numu) or `1UL<<3` (NC) | Omitted; caller uses returned `flag_numu_cc` bool | Intentional design change; no functional impact |
| Unused variables | `dir_beam`, `dir_drift`, `dir_vertical`, `flag_print`, `max_direct_length`, `max_medium_dQ_dx` | Removed | Cleanup |
| Muon extension (flag 3, >50 cm) | `find_cont_muon_segment(sg, pair_vertices.first)` / `.second` | Same via `find_vertices(graph, sg)` then `find_cont_muon_segment(graph, sg, pv1)` / `pv2` | Equivalent |

### 1.2 `count_daughters()` (both overloads)

**Verdict: Equivalent.**

| Aspect | Prototype | Toolkit | Assessment |
|---|---|---|---|
| Distance to main vertex | `sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))` | `ray_length(Ray{vtx_fit_pt(v1), mv_pt})` | Equivalent |
| Vertex position | `get_fit_pt()` directly | `vtx_fit_pt()` which falls back to `wcpt().point` if fit invalid | Toolkit is safer |
| BFS via `calculate_num_daughter_tracks` | Uses pointer-keyed map iteration | Uses Boost graph iteration | Equivalent; toolkit deterministic |
| Unused `muon_dir` TVector3 | Computed in both overloads | Removed | Cleanup |
| Shower overload: last segment | `get_last_segment_vertex_long_muon(segments_in_long_muon)` | Same | Identical |

### 1.3 `find_cont_muon_segment()`

**Verdict: Equivalent + bug fix.**

| Aspect | Prototype | Toolkit | Assessment |
|---|---|---|---|
| Direction vectors | `sg->cal_dir_3vector(vtx->get_fit_pt(), 15cm)` using TVector3 | `segment_cal_dir_3vector(sg, vtx_point, 15cm)` using WireCell::Vector | Equivalent |
| Angle calculation | `(pi - dir1.Angle(dir2)) / pi * 180` | `(pi - acos(clamp(dot/mag, -1, 1))) / pi * 180` | Toolkit adds `std::clamp` before `acos()` (line 963, 976) |
| Zero-magnitude guard | None | `if (dir1.magnitude() == 0 \|\| dir2.magnitude() == 0) continue` (line 958) | Toolkit prevents division by zero |
| Angle thresholds | `angle < 10 \|\| angle1 < 10 \|\| sg_length < 6cm && (angle < 15 \|\| angle1 < 15)` | Same with explicit parentheses | Equivalent; toolkit clarifies precedence |
| dQ/dx threshold | `ratio < 1.3 \|\| flag_ignore_dQ_dx` | Same | Identical |
| Selection criterion | `length * cos(angle/180*pi)` (max projected length) | Same | Identical |
| 50 cm direction check | `dir3 = sg->cal_dir_3vector(vtx, 50cm)` | Same | Identical |

### 1.4 `cal_numu_bdts_xgboost()` and sub-BDTs

**Verdict: Equivalent.**

| Aspect | Prototype | Toolkit | Assessment |
|---|---|---|---|
| Sub-BDT default values | `cal_numu_1_bdt(-0.4)`, `cal_numu_2_bdt(-0.1)`, `cal_numu_3_bdt(-0.2)`, `cal_cosmict_10_bdt(0.7)` | Same | Identical |
| TMVA variable order | 65+ variables in specific order | Same order verified line-by-line | Identical |
| `numu_cc_3_track_length` alias | `reader.AddVariable("numu_cc_3_track_length", &tagger_info.numu_cc_3_acc_track_length)` | Same | Identical |
| NaN guards | `cosmict_4_angle_beam`, `cosmict_7_angle_beam`, `cosmict_7_theta`, `cosmict_7_phi` | Same | Identical |
| Infinity guard | `numu_cc_1_dQ_dx_cut`: if `isinf` set to 10 | Same | Identical |
| `numu_cc_flag_1` loaded but not used as BDT input | Loaded into local but not added to reader | `(void)ti.numu_cc_flag_1.at(i)` | Consistent |
| `numu_3_score` | Computed but not added to final xgboost reader | Same | Consistent |
| Log-likelihood conversion | `TMath::Log10((1+val1)/(1-val1))` | `std::log10((1.0+val1)/(1.0-val1))` | Equivalent |
| Weight file path | Hardcoded `"input_data_files/numu_scalars_scores_0923.xml"` | Configurable via `m_numu_xgboost_xml` | Toolkit improvement |
| Return type | `float` return value | `void`; writes `ti.numu_score` | Equivalent |

### 1.5 `match_isFC` / `cluster_fc_check()`

**Verdict: Equivalent with multi-APA extension.**

| Aspect | Prototype | Toolkit | Assessment |
|---|---|---|---|
| FC determination | `fid->check_fully_contained(bundle, offset_x, ct_point_cloud, old_new_cluster_map, &fc_breakdown)` | `Facade::cluster_fc_check(*main_cluster, m_dv)` | Same two-round boundary check logic |
| Result type | `bool match_isFC` | `float match_isFC` (1.0f or 0.0f) | Compatible; TMVA needs float |
| Round 1 | Steiner boundary `flag_cosmic=true` | Same | Identical |
| Round 2 | Steiner boundary `flag_cosmic=false` (if round 1 finds no exit) | Same | Identical |
| Boundary checks | Direct fiducial, signal-processing (wire-collinear), dead-volume (transverse) | Same | Identical |
| Wire-plane geometry | Single APA assumed | Per-point APA/face via `dv->contained_by(p1)` and per-WPID lookups | Multi-APA extension |
| Two-endpoint protection | If both endpoints exit, re-check with direct fiducial | Same | Identical |
| Set in toolkit | N/A | `TaggerCheckNeutrino.cxx:382` | Correct |

---

## 2. Bugs Found

### BUG-1 (Medium): Unclamped `acos()` in `cluster_fc_check()`

**File**: `clus/src/Clustering_Util.cxx`, lines 143-159.

Six `acos()` calls compute projected wire-plane angles without clamping the
argument to [-1, 1]:

```cpp
double angle1 = acos(dir_1.dot(U_dir) / (dir_1.magnitude() * U_dir.magnitude()));
```

**Risk**: If `dir_1` has zero magnitude (direction is purely along the drift
axis, i.e. `dir.y() == 0 && dir.z() == 0`), this produces `0/0 = NaN`.
Even with non-zero magnitude, floating-point rounding can produce arguments
slightly outside [-1, 1], resulting in NaN from `acos()`.

Note that the toolkit already uses `std::clamp` correctly in
`find_cont_muon_segment()` (`NeutrinoVertexFinder.cxx:963`), so this is an
inconsistency within the toolkit.

**Fix**: Add a zero-magnitude guard for `dir_1` and use `std::clamp` before
each `acos()`:

```cpp
if (dir_1.magnitude() == 0) continue;  // purely drift-aligned, skip wire-angle checks

double angle1 = acos(std::clamp(dir_1.dot(U_dir) / (dir_1.magnitude() * U_dir.magnitude()), -1.0, 1.0));
```

Apply to all six `acos()` calls at lines 143, 146, 149, 152, 155, 158.

### BUG-2 (Low): Unclamped `acos()` for `main_angle`

**File**: `clus/src/Clustering_Util.cxx`, line 171.

```cpp
double main_angle = acos(dir_vec.dot(main_dir) /
                        (dir_vec.magnitude() * main_dir.magnitude()));
```

Same risk as BUG-1 but lower likelihood since `dir_vec` comes from
`vhough_transform` (unlikely to be zero).

**Fix**: Add `std::clamp`:

```cpp
double main_angle = acos(std::clamp(dir_vec.dot(main_dir) /
                        (dir_vec.magnitude() * main_dir.magnitude()), -1.0, 1.0));
```

---

## 3. Efficiency

### EFF-1 (High): TMVA::Reader recreation per event

**File**: `root/src/UbooneNumuBDTScorer.cxx`, lines 104-401.

Each call to `cal_numu_bdts_xgboost()` creates **5 TMVA::Reader** objects and
calls `BookMVA()` 5 times. `BookMVA` parses XML weight files and reconstructs
decision trees, which is expensive (potentially milliseconds per call).

**Current flow** (per event):
```
cal_numu_bdts_xgboost()
  -> cal_numu_1_bdt()    [creates Reader, BookMVA, EvaluateMVA, destroys Reader]
  -> cal_numu_2_bdt()    [same]
  -> cal_numu_3_bdt()    [same]
  -> cal_cosmict_10_bdt()[same]
  -> main xgboost        [creates Reader, BookMVA, EvaluateMVA, destroys Reader]
```

**Recommendation**: Cache the 5 `TMVA::Reader` objects as `mutable` class
members (or use lazy initialization). Call `BookMVA()` once at construction
time or on first invocation. The `EvaluateMVA()` calls update the bound
`float*` variables and re-score, which is the fast path.

This would also require the per-iteration local variables in `cal_numu_1_bdt`,
`cal_numu_2_bdt`, and `cal_cosmict_10_bdt` to become class members so their
addresses remain stable across calls.

---

## 4. Determinism

**Verdict: Fully addressed in toolkit.**

| Container | Prototype | Toolkit |
|---|---|---|
| Vertex sets | `std::set<ProtoVertex*>` (pointer-ordered) | `IndexedVertexSet` with `VertexIndexCmp` (graph-index ordered) |
| Segment sets | `std::set<ProtoSegment*>` (pointer-ordered) | `IndexedSegmentSet` with `SegmentIndexCmp` (graph-index ordered) |
| Shower sets | `std::set<WCShower*>` (pointer-ordered) | `IndexedShowerSet` with `ShowerIndexCmp` (shower-id ordered) |
| Vertex-to-segment map | `std::map<ProtoVertex*, std::set<ProtoSegment*>>` | `boost::out_edges(vd, graph)` (BGL insertion-order edge list) |
| Segment-to-vertex map | `std::map<ProtoSegment*, std::pair<ProtoVertex*, ProtoVertex*>>` | `boost::source()`/`boost::target()` from edge descriptor |

All pointer-based ordering in the prototype has been replaced with stable
index-based comparators. The BGL `adjacency_list` with `vecS` vertex/edge
storage uses stable numeric indices, ensuring deterministic iteration order
independent of memory layout.

No remaining pointer-based containers were found in any of the reviewed
functions or their callees (`calculate_num_daughter_tracks`,
`seg_at_main_vertex`, `get_last_segment_vertex_long_muon`).

---

## 5. Multi-APA/Face Handling

| Function | APA-sensitive? | Assessment |
|---|---|---|
| `numu_tagger()` | No | Operates on graph topology (segments, vertices, showers). APA-agnostic by design. |
| `count_daughters()` | No | Graph BFS traversal. APA-agnostic. |
| `find_cont_muon_segment()` | No | Uses segment directions and lengths, not raw wire info. APA-agnostic. |
| `cal_numu_bdts_xgboost()` | Indirect | Consumes pre-computed features. APA-correctness depends on upstream filling of `TaggerInfo`. |
| `cluster_fc_check()` | **Yes** | Correctly handles multi-APA: per-point APA lookup via `dv->contained_by(p1)`, per-WPID wire-plane direction maps. |

`cluster_fc_check()` is the only function that directly interacts with
detector geometry. Its multi-APA handling works as follows:

1. `compute_wireplane_params(wpids, dv, ...)` precomputes drift direction and
   U/V/W wire angles for each `WirePlaneId` present in the cluster's grouping.
2. For each candidate boundary point, `dv->contained_by(p1)` determines which
   APA/face the point belongs to.
3. The corresponding wire-plane directions are looked up from the precomputed
   maps before computing projected angles.
4. `fiducial_utils->inside_fiducial_volume()`, `check_signal_processing()`,
   and `check_dead_volume()` are called with the correct per-APA context.

This is a proper generalization of the prototype's single-APA logic.

---

## 6. Summary of Toolkit Improvements Over Prototype

1. **Deterministic ordering**: All pointer-based containers replaced with
   index-based comparators.
2. **`std::clamp` before `acos()`**: Added in `find_cont_muon_segment()`
   (NeutrinoVertexFinder.cxx:963, 976). Prevents NaN from FP rounding.
3. **Zero-magnitude guard**: Added in `find_cont_muon_segment()` (line 958).
4. **Safer vertex position**: `vtx_fit_pt()` falls back to `wcpt().point` if
   fit is invalid, avoiding use of uninitialized data.
5. **Explicit operator precedence**: Parentheses added in all compound
   conditions where C++ precedence could mislead readers.
6. **Dead code removal**: Unused `muon_dir`, `dir_beam`, `dir_drift`,
   `dir_vertical`, `flag_print`, `max_direct_length`, `max_medium_dQ_dx`
   removed.
7. **Configurable BDT weights**: Hardcoded XML paths replaced with
   configurable `m_*_xml` members resolved via `Persist::resolve()`.
8. **Multi-APA support**: `cluster_fc_check()` generalized to handle
   multiple APAs/faces via per-point geometry lookups.

---

## 7. Recommendations

| Priority | Item | Action | Status |
|---|---|---|---|
| **P1** | BUG-1: unclamped `acos()` in `cluster_fc_check()` | Add `std::clamp(-1.0, 1.0)` to 6 `acos()` calls at `Clustering_Util.cxx:143-158`; add zero-magnitude guard for `dir_1` | **Fixed** (commit aeb7980e) |
| **P1** | BUG-2: unclamped `acos()` for `main_angle` | Add `std::clamp` at `Clustering_Util.cxx:171` | **Fixed** (commit aeb7980e) |
| **P2** | EFF-1: TMVA::Reader recreation | Cache 5 Reader objects as class members in `UbooneNumuBDTScorer`; call `BookMVA` once | Open |
| **P3** | Unused variable suppressions | Remove `max_angle`, `max_ratio`, `max_ratio1_length` computations from `find_cont_muon_segment()` (lines 925-926, 993-994, 1000-1001) or document their retention | Open |
