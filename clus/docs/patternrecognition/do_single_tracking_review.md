# do_single_tracking Port Review

**Reviewed:** 2026-04-10  
**Fixes applied (2026-04-10):** §2.1 (comment added), §2.2 (FIXED), §2.3 (FIXED), §2.7 (FIXED), §2.8 (FIXED), §2.10 (FIXED), §2.11 (FIXED), §2.12 (FIXED), §2.13 (FIXED), §2.14 (FIXED), §2.15 (FIXED), §2.4/§5.1 (FIXED — multi-face blob+Steiner union in `form_point_association`), §5.2 (FIXED — option c: debug log on `contained_by()==(-1,-1)`), §3.D1 (FIXED — `vertex_index_map`/`segment_point_index_map` stable-int keys). §2.5, §2.6, §2.16 deferred per user request. §2.9 deferred (complex, requires channel→wire reverse lookup). §3.D2-D10 deferred (no stable ID on `Facade::Cluster`/`Facade::Blob`). S4 parameter parity verified — all defaults match prototype.  
**Entry point:** `TrackFitting::do_single_tracking` — `clus/src/TrackFitting.cxx:7985`  
**Prototype entry:** `WCPPID::PR3DCluster::do_tracking` — `prototype-dev/wire-cell/pid/src/PR3DCluster.cxx:33`  
**Scope:** Functional equivalence, bugs, efficiency, determinism, multi-APA/face correctness.  
**Method:** Three parallel deep-review sub-agents compared all helper functions against the
single-TPC prototype (`WCPPID::PR3DCluster`) and the prototype's `_multi_*` variants.  
All five review goals addressed:
1. Functional identity with prototype.
2. Bugs in toolkit.
3. Algorithmic efficiency.
4. Determinism (pointer-keyed containers).
5. Multi-APA / multi-face correctness.

**Key reference notes:**
- Toolkit blob wire ranges are half-open `[min, max)`. The `<` vs prototype's `<=` is intentional, not a bug.
- Prototype had separate `PR3DCluster_*` (single) and `PR3DCluster_multi_*` (multi-segment) variants. The toolkit merges them.
- `do_rough_path` / `adjust_rough_path` are ported but relocated — `PatternAlgorithms::do_rough_path` in `NeutrinoPatternBase.cxx:78` and `adjust_rough_path` in `TaggerCheckSTM.cxx:572`. Not in `TrackFitting`.
- `update_data_dQ_dx_fit` is ported as `TrackFitting::update_dQ_dx_data()` (`TrackFitting.cxx:5148`).

---

## 1. Function Map (prototype → toolkit)

| Prototype function | Prototype file:line | Toolkit function | Toolkit file:line | Notes |
|---|---|---|---|---|
| `do_tracking` | `PR3DCluster.cxx:33` | `do_single_tracking` | `TrackFitting.cxx:7985` | Entry point |
| `prepare_data` | `PR3DCluster_trajectory_fit.h:1804` | `prepare_data` | `TrackFitting.cxx:690` | Feature gap — see §2.8, §2.9 |
| `organize_wcps_path` | `PR3DCluster_trajectory_fit.h:2079` | `organize_orig_path` | `TrackFitting.cxx:1696` | Renamed |
| `organize_ps_path` | `PR3DCluster_trajectory_fit.h:1994` | `organize_ps_path` | `TrackFitting.cxx:1881` | Ported cleanly |
| `examine_end_ps` | `PR3DCluster_trajectory_fit.h:1933` | `examine_end_ps_vec` | `TrackFitting.cxx:1789` | See §2.7 |
| `form_map` | `PR3DCluster_trajectory_fit.h:1694` | `form_map` | `TrackFitting.cxx:3282` | Structural diff — see §2.4, §2.5 |
| `form_point_association` | `PR3DCluster_trajectory_fit.h:1197` | `form_point_association` | `TrackFitting.cxx:1966` | Significant drift — see §2.4–2.6 |
| `examine_point_association` | `PR3DCluster_trajectory_fit.h:804` | `examine_point_association` | `TrackFitting.cxx:2665` | Intentional deviation — see §2.17 |
| `trajectory_fit` | `PR3DCluster_trajectory_fit.h:39` | `trajectory_fit` | `TrackFitting.cxx:3973` | Frame-mixing bug — see §2.3 |
| `skip_trajectory_point` | `PR3DCluster_trajectory_fit.h:479` | `skip_trajectory_point` | `TrackFitting.cxx:4711` | Two bugs — see §2.1, §2.2 |
| `dQ_dx_fit` | `PR3DCluster_dQ_dx_fit.h:368` | `dQ_dx_fit` | `TrackFitting.cxx:6686` | See §2.11, §2.13–2.15 |
| `dQ_dx_fill` | `PR3DCluster_dQ_dx_fit.h:266` | `dQ_dx_fill` | `TrackFitting.cxx:5482` | Minor gap — see §4.3 |
| `cal_gaus_integral_seg` | `PR3DCluster_dQ_dx_fit.h:143` | `cal_gaus_integral_seg` | `TrackFitting.cxx:5131` | Div-by-zero — see §2.12 |
| `cal_gaus_integral` | `PR3DCluster_dQ_dx_fit.h:160` | `cal_gaus_integral` | `TrackFitting.cxx:4984` | See §2.14 |
| `cal_compact_matrix` | `PR3DCluster_dQ_dx_fit.h:3` | `calculate_compact_matrix` | `TrackFitting.cxx:5329` | Ported, added zero-guard |
| `cal_compact_matrix_multi` | `PR3DCluster_multi_dQ_dx_fit.h:967` | `calculate_compact_matrix_multi` | `TrackFitting.cxx:5188` | Determinism bug — see §3.D1 |
| `update_data_dQ_dx_fit` | `PR3DCluster_dQ_dx_fit.h:1047` | `update_dQ_dx_data` | `TrackFitting.cxx:5148` | Renamed + inverted loop |
| `do_rough_path` | `PR3DCluster_crawl.h:3` | `PatternAlgorithms::do_rough_path` | `NeutrinoPatternBase.cxx:78` | Relocated |
| `adjust_rough_path` | `PR3DCluster_crawl.h:22` | (local in TaggerCheckSTM) | `TaggerCheckSTM.cxx:572` | Relocated |
| `collect_charge_trajectory` | `PR3DCluster.cxx:1128` | **MISSING** | — | See §2.16 |

---

## 2. Bugs (S1 — Correctness)

### 2.1 `skip_trajectory_point` uses position index `i` instead of form_map count index — NOT A BUG (clarifying comment added)

**File:** `TrackFitting.cxx:4944`  
**Prototype:** `PR3DCluster_trajectory_fit.h:479` — takes separate `index` and `i` arguments; the `index` argument is the *count* assigned during `form_map`, while `i` is the position in `fine_tracking_path`.

After `trajectory_fit` removes or reorders points, `i` in `skip_trajectory_point` no longer aligns with the `m_3d_to_2d` key that was assigned during `form_map`. The dead-plane detection (`m_3d_to_2d.find(i)` at `:4944`) silently looks up the wrong point and produces the wrong dead-plane count.

**Fix:** Thread the `form_map` count index through `skip_trajectory_point` as a separate parameter, or pre-resolve `dead_plane_count` during `form_map` and store it alongside `pss_vec`.

---

### 2.2 `skip_trajectory_point` uses single `apa_face` for both points being compared — FIXED

**File:** `TrackFitting.cxx:4713-4714`

The function receives a single `apa_face` and projects both the current fit point `p` AND a neighbouring fit point under that one face. If the two points lie in different (apa, face) regions (possible for a track crossing an APA boundary), their wire/time arithmetic is mixed under the wrong geometry frame.

```cpp
// Current — both points projected under the same face:
int t1 = std::round((offset_t + slope_x * p_raw.x())/cur_ntime_ticks) * cur_ntime_ticks;
// ...
int t2 = std::round((offset_t + slope_x * prev_p_raw.x())/cur_ntime_ticks) * cur_ntime_ticks;
```

**Fix:** Resolve per-point `(apa, face)` via `m_dv->contained_by(p)` and `m_dv->contained_by(prev_point)` internally, look up `wpid_offsets`/`wpid_slopes` per-point, and project each point under its own frame. Skip the comparison when the two points live in different faces (distance calculation becomes ill-posed across APA boundaries).

---

### 2.3 `trajectory_fit` Gaussian weighting mixes coordinate frames — FIXED

**Files:** `TrackFitting.cxx:4055-4083`, `TrackFitting.cxx:3706-3718` (multi variant)

In the method-2 Gaussian division loop, the 3D point `p` is backward-transformed under its OWN (apa, face) via `m_dv->contained_by(...)` at line `4055`:

```cpp
auto test_wpid = m_dv->contained_by(pss_vec[idx].first);
auto p_raw = transform->backward(pss_vec[idx].first, cluster_t0, test_wpid.face(), test_wpid.apa());
```

But `central_t` and `central_ch` are then computed using `slope_x` and `offset_t` from the OUTER loop's `apa_face` (the pixel's (apa, face)):

```cpp
central_t = slope_x * p_raw.x() + offset_t;   // offset_t from outer (apa,face)!
```

For a 3D trajectory point near an APA boundary, when processing pixels that belong to a different face than the 3D point's own face, this multiplies one face's drift-corrected x by another face's slope and adds another face's offset — an incorrect mixed-frame projection.

**Fix:** Inside the per-2D-pixel loop, re-transform the 3D point under the PIXEL'S (apa, face) (from the outer `apa_face` loop variable). Use that face's `slope_x`, `offset_t`, `slope_y{u,v,w}`, etc. to compute `central_t` and `central_ch{u,v,w}`.

---

### 2.4 `form_point_association` blob filter is single-face only — FIXED

**File:** `TrackFitting.cxx:2038` (old filter removed)

The face filter `if (blob_wpid.apa()==apa && blob_wpid.face()==face)` silently dropped all blobs whose face differed from the closest cluster point's face, preventing charge collection from neighbouring faces for boundary-crossing tracks.

**Fix applied:** Removed the face filter. Nearby blobs from all faces are now collected into `nearby_blobs_set`. They are grouped by `(apa, face)` into `blobs_by_face` (a `std::map<pair<int,int>, vector<const Facade::Blob*>>`). Each face group is processed separately: the trajectory point `p` is re-projected into that face using `transform->backward(p, cluster_t0, face2, apa2)` and `convert_3Dpoint_time_ch`, and geometry (`angle_u/v/w`, `pitch_u/v/w`, `time_tick_width`) is looked up from `wpid_params`/`wpid_geoms` using the blob's own wpid. Output `Coord2D` records carry the blob's own `(apa2, face2)`.

Similarly, the Steiner branch now uses the Steiner closest-point's own `(st_apa, st_face)` from `closest_point_wpid` for all projections and geometry lookups — previously the outer-scope `apa`/`face` (from the primary blob face) were used unconditionally.

---

### 2.5 `form_point_association` Steiner branch collapses blob wire ranges to `[wire, wire+1)`

**File:** `TrackFitting.cxx:2287-2304`

For Steiner-path vertices, the toolkit inserts only a single wire span `[wire, wire+1)` into `map_time_wires`:

```cpp
int umin = vertex_wire_u, umax = vertex_wire_u+1;  // TrackFitting.cxx:2287-2289
```

The prototype (`PR3DCluster_trajectory_fit.h:1484-1660`) has two separate contribution sources for Steiner vertices: (i) the full wire span of each attached blob (`uwires.front()->index()` … `uwires.back()->index()`), and (ii) the raw-point `±1` neighbourhood formula. The toolkit only replicates the raw-point stencil and loses all blob-span information for Steiner vertices.

**Fix:** Add a blob-iteration branch for Steiner vertices that mirrors the non-Steiner blob loop (`:2060-2090`), using `blob->u_wire_index_min()` / `u_wire_index_max()` to fill `map_time_wires` with the full blob span.

---

### 2.6 `form_point_association` fallback projection has a narrower wire diamond than the prototype

**File:** `TrackFitting.cxx:2481`

```cpp
int wire_cut = time_cut / cur_ntime_ticks;   // wire half-width ÷ tick width (≈4)
```

This makes the wire half-span roughly 4× smaller than the prototype's fallback (`PR3DCluster_trajectory_fit.h:1669-1688`), which uses a full `[-time_cut, time_cut]` × `[-time_cut, time_cut]` L1 diamond where both axes have the same half-span (`time_cut` ≈ 5 time-slices).

**Fix:** Decouple wire half-width from time half-width. The wire loop should run from `fb_wire - time_cut` to `fb_wire + time_cut` (same radius as the prototype), not `fb_wire ± wire_cut`.

---

### 2.7 `examine_end_ps_vec` silently pops face-invalid points and may return an outside-detector point — FIXED

**File:** `TrackFitting.cxx:1804-1834`

When `contained_by(temp_start)` returns `face=-1`, the loop pops the point without checking `is_good_point`. If every point in `ps_list` has `face=-1`, the loop will exhaust the list, and the final `ps_list.push_front(temp_start)` at line `:1834` re-inserts the last-popped invalid point into the output. The caller then receives a point that is outside any face.

**Fix:** Distinguish between "face invalid" (point outside detector volume) and "bad charge" (good geometry, insufficient hits). For face-invalid points either: (a) break immediately and return the accumulated good points so far, or (b) accept the point unconditionally (like the prototype, which had no face concept) and leave the is_good_point test for in-face points only.

---

### 2.8 `prepare_data` dead-blob loop ignores multi-slice blobs — FIXED

**File:** `TrackFitting.cxx:811-847`

The dead-blob loop uses only `blob->slice_index_min()` for the time coordinate and never steps through the full slice range:

```cpp
int time_slice = blob->slice_index_min();   // only min — max is never used
// ... one entry per blob regardless of how many slices it spans
```

A dead blob that spans slices `[min, max)` should contribute one charge entry per slice. As written, only the first slice is covered.

**Fix:** Add an inner loop `for (int time_slice = blob->slice_index_min(); time_slice < blob->slice_index_max(); ++time_slice)` and move the per-wire entries inside it. Re-compute `charge` and `charge_err` per slice if the dead charge is distributed (simplest: uniform `blob_charge / (num_wires * (max-min))`).

---

### 2.9 `prepare_data` missing dead-channel neighbour error inflation — DEFERRED

**Prototype:** `PR3DCluster_trajectory_fit.h:1815-1857`

The prototype inflates `charge_err` for live pixels that are adjacent to dead channels:
- Induction planes: `charge_err *= 5`
- Collection plane: `charge_err *= 2.5`

This gives less weight to measurements near dead regions during the fit. The toolkit's `prepare_data` contains no equivalent logic.

**Fix:** After populating `m_charge_data`, iterate the dead-blob pixel set (known from the blob loop at `:778-847`). For each dead `(time, wire)` entry, check `(time ± 1, wire)`, `(time, wire ± 1)` in `m_charge_data`; if found and non-dead, inflate their `charge_err` by the appropriate factor. Alternatively, apply the inflation at row-insertion time inside `trajectory_fit` / `dQ_dx_fit` using `m_charge_data[key].flag == 0` as the dead indicator.

---

### 2.10 `dQ_dx_multi_fit` silent index-0 collapse via `operator[]` — FIXED

**File:** `TrackFitting.cxx:5889`

```cpp
pt_idx[i] = segment_point_index_map[std::make_pair(segment, i)];  // operator[]
```

`std::map::operator[]` inserts a zero-initialised value on a missing key. If `i` exceeds the range populated for `segment` in the earlier setup loop, every out-of-range point silently gets `pt_idx[i] = 0`, collapsing multiple trajectory indices to row 0 of the Eigen matrix and corrupting the fit.

**Fix:** Use `segment_point_index_map.at(std::make_pair(segment, i))` (throws on missing key — exposes the bug) or `find` + assert. Same applies at `:6499` and `:6673`.

---

### 2.11 `dQ_dx_multi_fit` reduced-χ² not guarded against zero prediction — FIXED

**File:** `TrackFitting.cxx:6631-6642`

`dQ_dx_fit` (single-track) correctly guards at `:7415`:

```cpp
if (pred_data_u_2D(it.row()) > 0)  // guard present
```

The multi variant at `:6631-6642` has no such guard. If `pred_data_*(row) == 0` (degenerate column), the division produces `Inf` or `NaN` in `reduced_chi2`.

**Fix:** Mirror the guard from `dQ_dx_fit:7415` into the same positions in `dQ_dx_multi_fit`.

---

### 2.12 `cal_gaus_integral_seg` divide-by-zero on degenerate zero-length segment — FIXED

**File:** `TrackFitting.cxx:5142`

```cpp
result /= result1;   // result1 can be zero
```

`result1` is the sum of distance-based weights. For a degenerate zero-length segment (current == previous == next fit point), all `weights[i]` become zero, `result1 = 0`, and division produces `NaN`, which propagates into `RU.insert(row, col)` of the Eigen matrix.

**Fix:**
```cpp
if (result1 > 0) result /= result1;
else result = 0.0;
```

---

### 2.13 `dQ_dx_fit` crashes on empty cluster — FIXED

**File:** `TrackFitting.cxx:6946`

```cpp
auto first_blob = cluster->children()[0];   // crash if children() is empty
```

**Fix:** Early-return if `cluster->children().empty()`. Since `cur_ntime_ticks` is used throughout the function, a `cur_ntime_ticks = 1` fallback should be assigned in that case or propagated differently.

---

### 2.14 `cal_gaus_integral` time-integration width uses a single blob's slice count — FIXED (empty-cluster guard added; per-blob cur_ntime_ticks deferred)

**File:** `TrackFitting.cxx:6946` (feeds into `:4984`)

`cur_ntime_ticks` is read once from `cluster->children()[0]->slice_index_max() - slice_index_min()`. If other blobs in the cluster have a different slice width (non-uniform rebinning), `cal_gaus_integral` integrates them over the wrong time window.

Additionally, within `cal_gaus_integral` at `:5012`, the wire dimension hardcodes `0.5` for the half-bin width while the time dimension uses `0.5 * cur_ntime_ticks`. This asymmetry is correct (wires are always 1-pitch wide) but should be commented.

**Fix:** Compute `cur_ntime_ticks` per-blob at the call site inside the blob loop in `dQ_dx_fit`, rather than once for the whole cluster.

---

### 2.15 Missing `wpid_offsets.find(...) == end()` guards in hot paths — FIXED

**Files:** `TrackFitting.cxx:4021, 4184, 4248, 4316, 4524`

`wpid_offsets[wpid]` and `wpid_slopes[wpid]` are dereferenced without checking the return of `find()`. The multi variant already guards at `:3779-3785`:

```cpp
if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) continue;
```

The single-track `trajectory_fit` and related loops do not have these guards. A `contained_by(...)` returning a wpid not in `wpid_offsets` (possible for boundary points) will segfault.

**Fix:** Add `if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) continue;` at each of the listed sites.

---

### 2.16 `collect_charge_trajectory` not ported

**Prototype:** `PR3DCluster.cxx:1128`

The prototype's `collect_charge_trajectory` sweeps a tube around the fine tracking path and gathers all live 2D charges within a `dis_cut`-radius into `collected_charge_map`. The toolkit has no equivalent. The `prepare_data` / `fill_global_rb_map` / `update_dQ_dx_data` pipeline is a different approach (blob-based, not tube-based).

**Assessment needed:** Determine whether any downstream consumer (NeutrinoID energy reco, kinematic reconstruction) reads from `collected_charge_map`. If so, this is a missing feature that must be ported. If the `global_rb_map`-based approach fully replaces it, document that explicitly.

---

### 2.17 `examine_point_association` — intentional but silent deviations from prototype

**File:** `TrackFitting.cxx:2796, 2836, 2876, 2917, 2958, 2999`

Two intentional changes that are improvements but break prototype-parity:

**(a) Time RMS normalisation:** The toolkit divides time residuals by `cur_ntime_ticks` before computing the outlier RMS:
```cpp
pow((time - ave_pos.second)/cur_ntime_ticks, 2)   // TrackFitting.cxx:2796
```
The prototype does not normalise (`pow(it->second - ave_pos.second, 2)`). This makes the RMS dimension-consistent (wire units) but changes the effective threshold `0.75*rms`.

**(b) Outlier wire-index fix:** The prototype incorrectly compares V-plane or U-plane outliers against `temp_results.at(3)` (which is the W wire index) due to a copy-paste error (prototype `PR3DCluster_trajectory_fit.h:944, 979, 1011, 1043, 1077, 1108, 1142, 1172`). The toolkit uses the correct per-plane expected wire.

**Action:** Document both deviations explicitly here (done) and verify with an end-to-end comparison plot that the change in (a) does not degrade track quality relative to the prototype.

---

## 3. Determinism Issues (S2)

The following pointer-keyed containers leak non-deterministic iteration order into observable state (Eigen row/column assignments, merged charge values). Listed in priority order.

### 3.D1 `vertex_index_map` and `segment_point_index_map` — CRITICAL — FIXED

**File:** `TrackFitting.cxx:5659-5660, 6217, 6508, 6648`

```cpp
std::map<std::shared_ptr<PR::Vertex>, int> vertex_index_map;     // 5659
std::map<std::pair<std::shared_ptr<PR::Segment>, int>, int> segment_point_index_map;  // 5660
```

`std::map` with `shared_ptr` keys orders by raw pointer address (default `std::less<shared_ptr<T>>` compares the stored pointer). The iteration at `:6217, 6508, 6648` over `vertex_index_map` drives:

```cpp
connected_vec[vertex_idx].push_back(fits[1].index);    // 6521
connected_vec[vertex_idx].push_back(fits[fits.size()-2].index);  // 6523
```

`connected_vec[i]` is a `std::vector`; the ORDER of its elements (set by `push_back`) determines what `calculate_compact_matrix_multi` reads at `connected_vec[i][0]` vs `connected_vec[i][1]`. Different runs with different ASLR / allocator state will construct different regularizer topologies and produce different fit outputs.

**Fix:** Key both maps by stable integer IDs:
- `vertex_index_map`: replace `shared_ptr<PR::Vertex>` key with `int` (the vertex's `fit().index`). Build the map via `vertex->fit().index → vertex`.
- `segment_point_index_map`: replace `shared_ptr<PR::Segment>` key with `int` (the segment's `get_graph_index()`). Build the map via `{graph_index, i} → matrix_index`.

---

### 3.D2 `m_clusters` and `m_loaded_clusters` — pointer-ordered sets — DEFERRED (no stable ID on Facade::Cluster)

**File:** `TrackFitting.h:569-570`

```cpp
std::set<Facade::Cluster*> m_clusters;
std::set<Facade::Cluster*> m_loaded_clusters;
```

Iterated in `prepare_data` and related functions. The iteration order drives the order in which clusters' charge data is merged into `m_charge_data`. When two clusters share a `(time, wire, plane, apa, face)` pixel, the merge at `:838-840` (charge accumulation for `flag==0`) depends on which cluster is processed first.

**Fix:** Replace with `std::set<Facade::Cluster*, ClusterIdLess>` where `ClusterIdLess` compares by `cluster->ident()`. Or use a sorted `std::vector<Facade::Cluster*>` populated once and kept sorted by `ident()`.

---

### 3.D3 `m_cluster_charge_data` — pointer-keyed unordered_map

**File:** `TrackFitting.h:582`

```cpp
std::unordered_map<Facade::Cluster*, std::unordered_map<CoordReadout,...>> m_cluster_charge_data;
```

Iterated at `TrackFitting.cxx:5589-5595`. Each cluster's charge data is pushed into `m_charge_data` in pointer-hash order, which depends on ASLR.

**Fix:** Convert to `std::map<int, std::map<CoordReadout,...>>` keyed by `cluster->ident()`. The inner `CoordReadout` map is already value-keyed and deterministic.

---

### 3.D4 `m_cluster_edges` — pointer-keyed unordered_map

**File:** `TrackFitting.h:575`

```cpp
std::unordered_map<Facade::Cluster*, std::vector<PR::edge_descriptor>> m_cluster_edges;
```

**Fix:** Key by `cluster->ident()` (int).

---

### 3.D5 `m_blobs` — pointer-ordered set

**File:** `TrackFitting.h:584`

```cpp
std::set<Facade::Blob*> m_blobs;
```

**Fix:** Use `std::set<Facade::Blob*, BlobIdLess>` with `BlobIdLess` comparing by `blob->ident()`.

---

### 3.D6 `global_rb_map` — unordered_set of Blob pointers as inner value

**File:** `TrackFitting.h:644`

```cpp
std::unordered_map<CoordReadout, std::unordered_set<Facade::Blob*>, CoordReadoutHash> global_rb_map;
```

The inner `unordered_set<Blob*>` is iterated at `:5166-5171` and `:1040`. Currently, the `is_shared` boolean and the `track_blobs_set` membership test are order-invariant (boolean short-circuit and `find()` respectively). However, this is fragile — any future change to the loop body can break determinism silently.

**Fix:** Replace inner `unordered_set<Facade::Blob*>` with `std::set<Facade::Blob*, BlobIdLess>`.

---

### 3.D7 `new_clusters` in `prepare_data` — pointer-ordered local set

**File:** `TrackFitting.cxx:696`

```cpp
std::set<Facade::Cluster*> new_clusters;
```

Drives the `m_charge_data` merge loop. See 3.D2 for fix pattern.

---

### 3.D8 `nearby_blobs_set` in `form_point_association`

**File:** `TrackFitting.cxx:2033`

```cpp
std::unordered_set<const Facade::Blob*> nearby_blobs_set;
```

Downstream consumers (`std::max` for time-slice bounds, `std::set<Coord2D>` for pixel insertion) are currently order-invariant, so the immediate output is deterministic. However, if §2.4's fix removes the face filter and the loop emits pixels conditionally, order-dependency will be introduced.

**Fix (proactive):** Change to `std::set<const Facade::Blob*, BlobIdLess>`.

---

### 3.D9 `track_blobs_set` in `update_dQ_dx_data`

**File:** `TrackFitting.cxx:5150`

```cpp
std::set<Facade::Blob*> track_blobs_set;   // pointer-ordered, used for find()
```

Only used for `find()` membership tests — the result is order-invariant. Still flagged per the policy of eliminating all pointer-keyed ordered containers.

**Fix:** `std::set<Facade::Blob*, BlobIdLess>`.

---

### 3.D10 `m_charge_data` / `m_orig_charge_data` — unordered_map with value key (weaker)

**File:** `TrackFitting.h:637-638`

```cpp
std::unordered_map<CoordReadout, ChargeMeasurement, CoordReadoutHash> m_charge_data, m_orig_charge_data;
```

`CoordReadout` is a value type with `operator<` defined; its hash in `CoordReadoutHash` should be stable. Iteration order is still implementation-defined (not sorted), but because all downstream consumers of `m_charge_data` either use `find()` (order-invariant) or copy keys into a sorted `std::map<CoordReadout,...>` before iterating, the current state is acceptable. Document this dependency.

**Recommendation:** Leave as-is for now but add a comment; if any new code iterates `m_charge_data` directly, it must first sort the keys.

---

## 4. Efficiency Wins (S3)

### 4.1 Hoist `m_pcts->pc_transform(...)` out of hot loops

**Files:** `TrackFitting.cxx:4109` (inside `trajectory_fit` per-point loop), `TrackFitting.cxx:7032, 7063` (inside `dQ_dx_fit`), and `cal_gaus_integral_seg` per-sample.

The call chain `m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()))` is repeated for every 3D trajectory point inside the fitting loops. The transform is constant for the lifetime of a given cluster.

**Fix:** Cache `auto transform = ...` per cluster BEFORE the per-point loop. For segments spanning one cluster (common case), this replaces N map lookups with 1.

---

### 4.2 Hoist `wpid_offsets.find(wpid)` and `wpid_slopes.find(wpid)` out of per-pixel loops

**Files:** `TrackFitting.cxx:4184, 4249, 4317, 1994-1995`

Per 2D pixel, `find()` on `wpid_offsets` and `wpid_slopes` (both `std::map<WirePlaneId,...>`) is repeated. Many pixels for the same trajectory point share one `(apa, face)` and therefore the same wpid.

**Fix:** Group the 2D pixels by `(apa, face)` before the inner loop, look up offsets/slopes once per group.

---

### 4.3 Unify `dQ_dx_fit` with `dQ_dx_multi_fit` fast-path for wire-indexed lookup

**Files:** `TrackFitting.cxx:7157-7244` (single) vs `TrackFitting.cxx:6104-6181` (multi)

`dQ_dx_multi_fit` already uses a `wire_idx_{U,V,W}` map to iterate only wires within a search radius (an O(log N + k) indexed scan). `dQ_dx_fit` retains the prototype's O(N_points × N_measurements) nested iteration over the full `map_U_charge_2D`. For large clusters this is a significant bottleneck.

**Fix:** Port the `wire_idx_*` indexed fast-path from `dQ_dx_multi_fit` into `dQ_dx_fit`. Both functions share the same `m_charge_data` source; the infrastructure already exists.

---

### 4.4 `form_point_association` — cache backward-transform per vertex

**File:** `TrackFitting.cxx:2277-2280`

For Steiner-path vertices, `transform->backward(vertex, ...)` plus three `convert_3Dpoint_time_ch(...)` are called per vertex per plane. These are three independent linear projections that could be computed once per vertex and reused.

**Fix:** Compute and cache `(u_raw, v_raw, w_raw, t_raw)` once per vertex before the plane loop.

---

### 4.5 `calculate_compact_matrix(_multi)` — replace `std::map<pair<int,int>, double>` with dense structure

**File:** `TrackFitting.cxx:5341-5372` (single), `TrackFitting.cxx:5188-5326` (multi)

`map_pair_values` is a `std::map<std::pair<int,int>, double>` used as a dense overlap-weight table. For typical trajectory sizes the keys are dense integer pairs in `[0, N) × [0, N)`. An `std::vector<std::vector<double>>` (or flat `std::vector<double>` with row-major indexing) indexed directly would be ~10× faster on random access.

---

### 4.6 `ptss = saved_pts` copy → move — FIXED

**File:** `TrackFitting.cxx:3358`

```cpp
ptss = saved_pts;   // full vector copy
```

**Fix:** `std::swap(ptss, saved_pts);` avoids the allocation.

---

### 4.7 `m_3d_to_2d[i]` silent insertion → `m_3d_to_2d.at(i)`

**File:** `TrackFitting.cxx:4105`

`operator[]` on `std::map<int, Point3DInfo>` silently inserts a zero-initialised entry when `i` is missing, masking logic bugs. Switch to `at()` (throws on miss, making bugs visible). Related to bug §2.1.

---

### 4.8 `cal_gaus_integral` erf caching

**File:** `TrackFitting.cxx:4984-5129`

`cal_gaus_integral` is called inside the innermost dQ/dx fit loop (once per 3D trajectory point per neighbouring 2D pixel, with 10 sample points each). Each call computes 4–16 `erf()` evaluations. Since `t_sigma`, `w_sigma` change slowly along the trajectory, tabulating `erf((bin ± 0.5 * cur_ntime_ticks - t_center) / (sqrt(2)*t_sigma))` per trajectory point (over a small wire/time window) could reduce erf evaluations by ~10×.

**Note:** This is a significant algorithmic change; prototype feasibility should be assessed first on a performance profile.

---

## 5. Multi-APA / Multi-Face Correctness (S5)

### 5.1 Charge aggregation should union across all (apa, face) at APA boundaries — PARTIALLY FIXED

Multiple sites filter 2D pixel associations to a single `(apa, face)`. For tracks that cross an APA boundary, the correct behaviour is to collect charge from ALL faces whose blobs the trajectory point intersects.

Sites addressed:
- `form_point_association:2038` — **FIXED** (see §2.4): face filter removed; blobs now grouped by face and each face projected independently.
- Steiner branch in `form_point_association` — **FIXED**: now uses `closest_point_wpid.apa()/face()` instead of outer-scope `apa`/`face`.

Sites still open (deferred):
- `dQ_dx_fit` — `precomp_UT[{apa,face}]` keyed by single face; cross-face charge invisible (deferred).
- `dQ_dx_multi_fit:5933` — points with `contained_by() == (-1,-1)` now log a debug message and continue (§5.2 option c applied).

### 5.2 Policy when `contained_by(p)` returns `(-1, -1)` — FIXED (option c)

At `dQ_dx_multi_fit` and `dQ_dx_fit`, when `contained_by(p)` returns `(-1,-1)`:

**Fix applied (option c):** Both sites now emit a `SPDLOG_LOGGER_DEBUG` message before continuing:
```cpp
if (test_wpid.apa() == -1 || test_wpid.face() == -1) {
    SPDLOG_LOGGER_DEBUG(s_log, "dQ_dx_multi_fit: trajectory point at ({:.2f},{:.2f},{:.2f}) "
        "has face=-1; regulariser-only constraint for this row", ...);
    continue;
}
```
The Eigen row remains regulariser-only for boundary points; the debug log allows monitoring of how often this occurs. Options (a) snap-to-nearest-face and (b) copy-previous-solution remain available for future improvement if tests show significant impact.

### 5.3 `skip_trajectory_point` is single-face per call

See §2.2. For multi-APA tracks, this function must resolve (apa, face) per point internally.

### 5.4 `examine_end_ps_vec` — correctly handles multi-face

**File:** `TrackFitting.cxx:1789-1878`

`examine_end_ps_vec` correctly calls `m_dv->contained_by(temp_start)` per test point and dispatches `m_grouping->is_good_point(raw, apa, face, ...)` with the per-point `(apa, face)`. This is the correct pattern — other functions should follow this model.

### 5.5 `examine_point_association` — fallback pixels may land in wrong face

**File:** `TrackFitting.cxx:2892-3020`

In the single-dead-plane fallback branches, the toolkit inserts pixels via:
```cpp
insert(Coord2D(apa, face, time_u, wire_u, channel_u, kUlayer))
```
using `(apa, face)` from `m_dv->contained_by(p)` at `:2676`. If the pre-existing associated pixels in `temp_2dut` were formed under a different `(apa, face)` (e.g., from `form_point_association` under the closest-cluster-point's face), the fallback pixels could be in a different face from the rest of the association. Downstream code that assumes all pixels for a given trajectory point share one (apa, face) will see mixed data.

**Fix:** Once §2.4's face-filter fix is applied, the pixel set will already be properly multi-face; the fallback should either be disabled or keyed to the dominant face in the existing pixel set.

---

## 6. Parameter Parity Checklist (S4)

Verify that `TrackFittingParams` defaults match the prototype's hardcoded values. Mismatches silently change fit results.

| `m_params` field | Prototype value | Prototype source | Verified? |
|---|---|---|---|
| `lambda` | `0.0005` | `PR3DCluster_dQ_dx_fit.h:933` | MATCH |
| `overlap_th` | `0.5` (i.e., `>1` ≡ `> 2×0.5`) | `:892-914` | MATCH |
| `default_dQ_dx` | `5000` | `PR3DCluster_dQ_dx_fit.h:358` | MATCH |
| `share_charge_err` | `8000` | `PR3DCluster_dQ_dx_fit.h:1071` | MATCH |
| `min_drift_time` | `≈50 µs` (50 µs floor) | `PR3DCluster_dQ_dx_fit.h:608` | MATCH |
| `default_charge_th` | `100` | `PR3DCluster_trajectory_fit.h:202` | MATCH |
| `default_charge_err` | `1000` | `PR3DCluster_trajectory_fit.h:203` | MATCH |
| `skip_ratio_cut` | `0.97` | `PR3DCluster_trajectory_fit.h:713` | MATCH |
| `skip_ratio_1_cut` | `0.75` | `:723` | MATCH |
| `skip_default_ratio_1` | `0.25` | `:791` | MATCH |
| `skip_angle_cut_1` | `160°` | `:745` | MATCH |
| `skip_angle_cut_2` | `90°` | `:781` | MATCH |
| `skip_angle_cut_3` | `45°` | `:733` | MATCH |
| `skip_dis_cut` | `0.5 cm` | `:786` | MATCH |
| `area_ratio1` | `1.8 mm·c` | prototype area-smoothing loop | MATCH |
| `area_ratio2` | `1.7` | same | MATCH |
| `dead_ind_weight`, `dead_col_weight`, `close_ind_weight`, `close_col_weight` | hardcoded | `PR3DCluster_dQ_dx_fit.h` multiple | MATCH (0.3, 0.9, 0.15, 0.45) |

---

## 7. Execution Order

The following order minimises rework: correctness first, then determinism, then efficiency.

| Phase | Work | Severity |
|---|---|---|
| 1 (this doc) | Audit document — team review | — |
| 2 | Verify parameter parity table (§6) | S4 |
| 3 | §2.1–2.3: `skip_trajectory_point` index bug, face-per-point, frame mixing in `trajectory_fit` | S1 critical |
| 4 | §2.4–2.7: `form_point_association` face filter, Steiner wire span, fallback diamond, `examine_end_ps_vec` | S1 multi-APA |
| 5 | §3.D1–3.D9: pointer-keyed container elimination (especially `vertex_index_map`) | S2 determinism |
| 6 | §2.8–2.9: `prepare_data` multi-slice dead blobs and error inflation | S1 feature |
| 7 | §2.10–2.15: dQ/dx defensive guards (zero prediction, empty cluster, NaN, wpid guards) | S1 crashes |
| 8 | §2.16: assess `collect_charge_trajectory` | S1/feature |
| 9 | §4.1–4.8: efficiency wins | S3 |

Each phase is a separate commit. Each phase must build cleanly before moving to the next.

---

## 8. Verification

1. **Determinism sweep.** Run the same event three times with ASLR on. The serialized output (`fine_tracking_path`, `dQ`, `dx`, `reduced_chi2`) must match bit-for-bit. Record a baseline BEFORE Phase 5 and confirm the same baseline AFTER.

2. **Single-face regression.** Pick a single-segment, single-face track. Compare `fine_tracking_path` to the prototype `do_tracking` output within 0.1 mm and `dQ/dx` within 1 %. Document all intentional deviations (§2.17).

3. **Multi-APA event.** Pick or construct an event with a track crossing an APA boundary. After Phase 4: verify `form_point_association` accumulates 2D pixels from both faces for boundary points; verify `reduced_chi2` is finite everywhere.

4. **Boundary-point policy.** After §5.2 decision is implemented, confirm that the fraction of `contained_by() == (-1,-1)` points across a representative sample is below an agreed threshold (e.g. < 1% of trajectory points).

5. **Regression suite.** Run `ctest -R clus` and the `uboone-mabc` jsonnet pipeline (if available) after each phase.

6. **Parameter parity.** After Phase 2, diff `m_params` defaults against the prototype hardcoded values using the table in §6.
