# Examination: UbooneMagnifyTrackingVisitor

> **Bug Fix Status**: Bug 8 (dtheta out-of-bounds crash) fixed. Bugs 1-7, 9-10 deferred — low severity or design decisions needing domain confirmation.
>
> **Efficiency Fix Status**: Efficiency 1 (redundant graph iteration) fixed — builds per-cluster edge/vertex maps in single O(E+V) pass. Efficiencies 2, 4-5 deferred — low impact.

Files examined:
- `root/src/UbooneMagnifyTrackingVisitor.cxx`
- `root/inc/WireCellRoot/UbooneMagnifyTrackingVisitor.h`
- `root/apps/wire-cell-uboone-magnify-tracking-convert.cxx`

---

## 1. Potential Bugs

### Bug 1: `reco_cluster_id` field set but never written to TTree
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 357, 397
- **Description:** The `WCPointTree` struct contains a `reco_cluster_id` field (header line 40) which is assigned at lines 357 and 397 (`point_tree.reco_cluster_id = cluster->get_cluster_id()`), but no TTree branch is ever created for this field. The `cluster_id` branch at line 299 uses `reco_mother_cluster_id` instead. This means per-point cluster identity information (which cluster a vertex or segment actually belongs to) is computed but silently discarded.
- **Severity:** Low
- **Why:** The data is still available via `reco_ndf` (which is set to cluster_id) and `reco_proto_cluster_id`, so downstream consumers can reconstruct this. However, it indicates either dead code or a missing branch.

### Bug 2: `reco_ndf` abused to carry cluster_id, type mismatch with branch spec
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 289, 360, 398
- **Description:** The branch `ndf` is declared as `"ndf/D"` (double), but it is filled with `cluster->get_cluster_id()` (an integer). The field `reco_ndf` in WCPointTree is `double`. This is an intentional encoding for the downstream app (`wire-cell-uboone-magnify-tracking-convert.cxx`, line 331: `std::round(ndf)`) which uses `ndf` as a cluster grouping key via `std::round()`. While it works due to the `std::round`, overloading a semantically named field (`ndf` = number of degrees of freedom) with cluster_id is confusing and fragile. If cluster IDs ever exceed 2^53, precision is lost.
- **Severity:** Low (functional but misleading)
- **Why:** The downstream app relies on `std::round(ndf)` to recover the integer cluster_id. This contract is implicit and undocumented.

### Bug 3: `real_cluster_id` and `sub_cluster_id` branches share the same address
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 302-303
- **Description:** Both branches point to `&point_tree.reco_proto_cluster_id`. They will always contain identical values. The comment says "Keep compatibility with legacy format" but the downstream app only reads `sub_cluster_id` (line 239 of the convert app). If any downstream consumer expects these to differ, it would get wrong results.
- **Severity:** Low
- **Why:** Documented as intentional for legacy compatibility. Not a bug per se, but worth noting.

### Bug 4: `mother_cluster_id` may remain -1 if no cluster has `main_cluster` flag
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 336-342
- **Description:** The code searches `all_clusters` for one with the `main_cluster` flag and uses that as `mother_cluster_id`. If no cluster has this flag set, `mother_cluster_id` stays -1, and all output entries will have `cluster_id = -1`. This silently produces potentially misleading output.
- **Severity:** Medium
- **Why:** Downstream code may filter or group by `cluster_id`. A value of -1 could either be misinterpreted or cause silent data loss. At minimum a warning should be logged.

### Bug 5: Residual range logic may be inverted for `dirsign == 1`
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 429-436
- **Description:** When `dirsign == 1` (forward direction), the residual range is computed as `L.back() - L[fits.size()-1-i]` via a reverse-indexed loop. When `dirsign == -1` (reverse), `rr_vec = L` directly. The convention for residual range is "distance remaining to the end of the track." For a forward-directed segment, the end is the last point, so residual range should decrease along the path -- which `L.back() - L[i]` would give. The current code `L.back() - L[fits.size()-1-i]` effectively assigns `rr_vec[fits.size()-1-i] = L.back() - L[fits.size()-1-i]`, which when iterated i=0..N-1 correctly gives `rr_vec[N-1]=0, rr_vec[N-2]=L[N-1]-L[N-2], ...`. This is actually equivalent to `rr_vec[j] = L.back() - L[j]` for all j, which is correct. So the forward case is correct. For `dirsign == -1`, the "end" is the first point (index 0), and `rr_vec = L` gives `rr_vec[0]=0`, increasing along the path -- this is correct for residual range measured from the start. However, the `dirsign == 0` (unknown) case falls through to the same `rr_vec = L` as `dirsign == -1`, which is an arbitrary choice. No warning is issued.
- **Severity:** Low
- **Why:** The `dirsign == 0` fallback is a reasonable default but should be documented or logged.

### Bug 6: Vertex endpoint residual-range markers may be applied to wrong ends
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 440-450
- **Description:** `find_vertices` returns `(start_vtx, end_vtx)` ordered by proximity to the segment's initial wcpoint. The code marks `rr_vec.front()` as -1 if `start_vtx` has degree > 1, and `rr_vec.back()` as -1 if `end_vtx` has degree > 1. This assumes the fits vector is ordered the same way as wcpts (start to end). If fits were ever reordered (e.g., by direction), the front/back association could be wrong. Currently this appears safe since fits and wcpts maintain consistent ordering, but the assumption is implicit.
- **Severity:** Low
- **Why:** Relies on an implicit ordering guarantee.

### Bug 7: App argument parsing has no bounds checking
- **File:** `wire-cell-uboone-magnify-tracking-convert.cxx`, lines 68-89
- **Description:** The argument parser `argv[i][1]` is accessed without checking `strlen(argv[i]) >= 2`. If an argument like `-` (single dash) is passed, `argv[i][1]` is the null terminator and `&argv[i][2]` points past the string. Additionally, there is no check for the `-` prefix itself; any two-character input would be parsed. No `default` case handles unknown flags.
- **Severity:** Medium
- **Why:** Can cause undefined behavior with malformed arguments. In practice, the app is likely only called from scripts with correct arguments.

### Bug 8: App dtheta index-out-of-bounds when single-point clusters exist
- **File:** `wire-cell-uboone-magnify-tracking-convert.cxx`, lines 537-539
- **Description:** For clusters with only 1 point (the `else` branch at line 537), the code does `dtheta->at(k).push_back(0)`. However, `dtheta` was never given a new sub-vector for this cluster -- the `push_back` for `dtheta` only happens in the `if (x2->at(k).size()>1)` branch (line 495). So `dtheta->at(k)` will access an index that was never created, causing an out-of-bounds access or crash.
- **Severity:** High
- **Why:** This will crash or produce undefined behavior whenever a single-point cluster exists. The `dtheta` vector lacks entries for clusters with size <= 1.

### Bug 9: App memory leaks -- heap-allocated vectors never deleted
- **File:** `wire-cell-uboone-magnify-tracking-convert.cxx`, lines 117-287
- **Description:** Many vectors are allocated with `new` (e.g., `vx`, `vy`, `vz`, `x`, `y`, `z`, `Q`, `max_dis`, `x2`, `y2`, etc.) but never `delete`d. The `file1` TFile is also never closed or deleted.
- **Severity:** Low
- **Why:** Since the app exits shortly after, the OS reclaims memory. But it is poor practice and could matter if the code is ever refactored into a library function.

### Bug 10: App map key uses floating-point discretization that can collide
- **File:** `wire-cell-uboone-magnify-tracking-convert.cxx`, lines 401, 477-480
- **Description:** The `map_point_index` key is created by dividing coordinates by `0.01*units::mm` and (on line 401 using implicit double-to-int truncation via `std::make_tuple` with `double` keys, but on lines 477-480 using explicit `int()` casts). The key types are inconsistent: line 401 stores `double` values (no `int()` cast) while lines 477-480 use `int()` casts for lookup. If the map key type is `tuple<int,int,int>` this means line 401 implicitly truncates via the tuple constructor, which should work, but the inconsistency is fragile. More importantly, two distinct points that round to the same grid cell will collide, causing the wrong index to be returned.
- **Severity:** Medium
- **Why:** Grid collisions silently corrupt the truth-reco charge association. The 0.01mm grid is fine enough that collisions are rare but not impossible for very close points.

---

## 2. Efficiency Improvements

### Efficiency 1: Redundant graph iteration for collecting clusters
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 318-333
- **Description:** The code iterates over all graph edges and all graph vertices to collect `all_clusters` into a `std::set`. Then in the main processing loop (lines 345-472), it iterates over `all_clusters`, and for each cluster, re-iterates over ALL edges and ALL vertices to find those belonging to the current cluster. This is O(C * (E + V)) instead of O(E + V).
- **Suggested improvement:** Build a `map<Facade::Cluster*, vector<edge_descriptor>>` and `map<Facade::Cluster*, vector<vertex_descriptor>>` in a single pass over edges and vertices. Then iterate the map entries.
- **Expected impact:** Medium -- for graphs with many clusters and many edges/vertices, this eliminates redundant iterations.

### Efficiency 2: Multiple lookups into `cluster_*` maps in write_proj_data
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 232-239
- **Description:** When building output vectors, the code looks up `cluster_time_slices[cid]`, `cluster_charges[cid]`, etc. separately for each cluster. These are separate map lookups but always for the same key.
- **Suggested improvement:** Use a single struct per cluster containing all vectors, stored in one map. This reduces from 5 map lookups to 1 per cluster.
- **Expected impact:** Low -- the number of clusters is typically small.

### Efficiency 3: Repeated `nticks_map.at(apa).at(face)` lookup per entry in write_proj_data
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, line 206
- **Description:** `nticks_map.at(apa).at(face)` is looked up inside the inner loop for every `(wire, time)` entry, even though `apa` and `face` only change in the outer loop.
- **Suggested improvement:** Move the lookup to the outer loop body (after line 201, before the inner `for`). The lookup is already at line 206 which is inside the outer loop but outside the inner loop, so this is actually fine. No change needed upon closer inspection.
- **Expected impact:** None (already correctly placed).

### Efficiency 4: App uses `new` for all vectors instead of stack allocation
- **File:** `wire-cell-uboone-magnify-tracking-convert.cxx`, lines 117-287
- **Description:** All vectors are heap-allocated with `new` and passed by pointer. Stack allocation or `std::unique_ptr` would be simpler and avoid the memory leak.
- **Suggested improvement:** Use stack-allocated `std::vector` objects and pass addresses.
- **Expected impact:** Low -- minor cleanup, no performance difference.

### Efficiency 5: Residual range vector allocated then potentially fully overwritten
- **File:** `UbooneMagnifyTrackingVisitor.cxx`, lines 412-413
- **Description:** `rr_vec` is initialized to all zeros, then `L` is computed, then `rr_vec` is either computed from L (dirsign==1) or assigned `rr_vec = L` (dirsign==-1 or 0). The initial zero-fill is wasted work.
- **Suggested improvement:** Defer initialization or use move semantics.
- **Expected impact:** Low -- the vector sizes are small (number of fit points per segment).

---

## 3. Algorithm and Code Explanation

### General Purpose

`UbooneMagnifyTrackingVisitor` is a ROOT file writer that runs as a visitor in the Wire-Cell Toolkit's multi-algorithm blob clustering (MABC) pipeline. It serializes track fitting and reconstruction results into a ROOT file format compatible with the MicroBooNE "Magnify" event display and analysis tools. The companion app `wire-cell-uboone-magnify-tracking-convert` reads this output and combines it with truth-level Monte Carlo information for track reconstruction performance evaluation.

### Key Data Structures

- **`WCPointTree`** (header lines 23-43): A POD-like struct that holds all per-point reconstruction variables. Acts as a reusable buffer whose address is bound to TTree branches; fields are updated and `Fill()` is called for each point.

- **`Clus::Facade::Grouping`**: The top-level clustering container that provides access to dead channels, track fitting results, and detector geometry (nticks_per_slice, anodes, detector volumes).

- **`PR::Graph`**: A Boost graph where vertices hold `PR::Vertex` objects and edges hold `PR::Segment` objects. Each segment contains a vector of `Fit` structs (3D point, charge, projections) and a vector of `WCPoint` structs (original measurement points).

- **`TrackFitting::FittedCharge2D`**: Stores 2D fitted charge results keyed by (APA, face, plane) and (wire, time), associating measured/predicted charges with the clusters that contributed to them.

### Algorithm Walkthrough

#### 1. `visit()` (line 66-102) -- Entry Point
- Finds the named grouping (default "live") from the ensemble.
- Sets anode planes and detector volumes on the grouping (needed for coordinate lookups).
- Opens the output ROOT file in RECREATE mode.
- Calls four write methods in sequence, then writes an empty `T_proj` tree for format compatibility.
- Closes and deletes the ROOT file.

#### 2. `write_bad_channels()` (lines 104-151)
- Creates `T_bad_ch` tree with per-channel dead-channel info.
- Iterates over all unique (APA, face) pairs from the grouping's wire-plane IDs.
- For each of the 3 planes, queries `get_all_dead_chs()` to obtain a list of (channel, time_range) pairs.
- Each dead channel entry is written as a row: channel ID, plane index, start/end time, plus run/event metadata.
- Exceptions from missing APA/face/plane combinations are caught and warned (robustness for multi-APA detectors).

#### 3. `write_trun()` (lines 153-173)
- Creates `Trun` tree with a single entry containing run/subrun/event numbers and the dQ/dx calibration parameters (scale and offset).
- These calibration parameters allow downstream consumers to reverse the charge transformation.

#### 4. `write_proj_data()` (lines 175-252)
- Retrieves `TrackFitting` from the grouping and accesses its `fitted_charge_2d` map.
- The fitted charge data is keyed by (APA, face, plane) -> (wire, time) -> FittedCharge2D.
- For each entry, computes:
  - **Channel number** using MicroBooNE convention: U=wire, V=2400+wire, W=4800+wire.
  - **Time slice** by dividing the raw time by `nticks_per_slice`.
- Groups all entries by cluster_id (from `fc.clusters`).
- Note: a single (wire, time) cell can belong to multiple clusters (via `fc.clusters` set), so the same measurement may appear under multiple cluster IDs.
- Writes a single-entry `T_proj_data` tree where each branch is a `vector<vector<int>>` -- outer index is cluster, inner index is the measurements for that cluster.

#### 5. `write_t_rec_data()` (lines 254-475)
This is the most complex method. It writes the `T_rec_charge` tree containing 3D reconstructed points with charge and fit quality.

**Phase 1: Setup (lines 258-312)**
- Gets the PR graph from TrackFitting.
- Determines `nticks_per_slice` from the grouping (uses the first available APA/face value -- appropriate for uBooNE's single-APA geometry).
- Creates the TTree with 19 branches bound to `WCPointTree` fields.
- Notable: `real_cluster_id` and `sub_cluster_id` both point to `reco_proto_cluster_id` for legacy compatibility.

**Phase 2: Collect clusters and find main cluster (lines 315-342)**
- Scans all edges (segments) and vertices in the graph to collect the set of all referenced clusters.
- Identifies the "mother" cluster (the one with `main_cluster` flag) whose ID becomes the `cluster_id` branch value for all points.

**Phase 3: Process vertices (lines 350-385)**
- Unless `flag_skip_vertex` is set, iterates all graph vertices.
- For each vertex belonging to the current cluster, writes its fitted position (converted to cm), charge (scaled by `dQdx_scale` + `dQdx_offset`), wire projections (pu, pv, pw), time projection (pt / nticks_per_slice), and reduced chi2.
- Vertices are marked with `flag_vertex = 1` and `rr = -1` (no residual range concept for vertices).

**Phase 4: Process segments (lines 388-471)**
- For each segment belonging to the current cluster:
  - Computes a `proto_cluster_id` encoding: `cluster_id * 1000 + graph_index` (allows distinguishing segments within a cluster, up to 1000 segments per cluster).
  - Determines track/shower classification from segment flags.
  - Extracts particle ID from segment's particle info if available, else defaults to 1 (shower) or 4 (track, likely proton PDG convention in this context).
  - **Residual range calculation** (lines 412-451):
    - Computes cumulative path length `L[i]` along the segment.
    - If `dirsign == 1` (forward): residual range = total_length - L[i], so it decreases to 0 at the endpoint.
    - If `dirsign == -1` or `0`: residual range = L[i], increasing from 0 at the start.
    - At branching vertices (degree > 1), the endpoint residual range is set to -1 as a sentinel, indicating the track continues beyond this segment.
  - Writes one TTree entry per fit point with position, charge, projections, and residual range.

### Downstream App: wire-cell-uboone-magnify-tracking-convert

The app reads the `T_rec_charge` tree produced by the visitor and optionally pairs it with Monte Carlo truth from a separate file. Key operations:

- **Cluster grouping**: Uses the `ndf` branch (which holds cluster_id as double) to group consecutive points into clusters. A new cluster starts whenever `std::round(ndf)` changes.
- **dQ/dx reversal**: Reads `dQdx_scale` and `dQdx_offset` from `Trun`, then reverses the transformation: `(dQ1 - dQdx_offset) / dQdx_scale` to recover original charge.
- **Truth matching** (MC mode): Builds a KD-tree of truth points, finds the closest truth point for each reco point, and computes distance/angle metrics for track comparison.
- **Output**: A combined ROOT file with `T_rec` (reorganized reco+truth data), `T_true` (truth tracks), `T_bad_ch` (cloned), `T_proj_data` (cloned), and `T_proj` (cloned).

### Important Edge Cases Handled
- Empty groupings (line 69-72): early return with debug message.
- Missing TrackFitting (lines 178-181, 259-262): early return with warning.
- Empty fitted charge data (lines 184-187): early return.
- Dead channel query failures (lines 142-147): caught exceptions, logged warnings.
- Segments with empty fits or wcpts (line 395): skipped.
- Branching vertices marked with sentinel rr=-1 (lines 443-450).
- Unknown direction (dirsign==0) falls through to same treatment as reverse.
- In the app, missing truth file is handled gracefully (lines 157-159).
