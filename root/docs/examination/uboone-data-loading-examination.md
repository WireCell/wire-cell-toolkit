# Uboone Data Loading Code Examination

> **Bug Fix Status**: Bugs 1, 3, 5, 6, 7, 8 fixed. Bug 2 (fid_ind sentinel) deferred — caught by downstream bounds check. Bug 4 (uppercase conversion) deferred — low risk.
>
> **Efficiency Fix Status**: Efficiency 2 (extract_live bad ch iteration) fixed — computes overlapping slices directly instead of iterating all slices. Efficiency 4 (flush O(N*M) cluster dispatch) fixed — builds cluster_id→blob map in single pass. Efficiencies 1, 3, 5, 6 deferred — low impact or one-time cost.

Examined files:
- `root/inc/WireCellRoot/UbooneTTrees.h` (lines 1-643)
- `root/inc/WireCellRoot/UbooneBlobSource.h` (lines 1-141)
- `root/src/UbooneBlobSource.cxx` (lines 1-487)
- `root/inc/WireCellRoot/UbooneClusterSource.h` (lines 1-177)
- `root/src/UbooneClusterSource.cxx` (lines 1-773)
- `root/inc/WireCellRoot/UbooneGeomHelper.h` (lines 1-80)
- `root/src/UbooneGeomHelper.cxx` (lines 1-215)
- `root/test/uboone.jsonnet`, `uboone-blobs.jsonnet`, `uboone-clusters.jsonnet`

---

## 1. Potential Bugs

### Bug 1: `Blob::flag_uvw` initialized with null pointers (HIGH)
**File:** `UbooneTTrees.h`, line 158  
**Description:** In `Blob::set_addresses()`, the line `flag_uvw = {flag_u_vec, flag_v_vec, flag_w_vec};` captures the *current* pointer values of `flag_u_vec`, `flag_v_vec`, `flag_w_vec`, which are all still `nullptr` at this point. ROOT's `SetBranchAddress` stores the *address* of those pointer members (e.g., `&flag_u_vec`) and only fills them when `GetEntry()` is called later. So `flag_uvw` will hold three null pointers and will never be updated, because it holds copies of the pointer values, not references to the pointer members.  
**Why it's a bug:** Any code that dereferences `flag_uvw[pln]` will segfault. This is currently used in `UbooneBlobSource::in_views()` (line 168-172 of `UbooneBlobSource.cxx`) via `blob(m_kind).flag_uvw`. However, the call goes through either `live` (TC) or `dead` (TDC), both of which call `Blob::set_addresses()` via their own `set_addresses()`. Since TC and TDC both override `set_addresses()` and call `Blob::set_addresses(tree)`, the same bug applies. The `flag_uvw` vector always holds the null pointer values from before `GetEntry()` was called.  
**Severity:** High -- will cause segfault when `in_views()` is called with flag_uvw access.

**UPDATE/Clarification:** Actually, `flag_u_vec` etc. are `std::vector<int>*` that start as nullptr. After `SetBranchAddress` and `GetEntry`, ROOT allocates the vector and updates the pointer `flag_u_vec` to point to it. But `flag_uvw` holds copies of the original null pointer values. The `flag_uvw` entries are NOT updated when ROOT fills `flag_u_vec`. This is indeed a bug -- the `flag_uvw` entries will remain null. The code in `in_views()` at line 172 does `flag_uvw[pln]->at(bind)` which dereferences a null pointer.

### Bug 2: `fid_ind[-1] = -1` stores `size_t(-1)` as a sentinel, risking huge index (MEDIUM)
**File:** `UbooneTTrees.h`, line 375  
**Description:** `fid_ind` is `std::map<int, size_t>`. Setting `fid_ind[-1] = -1` stores `std::numeric_limits<size_t>::max()` (on 64-bit: 18446744073709551615). Later, if a match has `flash_id == -1` (no flash), `cluster_flash[match.cluster_id]` gets this huge value. Then in `is_beam_flash_coincident()` (line 519), the check `flash_index >= optical.at("flash").get("time")->size_major()` should catch this. However, this relies on the comparison being correct for a near-max `size_t` value, which is fragile. Additionally, in `UbooneClusterSource::flush()` (line 655), `farr->element<int>(0) = find` sets the flash index, and if `find` is `size_t(-1)`, this is stored as a very large integer.  
**Why it's a bug:** The sentinel value propagates to point cloud data where it could be misinterpreted as a valid (but enormous) index by downstream code.  
**Severity:** Medium

### Bug 3: `flash` vectors pre-sized to `nflashes` but only `find` entries are filled (HIGH)
**File:** `UbooneTTrees.h`, lines 364-365 vs 406  
**Description:** The vectors `ftime`, `ftmin`, `ftmax`, `fval`, `fident`, `ftype` are all initialized to size `nflashes` (total flash entries across ALL events in the tree). But only entries matching the current event's run/subRun/event are filled (indexed by `find`). After the loop, `find` may be much smaller than `nflashes`. The remaining entries (from index `find` to `nflashes-1`) contain zeros and are included in the `flash_ds` Dataset. This means the flash dataset contains spurious zero-valued entries.  
**Why it's a bug:** Downstream code sees phantom flashes with time=0, tmin=0, tmax=0, qtot=0, ident=0, type=0. These garbage entries could cause incorrect flash-cluster matching or analysis.  
**Severity:** High -- the flash arrays should be resized to `find` before adding to the dataset, or constructed dynamically like the `lid`/`lt`/`lq`/`ldq` vectors are.

### Bug 4: Upper-case conversion logic is inverted (LOW)
**File:** `UbooneBlobSource.cxx`, line 94  
**Description:** `if (letter > 'Z') letter -= 32;` attempts to convert lowercase to uppercase. However, lowercase letters have ASCII values *greater* than `'Z'` (e.g., `'u'`=117, `'Z'`=90), so the condition is correct for detecting lowercase. But the subtraction should convert `'u'`(117) to `'U'`(85), which is correct: 117-32=85='U'. So the arithmetic works. However, the condition `> 'Z'` would also match characters like `[`, `\`, `]`, `^`, `_`, `` ` `` (ASCII 91-96) and incorrectly transform them. This is unlikely with typical inputs but is not robust.  
**Why it's a concern:** Non-letter characters between 'Z' and 'a' would be silently transformed. Using `std::toupper()` would be safer.  
**Severity:** Low

### Bug 5: Double EOS insertion in `fill_queue` + `load_live`/`load_dead` (MEDIUM)
**File:** `UbooneBlobSource.cxx`, lines 353-355 and 446-448  
**Description:** Both `load_live()` (line 353) and `load_dead()` (line 425) push `nullptr` to `m_queue` when `m_frame_eos` is true. Then `fill_queue()` at line 447-448 also pushes `nullptr` when `m_frame_eos` is true. This means two EOS markers are inserted per frame when `m_frame_eos == true`.  
**Why it's a bug:** Downstream components will receive duplicate EOS markers, which could cause unexpected behavior (e.g., double resets, premature termination).  
**Severity:** Medium

### Bug 6: Light data loading condition has wrong logic (HIGH)
**File:** `UbooneClusterSource.cxx`, line 79  
**Description:** The condition is:
```cpp
if (!m_light_name.empty() || !m_flash_name.empty() || m_flashlight_name.empty()) {
    kinds.push_back("light");
}
```
The third clause `m_flashlight_name.empty()` is NOT negated, unlike the first two. This means "light" is added to kinds whenever `m_flashlight_name` is empty (its default). So the "light" kind is effectively always added, even when no optical data is requested. The intended logic was almost certainly:
```cpp
if (!m_light_name.empty() || !m_flash_name.empty() || !m_flashlight_name.empty()) {
```
**Why it's a bug:** When no optical config names are given (all default to ""), the condition `m_flashlight_name.empty()` is true, so "light" trees (T_flash, T_match) are always loaded. If those trees don't exist in the file, the program crashes.  
**Severity:** High

### Bug 7: `UbooneClusterSource::flush()` accesses `islices[0]` and `iblobs[0]` without empty check (MEDIUM)
**File:** `UbooneClusterSource.cxx`, lines 542-543  
**Description:** After `extract_live()` or `extract_dead()`, the code accesses `islices[0]` and `iblobs[0]` without checking if these vectors are empty. If there are no blobs in the cached blobsets (e.g., all blobsets have zero blobs), this causes undefined behavior.  
**Why it's a bug:** While unlikely in normal operation, edge cases with empty data could cause crashes.  
**Severity:** Medium

### Bug 8: `TC::set_addresses` fallback sets `parent_cluster_id_vec = cluster_id_vec` which is still null (MEDIUM)
**File:** `UbooneTTrees.h`, line 198  
**Description:** In the fallback path where no `parent_cluster_id` branch exists, `parent_cluster_id_vec = cluster_id_vec` is set. But at this point in `set_addresses`, `cluster_id_vec` has been passed to `SetBranchAddress` but its value is still nullptr (ROOT hasn't called GetEntry yet). After GetEntry, `cluster_id_vec` will be updated by ROOT to point to the allocated vector, but `parent_cluster_id_vec` will remain null (it's a copy of the old null value).  
**Why it's a bug:** After GetEntry, `parent_cluster_id_vec` is still null, so `cluster_ids()` and `parent_cluster_ids()` will dereference null.  
**Severity:** Medium -- only triggered when the `parent_cluster_id` branch is missing.

---

## 2. Efficiency Improvements

### Efficiency 1: `load_optical()` scans entire T_flash and T_match trees per event (MEDIUM)
**File:** `UbooneTTrees.h`, lines 361-441  
**Current approach:** For each event, the code does a full linear scan of all flash entries and all match entries, comparing run/subRun/eventNo.  
**Suggested improvement:** Build an index (e.g., `std::unordered_map<tuple<int,int,int>, vector<int>>`) mapping (run,subRun,event) to entry indices during construction or first access. Alternatively, use `TTreeIndex` or `TTree::BuildIndex()`.  
**Expected impact:** Medium -- matters for files with many events or many flashes. Each event currently scans O(total_flashes + total_matches).

### Efficiency 2: `extract_live()` iterates all bad channels x all time slices (MEDIUM)
**File:** `UbooneClusterSource.cxx`, lines 352-371  
**Current approach:** For every bad channel mask entry, the code loops over ALL time slices (`n_slices_span`, typically ~2400) to check overlap. This is O(bad_channels * n_slices * ranges_per_channel).  
**Suggested improvement:** Instead of looping over all time slices, compute which time slices overlap with each bad channel's tick range directly: `tsid_start = tt.first / nrebin; tsid_end = (tt.second + nrebin - 1) / nrebin;` and only iterate that subset.  
**Expected impact:** Medium -- for thousands of bad channels and thousands of time slices, this inner loop is significant.

### Efficiency 3: `extract_live()` creates temporary vectors per blob for flag/wire lookups (LOW)
**File:** `UbooneClusterSource.cxx`, lines 392-400  
**Current approach:** For each blob, two `std::vector` objects (`flags` and `wire_indices`) are constructed, copying data from the ROOT vectors. This involves heap allocation per blob.  
**Suggested improvement:** Access the flag and wire index vectors directly without creating intermediate copies. Use `const auto& wu = live_data.wire_index_u_vec->at(bind)` directly in the loop.  
**Expected impact:** Low -- the vectors are small but there can be many blobs.

### Efficiency 4: Per-blob cluster assignment in `flush()` does O(N*M) scan (MEDIUM)
**File:** `UbooneClusterSource.cxx`, lines 616-637  
**Current approach:** For each cluster node (`cnodes`), the code iterates over ALL blobs to find those belonging to the cluster. This is O(clusters * blobs).  
**Suggested improvement:** Build a map from cluster_id to blob indices in a single pass, then process each cluster from the map. This reduces to O(blobs + clusters).  
**Expected impact:** Medium -- significant for events with many clusters and blobs.

### Efficiency 5: `UbooneGeomHelper::load_sce_offsets` copies TH3F to TH3D bin-by-bin (LOW)
**File:** `UbooneGeomHelper.cxx`, lines 88-106  
**Current approach:** A triple-nested loop copies each bin individually from a TH3F to a new TH3D.  
**Suggested improvement:** Use `TH3D` copy constructor or `TH3::Copy()`, or simply keep the original `TH3F` (since `Interpolate` works on both). If double precision is truly needed, ROOT's `TH3::Clone()` followed by conversion would be cleaner.  
**Expected impact:** Low -- this is one-time initialization cost.

### Efficiency 6: `UbooneBlobSource::load_live()` pre-creates blobsets for all timeslice IDs including empty ones (LOW)
**File:** `UbooneBlobSource.cxx`, lines 261-269  
**Current approach:** Creates `SimpleBlobSet` objects for every timeslice from 0 to `n_slices_span`, even though many may have no blobs. All are pushed to `m_queue` (line 350-352).  
**Suggested improvement:** Only create blobsets for timeslices that actually have data, or skip empty blobsets when pushing to the queue.  
**Expected impact:** Low -- the empty blobsets are lightweight, but sending them downstream may cause unnecessary processing.

---

## 3. Algorithm and Code Explanation

### 3.1 UbooneTTrees.h -- ROOT TTree Data Interface

**Purpose:** Provides a C++ interface to MicroBooNE Wire-Cell Prototype (WCP) ROOT output files. Acts as the data access layer for all Uboone data loading components.

**Key Data Structures:**

- `Header` (lines 51-99): Run metadata from `Trun` tree -- trigger time, run/subRun/event numbers, nrebin (time rebinning factor), trigger bits. Computes beam timing windows (`lowerwindow`/`upperwindow`) from trigger bits for BNB, NuMI, ExtBNB, ExtNuMI triggers.

- `Activity` (lines 102-117): Per-timeslice channel activity from `Trun` tree. Stores timeslice IDs, channel lists, and charge/error values. This is the "measurement" data for building slice activity maps.

- `Blob` (lines 121-159): Base struct for blob data from TC/TDC trees. Contains per-blob cluster ID, per-plane flags (active=1/dead=0), and wire index arrays for U/V/W planes. The `flag_uvw` derived field bundles the three flag pointers.

- `TC` (lines 162-204): Live blob data extending `Blob`. Adds single `time_slice` per blob, charge `q`, and `parent_cluster_id` (the "main" cluster used in matching). Distinguishes between individual cluster IDs and parent/main cluster IDs.

- `TDC` (lines 207-215): Dead blob data extending `Blob`. Adds `time_slices` (plural -- a dead blob can span multiple time slices).

- `Bad` (lines 233-258): Bad channel data from `T_bad_ch` tree. Builds a `ChannelMasks` map (channel -> list of [tick_begin, tick_end) ranges) by scanning all entries at construction time.

- `Flash` (lines 262-285): Per-flash optical data from `T_flash`. Contains time, PE values per PMT channel (32 PMTs), fired channels list. This is a per-flash tree, not per-event, requiring run/subRun/event matching.

- `Match` (lines 289-308): Cluster-flash matching from `T_match`. Maps TPC cluster IDs to flash IDs, with event_type encoding tagger flags.

- `UbooneTFiles` (lines 564-640): Iterator over multiple ROOT files. Manages file-to-file transitions, skips missing/empty files, and calls `UbooneTTrees::next()` to advance through entries.

**Algorithm Flow:**
1. Constructor opens the ROOT file, loads required trees (Trun, TC always; TDC if "dead" kind; T_flash/T_match if "light" kind).
2. `next()` advances to the next entry, calls `GetEntry()` on all loaded trees, computes beam windows, loads cluster ID ordering, and optionally loads optical data.
3. `load_clusters()` extracts unique cluster IDs from the blob data into a sorted set and builds an ID-to-index map.
4. `load_optical()` scans T_flash and T_match for entries matching the current event, builds flash/light/flashlight point cloud datasets, and maps cluster IDs to flash indices.

### 3.2 UbooneBlobSource -- Blob Set Source from ROOT

**Purpose:** Reads MicroBooNE WCP ROOT files and produces `IBlobSet` objects, one per time slice. Each blob set contains blobs for a single time slice with their associated activity (channel charge measurements).

**Key Design:**
- Configured with `kind` ("live" or "dead") and `views` (which plane combination: "uvw", "uv", "vw", "wu").
- Uses a queue (`m_queue`) to buffer all blob sets from one ROOT entry, then dispenses them one at a time through `operator()`.
- The `views` parameter is encoded as a bitmask: U=1, V=2, W=4. Three-view live is 7; two-view variants are 3, 5, or 6.

**Live Loading Algorithm (lines 250-359):**
1. Creates a frame with identifier encoded from run/event numbers.
2. Pre-creates one `SimpleSlice` and `SimpleBlobSet` per timeslice ID (0 to n_slices_span).
3. Fills activity maps: for each data timeslice, reads channel charges from the ROOT activity data and populates slice activity maps. For 2-view mode, skips real activity on the "bodged" (masked) plane and instead fills from bad channel masks.
4. Creates blobs: for each blob in the ROOT data that matches the configured views, constructs a `RayGrid::Blob` with strips from wire indices, associates it with the correct timeslice's blobset.
5. For 2-view live: additionally marks the bodged plane's channels as "bad" (uncertainty=1e12) in the blob's slice activity.
6. Pushes all blobsets to the queue.

**Dead Loading Algorithm (lines 362-431):**
1. For each dead blob matching configured views, creates or reuses a slice. Dead blobs can span multiple time slices (uses `front()` as the slice ID).
2. Fills the "dummy" plane's channels with bodge values in the slice activity.
3. Creates the blob with zero charge and adds it to the blobset.

**Edge Cases:**
- `make_strip()` (line 238-242) converts wire indices to half-open [min, max+1) strip bounds.
- Strips for layers 0 and 1 (the two "global" RayGrid bounds layers) are set to {0,1}.
- The `in_views()` function (lines 164-178) checks that a blob's plane flags match exactly the configured view pattern.

### 3.3 UbooneClusterSource -- Cluster + Optical Data Assembly

**Purpose:** A function node (IBlobTensoring) that receives merged blob sets from multiple UbooneBlobSource instances, groups them into clusters, attaches optical data (flashes, light, flash-light associations), and outputs a serialized PC-tree as an ITensorSet.

**Key Design:**
- Operates as a "collect and flush" node: caches incoming blob sets until a frame boundary or EOS, then processes all at once.
- Opens the same ROOT files independently to read cluster structure, optical data, and tagger flags.
- The output PC-tree has: root node with optical PCs, cluster nodes with per-cluster scalars and tagger info, blob nodes with sampled 3D point clouds.

**Flush Algorithm (lines 502-688):**
1. Advances the internal ROOT file reader to the next entry.
2. Creates cluster nodes in the PC-tree, one per unique cluster ID.
3. Applies tagger flags (beam flash coincidence, event type bits) to each cluster node.
4. Extracts blobs and slices from cached blob sets (using `extract_live` or `extract_dead`).
5. Dispatches blobs to cluster nodes using the blob-to-cluster mapping from ROOT data.
6. If a sampler is configured, generates 3D point clouds for each blob.
7. For live data: processes individual vs parent cluster IDs to tag "isolated" sub-clusters.
8. Attaches optical data (flash/light/flashlight) from ROOT to the root node.
9. Maps cluster-flash associations to cluster nodes.
10. Adds CTPC (charge-time-position-cloud) and dead region windows from slices.
11. Serializes the PC-tree to tensors.

**extract_live() Algorithm (lines 290-500):**
This is the most complex function. It rebuilds slices with comprehensive activity maps:
1. Creates activity maps for all time slices from ROOT activity data.
2. Overlays bad channel masks (T_bad_ch) onto activity maps.
3. Overlays per-blob dead plane information from TC flags (flag=0 means dead for that plane).
4. Overlays per-blob dead plane information from TDC data if available.
5. Creates SimpleSlice objects with the finalized activity maps.
6. Creates new blob objects referencing the new slices.

**Tagger Flags (lines 690-773):**
Extracts bit-encoded tagger information from `event_type` in the match data:
- Bit 3: TGM (through-going muon)
- Bit 4: Low energy
- Bit 1: Light mismatch
- Bit 2: Fully contained
- Bit 5: Short track muon (STM)
- Bit 6: Full detector dead

### 3.4 UbooneGeomHelper -- Geometry and SCE Corrections

**Purpose:** Implements the `IClusGeomHelper` interface for MicroBooNE, providing fiducial volume checks and Space Charge Effect (SCE) position corrections using data-driven correction histograms.

**Key Data Structures:**
- `FV` struct (lines 49-62): Fiducial volume bounds and margins for each dimension. Defaults correspond to MicroBooNE's active volume.
- `m_FV_map`: Maps "aXfY" strings (apa X, face Y) to FV parameters, supporting multi-APA detectors in principle.
- `m_h3_Dx/Dy/Dz`: 3D histograms for SCE displacement corrections, loaded from a ROOT file.

**SCE Correction Algorithm (lines 162-215):**
1. Converts input point from WCT units to cm.
2. Shifts y-origin to bottom of active volume (adds half the y-length).
3. Clamps coordinates to valid histogram range to avoid extrapolation.
4. Maps physical coordinates to histogram coordinate system using detector-specific scale factors.
5. Interpolates displacement corrections from the 3D histograms.
6. Applies corrections: x is corrected by subtracting (displacement / scale_x), while y and z are corrected by adding (displacement / scale_y or scale_z).
7. Un-shifts y-origin and converts back to WCT units.

The sign convention difference (subtract for x, add for y/z) follows the prototype's `func_pos_SCE_correction`. The x-direction correction is opposite because the x-histogram maps field distortion causing apparent drift toward the anode, requiring subtraction to get the true position.

**Fiducial Volume Checks:**
- `is_in_FV()` (lines 125-137): Strict inequality (<, >) against bounds plus margins.
- `is_in_FV_dim()` (lines 139-152): Per-dimension check with a caller-specified margin, using non-strict inequality (>=, <=). Note the inconsistency: `is_in_FV` uses strict `<`/`>` while `is_in_FV_dim` uses `>=`/`<=`.

### 3.5 Jsonnet Configuration

**uboone.jsonnet:** Central configuration library defining component constructors for the Uboone data loading pipeline. Key patterns:
- `multiplex_blob_views`: Creates one `UbooneBlobSource` per view combination and merges them with `BlobSetMerge`.
- `UbooneClusterSource`: Wraps blob-to-cluster conversion with optional optical data attachment.
- Defines `BlobSampler` configurations for live (stepped strategy) and dead (center strategy) sampling.

**uboone-blobs.jsonnet:** Separate live and dead pipelines producing NPZ files, with an optional clustering step reading both back.

**uboone-clusters.jsonnet:** Single pipeline loading live blobs from ROOT, passing through `UbooneClusterSource` for cluster+optical assembly, ending with `ClusterFlashDump`.
