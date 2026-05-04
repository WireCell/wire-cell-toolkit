# Improve_PR3DCluster Port Review

**Reviewed:** 2026-04-10
**Functions:** `Improve_PR3DCluster_2`, `Improve_PR3DCluster_1`, `Improve_PR3DCluster` (single-cluster and dual-cluster overloads) and all callees.
**Scope:** Functional equivalence, bugs, efficiency, determinism, multi-TPC/APA/face correctness.

---

## 1. Function Map (prototype → toolkit)

| Prototype (`prototype_base/pid/src/ImprovePR3DCluster.cxx`) | Toolkit (`clus/src/`) | Status |
|---|---|---|
| `Improve_PR3DCluster_2` @ L9 — top-level retiling driver | `ImproveCluster_2::mutate` — `improvecluster_2.cxx:77` | Ported |
| `Improve_PR3DCluster_1` @ L128 — dead + good channel fill | `ImproveCluster_1::mutate` — `improvecluster_1.cxx:87` | Ported (see §2.1) |
| `Improve_PR3DCluster` (single) @ L674 — path-tube + retile | folded into `get_activity_improved` + `hack_activity_improved` | Ported (see §2.2) |
| `Improve_PR3DCluster` (dual) @ L1351 — merge cluster1 path + cluster2 charges | two `hack_activity_improved` calls in `ImproveCluster_2::mutate:185,190` | Ported (see §2.3) |
| helper: `WCPPID::calc_sampling_points` | `Aux::sample_live` in both `mutate` loops | Ported differently — see §2.4 |

Sub-function map for `ImproveCluster_1` helper methods:

| Prototype stage | Toolkit helper | File:Line |
|---|---|---|
| Blob-wire harvest (`_1` L147–206) | `get_activity_improved` step 1 | `improvecluster_1.cxx:331–371` |
| Dead-channel fill (`_1` L208–315) | `get_activity_improved` step 2 | `:379–413` |
| Good-channel overlap (`_1` L340–427) | `get_activity_improved` step 3 | `:417–455` |
| `LowmemTiling::init_good_cells_with_charge` | `make_iblobs_improved` | `:846–944` |
| Connected-component filter (`_1` L499–639) | `remove_bad_blobs` | `:712–844` |
| Path-tube insertion (single L960–1097; dual L1632–1901) | `hack_activity_improved` | `:541–709` |
| `establish_same_mcell_steiner_edges` + Dijkstra + `remove` | `Steiner::Grapher`: `establish_same_blob_steiner_edges`, `graph_algorithms.shortest_path`, `remove_same_blob_steiner_edges` | `improvecluster_2.cxx:107,113,119,148,154,159` |

---

## 2. Functional Equivalence

### 2.1 `Improve_PR3DCluster_1` ↔ `ImproveCluster_1::mutate`

The prototype `Improve_PR3DCluster_1` does **not** contain any path-tube logic — it
only fills dead-channel and good-channel regions. The toolkit's corresponding step in
`ImproveCluster_1::mutate` correctly matches this by having the `hack_activity_improved`
call commented out (`improvecluster_1.cxx:139`). **This is correct behaviour**, but
the commented-out line is misleading. It should be deleted with an explanatory comment:

```cpp
// hack_activity_improved is intentionally NOT called here.
// The prototype's Improve_PR3DCluster_1 (L128-672) adds dead/good channels only.
// Path-tube hacking happens only in the _2-level via ImproveCluster_2::mutate.
```

### 2.2 `Improve_PR3DCluster` (single-cluster) ↔ `get_activity_improved` + `hack_activity_improved`

The single-cluster overload contributes:
- Blob-wire harvest (identical in structure — see §5 Bug #1 re: missing face filter)
- Dead-channel fill with 20 cm cut: **matches** — `dis_cut = 20 * units::cm` at `:377`
- Good-channel overlap with 20 cm cut: **matches** at `:417–455`
- Path-tube insertion: see §2.5 for the quantitative differences

### 2.3 `Improve_PR3DCluster` (dual-cluster) ↔ `ImproveCluster_2::mutate` two hack calls

In the prototype the dual-cluster overload (L1351):
1. Harvests wires from **cluster2** (the path-extended temp cluster).
2. Does dead/good channel fill using **cluster1**'s point cloud as the proximity reference.
3. Adds path tubes for both **cluster2**'s path and **cluster1**'s path.

In the toolkit (`improvecluster_2.cxx`):
1. `get_activity_improved(*orig_cluster, ...)` harvests from the **orig_cluster** (≡ prototype cluster1). **Divergence from prototype.** The prototype dual-overload harvests from `cluster2` (the temp cluster); the toolkit starts from `orig_cluster`. This is a semantic difference — needs owner decision on intent. The rest of the flow (dead/good fill against `orig_cluster`, two hack passes) does correctly mirror the prototype dual-overload.
2. `hack_activity_improved(*orig_cluster, ..., orig_path_point_indices, ...)` at `:185` — adds tube from orig path. **Matches** prototype L1773–1901.
3. `hack_activity_improved(temp_cluster, ..., temp_path_point_indices, ...)` at `:190` — adds tube from temp path. **Matches** prototype L1632–1768.

**Summary:** The toolkit `get_activity_improved` seeds from `orig_cluster` (cluster1 geometry) rather than from `temp_cluster` (cluster2 geometry). Whether this is intentional should be confirmed against the prototype. In practice the dead/good channel fill in both runs uses `orig_cluster`'s kd2d, so the dominant effect is the same, but the initial wire harvest differs.

### 2.4 Charge semantics in `make_iblobs_improved`

The prototype uses `LowmemTiling::init_good_cells_with_charge` which stores real
charge values on each wire of the new `SlimMergeGeomCell`. Cells that were added from
dead-channel or path-tube regions are given `charge=0, charge_err=0`.

The toolkit uses `1e-3` as a sentinel for "dead/forced" channels in `measures[plane][wire]`
and then converts it to `ISlice::value_t(0.0, 1e12)` in `make_iblobs_improved:908–909`.
This gives those channels a non-zero entry in the tiling activity (so they contribute a
blob strip) but a very large uncertainty. This is a deliberate improvement over the
prototype (dead channels cannot be used for charge estimation) and is **correct**. However,
it should be documented as an intentional divergence.

### 2.5 Path-tube geometry: quantitative differences

| Parameter | Prototype PID (`_2`, single, dual) | Prototype 2dtoy | Toolkit (`hack_activity_improved`) |
|---|---|---|---|
| Coverage weight: self (delta=0) | 2 | 1 | 2 |
| Coverage weight: left neighbor (delta=-1) | 1 | 2 | 1 |
| Coverage weight: right neighbor (delta=+1) | 1 | 1 | 1 |
| Coverage threshold (per plane) | `>=2` | `>0` | `>=2` |
| Coverage total threshold | n/a | `sum>=6` | n/a |
| Tube footprint | L2 disk: `dw²+dt²<=9` | L1 diamond: `|dw|+|dt|<=3` | L2 disk: `dw²+dt²<=9` |
| Added charge value | `0` | `0` | `1e-3` (sentinel) |

The toolkit uses the **PID variant** thresholds and geometry (L2 disk, per-plane ≥2),
not the 2dtoy variant (L1 diamond, total sum ≥6). This is the correct source-of-truth
since `ImproveCluster_2` is the port of the PID `Improve_PR3DCluster_2`.

**Coverage weights match the prototype.** The earlier draft of this table incorrectly
stated `delta=-1` has weight 2 in the prototype. Reading `ImprovePR3DCluster.cxx:1011`
shows `nu += 2` for `results.at(1)` (the self wire) and `nu++` for `results.at(1)-1`
(the left neighbor), i.e. self=2, left=1, right=1 — identical to the toolkit.
No divergence here; no further action required.

**Prototype ref:** `ImprovePR3DCluster.cxx:1002–1021` (coverage weights) and `:1047–1071` (tube).
**Toolkit ref:** `improvecluster_1.cxx:577–595` (coverage weights) and `:644–653` (tube).

---

## 3. Bugs

### Bug #1 — Missing face filter in `get_activity_improved` blob-harvest loop (CRITICAL)

**File:** `improvecluster_1.cxx:331–371`
**Severity:** Bug — incorrect result for any cluster spanning multiple APA faces.

The blob-harvest loop iterates `cluster.children()` without filtering by `(apa, face)`:

```cpp
// improvecluster_1.cxx:331-333
auto children = cluster.children();
for (auto child : children) {
    auto blob = child->value().facade<Blob>();
```

Compare with the base class `RetileCluster::get_activity` which explicitly checks:

```cpp
// retile_cluster.cxx:151 (reference implementation)
if (blob_wpid.apa() != apa || blob_wpid.face() != face) continue;
```

For a single-face cluster (e.g., MicroBooNE) this is harmless. For a multi-face cluster,
blobs from face B will be inserted into `u_time_chs` / `v_time_chs` / `w_time_chs` for
face A, and those wire indices will then be written into `measures[plane+2][wire_idx]`
which is sized to `m_plane_infos.at(apa).at(face)[plane].total_wires`. If a face-B blob
has a wire index that exceeds face-A's wire count, the subsequent `m[ch] = charge` at
`:482` will write out of bounds. If the counts are the same (symmetric faces), it will
silently map physically incorrect wires.

**Fix:** Add the same WPID check at the top of the blob-harvest loop:

```cpp
auto blob_wpid = blob->wpid();
if (blob_wpid.apa() != apa || blob_wpid.face() != face) continue;
```

### Bug #2 — Unguarded `begin()` dereference on potentially empty `map_slices_measures`

**Files:** `improvecluster_1.cxx:222`, `improvecluster_2.cxx:241`
**Severity:** Crash (undefined behaviour) on any face with zero tiled blobs.

```cpp
// improvecluster_1.cxx:222
int tick_span = map_slices_measures.begin()->first.second - map_slices_measures.begin()->first.first;
```

If `get_activity_improved` produces no slices for a given `(apa, face)` — which can
happen for a multi-face cluster where one face has no blobs within the cluster's bounding
box — `map_slices_measures.begin()` equals `end()` and dereferencing it is UB.

The same pattern also exists in the base class at `retile_cluster.cxx:672`.

**Fix:** Add an early-continue before this line:

```cpp
if (map_slices_measures.empty()) continue;
int tick_span = map_slices_measures.begin()->first.second - map_slices_measures.begin()->first.first;
```

Also add the same guard before the call to `hack_activity_improved` at `:581` which
makes the same assumption (reads `begin()->first` for `tick_span`).

### Bug #3 — `tick_span` derived from last blob, not from face metadata

**File:** `improvecluster_1.cxx:328,339`
**Severity:** Silent incorrect behavior if blobs in a face have unequal tick spans.

```cpp
// improvecluster_1.cxx:328-339
int tick_span = 1;
...
tick_span = time_slice_max - time_slice_min;  // overwritten per-blob
```

`tick_span` is overwritten for every blob; only the last blob's value survives.
If a cluster somehow contains blobs with different tick spans (or the cluster was built
from two different frame configurations), the iteration over dead/good channels at
`:382,394,405` uses the wrong step size for all except the last blob's face.

**Fix:** Derive `tick_span` from the grouping's per-face tick metadata rather than
from blob geometry:

```cpp
const double tick_s = m_grouping->get_tick().at(apa).at(face);
const int tick_span = static_cast<int>(std::round(tick_s));
```

Or, equivalently, assert that all blobs agree and log an error if they do not.

### Bug #4 — Null-dereference risk in `CreateSteinerGraph::visit` log line

**File:** `CreateSteinerGraph.cxx:91`
**Severity:** Crash if no cluster has `Flags::main_cluster` set.

```cpp
// CreateSteinerGraph.cxx:91 — before the null-check at L94
SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: {} clusters with beam_flash flag. {}", filtered_clusters.size(), main_cluster->ident());
```

`main_cluster` is `nullptr` if no cluster in the grouping carries `Flags::main_cluster`.
The guard `if (main_cluster != nullptr)` is at L94, but the log on L91 fires
unconditionally before that check.

**Fix:** Move the log after the guard, or make it conditional:

```cpp
SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: {} clusters with beam_flash flag. main={}",
    filtered_clusters.size(), main_cluster ? main_cluster->ident() : -1);
```

### Bug #5 — Early `return` in main-cluster failure skips associated-cluster processing

**File:** `CreateSteinerGraph.cxx:164–170`
**Severity:** Silent data loss — associated clusters are not processed when main cluster fails.

```cpp
// CreateSteinerGraph.cxx:164-170
if (!new_cluster.has_graph("steiner_graph")) {
    SPDLOG_LOGGER_WARN(...);
    grouping.destroy_child(new_cluster_ptr, true);
    return;  // <-- returns from visit() entirely
}
```

The associated-cluster loop is at L204–276. When this `return` fires, all associated
clusters for this event are skipped. The associated loop correctly uses `continue` for
its own failure at L254. The main-cluster branch should behave consistently:

```cpp
// Replace return with goto or restructure:
grouping.destroy_child(new_cluster_ptr, true);
goto process_assoc_clusters;  // or restructure into helper
```

A cleaner fix is to extract both branches into a `process_cluster()` helper and call it
from a single loop over `filtered_clusters`, with an `is_main` flag for the extra
kd/steiner probe at L187–190.

### Bug #6 — `remove_bad_blobs` representative-vertex non-determinism

**File:** `improvecluster_1.cxx:712–844`
**Severity:** Non-deterministic results across runs (see also §5).

The vertex IDs assigned to blobs (L725–731) depend on the iteration order of
`new_time_blob_map`'s inner `set<const Blob*>`:

```cpp
// improvecluster_1.cxx:725-731
for (const auto& [time_slice, new_blobs] : new_time_blob_map) {
    for (const Blob* blob : new_blobs) {   // pointer-ordered set
        map_index_blob[index] = blob;
        map_blob_index[blob]  = index;
        all_new_blobs.push_back(blob);
        index++;
    }
}
```

When `num_components > 1`, the "representative" for each component is the blob that
first receives that component's ID (the first vertex in that component encountered in
component-vector order). Whether a component is marked "good" depends entirely on which
blob is chosen as its representative (L783–829). Different heap layouts → different
pointer ordering → different representatives → potentially different physics outcomes.

**Fix:** Sort blobs deterministically before assigning vertex IDs (see §5).

### Bug #7 — `ImproveCluster_2` seeds `get_activity_improved` from `orig_cluster`, not `temp_cluster`

**File:** `improvecluster_2.cxx:180`
**Severity:** Confirmed functional divergence from prototype (see §2.3 for context).

In the prototype dual-cluster overload (`ImprovePR3DCluster.cxx:1369`), the initial
wire-harvest iterates `cluster2->get_mcells()` (the dead-channel-extended temp cluster),
while proximity checks use `cluster1->get_point_cloud()`. In the toolkit,
`get_activity_improved(*orig_cluster, ...)` does both harvest and proximity from
`orig_cluster`. This means the initial wire inventory is smaller (it doesn't include the
dead-channel extensions already added by `ImproveCluster_1`), which may reduce the
coverage of the final tiled cluster.

**Status:** Confirmed divergence. Fix requires splitting `get_activity_improved` into a
harvest step (operating on `temp_cluster`) and a proximity step (operating on
`orig_cluster`'s kd2d). Owner decision required before implementing.

---

## 4. Efficiency Improvements

### 4.1 `ImproveCluster_2` runs the per-face activity build twice

`ImproveCluster_2::mutate` calls `ImproveCluster_1::mutate(node)` at `:129`, which
internally runs a full per-face `get_activity_improved` + `make_iblobs_improved` +
`remove_bad_blobs` loop (N passes for N faces). Then `ImproveCluster_2::mutate` runs
its own per-face loop (another N passes) with `get_activity_improved` +
`hack_activity_improved` × 2 + `make_iblobs_improved` + `remove_bad_blobs`. The first
loop's activity data is discarded (only the blob geometry of `temp_cluster` is used).
For clusters with K faces and M time slices this is 2×K dead/good-channel queries and
2×K KD lookups that could be reduced to K.

**Proposed fix:** Expose a variant of `ImproveCluster_1::mutate` that returns
the `map_slices_measures` per face alongside the cluster node, so `ImproveCluster_2`
can reuse it instead of rebuilding.

### 4.2 ~70 lines of copy-paste in `CreateSteinerGraph::visit`

Lines L94–199 (main cluster) and L204–276 (associated clusters) are near-identical
pipelines: `mutate` → `find_graph` → `Grapher` → `establish_same_blob_steiner_edges`
→ `get_two_boundary_wcps` → `shortest_path` → `remove_same_blob_steiner_edges`
→ `create_steiner_tree` → safety-check → transfer pc/graph → `destroy_child`.

The only difference is that the main-cluster branch also calls `kd_steiner_knn` (L189)
and stores the `new_cluster`'s Steiner output on `main_cluster`. Any bug fixed in one
branch must currently be fixed in the other.

**Proposed fix:** Factor into a helper:

```cpp
bool process_cluster_steiner(Cluster* cluster, Cluster* ref_cluster,
                              Grouping& grouping, bool is_main);
```

### 4.3 `wpids_blob()` always rebuilds a `std::set` at every call site

Every caller of `wpids_blob()` does:

```cpp
auto wpids = orig_cluster->wpids_blob();  // returns vector
std::set<WirePlaneId> wpid_set(wpids.begin(), wpids.end());
```

(`improvecluster_1.cxx:101–102`, `improvecluster_2.cxx:165–166`, `retile_cluster.cxx:619–621`)

**Proposed fix:** Add `Cluster::wpids_blob_set()` returning a `std::set<WirePlaneId>` directly.

### 4.4 Connected-component filter uses Boost BGL where union-find suffices

`remove_bad_blobs` constructs a full `boost::adjacency_list` (`improvecluster_1.cxx:741–762`)
just to call `connected_components`. For the typical N (~10–200 blobs), a simple union-find
over `int` IDs would be faster, require no Boost headers in this translation unit, and
would naturally accept deterministic vertex IDs.

### 4.5 KD query inside the dead-channel loop fires once per (channel, time-slice) pair

For a dead wire spanning many time slices, the inner loop at `:381–390` issues one KD
query per time slice. The KD query result varies by time slice only through the 2D point
projection `(x_pos, y_pos)`. Batching by wire channel (query once per channel with the
midpoint time slice and use a slightly relaxed radius) could halve the number of queries
for long dead regions.

### 4.6 Triplicated U/V/W blocks in `get_activity_improved`

Steps 2 and 3 of `get_activity_improved` contain three identical dead-channel and
good-channel blocks (`improvecluster_1.cxx:380–413` and `:421–455`), one per plane.
Refactoring into a single loop over `plane ∈ {0,1,2}` would reduce the code by 2/3
and prevent future drift. A `dead_chs_ranges[3]` array and `time_chs[3]` array of
pointers would suffice:

```cpp
const auto* dead_ranges[3] = {&dead_uchs_range, &dead_vchs_range, &dead_wchs_range};
std::map<int, std::set<int>>* time_chs[3] = {&u_time_chs, &v_time_chs, &w_time_chs};
for (int pl = 0; pl < 3; ++pl) {
    for (const auto& [start, end] : *dead_ranges[pl]) { ... }
}
```

---

## 5. Determinism: Pointer-Keyed Containers

| Container | Location | Key type | Iteration? | Risk |
|---|---|---|---|---|
| `map_blob_index: std::map<const Blob*, int>` | `improvecluster_1.cxx:721` | pointer | lookup-only | Low — only `.at()` calls |
| `map_index_blob: std::map<int, const Blob*>` | `:720` | int | yes | Safe |
| `new_blobs: std::set<const Blob*>` (inner of `time_blob_map`) | `:725` | pointer | **yes** — vertex ID assignment | **High — drives CC labels** |
| `blobs_to_remove: std::set<const Blob*>` (return value) | `:768` | pointer | **yes** — callers iterate to remove | **Medium** |
| `slice_activity: map<IChannel::pointer, value_t>` | `:909` | pointer | by `RayGrid::make_activities` | Confirm downstream |

### 5.1 Primary issue: vertex-ID assignment in `remove_bad_blobs`

The loop at `improvecluster_1.cxx:725–731` iterates `new_time_blob_map`'s inner
`std::set<const Blob*>`. Because `std::set<T*>` sorts by pointer value (heap-dependent),
the vertex IDs differ between runs. Since the CC representative for each component is
the vertex with the smallest ID in that component, different runs may designate different
blobs as representatives and thus different components as "good".

**Fix:** Sort blobs deterministically before assigning vertex IDs. Use `blob->ident()` (a
stable integer) as the sort key:

```cpp
std::vector<const Blob*> all_blobs_sorted;
for (const auto& [ts, bset] : new_time_blob_map)
    for (const Blob* b : bset)
        all_blobs_sorted.push_back(b);
std::sort(all_blobs_sorted.begin(), all_blobs_sorted.end(),
    [](const Blob* a, const Blob* b){ return a->ident() < b->ident(); });
for (int i = 0; i < (int)all_blobs_sorted.size(); ++i) {
    map_index_blob[i] = all_blobs_sorted[i];
    map_blob_index[all_blobs_sorted[i]] = i;
}
```

Also change the return type from `std::set<const Blob*>` to `std::vector<const Blob*>`
(or at minimum sort the returned set's iteration in the caller) to ensure removal order
is deterministic.

### 5.2 `slice_activity` keyed by `IChannel::pointer`

`make_iblobs_improved` populates `SimpleSlice::slice_activity` using `IChannel::pointer`
keys (`improvecluster_1.cxx:905–915`). If `RayGrid::make_activities` iterates the
`slice_activity` map internally (rather than using the `measures` vectors), iteration
order is pointer-dependent. **Action:** Confirm whether `make_activities` reads from
`measures` (deterministic) or from `slice_activity` (non-deterministic). If the latter,
key by channel ident or wire index.

### 5.3 `std::set<WirePlaneId>` iteration order

The per-face loop iterates `wpid_set` (`improvecluster_1.cxx:124`). `WirePlaneId`
has a well-defined `operator<` (comparison on packed bits), so iteration is
**deterministic across runs for the same detector configuration**. This is safe as-is.
Worth adding a comment to confirm this is deliberate.

---

## 6. Multi-Face / Multi-APA Correctness

### 6.1 Per-face loop is present and structurally correct

Both `ImproveCluster_1::mutate` (`:124–230`) and `ImproveCluster_2::mutate` (`:171–249`)
loop over a `std::set<WirePlaneId>` derived from `orig_cluster->wpids_blob()`, covering
all `(apa, face)` pairs the cluster spans. Each helper is called with explicit `apa, face`
parameters, and every grouping API (`get_overlap_dead_chs`, `get_overlap_good_ch_charge`,
`convert_time_wire_2Dpoint`, `convert_3Dpoint_time_ch`, `kd2d`, `time_blob_map`, `get_tick`)
is called with `(apa, face)`. This is the correct multi-face generalization of the
prototype's single-face implicit design. **No structural gaps found.**

### 6.2 Bug #1 (face filter) is the main multi-face correctness problem

As described in §3 Bug #1, the blob-harvest loop in `get_activity_improved` processes
all cluster blobs regardless of face. This is the only place in the per-face code path
where the face gate is missing. All dead-channel and good-channel queries already carry
explicit `(apa, face)` parameters and are correct.

### 6.3 Path-tube interpolation across face boundaries

`hack_activity_improved` interpolates the path between consecutive path points and
assigns each interpolated point a `WirePlaneId` using `get_wireplaneid(p1, wpid_p, wpid2,
m_dv)` (`improvecluster_1.cxx:566`). For path segments that cross an (apa=0,face=0) →
(apa=0,face=1) boundary, the interpolated midpoints receive one face's WPID or the
other depending on `get_wireplaneid`'s logic. The per-face gate at `:588` then either
includes or excludes these midpoints from the current face's tube.

This is the same cross-face interpolation used by the base class `RetileCluster::hack_activity`
and was reviewed/confirmed acceptable in the Steiner graph review. It is worth verifying
on a PDHD 4-face cluster that the tube around a cross-face track segment is not
discontinuous.

### 6.4 `tick_span` per-face correctness

`m_grouping->get_tick().at(apa).at(face)` is used per-face in `make_iblobs_improved:854`.
The blob-harvest `tick_span` (Bug #3) is the only place tick is derived from blob geometry
rather than from the face metadata. Once Bug #3 is fixed, tick handling will be fully
per-face.

### 6.5 Wire-index vs channel-index in dead/good-channel queries

`get_overlap_dead_chs` and `get_overlap_good_ch_charge` are queried with `(min_ch, max_ch)`
bounds derived from `cluster.get_uvwt_min(apa, face)` (`:269–271`). In MicroBooNE,
wire index and channel index coincide. In detectors with non-trivial wire-to-channel
mappings (e.g., wrapped wires), these may differ. The inline comments at `:383`, `:395`,
`:406`, `:424`, `:436`, `:448` explicitly flag this:

```
// this should be wire index ... , so a mismatch between MicroBooNE and other detectors,
// need to be fixed at some points ...
```

Until this is resolved, the dead/good-channel fill will give incorrect results for any
detector where wire indices and channel IDs differ in the queried range. For single-face
MicroBooNE-style detectors this is harmless. For PDHD this needs verification.

### 6.6 `CreateSteinerGraph.cxx` is face-agnostic and delegates correctly

`CreateSteinerGraph::visit` treats clusters as opaque objects and delegates all per-face
logic to `IPCTreeMutate::mutate()`. The visit loop at L204–276 iterates all associated
clusters regardless of which face they belong to. There is no face-specific logic in
`CreateSteinerGraph.cxx` itself, which is the correct design. **No multi-face issues
found at this level.**

### 6.7 `SteinerFunctions::improve_grapher*` stubs

All four functions in `SteinerFunctions.cxx` (`improve_grapher`, `improve_grapher_1`,
`improve_grapher_2`, and the two-grapher overload of `improve_grapher`) raise
`LogicError("not implemented")`. A grep for callers in the toolkit shows zero active
call sites — they are dead code. These stubs should either be removed or retained with
a comment explaining they are reserved for a future port of the analogous prototype
path (which is itself unused in the prototype's NeutrinoID pipeline).

### 6.8 No jsonnet currently wires `ImproveCluster_2` into `CreateSteinerGraph`

`cfg/pgrapher/common/clus.jsonnet` defines `improve_cluster_2(...)` at line 414 but no
existing top-level configuration passes it as the `retiler` to `steiner(...)`. Until
an integration config is added, multi-face behavior of `ImproveCluster_2` cannot be
exercised in a full pipeline run. A multi-face smoke test (protoDUNE-HD or SBND-style
config) should be added.

---

## 7. Minor / Cosmetic Notes

- `ImproveCluster_1::configure` (`:30–75`) calls `RetileCluster::configure(cfg)` and
  then repeats the same `NeedDV::configure`, `NeedPCTS::configure`, sampler-parsing,
  and anode-parsing logic verbatim. The duplication is harmless now but will silently
  diverge if `RetileCluster::configure` is updated. Delete the redundant block and
  rely on the base-class call only.

- Stale comment in `retile_cluster.h:107–110`: `"fixme: this restricts the retiling to
  single-anode-face clusters"`. The code now correctly handles multiple faces via
  `wpid_set` iteration. Remove the comment.

- Comments in `improvecluster_1.cxx:465,491,516`: `"what to do the first two views???"`.
  These layers (`measures[0]` and `measures[1]`) receive `{1}` following the same
  convention as `RetileCluster::make_iblobs`. The comment should be replaced with
  an explanation of what layers 0 and 1 represent in the RayGrid scheme.

- Large blocks of commented-out debug `std::cout` and `SPDLOG_DEBUG` calls scattered
  throughout `improvecluster_1.cxx` and `improvecluster_2.cxx`. These should be cleaned
  up or promoted to proper `SPDLOG_LOGGER_TRACE` calls.

- `improvecluster_2.cxx:143`: `//temp_steiner_grapher.get_graph("basic_pid");` — this
  commented-out line hints at a rework. The active code at `:141` calls
  `temp_cluster.find_graph("ctpc_ref_pid", *orig_cluster, ...)` directly on the cluster
  (not via the local `temp_steiner_grapher`). The `temp_steiner_grapher` is then used
  only for `establish_same_blob_steiner_edges` and `remove_same_blob_steiner_edges`.
  This is correct but the commented-out line is confusing; delete it.

---

## 8. Prioritized Action List

### P0 — Correctness bugs

1. **Add face filter to `get_activity_improved` blob-harvest loop** (Bug #1,
   `improvecluster_1.cxx:331`). Affects all multi-face clusters; potential out-of-bounds
   write.
2. **Guard `map_slices_measures.begin()` against empty map** (Bug #2,
   `improvecluster_1.cxx:222`, `improvecluster_2.cxx:241`). Crash on any face
   producing zero slices.
3. **Fix null-dereference in `CreateSteinerGraph.cxx:91` log line** (Bug #4). Crash
   on any event with `beam_flash` clusters but no `main_cluster`.
4. **Fix early `return` in `CreateSteinerGraph.cxx:169`** (Bug #5). Data loss — all
   associated clusters are silently skipped.

### P1 — Determinism

5. **Sort blobs deterministically before assigning vertex IDs in `remove_bad_blobs`**
   (§5.1). Use `blob->ident()` as sort key; change return type to `std::vector`.
6. **Verify `slice_activity` pointer-key iteration path** (§5.2). Confirm
   `make_activities` reads only from `measures`, not from `slice_activity`.

### P2 — Functional alignment / review decisions

7. **Fix `get_activity_improved` seeding in `ImproveCluster_2`** (Bug #7, §2.3).
   Confirmed divergence: toolkit seeds from `orig_cluster`, prototype seeds from
   `cluster2` (temp_cluster). Fix requires splitting harvest from proximity in
   `get_activity_improved`. Owner decision required.
8. ~~**Confirm coverage-weight change** (§2.5)~~ — **Resolved**: weights match the
   prototype exactly (self=2, left=1, right=1). The earlier draft table was wrong.
   No code change needed.
9. ~~**Replace commented-out `hack_activity_improved` line in `ImproveCluster_1::mutate`**~~
   — **Done** (explanatory comment added).
10. **Resolve wire-index vs channel-ID ambiguity** noted at `improvecluster_1.cxx`
    (§6.5). Required for correctness on wrapped-wire or non-trivial channel-map detectors.

### P3 — Efficiency

11. ~~**Factor out ~70 duplicated lines** in `CreateSteinerGraph::visit` into a helper~~
    — **Done** (`process_cluster_steiner` lambda, `CreateSteinerGraph.cxx`).
12. ~~**Add `Cluster::wpids_blob_set()`**~~ — **Done** (`Facade_Cluster.h/cxx`; all
    four call sites updated).
13. ~~**Consider removing or implementing `SteinerFunctions::improve_grapher*` stubs**~~
    — **Done**: stubs retained as documentation with explanatory comments replacing the
    `raise<LogicError>` bodies (`SteinerFunctions.cxx/h`).
14. ~~**Address `tick_span` derivation from last blob** (Bug #3)~~
    — **Done**: now uses `m_grouping->get_nticks_per_slice().at(apa).at(face)`
    (`improvecluster_1.cxx`).
15. ~~**Refactor triplicated U/V/W blocks** in `get_activity_improved` into a plane loop~~
    — **Done** (`improvecluster_1.cxx`: steps 2, 3, and 4 each collapsed into a
    single `for (int pl = 0; pl < 3; ++pl)` loop).

### P4 — Testing

16. ~~**Add an integration config that wires `ImproveCluster_2` into `CreateSteinerGraph`**
    (§6.8)~~ — **Done**: `clus/test/test-porting/pdhd/clus.jsonnet` now includes
    `cm.steiner(retiler=improve_cluster_2)` with 4-APA × 2-face samplers and outputs
    `steiner_pc` via `bee_points_sets`.
17. **Bit-reproducibility check:** Run `ImproveCluster_2` twice on the same input;
    confirm blob counts and `steiner_graph` vertex/edge counts are identical after
    the determinism fixes above.

---

## 9. Verification Plan

1. **Build:** `ninja -C build clus` — must compile without errors after each fix batch.
2. **P0 bug checks:** Run on a MicroBooNE-equivalent single-face sample; confirm blob
   counts and Steiner graph topology are unchanged after Bug #1 fix (multi-face path was
   inactive; adding the filter must be a no-op for single-face clusters).
3. **Determinism check:** Run `ImproveCluster_2` twice on the same input file; diff the
   resulting `steiner_graph` edge list. Expect zero differences after §5 fixes.
4. **Multi-face smoke test:** Run on a protoDUNE-HD 4-face sample; confirm no crashes,
   no empty steiner graphs where blobs exist, and that `time_blob_map` entries across
   multiple faces are populated.
5. **Prototype parity:** On a uboone-equivalent sample, compare `new_cluster` blob
   count and time-slice coverage map between the prototype's `Improve_PR3DCluster_2`
   output and the toolkit's `ImproveCluster_2::mutate` output for the same input cluster.
   Allow small differences due to the intentional charge-sentinel change (§2.4); blob
   counts should agree within 1–2%.
