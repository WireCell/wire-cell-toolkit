# Flash-time (T0) grouped clustering — design note

**Status:** implemented. This note describes how the `cm.*` clustering visitors group clusters by
matched-flash time before merging, so that combined-stage clustering happens *within* a flash group
instead of across the whole event. The flash-aware behavior is gated by `use_flash_t0` (default OFF,
so per-APA / production output is bit-identical) and enabled in the SBND combined pipeline.

All `file:line` references are against the tree at the time of writing (branch `apply-pointcloud`).

---

## 1. Motivation

The SBND chain runs clustering in two stages:

- **Per-APA** clustering — one drift volume at a time, in local (uncorrected) `x`.
- **Combined / all-APA** clustering — after `PointTreeMerging`, a second `MultiAlgBlobClustering`
  runs on the global grouping. Its **first** step is `cm.switch_scope()`
  (`ClusteringSwitchScope`), which applies the per-cluster T0 correction and switches the active
  coordinate scope to the T0-corrected `x_t0cor` (`clustering_switch_scope.cxx`, output scope
  `{"3d", {"x_t0cor","y","z"}}` in `PCTransforms.cxx:145-149`).

In the combined stage every cluster's `x` has been shifted by *its own* matched-flash time
`x_corr = x_raw − dirx · cluster_t0 · drift_speed` (`PCTransforms.cxx:58-64,83-98`). Two clusters
matched to **different** flashes are therefore shifted by **different** amounts, and their corrected
geometry is no longer directly comparable. Merging across flashes in this scope can join tracks that
are not physically coincident.

The desired behavior:

- Divide the clusters into **flash-time groups** — clusters whose matched flash times differ by less
  than a window (default **80 ns**) belong to one group.
- Run the existing clustering logic **within** each group only.
- When **no** flash information is present, behave **exactly as today** — per-APA, no-flash output
  must stay bit-identical.

Functions in scope: `extend` (and `extend_loop`, which wraps it), `regular` (×2), `parallel_prolong`,
`close`, and `examine_bundles`. The first five share one merge pattern; `examine_bundles` is a special
case described in §6.

---

## 2. Where T0 lives, and the correct grouping key

QL charge–light matching stamps each matched cluster with two scalars
(`match/src/QLMatching.cxx:798,846`):

| scalar | value | accessor |
|---|---|---|
| `cluster_t0` | flash time in ns (`flash->get_time()*units::ns`) | `Cluster::get_cluster_t0()` |
| `flash` | `flash->get_flash_id()` (per-APA flash ident) | `Cluster::get_flash()` → index into the grouping `flash` PC |
| `matched_flash_gid` | `anode_ident*1e6 + global_flash_idx` — **globally unique** | `get_scalar("matched_flash_gid", …)` |

Unmatched clusters are stamped `cluster_t0 = -1e12` and `flash = -1`
(`match/src/QLMatching.cxx:315-316`). Note `get_cluster_t0()` returns `0` when the scalar is *absent*
(`clus/src/Facade_Cluster.cxx:190`) — so "no flash" has two possible t0 sentinels (`-1e12` and `0`).
**Use the integer `flash` scalar (`>= 0`) as the matched/unmatched discriminator**, and use
`cluster_t0` only to group the matched ones.

**Group on `cluster_t0` (time), never on the `flash` index.** The `flash` scalar is
`get_flash_id()` = `flash.ident()` (`match/src/Opflash.cxx:41`), which **can repeat across APAs** —
the code comment at `match/src/QLMatching.cxx:851` says exactly this, and it is the reason a separate
globally-unique `matched_flash_gid` exists (`QLMatching.cxx:807-808`). The same physical flash split
across APA0/APA1 carries colliding `flash` indices but ~equal time, so a time window of 80 ns is what
correctly merges cross-APA coincident clusters. Grouping by the `flash` index would wrongly separate
them.

For preserving flash provenance through a merge (§6), keep **`matched_flash_gid`** (globally unique),
not the `flash` index.

---

## 3. Unmatched clusters are already filtered out before these functions run

In the combined stage `switch_scope` runs first. For an unmatched cluster (`cluster_t0 = -1e12`):

1. The T0 x-shift is enormous, so every blob's corrected point lands outside the active volume
   (`PCTransforms.cxx:74-81,123-143`, `IDetectorVolumes::contained_by`).
2. `add_corrected_points` returns an all-fail per-blob filter (`Facade_Cluster.cxx:214-243`).
3. `clustering_switch_scope` separates the cluster into `id==0` (all blobs failed) and `id==1`
   (some blob passed) and sets `set_scope_filter(scope, id == 1)`
   (`clus/src/clustering_switch_scope.cxx:83-89`). The unmatched cluster lands in `id==0` with
   `scope_filter = false`.

Every downstream `cm.*` function skips clusters that fail `get_scope_filter(scope)` (e.g.
`clustering_regular.cxx:492`). **So in the combined stage these functions effectively only see
matched clusters with a valid `cluster_t0`** — the mixed matched/unmatched case largely does not
arise. The design below still adds a defensive rule: a cluster without a valid flash gets its own
singleton group and is never merged while flash-aware mode is active (see §5).

### T0 survives the `switch_scope` split — observed

`switch_scope` is a **divide** (it destroys each input cluster and recreates one or two sub-clusters
via `Grouping::separate()`), and it runs *after* QL matching, so the matched-flash scalars are
already on the clusters when it splits them. The scalars (`cluster_t0`, `flash`, `matched_flash_gid`)
live in the `cluster_scalar` local PC, which `Cluster::from()` copies to every child inside
`separate()` (`Facade_Grouping.cxx:163`, `Facade_Cluster.cxx:164-180`) — so both halves inherit the
parent's flash association.

This was verified empirically (not just by code reading): a temporary probe logging
`parent → sub` `(cluster_t0, flash)` across the `separate()` in `clustering_switch_scope.cxx`, on
SBND mc events idx 1–5, showed every sub-cluster carrying its parent's exact `cluster_t0` and `flash`
(47/47 single-sub cases on evt2 matched with zero mismatches; genuine **two-way** splits on evts
1/3/4/5 had *both* `id==0` and `id==1` halves carry the identical non-default values, e.g. evt4
parent ident=4 `t0=-1350658.12 flash=5` → both subs identical). The probe was removed after the run.

---

## 4. The common merge pattern (the five pairwise functions)

`clustering_extend`, `clustering_regular`, `clustering_parallel_prolong`, and `clustering_close` all
share the identical shape (line numbers from `clustering_regular.cxx`):

1. `auto live_clusters = live_grouping.children();`  (`:472`)
2. `sort_clusters(live_clusters);`  (`:486`)
3. O(N²) pairwise loops; for each accepted pair `boost::add_edge(...)` into a
   `cluster_connectivity_graph_t`, guarded per-cluster by `get_scope_filter(scope)`  (`:490-504`)
4. `merge_clusters(g, live_grouping);`  (`:507`) — connected-components → one merged cluster per
   component (`ClusteringFuncs.cxx:48-120`)

`ClusteringExtendLoop::visit` simply calls `clustering_extend(...)` 3–4× in a loop
(`clustering_extend.cxx:91-105`). **A guard placed inside `clustering_extend`'s pairwise loops covers
`extend_loop` automatically** — no separate change is needed there.

`examine_bundles` does **not** follow this pattern today: it only *labels* per-blob connectivity
(`clustering_examine_bundles.cxx:98,168`), it does not merge. See §6.

---

## 5. Proposed change to the five pairwise functions — one group guard

### Shared helper

Add a helper (in `ClusteringFuncs.{h,cxx}`) that assigns each cluster an integer **t0-group id**:

- Take the matched clusters (`get_scalar("flash",-1) >= 0`), sort by `cluster_t0`, walk the sorted
  list, and start a new group whenever the gap to the previous `cluster_t0` exceeds the window.
- A cluster **without** a valid flash gets a **unique singleton id**, so it can never share a group
  with any other cluster.
- Return a map `Cluster* → group_id` covering *every* cluster in the list.

### Config

Each affected visitor gains two parameters, surfaced through the jsonnet factory:

| param | type | default | meaning |
|---|---|---|---|
| `use_flash_t0` | bool | `false` | master toggle; `false` ⇒ today's behavior exactly |
| `flash_t0_window` | double | `80*units::ns` | grouping window on `cluster_t0` |

**Default OFF preserves production bit-identicality** (project convention: new reco variants must be
jsonnet-toggleable and default off so existing configs are bit-identical). The factory functions in
`cfg/pgrapher/common/clus.jsonnet` (`extend`, `regular`, `parallel_prolong`, `close`, `extend_loop`,
`examine_bundles`) gain the two params; they are enabled in the SBND **combined** `cm_pipeline`
(`cfg/pgrapher/experiment/sbnd/clus.jsonnet`).

### The guard

Compute the group-id map once, right after `live_grouping.children()`. In each inner pairwise test,
add a single guard:

```cpp
if (use_flash_t0 && group[c1] != group[c2]) continue;
```

Place it so it does **not** bypass the `used_clusters` early-termination bookkeeping in `extend` and
`close` (i.e. skip the edge, not the loop accounting).

**The guard must be unconditional on group id** (every cluster has a group, including no-flash
singletons). A weaker "skip only when *both* clusters have a flash and differ" guard leaves a
**bridging hole**: because the final merge is connected-components, a no-flash cluster `N` could edge
`A` (group 1) and `B` (group 2), and `connected_components` would then merge `A–N–B` across groups.
Giving `N` its own singleton group makes every cross-group edge impossible.

Result: cross-group edges are never added, so `merge_clusters` (unchanged) can never merge across
flash groups. Within a group, behavior is identical to today; with `use_flash_t0=false`, the whole
function is identical to today.

---

## 6. `examine_bundles` with T0 — a merge pre-step (as implemented)

Today `clustering_examine_bundles` loops over each cluster *individually*, runs
`Cluster::connected_blobs()` to partition that cluster's blobs, marks the main sub-component `-1`, and
stores `("isolated","perblob")` (`clustering_examine_bundles.cxx:84-170`); it reads no flash/T0 and
merges nothing.

The implemented T0 behavior (when `use_flash_t0` is on): **prepend a merge pre-step** that collapses
each matched-flash-time group of clusters into a single cluster, then run the **existing per-cluster
labeling loop unchanged**. The end state is one cluster per flash group carrying a fresh
`("isolated","perblob")` array with the main sub-component `-1` — exactly the shape
`clustering_isolated` produces, but grouped by flash time instead of geometry.

### Mechanism

1. `assign_flash_t0_groups(live_clusters, flash_t0_window)` (in `ClusteringFuncs.cxx`) assigns each
   cluster a flash-time group id; clusters with no valid flash get unique singleton ids.
2. Add an edge between every pair of in-scope clusters that share a group id, then
   `merge_clusters(g, live_grouping)` → one merged cluster per group. (No `savecc`/`perblob` is written
   by this merge; the labeling loop in step 3 writes the final `perblob`.)
3. Re-fetch `live_grouping.children()` and run the **existing** loop: `connected_blobs()` partitions
   each merged cluster, the longest sub-component is marked `-1`, and `put_pcarray("isolated","perblob")`
   records it.

With `use_flash_t0=false` the pre-step is skipped and the function is byte-identical to today.

### No split, no restore

An earlier draft of this note proposed a merge → separate → restore round-trip. That was only needed
if a **downstream** separator (the old `clustering_isolated` / `ClusteringRecoveringBundle`) re-split
the merged cluster and had to restore each fragment's flash. **The combined SBND pipeline drops those
downstream steps** (its tail is now just `examine_bundles`), so the merged cluster is terminal — there
is no later separation and therefore **no per-fragment flash restore is required**. "Select the main
cluster" is simply the `-1` main-marking inside `examine_bundles`; no separate node and no satellite
deletion (satellites become sub-components of the one merged cluster, recorded in `perblob`).

### Flash bookkeeping is fixed once, in `merge_clusters` (covers all clus call sites)

`Cluster::from` copies a scalar **only if the destination lacks it** (`Facade_Cluster.cxx:159-181`), so
a naive `merge_clusters` keeps the **first-encountered** member's `cluster_t0`/`flash`/`matched_flash_gid`
(arbitrary `std::set` order). This is a problem for **every** merging visitor, not just
`examine_bundles`. The fix lives in the shared primitive `merge_clusters` (`ClusteringFuncs.cxx`):
after each component is merged, if **any** contributing member has a valid `flash` (`>= 0`), set the
merged cluster's `cluster_t0`/`flash`/`matched_flash_gid` to those of the **longest** contributing
member. When no member carries a flash (the entire per-APA stage, where matching has not yet run) the
block is a no-op, so per-APA output stays bit-identical. This single change makes `extend`, `regular`,
`parallel_prolong`, `close`, `isolated`, `connect1`, `examine_bundles`, … all flash-correct.

This matters because `cluster_t0` is read downstream by `connect_graph_ctpc`/`connect_graph_relaxed`,
`NeutrinoStructureExaminer`, `MyFCN`, `FiducialUtils`, `TrackFitting`, `NeutrinoVertexFinder`,
`protect_overclustering` (`grep get_cluster_t0 clus/src/*.cxx`); the representative must be a real
member's t0, never 0.

`Grouping::separate()` needs no change: it already propagates the parent's (single) flash to each split
child via `c->from(*cluster)` (`Facade_Grouping.cxx`), which is correct for the single-flash clusters
that get separated in this pipeline (e.g. by `switch_scope`). A per-fragment restore would only be
needed if a *multi-flash* merged cluster were later re-split — which this pipeline does not do.

---

## 7. Configuration summary

| function | jsonnet (`cfg/pgrapher/common/clus.jsonnet`) | C++ visitor | new params |
|---|---|---|---|
| `extend` | `extend(...)` | `ClusteringExtend` | `use_flash_t0`, `flash_t0_window` |
| `extend_loop` | `extend_loop(...)` | `ClusteringExtendLoop` | (covered via `clustering_extend`) |
| `regular` | `regular(...)` | `ClusteringRegular` | `use_flash_t0`, `flash_t0_window` |
| `parallel_prolong` | `parallel_prolong(...)` | `ClusteringParallelProlong` | `use_flash_t0`, `flash_t0_window` |
| `close` | `close(...)` | `ClusteringClose` | `use_flash_t0`, `flash_t0_window` |
| `examine_bundles` | `examine_bundles(...)` | `ClusteringExamineBundles` | `use_flash_t0`, `flash_t0_window` |

All default `use_flash_t0=false`. Enable in the SBND combined `cm_pipeline`
(`cfg/pgrapher/experiment/sbnd/clus.jsonnet`); per-APA stays off (no flash there today) unless a flash
chain is later wired into the per-APA stage.

---

## 8. Open items / risks

- **Stale per-APA `perblob`** needs no wipe node: combined step 1 `switch_scope` loops every cluster
  (`clustering_switch_scope.cxx:64`) and calls `separate(...,remove=true)` unconditionally (`:83`);
  `NaryTree::separate` processes all non-negative group ids (`NaryTree.h:269`), so even a fully-in-FV
  cluster is destroyed and recreated via `Cluster::from()`, which copies only `cluster_scalar`
  (`Facade_Cluster.cxx:159-181`). The per-APA `("isolated","perblob")` array is therefore dropped at
  step 1 and never reaches the combined `examine_bundles`, which writes a fresh one.
- **Bit-identicality (toggle off):** verify the **per-APA** stage output is byte-for-byte unchanged
  (all five pairwise functions + `merge_clusters` are no-ops without a flash). This is the regression
  gate. The combined stage intentionally changes (its tail is replaced) and is not gated by the flag.
- **Correctness (toggle on):** verify no merged cluster spans a `cluster_t0` range > `flash_t0_window`,
  each flash group collapses to one cluster carrying `("isolated","perblob")` with exactly one main
  sub-component (`-1`), and the merged cluster's `cluster_t0` equals the longest contributing member's.
- **Greedy-chain edge case:** if real flashes were ever spaced < 80 ns apart, greedy chaining could
  link them into one group. In SBND practice matched clusters cluster tightly at discrete flash times,
  so this is not expected; the window is a tolerance for the same flash reconstructed slightly
  differently across APAs, not a way to bridge distinct flashes.
