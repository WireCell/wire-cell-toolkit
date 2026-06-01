# Flash-time (T0) grouped clustering — design note

**Status:** design / proposal. No code is changed by this document. It describes how to make the
`cm.*` clustering visitors group clusters by matched-flash time before merging, so that clustering
happens *within* a flash group instead of across the whole event.

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

## 6. `examine_bundles` with T0 — a new merge path

Today `clustering_examine_bundles` only labels per-blob connected components via
`Cluster::connected_blobs()` and stores `("isolated","perblob")` (`:98,168`); it reads no flash/T0
and merges nothing.

The new behavior (when `use_flash_t0` is on and flashes are present): form **cluster groups by T0
coincidence** and **merge each group into one cluster**.

### Mechanism (borrowed from `isolated`, criterion replaced)

`clustering_isolated.cxx` is the template for the *mechanism*: build candidate pairs → union-find →
add edges into a graph → `merge_clusters(g, live_grouping, "isolated")` (`:478`). But its grouping
**criterion is purely geometric** (closest-distance and PCA cuts, `:228,251,314-359`). For T0 mode,
**borrow the union-find + `merge_clusters` machinery, and replace the geometric criterion with T0
coincidence**: add an edge between every pair of clusters that share a t0-group id, then call

```cpp
merge_clusters(g, live_grouping, "isolated", "perblob");
```

so the merged cluster records, per blob, which source cluster it came from (`ClusteringFuncs.cxx:74,
101-110`). Keep the legacy per-blob *labeling* path when `use_flash_t0` is off.

### The flash-bookkeeping conflict, and the merge → separate round-trip

The requirement is stronger than just "pick a good representative". A merged cluster is not the end
of the line — it can be **separated again** downstream (the `("isolated","perblob")` label written by
`examine_bundles` is consumed by `ClusteringRecoveringBundle`, which calls
`grouping.separate(cluster, cc_vec)` to split a beam-flash cluster back into its components,
`ClusteringRecoveringBundle.cxx:106,148-149,163`). The flash information of each original piece must
**survive the merge and be restored on the later separation** — it must not be lost when the clusters
are combined.

Two facts make naive merging lossy:

1. **Merge keeps only one flash.** `merge_clusters` iterates a component's members in
   `std::set<descriptor>` order, calling `fresh_cluster.from(*live)` (`ClusteringFuncs.cxx:87-94`).
   `Cluster::from` copies a scalar **only if the destination does not already have it**
   (`Facade_Cluster.cxx:173-179`, the `if (arr && !arr1)` guard). So the merged cluster silently keeps
   the **first-encountered** member's `cluster_t0` / `flash` / `matched_flash_gid` — "first" is
   arbitrary `std::set` order, **not** the main/longest cluster. The other members' flash is dropped.
2. **Separate re-stamps every child with the parent's flash.** `Grouping::separate` calls
   `c->from(*cluster)` on every split child (`Facade_Grouping.cxx:160-163`). So after a separation each
   piece *inherits the merged cluster's single flash* and its own original flash is gone — unless it
   was stored somewhere blob-local and explicitly restored.

This matters beyond the merged cluster itself: `cluster_t0` is read downstream by
`connect_graph_ctpc` / `connect_graph_relaxed`, `NeutrinoStructureExaminer`, `MyFCN`, `FiducialUtils`,
`TrackFitting`, `NeutrinoVertexFinder`, `protect_overclustering`
(`grep get_cluster_t0 clus/src/*.cxx`).

**Recommended solution — save the flash per blob, alongside `perblob`; restore on separate.**

This reuses the established `("isolated","perblob")` idiom (`examine_bundles` writes it, `retile` and
`RecoveringBundle` read it) and adds a *parallel* per-blob array carrying flash provenance:

1. **At merge time** — when `merge_clusters` (savecc mode) walks the component members and fills the
   per-blob component-id array `("isolated","perblob")` (`ClusteringFuncs.cxx:101-110`), also fill a
   **parallel per-blob array `("flash_gid","perblob")`** (optionally `("flash_t0","perblob")`) with
   each source cluster's `matched_flash_gid` (and `cluster_t0`), captured from `*live` *before* it is
   destroyed (`:106`). Store the **gid** (globally unique), not the `flash` index (which can repeat
   across APAs, §2). Both arrays are blob-indexed in the same order, so blob `j`'s component and its
   original flash stay aligned even if connectivity is later recomputed.
2. **Representative on the merged ("combined") cluster** — set its `cluster_t0` / `flash` /
   `matched_flash_gid` to the **main cluster's** flash, where the main cluster is the component marked
   `-1` in the `perblob` array (the max-overlap / longest component, exactly the `-1` convention
   `examine_bundles` and `RecoveringBundle` already compute, `clustering_examine_bundles.cxx:103-166`,
   `ClusteringRecoveringBundle.cxx:112-146`). This is what the user asked for: "the flash information
   should be the main cluster's flash information when we created them." Set it explicitly after the
   merge so it does not depend on `from()`'s arbitrary first-wins order.
3. **On separate** — after `grouping.separate(cluster, cc_vec)`, each split child currently inherits
   the merged flash via `from()` (fact 2 above). Add a restore step (the natural home is
   `ClusteringRecoveringBundle::process_cluster`, right after the `separate` at
   `ClusteringRecoveringBundle.cxx:163-171`, where it already re-idents and flags the pieces): for each
   returned sub-cluster, look up the `("flash_gid","perblob")` value of the blobs that went into it and
   **override** its `cluster_t0` / `flash` / `matched_flash_gid` with that original flash. The main
   piece keeps the main flash; every other piece is restored to its own.

This is lossless: the merged cluster reports the main flash, and any later separation faithfully
restores each fragment's original flash from the per-blob record.

**Alternatives (noted, not recommended):**

- Blob-local scalar instead of a cluster-level `perblob` array: store the flash gid in each blob's own
  local PC so it travels with the blob through *any* `separate()` automatically (no read-before-split
  needed). More general, but departs from the existing `perblob` idiom; revisit only if a separation
  path that bypasses the `perblob` array appears.
- No tolerance: only merge clusters with *identical* `matched_flash_gid` (no 80 ns window) — sidesteps
  the conflict but defeats the cross-APA coincidence goal (§2).

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

- **Representative flash on the merged cluster** = the main component's (the `-1` blob group). Set it
  explicitly after the merge; must never fall back to 0.
- **Round-trip restore:** verify that merge → `RecoveringBundle` separate restores each fragment's
  original flash from `("flash_gid","perblob")` (i.e. a piece does **not** keep the merged/main flash
  it inherited via `from()`).
- **Bit-identicality (toggle off):** verify the combined-stage output is byte-for-byte unchanged when
  `use_flash_t0=false`, on a production SBND config.
- **Correctness (toggle on):** verify no merged cluster spans two `matched_flash_gid` groups, and that
  the parallel per-blob `("flash_gid","perblob")` array is populated after `examine_bundles` and stays
  blob-index-aligned with `("isolated","perblob")`.
- **Greedy-chain edge case:** if real flashes were ever spaced < 80 ns apart, greedy chaining could
  link them into one group. In SBND practice matched clusters cluster tightly at discrete flash times,
  so this is not expected; the window is a tolerance for the same flash reconstructed slightly
  differently across APAs, not a way to bridge distinct flashes.
