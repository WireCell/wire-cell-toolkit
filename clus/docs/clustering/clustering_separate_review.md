# ClusteringSeparate / clustering_separate() — Code Review

**Date**: 2026-04-10  
**Reviewed by**: Claude Code  
**Toolkit file**: `clus/src/clustering_separate.cxx`  
**Prototype reference**: `prototype_base/2dtoy/src/ToyClustering_separate.h`  
**Functions in scope**: `clustering_separate()` (driver), `Separate_1()`, `Separate_2()`,
`JudgeSeparateDec_1()`, `JudgeSeparateDec_2()`, `Facade::get_uvwt_range(cluster,b2id,id)`,
`Facade::get_length(cluster,b2id,id)` (`clus/src/Facade_Cluster.cxx:1788-1843`).

---

## 1. Logic fidelity

### 1.1 `sep_clusters[0]` / `sep_clusters[k≥2]` — NOT a bug (verified)

The driver only reads `sep_clusters.at(1)` for recursion, ignoring index 0 and k≥2.  
This looks like missing logic but is **correct**: `Separate_1` already inserts
`clusters_step0[0]`, `cluster2`, and each `saved_cluster[k]` into `live_grouping` via
`grouping->separate()` / `grouping->make_child()` before returning.  Those nodes already
live in the tree.  The driver only needs `sep_clusters.at(1)` to decide whether to recurse.
The prototype did the same conceptually (only `cluster2` was recursed on); the accounting
difference is an ownership model change, not a semantic difference.

### 1.2 `JudgeSeparateDec_1::temp_angle1` formulation (deliberate divergence)

**Prototype**: `asin(num_time_slices * tick_width / length)` — drift-axis extent from
time-slice count.  
**Toolkit**: `asin(|x_earliest - x_latest| / length)` — x-extent of 3D point extremes.

These are numerically different for sparse clusters.  Decision (confirmed by user): keep
the toolkit formulation.  The old `CHECKME` comment has been replaced with an explanation.

### 1.3 PDHD early return in `JudgeSeparateDec_1` (deliberate extension)

```cpp
if (angle1 < 5 && ratio2 < 0.05 && length > 300*units::cm) return false;
```

Not in prototype.  Guards against over-separating long, drift-aligned, very thin tracks
(through-going cosmics in PDHD where the drift volume is much longer than MicroBooNE).
Added a block comment explaining the rationale.

### 1.4 `JudgeSeparateDec_2` global FV envelope (correct by design)

Uses `wpid_all(0)` to query global detector bounds.  This is intentional: the function
detects whether a cluster exits the *cryostat*, not an individual TPC.  For a multi-TPC
detector (PDHD), a cosmic traversing one drift volume but staying within the outer
cryostat envelope should NOT be flagged as "exiting".  Added a detailed comment.

### 1.5 `JudgeSeparateDec_2` z-boundary convention (confirmed correct)

Metadata `FV_zmin = 15 cm` (inner), `FV_zmin_margin = 3 cm` (so outer = 12 cm).  
Matches prototype's hardcoded `zmin=15, zmin_outer=12`.  
Same for `FV_xmin=1, FV_xmin_margin=2 → outer=-1`.  All conventions verified correct.

---

## 2. Bugs fixed

### B1 — `num_outx_points` sign error (HIGH priority)

**Location**: `clustering_separate.cxx` — 6 occurrences in `JudgeSeparateDec_2`
(hy, ly, hz, lz, hx, lx directional blocks).

**Bug**: `pt.x() > det_FV_xmax - det_FV_xmax_margin` — subtracts the margin instead of
adding it.  With `FV_xmax=255 cm` and `FV_xmax_margin=2 cm`, the buggy condition triggers
at `x > 253 cm` instead of the correct `x > 257 cm`.  The `flag_outx` test at line 410
correctly uses `+`, so this was an isolated inconsistency.

**Effect**: `num_outx_points` counter fires ≈2cm too early on the upper-x side, causing
clusters that are just inside the outer boundary to be misclassified as "exiting", biasing
the cosmic separation decision.

**Fix**: Changed all 6 occurrences to `> det_FV_xmax + det_FV_xmax_margin`.

### B2 — `clusters_step0[0]` crash when no path blobs (HIGH priority)

**Location**: `clustering_separate.cxx:1373` (inside `Separate_1`).

**Bug**: The `find(0)` guard was commented out.  If every blob is assigned group id `1`
or `-1` (no blob survives on the path), `operator[]` default-inserts a null `Cluster*` for
key `0`, and the subsequent `->nchildren()` call segfaults.

**Fix**: Restored `if (clusters_step0.find(0) != clusters_step0.end())` and removed the
redundant `->nchildren() > 0` check (guaranteed non-empty by `grouping->separate`
semantics — only existing groups are inserted).

### B3 — Missing `test_wpid.apa() != -1` guard in point-classification loop (MEDIUM)

**Location**: `clustering_separate.cxx:1210` (inside `Separate_1`).

**Bug**: The path-building loop (line ~1200) guards `test_wpid.apa() != -1` before
`af_temp_cloud.at(...)`.  The cluster-point classification loop (line ~1210) had no such
guard; a point not mapped to any APA would crash with `at(-1)`.

**Fix**: Added `if (test_wpid.apa() == -1) continue;` at the top of the classification loop.

### B4 — Dead scope_transform variable (LOW)

**Location**: `clustering_separate.cxx:1353` (inside `Separate_1`).

**Bug**: `auto scope_transform = cluster->get_scope_transform(scope);` was captured before
`grouping->separate(cluster, ...)` destroys `cluster`, making the value unusable.  It was
never referenced again (the scope_transform was re-fetched on line ~1402 from
`clusters_step0[0]`).

**Fix**: Removed the dead assignment.

### B5 — Empty `cluster2` has no scope_transform (MEDIUM)

**Location**: `clustering_separate.cxx:1519-1522` (inside `Separate_1`).

**Bug**: `cluster2.set_scope_transform(...)` was only called when `to_be_merged_clusters`
is non-empty.  If everything ended up in `saved_clusters`, `cluster2` was pushed to
`final_clusters` with no scope_transform, causing downstream scope lookups to fail.

**Fix**: Added `else cluster2.set_scope_transform(scope, clusters_step0[0]->get_scope_transform(scope));`
as a fallback.

### B6 — Variable shadowing in `get_uvwt_range(cluster, b2id, id)` (LOW)

**Location**: `clus/src/Facade_Cluster.cxx:1798-1809`.

**Bug**: Outer loop `for (size_t i = ...)` iterates `b2id`.  Four inner loops
`for (int i = blob->xxx_min(); ...)` shadow it.  Currently correct by luck (outer `i` is
fully captured before inner loops start) but fragile.

**Fix**: Renamed inner loop counters to `wi` (wire) and `ti` (tick).

### B7 — `hx_points` loop reads `lx_points.at(j)` (MEDIUM)

**Location**: `clustering_separate.cxx:676-693` (inside `JudgeSeparateDec_2`).

**Bug**: Prototype bug faithfully ported.  The `hx_points` independent-surface
classification block reads `lx_points.at(j)` instead of `hx_points.at(j)`.  When
`hx_points.size() > lx_points.size()`, this is an out-of-bounds access.

**Fix**: Replaced all `lx_points.at(j)` references in the `hx_points` loop with
`hx_points.at(j)`.  Note: this is a deliberate divergence from the prototype.  It
changes which surface code is inserted for the highest-x cluster extremum.

---

## 3. Efficiency improvements

### C1 — Cache `get_length()`, `get_pca()`, `JudgeSeparateDec_1/2` per iteration

**Location**: outer loop in `clustering_separate()`, lines ~113-295.

**Before**: `cluster->get_length()` called 6+ times, `cluster->get_pca()` called 4-5
times, `JudgeSeparateDec_1` called twice on the same unmodified cluster.

**Fix**: 
- `const double cluster_length = cluster->get_length();` once at top of iteration.
- `const auto& pca = cluster->get_pca(); const double pca_ratio1 = pca.values[1]/pca.values[0];` once inside the `flag_top` block.
- `bool flag_dec1 = JudgeSeparateDec_1(cluster, ...);` cached once; reused in the
  `if (flag_proceed)` check below.

### C2 — Fuse 6 directional sweeps in `JudgeSeparateDec_2` (deferred)

The hy/ly/hz/lz/hx/lx extremum-classification blocks each make a full pass over
`boundary_points`.  They can be fused into one O(N) pass.  Not done in this review
(invasive restructuring); left as a future optimization.

### C3 — Hull caching (deferred)

`cluster->get_hull()` is recomputed on each call to `JudgeSeparateDec_2`.  If Cluster
grows a hull cache, repeated calls in the recursion (4 levels) would be free.

### C4 — `Simple3DPointCloud` lazy kd-tree (confirmed correct, no fix needed)

The `kd()` method in `Facade_Util.h:143` initializes `m_kd` lazily on first call.  The
prototype explicitly called `build_kdtree_index()` before pairwise distance queries;
the toolkit builds it automatically.  No regression.

---

## 4. Determinism

### D1 — Sort tiebreak using `cluster->ident()` instead of pointer address

**Location**: `clustering_separate.cxx:70-74`.

**Before**: `return c1 < c2;` — raw pointer comparison, non-deterministic across runs.

**Fix**: `return c1->ident() < c2->ident();`

`ident()` is read from the cluster's scalar PC array (set via `Grouping::enumerate_idents()`
or inherited from parent during `separate()`).  It provides a stable, run-independent
tiebreak for equal-length clusters.

Note: idents may not be globally unique (clusters produced from the same parent by
`separate()` share the same ident).  In the highly unlikely case of same-ident same-length
clusters, the ordering remains undefined but is at least better than a pointer.

### D2 — `BlobSet` uses `BlobLess` comparator (confirmed deterministic)

`BlobSet = std::set<const Blob*, BlobLess>` where `BlobLess` delegates to `blob_less(a,b)`
which compares on stable blob properties.  Not pointer-ordered.  No fix needed.

### D3 — Pointer-keyed maps used only for lookups (documented)

`std::map<const Blob*, int> mcell_np_map / mcell_np_map1 / mcell_index_map` appear in
`Separate_1` and `Separate_2`.  All are accessed only via `operator[]` or `find()` — never
iterated — so pointer-order non-determinism is irrelevant.  Added `// iteration-order safe:
lookups only` comments.

---

## 5. Multi-APA / multi-face support

### M1 — `JudgeSeparateDec_2` global FV bounds (correct by design — see §1.4)

### M2 — `Separate_2` per-(apa,face) time-slice grouping (correct)

`af_time_slices[apa][face]` partitions blobs by face before building adjacency edges.
`get_nticks_per_slice().at(apa).at(face)` used for per-TPC tick granularity.  Correct.

### M3 — `Separate_1` path-cloud widens for multi-face clusters (E1 fix)

**Before**: `af_temp_cloud[apa][face]` stores path points per face; each cluster point is
tested only against its own face's 2D projection.  A path crossing an APA boundary leaves
cluster points near the boundary unmatched (they see an empty projection in their face).

**Fix (E1)**: Added `const bool is_multi_face = (cluster->wpids().size() > 1)`.  When
true, the 2D distance for each cluster point is computed as the minimum across ALL face
projections, using a `min_2d_dis` lambda.  Dead-channel lookups remain face-local.

The single-face case (most common) is unchanged.  The multi-face path takes O(nfaces ×
npoints) instead of O(npoints) but is only active when needed.

### M4 — `wire_indices()[plane][j]` is face-scoped (confirmed)

`wire_indices()` returns a scoped view keyed by `m_scope_wire_index`, and
`get_dead_winds(apa, face, plane)` uses the same wire-index convention for that face.
Lookups at `dead_u_index.find(winds[0][j])` are consistent.  No fix needed.

### M5 — `get_length(cluster, b2id, id)` sums across faces (by design)

When a subcluster spans multiple (apa,face) pairs, `get_length` sums partial lengths per
wpid.  This matches `Cluster::get_length()` semantics.

---

## Summary of changes

| # | File | Change | Priority |
|---|------|---------|----------|
| B1 | clustering_separate.cxx | `num_outx_points` sign fix (6 sites) | HIGH |
| B2 | clustering_separate.cxx | Restore `clusters_step0.find(0)` guard | HIGH |
| B3 | clustering_separate.cxx | Add `test_wpid.apa() != -1` guard in point loop | MEDIUM |
| B4 | clustering_separate.cxx | Remove dead `scope_transform` variable | LOW |
| B5 | clustering_separate.cxx | Fix `cluster2` scope_transform fallback | MEDIUM |
| B6 | Facade_Cluster.cxx | Rename inner loop vars `i` → `wi`/`ti` | LOW |
| B7 | clustering_separate.cxx | Fix `hx_points` loop prototype bug | MEDIUM |
| C1 | clustering_separate.cxx | Cache `get_length`/`get_pca`/`JudgeSeparateDec_1` | PERF |
| D1 | clustering_separate.cxx | Deterministic sort tiebreak via `ident()` | MEDIUM |
| D3 | clustering_separate.cxx | Add `// lookups only` comments on pointer maps | LOW |
| E1 | clustering_separate.cxx | Multi-face min-distance widening in `Separate_1` | MEDIUM |
| E2 | clustering_separate.cxx | Clarify global FV comment, remove dead commented block | LOW |

---

## Remaining known limitations

- **C2/C3**: 6-sweep fusion and hull caching not implemented (deferred).
- **D1 residual**: Equal-ident same-length clusters (very rare) still have undefined order.
- **E1 caveat**: The `min_2d_dis` helper queries all face 2D clouds using their own wire
  angles.  For a point in face A, comparing its 3D position against face B's 2D projection
  (which uses face B's wire angles) is only physically meaningful if the path genuinely
  passes nearby in 3D — the 2D distance is a projection, not the true 3D distance.  This
  approach is conservative (will sometimes classify a cluster point as "near-path" when
  the path is actually far), which is acceptable for this use case.
- **BDT scoring**: Not in scope for this review.
