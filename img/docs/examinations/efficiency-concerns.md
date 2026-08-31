# Efficiency Concerns

## Algorithmic Complexity

### 1. InSliceDeghosting: O(n_2view * n_3view * n_channels) per slice

**File**: `src/InSliceDeghosting.cxx:362-516`

The local deghosting has three nested loop levels per time slice:
- Outer: iterate over 2-view blobs
- Middle: iterate over 3-view blobs (for adjacency check, lines 390-399)
  or over measures and their connected blobs (lines 481-511)
- Inner: per-plane channel iteration for overlap/adjacency

For slices with many blobs (e.g., busy events), this becomes:
`O(n_2view * n_3view)` for adjacency checks +
`O(n_2view * n_meas * n_blobs_per_meas)` for overlap scoring.

**Status**: NOT FIXED. The algorithm inherently requires pairwise comparisons
between 2-view and 3-view blobs. Reducing this would require fundamentally
different spatial indexing (e.g., grid-based lookup).

### 2. ProjectionDeghosting: O(n_clusters^2) per plane

**File**: `src/ProjectionDeghosting.cxx`

For each wire plane, each pair of 2D clusters is compared via `judge_coverage()`.
With `n` clusters per plane, this is O(n^2) coverage checks, each involving
sparse matrix operations.

**Status**: NOT FIXED. Pairwise coverage comparison is fundamental to the algorithm.

### 3. CMMModifier: O(n^3) nested loops

**File**: `src/CMMModifier.cxx:189-252`

Triple nested loop structure for channel mask processing:
- Outer: channels in mask map
- Middle: adjacent channel search (ch +/- 1, +/- 2)
- Inner: time range alignment

For large detectors with many channels, this compounds.

**Status**: NOT FIXED. Would require redesigning the adjacency search strategy.

### 4. BlobDepoFill: O(n_slices * n_depos * n_wires * n_blobs)

**File**: `src/BlobDepoFill.cxx:293-366`

For each slice, each depo, each wire, each blob -- check intersection.
This is effectively O(n^4) in the worst case. Could be accelerated with
spatial indexing (e.g., R-tree) to avoid brute-force intersection checks.

**Status**: NOT FIXED. Requires spatial indexing data structure (major refactor).

### 5. FrameQualityTagging: Unbounded search loops

**File**: `src/FrameQualityTagging.cxx:319-363`

Per-channel time expansion phase: while loops search for signal edges
without bounded iteration. For noisy data with irregular patterns, these
could iterate many times per channel.

**Status**: NOT FIXED. Would require understanding search intent and adding bounds.

---

## Redundant Computations

### 6. Projection2D: Multiple sparse matrix passes

**File**: `src/Projection2D.cxx:469-480`

The original code called `loop_count()` and `loop_sum()` 7 separate times on
3 sparse matrices, performing 7 full matrix iterations when 3 suffice.

**Status**: FIXED. Replaced 7 separate `loop_count`/`loop_sum` calls in
`judge_coverage_alt()` with 3 single-pass loops (one per matrix: ref, tar,
inter_proj). Each pass computes all needed statistics (live count, dead count,
charge sum) simultaneously. Reduces sparse matrix iterations from 7 to 3.

### 7. FrameQualityTagging: Repeated rebinning

**File**: `src/FrameQualityTagging.cxx:420-434`

The global fire rate calculation re-performs the same rebinning and content
summation that was already done during the per-plane analysis phase.

**Status**: NOT FIXED. Would require caching rebinned results across phases,
significant restructuring of the function.

### 8. BlobSetReframer: Channel cache misses trigger full rebuild

**File**: `src/BlobSetReframer.cxx:86-101`

Channel lookup cache was populated on demand. Each cache miss re-iterated all
planes to find the channel. Multiple misses caused repeated full plane scans.

**Status**: FIXED. Replaced per-miss cache repopulation with `ensure_cache()`
that pre-populates the entire chid->IChannel map once per face on first use.
Subsequent lookups are O(1) with no risk of redundant plane iteration.

### 9. ChargeSolving: dump_cg called unconditionally

**File**: `src/ChargeSolving.cxx:268`

`dump_cg()` iterated the entire cluster graph counting vertices and extracting
blob values, even when debug logging was disabled.

**Status**: FIXED (in bug fix pass). Guarded with `log->level() <= spdlog::level::debug`.

---

## Unnecessary Copies and Allocations

### 10. CSGraph: Expensive matrix inverse + Cholesky

**File**: `src/CSGraph.cxx:160`

The original code computed the full matrix inverse `mcov.inverse()` (O(n^3))
then performed Cholesky decomposition of the result (another O(n^3)).

**Status**: FIXED. Replaced with direct Cholesky of `mcov` followed by triangular
solve (`L.solve()`). The whitening transform U = L^{-1} is computed implicitly
via `L.solve(x)` instead of forming U explicitly. This:
- Eliminates the O(n^3) explicit inverse
- Uses more numerically stable triangular solve
- Produces mathematically equivalent whitening (U^T * U = mcov^{-1})

### 11. Projection2D: Dense conversion for file I/O

**File**: `src/Projection2D.cxx:297`

```cpp
Eigen::MatrixXf dense_m(proj2d.m_proj);
```

Converts entire sparse matrix to dense format for file writing.
For large (nchan x nslice) matrices this is O(n*m) memory even when
the matrix is highly sparse. Should write in sparse format.

**Status**: NOT FIXED. The `write()` function is a diagnostic/debug utility
(not a hot path), and `WireCell::Stream::write` expects dense array format.
Changing this would require Stream API modifications.

### 12. GlobalGeomClustering: Debug-only graph traversals

**File**: `src/GlobalGeomClustering.cxx`

Debug dump and connected_components calls ran unconditionally, traversing the
full graph even at production log levels.

**Status**: FIXED. Guarded all debug-only code (`dump_cg`, `connected_components`
for cluster counting, debug log messages) behind `log->level() <= spdlog::level::debug`.
This eliminates two full `connected_components` computations and four `dump_cg`
graph traversals in production runs.

### 13. BlobGrouping: Per-plane subgraph construction

**File**: `src/BlobGrouping.cxx:55-109`

Builds 3 separate subgraphs (one per plane) per slice by iterating
all adjacent vertices multiple times. Could collect edges in a single
pass and build subgraphs afterward.

**Status**: NOT FIXED. The current code already caches unique channels per plane
(`uniq_chans`) and the traversal pattern (b->w->c) is inherent to the graph
structure. Savings would be marginal.

### 14. MaskSlice: Linear search for plane membership

**File**: `src/MaskSlice.cxx:303-304`

Used `std::find` on vectors for plane membership checks in the per-trace hot loop.
With 3 planes the vectors are small, but the pattern is inefficient.

**Status**: FIXED. Precomputed `std::unordered_set<int>` from `m_active_planes`
at the start of `slice()`. Replaced `std::find()` on vector with O(1) set lookup.

---

## Memory and Allocation Patterns

### 15. set_intersection creates temporary vectors

**File**: `src/InSliceDeghosting.cxx:259-265`

`std::set_intersection` allocated a temporary `std::vector<int>` just to count
the common elements. This function is called in a tight inner loop (per blob pair).

**Status**: FIXED. Replaced with a manual merge-style counting loop that walks
both sorted sets simultaneously, counting matches without any allocation. The
algorithm is the same O(n+m) complexity as `set_intersection` but avoids the
heap allocation and copy overhead.

### 16. Unordered containers without reserve

Multiple files use `std::unordered_map` and `std::unordered_set` in loops
without calling `reserve()`. For known-size containers (e.g., number of
blobs, number of channels), pre-reserving avoids rehashing.

**Status**: NOT FIXED. The impact is typically small since these containers
grow incrementally and most instances have small sizes. The overhead of
rehashing is amortized O(1) per insertion.

---

## Summary

| # | Issue | Status | Impact |
|---|-------|--------|--------|
| 6 | Single-pass sparse matrix stats | FIXED | Reduces 7 matrix iterations to 3 |
| 8 | BlobSetReframer cache pre-population | FIXED | Eliminates redundant plane scans |
| 9 | ChargeSolving dump_cg guard | FIXED | Eliminates graph traversal at production log levels |
| 10 | CSGraph direct Cholesky solve | FIXED | Eliminates O(n^3) explicit matrix inverse |
| 12 | GlobalGeomClustering debug guards | FIXED | Eliminates 2 connected_components + 4 dump_cg in production |
| 14 | MaskSlice set-based plane lookup | FIXED | O(1) vs O(n) per-trace plane check |
| 15 | Count-only wire overlap | FIXED | Eliminates vector allocation per blob pair |
| 1-5,7 | Algorithmic complexity | NOT FIXED | Require fundamental algorithm redesign |
| 11,13,16 | Minor allocations | NOT FIXED | Low impact or requires API changes |
