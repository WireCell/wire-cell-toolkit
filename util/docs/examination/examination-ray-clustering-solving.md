# Ray Clustering & Solving - Code Examination

## Overview

The Ray Clustering and Ray Solving subsystems provide blob association and charge solving in the Wire-Cell Toolkit's RayGrid framework. Ray Clustering projects blobs onto 1-D grid arrays per layer, detects mutual overlap across multiple layers via recursive descent, tests containment relationships, and associates blobs between adjacent time slices. Ray Solving takes the resulting blob-measurement association graph, discovers connected components via Boost `connected_components`, and solves the linear system `m = G * s` (measurements = geometry-matrix times sources) using LASSO regularization to extract per-blob charge estimates.

## Files Examined

- `util/inc/WireCellUtil/RayClustering.h`
- `util/src/RayClustering.cxx`
- `util/inc/WireCellUtil/RaySolving.h`
- `util/src/RaySolving.cxx`

## Algorithm Description

### Blob Projection to 1-D Arrays

`projection()` (RayClustering.cxx:6-18) builds a sparse-array representation of blob occupancy along one layer. For each blob, it reads the blob's strip bounds in the given layer (`bounds.first` to `bounds.second`, half-open) and inserts the blob reference into every grid-index bucket in that range. The result (`blobproj_t`) is a `vector<blobvec_t>` indexed by grid index, dynamically resized to accommodate the maximum grid index seen.

### Blob Selection via Range Query

`select()` (RayClustering.cxx:21-34) collects all unique blobs that occupy any grid index within a half-open range `[first, second)`. It clamps the upper bound to the projection size and accumulates blobs into an `unordered_set` (using `blobref_hash`) for deduplication, then returns the unique set as a vector.

### Recursive Overlap Detection Across Layers

`overlap()` (RayClustering.cxx:36-57) recursively finds blobs that overlap a reference blob across multiple layers. Starting at the highest layer, it selects candidate blobs whose strip range (plus tolerance) overlaps the reference blob's strip in that layer, then re-projects those candidates onto the next-lower layer and recurses. Recursion terminates at layer 2 (the lowest wire-plane layer index) or when no candidates remain. This effectively performs a multi-layer intersection test.

### Surrounding Containment Test

`surrounding()` (RayClustering.cxx:69-92) checks if one blob's strips entirely contain another blob's strips across all layers, or vice versa. For each layer it checks whether a's bounds contain b's bounds and vice versa. Returns true only if one direction holds for every layer simultaneously.

### Blob Association

`associate()` (RayClustering.cxx:94-108) is the high-level API that pairs blobs from two time slices. It projects the second set's blobs onto the highest layer, then for each blob in the first set, finds overlapping blobs via the recursive `overlap()` call and invokes a user-supplied callback for each pair.

### Graph-Based Grouping (Boost connected_components)

`Grouping` (RaySolving.h:21-53, RaySolving.cxx:6-55) builds a tripartite graph of measurement nodes ('m'), wire nodes ('w'), and source nodes ('s'). Measurements and sources are connected to shared wire nodes, so that blobs sharing wires end up in the same connected component. `clusters()` calls `boost::connected_components` to partition the graph, then returns cluster sets containing only 'm' and 's' nodes (wire nodes are filtered out).

### LASSO-Based Charge Solving

`Solving` (RaySolving.h:55-88, RaySolving.cxx:57-184) takes grouped clusters and solves `m = G * s`. In `add()` (line 86-116), each cluster's measurement nodes are aggregated into a single measurement node (summed value, averaged weight), and source nodes are linked to it. In `solve()` (line 156-184), connected components are discovered again, and for each independent sub-problem, `solve_one()` (line 118-154) constructs a dense geometry matrix G (binary 0/1 entries indicating source-to-measurement edges), a measurement vector, and calls `Ress::solve` with LASSO parameters to obtain per-source charge values.

## Potential Bugs

### [BUG-1] blobref_hash casts pointer to int* -- undefined behavior on 64-bit

- **File**: `util/inc/WireCellUtil/RayClustering.h:18`
- **Severity**: High
- **Description**: The hash functor casts `&*blobref` (a `const Blob*`) to `int*` before hashing via `std::hash<int*>()`. While the intent is to hash the address, the intermediate cast to `int*` is a `reinterpret_cast` through C-style syntax from a `const` pointer to a non-const `int*`, discarding const-qualification. More importantly, on 64-bit systems `int*` and `const Blob*` have the same size so no truncation occurs in practice, but the cast is semantically incorrect -- it should use `std::hash<const void*>()` or `std::hash<size_t>()` with `reinterpret_cast<size_t>(...)` instead. If a platform had different pointer representations for different types, this could produce incorrect hashes.
- **Impact**: In practice on standard 64-bit platforms, the hash likely works correctly since all pointers are 8 bytes. However, the code is technically undefined behavior (violating strict aliasing) and would produce incorrect hashes on any platform where pointer representation differs by pointee type. The `const`-cast is also a const-correctness violation.

### [BUG-2] Tolerance bounds clamping mixes signed and unsigned types

- **File**: `util/src/RayClustering.cxx:40-41`
- **Severity**: Medium
- **Description**: The lower bound uses `std::max(0, strip.bounds.first - tolerance)` where both arguments are `int` (signed). The upper bound uses `std::min(proj.size(), (size_t)strip.bounds.second + tolerance)`. Here `strip.bounds.second` is `grid_index_t` (int, signed) and is cast to `size_t` (unsigned). If `strip.bounds.second` were negative (an error state or uninitialized), the cast to `size_t` would wrap to a very large number, and `std::min` with `proj.size()` would clamp it, but the resulting range could still be nonsensically large. Additionally, the result of the upper-bound `std::min` is `size_t`, but it is passed into `select()` which expects `grid_range_t` = `pair<int,int>`. When `hbound` (size_t) exceeds `INT_MAX`, the implicit narrowing conversion to `int` is undefined behavior.
- **Impact**: Could cause out-of-range access or nonsensical iteration ranges if blob strip bounds are negative or extremely large. In typical use the values are small non-negative integers so this does not manifest.

### [BUG-3] No error handling on LASSO solve result

- **File**: `util/src/RaySolving.cxx:146-153`
- **Severity**: Medium
- **Description**: `Ress::solve()` is called with LASSO parameters but the return value is used unconditionally. If the solver fails to converge (exceeds `max_iter`), the returned vector may contain non-converged values. There is no check for convergence, no residual validation, and no logging of potential solver failure. The `Ress::Params` struct has `max_iter = 100000` and `tolerance = 1e-3`, but there is no mechanism to detect or report failure.
- **Impact**: Non-converged or numerically unstable solutions would be silently accepted and stored in the solution map, leading to incorrect charge estimates propagated downstream without any warning.

### [BUG-4] Measurement weight averaging may be physically incorrect

- **File**: `util/src/RaySolving.cxx:111`
- **Severity**: Low
- **Description**: When multiple measurement nodes from a Grouping cluster are merged into a single Solving measurement node, the values are summed (`total_value`) but the weight is averaged (`total_weight / nms`). If measurements represent independent observations, their combined weight should arguably be the sum (or combined via inverse-variance weighting), not the arithmetic average. Averaging reduces the effective weight as more measurements are combined.
- **Impact**: The resulting measurement weight may underweight well-constrained measurements in the LASSO solver, leading to sub-optimal charge solutions. The `fixme` comment on line 143 also notes that measurement weights are not used in the geometry matrix construction.

### [BUG-5] measurement_node uses magic constant for ident

- **File**: `util/src/RaySolving.cxx:80`
- **Severity**: Low
- **Description**: `measurement_node()` sets `ident = 0xdeadbeef` for all measurement nodes. While measurement idents are not used downstream (only source idents appear in the solution map), this magic constant could cause confusion during debugging and means measurement nodes are not individually identifiable.
- **Impact**: Minimal functional impact, but reduces debuggability and traceability.

## Efficiency Considerations

### [EFF-1] Dense matrix for sparse bipartite graph

- **File**: `util/src/RaySolving.cxx:125`
- **Severity**: Medium
- **Description**: `solve_one()` constructs a fully dense `Eigen::MatrixXd` of size `measures x sources` (line 125: `Ress::matrix_t::Zero(measures.size(), sources.size())`), even though the geometry matrix G is binary and typically very sparse -- each source is connected to only a few measurements. For a sub-problem with N sources and M measurements, this allocates O(N*M) memory and the LASSO solver performs dense matrix operations on it.
- **Suggestion**: Use Eigen sparse matrix types (`Eigen::SparseMatrix<double>`) and a sparse-compatible solver. For sub-problems where the matrix is small (common case), the overhead is negligible, but for large connected components with hundreds of sources and measurements, sparse representation would significantly reduce memory and computation.

### [EFF-2] Redundant projections not cached across overlap calls

- **File**: `util/src/RayClustering.cxx:55-56`
- **Severity**: Medium
- **Description**: In the recursive `overlap()` function, each recursion level calls `projection(blobs, layer)` to re-project the filtered blob set onto the next layer. In `associate()` (line 94-108), `overlap()` is called once per blob in set `one`, and each call performs its own recursive projection chain. If many blobs in set `one` produce similar candidate sets, the same projections are recomputed many times.
- **Suggestion**: Cache projections by layer for the second blob set, or restructure the algorithm to project all blobs at each layer once and perform intersection lookups. Alternatively, build a multi-layer index up front.

### [EFF-3] O(total_width) sparse array in projection

- **File**: `util/src/RayClustering.cxx:8-18`
- **Severity**: Low
- **Description**: The `projection()` function builds a `vector<blobvec_t>` that is indexed by absolute grid index and dynamically resized to `max_grid_index + 1`. If grid indices span a wide range but blobs only occupy a small subset, most buckets are empty vectors, wasting memory. The `resize()` call on line 13 may also trigger multiple reallocations as the vector grows.
- **Suggestion**: Use an `unordered_map<grid_index_t, blobvec_t>` for truly sparse occupancy, or pre-compute the maximum grid index across all blobs to allocate once. In practice, grid indices are dense and contiguous so this is a minor concern.

### [EFF-4] Double connected_components call in Solving pipeline

- **File**: `util/src/RaySolving.cxx:40` and `util/src/RaySolving.cxx:159`
- **Severity**: Low
- **Description**: Connected components are computed twice: once in `Grouping::clusters()` to partition nodes into clusters, and again in `Solving::solve()` to partition the solving graph into independent sub-problems. While these operate on different graphs, the second decomposition could potentially be avoided if the cluster structure from Grouping were preserved through the `add()` step.
- **Suggestion**: Propagate cluster identity through `add()` to avoid recomputing connected components on the solving graph.

### [EFF-5] blobref_hash uses double indirection

- **File**: `util/inc/WireCellUtil/RayClustering.h:18`
- **Severity**: Low
- **Description**: The hash functor dereferences the iterator (`&*blobref`) to get the element address, then hashes that pointer. A simpler approach would hash the iterator's internal pointer or index directly, avoiding the dereference.
- **Suggestion**: Since `blobref_t` is a `vector::const_iterator`, a more direct hash could use the distance from a known base or simply `std::hash<const Blob*>()(&*blobref)` with proper typing.

## Summary

| ID    | Category   | Severity | Title                                                     |
|-------|------------|----------|-----------------------------------------------------------|
| BUG-1 | Bug        | High     | blobref_hash casts pointer to int* -- UB on 64-bit        |
| BUG-2 | Bug        | Medium   | Tolerance bounds clamping mixes signed and unsigned types  |
| BUG-3 | Bug        | Medium   | No error handling on LASSO solve result                    |
| BUG-4 | Bug        | Low      | Measurement weight averaging may be physically incorrect   |
| BUG-5 | Bug        | Low      | measurement_node uses magic constant for ident             |
| EFF-1 | Efficiency | Medium   | Dense matrix for sparse bipartite graph                    |
| EFF-2 | Efficiency | Medium   | Redundant projections not cached across overlap calls      |
| EFF-3 | Efficiency | Low      | O(total_width) sparse array in projection                  |
| EFF-4 | Efficiency | Low      | Double connected_components call in Solving pipeline       |
| EFF-5 | Efficiency | Low      | blobref_hash uses double indirection                       |
