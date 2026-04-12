# Ray Grid & Tiling - Code Examination

## Overview

This code implements a 2D ray grid coordinate system and blob-finding tiling algorithm for wire-cell reconstruction. The ray grid provides geometric primitives for working with parallel wire layers projected onto a common plane, while the tiling code builds "blobs" -- polygonal regions defined by the intersection of active wire strips across multiple layers. Supporting utilities generate symmetric ray pairs for standard detector geometries, simulate simple measurements, and compute inter-blob distances.

## Files Examined

- `util/inc/WireCellUtil/RayGrid.h` -- Ray grid coordinate system declarations
- `util/src/RayGrid.cxx` -- Ray grid coordinate system implementation
- `util/inc/WireCellUtil/RayTiling.h` -- Tiling, blob, strip, activity declarations
- `util/src/RayTiling.cxx` -- Tiling, blob, strip, activity implementation
- `util/inc/WireCellUtil/RayHelpers.h` -- Helper declarations for symmetric ray pairs and simulation
- `util/src/RayHelpers.cxx` -- Helper implementation
- `util/inc/WireCellUtil/RayTools.h` -- Relative distance utility declarations
- `util/src/RayTools.cxx` -- Relative distance utility implementation

## Algorithm Description

### Ray Grid Coordinate System (RayGrid.h / RayGrid.cxx)

The `Coordinates` class builds a complete geometric framework from a set of ray pairs, one pair per layer. Each ray pair defines two parallel rays whose separation determines the pitch direction and magnitude for that layer. The construction proceeds in three phases:

1. **Per-layer quantities** (lines 54-70 of RayGrid.cxx): For each layer, computes the pitch vector (via `ray_pitch` on the two parallel rays), its magnitude (`m_pitch_mag`), its unit direction (`m_pitch_dir`), and the center of ray 0 (`m_center`). All vectors are projected into the 2D plane by zeroing the normal axis.

2. **Pair-wise quantities** (lines 73-118): For each ordered pair of layers (l,m), computes the zero-crossing point `m_zero_crossing(l,m)` -- where ray 0 of layer l meets ray 0 of layer m -- and the ray-jump vector `m_ray_jump(l,m)` -- a vector along ray 0 of layer l from the crossing of ray 0 of layer m to the crossing of ray 1 of layer m. The diagonal elements store the center point (zero crossing) and the ray unit direction (ray jump).

3. **Triple-layer coefficients** (lines 122-152): Pre-computes two rank-3 tensors `m_a[l][m][n]` and `m_b[l][m][n]` that enable fast pitch-location lookups. Given a crossing of ray i in layer l and ray j in layer m, the pitch coordinate in layer n is: `j * a[l][m][n] + i * a[m][l][n] + b[l][m][n]`. The `a` tensor stores projections of ray-jump vectors onto pitch directions; `b` stores the offset between zero-crossing and center projections.

### Ring Points Sorting (RayGrid.cxx lines 177-207)

The `ring_points()` method takes a set of crossings (corner pairs), computes their 2D positions, finds the centroid, and sorts them by angle (via `atan2`) around that centroid. This produces a convex-hull-like ordering suitable for polygon rendering. The angle is computed using the Y=dir[2] and X=dir[1] components (Z and Y in detector coordinates), matching the 2D plane when normal_axis=0.

### Activity and Strip Construction (RayTiling.h / RayTiling.cxx)

An `Activity` represents a span of measured values along a single layer's pitch direction, starting at an integer offset. The constructor (lines 34-53 of RayTiling.cxx) trims leading and trailing sub-threshold values. The `active_ranges()` method (lines 91-114) scans for contiguous above-threshold regions; `make_strips()` converts these into `Strip` objects with half-open `[first, second)` grid-index bounds.

### Blob Construction via Strip Intersection (RayTiling.cxx)

The `Blob::add()` method (lines 164-237) incrementally builds a blob by intersecting strips:

- **First strip**: Stored directly, no corners.
- **Second strip**: Four corners from all pairwise edge crossings of the two strips.
- **Third+ strip**: Existing corners are tested for containment in the new strip. New corners (formed by crossing the new strip's edges with each existing strip's edges) are tested for containment in all other existing strips. The `in_strip()` helper (lines 142-162) applies a "nudge" to combat floating-point edge cases: the relative pitch is shifted toward the blob center by `nudge` fraction of a pitch before floor-truncation.

### Tiling and Projection (RayTiling.cxx)

The `Tiling` class orchestrates blob building:
- `operator()(Activity)` (lines 275-284): Creates one initial blob per active strip.
- `operator()(blobs_t, Activity)` (lines 359-380): Refines existing blobs with a new layer. Each blob is projected into the activity's layer to find the relevant pitch range; the activity is trimmed to that range, split into strips, and each strip is added to a copy of the blob.
- `projection()` (lines 286-330): Finds the pitch range covered by blob corners in the activity's layer, clips to the activity's own range, and returns the subspan.

### Pruning (RayTiling.cxx lines 390-444)

Post-construction, `prune()` tightens strip bounds to match the actual corner-defined region. For each wire-plane layer (index >= 2), it collects pitch projections of all corners (both direct grid indices and computed pitch-relative values for non-participating layers), then snaps the strip bounds inward using a nudge tolerance to handle near-integer values.

### Symmetric Ray Pair Generation (RayHelpers.cxx)

`symmetric_raypairs()` generates a standard 5-layer ray grid: horizontal bounds, vertical bounds, two induction planes at +/- angle, and a vertical collection plane. The induction wire positions are computed by projecting pitch jumps onto the Y and Z axes.

### Measurement Simulation (RayHelpers.cxx lines 62-99)

`make_measures()` takes a set of 3D points and bins them into per-layer histograms based on the pitch index of each point relative to each layer's center and pitch direction. Layers 0-1 (bounds layers) only accept points with pitch_index == 0.

### Relative Distance (RayTools.cxx)

`relative_distance()` computes a normalized distance between two strips (or two sets of strips). For a single pair of grid ranges, it divides the absolute difference of centers by the sum of widths. For strip vectors, it returns the mean of squared per-layer distances (an L2-like norm).

## Potential Bugs

### [BUG-1] active_ranges() uses literal 0.0 instead of m_threshold for exit test
- **File**: `util/src/RayTiling.cxx:103`
- **Severity**: Medium
- **Description**: The `active_ranges()` method uses `*it > m_threshold` (line 98) to detect entry into an active region, but uses `*it <= 0.0` (line 103) to detect exit. If a non-zero threshold is set (e.g., 0.5), a value like 0.3 would not trigger the exit condition even though it is below threshold. The exit test should be `*it <= m_threshold` for consistency with the entry test.
- **Impact**: With non-zero thresholds, active ranges may extend past the intended boundary, including sub-threshold bins. Strips would be wider than intended, potentially producing incorrect blobs.

### [BUG-2] ring_points() atan2 uses hardcoded axis indices [2],[1]
- **File**: `util/src/RayGrid.cxx:194`
- **Severity**: Low
- **Description**: The `ring_points()` method computes angles as `atan2(dir[2], dir[1])`, which uses the Z and Y components. This is correct only when `normal_axis=0` (X axis is the normal). If a different normal axis is used (e.g., normal_axis=2), the atan2 would use the wrong plane components, producing incorrect angular sorting. The axis choice should be derived from the `normal_axis` parameter used during construction.
- **Impact**: Incorrect polygon vertex ordering if the ray grid is constructed with a non-default normal axis. In practice, the toolkit always uses normal_axis=0, so this is unlikely to cause issues in current usage.

### [BUG-3] Nudge only applied to layer >= 2 in in_strip()
- **File**: `util/src/RayTiling.cxx:149`
- **Severity**: Low
- **Description**: The `in_strip()` function applies the floating-point nudge only for `strip.layer >= 2`, explicitly skipping layers 0 and 1 (the bounds layers). While this is likely intentional (bounds layers have only one pitch bin), the hardcoded layer index creates an implicit coupling to the convention that layers 0 and 1 are always bounds layers. If a Coordinates is constructed with a different layer ordering, this guard would be incorrect.
- **Impact**: Minimal in current usage since the convention is consistently followed. Could cause subtle edge-case failures if the layer ordering convention changes.

### [BUG-4] Division by zero possible in relative_distance()
- **File**: `util/src/RayTools.cxx:8`
- **Severity**: Medium
- **Description**: The function computes `num/den` where `den = std::abs(a.first-a.second) + std::abs(b.first-b.second)`. If both strips have zero width (first==second), den is zero, causing a division-by-zero resulting in NaN or infinity. While `Blob::valid()` rejects zero-width strips, there is no guard here and the function is public.
- **Impact**: Undefined floating-point behavior if called with two zero-width strips. Could propagate NaN through downstream calculations.

### [BUG-5] abs() applied to integer types in relative_distance() may use integer abs
- **File**: `util/src/RayTools.cxx:7-8`
- **Severity**: Medium
- **Description**: `grid_range_t` is `std::pair<grid_index_t, grid_index_t>` where `grid_index_t` is `int`. The expressions `a.first+a.second` and `a.first-a.second` are integer arithmetic. `std::abs()` on integers returns an integer. The result is then implicitly converted to `double` for the division. While the integer `abs` gives the same magnitude as `fabs` would on the converted value, the intermediate subtraction `(a.first+a.second) - (b.first+b.second)` is computed entirely in `int` and could overflow for very large grid indices (unlikely in practice). More importantly, the integer division semantics could surprise maintainers.
- **Impact**: Low practical impact since grid indices are typically small, but the code mixes integer and floating-point semantics in a way that could lead to confusion or subtle bugs with large indices.

### [BUG-6] prune() assumes strip ordering matches layer indices
- **File**: `util/src/RayTiling.cxx:440`
- **Severity**: Low
- **Description**: The `prune()` function accesses `strips[layer]` using the layer index as an array index (line 440: `strips[layer].bounds`). This assumes strips are stored in layer-index order (strip at position 0 is for layer 0, etc.). This is indeed how `make_blobs()` constructs them (activities are processed in order), but the assumption is not enforced and could break if blobs are constructed differently.
- **Impact**: Incorrect strip pruning if strips are not ordered by layer index. Currently safe because `make_blobs()` maintains the ordering.

### [BUG-7] make_measures() truncates pitch index via integer division
- **File**: `util/src/RayHelpers.cxx:76`
- **Severity**: Low
- **Description**: The pitch index is computed as `const int pit_ind = pit.dot(rel) / pitch_mags[ilayer]` which truncates toward zero (C++ integer conversion), not floor. For negative pitch values (which are filtered out on line 77), this is moot. But the truncation behavior differs from the `pitch_index()` method in Coordinates which uses `std::floor()`. For positive values the two are equivalent, but the inconsistency could cause confusion.
- **Impact**: No practical impact since negative indices are filtered, but the inconsistency with `pitch_index()` is a maintenance concern.

## Efficiency Considerations

### [EFF-1] Full angle sort in ring_points() vs cheaper radial sort
- **File**: `util/src/RayGrid.cxx:190-206`
- **Description**: The `ring_points()` method computes `atan2()` for every corner point, then sorts by angle using an indirect index array. The `atan2` call is relatively expensive (transcendental function). For convex polygon vertex ordering, a cheaper alternative exists: since all points lie on a convex hull (they are intersections of strip boundaries), one could use a cross-product-based comparison that avoids transcendental functions entirely.
- **Suggestion**: Replace `atan2`-based sorting with a cross-product-based angular comparison: `(a-center).cross(b-center)[normal_axis] > 0`. This avoids transcendental function calls while producing the same ordering. Alternatively, if the number of corners is always small (typically 4-8), the performance difference is negligible and clarity may be preferred.

### [EFF-2] Floating-point accumulation in center computation for blob corners
- **File**: `util/src/RayTiling.cxx:183-186`
- **Description**: When computing `center_in_new`, the code accumulates `pitch_location()` results for every existing corner by simple summation. This is also done for each old strip in `center_in_old` (lines 200-206). Each call to `pitch_location()` involves tensor lookups and a multiply-add. For blobs with many corners and many strips, this is O(corners * strips) work repeated for each new strip addition.
- **Suggestion**: Since pitch_location is linear in the grid indices, the center could be computed once from the mean grid indices of corners rather than re-summing pitch_location for each strip. However, the number of corners and strips is typically small (< 10), so this is a minor concern.

### [EFF-3] Redundant empty check in Tiling::projection()
- **File**: `util/src/RayTiling.cxx:307-309`
- **Description**: After populating `pitches` from `blob.corners()` (which was already checked for emptiness on line 300), the code checks `if (pitches.empty())` again on line 307. Since the loop on line 303 iterates over `blob.corners()` (same source checked on line 300), if corners is non-empty, pitches will be non-empty. This is a dead-code guard.
- **Suggestion**: Remove the redundant check on line 307, or consolidate both checks.

### [EFF-4] Blob copy in tiling refinement
- **File**: `util/src/RayTiling.cxx:370`
- **Description**: In `Tiling::operator()(blobs_t, Activity)`, each iteration copies the entire blob (`Blob newblob = blob`) before adding a new strip. This copies all strips and corners vectors. When there are many strips per blob, this copying overhead accumulates.
- **Suggestion**: If the strip does not produce valid corners, the copy was wasted. Consider checking the projection validity before copying, or using move semantics when only one strip remains.

### [EFF-5] Repeated strip-corner containment testing is O(strips^2 * corners)
- **File**: `util/src/RayTiling.cxx:209-234`
- **Description**: When adding the Nth strip to a blob, `Blob::add()` tests each of the 4 new candidate corners against all N-1 existing strips (lines 214-229). Combined with the surviving-corners test (lines 189-196), the total work per `add()` call is O(existing_corners + 4*(N-1)). Over all N strips, this is O(N^2) in total strip count. For the typical 5-layer case this is fine, but it could become costly if many layers are used.
- **Suggestion**: For the standard 3-5 wire plane case, this is adequate. If scaling to many layers, consider spatial indexing of strip boundaries.

### [EFF-6] symmetric_raypairs() recomputes pitch jump projections
- **File**: `util/src/RayHelpers.cxx:21-26`
- **Description**: The pitch-jump projection for each induction plane is computed twice: once for 0.5*pitch_mag (ray 0) and once for 1.5*pitch_mag (ray 1). The projection formula `pjumpu.dot(pjumpu) / (axis.dot(pjumpu))` scales linearly with pitch_mag magnitude, so the second computation could be derived from the first by scaling.
- **Suggestion**: Minor; compute the base projection once and scale for each ray. The function is called once during initialization so performance is not a concern.

## Summary

| ID | Type | Severity | Title |
|--------|------|----------|-------|
| BUG-1 | Bug | Medium | active_ranges() uses literal 0.0 instead of m_threshold for exit test |
| BUG-2 | Bug | Low | ring_points() atan2 uses hardcoded axis indices |
| BUG-3 | Bug | Low | Nudge only applied to layer >= 2 in in_strip() |
| BUG-4 | Bug | Medium | Division by zero possible in relative_distance() |
| BUG-5 | Bug | Medium | abs() on integer types in relative_distance() |
| BUG-6 | Bug | Low | prune() assumes strip ordering matches layer indices |
| BUG-7 | Bug | Low | make_measures() truncates pitch index via integer division |
| EFF-1 | Efficiency | Low | Full angle sort vs cheaper radial sort in ring_points() |
| EFF-2 | Efficiency | Low | Floating-point accumulation in center computation |
| EFF-3 | Efficiency | Low | Redundant empty check in projection() |
| EFF-4 | Efficiency | Low | Blob copy in tiling refinement |
| EFF-5 | Efficiency | Low | O(strips^2 * corners) containment testing |
| EFF-6 | Efficiency | Low | Recomputed pitch jump projections in symmetric_raypairs() |
