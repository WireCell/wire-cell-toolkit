# Spatial & Geometry - Code Examination

## Overview

This group of code provides spatial and geometric primitives for the Wire-Cell Toolkit: axis-aligned bounding boxes with intersection/union/containment queries, ray-box intersection algorithms, wire plane pitch/impact binning and coordinate transforms (Pimpos), and a 2D interval tree for rectangular region tracking (Rectangles). Together they support detector geometry calculations, spatial queries, and pixelization of wire plane data.

## Files Examined

- `util/inc/WireCellUtil/BoundingBox.h`
- `util/src/BoundingBox.cxx`
- `util/inc/WireCellUtil/Intersection.h`
- `util/src/Intersection.cxx`
- `util/inc/WireCellUtil/Pimpos.h`
- `util/src/Pimpos.cxx`
- `util/inc/WireCellUtil/Rectangles.h`

## Algorithm Description

### BoundingBox

An axis-aligned bounding box (AABB) represented as a `Ray` (pair of `Point`s) storing the min-corner (`m_bounds.first`) and max-corner (`m_bounds.second`). The box is lazily initialized via an `m_initialized` flag.

Key operations:
- **Expansion** (`operator()`): Grows the box to encompass a new point by taking per-axis min/max. First point initializes the box.
- **Inside test**: Per-axis range check, skipping zero-width dimensions.
- **Union** (`unite`): Per-axis min of mins, max of maxes.
- **Intersection** (`intersect(BoundingBox)`): Per-axis overlap calculation; returns empty box if any axis has no overlap.
- **Distance**: Computes Euclidean distance from a point to the nearest box surface (0 if inside), using per-axis clamped differences.
- **Ray intersection** (`intersect(Ray)`): Tests all 6 faces using `plane_split`/`plane_intersection` from Point.h, collects up to 2 intersection points, then orders them along the ray direction.
- **Cropping** (`crop`): Intersects a finite segment with the box, preserving original endpoints that lie inside.
- **Padding**: `pad_rel` scales padding by the diagonal vector; `pad_abs` scales by the diagonal unit vector.
- **Axis distances**: Returns signed distances from a point to the two box walls along a given axis, but only if the point's other two coordinates are within the box.

### Intersection (ray-box)

Provides `box_intersection` in two forms:

1. **Per-axis** (`box_intersection(axis0, ...)`): For a given axis, computes parametric ray-plane intersections with the two box faces orthogonal to that axis. For each face, computes `scale = (intercept - point[axis]) / dir[axis]`, then checks if the other two coordinates of the hit point lie within the box bounds. Returns a bitmask (0-3) indicating which faces were hit. Results are ordered to be parallel to the ray direction.

2. **Full box** (`box_intersection(bounds, point, dir, hits)`): Calls the per-axis version for all 3 axes, collects unique hit points in a `PointSet`, and returns the first two distinct hits. Uses the `PointSet` (a `std::set` with `ComparePoints` comparator using 1e-10 tolerance) to deduplicate corner hits.

### Pimpos (Pitch-Impact-Position)

Encapsulates geometry and binning for a plane of parallel, equidistant wires:

- **Coordinate system**: Three orthogonal axes -- axis0 = normal (wire x pitch), axis1 = wire direction, axis2 = pitch direction. The `transform` method projects a 3D point into this local frame.
- **Wire region binning** (`m_regionbins`): `nwires` bins, each centered on a wire, spanning half a pitch on either side. Total range is `[minwirepitch - pitch/2, maxwirepitch + pitch/2]`.
- **Impact binning** (`m_impactbins`): `nwires * nbins_per_wire` bins over the same range. Impact positions are bin edges and are in-phase with wire centers.
- **Closest wire/impact** (`closest`): Finds the nearest wire via `region_binning().bin(pitch)`, then computes the relative impact index by rounding the remainder divided by the impact bin size.
- **Wire impact index** (`wire_impact`): Returns `(0.5 + wireind) * nimpbins_per_wire`, the impact index coincident with a wire.
- **Reflection** (`reflect`): Reflects an impact index through a wire's impact index.

### Rectangles (2D Interval Tree)

A template class providing 2D interval tree functionality by nesting Boost ICL `interval_map`s:

- **Structure**: An outer `interval_map<XKey, ymap_t>` maps X-intervals to inner `interval_map<YKey, Set>` maps. When rectangles overlap, their associated value sets are merged (union via `operator+=` on the set type).
- **Add**: Inserts a half-open rectangle `[x1,x2) x [y1,y2)` with an associated value.
- **Regions**: Iterates all disjoint (non-overlapping) sub-rectangles produced by the ICL splitting, each with the combined set of all original rectangles that cover it.
- **Intersection query**: Restricts the X-map to a query X-interval using `& xi`, then restricts each resulting Y-map to the query Y-interval, returning only regions within the query rectangle.
- **Pixelization** (`pixelize` free function): Maps continuous-domain rectangles to discrete pixel indices using two `Binning` objects, clamping to valid bin ranges.

## Potential Bugs

### [BUG-1] PointSet deduplication uses fixed 1e-10 tolerance in box_intersection

- **File**: `util/src/Intersection.cxx`:64-96 (and `util/src/Point.cxx`:22-27)
- **Severity**: Medium
- **Description**: The full-box `box_intersection` function collects hit points into a `PointSet`, which uses `ComparePoints` with a hardcoded 1e-10 magnitude tolerance for equality. This tolerance is not related to the physical scale of the problem. For very small geometries (sub-nanometer coordinates), genuinely distinct intersection points may be falsely deduplicated. For very large geometries, points that should be considered the same (e.g., corner hits from two different axis projections) may not be merged due to floating-point drift exceeding 1e-10.
- **Impact**: Could return 1 hit instead of 2 (missing an intersection), or fail to merge corner hits (returning duplicate near-identical points consuming both hit slots, potentially masking the true second intersection). The bitmask returned would be incorrect.

### [BUG-2] BoundingBox::intersect(BoundingBox) does not set m_initialized for all dimensions in branch

- **File**: `util/src/BoundingBox.cxx`:224-228
- **Severity**: Low
- **Description**: The `intersect(BoundingBox)` method has an if/else branch on `result.m_initialized` (lines 225-228) that does exactly the same thing in both branches -- setting `result.m_bounds.first[ind]` and `result.m_bounds.second[ind]`. The `m_initialized` flag is only set to `true` after the loop (line 236), so the `if (!result.m_initialized)` branch is always taken. The code works correctly but the if/else is dead logic, suggesting a possible copy-paste error where the first iteration was meant to do something different (e.g., set all three axes of `.first` and `.second` for uninitialized dimensions to avoid stale zero values from the default `Ray`).
- **Impact**: Functionally benign since m_bounds starts as default Ray (origin) and all three axes are written in the loop. But the dead branch is misleading.

### [BUG-3] BoundingBox::intersect(Ray) may produce wrong results when line is parallel to a face

- **File**: `util/src/BoundingBox.cxx`:306-352
- **Severity**: Low
- **Description**: The `intersect(Ray)` method delegates to `plane_split` and `plane_intersection` for each of the 6 faces. When the ray lies exactly in a face plane (parallel to the face and coincident with it), `plane_split` may return false (the ray does not "split" the plane), so valid intersection segments along that face would be missed. This is an edge case for rays running along a box face.
- **Impact**: Edge case only. Rays exactly on a face boundary are unlikely in practice but could produce an empty result when one intersection point should be detected.

### [BUG-4] C-style round() used instead of std::round() in Pimpos::closest

- **File**: `util/src/Pimpos.cxx`:39
- **Severity**: Low
- **Description**: The expression `int(round(remainder / m_impactbins.binsize()))` uses unqualified `round()` from `<math.h>` (C library), which returns `double`. While functionally equivalent to `std::round()` in most implementations, the C-style `round` is technically not guaranteed to be in scope without `<cmath>` being included, and mixing C and C++ math functions is not best practice. Additionally, the cast `int(...)` truncates toward zero, which for `round`'s return value is correct, but a direct `std::lround()` would be more explicit and avoid any potential for intermediate floating-point representation issues.
- **Impact**: No practical bug in current compilers, but non-idiomatic.

### [BUG-5] No removal operation in Rectangles

- **File**: `util/inc/WireCellUtil/Rectangles.h`:76-166
- **Severity**: Low
- **Description**: The `Rectangles` class provides only `add` (and `operator+=`) but no `remove` or `erase` operation. While Boost ICL `interval_map` supports subtraction (`operator-=`), this is not exposed. Once a rectangle is added, there is no way to remove it short of rebuilding the entire structure.
- **Impact**: If callers ever need to remove rectangles, they must reconstruct the Rectangles object from scratch. This is a design limitation rather than a bug.

### [BUG-6] pixelize may produce empty intervals when bin mapping collapses range

- **File**: `util/inc/WireCellUtil/Rectangles.h`:189-190
- **Severity**: Medium
- **Description**: In the `pixelize` function, the lower and upper bounds of an interval are both clamped to `[0, nxbins-1]` (and similarly for Y). If the lower bound maps to a bin index >= the upper bound's bin index (e.g., a very narrow rectangle that fits within a single bin, where `lower(xi)` and `upper(xi)` map to the same bin), the resulting `right_open(n, n)` interval is empty. Furthermore, `upper(xi)` for a half-open input interval returns the exclusive upper bound, which `xbins.bin()` may map to a bin one past the intended range. The function does not check for or skip empty intervals before calling `ret.add()`.
- **Impact**: Could produce empty intervals that are silently added to the Rectangles structure, or off-by-one bin assignments for the upper bound of intervals. Narrow rectangles smaller than one pixel may be lost entirely.

### [BUG-7] pixelize takes xbins by value instead of by reference

- **File**: `util/inc/WireCellUtil/Rectangles.h`:179
- **Severity**: Low
- **Description**: The `pixelize` function signature is `pixelize(const Rec& inrec, const Binning xbins, const Binning& ybins)`. The `xbins` parameter is passed by value (missing `&`) while `ybins` is passed by const reference. This is likely an oversight.
- **Impact**: Unnecessary copy of the `Binning` object on every call. No correctness issue.

### [BUG-8] BoundingBox::pad_abs produces non-uniform padding for non-cubic boxes

- **File**: `util/src/BoundingBox.cxx`:90-95
- **Severity**: Medium
- **Description**: `pad_abs(double distance)` computes `ray_unit(m_bounds) * distance`, where `ray_unit` gives the unit vector along the diagonal of the box. For a non-cubic box, this unit vector has unequal components, so the padding is not uniform along each axis. For example, a long thin box would get most of its padding along the long axis and very little along the short axes. The name "pad_abs" and the parameter name "distance" suggest uniform absolute padding, but the actual behavior is proportional to the diagonal direction.
- **Impact**: Callers expecting uniform padding in all directions will get incorrect results. A flat box (e.g., one axis near zero) will get almost no padding on that axis.

## Efficiency Considerations

### [EFF-1] BoundingBox::intersect(Ray) collects all face intersections into a vector

- **File**: `util/src/BoundingBox.cxx`:309
- **Severity**: Low
- **Description**: The method tests all 6 faces and pushes intersection points into a `std::vector<Vector>`. Since at most 2 intersections are possible for a convex box, the early-exit optimization (stop after finding 2) is not implemented.
- **Suggestion**: Return early once `intersections.size() == 2` to avoid unnecessary plane intersection tests.

### [EFF-2] box_intersection per-axis normalizes bounds on every call

- **File**: `util/src/Intersection.cxx`:14-19
- **Severity**: Low
- **Description**: Each call to the per-axis `box_intersection` copies the bounds points and sorts each axis via `std::swap`. When called from the full-box version (3 times in a loop), this normalization is repeated 3 times.
- **Suggestion**: Normalize bounds once in the caller and pass pre-normalized bounds, or cache the normalized bounds.

### [EFF-3] Rectangles::regions() creates a new vector on every call

- **File**: `util/inc/WireCellUtil/Rectangles.h`:135-144
- **Severity**: Low
- **Description**: The `regions()` method allocates and returns a new `std::vector<region_t>` every time. If called frequently (e.g., in a loop), this incurs repeated heap allocations.
- **Suggestion**: Acceptable for current usage patterns. Could add an output iterator version or cache if profiling shows it as a hotspot.

### [EFF-4] BoundingBox::intersect(BoundingBox) and overlaps() both normalize min/max per axis

- **File**: `util/src/BoundingBox.cxx`:209-213, 249-252
- **Severity**: Low
- **Description**: Both `intersect(BoundingBox)` and `overlaps()` call `std::min`/`std::max` on each axis of both boxes to handle potentially un-normalized bounds. Since `operator()` always maintains `first <= second` per axis, this normalization is redundant for well-formed BoundingBox objects.
- **Suggestion**: If the class invariant guarantees normalized bounds, these min/max calls can be removed. Alternatively, document and enforce the invariant.

## Summary

| ID | Type | Severity | Title |
|----|------|----------|-------|
| BUG-1 | Bug | Medium | PointSet deduplication uses fixed 1e-10 tolerance in box_intersection |
| BUG-2 | Bug | Low | Dead if/else branch in BoundingBox::intersect(BoundingBox) |
| BUG-3 | Bug | Low | Ray-face intersection missed when ray lies in face plane |
| BUG-4 | Bug | Low | C-style round() instead of std::round() in Pimpos::closest |
| BUG-5 | Bug | Low | No removal operation in Rectangles |
| BUG-6 | Bug | Medium | pixelize may produce empty intervals or off-by-one bin mapping |
| BUG-7 | Bug | Low | pixelize takes xbins by value instead of by reference |
| BUG-8 | Bug | Medium | pad_abs produces non-uniform padding for non-cubic boxes |
| EFF-1 | Efficiency | Low | No early exit in BoundingBox::intersect(Ray) after finding 2 hits |
| EFF-2 | Efficiency | Low | Repeated bounds normalization in per-axis box_intersection |
| EFF-3 | Efficiency | Low | regions() allocates new vector on every call |
| EFF-4 | Efficiency | Low | Redundant min/max normalization in intersect/overlaps |
