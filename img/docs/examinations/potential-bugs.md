# Potential Bugs

Findings from code examination. Organized by severity.
"Confirmed" means the bug is clearly present in the code; "Suspect" means
the behavior may be intentional but appears incorrect or fragile.

---

## HIGH Severity

### 1. Bit-clear logic error in InSliceDeghosting

**File**: `src/InSliceDeghosting.cxx:66`

`pack &= (0 << p)` always clears ALL bits to 0. `(0 << p)` is always 0 regardless
of `p`, so `pack &= 0` clears the entire packed integer.

**Impact**: Every call to `tag(m, k, p, false)` destroys all tag bits for that blob,
not just the intended one. This means clearing any single quality flag (GOOD, BAD,
POTENTIAL_GOOD, etc.) inadvertently clears all flags. This could cause blobs to lose
their GOOD or POTENTIAL_GOOD status when only TO_BE_REMOVED was being cleared.

**Status**: FIXED. Changed to `pack &= ~(1 << p)`.

---

### 2. Division by zero in calculate_wire_overlap

**File**: `src/InSliceDeghosting.cxx:264`

No check that `wires1.size() != 0` before dividing. If a blob has no wires on a
given plane (possible for 2-view blobs on the dead plane), this divides by zero.

**Impact**: Floating-point division by zero produces infinity, which would pass
the `>= m_deghost_th` comparison and potentially cause incorrect ghost removal.

**Status**: FIXED. Added early return of 0.0 when `wires1` is empty.

---

### 3. Unchecked matrix inverse in CSGraph

**File**: `src/CSGraph.cxx:160`

If `mcov` is singular or near-singular (e.g., two measurements with identical
uncertainty = 0, or numerical issues), `mcov.inverse()` produces NaN/Inf values.
The subsequent Cholesky decomposition (`LLT`) of an invalid matrix fails silently
in Eigen (no exception by default), producing garbage whitened matrices.

**Impact**: Silently corrupted charge solutions. The `mcov.sum() == 0.0` early
return (line 88) only catches the all-zero case, not near-singular cases.

**Status**: FIXED. Added `llt.info() != Eigen::Success` check after decomposition;
logs warning and returns early on failure.

---

### 4. Single-blob scaling inconsistency in CSGraph

**File**: `src/CSGraph.cxx:94-104`

For the single-blob special case, the charge was averaged over measurements without
applying `params.scale`. The multi-blob path applies `solution[ind] * params.scale`
(line 200). This created an inconsistency where 1-blob subgraphs get different
scaling than multi-blob subgraphs.

Also, `value_t val` was not explicitly initialized, risking accumulation from
uninitialized state.

**Status**: FIXED. Applied `params.scale` to the single-blob result and initialized
`val` to `{0, 0}`.

---

### 5. Unreachable code in judge_coverage

**File**: `src/Projection2D.cxx:420`

A `return OTHER` statement after the closing brace of the if-else chain that already
returns on all branches. Dead code left from incomplete refactoring.

**Status**: FIXED. Removed unreachable `return OTHER`.

---

## MEDIUM Severity

### 6. Division by zero in Projection2D judge_coverage_alt

**File**: `src/Projection2D.cxx:517-520`

`small_charge` and `small_counts` can theoretically be 0 in the division
`(1 - common_charge / small_charge)`. While earlier returns likely prevent this
in practice, the logic was fragile and non-obvious.

**Status**: FIXED. Added explicit `small_charge == 0 || small_counts == 0` guard
returning OTHER before the division.

---

### 7. Hard-coded 3-plane assumption in BlobGrouping

**File**: `src/BlobGrouping.cxx:52`

```cpp
std::vector<bcdesc::graph_t> bcs(3); // fixme: hard-code 3 planes
```

Self-acknowledged bug. The code assumes exactly 3 wire planes. Would break for
2-plane or 4+ plane detector geometries.

**Impact**: Incorrect behavior or crash for non-standard geometries.

**Status**: NOT FIXED (per request). Pre-existing known limitation.

---

### 8. Incomplete DeadLiveMerging implementation

**File**: `src/DeadLiveMerging.cxx:82-83`

Comments indicate "dummy implementation" and "FIXME: which ident to use?".
The component merges graphs but doesn't perform actual dead/live tagging or
meaningful edge creation.

**Status**: NOT FIXED (per request). Pre-existing incomplete implementation.

---

### 9. Format string type mismatch in GeomClusteringUtil

**File**: `src/GeomClusteringUtil.cxx:26-28`

Used `%d` (integer format) for `double` values. Also, the error message on line 28
said "tmax" but printed `tmin`.

**Status**: FIXED. Changed `%d` to `%f` in both format strings, and fixed `tmin` to
`tmax` in the second error message.

---

### 10. Sentinel value comparison with floating point

**File**: `src/Projection2D.cxx:228-237`

Used floating-point literal `1e12` as sentinel and compared with `==`.
Also conflated "zero charge" with "no measurement".

**Status**: FIXED. Changed to `std::numeric_limits<double>::max()` and used
`sum_n == 0` as the "no data" condition instead of comparing against sentinel.

---

### 11. Hardcoded TODO threshold in blob_weight_uboone

**File**: `src/ChargeSolving.cxx:114-115`

The 300 charge threshold for connectivity-based weighting was hardcoded with a
TODO comment. This is the same value as `m_good_blob_charge_th` in InSliceDeghosting
but was not exposed to configuration.

**Status**: FIXED. Added `good_blob_charge_th` configuration parameter to
ChargeSolving (header + configure + default_configuration). The weighting functions
now receive the threshold as a parameter instead of using a hardcoded value.

---

### 12. BlobDeclustering passes nullptr slice

**File**: `src/BlobDeclustering.cxx:51`

Passed `nullptr` as slice to SimpleBlobSet constructor. Downstream consumers
that access the slice could get a null pointer dereference.

**Status**: FIXED. Now uses the first blob's slice if available, nullptr only when
the blob vector is empty.

---

## LOW Severity

### 13. Inconsistent dead pixel detection thresholds

**File**: `src/Projection2D.cxx`

- `get_projection()` uses `uncer_cut` (default 1e11) to detect dead channels (line 209)
- `judge_coverage()` uses `-uncer_cut` as threshold (line 369)
- `judge_coverage_alt()` uses `(-1) * uncer_cut` for dead pixel counting (line 472)
- Dead pixels are marked with `dead_default_charge = -1e12`

The `judge_coverage()` and `judge_coverage_alt()` use different criteria for "live"
pixels: the former counts anything > `-uncer_cut` as live (including zero charge),
while the latter only counts `x > 0` as live. This is likely intentional (the "alt"
method is more strict) but the difference is not documented.

**Status**: NOT FIXED. This appears to be intentional design for two different
coverage judgment modes, but warrants documentation.

---

### 14. Hash collision risk in pair_hash

**File**: `src/Projection2D.cxx:126-134`

XOR is a weak hash combiner -- `hash(a,b) == hash(b,a)` and collisions are
likely when h1 and h2 are similar. For (channel, slice) pairs this could cause
performance degradation in unordered containers.

**Status**: FIXED. Replaced `h1 ^ h2` with boost::hash_combine-style mixing:
`h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2)`.

---

### 15. Float-to-int precision loss in ChargeSolving

**File**: `src/ChargeSolving.cxx:104`

Cast floating-point time to int via `(int)csg[...].islice->start()`. If slice
start time is not an exact integer, precision loss occurs and the `time > cent_time`
/ `time < cent_time` comparisons operate on truncated values.

**Status**: FIXED. Removed the `(int)` casts, keeping the native floating-point
type for both `cent_time` and `time`. The comparisons now use full precision.

---

### 16. Excessive debug logging in production paths

**Files**: `src/ChargeSolving.cxx`

`dump_cg()` and `dump_sg()` iterate entire graphs to count vertices even when
logging is at info level, since the functions were called unconditionally.

**Status**: FIXED. Guarded `dump_cg()` with `log->level() <= spdlog::level::debug`
and `dump_sg()` loop with `log->level() <= spdlog::level::trace`.

---

### 17. FrameQualityTagging config default error

**File**: `src/FrameQualityTagging.cxx:83`

```cpp
m_n_fire_cut2 = get<int>(cfg, "n_fire_cut2", m_n_cover_cut2);  // BUG: wrong default
```

Copy-paste error: the default value for `m_n_fire_cut2` used `m_n_cover_cut2`
instead of `m_n_fire_cut2`. When the config key is absent, this causes
`m_n_fire_cut2` to get the value of `m_n_cover_cut2` instead of its own default.

**Status**: FIXED. Changed default to `m_n_fire_cut2`.

---

### 18. ShadowGhosting assumes non-empty string

**File**: `src/ShadowGhosting.cxx:49`

Accesses `m_shadow_type[0]` without checking string is non-empty. If configured
with an empty string, this is undefined behavior.

**Status**: FIXED. Added empty string check with error log and early return.
