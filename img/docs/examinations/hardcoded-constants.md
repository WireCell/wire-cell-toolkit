# Hardcoded Constants and Magic Numbers

This document catalogs hardcoded constants that may need attention when
porting to new detectors (e.g., PDHD) or when tuning algorithm performance.

---

## UBooNE-Specific Thresholds

These values were tuned for the MicroBooNE detector and may not generalize.

### MaskSlice default thresholds

**File**: `inc/WireCellImg/MaskSlice.h:77-79`

```cpp
m_nthreshold = {3.6, 3.6, 3.6};                    // per-plane multiplier
m_default_threshold = {5.87819e+02 * 4.0,           // U-plane: 2351.3
                       8.36644e+02 * 4.0,            // V-plane: 3346.6
                       5.67974e+02 * 4.0};           // W-plane: 2271.9
```

The `default_threshold` values are UBooNE RMS values multiplied by 4.
These are configurable but the defaults are detector-specific.

### FrameQualityTagging thresholds

**File**: `src/FrameQualityTagging.cxx:49-68`

| Constant | Value | Purpose |
|----------|-------|---------|
| `m_nrebin` | 4 | Rebinning factor |
| `m_length_cut` | 12 | Minimum region size |
| `m_time_cut` | 12 | Region connection threshold |
| `m_ch_threshold` | 100 | Channel coverage threshold |
| `m_n_cover_cut1` | 12 | Coverage cut level 1 |
| `m_n_fire_cut1` | 14 | Fire rate cut level 1 |
| `m_n_cover_cut2` | 6 | Coverage cut level 2 |
| `m_n_fire_cut2` | 14 | Fire rate cut level 2 |
| `m_n_cover_cut3` | 3 | Coverage cut level 3 |
| `m_n_fire_cut3` | 6 | Fire rate cut level 3 |
| `m_fire_threshold` | 0.22 | Multi-plane coincidence fraction |
| `m_global_threshold` | 0.048 | Global busy-ness limit |
| `m_min_time` | 3180 | Time window start (dataset-specific) |
| `m_max_time` | 7870 | Time window end (dataset-specific) |

The `m_min_time` and `m_max_time` are particularly concerning as they
define a time window specific to a particular dataset or run configuration.
All values are configurable but the defaults are UBooNE-tuned.

---

## Charge Thresholds

### Good blob charge threshold: 300

**Files**:
- `inc/WireCellImg/InSliceDeghosting.h:43`: `m_good_blob_charge_th{300.}`
- `src/ChargeSolving.cxx:115`: `if (iblob->value() < 300) continue;`

The value 300 appears in two places:
1. InSliceDeghosting: configurable via `good_blob_charge_th`
2. ChargeSolving `blob_weight_uboone`: **hardcoded with TODO comment**

These should be the same value and both should be configurable.

---

## LASSO Solver Constants

### Lambda and tolerance formulas

**File**: `src/CSGraph.cxx:144-145`

```cpp
double lambda = 3. / total_wire_charge / 2. * params.scale;
double tolerance = total_wire_charge / 3. / params.scale / R_mat.cols() * 0.005;
```

Breaking down:
- `lambda = 1.5 * scale / total_wire_charge` (L1 penalty strength)
- `tolerance = 0.005 * total_wire_charge / (3 * scale * n_blobs)` (convergence)

The factors 3, 2, and 0.005 are unexplained. These control the balance between
fitting the data and regularizing (shrinking) blob charges to zero.

### Max iterations

**File**: `src/CSGraph.cxx:146`

```cpp
rparams = Ress::Params{Ress::lasso, lambda, 100000, tolerance, true, false};
```

Max iterations = 100000 is hardcoded. No check whether convergence was achieved.

### Blob weight constants

**File**: `src/ChargeSolving.cxx:124-130`

```cpp
double weight = 9.;
if (next_con) { weight /= 3.; }
if (prev_con) { weight /= 3.; }
```

Base weight 9 with divisor 3 per connection. The choice of 9 and 3 (= 3^2 and 3^1)
creates a geometric progression: 9, 3, 1. Not documented why these values.

Also in BlobSolving.cxx:44-46:
```cpp
default_weight = 9;  // same value
reduction = 3;       // same divisor
homer = 2;           // "Max Power" -- apparently a joke/placeholder name
```

---

## Sentinel and Marker Values

### Dead pixel markers

| Constant | Value | File | Purpose |
|----------|-------|------|---------|
| `dead_default_charge` | -1e12 | Projection2D.h:61 | Charge for dead channels |
| `uncer_cut` | 1e11 | Projection2D.h:60 | Threshold for dead detection |
| `dummy_error` | 1e12 | MaskSlice.h:73 | Error for dummy planes |
| `masked_error` | 1e12 | MaskSlice.h:75 | Error for masked planes |
| `min_charge init` | 1e12 | Projection2D.cxx:228 | Sentinel for "no data" |

The use of 1e12 / -1e12 as sentinels is a recurring pattern. These values
serve as "effectively infinite" markers but are fragile if actual data
approaches these magnitudes (unlikely in practice but architecturally poor).

---

## Geometry Constants

### 3-plane assumption

**Files**: Multiple

| File | Line | Constant |
|------|------|----------|
| BlobGrouping.cxx | 52 | `bcs(3)` -- explicitly noted as bug |
| GridTiling.cxx | 147 | `iplane >= 3` check |
| Projection2D.cxx | 198-200 | `kUlayer, kVlayer, kWlayer` init |
| MaskSlice.cxx | 84-86 | `active_planes[0..2]` |
| FrameQualityTagging.cxx | various | Per-plane arrays of size 3 |
| InSliceDeghosting.cxx | 430 | `max(kUlayer, kVlayer, kWlayer)` |

While configurable in some places, the default and internal logic assumes
exactly 3 wire planes throughout.

### RayGrid bounding layers

**Files**: GridTiling.cxx:79, BlobSetReframer.cxx:52, 129

```cpp
nbounds_layers = 2;  // number of bounding box layers in RayGrid
strip.layer < 2;     // skip bounding layers
```

Assumes RayGrid uses exactly 2 bounding layers. This is a RayGrid convention
but should be queried from the RayGrid API rather than hardcoded.

### Dead time offset

**File**: `src/GeomClusteringUtil.cxx:34`

```cpp
static bool adjacent_dead(..., const double offset = 4*500*units::us)
```

Hard-coded dead-time tolerance: 4 ticks at 500 us/tick = 2 ms.
Not configurable. The TODO comment (line 33) says "confirm with Xin".

### Geometric clustering tolerances

**File**: `src/GeomClusteringUtil.cxx:56-66`

| Policy | max_rel_diff | gap_tol |
|--------|-------------|---------|
| uboone | 2 | {1:2, 2:1} |
| uboone_local | 2 | {1:2, 2:2} |
| simple | 1 | {1:0} |

These are hardcoded in the source. The policy names ("uboone", "uboone_local")
suggest detector-specific tuning.

---

## Deghosting Thresholds

### InSliceDeghosting defaults

**File**: `inc/WireCellImg/InSliceDeghosting.h:43-47`

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `m_good_blob_charge_th` | 300. | Charge threshold for "good" blobs |
| `m_deghost_th` | 0.75 | Wire overlap ratio threshold |
| `m_deghost_th1` | 0.5 | Secondary deghosting threshold |

### Projection2D coverage thresholds

**File**: `src/Projection2D.cxx:400-403`

```cpp
return x < -0.01;  // hardcoded tolerance for negative mask diff
return x > 0.01;   // hardcoded tolerance for positive mask diff
```

These 0.01 thresholds determine whether two projections differ. Not configurable.

---

## Summary: Constants Requiring Review for PDHD

| Priority | Constant | Current Value | Why Review |
|----------|----------|--------------|------------|
| HIGH | MaskSlice default_threshold | UBooNE RMS x 4 | Detector-specific |
| HIGH | FrameQualityTagging time window | 3180-7870 | Dataset-specific |
| HIGH | blob charge threshold 300 | Hardcoded in ChargeSolving | Should be configurable |
| MEDIUM | LASSO lambda/tolerance formulas | 3/2, 0.005 | Unexplained physics |
| MEDIUM | Blob weights 9/3 | Hardcoded | Unexplained |
| MEDIUM | Dead time offset 4*500us | Hardcoded | Detector-specific |
| LOW | Sentinel values 1e12 | Convention | Fragile pattern |
| LOW | 3-plane assumption | Throughout | Limits generality |
