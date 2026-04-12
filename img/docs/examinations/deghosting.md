# Deghosting Algorithms

## Background: What are Ghosts?

In a 3-plane wire TPC, each plane provides a 1D projection of the charge
distribution. A blob is the 3D intersection of active wire strips from all planes.
"Ghosts" are spurious blobs that arise when wires from unrelated true deposits
happen to cross, creating false intersections. Deghosting removes these artifacts.

## 1. InSliceDeghosting (`src/InSliceDeghosting.cxx`, 39K)

### Algorithm Overview

Operates within each time slice to identify and remove ghost blobs. Three phases:

#### Phase 1: Blob Quality Identification (`blob_quality_ident`, lines 269-308)

For each blob in the cluster graph:
1. If blob charge > `m_good_blob_charge_th` (default 300): tag as GOOD + POTENTIAL_GOOD
2. Check b-b connected neighbors: if any neighbor has charge > threshold, tag
   current blob as POTENTIAL_GOOD
3. This identifies blobs that are likely real based on charge and connectivity

#### Phase 2: Local Deghosting (`local_deghosting1`, lines 309-518)

Per time slice, performs wire-score-based ghost identification:

1. **Group blobs by view count** (line 336):
   - `view_groups[3]`: blobs with all 3 planes live (higher confidence)
   - `view_groups[2]`: blobs with 2 planes live (more likely ghosts)

2. **Build wire score map** (lines 346-358):
   For 3-view blobs, count how many blobs use each wire. High-score wires are
   shared by many blobs (ambiguous). Low-score wires belong to few blobs (more unique).
   Also accumulate scores for 2-view blob wires on live planes (lines 362-384).

3. **Find "cannot remove" blobs** (lines 386-400):
   A 2-view blob is protected if it is `adjacent()` to at least 2 POTENTIAL_GOOD
   3-view blobs. The adjacency test (lines 145-198) checks per-plane wire overlap:
   - Score 2 if channels overlap
   - Score 1 if channels are adjacent (differ by 1)
   - Score 0 if neither (immediately fails for that plane)
   - Total score >= 5 across planes means adjacent

4. **Score-based removal** (lines 402-516):
   For each 2-view blob:
   - Compute `blob_high_score_map`: max across planes of (average 1/wire_score)
   - For each overlapping higher-scoring blob, check:
     - `overlap_ratio >= m_deghost_th` (default 0.75)
     - `current_q2 > current_q1 * m_deghost_th`
   - If 2 such blobs found AND blob not in `cannot_remove`: tag as TO_BE_REMOVED

5. **3-view blob scores** (lines 433-461):
   Same wire-score computation for 3-view blobs, but POTENTIAL_GOOD blobs
   automatically get score 1 (lines 460).

#### Phase 3: Graph Filtering

Remove TO_BE_REMOVED blobs from the cluster graph. Optionally re-cluster remaining
blobs using `geom_clustering()` with the configured policy.

### Bit-tag System

Uses bit-packed integer tags per blob (enum BLOB_QUALITY_BITPOS):
- Bit 0: GOOD
- Bit 1: BAD
- Bit 2: POTENTIAL_GOOD
- Bit 3: POTENTIAL_BAD
- Bit 4: TO_BE_REMOVED

Tags set/cleared via template `tag()` function (lines 57-76).

---

## 2. Projection2D (`src/Projection2D.cxx`, 19K)

### Algorithm Overview

Builds 2D sparse matrix projections of 3D blob clusters onto wire planes.

#### `get_projection()` (lines 136-266)

For a group of blob vertices:
1. Map slices to blobs (`map_s2vb`) and blobs to channels (`map_b2c`)
2. For each (slice, blob, channel):
   - Get channel charge from slice activity
   - If uncertainty > `uncer_cut` (default 1e11): mark as dead (`dead_default_charge = -1e12`)
   - Otherwise accumulate charge per plane
   - Build Eigen triplets: `{channel_index, time_slice_index, charge}`
3. Construct sparse matrices per plane via `setFromTriplets`
4. Also compute: estimated_minimum_charge (min across planes per blob, summed),
   estimated_total_charge, number_blobs, number_slices

#### `get_geom_clusters()` (lines 55-123)

Groups blobs by geometric connectivity:
1. Extract blob-only subgraph from cluster graph
2. Run `boost::connected_components` to find connected blob groups
3. Return map: group_id -> set of blob vertex descriptors

#### Coverage Judgment (`judge_coverage`, lines 365-448)

Compares two projections (ref vs tar) to determine subset relationships:
1. Create binary masks: pixel is "live" if charge > `-uncer_cut`
2. Check: both empty, one empty, or compare `ref_mask - tar_mask`
3. Uses threshold 0.01 (hardcoded, lines 400-403) for positive/negative detection
4. Returns: REF_COVERS_TAR, TAR_COVERS_REF, REF_EQ_TAR, BOTH_EMPTY, or OTHER

#### Alternative Coverage (`judge_coverage_alt`, lines 454-527)

More nuanced coverage with charge and count ratios:
1. Count live pixels and dead pixels in ref and tar
2. Count intersection pixels
3. If no intersection but both non-empty: OTHER
4. If intersection equals one set entirely: direct coverage
5. Otherwise: fractional comparison with configurable cut values:
   ```
   (1 - common_charge/small_charge) < min(cut[0] * (small+dead)/small, cut[1])
   AND
   (1 - common_counts/small_counts) < min(cut[2] * (small+dead)/small, cut[3])
   ```

---

## 3. ProjectionDeghosting (`src/ProjectionDeghosting.cxx`, 21K)

### Algorithm Overview

Global deghosting using 2D projection comparisons across clusters.

#### Main Flow

1. **Get geometric clusters**: Group blobs by connectivity
2. **Build projections**: For each cluster, compute 2D projections on all 3 planes
3. **Per-plane coverage analysis**: For each plane, compare projections between
   clusters using `judge_coverage()` or `judge_coverage_alt()`
4. **Remove covered clusters**: If a cluster's projection is covered by another's
   on enough planes, mark for removal
5. **Global cuts**: Apply chi-squared-like metrics combining time slices and
   charge per blob to filter remaining clusters

#### Configuration

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `nchan` | varies | Number of channels (matrix dimension) |
| `nslice` | varies | Number of time slices (matrix dimension) |
| `uncer_cut` | 1e11 | Threshold for dead pixel detection |
| `dead_default_charge` | -1e12 | Charge assigned to dead pixels |
| `global_deghosting_cut_nparas` | varies | Number of parameters per flag level |
| `global_deghosting_cut_values` | varies | Cut values for 3-view and 2-view clusters |
| `judge_alt_cut_values` | [4 values] | Tolerance cuts for alternative coverage |

---

## 4. ShadowGhosting (`src/ShadowGhosting.cxx`)

### Status: Incomplete / Pass-through

Currently creates blob and cluster shadow graphs via `BlobShadow::shadow()` and
`ClusterShadow::shadow()` but does not use them -- the input cluster is passed
through unchanged (line 39-40). The shadow computation appears to be scaffolding
for future implementation.

---

## Deghosting Strategy Summary

| Component | Scope | Method | When Used |
|-----------|-------|--------|-----------|
| InSliceDeghosting | Per-slice | Wire-score + adjacency | After initial clustering |
| ProjectionDeghosting | Global | 2D projection coverage | After charge solving |
| ShadowGhosting | Global | Cross-view shadows | Not yet functional |

The deghosting pipeline typically runs InSliceDeghosting first (local cleanup),
then ProjectionDeghosting (global cleanup after charge solving provides better
charge estimates).
