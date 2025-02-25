# Understanding InSliceDeghosting in Wire-Cell Toolkit

## Overview

The `InSliceDeghosting` class is a component of the Wire-Cell Toolkit designed to identify and remove "ghost" signals in Time Projection Chamber (TPC) detectors. This document explains the purpose, methodology, and implementation of this class.

## What is Ghosting?

In wire chamber detectors like Liquid Argon Time Projection Chambers (LArTPCs), multiple wire planes are used to detect charge signals and reconstruct 3D positions. Each plane provides a 2D projection of the charge deposits. When multiple unrelated charge deposits occur in different locations, their wire plane projections can create false intersections, leading to "ghost" signals.

![Ghost Signals Illustration](https://raw.githubusercontent.com/WireCell/wire-cell-graphics/master/ghosting.png)

*Note: Image shows how multiple real charge deposits (green) can create false intersections (red) in the wire plane readouts.*

## Class Purpose

`InSliceDeghosting` performs "in-slice de-ghosting," which:

1. Identifies high-confidence "good" blobs based on charge thresholds
2. Assigns quality tags to blobs (good, bad, potential_good, potential_bad, to_be_removed)
3. Uses geometric consistency to determine which ambiguous blobs are likely ghosts
4. Removes identified ghost blobs from the cluster
5. Creates new blob-to-blob connections between remaining blobs

## Key Concepts

### Wire-Cell Detector Geometry

- **AnodePlane**: Represents the entire detection system
- **Face**: A region of the anode plane, typically with multiple wire planes
- **Wire Planes**: Each face has multiple wire planes (typically 3 in U, V, W orientations)
- **Blobs**: Detector response regions representing charge deposits
- **Clusters**: Groups of related blobs

### Quality Tags

The class uses a bitwise tagging system to mark blobs:

```cpp
enum blob_quality_bitpos {
    good,              // High confidence real blob
    bad,               // High confidence ghost blob
    potential_good,    // Likely real blob
    potential_bad,     // Likely ghost blob
    to_be_removed      // Marked for removal
};
```

## Algorithm Steps

The de-ghosting algorithm follows these general steps:

1. **Initial Quality Identification** (`blob_quality_ident`):
   - Tags blobs with charge above threshold as `good` and `potential_good`
   - Considers blob-to-blob connections to enhance quality assessment

2. **Local De-ghosting** (`local_deghosting` or `local_deghosting1`):
   - Groups blobs by number of active wire planes they appear in
   - Uses wire channel consistency across planes to identify ghosts
   - Tags low-confidence blobs for removal

3. **Geometric Clustering**:
   - Removes tagged ghost blobs
   - Creates new blob-to-blob edges within quality groups
   - Uses ray grid overlap logic to determine proximity

## Configuration Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `dryrun` | bool | If true, outputs original clusters without changes |
| `good_blob_charge_th` | double | Charge threshold for identifying "good" blobs |
| `clustering_policy` | string | Policy used for geometric clustering ("uboone", "simple", etc.) |
| `config_round` | int | Which algorithm variant to use (1, 2, or 3) |
| `deghost_th` | float | Threshold for de-ghosting decision logic |
| `deghost_th1` | float | Alternative threshold for second de-ghosting algorithm |

## Usage Example

```jsonnet
local wc = import "wirecell.jsonnet";

local deghosting = {
  type: "InSliceDeghosting",
  data: {
    good_blob_charge_th: 300.0,
    clustering_policy: "uboone",
    config_round: 1,
    deghost_th: 0.75,
    deghost_th1: 0.5,
  }
};
```

## Implementation Details

### Wire Plane Consistency

The algorithm uses wire channel crossings to determine when a blob is likely a ghost:

```
For each two-wire-plane blob:
  Check if adjacent to any good three-wire-plane blob
  If not adjacent to at least two good blobs:
    Calculate "score" based on wire channel consistency
    If score below threshold, mark as ghost
```

### Face and Wire Planes

The algorithm extensively uses the "face" concept:

1. Blobs store which face they belong to
2. Each face has wire planes (typically U, V, W)
3. Channels are organized by plane within each face
4. Spatial consistency is checked using the face's coordinate system

## Performance Considerations

- The algorithm uses bitwise operations for efficient tagging
- Connected component analysis identifies blob groups
- Geometric clustering operations can be computationally intensive
- The algorithm can be configured for different detector geometries

## References

- [Wire-Cell Documentation](https://wirecell.github.io/)
- [Cluster Shadow Documentation](https://github.com/wirecell/wire-cell-toolkit/blob/master/aux/docs/cluster-shadow.org)
- [LArTPC Reconstruction Techniques](https://arxiv.org/abs/1804.02583)