# Understanding ProjectionDeghosting in Wire-Cell Toolkit

## Introduction

The `ProjectionDeghosting` class is a critical component in the Wire-Cell Toolkit, designed to address a fundamental challenge in 3D reconstruction for Liquid Argon Time Projection Chambers (LArTPCs). This document explains the purpose, functionality, and implementation details of this class to help developers and physicists understand how it contributes to the overall reconstruction chain.

## The Ghosting Problem

In LArTPCs, charged particles create ionization tracks in liquid argon. These ionization electrons drift to wire planes that are oriented at different angles. Each wire plane provides a 2D projection of the particle trajectory from a different perspective. When attempting to reconstruct 3D positions from these 2D projections, ambiguities can arise:

- Multiple possible 3D points may be consistent with the same set of wire signals
- These ambiguities create "ghost" hits - false intersections that don't correspond to actual energy deposits
- Ghosts can significantly degrade reconstruction quality and create spurious features in event displays

![Ghosting Illustration](https://i.imgur.com/Lm1YLFD.png)
*Conceptual illustration: Left - true particle trajectory. Right - ambiguous reconstructions including ghosts.*

## How ProjectionDeghosting Works

The `ProjectionDeghosting` class implements a sophisticated algorithm to identify and remove ghost hits by analyzing the 2D projections of reconstructed 3D clusters. Its operation can be broken down into these key stages:

### 1. Blob and Cluster Shadows

The algorithm first constructs specialized graph representations:

- **BlobShadow Graph**: Represents relationships between "blobs" (3D reconstructed charge deposits)
- **ClusterShadow Graph**: Represents connections between clusters of blobs
  
These graphs capture how blobs and clusters relate to each other across different wire plane views.

### 2. 2D Projections

For each cluster, the algorithm creates 2D projections for each wire plane layer:

```cpp
projection2d::layerprojection2dmap& proj_cluster = get_projection(
    id2lproj, cs_cluster, in_graph, b_cluster, m_nchan, m_nslice, m_uncer_cut, m_dead_default_charge);
```

These projections represent how each 3D cluster would appear when viewed from each wire plane's perspective.

### 3. Coverage Analysis

The algorithm then analyzes the "coverage" between different projections:

```cpp
coverage_alt = projection2d::judge_coverage_alt(proj2d_comp_3dclus, 
    proj2d_clust3d, m_judge_alt_cut_values, m_uncer_cut);
```

Coverage analysis examines how projections from different clusters overlap or contain each other. Specific patterns of coverage provide strong indicators of whether a cluster is real or a ghost.

### 4. Ghost Identification

Based on the coverage analysis, charge distribution, and other metrics, the algorithm makes decisions about which clusters are likely to be ghosts:

```cpp
if(sqrt(pow(n_timeslices / m_global_deghosting_cut_values.at(0), 2) + 
   pow(min_charge / n_blobs / m_global_deghosting_cut_values.at(1), 2)) < 1 || 
   min_charge / n_blobs / m_global_deghosting_cut_values.at(2) < 1.) {
    saved = 0;  // Mark as ghost
} else {
    saved = 1;  // Keep
}
```

The algorithm applies configurable thresholds to make these determinations, considering:
- Time extent of the cluster
- Charge distribution
- Projection coverage patterns
- Blob count and density

### 5. Cluster Pruning

Finally, the algorithm removes clusters identified as ghosts:

```cpp
auto out_graph = remove_blobs(in_graph, tagged_bs, true);
```

This produces a refined cluster graph with ghost hits removed, resulting in a more accurate 3D reconstruction.

## Key Concepts

### AnodePlane and Face

In Wire-Cell:
- An `AnodePlane` represents a physical detector component with multiple wire planes
- Each `AnodePlane` can have one or more "faces" (typically two, one on each side)
- A "face" (`ianodeface`) is a collection of wire planes that share the same sensitive volume

The "face" concept is crucial for deghosting because:
1. Each face contains wire planes at different angles (typically U, V, and W views)
2. Particles passing through the detector create signals on wires from each plane of the same face
3. Ghosts can only form between wires from the same face

The algorithm tracks which face each blob belongs to:
```cpp
auto iblob = get<cluster_node_t::blob_t>(cg[vtx].ptr);
face2blobs[iblob->face()].push_back(vdesc);
```

### Projections and Coverage

The algorithm uses several types of coverage relationships between projections:
- **ref_covers_tar**: One projection fully covers another
- **tar_covers_ref**: One projection is fully covered by another
- **ref_eq_tar**: Projections are equivalent
- **both_empty**: Both projections are empty
- **other**: Other relationship patterns

These coverage relationships form distinctive patterns that help distinguish real tracks from ghost hits.

## Configuration Parameters

The `ProjectionDeghosting` class uses several configurable parameters:

```cpp
cfg["verbose"] = m_verbose;
cfg["nchan"] = (unsigned int) m_nchan;
cfg["nslice"] = (unsigned int) (m_nslice);
cfg["dryrun"] = m_dryrun;
cfg["global_deghosting_cut_nparas"] = m_global_deghosting_cut_nparas;
```

Key parameters include:
- **nchan/nslice**: Dimensions of the detector readout
- **global_deghosting_cut_values**: Values controlling ghost identification thresholds
- **judge_alt_cut_values**: Parameters for the coverage judgment algorithm
- **uncer_cut**: Threshold for handling measurement uncertainties
- **dead_default_charge**: Value to use for dead/inactive channels

## Integration with Wire-Cell

The `ProjectionDeghosting` class integrates with the broader Wire-Cell reconstruction chain:

1. It implements the `iclusterfilter` interface, allowing it to be used as a filter in reconstruction pipelines
2. It takes cluster graphs as input and produces filtered cluster graphs as output
3. It can be configured through Wire-Cell's JSON configuration system
4. It's typically positioned in the reconstruction chain after initial clustering but before 3D point extraction

## Performance Considerations

The effectiveness of projection deghosting depends on several factors:

- **Detector geometry**: The angles between wire planes affect ghosting patterns
- **Signal quality**: Cleaner signals lead to better disambiguation
- **Tuning parameters**: The thresholds need to be tuned for specific detectors
- **Computational cost**: The algorithm performs significant graph operations and may be computationally intensive

## Example Usage

In a Wire-Cell configuration:

```json
{
  "configs": [
    {
      "data": {
        "verbose": false,
        "nchan": 8256,
        "nslice": 9592,
        "global_deghosting_cut_nparas": 3,
        "global_deghosting_cut_values": [3.0, 3000.0, 2000.0, 8.0, 8000.0, 4000.0, 8.0, 8000.0, 6000.0],
        "judge_alt_cut_values": [0.05, 0.33, 0.15, 0.33]
      },
      "name": "ProjectionDeghosting"
    }
  ]
}
```

## Summary

The `ProjectionDeghosting` class provides a sophisticated solution to the ghosting problem in LArTPC reconstruction. By analyzing 2D projections of 3D clusters and applying geometric constraints, it can effectively distinguish real particle trajectories from ghost hits, leading to significantly improved 3D reconstruction quality.

Understanding this class is essential for anyone working with Wire-Cell reconstruction, particularly when dealing with complex event topologies where ghosting can severely impact reconstruction performance.

## References

1. Wire-Cell Toolkit documentation
2. LArTPC reconstruction techniques
3. MicroBooNE and DUNE reconstruction papers
4. Wire-Cell developer notes on clustering algorithms