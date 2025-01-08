# clustering_separate Function Analysis

## Function Signature
```cpp
void clustering_separate(
    Grouping& live_grouping,
    std::map<int, std::pair<double, double>>& dead_u_index,
    std::map<int, std::pair<double, double>>& dead_v_index,
    std::map<int, std::pair<double, double>>& dead_w_index,
    const bool use_ctpc
)
```

## Purpose
Primary function for analyzing and separating potentially merged or crossed particle tracks based on geometric and topological criteria.

## Core Algorithm Flow

### 1. Initial Setup
```cpp
// Get and sort clusters
std::vector<Cluster*> live_clusters = live_grouping.children();
std::sort(live_clusters.begin(), live_clusters.end(), 
    [](const Cluster* cluster1, const Cluster* cluster2) {
        return cluster1->get_length() > cluster2->get_length();
    });

// Define reference directions
geo_point_t drift_dir(1, 0, 0);
geo_point_t beam_dir(0, 0, 1);
geo_point_t vertical_dir(0, 1, 0);

// Get time slice parameters
const auto& mp = live_grouping.get_params();
double live_time_slice_width = mp.nticks_live_slice * mp.tick_drift;

// Storage for processing results
std::vector<Cluster*> new_clusters;
std::vector<Cluster*> del_clusters;
```

### 2. Main Processing Loop
```cpp
for (size_t i = 0; i != live_clusters.size(); i++) {
    Cluster* cluster = live_clusters.at(i);
    
    // Only process large clusters
    if (cluster->get_length() > 100 * units::cm) {
        // Processing logic here
    }
}
```

### 3. Cluster Analysis Decision Tree
```cpp
if (cluster->get_length() > 100 * units::cm) {
    std::vector<geo_point_t> boundary_points;
    std::vector<geo_point_t> independent_points;

    // First decision criterion
    bool flag_proceed = JudgeSeparateDec_2(cluster, drift_dir, 
                                          boundary_points, 
                                          independent_points, 
                                          cluster->get_length());

    // Secondary analysis if first fails
    if (!flag_proceed && 
        cluster->get_length() > 100 * units::cm &&
        JudgeSeparateDec_1(cluster, drift_dir, cluster->get_length(), live_time_slice_width) &&
        independent_points.size() > 0) {
        // Additional analysis
    }
}
```

### 4. Topology Analysis
```cpp
bool flag_top = false;
for (size_t j = 0; j != independent_points.size(); j++) {
    if (independent_points.at(j).y() > 101.5 * units::cm) {
        flag_top = true;
        break;
    }
}

// Get main direction from PCA
geo_point_t main_dir(cluster->get_pca_axis(0).x(),
                     cluster->get_pca_axis(0).y(),
                     cluster->get_pca_axis(0).z());
```

### 5. Position-Based Analysis
For top region clusters (flag_top true):
```cpp
if (flag_top) {
    if (fabs(main_dir.angle(beam_dir) - 3.1415926/2.)/3.1415926*180. < 16 ||
        fabs(main_dir.angle(beam_dir) - 3.1415926/2.)/3.1415926*180. < 33 && 
            cluster->get_length() > 160*units::cm ||
        fabs(main_dir.angle(beam_dir) - 3.1415926/2.)/3.1415926*180. < 40 && 
            cluster->get_length() > 260*units::cm ||
        // Additional angle/length conditions...
        ) {
        flag_proceed = true;
    }
    else {
        // Secondary analysis for specific angles
        if (fabs(main_dir.angle(beam_dir) - 3.1415926/2.)/3.1415926*180. < 40 &&
            cluster->get_pca_value(1) > 0.2 * cluster->get_pca_value(0)) {
            // Try separation with larger distance cut
            std::vector<Cluster*> temp_sep_clusters = Separate_2(cluster, 10*units::cm);
            // Check results
        }
    }
}
```

For non-top region clusters:
```cpp
else {
    if (fabs(main_dir.angle(beam_dir) - 3.1415926/2.)/3.1415926*180. < 4 &&
            cluster->get_length() > 170*units::cm ||
        fabs(main_dir.angle(beam_dir) - 3.1415926/2.)/3.1415926*180. < 25 &&
            cluster->get_length() > 210*units::cm ||
        // Additional criteria...
        ) {
        flag_proceed = true;
    }
}
```


# Separation Execution Analysis

## Entry Point
```cpp
if (flag_proceed) {
    if (JudgeSeparateDec_1(cluster, drift_dir, cluster->get_length(), live_time_slice_width)) {
        // Main separation path
    }
    else if (cluster->get_length() < 6*units::m) {
        // Alternative separation path for shorter clusters
    }
}
```

## 1. Main Separation Path

### Initial Separation
```cpp
// Primary separation attempt
std::vector<Cluster*> sep_clusters = Separate_1(use_ctpc, cluster, 
    boundary_points, independent_points,
    dead_u_index, dead_v_index, dead_w_index,
    cluster->get_length());

// Process first separated cluster
Cluster* cluster1 = sep_clusters.at(0);
new_clusters.push_back(cluster1);
del_clusters.push_back(cluster);
```

### Multiple Cluster Handling
```cpp
if (sep_clusters.size() >= 2) {  // If more than one cluster produced
    // Add additional clusters (beyond first two)
    for (size_t k = 2; k < sep_clusters.size(); k++) {
        new_clusters.push_back(sep_clusters.at(k));
    }

    // Process second cluster
    std::vector<Cluster*> temp_del_clusters;
    Cluster* cluster2 = sep_clusters.at(1);
    double length_1 = cluster2->get_length();
    Cluster* final_sep_cluster = cluster2;
```

### Recursive Separation Level 1
```cpp
if (length_1 > 100 * units::cm) {
    boundary_points.clear();
    independent_points.clear();

    if (JudgeSeparateDec_1(cluster2, drift_dir, length_1, live_time_slice_width) &&
        JudgeSeparateDec_2(cluster2, drift_dir, boundary_points, independent_points, length_1)) {
        
        // Second level separation
        std::vector<Cluster*> sep_clusters = Separate_1(use_ctpc, cluster2, 
            boundary_points, independent_points,
            dead_u_index, dead_v_index, dead_w_index, length_1);
        
        Cluster* cluster3 = sep_clusters.at(0);
        new_clusters.push_back(cluster3);
        temp_del_clusters.push_back(cluster2);
```

### Recursive Separation Level 2
```cpp
        if (sep_clusters.size() >= 2) {
            // Add additional clusters from second separation
            for (size_t k = 2; k < sep_clusters.size(); k++) {
                new_clusters.push_back(sep_clusters.at(k));
            }

            Cluster* cluster4 = sep_clusters.at(1);
            final_sep_cluster = cluster4;
            length_1 = cluster4->get_length();

            // Third level check
            if (length_1 > 100 * units::cm) {
                boundary_points.clear();
                independent_points.clear();
```

### Recursive Separation Level 3
```cpp
                if (JudgeSeparateDec_1(cluster4, drift_dir, length_1, live_time_slice_width) &&
                    JudgeSeparateDec_2(cluster4, drift_dir, boundary_points, independent_points, length_1)) {
                    
                    std::vector<Cluster*> sep_clusters = Separate_1(use_ctpc, cluster4, 
                        boundary_points, independent_points,
                        dead_u_index, dead_v_index, dead_w_index, length_1);

                    Cluster* cluster5 = sep_clusters.at(0);
                    new_clusters.push_back(cluster5);
                    temp_del_clusters.push_back(cluster4);

                    if (sep_clusters.size() >= 2) {
                        for (size_t k = 2; k < sep_clusters.size(); k++) {
                            new_clusters.push_back(sep_clusters.at(k));
                        }
                        Cluster* cluster6 = sep_clusters.at(1);
                        final_sep_cluster = cluster6;
                    }
                    else {
                        final_sep_cluster = 0;
                    }
                }
            }
        }
```

### Final Cluster Processing
```cpp
if (final_sep_cluster != 0) {
    length_1 = final_sep_cluster->get_length();

    if (length_1 > 60 * units::cm) {
        boundary_points.clear();
        independent_points.clear();
        
        // Final separation attempt for long clusters
        if (JudgeSeparateDec_1(final_sep_cluster, drift_dir, length_1, live_time_slice_width) &&
            JudgeSeparateDec_2(final_sep_cluster, drift_dir, boundary_points, independent_points, length_1) &&
            independent_points.size() > 0) {
            
            std::vector<Cluster*> sep_clusters = Separate_1(use_ctpc, final_sep_cluster,
                boundary_points, independent_points,
                dead_u_index, dead_v_index, dead_w_index, length_1);

            Cluster* cluster5 = sep_clusters.at(0);
            new_clusters.push_back(cluster5);
            temp_del_clusters.push_back(final_sep_cluster);

            if (sep_clusters.size() >= 2) {
                // Process additional clusters
                for (size_t k = 2; k < sep_clusters.size(); k++) {
                    new_clusters.push_back(sep_clusters.at(k));
                }
                final_sep_cluster = sep_clusters.at(1);
            }
            else {
                final_sep_cluster = 0;
            }
        }
    }
```

### Final Separation Using Separate_2
```cpp
    if (final_sep_cluster != 0) {
        // Use simpler separation algorithm for final pass
        std::vector<Cluster*> final_sep_clusters = Separate_2(final_sep_cluster);
        for (auto it = final_sep_clusters.begin(); it != final_sep_clusters.end(); it++) {
            new_clusters.push_back(*it);
        }
        temp_del_clusters.push_back(final_sep_cluster);
    }
}
```

## 2. Alternative Path for Shorter Clusters
```cpp
else if (cluster->get_length() < 6*units::m) {
    std::vector<Cluster*> sep_clusters = Separate_1(use_ctpc, cluster,
        boundary_points, independent_points,
        dead_u_index, dead_v_index, dead_w_index,
        cluster->get_length());

    Cluster* cluster1 = sep_clusters.at(0);
    new_clusters.push_back(cluster1);
    del_clusters.push_back(cluster);

    if (sep_clusters.size() >= 2) {
        // Process additional clusters similar to main path
        // but with simplified logic
    }
}
```

## Key Features of Separation Execution

1. **Recursive Structure**
   - Up to 3 levels of recursive separation
   - Each level handles progressively smaller clusters
   - Different criteria at each level

2. **Length-Based Processing**
   - Primary threshold: 100 cm
   - Secondary threshold: 60 cm
   - Maximum length: 6 meters

3. **Memory Management**
   - Tracks clusters to be deleted
   - Manages temporary clusters
   - Maintains cluster hierarchy

4. **Protection Mechanisms**
   - Multiple validation checks
   - Size-based restrictions
   - Geometric criteria at each level

5. **Alternative Processing**
   - Special handling for shorter clusters
   - Simplified logic for certain cases
   - Final cleanup using Separate_2

This detailed separation execution shows how the function handles complex cluster configurations through multiple levels of recursive separation while maintaining cluster integrity and proper memory management.




### 7. Result Management
```cpp
// Update cluster collections
for (auto it = new_clusters.begin(); it != new_clusters.end(); it++) {
    Cluster* ncluster = (*it);
    live_clusters.push_back(ncluster);
}

// Remove processed clusters
for (auto it = del_clusters.begin(); it != del_clusters.end(); it++) {
    Cluster* ocluster = (*it);
    live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
}
```

## Key Decision Points

1. **Initial Size Filter**: Only processes clusters > 100 cm
2. **Top/Bottom Region**: Different criteria based on vertical position
3. **Angle Analysis**: Multiple angle thresholds relative to beam direction
4. **Length Considerations**: Various length thresholds for different conditions
5. **PCA-based Decisions**: Uses cluster shape characteristics

## Processing Thresholds

- **Minimum Cluster Size**: 100 cm
- **Top Region Boundary**: 101.5 cm
- **Length Thresholds**: Various (160 cm, 210 cm, 260 cm, etc.)
- **Angle Thresholds**: Multiple values (16°, 33°, 40°, etc.)
- **Maximum Length**: 6 meters

## Function Outcomes

1. **Cluster Separation**: Creates new separated clusters
2. **Cluster Cleanup**: Removes processed clusters
3. **Collection Updates**: Maintains live cluster collection
4. **Dead Region Handling**: Considers dead wire regions in separation

This function serves as the main coordinator for cluster separation, making decisions based on geometric and topological characteristics while maintaining the integrity of the detector's physical constraints.