# Simple3DPointCloud Class Analysis

## Core Purpose
The Simple3DPointCloud class is designed to manage and perform operations on 3D point cloud data with efficient spatial queries using a k-d tree data structure. It's part of the WireCell::PointCloud::Facade namespace and appears to be used in particle physics detector data processing.

## Class Structure

### Key Member Variables
1. `m_points` (points_type)
   - A 3-dimensional vector structure storing point coordinates
   - Organized as three separate arrays for x, y, z coordinates
   - Uses a columnar data layout for better memory efficiency

2. `m_kd` (std::unique_ptr<nfkd_t>)
   - A lazy-initialized k-d tree for spatial queries
   - Uses NFKDVec::Tree with dynamic indexing
   - Mutable to allow lazy initialization in const methods

### Important Type Definitions
```cpp
using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexDynamic>;
using points_type = nfkd_t::points_type;
using results_type = nfkd_t::results_type;
using point_type = std::vector<double>;
```

## Key Algorithms and Operations

### 1. Point Addition (add method)
```cpp
void add(const point_type& new_pt) {
    // Validate 3D point
    if (new_pt.size() != 3) {
        raise<ValueError>("points must be 3D");
    }
    // Add coordinates to respective arrays
    for (size_t ind=0; ind<3; ++ind) {
        points()[ind].push_back(new_pt[ind]);
    }
    // Update k-d tree
    kd().append({{new_pt[0]}, {new_pt[1]}, {new_pt[2]}});
}
```

### 2. K-D Tree Management
- Lazy initialization pattern used for the k-d tree
- Two versions (const and non-const) of the kd() method
- Tree can be rebuilt on demand using the rebuild parameter

### 3. Nearest Neighbor Search Operations

#### A. Closest Point Index Query
```cpp
results_type get_closest_index(const geo_point_t& p, const size_t N) const {
    return kd().knn(N, p);  // K-nearest neighbors search
}
```

#### B. Closest Point Search
```cpp
std::pair<size_t, geo_point_t> get_closest_wcpoint(const geo_point_t& p) const {
    const auto knn_res = kd().knn(1, p);
    // Returns index and actual point coordinates
}
```

#### C. Directional Search (get_closest_point_along_vec)
Searches for points along a specified direction with constraints:
- Starting point (p_test1)
- Direction vector (dir)
- Search distance (test_dis)
- Step size (dis_step)
- Angular cutoff (angle_cut)
- Distance cutoff (dis_cut)

## Key Features

1. **Efficiency**
   - Columnar data storage for better memory access patterns
   - Lazy initialization of k-d tree
   - Efficient spatial queries using k-d tree structure

2. **Flexibility**
   - Supports both exact and approximate nearest neighbor searches
   - Allows directional searches with constraints
   - Supports multiple point query methods

3. **Safety**
   - Input validation for 3D points
   - Error handling for invalid operations
   - Const-correctness for thread safety

4. **Memory Management**
   - Smart pointer usage for k-d tree
   - Automatic cleanup through RAII
   - Efficient memory usage through columnar storage

## Common Use Patterns

1. **Building Point Cloud**
```cpp
Simple3DPointCloud cloud;
cloud.add({x, y, z});  // Add points one at a time
```

2. **Spatial Queries**
```cpp
// Find N nearest neighbors
auto neighbors = cloud.get_closest_index(point, N);

// Find single nearest point
auto [index, point] = cloud.get_closest_wcpoint(query_point);
```

3. **Directional Searches**
```cpp
auto [index, distance] = cloud.get_closest_point_along_vec(
    start_point, direction, 
    search_distance, step_size,
    angle_cutoff, distance_cutoff
);
```

## Implementation Notes

1. The class uses a columnar data structure instead of an array of structs, which can provide better cache performance for certain operations.

2. The k-d tree is lazily initialized and can be rebuilt on demand, allowing for efficient updates when the point cloud changes.

3. The class provides both exact (knn) and constrained (directional) search capabilities.

4. Error handling is implemented for invalid inputs and operations.

5. The class supports integration with other components through its template-based closest points comparison functionality.

# Multi2DPointCloud Class Analysis

## Core Purpose
The Multi2DPointCloud class is designed to manage and perform operations on 2D projections of 3D point clouds in three different planes (u, v, w), commonly used in particle physics wire chamber detectors. It maintains separate 2D point clouds for each plane and provides efficient spatial queries using k-d trees.

## Class Structure

### Key Member Variables
1. `m_points[3]` (points_type[3])
   - Array of 3 two-dimensional vector structures storing point coordinates
   - Each points_type stores x and y coordinates for a plane
   - Uses columnar data layout for efficiency

2. `m_kd[3]` (std::unique_ptr<nfkd_t>[3])
   - Array of 3 lazy-initialized k-d trees for spatial queries
   - One tree per plane (u, v, w)
   - Mutable to allow lazy initialization in const methods

3. `angle_uvw[3]` (double[3])
   - Array storing the angles for each plane (u, v, w)
   - Used for projecting 3D points onto 2D planes

### Important Type Definitions
```cpp
using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexDynamic>;
using coordinates_type = nfkd_t::coordinates_type;
using points_type = nfkd_t::points_type;
using results_type = nfkd_t::results_type;
using point_type = std::vector<double>;
```

## Key Algorithms and Operations

### 1. Initialization
```cpp
Multi2DPointCloud(double angle_u, double angle_v, double angle_w) 
    : angle_uvw{angle_u, angle_v, angle_w} {
    for (size_t plane = 0; plane < 3; ++plane) {
        points(plane).resize(2);  // 2D points need only x,y coordinates
    }
}
```

### 2. Point Addition and Projection
```cpp
void add(const geo_point_t& new_pt) {
    for (size_t plane = 0; plane < 3; ++plane) {
        // Project 3D point onto 2D plane using rotation matrix
        double x = new_pt[0];  // x coordinate remains unchanged
        // Calculate y coordinate using rotation matrix
        double y = cos(angle_uvw[plane]) * new_pt[2] - 
                  sin(angle_uvw[plane]) * new_pt[1];
        
        // Store projected coordinates
        points(plane)[0].push_back(x);
        points(plane)[1].push_back(y);
        // Update k-d tree for the plane
        kd(plane).append({{x}, {y}});
    }
}
```

### 3. K-D Tree Management
Maintains separate k-d trees for each plane with lazy initialization:
```cpp
const nfkd_t& kd(const size_t plane, const bool rebuild=false) const {
    if (rebuild) m_kd[plane] = nullptr;
    if (m_kd[plane]) return *m_kd[plane];
    m_kd[plane] = std::make_unique<nfkd_t>(2);
    return *m_kd[plane];
}
```

### 4. Spatial Query Operations

#### A. Closest Point Distance Query
```cpp
std::pair<int, double> get_closest_2d_dis(const geo_point_t& p, size_t plane) const {
    // Project 3D point to 2D
    double x = p[0];
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    
    // Perform k-nearest neighbor search
    const auto& res = kd(plane).knn(1, {x, y});
    
    // Return index and distance
    if (res.size() == 1)
        return std::make_pair(res[0].first, sqrt(res[0].second));
    else
        return std::make_pair(-1, 1e9);
}
```

#### B. Radius Search
```cpp
std::vector<std::pair<size_t, double>> get_closest_2d_index_radius(
    const geo_point_t& p, const double radius, size_t plane) const {
    // Project point to 2D
    double x = p[0];
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    
    // Search within radius
    const auto& res = kd(plane).radius(radius * radius, {x, y});
    
    // Convert results
    std::vector<std::pair<size_t, double>> ret;
    for (const auto& r : res) {
        ret.push_back(std::make_pair(r.first, sqrt(r.second)));
    }
    return ret;
}
```

#### C. K-Nearest Neighbors Search
```cpp
std::vector<std::pair<size_t, double>> get_closest_2d_index_knn(
    const geo_point_t& p, const int N, size_t plane) const {
    // Project point to 2D
    double x = p[0];
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    
    // Find N nearest neighbors
    const auto& res = kd(plane).knn(N, {x, y});
    
    // Convert results
    std::vector<std::pair<size_t, double>> ret;
    for (const auto& r : res) {
        ret.push_back(std::make_pair(r.first, r.second));
    }
    return ret;
}
```

## Key Features

1. **Multi-Plane Management**
   - Maintains separate 2D projections for each plane
   - Independent k-d trees for efficient spatial queries in each plane
   - Consistent projection transformations

2. **Efficient Projection**
   - Uses rotation matrices for 3D to 2D projection
   - Preserves x-coordinate while transforming y-z plane
   - Optimized for wire chamber geometry

3. **Flexible Querying**
   - Supports nearest neighbor searches
   - Implements radius-based searches
   - Provides k-nearest neighbors queries
   - Returns both indices and distances

4. **Memory Efficiency**
   - Columnar data storage
   - Lazy k-d tree initialization
   - Smart pointer management

## Common Use Patterns

1. **Initialization and Point Addition**
```cpp
Multi2DPointCloud cloud(angle_u, angle_v, angle_w);
cloud.add(point3d);  // Automatically projects to all planes
```

2. **Spatial Queries Per Plane**
```cpp
// Find closest point in a plane
auto [index, distance] = cloud.get_closest_2d_dis(point3d, plane);

// Find points within radius
auto neighbors = cloud.get_closest_2d_index_radius(point3d, radius, plane);

// Find k nearest neighbors
auto knn = cloud.get_closest_2d_index_knn(point3d, k, plane);
```

## Implementation Notes

1. The class uses a projection scheme specific to wire chamber detector geometry where each plane represents a different wire orientation.

2. The projection preserves the x-coordinate and transforms the y-z coordinates based on the plane angle.

3. Each plane maintains its own k-d tree for efficient spatial queries in the projected 2D space.

4. The class provides comprehensive error handling and boundary checking.

5. The implementation is optimized for the specific needs of particle physics detector data processing.


# DynamicPointCloud Class Analysis

## Core Purpose
The DynamicPointCloud class is a sophisticated data structure that combines 3D point cloud management with 2D projections, specifically designed for particle physics detector track analysis. It maintains both 3D and 2D representations of points while associating them with clusters and blobs, making it particularly useful for particle track reconstruction.

## Class Structure

### Key Member Variables
1. `m_pc2d` (Multi2DPointCloud)
   - Handles 2D projections in three planes (u, v, w)
   - Manages wire plane geometry and projections

2. `m_pc3d` (Simple3DPointCloud)
   - Stores and manages the full 3D point cloud data
   - Provides 3D spatial queries

3. `m_winds[3]` (std::vector<int>[3])
   - Stores wire indices for each plane (u, v, w)
   - Used for detector readout mapping

4. `m_clusters` (std::vector<const Cluster*>)
   - Stores pointers to associated cluster objects
   - Maps points to their parent clusters

5. `m_blobs` (std::vector<const Blob*>)
   - Stores pointers to associated blob objects
   - Contains charge deposit information

### Important Type Definitions
```cpp
using points3d_type = Simple3DPointCloud::points_type;
using points2d_type = Multi2DPointCloud::points_type;
using point_type = std::vector<double>;
```

## Key Algorithms and Operations

### 1. Point Addition Methods

#### A. Standard Point Addition
```cpp
void add_points(const Cluster* cluster, const int flag=0, 
                const double step = 0.6*units::cm) {
    size_t current_size = get_num_points();
    const auto& winds = cluster->wire_indices();

    if (flag == 0) {
        // Add actual points from cluster
        for (size_t i = 0; i != cluster->npoints(); i++) {
            // Store cluster reference
            m_clusters.push_back(cluster);
            
            // Add 3D point
            m_pc3d.add({cluster->point3d(i).x(), 
                       cluster->point3d(i).y(), 
                       cluster->point3d(i).z()});
            
            // Add 2D projections
            m_pc2d.add(cluster->point3d(i));
            
            // Store wire indices
            for (size_t plane = 0; plane < 3; ++plane) {
                m_winds[plane].push_back(winds[plane][i]);
            }
            
            // Store blob reference
            m_blobs.push_back(cluster->blob_with_point(i));
        }
    }
    else {
        // Add skeleton points with interpolation
        const std::list<size_t>& path_wcps = cluster->get_path_wcps();
        
        // Interpolate points along path
        geo_point_t prev_wcp = cluster->point3d(path_wcps.front());
        for (auto it = path_wcps.begin(); it != path_wcps.end(); it++) {
            geo_point_t test_point = cluster->point3d(*it);
            double dis = (test_point - prev_wcp).magnitude();
            
            if (dis <= step) {
                // Add point directly if close enough
                points.push_back(test_point);
            }
            else {
                // Interpolate points along segment
                int num_points = int(dis / step) + 1;
                for (int k = 0; k != num_points; k++) {
                    double t = (k + 1.0) / num_points;
                    geo_point_t current_pt = prev_wcp + 
                        (test_point - prev_wcp) * t;
                    points.push_back(current_pt);
                }
            }
            prev_wcp = test_point;
        }
    }
}
```

#### B. Directional Point Addition
```cpp
void add_points(const Cluster* cluster, const geo_point_t& p_test,
                const geo_point_t& dir_unmorm, const double range,
                const double step, const double angle) {
    geo_point_t dir = dir_unmorm.norm();
    int num_points = int(range / step) + 1;
    
    for (int k = 0; k != num_points; k++) {
        // Calculate distance cut based on angle
        double dis_cut = std::min(
            std::max(2.4 * units::cm, 
                    k * step * sin(angle / 180. * 3.1415926)), 
            13 * units::cm);
        
        // Add point with calculated position
        m_clusters.push_back(cluster);
        m_blobs.push_back(nullptr);
        
        geo_point_t new_point = p_test + dir * k * step;
        m_pc3d.add(new_point);
        m_pc2d.add(new_point);
        
        // Store distance cut as wire index
        for (int plane = 0; plane < 3; ++plane) {
            m_winds[plane].push_back(int(dis_cut));
        }
    }
}
```

### 2. Spatial Query Operations

#### A. 2D Point Information Retrieval
```cpp
std::vector<std::tuple<double, const Cluster*, size_t>> 
get_2d_points_info(const geo_point_t& p, const double radius, 
                   const int plane) {
    // Get points within radius in specified plane
    auto results = m_pc2d.get_closest_2d_index_radius(p, radius, plane);
    
    // Build result tuples with distance, cluster, and index
    std::vector<std::tuple<double, const Cluster*, size_t>> return_results;
    for (const auto& [index, distance] : results) {
        return_results.push_back(std::make_tuple(
            distance,
            m_clusters.at(index),
            index
        ));
    }
    return return_results;
}
```

### 3. Advanced Analysis Methods

#### A. Hough Transform
```cpp
std::pair<double, double> hough_transform(const geo_point_t& origin, 
                                        const double dis) const {
    // Collect points within distance
    std::vector<geo_point_t> pts;
    std::vector<const Blob*> blobs;
    auto results = m_pc3d.kd().radius(dis * dis, origin);
    
    // Build histogram in parameter space
    auto hist = make_histogram(...);
    
    for (const auto& [point_index, _] : results) {
        const auto* blob = m_blobs[point_index];
        auto charge = blob->charge();
        if (charge <= 0) continue;
        
        const auto& pt = m_pc3d.point(point_index);
        const Vector dir = (pt - origin).norm();
        const double r = (pt - origin).magnitude();
        
        // Calculate parameters and fill histogram
        const double p1 = theta_param(dir);
        const double p2 = phi_param(dir);
        
        // Weight based on distance and charge
        double weight = charge / blob->npoints();
        if (r >= 10 * units::cm) {
            weight *= pow(10 * units::cm / r, 2);
        }
        
        hist(p1, p2, weight);
    }
    
    // Find maximum bin
    auto max_bin = find_maximum_bin(hist);
    return {max_bin.center(0), max_bin.center(1)};
}
```

## Key Features

1. **Dual Representation**
   - Maintains both 3D and 2D point clouds
   - Efficiently handles projections and transformations
   - Preserves relationships between different views

2. **Cluster Association**
   - Associates points with physics clusters
   - Maintains blob information for charge analysis
   - Supports track reconstruction

3. **Flexible Point Addition**
   - Supports direct point addition
   - Provides path interpolation
   - Allows directional point generation

4. **Advanced Analysis**
   - Implements Hough transform for track finding
   - Supports various spatial queries
   - Handles charge-weighted calculations

## Common Use Patterns

1. **Cluster Processing**
```cpp
DynamicPointCloud cloud(angle_u, angle_v, angle_w);
cloud.add_points(cluster);  // Add all points from cluster
```

2. **Track Finding**
```cpp
// Find track direction using Hough transform
auto [theta, phi] = cloud.hough_transform(origin, search_radius);
```

3. **Spatial Analysis**
```cpp
// Find nearby points in 2D projection
auto points = cloud.get_2d_points_info(point, radius, plane);
```

## Implementation Notes

1. The class combines Simple3DPointCloud and Multi2DPointCloud for comprehensive spatial analysis.

2. Point addition methods handle both direct points and interpolated paths.

3. The Hough transform implementation is optimized for particle track finding.

4. Wire indices are used to map points to detector readout channels.

5. The implementation supports charge-weighted analysis for better track reconstruction.

