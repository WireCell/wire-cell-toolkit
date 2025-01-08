# Facade_Cluster Class Documentation

## Overview
The `Facade_Cluster` class is part of the WireCell Point Cloud namespace and provides a high-level interface for working with clusters of 3D points and blobs in particle detector data analysis. It acts as a facade over a PC (Point Cloud) tree, giving semantics to what would otherwise be simple nodes.

## Key Features

### Point Cloud Management
- Supports both 2D and 3D point clouds
- Provides k-d tree functionality for efficient spatial queries
- Manages collections of blobs (groups of points) within the cluster

### Spatial Operations
- Distance calculations between points and clusters
- Finding nearest neighbors and points within a radius
- Calculation of geometric properties (center, PCA axes)
- Convex hull computation
- Hough transformations for direction finding

### Graph Operations
- Supports creation and manipulation of cluster graphs
- Implements Dijkstra's shortest path algorithm
- Provides connectivity analysis between blobs

### Time-Based Analysis
- Manages time-slice based organization of blobs
- Supports time-ordered operations and queries
- Tracks temporal relationships between points and blobs

## Complete API Reference

### Core Class Management

#### Constructor & Destructor
```cpp
Cluster();
virtual ~Cluster();
```
- Default constructor initializes an empty cluster
- Virtual destructor ensures proper cleanup of derived classes

#### Grouping Access
```cpp
Grouping* grouping();
const Grouping* grouping() const;
```
- Returns pointer to the parent grouping that contains this cluster
- Used to access parameters and higher-level organization
- Const version available for read-only access

### Point Cloud Operations

#### 3D Point Cloud Management
```cpp
const sv3d_t& sv3d() const;
```
- Returns the scoped view for the "3d" point cloud (x,y,z coordinates)
- Used internally for point cloud management
- Provides access to underlying data structure

```cpp
const kd3d_t& kd3d() const;
```
- Returns the k-d tree for "3d" point cloud
- May trigger k-d tree building if not already constructed
- Essential for efficient spatial queries

```cpp
const points_type& points() const;
```
- Returns full array of point coordinates
- Recommended for bulk point access
- More efficient than accessing individual points repeatedly

#### 2D Point Cloud Management
```cpp
const sv2d_t& sv2d(const size_t plane) const;
```
- Returns 2D scoped view for specified plane (0=U, 1=V, 2=W)
- Provides access to 2D projections for wire plane analysis
- Essential for wire-plane specific operations

```cpp
const kd2d_t& kd2d(const size_t plane) const;
```
- Returns 2D k-d tree for specified plane
- Enables efficient spatial queries in 2D projections
- Used for wire-plane specific distance calculations

### Point Access and Queries

#### Point Information
```cpp
geo_point_t point3d(size_t point_index) const;
geo_point_t point(size_t point_index) const;
```
- Returns 3D point at given k-d tree index
- point() is alias for point3d() to match Simple3DPointCloud interface
- Use with caution in tight loops - prefer bulk access through points()

```cpp
int npoints() const;
```
- Returns total number of points in cluster
- Uses cached value when available
- Thread-safe for read access

```cpp
size_t nbpoints() const;
```
- Returns total number of points according to sum of Blob::nbpoints()
- Equivalent to WCP's get_num_points()
- May differ from npoints() due to different counting methods

#### Spatial Queries
```cpp
kd_results_t kd_radius(double radius, const geo_point_t& query_point) const;
```
- Performs radius search in k-d tree
- Returns points within specified radius of query point
- Note: radius parameter is linear distance, not squared

```cpp
kd_results_t kd_knn(int nnearest, const geo_point_t& query_point) const;
```
- Performs k-nearest neighbor search
- Returns the specified number of closest points to query point
- Efficient for finding local point neighborhoods

```cpp
std::vector<geo_point_t> kd_points(const kd_results_t& res);
std::vector<geo_point_t> kd_points(const kd_results_t& res) const;
```
- Converts k-d tree query results to vector of points
- Both mutable and const versions available
- Useful for processing query results

#### Closest Point Finding
```cpp
std::pair<geo_point_t, const Blob*> get_closest_point_blob(const geo_point_t& point) const;
```
- Returns closest point and its containing blob to given point
- Useful for associating external points with cluster structure
- Returns {point, nullptr} if no points found

```cpp
std::pair<size_t, geo_point_t> get_closest_wcpoint(const geo_point_t& p) const;
```
- Returns index and coordinates of closest point to given point
- Equivalent to WCP's get_closest_wcpoint functionality
- Returns {-1, nullptr} if no points found

```cpp
size_t get_closest_point_index(const geo_point_t& point) const;
```
- Returns k-d tree index of closest point to given point
- Throws ValueError if cluster is empty
- More efficient than get_closest_wcpoint when only index is needed

```cpp
double get_closest_dis(const geo_point_t& point) const;
```
- Returns distance to closest point from given point
- Throws ValueError if cluster is empty
- Efficient when only distance is needed

```cpp
std::vector<size_t> get_closest_2d_index(
    const geo_point_t& p, 
    const double search_radius, 
    const int plane) const;
```
- Finds indices of points within search radius in 2D projection
- Operates on specified wire plane (0=U, 1=V, 2=W)
- Essential for wire-plane specific analysis

```cpp
template <typename PCType>
std::tuple<int, int, double> get_closest_points(const PCType& two) const;
```
- Finds closest points between this cluster and another point cloud
- Returns tuple containing:
  - Index of closest point in this cluster
  - Index of closest point in the other point cloud
  - Distance between these points
- Uses iterative search to find global minimum distance
- Checks multiple starting positions to avoid local minima
- Template allows comparison with any point cloud type implementing required interface

```cpp
std::pair<geo_point_t, double> get_closest_point_along_vec(
    geo_point_t& p_test,  // Starting point
    geo_point_t dir,      // Direction vector
    double test_dis,      // Maximum search distance
    double dis_step,      // Step size for search
    double angle_cut,     // Maximum allowed angle deviation
    double dis_cut        // Maximum allowed distance from line
) const;
```
- Searches for closest point along a specified direction vector
- Parameters:
  - p_test: Starting point for the search
  - dir: Direction to search along
  - test_dis: Maximum distance to search
  - dis_step: Distance between test points along direction
  - angle_cut: Maximum allowed angular deviation (in degrees)
  - dis_cut: Maximum allowed perpendicular distance from search line
- Returns:
  - Closest point found meeting criteria
  - Distance from starting point
- Used for tracking trajectory analysis
- Important for finding continuation points in tracks
- Includes angular constraints to ensure consistent direction

#### Point Counting and Analysis
```cpp
int nnearby(const geo_point_t& point, double radius) const;
```
- Counts points within radius of given point
- Uses linear distance measure (not squared)
- Efficient for density estimation

```cpp
std::pair<int, int> ndipole(
    const geo_point_t& point, 
    const geo_point_t& dir, 
    const double dis=-1) const;
```
- Returns count of points in forward/backward direction from point
- Optional distance cutoff (disabled if negative)
- Useful for analyzing point cloud directionality

### Blob Management

#### Blob Access
```cpp
std::vector<Blob*> kd_blobs();
std::vector<const Blob*> kd_blobs() const;
```
- Returns all blobs in k-d tree order
- Order differs from children() order and sort_blobs() order
- Both mutable and const versions available

```cpp
Blob* blob_with_point(size_t point_index);
const Blob* blob_with_point(size_t point_index) const;
```
- Returns blob containing the point at given k-d tree index
- Essential for mapping between points and their containing blobs
- Both mutable and const versions available

```cpp
std::vector<Blob*> blobs_with_points(const kd_results_t& res);
std::vector<const Blob*> blobs_with_points(const kd_results_t& res) const;
```
- Returns blobs containing points from k-d tree query results
- Maintains order of input results
- Useful for processing spatial query results

```cpp
std::vector<int> get_blob_indices(const Blob*) const;
```
- Returns vector of point indices belonging to specified blob
- Uses lazy initialization for efficiency
- Important for blob-point relationship mapping

#### Blob Information
```cpp
void print_blobs_info() const;
```
- Prints detailed information about all blobs in cluster
- Outputs wire index ranges (U, V, W) and time slice indices
- Useful for debugging and analysis

```cpp
const Blob* get_first_blob() const;
```
- Returns blob at earliest time
- Throws ValueError if cluster is empty
- Based on time_blob_map ordering

```cpp
const Blob* get_last_blob() const;
```
- Returns blob at latest time
- Throws ValueError if cluster is empty
- Based on time_blob_map ordering

```cpp
size_t get_num_time_slices() const;
```
- Returns number of unique time slices in cluster
- Based on time_blob_map size
- Important for temporal analysis

#### Blob Relationships
```cpp
std::vector<const Blob*> is_connected(const Cluster& c, const int offset) const;
```
- Determines connectivity between this cluster and another
- Offset parameter controls connection tolerance
- Returns vector of connecting blobs

```cpp
const_blob_point_map_t get_closest_blob(
    const geo_point_t& point, 
    double radius) const;
```
- Returns map of blobs and their closest points within radius
- Each returned blob has at least one point within radius
- Useful for analyzing local blob structure

### Geometric Analysis

#### Position and Direction
```cpp
geo_point_t calc_ave_pos(const geo_point_t& origin, const double dis) const;
```
- Calculates charge-weighted average position of nearby points
- Uses points within specified distance of origin
- Important for smoothing and local averaging

```cpp
std::pair<geo_point_t,geo_point_t> get_two_extreme_points() const;
```
- Finds most distant pair of points in cluster
- Uses local averaging for stability
- Important for determining cluster extent

```cpp
std::pair<geo_point_t, geo_point_t> get_highest_lowest_points(size_t axis = 1) const;
```
- Returns points at extremes of given Cartesian axis
- Default is Y-axis (axis=1)
- Returns points in descending order

```cpp
std::pair<geo_point_t, geo_point_t> get_earliest_latest_points() const;
```
- Returns points at extremes of X-axis (drift direction)
- Returns points in ascending order
- Important for drift time analysis

```cpp
std::pair<geo_point_t, geo_point_t> get_front_back_points() const;
```
- Returns points at extremes of Z-axis
- Used for longitudinal extent analysis
- Important for track reconstruction

#### Shape Analysis
```cpp
std::vector<geo_point_t> get_hull() const;
```
- Computes convex hull of cluster points
- Uses QuickHull algorithm
- Important for shape analysis and visualization

```cpp
double get_length() const;
```
- Returns geometric size of cluster
- Based on transverse extents and time
- Uses cached value for efficiency

```cpp
std::tuple<int, int, int, int> get_uvwt_range() const;
std::tuple<int, int, int, int> get_uvwt_min() const;
std::tuple<int, int, int, int> get_uvwt_max() const;
```
- Return wire indices and time ranges
- Important for detector coordinate analysis
- Provide min/max values for U, V, W coordinates and time

```cpp
geo_point_t get_center() const;
```
- Returns the geometric center (centroid) of the cluster
- Calculated as average position of all points
- Uses charge weighting if available
- Cached for efficiency after first calculation
- Important for PCA and other geometric calculations
- Triggers PCA calculation if not already performed

```cpp
geo_vector_t get_pca_axis(int axis) const;
```
- Returns principal component axis vector for specified component
- Parameters:
  - axis: Index of principal component (0, 1, or 2)
    - 0: Primary (longest) axis
    - 1: Secondary axis
    - 2: Tertiary (shortest) axis
- Returns normalized direction vector
- Triggers PCA calculation if not already performed
- Important for determining cluster orientation
- Throws IndexError if axis is invalid

```cpp
double get_pca_value(int axis) const;
```
- Returns eigenvalue for specified principal component
- Parameters:
  - axis: Index of principal component (0, 1, or 2)
    - 0: Largest eigenvalue (most variance)
    - 1: Middle eigenvalue
    - 2: Smallest eigenvalue
- Eigenvalues indicate spread of points along each axis
- Useful for shape analysis:
  - Large ratio between values indicates elongated structure
  - Similar values indicate spherical structure
- Triggers PCA calculation if not already performed
- Throws IndexError if axis is invalid

These functions are particularly important for:
- Track finding and reconstruction
- Cluster shape analysis
- Trajectory determination
- Pattern recognition
- Quality assessment of clustering

Note on PCA Implementation:
- PCA calculation is performed lazily (only when needed)
- Results are cached for efficiency
- Uses Eigen library for eigenvalue decomposition
- Considers charge weighting when available
- Thread-safe for const access
- Invalidated when cluster structure changes

#### Direction Finding
```cpp
std::pair<double, double> hough_transform(
    const geo_point_t& point, 
    const double radius,
    HoughParamSpace param_space = HoughParamSpace::theta_phi,
    std::shared_ptr<const Simple3DPointCloud> s3dpc = nullptr,
    const std::vector<size_t>& global_indices = {}) const;
```
- Performs Hough transform for direction finding
- Supports different parameter spaces (theta-phi or costheta-phi)
- Essential for track direction determination
- Optional external point cloud support

```cpp
geo_vector_t vhough_transform(...);  // Same parameters as hough_transform
```
- Converts Hough transform results to directional vector
- More convenient than raw Hough parameters
- Returns normalized direction vector

### Graph Operations

#### Graph Construction and Management
```cpp
void Create_graph(const bool use_ctpc = true) const;
```
- Creates graph representation of cluster
- Optional CTPC (Continuous Track Point Cloud) usage
- Foundation for path finding operations

```cpp
void Establish_close_connected_graph() const;
```
- Creates edges between points within blobs and overlapping blobs
- Uses distance-based cuts for edge creation
- Important for initial graph structure

```cpp
void Connect_graph(const bool use_ctpc = false) const;
```
- Connects graph components with additional edges
- Optional CTPC usage for validation
- Completes graph connectivity

#### Path Finding
```cpp
void dijkstra_shortest_paths(const size_t pt_idx, const bool use_ctpc = true) const;
```
- Computes shortest paths from given start point
- Uses Dijkstra's algorithm
- Essential for track path finding

```cpp
void cal_shortest_path(const size_t dest_wcp_index) const;
```
- Calculates specific shortest path to destination
- Uses results from dijkstra_shortest_paths
- Updates internal path storage

```cpp
const std::list<size_t>& get_path_wcps() const;
```
- Returns current path point indices
- Available after path finding operations
- Represents ordered sequence of points

#### Path Analysis
```cpp
geo_point_t get_furthest_wcpoint(
    geo_point_t old_wcp, 
    geo_point_t dir,
    const double step = 5*units::cm,
    const int allowed_nstep = 12) const;
```
- Finds furthest point along specified direction
- Uses step size and maximum steps for search
- Important for track extension

```cpp
void adjust_wcpoints_parallel(size_t& start_idx, size_t& end_idx) const;
```
- Adjusts endpoint positions to align with point cloud
- Updates indices in place
- Important for track endpoint refinement

```cpp
bool construct_skeleton(const bool use_ctpc);
```
- Builds skeletal representation of cluster
- Returns false if skeleton already exists
- Important for structural analysis

### Time-Based Operations
```cpp
const time_blob_map_t& time_blob_map() const;
```
- Returns mapping of time slices to blob sets
- Uses lazy initialization
- Essential for temporal analysis

```cpp
std::unordered_map<int, Cluster*> examine_x_boundary(
    const double low_limit = -1*units::cm,
    const double high_limit = 257*units::cm);
```
- Analyzes cluster for boundary crossing
- Can split cluster at boundaries
- Returns map of resulting clusters

### Quality Assessment
```cpp
bool judge_vertex(
    geo_point_t& p_test,
    const double asy_cut = 1./3.,
    const double occupied_cut = 0.85);
```
- Evaluates if point is a vertex
- Uses asymmetry and occupancy criteria
- Updates point position during evaluation

```cpp
bool sanity(Log::logptr_t log = nullptr) const;
```
- Verifies internal consistency of cluster
- Optional logging of issues
- Important for debugging and validation

### Utility Functions
```cpp
size_t hash() const;
```
- Generates hash value representing cluster content
- Based on length and blob hashes
- Useful for comparison and caching

```cpp
static bool cluster_less(const Cluster* a, const Cluster* b);
static void sort_clusters(std::vector<const Cluster*>& clusters);
static void sort_clusters(std::vector<Cluster*>& clusters);
```
- Comparison and sorting functions for clusters
- Based on multiple criteria (length, points, coordinates)
- Both const and non-const versions available

## Important Notes

1. **Performance Considerations**
   - Many operations use lazy evaluation and caching
   - K-d trees are built on-demand
   - Complex operations should be used judiciously in tight loops

2. **Thread Safety**
   - Most const methods are thread-safe
   - Caching operations may not be thread-safe
   - Graph operations should be synchronized if used in multi-threaded context

3. **Memory Management**
   - The class manages various internal caches and data structures
   - Users should be aware of potential memory usage with large point clouds

## Dependencies
- Boost Graph Library
- WireCell utilities and interfaces
- Point Cloud data structures
- K-d tree implementations


## Best Practices

1. **Efficient Point Access**
   - Use bulk point access methods when possible
   - Avoid repeated single point queries in loops
   - Leverage k-d tree functionality for spatial queries

2. **Graph Operations**
   - Create graphs only when needed
   - Reuse path finding results when possible
   - Consider using cached results for repeated operations

3. **Memory Optimization**
   - Clear caches if memory becomes a concern
   - Use appropriate container types for point storage
   - Monitor memory usage with large datasets

## Common Pitfalls

1. **Performance Issues**
   - Avoid repeated construction of k-d trees
   - Don't perform point-by-point operations when bulk operations are available
   - Be careful with large-scale graph operations

2. **Accuracy Considerations**
   - Be aware of floating-point precision in spatial calculations
   - Consider distance metrics carefully in spatial queries
   - Validate results when working with edge cases

3. **Resource Management**
   - Don't assume caches are always valid
   - Clear unnecessary data when working with memory constraints
   - Be mindful of graph construction costs

## Contributing
When extending or modifying the Facade_Cluster class:
- Maintain const correctness
- Update caching mechanisms appropriately
- Document performance implications
- Add appropriate test cases
- Follow existing coding style and conventions

