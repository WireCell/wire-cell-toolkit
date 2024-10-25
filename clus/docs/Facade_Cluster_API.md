# Facade_Cluster Class Documentation

[Previous overview sections remain the same...]

## Complete API Reference

### Constructor & Destructor
```cpp
Cluster();
virtual ~Cluster();
```
- Default constructor initializes an empty cluster
- Virtual destructor ensures proper cleanup of derived classes

### Grouping Access
```cpp
Grouping* grouping();
const Grouping* grouping() const;
```
- Returns pointer to the parent grouping that contains this cluster
- Const version available for read-only access

### Point Cloud Views and Access

#### 3D Point Cloud Operations
```cpp
const sv3d_t& sv3d() const;
```
- Returns the scoped view for the "3d" point cloud (x,y,z coordinates)
- Used internally for point cloud management

```cpp
const kd3d_t& kd3d() const;
```
- Returns the k-d tree for "3d" point cloud
- May trigger k-d tree building if not already constructed

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

```cpp
std::vector<geo_point_t> kd_points(const kd_results_t& res);
std::vector<geo_point_t> kd_points(const kd_results_t& res) const;
```
- Converts k-d tree query results to vector of points
- Both mutable and const versions available

### Blob Management

```cpp
std::vector<Blob*> kd_blobs();
std::vector<const Blob*> kd_blobs() const;
```
- Returns all blobs in k-d tree order
- Order differs from children() order and sort_blobs() order

```cpp
Blob* blob_with_point(size_t point_index);
const Blob* blob_with_point(size_t point_index) const;
```
- Returns blob containing the point at given k-d tree index
- Both mutable and const versions available

```cpp
std::vector<Blob*> blobs_with_points(const kd_results_t& res);
std::vector<const Blob*> blobs_with_points(const kd_results_t& res) const;
```
- Returns blobs containing points from k-d tree query results
- Both mutable and const versions available

### Point Access and Information

```cpp
geo_point_t point3d(size_t point_index) const;
geo_point_t point(size_t point_index) const;
```
- Returns 3D point at given k-d tree index
- point() is alias for point3d() to match Simple3DPointCloud interface

```cpp
const points_type& points() const;
```
- Returns full array of point coordinates
- Recommended for bulk point access

```cpp
int npoints() const;
```
- Returns total number of points in cluster
- Uses cached value when available

```cpp
size_t nbpoints() const;
```
- Returns total number of points according to sum of Blob::nbpoints()
- Equivalent to WCP's get_num_points()

### Spatial Queries and Analysis

```cpp
geo_point_t calc_ave_pos(const geo_point_t& origin, const double dis) const;
```
- Calculates charge-weighted average position of points within distance of origin

```cpp
std::pair<geo_point_t, const Blob*> get_closest_point_blob(const geo_point_t& point) const;
```
- Returns closest point and its containing blob to given point

```cpp
std::pair<size_t, geo_point_t> get_closest_wcpoint(const geo_point_t& p) const;
```
- Returns index and coordinates of closest point to given point

```cpp
size_t get_closest_point_index(const geo_point_t& point) const;
```
- Returns k-d tree index of closest point to given point

```cpp
double get_closest_dis(const geo_point_t& point) const;
```
- Returns distance to closest point from given point

```cpp
int nnearby(const geo_point_t& point, double radius) const;
```
- Counts points within radius of given point
- Uses linear distance measure

```cpp
std::pair<int, int> ndipole(const geo_point_t& point, const geo_point_t& dir, const double dis=-1) const;
```
- Returns count of points in forward/backward direction from point
- Optional distance cutoff (disabled if negative)

### Geometric Analysis

```cpp
std::pair<geo_point_t,geo_point_t> get_two_extreme_points() const;
```
- Finds most distant pair of points along principal axes
- Adjusts positions using local averaging

```cpp
std::pair<geo_point_t, geo_point_t> get_highest_lowest_points(size_t axis = 1) const;
```
- Returns points at extremes of given Cartesian axis
- Default is Y-axis (axis=1)
- Returns points in descending order

```cpp
std::pair<geo_point_t, geo_point_t> get_earliest_latest_points() const;
```
- Returns points at extremes of X-axis
- Returns points in ascending order

```cpp
std::pair<geo_point_t, geo_point_t> get_front_back_points() const;
```
- Returns points at extremes of Z-axis

### Path Finding and Trajectories

```cpp
geo_point_t get_furthest_wcpoint(geo_point_t old_wcp, geo_point_t dir, 
                                const double step = 5*units::cm, 
                                const int allowed_nstep = 12) const;
```
- Finds furthest point along given direction
- Uses step size and maximum number of steps for search

```cpp
void adjust_wcpoints_parallel(size_t& start_idx, size_t& end_idx) const;
```
- Adjusts start and end point positions to align with point cloud

```cpp
bool construct_skeleton(const bool use_ctpc);
```
- Builds skeletal representation of cluster
- Returns false if skeleton already exists

### Graph Operations

```cpp
void Create_graph(const bool use_ctpc = true) const;
```
- Creates graph representation of cluster
- Optional CTPC (Continuous Track Point Cloud) usage

```cpp
void Establish_close_connected_graph() const;
```
- Creates edges between points within blobs and between overlapping blobs
- Uses distance-based cuts

```cpp
void Connect_graph(const bool use_ctpc = false) const;
```
- Connects graph components with additional edges
- Optional CTPC usage

```cpp
void dijkstra_shortest_paths(const size_t pt_idx, const bool use_ctpc = true) const;
```
- Computes shortest paths from given point using Dijkstra's algorithm

```cpp
void cal_shortest_path(const size_t dest_wcp_index) const;
```
- Calculates specific shortest path to destination point

### Time-Based Operations

```cpp
const time_blob_map_t& time_blob_map() const;
```
- Returns mapping of time slices to blob sets
- Lazy initialization

```cpp
const Blob* get_first_blob() const;
```
- Returns blob at earliest time
- Throws ValueError if cluster is empty

```cpp
const Blob* get_last_blob() const;
```
- Returns blob at latest time
- Throws ValueError if cluster is empty

```cpp
size_t get_num_time_slices() const;
```
- Returns number of unique time slices in cluster

### Boundary Analysis

```cpp
std::unordered_map<int, Cluster*> examine_x_boundary(
    const double low_limit = -1*units::cm, 
    const double high_limit = 257*units::cm);
```
- Analyzes cluster crossing boundaries
- Returns map of broken clusters if splitting occurs

### Quality Assessment

```cpp
bool judge_vertex(geo_point_t& p_test, 
                 const double asy_cut = 1./3., 
                 const double occupied_cut = 0.85);
```
- Evaluates if given point is a vertex based on asymmetry and occupancy
- Updates point position during evaluation

### Direction Finding

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

```cpp
geo_vector_t vhough_transform(
    const geo_point_t& point, 
    const double radius,
    HoughParamSpace param_space = HoughParamSpace::theta_phi,
    std::shared_ptr<const Simple3DPointCloud> = nullptr,
    const std::vector<size_t>& global_indices = {}) const;
```
- Converts Hough transform results to directional vector

### Utility Functions

```cpp
size_t hash() const;
```
- Generates hash value representing cluster content

```cpp
bool sanity(Log::logptr_t log = nullptr) const;
```
- Verifies internal consistency of cluster
- Optional logging of issues

```cpp
std::vector<geo_point_t> get_hull() const;
```
- Computes convex hull of cluster points

```cpp
std::tuple<int, int, int, int> get_uvwt_range() const;
std::tuple<int, int, int, int> get_uvwt_min() const;
std::tuple<int, int, int, int> get_uvwt_max() const;
```
- Returns various coordinate ranges for cluster
- Provides minimum and maximum values for U, V, W coordinates and time

### Static Utility Functions

```cpp
static bool cluster_less(const Cluster* a, const Cluster* b);
```
- Comparison function for ordering clusters
- Used in sorting operations

```cpp
static void sort_clusters(std::vector<const Cluster*>& clusters);
static void sort_clusters(std::vector<Cluster*>& clusters);
```
- Sorts clusters in descending order
- Both const and non-const versions available
