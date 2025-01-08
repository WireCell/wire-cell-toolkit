

# Detailed Algorithm Explanations for Facade_Cluster

## 1. Graph Construction Algorithms

### 1.1 Establish_close_connected_graph
This algorithm creates the initial graph structure connecting points within and between blobs.

```python
Algorithm: Establish_close_connected_graph

Input: Collection of blobs and their points
Output: Connected graph representing close points

1. Initialize data structures:
   - Create maps for wire indices (U, V, W planes)
   - Initialize empty graph

2. For each blob:
   - Map points to wire indices
   - Create wire-index to point mappings for each plane
   
3. Create in-blob connections:
   For each point p1 in blob:
       For each point p2 in same blob:
           if are_connected(p1, p2):
               distance = calculate_distance(p1, p2)
               add_edge(p1, p2, distance)

4. Create between-blob connections:
   For each time slice:
       For each blob pair in time slice:
           if blobs_overlap(blob1, blob2):
               connect_overlapping_blobs(blob1, blob2)

5. Apply connection criteria:
   - Maximum wire interval check
   - Distance threshold check
   - Angular separation check
```

Key Features:
- Wire plane-based connectivity check
- Distance-based edge creation
- Time slice consideration
- Angular constraints

### 1.2 Connect_graph
This algorithm enhances connectivity between graph components.

```python
Algorithm: Connect_graph

Input: Initial graph from Establish_close_connected_graph
Output: Fully connected graph with additional edges

1. Find connected components:
   components = find_connected_components(graph)

2. For each component pair:
   2.1 Create point clouds for each component
   2.2 Find closest points between components:
       - Use k-d tree search
       - Consider multiple starting points
       - Apply Hough transform for direction

3. For each potential connection:
   3.1 Validate path:
       - Check point density
       - Verify trajectory smoothness
       - Consider detector geometry
   
   3.2 If path valid:
       - Add new edges
       - Update connectivity

4. Optional CTPC validation:
   If use_ctpc:
       - Validate paths through detector
       - Remove invalid connections
```

## 2. Path Finding Algorithms

### 2.1 Dijkstra Implementation
Customized Dijkstra's algorithm for track finding.

```python
Algorithm: dijkstra_shortest_paths

Input: 
- Starting point index
- Graph structure
- Optional CTPC validation

Output: 
- Distance to all points
- Parent pointers for path reconstruction

1. Initialize:
   - Set all distances to infinity
   - Set source distance to 0
   - Create priority queue Q

2. Custom distance metric:
   distance_metric(v1, v2):
       base_distance = euclidean_distance(v1, v2)
       if use_ctpc:
           quality = evaluate_path_quality(v1, v2)
           return base_distance * quality_factor(quality)
       return base_distance

3. Main loop:
   While Q not empty:
       u = Q.extract_min()
       For each neighbor v of u:
           alt = distance[u] + distance_metric(u, v)
           if alt < distance[v]:
               distance[v] = alt
               parent[v] = u
               Q.decrease_key(v, alt)

4. Path reconstruction:
   reconstruct_path(parent, target):
       path = empty_list
       while target ≠ source:
           path.prepend(target)
           target = parent[target]
       path.prepend(source)
       return path
```

## 3. Spatial Analysis Algorithms

### 3.1 PCA Implementation
Principal Component Analysis for cluster orientation.

```python
Algorithm: Calculate_PCA

Input: Point cloud with optional charge weights
Output: Principal axes and eigenvalues

1. Calculate center:
   center = weighted_mean(points, charges)

2. Build covariance matrix:
   For each point p:
       p_centered = p - center
       For i in [0,1,2]:
           For j in [i,2]:
               cov[i,j] += weight * p_centered[i] * p_centered[j]
               if i != j: cov[j,i] = cov[i,j]

3. Eigendecomposition:
   - Use Eigen library for decomposition
   - Sort eigenvalues in descending order
   - Normalize eigenvectors

4. Cache results:
   - Store center point
   - Store principal axes
   - Store eigenvalues

5. Quality checks:
   - Verify orthogonality
   - Check eigenvalue ratios
   - Validate axis directions
```

### 3.2 Hough Transform for Direction Finding

```python
Algorithm: vhough_transform

Input:
- Reference point
- Search radius
- Parameter space type (theta-phi or costheta-phi)

Output: Direction vector

1. Initialize parameter space:
   - Create 2D histogram (180×360 bins)
   - Define parameter ranges based on space type

2. Collect points:
   points = find_points_in_radius(reference, radius)

3. For each point:
   3.1 Calculate direction vector to reference
   3.2 Convert to chosen parameter space:
       If theta-phi:
           theta = acos(dir.z)
           phi = atan2(dir.y, dir.x)
       If costheta-phi:
           costheta = dir.z
           phi = atan2(dir.y, dir.x)
   
   3.3 Weight contribution:
       weight = calculate_weight(point)
       histogram[theta_bin][phi_bin] += weight

4. Find maximum:
   - Locate histogram maximum
   - Convert parameters back to vector
   - Normalize result

5. Optional refinement:
   - Local maximum fitting
   - Multiple peak detection
   - Angular uncertainty estimation
```

## 4. Connectivity Analysis

### 4.1 Blob Connection Algorithm

```python
Algorithm: is_connected

Input: 
- Two clusters
- Offset parameter for tolerance

Output: Vector of connecting blobs

1. Time slice mapping:
   - Create time-ordered blob maps for both clusters
   - Find overlapping time ranges

2. For each time overlap:
   2.1 Check blob pairs:
       For blob1 in cluster1[time]:
           For blob2 in cluster2[time]:
               if blobs_overlap(blob1, blob2, offset):
                   add_to_connections(blob1)

3. Overlap check:
   blob_overlap(blob1, blob2, offset):
       - Check wire index ranges
       - Apply offset tolerance
       - Verify spatial proximity
       - Consider charge distributions

4. Return connecting blobs:
   - Sort by connection strength
   - Filter weak connections
   - Validate consistency
```

## 5. Point Cloud Operations

### 5.1 Closest Point Search Algorithm

```python
Algorithm: get_closest_points

Input: Two clusters
Output: Closest point indices and distance

1. Initialize search:
   min_distance = infinity
   best_points = null

2. Multi-start search:
   2.1 Try endpoints:
       - Start-Start
       - Start-End
       - End-Start
       - End-End

   2.2 For each start pair:
       current = start_pair
       while not converged:
           p1 = find_closest(cluster1, p2)
           p2 = find_closest(cluster2, p1)
           if no_improvement:
               break
           update_best_if_better(p1, p2)

3. Local refinement:
   - Use k-d tree for efficient nearby point search
   - Consider charge weights
   - Apply geometric constraints

4. Validation:
   - Check distance consistency
   - Verify path feasibility
   - Consider detector constraints
```

Would you like me to:
1. Add more specific implementation details for any algorithm?
2. Explain the mathematical foundations of certain algorithms?
3. Provide complexity analysis for the algorithms?
4. Show more example use cases?