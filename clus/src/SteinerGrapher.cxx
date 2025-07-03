#include "WireCellUtil/Exceptions.h"
#include "SteinerGrapher.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellClus/Graphs.h"
#include <algorithm>
#include <set>
#include <map>
#include <boost/graph/copy.hpp>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_graph()
{
    //  auto cell_points_map = form_cell_points_map();
    // std::cout << "Xin: " << " " << cell_points_map.size() << std::endl;

    return fake_steiner_graph(); // see SteinerGrapher_face.cxx
}


Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_tree(
    const Facade::Cluster* reference_cluster,
    const std::vector<size_t>& path_point_indices,
    bool disable_dead_mix_cell,
    const std::string& steiner_pc_name)
{
    log->debug("create_steiner_tree: starting with reference_cluster={}, path_size={}", 
               (reference_cluster ? "provided" : "null"), path_point_indices.size());

    // Phase 1: Find initial steiner terminals
    vertex_set steiner_terminals = find_steiner_terminals(disable_dead_mix_cell);
    log->debug("create_steiner_tree: found {} initial steiner terminals", steiner_terminals.size());

    if (steiner_terminals.empty()) {
        log->warn("create_steiner_tree: no steiner terminals found, returning empty graph");
        return graph_type(0);
    }

    // Phase 2: Apply reference cluster spatial filtering
    if (reference_cluster) {
        vertex_set original_size = steiner_terminals;
        steiner_terminals = filter_by_reference_cluster(steiner_terminals, reference_cluster);
        log->debug("create_steiner_tree: reference cluster filtering: {} -> {} terminals", 
                   original_size.size(), steiner_terminals.size());
    }

    // Phase 3: Apply path-based filtering if path is provided
    if (!path_point_indices.empty()) {
        vertex_set pre_path_size = steiner_terminals;
        steiner_terminals = filter_by_path_constraints(steiner_terminals, path_point_indices);
        log->debug("create_steiner_tree: path filtering: {} -> {} terminals", 
                   pre_path_size.size(), steiner_terminals.size());
    }

    // Phase 4: Add extreme points
    vertex_set extreme_points = get_extreme_points_for_reference(reference_cluster);
    steiner_terminals.insert(extreme_points.begin(), extreme_points.end());
    log->debug("create_steiner_tree: added {} extreme points, total terminals: {}", 
               extreme_points.size(), steiner_terminals.size());

    if (steiner_terminals.empty()) {
        log->warn("create_steiner_tree: no terminals remain after filtering, returning empty graph");
        return graph_type(0);
    }

    // Phase 5: Create subset point cloud for steiner points
    auto steiner_pc = create_steiner_subset_pc(steiner_terminals);
    put_point_cloud(std::move(steiner_pc), steiner_pc_name);
    log->debug("create_steiner_tree: created steiner subset point cloud '{}'", steiner_pc_name);

    // Phase 6: Build Steiner tree on the subset
    const auto& base_graph = get_graph("basic");
    auto& graph_algo = m_cluster.graph_algorithms("basic");

    // Create filtered graph with only steiner vertices
    auto filtered_graph = graph_algo.reduce(steiner_terminals);
    
    // Convert filtered graph to proper graph type using boost::copy_graph
    graph_type steiner_graph_input;
    boost::copy_graph(filtered_graph, steiner_graph_input);
    
    // Convert terminal set to vector for algorithms
    std::vector<vertex_type> terminal_vector(steiner_terminals.begin(), steiner_terminals.end());

    // Apply Voronoi tessellation and create Steiner graph
    auto vor = Graphs::Weighted::voronoi(steiner_graph_input, terminal_vector);
    auto steiner_result = Graphs::Weighted::steiner_graph(steiner_graph_input, vor);

    log->debug("create_steiner_tree: created steiner graph with {} vertices, {} edges", 
               boost::num_vertices(steiner_result), boost::num_edges(steiner_result));

    return steiner_result;
}

// ========================================
// Helper Method Implementations
// ========================================

Steiner::Grapher::vertex_set Steiner::Grapher::filter_by_reference_cluster(
    const vertex_set& terminals,
    const Facade::Cluster* reference_cluster) const
{
    if (!reference_cluster) {
        return terminals;
    }

    vertex_set filtered_terminals;
    
    // Get reference cluster's time-blob mapping
    const auto& ref_time_blob_map = reference_cluster->time_blob_map();
    
    if (ref_time_blob_map.empty()) {
        log->debug("filter_by_reference_cluster: reference cluster has empty time_blob_map");
        return terminals;
    }

    // Filter terminals based on spatial relationship with reference cluster
    for (auto terminal_idx : terminals) {
        if (is_point_spatially_related_to_reference(terminal_idx, ref_time_blob_map)) {
            filtered_terminals.insert(terminal_idx);
        }
    }

    return filtered_terminals;
}

Steiner::Grapher::vertex_set Steiner::Grapher::filter_by_path_constraints(
    const vertex_set& terminals,
    const std::vector<size_t>& path_point_indices) const
{
    if (path_point_indices.empty()) {
        return terminals;
    }

    // Create path point cloud for distance calculations
    auto path_pc = create_path_point_cloud(path_point_indices);
    
    vertex_set filtered_terminals;
    
    // Distance thresholds from prototype
    const double distance_3d_threshold = 6.0 * units::cm;
    const double distance_2d_threshold = 1.8 * units::cm;

    for (auto terminal_idx : terminals) {
        Point point = m_cluster.point3d(terminal_idx);
        
        // Calculate distances similar to prototype's logic
        double dis_3d = calculate_closest_distance_3d(point, path_pc);
        auto dis_2d = calculate_closest_distances_2d(point, path_pc);
        
        // Apply prototype's filtering logic:
        // Remove if close in 2D projections but far in 3D
        bool close_in_2d = (dis_2d[0] < distance_2d_threshold && dis_2d[1] < distance_2d_threshold) ||
                          (dis_2d[0] < distance_2d_threshold && dis_2d[2] < distance_2d_threshold) ||
                          (dis_2d[1] < distance_2d_threshold && dis_2d[2] < distance_2d_threshold);
        
        bool should_remove = close_in_2d && (dis_3d > distance_3d_threshold);
        
        if (!should_remove) {
            filtered_terminals.insert(terminal_idx);
        }
    }

    return filtered_terminals;
}

Steiner::Grapher::vertex_set Steiner::Grapher::get_extreme_points_for_reference(
    const Facade::Cluster* reference_cluster) const
{
    vertex_set extreme_points;
    
    // Use cluster's existing extreme point calculation
    // If reference cluster is provided, we should get extreme points that are
    // spatially related to it. For now, use all extreme points and filter later.
    
    try {
        // Get extreme points using cluster's method
        // This maps to prototype's get_extreme_wcps() functionality
        auto extreme_point_groups = m_cluster.get_extreme_wcps(reference_cluster);
        
        // Convert extreme points to vertex indices
        for (const auto& point_group : extreme_point_groups) {
            for (const auto& point : point_group) {
                // Find the vertex index for this point
                // This requires mapping 3D point back to vertex index
                auto closest_idx = find_closest_vertex_to_point(point);
                if (closest_idx != SIZE_MAX) {
                    extreme_points.insert(closest_idx);
                }
            }
        }
    } catch (const std::exception& e) {
        log->warn("get_extreme_points_for_reference: failed to get extreme points: {}", e.what());
    }

    return extreme_points;
}

PointCloud::Dataset Steiner::Grapher::create_path_point_cloud(
    const std::vector<size_t>& path_indices) const
{
    if (path_indices.empty()) {
        return PointCloud::Dataset();
    }

    const double step_dis = 0.6 * units::cm;
    std::vector<Point> interpolated_points;
    std::vector<double> x_coords, y_coords, z_coords;
    
    Point prev_point = m_cluster.point3d(path_indices[0]);
    interpolated_points.push_back(prev_point);

    for (size_t i = 1; i < path_indices.size(); ++i) {
        Point current_point = m_cluster.point3d(path_indices[i]);
        double distance = (current_point - prev_point).magnitude();
        
        if (distance <= step_dis) {
            interpolated_points.push_back(current_point);
        } else {
            // Interpolate points at step_dis intervals (matching prototype logic)
            int num_steps = static_cast<int>(distance / step_dis);
            for (int step = 1; step <= num_steps; ++step) {
                double fraction = static_cast<double>(step) / num_steps;
                Point interpolated = prev_point + (current_point - prev_point) * fraction;
                interpolated_points.push_back(interpolated);
            }
            interpolated_points.push_back(current_point);
        }
        prev_point = current_point;
    }

    return points_to_dataset(interpolated_points);
}

PointCloud::Dataset Steiner::Grapher::create_steiner_subset_pc(
    const vertex_set& steiner_indices) const
{
    // Get the original point cloud from the cluster's scoped view
    const auto& sv = m_cluster.sv3d();
    const auto original_pc = sv.flat_pc("3d");
    
    // Convert set to vector for subset operation
    std::vector<size_t> indices_vector(steiner_indices.begin(), steiner_indices.end());
    
    // Create subset point cloud containing only steiner points
    return original_pc.subset(indices_vector);
}

bool Steiner::Grapher::is_point_spatially_related_to_reference(
    size_t point_idx,
    const Facade::Cluster::time_blob_map_t& ref_time_blob_map) const
{
    // Delegate to the cluster's existing method which implements the proper logic
    // for checking spatial relationships with the complex time_blob_map structure
    return m_cluster.is_point_spatially_related_to_time_blobs(point_idx, ref_time_blob_map);
}

double Steiner::Grapher::calculate_closest_distance_3d(
    const Point& point,
    const PointCloud::Dataset& path_pc) const
{
    if (path_pc.size_major() == 0) {
        return std::numeric_limits<double>::max();
    }

    const auto& x_coords = path_pc.get("x")->elements<double>();
    const auto& y_coords = path_pc.get("y")->elements<double>();
    const auto& z_coords = path_pc.get("z")->elements<double>();

    double min_distance = std::numeric_limits<double>::max();
    
    for (size_t i = 0; i < x_coords.size(); ++i) {
        Point path_point(x_coords[i], y_coords[i], z_coords[i]);
        double distance = (point - path_point).magnitude();
        min_distance = std::min(min_distance, distance);
    }

    return min_distance;
}

std::array<double, 3> Steiner::Grapher::calculate_closest_distances_2d(
    const Point& point,
    const PointCloud::Dataset& path_pc) const
{
    std::array<double, 3> min_distances = {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()
    };

    if (path_pc.size_major() == 0) {
        return min_distances;
    }

    const auto& x_coords = path_pc.get("x")->elements<double>();
    const auto& y_coords = path_pc.get("y")->elements<double>();
    const auto& z_coords = path_pc.get("z")->elements<double>();

    // Calculate wire plane angles (approximation)
    const double angle_u = 1.047;  // ~60 degrees
    const double angle_v = -1.047; // ~-60 degrees  
    const double angle_w = 0.0;    // 0 degrees

    for (size_t i = 0; i < x_coords.size(); ++i) {
        Point path_point(x_coords[i], y_coords[i], z_coords[i]);
        
        // U projection distance
        double u_proj_point = std::cos(angle_u) * point.z() - std::sin(angle_u) * point.y();
        double u_proj_path = std::cos(angle_u) * path_point.z() - std::sin(angle_u) * path_point.y();
        double u_distance = std::abs(u_proj_point - u_proj_path);
        min_distances[0] = std::min(min_distances[0], u_distance);
        
        // V projection distance
        double v_proj_point = std::cos(angle_v) * point.z() - std::sin(angle_v) * point.y();
        double v_proj_path = std::cos(angle_v) * path_point.z() - std::sin(angle_v) * path_point.y();
        double v_distance = std::abs(v_proj_point - v_proj_path);
        min_distances[1] = std::min(min_distances[1], v_distance);
        
        // W projection distance
        double w_proj_point = std::cos(angle_w) * point.z() - std::sin(angle_w) * point.y();
        double w_proj_path = std::cos(angle_w) * path_point.z() - std::sin(angle_w) * path_point.y();
        double w_distance = std::abs(w_proj_point - w_proj_path);
        min_distances[2] = std::min(min_distances[2], w_distance);
    }

    return min_distances;
}

PointCloud::Dataset Steiner::Grapher::points_to_dataset(const std::vector<Point>& points) const
{
    std::vector<double> x_coords, y_coords, z_coords;
    x_coords.reserve(points.size());
    y_coords.reserve(points.size());
    z_coords.reserve(points.size());

    for (const auto& point : points) {
        x_coords.push_back(point.x());
        y_coords.push_back(point.y());
        z_coords.push_back(point.z());
    }

    PointCloud::Dataset dataset;
    dataset.add("x", PointCloud::Array(x_coords));
    dataset.add("y", PointCloud::Array(y_coords));
    dataset.add("z", PointCloud::Array(z_coords));

    return dataset;
}

// ========================================
// Additional Helper Methods (to be implemented)
// ========================================

// These methods need cluster-specific implementation:

size_t Steiner::Grapher::find_closest_vertex_to_point(const Point& point) const
{
    // Find the vertex index closest to the given 3D point
    double min_distance = std::numeric_limits<double>::max();
    size_t closest_idx = SIZE_MAX;
    
    for (size_t i = 0; i < m_cluster.npoints(); ++i) {
        Point vertex_point = m_cluster.point3d(i);
        double distance = (point - vertex_point).magnitude();
        if (distance < min_distance) {
            min_distance = distance;
            closest_idx = i;
        }
    }
    
    return closest_idx;
}

bool Steiner::Grapher::is_point_near_blob(const Point& point, const Facade::Blob* blob) const
{
    if (!blob) return false;
    
    // Check if point is within blob's spatial extent
    // This is an approximation - you may need to implement more sophisticated logic
    Point blob_center = blob->center_pos();
    double distance = (point - blob_center).magnitude();
    double blob_radius = 1*units::cm ;//blob->extent();  // or some characteristic size
    
    // Use tolerance similar to prototype's wire index ranges
    const double spatial_tolerance = 2.0 * units::cm;
    
    return distance <= (blob_radius + spatial_tolerance);
}



// Overloaded version for specific blobs (equivalent to prototype's mcells parameter)
Steiner::Grapher::vertex_set Steiner::Grapher::find_peak_point_indices(
    const std::vector<const Facade::Blob*>& target_blobs, 
    bool disable_dead_mix_cell)
{
    vertex_set peak_points;
    
    if (target_blobs.empty()) {
        return peak_points;
    }
    
    // Get the blob-to-points mapping
    auto cell_points_map = form_cell_points_map();
    
    // Collect indices only from target blobs
    vertex_set all_indices;
    for (const auto* blob : target_blobs) {
        auto it = cell_points_map.find(blob);
        if (it != cell_points_map.end()) {
            all_indices.insert(it->second.begin(), it->second.end());
        }
    }
    
    if (all_indices.empty()) {
        return peak_points;
    }
    
    // Rest of the implementation is similar to the above version
    // but only operates on the filtered point indices
    
    // Calculate charges and find candidates
    std::map<size_t, double> map_index_charge;
    std::set<std::pair<double, size_t>, std::greater<std::pair<double, size_t>>> candidates_set;
    
    const double charge_threshold = 4000.0;
    
    for (size_t point_idx : all_indices) {
        auto [charge_quality, charge] = m_cluster.calc_charge_wcp(point_idx, 0.0, disable_dead_mix_cell);
        
        map_index_charge[point_idx] = charge;
        
        if (charge > charge_threshold && charge_quality) {
            candidates_set.insert(std::make_pair(charge, point_idx));
        }
    }
    
    if (candidates_set.empty()) {
        return peak_points;
    }
    
    // Convert candidates to vector for graph operations
    std::vector<size_t> vec_peak_indices;
    vec_peak_indices.reserve(candidates_set.size());
    for (const auto& [charge, point_idx] : candidates_set) {
        vec_peak_indices.push_back(point_idx);
    }
    
    // Find connected components among peak candidates
    // Note: This is a simplified version - the prototype uses boost::connected_components
    // but we need to adapt this to work with the toolkit's graph structure
    
    std::vector<int> component(vec_peak_indices.size());
    
    // For now, implement a simple connected components algorithm
    // Each point starts as its own component
    std::iota(component.begin(), component.end(), 0);
    
    // Check adjacency in the graph and merge components
    // This is a simplified version - you may need to adapt based on the actual graph structure
    for (size_t i = 0; i < vec_peak_indices.size(); ++i) {
        for (size_t j = i + 1; j < vec_peak_indices.size(); ++j) {
            size_t idx1 = vec_peak_indices[i];
            size_t idx2 = vec_peak_indices[j];
            
            // Check if points are connected in graph
            // Note: This needs to be adapted to your specific graph representation
            // For now, use spatial proximity as a simple heuristic
            auto p1 = m_cluster.point3d(idx1);
            auto p2 = m_cluster.point3d(idx2);
            double distance = sqrt(pow(p1.x() - p2.x(), 2) + 
                                 pow(p1.y() - p2.y(), 2) + 
                                 pow(p1.z() - p2.z(), 2));
            
            const double connection_threshold = 3.0; // Adjust as needed
            if (distance < connection_threshold) {
                // Merge components
                int comp1 = component[i];
                int comp2 = component[j];
                if (comp1 != comp2) {
                    // Replace all comp2 with comp1
                    for (size_t k = 0; k < component.size(); ++k) {
                        if (component[k] == comp2) {
                            component[k] = comp1;
                        }
                    }
                }
            }
        }
    }
    
    // Find unique components
    std::set<int> unique_components(component.begin(), component.end());
    int num_components = unique_components.size();
    
    if (num_components == 0) {
        return peak_points;
    }
    
    // For each component, find the point closest to center of mass
    std::vector<double> min_dis(num_components, 1e9);
    std::vector<size_t> min_index(num_components, 0);
    std::vector<int> ncounts(num_components, 0);
    std::vector<WireCell::Point> centers(num_components, WireCell::Point(0, 0, 0));
    
    // Map component IDs to array indices
    std::map<int, int> comp_to_idx;
    int idx = 0;
    for (int comp : unique_components) {
        comp_to_idx[comp] = idx++;
    }
    
    // Calculate center of mass for each component
    for (size_t i = 0; i < component.size(); ++i) {
        int comp_idx = comp_to_idx[component[i]];
        ncounts[comp_idx]++;
        
        auto point = m_cluster.point3d(vec_peak_indices[i]);
        centers[comp_idx] = centers[comp_idx] + WireCell::Vector(point.x(), point.y(), point.z());
    }
    
    // Average the positions
    for (int i = 0; i < num_components; ++i) {
        if (ncounts[i] > 0) {
            centers[i] = centers[i] * (1.0 / ncounts[i]);
        }
    }
    
    // Find closest point to center for each component
    for (size_t i = 0; i < component.size(); ++i) {
        int comp_idx = comp_to_idx[component[i]];
        auto point = m_cluster.point3d(vec_peak_indices[i]);
        
        double dis = pow(centers[comp_idx].x() - point.x(), 2) +
                    pow(centers[comp_idx].y() - point.y(), 2) +
                    pow(centers[comp_idx].z() - point.z(), 2);
        
        if (dis < min_dis[comp_idx]) {
            min_dis[comp_idx] = dis;
            min_index[comp_idx] = vec_peak_indices[i];
        }
    }
    
    // Add the representative points to the result
    for (int i = 0; i < num_components; ++i) {
        peak_points.insert(min_index[i]);
    }
    
    return peak_points;
}


Steiner::Grapher::vertex_set Steiner::Grapher::find_steiner_terminals(bool disable_dead_mix_cell)
{
    vertex_set steiner_terminals;
    
    // Get the blob-to-points mapping 
    auto cell_points_map = form_cell_points_map();
    
    if (cell_points_map.empty()) {
        return steiner_terminals;
    }
    
    // Process each blob individually (following prototype pattern)
    for (const auto& [blob, point_indices] : cell_points_map) {
        // Create a single-blob vector for processing
        std::vector<const Blob*> single_blob = {blob};
        
        // Find peak points for this specific blob
        auto blob_peaks = find_peak_point_indices(single_blob, disable_dead_mix_cell);
        
        // Add to overall steiner terminals set
        steiner_terminals.insert(blob_peaks.begin(), blob_peaks.end());
    }
    
    return steiner_terminals;
}



Steiner::Grapher::blob_vertex_map Steiner::Grapher::form_cell_points_map()
{
    blob_vertex_map cell_points;
    
    // Get the 3D scoped view with x,y,z coordinates
    const auto& sv = m_cluster.sv3d();   // This is default scoped view with 3D coordinates ... 
    const auto& nodes = sv.nodes();   // These are the blob nodes
    const auto& skd = sv.kd();        // K-d tree with point data
    
    // Check if we have valid data
    if (nodes.empty() || skd.npoints() == 0) {
        return cell_points;  // Return empty map if no data
    }
    
    // The major indices tell you which blob each point belongs to
    const auto& majs = skd.major_indices();  // blob index for each point
    
    // Iterate through all points and assign them to their respective blobs
    for (size_t point_idx = 0; point_idx < skd.npoints(); ++point_idx) {
        size_t blob_idx = majs[point_idx];
        
        // Bounds check
        if (blob_idx >= nodes.size()) {
            continue;  // Skip invalid blob indices
        }
        
        // Get the blob facade from the node
        const auto* blob = nodes[blob_idx]->value.facade<Blob>();
        if (!blob) {
            continue;  // Skip if facade creation failed
        }
        
        // Initialize the set for this blob if it doesn't exist
        if (cell_points.find(blob) == cell_points.end()) {
            cell_points[blob] = vertex_set();
        }
        
        // Add this point index to the blob's set
        cell_points[blob].insert(point_idx);
    }
    
    return cell_points;
}

void Steiner::Grapher::establish_same_blob_steiner_edges(const std::string& graph_name, 
                                                        bool disable_dead_mix_cell, int flag)
{
    if (!m_cluster.has_graph(graph_name)) {
        log->error("Graph '{}' does not exist in cluster", graph_name);
        return;
    }

    auto& graph = m_cluster.get_graph(graph_name);
    edge_set added_edges;

    // Step 1: Find Steiner terminals using the existing implementation
    vertex_set steiner_terminals = find_steiner_terminals(disable_dead_mix_cell);
    
    log->debug("Found {} Steiner terminals for same-blob edge establishment", steiner_terminals.size());

    // Step 2: Get the blob-to-points mapping (equivalent to map_mcell_all_indices in prototype)
    auto cell_points_map = form_cell_points_map();
    
    if (cell_points_map.empty()) {
        log->warn("No blob-to-points mapping available for Steiner edge establishment");
        return;
    }

    log->debug("Processing {} blobs for same-blob edges", cell_points_map.size());

    // Step 3: For each blob, add edges between all pairs of points (following prototype logic)
    for (const auto& [blob, point_indices] : cell_points_map) {
        if (point_indices.size() < 2) {
            continue; // Need at least 2 points to make edges
        }

        // // Skip dead cells if requested
        // if (disable_dead_mix_cell && blob->is_dead()) {
        //     continue;
        // }

        log->debug("Processing blob with {} points", point_indices.size());

        // Convert set to vector for easier iteration
        std::vector<vertex_type> points_vec(point_indices.begin(), point_indices.end());

        // Add edges between all pairs of points in the same blob (following prototype)
        for (size_t i = 0; i < points_vec.size(); ++i) {
            vertex_type index1 = points_vec[i];
            bool flag_index1 = (steiner_terminals.find(index1) != steiner_terminals.end());
            
            for (size_t j = i + 1; j < points_vec.size(); ++j) {
                vertex_type index2 = points_vec[j];
                bool flag_index2 = (steiner_terminals.find(index2) != steiner_terminals.end());

                // Calculate base distance between points
                double distance = calculate_distance(index1, index2);
                
                // Determine edge weight based on terminal status (following prototype logic)
                double edge_weight = 0.0;
                bool add_edge = false;
                
                if (flag_index1 && flag_index2) {
                    // Both are steiner terminals: weight = distance * 0.8
                    edge_weight = distance * 0.8;
                    add_edge = true;
                } else if (flag_index1 || flag_index2) {
                    // One is steiner terminal: weight = distance * 0.9  
                    edge_weight = distance * 0.9;
                    add_edge = true;
                }
                // If neither is a steiner terminal, don't add edge (add_edge stays false)

                if (add_edge) {
                    // Add edge with calculated weight
                    auto [edge, success] = boost::add_edge(index1, index2, edge_weight, graph);
                    if (success) {
                        added_edges.insert(edge);
                        log->debug("Added same-blob edge: {} -- {} (distance: {:.3f} cm, weight: {:.3f}, terminals: {}/{})", 
                                  index1, index2, distance / units::cm, edge_weight / units::cm,
                                  flag_index1 ? "T" : "N", flag_index2 ? "T" : "N");
                    }
                }
            }
        }
    }

    // Store the added edges for later removal
    store_added_edges(graph_name, added_edges);

    // Invalidate any cached GraphAlgorithms that use this graph
    invalidate_graph_algorithms_cache(graph_name);

    log->info("Added {} same-blob edges to graph '{}' from {} total points ({} steiner terminals)", 
             added_edges.size(), graph_name, 
             std::accumulate(cell_points_map.begin(), cell_points_map.end(), 0,
                           [](int sum, const auto& pair) { return sum + pair.second.size(); }),
             steiner_terminals.size());
}


void Steiner::Grapher::remove_same_blob_steiner_edges(const std::string& graph_name)
{
    if (!m_cluster.has_graph(graph_name)) {
        log->warn("Graph '{}' does not exist, cannot remove edges", graph_name);
        return;
    }

    auto it = m_added_edges_by_graph.find(graph_name);
    if (it == m_added_edges_by_graph.end() || it->second.empty()) {
        log->debug("No edges to remove for graph '{}'", graph_name);
        return;
    }

    auto& graph = m_cluster.get_graph(graph_name);
    const auto& edges_to_remove = it->second;
    size_t removed_count = 0;

    // Remove the edges
    for (const auto& edge : edges_to_remove) {
        boost::remove_edge(edge, graph);
        ++removed_count;
    }

    // Clear the tracking for this graph
    it->second.clear();

    // Invalidate any cached GraphAlgorithms that use this graph
    invalidate_graph_algorithms_cache(graph_name);

    log->info("Removed {} same-blob Steiner edges from graph '{}'", removed_count, graph_name);
}

void Steiner::Grapher::invalidate_graph_algorithms_cache(const std::string& graph_name)
{
    // Use the new public method we'll add to Cluster
    m_cluster.clear_graph_algorithms_cache(graph_name);
}

void Steiner::Grapher::store_added_edges(const std::string& graph_name, const edge_set& edges)
{
    // Add to existing set if graph already has tracked edges
    auto& tracked_edges = m_added_edges_by_graph[graph_name];
    tracked_edges.insert(edges.begin(), edges.end());
}

bool Steiner::Grapher::same_blob(vertex_type v1, vertex_type v2) const
{
    const auto* blob1 = get_blob_for_vertex(v1);
    const auto* blob2 = get_blob_for_vertex(v2);
    return (blob1 && blob2 && blob1 == blob2);
}

double Steiner::Grapher::calculate_distance(vertex_type v1, vertex_type v2) const
{
     // Use cluster's point3d method to get 3D coordinates  
    // This is the standard way to access point coordinates in the toolkit
    auto point1 = m_cluster.point3d(v1);
    auto point2 = m_cluster.point3d(v2);
    
    // Calculate Euclidean distance
    double dx = point2.x() - point1.x();
    double dy = point2.y() - point1.y();
    double dz = point2.z() - point1.z();
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

const Facade::Blob* Steiner::Grapher::get_blob_for_vertex(vertex_type vertex) const
{
    const auto& sv = m_cluster.sv3d();
    const auto& nodes = sv.nodes();
    const auto& skd = sv.kd();
    
    if (vertex >= skd.npoints()) {
        return nullptr;
    }
    
    const auto& majs = skd.major_indices();
    size_t blob_idx = majs[vertex];
    
    if (blob_idx >= nodes.size()) {
        return nullptr;
    }
    
    return nodes[blob_idx]->value.facade<Blob>();
}