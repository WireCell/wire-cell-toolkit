#include "WireCellUtil/Exceptions.h"
#include "SteinerGrapher.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_graph()
{
    return fake_steiner_graph(); // see SteinerGrapher_face.cxx
}






Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_tree(/*what type for point_cloud_steiner?*/)
{
    Steiner::Grapher::graph_type graph_steiner;
    raise<LogicError>("not implemented yet");
    return graph_steiner;
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