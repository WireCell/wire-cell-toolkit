#include "WireCellUtil/Exceptions.h"
#include "SteinerGrapher.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_graph()
{
    return fake_steiner_graph(); // see SteinerGrapher_face.cxx
}


void Steiner::Grapher::calc_sampling_points(/*, ...*/)
{
    raise<LogicError>("not implemented yet");

    // Xin, this either needs access to IBlobs in order to sample from scratch
    // or it needs to take points from previously sampled clusters.
}





void Steiner::Grapher::establish_same_blob_steiner_edges(graph_type& graph, bool disable_dead_mix_cell, int flag)
{
    raise<LogicError>("not implemented yet");

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