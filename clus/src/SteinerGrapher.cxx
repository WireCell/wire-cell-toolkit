#include "WireCellUtil/Exceptions.h"
#include "SteinerGrapher.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellClus/Graphs.h"
#include "WireCellClus/DynamicPointCloud.h"
#include <algorithm>
#include <set>
#include <unordered_set>
#include <map>
#include <chrono>
#include <boost/graph/copy.hpp>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;


void Steiner::Grapher::create_steiner_tree(
    const Facade::Cluster* reference_cluster, // may not be the same as m_cluster
    const std::vector<size_t>& path_point_indices, // of m_cluster
    const std::string& graph_name,
    const std::string& steiner_graph_name,
    bool disable_dead_mix_cell,
    const std::string& steiner_pc_name)
{
    SPDLOG_LOGGER_DEBUG(log, "create_steiner_tree: starting with reference_cluster={}, path_size={}", 
               (reference_cluster ? "provided" : "null"), path_point_indices.size());

    // Phase 1: Find initial steiner terminals
    vertex_set steiner_terminals = find_steiner_terminals(graph_name, disable_dead_mix_cell);
    SPDLOG_LOGGER_DEBUG(log, "create_steiner_tree: found {} initial steiner terminals", steiner_terminals.size());

    // std::cout << "Test1: " << steiner_terminals.size() << std::endl;

    if (steiner_terminals.size() < 2) {
        SPDLOG_LOGGER_WARN(log, "create_steiner_tree: only {} steiner terminal(s) found (need >=2), returning empty graph",
                           steiner_terminals.size());
        return;
    }

    // Phase 2: Apply reference cluster spatial filtering
    if (reference_cluster) {
        vertex_set original_size = steiner_terminals;
        steiner_terminals = filter_by_reference_cluster(steiner_terminals, reference_cluster);
        SPDLOG_LOGGER_DEBUG(log, "create_steiner_tree: reference cluster filtering: {} -> {} terminals", 
                   original_size.size(), steiner_terminals.size());
    }

    // std::cout << "Test2: " << steiner_terminals.size() << std::endl;


    // Phase 3: Apply path-based filtering if path is provided
    if (!path_point_indices.empty()) {
        vertex_set pre_path_size = steiner_terminals;
        steiner_terminals = filter_by_path_constraints(steiner_terminals, path_point_indices);
        SPDLOG_LOGGER_DEBUG(log, "create_steiner_tree: path filtering: {} -> {} terminals", 
                   pre_path_size.size(), steiner_terminals.size());
    }

    // std::cout << "Test3: " << steiner_terminals.size() << std::endl;


    // Phase 4: Add extreme points
    vertex_set extreme_points = get_extreme_points_for_reference(reference_cluster);
    steiner_terminals.insert(extreme_points.begin(), extreme_points.end());
    SPDLOG_LOGGER_DEBUG(log, "create_steiner_tree: added {} extreme points, total terminals: {}", 
               extreme_points.size(), steiner_terminals.size());

    // std::cout << "Test4: " << steiner_terminals.size() << std::endl;


    if (steiner_terminals.size() < 2) {
        SPDLOG_LOGGER_WARN(log, "create_steiner_tree: only {} terminal(s) remain after filtering (need >=2), returning empty graph",
                           steiner_terminals.size());
        return;
    }

    const auto& base_graph = get_graph(graph_name);

   const auto& original_pc = get_point_cloud("default"); // default pc ...
    
    // Configure charge weighting to match prototype values
    Graphs::Weighted::ChargeWeightingConfig charge_config;
    charge_config.Q0 = 10000.0;           // From prototype
    charge_config.factor1 = 0.8;          // From prototype  
    charge_config.factor2 = 0.4;          // From prototype
    charge_config.enable_weighting = true; // Enable charge weighting
    
    // Use the enhanced approach with cluster reference for charge calculation
    auto steiner_result = Graphs::Weighted::create_enhanced_steiner_graph(
        base_graph, steiner_terminals, original_pc, m_cluster, charge_config);

    // std::cout << "Test5: " <<  " Graph vertices: " << boost::num_vertices(steiner_result.graph) << ", edges: " << boost::num_edges(steiner_result.graph) << std::endl;

    // just run this once the steiner graph is created
    Graphs::Weighted::establish_same_blob_steiner_edges_steiner_graph(steiner_result, m_cluster);

    // std::cout << "Test5: " <<  " Graph vertices: " << boost::num_vertices(steiner_result.graph) << ", edges: " << boost::num_edges(steiner_result.graph) << std::endl;

    // Phase 6: Store results for later access
    m_flag_steiner_terminal = steiner_result.flag_steiner_terminal;
    m_old_to_new_index = steiner_result.old_to_new_index;
    m_new_to_old_index = steiner_result.new_to_old_index;
    
    // Store the subset point cloud
    if (!steiner_pc_name.empty()) {
        put_point_cloud(std::move(steiner_result.point_cloud), steiner_pc_name);
        SPDLOG_LOGGER_DEBUG(log, "create_steiner_tree: created steiner subset point cloud '{}'", steiner_pc_name);
    }

    // can I do this to store the graph ???
    m_cluster.give_graph(steiner_graph_name, std::move(steiner_result.graph));

    
    SPDLOG_LOGGER_DEBUG(log, "create_steiner_tree: created reduced steiner graph with {} vertices (was {}), {} edges", 
               boost::num_vertices(steiner_result.graph), boost::num_vertices(base_graph),
               boost::num_edges(steiner_result.graph));


    // return steiner_result.graph;

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
    const auto& ref_time_blob_map = reference_cluster->time_blob_map(); // this one has the time blob map ...
    
    if (ref_time_blob_map.empty()) {
        SPDLOG_LOGGER_DEBUG(log, "filter_by_reference_cluster: reference cluster has empty time_blob_map");
        return terminals;
    }

    // Filter terminals based on spatial relationship with reference cluster
    for (auto terminal_idx : terminals) {
        //std::cout << "Test: " << terminal_idx << " " << ref_time_blob_map.size() << std::endl;
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

    const auto wpids_set = m_cluster.wpids_blob_set();
    std::vector<WireCell::WirePlaneId> wpids(wpids_set.begin(), wpids_set.end());

    std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;

    // Access the detector volumes from the config
    IDetectorVolumes::pointer dv = m_config.dv;

    
    for (const auto& wpid : wpids) {
        int apa = wpid.apa();
        int face = wpid.face();

        // Create wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);
     
        // Get drift direction based on face orientation
        int face_dirx = dv->face_dirx(wpid_u);
        geo_point_t drift_dir(face_dirx, 0, 0);
        
        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        wpid_params[wpid] = std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
    }
    
    auto path_point_cloud = std::make_shared<DynamicPointCloud>(wpid_params);
    const double step_dis = 0.6 * units::cm;

    path_point_cloud->add_points(make_points_cluster_skeleton(&m_cluster, dv, wpid_params, path_point_indices, false, step_dis)); // check path_indices & m_cluster

    // std::cout << path_point_cloud->get_points().size() << " points in path point cloud" << std::endl;

   
    
    vertex_set filtered_terminals;
    
    // Distance thresholds from prototype
    const double distance_3d_threshold = 6.0 * units::cm;
    const double distance_2d_threshold = 1.8 * units::cm;

    for (auto terminal_idx : terminals) {
        Point point = m_cluster.point3d(terminal_idx);
        auto wpid = m_cluster.wpid(terminal_idx);

        // Calculate distances similar to prototype's logic

        auto result_3d = path_point_cloud->kd3d().knn(1, point);
        double dis_3d = sqrt(result_3d[0].second);

        auto dis_2d_u =   path_point_cloud->get_closest_2d_point_info(point, 0, wpid.face(), wpid.apa());
        auto dis_2d_v =   path_point_cloud->get_closest_2d_point_info(point, 1, wpid.face(), wpid.apa());
        auto dis_2d_w =   path_point_cloud->get_closest_2d_point_info(point, 2, wpid.face(), wpid.apa());

        double dis_2d[3]={std::get<0>(dis_2d_u), std::get<0>(dis_2d_v), std::get<0>(dis_2d_w)};

        
        // Apply prototype's filtering logic:
        // Remove if close in 2D projections but far in 3D
        bool close_in_2d = (dis_2d[0] < distance_2d_threshold && dis_2d[1] < distance_2d_threshold) ||
                          (dis_2d[0] < distance_2d_threshold && dis_2d[2] < distance_2d_threshold) ||
                          (dis_2d[1] < distance_2d_threshold && dis_2d[2] < distance_2d_threshold);
        bool should_remove = close_in_2d && (dis_3d > distance_3d_threshold);
        
        // std::cout << "Test1: " << point << " " << dis_3d << " " << dis_2d[0] << " " << dis_2d[1] << " " << dis_2d[2] << std::endl;

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
        SPDLOG_LOGGER_WARN(log, "get_extreme_points_for_reference: failed to get extreme points: {}", e.what());
    }

    return extreme_points;
}





bool Steiner::Grapher::is_point_spatially_related_to_reference(
    size_t point_idx,
    const Facade::Cluster::time_blob_map_t& ref_time_blob_map) const
{
    // Delegate to the cluster's existing method which implements the proper logic
    // for checking spatial relationships with the complex time_blob_map structure
    return m_cluster.is_point_spatially_related_to_time_blobs(point_idx, ref_time_blob_map, true);
}


// These methods need cluster-specific implementation:
size_t Steiner::Grapher::find_closest_vertex_to_point(const Point& point) const
{
    // Find the vertex index closest to the given 3D point
    // check scope ???
    auto closest_idx = m_cluster.get_closest_point_index(point);
    
    return closest_idx;
}



// Overloaded version for specific blobs (equivalent to prototype's mcells parameter)
Steiner::Grapher::vertex_set Steiner::Grapher::find_peak_point_indices(
    const std::vector<const Facade::Blob*>& target_blobs, const std::string& graph_name,
    bool disable_dead_mix_cell, int nlevel)
{
    if (target_blobs.empty()) {
        return {};
    }

    // Build a reverse lookup: Blob* -> node index, so we can query the
    // node-index-keyed cell_points_map using a Blob* from the caller.
    const auto& sv = m_cluster.sv3d();
    const auto& nodes = sv.nodes();
    std::map<const Facade::Blob*, size_t> blob_to_node_idx;
    for (size_t i = 0; i < nodes.size(); ++i) {
        const auto* b = nodes[i]->value.facade<Blob>();
        if (b) blob_to_node_idx[b] = i;
    }

    auto cell_points_map = form_cell_points_map();
    vertex_set all_indices;
    for (const auto* blob : target_blobs) {
        auto rev_it = blob_to_node_idx.find(blob);
        if (rev_it == blob_to_node_idx.end()) continue;
        auto it = cell_points_map.find(rev_it->second);
        if (it != cell_points_map.end()) {
            all_indices.insert(it->second.begin(), it->second.end());
        }
    }
    return find_peak_point_indices(all_indices, graph_name, disable_dead_mix_cell, nlevel);
}

Steiner::Grapher::vertex_set Steiner::Grapher::find_peak_point_indices(
    const vertex_set& all_indices, const std::string& graph_name,
    bool disable_dead_mix_cell, int nlevel)
{
    vertex_set peak_points;

    if (all_indices.empty()) {
        return peak_points;
    }

    // Calculate charges and find candidates
    std::map<size_t, double> map_index_charge;
    std::set<std::pair<double, size_t>, std::greater<std::pair<double, size_t>>> candidates_set;
    
    const double charge_threshold = 4000.0;
    
    for (size_t point_idx : all_indices) {
        auto [charge_quality, charge] = m_cluster.calc_charge_wcp(point_idx, charge_threshold, disable_dead_mix_cell);
        
        map_index_charge[point_idx] = charge;
        
        if (charge > charge_threshold && charge_quality) {
            candidates_set.insert(std::make_pair(charge, point_idx));
        }
    }
    
    // std::cout << "Xin2: candidates_set size: " << candidates_set.size() << std::endl;

    if (candidates_set.empty()) {
        return peak_points;
    }
    
   // Get access to the underlying boost graph
    // NOTE: This assumes there's a way to get the graph - you'll need to implement
    // get_graph() method or similar in your Cluster interface
    const auto& graph = m_cluster.get_graph(graph_name); // You need to implement this
    
    std::set<size_t> peak_indices;
    std::set<size_t> non_peak_indices;
    
    // Core algorithm from prototype: process candidates in order of decreasing charge
    for (const auto& [current_charge, current_index] : candidates_set) {
        
        // Find all vertices within nlevel hops using graph traversal (prototype logic)
        std::set<size_t> total_vertices_found;
        total_vertices_found.insert(current_index);
        
        // Breadth-first exploration for nlevel steps
        std::set<size_t> vertices_to_be_examined;
        vertices_to_be_examined.insert(current_index);
        
        for (int level = 0; level < nlevel; ++level) {
            std::set<size_t> vertices_saved_for_next;
            
            for (size_t temp_current_index : vertices_to_be_examined) {
                // Get adjacent vertices using boost graph interface
                auto [neighbors_begin, neighbors_end] = boost::adjacent_vertices(temp_current_index, graph);
                
                for (auto neighbor_it = neighbors_begin; neighbor_it != neighbors_end; ++neighbor_it) {
                    size_t neighbor_index = *neighbor_it;
                    
                    if (total_vertices_found.find(neighbor_index) == total_vertices_found.end()) {
                        total_vertices_found.insert(neighbor_index);
                        vertices_saved_for_next.insert(neighbor_index);
                    }
                }
            }
            vertices_to_be_examined = vertices_saved_for_next;
        }
        total_vertices_found.erase(current_index);
        
        // Peak selection logic (following prototype)
        if (peak_indices.empty()) {
            // First candidate becomes a peak
            peak_indices.insert(current_index);
            
            // Mark neighbors with lower charge as non-peaks
            for (size_t neighbor_idx : total_vertices_found) {
                if (map_index_charge.find(neighbor_idx) == map_index_charge.end()) continue;
                
                if (current_charge > map_index_charge[neighbor_idx]) {
                    non_peak_indices.insert(neighbor_idx);
                }
            }
        } else {
            // Skip if already classified
            if (peak_indices.find(current_index) != peak_indices.end() ||
                non_peak_indices.find(current_index) != non_peak_indices.end()) {
                continue;
            }
            
            bool flag_insert = true;
            
            // Check against neighbors
            for (size_t neighbor_idx : total_vertices_found) {
                if (map_index_charge.find(neighbor_idx) == map_index_charge.end()) continue;
                
                if (current_charge > map_index_charge[neighbor_idx]) {
                    non_peak_indices.insert(neighbor_idx);
                } else if (current_charge < map_index_charge[neighbor_idx]) {
                    flag_insert = false;
                    break;
                }
            }
            
            if (flag_insert) {
                peak_indices.insert(current_index);
            }
        }
    }
    
    // Connected components analysis to merge nearby peaks (prototype logic)
    if (peak_indices.size() > 1) {
        std::vector<size_t> vec_peak_indices(peak_indices.begin(), peak_indices.end());
        peak_indices.clear();
        
        const size_t N = vec_peak_indices.size();
        
        // Create temporary graph for peak connectivity
        boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
            boost::no_property, boost::property<boost::edge_weight_t, double>>
            temp_graph(N);
        
        // Check connectivity in original graph and replicate in temp graph.
        // Use k = j+1 to visit each undirected pair once (halves boost::edge() calls).
        for (size_t j = 0; j < N; ++j) {
            for (size_t k = j + 1; k < N; ++k) {
                if (boost::edge(vec_peak_indices[j], vec_peak_indices[k], graph).second) {
                    boost::add_edge(j, k, temp_graph);
                }
            }
        }
        
        // Find connected components
        std::vector<int> component(boost::num_vertices(temp_graph));
        const int num = boost::connected_components(temp_graph, &component[0]);
        
        // For each component, find the point closest to center of mass
        std::vector<double> min_dis(num, 1e9);
        std::vector<WireCell::Point> points(num, WireCell::Point(0, 0, 0));
        std::vector<size_t> min_index(num, 0);
        std::vector<int> ncounts(num, 0);
        
        // Calculate center of mass for each component
        for (size_t i = 0; i < component.size(); ++i) {
            ncounts[component[i]]++;
            auto point = m_cluster.point3d(vec_peak_indices[i]);
            points[component[i]] = points[component[i]] + WireCell::Vector(point.x(), point.y(), point.z());
        }
        
        // Average the positions
        for (int i = 0; i < num; ++i) {
            if (ncounts[i] > 0) {
                points[i] = points[i] * (1.0 / ncounts[i]);
            }
        }
        
        // Find closest point to center for each component
        for (size_t i = 0; i < component.size(); ++i) {
            auto point = m_cluster.point3d(vec_peak_indices[i]);
            
            double dis = pow(points[component[i]].x() - point.x(), 2) +
                        pow(points[component[i]].y() - point.y(), 2) +
                        pow(points[component[i]].z() - point.z(), 2);
            
            if (dis < min_dis[component[i]]) {
                min_dis[component[i]] = dis;
                min_index[component[i]] = vec_peak_indices[i];
            }
        }
        
        // Add the representative points to the result
        for (int i = 0; i < num; ++i) {
            peak_indices.insert(min_index[i]);
        }
    }
    
    return peak_indices;
}


Steiner::Grapher::vertex_set Steiner::Grapher::find_steiner_terminals(const std::string& graph_name, bool disable_dead_mix_cell)
{
    auto cell_points_map = form_cell_points_map();
    return find_steiner_terminals(graph_name, disable_dead_mix_cell, cell_points_map);
}

Steiner::Grapher::vertex_set Steiner::Grapher::find_steiner_terminals(const std::string& graph_name, bool disable_dead_mix_cell, const blob_vertex_map& cell_points_map)
{
    vertex_set steiner_terminals;
    
    if (cell_points_map.empty()) {
        return steiner_terminals;
    }
    
    // Process each blob individually (following prototype pattern).
    // Iterating by blob node index (size_t key) gives deterministic order.
    for (const auto& [blob_node_idx, point_indices] : cell_points_map) {
        (void)blob_node_idx;  // key used only for ordering; value suffices here
        auto blob_peaks = find_peak_point_indices(point_indices, graph_name, disable_dead_mix_cell);
        steiner_terminals.insert(blob_peaks.begin(), blob_peaks.end());
    }

    return steiner_terminals;
}



Steiner::Grapher::blob_vertex_map Steiner::Grapher::form_cell_points_map()
{
    blob_vertex_map cell_points;

    // Get the 3D scoped view with x,y,z coordinates
    const auto& sv = m_cluster.sv3d();
    const auto& nodes = sv.nodes();   // These are the blob nodes
    const auto& skd = sv.kd();        // K-d tree with point data

    // Check if we have valid data
    if (nodes.empty() || skd.npoints() == 0) {
        return cell_points;
    }

    // The major indices tell you which blob (node index) each point belongs to.
    // We key the map on the node index (a stable size_t from the traversal-ordered
    // nodes vector) rather than the Blob* address to guarantee deterministic
    // iteration order across runs.
    const auto& majs = skd.major_indices();

    for (size_t point_idx = 0; point_idx < skd.npoints(); ++point_idx) {
        size_t blob_node_idx = majs[point_idx];
        if (blob_node_idx >= nodes.size()) continue;

        // Verify the facade is valid (guards downstream null dereferences)
        const auto* blob = nodes[blob_node_idx]->value.facade<Blob>();
        if (!blob) continue;

        // operator[] default-constructs the vertex_set on first access;
        // no need for a separate find() guard.
        cell_points[blob_node_idx].insert(point_idx);
    }

    return cell_points;
}

void Steiner::Grapher::establish_same_blob_steiner_edges(const std::string& graph_name, 
                                                        bool disable_dead_mix_cell)
{
    using Clock = std::chrono::steady_clock;
    // using MS = std::chrono::duration<double, std::milli>;
    auto t0 = Clock::now();

    if (!m_cluster.has_graph(graph_name)) {
        SPDLOG_LOGGER_ERROR(log, "Graph '{}' does not exist in cluster", graph_name);
        return;
    }

    auto& graph = m_cluster.get_graph(graph_name);
    std::vector<edge_type> added_edges;

    // Compute blob->points map once and reuse it for both terminal finding and edge adding
    auto cell_points_map = form_cell_points_map();
    // if (m_perf) std::cout << "establish_same_blob_steiner_edges[" << graph_name << "] timing: form_cell_points_map took " << MS(Clock::now()-t0).count() << " ms  blobs=" << cell_points_map.size() << std::endl;

    if (cell_points_map.empty()) {
        SPDLOG_LOGGER_WARN(log, "No blob-to-points mapping available for Steiner edge establishment");
        return;
    }

    // Step 1: Find Steiner terminals, passing the already-computed map to avoid
    // calling form_cell_points_map() a second time inside find_steiner_terminals().
    t0 = Clock::now();
    vertex_set steiner_terminals = find_steiner_terminals(graph_name, disable_dead_mix_cell, cell_points_map);
    // if (m_perf) std::cout << "establish_same_blob_steiner_edges[" << graph_name << "] timing: find_steiner_terminals took " << MS(Clock::now()-t0).count() << " ms  terminals=" << steiner_terminals.size() << std::endl;
    
    SPDLOG_LOGGER_DEBUG(log, "Found {} Steiner terminals for same-blob edge establishment", steiner_terminals.size());
    // SPDLOG_LOGGER_DEBUG(log, "Processing {} blobs for same-blob edges", cell_points_map.size());

    // Build an O(1)-lookup set from the terminals (std::set lookup is O(log N)).
    std::unordered_set<vertex_type> terminal_set(steiner_terminals.begin(), steiner_terminals.end());

    // Step 2: For each blob, add edges between all pairs of points that include
    // at least one Steiner terminal.  Cache 3D coordinates per blob to avoid
    // repeated point3d() lookups inside the O(M^2) inner loop.
    t0 = Clock::now();
    // Iterate blobs in node-index order (deterministic) and add intra-blob edges.
    for (const auto& [blob_node_idx, point_indices] : cell_points_map) {
        (void)blob_node_idx;
        if (point_indices.size() < 2) continue;

        // Convert set to vector for O(1) index access
        std::vector<vertex_type> points_vec(point_indices.begin(), point_indices.end());
        const size_t M = points_vec.size();

        // --- Early-exit: skip blobs with zero terminals (avoids O(M²) loop entirely) ---
        size_t n_terminals = 0;
        size_t sole_term_k = M;  // index of the single terminal if n_terminals==1
        for (size_t k = 0; k < M; ++k) {
            if (terminal_set.count(points_vec[k])) {
                ++n_terminals;
                sole_term_k = k;
                if (n_terminals > 1) break;
            }
        }
        if (n_terminals == 0) continue;

        // Cache 3D coordinates — avoids 2 × M² point3d() calls in the loop below
        std::vector<Point> pts(M);
        for (size_t k = 0; k < M; ++k) {
            pts[k] = m_cluster.point3d(points_vec[k]);
        }

        // --- Single-terminal fast path: only M-1 pairs possible ---
        if (n_terminals == 1) {
            vertex_type term_v = points_vec[sole_term_k];
            const auto& p_term = pts[sole_term_k];
            for (size_t j = 0; j < M; ++j) {
                if (j == sole_term_k) continue;
                vertex_type other_v = points_vec[j];
                if (!boost::edge(term_v, other_v, graph).second) {
                    const auto& p_other = pts[j];
                    double dx = p_other.x()-p_term.x(), dy = p_other.y()-p_term.y(), dz = p_other.z()-p_term.z();
                    double w = 0.9 * std::sqrt(dx*dx + dy*dy + dz*dz);
                    auto [e, ok] = boost::add_edge(term_v, other_v, w, graph);
                    if (ok) added_edges.push_back(e);
                }
            }
            continue;
        }

        // --- General path: multiple terminals, O(M²) pairs ---
        for (size_t i = 0; i < M; ++i) {
            vertex_type index1 = points_vec[i];
            bool flag_index1 = (terminal_set.count(index1) > 0);

            for (size_t j = i + 1; j < M; ++j) {
                vertex_type index2 = points_vec[j];
                bool flag_index2 = (terminal_set.count(index2) > 0);

                double factor;
                if (flag_index1 && flag_index2)       factor = 0.8;
                else if (flag_index1 || flag_index2)  factor = 0.9;
                else continue;  // neither is a terminal — skip

                if (!boost::edge(index1, index2, graph).second) {
                    const auto& p1 = pts[i];
                    const auto& p2 = pts[j];
                    double dx = p2.x()-p1.x(), dy = p2.y()-p1.y(), dz = p2.z()-p1.z();
                    double w = factor * std::sqrt(dx*dx + dy*dy + dz*dz);
                    auto [e, ok] = boost::add_edge(index1, index2, w, graph);
                    if (ok) added_edges.push_back(e);
                }
            }
        }
    }
    // if (m_perf) std::cout << "establish_same_blob_steiner_edges[" << graph_name << "] timing: edge loop took " << MS(Clock::now()-t0).count() << " ms  pairs_checked=" << total_pairs_checked << "  edges_added=" << total_edges_added << std::endl;

    // Store the added edges for later removal
    t0 = Clock::now();
    store_added_edges(graph_name, added_edges);
    invalidate_graph_algorithms_cache(graph_name);
    // if (m_perf) std::cout << "establish_same_blob_steiner_edges[" << graph_name << "] timing: store+invalidate took " << MS(Clock::now()-t0).count() << " ms" << std::endl;

    // SPDLOG_LOGGER_INFO(log, "Added {} same-blob edges to graph '{}' from {} total points ({} steiner terminals)", 
    //          added_edges.size(), graph_name, 
    //          std::accumulate(cell_points_map.begin(), cell_points_map.end(), 0,
    //                        [](int sum, const auto& pair) { return sum + pair.second.size(); }),
    //          steiner_terminals.size());
}


void Steiner::Grapher::remove_same_blob_steiner_edges(const std::string& graph_name)
{
    if (!m_cluster.has_graph(graph_name)) {
        SPDLOG_LOGGER_WARN(log, "Graph '{}' does not exist, cannot remove edges", graph_name);
        return;
    }

    auto it = m_added_edges_by_graph.find(graph_name);
    if (it == m_added_edges_by_graph.end() || it->second.empty()) {
        SPDLOG_LOGGER_DEBUG(log, "No edges to remove for graph '{}'", graph_name);
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

    (void)removed_count;
    // SPDLOG_LOGGER_INFO(log, "Removed {} same-blob Steiner edges from graph '{}'", removed_count, graph_name);
}

void Steiner::Grapher::invalidate_graph_algorithms_cache(const std::string& graph_name)
{
    // Use the new public method we'll add to Cluster
    m_cluster.clear_graph_algorithms_cache(graph_name);
}

void Steiner::Grapher::store_added_edges(const std::string& graph_name, const std::vector<edge_type>& edges)
{
    auto& tracked_edges = m_added_edges_by_graph[graph_name];
    tracked_edges.insert(tracked_edges.end(), edges.begin(), edges.end());
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



namespace WireCell::Clus::Graphs::Weighted{


double calculate_charge_weighted_distance(
    double geometric_distance,
    double charge_source,
    double charge_target, 
    const ChargeWeightingConfig& config)
{
    if (!config.enable_weighting) {
        return geometric_distance;
    }
    
    // Apply prototype charge weighting formula
    double weight_factor = config.factor1 + config.factor2 * 
        (0.5 * config.Q0 / (charge_source + config.Q0) + 
         0.5 * config.Q0 / (charge_target + config.Q0));
    
    return geometric_distance * weight_factor;
}

std::map<Weighted::vertex_type, double> calculate_vertex_charges(
    const vertex_set& vertices, 
    const PointCloud::Dataset& pc,
    const WireCell::Clus::Facade::Cluster& cluster,
    double charge_cut = 4000.0,
    bool disable_dead_mix_cell = true)
{
   std::map<vertex_type, double> charges;
    
    if (pc.size_major() == 0) {
        return charges;
    }
    
    // Calculate charge for each vertex using the existing cluster method
    for (auto vtx : vertices) {
        if (vtx < pc.size_major()) {
            // Use the existing calc_charge_wcp method from Facade::Cluster
            auto charge_result = cluster.calc_charge_wcp(vtx, charge_cut, disable_dead_mix_cell);
            charges[vtx] = charge_result.second;
        }
    }
    
    return charges;
}

void establish_same_blob_steiner_edges_steiner_graph(EnhancedSteinerResult& result, const WireCell::Clus::Facade::Cluster& cluster) {
    // Key by blob node index (size_t, traversal-ordered) for deterministic iteration.
    using blob_vertex_map = std::map<size_t, std::set<size_t>>;
    blob_vertex_map cell_points_map;

    const auto& sv = cluster.sv3d();
    const auto& nodes = sv.nodes();
    const auto& skd = sv.kd();
    const auto& majs = skd.major_indices();

    for (const auto& [old_index, new_index] : result.old_to_new_index) {
        size_t blob_node_idx = majs[old_index];
        if (blob_node_idx >= nodes.size()) continue;
        // Null guard: skip if blob facade is invalid
        const auto* blob = nodes[blob_node_idx]->value.facade<Facade::Blob>();
        if (!blob) continue;
        cell_points_map[blob_node_idx].insert(new_index);
    }

    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_arr = result.point_cloud.get(coords.at(0))->elements<double>();
    const auto& y_arr = result.point_cloud.get(coords.at(1))->elements<double>();
    const auto& z_arr = result.point_cloud.get(coords.at(2))->elements<double>();

    for (const auto& [blob_node_idx, point_indices] : cell_points_map) {
        (void)blob_node_idx;
        std::vector<vertex_type> points_vec(point_indices.begin(), point_indices.end());

        for (size_t i = 0; i < points_vec.size(); ++i) {
            vertex_type index1 = points_vec[i];
            Point point1(x_arr[index1], y_arr[index1], z_arr[index1]);
            for (size_t j = i + 1; j < points_vec.size(); ++j) {
                vertex_type index2 = points_vec[j];
                Point point2(x_arr[index2], y_arr[index2], z_arr[index2]);

                if (result.flag_steiner_terminal[index1] || result.flag_steiner_terminal[index2]) {
                    double distance = (point1-point2).magnitude();
                    if (!boost::edge(index1, index2, result.graph).second) {
                        boost::add_edge(index1, index2, distance, result.graph);
                    }
                }
            }
        }
    }
}


// Updated enhanced steiner graph function to use the new charge calculation
EnhancedSteinerResult create_enhanced_steiner_graph(
            const graph_type& base_graph,
            const vertex_set& terminal_vertices,
            const PointCloud::Dataset& original_pc,
            const WireCell::Clus::Facade::Cluster& cluster, // Added cluster parameter
            const ChargeWeightingConfig& charge_config,
            bool disable_dead_mix_cell
        )
{
    using namespace WireCell::Clus::Graphs::Weighted;
    
    EnhancedSteinerResult result;
    
    // Step 1: Create Voronoi tessellation
    std::vector<vertex_type> terminal_vector(terminal_vertices.begin(), terminal_vertices.end());
    auto vor = voronoi(base_graph, terminal_vector);
    
    // Step 2: Build complete terminal distance map (matches prototype map_saved_edge)
    auto edge_weight = get(boost::edge_weight, base_graph);
    std::map<vertex_pair, std::pair<double, edge_type>> map_saved_edge;
    std::vector<edge_type> all_terminal_connecting_edges;
    
    // Find best edges between all terminal pairs (matches prototype logic exactly)
    auto [edge_iter, edge_end] = boost::edges(base_graph);
    for (auto fine_edge : boost::make_iterator_range(edge_iter, edge_end)) {
        const vertex_type fine_tail = boost::source(fine_edge, base_graph);
        const vertex_type fine_head = boost::target(fine_edge, base_graph);
        const double fine_distance = edge_weight[fine_edge];
        
        const vertex_type term_tail = vor.terminal[fine_tail];
        const vertex_type term_head = vor.terminal[fine_head];
        
        // Skip edges within same terminal region
        if (term_tail == term_head) {
            continue;
        }
        
        // Calculate total distance: path_to_terminal + edge + path_to_terminal
        const double total_distance = vor.distance[fine_tail] + fine_distance + vor.distance[fine_head];
        const vertex_pair term_vp = make_vertex_pair(term_tail, term_head);
        
        // Check if this is the best edge for this terminal pair (matches prototype logic)
        auto it = map_saved_edge.find(term_vp);
        if (it == map_saved_edge.end()) {
            // Try reverse pair
            vertex_pair reverse_vp = make_vertex_pair(term_head, term_tail);
            auto reverse_it = map_saved_edge.find(reverse_vp);
            if (reverse_it == map_saved_edge.end()) {
                // First edge for this terminal pair
                map_saved_edge[term_vp] = std::make_pair(total_distance, fine_edge);
            } else if (total_distance < reverse_it->second.first) {
                // Better than existing reverse pair
                map_saved_edge.erase(reverse_it);
                map_saved_edge[term_vp] = std::make_pair(total_distance, fine_edge);
            }
        } else if (total_distance < it->second.first) {
            // Better than existing edge for this pair
            it->second = std::make_pair(total_distance, fine_edge);
        }
    }
    
    // Step 3: Extract all selected terminal connecting edges (matches prototype terminal_edge)
    for (const auto& [term_pair, edge_info] : map_saved_edge) {
        all_terminal_connecting_edges.push_back(edge_info.second);
    }
    
    // std::cout << "Terminal connecting edges: " << total.size() << std::endl;
    
    // Step 4: Build complete edge set by including all paths (matches prototype unique_edges logic).
    // Store as (vertex_pair, edge_type) so we can deduplicate and sort by stable vertex-index
    // values rather than by edge_type (which embeds a pointer — non-deterministic across runs).
    vertex_set selected_vertices;
    std::vector<std::pair<vertex_pair, edge_type>> tree_edge_pairs;

    for (auto edge : all_terminal_connecting_edges) {
        auto vs = boost::source(edge, base_graph);
        auto vt = boost::target(edge, base_graph);
        tree_edge_pairs.push_back({make_vertex_pair(vs, vt), edge});

        for (auto endpoint : {vs, vt}) {
            vertex_type current_vtx = endpoint;
            while (vor.terminal[current_vtx] != current_vtx) {
                auto path_edge = vor.last_edge[current_vtx];
                auto ps = boost::source(path_edge, base_graph);
                auto pt = boost::target(path_edge, base_graph);
                tree_edge_pairs.push_back({make_vertex_pair(ps, pt), path_edge});
                current_vtx = ps;
            }
        }
    }

    // Step 5: Sort by vertex_pair (deterministic size_t comparison) and deduplicate.
    std::sort(tree_edge_pairs.begin(), tree_edge_pairs.end(),
        [](const auto& a, const auto& b){ return a.first < b.first; });
    tree_edge_pairs.erase(
        std::unique(tree_edge_pairs.begin(), tree_edge_pairs.end(),
            [](const auto& a, const auto& b){ return a.first == b.first; }),
        tree_edge_pairs.end());

    // Collect all vertices from the unique edges
    for (const auto& [vp, edge] : tree_edge_pairs) {
        selected_vertices.insert(vp.first);
        selected_vertices.insert(vp.second);
    }
    
    // Step 6: Create index mappings (same as before)
    std::vector<vertex_type> selected_vector(selected_vertices.begin(), selected_vertices.end());
    for (size_t i = 0; i < selected_vector.size(); ++i) {
        result.old_to_new_index[selected_vector[i]] = i;
        result.new_to_old_index[i] = selected_vector[i];
    }
    
    // Step 7: Create flag_steiner_terminal (same as before)
    result.flag_steiner_terminal.resize(selected_vector.size());
    for (size_t i = 0; i < selected_vector.size(); ++i) {
        vertex_type old_idx = selected_vector[i];
        result.flag_steiner_terminal[i] = (terminal_vertices.find(old_idx) != terminal_vertices.end());
    }
    
    // Step 8: Calculate charges for vertices (same as before)
    if (original_pc.size_major() > 0 && charge_config.enable_weighting) {
        result.vertex_charges = calculate_vertex_charges(
            selected_vertices, 
            original_pc, 
            cluster,
            4000.0,  // charge_cut from prototype 
            disable_dead_mix_cell
        );
    }
    
    // Step 9: Create subset point cloud (same as before)
    if (original_pc.size_major() > 0) {
        std::vector<size_t> subset_indices(selected_vertices.begin(), selected_vertices.end());
        result.point_cloud = original_pc.subset(subset_indices);

        // Add the flag_steiner_terminal boolean array to point cloud
        // Store as int to match all read sites that use elements<int>()
        std::vector<int> steiner_flags_int(result.flag_steiner_terminal.begin(), result.flag_steiner_terminal.end());
        PointCloud::Array steiner_flag_array(steiner_flags_int);
        result.point_cloud.add("flag_steiner_terminal", std::move(steiner_flag_array));
    }
    
    // Step 10: Create reduced graph with ALL unique edges.
    // Iterate tree_edge_pairs (sorted by vertex_pair) so edge insertion order is deterministic.
    result.graph = graph_type(selected_vector.size());

    for (const auto& [vp, edge] : tree_edge_pairs) {
        vertex_type old_source = vp.first;
        vertex_type old_target = vp.second;

        if (result.old_to_new_index.find(old_source) == result.old_to_new_index.end() ||
            result.old_to_new_index.find(old_target) == result.old_to_new_index.end()) {
            continue;
        }

        vertex_type new_source = result.old_to_new_index[old_source];
        vertex_type new_target = result.old_to_new_index[old_target];

        double geometric_distance = edge_weight[edge];

        double final_distance = geometric_distance;
        if (charge_config.enable_weighting && !result.vertex_charges.empty()) {
            double charge_source = result.vertex_charges.count(old_source) ?
                                  result.vertex_charges[old_source] : 0.0;
            double charge_target = result.vertex_charges.count(old_target) ?
                                  result.vertex_charges[old_target] : 0.0;
            // Prototype formula: dis * (factor1 + factor2 * (0.5*Q0/(Qs+Q0) + 0.5*Q0/(Qt+Q0)))
            double weight_factor = charge_config.factor1 + charge_config.factor2 *
                (0.5 * charge_config.Q0 / (charge_source + charge_config.Q0) +
                 0.5 * charge_config.Q0 / (charge_target + charge_config.Q0));
            final_distance = geometric_distance * weight_factor;
        }
        if (!boost::edge(new_source, new_target, result.graph).second) {
            boost::add_edge(new_source, new_target, final_distance, result.graph);
        }
    }
    
    // std::cout << "Final graph - vertices: " << boost::num_vertices(result.graph) 
    //           << ", edges: " << boost::num_edges(result.graph) << std::endl;
    
    result.steiner_terminal_indices = terminal_vertices;
    return result;
}

}
