#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
// #include "WireCellClus/Graphs/Weighted.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

// Edge property tag for Boost Graph
struct edge_base_t {
    typedef boost::edge_property_tag kind;
};

// Helper struct to track segment candidates
struct Res_proto_segment {
    int group_num;
    int number_points;
    size_t special_A;
    size_t special_B;
    double length;
    int number_not_faked;
    double max_dis_u;
    double max_dis_v;
    double max_dis_w;
};

void PatternAlgorithms::find_other_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_break_track, double search_range, double scaling_2d)
{
    // Get steiner point cloud data
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    const auto& wpid_array = steiner_pc.get("wpid")->elements<WirePlaneId>();
    
    const size_t N = x_coords.size();
    if (N == 0) return;
    
    // Step 1: Tag points near existing segments
    std::vector<bool> flag_tagged(N, false);
    // int num_tagged = 0;
    
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    
    // Get all existing segments in this cluster
    std::set<SegmentPtr> existing_segments = find_cluster_segments(graph, cluster);
    
    for (size_t i = 0; i < N; i++) {
        Facade::geo_point_t p(x_coords[i], y_coords[i], z_coords[i]);
        double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
        double min_3d_dis = 1e9;
        
        WirePlaneId wpid = wpid_array[i];
        int apa = wpid.apa();
        int face = wpid.face();
        
        // Check distances to existing segments
        for (auto seg : existing_segments) {
            // Get closest 3D point
            auto closest_result = segment_get_closest_point(seg, p, "fit");
            double dis_3d = closest_result.first;  // distance is already computed
            
            if (dis_3d < min_3d_dis) min_3d_dis = dis_3d;
            
            if (dis_3d < search_range) {
                flag_tagged[i] = true;
                // num_tagged++;
                break;
            }
            
            // Get 2D distances
            auto closest_2d = segment_get_closest_2d_distances(seg, p, apa, face, "fit");
            double dis_u = std::get<0>(closest_2d);
            double dis_v = std::get<1>(closest_2d);
            double dis_w = std::get<2>(closest_2d);
            
            if (dis_u < min_dis_u) min_dis_u = dis_u;
            if (dis_v < min_dis_v) min_dis_v = dis_v;
            if (dis_w < min_dis_w) min_dis_w = dis_w;
        }
        
        // Additional tagging based on 2D projections and dead channels
        if (!flag_tagged[i]) {
            auto p_raw = transform->backward(p, cluster_t0, face, apa);
            
            bool u_ok = (min_dis_u < scaling_2d * search_range || 
                        cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0));
            bool v_ok = (min_dis_v < scaling_2d * search_range || 
                        cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1));
            bool w_ok = (min_dis_w < scaling_2d * search_range || 
                        cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2));
            
            if (u_ok && v_ok && w_ok) {
                flag_tagged[i] = true;
            }
        }
    }
    
    // Step 2: Get terminal vertices
    const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();
    std::vector<size_t> terminals;
    std::map<size_t, size_t> map_oindex_tindex;
    
    for (size_t i = 0; i < flag_steiner_terminal.size(); i++) {
        if (flag_steiner_terminal[i]) {
            map_oindex_tindex[i] = terminals.size();
            terminals.push_back(i);
        }
    }
    
    if (terminals.empty()) return;
    
    // Step 3: Compute Voronoi diagram
    const auto& steiner_graph = cluster.get_graph("steiner_graph");
    using namespace Graphs::Weighted;
    auto vor = voronoi(steiner_graph, terminals);
    
    // Step 4: Build terminal graph with MST
    using Base = boost::property<edge_base_t, edge_type>;
    using WeightProperty = boost::property<boost::edge_weight_t, double, Base>;
    using TerminalGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                                 boost::no_property, WeightProperty>;
    
    TerminalGraph terminal_graph(N);
    std::map<std::pair<size_t, size_t>, std::pair<double, edge_type>> map_saved_edge;
    
    auto edge_weight = get(boost::edge_weight, steiner_graph);
    
    for (auto w : boost::make_iterator_range(edges(steiner_graph))) {
        size_t nearest_to_source = vor.terminal[source(w, steiner_graph)];
        size_t nearest_to_target = vor.terminal[target(w, steiner_graph)];
        
        if (nearest_to_source != nearest_to_target) {
            double weight = vor.distance[source(w, steiner_graph)] + 
                          vor.distance[target(w, steiner_graph)] + 
                          edge_weight[w];
            
            auto edge_pair1 = std::make_pair(nearest_to_source, nearest_to_target);
            auto edge_pair2 = std::make_pair(nearest_to_target, nearest_to_source);
            
            auto it1 = map_saved_edge.find(edge_pair1);
            auto it2 = map_saved_edge.find(edge_pair2);
            
            if (it1 != map_saved_edge.end()) {
                if (weight < it1->second.first) {
                    it1->second = std::make_pair(weight, w);
                }
            } else if (it2 != map_saved_edge.end()) {
                if (weight < it2->second.first) {
                    it2->second = std::make_pair(weight, w);
                }
            } else {
                map_saved_edge[edge_pair1] = std::make_pair(weight, w);
            }
        }
    }
    
    // Add edges to terminal graph
    for (const auto& [edge_pair, weight_info] : map_saved_edge) {
        boost::add_edge(edge_pair.first, edge_pair.second,
                       WeightProperty(weight_info.first, Base(weight_info.second)),
                       terminal_graph);
    }
    
    // Step 5: Find minimum spanning tree
    std::vector<boost::graph_traits<TerminalGraph>::edge_descriptor> mst_edges;
    boost::kruskal_minimum_spanning_tree(terminal_graph, std::back_inserter(mst_edges));
    
    // Step 6: Create cluster graph based on tagging
    TerminalGraph terminal_graph_cluster(terminals.size());
    std::map<size_t, std::set<size_t>> map_connection;
    
    for (const auto& edge : mst_edges) {
        size_t source_idx = boost::source(edge, terminal_graph);
        size_t target_idx = boost::target(edge, terminal_graph);
        
        if (flag_tagged[source_idx] == flag_tagged[target_idx]) {
            boost::add_edge(map_oindex_tindex[source_idx], 
                          map_oindex_tindex[target_idx], 
                          terminal_graph_cluster);
        } else {
            if (map_connection.find(source_idx) == map_connection.end()) {
                std::set<size_t> temp_results;
                temp_results.insert(target_idx);
                map_connection[source_idx] = temp_results;
            } else {
                map_connection[source_idx].insert(target_idx);
            }
            
            if (map_connection.find(target_idx) == map_connection.end()) {
                std::set<size_t> temp_results;
                temp_results.insert(source_idx);
                map_connection[target_idx] = temp_results;
            } else {
                map_connection[target_idx].insert(source_idx);
            }
        }
    }
    
    // Step 7: Find connected components
    std::vector<int> component(boost::num_vertices(terminal_graph_cluster));
    const int num_components = boost::connected_components(terminal_graph_cluster, &component[0]);
    
    std::vector<int> ncounts(num_components, 0);
    std::vector<std::vector<size_t>> sep_clusters(num_components);
    
    for (size_t i = 0; i < component.size(); ++i) {
        ncounts[component[i]]++;
        sep_clusters[component[i]].push_back(terminals[i]);
    }
    
    // Step 8: Analyze each cluster and filter
    std::vector<Res_proto_segment> temp_segments(num_components);
    std::set<int> remaining_segments;
    
    for (int i = 0; i < num_components; i++) {
        // Skip if inside original track
        if (flag_tagged[sep_clusters[i].front()]) {
            continue;
        }
        
        remaining_segments.insert(i);
        temp_segments[i].group_num = i;
        temp_segments[i].number_points = ncounts[i];
        
        // Find connection point (special_A)
        size_t special_A = SIZE_MAX;
        std::vector<size_t> candidates_special_A;
        
        for (int j = 0; j < ncounts[i]; j++) {
            if (map_connection.find(sep_clusters[i][j]) != map_connection.end()) {
                candidates_special_A.push_back(sep_clusters[i][j]);
            }
        }
        
        if (!candidates_special_A.empty()) {
            special_A = candidates_special_A.front();
        }
        
        // Find furthest point (special_B)
        size_t special_B = special_A;
        double min_dis = 0;
        int number_not_faked = 0;
        double max_dis_u = 0, max_dis_v = 0, max_dis_w = 0;
        
        for (int j = 0; j < ncounts[i]; j++) {
            size_t idx = sep_clusters[i][j];
            double dis = std::sqrt(
                std::pow(x_coords[idx] - x_coords[special_A], 2) +
                std::pow(y_coords[idx] - y_coords[special_A], 2) +
                std::pow(z_coords[idx] - z_coords[special_A], 2));
            
            if (dis > min_dis) {
                min_dis = dis;
                special_B = idx;
            }
            
            // Check if point is fake (too close to existing segments)
            Facade::geo_point_t p(x_coords[idx], y_coords[idx], z_coords[idx]);
            double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
            
            WirePlaneId wpid = wpid_array[idx];
            int apa = wpid.apa();
            int face = wpid.face();
            
            for (auto seg : existing_segments) {
                auto closest_2d = segment_get_closest_2d_distances(seg, p, apa, face, "fit");
                double dis_u = std::get<0>(closest_2d);
                double dis_v = std::get<1>(closest_2d);
                double dis_w = std::get<2>(closest_2d);
                
                if (dis_u < min_dis_u) min_dis_u = dis_u;
                if (dis_v < min_dis_v) min_dis_v = dis_v;
                if (dis_w < min_dis_w) min_dis_w = dis_w;
            }
            
            auto p_raw = transform->backward(p, cluster_t0, face, apa);
            
            int flag_num = 0;
            if (min_dis_u > scaling_2d * search_range && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0)) flag_num++;
            if (min_dis_v > scaling_2d * search_range && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1)) flag_num++;
            if (min_dis_w > scaling_2d * search_range && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2)) flag_num++;
            
            if (min_dis_u > max_dis_u && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0)) max_dis_u = min_dis_u;
            if (min_dis_v > max_dis_v && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1)) max_dis_v = min_dis_v;
            if (min_dis_w > max_dis_w && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2)) max_dis_w = min_dis_w;
            
            if (flag_num >= 2) number_not_faked++;
        }
        
        double length = std::sqrt(
            std::pow(x_coords[special_A] - x_coords[special_B], 2) +
            std::pow(y_coords[special_A] - y_coords[special_B], 2) +
            std::pow(z_coords[special_A] - z_coords[special_B], 2));
        
        // Adjust special_A if length is too short
        if (length < 3 * units::cm && special_A != SIZE_MAX) {
            size_t save_index = special_A;
            double save_dis = 1e9;
            for (auto it1 = map_connection[special_A].begin(); 
                 it1 != map_connection[special_A].end(); it1++) {
                double temp_dis = std::sqrt(
                    std::pow(x_coords[special_A] - x_coords[*it1], 2) +
                    std::pow(y_coords[special_A] - y_coords[*it1], 2) +
                    std::pow(z_coords[special_A] - z_coords[*it1], 2));
                if (temp_dis < save_dis) {
                    save_index = *it1;
                    save_dis = temp_dis;
                }
            }
            special_A = save_index;
            length = std::sqrt(
                std::pow(x_coords[special_A] - x_coords[special_B], 2) +
                std::pow(y_coords[special_A] - y_coords[special_B], 2) +
                std::pow(z_coords[special_A] - z_coords[special_B], 2));
        }
        
        temp_segments[i].special_A = special_A;
        temp_segments[i].special_B = special_B;
        temp_segments[i].length = length;
        temp_segments[i].number_not_faked = number_not_faked;
        temp_segments[i].max_dis_u = max_dis_u;
        temp_segments[i].max_dis_v = max_dis_v;
        temp_segments[i].max_dis_w = max_dis_w;
        
        // Apply quality cuts
        if ((temp_segments[i].number_points == 1) ||
            (number_not_faked == 0 &&
             ((length < 3.5 * units::cm) ||
              (((number_not_faked < 0.25 * temp_segments[i].number_points) ||
               (number_not_faked < 0.4 * temp_segments[i].number_points && length < 7 * units::cm)) &&
              max_dis_u / units::cm < 3 && max_dis_v / units::cm < 3 && max_dis_w / units::cm < 3 &&
              max_dis_u + max_dis_v + max_dis_w < 6 * units::cm)))) {
            remaining_segments.erase(i);
        }
    }
    
    // Step 9: Process remaining segments in order of quality
    std::vector<SegmentPtr> new_segments_for_tracking;
    
    while (!remaining_segments.empty()) {
        // Find the best segment (most non-faked points, then longest)
        double max_number_not_faked = 0;
        double max_length = 0;
        int max_length_cluster = -1;
        
        for (auto it = remaining_segments.begin(); it != remaining_segments.end(); it++) {
            if (temp_segments[*it].number_not_faked > max_number_not_faked) {
                max_length_cluster = *it;
                max_number_not_faked = temp_segments[*it].number_not_faked;
                max_length = temp_segments[*it].length;
            } else if (temp_segments[*it].number_not_faked == max_number_not_faked) {
                if (temp_segments[*it].length > max_length) {
                    max_length_cluster = *it;
                    max_number_not_faked = temp_segments[*it].number_not_faked;
                    max_length = temp_segments[*it].length;
                }
            }
        }
        
        if (max_length_cluster == -1) break;
        
        remaining_segments.erase(max_length_cluster);
        
        // Create new segment
        size_t special_A = temp_segments[max_length_cluster].special_A;
        size_t special_B = temp_segments[max_length_cluster].special_B;
        
        if (special_A == SIZE_MAX || special_B == SIZE_MAX) continue;
        
        // Use Dijkstra to find path
        auto path_indices = cluster.graph_algorithms("steiner_graph").shortest_path(special_A, special_B);
        
        std::vector<Facade::geo_point_t> path_points;
        for (size_t idx : path_indices) {
            path_points.emplace_back(x_coords[idx], y_coords[idx], z_coords[idx]);
        }
        
        if (path_points.size() <= 1) continue;
        
        // Create vertices
        VertexPtr v1 = make_vertex(graph);
        v1->wcpt().point = path_points.front();
        v1->cluster(&cluster);
        
        VertexPtr v2 = make_vertex(graph);
        v2->wcpt().point = path_points.back();
        v2->cluster(&cluster);
        
        // Create segment
        auto new_seg = create_segment_for_cluster(cluster, dv, path_points);
        if (!new_seg) {
            remove_vertex(graph, v1);
            remove_vertex(graph, v2);
            continue;
        }
        
        add_segment(graph, new_seg, v1, v2);
        new_segments_for_tracking.push_back(new_seg);
        existing_segments.insert(new_seg);
        
        // Re-evaluate remaining segments
        std::set<int> tmp_del_set;
        for (auto it = remaining_segments.begin(); it != remaining_segments.end(); it++) {
            temp_segments[*it].number_not_faked = 0;
            temp_segments[*it].max_dis_u = 0;
            temp_segments[*it].max_dis_v = 0;
            temp_segments[*it].max_dis_w = 0;
            
            for (int j = 0; j < ncounts[*it]; j++) {
                size_t idx = sep_clusters[*it][j];
                Facade::geo_point_t p(x_coords[idx], y_coords[idx], z_coords[idx]);
                double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
                
                WirePlaneId wpid = wpid_array[idx];
                int apa = wpid.apa();
                int face = wpid.face();
                
                for (auto seg : existing_segments) {
                    auto closest_2d = segment_get_closest_2d_distances(seg, p, apa, face, "fit");
                    double dis_u = std::get<0>(closest_2d);
                    double dis_v = std::get<1>(closest_2d);
                    double dis_w = std::get<2>(closest_2d);
                    
                    if (dis_u < min_dis_u) min_dis_u = dis_u;
                    if (dis_v < min_dis_v) min_dis_v = dis_v;
                    if (dis_w < min_dis_w) min_dis_w = dis_w;
                }
                
                auto p_raw = transform->backward(p, cluster_t0, face, apa);
                
                int flag_num = 0;
                if (min_dis_u > scaling_2d * search_range && 
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0)) flag_num++;
                if (min_dis_v > scaling_2d * search_range && 
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1)) flag_num++;
                if (min_dis_w > scaling_2d * search_range && 
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2)) flag_num++;
                
                if (flag_num >= 2) temp_segments[*it].number_not_faked++;
                
                if (min_dis_u > temp_segments[*it].max_dis_u && 
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0)) 
                    temp_segments[*it].max_dis_u = min_dis_u;
                if (min_dis_v > temp_segments[*it].max_dis_v && 
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1)) 
                    temp_segments[*it].max_dis_v = min_dis_v;
                if (min_dis_w > temp_segments[*it].max_dis_w && 
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2)) 
                    temp_segments[*it].max_dis_w = min_dis_w;
            }
            
            // Apply quality cuts again
            if ((temp_segments[*it].number_points == 1) ||
                (temp_segments[*it].number_not_faked == 0 &&
                 ((temp_segments[*it].length < 3.5 * units::cm) || (
                  ((temp_segments[*it].number_not_faked < 0.25 * temp_segments[*it].number_points) ||
                   (temp_segments[*it].number_not_faked < 0.4 * temp_segments[*it].number_points && 
                    temp_segments[*it].length < 7 * units::cm)) &&
                  temp_segments[*it].max_dis_u / units::cm < 3 && 
                  temp_segments[*it].max_dis_v / units::cm < 3 && 
                  temp_segments[*it].max_dis_w / units::cm < 3 &&
                  temp_segments[*it].max_dis_u + temp_segments[*it].max_dis_v + 
                  temp_segments[*it].max_dis_w < 6 * units::cm)))) {
                tmp_del_set.insert(*it);
            }
        }
        
        for (auto it = tmp_del_set.begin(); it != tmp_del_set.end(); it++) {
            remaining_segments.erase(*it);
        }
    }
    
    // Step 10: Perform tracking on new segments
    if (!new_segments_for_tracking.empty()) {
        for (auto seg : new_segments_for_tracking) {
            track_fitter.add_segment(seg);
        }
        track_fitter.do_multi_tracking(true, true, true);
        
        // Optionally break long segments if requested
        if (flag_break_track) {
            std::vector<SegmentPtr> segments_to_break(new_segments_for_tracking.begin(), 
                                                       new_segments_for_tracking.end());
            break_segments(graph, track_fitter, dv, segments_to_break);
        }
    }
}
