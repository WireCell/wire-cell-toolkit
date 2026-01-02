#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/FiducialUtils.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

bool WireCell::Clus::PR::PatternAlgorithms::search_for_vertex_activities(Graph& graph, VertexPtr vertex, std::set<SegmentPtr>& segments_set, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double search_range){
    // Get steiner point cloud and terminal flags
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();
    
    // Get transform and grouping for point validation
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) return false;
    
    // Get vertex position
    WireCell::Point vtx_point = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
    
    // Collect directions from existing segments
    std::vector<WireCell::Vector> saved_dirs;
    for (auto seg : segments_set) {
        WireCell::Vector dir = vertex_segment_get_dir(vertex, seg, graph, 5*units::cm);
        if (dir.magnitude() != 0) {
            saved_dirs.push_back(dir);
        }
    }
    
    // Get candidate points within search range
    auto candidate_results = cluster.kd_steiner_radius(search_range, vtx_point, "steiner_pc");
    
    double max_dis = 0;
    size_t max_idx = 0;
    bool found = false;
    
    // First round: look for points with good angular separation and charge
    for (const auto& [idx, dist_sq] : candidate_results) {
        if (!flag_steiner_terminal[idx]) continue;
        
        // Skip if this is the vertex itself
        if (ray_length(Ray{WireCell::Point(x_coords[idx], y_coords[idx], z_coords[idx]), vertex->wcpt().point}) < 0.01*units::cm) continue;
        
        // double dis = std::sqrt(dist_sq);
        WireCell::Point test_p(x_coords[idx], y_coords[idx], z_coords[idx]);
        
        // Find minimum distance to all segments
        double min_dis = 1e9;
        double min_dis_u = 1e9;
        double min_dis_v = 1e9;
        double min_dis_w = 1e9;
        
        auto test_wpid = dv->contained_by(test_p);
        if (test_wpid.face() == -1 || test_wpid.apa() == -1) continue;
        
        for (auto it = boost::edges(graph).first; it != boost::edges(graph).second; ++it) {
            SegmentPtr sg = graph[*it].segment;
            if (!sg || sg->cluster() != &cluster) continue;
            
            auto [dist_3d, closest_pt] = segment_get_closest_point(sg, test_p, "fit");
            if (dist_3d < min_dis) min_dis = dist_3d;
            
            auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(sg, test_p, test_wpid.apa(), test_wpid.face(), "fit");
            if (dist_u < min_dis_u) min_dis_u = dist_u;
            if (dist_v < min_dis_v) min_dis_v = dist_v;
            if (dist_w < min_dis_w) min_dis_w = dist_w;
        }
        
        if (min_dis > 0.6*units::cm && min_dis_u + min_dis_v + min_dis_w > 1.2*units::cm) {
            WireCell::Vector dir(test_p.x() - vtx_point.x(), test_p.y() - vtx_point.y(), test_p.z() - vtx_point.z());
            double sum_angle = 0;
            double min_angle = 1e9;
            
            for (size_t j = 0; j < saved_dirs.size(); j++) {
                double angle = std::acos(dir.dot(saved_dirs[j]) / (dir.magnitude() * saved_dirs[j].magnitude())) / 3.1415926 * 180.0;
                sum_angle += angle;
                if (angle < min_angle) min_angle = angle;
            }
            
            // Get average charge
            double sum_charge = 0;
            int ncount = 0;
            auto test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
            
            for (int plane = 0; plane < 3; plane++) {
                if (!grouping->get_closest_dead_chs(test_p_raw, 1, test_wpid.apa(), test_wpid.face(), plane)) {
                    sum_charge += grouping->get_ave_charge(test_p_raw, 0.3*units::cm, plane, test_wpid.apa(), test_wpid.face());
                    ncount++;
                }
            }
            if (ncount != 0) sum_charge /= ncount;
            
            if ((sum_angle) * (sum_charge + 1e-9) > max_dis) {
                max_dis = (sum_angle) * (sum_charge + 1e-9);
                max_idx = idx;
                found = true;
            }
        }
    }
    
    // Second round: if nothing found, use relaxed criteria
    if (max_dis == 0) {
        for (const auto& [idx, dist_sq] : candidate_results) {
            if (!flag_steiner_terminal[idx]) continue;
            
            // Skip if this is the vertex itself
            if (ray_length(Ray{WireCell::Point(x_coords[idx], y_coords[idx], z_coords[idx]), vertex->wcpt().point}) < 0.01*units::cm) continue;
            
            // double dis = std::sqrt(dist_sq);
            WireCell::Point test_p(x_coords[idx], y_coords[idx], z_coords[idx]);
            
            // Find minimum distance to all segments
            double min_dis = 1e9;
            double min_dis_u = 1e9;
            double min_dis_v = 1e9;
            double min_dis_w = 1e9;
            
            auto test_wpid = dv->contained_by(test_p);
            if (test_wpid.face() == -1 || test_wpid.apa() == -1) continue;
            
            for (auto it = boost::edges(graph).first; it != boost::edges(graph).second; ++it) {
                SegmentPtr sg = graph[*it].segment;
                if (!sg || sg->cluster() != &cluster) continue;
                
                auto [dist_3d, closest_pt] = segment_get_closest_point(sg, test_p, "fit");
                if (dist_3d < min_dis) min_dis = dist_3d;
                
                auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(sg, test_p, test_wpid.apa(), test_wpid.face(), "fit");
                if (dist_u < min_dis_u) min_dis_u = dist_u;
                if (dist_v < min_dis_v) min_dis_v = dist_v;
                if (dist_w < min_dis_w) min_dis_w = dist_w;
            }
            
            if (min_dis > 0.36*units::cm && min_dis_u + min_dis_v + min_dis_w > 0.8*units::cm) {
                // Get average charge
                double sum_charge = 0;
                int ncount = 0;
                auto test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                
                for (int plane = 0; plane < 3; plane++) {
                    if (!grouping->get_closest_dead_chs(test_p_raw, 1, test_wpid.apa(), test_wpid.face(), plane)) {
                        sum_charge += grouping->get_ave_charge(test_p_raw, 0.3*units::cm, plane, test_wpid.apa(), test_wpid.face());
                        ncount++;
                    }
                }
                if (ncount != 0) sum_charge /= ncount;
                
                if (min_dis + (min_dis_u + min_dis_v + min_dis_w) / std::sqrt(3.0) > max_dis && sum_charge > 20000) {
                    max_dis = min_dis + (min_dis_u + min_dis_v + min_dis_w) / std::sqrt(3.0);
                    max_idx = idx;
                    found = true;
                }
            }
        }
    }
    
    // If a good candidate was found, create new vertex and segment
    if (found && max_dis != 0) {
        WireCell::Point max_point(x_coords[max_idx], y_coords[max_idx], z_coords[max_idx]);
        
        // Create new vertex at the found point
        auto v1 = make_vertex(graph);
        WCPoint new_wcp;
        new_wcp.point = max_point;
        v1->wcpt(new_wcp).cluster(&cluster);
        
        // Build path from vertex to new vertex using steiner point cloud
        WireCell::Point mid_p(
            (vertex->wcpt().point.x() + max_point.x()) / 2.0,
            (vertex->wcpt().point.y() + max_point.y()) / 2.0,
            (vertex->wcpt().point.z() + max_point.z()) / 2.0
        );
        
        auto [mid_idx, mid_pt] = cluster.get_closest_wcpoint(mid_p);
        
        std::list<WireCell::Point> wcp_list;
        wcp_list.push_back(vertex->wcpt().point);
        
        if (ray_length(Ray{mid_pt, wcp_list.back()}) > 0.01*units::cm) {
            wcp_list.push_back(mid_pt);
        }
        
        if (ray_length(Ray{max_point, wcp_list.back()}) > 0.01*units::cm) {
            wcp_list.push_back(max_point);
        }
        
        if (wcp_list.size() > 1) {
            std::cout << "Cluster: " << cluster.ident() << " Vertex Activity Found at " << mid_p << std::endl;
            
            // Convert to vector for segment creation
            std::vector<WireCell::Point> path_points(wcp_list.begin(), wcp_list.end());
            
            // Create new segment
            auto sg1 = create_segment_for_cluster(cluster, dv, path_points, 0);
            if (sg1) {
                add_segment(graph, sg1, v1, vertex);
                return true;
            }
        }
    }
    
    return false;
}

std::tuple<bool, int, int> PatternAlgorithms::examine_main_vertex_candidate(Graph& graph, VertexPtr vertex){
    bool flag_in = false;
    int ntracks = 0;
    int nshowers = 0;
    SegmentPtr shower_cand = nullptr;
    SegmentPtr track_cand = nullptr;
    
    // Get all segments connected to this vertex
    if (!vertex || !vertex->descriptor_valid()) {
        return std::make_tuple(flag_in, ntracks, nshowers);
    }
    
    auto vd = vertex->get_descriptor();
    auto edge_range = boost::out_edges(vd, graph);
    
    for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg) continue;
        
        // Check if segment is a shower (has kShowerTrajectory or kShowerTopology flags)
        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                        sg->flags_any(SegmentFlags::kShowerTopology);
        
        if (is_shower) {
            nshowers++;
            shower_cand = sg;
        } else {
            ntracks++;
            track_cand = sg;
        }
        
        // Determine which end of segment connects to vertex
        const auto& wcps = sg->wcpts();
        if (wcps.empty()) continue;
        
        bool flag_start = (ray_length(Ray{wcps.front().point, vertex->wcpt().point}) <
                          ray_length(Ray{wcps.back().point, vertex->wcpt().point}));
        
        // Check if segment is pointing IN to the vertex (strong direction)
        int dir_sign = sg->dirsign();
        bool is_dir_weak = sg->dir_weak();
        
        if (flag_start && dir_sign == -1 && !is_dir_weak) {
            flag_in = true;
            break;
        } else if (!flag_start && dir_sign == 1 && !is_dir_weak) {
            flag_in = true;
            break;
        }
    }
    
    // Check Michel electron case: 2 segments (1 track + 1 shower)
    int num_segments = 0;
    auto edge_range2 = boost::out_edges(vd, graph);
    for (auto eit = edge_range2.first; eit != edge_range2.second; ++eit) {
        if (graph[*eit].segment) num_segments++;
    }
    
    if (num_segments == 2 && ntracks == 1 && nshowers == 1 && track_cand && shower_cand) {
        // Calculate the number of daughter showers
        auto pair_result = calculate_num_daughter_showers(graph, vertex, shower_cand);
        
        if (pair_result.first <= 3 && pair_result.second < 30 * units::cm) {
            const auto& track_wcps = track_cand->wcpts();
            if (!track_wcps.empty()) {
                bool flag_start = (ray_length(Ray{track_wcps.front().point, vertex->wcpt().point}) <
                                  ray_length(Ray{track_wcps.back().point, vertex->wcpt().point}));
                
                int track_dir = track_cand->dirsign();
                if ((flag_start && track_dir == -1) || (!flag_start && track_dir == 1)) {
                    flag_in = true;
                }
            }
        }
    }
    
    return std::make_tuple(flag_in, ntracks, nshowers);
}

VertexPtr PatternAlgorithms::compare_main_vertices_all_showers(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    if (vertex_candidates.empty()) return nullptr;
    
    VertexPtr temp_main_vertex = vertex_candidates.front();
    
    // Collect all points from segments and vertices in the cluster
    std::vector<Facade::geo_point_t> pts;
    
    // Collect points from segments
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        
        const auto& wcpts = sg->wcpts();
        if (wcpts.size() <= 2) continue;
        
        for (size_t i = 1; i + 1 < wcpts.size(); i++) {
            pts.push_back(wcpts[i].point);
        }
    }
    
    // Collect points from vertices
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;
        pts.push_back(vtx->wcpt().point);
    }
    
    if (pts.size() <= 3) {
        return temp_main_vertex;
    }
    
    // Calculate PCA main axis
    auto pair_result = calc_PCA_main_axis(pts);
    Facade::geo_vector_t dir = pair_result.second;
    Facade::geo_point_t center = pair_result.first;
    
    // Find min and max vertices along the main axis
    double min_val = 1e9, max_val = -1e9;
    VertexPtr min_vtx = nullptr, max_vtx = nullptr;
    
    for (auto vtx : vertex_candidates) {
        double val = (vtx->wcpt().point.x() - center.x()) * dir.x() + 
                    (vtx->wcpt().point.y() - center.y()) * dir.y() + 
                    (vtx->wcpt().point.z() - center.z()) * dir.z();
        
        // Adjust for single short segment vertices
        if (vtx->descriptor_valid()) {
            auto vd = vtx->get_descriptor();
            auto edge_range = boost::out_edges(vd, graph);
            int num_segs = 0;
            double seg_length = 0;
            for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
                if (graph[*e_it].segment) {
                    num_segs++;
                    seg_length = segment_track_length(graph[*e_it].segment);
                }
            }
            
            if (num_segs == 1 && seg_length < 1 * units::cm) {
                if (val > 0) val -= 0.5 * units::cm;
                else if (val < 0) val += 0.5 * units::cm;
            }
        }
        
        if (val > max_val) {
            max_val = val;
            max_vtx = vtx;
        }
        if (val < min_val) {
            min_val = val;
            min_vtx = vtx;
        }
    }
    
    if (!min_vtx || !max_vtx || min_vtx == max_vtx) {
        return temp_main_vertex;
    }
    
    // Check if steiner point cloud exists
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    if (steiner_pc.size() < 3) {
        // Pick forward vertex based on z coordinate
        if (max_vtx->wcpt().point.z() < min_vtx->wcpt().point.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
        return temp_main_vertex;
    }
    
    // Find path between min and max vertices using steiner graph
    auto path_points = do_rough_path(cluster, max_vtx->wcpt().point, min_vtx->wcpt().point);
    
    if (path_points.size() <= 2) {
        // Pick forward vertex based on z coordinate
        if (max_vtx->wcpt().point.z() < min_vtx->wcpt().point.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
        return temp_main_vertex;
    }
    
    // Create temporary local graph for fitting
    auto local_graph = std::make_shared<PR::Graph>();
    
    // Create temporary vertices
    auto tmp_v1 = make_vertex(*local_graph);
    tmp_v1->wcpt().point = path_points.front();
    tmp_v1->cluster(&cluster);
    
    auto tmp_v2 = make_vertex(*local_graph);
    tmp_v2->wcpt().point = path_points.back();
    tmp_v2->cluster(&cluster);
    
    // Create temporary segment
    auto tmp_sg = create_segment_for_cluster(cluster, dv, path_points);
    if (!tmp_sg) {
        if (max_vtx->wcpt().point.z() < min_vtx->wcpt().point.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
        return temp_main_vertex;
    }
    
    // Add segment to local graph
    add_segment(*local_graph, tmp_sg, tmp_v1, tmp_v2);
    
    // Create local fitter with same configuration as input fitter
    TrackFitting local_fitter(TrackFitting::FittingType::Multiple);
    local_fitter.set_parameters(track_fitter.get_parameters());
    local_fitter.add_graph(local_graph);
    
    // Do fitting on local graph
    local_fitter.do_multi_tracking(true, true, false);
    
    // Create fit point cloud
    create_segment_fit_point_cloud(tmp_sg, dv, "fit");
    
    // Associate points from cluster to segment
    clustering_points_segments({tmp_sg}, dv, "associate_points", 0.5*units::cm, 3.0);
    
    // Determine shower direction
    segment_determine_shower_direction(tmp_sg, particle_data, recomb_model, "associate_points");
    
    double tmp_sg_length = segment_track_length(tmp_sg);
    int tmp_sg_dir = tmp_sg->dirsign();
    
    // Decide which vertex should be the main vertex based on direction
    if (tmp_sg_dir == 1) {
        temp_main_vertex = max_vtx;
    } else if (tmp_sg_dir == -1) {
        temp_main_vertex = min_vtx;
    } else {
        // No clear direction, pick forward vertex
        if (max_vtx->wcpt().point.z() < min_vtx->wcpt().point.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
    }
    
    // For large showers, always pick forward vertex
    if (tmp_sg_length > 80 * units::cm && 
        std::abs(max_vtx->wcpt().point.z() - min_vtx->wcpt().point.z()) > 40 * units::cm) {
        if (max_vtx->wcpt().point.z() < min_vtx->wcpt().point.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
    }
    
    // Local graph and temporary elements automatically cleaned up when going out of scope
    
    return temp_main_vertex;
}

float PatternAlgorithms::calc_conflict_maps(Graph& graph, VertexPtr vertex){
    // Assume temp_vertex is true neutrino vertex, calculate conflicts in the system
    float num_conflicts = 0;
    
    // Map segment to its direction (start vertex -> end vertex)
    std::map<SegmentPtr, std::pair<VertexPtr, VertexPtr>> map_seg_dir;
    std::set<VertexPtr> used_vertices;
    
    if (!vertex || !vertex->descriptor_valid()) return num_conflicts;
    
    // Start from the assumed neutrino vertex
    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;
    auto vd = vertex->get_descriptor();
    auto edge_range = boost::out_edges(vd, graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        SegmentPtr seg = graph[*e_it].segment;
        if (seg) {
            segments_to_be_examined.push_back(std::make_pair(vertex, seg));
        }
    }
    used_vertices.insert(vertex);
    
    // Propagate through the graph and build direction map
    while (!segments_to_be_examined.empty()) {
        std::vector<std::pair<VertexPtr, SegmentPtr>> temp_segments;
        
        for (const auto& [prev_vtx, current_sg] : segments_to_be_examined) {
            // Skip if already examined
            if (map_seg_dir.find(current_sg) != map_seg_dir.end()) continue;
            
            // Find the other vertex of this segment
            VertexPtr curr_vertex = find_other_vertex(graph, current_sg, prev_vtx);
            if (!curr_vertex) continue;
            
            // Record the direction: prev_vtx -> curr_vertex
            map_seg_dir[current_sg] = std::make_pair(prev_vtx, curr_vertex);
            
            // Skip if we've already processed this vertex
            if (used_vertices.find(curr_vertex) != used_vertices.end()) continue;
            
            // Add all segments connected to curr_vertex for examination
            if (curr_vertex->descriptor_valid()) {
                auto curr_vd = curr_vertex->get_descriptor();
                auto curr_edge_range = boost::out_edges(curr_vd, graph);
                for (auto e_it = curr_edge_range.first; e_it != curr_edge_range.second; ++e_it) {
                    SegmentPtr seg = graph[*e_it].segment;
                    if (seg) {
                        temp_segments.push_back(std::make_pair(curr_vertex, seg));
                    }
                }
            }
            used_vertices.insert(curr_vertex);
        }
        segments_to_be_examined = temp_segments;
    }
    
    // Check segments for direction conflicts
    for (const auto& [sg, vtx_pair] : map_seg_dir) {
        VertexPtr start_vtx = vtx_pair.first;
        
        // Determine which end of segment connects to start_vtx
        const auto& wcps = sg->wcpts();
        if (wcps.empty()) continue;
        
        bool flag_start = (ray_length(Ray{wcps.front().point, start_vtx->wcpt().point}) <
                          ray_length(Ray{wcps.back().point, start_vtx->wcpt().point}));
        
        // Check if segment has direction and is long enough to matter
        int dir_sign = sg->dirsign();
        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                        sg->flags_any(SegmentFlags::kShowerTopology);
        double sg_length = segment_track_length(sg);
        
        if (dir_sign != 0 && ((is_shower && sg_length > 5*units::cm) || !is_shower)) {
            // Check if direction conflicts with topology
            if ((flag_start && dir_sign == -1) || (!flag_start && dir_sign == 1)) {
                if (!sg->dir_weak()) {
                    num_conflicts += 1.0;
                } else {
                    num_conflicts += 0.5;
                }
            }
        }
    }
    
    // Beam direction (along z-axis)
    Facade::geo_vector_t dir_beam(0, 0, 1);
    
    // Check vertices for topology conflicts
    for (VertexPtr vtx : used_vertices) {
        if (!vtx->descriptor_valid()) continue;
        
        auto vtx_vd = vtx->get_descriptor();
        auto vtx_edge_range = boost::out_edges(vtx_vd, graph);
        
        // Count number of segments
        int num_segments = 0;
        for (auto e_it = vtx_edge_range.first; e_it != vtx_edge_range.second; ++e_it) {
            if (graph[*e_it].segment) num_segments++;
        }
        if (num_segments <= 1) continue;
        
        int n_in = 0;
        int n_in_shower = 0;
        int n_out_tracks = 0;
        int n_out_showers = 0;
        
        std::map<SegmentPtr, Facade::geo_vector_t> map_in_segment_dirs;
        std::map<SegmentPtr, Facade::geo_vector_t> map_out_segment_dirs;
        
        // Analyze each segment connected to this vertex
        for (auto e_it = vtx_edge_range.first; e_it != vtx_edge_range.second; ++e_it) {
            SegmentPtr sg = graph[*e_it].segment;
            if (!sg || map_seg_dir.find(sg) == map_seg_dir.end()) continue;
            
            VertexPtr start_vtx = map_seg_dir[sg].first;
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                            sg->flags_any(SegmentFlags::kShowerTopology);
            // bool is_shower_traj = sg->flags_any(SegmentFlags::kShowerTrajectory);
            
            if (vtx != start_vtx) {
                // Segment is incoming to this vertex
                n_in++;
                if (is_shower) n_in_shower++;
                map_in_segment_dirs[sg] = segment_cal_dir_3vector(sg, vtx->wcpt().point, 10*units::cm);
            } else {
                // Segment is outgoing from this vertex
                if (!is_shower) {
                    n_out_tracks++;
                } else {
                    n_out_showers++;
                }
                map_out_segment_dirs[sg] = segment_cal_dir_3vector(sg, vtx->wcpt().point, 10*units::cm);
            }
        }
        
        // Check angles between incoming and outgoing segments
        if (!map_in_segment_dirs.empty() && !map_out_segment_dirs.empty()) {
            double max_angle = -1;
            SegmentPtr sg1 = nullptr;
            SegmentPtr sg2 = nullptr;
            
            for (const auto& [in_sg, in_dir] : map_in_segment_dirs) {
                for (const auto& [out_sg, out_dir] : map_out_segment_dirs) {
                    double angle = std::acos(std::clamp(in_dir.dot(out_dir) / 
                                   (in_dir.magnitude() * out_dir.magnitude()), -1.0, 1.0)) 
                                   * 180.0 / M_PI;
                    if (angle > max_angle) {
                        max_angle = angle;
                        sg1 = in_sg;
                        sg2 = out_sg;
                    }
                }
            }
            
            if (sg1 && sg2) {
                bool flag_check = true;
                bool is_sg2_shower_traj = sg2->flags_any(SegmentFlags::kShowerTrajectory);
                bool is_sg1_shower = sg1->flags_any(SegmentFlags::kShowerTrajectory) || 
                                    sg1->flags_any(SegmentFlags::kShowerTopology);
                bool is_sg2_shower = sg2->flags_any(SegmentFlags::kShowerTrajectory) || 
                                    sg2->flags_any(SegmentFlags::kShowerTopology);
                
                // Skip check for shower trajectories or both showers
                if (is_sg2_shower_traj || (is_sg1_shower && is_sg2_shower)) {
                    flag_check = false;
                }
                
                double angle_beam = std::acos(std::clamp(map_in_segment_dirs[sg1].dot(dir_beam) / 
                                    map_in_segment_dirs[sg1].magnitude(), -1.0, 1.0)) 
                                    * 180.0 / M_PI;
                
                if (max_angle >= 0 && flag_check) {
                    if (max_angle < 35) {
                        num_conflicts += 5.0;
                    } else if (max_angle < 70) {
                        num_conflicts += 3.0;
                    } else if (max_angle < 85) {
                        num_conflicts += 1.0;
                    } else if (max_angle < 110) {
                        num_conflicts += 0.25;
                    }
                    
                    // Additional penalty for backward-going particles
                    if (angle_beam < 60 && max_angle < 110) {
                        num_conflicts += 1.0;
                    } else if (angle_beam < 45 && max_angle < 70) {
                        num_conflicts += 3.0;
                    }
                }
            }
        }
        
        // Penalize multiple incoming particles
        if (n_in > 1) {
            if (n_in != n_in_shower) {
                num_conflicts += (n_in - 1);
            } else {
                num_conflicts += (n_in - 1) / 2.0;
            }
        }
        
        // Penalize showers in with tracks out (suspicious topology)
        if (n_in_shower > 0 && n_out_tracks > 0) {
            num_conflicts += std::min(n_in_shower, n_out_tracks);
        }
        (void)n_out_showers; // to avoid unused variable warning
    }
    
    return num_conflicts;
}

VertexPtr PatternAlgorithms::compare_main_vertices(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates){
    if (vertex_candidates.empty()) return nullptr;
    
    std::map<VertexPtr, double> map_vertex_num;
    for (auto vtx : vertex_candidates) {
        map_vertex_num[vtx] = 0;
    }
    
    // Find the longest muon candidate
    SegmentPtr max_length_muon = nullptr;
    double max_length = 0;
    
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        
        // Skip showers
        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                        sg->flags_any(SegmentFlags::kShowerTopology);
        if (is_shower) continue;
        
        // Skip protons
        if (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 2212) continue;
        
        double length = segment_track_length(sg);
        if (length > max_length) {
            max_length = length;
            max_length_muon = sg;
        }
    }
    
    // Analyze proton topology for each vertex candidate
    for (auto vtx : vertex_candidates) {
        if (!vtx->descriptor_valid()) continue;
        
        int n_proton_in = 0;
        int n_proton_out = 0;
        
        auto vd = vtx->get_descriptor();
        auto edge_range = boost::out_edges(vd, graph);
        
        for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
            SegmentPtr sg = graph[*e_it].segment;
            if (!sg) continue;
            
            bool is_proton = sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 2212;
            if (!is_proton) continue;
            
            int dir_sign = sg->dirsign();
            bool is_weak = sg->dir_weak();
            
            if ((is_weak || dir_sign == 0)) {
                VertexPtr other_vertex = find_other_vertex(graph, sg, vtx);
                if (!other_vertex || !other_vertex->descriptor_valid()) continue;
                
                auto other_vd = other_vertex->get_descriptor();
                auto other_edge_range = boost::out_edges(other_vd, graph);
                
                int num_segs = 0;
                for (auto oe_it = other_edge_range.first; oe_it != other_edge_range.second; ++oe_it) {
                    if (graph[*oe_it].segment) num_segs++;
                }
                
                if (num_segs > 1) {
                    for (auto oe_it = other_edge_range.first; oe_it != other_edge_range.second; ++oe_it) {
                        SegmentPtr other_sg = graph[*oe_it].segment;
                        if (!other_sg) continue;
                        
                        const auto& wcps = other_sg->wcpts();
                        if (wcps.empty()) continue;
                        
                        bool flag_start = (ray_length(Ray{wcps.front().point, other_vertex->wcpt().point}) <
                                          ray_length(Ray{wcps.back().point, other_vertex->wcpt().point}));
                        
                        int other_dir = other_sg->dirsign();
                        bool other_weak = other_sg->dir_weak();
                        bool is_other_proton = other_sg->has_particle_info() && 
                                              std::abs(other_sg->particle_info()->pdg()) == 2212;
                        
                        if (!other_weak && is_other_proton) {
                            if ((flag_start && other_dir == 1) || (!flag_start && other_dir == -1)) {
                                n_proton_out++;
                            }
                            if ((flag_start && other_dir == -1) || (!flag_start && other_dir == 1)) {
                                n_proton_in++;
                            }
                        }
                        
                        if ((other_weak || other_dir == 0) && is_other_proton) {
                            n_proton_in++;
                        }
                    }
                }
            }
        }
        
        // Score proton topology
        if (n_proton_in > n_proton_out) {
            map_vertex_num[vtx] -= (n_proton_in - n_proton_out) / 4.0;
        } else {
            map_vertex_num[vtx] -= (n_proton_in - n_proton_out) / 4.0 - (n_proton_in + n_proton_out) / 8.0;
        }
    }
    
    // Score based on z position (prefer forward/upstream vertices)
    double min_z = 1e9;
    for (auto vtx : vertex_candidates) {
        if (vtx->wcpt().point.z() < min_z) min_z = vtx->wcpt().point.z();
    }
    
    for (auto vtx : vertex_candidates) {
        if (!vtx->descriptor_valid()) continue;
        
        // Position penalty
        map_vertex_num[vtx] -= (vtx->wcpt().point.z() - min_z) / (200 * units::cm);
        
        // Score based on connected segments
        auto vd = vtx->get_descriptor();
        auto edge_range = boost::out_edges(vd, graph);
        
        for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
            SegmentPtr sg = graph[*e_it].segment;
            if (!sg) continue;
            
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                            sg->flags_any(SegmentFlags::kShowerTopology);
            // bool is_shower_traj = sg->flags_any(SegmentFlags::kShowerTrajectory);
            
            if (is_shower) {
                map_vertex_num[vtx] += 1.0 / 4.0 / 2.0; // number of showers
                
                auto pair_results = calculate_num_daughter_showers(graph, vtx, sg);
                if (pair_results.second > 45 * units::cm) {
                    map_vertex_num[vtx] += 1.0 / 4.0 / 2.0;
                }
            } else {
                map_vertex_num[vtx] += 1.0 / 4.0; // number of tracks
            }
            
            int dir_sign = sg->dirsign();
            bool is_weak = sg->dir_weak();
            bool is_proton = sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 2212;
            
            if (is_proton && dir_sign != 0 && !is_weak) {
                map_vertex_num[vtx] += 1.0 / 4.0; // has a clear proton
            } else if (dir_sign != 0 && !is_shower) {
                map_vertex_num[vtx] += 1.0 / 4.0 / 2.0; // has direction with track
            }
            
            if (max_length > 35 * units::cm && sg == max_length_muon) {
                map_vertex_num[vtx] += 1.0 / 4.0 / 2.0; // long muon adds weight
            }
        }
    }
    
    // Score based on fiducial volume (removed offset_x in WCP, need to validate)
    auto fiducial_utils = cluster.grouping()->get_fiducialutils();
    if (fiducial_utils) {
        for (auto vtx : vertex_candidates) {
            if (fiducial_utils->inside_fiducial_volume(vtx->wcpt().point)) {
                map_vertex_num[vtx] += 0.5; // good - inside fiducial volume
            }
        }
    }
    
    // Score based on topology conflicts
    for (auto vtx : vertex_candidates) {
        double num_conflicts = calc_conflict_maps(graph, vtx);
        map_vertex_num[vtx] -= num_conflicts / 4.0;
    }
    
    // Find the vertex with maximum score
    double max_val = -1e9;
    VertexPtr max_vertex = nullptr;
    
    for (auto vtx : vertex_candidates) {
        if (map_vertex_num[vtx] > max_val) {
            max_val = map_vertex_num[vtx];
            max_vertex = vtx;
        }
    }
    
    return max_vertex;
}


std::pair<SegmentPtr, VertexPtr> PatternAlgorithms::find_cont_muon_segment(Graph &graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx){
    SegmentPtr sg1 = nullptr;
    VertexPtr vtx1 = nullptr;
    
    double max_length = 0;
    double max_angle = 0;
    double max_ratio = 0;
    
    bool flag_cont = false;
    
    double max_ratio1 = 0;
    double max_ratio1_length = 0;
    
    double sg_length = segment_track_length(sg);
    
    // Get vertex point
    WireCell::Point vtx_point = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
    
    if (!vtx->descriptor_valid()) {
        return std::make_pair(sg1, vtx1);
    }
    
    // Iterate through all segments connected to this vertex
    auto vd = vtx->get_descriptor();
    auto edge_range = boost::out_edges(vd, graph);
    
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        SegmentPtr sg2 = graph[*e_it].segment;
        if (!sg2 || sg2 == sg) continue;
        
        // Find the other vertex of sg2
        VertexPtr vtx2 = find_other_vertex(graph, sg2, vtx);
        if (!vtx2) continue;
        
        // Calculate direction vectors at 15cm from vertex
        Facade::geo_vector_t dir1 = segment_cal_dir_3vector(sg, vtx_point, 15*units::cm);
        Facade::geo_vector_t dir2 = segment_cal_dir_3vector(sg2, vtx_point, 15*units::cm);
        
        if (dir1.magnitude() == 0 || dir2.magnitude() == 0) continue;
        
        double length = segment_track_length(sg2);
        
        // Calculate angle (180° - angle between directions)
        double cos_angle = std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0);
        double angle = (M_PI - std::acos(cos_angle)) / M_PI * 180.0;
        
        // Calculate dQ/dx ratio
        double ratio = segment_median_dQ_dx(sg2) / (43e3 / units::cm);
        
        // For longer segments, also check angle at 50cm
        double angle1 = angle;
        if (length > 50*units::cm) {
            Facade::geo_vector_t dir3 = segment_cal_dir_3vector(sg, vtx_point, 50*units::cm);
            Facade::geo_vector_t dir4 = segment_cal_dir_3vector(sg2, vtx_point, 50*units::cm);
            
            if (dir3.magnitude() > 0 && dir4.magnitude() > 0) {
                double cos_angle1 = std::clamp(dir3.dot(dir4) / (dir3.magnitude() * dir4.magnitude()), -1.0, 1.0);
                angle1 = (M_PI - std::acos(cos_angle1)) / M_PI * 180.0;
            }
        }
        
        // Check if this segment qualifies as a continuation
        bool angle_ok = (angle < 10.0 || angle1 < 10.0 || 
                        (sg_length < 6*units::cm && (angle < 15.0 || angle1 < 15.0)));
        bool ratio_ok = (ratio < 1.3 || flag_ignore_dQ_dx);
        
        if (angle_ok && ratio_ok) {
            flag_cont = true;
            
            // Select segment with maximum projected length
            double projected_length = length * std::cos(angle / 180.0 * M_PI);
            if (projected_length > max_length) {
                max_length = projected_length;
                max_angle = angle;
                max_ratio = ratio;
                sg1 = sg2;
                vtx1 = vtx2;
            }
        } else {
            // Track maximum dQ/dx ratio among non-qualifying segments
            if (ratio > max_ratio1) {
                max_ratio1 = ratio;
                max_ratio1_length = length;
            }
        }
    }

    (void)max_angle;
    (void)max_ratio;
    (void)max_ratio1_length;
    
    if (flag_cont) {
        return std::make_pair(sg1, vtx1);
    } else {
        return std::make_pair(nullptr, nullptr);
    }
}
