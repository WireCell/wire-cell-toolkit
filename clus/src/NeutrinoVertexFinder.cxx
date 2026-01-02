#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

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
    bool has_direction = segment_determine_shower_direction(tmp_sg, particle_data, recomb_model, "associate_points");
    
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