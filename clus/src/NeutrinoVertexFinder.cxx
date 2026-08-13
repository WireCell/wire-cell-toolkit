#include "WireCellClus/NeutrinoPatternBase.h"
#include <algorithm>
#include <vector>
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/MyFCN.h"

#include "WireCellAux/Logger.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/BuildConfig.h"

#ifdef HAVE_PYTHON_INC
#include "WCPPyUtil/SCN_Vertex.h"
#endif

#include <chrono>
#include <cmath>
#include <unordered_set>
#include <cstdlib>
#include <cstdio>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// doc sbnd_xin/docs/pr/40: WCT_PID_WRITE_DEBUG (see PRSegment.cxx) logs
// direct set_pdg() mutations on the pointee, which bypass the
// Segment::particle_info() setter hook entirely.
static inline void pr40_probe_setpdg(SegmentPtr sg, int new_pdg, const char* site) {
    static const bool dbg = std::getenv("WCT_PID_WRITE_DEBUG") != nullptr;
    if (!dbg || !sg) return;
    if (std::abs(new_pdg) != 11) return;
    const int cid = sg->cluster() ? sg->cluster()->get_cluster_id() : -1;
    std::fprintf(stderr, "PID_WRITE_DEBUG set_pdg id=%d clus=%d gidx=%zu pdg -> %d  at %s\n",
                 sg->id(), cid, sg->get_graph_index(), new_pdg, site);
}

// doc sbnd_xin/docs/pr/74 round 2: improve_vertex's topology re-exam can
// clear a stage-3 kShowerTopology flag with no PID-trace fingerprint (the
// low-score branch below writes no pdg), which is why 90055's "second
// demotion" was unattributable in round 1.  Env-gated, default off.
static inline void pr74_probe_topo_reexam(SegmentPtr sg, const char* what) {
    static const bool dbg = std::getenv("WCT_SHOWER_TOPO_DEBUG") != nullptr;
    if (!dbg || !sg) return;
    const int cid = sg->cluster() ? sg->cluster()->get_cluster_id() : -1;
    std::fprintf(stderr, "TOPO_REEXAM id=%d clus=%d gidx=%zu %s pdg=%d score=%.3f\n",
                 sg->id(), cid, sg->get_graph_index(), what,
                 (sg->particle_info() ? sg->particle_info()->pdg() : 0),
                 sg->particle_score());
}

// doc sbnd_xin/docs/pr/32 §11 F1 -- the prototype's ProtoVertex::get_fit_pt().
//
// use_fit == false reproduces today's raw `wcpt()` read exactly, so every call
// site below is byte-identical with the knob off.  use_fit == true is the
// prototype's convention, and the `fit().valid()` fallback stands in for
// ProtoVertex's constructor initialising fit_pt = wcpt (ProtoVertex.cxx:17-19),
// which PR::Vertex does not do -- without it a never-fitted vertex would put
// (0,0,0) into a score ladder.
// The counter block is unconditional -- it runs in the knob-off arm too, which
// is what measures how often the divergence is reachable at all.  Incrementing
// an atomic cannot change reconstruction output.
static inline WireCell::Point pr32_vtx_pt(const VertexPtr& v, bool use_fit)
{
    g_pr32_audit.f1_reads.fetch_add(1, std::memory_order_relaxed);
    if (v->fit().valid()) {
        g_pr32_audit.f1_fit_valid.fetch_add(1, std::memory_order_relaxed);
        const uint64_t um = static_cast<uint64_t>(
            (v->fit().point - v->wcpt().point).magnitude() / WireCell::units::micrometer);
        g_pr32_audit.f1_moved_um_sum.fetch_add(um, std::memory_order_relaxed);
        uint64_t prev = g_pr32_audit.f1_moved_um_max.load(std::memory_order_relaxed);
        while (um > prev &&
               !g_pr32_audit.f1_moved_um_max.compare_exchange_weak(prev, um,
                                                                   std::memory_order_relaxed)) {}
        if (use_fit) return v->fit().point;
    }
    return v->wcpt().point;
}

// doc sbnd_xin/docs/pr/32 §11 -- per-event diagnostics for the four findings.
namespace WireCell::Clus::PR {
    Pr32AuditCounters g_pr32_audit;
}

// doc sbnd_xin/docs/pr/32 §11 F2(a) -- the improve_vertex shower-trajectory
// outer gate, at both of its sites.
//
// parity == false is today's path: recompute at (10 cm, m_mip_dqdx) and use
// that.  The recomputation is not free of side effects -- it can SET
// kShowerTrajectory on a segment whose stored label said otherwise -- which is
// itself part of the divergence.
//
// parity == true is the prototype's: read the stored label
// (NeutrinoID_improve_vertex.h:248 and :287) and do not touch it.  That is only
// meaningful once PR::g_shower_traj_refresh_flag lets the label be cleared;
// without it the stored bit is a monotone OR over history.
//
// The disagreement counter fills on the legacy path only, because that is the
// path where both answers are already in hand.  Asking for the fresh answer
// under parity would re-introduce the mutation the knob exists to remove.
static inline bool pr32_shower_traj_gate(const SegmentPtr& sg, bool parity, double mip_dqdx)
{
    g_pr32_audit.f2_gate_calls.fetch_add(1, std::memory_order_relaxed);
    const bool stored = sg->flags_any(SegmentFlags::kShowerTrajectory);
    if (parity) return stored;
    const bool fresh = segment_is_shower_trajectory(sg, 10*WireCell::units::cm, mip_dqdx);
    if (stored != fresh) {
        g_pr32_audit.f2_gate_disagree.fetch_add(1, std::memory_order_relaxed);
    }
    return fresh;
}

// Temporary determinism-debug helper (WCT_DET_DEBUG=1): content hash of the
// PR graph including wcpt/fit coordinates, to localize run-to-run divergence.
namespace {
    inline uint64_t detv_fnv(uint64_t h, const void* p, size_t n) {
        const unsigned char* b = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ULL; }
        return h;
    }
    void detv_dump(const char* tag, WireCell::Clus::PR::Graph& g) {
        static const bool on = (std::getenv("WCT_DET_DEBUG") != nullptr);
        if (!on) return;
        uint64_t h = 1469598103934665603ULL;
        size_t nvtx = 0, nseg = 0;
        for (auto nd : WireCell::Clus::PR::ordered_nodes(g)) {
            auto vtx = g[nd].vertex;
            if (!vtx) continue;
            ++nvtx;
            const auto& p = vtx->wcpt().point;
            double xyz[3] = {p.x(), p.y(), p.z()};
            h = detv_fnv(h, xyz, sizeof(xyz));
        }
        for (auto ed : WireCell::Clus::PR::ordered_edges(g)) {
            auto seg = g[ed].segment;
            if (!seg) continue;
            ++nseg;
            for (const auto& w : seg->wcpts()) {
                double xyz[3] = {w.point.x(), w.point.y(), w.point.z()};
                h = detv_fnv(h, xyz, sizeof(xyz));
            }
        }
        fprintf(stderr, "WCT_DETV %-28s nvtx=%zu nseg=%zu h=%016lx\n", tag, nvtx, nseg, (unsigned long)h);
    }
}

bool WireCell::Clus::PR::PatternAlgorithms::search_for_vertex_activities(Graph& graph, VertexPtr vertex, std::vector<SegmentPtr>& segments_set, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double search_range){
    s_log->trace("search_for_vertex_activities: cluster {} vertex ({:.2f}, {:.2f}, {:.2f}) search_range={:.2f} cm nsegs={}",
        cluster.ident(),
        vertex->wcpt().point.x()/units::cm, vertex->wcpt().point.y()/units::cm, vertex->wcpt().point.z()/units::cm,
        search_range/units::cm, segments_set.size());

    if (!cluster.has_pc("steiner_pc")) return false;

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

        for (auto e : ordered_edges(graph)) {
            SegmentPtr sg = graph[e].segment;
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
                    sum_charge += grouping->get_ave_charge(test_p_raw, test_wpid.apa(), test_wpid.face(), plane, 0.3*units::cm);
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
            
            for (auto e : ordered_edges(graph)) {
                SegmentPtr sg = graph[e].segment;
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
                        sum_charge += grouping->get_ave_charge(test_p_raw, test_wpid.apa(), test_wpid.face(), plane, 0.3*units::cm);
                        ncount++;
                    }
                }
                if (ncount != 0) sum_charge /= ncount;
                
                // prototype (NeutrinoID.h:1172): sum_charge > 20000 e = 20/43
                // of the MIP-median dQ/dx integrated over 1 cm.  Expressed on
                // m_mip_dqdx_median so it scales with the detector; exactly
                // 20000 at the uBooNE default (43000*20 = 860000, /43 = 20000,
                // FP-exact).  docs/pr/2 sec 7.3.
                if (min_dis + (min_dis_u + min_dis_v + min_dis_w) / std::sqrt(3.0) > max_dis &&
                    sum_charge > m_mip_dqdx_median * units::cm * 20.0 / 43.0) {
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
            s_log->trace("search_for_vertex_activities: cluster {} vertex activity found at ({:.2f}, {:.2f}, {:.2f}), max_dis={:.4f}",
                cluster.ident(), max_point.x()/units::cm, max_point.y()/units::cm, max_point.z()/units::cm, max_dis);

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
        s_log->trace("examine_main_vertex_candidate: invalid vertex, returning early");
        return std::make_tuple(flag_in, ntracks, nshowers);
    }
    s_log->trace("examine_main_vertex_candidate: vertex ({:.2f}, {:.2f}, {:.2f})",
        vertex->wcpt().point.x()/units::cm, vertex->wcpt().point.y()/units::cm, vertex->wcpt().point.z()/units::cm);
    
    auto vd = vertex->get_descriptor();

    // Stable edge-index order: the early 'break' below leaves ntracks/nshowers
    // as partial counts, so the iteration order escapes into the result.
    for (auto edesc : sorted_out_edges(vd, graph)) {
        SegmentPtr sg = graph[edesc].segment;
        if (!sg) continue;
        
        // matches prototype get_flag_shower() = kShowerTrajectory || kShowerTopology || particle_type==11
        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                        sg->flags_any(SegmentFlags::kShowerTopology) ||
                        (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
        
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
        bool is_dir_weak = seg_dir_weak(sg);
        
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
    const auto edge_range2 = sorted_out_edges(vd, graph);
    for (auto eit : edge_range2) {
        if (graph[eit].segment) num_segments++;
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
    s_log->trace("compare_main_vertices_all_showers: cluster {} ncandidates={}", cluster.ident(), vertex_candidates.size());
    if (vertex_candidates.empty()) return nullptr;

    VertexPtr temp_main_vertex = vertex_candidates.front();
    
    // Collect all points from segments and vertices in the cluster
    std::vector<Facade::geo_point_t> pts;
    
    // Collect points from segments.  Stable edge/node-index order: pts feeds
    // FP accumulation downstream, so pointer-order iteration would leak.
    for (auto ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;

        const auto& wcpts = sg->wcpts();
        if (wcpts.size() <= 2) continue;

        for (size_t i = 1; i + 1 < wcpts.size(); i++) {
            pts.push_back(wcpts[i].point);
        }
    }

    // Collect points from vertices.
    // doc pr/32 §11 F1: prototype :1468 pushes get_fit_pt().
    for (auto nd : ordered_nodes(graph)) {
        VertexPtr vtx = graph[nd].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;
        pts.push_back(pr32_vtx_pt(vtx, m_vertex_dir_use_fit_point));
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
        // doc pr/32 §11 F1: prototype :1480 projects get_fit_pt().
        const WireCell::Point vpt = pr32_vtx_pt(vtx, m_vertex_dir_use_fit_point);
        double val = (vpt.x() - center.x()) * dir.x() +
                    (vpt.y() - center.y()) * dir.y() +
                    (vpt.z() - center.z()) * dir.z();
        
        // Adjust for single short segment vertices
        if (vtx->descriptor_valid()) {
            auto vd = vtx->get_descriptor();
            const auto edge_range = sorted_out_edges(vd, graph);
            int num_segs = 0;
            double seg_length = 0;
            for (auto e_it : edge_range) {
                if (graph[e_it].segment) {
                    num_segs++;
                    seg_length = segment_track_length(graph[e_it].segment);
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

    // doc pr/32 §11 F1: every remaining use of these two vertices' positions --
    // the four forward-z tie-breaks (prototype :1542, :1551, :1567, :1574) and
    // the Steiner path endpoints (prototype :1504-1505) -- reads get_fit_pt()
    // on the prototype side.  do_rough_path snaps to the nearest steiner_pc
    // point internally (NeutrinoPatternBase.cxx:97-104), which is exactly the
    // prototype's get_closest_wcpoint(get_fit_pt()), so the fit point is the
    // faithful argument there too.  Named locals: do_rough_path takes
    // non-const Point&.
    WireCell::Point max_pt = pr32_vtx_pt(max_vtx, m_vertex_dir_use_fit_point);
    WireCell::Point min_pt = pr32_vtx_pt(min_vtx, m_vertex_dir_use_fit_point);

    // Check if steiner point cloud exists and has enough points
    if (!cluster.has_pc("steiner_pc")) {
        // Fall through to backup vertex selection below
        // Pick forward vertex based on z coordinate
        if (max_pt.z() < min_pt.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
        return temp_main_vertex;
    }

    // Find path between min and max vertices using steiner graph
    auto path_points = do_rough_path(cluster, max_pt, min_pt);

    if (path_points.size() <= 2) {
        // Pick forward vertex based on z coordinate
        if (max_pt.z() < min_pt.z()) {
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
        if (max_pt.z() < min_pt.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
        return temp_main_vertex;
    }
    
    // Add segment to local graph
    add_segment(*local_graph, tmp_sg, tmp_v1, tmp_v2);
    
    // Create local fitter that inherits pre-built geometry and cluster charge data
    // from the parent fitter, avoiding expensive BuildGeometry() + prepare_data() calls.
    TrackFitting local_fitter(TrackFitting::FittingType::Multiple);
    local_fitter.inherit_from(track_fitter, &cluster);
    local_fitter.add_graph(local_graph);
    
    // Do fitting on local graph
    // 3-argument form: flag_exclusion takes its default (false), and that is
    // correct here regardless of m_fit_exclusion (doc pr/30 §11 P1).
    // local_graph holds exactly ONE segment, so update_association -- which
    // strips the 2-D cells belonging to OTHER segments -- has nothing to
    // strip.  Knobbing this site would be a no-op that reads as a divergence.
    local_fitter.do_multi_tracking(true, true, false);
    
    // Create fit point cloud
    create_segment_fit_point_cloud(tmp_sg, dv, "fit");
    
    // Associate points from cluster to segment
    clustering_points_segments({tmp_sg}, dv, "associate_points", 0.5*units::cm, 3.0);
    
    // Determine shower direction
    segment_determine_shower_direction(tmp_sg, particle_data, recomb_model, "associate_points", m_mip_dqdx_median, 0.4*units::cm, m_mip_dqdx, m_dir_track_median_local);
    
    double tmp_sg_length = segment_track_length(tmp_sg);
    int tmp_sg_dir = tmp_sg->dirsign();
    
    // Decide which vertex should be the main vertex based on direction
    if (tmp_sg_dir == 1) {
        temp_main_vertex = max_vtx;
    } else if (tmp_sg_dir == -1) {
        temp_main_vertex = min_vtx;
    } else {
        // No clear direction, pick forward vertex
        if (max_pt.z() < min_pt.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
    }

    // For large showers, always pick forward vertex
    if (tmp_sg_length > 80 * units::cm &&
        std::abs(max_pt.z() - min_pt.z()) > 40 * units::cm) {
        if (max_pt.z() < min_pt.z()) {
            temp_main_vertex = max_vtx;
        } else {
            temp_main_vertex = min_vtx;
        }
    }
    
    // Local graph and temporary elements automatically cleaned up when going out of scope
    s_log->trace("compare_main_vertices_all_showers: selected vertex ({:.2f}, {:.2f}, {:.2f}) sg_length={:.2f}cm sg_dir={}",
        temp_main_vertex->wcpt().point.x()/units::cm, temp_main_vertex->wcpt().point.y()/units::cm,
        temp_main_vertex->wcpt().point.z()/units::cm, tmp_sg_length/units::cm, tmp_sg_dir);

    return temp_main_vertex;
}

float PatternAlgorithms::calc_conflict_maps(Graph& graph, VertexPtr vertex){
    // Assume temp_vertex is true neutrino vertex, calculate conflicts in the system
    float num_conflicts = 0;
    
    // Map segment to its direction (start vertex -> end vertex)
    std::map<SegmentPtr, std::pair<VertexPtr, VertexPtr>> map_seg_dir;
    std::set<VertexPtr> used_vertices;
    
    if (!vertex || !vertex->descriptor_valid()) return num_conflicts;
    
    // Start from the assumed neutrino vertex.  Stable edge-index order: on a
    // graph with cycles the BFS first-writer wins in map_seg_dir, so frontier
    // order affects the recorded directions (and the conflict score).
    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;
    auto vd = vertex->get_descriptor();
    for (auto edesc : sorted_out_edges(vd, graph)) {
        SegmentPtr seg = graph[edesc].segment;
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
            // (stable edge-index order — see seed loop above).
            if (curr_vertex->descriptor_valid()) {
                auto curr_vd = curr_vertex->get_descriptor();
                for (auto edesc : sorted_out_edges(curr_vd, graph)) {
                    SegmentPtr seg = graph[edesc].segment;
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
                        sg->flags_any(SegmentFlags::kShowerTopology) ||
                        (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
        double sg_length = segment_track_length(sg);
        
        if (dir_sign != 0 && ((is_shower && sg_length > 5*units::cm) || !is_shower)) {
            // Check if direction conflicts with topology
            if ((flag_start && dir_sign == -1) || (!flag_start && dir_sign == 1)) {
                if (!seg_dir_weak(sg)) {
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
        const auto vtx_edge_range = sorted_out_edges(vtx_vd, graph);
        
        // Count number of segments
        int num_segments = 0;
        for (auto e_it : vtx_edge_range) {
            if (graph[e_it].segment) num_segments++;
        }
        if (num_segments <= 1) continue;
        
        int n_in = 0;
        int n_in_shower = 0;
        int n_out_tracks = 0;
        int n_out_showers = 0;
        
        // Index-ordered, not address-ordered: the (in,out) search below keeps the
        // FIRST pair reaching max_angle (strict `>` at :697), so an exact tie --
        // two outgoing legs at the same angle to one incoming leg -- was resolved
        // by heap layout.  Both maps are insert-and-iterate only; the operator[]
        // at :719 uses a key taken from this same map, so switching the key from
        // address identity to graph-index identity does not change what it finds
        // (every key here came from a live graph edge, so its index is set and
        // unique -- never the SIZE_MAX sentinel, doc pr/28 sec 14.7).
        std::map<SegmentPtr, Facade::geo_vector_t, PR::SegmentIndexCmp> map_in_segment_dirs;
        std::map<SegmentPtr, Facade::geo_vector_t, PR::SegmentIndexCmp> map_out_segment_dirs;
        
        // doc pr/32 §11 F1: the point BOTH direction vectors are measured from.
        // segment_cal_dir_3vector takes a non-const Point&, so this must be a
        // named local; it is the same vertex for the whole loop.
        WireCell::Point vtx_dir_pt = pr32_vtx_pt(vtx, m_vertex_dir_use_fit_point);

        // Analyze each segment connected to this vertex
        for (auto e_it : vtx_edge_range) {
            SegmentPtr sg = graph[e_it].segment;
            if (!sg || map_seg_dir.find(sg) == map_seg_dir.end()) continue;

            VertexPtr start_vtx = map_seg_dir[sg].first;
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                            sg->flags_any(SegmentFlags::kShowerTopology) ||
                            (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
            // bool is_shower_traj = sg->flags_any(SegmentFlags::kShowerTrajectory);

            if (vtx != start_vtx) {
                // Segment is incoming to this vertex
                n_in++;
                if (is_shower) n_in_shower++;
                // doc pr/32 §11 F1: prototype NeutrinoID_track_shower.h:1804
                // measures from get_fit_pt(), not the Steiner-snapped wcpt.
                map_in_segment_dirs[sg] = segment_cal_dir_3vector(sg, vtx_dir_pt, 10*units::cm);
            } else {
                // Segment is outgoing from this vertex
                if (!is_shower) {
                    n_out_tracks++;
                } else {
                    n_out_showers++;
                }
                // doc pr/32 §11 F1: prototype :1808, same substitution.
                map_out_segment_dirs[sg] = segment_cal_dir_3vector(sg, vtx_dir_pt, 10*units::cm);
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
                                    sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                    (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11);
                bool is_sg2_shower = sg2->flags_any(SegmentFlags::kShowerTrajectory) ||
                                    sg2->flags_any(SegmentFlags::kShowerTopology) ||
                                    (sg2->has_particle_info() && std::abs(sg2->particle_info()->pdg()) == 11);
                
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
    s_log->trace("compare_main_vertices: cluster {} ncandidates={}", cluster.ident(), vertex_candidates.size());
    if (vertex_candidates.empty()) return nullptr;

    // doc pr/32 §11 F3 (was P7): two of the six scoring blocks guard on
    // descriptor_valid() and four do not -- the min_z scan, the fiducial +0.5,
    // the conflict penalty and the argmax.  An invalid-descriptor candidate
    // therefore reaches the argmax carrying `0 + (0.5 if in FV) - conflicts/4`,
    // and real candidates routinely score negative, so it can win.  Filtering
    // once here makes the whole scorer see one candidate set instead of
    // patching four more guards in.  The caller's vector is NOT mutated.
    //
    // The drop is COUNTED unconditionally: if it is 0 on every event then P7 is
    // confirmed dead code and this knob can be retired, which is a measurement
    // rather than the control-flow argument doc pr/32 §10.4 had to settle for.
    std::vector<VertexPtr> filtered;
    for (auto vtx : vertex_candidates) {
        if (vtx->descriptor_valid()) filtered.push_back(vtx);
        else g_pr32_audit.f3_dropped.fetch_add(1, std::memory_order_relaxed);
    }
    g_pr32_audit.f3_candidates.fetch_add(vertex_candidates.size(), std::memory_order_relaxed);
    std::vector<VertexPtr>& scored =
        (m_main_vertex_require_descriptor ? filtered : vertex_candidates);
    if (scored.empty()) return nullptr;

    std::map<VertexPtr, double> map_vertex_num;
    for (auto vtx : scored) {
        map_vertex_num[vtx] = 0;
    }
    
    // Find the longest muon candidate
    SegmentPtr max_length_muon = nullptr;
    double max_length = 0;
    
    // ordered_edges: max_length_muon is a tie-broken max over all segments.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;

        // Skip showers
        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                        sg->flags_any(SegmentFlags::kShowerTopology) ||
                        (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
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
    for (auto vtx : scored) {
        if (!vtx->descriptor_valid()) continue;
        
        int n_proton_in = 0;
        int n_proton_out = 0;
        
        auto vd = vtx->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);
        
        for (auto e_it : edge_range) {
            SegmentPtr sg = graph[e_it].segment;
            if (!sg) continue;
            
            bool is_proton = sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 2212;
            if (!is_proton) continue;
            
            int dir_sign = sg->dirsign();
            bool is_weak = seg_dir_weak(sg);
            
            if ((is_weak || dir_sign == 0)) {
                VertexPtr other_vertex = find_other_vertex(graph, sg, vtx);
                if (!other_vertex || !other_vertex->descriptor_valid()) continue;
                
                auto other_vd = other_vertex->get_descriptor();
                const auto other_edge_range = sorted_out_edges(other_vd, graph);
                
                int num_segs = 0;
                for (auto oe_it : other_edge_range) {
                    if (graph[oe_it].segment) num_segs++;
                }
                
                if (num_segs > 1) {
                    for (auto oe_it : other_edge_range) {
                        SegmentPtr other_sg = graph[oe_it].segment;
                        if (!other_sg) continue;
                        
                        const auto& wcps = other_sg->wcpts();
                        if (wcps.empty()) continue;
                        
                        bool flag_start = (ray_length(Ray{wcps.front().point, other_vertex->wcpt().point}) <
                                          ray_length(Ray{wcps.back().point, other_vertex->wcpt().point}));
                        
                        int other_dir = other_sg->dirsign();
                        bool other_weak = seg_dir_weak(other_sg);
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
    // Use fit position (improve_vertex has already run in this branch)
    auto vtx_pt = [](VertexPtr v) -> WireCell::Point {
        return v->fit().valid() ? v->fit().point : v->wcpt().point;
    };
    double min_z = 1e9;
    for (auto vtx : scored) {
        if (vtx_pt(vtx).z() < min_z) min_z = vtx_pt(vtx).z();
    }

    for (auto vtx : scored) {
        if (!vtx->descriptor_valid()) continue;

        // Position penalty (uBooNE 200 cm; m_vertex_z_prior_scale, doc pr/2 2e(iv))
        map_vertex_num[vtx] -= (vtx_pt(vtx).z() - min_z) / m_vertex_z_prior_scale;
        
        // Score based on connected segments
        auto vd = vtx->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);
        
        for (auto e_it : edge_range) {
            SegmentPtr sg = graph[e_it].segment;
            if (!sg) continue;
            
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                            sg->flags_any(SegmentFlags::kShowerTopology) ||
                            (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
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
            bool is_weak = seg_dir_weak(sg);
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
    
    // Score based on fiducial volume — use fit position (improve_vertex has already run)
    auto fiducial_utils = cluster.grouping()->get_fiducialutils();
    if (fiducial_utils) {
        for (auto vtx : scored) {
            if (fiducial_utils->inside_fiducial_volume(vtx_pt(vtx))) {
                map_vertex_num[vtx] += 0.5; // good - inside fiducial volume
            }
        }
    }
    
    // Score based on topology conflicts
    for (auto vtx : scored) {
        double num_conflicts = calc_conflict_maps(graph, vtx);
        map_vertex_num[vtx] -= num_conflicts / 4.0;
    }
    
    // Find the vertex with maximum score
    double max_val = -1e9;
    VertexPtr max_vertex = nullptr;
    
    for (auto vtx : scored) {
        if (map_vertex_num[vtx] > max_val) {
            max_val = map_vertex_num[vtx];
            max_vertex = vtx;
        }
    }

    if (max_vertex) {
        s_log->trace("compare_main_vertices: selected vertex ({:.2f}, {:.2f}, {:.2f}) score={:.3f}",
            max_vertex->wcpt().point.x()/units::cm, max_vertex->wcpt().point.y()/units::cm,
            max_vertex->wcpt().point.z()/units::cm, max_val);
    } else {
        s_log->trace("compare_main_vertices: no vertex selected");
    }
    return max_vertex;
}


std::pair<SegmentPtr, VertexPtr> PatternAlgorithms::find_cont_muon_segment(Graph &graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx, bool flag_stub_bridge){
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
    const auto edge_range = sorted_out_edges(vd, graph);
    
    for (auto e_it : edge_range) {
        SegmentPtr sg2 = graph[e_it].segment;
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
        double ratio = segment_median_dQ_dx(sg2) / m_mip_dqdx_median;
        
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

        // doc sbnd_xin/docs/pr/46 stub bridge: a short (< 6 cm) incoming
        // segment carries no reliable direction, so a broken long muon's
        // first piece fails the 15 deg test against its own continuation
        // (evt 55595: 30.5-35.4 deg vs the 192.9 cm MIP muon).  Accept a
        // long (> 35 cm) MIP candidate up to 45 deg, unless the junction
        // has another substantial track-like arm (> 10 cm, not
        // shower-flagged, pdg != +-11): multiple substantial outgoing
        // tracks mean a genuine hadronic vertex (evt 66118 at 70-82 deg is
        // additionally excluded by the angle cap).  Delta-ray electrons and
        // tiny fragments do not veto.  See m_long_muon_stub_bridge.
        if (flag_stub_bridge && !angle_ok && ratio_ok &&
            sg_length < 6*units::cm && length > 35*units::cm &&
            (angle < 45.0 || angle1 < 45.0)) {
            bool has_other_track = false;
            for (auto e2_it : edge_range) {
                SegmentPtr sg3 = graph[e2_it].segment;
                if (!sg3 || sg3 == sg || sg3 == sg2) continue;
                if (sg3->flags_any(SegmentFlags::kShowerTrajectory) ||
                    sg3->flags_any(SegmentFlags::kShowerTopology)) continue;
                if (sg3->has_particle_info() &&
                    std::abs(sg3->particle_info()->pdg()) == 11) continue;
                if (segment_track_length(sg3) > 10*units::cm) {
                    has_other_track = true;
                    break;
                }
            }
            if (!has_other_track) angle_ok = true;
        }

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

bool PatternAlgorithms::examine_direction(Graph& graph, VertexPtr vertex, VertexPtr main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_final){
    if (!vertex || !vertex->cluster()) return false;

    Facade::Cluster& cluster = *vertex->cluster();
    s_log->trace("examine_direction: cluster {} vtx ({:.2f},{:.2f},{:.2f}) is_main_vtx={} flag_final={}",
        cluster.ident(),
        vertex->wcpt().point.x()/units::cm, vertex->wcpt().point.y()/units::cm, vertex->wcpt().point.z()/units::cm,
        (vertex == main_vertex), flag_final);
    
    // Calculate cluster statistics
    double max_vtx_length = 0;
    double min_vtx_length = 1e9;
    int num_total_segments = 0;
    bool flag_only_showers = true;
    
    // Examine vertex segments to determine characteristics
    if (vertex->descriptor_valid()) {
        auto vd = vertex->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);
        for (auto e_it : edge_range) {
            SegmentPtr seg = graph[e_it].segment;
            if (!seg) continue;
            
            double length = segment_track_length(seg);
            if (length > max_vtx_length) max_vtx_length = length;
            if (length < min_vtx_length) min_vtx_length = length;
        }
    }
    
    // Check all vertices in the cluster, in stable node-index order.
    for (const auto& vd_it : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd_it].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;

        auto results = examine_main_vertex_candidate(graph, vtx);
        bool flag_in = std::get<0>(results);
        int ntracks = std::get<1>(results);
        // int nshowers = std::get<2>(results);
        
        if (!flag_in && ntracks > 0) {
            flag_only_showers = false;
        }
    }
    
    // Count total segments in cluster
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr seg = graph[ed].segment;
        if (seg && seg->cluster() == &cluster) {
            num_total_segments++;
        }
    }
    
    // Determine if only showers based on topology
    if (vertex->descriptor_valid()) {
        auto vd = vertex->get_descriptor();
        int num_vertex_segments = 0;
        const auto edge_range = sorted_out_edges(vd, graph);
        for (auto e_it : edge_range) {
            if (graph[e_it].segment) num_vertex_segments++;
        }
        
        if ((num_vertex_segments == 2 && (max_vtx_length > 30*units::cm || min_vtx_length > 15*units::cm)) ||
            (num_vertex_segments > 2 && num_total_segments > 4) ||
            (num_vertex_segments > 3)) {
            flag_only_showers = false;
        }
    }
    
    s_log->trace("examine_direction: cluster {} topology: num_total_segs={} max_vtx_len={:.2f}cm min_vtx_len={:.2f}cm flag_only_showers={}",
        cluster.ident(), num_total_segments,
        max_vtx_length/units::cm, min_vtx_length/units::cm, flag_only_showers);

    // Beam direction (along z-axis)
    Facade::geo_vector_t drift_dir(1, 0, 0);
    
    // Track used vertices and segments
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    // Start propagation from the main vertex.  Stable edge-index order: the
    // BFS assigns each segment's direction/particle type on first discovery,
    // so frontier order is a physics decision on any multi-path graph.
    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;
    if (vertex->descriptor_valid()) {
        auto vd = vertex->get_descriptor();
        for (auto edesc : sorted_out_edges(vd, graph)) {
            SegmentPtr seg = graph[edesc].segment;
            if (seg) {
                segments_to_be_examined.push_back(std::make_pair(vertex, seg));
            }
        }
    }
    used_vertices.insert(vertex);
    
    // Propagate through the graph setting directions and particle types
    while (!segments_to_be_examined.empty()) {
        std::vector<std::pair<VertexPtr, SegmentPtr>> temp_segments;
        
        for (const auto& [prev_vtx, current_sg] : segments_to_be_examined) {
            if (!prev_vtx->descriptor_valid()) continue;
            
            // std::cout << "examine segment " << current_sg->get_graph_index() << std::endl;

            // Check for incoming showers
            bool flag_shower_in = false;
            std::vector<SegmentPtr> in_showers;
            
            auto prev_vd = prev_vtx->get_descriptor();
            const auto prev_edge_range = sorted_out_edges(prev_vd, graph);
            for (auto e_it : prev_edge_range) {
                SegmentPtr sg = graph[e_it].segment;
                if (!sg || sg == current_sg) continue;
                // Only count a segment as "incoming shower" if it was already processed
                // in an earlier BFS level. Sibling segments in the same BFS level may
                // not have their direction set yet, causing order-dependent results.
                if (used_segments.find(sg) == used_segments.end()) continue;

                const auto& wcps = sg->wcpts();
                if (wcps.empty()) continue;

                bool flag_start = (ray_length(Ray{wcps.front().point, prev_vtx->wcpt().point}) <
                                  ray_length(Ray{wcps.back().point, prev_vtx->wcpt().point}));

                int dir_sign = sg->dirsign();
                if ((flag_start && dir_sign == -1) || (!flag_start && dir_sign == 1)) {
                    if (sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                        sg->flags_any(SegmentFlags::kShowerTopology) ||
                        (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11)) {
                        flag_shower_in = true;
                        in_showers.push_back(sg);
                        break;
                    }
                }
            }
            
            if (used_segments.find(current_sg) != used_segments.end()) continue;

            double length = segment_track_length(current_sg);
            bool is_shower = current_sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                            current_sg->flags_any(SegmentFlags::kShowerTopology) ||
                            (current_sg->has_particle_info() && std::abs(current_sg->particle_info()->pdg()) == 11);

            s_log->trace("examine_direction: cluster {} processing seg: length={:.2f}cm is_shower={} flag_shower_in={} cur_dirsign={} dir_weak={}",
                cluster.ident(), length/units::cm, is_shower, flag_shower_in,
                current_sg->dirsign(), seg_dir_weak(current_sg));

            // Determine segment direction
            if (current_sg->dirsign() == 0 || seg_dir_weak(current_sg) || is_shower || flag_final) {
                const auto& wcps = current_sg->wcpts();
                s_log->trace("examine_direction: processing seg nfits={} nwcpts={} dirsign_before={}"
                    " is_shower={} shower_topo={} shower_traj={}",
                    current_sg->fits().size(), wcps.size(), current_sg->dirsign(),
                    is_shower ? 1 : 0,
                    current_sg->flags_any(SegmentFlags::kShowerTopology) ? 1 : 0,
                    current_sg->flags_any(SegmentFlags::kShowerTrajectory) ? 1 : 0);
                if (!wcps.empty()) {
                    bool flag_start = (ray_length(Ray{wcps.front().point, prev_vtx->wcpt().point}) <
                                      ray_length(Ray{wcps.back().point, prev_vtx->wcpt().point}));

                    // Set direction
                    if (flag_start) {
                        current_sg->dirsign(1);
                    } else {
                        current_sg->dirsign(-1);
                    }
                    s_log->trace("examine_direction:   → dirsign set to {}", current_sg->dirsign());
                    
                    // Determine particle type
                    if (flag_shower_in && current_sg->dirsign() == 0 && !is_shower) {
                        auto four_momentum = segment_cal_4mom(current_sg, 11, particle_data, recomb_model, m_mip_dqdx);
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11), four_momentum);
                        current_sg->particle_info(pinfo);
                    } else if (flag_shower_in && length < 2.0*units::cm && !is_shower) {
                        auto four_momentum = segment_cal_4mom(current_sg, 11, particle_data, recomb_model, m_mip_dqdx);
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11), four_momentum);
                        current_sg->particle_info(pinfo);
                    } else if (flag_shower_in) {
                        // Matches prototype: flag_shower_in && (|pdg|==13 || pdg==0).
                        // no-particle-info means pdg defaults to 0, which also qualifies.
                        int cur_pdg = current_sg->has_particle_info() ? current_sg->particle_info()->pdg() : 0;
                        if (std::abs(cur_pdg) == 13 || cur_pdg == 0) {
                            // doc sbnd_xin/docs/pr/74 round 2 P1 -- see
                            // m_shower_in_cascade_guard's docstring.  Knob off
                            // => the veto is never evaluated => byte-identical.
                            if (m_shower_in_cascade_guard &&
                                segment_shower_in_cascade_vetoed(current_sg, m_mip_dqdx_median,
                                                                 m_shower_in_max_len, m_shower_in_mip_hi)) {
                                SPDLOG_LOGGER_DEBUG(s_log,
                                    "pr74 shower_in_cascade_guard: veto e- relabel gidx={} L={:.1f}cm pdg={}",
                                    current_sg->get_graph_index(), length/units::cm, cur_pdg);
                            } else {
                                auto four_momentum = segment_cal_4mom(current_sg, 11, particle_data, recomb_model, m_mip_dqdx);
                                auto pinfo = std::make_shared<Aux::ParticleInfo>(11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11), four_momentum);
                                current_sg->particle_info(pinfo);
                            }
                        }
                    } else {
                        auto pair_result = calculate_num_daughter_showers(graph, prev_vtx, current_sg);
                        auto pair_result1 = calculate_num_daughter_showers(graph, prev_vtx, current_sg, false);
                        int num_daughter_showers = pair_result.first;
                        double length_daughter_showers = pair_result.second;
                        
                        // Check if should be electron based on daughter showers
                        int current_pdg = current_sg->has_particle_info() ? current_sg->particle_info()->pdg() : 0;
                        if (current_pdg != 11 && 
                            (num_daughter_showers >= 4 || 
                             (length_daughter_showers > 50*units::cm && num_daughter_showers >= 2)) &&
                            pair_result.second > pair_result1.second - length - pair_result.second) {
                            
                            // Check angles with connected segments
                            bool flag_change = false;
                            VertexPtr next_vertex = find_other_vertex(graph, current_sg, prev_vtx);
                            if (next_vertex && next_vertex->descriptor_valid()) {
                                Facade::geo_vector_t tmp_dir1 = segment_cal_dir_3vector(current_sg, next_vertex->wcpt().point, 15*units::cm);
                                
                                auto next_vd = next_vertex->get_descriptor();
                                const auto next_edge_range = sorted_out_edges(next_vd, graph);
                                for (auto ne_it : next_edge_range) {
                                    SegmentPtr other_sg = graph[ne_it].segment;
                                    if (!other_sg || other_sg == current_sg) continue;
                                    
                                    Facade::geo_vector_t tmp_dir2 = segment_cal_dir_3vector(other_sg, next_vertex->wcpt().point, 15*units::cm);
                                    
                                    if (tmp_dir1.magnitude() > 0 && tmp_dir2.magnitude() > 0) {
                                        double angle = std::acos(std::clamp(tmp_dir1.dot(tmp_dir2) / 
                                                      (tmp_dir1.magnitude() * tmp_dir2.magnitude()), -1.0, 1.0)) 
                                                      * 180.0 / M_PI;
                                        double angle_drift1 = std::acos(std::clamp(drift_dir.dot(tmp_dir1) / 
                                                             (drift_dir.magnitude() * tmp_dir1.magnitude()), -1.0, 1.0)) 
                                                             * 180.0 / M_PI;
                                        double angle_drift2 = std::acos(std::clamp(drift_dir.dot(tmp_dir2) / 
                                                             (drift_dir.magnitude() * tmp_dir2.magnitude()), -1.0, 1.0)) 
                                                             * 180.0 / M_PI;
                                        double other_length = segment_track_length(other_sg);
                                        
                                        if (angle > 155 ||
                                            (angle > 135 && std::abs(angle_drift1 - 90) < 10 && std::abs(angle_drift2 - 90) < 10) ||
                                            (angle > 135 && other_length < 6*units::cm)) {
                                            flag_change = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            
                            if (flag_change) {
                                auto four_momentum = segment_cal_4mom(current_sg, 11, particle_data, recomb_model, m_mip_dqdx);
                                auto pinfo = std::make_shared<Aux::ParticleInfo>(11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11), four_momentum);
                                current_sg->particle_info(pinfo);
                            }
                        } else if (current_pdg == 11 && num_daughter_showers <= 2 && !flag_shower_in &&
                                  !current_sg->flags_any(SegmentFlags::kShowerTopology) &&
                                  !current_sg->flags_any(SegmentFlags::kShowerTrajectory) &&
                                  length > 10*units::cm && !flag_only_showers) {

                            // Matches prototype: score <= 100 guard prevents reclassification
                            // when the classifier is very confident it is an electron
                            if (current_sg->particle_score() <= 100) {
                                double direct_length = segment_track_direct_length(current_sg);
                                if (direct_length >= 34*units::cm ||
                                    (direct_length < 34*units::cm && direct_length > 0.93 * length)) {
                                    auto four_momentum = segment_cal_4mom(current_sg, 13, particle_data, recomb_model, m_mip_dqdx);
                                    auto pinfo = std::make_shared<Aux::ParticleInfo>(13, particle_data->get_particle_mass(13), particle_data->pdg_to_name(13), four_momentum);
                                    current_sg->particle_info(pinfo);
                                }
                            }
                        } else if (current_pdg == 11 && current_sg->flags_any(SegmentFlags::kShowerTrajectory) &&
                                  num_daughter_showers == 1 && !flag_only_showers) {
                            auto pair_result1 = calculate_num_daughter_showers(graph, prev_vtx, current_sg, false);
                            if (pair_result1.second > 3*length && pair_result1.second - length > 12*units::cm) {
                                current_sg->unset_flags(SegmentFlags::kShowerTrajectory);
                                auto four_momentum = segment_cal_4mom(current_sg, 13, particle_data, recomb_model, m_mip_dqdx);
                                auto pinfo = std::make_shared<Aux::ParticleInfo>(13, particle_data->get_particle_mass(13), particle_data->pdg_to_name(13), four_momentum);
                                current_sg->particle_info(pinfo);
                            }
                        }
                    }
                    
                    // Default particle type assignment for undetermined particles from main vertex
                    if (vertex == main_vertex) {
                        int current_pdg = current_sg->has_particle_info() ? current_sg->particle_info()->pdg() : 0;
                        if (current_pdg == 0 && !is_shower) {
                            if (flag_only_showers) {
                                auto four_momentum = segment_cal_4mom(current_sg, 11, particle_data, recomb_model, m_mip_dqdx);
                                auto pinfo = std::make_shared<Aux::ParticleInfo>(11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11), four_momentum);
                                current_sg->particle_info(pinfo);
                            } else {
                                double dqdx_ratio = segment_median_dQ_dx(current_sg) / m_mip_dqdx_median;
                                s_log->trace("examine_direction: cluster {} segment_graph_index={} dqdx_ratio={:.3f} length={:.2f}cm nfits={}",
                                    cluster.ident(), current_sg->get_graph_index(), dqdx_ratio, length/units::cm, current_sg->fits().size());
                                // std::cout << "examine_direction: A segment graph index " << current_sg->get_graph_index() << " dqdx_ratio=" << dqdx_ratio << " length=" << length/units::cm << "cm nfits=" << current_sg->fits().size() << std::endl;
                                if (dqdx_ratio > 1.4) {
                                    auto four_momentum = segment_cal_4mom(current_sg, 2212, particle_data, recomb_model, m_mip_dqdx);
                                    auto pinfo = std::make_shared<Aux::ParticleInfo>(2212, particle_data->get_particle_mass(2212), particle_data->pdg_to_name(2212), four_momentum);
                                    current_sg->particle_info(pinfo);
                                } else {
                                    auto four_momentum = segment_cal_4mom(current_sg, 13, particle_data, recomb_model, m_mip_dqdx);
                                    auto pinfo = std::make_shared<Aux::ParticleInfo>(13, particle_data->get_particle_mass(13), particle_data->pdg_to_name(13), four_momentum);
                                    current_sg->particle_info(pinfo);
                                }
                            }
                        }
                    }
                    
                    current_sg->dir_weak(true);
                    s_log->trace("examine_direction: cluster {} seg assigned: dirsign={} pdg={} length={:.2f}cm",
                        cluster.ident(), current_sg->dirsign(),
                        current_sg->has_particle_info() ? current_sg->particle_info()->pdg() : 0,
                        length/units::cm);
                }
            } else if (current_sg->dirsign() != 0 && !seg_dir_weak(current_sg)) {
                // Strong direction already set
                auto pair_result = calculate_num_daughter_showers(graph, prev_vtx, current_sg);
                int num_daughter_showers = pair_result.first;
                
                int current_pdg = current_sg->has_particle_info() ? current_sg->particle_info()->pdg() : 0;
                if (current_pdg == 2212 && flag_shower_in && num_daughter_showers == 0) {
                    for (auto in_shower : in_showers) {
                        double dqdx_ratio = segment_median_dQ_dx(in_shower) / m_mip_dqdx_median;

                        // std::cout << "examine_direction: B segment graph index " << current_sg->get_graph_index() << " dqdx_ratio=" << dqdx_ratio << " length=" << length/units::cm << "cm nfits=" << current_sg->fits().size() << std::endl;

                        if (dqdx_ratio > 1.3) {
                            auto four_momentum = segment_cal_4mom(in_shower, 2212, particle_data, recomb_model, m_mip_dqdx);
                            auto pinfo = std::make_shared<Aux::ParticleInfo>(2212, particle_data->get_particle_mass(2212), particle_data->pdg_to_name(2212), four_momentum);
                            in_shower->particle_info(pinfo);
                        } else {
                            auto four_momentum = segment_cal_4mom(in_shower, 211, particle_data, recomb_model, m_mip_dqdx);
                            auto pinfo = std::make_shared<Aux::ParticleInfo>(211, particle_data->get_particle_mass(211), particle_data->pdg_to_name(211), four_momentum);
                            in_shower->particle_info(pinfo);
                        }
                        in_shower->unset_flags(SegmentFlags::kShowerTrajectory);
                        pr74_probe_topo_flag(in_shower, "unset", "NeutrinoVertexFinder.cxx:examine_direction-B");
                        in_shower->unset_flags(SegmentFlags::kShowerTopology);
                    }
                }
            }
            
            used_segments.insert(current_sg);
            
            // Find next vertex and add its segments to examination list
            VertexPtr curr_vertex = find_other_vertex(graph, current_sg, prev_vtx);
            if (!curr_vertex || used_vertices.find(curr_vertex) != used_vertices.end()) continue;
            
            if (curr_vertex->descriptor_valid()) {
                auto curr_vd = curr_vertex->get_descriptor();
                // Stable edge-index order: BFS discovery order assigns
                // direction/particle type.
                for (auto edesc : sorted_out_edges(curr_vd, graph)) {
                    SegmentPtr seg = graph[edesc].segment;
                    if (seg) {
                        temp_segments.push_back(std::make_pair(curr_vertex, seg));
                    }
                }
            }
            used_vertices.insert(curr_vertex);
        }
        segments_to_be_examined = temp_segments;
    }
    
    // Find long muon candidates
    bool flag_fill_long_muon = true;
    for (auto seg : segments_in_long_muon) {
        if (seg->cluster() == &cluster) {
            flag_fill_long_muon = false;
            break;
        }
    }
    
    if (flag_fill_long_muon && vertex->descriptor_valid()) {
        auto vd = vertex->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);
        
        for (auto e_it : edge_range) {
            SegmentPtr sg = graph[e_it].segment;
            if (!sg) continue;
            
            double dqdx_ratio = segment_median_dQ_dx(sg) / m_mip_dqdx_median;
            if (dqdx_ratio > 1.3) continue;
            
            VertexPtr vtx = find_other_vertex(graph, sg, vertex);
            if (!vtx) continue;
            
            std::vector<SegmentPtr> acc_segments;
            std::vector<VertexPtr> acc_vertices;
            acc_segments.push_back(sg);
            acc_vertices.push_back(vtx);
            
            auto results = find_cont_muon_segment(graph, sg, vtx, false, m_long_muon_stub_bridge);
            while (results.first != nullptr) {
                acc_segments.push_back(results.first);
                acc_vertices.push_back(results.second);
                results = find_cont_muon_segment(graph, results.first, results.second, false, m_long_muon_stub_bridge);
            }
            
            double total_length = 0, max_length = 0;
            for (auto acc_seg : acc_segments) {
                double length = segment_track_length(acc_seg);
                total_length += length;
                if (length > max_length) max_length = length;
            }
            
            if (total_length > 45*units::cm && max_length > 35*units::cm && acc_segments.size() > 1) {
                s_log->trace("examine_direction: cluster {} found long muon chain: total_length={:.2f}cm max_length={:.2f}cm nsegs={}",
                    cluster.ident(), total_length/units::cm, max_length/units::cm, acc_segments.size());
                for (auto acc_seg : acc_segments) {
                    auto four_momentum = segment_cal_4mom(acc_seg, 13, particle_data, recomb_model, m_mip_dqdx);
                    auto pinfo = std::make_shared<Aux::ParticleInfo>(13, particle_data->get_particle_mass(13), particle_data->pdg_to_name(13), four_momentum);
                    acc_seg->particle_info(pinfo);
                    acc_seg->unset_flags(SegmentFlags::kShowerTrajectory);
                    pr74_probe_topo_flag(acc_seg, "unset", "NeutrinoVertexFinder.cxx:examine_direction-longmu");
                    acc_seg->unset_flags(SegmentFlags::kShowerTopology);
                    segments_in_long_muon.insert(acc_seg);
                }
                for (auto acc_vtx : acc_vertices) {
                    vertices_in_long_muon.insert(acc_vtx);
                }
            }
        }
    }
    
    // Find muon candidate and make others pions
    if (vertex->descriptor_valid()) {
        SegmentPtr muon_sg = nullptr;
        double muon_length = 0;
        std::vector<SegmentPtr> pion_sgs;
        // doc pr/43 round 2 K1 bookkeeping: candidates rejected only by the
        // chain-walk proton veto, and the winner legacy selection would have
        // picked (the demote-all guard's fallback).
        std::vector<SegmentPtr> chain_vetoed_sgs;
        SegmentPtr legacy_muon_sg = nullptr;
        double legacy_muon_length = 0;
        
        auto vd = vertex->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);
        
        for (auto e_it : edge_range) {
            SegmentPtr sg = graph[e_it].segment;
            if (!sg) continue;
            
            int pdg = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;
            if (std::abs(pdg) == 13) {
                if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()) {
                    // doc sbnd_xin/docs/pr/43 round 2 K2 -- legacy skips the
                    // long-muon chain entirely, so it neither competes for
                    // nor claims the vertex muon slot and a second, shorter
                    // pdg-13 arm survives as "the" muon (evt 56463: 14006
                    // stays mu- next to the 411 cm chain 14005+14007).
                    // Knob on: the chain claims the slot with its summed
                    // length; it is never pushed to pion_sgs, so it is never
                    // demoted itself, and other pdg-13 arms demote below.
                    if (m_single_muon_long_muon_claim) {
                        double chain_length = 0;
                        for (auto lm_sg : segments_in_long_muon) {
                            if (lm_sg && lm_sg->cluster() == &cluster)
                                chain_length += segment_track_length(lm_sg);
                        }
                        if (chain_length > muon_length) {
                            muon_length = chain_length;
                            muon_sg = sg;
                        }
                    }
                    continue;
                }

                VertexPtr other_vertex = find_other_vertex(graph, sg, vertex);
                if (!other_vertex || !other_vertex->descriptor_valid()) continue;

                int n_proton = 0;
                auto other_vd = other_vertex->get_descriptor();
                const auto other_edge_range = sorted_out_edges(other_vd, graph);
                for (auto oe_it : other_edge_range) {
                    SegmentPtr other_sg = graph[oe_it].segment;
                    if (other_sg && other_sg->has_particle_info() &&
                        std::abs(other_sg->particle_info()->pdg()) == 2212) {
                        n_proton++;
                    }
                }

                double sg_length = segment_track_length(sg);
                // doc sbnd_xin/docs/pr/43 round 2 K1 -- the 1-hop n_proton
                // veto misses a proton behind a short degree-2 stub (evt
                // 54351: 17007 -> 17005 -> proton 17011).  Knob on: extend
                // the veto through the bounded continuation chain.  The
                // demote-all guard lives after the loop: a chain-vetoed
                // winner is only rejected if some chain-proton-free
                // candidate exists.
                bool chain_proton = false;
                if (m_single_muon_proton_chain_veto && n_proton == 0 && vertex->descriptor_valid()) {
                    chain_proton = segment_chain_has_proton(graph, sg, vertex, m_mip_dqdx_median);
                }
                if (sg_length > muon_length && n_proton == 0 && !chain_proton) {
                    muon_length = sg_length;
                    muon_sg = sg;
                }
                if (sg_length > legacy_muon_length && n_proton == 0) {
                    legacy_muon_length = sg_length;
                    legacy_muon_sg = sg;
                }
                if (chain_proton) chain_vetoed_sgs.push_back(sg);
                pion_sgs.push_back(sg);
            } else if (pdg == 0) {
                VertexPtr other_vertex = find_other_vertex(graph, sg, vertex);
                if (!other_vertex || !other_vertex->descriptor_valid()) continue;
                
                int n_proton = 0;
                auto other_vd = other_vertex->get_descriptor();
                const auto other_edge_range = sorted_out_edges(other_vd, graph);
                for (auto oe_it : other_edge_range) {
                    SegmentPtr other_sg = graph[oe_it].segment;
                    if (other_sg && other_sg->has_particle_info() && 
                        std::abs(other_sg->particle_info()->pdg()) == 2212) {
                        n_proton++;
                    }
                }
                
                if (n_proton > 0) {
                    double dqdx_ratio = segment_median_dQ_dx(sg) / m_mip_dqdx_median;

                    // std::cout << "examine_direction: C segment graph index " << sg->get_graph_index() << " dqdx_ratio=" << dqdx_ratio << " nfits=" << sg->fits().size() << std::endl;

                    if (dqdx_ratio > 1.3) {
                        auto four_momentum = segment_cal_4mom(sg, 2212, particle_data, recomb_model, m_mip_dqdx);
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(2212, particle_data->get_particle_mass(2212), particle_data->pdg_to_name(2212), four_momentum);
                        sg->particle_info(pinfo);
                    } else {
                        auto four_momentum = segment_cal_4mom(sg, 211, particle_data, recomb_model, m_mip_dqdx);
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(211, particle_data->get_particle_mass(211), particle_data->pdg_to_name(211), four_momentum);
                        sg->particle_info(pinfo);
                    }
                }
            }
        }
        
        s_log->trace("examine_direction: cluster {} muon/pion selection: muon_length={:.2f}cm n_muon_candidates={} has_muon={}",
            cluster.ident(), muon_length/units::cm, (int)pion_sgs.size(), (muon_sg != nullptr));

        // doc pr/43 round 2 K1 demote-all guard: if the chain veto rejected
        // every otherwise-qualified candidate, fall back to the legacy
        // winner rather than demoting all arms (a topology legacy would not
        // have produced).
        if (m_single_muon_proton_chain_veto && !muon_sg && legacy_muon_sg) {
            muon_sg = legacy_muon_sg;
            muon_length = legacy_muon_length;
        }

        // Convert non-muon candidates to pions
        for (auto pion_sg : pion_sgs) {
            if (pion_sg == muon_sg) continue;
            auto four_momentum = segment_cal_4mom(pion_sg, 211, particle_data, recomb_model, m_mip_dqdx);
            auto pinfo = std::make_shared<Aux::ParticleInfo>(211, particle_data->get_particle_mass(211), particle_data->pdg_to_name(211), four_momentum);
            pion_sg->particle_info(pinfo);
        }

        // doc pr/43 round 2 K1 -- a chain-vetoed candidate's own
        // continuation stubs (evt 54351: 17005, 2.7 cm, pdg 13 between the
        // demoted 17007 and the proton) must read pion with it; demoting
        // only the head leaves an inconsistent muon stub mid-chain.
        if (m_single_muon_proton_chain_veto) {
            for (auto vet_sg : chain_vetoed_sgs) {
                if (vet_sg == muon_sg) continue;
                for (auto chain_sg : segment_chain_continuation(graph, vet_sg, vertex)) {
                    if (!chain_sg->has_particle_info() ||
                        std::abs(chain_sg->particle_info()->pdg()) != 13) continue;
                    auto four_momentum = segment_cal_4mom(chain_sg, 211, particle_data, recomb_model, m_mip_dqdx);
                    auto pinfo = std::make_shared<Aux::ParticleInfo>(211, particle_data->get_particle_mass(211), particle_data->pdg_to_name(211), four_momentum);
                    chain_sg->particle_info(pinfo);
                }
            }
        }
    }
    
    // Find Michel electrons
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        
        // Check if segment has particle info with mass but no 4-momentum yet, and is not shower topology
        bool has_4mom = sg->has_particle_info();
        if (has_4mom && sg->particle_info()->mass() > 0 && 
            !sg->flags_any(SegmentFlags::kShowerTopology)) {
            
            if (!seg_dir_weak(sg)) {
                // Strong direction - calculate 4-momentum
                int pdg = sg->particle_info()->pdg();
                auto four_momentum = segment_cal_4mom(sg, pdg, particle_data, recomb_model, m_mip_dqdx);
                auto pinfo = std::make_shared<Aux::ParticleInfo>(pdg, particle_data->get_particle_mass(pdg), particle_data->pdg_to_name(pdg), four_momentum);
                sg->particle_info(pinfo);
            } else {
                // Weak direction - need to check endpoint conditions
                // Find the two vertices of this segment
                VertexPtr start_v = nullptr, end_v = nullptr;
                
                for (const auto& vd_it : ordered_nodes(graph)) {
                    VertexPtr vtx = graph[vd_it].vertex;
                    if (!vtx || !vtx->descriptor_valid()) continue;
                    
                    // Check if this vertex is connected to our segment
                    auto vtx_vd = vtx->get_descriptor();
                    const auto vtx_edge_range = sorted_out_edges(vtx_vd, graph);
                    for (auto ve_it : vtx_edge_range) {
                        if (graph[ve_it].segment == sg) {
                            // This vertex is connected to our segment
                            const auto& wcps = sg->wcpts();
                            if (!wcps.empty()) {
                                if (ray_length(Ray{wcps.front().point, vtx->wcpt().point}) < 0.01*units::cm) {
                                    start_v = vtx;
                                } else if (ray_length(Ray{wcps.back().point, vtx->wcpt().point}) < 0.01*units::cm) {
                                    end_v = vtx;
                                }
                            }
                            break;
                        }
                    }
                    if (start_v && end_v) break;
                }
                
                if (!start_v || !end_v) continue;
                
                int dir_sign = sg->dirsign();
                auto fiducial_utils = cluster.grouping()->get_fiducialutils();
                
                // Count segments at start and end vertices
                int num_segs_start = 0, num_segs_end = 0;
                if (start_v->descriptor_valid()) {
                    auto start_vd = start_v->get_descriptor();
                    const auto start_edge_range = sorted_out_edges(start_vd, graph);
                    for (auto se_it : start_edge_range) {
                        if (graph[se_it].segment) num_segs_start++;
                    }
                }
                if (end_v->descriptor_valid()) {
                    auto end_vd = end_v->get_descriptor();
                    const auto end_edge_range = sorted_out_edges(end_vd, graph);
                    for (auto ee_it : end_edge_range) {
                        if (graph[ee_it].segment) num_segs_end++;
                    }
                }
                
                // Check if endpoint is in fiducial volume or is Michel electron candidate
                WireCell::Point end_pt = end_v->fit().valid() ? end_v->fit().point : end_v->wcpt().point;
                WireCell::Point start_pt = start_v->fit().valid() ? start_v->fit().point : start_v->wcpt().point;
                
                bool should_calc = false;
                bool flag_Michel_triggered = false;

                // Case 1: Direction is outward and endpoint is isolated in fiducial volume
                if (dir_sign == 1 && num_segs_end == 1 && fiducial_utils &&
                    fiducial_utils->inside_fiducial_volume(end_pt)) {
                    should_calc = true;
                } else if (dir_sign == -1 && num_segs_start == 1 && fiducial_utils &&
                          fiducial_utils->inside_fiducial_volume(start_pt)) {
                    should_calc = true;
                }
                // Case 2: Check for Michel electron topology at end vertex (2 segments, one is shower)
                else if (num_segs_end == 2 && end_v->descriptor_valid()) {
                    bool flag_Michel = false;
                    auto end_vd = end_v->get_descriptor();
                    const auto end_edge_range = sorted_out_edges(end_vd, graph);
                    for (auto ee_it : end_edge_range) {
                        SegmentPtr other_sg = graph[ee_it].segment;
                        if (!other_sg || other_sg == sg) continue;
                        // Check if other segment is a shower (kShowerTrajectory flag or electron PDG)
                        if (other_sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                            (other_sg->has_particle_info() && std::abs(other_sg->particle_info()->pdg()) == 11)) {
                            flag_Michel = true;
                            break;
                        }
                    }
                    if (flag_Michel) { should_calc = true; flag_Michel_triggered = true; }
                }
                // Case 3: Check for Michel electron topology at start vertex (2 segments, one is shower)
                else if (num_segs_start == 2 && start_v->descriptor_valid()) {
                    bool flag_Michel = false;
                    auto start_vd = start_v->get_descriptor();
                    const auto start_edge_range = sorted_out_edges(start_vd, graph);
                    for (auto se_it : start_edge_range) {
                        SegmentPtr other_sg = graph[se_it].segment;
                        if (!other_sg || other_sg == sg) continue;
                        // Check if other segment is a shower (kShowerTrajectory flag or electron PDG)
                        if (other_sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                            (other_sg->has_particle_info() && std::abs(other_sg->particle_info()->pdg()) == 11)) {
                            flag_Michel = true;
                            break;
                        }
                    }
                    if (flag_Michel) { should_calc = true; flag_Michel_triggered = true; }
                }

                if (should_calc) {
                    int pdg = sg->particle_info()->pdg();
                    // A stopped proton cannot produce a Michel electron; if Michel topology is
                    // present at the endpoint, this track is more likely a misidentified muon.
                    if (flag_Michel_triggered && pdg == 2212) pdg = 13;
                    auto four_momentum = segment_cal_4mom(sg, pdg, particle_data, recomb_model, m_mip_dqdx);
                    auto pinfo = std::make_shared<Aux::ParticleInfo>(pdg, particle_data->get_particle_mass(pdg), particle_data->pdg_to_name(pdg), four_momentum);
                    sg->particle_info(pinfo);
                }
            }
        }
    }
    
    // Final pass: ensure every shower segment has particle_info set to electron.
    // Mirrors prototype get_particle_type() which returns 11 for any shower regardless
    // of whether determine_dir_shower_trajectory was called.
    set_default_shower_particle_info(graph, cluster, particle_data, recomb_model, main_vertex);

    // doc sbnd_xin/docs/pr/40 round 4 F8 -- same call site (last word before
    // shower_clustering_with_nv), same per-cluster main_vertex.  No-op when
    // m_muon_multi_proton_pion is false.
    override_muon_multi_proton_pion(graph, cluster, particle_data, recomb_model, main_vertex);

    // doc sbnd_xin/docs/pr/40 round 6 F14 -- same call-site convention as F8.
    // No-op when m_michel_stem_muon_rescue is false.
    override_michel_stem_muon(graph, cluster, particle_data, recomb_model, main_vertex);

    return examine_maps(graph, cluster);
}

bool PatternAlgorithms::eliminate_short_vertex_activities(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, std::set<SegmentPtr>& existing_segments, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_updated = false;
    bool flag_continue = true;
        
    while (flag_continue) {
        flag_continue = false;
        std::set<SegmentPtr> to_be_removed_segments;
        std::set<VertexPtr> to_be_removed_vertices;
        
        // Iterate through all edges (segments) in the graph.  Stable
        // edge-index order: the first removable segment wins (remove +
        // restart), so iteration order is a topology decision.
        for (const auto& ed : ordered_edges(graph)) {
            SegmentPtr sg = graph[ed].segment;
            if (!sg || sg->cluster() != &cluster) continue;
            if (existing_segments.find(sg) != existing_segments.end()) continue;

            // Get the two vertices of this segment
            auto [v1, v2] = find_vertices(graph, sg);
            if (!v1 || !v2) continue;
            
            // Count segments at each vertex
            int num_segs_v1 = 0, num_segs_v2 = 0;
            if (v1->descriptor_valid()) {
                auto v1d = v1->get_descriptor();
                const auto v1_edge_range = sorted_out_edges(v1d, graph);
                for (auto ve_it : v1_edge_range) {
                    if (graph[ve_it].segment) num_segs_v1++;
                }
            }
            if (v2->descriptor_valid()) {
                auto v2d = v2->get_descriptor();
                const auto v2_edge_range = sorted_out_edges(v2d, graph);
                for (auto ve_it : v2_edge_range) {
                    if (graph[ve_it].segment) num_segs_v2++;
                }
            }
            
            double length = segment_track_direct_length(sg);
            
            // Check Case 1: v1 has 1 segment, v2 has >=3 segments
            if (num_segs_v1 == 1 && num_segs_v2 >= 3) {
                if (length < 0.36*units::cm) {
                    to_be_removed_segments.insert(sg);
                    to_be_removed_vertices.insert(v1);
                    flag_continue = true;
                    break;
                } else if (length < 0.5*units::cm && num_segs_v2 > 3) {
                    to_be_removed_segments.insert(sg);
                    to_be_removed_vertices.insert(v1);
                    flag_continue = true;
                    break;
                }
            }
            // Check Case 2: v2 has 1 segment, v1 has >=3 segments
            else if (num_segs_v2 == 1 && num_segs_v1 >= 3) {
                if (length < 0.36*units::cm) {
                    to_be_removed_segments.insert(sg);
                    to_be_removed_vertices.insert(v2);
                    flag_continue = true;
                    break;
                } else if (length < 0.5*units::cm && num_segs_v1 > 3) {
                    to_be_removed_segments.insert(sg);
                    to_be_removed_vertices.insert(v2);
                    flag_continue = true;
                    break;
                }
            }
            
            // Check Case 3: Very short segments (< 0.1 cm) connected to main_vertex
            // Prototype always removes v2, relying on set ordering so main_vertex is v1.
            // Toolkit ordering is non-deterministic, so explicitly remove the non-main vertex.
            if (!flag_continue) {
                if ((v1 == main_vertex && num_segs_v1 > 1) || (v2 == main_vertex && num_segs_v2 > 1)) {
                    if (length < 0.1*units::cm) {
                        to_be_removed_segments.insert(sg);
                        VertexPtr to_remove = (v1 == main_vertex) ? v2 : v1;
                        to_be_removed_vertices.insert(to_remove);
                        flag_continue = true;
                        break;
                    }
                }
            }
            
            // Check Case 4: Isolated vertex close to another segment
            if (!flag_continue) {
                WireCell::Point v1_pt = v1->fit().valid() ? v1->fit().point : v1->wcpt().point;
                WireCell::Point v2_pt = v2->fit().valid() ? v2->fit().point : v2->wcpt().point;
                // auto v1_wpid = dv->contained_by(v1_pt);
                // auto v2_wpid = dv->contained_by(v2_pt);
                
                if (num_segs_v1 == 1 && num_segs_v2 > 1 && v2->descriptor_valid()) {
                    auto v2d = v2->get_descriptor();
                    const auto v2_edge_range = sorted_out_edges(v2d, graph);
                    for (auto ve_it : v2_edge_range) {
                        SegmentPtr sg1 = graph[ve_it].segment;
                        if (!sg1 || sg1 == sg) continue;
                        
                        auto [dis, closest_pt] = segment_get_closest_point(sg1, v1_pt, "fit");
                        double seg_length = segment_track_length(sg);
                        
                        if (dis < 0.36*units::cm) {
                            to_be_removed_segments.insert(sg);
                            to_be_removed_vertices.insert(v1);
                            flag_continue = true;
                            break;
                        } else if ((v2 == main_vertex && dis < 0.45*units::cm && seg_length < 0.45*units::cm)) {
                            to_be_removed_segments.insert(sg);
                            to_be_removed_vertices.insert(v1);
                            flag_continue = true;
                            break;
                        }
                    }
                } else if (num_segs_v2 == 1 && num_segs_v1 > 1 && v1->descriptor_valid()) {
                    auto v1d = v1->get_descriptor();
                    const auto v1_edge_range = sorted_out_edges(v1d, graph);
                    for (auto ve_it : v1_edge_range) {
                        SegmentPtr sg1 = graph[ve_it].segment;
                        if (!sg1 || sg1 == sg) continue;
                        
                        auto [dis, closest_pt] = segment_get_closest_point(sg1, v2_pt, "fit");
                        double seg_length = segment_track_length(sg);
                        
                        if (dis < 0.36*units::cm) {
                            to_be_removed_segments.insert(sg);
                            to_be_removed_vertices.insert(v2);
                            flag_continue = true;
                            break;
                        } else if ((v1 == main_vertex && dis < 0.45*units::cm && seg_length < 0.45*units::cm)) {
                            to_be_removed_segments.insert(sg);
                            to_be_removed_vertices.insert(v2);
                            flag_continue = true;
                            break;
                        }
                    }
                }
            }
            
            // Check Case 5: Segment not in existing_segments and all points close to existing segments
            if (!flag_continue && existing_segments.find(sg) == existing_segments.end() && length > 0.45*units::cm) {
                const auto& wcpts = sg->wcpts();
                int n_good = 0;

                // Hoisted out of the point x existing-segment double loop: the
                // "is this segment still in the graph" test is a pure existence
                // check, so it can be answered once from a snapshot instead of
                // re-walking every edge per (point, existing segment) pair.
                // POINTER identity, like the scan it replaces -- see the note
                // on this function in NeutrinoPatternBase.h.  Unordered because
                // it is only ever count()-ed, never iterated.
                std::unordered_set<SegmentPtr> segs_in_graph;
                for (const auto& ed2 : ordered_edges(graph)) {
                    if (graph[ed2].segment) segs_in_graph.insert(graph[ed2].segment);
                }

                for (size_t i = 0; i < wcpts.size(); i++) {
                    WireCell::Point pt = wcpts[i].point;
                    auto wpid = dv->contained_by(pt);
                    if (wpid.face() == -1 || wpid.apa() == -1) continue;

                    double dis_u = 1e9, dis_v = 1e9, dis_w = 1e9;

                    for (auto existing_sg : existing_segments) {
                        if (!segs_in_graph.count(existing_sg)) continue;

                        auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(existing_sg, pt, wpid.apa(), wpid.face(), "fit");
                        if (dist_u < dis_u) dis_u = dist_u;
                        if (dist_v < dis_v) dis_v = dist_v;
                        if (dist_w < dis_w) dis_w = dist_w;
                    }
                    
                    if ((dis_u > 0.45*units::cm || dis_v > 0.45*units::cm || dis_w > 0.45*units::cm)) {
                        n_good++;
                    }
                }
                
                if (n_good == 0) {
                    to_be_removed_segments.insert(sg);
                    if (num_segs_v1 == 1) to_be_removed_vertices.insert(v1);
                    if (num_segs_v2 == 1) to_be_removed_vertices.insert(v2);
                }
            }
            
            if (flag_continue) break;
        }
        
        // Remove segments and vertices
        for (auto sg : to_be_removed_segments) {
            flag_updated = true;
            remove_segment(graph, sg);
        }
        for (auto vtx : to_be_removed_vertices) {
            remove_vertex(graph, vtx);
        }
    }
    
    return flag_updated;
}


bool PatternAlgorithms::fit_vertex(Facade::Cluster& cluster, VertexPtr vertex, VertexPtr main_vertex, std::vector<SegmentPtr>& sg_set, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Allow to move 1.5 cm - create MyFCN object with constraint parameters
    MyFCN fcn(vertex, true, 0.43*units::cm, 1.5*units::cm, 0.9*units::cm, 6*units::cm);
    
    // std::cout << "MyFCN: " << sg_set.size() << " segments to be fitted for vertex " << vertex->fit_index() << std::endl;

    // Add all segments to the fitting.  With m_fit_vertex_min_seg_length set
    // (doc sbnd_xin/docs/pr/9 sec. 11 F3c; 0 = legacy include-all), segments
    // shorter than the cut are excluded from the position fit so that
    // vertex-activity stubs cannot drag the vertex -- unless fewer than 2
    // segments survive, in which case the full set is kept.
    std::vector<SegmentPtr> fit_segments(sg_set.begin(), sg_set.end());
    if (m_fit_vertex_min_seg_length > 0) {
        // Measure in WCPT space: the fit cloud of a freshly created
        // vertex-activity stub can spread well past its 0.6 cm wcpt path
        // (evt 172230: fit span > 1 cm before any refit), so fit-based
        // lengths cannot identify the stub.  The wcpt path is the segment's
        // graph-topology extent.
        auto wcpt_path_length = [](SegmentPtr sg) {
            const auto& wcpts = sg->wcpts();
            double len = 0;
            for (size_t i = 0; i + 1 < wcpts.size(); ++i) {
                len += (wcpts[i+1].point - wcpts[i].point).magnitude();
            }
            return len;
        };
        std::vector<SegmentPtr> long_segments;
        for (auto sg : fit_segments) {
            const double wlen = wcpt_path_length(sg);
            if (sg_set.size() > 2) {
                s_log->trace("fit_vertex: cluster {} candidate segment wcpt_len={:.2f} cm fit_len={:.2f} cm (cut {:.2f} cm)",
                    cluster.ident(), wlen/units::cm,
                    segment_track_direct_length(sg)/units::cm,
                    m_fit_vertex_min_seg_length/units::cm);
            }
            if (wlen >= m_fit_vertex_min_seg_length) {
                long_segments.push_back(sg);
            }
        }
        if (long_segments.size() < fit_segments.size()) {
            if (long_segments.size() <= 2) {
                // All extra legs are below the cut: this is effectively a
                // two-leg (or fewer) vertex whose position the plain two-leg
                // fit already determined.  Refitting on the re-tracked fit
                // clouds reproduces the stub-induced drag even with the stub
                // excluded from the MyFCN accumulation (evt 172230: identical
                // 2.4 cm displacement -- the stub corrupts the fit through
                // do_multi_tracking charge competition, not through its own
                // PCA term), so skip the fit outright.
                s_log->trace("fit_vertex: cluster {} skipping vertex fit: only {} segment(s) >= {:.2f} cm (of {})",
                    cluster.ident(), long_segments.size(),
                    m_fit_vertex_min_seg_length/units::cm, fit_segments.size());
                return false;
            }
            s_log->trace("fit_vertex: cluster {} excluding {} short segment(s) (< {:.2f} cm) from the vertex fit",
                cluster.ident(), fit_segments.size() - long_segments.size(),
                m_fit_vertex_min_seg_length/units::cm);
            fit_segments = long_segments;
        }
    }
    // doc sbnd_xin/docs/pr/51 round 7: robust per-leg direction windows
    // (design at MyFCN.h RobustParams; knob block at NeutrinoPatternBase.h
    // m_mvfit_robust).  Gate on vertex == main_vertex directly --
    // enforce_two_track_fit is only set AFTER the AddSegment loop below, so
    // it cannot be consulted for this.
    const bool robust_here = m_mvfit_robust && (!m_mvfit_main_only || vertex == main_vertex);
    if (robust_here) {
        MyFCN::RobustParams rp;
        rp.on = true;
        rp.min_len = m_mvfit_min_len;
        rp.rin_margin = m_mvfit_rin_margin;
        rp.rout_frac = m_mvfit_rout_frac;
        rp.rout_min = m_mvfit_rout_min;
        rp.rout_max = m_mvfit_rout_max;
        rp.angle = m_mvfit_angle;
        rp.min_pts = m_mvfit_min_pts;
        rp.min_aniso = m_mvfit_min_aniso;
        rp.prior_range = m_mvfit_prior_range;
        fcn.set_robust(rp);
    }
    for (auto it = fit_segments.begin(); it != fit_segments.end(); it++) {
        fcn.AddSegment(*it);
    }

    // If this is the main vertex, enforce two track fit
    if (vertex == main_vertex) fcn.set_enforce_two_track_fit(true);
    
    // Perform vertex fitting
    std::pair<bool, Facade::geo_point_t> results = fcn.FitVertex();
    
    // Get grouping for charge calculation
    auto grouping = cluster.grouping();
    if (!grouping) {
        // && short-circuits exactly like the previous `if (results.first)`:
        // UpdateInfo is not called when the fit failed.  Its return now
        // propagates so a guard-aborted (no-op) update is not reported as a
        // successful fit -- see MyFCN.h and docs/pr/28 sec. 3.3.
        return results.first && fcn.UpdateInfo(results.second, cluster, track_fitter, dv);
    }
    
    // Get transform for coordinate conversion
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    // double cluster_t0 = cluster.get_cluster_t0();
    
    // Get old and new vertex positions
    Facade::geo_point_t old_pos = vertex->fit().point;
    Facade::geo_point_t new_pos = results.second;
    
    // Get APA/face for old and new positions
    auto old_wpid = dv->contained_by(old_pos);
    auto new_wpid = dv->contained_by(new_pos);
    
    // Calculate average charge at old and new positions
    double old_charge = 0;
    double new_charge = 0;
    
    // Sampling radius: the prototype (NeutrinoID_improve_vertex.h:22-23) calls
    // get_ave_3d_charge(pt) with no radius, taking the ToyCTPointCloud.h:35
    // default of 0.3 cm; Facade_Grouping.h:315 carries the same 0.3 cm default.
    // This passed 0.6 cm explicitly (since 6b66163e, uncommented), sampling 8x
    // the volume and averaging in charge from well outside the vertex.  Left
    // implicit so the two defaults cannot drift apart.  docs/pr/28 sec. 3.3.
    if (old_wpid.apa() != -1 && old_wpid.face() != -1) {
        old_charge = grouping->get_ave_3d_charge(old_pos, old_wpid.apa(), old_wpid.face());
    }

    if (new_wpid.apa() != -1 && new_wpid.face() != -1) {
        new_charge = grouping->get_ave_3d_charge(new_pos, new_wpid.apa(), new_wpid.face());
    }
    
    // Check charge conditions - if new position has much lower charge, keep old position.
    // prototype (NeutrinoID_improve_vertex.h:26,28): 5000 / 8000 e = 5/43 and
    // 8/43 of the MIP-median dQ/dx integrated over 1 cm.  Expressed on
    // m_mip_dqdx_median so they scale with the detector; exactly 5000 / 8000
    // at the uBooNE default (215000/43 and 344000/43, FP-exact).  docs/pr/2 sec 7.3.
    if (new_charge < m_mip_dqdx_median * units::cm * 5.0 / 43.0 && new_charge < 0.4*old_charge) {
        results.second = old_pos;
    } else if (new_charge < m_mip_dqdx_median * units::cm * 8.0 / 43.0 && new_charge < 0.6*old_charge) {
        // Reduce the strength - keep old position
        results.second = old_pos;
        new_charge = old_charge;
    }

    // doc sbnd_xin/docs/pr/51 round 7: a charge-vetoed robust fit must
    // re-seat exactly as production would -- otherwise UpdateInfo below would
    // straighten the path out to the substituted outer centers (~6 cm) around
    // an unmoved vertex, a larger graph rewrite than production's for a fit
    // that was rejected.  Exact FP equality detects the veto: both veto
    // branches assign old_pos by copy (an unvetoed fit landing exactly on
    // old_pos restores too, which then IS production behavior).
    if (fcn.robust_substituted() > 0 &&
        results.second.x() == old_pos.x() &&
        results.second.y() == old_pos.y() &&
        results.second.z() == old_pos.z()) {
        fcn.restore_production_pca();
    }

    // Update vertex and segment information with fitted position.  The && is
    // the former `if (results.first)` guard; UpdateInfo's return propagates so
    // that a guard-aborted no-op is not reported as a fit (docs/pr/28 sec 3.3).
    return results.first && fcn.UpdateInfo(results.second, cluster, track_fitter, dv);
}


// doc sbnd_xin/docs/pr/50: main-vertex kink-consistency snap.  Design and
// guards documented at the m_vertex_kink_snap member block
// (NeutrinoPatternBase.h).  Toolkit-only; no prototype counterpart.
bool PatternAlgorithms::snap_main_vertex_to_kink(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    if (!m_vertex_kink_snap) return false;
    if (!main_vertex || main_vertex->cluster() != &cluster) return false;
    if (!main_vertex->descriptor_valid()) return false;
    // G1: never steal a protected junction (pr/48 TEB breaks etc.).
    if (main_vertex->flags_any(VertexFlags::kProtectedBreak)) return false;

    const WireCell::Point vpt = main_vertex->wcpt().point;
    const auto vd = main_vertex->get_descriptor();

    // Incident segments in stable edge-index order.
    std::vector<SegmentPtr> incident;
    for (auto edesc : sorted_out_edges(vd, graph)) {
        SegmentPtr sg = graph[edesc].segment;
        if (sg && sg->cluster() == &cluster) incident.push_back(sg);
    }
    if (incident.size() < 2) return false;  // G2: pass-through needs two arms

    // Oriented wcpt path (index 0 at the main vertex) + near-vertex outward
    // chord per segment.  A segment not touching the main vertex per
    // find_vertices (shouldn't happen for an incident edge) is skipped.
    struct Arm {
        SegmentPtr seg;
        std::vector<WireCell::Point> pts;
        WireCell::Vector chord;   // unit, vertex -> ~3 cm out
        bool reversed{false};
    };
    std::vector<Arm> arms;
    for (auto& sg : incident) {
        auto [vf, vb] = find_vertices(graph, sg);
        Arm a; a.seg = sg;
        const auto& wcpts = sg->wcpts();
        if (wcpts.size() < 2) continue;
        if (vf == main_vertex) {
            for (const auto& w : wcpts) a.pts.push_back(w.point);
        } else if (vb == main_vertex) {
            a.reversed = true;
            for (auto it = wcpts.rbegin(); it != wcpts.rend(); ++it) a.pts.push_back(it->point);
        } else {
            continue;
        }
        double acc = 0;
        size_t k = 1;
        for (; k + 1 < a.pts.size(); k++) {
            acc += (a.pts[k] - a.pts[k-1]).magnitude();
            if (acc >= 3*units::cm) break;
        }
        a.chord = (a.pts[k] - a.pts[0]).norm();
        arms.push_back(std::move(a));
    }
    if (arms.size() < 2) return false;

    // Trace-level geometry dump for tuning (one line per near-vertex path
    // point; trace builds only).
    if (s_log->should_log(spdlog::level::trace)) {
        for (size_t ia = 0; ia < arms.size(); ia++) {
            double acc = 0;
            for (size_t k = 0; k < arms[ia].pts.size() && acc < 8*units::cm; k++) {
                if (k) acc += (arms[ia].pts[k] - arms[ia].pts[k-1]).magnitude();
                s_log->trace("snap_path: arm {} k={} arc={:.2f} cm p=({:.2f},{:.2f},{:.2f})",
                    ia, k, acc/units::cm,
                    arms[ia].pts[k].x()/units::cm, arms[ia].pts[k].y()/units::cm, arms[ia].pts[k].z()/units::cm);
            }
        }
    }

    // Bragg-hot test per arm (computed lazily, cached): median dQ/dx of the
    // arm's fits within 3 cm of the vertex end, as a multiple of the MIP
    // median.  -1 = not enough samples (test abstains).
    std::vector<double> hot_cache(arms.size(), -2);
    auto arm_hot_ratio = [&](size_t ib) -> double {
        if (hot_cache[ib] > -2) return hot_cache[ib];
        double out = -1;
        const auto& fits = arms[ib].seg->fits();
        std::vector<double> dqdx;
        if (fits.size() >= 2) {
            double acc = 0;
            if (!arms[ib].reversed) {
                for (size_t k = 0; k < fits.size(); k++) {
                    if (k) acc += (fits[k].point - fits[k-1].point).magnitude();
                    if (acc > 3*units::cm) break;
                    if (fits[k].dx > 0) dqdx.push_back(fits[k].dQ / fits[k].dx);
                }
            } else {
                for (size_t k = fits.size(); k-- > 0;) {
                    if (k + 1 < fits.size()) acc += (fits[k].point - fits[k+1].point).magnitude();
                    if (acc > 3*units::cm) break;
                    if (fits[k].dx > 0) dqdx.push_back(fits[k].dQ / fits[k].dx);
                }
            }
        }
        if (dqdx.size() >= 2 && m_mip_dqdx_median > 0) {
            std::nth_element(dqdx.begin(), dqdx.begin() + dqdx.size()/2, dqdx.end());
            out = dqdx[dqdx.size()/2] / m_mip_dqdx_median;
        }
        hot_cache[ib] = out;
        return out;
    };

    // Scan every arm, every candidate turn (the strongest turn is not
    // necessarily the right corner -- in 172230 a secondary wiggle at
    // 4.9 cm out-turned the true corner at 2.4 cm by 3 deg).  Selection
    // among survivors is by PASS-THROUGH EVIDENCE: smallest bendV, then
    // larger turn, then smaller arc -- all deterministic scalars.
    int best_arm = -1;
    VertexKinkScanResult best;
    double best_bendV = 181;
    for (size_t ia = 0; ia < arms.size(); ia++) {
        const auto cands = path_scan_vertex_kink(arms[ia].pts, m_vks_min_dis, m_vks_radius,
                                                 m_vks_skirt, m_vks_baseline, m_vks_angle, m_vks_min_arm);
        for (const auto& res : cands) {
            // G3 pass-through: the stub V->K must continue another incident
            // arm straight through V.  bendV = deviation from straight.
            const auto u = (arms[ia].pts[res.idx] - vpt).norm();
            double bendV = 181;
            int opp = -1;
            for (size_t ib = 0; ib < arms.size(); ib++) {
                if (ib == ia) continue;
                const double ang = std::acos(std::max(-1.0, std::min(1.0, u.dot(arms[ib].chord)))) / M_PI * 180.0;
                const double dev = 180.0 - ang;
                if (dev < bendV) { bendV = dev; opp = static_cast<int>(ib); }
            }
            if (opp < 0 || bendV > m_vks_collinear) {
                s_log->trace("snap_main_vertex_to_kink: cluster {} arm {} turn={:.1f} deg at {:.2f} cm K=({:.2f},{:.2f},{:.2f}) declined (bendV={:.1f} > collinear {})",
                    cluster.ident(), ia, res.turn_deg, res.arc/units::cm,
                    arms[ia].pts[res.idx].x()/units::cm, arms[ia].pts[res.idx].y()/units::cm, arms[ia].pts[res.idx].z()/units::cm,
                    bendV, m_vks_collinear);
                continue;
            }
            // G4 strength margin: the interior turn must clearly beat the
            // residual bend at the vertex.
            if (res.turn_deg < bendV + m_vks_margin) {
                s_log->trace("snap_main_vertex_to_kink: cluster {} turn={:.1f} deg at {:.2f} cm declined (margin vs bendV={:.1f})",
                    cluster.ident(), res.turn_deg, res.arc/units::cm, bendV);
                continue;
            }
            // G5 fit-miss: the snap only makes sense when the FITTED
            // trajectory does not already honor the image corner -- a
            // genuinely modeled kink (real secondary scatter, real junction)
            // has fit points passing through it, while the 172230 failure's
            // fit rounds the corner and misses it.  min distance from K to
            // any incident arm's fit points; below m_vks_fit_miss => the
            // fit already goes there => decline.  (A Bragg-hot veto was
            // tried and MISFIRES on the failure class itself: with the
            // vertex misplaced onto the muon, the arm near V reads the
            // proton's charge, 1.91x MIP in 172230 -- local dQ/dx cannot
            // tell a real hot junction from misassigned charge.
            // m_vks_hot_ratio survives as an optional extra veto, default
            // 0 = off.)
            {
                const auto& K = arms[ia].pts[res.idx];
                double fit_miss = 1e9;
                for (const auto& a2 : arms) {
                    for (const auto& f : a2.seg->fits()) {
                        const double d = (f.point - K).magnitude();
                        if (d < fit_miss) fit_miss = d;
                    }
                }
                s_log->trace("snap_main_vertex_to_kink: cluster {} arm {} cand turn={:.1f} deg at {:.2f} cm K=({:.2f},{:.2f},{:.2f}) bendV={:.1f} fit_miss={:.2f} cm",
                    cluster.ident(), ia, res.turn_deg, res.arc/units::cm,
                    K.x()/units::cm, K.y()/units::cm, K.z()/units::cm, bendV, fit_miss/units::cm);
                if (fit_miss < m_vks_fit_miss) {
                    continue;
                }
            }
            if (m_vks_hot_ratio > 0) {
                const double hr = arm_hot_ratio(opp);
                if (hr > m_vks_hot_ratio) {
                    s_log->trace("snap_main_vertex_to_kink: cluster {} turn={:.1f} deg declined (opposite arm Bragg-hot, med/mip={:.2f})",
                        cluster.ident(), res.turn_deg, hr);
                    continue;
                }
            }
            // G6: an existing vertex of this cluster already sits at the turn.
            {
                const auto& K = arms[ia].pts[res.idx];
                bool taken = false;
                for (const auto& nd : ordered_nodes(graph)) {
                    VertexPtr ov = graph[nd].vertex;
                    if (!ov || ov->cluster() != &cluster || ov == main_vertex) continue;
                    if ((ov->wcpt().point - K).magnitude() < m_vks_min_dis) { taken = true; break; }
                }
                if (taken) {
                    s_log->trace("snap_main_vertex_to_kink: cluster {} turn={:.1f} deg at {:.2f} cm declined (vertex already there)",
                        cluster.ident(), res.turn_deg, res.arc/units::cm);
                    continue;
                }
            }
            const bool better = (best_arm < 0)
                || (bendV < best_bendV)
                || (bendV == best_bendV && res.turn_deg > best.turn_deg)
                || (bendV == best_bendV && res.turn_deg == best.turn_deg && res.arc < best.arc);
            if (better) {
                best_arm = static_cast<int>(ia);
                best = res;
                best_bendV = bendV;
            }
        }
    }
    if (best_arm < 0) return false;

    const WireCell::Point K = arms[best_arm].pts[best.idx];
    auto [ok, segs, vtx_new] = break_segment(graph, arms[best_arm].seg, K, particle_data, recomb_model, dv);
    if (!ok || !vtx_new) {
        SPDLOG_LOGGER_DEBUG(s_log, "snap_main_vertex_to_kink: cluster {} break_segment declined at ({:.2f},{:.2f},{:.2f})",
            cluster.ident(), K.x()/units::cm, K.y()/units::cm, K.z()/units::cm);
        return false;
    }
    // break_segment does not set the cluster: without it the new vertex is
    // invisible to every candidate loop (same gotcha as break_two_end_dqdx).
    vtx_new->cluster(&cluster);
    // Shield against examine_vertices / ES2 / ES3 re-merging the deliberate
    // break (the pr/48 kProtectedBreak contract).
    vtx_new->set_flags(VertexFlags::kProtectedBreak);
    (void)segs;

    SPDLOG_LOGGER_DEBUG(s_log, "snap_main_vertex_to_kink: SNAP cluster {} old=({:.2f},{:.2f},{:.2f}) new=({:.2f},{:.2f},{:.2f}) turn={:.1f} deg arc={:.2f} cm bendV={:.1f} deg",
        cluster.ident(),
        vpt.x()/units::cm, vpt.y()/units::cm, vpt.z()/units::cm,
        K.x()/units::cm, K.y()/units::cm, K.z()/units::cm,
        best.turn_deg, best.arc/units::cm, best_bendV);

    main_vertex = vtx_new;
    // One full refit so the trajectory re-anchors on the corner before
    // improve_vertex polishes (same argument pattern as its refits).
    track_fitter.do_multi_tracking(true, true, true, m_fit_exclusion, false, &cluster);
    return true;
}


void PatternAlgorithms::improve_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_search_vertex_activity, bool flag_final_vertex){
    s_log->trace("improve_vertex: cluster {} flag_search_vertex_activity={} flag_final_vertex={}", cluster.ident(), flag_search_vertex_activity, flag_final_vertex);

    IndexedVertexSet fitted_vertices;  // order by stable graph index, not pointer address
    std::set<SegmentPtr> existing_segments;   // pointer-keyed on purpose, see header
    
    // Check if all segments are showers, no need to fit vertex with only two legs
    bool flag_skip_two_legs = false;
    {
        int ntracks = 0;
        for (const auto& ed : ordered_edges(graph)) {
            SegmentPtr sg = graph[ed].segment;
            if (!sg || sg->cluster() != &cluster) continue;
            existing_segments.insert(sg);
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology) ||
                             (sg->particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
            if (!is_shower) ntracks++;
        }
        if (ntracks == 0)
            flag_skip_two_legs = true; // all showers
    }
    

    bool flag_found_vertex_activities = false;
    
    // Search for vertex activities
    if (flag_search_vertex_activity) {
    detv_dump("iv:entry", graph);
        if (examine_structure_4(main_vertex, flag_final_vertex, graph, cluster, track_fitter, dv)) {
            flag_found_vertex_activities = true;
            track_fitter.do_multi_tracking(true, true, true, m_fit_exclusion, false, &cluster);
            detv_dump("iv:examine_structure_4a", graph);
        }
    }
    
    bool flag_update_fit = false;
    
    // Find and fit vertices
    for (auto v : ordered_nodes(graph)) {
        VertexPtr vtx = graph[v].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;

        // Get segments connected to this vertex.  Order matters: it feeds
        // MyFCN's normal-equation accumulation in fit_vertex, so use the
        // stable edge-index order, not boost::out_edges pointer order.
        std::vector<SegmentPtr> vertex_segments;
        for (auto edesc : sorted_out_edges(v, graph)) {
            SegmentPtr sg = graph[edesc].segment;
            if (sg && sg->cluster() == &cluster) {
                vertex_segments.push_back(sg);
            }
        }

        if (vertex_segments.size() <= 2 && vtx != main_vertex) continue;
        
        int ntracks = 0, nshowers = 0;
        int n_long_muons = 0;
        for (auto sg : vertex_segments) {
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology) ||
                             (sg->particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
            if (is_shower) nshowers++;
            else ntracks++;
            if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()) n_long_muons++;
        }
        
        if (ntracks == 0 && vtx != main_vertex) continue;
        if (flag_skip_two_legs && vertex_segments.size() <= 2) continue;
        
        auto wcp_save = vtx->wcpt();

        s_log->trace("improve_vertex: cluster {} fitting vertex ({:.2f}, {:.2f}, {:.2f}) nsegs={}",
            cluster.ident(),
            vtx->wcpt().point.x()/units::cm, vtx->wcpt().point.y()/units::cm, vtx->wcpt().point.z()/units::cm,
            vertex_segments.size());
        bool flag_update = fit_vertex(cluster, vtx, main_vertex, vertex_segments, track_fitter, dv);
        if (flag_update) fitted_vertices.insert(vtx);
        if (flag_update) {
            flag_update_fit = true;

            double tmp_dis = std::sqrt(std::pow(wcp_save.point.x() - vtx->wcpt().point.x(), 2) +
                                      std::pow(wcp_save.point.y() - vtx->wcpt().point.y(), 2) +
                                      std::pow(wcp_save.point.z() - vtx->wcpt().point.z(), 2));
            s_log->trace("improve_vertex: cluster {} fit_vertex done, vertex moved {:.3f} cm -> ({:.2f}, {:.2f}, {:.2f})",
                cluster.ident(), tmp_dis/units::cm,
                vtx->wcpt().point.x()/units::cm, vtx->wcpt().point.y()/units::cm, vtx->wcpt().point.z()/units::cm);

            if (tmp_dis > 0.5*units::cm) { // if the vertex moved far, refit
                track_fitter.do_multi_tracking(true, true, true, m_fit_exclusion, false, &cluster);
                fit_vertex(cluster, vtx, main_vertex, vertex_segments, track_fitter, dv);
                s_log->trace("improve_vertex: cluster {} second fit_vertex done -> ({:.2f}, {:.2f}, {:.2f})",
                    cluster.ident(),
                    vtx->wcpt().point.x()/units::cm, vtx->wcpt().point.y()/units::cm, vtx->wcpt().point.z()/units::cm);
            }
        } else {
            s_log->trace("improve_vertex: cluster {} fit_vertex made no update", cluster.ident());
        }

        (void)nshowers;
        (void)n_long_muons;
    }
    
    if (flag_update_fit) {
        // Do the overall fit again
        track_fitter.do_multi_tracking(true, true, true, m_fit_exclusion, false, &cluster);
        detv_dump("iv:post_fitloop_dmt", graph);

        bool flag_keep_main_vertex = false;
        Facade::geo_point_t main_vtx_pt;
        if (main_vertex != nullptr) {
            main_vtx_pt = main_vertex->fit().point;
            flag_keep_main_vertex = true;
        }
        
        examine_vertices(graph, cluster, track_fitter, dv, main_vertex);
        detv_dump("iv:examine_vertices_a", graph);
        
        if (flag_keep_main_vertex) {
            // Check if main_vertex still exists in graph
            bool found_main_vertex = false;
            for (auto v : ordered_nodes(graph)) {
                if (graph[v].vertex == main_vertex) {
                    found_main_vertex = true;
                    break;
                }
            }

            if (!found_main_vertex) {
                double min_dis = 1e9;
                for (auto v : ordered_nodes(graph)) {
                    VertexPtr vtx = graph[v].vertex;
                    if (!vtx || vtx->cluster() != &cluster) continue;
                    double dis = std::sqrt(std::pow(vtx->fit().point.x() - main_vtx_pt.x(), 2) + 
                                          std::pow(vtx->fit().point.y() - main_vtx_pt.y(), 2) + 
                                          std::pow(vtx->fit().point.z() - main_vtx_pt.z(), 2));
                    if (dis < min_dis) {
                        min_dis = dis;
                        main_vertex = vtx;
                    }
                }
            }
        }
    }
    
    std::vector<VertexPtr> refit_vertices;
    flag_update_fit = false;
    
    if (flag_search_vertex_activity) {
        // Search for vertex activities again
        if (!flag_found_vertex_activities) {
            if (examine_structure_4(main_vertex, flag_final_vertex, graph, cluster, track_fitter, dv)) {
                flag_found_vertex_activities = true;
                
                 detv_dump("iv:examine_structure_4b", graph);
                // Get segments connected to main_vertex
                std::vector<SegmentPtr> main_vertex_segments;
                for (auto v : ordered_nodes(graph)) {
                    if (graph[v].vertex == main_vertex) {
                        // Stable order.  Byte-identical here (the vector feeds
                        // only the size()==3 test below), but improve_vertex
                        // now has no unordered edge iteration left.
                        for (auto edesc : sorted_out_edges(v, graph)) {
                            SegmentPtr sg = graph[edesc].segment;
                            if (sg && sg->cluster() == &cluster) {
                                main_vertex_segments.push_back(sg);
                            }
                        }
                        break;
                    }
                }
                
                if (main_vertex_segments.size() == 3) refit_vertices.push_back(main_vertex);
                flag_update_fit = true;
            }
        }
        
        for (auto v : ordered_nodes(graph)) {
            VertexPtr vtx = graph[v].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;

            // Get segments connected to this vertex.  Stable edge-index order:
            // this vector feeds search_for_vertex_activities, whose decisions
            // must not depend on boost pointer order.
            std::vector<SegmentPtr> vertex_segments;
            for (auto edesc : sorted_out_edges(v, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (sg && sg->cluster() == &cluster) {
                    vertex_segments.push_back(sg);
                }
            }

            if (vertex_segments.size() <= 2 && vtx != main_vertex) continue;

            int ntracks = 0, nshowers = 0;
            for (auto sg : vertex_segments) {
                bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology) ||
                                 (sg->particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
                if (is_shower) nshowers++;
                else ntracks++;
            }

            if (ntracks == 0 && vtx != main_vertex) continue;
            if (vertices_in_long_muon.find(vtx) != vertices_in_long_muon.end()) continue;
            if (vtx == main_vertex && flag_found_vertex_activities) continue;
            if (flag_skip_two_legs && vertex_segments.size() <= 2) continue;

            double search_range = 1.5*units::cm;
            if (vertex_segments.size() == 1) search_range = 3.0*units::cm;

            bool flag_update = search_for_vertex_activities(graph, vtx, vertex_segments, cluster, track_fitter, dv, search_range);
            detv_dump("iv:search_vtx_act", graph);

            if (flag_update) {
                // Get updated segments.  Stable order, matching the first loop
                // at :2247: this set feeds only the size()==3 test below, so
                // the change is byte-identical, but an iterated unordered edge
                // set must not be left on a path that could later feed output.
                vertex_segments.clear();
                for (auto edesc : sorted_out_edges(v, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (sg && sg->cluster() == &cluster) {
                        vertex_segments.push_back(sg);
                    }
                }
                if (vertex_segments.size() == 3) refit_vertices.push_back(vtx);
                flag_update_fit = true;
            }
            (void)nshowers;
        }
        
        if (flag_update_fit) {
            // Do the overall fit again
            track_fitter.do_multi_tracking(true, true, true, m_fit_exclusion, false, &cluster);
            detv_dump("iv:post_search_dmt", graph);
            flag_update_fit = false;

            // Redo the fit
            for (auto vtx : refit_vertices) {
                std::vector<SegmentPtr> vertex_segments;
                for (auto v : ordered_nodes(graph)) {
                    if (graph[v].vertex == vtx) {
                        // Stable order: feeds MyFCN accumulation via fit_vertex.
                        for (auto edesc : sorted_out_edges(v, graph)) {
                            SegmentPtr sg = graph[edesc].segment;
                            if (sg && sg->cluster() == &cluster) {
                                vertex_segments.push_back(sg);
                            }
                        }
                        break;
                    }
                }

                s_log->trace("improve_vertex: cluster {} refit pass, fitting vertex ({:.2f}, {:.2f}, {:.2f}) nsegs={}",
                    cluster.ident(),
                    vtx->wcpt().point.x()/units::cm, vtx->wcpt().point.y()/units::cm, vtx->wcpt().point.z()/units::cm,
                    vertex_segments.size());
                bool flag_update = fit_vertex(cluster, vtx, main_vertex, vertex_segments, track_fitter, dv);
                if (flag_update) fitted_vertices.insert(vtx);
                if (flag_update) {
                    flag_update_fit = true;
                    s_log->trace("improve_vertex: cluster {} refit fit_vertex updated vertex -> ({:.2f}, {:.2f}, {:.2f})",
                        cluster.ident(),
                        vtx->wcpt().point.x()/units::cm, vtx->wcpt().point.y()/units::cm, vtx->wcpt().point.z()/units::cm);
                }
            }
            if (flag_update_fit) {
                track_fitter.do_multi_tracking(true, true, true, m_fit_exclusion, false, &cluster);
            }
        }
        
        // Eliminate short tracks
        if (eliminate_short_vertex_activities(graph, cluster, main_vertex, existing_segments, track_fitter, dv)) {
            track_fitter.do_multi_tracking(true, true, true, m_fit_exclusion, false, &cluster);
            detv_dump("iv:final_dmt", graph);
        }
        
        // Determine directions for segments (ordered_edges: this mutates
        // per-segment direction/PID, so iteration order is an input).
        for (const auto& ed : ordered_edges(graph)) {
            SegmentPtr sg1 = graph[ed].segment;
            if (!sg1 || sg1->cluster() != &cluster) continue;

            if (!sg1->particle_info()) {
                segment_is_shower_topology(sg1, false, m_mip_dqdx_median, m_shower_topo_demote_len, m_shower_topo_reset, m_shower_topo_dqdx_guard);

                VertexPtr start_v = nullptr, end_v = nullptr;
                auto source_v = boost::source(ed, graph);
                auto target_v = boost::target(ed, graph);
                
                auto& wcpts = sg1->wcpts();
                if (!wcpts.empty()) {
                    if ((graph[source_v].vertex->wcpt().point - wcpts.front().point).magnitude() < 0.01*units::cm) {
                        start_v = graph[source_v].vertex;
                        end_v = graph[target_v].vertex;
                    } else {
                        end_v = graph[source_v].vertex;
                        start_v = graph[target_v].vertex;
                    }
                }
                
                if (start_v && end_v) {
                    int start_n = boost::out_degree(source_v, graph);
                    int end_n = boost::out_degree(target_v, graph);
                    
                    // doc pr/32 §11 F2(a): prototype
                    // NeutrinoID_improve_vertex.h:248 READS the stored flag
                    // here.  The toolkit recomputes, which also has the side
                    // effect of SETTING kShowerTrajectory on a segment the
                    // stored label said was not one.
                    if (pr32_shower_traj_gate(sg1, m_shower_traj_recheck_parity, m_mip_dqdx)) {
                        // Trajectory shower
                        segment_determine_shower_direction_trajectory(sg1, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, false, track_pid_options());
                    } else {
                        segment_determine_dir_track(sg1, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, false, track_pid_options());
                    }
                }
            }
        }
    } else { // flag_search_vertex_activity
        for (auto vtx : fitted_vertices) {
            // Find vertex descriptor
            for (auto v : ordered_nodes(graph)) {
                if (graph[v].vertex != vtx) continue;
                
                // Stable edge-index order.  This loop MUTATES segments --
                // segment_is_shower_topology sets kShowerTopology and
                // segment_determine_dir_track writes dirsign/particle_info --
                // and boost::out_edges iterates a setS edge list in *pointer*
                // order (PRTrajectoryView.h:151), so the results were
                // run-to-run nondeterministic.  docs/pr/28 sec 3.3 / sec 9.
                for (auto edesc : sorted_out_edges(v, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (!sg || sg->cluster() != &cluster) continue;

                    if (!sg->particle_info()) segment_is_shower_topology(sg, false, m_mip_dqdx_median, m_shower_topo_demote_len, m_shower_topo_reset, m_shower_topo_dqdx_guard);

                    VertexPtr start_v = nullptr, end_v = nullptr;
                    auto source_v = boost::source(edesc, graph);
                    auto target_v = boost::target(edesc, graph);
                    
                    auto& wcpts = sg->wcpts();
                    if (!wcpts.empty()) {
                        if ((graph[source_v].vertex->wcpt().point - wcpts.front().point).magnitude() < 0.01*units::cm) {
                            start_v = graph[source_v].vertex;
                            end_v = graph[target_v].vertex;
                        } else {
                            end_v = graph[source_v].vertex;
                            start_v = graph[target_v].vertex;
                        }
                    }
                    
                    bool flag_print = false;
                    bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology) ||
                                     (sg->particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
                    if (start_v && end_v && !is_shower) {
                        // Track
                        int start_n = boost::out_degree(source_v, graph);
                        int end_n = boost::out_degree(target_v, graph);
                        segment_determine_dir_track(sg, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, flag_print, track_pid_options());
                    }
                }
            }
        }
    }
    
    // Handle special cases for main_vertex segments
    if (main_vertex != nullptr && main_vertex->cluster() == &cluster) {
        // Find main_vertex descriptor
        for (auto v : ordered_nodes(graph)) {
            if (graph[v].vertex != main_vertex) continue;
            
            // Stable edge-index order -- same reason as the fitted_vertices
            // loop above: this one mutates harder still (unset_flags plus
            // segment_determine_dir_track on every branch).  docs/pr/28 sec 9.
            for (auto edesc : sorted_out_edges(v, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg || sg->cluster() != &cluster) continue;

                std::pair<int, double> pair_result = calculate_num_daughter_showers(graph, main_vertex, sg, false);
                
                double medium_dQdx = segment_median_dQ_dx(sg);
                // doc pr/32 §11 F2: outer gate = prototype
                // NeutrinoID_improve_vertex.h:287 (stored flag); inner test =
                // :288, `!sg->is_shower_trajectory()` at the 10 cm default with
                // the function's own 50000-analog scale.  With the gate ALSO
                // recomputing (today's path), an inner call at those same
                // parameters would negate the line above it and the body would
                // be dead -- which is why the toolkit's inner call runs at
                // 1.0 cm.  The two move together or not at all.
                if ((pair_result.first <= 2 || (medium_dQdx/m_mip_dqdx_median > 1.6 && pair_result.first <= 3)) && pr32_shower_traj_gate(sg, m_shower_traj_recheck_parity, m_mip_dqdx)) {
                    const bool inner_traj = m_shower_traj_recheck_parity
                        ? segment_is_shower_trajectory(sg, 10*units::cm,  m_mip_dqdx)
                        : segment_is_shower_trajectory(sg, 1.0*units::cm, m_mip_dqdx_median);
                    if (!inner_traj) {
                        g_pr32_audit.f2_body_runs.fetch_add(1, std::memory_order_relaxed);
                        VertexPtr start_v = nullptr, end_v = nullptr;
                        auto source_v = boost::source(edesc, graph);
                        auto target_v = boost::target(edesc, graph);
                        
                        auto& wcpts = sg->wcpts();
                        if (!wcpts.empty()) {
                            if ((graph[source_v].vertex->wcpt().point - wcpts.front().point).magnitude() < 0.01*units::cm) {
                                start_v = graph[source_v].vertex;
                                end_v = graph[target_v].vertex;
                            } else {
                                end_v = graph[source_v].vertex;
                                start_v = graph[target_v].vertex;
                            }
                        }
                        
                        if (start_v && end_v) {
                            int start_n = boost::out_degree(source_v, graph);
                            int end_n = boost::out_degree(target_v, graph);
                            segment_determine_dir_track(sg, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, false, track_pid_options());
                        }
                    }
                }
                
                // Examine topology case
                if (pair_result.first == 1 && segment_is_shower_topology(sg, false, m_mip_dqdx_median, m_shower_topo_demote_len, m_shower_topo_reset, m_shower_topo_dqdx_guard)) {
                    pr74_probe_topo_reexam(sg, "enter");
                    int dir_save = sg->dirsign();
                    
                    VertexPtr start_v = nullptr, end_v = nullptr;
                    auto source_v = boost::source(edesc, graph);
                    auto target_v = boost::target(edesc, graph);
                    
                    auto& wcpts = sg->wcpts();
                    if (!wcpts.empty()) {
                        if ((graph[source_v].vertex->wcpt().point - wcpts.front().point).magnitude() < 0.01*units::cm) {
                            start_v = graph[source_v].vertex;
                            end_v = graph[target_v].vertex;
                        } else {
                            end_v = graph[source_v].vertex;
                            start_v = graph[target_v].vertex;
                        }
                    }
                    
                    if (start_v && end_v) {
                        sg->unset_flags(SegmentFlags::kShowerTopology);
                        int start_n = boost::out_degree(source_v, graph);
                        int end_n = boost::out_degree(target_v, graph);
                        segment_determine_dir_track(sg, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, false, track_pid_options());
                        pr74_probe_topo_reexam(sg, "after-pid");

                        if ((sg->particle_info() && sg->particle_info()->pdg() == 2212 && sg->particle_score() < 0.09) ||
                            (sg->particle_info() && sg->particle_info()->pdg() == 13 && sg->particle_score() < 0.06)) {
                            sg->unset_flags(SegmentFlags::kShowerTopology);
                            pr74_probe_topo_reexam(sg, "re-demote(no-pdg-write)");
                        } else {
                            sg->set_flags(SegmentFlags::kShowerTopology);
                            if (!sg->particle_info()) {
                                sg->particle_info() = std::make_shared<Aux::ParticleInfo>();
                            }
                            sg->particle_info()->set_pdg(11);
                            pr40_probe_setpdg(sg, 11, "NeutrinoVertexFinder.cxx:topo-escape(M5)");
                            sg->particle_score(100);
                            sg->dirsign(dir_save);
                            sg->particle_info()->set_mass(particle_data->get_particle_mass(11));
                        }
                    }
                }
                
                if (flag_skip_two_legs && existing_segments.find(sg) == existing_segments.end()) {
                    VertexPtr start_v = nullptr, end_v = nullptr;
                    auto source_v = boost::source(edesc, graph);
                    auto target_v = boost::target(edesc, graph);
                    
                    auto& wcpts = sg->wcpts();
                    if (!wcpts.empty()) {
                        if ((graph[source_v].vertex->wcpt().point - wcpts.front().point).magnitude() < 0.01*units::cm) {
                            start_v = graph[source_v].vertex;
                            end_v = graph[target_v].vertex;
                        } else {
                            end_v = graph[source_v].vertex;
                            start_v = graph[target_v].vertex;
                        }
                    }
                    
                    if (start_v && end_v) {
                        int start_n = boost::out_degree(source_v, graph);
                        int end_n = boost::out_degree(target_v, graph);
                        segment_determine_dir_track(sg, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, false, track_pid_options());
                        
                        if ((!sg->particle_info() || seg_dir_weak(sg)) && medium_dQdx/m_mip_dqdx_median < 1.3) {
                            if (!sg->particle_info()) {
                                sg->particle_info() = std::make_shared<Aux::ParticleInfo>();
                            }
                            sg->particle_info()->set_pdg(11);
                            pr40_probe_setpdg(sg, 11, "NeutrinoVertexFinder.cxx:flag_skip_two_legs");
                            sg->particle_score(100);
                            sg->particle_info()->set_mass(particle_data->get_particle_mass(11));
                        }
                    }
                }
            }
            break;
        }
    }

    s_log->trace("improve_vertex: cluster {} done", cluster.ident());
}

void PatternAlgorithms::determine_main_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    s_log->trace("determine_main_vertex: cluster {}", cluster.ident());

    // Find the main vertex - check if we have only showers
    bool flag_save_only_showers = true;
    for (auto v : ordered_nodes(graph)) {
        VertexPtr vtx = graph[v].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;

        auto results = examine_main_vertex_candidate(graph, vtx);
        bool flag_in = std::get<0>(results);
        int ntracks = std::get<1>(results);
        // int nshowers = std::get<2>(results);

        if (!flag_in) {
            if (ntracks > 0) {
                flag_save_only_showers = false;
                break;
            }
        }
    }
    MS t_scan_only_showers(Clock::now() - t0); t0 = Clock::now();

    s_log->trace("determine_main_vertex: cluster {} flag_save_only_showers={}", cluster.ident(), flag_save_only_showers);

    // Improve vertex if not only showers and cluster is main cluster
    if (!flag_save_only_showers && cluster.get_flag(Facade::Flags::main_cluster)) {
        s_log->trace("determine_main_vertex: cluster {} calling improve_vertex + fix_maps_shower_in_track_out", cluster.ident());
        improve_vertex(graph, cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, track_fitter, dv, particle_data, recomb_model, false);
        // Fix maps with shower in and track out
        fix_maps_shower_in_track_out(graph, cluster);
    }
    MS t_improve_vertex(Clock::now() - t0); t0 = Clock::now();

    // Build map of vertex candidates and their track/shower counts
    std::map<VertexPtr, std::pair<int, int>> map_vertex_track_shower;
    std::vector<VertexPtr> main_vertex_candidates;

    for (auto v : ordered_nodes(graph)) {
        VertexPtr vtx = graph[v].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;

        auto results = examine_main_vertex_candidate(graph, vtx);
        bool flag_in = std::get<0>(results);
        int ntracks = std::get<1>(results);
        int nshowers = std::get<2>(results);

        if (!flag_in) {
            map_vertex_track_shower[vtx] = std::make_pair(ntracks, nshowers);
        }
    }
    MS t_build_candidate_map(Clock::now() - t0); t0 = Clock::now();

    // Select main vertex candidates based on topology
    if (flag_save_only_showers) {
        // For all showers case, add vertices with 1 segment first
        for (auto v : ordered_nodes(graph)) {
            VertexPtr vtx = graph[v].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;

            int num_segs = 0;
            if (vtx->descriptor_valid()) {
                auto vd = vtx->get_descriptor();
                const auto edge_range = sorted_out_edges(vd, graph);
                for (auto e_it : edge_range) {
                    if (graph[e_it].segment) num_segs++;
                }
            }

            if (num_segs == 1) {
                main_vertex_candidates.push_back(vtx);
            }
        }

        // Add remaining candidates in insertion order
        for (auto v : ordered_nodes(graph)) {
            VertexPtr vtx = graph[v].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            if (map_vertex_track_shower.find(vtx) == map_vertex_track_shower.end()) continue;
            if (std::find(main_vertex_candidates.begin(), main_vertex_candidates.end(), vtx) == main_vertex_candidates.end()) {
                main_vertex_candidates.push_back(vtx);
            }
        }
    } else {
        // For mixed case, only add vertices with tracks in insertion order
        for (auto v : ordered_nodes(graph)) {
            VertexPtr vtx = graph[v].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            auto it = map_vertex_track_shower.find(vtx);
            if (it != map_vertex_track_shower.end() && it->second.first > 0) {
                main_vertex_candidates.push_back(vtx);
            }
        }
    }
    MS t_select_candidates(Clock::now() - t0); t0 = Clock::now();

    // doc pr/32 §11 F4 (was P12): the prototype records the per-cluster
    // candidate list HERE, after both branches have pushed and before any
    // filtering (NeutrinoID_track_shower.h:1332), and exposes it via
    // get_map_cluster_candidate_vertices() (NeutrinoID.h:1720).  The toolkit
    // has no equivalent, which is why doc pr/27 §6 records "candidate list,
    // per-cluster and global -> nothing".
    //
    // Ported as a per-vertex flag rather than a new container: it needs no
    // plumbing, it travels with the vertex, and the per-cluster list is
    // recoverable as "the flagged vertices of that cluster".  DIAGNOSTIC ONLY
    // -- nothing in clus reads kMainCandidate except PrDisplayDump, exactly as
    // in the prototype where every consumer is an app-level tree filler.
    if (m_main_vertex_candidate_flag) {
        for (auto vtx : main_vertex_candidates) {
            vtx->set_flags(PR::VertexFlags::kMainCandidate);
            g_pr32_audit.f4_flagged.fetch_add(1, std::memory_order_relaxed);
        }
    }

    s_log->trace("determine_main_vertex: cluster {} ncandidates={}", cluster.ident(), main_vertex_candidates.size());

    // Determine main vertex based on candidates
    if (flag_save_only_showers) {
        if (main_vertex_candidates.size() > 0) {
            s_log->trace("determine_main_vertex: cluster {} all-showers path, calling compare_main_vertices_all_showers", cluster.ident());
            SPDLOG_LOGGER_TRACE(s_log, "determine_main_vertex: cluster {} all-showers, {} candidates",
                cluster.get_cluster_id(), main_vertex_candidates.size());
            main_vertex = compare_main_vertices_all_showers(graph, cluster, main_vertex_candidates, track_fitter, dv, particle_data, recomb_model);
        } else {
            s_log->trace("determine_main_vertex: cluster {} all-showers but no candidates, early return", cluster.ident());
            if (m_perf) {
                MS t_total_ms(Clock::now() - t_total);
                SPDLOG_LOGGER_TRACE(s_log,
                    "determine_main_vertex timing (early return, all-showers no candidates): "
                    "scan_only_showers={:.3f}ms improve_vertex={:.3f}ms "
                    "build_candidate_map={:.3f}ms select_candidates={:.3f}ms TOTAL={:.3f}ms",
                    t_scan_only_showers.count(), t_improve_vertex.count(),
                    t_build_candidate_map.count(), t_select_candidates.count(), t_total_ms.count());
            }
            return;
        }
    } else {
        // Examine main vertex candidates to filter and identify back-to-back tracks
        s_log->trace("determine_main_vertex: cluster {} calling examine_main_vertices_local", cluster.ident());
        examine_main_vertices_local(graph, main_vertex_candidates, particle_data, recomb_model);
        s_log->trace("determine_main_vertex: cluster {} after examine_main_vertices_local, ncandidates={}", cluster.ident(), main_vertex_candidates.size());

        for (auto vtx : main_vertex_candidates) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            std::string seg_list;
            if (vtx->descriptor_valid()) {
                auto vd = vtx->get_descriptor();
                const auto edge_range = sorted_out_edges(vd, graph);
                for (auto e_it : edge_range) {
                    SegmentPtr sg = graph[e_it].segment;
                    if (sg) seg_list += std::to_string(sg->id()) + ", ";
                }
            }
            SPDLOG_LOGGER_TRACE(s_log, "determine_main_vertex: cluster {} candidate pos=({:.1f},{:.1f},{:.1f}) connecting to: {}",
                cluster.get_cluster_id(), vtx_pt.x()/units::cm, vtx_pt.y()/units::cm, vtx_pt.z()/units::cm, seg_list);
        }

        if (main_vertex_candidates.size() == 1) {
            s_log->trace("determine_main_vertex: cluster {} single candidate, selecting directly", cluster.ident());
            main_vertex = main_vertex_candidates.front();
        } else if (main_vertex_candidates.size() > 1) {
            s_log->trace("determine_main_vertex: cluster {} multiple candidates, calling compare_main_vertices", cluster.ident());
            main_vertex = compare_main_vertices(graph, cluster, main_vertex_candidates);
        } else {
            if (m_perf) {
                MS t_total_ms(Clock::now() - t_total);
                SPDLOG_LOGGER_TRACE(s_log,
                    "determine_main_vertex timing (early return, no candidates): "
                    "scan_only_showers={:.3f}ms improve_vertex={:.3f}ms "
                    "build_candidate_map={:.3f}ms select_candidates={:.3f}ms TOTAL={:.3f}ms",
                    t_scan_only_showers.count(), t_improve_vertex.count(),
                    t_build_candidate_map.count(), t_select_candidates.count(), t_total_ms.count());
            }
            return;
        }
    }
    MS t_select_main_vertex(Clock::now() - t0); t0 = Clock::now();

    // Examine structure for non-shower cases
    if (!flag_save_only_showers) {
        examine_structure_final(graph, main_vertex, cluster, track_fitter, dv);
    }
    MS t_examine_structure_final(Clock::now() - t0); t0 = Clock::now();

    // Examine directions
    bool flag_check = examine_direction(graph, main_vertex, main_vertex, vertices_in_long_muon, segments_in_long_muon, particle_data, recomb_model, false);
    if (!flag_check) {
        SPDLOG_LOGGER_TRACE(s_log, "determine_main_vertex: cluster {} inconsistency for track directions",
            cluster.get_cluster_id());
    }
    MS t_examine_direction(Clock::now() - t0);

    if (m_perf) {
        MS t_total_ms(Clock::now() - t_total);
        SPDLOG_LOGGER_TRACE(s_log,
            "determine_main_vertex timing: "
            "scan_only_showers={:.3f}ms improve_vertex={:.3f}ms "
            "build_candidate_map={:.3f}ms select_candidates={:.3f}ms "
            "select_main_vertex={:.3f}ms examine_structure_final={:.3f}ms "
            "examine_direction={:.3f}ms ",
            t_scan_only_showers.count(), t_improve_vertex.count(),
            t_build_candidate_map.count(), t_select_candidates.count(),
            t_select_main_vertex.count(), t_examine_structure_final.count(),
            t_examine_direction.count());
        SPDLOG_LOGGER_TRACE(s_log, "determine_main_vertex timing: TOTAL={:.3f}ms", t_total_ms.count());
    }

    {
        WireCell::Point vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
        std::string seg_list;
        if (main_vertex->descriptor_valid()) {
            auto vd = main_vertex->get_descriptor();
            const auto edge_range = sorted_out_edges(vd, graph);
            for (auto e_it : edge_range) {
                SegmentPtr sg = graph[e_it].segment;
                if (sg) seg_list += std::to_string(sg->id()) + ", ";
            }
        }
        SPDLOG_LOGGER_TRACE(s_log, "determine_main_vertex: cluster {} main vertex pos=({:.1f},{:.1f},{:.1f}) connecting to: {}",
            cluster.get_cluster_id(), vtx_pt.x()/units::cm, vtx_pt.y()/units::cm, vtx_pt.z()/units::cm, seg_list);
    }

    if (main_vertex) {
        s_log->trace("determine_main_vertex: cluster {} done, main_vertex ({:.2f}, {:.2f}, {:.2f})",
            cluster.ident(),
            main_vertex->wcpt().point.x()/units::cm, main_vertex->wcpt().point.y()/units::cm,
            main_vertex->wcpt().point.z()/units::cm);
    } else {
        s_log->trace("determine_main_vertex: cluster {} done, main_vertex is null", cluster.ident());
    }

   
}

void PatternAlgorithms::change_daughter_type(Graph& graph, VertexPtr vertex, SegmentPtr segment, int particle_type, double mass, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    // Find the other vertex of this segment
    VertexPtr other_vtx = find_other_vertex(graph, segment, vertex);
    if (!other_vtx) return;
    
    // Get vertex point
    WireCell::Point other_vtx_pt = other_vtx->fit().valid() ? other_vtx->fit().point : other_vtx->wcpt().point;
    
    // Calculate direction from the other vertex
    Facade::geo_vector_t dir1 = segment_cal_dir_3vector(segment, other_vtx_pt, 15*units::cm);
    
    // Check if other vertex has valid descriptor
    if (!other_vtx->descriptor_valid()) return;
    
    // Iterate through all segments connected to the other vertex
    auto other_vd = other_vtx->get_descriptor();
    const auto edge_range = sorted_out_edges(other_vd, graph);
    
    for (auto e_it : edge_range) {
        SegmentPtr sg1 = graph[e_it].segment;
        if (!sg1 || sg1 == segment) continue;
        
        // Skip if already the same particle type
        int current_pdg = sg1->has_particle_info() ? sg1->particle_info()->pdg() : 0;
        if (current_pdg == particle_type) continue;
        
        // Skip if shower trajectory or has strong direction
        if (sg1->flags_any(SegmentFlags::kShowerTrajectory)) continue;
        if (!seg_dir_weak(sg1) && sg1->dirsign() != 0) continue;
        
        // Check shower topology case: long segments with no direction.
        // Matches prototype: if large-shower-topology → try 170°, then fall through to 165° check.
        // All other segments (small-shower-topology, non-topology) → skip.
        if (sg1->flags_any(SegmentFlags::kShowerTopology) &&
            sg1->dirsign() == 0 &&
            segment_track_length(sg1) > 40*units::cm) {

            // Calculate direction at 40cm
            Facade::geo_vector_t dir2 = segment_cal_dir_3vector(sg1, other_vtx_pt, 40*units::cm);

            if (dir1.magnitude() > 0 && dir2.magnitude() > 0) {
                double cos_angle = std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0);
                double angle = std::acos(cos_angle) / M_PI * 180.0;

                if (angle > 170) {
                    // Change particle type and mass
                    if (!sg1->particle_info()) {
                        sg1->particle_info() = std::make_shared<Aux::ParticleInfo>();
                    }
                    sg1->particle_info()->set_pdg(particle_type);
                    pr40_probe_setpdg(sg1, particle_type, "NeutrinoVertexFinder.cxx:change_daughter_type(170deg)");
                    sg1->particle_info()->set_mass(mass);
                    pr74_probe_topo_flag(sg1, "unset", "NeutrinoVertexFinder.cxx:change_daughter_type");
                    sg1->unset_flags(SegmentFlags::kShowerTopology);

                    // Recursively propagate changes
                    change_daughter_type(graph, other_vtx, sg1, particle_type, mass, particle_data, recomb_model);
                    VertexPtr sg1_other_vtx = find_other_vertex(graph, sg1, other_vtx);
                    if (sg1_other_vtx) {
                        change_daughter_type(graph, sg1_other_vtx, sg1, particle_type, mass, particle_data, recomb_model);
                    }
                }
            }
            // Fall through to general 165° check below (matches prototype)
        } else {
            // Not a large-shower-topology segment → skip (matches prototype's else continue)
            continue;
        }

        // General check: only reached for large-shower-topology segments (length > 40cm implies > 10cm)
        if (segment_track_length(sg1) > 10*units::cm) {
            // Calculate direction at 15cm
            Facade::geo_vector_t dir2 = segment_cal_dir_3vector(sg1, other_vtx_pt, 15*units::cm);
            
            if (dir1.magnitude() > 0 && dir2.magnitude() > 0) {
                double cos_angle = std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0);
                double angle = std::acos(cos_angle) / M_PI * 180.0;
                
                if (angle > 165) {
                    // Change particle type and mass
                    if (!sg1->particle_info()) {
                        sg1->particle_info() = std::make_shared<Aux::ParticleInfo>();
                    }
                    sg1->particle_info()->set_pdg(particle_type);
                    pr40_probe_setpdg(sg1, particle_type, "NeutrinoVertexFinder.cxx:change_daughter_type(165deg)");
                    sg1->particle_info()->set_mass(mass);

                    // Recursively propagate changes
                    change_daughter_type(graph, other_vtx, sg1, particle_type, mass, particle_data, recomb_model);
                    VertexPtr sg1_other_vtx = find_other_vertex(graph, sg1, other_vtx);
                    if (sg1_other_vtx) {
                        change_daughter_type(graph, sg1_other_vtx, sg1, particle_type, mass, particle_data, recomb_model);
                    }
                }
            }
        }
    }
}

void PatternAlgorithms::examine_main_vertices_local(Graph& graph, std::vector<VertexPtr>& vertices, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    s_log->trace("examine_main_vertices_local: nvertices={}", vertices.size());
    if (vertices.size() == 1) return;

    double max_length = 0;
    IndexedVertexSet tmp_vertices;  // order by stable graph index, not pointer address

    for (auto vtx : vertices) {
        if (!vtx || !vtx->descriptor_valid()) continue;

        // Count segments connected to this vertex
        int num_segs = 0;
        auto vd = vtx->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);
        for (auto e_it : edge_range) {
            if (graph[e_it].segment) num_segs++;
        }
        
        // If only 1 segment, add to tmp_vertices
        if (num_segs == 1) {
            tmp_vertices.insert(vtx);
        } else {
            // Check pairs of segments for back-to-back tracks
            std::set<SegmentPtr> used_segments;
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            
            // Collect all segments and check pairs
            std::vector<SegmentPtr> vertex_segments;
            for (auto e_it : edge_range) {
                SegmentPtr sg = graph[e_it].segment;
                if (sg) vertex_segments.push_back(sg);
            }
            
            for (size_t i = 0; i < vertex_segments.size(); i++) {
                SegmentPtr sg1 = vertex_segments[i];
                double length1 = segment_track_length(sg1);
                if (length1 < 10*units::cm) continue;
                
                Facade::geo_vector_t dir1 = segment_cal_dir_3vector(sg1, vtx_pt, 15*units::cm);
                Facade::geo_vector_t dir3 = segment_cal_dir_3vector(sg1, vtx_pt, 30*units::cm);
                
                if (length1 > max_length) max_length = length1;
                
                for (size_t j = i + 1; j < vertex_segments.size(); j++) {
                    SegmentPtr sg2 = vertex_segments[j];
                    double length2 = segment_track_length(sg2);
                    if (length2 < 10*units::cm) continue;
                    
                    Facade::geo_vector_t dir2 = segment_cal_dir_3vector(sg2, vtx_pt, 15*units::cm);
                    Facade::geo_vector_t dir4 = segment_cal_dir_3vector(sg2, vtx_pt, 30*units::cm);
                    
                    double angle1 = 0, angle2 = 0;
                    if (dir1.magnitude() > 0 && dir2.magnitude() > 0) {
                        double cos_angle = std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0);
                        angle1 = std::acos(cos_angle) / M_PI * 180.0;
                    }
                    if (dir3.magnitude() > 0 && dir4.magnitude() > 0) {
                        double cos_angle = std::clamp(dir3.dot(dir4) / (dir3.magnitude() * dir4.magnitude()), -1.0, 1.0);
                        angle2 = std::acos(cos_angle) / M_PI * 180.0;
                    }
                    
                    int pdg1 = sg1->has_particle_info() ? sg1->particle_info()->pdg() : 0;
                    int pdg2 = sg2->has_particle_info() ? sg2->particle_info()->pdg() : 0;
                    
                    // Check for back-to-back muon tracks
                    if ((angle1 > 165 || angle2 > 165) && 
                        (pdg1 == 13 || pdg2 == 13) && 
                        (length1 > 30*units::cm || length2 > 30*units::cm)) {
                        used_segments.insert(sg1);
                        used_segments.insert(sg2);
                    }
                    // Check for back-to-back proton tracks
                    else if ((angle1 > 170 || angle2 > 170) &&
                            ((pdg1 == 2212 && (pdg2 == 0 || pdg2 == 2212)) ||
                             (pdg2 == 2212 && (pdg1 == 0 || pdg1 == 2212))) &&
                            (length1 > 20*units::cm && length2 > 20*units::cm)) {
                        used_segments.insert(sg1);
                        used_segments.insert(sg2);
                    }
                }
            }
            
            // If we found back-to-back tracks, check remaining segments
            if (used_segments.size() > 0) {
                bool flag_skip = true;
                
                for (auto e_it : edge_range) {
                    SegmentPtr sg1 = graph[e_it].segment;
                    if (!sg1) continue;
                    if (used_segments.find(sg1) != used_segments.end()) continue;
                    
                    double length = segment_track_length(sg1);
                    bool is_shower = sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                    sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                    (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11);

                    if (is_shower) {
                        // Check shower significance
                        auto pair_result = calculate_num_daughter_showers(graph, vtx, sg1, false);
                        if (pair_result.second > 35*units::cm) {
                            flag_skip = false;
                            break;
                        }
                    } else {
                        // Check track significance
                        if (!seg_dir_weak(sg1) && length > 6*units::cm) {
                            flag_skip = false;
                            break;
                        }
                    }
                }
                
                if (!flag_skip) {
                    tmp_vertices.insert(vtx);
                } else {
                    // Change particle types to muons for back-to-back tracks
                    double muon_mass = particle_data->get_particle_mass(13);
                    
                    // used_segments is a std::set<SegmentPtr>, i.e. ADDRESS
                    // ordered, and this body mutates: it rewrites particle_info
                    // (pdg/mass), clears kShowerTopology, and calls
                    // change_daughter_type, which propagates the type change
                    // through the graph.  A later iteration then reads a
                    // current_pdg (:2953 below) that an earlier one may have
                    // written, so the outcome depended on heap layout.  Iterate
                    // a copy sorted by graph index instead.
                    //
                    // The SET ITSELF stays pointer-keyed on purpose: its
                    // find() at :2905 must remain IDENTITY, and an index-ordered
                    // set would compare by Segment::get_graph_index(), which is
                    // not guaranteed unique across live Segment objects (doc
                    // pr/28 §11.7, §14.7 -- SIZE_MAX is shared by every segment
                    // not yet added to a graph).  Sorting only at the point of
                    // iteration changes the order and nothing else.  The key is
                    // unique here because every member came from vertex_segments
                    // above, i.e. from live graph edges.
                    std::vector<SegmentPtr> ordered_used_segments(used_segments.begin(),
                                                                  used_segments.end());
                    std::sort(ordered_used_segments.begin(), ordered_used_segments.end(),
                              [](const SegmentPtr& a, const SegmentPtr& b) {
                                  return a->get_graph_index() < b->get_graph_index();
                              });
                    SPDLOG_LOGGER_TRACE(s_log,
                        "improve_vertex: back-to-back retype over {} segment(s) in graph-index order "
                        "(order matters only when this is >1)",
                        ordered_used_segments.size());

                    for (auto sg1 : ordered_used_segments) {
                        // Skip shower trajectory
                        if (sg1->flags_any(SegmentFlags::kShowerTrajectory)) continue;
                        
                        // Handle shower topology
                        if (sg1->flags_any(SegmentFlags::kShowerTopology)) {
                            if (segment_track_length(sg1) > 40*units::cm && sg1->dirsign() == 0) {
                                if (!sg1->particle_info()) {
                                    sg1->particle_info() = std::make_shared<Aux::ParticleInfo>();
                                }
                                sg1->particle_info()->set_pdg(13);
                                sg1->particle_info()->set_mass(muon_mass);
                                pr74_probe_topo_flag(sg1, "unset", "NeutrinoVertexFinder.cxx:examine_main_vertices_local");
                                sg1->unset_flags(SegmentFlags::kShowerTopology);
                                
                                change_daughter_type(graph, vtx, sg1, 13, muon_mass, particle_data, recomb_model);
                                VertexPtr other_vtx = find_other_vertex(graph, sg1, vtx);
                                if (other_vtx) {
                                    change_daughter_type(graph, other_vtx, sg1, 13, muon_mass, particle_data, recomb_model);
                                }
                            } else {
                                continue;
                            }
                        }
                        // Handle non-muon particles
                        else {
                            int current_pdg = sg1->has_particle_info() ? sg1->particle_info()->pdg() : 0;
                            if (current_pdg != 13) {
                                if (!sg1->particle_info()) {
                                    sg1->particle_info() = std::make_shared<Aux::ParticleInfo>();
                                }
                                sg1->particle_info()->set_pdg(13);
                                sg1->particle_info()->set_mass(muon_mass);
                                
                                change_daughter_type(graph, vtx, sg1, 13, muon_mass, particle_data, recomb_model);
                                VertexPtr other_vtx = find_other_vertex(graph, sg1, vtx);
                                if (other_vtx) {
                                    change_daughter_type(graph, other_vtx, sg1, 13, muon_mass, particle_data, recomb_model);
                                }
                            }
                        }
                        
                        // Find continuation muon segments and add final vertex
                        std::vector<VertexPtr> acc_vertices;
                        auto results = find_cont_muon_segment(graph, sg1, vtx);
                        while (results.first != nullptr) {
                            acc_vertices.push_back(results.second);
                            results = find_cont_muon_segment(graph, results.first, results.second);
                        }
                        
                        if (acc_vertices.size() > 0 && acc_vertices.back() != nullptr) {
                            tmp_vertices.insert(acc_vertices.back());
                        }
                    }
                }
            } else {
                // No back-to-back tracks found, keep this vertex
                tmp_vertices.insert(vtx);
            }
        }
    }
    
    // Update vertices collection
    if (tmp_vertices.size() == 0) {
        s_log->trace("examine_main_vertices_local: no vertices survived filtering, returning unchanged");
        return;
    }

    vertices.clear();
    vertices.resize(tmp_vertices.size());
    std::copy(tmp_vertices.begin(), tmp_vertices.end(), vertices.begin());
    s_log->trace("examine_main_vertices_local: done, nvertices_out={}", vertices.size());
}


VertexPtr PatternAlgorithms::compare_main_vertices_global(Graph& graph, std::vector<VertexPtr>& vertex_candidates, Facade::Cluster& main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    if (vertex_candidates.empty()) return nullptr;

    // Sort candidates by cluster_id for deterministic ordering independent of pointer address
    std::sort(vertex_candidates.begin(), vertex_candidates.end(),
              [](const VertexPtr& a, const VertexPtr& b) {
                  int aid = (a && a->cluster()) ? a->cluster()->get_cluster_id() : -1;
                  int bid = (b && b->cluster()) ? b->cluster()->get_cluster_id() : -1;
                  return aid < bid;
              });

    // Initialize scoring map
    std::map<VertexPtr, double> map_vertex_num;
    for (auto vtx : vertex_candidates) {
        map_vertex_num[vtx] = 0;
    }
    
    // Score based on z position (prefer earlier/upstream vertices)
    double min_z = 1e9;
    for (auto vtx : vertex_candidates) {
        WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
        if (vtx_pt.z() < min_z) min_z = vtx_pt.z();
    }
    
    for (auto vtx : vertex_candidates) {
        WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
        // uBooNE 200 cm; m_vertex_z_prior_scale, doc pr/2 2e(iv)
        map_vertex_num[vtx] -= (vtx_pt.z() - min_z) / m_vertex_z_prior_scale;
        
        // Score based on segments connected to this vertex
        if (vtx->descriptor_valid()) {
            auto vd = vtx->get_descriptor();
            const auto edge_range = sorted_out_edges(vd, graph);
            
            for (auto e_it : edge_range) {
                SegmentPtr sg = graph[e_it].segment;
                if (!sg) continue;
                
                bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                sg->flags_any(SegmentFlags::kShowerTopology) ||
                                (sg->has_particle_info() && sg->particle_info() && std::abs(sg->particle_info()->pdg()) == 11);

                if (is_shower) {
                    map_vertex_num[vtx] += 1.0 / 4.0 / 2.0;  // showers count less
                } else {
                    map_vertex_num[vtx] += 1.0 / 4.0;  // tracks
                }
                
                // Bonus for clear protons or tracks with direction
                int pdg = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;
                int dirsign = sg->dirsign();
                bool is_dir_weak = seg_dir_weak(sg);
                
                if (pdg == 2212 && dirsign != 0 && !is_dir_weak) {
                    map_vertex_num[vtx] += 1.0 / 4.0;  // clear proton
                } else if (dirsign != 0 && !is_shower) {
                    map_vertex_num[vtx] += 1.0 / 4.0 / 2.0;  // track with direction
                }
            }
        }
        
        // Bonus if vertex is in main cluster
        if (vtx->cluster() == &main_cluster) {
            map_vertex_num[vtx] += 0.25;
        }
        
        SPDLOG_LOGGER_TRACE(s_log, "compare_main_vertices_global: cluster {} score_A={:.4f} z_norm={:.4f}",
            vtx->cluster() ? vtx->cluster()->get_cluster_id() : -1,
            map_vertex_num[vtx], (vtx_pt.z() - min_z) / m_vertex_z_prior_scale);
    }

    // Score based on fiducial volume
    auto grouping = main_cluster.grouping();
    auto fiducial_utils = grouping ? grouping->get_fiducialutils() : nullptr;
    
    for (auto vtx : vertex_candidates) {
        WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
        
        bool in_fv = fiducial_utils && fiducial_utils->inside_fiducial_volume(vtx_pt);
        if (in_fv || vtx->cluster() == &main_cluster) {
            map_vertex_num[vtx] += 0.5;
        }
        
        SPDLOG_LOGGER_TRACE(s_log, "compare_main_vertices_global: cluster {} score_B={:.4f} in_fv={}",
            vtx->cluster() ? vtx->cluster()->get_cluster_id() : -1,
            map_vertex_num[vtx], in_fv);
    }

    // Calculate direction for each vertex
    std::map<VertexPtr, Facade::geo_vector_t> map_vertex_dir;
    for (auto vtx : vertex_candidates) {
        map_vertex_dir[vtx] = vertex_get_dir(vtx, graph, 5 * units::cm);
    }
    
    // Score based on whether other vertices point toward this vertex
    for (auto vtx : vertex_candidates) {
        WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
        double delta = 0;
        
        for (auto vtx1 : vertex_candidates) {
            if (vtx1 == vtx) continue;
            
            WireCell::Point vtx1_pt = vtx1->fit().valid() ? vtx1->fit().point : vtx1->wcpt().point;
            
            // Direction from vtx to vtx1
            Facade::geo_vector_t dir(vtx1_pt.x() - vtx_pt.x(),
                                     vtx1_pt.y() - vtx_pt.y(),
                                     vtx1_pt.z() - vtx_pt.z());
            
            Facade::geo_vector_t dir1 = map_vertex_dir[vtx1];
            
            if (dir.magnitude() > 0 && dir1.magnitude() > 0) {
                double cos_angle = std::clamp(dir.dot(dir1) / (dir.magnitude() * dir1.magnitude()), -1.0, 1.0);
                double angle = std::acos(cos_angle) / M_PI * 180.0;
                
                if (angle < 15) {
                    map_vertex_num[vtx] += 0.25;
                    delta++;
                } else if (angle < 30) {
                    map_vertex_num[vtx] += 0.25 / 2.0;
                    delta++;
                }
            }
        }
        
        // Penalize isolated vertices not in main cluster
        if (delta == 0) {
            double total_length = 0;
            int num_tracks = 0;
            
            // ordered_edges: total_length is an FP accumulation.
            for (const auto& ed : ordered_edges(graph)) {
                SegmentPtr seg = graph[ed].segment;
                if (!seg || seg->cluster() != vtx->cluster()) continue;

                total_length += segment_track_length(seg);
                num_tracks++;
            }
            
            if (vtx->cluster() != &main_cluster && total_length < 6 * units::cm) {
                map_vertex_num[vtx] -= 0.25 * num_tracks;
            }
        }
        
        SPDLOG_LOGGER_TRACE(s_log, "compare_main_vertices_global: cluster {} score_E={:.4f}",
            vtx->cluster() ? vtx->cluster()->get_cluster_id() : -1,
            map_vertex_num[vtx]);
    }

    // Find vertex with maximum score
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

Facade::Cluster* PatternAlgorithms::check_switch_main_cluster(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    if (!main_cluster) return main_cluster;
    
    bool flag_all_showers = false;
    
    VertexPtr temp_main_vertex = nullptr;
    
    // Check if main cluster has a main vertex
    if (map_cluster_main_vertices.find(main_cluster) != map_cluster_main_vertices.end()) {
        temp_main_vertex = map_cluster_main_vertices[main_cluster];
        int n_showers = 0;
        int n_total = 0;
        
        // Count showers connected to this vertex
        if (temp_main_vertex && temp_main_vertex->descriptor_valid()) {
            auto vd = temp_main_vertex->get_descriptor();
            const auto edge_range = sorted_out_edges(vd, graph);
            
            for (auto e_it : edge_range) {
                SegmentPtr seg = graph[e_it].segment;
                if (!seg) continue;
                
                n_total++;
                bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                seg->flags_any(SegmentFlags::kShowerTopology) ||
                                (seg->has_particle_info() && seg->particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
                if (is_shower) n_showers++;
            }
        }

        if (n_total > 0 && n_showers == n_total) {
            flag_all_showers = true;
        }
    } else {
        flag_all_showers = true;
    }
    
    // If all showers, consider switching main cluster
    if (flag_all_showers) {
        // Collect all vertex candidates, sorted by cluster_id for determinism
        std::vector<VertexPtr> vertex_candidates;
        for (auto& [cluster, vertex] : map_cluster_main_vertices) {
            if (vertex) {
                vertex_candidates.push_back(vertex);
            }
        }
        std::sort(vertex_candidates.begin(), vertex_candidates.end(),
                  [](const VertexPtr& a, const VertexPtr& b) {
                      int aid = (a && a->cluster()) ? a->cluster()->get_cluster_id() : -1;
                      int bid = (b && b->cluster()) ? b->cluster()->get_cluster_id() : -1;
                      return aid < bid;
                  });

        for (auto vtx : vertex_candidates) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            int cluster_id = vtx->cluster() ? vtx->cluster()->get_cluster_id() : -1;
            std::string seg_list;
            if (vtx->descriptor_valid()) {
                auto vd = vtx->get_descriptor();
                const auto edge_range = sorted_out_edges(vd, graph);
                for (auto e_it : edge_range) {
                    SegmentPtr seg = graph[e_it].segment;
                    if (seg) seg_list += std::to_string(seg->id()) + ", ";
                }
            }
            SPDLOG_LOGGER_TRACE(s_log, "check_switch_main_cluster: candidate cluster {} pos=({:.1f},{:.1f},{:.1f}) connecting to: {}",
                cluster_id, vtx_pt.x()/units::cm, vtx_pt.y()/units::cm, vtx_pt.z()/units::cm, seg_list);
        }
        
        // Compare all vertex candidates to find the best one
        VertexPtr temp_main_vertex_1 = nullptr;
        if (!vertex_candidates.empty()) {
            temp_main_vertex_1 = compare_main_vertices_global(graph, vertex_candidates, *main_cluster, track_fitter, dv);
        }
        
        // Check if we should switch
        if (temp_main_vertex_1 && temp_main_vertex_1 != temp_main_vertex) {
            int old_id = temp_main_vertex ? (temp_main_vertex->cluster() ? temp_main_vertex->cluster()->get_cluster_id() : -1) : -1;
            int new_id = temp_main_vertex_1->cluster() ? temp_main_vertex_1->cluster()->get_cluster_id() : -1;
            SPDLOG_LOGGER_TRACE(s_log, "check_switch_main_cluster: switch main cluster {} -> {}", old_id, new_id);
            
            // Find which cluster this vertex belongs to and swap
            for (auto& [cluster, vertex] : map_cluster_main_vertices) {
                if (vertex == temp_main_vertex_1 && cluster != main_cluster) {
                    main_cluster = swap_main_cluster(*cluster, *main_cluster, other_clusters);
                    break;
                }
            }
        }
    }
    
    return main_cluster;
}

Facade::Cluster* PatternAlgorithms::check_switch_main_cluster_2(Graph& graph, VertexPtr temp_main_vertex, Facade::Cluster* max_length_cluster, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters){
    if (!temp_main_vertex || !max_length_cluster || !main_cluster) return main_cluster;
    
    bool flag_switch = false;
    
    // Count showers connected to this vertex
    int n_showers = 0;
    int n_total = 0;
    
    if (temp_main_vertex->descriptor_valid()) {
        auto vd = temp_main_vertex->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);
        
        for (auto e_it : edge_range) {
            SegmentPtr seg = graph[e_it].segment;
            if (!seg) continue;
            
            n_total++;
            bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) ||
                            seg->flags_any(SegmentFlags::kShowerTopology) ||
                            (seg->has_particle_info() && seg->particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
            if (is_shower) n_showers++;
        }
    }
    
    // If all segments are showers, consider switching
    if (n_total > 0 && n_showers == n_total) {
        flag_switch = true;
    }
    
    if (flag_switch) {
        SPDLOG_LOGGER_TRACE(s_log, "check_switch_main_cluster_2: switch main cluster {} -> {}",
            main_cluster->get_cluster_id(), max_length_cluster->get_cluster_id());
        main_cluster = swap_main_cluster(*max_length_cluster, *main_cluster, other_clusters);
    }
    
    return main_cluster;
}

bool PatternAlgorithms::determine_overall_main_vertex_DL(
    Graph& graph,
    ClusterVertexMap& map_cluster_main_vertices,
    Facade::Cluster*& main_cluster,
    std::vector<Facade::Cluster*>& other_clusters,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    TrackFitting& track_fitter,
    IDetectorVolumes::pointer dv,
    const Clus::ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    const std::string& dl_weights,
    double dl_vtx_cut,
    double dQdx_scale,
    double dQdx_offset,
    bool flag_rerank,
    int dl_vtx_top_k,
    double dl_vtx_min_accept_score,
    double dl_vtx_score_scale,
    bool dl_vtx_swap_guard)
{
    bool flag_change = false;

#ifdef HAVE_PYTHON_INC
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    // Collect point cloud (x,y,z,q) from all vertices and segment interior points.
    // Use every vertex/segment in the graph (equivalent to prototype's map_vertex_segments
    // and map_segment_vertices which cover all clusters in the event).
    std::vector<std::vector<float>> vec_xyzq(4);
    std::vector<VertexPtr> cand_vertices;

    // Stable node-index order: cand_vertices drives argmin vertex picks and
    // vec_xyzq is the DL model input point cloud.
    for (const auto& nd : ordered_nodes(graph)) {
        auto vtx = graph[nd].vertex;
        if (!vtx) continue;
        cand_vertices.push_back(vtx);
        auto pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
        vec_xyzq[0].push_back(static_cast<float>(pt.x() / units::cm));
        vec_xyzq[1].push_back(static_cast<float>(pt.y() / units::cm));
        vec_xyzq[2].push_back(static_cast<float>(pt.z() / units::cm));
        double dQ = vtx->fit().valid() ? vtx->fit().dQ : 0.0;
        vec_xyzq[3].push_back(static_cast<float>(dQ * dQdx_scale + dQdx_offset));
    }

    // Segment interior points: skip endpoints (those are vertex positions)
    for (auto e : ordered_edges(graph)) {
        SegmentPtr sg = graph[e].segment;
        if (!sg) continue;
        const auto& fits = sg->fits();
        for (size_t i = 1; i + 1 < fits.size(); ++i) {
            const auto& fit = fits[i];
            vec_xyzq[0].push_back(static_cast<float>(fit.point.x() / units::cm));
            vec_xyzq[1].push_back(static_cast<float>(fit.point.y() / units::cm));
            vec_xyzq[2].push_back(static_cast<float>(fit.point.z() / units::cm));
            vec_xyzq[3].push_back(static_cast<float>(fit.dQ * dQdx_scale + dQdx_offset));
        }
    }
    MS t_collect_pc(Clock::now() - t0); t0 = Clock::now();

    if (vec_xyzq[0].empty()) return false;

    try {
        // Request top_k voxels from Python: legacy mode (flag_rerank==false) uses top_k=1
        // which returns the same 3-float argmax payload as before (byte-for-byte identical).
        // Rerank mode requests dl_vtx_top_k voxels and receives 4*K floats [x,y,z,score per voxel].
        int top_k_arg = flag_rerank ? std::max(1, dl_vtx_top_k) : 1;
        auto dnn_vtx = WCPPyUtil::SCN_Vertex("SCN_Vertex", "SCN_Vertex", dl_weights, vec_xyzq, "float32", false, top_k_arg);
        MS t_scn_inference(Clock::now() - t0); t0 = Clock::now();

        // -----------------------------------------------------------------------
        // Shared output variables: determined by either the legacy or rerank path
        // -----------------------------------------------------------------------
        double min_dis = 1e9;
        VertexPtr min_vertex = nullptr;
        bool flag_pass = false;

        if (!flag_rerank) {
            // ===================================================================
            // LEGACY PATH: top-1 argmax + nearest-neighbor snap + hard gates
            // (identical to the original implementation; behaviour is unchanged)
            // ===================================================================
            if (dnn_vtx.size() != 3) {
                SPDLOG_LOGGER_WARN(s_log, "determine_overall_main_vertex_DL: unexpected DNN output size {}", dnn_vtx.size());
                return false;
            }

            double x_reg = dnn_vtx[0] * units::cm;
            double y_reg = dnn_vtx[1] * units::cm;
            double z_reg = dnn_vtx[2] * units::cm;
            SPDLOG_LOGGER_TRACE(s_log, "determine_overall_main_vertex_DL: DNN prediction: ({:.2f}, {:.2f}, {:.2f}) cm",
                                dnn_vtx[0], dnn_vtx[1], dnn_vtx[2]);

            // Find nearest candidate vertex to DL prediction
            for (auto vtx : cand_vertices) {
                auto pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                double dis = std::sqrt(std::pow(pt.x() - x_reg, 2) +
                                       std::pow(pt.y() - y_reg, 2) +
                                       std::pow(pt.z() - z_reg, 2));
                if (dis < min_dis) {
                    min_dis = dis;
                    min_vertex = vtx;
                }
            }

            if (!min_vertex) return false;

            // Direction sanity check: reject if ALL connected long tracks point away from vertex
            flag_pass = true;
            {
                int num_bad = 0;
                int num_tracks = 0;
                if (min_vertex->descriptor_valid()) {
                    auto vd = min_vertex->get_descriptor();
                    const auto edge_range = sorted_out_edges(vd, graph);
                    for (auto e_it : edge_range) {
                        SegmentPtr sg = graph[e_it].segment;
                        if (!sg) continue;
                        double length = segment_track_length(sg);
                        double medium_dqdx = segment_median_dQ_dx(sg);
                        double dQ_dx_cut = muon_dqdx_cut(length) * m_mip_dqdx_median;
                        auto [v1, v2] = find_vertices(graph, sg);
                        bool flag_start = (v1 == min_vertex);
                        if (length > 15 * units::cm && !seg_dir_weak(sg) && !flag_start && medium_dqdx > dQ_dx_cut) {
                            num_bad++;
                        }
                        num_tracks++;
                    }
                }
                if (num_bad > 0 && num_bad == num_tracks) flag_pass = false;
            }

            // Distance cut
            if (min_dis > dl_vtx_cut) flag_pass = false;

        } else {
            // ===================================================================
            // RERANK PATH: top-K voxels, snap each to nearest ProtoVertex,
            // score by composite heuristic, pick argmax.
            // ===================================================================
            if (dnn_vtx.size() < 4 || dnn_vtx.size() % 4 != 0) {
                SPDLOG_LOGGER_WARN(s_log, "determine_overall_main_vertex_DL: unexpected top-K payload size {}", dnn_vtx.size());
                return false;
            }
            int K = static_cast<int>(dnn_vtx.size() / 4);
            SPDLOG_LOGGER_TRACE(s_log, "determine_overall_main_vertex_DL: rerank mode, K={}", K);

            // --- Parse K (pred_pt, dl_score) from the payload ---
            // dl_score = sigmoid(vertex_class) - sigmoid(bg_class), range [-1, +1].
            // Values near +1: strong vertex signal. Near 0: model is uncertain (both
            // sigmoid outputs ~0.5). Near -1: strong background.  When all top-K scores
            // are small (e.g. <0.01), the DL term barely differentiates candidates and
            // the geometric terms below dominate the final ranking — which is intentional.
            struct DLVoxel { WireCell::Point pred_pt; double dl_score; };
            std::vector<DLVoxel> dl_voxels;
            dl_voxels.reserve(K);
            double dl_score_max = -1e9, dl_score_min = 1e9;
            for (int i = 0; i < K; ++i) {
                WireCell::Point p(dnn_vtx[4*i+0] * units::cm,
                                  dnn_vtx[4*i+1] * units::cm,
                                  dnn_vtx[4*i+2] * units::cm);
                double s = static_cast<double>(dnn_vtx[4*i+3]);
                dl_voxels.push_back({p, s});
                if (s > dl_score_max) dl_score_max = s;
                if (s < dl_score_min) dl_score_min = s;
            }
            SPDLOG_LOGGER_TRACE(s_log,
                "  DL top-{} scores in [{:.4f}, {:.4f}] (scale: [-1,+1]; near 0 = model uncertain)",
                K, dl_score_min, dl_score_max);
            {
                double dl_mean = 0;
                for (const auto& v : dl_voxels) dl_mean += v.dl_score;
                dl_mean /= static_cast<double>(dl_voxels.size());
                double dl_var = 0;
                for (const auto& v : dl_voxels) dl_var += (v.dl_score - dl_mean) * (v.dl_score - dl_mean);
                double dl_std = std::sqrt(dl_var / std::max<size_t>(1, dl_voxels.size() - 1));
                const char* regime = (dl_score_max > 0.1) ? "confident" : (dl_score_max < 0.02 ? "uncertain" : "intermediate");
                SPDLOG_LOGGER_TRACE(s_log, "  DL score stats: mean={:.4f} std={:.4f} regime={}", dl_mean, dl_std, regime);
            }

            // --- Snap each DL voxel to nearest ProtoVertex; deduplicate ---
            // For duplicate snaps (same vtx), keep the entry with the higher dl_score.
            struct SnappedCand { VertexPtr vtx; double snap_dis; double dl_score; int voxel_rank; };
            std::map<VertexPtr, SnappedCand> snap_map;
            for (int vi = 0; vi < static_cast<int>(dl_voxels.size()); ++vi) {
                const auto& vox = dl_voxels[vi];
                double best_dis = 1e9;
                VertexPtr best_vtx = nullptr;
                for (auto vtx : cand_vertices) {
                    auto pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                    double d = std::sqrt(std::pow(pt.x() - vox.pred_pt.x(), 2) +
                                        std::pow(pt.y() - vox.pred_pt.y(), 2) +
                                        std::pow(pt.z() - vox.pred_pt.z(), 2));
                    if (d < best_dis) { best_dis = d; best_vtx = vtx; }
                }
                if (!best_vtx) continue;
                auto it = snap_map.find(best_vtx);
                if (it == snap_map.end() || vox.dl_score > it->second.dl_score)
                    snap_map[best_vtx] = {best_vtx, best_dis, vox.dl_score, vi};
            }
            if (snap_map.empty()) return false;

            std::vector<SnappedCand> snapped;
            snapped.reserve(snap_map.size());
            for (auto& [vtx, sc] : snap_map) snapped.push_back(sc);
            // snap_map is keyed on VertexPtr with the default (ADDRESS) comparator,
            // so the walk above lands in heap order, and the scoring loop below
            // keeps the first candidate reaching best_score (strict `>`) -- an
            // exact tie was decided by allocation layout.  Sort the ITERATION, not
            // the container: snap_map.find() at :3534 must stay address identity
            // (an index-keyed map would compare on get_graph_index(), which is
            // SIZE_MAX for any vertex not yet added to a graph, doc pr/28 sec 14.7).
            std::sort(snapped.begin(), snapped.end(),
                      [](const SnappedCand& a, const SnappedCand& b) {
                          return a.vtx->get_graph_index() < b.vtx->get_graph_index();
                      });

            // --- Precompute cluster total track lengths (avoid O(N_edges) per candidate) ---
            std::map<Facade::Cluster*, double> cluster_total_len;
            {
                // ordered_edges: the += is an FP accumulation.
                for (const auto& ed : ordered_edges(graph)) {
                    SegmentPtr seg = graph[ed].segment;
                    if (seg && seg->cluster())
                        cluster_total_len[seg->cluster()] += segment_track_length(seg);
                }
            }

            // --- Print per-voxel snap summary (voxel → snapped cluster) ---
            // This bridges the DL score list to the geometric scoring below.
            // Geometric terms each contribute ~0.125-0.5; DL term contributes dl_score
            // (range [-1,+1]). When dl_score << 1, geometry dominates — this is expected.
            for (const auto& sc : snapped) {
                double L = 0.0;
                if (sc.vtx->cluster()) {
                    auto it = cluster_total_len.find(sc.vtx->cluster());
                    if (it != cluster_total_len.end()) L = it->second;
                }
                auto pt = sc.vtx->fit().valid() ? sc.vtx->fit().point : sc.vtx->wcpt().point;
                SPDLOG_LOGGER_TRACE(s_log,
                    "  DL voxel {} (dl_score={:.4f}) → cluster={} pos=({:.1f},{:.1f},{:.1f})cm "
                    "L={:.1f}cm snap={:.2f}cm{}",
                    sc.voxel_rank, sc.dl_score,
                    sc.vtx->cluster() ? sc.vtx->cluster()->get_cluster_id() : -1,
                    pt.x()/units::cm, pt.y()/units::cm, pt.z()/units::cm,
                    L/units::cm, sc.snap_dis/units::cm,
                    (sc.vtx->cluster() == main_cluster) ? " [main_cluster]" : "");
            }

            // --- Find min z over snapped candidates for forward-z penalty ---
            double min_z_set = 1e9;
            for (const auto& sc : snapped) {
                auto pt = sc.vtx->fit().valid() ? sc.vtx->fit().point : sc.vtx->wcpt().point;
                if (pt.z() < min_z_set) min_z_set = pt.z();
            }

            // --- Score each candidate ---
            // Empirical weights tuned on 36 annotated events (2026-04-15).
            // Previously scored: segs, ltrk, mult, flg_in, conf, ptbk — removed after empirical
            // analysis showed them disabled or actively harmful. The three active geometric signals
            // (main, clen, isol) dominate when DL is uncertain (scores ~0.005); DL dominates
            // when confident (scores >0.1). fwd_z kept as a vestigial tiebreaker, capped at 0.25.
            constexpr double W_MAIN    = 2.0;   // host cluster == main_cluster
            constexpr double W_CLEN    = 2.0;   // saturating bonus for long host clusters
            constexpr double W_CLEN_L  = 60.0;  // saturation length in cm
            constexpr double W_ISOL    = 2.0;   // short isolated non-main cluster penalty
            constexpr double W_ISOL_L  = 6.0;   // isolation length cutoff in cm
            constexpr double W_FV      = 0.5;   // inside fiducial volume
            constexpr double W_SNAP_L  = 5.0;   // soft snap penalty denominator in cm
            constexpr double W_SNAP_MAX = 2.0;  // soft snap penalty saturation
            constexpr double W_FWD_Z   = 0.25;  // vestigial upstream-z tiebreaker (max |penalty|)
            constexpr double W_FWD_L   = 400.0; // fwd_z normalization in cm

            double best_score = -1e9;
            VertexPtr best_vtx_rerank = nullptr;
            double best_snap_dis = 1e9;

            for (const auto& sc : snapped) {
                VertexPtr vtx = sc.vtx;
                auto vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;

                // doc sbnd_xin/docs/pr/51 (18255-506746): a confident DL
                // voxel (s_dl = raw_score x score_scale, here +576) swamps
                // every +-2 structural term, so the rerank can move the
                // vertex onto a different, non-flash-matched cluster (28 cm
                // off the first-round main cluster).  With the guard on,
                // cross-cluster candidates never enter the acceptance; if
                // none survive, flag_pass stays false and the traditional
                // fallback runs.  Default false => byte-identical.
                if (dl_vtx_swap_guard && vtx->cluster() && main_cluster && vtx->cluster() != main_cluster) {
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "dl_swap_guard: skipping cross-cluster DL candidate voxel {} "
                        "(dl_score={:.4f}) cluster={} != main {} pos=({:.1f},{:.1f},{:.1f})cm snap={:.2f}cm",
                        sc.voxel_rank, sc.dl_score,
                        vtx->cluster()->get_cluster_id(), main_cluster->get_cluster_id(),
                        vtx_pt.x()/units::cm, vtx_pt.y()/units::cm, vtx_pt.z()/units::cm,
                        sc.snap_dis/units::cm);
                    continue;
                }

                // (1) DL confidence — sigmoid-diff in [-1,+1], scaled by dl_vtx_score_scale.
                double s_dl = sc.dl_score * dl_vtx_score_scale;

                // (2) Soft snap penalty: 0 at perfect snap, saturates at -W_SNAP_MAX.
                //     Replaces the old hard "min_dis > 2*dl_vtx_cut" gate, which caused false
                //     rejections on long main-cluster candidates with snap distances of 5-8 cm.
                double s_snap = -std::min(W_SNAP_MAX, sc.snap_dis / (W_SNAP_L * units::cm));

                // (3) Vestigial forward-z tiebreaker: max |penalty| = W_FWD_Z = 0.25.
                //     Capped so it never swamps main/clen bonuses. Breaks upstream/downstream
                //     ties when two candidates land on the same long main cluster.
                double s_fwd_z = -W_FWD_Z * std::clamp(
                    (vtx_pt.z() - min_z_set) / (W_FWD_L * units::cm), 0.0, 1.0);

                // (4) Host cluster total track length: continuous, saturates at W_CLEN_L cm.
                double L_host = 0.0;
                if (vtx->cluster()) {
                    auto it = cluster_total_len.find(vtx->cluster());
                    if (it != cluster_total_len.end()) L_host = it->second;
                }
                double s_clen = W_CLEN * std::min(1.0, L_host / (W_CLEN_L * units::cm));

                // (5) Isolated-short cluster penalty.
                double s_isol = 0.0;
                if (L_host < W_ISOL_L * units::cm && vtx->cluster() != main_cluster)
                    s_isol = -W_ISOL;

                // (6) Main cluster bonus.
                double s_main = (vtx->cluster() == main_cluster) ? W_MAIN : 0.0;

                // (7) Fiducial volume bonus.
                double s_fv = 0.0;
                if (vtx->cluster()) {
                    auto grouping = vtx->cluster()->grouping();
                    auto fv = grouping ? grouping->get_fiducialutils() : nullptr;
                    if (fv && fv->inside_fiducial_volume(vtx_pt)) s_fv = W_FV;
                }

                double score = s_dl + s_snap + s_fwd_z + s_clen + s_isol + s_main + s_fv;

                SPDLOG_LOGGER_TRACE(s_log,
                    "DL rerank cand [voxel {}] cluster={} pos=({:.1f},{:.1f},{:.1f})cm "
                    "L={:.1f}cm snap={:.2f}cm | "
                    "dl={:+.4f} snap={:+.3f} fwd_z={:+.3f} clen={:+.3f} "
                    "isol={:+.3f} main={:+.3f} fv={:+.3f} | TOTAL={:+.3f}",
                    sc.voxel_rank,
                    vtx->cluster() ? vtx->cluster()->get_cluster_id() : -1,
                    vtx_pt.x()/units::cm, vtx_pt.y()/units::cm, vtx_pt.z()/units::cm,
                    L_host/units::cm, sc.snap_dis/units::cm,
                    s_dl, s_snap, s_fwd_z, s_clen, s_isol, s_main, s_fv, score);

                if (score > best_score) {
                    best_score      = score;
                    best_vtx_rerank = vtx;
                    best_snap_dis   = sc.snap_dis;
                }
            }

            min_vertex = best_vtx_rerank;
            min_dis    = best_snap_dis;

            // Accept if the winner clears the composite-score threshold.
            // The soft s_snap penalty already limits distant-snap acceptance;
            // no hard snap gate is needed (and it was causing false rejections
            // on long main-cluster candidates at 5-8 cm snap distance).
            flag_pass = (min_vertex != nullptr)
                     && (best_score >= dl_vtx_min_accept_score);

            if (flag_pass) {
                SPDLOG_LOGGER_TRACE(s_log,
                    "determine_overall_main_vertex_DL: rerank selected cluster={} "
                    "snap_dis={:.2f}cm composite_score={:.4f}",
                    min_vertex->cluster() ? min_vertex->cluster()->get_cluster_id() : -1,
                    min_dis / units::cm, best_score);
            } else {
                SPDLOG_LOGGER_TRACE(s_log,
                    "determine_overall_main_vertex_DL: rerank rejected "
                    "(best_score={:.4f} < threshold={:.4f}), staying with traditional vertex",
                    best_score, dl_vtx_min_accept_score);
            }
        }

        MS t_selection(Clock::now() - t0); t0 = Clock::now();

        MS t_examine_direction(MS::zero());
        MS t_proton_tagging(MS::zero());
        MS t_cleanup_long_muon(MS::zero());

        // doc pr/48: a two-end-break junction vertex carries direct dQ/dx
        // evidence -- two stopping-template arms, the very signature the
        // break was accepted on -- that the image-based DL vertex cannot
        // see; the back-to-back mid-track junction is precisely the
        // topology DL and angular methods miss (that is why the break pass
        // exists).  Never let the DL rerank move the main vertex OFF such a
        // vertex (51513/56211/57485: the traditional chain selects the
        // junction and the DL snap lands on a Bragg tip).  A DL choice that
        // AGREES with the protected vertex still passes.  Inert unless
        // m_two_end_break (the flag is only ever set by that pass).
        if (flag_pass && m_two_end_break && main_cluster) {
            auto trad_it = map_cluster_main_vertices.find(main_cluster);
            if (trad_it != map_cluster_main_vertices.end() && trad_it->second &&
                trad_it->second->flags_any(VertexFlags::kProtectedBreak) &&
                min_vertex != trad_it->second) {
                SPDLOG_LOGGER_TRACE(s_log,
                    "determine_overall_main_vertex_DL: keeping protected two-end-break "
                    "vertex over the DL choice (doc pr/48)");
                flag_pass = false;
            }
        }

        if (flag_pass) {
            flag_change = true;
            SPDLOG_LOGGER_TRACE(s_log,
                "determine_overall_main_vertex_DL: switching to DL vertex (dis={:.2f} cm)",
                min_dis / units::cm);

            // Switch main_cluster if DL vertex belongs to a different cluster
            if (main_cluster && min_vertex->cluster() && min_vertex->cluster() != main_cluster) {
                main_cluster = swap_main_cluster(*min_vertex->cluster(), *main_cluster, other_clusters);
            }

            // Record DL-chosen vertex as the neutrino vertex for the main cluster
            map_cluster_main_vertices[main_cluster] = min_vertex;

            VertexPtr main_vertex = min_vertex;

            // Re-examine track directions with new main vertex as reference
            examine_direction(graph, main_vertex, main_vertex,
                              vertices_in_long_muon, segments_in_long_muon,
                              particle_data, recomb_model, false);
            t_examine_direction = MS(Clock::now() - t0); t0 = Clock::now();

            // Proton tagging: short stubs with high dQ/dx near the new vertex are likely protons
            if (main_vertex->descriptor_valid()) {
                auto vd = main_vertex->get_descriptor();
                const auto edge_range = sorted_out_edges(vd, graph);
                for (auto e_it : edge_range) {
                    SegmentPtr sg = graph[e_it].segment;
                    if (!sg) continue;
                    auto pair_results = calculate_num_daughter_showers(graph, main_vertex, sg, false);
                    double length = segment_track_length(sg);
                    double median_dqdx = segment_median_dQ_dx(sg) / m_mip_dqdx_median;
                    if (pair_results.first == 1 && length < 1.5 * units::cm && median_dqdx > 1.6) {
                        auto four_momentum = segment_cal_4mom(sg, 2212, particle_data, recomb_model, m_mip_dqdx);
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(
                            2212, particle_data->get_particle_mass(2212),
                            particle_data->pdg_to_name(2212), four_momentum);
                        sg->particle_info(pinfo);
                    }
                }
            }
            t_proton_tagging = MS(Clock::now() - t0); t0 = Clock::now();

            // Long muon cleanup: remove segments/vertices that no longer belong to main cluster
            {
                std::set<SegmentPtr> tmp_segments;
                std::set<VertexPtr>  tmp_vertices;
                for (auto seg : segments_in_long_muon) {
                    if (seg && seg->cluster() != main_cluster) tmp_segments.insert(seg);
                }
                for (auto vtx : vertices_in_long_muon) {
                    if (vtx && vtx->cluster() != main_cluster) tmp_vertices.insert(vtx);
                }
                for (auto seg : tmp_segments) segments_in_long_muon.erase(seg);
                for (auto vtx : tmp_vertices) vertices_in_long_muon.erase(vtx);
            }
            t_cleanup_long_muon = MS(Clock::now() - t0);
        } else {
            SPDLOG_LOGGER_TRACE(s_log, "determine_overall_main_vertex_DL: staying with traditional vertex");
        }

        if (m_perf) {
            MS t_total_ms(Clock::now() - t_total);
            SPDLOG_LOGGER_TRACE(s_log,
                "determine_overall_main_vertex_DL timing: "
                "collect_pc={:.3f}ms scn_inference={:.3f}ms selection={:.3f}ms "
                "examine_direction={:.3f}ms "
                "proton_tagging={:.3f}ms cleanup_long_muon={:.3f}ms TOTAL={:.3f}ms",
                t_collect_pc.count(), t_scn_inference.count(), t_selection.count(),
                t_examine_direction.count(),
                t_proton_tagging.count(), t_cleanup_long_muon.count(), t_total_ms.count());
        }
    }
    catch (const std::exception& ex) {
        SPDLOG_LOGGER_WARN(s_log, "determine_overall_main_vertex_DL: DL vertex failed: {}", ex.what());
    }
#endif  // HAVE_PYTHON_INC

    return flag_change;
}

VertexPtr PatternAlgorithms::determine_overall_main_vertex(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_dev_chain){
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    if (!main_cluster) return nullptr;

    // Find cluster with maximum length.
    // Collect all unique candidates, sort by cluster_id for determinism, then scan.
    Facade::Cluster* max_length_cluster = nullptr;
    double max_length = 0;

    std::vector<Facade::Cluster*> all_candidate_clusters;
    for (auto& [cluster, vertex] : map_cluster_main_vertices) {
        if (cluster) all_candidate_clusters.push_back(cluster);
    }
    all_candidate_clusters.push_back(main_cluster);
    for (auto cluster : other_clusters) {
        if (cluster) all_candidate_clusters.push_back(cluster);
    }
    // Remove duplicates while preserving deterministic order
    std::sort(all_candidate_clusters.begin(), all_candidate_clusters.end(),
              [](Facade::Cluster* a, Facade::Cluster* b) {
                  return a->get_cluster_id() < b->get_cluster_id();
              });
    all_candidate_clusters.erase(std::unique(all_candidate_clusters.begin(),
                                             all_candidate_clusters.end()),
                                 all_candidate_clusters.end());
    for (auto cluster : all_candidate_clusters) {
        double length = cluster->get_length();
        if (length > max_length) {
            max_length = length;
            max_length_cluster = cluster;
        }
    }
    MS t_find_max_length(Clock::now() - t0); t0 = Clock::now();

    // Examine main vertices first
    examine_main_vertices(graph, map_cluster_main_vertices, main_cluster, other_clusters);
    MS t_examine_main_vertices(Clock::now() - t0); t0 = Clock::now();

    // Check for main cluster switch
    if (flag_dev_chain) {
        // Development chain: use compare_main_vertices_global to find the globally best vertex
        main_cluster = check_switch_main_cluster(graph, map_cluster_main_vertices, main_cluster, other_clusters,
                                                 track_fitter, dv);
    } else {
        // Frozen chain: only switch if a *different* cluster is significantly longer and main vertex is all showers
        if (max_length_cluster && main_cluster && max_length_cluster != main_cluster) {
            double main_length = main_cluster->get_length();
            if (max_length > main_length * 0.8) {
                VertexPtr temp_main_vertex = map_cluster_main_vertices.find(main_cluster) != map_cluster_main_vertices.end()
                                             ? map_cluster_main_vertices[main_cluster] : nullptr;
                if (temp_main_vertex) {
                    main_cluster = check_switch_main_cluster_2(graph, temp_main_vertex, max_length_cluster, main_cluster,
                                                               other_clusters);
                }
            }
        }
    }
    MS t_check_switch(Clock::now() - t0); t0 = Clock::now();

    // Get the main vertex
    VertexPtr main_vertex = nullptr;
    if (map_cluster_main_vertices.find(main_cluster) != map_cluster_main_vertices.end()) {
        main_vertex = map_cluster_main_vertices[main_cluster];
    }

    if (!main_vertex) {
        if (m_perf) {
            MS t_total_ms(Clock::now() - t_total);
            SPDLOG_LOGGER_TRACE(s_log,
                "determine_overall_main_vertex timing (early return, no main vertex): "
                "find_max_length={:.3f}ms examine_main_vertices={:.3f}ms "
                "check_switch={:.3f}ms TOTAL={:.3f}ms",
                t_find_max_length.count(), t_examine_main_vertices.count(),
                t_check_switch.count(), t_total_ms.count());
        }
        return nullptr;
    }

    // Examine tracks connected to main vertex - look for short high dQ/dx proton candidates
    if (main_vertex->descriptor_valid()) {
        auto vd = main_vertex->get_descriptor();
        const auto edge_range = sorted_out_edges(vd, graph);

        for (auto e_it : edge_range) {
            SegmentPtr sg = graph[e_it].segment;
            if (!sg) continue;

            auto pair_results = calculate_num_daughter_showers(graph, main_vertex, sg, false);
            double length = segment_track_length(sg);
            double median_dqdx = segment_median_dQ_dx(sg) / m_mip_dqdx_median;

            // Short segment with only 1 daughter and high dQ/dx -> likely proton
            if (pair_results.first == 1 && length < 1.5 * units::cm && median_dqdx > 1.6) {
                auto four_momentum = segment_cal_4mom(sg, 2212, particle_data, recomb_model, m_mip_dqdx);
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    2212, particle_data->get_particle_mass(2212),
                    particle_data->pdg_to_name(2212), four_momentum);
                sg->particle_info(pinfo);
            }
        }
    }
    MS t_proton_tagging(Clock::now() - t0); t0 = Clock::now();

    // Clean up long muons - remove segments/vertices not in main cluster
    {
        std::set<SegmentPtr> tmp_segments;
        std::set<VertexPtr> tmp_vertices;

        // Find segments not in main cluster
        for (auto seg : segments_in_long_muon) {
            if (seg && seg->cluster() != main_cluster) {
                tmp_segments.insert(seg);
            }
        }

        // Find vertices not in main cluster
        for (auto vtx : vertices_in_long_muon) {
            if (vtx && vtx->cluster() != main_cluster) {
                tmp_vertices.insert(vtx);
            }
        }

        // Remove them from the long muon sets
        for (auto seg : tmp_segments) {
            segments_in_long_muon.erase(seg);
        }

        for (auto vtx : tmp_vertices) {
            vertices_in_long_muon.erase(vtx);
        }
    }
    MS t_cleanup_long_muon(Clock::now() - t0);

    if (m_perf) {
        MS t_total_ms(Clock::now() - t_total);
        SPDLOG_LOGGER_TRACE(s_log,
            "determine_overall_main_vertex timing: "
            "find_max_length={:.3f}ms examine_main_vertices={:.3f}ms "
            "check_switch={:.3f}ms proton_tagging={:.3f}ms "
            "cleanup_long_muon={:.3f}ms TOTAL={:.3f}ms",
            t_find_max_length.count(), t_examine_main_vertices.count(),
            t_check_switch.count(), t_proton_tagging.count(),
            t_cleanup_long_muon.count(), t_total_ms.count());
    }

    return main_vertex;
}