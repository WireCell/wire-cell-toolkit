// NeutrinoTaggerSinglePhoton.cxx  — Subsection 1 of 4
//
// Ports from prototype:
//   NeutrinoID_singlephoton_tagger.h
//     bad_reconstruction_1_sp  (line 4170) → fills shw_sp_br2_* TaggerInfo fields
//     bad_reconstruction_sp    (line 3766) → fills shw_sp_br1_* TaggerInfo fields
//
// Namespace/class: WireCell::Clus::PR::PatternAlgorithms
//
// Design: all helpers are static file-local functions taking SpContext& ctx.
// The public entry point PatternAlgorithms::singlephoton_tagger() will be added
// in Subsection 4 and declared in NeutrinoPatternBase.h.
//
// Translation conventions (see neutrino_id_function_map.md):
//   map_vertex_segments[vtx]              → boost::out_edges / vtx_degree
//   map_vertex_segments[v].size()         → vtx_degree(v, ctx.graph)
//   sg->get_length()                      → segment_track_length(sg)
//   sg->get_direct_length()               → segment_track_direct_length(sg)
//   sg->get_medium_dQ_dx()               → segment_median_dQ_dx(sg)
//   sg->cal_dir_3vector(pt, dis)          → segment_cal_dir_3vector(sg, pt, dis)
//   shower->cal_dir_3vector(pt, dis)      → shower_cal_dir_3vector(*shower, pt, dis)
//   shower->get_start_vertex().first      → shower->get_start_vertex_and_type().first
//   shower->get_start_segment()           → shower->start_segment()
//   shower->fill_point_vec(pts, true)     → shower->fill_point_vector(pts, true)
//   sg->get_flag_shower_topology()        → sg->flags_any(SegmentFlags::kShowerTopology)
//   sg->get_flag_shower_trajectory()      → sg->flags_any(SegmentFlags::kShowerTrajectory)
//   sg->get_flag_avoid_muon_check()       → sg->flags_any(SegmentFlags::kAvoidMuonCheck)
//   find_other_vertex(sg, v)              → find_other_vertex(ctx.graph, sg, v)
//   find_vertices(sg)                     → find_vertices(ctx.graph, sg)
//   find_cont_muon_segment_nue(sg,v,f)    → find_cont_muon_segment_nue(ctx.graph, sg, v, f)
//   main_cluster->Calc_PCA(pts)
//     + get_PCA_axis(0)                   → ctx.self.calc_PCA_main_axis(pts).second
//   vertex wcpt-index front/back check   → seg_endpoint_near(sg, vtx_fit_pt(vtx))
//   tagger_info.X (flag_fill-gated)       → ti.X (unconditional)
//   flag_print output                     → dropped

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;
using namespace WireCell;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// ---------------------------------------------------------------------------
// File-local helpers
// ---------------------------------------------------------------------------

// Vertex fit point (with wcpt fallback).
static inline Point vtx_fit_pt(VertexPtr v) {
    return v->fit().valid() ? v->fit().point : v->wcpt().point;
}

// Number of graph edges (segments) at a vertex.
// Replaces prototype's map_vertex_segments[vtx].size().
static inline int vtx_degree(VertexPtr vtx, const Graph& graph) {
    if (!vtx || !vtx->descriptor_valid()) return 0;
    return static_cast<int>(boost::out_degree(vtx->get_descriptor(), graph));
}

// Best-estimate shower energy: kine_best if set, else kine_charge.
static inline double shower_energy(ShowerPtr shower) {
    return (shower->get_kine_best() != 0) ? shower->get_kine_best()
                                           : shower->get_kine_charge();
}

// True when segment is shower-like (trajectory, topology, or electron PDG).
// Replaces prototype's sg->get_flag_shower().
static inline bool seg_is_shower(SegmentPtr seg) {
    return seg->flags_any(SegmentFlags::kShowerTrajectory) ||
           seg->flags_any(SegmentFlags::kShowerTopology)   ||
           (seg->has_particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
}

// Endpoint of seg that is geometrically nearest to ref_pt.
// Replaces prototype's wcpt-index comparison for front/back of segment.
static Point seg_endpoint_near(SegmentPtr seg, const Point& ref_pt) {
    const auto& fits = seg->fits();
    Point front = fits.front().point;
    Point back  = fits.back().point;
    return (ray_length(Ray{ref_pt, front}) <= ray_length(Ray{ref_pt, back}))
           ? front : back;
}

// ---------------------------------------------------------------------------
// SpContext: file-local bundle of shared state for singlephoton_tagger helpers.
//
// Constructed once inside PatternAlgorithms::singlephoton_tagger() (Subsection 4)
// and passed by reference to every helper.
// ---------------------------------------------------------------------------
struct SpContext {
    PatternAlgorithms& self;                        // for member functions (calc_PCA_main_axis, etc.)
    Graph& graph;
    Facade::Cluster* main_cluster;                  // needed by bad_reconstruction_1_sp (PCA)
                                                    // and low_energy_michel_sp (check_direction)
    VertexPtr main_vertex;
    int apa{0}, face{0};                            // for SCE correction and point-cloud queries
    IndexedShowerSet& showers;
    VertexShowerSetMap& map_vertex_to_shower;
    ShowerIntMap& map_shower_pio_id;
    std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers;
    std::map<int, std::pair<double,int>>& map_pio_id_mass;
    IDetectorVolumes::pointer dv;
    IClusGeomHelper::pointer geom_helper;           // nullable; for SCE correction in entry point
};

// ===========================================================================
// bad_reconstruction_sp
//
// Determines whether a shower is a "bad reconstruction" (i.e. really a long
// muon-like track mis-identified as a shower).  Three sub-checks:
//   br1_1 : stem-length / topology check
//   br1_2 : longest muon-like track inside the shower (via find_cont_muon_segment_nue)
//   br1_3 : straight track near the far end of the start segment
//
// Fills TaggerInfo shw_sp_br1_* fields (unconditionally).
//
// Prototype: WCPPID::NeutrinoID::bad_reconstruction_sp()
//            NeutrinoID_singlephoton_tagger.h line 3766.
//
// The logic is identical to bad_reconstruction() in NeutrinoTaggerCosmic.cxx;
// only the TaggerInfo field names differ (shw_sp_br1_* vs br1_*).
// ===========================================================================
static bool bad_reconstruction_sp(SpContext& ctx, ShowerPtr shower, TaggerInfo& ti)
{
    bool flag_bad_shower_1 = false;
    bool flag_bad_shower_2 = false;
    bool flag_bad_shower_3 = false;

    double Eshower = shower_energy(shower);

    SegmentPtr sg = shower->start_segment();
    if (!sg) return false;

    auto [vtx, start_type] = shower->get_start_vertex_and_type();

    // Collect all segments/vertices in the shower.
    // Replaces prototype's shower->get_map_seg_vtxs() / get_map_vtx_segs().
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    // -------------------------------------------------------------------
    // Sub-check 1 (br1_1): stem characteristics.
    // Prototype lines 3789-3808.
    // -------------------------------------------------------------------
    {
        double sg_length = segment_track_length(sg);

        if (start_type == 1 && vtx_degree(vtx, ctx.graph) == 1 &&
            Eshower < 120 * units::MeV && (int)shower_segs.size() <= 3) {
            bool topo = sg->flags_any(SegmentFlags::kShowerTopology);
            bool traj = sg->flags_any(SegmentFlags::kShowerTrajectory);
            if (!topo && !traj && sg_length > 10 * units::cm)
                flag_bad_shower_1 = true;
        }
        if (sg_length > 80 * units::cm)  // stem too long → definitely bad
            flag_bad_shower_1 = true;

        ti.shw_sp_br1_1_flag              = !flag_bad_shower_1;
        ti.shw_sp_br1_1_shower_type       = start_type;
        ti.shw_sp_br1_1_vtx_n_segs        = vtx_degree(vtx, ctx.graph);
        ti.shw_sp_br1_1_energy            = static_cast<float>(Eshower / units::MeV);
        ti.shw_sp_br1_1_n_segs            = static_cast<float>(shower_segs.size());
        ti.shw_sp_br1_1_flag_sg_topology  = sg->flags_any(SegmentFlags::kShowerTopology);
        ti.shw_sp_br1_1_flag_sg_trajectory= sg->flags_any(SegmentFlags::kShowerTrajectory);
        ti.shw_sp_br1_1_sg_length         = static_cast<float>(sg_length / units::cm);
    }

    // -------------------------------------------------------------------
    // Sub-check 2 (br1_2): look for a long muon-like track inside the shower.
    // For each shower segment, try to extend it via find_cont_muon_segment_nue.
    // Prototype lines 3813-3974.
    // -------------------------------------------------------------------
    {
        double max_length      = 0;
        int    n_connected     = 0;
        int    n_connected1    = 0;
        double max_length_ratio= 0;

        for (SegmentPtr sg1 : shower_segs) {
            double length        = segment_track_length(sg1);
            double direct_length = segment_track_direct_length(sg1);
            bool   topo          = sg1->flags_any(SegmentFlags::kShowerTopology);
            bool   avoid_muon    = sg1->flags_any(SegmentFlags::kAvoidMuonCheck);

            if (avoid_muon) continue;
            if (topo && direct_length <= 0.94 * length) continue;

            auto [sv1, sv2] = find_vertices(ctx.graph, sg1);
            if (!sv1 || !sv2) continue;

            double tmp_length = length;
            int    tmp_nc1    = 0;

            if (sv1 != ctx.main_vertex) {
                auto [ext_sg, ext_vtx] = ctx.self.find_cont_muon_segment_nue(ctx.graph, sg1, sv1, true);
                if (ext_sg) {
                    bool ext_topo  = ext_sg->flags_any(SegmentFlags::kShowerTopology);
                    bool ext_avoid = ext_sg->flags_any(SegmentFlags::kAvoidMuonCheck);
                    double ext_dl  = segment_track_direct_length(ext_sg);
                    double ext_len = segment_track_length(ext_sg);
                    if (!ext_avoid && (!ext_topo || ext_dl > 0.94 * ext_len)) {
                        tmp_length += ext_len;
                        tmp_nc1    += vtx_degree(ext_vtx, ctx.graph) - 1;
                    }
                }
            }
            if (sv2 != ctx.main_vertex) {
                auto [ext_sg, ext_vtx] = ctx.self.find_cont_muon_segment_nue(ctx.graph, sg1, sv2, true);
                if (ext_sg) {
                    bool ext_topo  = ext_sg->flags_any(SegmentFlags::kShowerTopology);
                    bool ext_avoid = ext_sg->flags_any(SegmentFlags::kAvoidMuonCheck);
                    double ext_dl  = segment_track_direct_length(ext_sg);
                    double ext_len = segment_track_length(ext_sg);
                    if (!ext_avoid && (!ext_topo || ext_dl > 0.94 * ext_len)) {
                        tmp_length += ext_len;
                        tmp_nc1    += vtx_degree(ext_vtx, ctx.graph) - 1;
                    }
                }
            }

            // 6cm offset for topology segments or segments outside the start cluster
            double length_offset = 0;
            int start_cl = sg->cluster()  ? sg->cluster()->get_cluster_id()  : -1;
            int sg1_cl   = sg1->cluster() ? sg1->cluster()->get_cluster_id() : -1;
            if (topo || sg1_cl != start_cl) length_offset = 6 * units::cm;

            double eff_length = tmp_length - length_offset;
            if (eff_length > max_length) {
                max_length        = eff_length;
                max_length_ratio  = (length > 0) ? direct_length / length : 0;
                n_connected1      = tmp_nc1;
                n_connected       = 0;
                if (sv1 != ctx.main_vertex) n_connected += vtx_degree(sv1, ctx.graph) - 1;
                if (sv2 != ctx.main_vertex) n_connected += vtx_degree(sv2, ctx.graph) - 1;
            }
        }

        auto check_len2 = [&](int nc, double ml, double t0, double t1, double t2, double t3) {
            if (nc <= 1 && ml > t0) return true;
            if (nc == 2 && ml > t1) return true;
            if (nc == 3 && ml > t2) return true;
            if (ml > t3)            return true;
            return false;
        };

        if (Eshower < 200 * units::MeV) {
            flag_bad_shower_2 = check_len2(n_connected, max_length,
                38*units::cm, 42*units::cm, 46*units::cm, 50*units::cm);
        } else if (Eshower < 400 * units::MeV) {
            flag_bad_shower_2 = check_len2(n_connected, max_length,
                42*units::cm, 49*units::cm, 52*units::cm, 55*units::cm);
            if (n_connected + n_connected1 > 4 && max_length <= 72 * units::cm)
                flag_bad_shower_2 = false;
        } else if (Eshower < 600 * units::MeV) {
            flag_bad_shower_2 = check_len2(n_connected, max_length,
                45*units::cm, 48*units::cm, 54*units::cm, 62*units::cm);
        } else if (Eshower < 800 * units::MeV) {
            flag_bad_shower_2 = check_len2(n_connected, max_length,
                51*units::cm, 52*units::cm, 56*units::cm, 62*units::cm);
            if (flag_bad_shower_2) {
                if ((vtx_degree(ctx.main_vertex, ctx.graph) == 1 && max_length < 68 * units::cm) ||
                    (n_connected >= 6 && max_length < 76 * units::cm))
                    flag_bad_shower_2 = false;
            }
            if (shower->get_num_segments() >= 15 && max_length < 60 * units::cm)
                flag_bad_shower_2 = false;
        } else if (Eshower < 1500 * units::MeV) {
            flag_bad_shower_2 = check_len2(n_connected, max_length,
                55*units::cm, 60*units::cm, 65*units::cm, 75*units::cm);
        } else {
            flag_bad_shower_2 = check_len2(n_connected, max_length,
                55*units::cm, 65*units::cm, 70*units::cm, 75*units::cm);
        }

        if (Eshower > 1000 * units::MeV && flag_bad_shower_2 && max_length_ratio < 0.95)
            flag_bad_shower_2 = false;

        double total_len = shower->get_total_length();
        if (max_length > 0.75 * total_len && max_length > 35 * units::cm)
            flag_bad_shower_2 = true;

        ti.shw_sp_br1_2_flag            = !flag_bad_shower_2;
        ti.shw_sp_br1_2_energy          = static_cast<float>(Eshower / units::MeV);
        ti.shw_sp_br1_2_n_connected     = n_connected;
        ti.shw_sp_br1_2_max_length      = static_cast<float>(max_length / units::cm);
        ti.shw_sp_br1_2_n_connected_1   = n_connected1;
        ti.shw_sp_br1_2_vtx_n_segs      = vtx_degree(ctx.main_vertex, ctx.graph);
        ti.shw_sp_br1_2_n_shower_segs   = shower->get_num_segments();
        ti.shw_sp_br1_2_max_length_ratio= static_cast<float>(max_length_ratio);
        ti.shw_sp_br1_2_shower_length   = static_cast<float>(total_len / units::cm);
    }

    // -------------------------------------------------------------------
    // Sub-check 3 (br1_3): long straight track near the far end of the start
    // segment ("main length" test).
    // Prototype lines 3976-4161.
    // -------------------------------------------------------------------
    {
        double max_length  = 0;
        int    n_connected = 0;
        double main_length = segment_track_length(sg);

        VertexPtr other_vtx = find_other_vertex(ctx.graph, sg, vtx);
        Point     other_pt  = vtx_fit_pt(other_vtx);

        if (main_length > 10 * units::cm && other_vtx) {
            Vector dir1 = segment_cal_dir_3vector(sg, other_pt, 15 * units::cm);

            for (SegmentPtr sg1 : shower_segs) {
                if (sg1 == sg) continue;
                bool topo  = sg1->flags_any(SegmentFlags::kShowerTopology);
                bool traj  = sg1->flags_any(SegmentFlags::kShowerTrajectory);
                double sg1_len = segment_track_length(sg1);
                if (topo || traj || sg1_len < 10 * units::cm) continue;

                auto [pv1, pv2] = find_vertices(ctx.graph, sg1);
                if (!pv1 || !pv2) continue;

                Point  pt1  = vtx_fit_pt(pv1);
                Point  pt2  = vtx_fit_pt(pv2);
                double dis1 = ray_length(Ray{pt1, other_pt});
                double dis2 = ray_length(Ray{pt2, other_pt});

                double tmp_length1 = 0;
                int    tmp_nc      = 0;

                if (dis1 < 5 * units::cm) {
                    Vector dir2  = segment_cal_dir_3vector(sg1, pt1, 15 * units::cm);
                    double angle = dir1.angle(dir2) / M_PI * 180.0;
                    if (angle > 170) {
                        if (sg1_len + dis1 > tmp_length1) {
                            tmp_length1 = sg1_len + dis1;
                            tmp_nc = vtx_degree(pv2, ctx.graph) - 1;
                        }
                    }
                } else if (dis2 < 5 * units::cm) {
                    Vector dir2  = segment_cal_dir_3vector(sg1, pt2, 15 * units::cm);
                    double angle = dir1.angle(dir2) / M_PI * 180.0;
                    if (angle > 170) {
                        if (sg1_len + dis2 > tmp_length1) {
                            tmp_length1 = sg1_len + dis2;
                            tmp_nc = vtx_degree(pv1, ctx.graph) - 1;
                        }
                    }
                } else {
                    Point  close_pt;
                    Vector dir2;
                    if (dis1 < dis2) {
                        dir2     = segment_cal_dir_3vector(sg1, pt1, 15 * units::cm);
                        close_pt = pt1;
                    } else {
                        dir2     = segment_cal_dir_3vector(sg1, pt2, 15 * units::cm);
                        close_pt = pt2;
                    }
                    double angle = dir1.angle(dir2) / M_PI * 180.0;
                    if (angle > 165) {
                        for (SegmentPtr sg2 : shower_segs) {
                            if (sg2 == sg || sg2 == sg1) continue;
                            double sg2_len = segment_track_length(sg2);
                            if (sg2_len < 10 * units::cm) continue;

                            auto [pv1_2, pv2_2] = find_vertices(ctx.graph, sg2);
                            if (!pv1_2 || !pv2_2) continue;

                            double d3 = ray_length(Ray{vtx_fit_pt(pv1_2), other_pt});
                            double d4 = ray_length(Ray{vtx_fit_pt(pv2_2), other_pt});

                            double angle1 = 0;
                            if (d3 < 6 * units::cm) {
                                Point pt1_2 = vtx_fit_pt(pv1_2);
                                Vector dcheck = segment_cal_dir_3vector(sg2, pt1_2, 15 * units::cm);
                                angle1 = dir1.angle(dcheck) / M_PI * 180.0;
                            } else if (d4 < 6 * units::cm) {
                                Point pt2_2 = vtx_fit_pt(pv2_2);
                                Vector dcheck = segment_cal_dir_3vector(sg2, pt2_2, 15 * units::cm);
                                angle1 = dir1.angle(dcheck) / M_PI * 180.0;
                            }
                            if (angle1 > 170 && sg2_len > 0.75 * std::min(dis1, dis2)) {
                                if (sg1_len + sg2_len > tmp_length1) {
                                    tmp_length1 = sg1_len + sg2_len;
                                    tmp_nc = (dis2 < dis1) ? vtx_degree(pv1, ctx.graph) - 1
                                                           : vtx_degree(pv2, ctx.graph) - 1;
                                }
                            }
                        }
                    }
                }

                if (tmp_length1 + main_length > max_length) {
                    max_length  = tmp_length1 + main_length;
                    n_connected = tmp_nc;
                }
            }
        }

        auto check_len3 = [&](int nc, double ml, double t0, double t1, double t2, double t3) {
            if (nc <= 1 && ml > t0) return true;
            if (nc == 2 && ml > t1) return true;
            if (nc == 3 && ml > t2) return true;
            if (ml > t3)            return true;
            return false;
        };

        if (Eshower < 200 * units::MeV) {
            flag_bad_shower_3 = check_len3(n_connected, max_length,
                36*units::cm, 42*units::cm, 48*units::cm, 54*units::cm);
        } else if (Eshower < 400 * units::MeV) {
            flag_bad_shower_3 = check_len3(n_connected, max_length,
                45*units::cm, 42*units::cm, 42*units::cm, 50*units::cm);
        } else if (Eshower < 800 * units::MeV) {
            flag_bad_shower_3 = check_len3(n_connected, max_length,
                55*units::cm, 60*units::cm, 75*units::cm, 80*units::cm);
            if (shower->get_num_segments() > 20 && max_length < 90 * units::cm)
                flag_bad_shower_3 = false;
        } else if (Eshower < 1500 * units::MeV) {
            flag_bad_shower_3 = check_len3(n_connected, max_length,
                55*units::cm, 60*units::cm, 75*units::cm, 80*units::cm);
        } else {
            flag_bad_shower_3 = check_len3(n_connected, max_length,
                50*units::cm, 60*units::cm, 75*units::cm, 80*units::cm);
        }

        if (flag_bad_shower_3) {
            bool sg_topo = sg->flags_any(SegmentFlags::kShowerTopology);
            bool sg_traj = sg->flags_any(SegmentFlags::kShowerTrajectory);
            if ((!sg_topo || (sg_topo && Eshower < 200 * units::MeV)) &&
                !sg_traj && shower->get_num_main_segments() == 1) {
                if (max_length <= segment_track_length(sg))
                    flag_bad_shower_3 = false;
            } else {
                flag_bad_shower_3 = false;
            }
        }

        ti.shw_sp_br1_3_flag               = !flag_bad_shower_3;
        ti.shw_sp_br1_3_energy             = static_cast<float>(Eshower / units::MeV);
        ti.shw_sp_br1_3_n_connected_p      = n_connected;
        ti.shw_sp_br1_3_max_length_p       = static_cast<float>(max_length / units::cm);
        ti.shw_sp_br1_3_n_shower_segs      = shower->get_num_segments();
        ti.shw_sp_br1_3_flag_sg_topology   = sg->flags_any(SegmentFlags::kShowerTopology);
        ti.shw_sp_br1_3_flag_sg_trajectory = sg->flags_any(SegmentFlags::kShowerTrajectory);
        ti.shw_sp_br1_3_n_shower_main_segs = shower->get_num_main_segments();
        ti.shw_sp_br1_3_sg_length          = static_cast<float>(segment_track_length(sg) / units::cm);
    }

    bool flag_bad = flag_bad_shower_1 || flag_bad_shower_2 || flag_bad_shower_3;
    ti.shw_sp_br1_flag = !flag_bad;
    return flag_bad;
}

// ===========================================================================
// bad_reconstruction_1_sp
//
// PCA-based stem/shower-direction mismatch check.  Tests whether the shower's
// main PCA axis disagrees with the local stem direction; if so the shower is
// likely a mis-reconstructed track.
//
// Fills TaggerInfo shw_sp_br2_* fields (unconditionally).
//
// Prototype: WCPPID::NeutrinoID::bad_reconstruction_1_sp()
//            NeutrinoID_singlephoton_tagger.h line 4170.
//
// Translation notes:
//   main_cluster->Calc_PCA(pts) + get_PCA_axis(0) → ctx.self.calc_PCA_main_axis(pts).second
//   vertex wcpt-index front/back check            → seg_endpoint_near()
//   map_vertex_segments[other_vertex] iteration   → boost::out_edges loop
// ===========================================================================
static bool bad_reconstruction_1_sp(SpContext& ctx, ShowerPtr shower,
                                     bool flag_single_shower, int num_valid_tracks,
                                     TaggerInfo& ti)
{
    Vector dir_drift(1, 0, 0);
    bool flag_bad_shower = false;

    double Eshower = shower_energy(shower);

    SegmentPtr sg  = shower->start_segment();
    if (!sg) return false;
    VertexPtr vertex = shower->get_start_vertex_and_type().first;

    // Collect shower points (main-cluster only, flag_main=true).
    // Prototype: shower->fill_point_vec(tmp_pts, true)
    std::vector<Point> tmp_pts;
    shower->fill_point_vector(tmp_pts, /*flag_main=*/true);

    // Overall shower direction from the vertex end.
    // Prototype: dir_shower = shower->cal_dir_3vector(vertex_point, 100*units::cm)
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));
    Vector dir_shower  = shower_cal_dir_3vector(*shower, vertex_point, 100 * units::cm);

    // PCA of shower points → main axis.
    // Prototype: main_cluster->Calc_PCA(tmp_pts); dir1 = get_PCA_axis(0)
    Vector dir1;
    if (!tmp_pts.empty())
        dir1 = ctx.self.calc_PCA_main_axis(tmp_pts).second;

    double angle  = 0;  // angle between PCA axis and stem direction
    double angle1 = 0;  // shower drift-angle deviation (drift-perpendicular)
    double angle2 = std::fabs(dir1.angle(dir_drift) / M_PI * 180.0 - 90.0);
    double angle3 = 0;  // stem vs shower-direction angle

    // Determine stem direction from the vertex-side endpoint.
    const auto& sg_fits = sg->fits();
    Point sg_front = sg_fits.front().point;
    Point sg_back  = sg_fits.back().point;
    bool vertex_at_front = (ray_length(Ray{vtx_fit_pt(vertex), sg_front}) <=
                            ray_length(Ray{vtx_fit_pt(vertex), sg_back}));

    if (vertex_at_front) {
        Vector dir2 = segment_cal_dir_3vector(sg, sg_front, 5 * units::cm);
        Vector dir3 = shower_cal_dir_3vector(*shower, sg_front, 30 * units::cm);
        angle  = dir1.angle(dir2) / M_PI * 180.0;
        if (angle > 90) angle = 180.0 - angle;
        angle1 = std::fabs(dir3.angle(dir_drift) / M_PI * 180.0 - 90.0);
        angle3 = dir_shower.angle(dir2) / M_PI * 180.0;
    } else {
        Vector dir2 = segment_cal_dir_3vector(sg, sg_back, 5 * units::cm);
        Vector dir3 = shower_cal_dir_3vector(*shower, sg_back, 30 * units::cm);
        angle  = dir1.angle(dir2) / M_PI * 180.0;
        if (angle > 90) angle = 180.0 - angle;
        angle1 = std::fabs(dir3.angle(dir_drift) / M_PI * 180.0 - 90.0);
        angle3 = dir_shower.angle(dir2) / M_PI * 180.0;
    }

    // Find the maximum opening angle between the stem and any other segment
    // at the far-end vertex of the stem.
    double max_angle   = 0;
    VertexPtr other_vertex = find_other_vertex(ctx.graph, sg, vertex);
    if (other_vertex) {
        Point other_pt  = vtx_fit_pt(other_vertex);
        Vector dir_stem = segment_cal_dir_3vector(sg, other_pt, 10 * units::cm);
        auto [eit, eend] = boost::out_edges(other_vertex->get_descriptor(), ctx.graph);
        for (; eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            Vector dir_other = segment_cal_dir_3vector(sg1, other_pt, 10 * units::cm);
            double ang = dir_stem.angle(dir_other) / M_PI * 180.0;
            if (ang > max_angle) max_angle = ang;
        }
    }

    // Apply cuts: flag as bad if stem direction disagrees with shower/PCA axes.
    if (flag_single_shower || num_valid_tracks == 0) {
        if (Eshower > 1000 * units::MeV) {
            // no cut at very high energy
        } else if (Eshower > 500 * units::MeV) {
            if ((angle1 > 10 || angle2 > 10) && angle > 30) {
                if (angle3 > 3) flag_bad_shower = true;
            }
        } else {
            if (((angle > 25 && shower->get_num_main_segments() > 1) || angle > 30) &&
                (angle1 > 7.5 || angle2 > 7.5)) {
                flag_bad_shower = true;
            }
        }
    }

    // Additional cuts (prototype lines 4248-4252)
    if (angle > 40 && (angle1 > 7.5 || angle2 > 7.5) && max_angle < 100)
        flag_bad_shower = true;
    if (angle > 20 && (angle1 > 7.5 || angle2 > 7.5) &&
        segment_track_length(sg) > 21 * units::cm && Eshower < 600 * units::MeV &&
        sg->flags_any(SegmentFlags::kShowerTrajectory)) {
        flag_bad_shower = true;
    }

    ti.shw_sp_br2_flag               = !flag_bad_shower;
    ti.shw_sp_br2_flag_single_shower  = flag_single_shower;
    ti.shw_sp_br2_num_valid_tracks    = num_valid_tracks;
    ti.shw_sp_br2_energy              = static_cast<float>(Eshower / units::MeV);
    ti.shw_sp_br2_angle1              = static_cast<float>(angle1);
    ti.shw_sp_br2_angle2              = static_cast<float>(angle2);
    ti.shw_sp_br2_angle               = static_cast<float>(angle);
    ti.shw_sp_br2_angle3              = static_cast<float>(angle3);
    ti.shw_sp_br2_n_shower_main_segs  = shower->get_num_main_segments();
    ti.shw_sp_br2_max_angle           = static_cast<float>(max_angle);
    ti.shw_sp_br2_sg_length           = static_cast<float>(segment_track_length(sg) / units::cm);
    ti.shw_sp_br2_flag_sg_trajectory  = sg->flags_any(SegmentFlags::kShowerTrajectory);

    return flag_bad_shower;
}

// ===========================================================================
// bad_reconstruction_2_sp
//
// Eight sub-checks (br3_1 … br3_8) testing whether the shower is really a
// mis-reconstructed track based on segment topology, dQ/dx, and angular
// structure relative to the stem direction.
//
// Fills TaggerInfo shw_sp_br3_* fields (unconditionally).
//
// Prototype: WCPPID::NeutrinoID::bad_reconstruction_2_sp()
//            NeutrinoID_singlephoton_tagger.h line 3461.
//
// Logic is identical to bad_reconstruction_2() in NeutrinoTaggerNuE.cxx;
// only TaggerInfo field names differ (shw_sp_br3_* vs br3_*).
// ===========================================================================
static bool bad_reconstruction_2_sp(SpContext& ctx,
                                     VertexPtr vertex, ShowerPtr shower,
                                     TaggerInfo& ti)
{
    bool flag_bad1 = false, flag_bad2 = false;
    bool flag_bad3_save = false, flag_bad4 = false, flag_bad5 = false;
    bool flag_bad6_save = false, flag_bad7 = false, flag_bad8 = false;

    Vector drift_dir(1, 0, 0);
    double Eshower = shower_energy(shower);

    SegmentPtr sg            = shower->start_segment();
    double total_length      = shower->get_total_length();
    double total_main_length = sg->cluster() ? shower->get_total_length(sg->cluster()) : 0;
    double length            = segment_track_length(sg);
    double direct_length     = segment_track_direct_length(sg);

    // End-to-end direction of start segment.
    const auto& sg_fits  = sg->fits();
    Vector dir_two_end   = sg_fits.empty() ? Vector(0, 0, 0)
                           : (sg_fits.front().point - sg_fits.back().point);

    // -------------------------------------------------------------------
    // br3_1: straight low-energy shower.
    // Prototype lines 3493-3499.
    // -------------------------------------------------------------------
    if (Eshower < 100*units::MeV && shower->get_num_segments() == 1 &&
        !sg->flags_any(SegmentFlags::kShowerTrajectory) &&
        direct_length / length > 0.95) flag_bad1 = true;
    if (Eshower < 100*units::MeV && total_main_length/total_length > 0.95 &&
        length/total_length > 0.85 &&
        (direct_length/length > 0.95 ||
         std::fabs(dir_two_end.angle(drift_dir)/M_PI*180.0 - 90.0) < 5.0) &&
        sg->flags_any(SegmentFlags::kShowerTrajectory)) flag_bad1 = true;
    if (Eshower < 200*units::MeV && total_main_length/total_length > 0.96 &&
        length/total_length > 0.925 &&
        (direct_length/length > 0.95 ||
         (std::fabs(dir_two_end.angle(drift_dir)/M_PI*180.0 - 90.0) < 5.0 &&
          sg->flags_any(SegmentFlags::kShowerTrajectory))) &&
        length > 25*units::cm) flag_bad1 = true;
    if (Eshower < 100*units::MeV && total_main_length/total_length > 0.95 &&
        length/total_length > 0.95 && direct_length/length > 0.95 &&
        sg->flags_any(SegmentFlags::kShowerTopology)) flag_bad1 = true;

    ti.shw_sp_br3_1_energy            = Eshower / units::MeV;
    ti.shw_sp_br3_1_n_shower_segments = shower->get_num_segments();
    ti.shw_sp_br3_1_sg_flag_trajectory= sg->flags_any(SegmentFlags::kShowerTrajectory);
    ti.shw_sp_br3_1_sg_flag_topology  = sg->flags_any(SegmentFlags::kShowerTopology);
    ti.shw_sp_br3_1_sg_direct_length  = direct_length / units::cm;
    ti.shw_sp_br3_1_sg_length         = length / units::cm;
    ti.shw_sp_br3_1_total_main_length = total_main_length / units::cm;
    ti.shw_sp_br3_1_total_length      = total_length / units::cm;
    ti.shw_sp_br3_1_iso_angle         = std::fabs(dir_two_end.angle(drift_dir)/M_PI*180.0 - 90.0);
    ti.shw_sp_br3_1_flag              = !flag_bad1;

    // -------------------------------------------------------------------
    // br3_2: shower segment type composition + fiducial check.
    // Prototype lines 3520-3550.
    // -------------------------------------------------------------------
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    int n_ele = 0, n_other = 0;
    for (SegmentPtr sg1 : shower_segs) {
        if (sg1->cluster() != sg->cluster()) continue;
        double med   = segment_median_dQ_dx(sg1) / (43e3/units::cm);
        double ratio = segment_track_direct_length(sg1) / segment_track_length(sg1);
        if (sg1->flags_any(SegmentFlags::kShowerTopology) ||
            (sg1->flags_any(SegmentFlags::kShowerTrajectory) && med < 1.3) ||
            ratio < 0.92) ++n_ele;
        else if (med > 1.3 || ratio > 0.95) ++n_other;
    }

    VertexPtr other_vertex = find_other_vertex(ctx.graph, sg, vertex);

    FiducialUtilsPtr fiducial_utils;
    if (ctx.main_cluster && ctx.main_cluster->grouping())
        fiducial_utils = ctx.main_cluster->grouping()->get_fiducialutils();
    bool other_fid = fiducial_utils
                     ? fiducial_utils->inside_fiducial_volume(vtx_fit_pt(other_vertex)) : true;

    if (Eshower < 150*units::MeV && total_main_length/total_length > 0.95 &&
        ((n_ele == 0 && n_other > 0) ||
         (n_ele == 1 && n_ele < n_other && n_other <= 2))) flag_bad2 = true;
    if (Eshower < 150*units::MeV && total_main_length/total_length > 0.95 &&
        n_ele == 1 && n_other == 0 && !other_fid) flag_bad2 = true;

    ti.shw_sp_br3_2_n_ele             = n_ele;
    ti.shw_sp_br3_2_n_other           = n_other;
    ti.shw_sp_br3_2_energy            = Eshower / units::MeV;
    ti.shw_sp_br3_2_total_main_length = total_main_length / units::cm;
    ti.shw_sp_br3_2_total_length      = total_length / units::cm;
    ti.shw_sp_br3_2_other_fid         = other_fid;
    ti.shw_sp_br3_2_flag              = !flag_bad2;

    // -------------------------------------------------------------------
    // br3_3 / br3_4: backward segments in main cluster.
    // Prototype lines 3556-3617.
    // -------------------------------------------------------------------
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));
    Point other_point  = (ray_length(Ray{vtx_fit_pt(vertex), sg_fits.front().point}) <=
                          ray_length(Ray{vtx_fit_pt(vertex), sg_fits.back().point}))
                         ? sg_fits.back().point : sg_fits.front().point;

    Vector dir_stem = segment_cal_dir_3vector(sg, vertex_point, 15 * units::cm);

    double acc_length = 0, total_main_len2 = 0;
    for (SegmentPtr sg1 : shower_segs) {
        if (sg1->cluster() != sg->cluster()) continue;
        const auto& fits1 = sg1->fits();
        Point front1 = fits1.front().point, back1 = fits1.back().point;
        double d1 = ray_length(Ray{vertex_point, front1});
        double d2 = ray_length(Ray{vertex_point, back1});
        Vector dir1 = (d1 < d2) ? (back1 - front1) : (front1 - back1);
        double len1 = segment_track_length(sg1);
        bool flag_bad3 = false;
        if (dir1.magnitude() > 10*units::cm) {
            double angle = dir1.angle(dir_stem) / M_PI * 180.0;
            if (angle > 90)  acc_length += len1;
            if (angle > 150 && Eshower < 600*units::MeV) flag_bad3 = true;
            if (angle > 105 && len1 > 15*units::cm && Eshower < 600*units::MeV) flag_bad3 = true;

            ti.shw_sp_br3_3_v_energy.push_back(Eshower / units::MeV);
            ti.shw_sp_br3_3_v_angle.push_back(angle);
            ti.shw_sp_br3_3_v_dir_length.push_back(dir1.magnitude() / units::cm);
            ti.shw_sp_br3_3_v_length.push_back(len1 / units::cm);
            ti.shw_sp_br3_3_v_flag.push_back(!flag_bad3);
        }
        total_main_len2 += len1;
        if (flag_bad3) flag_bad3_save = true;
    }
    if (acc_length > 0.33 * total_main_len2 && Eshower < 600*units::MeV) flag_bad4 = true;

    ti.shw_sp_br3_4_acc_length   = acc_length / units::cm;
    ti.shw_sp_br3_4_total_length = total_main_len2 / units::cm;
    ti.shw_sp_br3_4_energy       = Eshower / units::MeV;
    ti.shw_sp_br3_4_flag         = !flag_bad4;

    // -------------------------------------------------------------------
    // br3_5: average position of non-stem main-cluster segments.
    // Prototype lines 3621-3674.
    // -------------------------------------------------------------------
    {
        Point  ave_p(0, 0, 0);
        int    num_p = 0, n_seg = 0;
        double side_total_length = 0;
        for (SegmentPtr sg1 : shower_segs) {
            if (sg1->cluster() != sg->cluster() || sg1 == sg) continue;
            for (const auto& fit : sg1->fits()) {
                ave_p = Point(ave_p.x() + fit.point.x(),
                              ave_p.y() + fit.point.y(),
                              ave_p.z() + fit.point.z());
                ++num_p;
            }
            ++n_seg;
            side_total_length += segment_track_length(sg1);
        }

        if (num_p > 0) {
            ave_p = Point(ave_p.x()/num_p, ave_p.y()/num_p, ave_p.z()/num_p);
            Vector dir1      = ave_p - other_point;
            bool avoid_check = sg->flags_any(SegmentFlags::kAvoidMuonCheck);

            if ((dir1.magnitude() > 3*units::cm || side_total_length > 6*units::cm) &&
                (!avoid_check || n_seg > 1) &&
                dir_stem.angle(dir1)/M_PI*180.0 > 60 &&
                length > 10*units::cm && Eshower < 250*units::MeV)
                flag_bad5 = true;
            // 7018_888_44410
            if (shower->get_num_main_segments() + 6 < shower->get_num_segments() &&
                sg->cluster() &&
                shower->get_total_length(sg->cluster()) < 0.7 * shower->get_total_length() &&
                Eshower < 250*units::MeV)
                flag_bad5 = false;

            ti.shw_sp_br3_5_v_dir_length.push_back(dir1.magnitude() / units::cm);
            ti.shw_sp_br3_5_v_total_length.push_back(side_total_length / units::cm);
            ti.shw_sp_br3_5_v_flag_avoid_muon_check.push_back(avoid_check);
            ti.shw_sp_br3_5_v_n_seg.push_back(n_seg);
            ti.shw_sp_br3_5_v_angle.push_back(dir_stem.angle(dir1) / M_PI * 180.0);
            ti.shw_sp_br3_5_v_sg_length.push_back(length / units::cm);
            ti.shw_sp_br3_5_v_energy.push_back(Eshower / units::MeV);
            ti.shw_sp_br3_5_v_n_main_segs.push_back(shower->get_num_main_segments());
            ti.shw_sp_br3_5_v_n_segs.push_back(shower->get_num_segments());
            ti.shw_sp_br3_5_v_shower_main_length.push_back(
                sg->cluster() ? shower->get_total_length(sg->cluster()) / units::cm : 0);
            ti.shw_sp_br3_5_v_shower_total_length.push_back(shower->get_total_length() / units::cm);
            ti.shw_sp_br3_5_v_flag.push_back(!flag_bad5);
        }
    }

    // -------------------------------------------------------------------
    // br3_6 / br3_7: segments at the far end of the stem.
    // Prototype lines 3680-3725.
    // -------------------------------------------------------------------
    other_vertex = find_other_vertex(ctx.graph, sg, vertex);
    double min_angle = 180;
    if (other_vertex && other_vertex->descriptor_valid()) {
        Point  ovp = vtx_fit_pt(other_vertex);
        size_t n_other_vtx_segs = boost::out_degree(other_vertex->get_descriptor(), ctx.graph);
        for (auto [eit, eend] = boost::out_edges(other_vertex->get_descriptor(), ctx.graph);
             eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            VertexPtr vtx_1 = find_other_vertex(ctx.graph, sg1, other_vertex);
            if (!vtx_1) continue;
            Vector dir1   = vtx_fit_pt(vtx_1) - ovp;
            double angle  = dir1.angle(dir_stem) / M_PI * 180.0;
            double angle1 = std::max(
                std::fabs(90.0 - dir_stem.angle(drift_dir)/M_PI*180.0),
                std::fabs(90.0 - dir1.angle(drift_dir)/M_PI*180.0));
            double sg1_len = segment_track_length(sg1);
            double sg1_dir = segment_track_direct_length(sg1);
            bool flag_bad6 = false;
            if (angle > 150 && angle1 > 10 &&
                !sg1->flags_any(SegmentFlags::kShowerTrajectory) &&
                sg1_dir/sg1_len > 0.9 &&
                sg1_len > 7.5*units::cm && n_other_vtx_segs <= 4 &&
                Eshower < 600*units::MeV) flag_bad6 = true;
            if (angle < min_angle && sg1_len > 6*units::cm) min_angle = angle;

            ti.shw_sp_br3_6_v_angle.push_back(angle);
            ti.shw_sp_br3_6_v_angle1.push_back(angle1);
            ti.shw_sp_br3_6_v_flag_shower_trajectory.push_back(
                sg1->flags_any(SegmentFlags::kShowerTrajectory));
            ti.shw_sp_br3_6_v_direct_length.push_back(sg1_dir / units::cm);
            ti.shw_sp_br3_6_v_length.push_back(sg1_len / units::cm);
            ti.shw_sp_br3_6_v_n_other_vtx_segs.push_back(n_other_vtx_segs);
            ti.shw_sp_br3_6_v_energy.push_back(Eshower / units::MeV);
            ti.shw_sp_br3_6_v_flag.push_back(!flag_bad6);
            if (flag_bad6) flag_bad6_save = true;
        }
    }

    double shower_main_len = vertex->cluster() ? shower->get_total_length(vertex->cluster()) : 0;
    if (Eshower < 200*units::MeV && min_angle > 60 &&
        length < 0.2 * shower_main_len) flag_bad7 = true;

    ti.shw_sp_br3_7_energy             = Eshower / units::MeV;
    ti.shw_sp_br3_7_min_angle          = min_angle;
    ti.shw_sp_br3_7_sg_length          = length / units::cm;
    ti.shw_sp_br3_7_shower_main_length = shower_main_len / units::cm;
    ti.shw_sp_br3_7_flag               = !flag_bad7;

    // -------------------------------------------------------------------
    // br3_8: sliding-window dQ/dx peak across main-cluster shower segments.
    // Prototype lines 3731-3754.
    // -------------------------------------------------------------------
    double max_dQ_dx = 0;
    for (SegmentPtr sg1 : shower_segs) {
        if (sg1->cluster() != vertex->cluster()) continue;
        int n = (int)sg1->fits().size();
        for (int i = 0; i < n - 5; ++i) {
            double med = segment_median_dQ_dx(sg1, i, i+5) / (43e3/units::cm);
            if (med > max_dQ_dx) max_dQ_dx = med;
        }
    }
    if (max_dQ_dx > 1.85 && Eshower < 150*units::MeV &&
        shower->get_num_main_segments() <= 2 &&
        vertex->cluster() &&
        shower->get_total_length(vertex->cluster()) > shower->get_total_length() * 0.8)
        flag_bad8 = true;

    ti.shw_sp_br3_8_max_dQ_dx          = max_dQ_dx;
    ti.shw_sp_br3_8_energy             = Eshower / units::MeV;
    ti.shw_sp_br3_8_n_main_segs        = shower->get_num_main_segments();
    ti.shw_sp_br3_8_shower_main_length = vertex->cluster()
                                         ? shower->get_total_length(vertex->cluster()) / units::cm : 0;
    ti.shw_sp_br3_8_shower_length      = shower->get_total_length() / units::cm;
    ti.shw_sp_br3_8_flag               = !flag_bad8;

    bool flag_bad = flag_bad1 || flag_bad2 || flag_bad3_save || flag_bad4 ||
                    flag_bad5 || flag_bad6_save || flag_bad7 || flag_bad8;
    ti.shw_sp_br3_flag = !flag_bad;
    return flag_bad;
}

// ===========================================================================
// bad_reconstruction_3_sp
//
// Two sub-checks (br4_1, br4_2):
//   br4_1: main-cluster fraction vs distance to closest off-cluster segment
//   br4_2: angular distribution of all shower hit points relative to shower dir
//
// Fills TaggerInfo shw_sp_br4_* fields (unconditionally).
//
// Prototype: WCPPID::NeutrinoID::bad_reconstruction_3_sp()
//            NeutrinoID_singlephoton_tagger.h line 3223.
//
// Logic is identical to bad_reconstruction_3() in NeutrinoTaggerNuE.cxx;
// only TaggerInfo field names differ (shw_sp_br4_* vs br4_*).
// ===========================================================================
static bool bad_reconstruction_3_sp(SpContext& ctx,
                                     VertexPtr vertex, ShowerPtr shower,
                                     TaggerInfo& ti)
{
    bool flag_bad1 = false, flag_bad2 = false;

    Vector drift_dir(1, 0, 0);
    double Eshower      = shower_energy(shower);
    double main_length  = vertex->cluster() ? shower->get_total_length(vertex->cluster()) : 0;
    double total_length = shower->get_total_length();

    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    // -------------------------------------------------------------------
    // br4_1: find the farthest shower vertex in the main cluster, then
    // measure distance to closest off-cluster shower segment.
    // Prototype lines 3245-3305.
    // -------------------------------------------------------------------
    double max_dis_v = 0;
    Point  max_p     = vtx_fit_pt(vertex);
    for (VertexPtr vtx1 : shower_vtxs) {
        if (vtx1->cluster() != vertex->cluster()) continue;
        double d = ray_length(Ray{vtx_fit_pt(vertex), vtx_fit_pt(vtx1)});
        if (d > max_dis_v) { max_dis_v = d; max_p = vtx_fit_pt(vtx1); }
    }

    double min_dis = 1e9;
    for (SegmentPtr sg1 : shower_segs) {
        if (sg1->cluster() == vertex->cluster()) continue;
        if (segment_track_length(sg1) < 6*units::cm) continue;
        double d = segment_get_closest_point(sg1, max_p).first;
        if (d < min_dis) min_dis = d;
    }

    double acc_close_length = 0;
    int    num_close        = 0;
    double min_dis1         = 1e9;
    for (SegmentPtr sg1 : shower_segs) {
        if (sg1->cluster() == vertex->cluster()) continue;
        double d = segment_get_closest_point(sg1, max_p).first;
        if (d < min_dis) {
            double len1 = segment_track_length(sg1);
            acc_close_length += len1;
            ++num_close;
            if (len1 > 3*units::cm && d < min_dis1) min_dis1 = d;
        }
    }
    if (min_dis1 > 1e8) min_dis1 = 0;

    SegmentPtr start_sg = shower->start_segment();
    if (acc_close_length > 10*units::cm ||
        (num_close >= 3 && acc_close_length > 4.5*units::cm) ||
        (start_sg && start_sg->flags_any(SegmentFlags::kAvoidMuonCheck)))
        min_dis = min_dis1;

    VertexPtr start_vtx = shower->get_start_vertex_and_type().first;
    size_t n_vtx_segs = (start_vtx && start_vtx->descriptor_valid())
                        ? boost::out_degree(start_vtx->get_descriptor(), ctx.graph) : 0;

    if (min_dis < 1e7) {
        if (main_length < 0.40 * total_length && min_dis > 40*units::cm) flag_bad1 = true;
        if (main_length < 0.25 * total_length && min_dis > 33*units::cm) flag_bad1 = true;
        if (main_length < 0.16 * total_length && min_dis > 23*units::cm) flag_bad1 = true;
        if (main_length < 0.10 * total_length && min_dis > 18*units::cm) flag_bad1 = true;
        if (main_length < 0.05 * total_length && min_dis >  8*units::cm) flag_bad1 = true;
        if (main_length < 8*units::cm && main_length < 0.1*total_length &&
            ((min_dis > 8*units::cm && Eshower < 300*units::MeV) || min_dis > 14*units::cm))
            flag_bad1 = true;

        if (flag_bad1 && start_sg && start_sg->flags_any(SegmentFlags::kAvoidMuonCheck) &&
            main_length > 12*units::cm && main_length > 0.1*total_length &&
            min_dis < 40*units::cm) flag_bad1 = false;
        if (n_vtx_segs == 1 &&
            ((main_length > 20*units::cm && min_dis < 40*units::cm && main_length > 0.1*total_length) ||
             (main_length > 15*units::cm && min_dis < 32*units::cm && main_length > 0.15*total_length)))
            flag_bad1 = false;
        if (flag_bad1 && main_length > 30*units::cm && shower->get_num_main_segments() >= 4)
            flag_bad1 = false;
    }

    ti.shw_sp_br4_1_shower_main_length    = main_length / units::cm;
    ti.shw_sp_br4_1_shower_total_length   = total_length / units::cm;
    ti.shw_sp_br4_1_min_dis               = min_dis / units::cm;
    ti.shw_sp_br4_1_energy                = Eshower / units::MeV;
    ti.shw_sp_br4_1_flag_avoid_muon_check = (start_sg && start_sg->flags_any(SegmentFlags::kAvoidMuonCheck));
    ti.shw_sp_br4_1_n_vtx_segs           = (int)n_vtx_segs;
    ti.shw_sp_br4_1_n_main_segs          = shower->get_num_main_segments();
    ti.shw_sp_br4_1_flag                 = !flag_bad1;

    // -------------------------------------------------------------------
    // br4_2: angular distribution of shower fit points relative to shower dir.
    // Prototype lines 3322-3448.
    // -------------------------------------------------------------------
    {
        SegmentPtr sg = shower->start_segment();
        Point vp = seg_endpoint_near(sg, vtx_fit_pt(vertex));

        Vector dir_sg = segment_cal_dir_3vector(sg, vp, 15*units::cm);
        Vector dir;
        if (segment_track_length(sg) > 12*units::cm)
            dir = segment_cal_dir_3vector(sg, vp, 15*units::cm);
        else
            dir = shower_cal_dir_3vector(*shower, vp, 15*units::cm);
        if (std::fabs(dir.angle(drift_dir) / M_PI * 180.0 - 90.0) < 10.0)
            dir = shower_cal_dir_3vector(*shower, vp, 25*units::cm);

        int ncount = 0, ncount1 = 0;
        int ncount_15 = 0, ncount1_15 = 0;
        int ncount_25 = 0, ncount1_25 = 0;
        int ncount_35 = 0, ncount1_35 = 0;
        int ncount_45 = 0, ncount1_45 = 0;

        auto count_dir = [&](const Vector& dir1, bool is_main) {
            double ang = dir1.angle(dir) / M_PI * 180.0;
            if (ang < 15) ++ncount1_15;
            if (ang < 25) ++ncount1_25;
            if (ang < 35) ++ncount1_35;
            if (ang < 45) ++ncount1_45;
            ++ncount1;
            if (is_main) {
                if (ang < 15) ++ncount_15;
                if (ang < 25) ++ncount_25;
                if (ang < 35) ++ncount_35;
                if (ang < 45) ++ncount_45;
                ++ncount;
            }
        };

        for (SegmentPtr sg1 : shower_segs) {
            bool is_main = (sg1->cluster() == vertex->cluster());
            const auto& fits1 = sg1->fits();
            // Interior fit points
            for (size_t i = 1; i + 1 < fits1.size(); ++i)
                count_dir(fits1[i].point - vp, is_main);
            // Endpoint vertices
            auto [v1, v2] = find_vertices(ctx.graph, sg1);
            for (VertexPtr ep : {v1, v2})
                if (ep) count_dir(vtx_fit_pt(ep) - vp, is_main);
        }

        if (ncount_45 < 0.7*ncount || ncount_25 < 0.6*ncount ||
            (ncount_25 < 0.8*ncount && ncount_15 < 0.3*ncount) ||
            (ncount_15 < 0.35*ncount && ncount_25 > 0.9*ncount && Eshower < 1000*units::MeV))
            flag_bad2 = true;

        double iso_angle   = std::fabs(dir.angle(drift_dir)    / M_PI * 180.0 - 90.0);
        double iso_angle1  = std::fabs(dir_sg.angle(drift_dir) / M_PI * 180.0 - 90.0);
        double sg_dir_angle= dir_sg.angle(dir) / M_PI * 180.0;

        bool cut_b =
            (ncount1_15 < 0.35*ncount1 && iso_angle > 15 &&
             ((ncount1_25 < 0.95*ncount1 && Eshower < 1000*units::MeV) || Eshower >= 1000*units::MeV)) ||
            (ncount1_15 < 0.2*ncount1 && ncount1_25 < 0.45*ncount1 && Eshower < 600*units::MeV) ||
            (sg_dir_angle > 25 && std::max(iso_angle, iso_angle1) > 8 &&
             ((ncount1_15 < 0.8*ncount1 && Eshower < 1000*units::MeV) || Eshower >= 1000*units::MeV)) ||
            (sg_dir_angle > 20 && std::max(iso_angle, iso_angle1) > 5 &&
             ncount1_15 < 0.5*ncount1);
        if (cut_b) flag_bad2 = true;

        if (ncount > 0) {
            ti.shw_sp_br4_2_ratio_45  = ncount_45  / (ncount + 1e-9);
            ti.shw_sp_br4_2_ratio_35  = ncount_35  / (ncount + 1e-9);
            ti.shw_sp_br4_2_ratio_25  = ncount_25  / (ncount + 1e-9);
            ti.shw_sp_br4_2_ratio_15  = ncount_15  / (ncount + 1e-9);
        } else {
            ti.shw_sp_br4_2_ratio_45  = ti.shw_sp_br4_2_ratio_35  = 1;
            ti.shw_sp_br4_2_ratio_25  = ti.shw_sp_br4_2_ratio_15  = 1;
        }
        ti.shw_sp_br4_2_energy = Eshower / units::MeV;
        if (ncount1 > 0) {
            ti.shw_sp_br4_2_ratio1_45 = ncount1_45 / (ncount1 + 1e-9);
            ti.shw_sp_br4_2_ratio1_35 = ncount1_35 / (ncount1 + 1e-9);
            ti.shw_sp_br4_2_ratio1_25 = ncount1_25 / (ncount1 + 1e-9);
            ti.shw_sp_br4_2_ratio1_15 = ncount1_15 / (ncount1 + 1e-9);
        } else {
            ti.shw_sp_br4_2_ratio1_45 = ti.shw_sp_br4_2_ratio1_35 = 1;
            ti.shw_sp_br4_2_ratio1_25 = ti.shw_sp_br4_2_ratio1_15 = 1;
        }
        ti.shw_sp_br4_2_iso_angle  = iso_angle;
        ti.shw_sp_br4_2_iso_angle1 = iso_angle1;
        ti.shw_sp_br4_2_angle      = sg_dir_angle;
        ti.shw_sp_br4_2_flag       = !flag_bad2;
    }

    bool flag_bad = flag_bad1 || flag_bad2;
    ti.shw_sp_br4_flag = !flag_bad;
    return flag_bad;
}

// ===========================================================================
// mip_identification_sp
//
// Classifies the shower stem as MIP-like (mip_id=1), photon-like (-1), or
// ambiguous (0) based on the near-vertex dQ/dx pattern.
// Returns mip_id and fills many TaggerInfo shw_sp_* fields.
//
// Prototype: WCPPID::NeutrinoID::mip_identification_sp()
//            NeutrinoID_singlephoton_tagger.h line 2102.
//
// Translation notes:
//   shower->get_stem_dQ_dx(vertex, sg, 20) → identical (toolkit takes VertexPtr,SegmentPtr)
//   map_vertex_segments[v].size()           → vtx_degree(v, ctx.graph)
//   map_vertex_segments[v] iteration        → boost::out_edges loop
//   shower->get_total_length(id)            → shower->get_total_length(cluster*)
//   sg->get_point_vec().size()              → sg->fits().size()
//   shower->fill_sets for min_dis           → iterate shower_segs / shower_vtxs
// ===========================================================================
static int mip_identification_sp(SpContext& ctx,
                                  VertexPtr vertex, SegmentPtr sg, ShowerPtr shower,
                                  bool flag_single_shower, bool flag_strong_check,
                                  TaggerInfo& ti)
{
    int mip_id = 1; // 1=good (MIP), -1=bad (photon), 0=not sure

    Vector dir_beam(0, 0, 1);
    Vector dir_drift(1, 0, 0);
    double Eshower = shower_energy(shower);

    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));

    // Stem direction: use start segment if long enough, else use shower direction
    Vector dir_shower;
    SegmentPtr start_sg = shower->start_segment();
    if (start_sg && segment_track_length(start_sg) > 12*units::cm) {
        Point vp = vertex_point;
        dir_shower = segment_cal_dir_3vector(start_sg, vp, 15*units::cm);
    } else {
        dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
    }
    if (std::fabs(dir_shower.angle(dir_drift) / M_PI * 180.0 - 90.0) < 10 ||
        Eshower > 800*units::MeV)
        dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 25*units::cm);
    if (dir_shower.magnitude() > 0) dir_shower = dir_shower.norm();

    double dQ_dx_cut = 1.45;
    if      (Eshower > 1200*units::MeV) dQ_dx_cut = 1.85;
    else if (Eshower > 1000*units::MeV) dQ_dx_cut = 1.6;
    else if (Eshower <  550*units::MeV) dQ_dx_cut = 1.3;
    if      (Eshower <  300*units::MeV) dQ_dx_cut = 1.3;

    std::vector<double> vec_dQ_dx = shower->get_stem_dQ_dx(vertex, sg, 20);

    std::vector<int> vec_threshold(vec_dQ_dx.size(), 0);
    for (size_t i = 0; i < vec_dQ_dx.size(); ++i)
        if (vec_dQ_dx[i] > dQ_dx_cut) vec_threshold[i] = 1;

    int n_end_reduction = 0;
    double prev_dQ = vec_dQ_dx.front();
    for (size_t i = 1; i < vec_dQ_dx.size(); ++i) {
        if (vec_dQ_dx[i] < prev_dQ) {
            n_end_reduction = i;
            prev_dQ = vec_dQ_dx[i];
            if (vec_dQ_dx[i] < dQ_dx_cut) break;
        }
    }

    int n_first_mip = 0;
    for (size_t i = 0; i < vec_dQ_dx.size(); ++i) {
        n_first_mip = i;
        if (vec_threshold[i] == 0) break;
    }

    int n_first_non_mip = n_first_mip;
    for (size_t i = n_first_non_mip; i < vec_dQ_dx.size(); ++i) {
        n_first_non_mip = i;
        if (vec_threshold[i] == 1) break;
    }

    int n_first_non_mip_1 = n_first_mip;
    for (size_t i = n_first_non_mip; i < vec_dQ_dx.size(); ++i) {
        n_first_non_mip_1 = i;
        if (vec_threshold[i] == 1 && i+1 < vec_dQ_dx.size())
            if (vec_threshold[i+1] == 1) break;
    }

    int n_first_non_mip_2 = n_first_mip;
    for (size_t i = n_first_non_mip; i < vec_dQ_dx.size(); ++i) {
        n_first_non_mip_2 = i;
        if (vec_threshold[i] == 1 && i+1 < vec_dQ_dx.size())
            if (vec_threshold[i+1] == 1 && i+2 < vec_dQ_dx.size())
                if (vec_threshold[i+2] == 1) break;
    }

    double lowest_dQ_dx  = 100; int n_lowest  = 0;
    double highest_dQ_dx = 0;   int n_highest = 0;
    int n_below_threshold = 0;
    int n_below_zero = 0;
    for (size_t i = n_first_mip; i < (size_t)n_first_non_mip_2; ++i) {
        if (vec_dQ_dx[i] < lowest_dQ_dx && (int)i <= 12) {
            lowest_dQ_dx = vec_dQ_dx[i];
            n_lowest = i;
        }
        if (vec_dQ_dx[i] > highest_dQ_dx) {
            highest_dQ_dx = vec_dQ_dx[i];
            n_highest = i;
        }
        if (vec_dQ_dx[i] < dQ_dx_cut) ++n_below_threshold;
        if (vec_dQ_dx[i] < 0) ++n_below_zero;
    }

    // Primary MIP classification
    if (n_first_non_mip_2 - n_first_mip >= 2 &&
        (n_first_mip <= 2 ||
         (n_first_mip <= n_end_reduction &&
          (n_first_mip <= 3 ||
           (n_first_mip <= 4 && n_first_non_mip_1 - n_first_mip > 5 && Eshower > 150*units::MeV) ||
           (n_first_mip <= 4 && Eshower > 600*units::MeV) ||
           (n_first_mip <= 5 && Eshower > 800*units::MeV) ||
           (n_first_mip <= 6 && Eshower > 1000*units::MeV) ||
           (n_first_mip <= 10 && Eshower > 1000*units::MeV && n_first_non_mip_1 - n_first_mip > 5) ||
           (n_first_mip <= 10 && Eshower > 1250*units::MeV)))))
        mip_id = 1;
    else
        mip_id = -1;

    double max_dQ_dx_sample = 0;
    for (size_t i = n_first_non_mip_2; i < n_first_non_mip_2 + 3; ++i) {
        if (i >= vec_dQ_dx.size()) break;
        if (vec_dQ_dx[i] > max_dQ_dx_sample) max_dQ_dx_sample = vec_dQ_dx[i];
    }

    // Refinement: 7013_63_3191 etc.
    if (mip_id == -1 && n_first_mip <= n_end_reduction && n_first_mip <= 5 &&
        ((((n_first_non_mip_2 - n_first_mip >= 8 && n_first_non_mip - n_first_mip >= 7) ||
           (n_first_non_mip_2 - n_first_mip >= 5 && max_dQ_dx_sample < 1.6)) &&
          std::max(vec_dQ_dx[0], vec_dQ_dx[1]) > 1.75) ||
         (n_first_non_mip_2 - n_first_mip >= 5 && std::max(vec_dQ_dx[0], vec_dQ_dx[1]) > 3.0)))
        mip_id = 0;

    // 6640_171_8560 + 7014_954_47722
    if (mip_id == -1 && n_first_mip <= n_end_reduction &&
        std::max(vec_dQ_dx[0], vec_dQ_dx[1]) < 1.45 &&
        (n_first_non_mip_2 - n_first_mip + n_end_reduction >= 12 &&
         n_below_threshold + n_end_reduction >= 10) &&
        n_first_non_mip_2 - n_first_mip >= 4 &&
        (n_end_reduction < 4 && Eshower < 100*units::MeV ||
         n_end_reduction < 7 && Eshower < 200*units::MeV && Eshower >= 100*units::MeV ||
         Eshower >= 200*units::MeV)) {
        if (flag_single_shower) mip_id = 0;
        else                    mip_id = 1;
    }

    // Strong check mode
    if (flag_strong_check) {
        if (!((n_first_mip <= 2 ||
               (n_first_mip <= n_end_reduction &&
                (n_first_mip <= 3 ||
                 (n_first_mip <= 4 && Eshower > 600*units::MeV) ||
                 (n_first_mip <= 5 && Eshower > 800*units::MeV) ||
                 (n_first_mip <= 6 && Eshower > 1000*units::MeV) ||
                 (n_first_mip <= 10 && Eshower > 1250*units::MeV)))) &&
              (n_first_non_mip_2 - n_first_mip > 3 ||
               (n_first_non_mip_2 - n_first_mip == 3 &&
                std::max(vec_dQ_dx[0], vec_dQ_dx[1]) > 3.3))))
            mip_id = -1;
        if (mip_id == -1 && n_first_mip <= n_end_reduction &&
            n_first_mip <= 5 && n_first_non_mip_2 - n_first_mip >= 7 &&
            std::max(vec_dQ_dx[0], vec_dQ_dx[1]) > 3.3)
            mip_id = 0;
    }

    // n_good_tracks at this vertex
    int n_good_tracks = 0;
    if (vertex && vertex->descriptor_valid()) {
        auto vd = vertex->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            if (!sg1->dir_weak() || segment_track_length(sg1) > 10*units::cm) ++n_good_tracks;
        }
    }
    if (Eshower < 600*units::MeV) {
        // 6043_4_243
        if (n_good_tracks > 1 && n_first_non_mip_2 <= 2) mip_id = -1;
    }

    // flag_all_above: all first 6 dQ/dx values > 1.2
    bool flag_all_above = true;
    for (size_t i = 0; i < vec_dQ_dx.size(); ++i) {
        if (vec_dQ_dx[i] < 1.2) { flag_all_above = false; break; }
        if ((int)i > 5) break;
    }

    int vtx_deg = vtx_degree(vertex, ctx.graph);

    // Energy-dependent single-shower / multi-track corrections
    if (mip_id == 1 && vtx_deg == 1 && Eshower < 500*units::MeV) {
        if (Eshower < 180*units::MeV || n_first_mip > 0 ||
            (vec_dQ_dx[0] > 1.15 && n_end_reduction >= n_first_mip && Eshower < 360*units::MeV))
            mip_id = 0;
        if (flag_single_shower && Eshower < 400*units::MeV && n_end_reduction > 0) mip_id = 0;
    } else if (mip_id == 1 && vtx_deg > 1 && Eshower < 300*units::MeV) {
        if (vec_dQ_dx.size() >= 3)
            if (vec_dQ_dx[1] < 0.6 || vec_dQ_dx[2] < 0.6) mip_id = 0;
    } else if (mip_id == 1 && vtx_deg > 1 && Eshower < 600*units::MeV) {
        Vector d15 = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
        if (d15.angle(dir_beam) / M_PI * 180.0 > 60 || n_first_non_mip_1 == 1) mip_id = 0;
        if (flag_all_above) mip_id = 0;
    } else if (mip_id == 1 && flag_single_shower && Eshower < 900*units::MeV) {
        if (n_first_mip != 0) mip_id = 0;
    }

    // Angular cuts on low-energy showers
    {
        Vector d15 = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
        double ang_beam = d15.angle(dir_beam) / M_PI * 180.0;
        double iso = std::fabs(dir_shower.angle(dir_drift) / M_PI * 180.0 - 90.0);
        if (Eshower < 300*units::MeV) {
            if (ang_beam > 40) {
                // 7018_926_46331
                if (((n_first_non_mip_2 - n_first_mip <= 3 && n_first_non_mip_2 <= 3) ||
                     n_first_non_mip_2 - n_first_mip <= 2) &&
                    n_first_mip <= 1 && max_dQ_dx_sample > 1.9) mip_id = -1;
                // 5337_192_9614
                if (flag_single_shower && n_first_mip >= 3 &&
                    n_first_non_mip - n_first_mip <= 1 &&
                    std::max(vec_dQ_dx[0], vec_dQ_dx[1]) < 2.7) mip_id = -1;
                // 7048_108_5419
                if (flag_single_shower && n_first_mip >= 2 &&
                    std::max(vec_dQ_dx[0], vec_dQ_dx[1]) < 2.7) mip_id = -1;
            }
            // 7021_586_29303
            if (ang_beam > 30 && Eshower < 200*units::MeV && flag_single_shower) {
                if (vec_dQ_dx[0] > 1.5 && n_first_mip > 0 &&
                    std::max(vec_dQ_dx[0], vec_dQ_dx[1]) < 2.7) mip_id = -1;
            }
        }

        // min_dQ_dx_5: minimum dQ/dx in first 6 bins
        double min_dQ_dx_5 = 1;
        // 7017_1631_81564
        if (flag_single_shower && Eshower < 500*units::MeV &&
            shower->get_total_length(sg->cluster()) > shower->get_total_length() * 0.95) {
            min_dQ_dx_5 = 1e9;
            for (size_t i = 0; i < vec_dQ_dx.size(); ++i) {
                if (vec_dQ_dx[i] < min_dQ_dx_5) min_dQ_dx_5 = vec_dQ_dx[i];
                if ((int)i > 5) break;
            }
            if (n_first_non_mip_2 - n_first_mip <= 2 && min_dQ_dx_5 > 1.3) mip_id = -1;
        }

        // dQ/dx shape cuts
        if (mip_id == 1) {
            if (n_below_threshold <= 5 &&
                (lowest_dQ_dx < 0.7 ||
                 (lowest_dQ_dx > 1.1 && iso < 15))) mip_id = 0;
        }
        // 7018_235_11772
        if (lowest_dQ_dx > 1.3 && iso < 15 && Eshower < 1000*units::MeV) mip_id = -1;
        // 7049_1241_62062
        if (lowest_dQ_dx < 0 && Eshower < 800*units::MeV && n_below_zero > 2) mip_id = -1;
        // 7025_380_19030
        if (lowest_dQ_dx < 0 && Eshower < 800*units::MeV && n_below_zero <= 2 &&
            highest_dQ_dx > 1.3) mip_id = -1;
        // 7054_1985_99267
        if (lowest_dQ_dx < 0 && n_lowest <= 1 && n_highest < n_lowest &&
            highest_dQ_dx < 0.9) mip_id = -1;
        // 7049_1033_51667
        if (lowest_dQ_dx < 0.6 && highest_dQ_dx < 0.8 && n_lowest <= 1 && n_highest <= 1 &&
            n_first_non_mip_2 - n_first_mip <= 2 && max_dQ_dx_sample > 1.8) mip_id = -1;
        // 7017_1508_75440
        if (lowest_dQ_dx < 0.6 && highest_dQ_dx > 1.3 && n_highest > 1 && n_highest < 4 &&
            Eshower < 1000*units::MeV && std::fabs(n_lowest - n_highest) > 1) mip_id = -1;
        // 7003_1754_87734 + anti 7026_54_2747
        if (lowest_dQ_dx < 0.9 && n_lowest <= 1 && highest_dQ_dx > 1.2 &&
            n_below_threshold <= 4 && Eshower < 1000*units::MeV && iso > 10 &&
            n_first_non_mip_2 < 5 && max_dQ_dx_sample > 1.9) mip_id = -1;
        // 7012_1370_68520 + anti 7010_451_22560
        if (n_lowest <= 2 && n_highest > n_lowest && lowest_dQ_dx > 1.1 && iso < 5 &&
            vtx_deg > 1) mip_id = -1;
        // 7055_147_7354
        if (n_lowest <= 3 && lowest_dQ_dx < 0.7 && highest_dQ_dx > 1.3 &&
            n_highest < n_lowest && iso < 5) mip_id = -1;
        // 7012_297_14884
        if (flag_single_shower && n_below_threshold <= 3 && highest_dQ_dx > 1.2 &&
            Eshower < 800*units::MeV && iso > 7.5) {
            mip_id = -1;
            // 7012_366_18344
            if (n_below_threshold == 3 && std::max(vec_dQ_dx[0], vec_dQ_dx[1]) > 3.5)
                mip_id = 0;
        }
        // 6936_165_8288
        if (Eshower < 800*units::MeV && lowest_dQ_dx < 0.2 && n_lowest <= 3 && iso > 15 &&
            segment_track_length(sg) < 5*units::cm) mip_id = -1;

        // Energy/shower-topology corrections from other showers
        double E_direct_total_energy = 0, E_direct_max_energy = 0;
        double E_indirect_total_energy = 0, E_indirect_max_energy = 0;
        int n_direct_showers = 0, n_indirect_showers = 0;

        for (ShowerPtr shower1 : ctx.showers) {
            SegmentPtr sg1 = shower1->start_segment();
            if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11) continue;
            if (shower1 == shower) continue;
            auto [vtx1, conn1] = shower1->get_start_vertex_and_type();
            double E1 = shower_energy(shower1);
            if (conn1 == 1) {
                E_direct_total_energy += E1;
                if (E1 > E_direct_max_energy) E_direct_max_energy = E1;
                if (E1 > 80*units::MeV) ++n_direct_showers;
            } else if (conn1 == 2) {
                E_indirect_total_energy += E1;
                if (E1 > E_indirect_max_energy) E_indirect_max_energy = E1;
                if (E1 > 80*units::MeV) ++n_indirect_showers;
            }
        }
        (void)E_direct_total_energy; (void)n_direct_showers;
        (void)E_indirect_total_energy; (void)n_indirect_showers;

        // 7049_1070_53534
        if (flag_single_shower && std::max(vec_dQ_dx[0], vec_dQ_dx[1]) > 1.6 &&
            std::max(vec_dQ_dx[0], vec_dQ_dx[1]) < 3.5 && Eshower < 350*units::MeV &&
            E_indirect_max_energy > 70*units::MeV) mip_id = -1;
        // 7012_1450_72525
        if (flag_single_shower && E_indirect_max_energy > 0.33 * Eshower && mip_id == 1)
            mip_id = 0;
        // 7023_28_1419
        if (mip_id == 0 && Eshower < 250*units::MeV &&
            sg->flags_any(SegmentFlags::kShowerTrajectory) &&
            segment_track_length(sg) < 5*units::cm) mip_id = -1;

        // min_dis: minimum distance from non-main shower vertices to main-cluster shower vertices
        double length1 = sg->cluster() ? shower->get_total_length(sg->cluster()) : 0;
        double length2 = shower->get_total_length();
        double min_dis = 1e9;
        {
            IndexedSegmentSet sh_segs;
            IndexedVertexSet  sh_vtxs;
            shower->fill_sets(sh_vtxs, sh_segs, false);

            Facade::Cluster* sg_cl = sg->cluster();
            std::vector<VertexPtr> main_cl_vtxs;
            for (VertexPtr v2 : sh_vtxs)
                if (v2->cluster() == sg_cl) main_cl_vtxs.push_back(v2);

            for (SegmentPtr sg1 : sh_segs) {
                if (sg1->cluster() == sg_cl) continue;
                if (segment_track_length(sg1) < 3*units::cm) continue;
                auto [va, vb] = find_vertices(ctx.graph, sg1);
                for (VertexPtr vtx1 : {va, vb}) {
                    if (!vtx1) continue;
                    Point p1 = vtx_fit_pt(vtx1);
                    for (VertexPtr vtx2 : main_cl_vtxs) {
                        double d = ray_length(Ray{p1, vtx_fit_pt(vtx2)});
                        if (d < min_dis) min_dis = d;
                    }
                }
            }
        }

        if (mip_id == 1) {
            // 7012_1646_82342
            if (length1 < 0.1 * length2 && length1 < 10*units::cm && min_dis > 8*units::cm)
                mip_id = 0;
        }

        // n_other_vertex: degree of the far-end vertex of sg
        int n_other_vertex = 0;
        VertexPtr other_vertex = find_other_vertex(ctx.graph, sg, vertex);
        if (other_vertex) {
            int other_deg = vtx_degree(other_vertex, ctx.graph);
            if (other_deg > 2 &&
                (int)sg->fits().size() <= n_first_mip + 1 && n_first_mip > 2)
                mip_id = -1;
            n_other_vertex = other_deg;
        }

        // Median dQ/dx quality check
        double medium_dQ_dx = 1;
        {
            std::vector<double> tmp = vec_dQ_dx;
            std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
            medium_dQ_dx = *std::next(tmp.begin(), tmp.size()/2);
        }
        if (medium_dQ_dx < 0.75 && Eshower < 150*units::MeV) mip_id = -1;

        // Pad to ≥20 entries then fill shw_sp_vec_dQ_dx_* (0–19)
        while (vec_dQ_dx.size() < 20) vec_dQ_dx.push_back(3);

        // Median dedx from first 7 dQ/dx values
        std::vector<float> dqdx7(vec_dQ_dx.begin(), vec_dQ_dx.begin() + 7);
        std::sort(dqdx7.begin(), dqdx7.end());
        size_t mid = dqdx7.size() / 2;
        float median_dqdx = (dqdx7.size() % 2 == 0)
                            ? (dqdx7[mid] + dqdx7[mid-1]) / 2.0f
                            : dqdx7[mid];
        float mean_dqdx = 0;
        for (float d : dqdx7) mean_dqdx += d;
        mean_dqdx /= (float)dqdx7.size();

        const float alpha = 1.0f, beta = 0.255f;
        float median_dedx = (std::exp(median_dqdx * 43e3 * 23.6e-6 * beta / 1.38f / 0.273f) - alpha)
                            / (beta / 1.38f / 0.273f);
        if (median_dedx < 0) median_dedx = 0;
        if (median_dedx > 50) median_dedx = 50;
        float mean_dedx   = (std::exp(mean_dqdx   * 43e3 * 23.6e-6 * beta / 1.38f / 0.273f) - alpha)
                            / (beta / 1.38f / 0.273f);
        if (mean_dedx < 0) mean_dedx = 0;
        if (mean_dedx > 50) mean_dedx = 50;

        // Fill TaggerInfo (unconditional in toolkit)
        if (min_dis > 1000*units::cm) min_dis = 1000*units::cm;
        ti.shw_sp_flag = (mip_id == 1) ? 0.0f : 1.0f;
        ti.shw_sp_energy            = Eshower / units::MeV;
        ti.shw_sp_max_dQ_dx_sample  = max_dQ_dx_sample;
        ti.shw_sp_vec_dQ_dx_0       = vec_dQ_dx[0];
        ti.shw_sp_vec_dQ_dx_1       = vec_dQ_dx[1];
        ti.shw_sp_n_below_threshold = n_below_threshold;
        ti.shw_sp_n_good_tracks     = n_good_tracks;
        ti.shw_sp_n_vertex          = vtx_deg;
        ti.shw_sp_angle_beam        = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm)
                                      .angle(dir_beam) / M_PI * 180.0;
        ti.shw_sp_flag_all_above    = flag_all_above;
        ti.shw_sp_length_main       = (sg->cluster() ? shower->get_total_length(sg->cluster()) : 0) / units::cm;
        ti.shw_sp_length_total      = shower->get_total_length() / units::cm;
        ti.shw_sp_min_dQ_dx_5       = min_dQ_dx_5;
        ti.shw_sp_lowest_dQ_dx      = lowest_dQ_dx;
        ti.shw_sp_iso_angle         = iso;
        ti.shw_sp_n_below_zero      = n_below_zero;
        ti.shw_sp_highest_dQ_dx     = highest_dQ_dx;
        ti.shw_sp_n_lowest          = n_lowest;
        ti.shw_sp_n_highest         = n_highest;
        ti.shw_sp_stem_length       = segment_track_length(sg) / units::cm;
        ti.shw_sp_E_indirect_max_energy = E_indirect_max_energy / units::MeV;
        ti.shw_sp_flag_stem_trajectory  = sg->flags_any(SegmentFlags::kShowerTrajectory) ? 1.0f : 0.0f;
        ti.shw_sp_min_dis           = min_dis / units::cm;
        ti.shw_sp_n_other_vertex    = n_other_vertex;
        ti.shw_sp_n_stem_size       = (float)sg->fits().size();
        ti.shw_sp_medium_dQ_dx      = medium_dQ_dx;
        ti.shw_sp_filled            = 1;
        ti.shw_sp_vec_dQ_dx_2  = vec_dQ_dx[2];
        ti.shw_sp_vec_dQ_dx_3  = vec_dQ_dx[3];
        ti.shw_sp_vec_dQ_dx_4  = vec_dQ_dx[4];
        ti.shw_sp_vec_dQ_dx_5  = vec_dQ_dx[5];
        ti.shw_sp_vec_dQ_dx_6  = vec_dQ_dx[6];
        ti.shw_sp_vec_dQ_dx_7  = vec_dQ_dx[7];
        ti.shw_sp_vec_dQ_dx_8  = vec_dQ_dx[8];
        ti.shw_sp_vec_dQ_dx_9  = vec_dQ_dx[9];
        ti.shw_sp_vec_dQ_dx_10 = vec_dQ_dx[10];
        ti.shw_sp_vec_dQ_dx_11 = vec_dQ_dx[11];
        ti.shw_sp_vec_dQ_dx_12 = vec_dQ_dx[12];
        ti.shw_sp_vec_dQ_dx_13 = vec_dQ_dx[13];
        ti.shw_sp_vec_dQ_dx_14 = vec_dQ_dx[14];
        ti.shw_sp_vec_dQ_dx_15 = vec_dQ_dx[15];
        ti.shw_sp_vec_dQ_dx_16 = vec_dQ_dx[16];
        ti.shw_sp_vec_dQ_dx_17 = vec_dQ_dx[17];
        ti.shw_sp_vec_dQ_dx_18 = vec_dQ_dx[18];
        ti.shw_sp_vec_dQ_dx_19 = vec_dQ_dx[19];
        ti.shw_sp_vec_median_dedx = median_dedx;
        ti.shw_sp_vec_mean_dedx   = mean_dedx;
    }

    return mip_id;
}

// ===========================================================================
// high_energy_overlapping_sp
//
// Two sub-checks for high-energy shower overlap:
//   hol_1: other showers/tracks at vertex with small opening angle or all-shower
//           topology suggest overlapping clusters.
//   hol_2: consecutive stem fit points within 0.6 cm of another segment at vtx
//           with a matching dQ/dx.
//
// Fills TaggerInfo shw_sp_hol_* fields (unconditionally).
//
// Prototype: WCPPID::NeutrinoID::high_energy_overlapping_sp()
//            NeutrinoID_singlephoton_tagger.h line 2596.
//
// Logic identical to high_energy_overlapping() in NeutrinoTaggerNuE.cxx;
// only TaggerInfo field names differ (shw_sp_hol_* vs hol_*).
// ===========================================================================
static bool high_energy_overlapping_sp(SpContext& ctx, ShowerPtr shower, TaggerInfo& ti)
{
    bool flag_overlap1 = false;
    bool flag_overlap2 = false;

    double Eshower = shower_energy(shower);

    auto [vtx, conn_type] = shower->get_start_vertex_and_type();
    SegmentPtr sg = shower->start_segment();

    auto vec_dQ_dx = shower->get_stem_dQ_dx(vtx, sg, 20);
    double max_dQ_dx = 0;
    for (size_t i = 0; i < vec_dQ_dx.size(); ++i) {
        if (vec_dQ_dx[i] > max_dQ_dx) max_dQ_dx = vec_dQ_dx[i];
        if ((int)i == 2) break;
    }

    const auto& sg_fits = sg->fits();
    bool flag_start = !sg_fits.empty() &&
        (ray_length(Ray{vtx_fit_pt(vtx), sg_fits.front().point}) <=
         ray_length(Ray{vtx_fit_pt(vtx), sg_fits.back().point}));
    Point vtx_point = flag_start ? sg_fits.front().point : sg_fits.back().point;

    if (conn_type == 1 && vtx && vtx->descriptor_valid()) {
        // ------------------------------------------------------------------
        // hol_1: n_valid_tracks and min_angle at vertex
        // 7012_1195_59764 + 7017_1158_57929
        // ------------------------------------------------------------------
        Point vp = vtx_point;
        Vector dir1 = segment_cal_dir_3vector(sg, vp, 15*units::cm);
        int    n_valid_tracks = 0;
        double min_angle      = 180;
        double min_length     = 0;
        bool   flag_all_showers = true;

        auto vd = vtx->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            bool is_pdg11 = sg1->has_particle_info() && sg1->particle_info()->pdg() == 11;
            bool is_weak_muon = sg1->has_particle_info() &&
                                sg1->particle_info()->pdg() == 13 &&
                                sg1->dir_weak() &&
                                segment_track_length(sg1) < 6*units::cm;
            if (is_pdg11 || is_weak_muon) {
                Point dp = vtx_point;
                Vector dir2 = segment_cal_dir_3vector(sg1, dp, 5*units::cm);
                if (dir2.magnitude() == 0) { flag_all_showers = false; continue; }
                double angle = dir1.angle(dir2) / M_PI * 180.0;
                if (angle < min_angle) {
                    min_angle = angle;
                    min_length = segment_track_length(sg1);
                }
            } else {
                flag_all_showers = false;
            }
            double norm_dQ = segment_median_dQ_dx(sg1) / (43e3 / units::cm);
            bool is_proton = sg1->has_particle_info() && sg1->particle_info()->pdg() == 2212;
            if ((!sg1->dir_weak() || is_proton || segment_track_length(sg1) > 20*units::cm) &&
                !seg_is_shower(sg1))
                ++n_valid_tracks;
            else if (norm_dQ > 2.0 && segment_track_length(sg1) > 1.8*units::cm)
                ++n_valid_tracks; // 7010_20_1012
        }

        int num_showers = 0;
        auto mv_it = ctx.map_vertex_to_shower.find(vtx);
        if (mv_it != ctx.map_vertex_to_shower.end()) {
            for (ShowerPtr shower1 : mv_it->second) {
                if (shower1 == shower) continue;
                SegmentPtr sg1 = shower1->start_segment();
                if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11)
                    continue;
                auto [vtx1, conn1] = shower1->get_start_vertex_and_type();
                double E1 = shower_energy(shower1);
                if (conn1 == 1 && E1 > 250*units::MeV) ++n_valid_tracks;
                if (Eshower > 60*units::MeV && conn1 < 3) ++num_showers;
            }
        }
        // 6689_127_6366
        if (max_dQ_dx > 3.6 && n_valid_tracks == 0) ++n_valid_tracks;
        if (max_dQ_dx > 2.8 && num_showers >= 2 && min_angle > 40) ++n_valid_tracks;

        if (n_valid_tracks == 0 && min_angle < 30 && Eshower < 1500*units::MeV)
            flag_overlap1 = true;
        if (n_valid_tracks == 0 && min_angle < 60 && flag_all_showers &&
            Eshower < 300*units::MeV && Eshower < 1500*units::MeV)
            flag_overlap1 = true;
        if (n_valid_tracks == 0 && min_angle < 60 && flag_all_showers &&
            Eshower < 800*units::MeV && min_length < 5*units::cm && Eshower < 1500*units::MeV)
            flag_overlap1 = true;

        ti.shw_sp_hol_1_n_valid_tracks  = n_valid_tracks;
        ti.shw_sp_hol_1_min_angle       = min_angle;
        ti.shw_sp_hol_1_energy          = Eshower / units::MeV;
        ti.shw_sp_hol_1_flag_all_shower = flag_all_showers;
        ti.shw_sp_hol_1_min_length      = min_length / units::cm;
        ti.shw_sp_hol_1_flag            = !flag_overlap1;

        // ------------------------------------------------------------------
        // hol_2: ncount consecutive close fit points + dQ/dx of min_sg
        // ------------------------------------------------------------------
        {
            Point dp2 = vtx_point;
            Vector dir1_8 = segment_cal_dir_3vector(sg, dp2, 8*units::cm);
            double min_ang2 = 180;
            SegmentPtr min_sg = nullptr;

            for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
                SegmentPtr sg1 = ctx.graph[*eit].segment;
                if (!sg1 || sg1 == sg) continue;
                Point dp3 = vtx_point;
                Vector dir2 = segment_cal_dir_3vector(sg1, dp3, 5*units::cm);
                if (dir2.magnitude() == 0) continue;
                double ang = dir1_8.angle(dir2) / M_PI * 180.0;
                if (ang < min_ang2) { min_ang2 = ang; min_sg = sg1; }
            }

            int ncount = 0;
            auto iterate_pts = [&](auto begin_it, auto end_it) {
                for (auto it1 = begin_it; it1 != end_it; ++it1) {
                    double min_dis = 1e9;
                    for (auto [eit2, eend2] = boost::out_edges(vd, ctx.graph);
                         eit2 != eend2; ++eit2) {
                        SegmentPtr sg2 = ctx.graph[*eit2].segment;
                        if (!sg2 || sg2 == sg) continue;
                        double dis = segment_get_closest_point(sg2, it1->point).first;
                        if (dis < min_dis) min_dis = dis;
                    }
                    if (min_dis < 0.6*units::cm) ++ncount;
                    else break;
                }
            };
            if (flag_start) iterate_pts(sg_fits.begin(),  sg_fits.end());
            else            iterate_pts(sg_fits.rbegin(), sg_fits.rend());

            double medium_dQ_dx = 0;
            if (min_sg) {
                const auto& min_fits = min_sg->fits();
                int n_min = (int)min_fits.size();
                bool min_front_near =
                    (ray_length(Ray{vtx_point, min_fits.front().point}) <=
                     ray_length(Ray{vtx_point, min_fits.back().point}));
                if (min_front_near)
                    medium_dQ_dx = segment_median_dQ_dx(min_sg, 0, ncount) / (43e3/units::cm);
                else
                    medium_dQ_dx = segment_median_dQ_dx(min_sg, n_min-1-ncount, n_min-1) / (43e3/units::cm);
            }

            if (min_ang2 < 15 && medium_dQ_dx > 0.95 && ncount > 5  && Eshower < 1500*units::MeV)
                flag_overlap2 = true;
            if (min_ang2 < 7.5 && medium_dQ_dx > 0.8  && ncount > 8  && Eshower < 1500*units::MeV)
                flag_overlap2 = true;
            if (min_ang2 < 5   && ncount > 12 && medium_dQ_dx > 0.5  && Eshower < 1500*units::MeV)
                flag_overlap2 = true;

            ti.shw_sp_hol_2_min_angle    = min_ang2;
            ti.shw_sp_hol_2_medium_dQ_dx = medium_dQ_dx;
            ti.shw_sp_hol_2_ncount       = ncount;
            ti.shw_sp_hol_2_energy       = Eshower / units::MeV;
            ti.shw_sp_hol_2_flag         = !flag_overlap2;
        }
    }

    bool flag_overlap = flag_overlap1 || flag_overlap2;
    ti.shw_sp_hol_flag = !flag_overlap;
    return flag_overlap;
}

// ===========================================================================
// low_energy_overlapping_sp
//
// Three complementary checks for low-energy shower overlap:
//   lol_1: two shower-internal segments at a main-cluster vertex form
//     a small opening angle.
//   lol_2: short collinear muon/weak segment at vertex parallel to stem.
//   lol_3: backward-going shower with no valid tracks, or isolated
//     shower cluster with too many outward-pointing hits.
//
// Fills TaggerInfo shw_sp_lol_* fields (unconditionally).
//
// Prototype: WCPPID::NeutrinoID::low_energy_overlapping_sp()
//            NeutrinoID_singlephoton_tagger.h line 2797.
//
// Logic identical to low_energy_overlapping() in NeutrinoTaggerNuE.cxx;
// only TaggerInfo field names differ (shw_sp_lol_* vs lol_*).
// ===========================================================================
static bool low_energy_overlapping_sp(SpContext& ctx, ShowerPtr shower, TaggerInfo& ti)
{
    bool flag_overlap_1_save = false;
    bool flag_overlap_2_save = false;
    bool flag_overlap_3      = false;

    Vector dir_beam(0, 0, 1);
    double Eshower = shower_energy(shower);

    auto [vtx, conn_type] = shower->get_start_vertex_and_type();
    SegmentPtr sg = shower->start_segment();
    Point vtx_point = seg_endpoint_near(sg, vtx_fit_pt(vtx));

    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    int nseg = 0;
    for (SegmentPtr sg1 : shower_segs)
        if (sg1->cluster() == sg->cluster()) ++nseg;

    Point vp_copy = vtx_point;
    Vector dir1_stem = segment_cal_dir_3vector(sg, vp_copy, 5*units::cm);
    double angle_beam = dir1_stem.angle(dir_beam) / M_PI * 180.0;

    int    n_valid_tracks = 0;
    double min_angle_vtx  = 180;
    size_t n_vtx_segs_global = (vtx && vtx->descriptor_valid())
                               ? boost::out_degree(vtx->get_descriptor(), ctx.graph) : 0;

    if (vtx && vtx->descriptor_valid()) {
        auto vd = vtx->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            Point dp = vtx_point;
            Vector dir2 = segment_cal_dir_3vector(sg1, dp, 5*units::cm);
            double tmp_angle = dir2.angle(dir1_stem) / M_PI * 180.0;
            if (tmp_angle < min_angle_vtx) min_angle_vtx = tmp_angle;
            bool is_proton = sg1->has_particle_info() && sg1->particle_info()->pdg() == 2212;
            if ((!sg1->dir_weak() || is_proton || segment_track_length(sg1) > 20*units::cm) &&
                !seg_is_shower(sg1))
                ++n_valid_tracks;
        }
    }

    // n_out / n_sum: shower hits outside 15° cone
    int n_sum = 0, n_out = 0;
    {
        Point vp2 = vtx_point;
        Vector dir1_15 = segment_cal_dir_3vector(sg, vp2, 15*units::cm);
        for (VertexPtr vtx1 : shower_vtxs) {
            Vector dir2 = vtx_fit_pt(vtx1) - vtx_point;
            ++n_sum;
            if (dir1_15.angle(dir2) / M_PI * 180.0 > 15) ++n_out;
        }
        for (SegmentPtr sg1 : shower_segs) {
            const auto& fits1 = sg1->fits();
            for (size_t i = 1; i + 1 < fits1.size(); ++i) {
                Vector dir2 = fits1[i].point - vtx_point;
                ++n_sum;
                if (dir1_15.angle(dir2) / M_PI * 180.0 > 15) ++n_out;
            }
        }
    }

    // ------------------------------------------------------------------
    // lol_1: shower-internal vertices with 2 connected shower segs
    //        forming a small opening angle.
    // ------------------------------------------------------------------
    for (VertexPtr vtx1 : shower_vtxs) {
        if (!vtx1->cluster() || vtx1->cluster() != sg->cluster()) continue;
        if (!vtx1->descriptor_valid()) continue;

        std::vector<SegmentPtr> vtx_ss;
        for (auto [eit, eend] = boost::out_edges(vtx1->get_descriptor(), ctx.graph);
             eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (sg1 && shower_segs.count(sg1)) vtx_ss.push_back(sg1);
        }
        if (vtx_ss.empty()) continue;

        Point vp1 = vtx_fit_pt(vtx1);
        Point d1p = vp1, d2p = vp1;
        Vector dv1 = segment_cal_dir_3vector(vtx_ss.front(), d1p, 5*units::cm);
        Vector dv2 = segment_cal_dir_3vector(vtx_ss.back(),  d2p, 5*units::cm);
        double open_angle = dv1.angle(dv2) / M_PI * 180.0;

        bool flag_ov1 = false;
        if (vtx_ss.size() == 2 &&
            open_angle < 36 && nseg == 2 && Eshower < 150*units::MeV &&
            n_vtx_segs_global == 1)
            flag_ov1 = true;

        ti.shw_sp_lol_1_v_energy.push_back(Eshower / units::MeV);
        ti.shw_sp_lol_1_v_vtx_n_segs.push_back(n_vtx_segs_global);
        ti.shw_sp_lol_1_v_nseg.push_back(nseg);
        ti.shw_sp_lol_1_v_angle.push_back(open_angle);
        ti.shw_sp_lol_1_v_flag.push_back(!flag_ov1);

        if (flag_ov1) flag_overlap_1_save = true;
    }

    // ------------------------------------------------------------------
    // lol_2: short collinear muon/weak segment at vertex.
    // ------------------------------------------------------------------
    if (vtx && vtx->descriptor_valid()) {
        auto vd = vtx->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            bool flag_ov2 = false;
            Point dp = vtx_point;
            Vector dir2 = segment_cal_dir_3vector(sg1, dp, 5*units::cm);
            double ang2 = dir2.angle(dir1_stem) / M_PI * 180.0;
            double len1 = segment_track_length(sg1);
            bool is_muon = sg1->has_particle_info() && sg1->particle_info()->pdg() == 13;
            double main_len = sg->cluster() ? shower->get_total_length(sg->cluster()) : 0;

            // 7017_1604_80242
            if (((len1 < 30*units::cm && ang2 < 10) ||
                 (len1 < 7.5*units::cm && ang2 < 17.5)) && is_muon &&
                n_vtx_segs_global > 1 && Eshower < 300*units::MeV && main_len < 20*units::cm)
                flag_ov2 = true;
            // 7020_249_12479
            if (sg1->dir_weak() && len1 < 8*units::cm && ang2 < 30 &&
                n_vtx_segs_global == 2 && Eshower < 400*units::MeV)
                flag_ov2 = true;

            ti.shw_sp_lol_2_v_flag.push_back(!flag_ov2);
            ti.shw_sp_lol_2_v_length.push_back(len1 / units::cm);
            ti.shw_sp_lol_2_v_angle.push_back(ang2);
            ti.shw_sp_lol_2_v_type.push_back(sg1->has_particle_info()
                                              ? (float)sg1->particle_info()->pdg() : 0.0f);
            ti.shw_sp_lol_2_v_vtx_n_segs.push_back(n_vtx_segs_global);
            ti.shw_sp_lol_2_v_energy.push_back(Eshower / units::MeV);
            ti.shw_sp_lol_2_v_shower_main_length.push_back(main_len / units::cm);
            ti.shw_sp_lol_2_v_flag_dir_weak.push_back(sg1->dir_weak());

            if (flag_ov2) flag_overlap_2_save = true;
        }
    }

    // ------------------------------------------------------------------
    // lol_3: backward/isolated low-energy shower.
    // ------------------------------------------------------------------
    double main_length = sg->cluster() ? shower->get_total_length(sg->cluster()) : 0;

    // 7010_194_9750
    if (angle_beam > 60 && n_valid_tracks == 0 && min_angle_vtx < 80 &&
        n_vtx_segs_global > 1 && Eshower < 300*units::MeV && main_length < 20*units::cm)
        flag_overlap_3 = true;
    if (n_vtx_segs_global == 1 && main_length < 15*units::cm &&
        Eshower > 30*units::MeV && Eshower < 250*units::MeV &&
        n_out > n_sum / 3)
        flag_overlap_3 = true;

    ti.shw_sp_lol_3_flag               = !flag_overlap_3;
    ti.shw_sp_lol_3_angle_beam         = angle_beam;
    ti.shw_sp_lol_3_min_angle          = min_angle_vtx;
    ti.shw_sp_lol_3_n_valid_tracks     = n_valid_tracks;
    ti.shw_sp_lol_3_vtx_n_segs         = n_vtx_segs_global;
    ti.shw_sp_lol_3_energy             = Eshower / units::MeV;
    ti.shw_sp_lol_3_shower_main_length = main_length / units::cm;
    ti.shw_sp_lol_3_n_sum              = n_sum;
    ti.shw_sp_lol_3_n_out              = n_out;

    bool flag_overlap = flag_overlap_1_save || flag_overlap_2_save || flag_overlap_3;
    ti.shw_sp_lol_flag = !flag_overlap;
    return flag_overlap;
}

// ===========================================================================
// pi0_identification_sp
//
// Two sub-checks for pi0 identification:
//   pio_1: shower belongs to a reconstructed pi0 pair (from map_shower_pio_id);
//           check whether the pair kinematics are consistent with pi0 decay.
//   pio_2: shower is NOT in the pi0 map; look for a back-to-back cluster
//           that could be the other photon.
//
// Fills TaggerInfo shw_sp_pio_* fields (unconditionally).
// Returns pi0_flag_pio = (shower IS in map_shower_pio_id), not flag_pi0.
//
// Prototype: WCPPID::NeutrinoID::pi0_identification_sp()
//            NeutrinoID_singlephoton_tagger.h line 2955.
//
// Logic closely follows pi0_identification() in NeutrinoTaggerNuE.cxx;
// differences: field prefix shw_sp_pio_*, separate pi0_flag_pio return,
// and ti.shw_sp_pio_flag / ti.shw_sp_pio_1_flag set explicitly.
// ===========================================================================
static bool pi0_identification_sp(SpContext& ctx,
                                   VertexPtr vertex, SegmentPtr sg, ShowerPtr shower,
                                   double threshold, TaggerInfo& ti)
{
    bool flag_pi0_1 = false;
    bool flag_pi0_2 = false;

    Point vertex_point = vtx_fit_pt(vertex);

    // Collect all vertices that already belong to known pi0-paired showers.
    IndexedVertexSet used_vertices;
    for (auto& [shower1, pio_id] : ctx.map_shower_pio_id) {
        IndexedVertexSet vtxs;
        IndexedSegmentSet segs;
        shower1->fill_sets(vtxs, segs, false);
        used_vertices.insert(vtxs.begin(), vtxs.end());
    }

    auto it = ctx.map_shower_pio_id.find(shower);
    bool pi0_flag_pio = (it != ctx.map_shower_pio_id.end());
    ti.shw_sp_pio_flag_pio = pi0_flag_pio ? 1.0f : 0.0f;

    if (it != ctx.map_shower_pio_id.end()) {
        // ----------------------------------------------------------------
        // pio_1: shower belongs to a reconstructed pi0 pair.
        // ----------------------------------------------------------------
        auto& tmp_pi0_showers = ctx.map_pio_id_showers[it->second];
        auto  mass_pair       = ctx.map_pio_id_mass[it->second];

        double Eshower_1 = tmp_pi0_showers.front()->get_kine_charge();
        double Eshower_2 = tmp_pi0_showers.back()->get_kine_charge();

        double dis1 = ray_length(Ray{tmp_pi0_showers.front()->get_start_point(), vertex_point});
        double dis2 = ray_length(Ray{tmp_pi0_showers.back()->get_start_point(),  vertex_point});

        ti.shw_sp_pio_1_mass      = mass_pair.first  / units::MeV;
        ti.shw_sp_pio_1_pio_type  = mass_pair.second;
        ti.shw_sp_pio_1_energy_1  = Eshower_1 / units::MeV;
        ti.shw_sp_pio_1_energy_2  = Eshower_2 / units::MeV;
        ti.shw_sp_pio_1_dis_1     = dis1 / units::cm;
        ti.shw_sp_pio_1_dis_2     = dis2 / units::cm;

        bool mass_ok_1 = (std::fabs(mass_pair.first - 135*units::MeV) < 35*units::MeV &&
                          mass_pair.second == 1);
        bool mass_ok_2 = (std::fabs(mass_pair.first - 135*units::MeV) < 60*units::MeV &&
                          mass_pair.second == 2);

        if (mass_ok_1 || mass_ok_2) {
            if (std::min(Eshower_1, Eshower_2) > 15*units::MeV &&
                std::fabs(Eshower_1 - Eshower_2) / (Eshower_1 + Eshower_2) < 0.87)
                flag_pi0_1 = true;
            // 6058_43_2166, 7017_364_18210
            if (std::min(Eshower_1, Eshower_2) > std::max(10*units::MeV, threshold) &&
                std::max(Eshower_1, Eshower_2) < 400*units::MeV)
                flag_pi0_1 = true;

            // Veto: asymmetric pair with large separation (7049_875_43775)
            if (flag_pi0_1) {
                bool veto_1 = (std::min(Eshower_1, Eshower_2) < 30*units::MeV &&
                               std::max(dis1, dis2) > 80*units::cm &&
                               std::fabs(Eshower_1 - Eshower_2) / (Eshower_1 + Eshower_2) > 0.87 &&
                               std::min(dis1, dis2) == 0);
                bool veto_2 = (std::min(Eshower_1, Eshower_2) < 30*units::MeV &&
                               std::max(dis1, dis2) > 120*units::cm &&
                               std::fabs(Eshower_1 - Eshower_2) / (Eshower_1 + Eshower_2) > 0.80 &&
                               std::min(dis1, dis2) == 0);
                if (veto_1 || veto_2) flag_pi0_1 = false;
            }
        }

    } else {
        // ----------------------------------------------------------------
        // pio_2: shower not in pi0 map — look for back-to-back cluster.
        // ----------------------------------------------------------------
        Point sg_vp = vertex_point;
        Vector dir1 = segment_cal_dir_3vector(sg, sg_vp, 12*units::cm);

        if (dir1.magnitude() > 0) {
            // Precompute total track length per cluster
            std::map<Facade::Cluster*, double> cluster_acc_length;
            for (auto [eit, eend] = boost::edges(ctx.graph); eit != eend; ++eit) {
                SegmentPtr sg1 = ctx.graph[*eit].segment;
                if (sg1 && sg1->cluster())
                    cluster_acc_length[sg1->cluster()] += segment_track_length(sg1);
            }

            for (const auto& vd : graph_nodes(ctx.graph)) {
                VertexPtr vtx1 = ctx.graph[vd].vertex;
                if (!vtx1) continue;
                if (vtx1->cluster() == vertex->cluster()) continue;
                if (used_vertices.count(vtx1)) continue;

                double acc_length = 0;
                auto cl_it = cluster_acc_length.find(vtx1->cluster());
                if (cl_it != cluster_acc_length.end()) acc_length = cl_it->second;

                Point vtx1_pt = vtx_fit_pt(vtx1);
                Vector dir2 = vtx1_pt - vertex_point;
                double dis2 = dir2.magnitude();
                if (dis2 <= 0) continue;

                double back_angle = 180.0 - dir1.angle(dir2) / M_PI * 180.0;

                if (dis2 < 36*units::cm && back_angle < 7.5 && acc_length > 0) {
                    flag_pi0_2 = true;
                    ti.shw_sp_pio_2_v_flag.push_back(0.0f);
                    ti.shw_sp_pio_2_v_dis2.push_back(dis2 / units::cm);
                    ti.shw_sp_pio_2_v_angle2.push_back(back_angle);
                    ti.shw_sp_pio_2_v_acc_length.push_back(acc_length / units::cm);
                } else {
                    ti.shw_sp_pio_2_v_flag.push_back(1.0f);
                    ti.shw_sp_pio_2_v_dis2.push_back(dis2 / units::cm);
                    ti.shw_sp_pio_2_v_angle2.push_back(back_angle);
                    ti.shw_sp_pio_2_v_acc_length.push_back(acc_length / units::cm);
                }
            }
        }
    }

    ti.shw_sp_pio_1_flag = !flag_pi0_1;
    ti.shw_sp_pio_flag   = (flag_pi0_1 || flag_pi0_2) ? 0.0f : 1.0f;

    return pi0_flag_pio;
}
