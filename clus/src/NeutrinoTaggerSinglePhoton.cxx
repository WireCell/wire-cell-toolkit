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
