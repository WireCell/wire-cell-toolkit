// NeutrinoTaggerCosmic.cxx
//
// Ports of:
//   WCPPID::NeutrinoID::bad_reconstruction()  (NeutrinoID_nue_tagger.h:3450)
//   WCPPID::NeutrinoID::cosmic_tagger()       (NeutrinoID_cosmic_tagger.h)
//
// Namespace/class: WireCell::Clus::PR::PatternAlgorithms
//
// Translation conventions (see neutrino_id_function_map.md):
//   main_vertex->get_fit_pt()         → vtx_pt(vtx)  [lambda defined locally]
//   map_vertex_segments[vtx]          → boost::out_edges(vtx->get_descriptor(), graph)
//   sg->get_flag_shower()             → seg->flags_any(kShowerTrajectory) || flags_any(kShowerTopology) || abs(pdg)==11
//   sg->get_length()                  → segment_track_length(sg)
//   sg->get_medium_dQ_dx(n1,n2)       → segment_median_dQ_dx(sg, n1, n2)
//   sg->cal_dir_3vector(pt, dis)      → segment_cal_dir_3vector(sg, pt, dis)
//   sg->get_particle_type()           → sg->particle_info()->pdg()
//   shower->get_start_vertex()        → shower->get_start_vertex_and_type()
//   shower->get_start_segment()       → shower->start_segment()
//   shower->cal_dir_3vector(pt, dis)  → shower_cal_dir_3vector(*shower, pt, dis)
//   find_other_vertex(sg, vtx)        → find_other_vertex(graph, sg, vtx)
//   map_cluster_length[cl]            → cl->get_length()
//   main_cluster->Calc_PCA(pts)
//     + get_PCA_axis(0)               → calc_PCA_main_axis(pts).second
//   fid->inside_fiducial_volume(p,..) → fiducial_utils->inside_fiducial_volume(p)

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"
#include <cmath>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;
using namespace WireCell;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// ---------------------------------------------------------------------------
// Helpers: spherical angles for a direction vector.
// D3Vector has no theta()/phi() methods; compute explicitly.
//   theta = polar angle from +z  (range [0, pi])
//   phi   = azimuthal angle      (range (-pi, pi])
// ---------------------------------------------------------------------------
static inline double vec_theta(const Vector& v) {
    double m = v.magnitude();
    return (m > 0) ? std::acos(v.z() / m) : 0.0;
}
static inline double vec_phi(const Vector& v) {
    return std::atan2(v.y(), v.x());
}

// ---------------------------------------------------------------------------
// Helper: vertex fit point (with wcpt fallback), matching pattern in
// NeutrinoVertexFinder.cxx.
// ---------------------------------------------------------------------------
static inline Point vtx_fit_pt(VertexPtr v) {
    if (!v) return Point{};
    return v->fit().valid() ? v->fit().point : v->wcpt().point;
}

// ---------------------------------------------------------------------------
// Helper: number of edges (segments) at a vertex in the main graph.
// Prototype: map_vertex_segments[vtx].size()
// ---------------------------------------------------------------------------
static inline int vtx_degree(VertexPtr vtx, const Graph& graph) {
    if (!vtx || !vtx->descriptor_valid()) return 0;
    return static_cast<int>(boost::out_degree(vtx->get_descriptor(), graph));
}

// ---------------------------------------------------------------------------
// Helper: is this segment a "shower" segment?
// Prototype: sg->get_flag_shower() = kShowerTrajectory || kShowerTopology || abs(pdg)==11
// ---------------------------------------------------------------------------
static inline bool seg_is_shower(SegmentPtr seg) {
    return seg->flags_any(SegmentFlags::kShowerTrajectory) ||
           seg->flags_any(SegmentFlags::kShowerTopology)   ||
           (seg->has_particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
}

// shower_energy helper was removed.
// PRShower::get_kine_best() already falls back to kenergy_charge when kenergy_best==0,
// so shower->get_kine_best() is the correct and sufficient call everywhere.

// ===========================================================================
// bad_reconstruction
//
// Determines if a shower is likely a mis-identified muon/track by checking:
//   flag_bad_shower_1 : stem too long, or un-shower-like stem with low energy
//   flag_bad_shower_2 : long muon-like track inside shower (using find_cont_muon_segment_nue)
//   flag_bad_shower_3 : long muon-like track near start segment (nearby segment geometry)
//
// Prototype: WCPPID::NeutrinoID::bad_reconstruction in NeutrinoID_nue_tagger.h:3450.
// flag_fill: if true, fill tagger_info BDT features; ti must be non-null in that case.
// ===========================================================================
bool PatternAlgorithms::bad_reconstruction(
    Graph& graph,
    VertexPtr main_vertex,
    ShowerPtr shower,
    bool flag_fill,
    TaggerInfo* ti)
{
    bool flag_bad_shower_1 = false;
    bool flag_bad_shower_2 = false;
    bool flag_bad_shower_3 = false;

    double Eshower = shower->get_kine_best();

    SegmentPtr sg = shower->start_segment();
    if (!sg) return false;

    auto [vtx, start_type] = shower->get_start_vertex_and_type();

    // -----------------------------------------------------------------------
    // Collect all segments and vertices in this shower.
    // Replaces prototype's map_seg_vtxs / map_vtx_segs.
    // -----------------------------------------------------------------------
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    // -----------------------------------------------------------------------
    // Sub-check 1: stem characteristics.
    // -----------------------------------------------------------------------
    {
        double sg_length = segment_track_length(sg);

        if (start_type == 1 && vtx_degree(vtx, graph) == 1 &&
            Eshower < 120 * units::MeV && (int)shower_segs.size() <= 3) {
            bool topo = sg->flags_any(SegmentFlags::kShowerTopology);
            bool traj = sg->flags_any(SegmentFlags::kShowerTrajectory);
            if (!topo && !traj && sg_length > 10 * units::cm)
                flag_bad_shower_1 = true;
        }
        // stem too long → definitely bad
        if (sg_length > 80 * units::cm)
            flag_bad_shower_1 = true;

        if (flag_fill && ti) {
            ti->br1_1_flag              = !flag_bad_shower_1;
            ti->br1_1_shower_type       = start_type;
            ti->br1_1_vtx_n_segs        = vtx_degree(vtx, graph);
            ti->br1_1_energy            = static_cast<float>(Eshower / units::MeV);
            ti->br1_1_n_segs            = static_cast<float>(shower_segs.size());
            ti->br1_1_flag_sg_topology  = sg->flags_any(SegmentFlags::kShowerTopology);
            ti->br1_1_flag_sg_trajectory= sg->flags_any(SegmentFlags::kShowerTrajectory);
            ti->br1_1_sg_length         = static_cast<float>(sg_length / units::cm);
        }
    }

    // -----------------------------------------------------------------------
    // Sub-check 2: look for a long muon-like track inside the shower.
    // For each shower segment sg1, try to extend it via find_cont_muon_segment_nue.
    // Prototype: NeutrinoID_nue_tagger.h ~3505-3610.
    // -----------------------------------------------------------------------
    {
        double max_length    = 0;
        int    n_connected   = 0;
        int    n_connected1  = 0;
        double max_length_ratio = 0;

        for (SegmentPtr sg1 : shower_segs) {
            double length        = segment_track_length(sg1);
            double direct_length = segment_track_direct_length(sg1);
            bool   topo          = sg1->flags_any(SegmentFlags::kShowerTopology);
            bool   avoid_muon    = sg1->flags_any(SegmentFlags::kAvoidMuonCheck);

            if (avoid_muon) continue;
            // Only process non-topology or nearly-straight segments
            if (topo && direct_length <= 0.94 * length) continue;

            auto [sv1, sv2] = find_vertices(graph, sg1);
            if (!sv1 || !sv2) continue;

            double tmp_length = length;
            int    tmp_nc1    = 0;

            // Try to extend through sv1 (if not main vertex)
            if (sv1 != main_vertex) {
                auto [ext_sg, ext_vtx] = find_cont_muon_segment_nue(graph, sg1, sv1, true);
                if (ext_sg) {
                    bool ext_topo  = ext_sg->flags_any(SegmentFlags::kShowerTopology);
                    bool ext_avoid = ext_sg->flags_any(SegmentFlags::kAvoidMuonCheck);
                    double ext_dl  = segment_track_direct_length(ext_sg);
                    double ext_len = segment_track_length(ext_sg);
                    if (!ext_avoid && (!ext_topo || ext_dl > 0.94 * ext_len)) {
                        tmp_length += ext_len;
                        tmp_nc1    += vtx_degree(ext_vtx, graph) - 1;
                    }
                }
            }
            // Try to extend through sv2 (if not main vertex)
            if (sv2 != main_vertex) {
                auto [ext_sg, ext_vtx] = find_cont_muon_segment_nue(graph, sg1, sv2, true);
                if (ext_sg) {
                    bool ext_topo  = ext_sg->flags_any(SegmentFlags::kShowerTopology);
                    bool ext_avoid = ext_sg->flags_any(SegmentFlags::kAvoidMuonCheck);
                    double ext_dl  = segment_track_direct_length(ext_sg);
                    double ext_len = segment_track_length(ext_sg);
                    if (!ext_avoid && (!ext_topo || ext_dl > 0.94 * ext_len)) {
                        tmp_length += ext_len;
                        tmp_nc1    += vtx_degree(ext_vtx, graph) - 1;
                    }
                }
            }

            // Apply 6cm offset for topology segments or cross-cluster segments
            double length_offset = 0;
            int start_cl = sg->cluster() ? sg->cluster()->get_cluster_id() : -1;
            int sg1_cl   = sg1->cluster() ? sg1->cluster()->get_cluster_id() : -1;
            if (topo || sg1_cl != start_cl) length_offset = 6 * units::cm;

            double eff_length = tmp_length - length_offset;
            if (eff_length > max_length) {
                max_length       = eff_length;
                max_length_ratio = (length > 0) ? direct_length / length : 0;
                n_connected1     = tmp_nc1;
                n_connected      = 0;
                if (sv1 != main_vertex) n_connected += vtx_degree(sv1, graph) - 1;
                if (sv2 != main_vertex) n_connected += vtx_degree(sv2, graph) - 1;
            }
        }

        // Apply length thresholds depending on energy
        auto check_len = [&](int nc, double ml, double t0, double t1, double t2, double t3) {
            if (nc <= 1 && ml > t0) return true;
            if (nc == 2 && ml > t1) return true;
            if (nc == 3 && ml > t2) return true;
            if (ml > t3)            return true;
            return false;
        };

        if (Eshower < 200 * units::MeV) {
            flag_bad_shower_2 = check_len(n_connected, max_length,
                38*units::cm, 42*units::cm, 46*units::cm, 50*units::cm);
        } else if (Eshower < 400 * units::MeV) {
            flag_bad_shower_2 = check_len(n_connected, max_length,
                42*units::cm, 49*units::cm, 52*units::cm, 55*units::cm);
            if (n_connected + n_connected1 > 4 && max_length <= 72 * units::cm)
                flag_bad_shower_2 = false;
        } else if (Eshower < 600 * units::MeV) {
            flag_bad_shower_2 = check_len(n_connected, max_length,
                45*units::cm, 48*units::cm, 54*units::cm, 62*units::cm);
        } else if (Eshower < 800 * units::MeV) {
            flag_bad_shower_2 = check_len(n_connected, max_length,
                51*units::cm, 52*units::cm, 56*units::cm, 62*units::cm);
            if (flag_bad_shower_2) {
                if ((vtx_degree(main_vertex, graph) == 1 && max_length < 68 * units::cm) ||
                    (n_connected >= 6 && max_length < 76 * units::cm))
                    flag_bad_shower_2 = false;
            }
            if (shower->get_num_segments() >= 15 && max_length < 60 * units::cm)
                flag_bad_shower_2 = false;
        } else if (Eshower < 1500 * units::MeV) {
            flag_bad_shower_2 = check_len(n_connected, max_length,
                55*units::cm, 60*units::cm, 65*units::cm, 75*units::cm);
        } else {
            flag_bad_shower_2 = check_len(n_connected, max_length,
                55*units::cm, 65*units::cm, 70*units::cm, 75*units::cm);
        }

        if (Eshower > 1000 * units::MeV && flag_bad_shower_2 && max_length_ratio < 0.95)
            flag_bad_shower_2 = false;

        double total_len = shower->get_total_length();
        if (max_length > 0.75 * total_len && max_length > 35 * units::cm)
            flag_bad_shower_2 = true;

        if (flag_fill && ti) {
            ti->br1_2_flag           = !flag_bad_shower_2;
            ti->br1_2_energy         = static_cast<float>(Eshower / units::MeV);
            ti->br1_2_n_connected    = n_connected;
            ti->br1_2_max_length     = static_cast<float>(max_length / units::cm);
            ti->br1_2_n_connected_1  = n_connected1;
            ti->br1_2_vtx_n_segs     = vtx_degree(main_vertex, graph);
            ti->br1_2_n_shower_segs  = shower->get_num_segments();
            ti->br1_2_max_length_ratio= static_cast<float>(max_length_ratio);
            ti->br1_2_shower_length  = static_cast<float>(total_len / units::cm);
        }
    }

    // -----------------------------------------------------------------------
    // Sub-check 3: look for a long straight track near the far end of the
    // start segment (the "main length" test).
    // Prototype: NeutrinoID_nue_tagger.h ~3620-3840.
    // -----------------------------------------------------------------------
    {
        double max_length   = 0;
        int    n_connected  = 0;
        double main_length  = segment_track_length(sg);

        VertexPtr other_vtx = find_other_vertex(graph, sg, vtx);

        if (main_length > 10 * units::cm && other_vtx) {
            Point  other_pt = vtx_fit_pt(other_vtx);
            Vector dir1 = segment_cal_dir_3vector(sg, other_pt, 15 * units::cm);

            // Examine all other shower segments for track-like continuations
            for (SegmentPtr sg1 : shower_segs) {
                if (sg1 == sg) continue;
                bool topo = sg1->flags_any(SegmentFlags::kShowerTopology);
                bool traj = sg1->flags_any(SegmentFlags::kShowerTrajectory);
                double sg1_len = segment_track_length(sg1);

                if (topo || traj || sg1_len < 10 * units::cm) continue;

                auto [pv1, pv2] = find_vertices(graph, sg1);
                if (!pv1 || !pv2) continue;

                Point pt1 = vtx_fit_pt(pv1);
                Point pt2 = vtx_fit_pt(pv2);

                double dis1 = ray_length(Ray{pt1, other_pt});
                double dis2 = ray_length(Ray{pt2, other_pt});

                double tmp_length1 = 0;
                int    tmp_nc = 0;

                if (dis1 < 5 * units::cm) {
                    Vector dir2  = segment_cal_dir_3vector(sg1, pt1, 15 * units::cm);
                    double angle = dir1.angle(dir2) / M_PI * 180.0;
                    if (angle > 170) {
                        if (sg1_len + dis1 > tmp_length1) {
                            tmp_length1 = sg1_len + dis1;
                            tmp_nc = vtx_degree(pv2, graph) - 1;
                        }
                    }
                } else if (dis2 < 5 * units::cm) {
                    Vector dir2  = segment_cal_dir_3vector(sg1, pt2, 15 * units::cm);
                    double angle = dir1.angle(dir2) / M_PI * 180.0;
                    if (angle > 170) {
                        if (sg1_len + dis2 > tmp_length1) {
                            tmp_length1 = sg1_len + dis2;
                            tmp_nc = vtx_degree(pv1, graph) - 1;
                        }
                    }
                } else {
                    // Neither end is close — check if angle matches and look for a
                    // third segment bridging the gap (prototype logic ~3720-3760).
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
                        // Check remaining shower segments for a bridge
                        for (SegmentPtr sg2 : shower_segs) {
                            if (sg2 == sg || sg2 == sg1) continue;
                            double sg2_len = segment_track_length(sg2);
                            if (sg2_len < 10 * units::cm) continue;

                            auto [pv1_2, pv2_2] = find_vertices(graph, sg2);
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
                                    tmp_nc = (dis2 < dis1) ? vtx_degree(pv1, graph) - 1
                                                           : vtx_degree(pv2, graph) - 1;
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

        // Apply length thresholds (same structure as sub-check 2)
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

        // Only truly bad if: stem has no shower topology/trajectory AND shower has 1 main seg
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

        if (flag_fill && ti) {
            ti->br1_3_flag              = !flag_bad_shower_3;
            ti->br1_3_energy            = static_cast<float>(Eshower / units::MeV);
            ti->br1_3_n_connected_p     = n_connected;
            ti->br1_3_max_length_p      = static_cast<float>(max_length / units::cm);
            ti->br1_3_n_shower_segs     = shower->get_num_segments();
            ti->br1_3_flag_sg_topology  = sg->flags_any(SegmentFlags::kShowerTopology);
            ti->br1_3_flag_sg_trajectory= sg->flags_any(SegmentFlags::kShowerTrajectory);
            ti->br1_3_n_shower_main_segs= shower->get_num_main_segments();
            ti->br1_3_sg_length         = static_cast<float>(segment_track_length(sg) / units::cm);
        }
    }

    bool flag_bad = flag_bad_shower_1 || flag_bad_shower_2 || flag_bad_shower_3;
    if (flag_fill && ti) ti->br1_flag = !flag_bad;
    return flag_bad;
}

// ===========================================================================
// cosmic_tagger
//
// Runs 10 cosmic-rejection checks on the event and fills the cosmict_* BDT
// features in TaggerInfo.  Returns true if any check fires (event is cosmic).
//
// Prototype: WCPPID::NeutrinoID::cosmic_tagger in NeutrinoID_cosmic_tagger.h.
//
// Parameters mirror the state that NeutrinoID carries as member variables:
//   graph                 – the PR graph containing all segments/vertices
//   main_vertex           – the main neutrino vertex
//   showers               – all reconstructed showers (from shower_clustering_with_nv)
//   map_segment_in_shower – maps any shower segment → its parent shower
//   map_vertex_to_shower  – maps a vertex → set of showers that START at it
//   segments_in_long_muon – segments that belong to a long muon shower chain
//   main_cluster          – the primary cluster (contains main_vertex)
//   all_clusters          – all clusters including other_clusters
//   dv                    – detector volumes (provides fiducial-volume test)
//   ti                    – TaggerInfo to fill
// ===========================================================================
bool PatternAlgorithms::cosmic_tagger(
    Graph&              graph,
    VertexPtr           main_vertex,
    IndexedShowerSet&   showers,
    ShowerSegmentMap&   map_segment_in_shower,
    VertexShowerSetMap& map_vertex_to_shower,
    IndexedSegmentSet&  segments_in_long_muon,
    Facade::Cluster*    main_cluster,
    std::vector<Facade::Cluster*>& all_clusters,
    IDetectorVolumes::pointer dv,
    TaggerInfo&         ti)
{
    // Direction reference vectors (beam = z, drift = x, vertical = y)
    const Vector dir_beam    (0, 0, 1);
    const Vector dir_drift   (1, 0, 0);
    const Vector dir_vertical(0, 1, 0);

    // get_last_segment_vertex_long_muon takes std::set<SegmentPtr> (no comparator),
    // but IndexedSegmentSet uses SegmentIndexCmp.  Build a plain set once.
    std::set<SegmentPtr> segments_in_long_muon_plain(
        segments_in_long_muon.begin(), segments_in_long_muon.end());

    bool flag_cosmic_1  = false;
    bool flag_cosmic_2  = false;
    bool flag_cosmic_3  = false;
    bool flag_cosmic_4  = false;
    bool flag_cosmic_5  = false;
    bool flag_cosmic_6  = false;
    bool flag_cosmic_7  = false;
    bool flag_cosmic_8  = false;
    bool flag_cosmic_9  = false;
    bool flag_cosmic_10 = false;
    bool flag_cosmic_10_save = false;

    // Fiducial volume utility (obtained from the grouping attached to main_cluster)
    FiducialUtilsPtr fiducial_utils;
    if (main_cluster && main_cluster->grouping())
        fiducial_utils = main_cluster->grouping()->get_fiducialutils();

    auto inside_fv = [&](const Point& p) -> bool {
        if (!fiducial_utils) return true; // no FV check available → conservative
        return fiducial_utils->inside_fiducial_volume(p);
    };

    // Helper: iterate all segments connected to a vertex in the graph
    auto segs_at_vtx = [&](VertexPtr vtx, auto callback) {
        if (!vtx || !vtx->descriptor_valid()) return;
        auto vd = vtx->get_descriptor();
        for (auto [eit, eit_end] = boost::out_edges(vd, graph); eit != eit_end; ++eit) {
            SegmentPtr seg = graph[*eit].segment;
            if (seg) callback(seg);
        }
    };

    // Main vertex position
    Point mv_pt = vtx_fit_pt(main_vertex);

    // -----------------------------------------------------------------------
    // Determine whether main vertex is inside FV.
    // Prototype: if distance between fit point and wcp point > 5cm, use wcp point.
    // -----------------------------------------------------------------------
    {
        Point mv_wcp_pt = main_vertex->wcpt().point;
        double tmp_dis  = ray_length(Ray{mv_pt, mv_wcp_pt});
        Point test_p    = (tmp_dis > 5 * units::cm) ? mv_wcp_pt : mv_pt;

        // Use tighter tolerance (STM-like -1.5cm) for the vertex outside-FV check.
        // The prototype calls fid->inside_fiducial_volume(test_p, offset_x, &stm_tol_vec)
        // with all tolerances set to -1.5 cm, meaning the fiducial volume is shrunk by 1.5cm.
        // The toolkit FiducialUtils::inside_fiducial_volume already accounts for offset;
        // the stm_tol_vec shrinkage is not yet parameterizable, so we use the default here.
        // TODO: if stm tolerance is needed, extend FiducialUtils API.
        if (fiducial_utils && !fiducial_utils->inside_fiducial_volume(test_p))
            flag_cosmic_1 = true;
    }

    // -----------------------------------------------------------------------
    // Section: single-muon cosmic test (flags 2-5)
    //
    // Find the longest muon-like segment (or long-muon shower) connected to
    // main_vertex and check its direction / endpoint topology.
    // Prototype: NeutrinoID_cosmic_tagger.h lines 43-265.
    // -----------------------------------------------------------------------
    {
        double max_length = 0;
        SegmentPtr muon     = nullptr;
        ShowerPtr  long_muon= nullptr;
        int valid_tracks    = 0;

        // Find the dominant muon/muon-shower candidate
        segs_at_vtx(main_vertex, [&](SegmentPtr sg) {
            double length      = segment_track_length(sg);
            double med_dqdx    = segment_median_dQ_dx(sg);
            double dqdx_cut    = 0.8866 + 0.9533 * std::pow(18 * units::cm / length, 0.4234);
            int    pdg         = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;
            bool   is_shower   = seg_is_shower(sg);

            bool muon_like = (pdg == 13) ||
                             (!is_shower && med_dqdx < dqdx_cut * 1.05 * 43e3 / units::cm && pdg != 211);

            if (muon_like) {
                if (segments_in_long_muon.count(sg)) {
                    // This segment belongs to a long muon shower chain
                    ShowerPtr lm = map_segment_in_shower.count(sg) ? map_segment_in_shower.at(sg) : nullptr;
                    if (lm) {
                        double lm_len = lm->get_total_track_length();
                        if (lm_len > max_length) {
                            max_length = lm_len;
                            long_muon  = lm;
                            muon       = nullptr;
                        }
                    }
                } else {
                    if (length > max_length) {
                        max_length = length;
                        muon       = sg;
                        long_muon  = nullptr;
                    }
                }
            }
        });

        // Count valid (non-muon) tracks and showers at main vertex
        segs_at_vtx(main_vertex, [&](SegmentPtr sg) {
            if (sg == muon) return;
            if (long_muon && sg == long_muon->start_segment()) return;

            double length = segment_track_length(sg);
            if (length > 2.5 * units::cm ||
                (length > 0.9 * units::cm && sg->has_particle_info() && sg->particle_info()->pdg() == 2212) ||
                !sg->dir_weak()) {
                if (length < 15 * units::cm && segment_median_dQ_dx(sg) / (43e3 / units::cm) < 0.75 && sg->dir_weak())
                    return;
                valid_tracks++;
            }
        });

        // Count high-energy showers at main vertex
        int connected_showers = 0;
        if (map_vertex_to_shower.count(main_vertex)) {
            for (ShowerPtr shower : map_vertex_to_shower.at(main_vertex)) {
                double Eshower = shower->get_kine_best();
                auto [sv, stype] = shower->get_start_vertex_and_type();
                if (stype > 2) continue;
                if (shower == long_muon) continue;
                if (Eshower > 150 * units::MeV && !bad_reconstruction(graph, main_vertex, shower))
                    valid_tracks++;
                if (stype == 1 && Eshower > 70 * units::MeV &&
                    !bad_reconstruction(graph, main_vertex, shower) &&
                    shower->start_segment() &&
                    shower->start_segment()->has_particle_info() &&
                    shower->start_segment()->particle_info()->pdg() != 13)
                    connected_showers++;
            }
        }

        // Flag 2: single segment muon going in wrong direction or exiting FV
        if (muon) {
            VertexPtr other_vtx = find_other_vertex(graph, muon, main_vertex);
            Point test_p1 = vtx_fit_pt(other_vtx);
            Vector dir    = segment_cal_dir_3vector(muon, mv_pt, 15 * units::cm);
            bool flag_inside  = inside_fv(test_p1);

            // dQ/dx at near and far ends (first/last 10 fit points)
            double dQ_dx_front = 0, dQ_dx_end = 0;
            const auto& wcps = muon->wcpts();
            bool start_is_main = !wcps.empty() &&
                ray_length(Ray{wcps.front().point, main_vertex->wcpt().point}) <
                ray_length(Ray{wcps.back().point,  main_vertex->wcpt().point});

            if (start_is_main) {
                dQ_dx_front = segment_median_dQ_dx(muon, 0, 10);
                int n = (int)muon->fits().size();
                dQ_dx_end = segment_median_dQ_dx(muon, n - 10, n);
            } else {
                int n = (int)muon->fits().size();
                dQ_dx_end   = segment_median_dQ_dx(muon, 0, 10);
                dQ_dx_front = segment_median_dQ_dx(muon, n - 10, n);
            }

            int    n_muon_tracks       = calculate_num_daughter_tracks(graph, main_vertex, muon, false, 3 * units::cm).first;
            double total_shower_length = calculate_num_daughter_showers(graph, main_vertex, muon).second;

            ti.cosmict_2_filled           = 1;
            ti.cosmict_2_particle_type    = muon->has_particle_info() ? muon->particle_info()->pdg() : 0;
            ti.cosmict_2_n_muon_tracks    = n_muon_tracks;
            ti.cosmict_2_total_shower_length = static_cast<float>(total_shower_length / units::cm);
            ti.cosmict_2_flag_inside      = flag_inside;
            ti.cosmict_2_angle_beam       = static_cast<float>(dir.angle(dir_beam) / M_PI * 180.0);
            ti.cosmict_2_flag_dir_weak    = muon->dir_weak();
            ti.cosmict_2_dQ_dx_end        = static_cast<float>(dQ_dx_end / (43e3 / units::cm));
            ti.cosmict_2_dQ_dx_front      = static_cast<float>(dQ_dx_front / (43e3 / units::cm));
            ti.cosmict_2_theta            = static_cast<float>(vec_theta(dir) / M_PI * 180.0);
            ti.cosmict_2_phi              = static_cast<float>(vec_phi(dir) / M_PI * 180.0);
            ti.cosmict_2_valid_tracks     = valid_tracks;

            int pdg = muon->has_particle_info() ? muon->particle_info()->pdg() : 0;
            if (pdg == 13 && n_muon_tracks <= 2 && total_shower_length < 40 * units::cm &&
                (((!flag_inside) && dir.angle(dir_beam) / M_PI * 180.0 > 40) ||
                 (flag_inside && ((muon->dir_weak() && !(dQ_dx_end > 1.4 * 43e3 / units::cm && dQ_dx_end > 1.2 * dQ_dx_front)) ||
                                  dir.angle(dir_beam) / M_PI * 180.0 > 60))) &&
                (vec_theta(dir) / M_PI * 180.0 >= 100.0 || std::fabs(std::fabs(vec_phi(dir) / M_PI * 180.0) - 90) <= 50) &&
                valid_tracks == 0)
                flag_cosmic_2 = true;

        } else if (long_muon) {
            // Flag 3: long muon shower chain going in wrong direction
            SegmentPtr start_sg = long_muon->start_segment();
            auto [last_sg, other_vtx] = long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon_plain);

            Point test_p1 = vtx_fit_pt(other_vtx);
            Vector dir    = shower_cal_dir_3vector(*long_muon, mv_pt, 30 * units::cm);
            bool flag_inside = inside_fv(test_p1);

            double dQ_dx_front = 0, dQ_dx_end = 0;
            if (start_sg) {
                const auto& wcps = start_sg->wcpts();
                bool start_is_main = !wcps.empty() &&
                    ray_length(Ray{wcps.front().point, main_vertex->wcpt().point}) <
                    ray_length(Ray{wcps.back().point,  main_vertex->wcpt().point});
                if (start_is_main) {
                    dQ_dx_front = segment_median_dQ_dx(start_sg, 0, 10);
                } else {
                    int n = (int)start_sg->fits().size();
                    dQ_dx_front = segment_median_dQ_dx(start_sg, n - 10, n);
                }
            }
            if (last_sg) {
                const auto& wcps = last_sg->wcpts();
                bool end_is_other = !wcps.empty() &&
                    ray_length(Ray{wcps.front().point, other_vtx->wcpt().point}) <
                    ray_length(Ray{wcps.back().point,  other_vtx->wcpt().point});
                if (end_is_other) {
                    dQ_dx_end = segment_median_dQ_dx(last_sg, 0, 10);
                } else {
                    int n = (int)last_sg->fits().size();
                    dQ_dx_end = segment_median_dQ_dx(last_sg, n - 10, n);
                }
            }

            ti.cosmict_3_filled        = 1;
            ti.cosmict_3_flag_inside   = flag_inside;
            ti.cosmict_3_angle_beam    = static_cast<float>(dir.angle(dir_beam) / M_PI * 180.0);
            ti.cosmict_3_flag_dir_weak = last_sg ? last_sg->dir_weak() : false;
            ti.cosmict_3_dQ_dx_front   = static_cast<float>(dQ_dx_front / (43e3 / units::cm));
            ti.cosmict_3_dQ_dx_end     = static_cast<float>(dQ_dx_end   / (43e3 / units::cm));
            ti.cosmict_3_theta         = static_cast<float>(vec_theta(dir) / M_PI * 180.0);
            ti.cosmict_3_phi           = static_cast<float>(vec_phi(dir)   / M_PI * 180.0);
            ti.cosmict_3_valid_tracks  = valid_tracks;

            bool dir_weak_end = last_sg && last_sg->dir_weak();
            if ((((!flag_inside) && dir.angle(dir_beam) / M_PI * 180.0 > 40) ||
                 (flag_inside && ((dir_weak_end && !(dQ_dx_end > 1.4 * 43e3 / units::cm && dQ_dx_end > 1.2 * dQ_dx_front)) ||
                                  dir.angle(dir_beam) / M_PI * 180.0 > 60))) &&
                (vec_theta(dir) / M_PI * 180.0 >= 100.0 || std::fabs(std::fabs(vec_phi(dir) / M_PI * 180.0) - 90) <= 50) &&
                valid_tracks == 0)
                flag_cosmic_3 = true;
        }

        // Flags 4/5: muon exiting detector upstream (angle > 100° from beam)
        if (muon) {
            VertexPtr other_vtx = find_other_vertex(graph, muon, main_vertex);
            Point test_p1 = vtx_fit_pt(other_vtx);
            Vector dir    = segment_cal_dir_3vector(muon, mv_pt, 15 * units::cm);
            bool flag_inside = inside_fv(test_p1);

            ti.cosmict_4_filled           = 1;
            ti.cosmict_4_flag_inside      = flag_inside;
            ti.cosmict_4_angle_beam       = static_cast<float>(dir.angle(dir_beam) / M_PI * 180.0);
            ti.cosmict_4_connected_showers = connected_showers;

            if (!flag_inside && dir.angle(dir_beam) / M_PI * 180.0 > 100 && connected_showers == 0)
                flag_cosmic_4 = true;

        } else if (long_muon) {
            auto [last_sg, other_vtx] = long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon_plain);
            Point test_p1 = vtx_fit_pt(other_vtx);
            Vector dir    = shower_cal_dir_3vector(*long_muon, mv_pt, 30 * units::cm);
            bool flag_inside = inside_fv(test_p1);

            ti.cosmict_5_filled           = 1;
            ti.cosmict_5_flag_inside      = flag_inside;
            ti.cosmict_5_angle_beam       = static_cast<float>(dir.angle(dir_beam) / M_PI * 180.0);
            ti.cosmict_5_connected_showers = connected_showers;

            if (!flag_inside && dir.angle(dir_beam) / M_PI * 180.0 > 100 && connected_showers == 0)
                flag_cosmic_5 = true;
        }
    }

    // -----------------------------------------------------------------------
    // Section: stopped-muon-with-Michel-electron test (flags 6-8).
    // Prototype: NeutrinoID_cosmic_tagger.h lines 268-588.
    // -----------------------------------------------------------------------
    {
        ShowerPtr  michel_ele = nullptr;
        double     michel_energy = 0;
        SegmentPtr muon     = nullptr;
        SegmentPtr muon_2nd = nullptr;
        ShowerPtr  long_muon= nullptr;
        double     max_length = 0;
        int        valid_tracks = 0;

        // Find muon + michel candidates at main vertex
        segs_at_vtx(main_vertex, [&](SegmentPtr sg) {
            bool is_shower = seg_is_shower(sg);
            double length  = segment_track_length(sg);
            int pdg        = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;

            if (is_shower) {
                // Check if this segment starts a shower that could be michel
                for (const ShowerPtr& shower : showers) {
                    if (shower->start_segment() != sg) continue;
                    double E = shower->get_kine_best();
                    if (E > michel_energy) {
                        michel_energy = E;
                        michel_ele    = shower;
                    }
                    break;
                }
            } else {
                if (pdg == 13) {
                    if (segments_in_long_muon.count(sg)) {
                        ShowerPtr lm = map_segment_in_shower.count(sg) ? map_segment_in_shower.at(sg) : nullptr;
                        if (lm) {
                            double lm_len = lm->get_total_track_length();
                            if (lm_len > max_length) {
                                max_length = lm_len;
                                long_muon  = lm;
                                muon       = nullptr;
                            }
                        }
                    } else {
                        if (length > max_length) {
                            max_length = length;
                            muon       = sg;
                            long_muon  = nullptr;
                        }
                    }
                }
            }
        });

        // Find the second-longest non-muon track
        double max_length2 = 0;
        segs_at_vtx(main_vertex, [&](SegmentPtr sg) {
            if (sg == muon) return;
            if (long_muon && sg == long_muon->start_segment()) return;
            if (michel_ele && sg == michel_ele->start_segment()) return;
            if (seg_is_shower(sg)) return;
            int pdg = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;
            if (pdg != 2212) {
                double length = segment_track_length(sg);
                if (length > max_length2) {
                    max_length2 = length;
                    muon_2nd = sg;
                }
            }
        });

        // Count valid remaining tracks
        segs_at_vtx(main_vertex, [&](SegmentPtr sg) {
            if (sg == muon) return;
            if (long_muon && sg == long_muon->start_segment()) return;
            if (michel_ele && sg == michel_ele->start_segment()) return;
            if (seg_is_shower(sg)) return;
            if (sg == muon_2nd) return;
            double length = segment_track_length(sg);
            if (length > 2.5 * units::cm ||
                (length > 0.9 * units::cm && sg->has_particle_info() && sg->particle_info()->pdg() == 2212) ||
                !sg->dir_weak()) {
                if (length < 15 * units::cm && segment_median_dQ_dx(sg) / (43e3 / units::cm) < 0.75 && sg->dir_weak())
                    return;
                valid_tracks++;
            }
        });

        // Count high-energy non-michel showers at main vertex
        if (map_vertex_to_shower.count(main_vertex)) {
            for (ShowerPtr shower : map_vertex_to_shower.at(main_vertex)) {
                if (shower == michel_ele || shower == long_muon) continue;
                double E = shower->get_kine_best();
                auto [sv, stype] = shower->get_start_vertex_and_type();
                if (stype > 2) continue;
                if (E > 150 * units::MeV && !bad_reconstruction(graph, main_vertex, shower))
                    valid_tracks++;
            }
        }

        // Refine valid_tracks using angular relationship with muon and muon_2nd
        if ((muon || long_muon) && michel_ele && muon_2nd) {
            double length2nd = segment_track_length(muon_2nd);
            if (length2nd > 2.5 * units::cm ||
                (length2nd > 0.9 * units::cm && muon_2nd->has_particle_info() && muon_2nd->particle_info()->pdg() == 2212) ||
                !muon_2nd->dir_weak()) {
                if (!(length2nd < 15 * units::cm && segment_median_dQ_dx(muon_2nd) / (43e3 / units::cm) < 0.75 && muon_2nd->dir_weak()))
                    valid_tracks++;
            }

            Vector dir1 = muon ? segment_cal_dir_3vector(muon, mv_pt, 15 * units::cm)
                               : shower_cal_dir_3vector(*long_muon, mv_pt, 30 * units::cm);
            Vector dir2 = shower_cal_dir_3vector(*michel_ele, mv_pt, 15 * units::cm);
            Vector dir3 = segment_cal_dir_3vector(muon_2nd, mv_pt, 15 * units::cm);

            double Emi = michel_energy;
            if (Emi < 25 * units::cm && dir1.angle(dir2) / M_PI * 180.0 > 170)
                valid_tracks--;
            else if (length2nd < 5 * units::cm && dir1.angle(dir3) / M_PI * 180.0 > 170)
                valid_tracks--;

            if (dir1.angle(dir3) / M_PI * 180.0 > 175 && valid_tracks <= 2) {
                if (!(length2nd < 5 * units::cm && dir1.angle(dir3) / M_PI * 180.0 > 170))
                    valid_tracks--;
                segs_at_vtx(main_vertex, [&](SegmentPtr sg1) {
                    if (sg1 == muon || sg1 == muon_2nd) return;
                    if (michel_ele && sg1 == michel_ele->start_segment()) return;
                    if (long_muon && sg1 == long_muon->start_segment()) return;
                    if (sg1->dir_weak() && segment_track_length(sg1) < 5 * units::cm)
                        valid_tracks--;
                });
            }

            if (valid_tracks < 0) valid_tracks = 0;
        }

        // Flags 6 and 7
        if ((muon || long_muon) &&
            ((michel_ele && valid_tracks == 0) || (valid_tracks == 0 && muon_2nd && !michel_ele))) {

            double dQ_dx_front = 0, dQ_dx_end = 0;
            Vector dir;
            bool   flag_inside   = false;
            bool   flag_weak_dir = false;
            double muon_length   = 0;
            int    n_muon_tracks = 0;
            double total_shower_length = 0;

            if (muon) {
                VertexPtr other_vtx = find_other_vertex(graph, muon, main_vertex);
                flag_inside   = inside_fv(vtx_fit_pt(other_vtx));
                dir           = segment_cal_dir_3vector(muon, mv_pt, 15 * units::cm);
                const auto& wcps = muon->wcpts();
                bool start_is_main = !wcps.empty() &&
                    ray_length(Ray{wcps.front().point, main_vertex->wcpt().point}) <
                    ray_length(Ray{wcps.back().point,  main_vertex->wcpt().point});
                if (start_is_main) {
                    dQ_dx_front = segment_median_dQ_dx(muon, 0, 10);
                    int n = (int)muon->fits().size();
                    dQ_dx_end = segment_median_dQ_dx(muon, n - 10, n);
                } else {
                    int n = (int)muon->fits().size();
                    dQ_dx_end   = segment_median_dQ_dx(muon, 0, 10);
                    dQ_dx_front = segment_median_dQ_dx(muon, n - 10, n);
                }
                flag_weak_dir    = muon->dir_weak();
                muon_length      = segment_track_length(muon);
                n_muon_tracks    = calculate_num_daughter_tracks(graph, main_vertex, muon, false, 3 * units::cm).first;
                total_shower_length = calculate_num_daughter_showers(graph, main_vertex, muon).second;

            } else if (long_muon) {
                SegmentPtr start_sg = long_muon->start_segment();
                auto [last_sg, other_vtx] = long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon_plain);
                flag_inside = inside_fv(vtx_fit_pt(other_vtx));
                dir = shower_cal_dir_3vector(*long_muon, mv_pt, 30 * units::cm);
                if (start_sg) {
                    const auto& wcps = start_sg->wcpts();
                    bool start_is_main = !wcps.empty() &&
                        ray_length(Ray{wcps.front().point, main_vertex->wcpt().point}) <
                        ray_length(Ray{wcps.back().point,  main_vertex->wcpt().point});
                    dQ_dx_front = start_is_main
                        ? segment_median_dQ_dx(start_sg, 0, 10)
                        : segment_median_dQ_dx(start_sg, (int)start_sg->fits().size() - 10, (int)start_sg->fits().size());
                }
                if (last_sg) {
                    const auto& wcps = last_sg->wcpts();
                    bool end_is_other = !wcps.empty() &&
                        ray_length(Ray{wcps.front().point, other_vtx->wcpt().point}) <
                        ray_length(Ray{wcps.back().point,  other_vtx->wcpt().point});
                    dQ_dx_end = end_is_other
                        ? segment_median_dQ_dx(last_sg, 0, 10)
                        : segment_median_dQ_dx(last_sg, (int)last_sg->fits().size() - 10, (int)last_sg->fits().size());
                    flag_weak_dir = last_sg->dir_weak();
                }
                muon_length   = long_muon->get_total_track_length();
                n_muon_tracks = 1;
                total_shower_length = 0;
            }

            // Flag 6: muon_2nd exits FV in opposite direction from muon
            bool flag_sec = false;
            if (michel_ele) {
                double Emi = michel_ele->get_kine_best();
                if (Emi < 70 * units::MeV) {
                    flag_sec = true;
                } else {
                    Vector dir2 = shower_cal_dir_3vector(*michel_ele, mv_pt, 30 * units::cm);
                    double tmp_angle = dir.angle(dir2) / M_PI * 180.0;
                    if ((muon_length > 120 * units::cm && tmp_angle > 155) ||
                        (muon_length > 20 * units::cm  && tmp_angle > 175) ||
                        (muon_length > 60 * units::cm  && tmp_angle > 167.5))
                        flag_sec = true;
                }
            } else if (muon_2nd) {
                if (muon_2nd->dir_weak() && segment_track_length(muon_2nd) < 8 * units::cm)
                    flag_sec = true;

                Vector dir2 = segment_cal_dir_3vector(muon_2nd, mv_pt, 15 * units::cm);
                VertexPtr other_vtx2 = find_other_vertex(graph, muon_2nd, main_vertex);

                ti.cosmict_6_filled       = 1;
                ti.cosmict_6_flag_dir_weak= muon_2nd->dir_weak();
                ti.cosmict_6_flag_inside  = inside_fv(vtx_fit_pt(other_vtx2));
                ti.cosmict_6_angle        = static_cast<float>(dir.angle(dir2) / M_PI * 180.0);

                if (muon_2nd->dir_weak() && !inside_fv(vtx_fit_pt(other_vtx2))) {
                    if (dir.angle(dir2) / M_PI * 180.0 > 170)
                        flag_cosmic_6 = true;
                }
            }

            // Flag 7: stopped muon + secondary (michel or nearly back-to-back track)
            ti.cosmict_7_filled          = 1;
            ti.cosmict_7_flag_sec        = flag_sec;
            ti.cosmict_7_n_muon_tracks   = n_muon_tracks;
            ti.cosmict_7_total_shower_length = static_cast<float>(total_shower_length / units::cm);
            ti.cosmict_7_flag_inside     = flag_inside;
            ti.cosmict_7_angle_beam      = static_cast<float>(dir.angle(dir_beam) / M_PI * 180.0);
            ti.cosmict_7_flag_dir_weak   = flag_weak_dir;
            ti.cosmict_7_dQ_dx_end       = static_cast<float>(dQ_dx_end   / (43e3 / units::cm));
            ti.cosmict_7_dQ_dx_front     = static_cast<float>(dQ_dx_front / (43e3 / units::cm));
            ti.cosmict_7_theta           = static_cast<float>(vec_theta(dir) / M_PI * 180.0);
            ti.cosmict_7_phi             = static_cast<float>(vec_phi(dir)   / M_PI * 180.0);

            if (flag_sec && n_muon_tracks <= 2 && total_shower_length < 40 * units::cm) {
                if ((((!flag_inside) && dir.angle(dir_beam) / M_PI * 180.0 > 40) ||
                     (flag_inside && ((flag_weak_dir && !(dQ_dx_end > 1.4 * 43e3 / units::cm && dQ_dx_end > 1.2 * dQ_dx_front)) ||
                                      dir.angle(dir_beam) / M_PI * 180.0 > 60))) &&
                    (vec_theta(dir) / M_PI * 180.0 >= 100.0 || std::fabs(std::fabs(vec_phi(dir) / M_PI * 180.0) - 90) <= 50))
                    flag_cosmic_7 = true;
            }

        } else if ((muon || long_muon) && valid_tracks == 1) {
            // Flag 8: muon + one other track, the other track exits FV backwards
            Vector dir;
            double muon_length;
            if (muon) {
                dir          = segment_cal_dir_3vector(muon, mv_pt, 15 * units::cm);
                muon_length  = segment_track_length(muon);
            } else {
                dir          = shower_cal_dir_3vector(*long_muon, mv_pt, 30 * units::cm);
                muon_length  = long_muon->get_total_track_length();
            }

            double acc_length = 0;
            bool   flag_out   = false;
            segs_at_vtx(main_vertex, [&](SegmentPtr sg1) {
                if (sg1 == muon) return;
                if (long_muon && sg1 == long_muon->start_segment()) return;
                Vector dir2 = segment_cal_dir_3vector(sg1, mv_pt, 15 * units::cm);
                if (dir.angle(dir2) / M_PI * 180.0 > 165) {
                    VertexPtr other_vtx = find_other_vertex(graph, sg1, main_vertex);
                    if (!inside_fv(vtx_fit_pt(other_vtx))) flag_out = true;
                } else {
                    acc_length += segment_track_length(sg1);
                }
            });

            ti.cosmict_8_filled      = 1;
            ti.cosmict_8_flag_out    = flag_out;
            ti.cosmict_8_muon_length = static_cast<float>(muon_length / units::cm);
            ti.cosmict_8_acc_length  = static_cast<float>(acc_length  / units::cm);

            if (flag_out && muon_length > 100 * units::cm && acc_length < 12 * units::cm)
                flag_cosmic_8 = true;
        }
    }

    // Cluster id of main_cluster — used in flag 9 and flag 10 blocks.
    const int main_cl_id = main_cluster ? main_cluster->get_cluster_id() : -1;

    // -----------------------------------------------------------------------
    // Section: global cluster-direction PCA analysis (flag 9).
    // Checks if the event looks like a collection of cosmic tracks by examining
    // the PCA direction of each cluster.
    // Prototype: NeutrinoID_cosmic_tagger.h lines 594-796.
    // -----------------------------------------------------------------------
    {
        // Collect per-cluster point clouds and lengths
        std::map<int, std::vector<Point>> map_cl_pts;
        std::map<int, Point>  map_cl_high_pt;
        std::map<int, int>    map_cl_shower_pts;
        std::map<int, double> map_cl_length;
        std::set<int>         big_cluster_ids;

        int    num_small_pieces   = 0;
        double acc_small_length   = 0;

        for (Facade::Cluster* cl : all_clusters) {
            if (!cl) continue;
            int    cl_id = cl->get_cluster_id();
            double len   = cl->get_length();
            map_cl_length[cl_id] = len;
            map_cl_shower_pts[cl_id] = 0;
            if (len > 3 * units::cm) {
                big_cluster_ids.insert(cl_id);
            } else {
                // Check if top part of detector (y > 50cm)
                const auto& pca = cl->get_pca();
                if (pca.center.y() > 50 * units::cm) {
                    num_small_pieces++;
                    acc_small_length += len;
                }
            }
        }

        // Collect segment points for big clusters
        for (auto [eit, eit_end] = boost::edges(graph); eit != eit_end; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg || !sg->cluster()) continue;
            int cl_id = sg->cluster()->get_cluster_id();
            if (!big_cluster_ids.count(cl_id)) continue;

            const auto& fits = sg->fits();
            if ((int)fits.size() <= 2) continue;

            if (!map_cl_high_pt.count(cl_id))
                map_cl_high_pt[cl_id] = fits.front().point;

            for (size_t i = 1; i + 1 < fits.size(); i++) {
                Point p = fits[i].point;
                map_cl_pts[cl_id].push_back(p);
                if (p.y() > map_cl_high_pt[cl_id].y())
                    map_cl_high_pt[cl_id] = p;
            }
            if (seg_is_shower(sg))
                map_cl_shower_pts[cl_id] += static_cast<int>(fits.size());
        }

        // Add vertex positions to cluster point clouds
        for (auto [vit, vit_end] = boost::vertices(graph); vit != vit_end; ++vit) {
            VertexPtr vtx = graph[*vit].vertex;
            if (!vtx || !vtx->cluster()) continue;
            int cl_id = vtx->cluster()->get_cluster_id();
            if (!big_cluster_ids.count(cl_id)) continue;
            Point p = vtx_fit_pt(vtx);
            if (map_cl_high_pt.count(cl_id) && p.y() > map_cl_high_pt[cl_id].y())
                map_cl_high_pt[cl_id] = p;
            map_cl_pts[cl_id].push_back(p);
        }

        // Compute PCA-based angles for each cluster
        int    num_cosmic      = 0;
        double acc_cosmic_len  = 0;
        double acc_total_len   = 0;
        double highest_y       = -100 * units::cm;
        double max_length      = 0;
        bool   flag_main_cluster = true;

        for (auto& [cl_id, pts] : map_cl_pts) {
            if (pts.empty()) continue;

            int    n_shower_pts = map_cl_shower_pts[cl_id];
            double cl_len       = map_cl_length[cl_id];
            double shower_frac  = (double)n_shower_pts / (double)pts.size();

            double angle_cosmic, angle_beam;

            // For the main cluster or large clusters, use centroid-based direction
            bool use_centroid = (cl_id == main_cl_id &&
                                 ((shower_frac > 0.7 && cl_len < 45 * units::cm) ||
                                  (shower_frac < 0.7 && cl_len > 40 * units::cm))) ||
                                cl_len > 60 * units::cm;

            if (use_centroid) {
                Vector sum(0, 0, 0);
                for (const Point& p : pts) {
                    sum += Vector(p.x() - mv_pt.x(), p.y() - mv_pt.y(), p.z() - mv_pt.z());
                }
                angle_beam   = sum.angle(dir_beam) / M_PI * 180.0;
                if (angle_beam > 90) angle_beam = 180 - angle_beam;
                angle_cosmic = 180.0 - sum.angle(dir_vertical) / M_PI * 180.0;
            } else {
                // PCA via calc_PCA_main_axis
                auto pca_axis = calc_PCA_main_axis(pts).second;
                angle_cosmic  = pca_axis.angle(dir_vertical) / M_PI * 180.0;
                angle_beam    = pca_axis.angle(dir_beam)     / M_PI * 180.0;
                if (angle_cosmic > 90) angle_cosmic = 180 - angle_cosmic;
                if (angle_beam   > 90) angle_beam   = 180 - angle_beam;
            }

            // Determine if main cluster is non-cosmic (neutrino-like track)
            if (cl_id == main_cl_id) {
                if (shower_frac < 0.3 && cl_len > 20 * units::cm && angle_cosmic > 40)
                    flag_main_cluster = false;
                else if (cl_len > 80 * units::cm && angle_cosmic > 25)
                    flag_main_cluster = false;
            }

            double highest_y_cl = map_cl_high_pt.count(cl_id) ? map_cl_high_pt[cl_id].y() : -999 * units::cm;
            if (highest_y_cl > highest_y) highest_y = highest_y_cl;

            // Count as cosmic if nearly vertical and not along beam
            if (shower_frac > 0.7) {
                if ((angle_cosmic < 30 && angle_beam > 30) ||
                    (angle_cosmic < 35 && angle_beam > 40)) {
                    acc_cosmic_len += cl_len;
                    num_cosmic++;
                    if (cl_len > max_length) max_length = cl_len;
                }
            } else {
                bool top_main = (cl_id == main_cl_id && highest_y_cl > 100 * units::cm);
                if (angle_cosmic < 20 || (angle_cosmic < 30 && top_main)) {
                    acc_cosmic_len += cl_len;
                    num_cosmic++;
                    if (cl_len > max_length) max_length = cl_len;
                }
            }

            acc_total_len += cl_len;
        }

        // Cosmic candidate if enough of the event is vertically oriented
        bool flagp_cosmic = false;
        if (((num_cosmic > 2  && acc_cosmic_len + acc_small_length > 0.55 * acc_total_len) ||
             (num_cosmic >= 2 && acc_cosmic_len + acc_small_length > 0.70 * acc_total_len) ||
             (num_cosmic == 1 && acc_cosmic_len + acc_small_length > 0.625 * acc_total_len && highest_y > 102 * units::cm)) &&
            mv_pt.y() > 0 && flag_main_cluster && highest_y > 80 * units::cm)
            flagp_cosmic = true;

        if (num_cosmic == 1 && acc_cosmic_len > 100 * units::cm)
            flagp_cosmic = false;

        if (flagp_cosmic) {
            // Cross-check: neutrino-like shower at main vertex reduces false-positive
            int    n_solid_tracks             = 0;
            int    n_direct_showers           = 0;
            double energy_direct_showers      = 0;
            int    n_main_showers             = 0;
            double energy_main_showers        = 0;
            int    n_indirect_showers         = 0;
            double energy_indirect_showers    = 0;

            segs_at_vtx(main_vertex, [&](SegmentPtr sg1) {
                if (seg_is_shower(sg1)) return;
                double length = segment_track_length(sg1);
                if ((sg1->dir_weak() && length > 10 * units::cm) || !sg1->dir_weak())
                    n_solid_tracks++;
            });

            for (const ShowerPtr& shower : showers) {
                if (shower->start_segment() &&
                    shower->start_segment()->has_particle_info() &&
                    shower->start_segment()->particle_info()->pdg() != 11) continue;
                double E = shower->get_kine_best();
                if (E <= 60 * units::MeV) continue;

                auto [sv, stype] = shower->get_start_vertex_and_type();
                if (sv == main_vertex && stype == 1) {
                    n_main_showers++;
                    energy_main_showers += E;
                }
                if (sv && sv->cluster() && main_cluster &&
                    sv->cluster()->get_cluster_id() == main_cl_id && stype == 1) {
                    n_direct_showers++;
                    energy_direct_showers += E;
                }
                if (stype > 1) {
                    n_indirect_showers++;
                    energy_indirect_showers += E;
                }
            }

            ti.cosmic_n_solid_tracks           = n_solid_tracks;
            ti.cosmic_energy_main_showers      = static_cast<float>(energy_main_showers    / units::MeV);
            ti.cosmic_energy_direct_showers    = static_cast<float>(energy_direct_showers  / units::MeV);
            ti.cosmic_energy_indirect_showers  = static_cast<float>(energy_indirect_showers/ units::MeV);
            ti.cosmic_n_direct_showers         = n_direct_showers;
            ti.cosmic_n_indirect_showers       = n_indirect_showers;
            ti.cosmic_n_main_showers           = n_main_showers;
            ti.cosmic_filled                   = 1;

            if (((n_solid_tracks > 0 && energy_main_showers > 80 * units::MeV) ||
                 energy_main_showers > 200 * units::MeV) &&
                (energy_indirect_showers < energy_main_showers * 0.6 ||
                 (n_solid_tracks > 0 && energy_indirect_showers < energy_main_showers))) {
                flagp_cosmic  = false;
                ti.cosmic_flag = true;
            } else {
                ti.cosmic_flag = false;
            }
        }

        if (flagp_cosmic)
            flag_cosmic_9 = true;
    }

    // -----------------------------------------------------------------------
    // Section: front-end vertex check (flag 10).
    // Flag tracks that are near z=0 (front face), not in FV, and pointing
    // along the beam direction.
    // Prototype: NeutrinoID_cosmic_tagger.h lines 799-835.
    // -----------------------------------------------------------------------
    {
        for (auto [vit, vit_end] = boost::vertices(graph); vit != vit_end; ++vit) {
            VertexPtr vtx = graph[*vit].vertex;
            if (!vtx) continue;
            if (!vtx->cluster() || !main_cluster) continue;
            if (vtx->cluster()->get_cluster_id() != main_cluster->get_cluster_id()) continue;

            Point vpt = vtx_fit_pt(vtx);
            if (!inside_fv(vpt) && vpt.z() < 15 * units::cm) {
                segs_at_vtx(vtx, [&](SegmentPtr sg) {
                    flag_cosmic_10 = false;
                    Vector dir     = segment_cal_dir_3vector(sg, vpt, 15 * units::cm);
                    double angle_beam = dir.angle(dir_beam) / M_PI * 180.0;
                    if (angle_beam > 90) angle_beam = 180 - angle_beam;
                    double length = segment_track_length(sg);

                    if (!seg_is_shower(sg) && sg->dir_weak() && angle_beam < 25 && length > 10 * units::cm)
                        flag_cosmic_10 = true;

                    ti.cosmict_10_flag_inside.push_back(inside_fv(vpt));
                    ti.cosmict_10_vtx_z.push_back(static_cast<float>(vpt.z() / units::cm));
                    ti.cosmict_10_flag_shower.push_back(seg_is_shower(sg));
                    ti.cosmict_10_flag_dir_weak.push_back(sg->dir_weak());
                    ti.cosmict_10_angle_beam.push_back(static_cast<float>(angle_beam));
                    ti.cosmict_10_length.push_back(static_cast<float>(length / units::cm));
                    ti.cosmict_flag_10.push_back(flag_cosmic_10);

                    if (flag_cosmic_10) flag_cosmic_10_save = true;
                });
            }
        }
    }

    // Collect all flags into tagger_info
    ti.cosmict_flag_1 = flag_cosmic_1;
    ti.cosmict_flag_2 = flag_cosmic_2;
    ti.cosmict_flag_3 = flag_cosmic_3;
    ti.cosmict_flag_4 = flag_cosmic_4;
    ti.cosmict_flag_5 = flag_cosmic_5;
    ti.cosmict_flag_6 = flag_cosmic_6;
    ti.cosmict_flag_7 = flag_cosmic_7;
    ti.cosmict_flag_8 = flag_cosmic_8;
    ti.cosmict_flag_9 = flag_cosmic_9;
    // ti.cosmict_flag_10 is filled per-entry above in the push_back loop

    bool flag_cosmic = flag_cosmic_1 || flag_cosmic_2 || flag_cosmic_3 ||
                       flag_cosmic_4 || flag_cosmic_5 || flag_cosmic_6 ||
                       flag_cosmic_7 || flag_cosmic_8 || flag_cosmic_9 ||
                       flag_cosmic_10_save;

    ti.cosmict_flag = flag_cosmic;
    return flag_cosmic;
}
