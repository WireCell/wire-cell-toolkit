// NeutrinoTaggerNuE.cxx
//
// Ports from prototype:
//   NeutrinoID_nue_tagger.h  — low_energy_michel, single_shower, angular_cut,
//                               track_overclustering, broken_muon_id, shower_to_wall,
//                               gap_identification, mip_quality, mip_identification,
//                               high_energy_overlapping, low_energy_overlapping,
//                               pi0_identification, single_shower_pio_tagger,
//                               bad_reconstruction_3, bad_reconstruction_2, bad_reconstruction_1
//   NeutrinoID_nue_functions.h — stem_direction, multiple_showers, other_showers,
//                                stem_length, vertex_inside_shower, compare_muon_energy
//
// Namespace/class: WireCell::Clus::PR::PatternAlgorithms
//
// Design: all helpers are static file-local functions that take NuEContext& ctx.
// The public entry point nue_tagger() (declared in NeutrinoPatternBase.h) takes
// individual parameters matching the existing numu_tagger() style, constructs
// NuEContext internally, and delegates to the helpers.
//
// Translation conventions (see neutrino_id_function_map.md for complete map):
//   map_vertex_segments[vtx]              → boost::out_edges(vtx->get_descriptor(), ctx.graph)
//   calculate_num_daughter_tracks(v,sg,f) → calculate_num_daughter_tracks(ctx.graph, v, sg, f, 0)
//   tagger_info.X = v  (flag_fill gate)   → ti.X = v  (always fill)
//   flag_print output                     → dropped
//   main_cluster->Calc_PCA(pts)
//     + get_PCA_axis(0)                   → calc_PCA_main_axis(pts).second
//   fid->inside_fiducial_volume(p,.)      → fiducial_utils->inside_fiducial_volume(p)
//   sg->cal_dir_3vector(pt, dis)          → segment_cal_dir_3vector(sg, pt, dis)
//   shower->cal_dir_3vector(pt, dis)      → shower_cal_dir_3vector(*shower, pt, dis)
//   shower->get_start_vertex().first      → shower->get_start_vertex_and_type().first
//   shower->get_start_segment()           → shower->start_segment()
//   sg->get_point_vec().front()           → sg->fits().front().point
//   sg->get_length()                      → segment_track_length(sg)
//   sg->get_length(n1, n2)               → segment_track_length(sg, 0, n1, n2)
//   sg->get_direct_length(n1, n2)        → segment_track_direct_length(sg, n1, n2)
//   sg->get_flag_avoid_muon_check()       → sg->flags_any(SegmentFlags::kAvoidMuonCheck)
//   sg->get_flag_shower()                 → seg_is_shower(sg)
//   sg->is_dir_weak()                     → sg->dir_weak()
//   sg->get_particle_type()               → sg->particle_info()->pdg()
//   sg->get_medium_dQ_dx()               → segment_median_dQ_dx(sg)
//   shower->get_total_length(cluster_id) → shower->get_total_length(sg->cluster())
//   TPCParams::get_muon_r2ke()->Eval(L)  → ctx.particle_data->get_range_function("muon")
//                                            ->scalar_function(L/cm) * MeV
//   vertex wcpt-index front/back check   → geometric proximity of sg endpoints to vertex
//   neutrino_type flag                    → dropped (only appeared in debug prints)

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
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

// Vertex fit point (with wcpt fallback).  Same as in other NuE tagger files.
static inline Point vtx_fit_pt(VertexPtr v) {
    return v->fit().valid() ? v->fit().point : v->wcpt().point;
}

// True if a segment is shower-like (trajectory, topology, or PDG == 11).
static inline bool seg_is_shower(SegmentPtr seg) {
    return seg->flags_any(SegmentFlags::kShowerTrajectory) ||
           seg->flags_any(SegmentFlags::kShowerTopology)   ||
           (seg->has_particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
}

// Determine which end of seg is nearest to a given point.
// Returns the endpoint (fit point) that is geometrically closest to ref_pt.
// Used to replace prototype's wcpt-index comparison for front/back of segment.
static Point seg_endpoint_near(SegmentPtr seg, const Point& ref_pt) {
    const auto& fits = seg->fits();
    Point front = fits.front().point;
    Point back  = fits.back().point;
    return (ray_length(Ray{ref_pt, front}) <= ray_length(Ray{ref_pt, back}))
           ? front : back;
}

// ---------------------------------------------------------------------------
// NuEContext: file-local bundle of all shared state for nue_tagger helpers.
//
// The public entry point PatternAlgorithms::nue_tagger() (declared in
// NeutrinoPatternBase.h) takes individual parameters and constructs this
// internally.  All helper functions take (NuEContext& ctx, ..., TaggerInfo& ti).
// ---------------------------------------------------------------------------
struct NuEContext {
    PatternAlgorithms& self;                            // for calling member functions
    Graph& graph;
    Facade::Cluster* main_cluster;
    VertexPtr main_vertex;
    int apa{0}, face{0};                            // for point-cloud queries — set by caller
    IndexedShowerSet& showers;
    VertexShowerSetMap& map_vertex_to_shower;
    IndexedShowerSet& pi0_showers;
    ShowerIntMap& map_shower_pio_id;
    std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers;
    std::map<int, std::pair<double,int>>& map_pio_id_mass;
    IDetectorVolumes::pointer dv;
    ParticleDataSet::pointer particle_data;
};

// ===========================================================================
// low_energy_michel
//
// Checks whether the shower looks like a low-energy Michel electron:
// too short, or charge dominated by shower topology rather than MIP dQ/dx.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::low_energy_michel()
// Fills: ti.lem_*
// ===========================================================================
static bool low_energy_michel(NuEContext& ctx, ShowerPtr shower, TaggerInfo& ti) {
    bool flag_bad = false;

    double E_dQdx   = shower->get_kine_dQdx();
    double E_charge = shower->get_kine_charge();

    // Collect shower internal topology so we can count branching vertices.
    // Prototype uses shower->get_map_vtx_segs() (shower-internal map).
    // Toolkit: fill_sets gives all shower vertices + segments, then filter
    // to the main cluster and count shower-internal vertex degree.
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    SegmentPtr       start_sg = shower->start_segment();
    Facade::Cluster* start_cl = start_sg ? start_sg->cluster() : nullptr;

    // n_3seg: number of vertices in the shower's main cluster that connect to
    //         ≥ 3 shower segments (i.e. branching points).
    int n_3seg = 0;
    for (VertexPtr vtx1 : shower_vtxs) {
        if (!vtx1->cluster() || vtx1->cluster() != start_cl) continue;
        if (!vtx1->descriptor_valid()) continue;
        int deg = 0;
        for (auto [eit, eend] = boost::out_edges(vtx1->get_descriptor(), ctx.graph);
             eit != eend; ++eit) {
            if (shower_segs.count(ctx.graph[*eit].segment)) ++deg;
        }
        if (deg >= 3) ++n_3seg;
    }

    double total_length = shower->get_total_length();
    double main_length  = start_cl ? shower->get_total_length(start_cl) : total_length;

    // Short-shower criterion (7003_1226_61350)
    if ((total_length < 25*units::cm && main_length > 0.75 * total_length && n_3seg == 0) ||
        (total_length < 18*units::cm && main_length > 0.75 * total_length && n_3seg > 0))
        flag_bad = true;

    // Low-charge MIP criterion (7004_1291_64560 + 7026_879_43995)
    if (E_charge < 100*units::MeV && E_dQdx < 0.7 * E_charge &&
        shower->get_num_segments() == shower->get_num_main_segments())
        flag_bad = true;

    ti.lem_shower_total_length  = total_length / units::cm;
    ti.lem_shower_main_length   = main_length / units::cm;
    ti.lem_n_3seg               = n_3seg;
    ti.lem_e_charge             = E_charge / units::MeV;
    ti.lem_e_dQdx               = E_dQdx / units::MeV;
    ti.lem_shower_num_segs      = shower->get_num_segments();
    ti.lem_shower_num_main_segs = shower->get_num_main_segments();
    ti.lem_flag                 = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// stem_length
//
// Rejects events where the shower stem segment is too long for a genuine
// electromagnetic shower at the given energy.
//
// Prototype: NeutrinoID_nue_functions.h, WCPPID::NeutrinoID::stem_length()
// Fills: ti.stem_len_*
// ===========================================================================
static bool stem_length(NuEContext& ctx, ShowerPtr shower, double energy, TaggerInfo& ti) {
    bool flag_bad = false;

    SegmentPtr sg     = shower->start_segment();
    VertexPtr  vertex = shower->get_start_vertex_and_type().first;

    // pair_result.first  = number of daughter tracks beyond sg
    // pair_result.second = total track length beyond sg
    auto pair_result = ctx.self.calculate_num_daughter_tracks(ctx.graph, vertex, sg, /*count_shower=*/true, 0);

    double sg_length = segment_track_length(sg);
    if (energy < 500*units::MeV && sg_length > 50*units::cm &&
        !sg->flags_any(SegmentFlags::kAvoidMuonCheck)) {
        flag_bad = true;
        // Exception: many daughters → the stem is a muon with secondary kinks
        if (pair_result.first > 6 && sg_length < 55*units::cm) flag_bad = false;
    }

    ti.stem_len_energy               = energy / units::MeV;
    ti.stem_len_length               = sg_length / units::cm;
    ti.stem_len_flag_avoid_muon_check = sg->flags_any(SegmentFlags::kAvoidMuonCheck);
    ti.stem_len_num_daughters        = pair_result.first;
    ti.stem_len_daughter_length      = pair_result.second / units::cm;
    ti.stem_len_flag                 = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// angular_cut
//
// Rejects events where the track material is predominantly backward relative
// to the shower direction (suggesting the "shower" is actually a hadronic
// interaction product going forward while tracks go backward), or where
// shower vertices lie outside the fiducial volume.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::angular_cut()
// Fills: ti.anc_*
// Note: the prototype parameter `angle` is the shower angle w.r.t. beam
//       (computed by the caller in nue_tagger).  A loop-local variable also
//       named `angle` in the prototype is renamed here to `seg_angle` to
//       avoid shadowing.
// ===========================================================================
static bool angular_cut(NuEContext& ctx, ShowerPtr shower,
                        double energy, double angle, TaggerInfo& ti) {
    bool flag_bad = false;

    VertexPtr  vertex = shower->get_start_vertex_and_type().first;
    SegmentPtr sg     = shower->start_segment();
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));

    Vector dir_beam(0, 0, 1);
    Vector dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 30*units::cm);

    double acc_forward_length  = 0;
    double acc_forward_length1 = 0;
    double acc_backward_length = 0;
    double max_angle  = 0;
    double max_length = 0;

    if (vertex && vertex->descriptor_valid()) {
        auto vd = vertex->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1) continue;
            Vector dir1     = segment_cal_dir_3vector(sg1, vertex_point, 15*units::cm);
            double seg_angle = dir1.angle(dir_beam) / M_PI * 180.0;

            // Accumulate track lengths in each hemisphere relative to beam
            auto pair_result = ctx.self.calculate_num_daughter_tracks(ctx.graph, vertex, sg1,
                                                                     /*count_shower=*/true, 0);
            if (seg_angle > 90) {
                acc_backward_length += pair_result.second;
            } else {
                acc_forward_length += pair_result.second;
                if (seg_angle < 85) acc_forward_length1 += pair_result.second;
            }

            // Track the segment most anti-aligned with the shower direction
            double seg_angle1 = dir1.angle(dir_shower) / M_PI * 180.0;
            if (seg_angle1 > max_angle) {
                max_angle  = seg_angle1;
                max_length = pair_result.second;
            }
        }
    }

    double shower_main_length  = sg->cluster() ? shower->get_total_length(sg->cluster()) : 0;
    double shower_total_length = shower->get_total_length();

    // Primary angular cut: backward-dominated track topology.
    // Inner sub-conditions use explicit parentheses around each && clause
    // to suppress -Wparentheses while preserving prototype logic exactly.
    bool cut_1 = (energy < 650*units::MeV && angle > 160);

    bool cut_2 = (energy < 650*units::MeV && angle > 135 &&
                  max_angle > 170 && max_length > 12*units::cm);

    bool cut_3 = (energy < 650*units::MeV && angle > 135 &&
                  ((acc_forward_length < 0.8 * acc_backward_length && acc_forward_length < 15*units::cm) ||
                   (acc_forward_length < 0.6 * acc_backward_length && acc_forward_length >= 15*units::cm &&
                       acc_backward_length >= 80*units::cm) ||
                   (acc_forward_length < 0.4 * acc_backward_length && acc_forward_length >= 15*units::cm)));

    bool cut_4 = (energy < 650*units::MeV &&
                  (acc_forward_length == 0 || acc_forward_length < 0.03 * acc_backward_length ||
                   acc_forward_length1 == 0 || acc_forward_length1 < 0.03 * acc_backward_length) &&
                  acc_backward_length > 0 &&
                  ((acc_backward_length - shower_main_length > acc_forward_length && angle > 90) ||
                   angle <= 90));

    bool cut_5 = (energy >= 650*units::MeV && angle > 90);

    if (cut_1 || cut_2 || cut_3 || cut_4 || cut_5) {
        flag_bad = true;
    }

    // Secondary cut: shower vertices outside fiducial volume
    bool flag_main_outside = false;
    FiducialUtilsPtr fiducial_utils;
    if (ctx.main_cluster && ctx.main_cluster->grouping())
        fiducial_utils = ctx.main_cluster->grouping()->get_fiducialutils();

    if (fiducial_utils) {
        IndexedSegmentSet fv_segs;
        IndexedVertexSet  fv_vtxs;
        shower->fill_sets(fv_vtxs, fv_segs, /*flag_exclude_start_segment=*/false);
        VertexPtr        start_vtx = shower->get_start_vertex_and_type().first;
        Facade::Cluster* start_cl  = sg->cluster();
        for (VertexPtr vtx1 : fv_vtxs) {
            if (vtx1 == start_vtx) continue;
            if (!vtx1->cluster() || vtx1->cluster() != start_cl) continue;
            if (!fiducial_utils->inside_fiducial_volume(vtx_fit_pt(vtx1)))
                flag_main_outside = true;
        }
    }

    if ((angle > 90 || energy < 300*units::MeV ||
         (angle > 60 && energy < 800*units::MeV)) && flag_main_outside) {
        flag_bad = true;
        // Exception: shower mostly outside and very long → may be OK at shallow angle
        if (shower_main_length < 0.5 * shower_total_length &&
            shower_total_length > 80*units::cm && angle < 90)
            flag_bad = false;
    }

    ti.anc_energy              = energy / units::MeV;
    ti.anc_angle               = angle;
    ti.anc_max_angle           = max_angle;
    ti.anc_max_length          = max_length / units::cm;
    ti.anc_acc_forward_length  = acc_forward_length / units::cm;
    ti.anc_acc_backward_length = acc_backward_length / units::cm;
    ti.anc_acc_forward_length1 = acc_forward_length1 / units::cm;
    ti.anc_shower_main_length  = shower_main_length / units::cm;
    ti.anc_shower_total_length = shower_total_length / units::cm;
    ti.anc_flag_main_outside   = flag_main_outside;
    ti.anc_flag                = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// compare_muon_energy
//
// Rejects events where the reconstructed muon energy (from an external
// muon-length estimate) is comparable to or larger than the shower energy,
// suggesting the candidate shower is the muon rather than an electron.
//
// Prototype: NeutrinoID_nue_functions.h, WCPPID::NeutrinoID::compare_muon_energy()
// Fills: ti.cme_*
// Note: the prototype's `neutrino_type` flag was used only in a debug print
//       and is dropped here.
// ===========================================================================
static bool compare_muon_energy(NuEContext& ctx, ShowerPtr shower,
                                double energy, double muon_length, TaggerInfo& ti) {
    bool flag_bad = false;

    Vector dir_drift(1, 0, 0);
    VertexPtr  vertex = shower->get_start_vertex_and_type().first;
    SegmentPtr sg     = shower->start_segment();
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));

    // Shower direction
    Vector dir_shower;
    if (segment_track_length(sg) > 12*units::cm) {
        dir_shower = segment_cal_dir_3vector(sg, vertex_point, 15*units::cm);
    } else {
        dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
    }
    if (std::fabs(dir_shower.angle(dir_drift) / M_PI * 180.0 - 90.0) < 10.0 ||
        energy > 800*units::MeV)
        dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 25*units::cm);
    dir_shower = dir_shower.norm();

    Vector dir_beam(0, 0, 1);

    // Muon kinetic energy from range (TPCParams::get_muon_r2ke() equivalent)
    auto muon_range_fn = ctx.particle_data->get_range_function("muon");
    double E_muon = muon_range_fn->scalar_function(muon_length / units::cm) * units::MeV;

    // Check muon-like segments at main_vertex: if any have more range than
    // the input muon_length, update E_muon
    if (ctx.main_vertex && ctx.main_vertex->descriptor_valid()) {
        auto vd = ctx.main_vertex->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || !sg1->has_particle_info()) continue;
            int pdg = sg1->particle_info()->pdg();
            if (pdg == 13 || pdg == 2212) {
                double length       = segment_track_length(sg1);
                double medium_dQ_dx = segment_median_dQ_dx(sg1);
                double dQ_dx_cut    = 0.8866 + 0.9533 * std::pow(18*units::cm / length, 0.4234);
                if (medium_dQ_dx < dQ_dx_cut * 43e3 / units::cm) {
                    double tmp_energy = muon_range_fn->scalar_function(length / units::cm) * units::MeV;
                    if (tmp_energy > E_muon) E_muon = tmp_energy;
                }
            }
        }
    }

    double tmp_shower_total_length = shower->get_total_length();

    if ((E_muon > energy && energy < 550*units::MeV) ||
        muon_length > tmp_shower_total_length ||
        muon_length > 80*units::cm ||
        (muon_length > 0.6 * tmp_shower_total_length && energy < 500*units::MeV))
        flag_bad = true;

    ti.cme_mu_energy  = E_muon / units::MeV;
    ti.cme_energy     = energy / units::MeV;
    ti.cme_mu_length  = muon_length / units::cm;
    ti.cme_length     = tmp_shower_total_length / units::cm;
    ti.cme_angle_beam = dir_beam.angle(dir_shower) / M_PI * 180.0;
    ti.cme_flag       = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// stem_direction
//
// Checks whether the shower stem segment direction is inconsistent with
// an electromagnetic shower: if the stem is significantly misaligned with
// the shower's PCA axis or the shower is nearly parallel to the drift
// direction, this can indicate a mis-reconstructed or hadronic event.
//
// Prototype: NeutrinoID_nue_functions.h, WCPPID::NeutrinoID::stem_direction()
// Fills: ti.stem_dir_*
//
// Translation note:
//   Prototype: main_cluster->Calc_PCA(shower_pts) + get_PCA_axis(0)
//   Toolkit:   calc_PCA_main_axis(shower_pts).second
//   Both compute the first PCA axis of the shower's point cloud.
//   The toolkit uses shower->fill_point_vector() to collect the same points.
// ===========================================================================
static bool stem_direction(NuEContext& ctx, ShowerPtr shower, double energy, TaggerInfo& ti) {
    bool flag_bad = false;

    Vector dir_drift(1, 0, 0);
    SegmentPtr sg = shower->start_segment();

    // PCA of shower points — mirrors prototype's main_cluster->Calc_PCA(tmp_pts)
    std::vector<Point> tmp_pts;
    shower->fill_point_vector(tmp_pts, /*flag_main=*/true);
    Vector dir1;
    if (!tmp_pts.empty())
        dir1 = ctx.self.calc_PCA_main_axis(tmp_pts).second;

    // Determine which end of sg is at main_vertex (geometric proximity)
    const auto& sg_fits  = sg->fits();
    Point sg_front       = sg_fits.front().point;
    Point sg_back        = sg_fits.back().point;
    Point mv_pt          = vtx_fit_pt(ctx.main_vertex);
    bool  front_is_mv    = (ray_length(Ray{mv_pt, sg_front}) <= ray_length(Ray{mv_pt, sg_back}));
    Point vertex_point   = front_is_mv ? sg_front : sg_back;

    Vector dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 100*units::cm);

    double angle  = 0;
    double angle1 = 0;
    double ratio  = 0;
    // angle2: how far the PCA axis deviates from the drift-perpendicular plane
    double angle2 = std::fabs(dir1.angle(dir_drift) / M_PI * 180.0 - 90.0);
    double angle3 = 0;

    if (front_is_mv) {
        Vector dir2 = segment_cal_dir_3vector(sg, sg_front, 5*units::cm);
        Vector dir3 = shower_cal_dir_3vector(*shower, sg_front, 30*units::cm);
        angle  = dir1.angle(dir2) / M_PI * 180.0;
        angle3 = dir2.angle(dir_shower) / M_PI * 180.0;
        if (angle > 90) angle = 180.0 - angle;
        angle1 = std::fabs(dir3.angle(dir_drift) / M_PI * 180.0 - 90.0);
        double len0_10 = segment_track_length(sg, 0, 0, 10);
        if (len0_10 > 0)
            ratio = segment_track_direct_length(sg, 0, 10) / len0_10;
    } else {
        Vector dir2 = segment_cal_dir_3vector(sg, sg_back, 5*units::cm);
        Vector dir3 = shower_cal_dir_3vector(*shower, sg_back, 30*units::cm);
        int num = (int)sg_fits.size() - 1;
        angle  = dir1.angle(dir2) / M_PI * 180.0;
        angle3 = dir2.angle(dir_shower) / M_PI * 180.0;
        if (angle > 90) angle = 180.0 - angle;
        angle1 = std::fabs(dir3.angle(dir_drift) / M_PI * 180.0 - 90.0);
        double len_nm10_n = segment_track_length(sg, 0, num - 10, num);
        if (len_nm10_n > 0)
            ratio = segment_track_direct_length(sg, num - 10, num) / len_nm10_n;
    }

    if (angle > 18) {
        if (energy > 1000*units::MeV) {
            // No stem-direction cut at very high energy
        } else if (energy > 500*units::MeV) {
            // High energy: cut on large misalignment + drift-parallel topology
            if (((angle1 > 12.5 || angle2 > 12.5) && angle > 25) ||
                ((angle1 > 10.0 || angle2 > 10.0) && angle > 32)) {
                if (angle3 > 3) flag_bad = true;  // 7006_293_14696
            }
        } else {
            if (angle > 25 && (angle1 > 7.5 || angle2 > 7.5)) {
                flag_bad = true;
            } else if ((angle1 > 7.5 || angle2 > 7.5) && ratio < 0.97) {
                flag_bad = true;
            }
        }
    }

    ti.stem_dir_angle  = angle;
    ti.stem_dir_energy = energy / units::MeV;
    ti.stem_dir_angle1 = angle1;
    ti.stem_dir_angle2 = angle2;
    ti.stem_dir_angle3 = angle3;
    ti.stem_dir_ratio  = ratio;
    ti.stem_dir_filled = 1.0f;
    ti.stem_dir_flag   = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// Forward declaration: bad_reconstruction (implemented in a later step).
// Called as a predicate from multiple_showers and other_showers; callers
// pass a local dummy TaggerInfo to avoid corrupting the real ti.
// ===========================================================================
static bool bad_reconstruction(NuEContext& ctx, ShowerPtr shower, TaggerInfo& ti);

// ===========================================================================
// pi0_identification
//
// Determines whether a shower is likely the photon from a pi0 decay.
// Two complementary checks:
//   flag_pi0_1: shower is already in the pi0 map (from shower_clustering_with_nv)
//               and the paired photon mass and asymmetry are consistent with pi0.
//   flag_pi0_2: no pi0 pair found, but a nearby cluster in the back-to-back
//               direction (relative to the shower stem) suggests a second photon.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::pi0_identification()
// Fills: ti.pio_*
//
// Note: when calling this function as a predicate from within other helpers
// (multiple_showers, other_showers), pass a local `TaggerInfo dummy{}` to
// avoid corrupting the real ti.  The prototype used flag_fill=false for this.
// ===========================================================================
static bool pi0_identification(NuEContext& ctx,
                               VertexPtr vertex, SegmentPtr sg, ShowerPtr shower,
                               double threshold, TaggerInfo& ti) {
    bool flag_pi0_1 = false;
    bool flag_pi0_2 = false;

    Point vertex_point = vtx_fit_pt(vertex);

    // Collect all vertices that already belong to known pi0-paired showers.
    // These are excluded from the flag_pi0_2 geometry search below.
    IndexedVertexSet used_vertices;
    for (auto& [shower1, pio_id] : ctx.map_shower_pio_id) {
        IndexedVertexSet vtxs;
        IndexedSegmentSet segs;
        shower1->fill_sets(vtxs, segs, /*flag_exclude_start_segment=*/false);
        used_vertices.insert(vtxs.begin(), vtxs.end());
    }

    auto it = ctx.map_shower_pio_id.find(shower);
    ti.pio_flag_pio = (it != ctx.map_shower_pio_id.end());

    if (it != ctx.map_shower_pio_id.end()) {
        // ----------------------------------------------------------------
        // flag_pi0_1 branch: shower is in a reconstructed pi0 pair.
        // ----------------------------------------------------------------
        auto& tmp_pi0_showers = ctx.map_pio_id_showers[it->second];
        auto  mass_pair       = ctx.map_pio_id_mass[it->second];

        double Eshower_1 = tmp_pi0_showers.front()->get_kine_charge();
        double Eshower_2 = tmp_pi0_showers.back()->get_kine_charge();

        double dis1 = ray_length(Ray{tmp_pi0_showers.front()->get_start_point(), vertex_point});
        double dis2 = ray_length(Ray{tmp_pi0_showers.back()->get_start_point(),  vertex_point});

        ti.pio_1_mass      = mass_pair.first  / units::MeV;
        ti.pio_1_pio_type  = mass_pair.second;
        ti.pio_1_energy_1  = Eshower_1 / units::MeV;
        ti.pio_1_energy_2  = Eshower_2 / units::MeV;
        ti.pio_1_dis_1     = dis1 / units::cm;
        ti.pio_1_dis_2     = dis2 / units::cm;

        bool mass_ok_1 = (std::fabs(mass_pair.first - 135*units::MeV) < 35*units::MeV &&
                          mass_pair.second == 1);
        bool mass_ok_2 = (std::fabs(mass_pair.first - 135*units::MeV) < 60*units::MeV &&
                          mass_pair.second == 2);

        if (mass_ok_1 || mass_ok_2) {
            // Symmetric photon pair
            if (std::min(Eshower_1, Eshower_2) > 15*units::MeV &&
                std::fabs(Eshower_1 - Eshower_2) / (Eshower_1 + Eshower_2) < 0.87)
                flag_pi0_1 = true;
            // Low-energy or balanced pair (6058_43_2166, 7017_364_18210)
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

        ti.pio_1_flag = !flag_pi0_1;

    } else {
        // ----------------------------------------------------------------
        // flag_pi0_2 branch: shower not in pi0 map.
        // Look for another cluster in the back-to-back direction that could
        // be the other photon from a pi0 decay.
        // ----------------------------------------------------------------
        Vector dir1 = segment_cal_dir_3vector(sg, vertex_point, 12*units::cm);

        if (dir1.magnitude() > 0) {
            // Precompute total track length per cluster (used to check acc_length > 0)
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

                Point   vtx1_pt = vtx_fit_pt(vtx1);
                Vector  dir2    = vtx1_pt - vertex_point;
                double  dis2    = dir2.magnitude();

                if (dis2 <= 0) continue;

                // Back-to-back angle: 180° - angle(dir1, dir2) < 7.5° means
                // vtx1 is nearly anti-parallel to the shower direction.
                double back_angle = 180.0 - dir1.angle(dir2) / M_PI * 180.0;

                if (dis2 < 36*units::cm && back_angle < 7.5 && acc_length > 0) {
                    flag_pi0_2 = true;
                    ti.pio_2_v_flag.push_back(0.0f);  // 0 = "this pair IS pi0"
                    ti.pio_2_v_dis2.push_back(dis2 / units::cm);
                    ti.pio_2_v_angle2.push_back(back_angle);
                    ti.pio_2_v_acc_length.push_back(acc_length / units::cm);
                } else {
                    ti.pio_2_v_flag.push_back(1.0f);  // 1 = "this pair is NOT pi0"
                    ti.pio_2_v_dis2.push_back(dis2 / units::cm);
                    ti.pio_2_v_angle2.push_back(back_angle);
                    ti.pio_2_v_acc_length.push_back(acc_length / units::cm);
                }
            }
        }
    }

    return flag_pi0_1 || flag_pi0_2;
}

// ===========================================================================
// single_shower
//
// Evaluates whether a single-shower topology passes geometric and dQ/dx
// quality cuts specific to the nue BDT.  The flag_single_shower parameter
// selects between two distinct cut sets:
//   true  — truly isolated single shower (no other vertex-connected tracks)
//   false — single shower but with vertex-connected tracks present
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::single_shower()
// Fills: ti.spt_*
//
// Translation notes:
//   shower->get_stem_dQ_dx() returns NORMALIZED dQ/dx (divided by 43e3/cm
//   internally), so comparison to 3.0, 2.4 is dimensionless — no extra factor.
//   map_vertex_segments[vertex].size() → boost::out_degree(vertex, graph).
// ===========================================================================
static bool single_shower(NuEContext& ctx, ShowerPtr shower,
                          bool flag_single_shower, TaggerInfo& ti) {
    bool flag_bad = false;

    Vector dir_beam(0, 0, 1);
    Vector dir_drift(1, 0, 0);
    Vector dir_vertical(0, 1, 0);

    VertexPtr  vertex = shower->get_start_vertex_and_type().first;
    SegmentPtr sg     = shower->start_segment();
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));

    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best()
                     : shower->get_kine_charge();

    // Normalized stem dQ/dx (first ≤3 fit points near the vertex)
    auto vec_dQ_dx = shower->get_stem_dQ_dx(vertex, sg, 20);
    double max_dQ_dx = 0;
    for (size_t i = 0; i < vec_dQ_dx.size(); ++i) {
        if (vec_dQ_dx[i] > max_dQ_dx) max_dQ_dx = vec_dQ_dx[i];
        if (i == 2) break;
    }

    // Primary shower direction
    Vector dir_shower;
    if (segment_track_length(sg) > 12*units::cm) {
        dir_shower = segment_cal_dir_3vector(sg, vertex_point, 15*units::cm);
    } else {
        dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
    }
    if (std::fabs(dir_shower.angle(dir_drift) / M_PI * 180.0 - 90.0) < 10.0 ||
        Eshower > 800*units::MeV)
        dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 25*units::cm);
    dir_shower = dir_shower.norm();

    // Secondary direction for angle_beam_1 / angle_drift_1 fills
    Vector dir_shower1 = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);

    double angle_beam     = dir_shower.angle(dir_beam)     / M_PI * 180.0;
    double angle_vertical = dir_vertical.angle(dir_shower) / M_PI * 180.0;
    double angle_drift    = std::fabs(M_PI/2.0 - dir_shower.angle(dir_drift)) / M_PI * 180.0;

    // Count valid tracks at the shower vertex (used for both cut paths)
    int    num_valid_tracks = 0;
    double max_length       = 0;

    if (vertex && vertex->descriptor_valid()) {
        auto vd = vertex->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;

            double medium_dQ_dx = segment_median_dQ_dx(sg1) / (43e3 / units::cm);
            double length       = segment_track_length(sg1);

            // 7022_110_5542: count non-shower-flagged segments with sufficient length/dQ_dx
            if (!seg_is_shower(sg1) &&
                (!sg1->dir_weak() || (sg1->dir_weak() && length > 4.2*units::cm) ||
                 (length > 0.6*units::cm && medium_dQ_dx > 3) ||
                 (length > 1.6*units::cm && medium_dQ_dx > 2.2))) {
                ++num_valid_tracks;
                if (length > max_length) max_length = length;
            } else {
                double dQ_dx_cut = 0.8866 + 0.9533 * std::pow(18*units::cm / length, 0.4234);
                if (medium_dQ_dx > dQ_dx_cut) {
                    ++num_valid_tracks;
                    if (length > max_length) max_length = length;
                }
            }
        }
    }

    double shower_main_length  = sg->cluster() ? shower->get_total_length(sg->cluster()) : 0;
    double shower_total_length = shower->get_total_length();
    size_t n_vtx_segs = vertex && vertex->descriptor_valid()
                        ? boost::out_degree(vertex->get_descriptor(), ctx.graph) : 0;

    if (flag_single_shower) {
        // 6572_18_948 + 7020_473_23679: shower almost entirely outside main cluster
        if (Eshower < 600*units::MeV &&
            ((shower_main_length < 0.1 * shower_total_length && angle_beam > 40) ||
             (angle_beam <= 40 && shower_main_length < 0.08 * shower_total_length)))
            flag_bad = true;
        // Nearly-vertical or nearly-horizontal shower with weak dQ/dx start
        if ((angle_vertical < 20 || angle_vertical > 160) &&
            angle_beam > 80 && max_dQ_dx < 3.0)
            flag_bad = true;
        // Drift-parallel shower with weak dQ/dx start (not forward-going)
        if ((angle_beam > 15 || dir_shower1.angle(dir_beam) / M_PI * 180.0 > 15) &&
            (angle_drift < 5 ||
             std::fabs(M_PI/2.0 - dir_shower1.angle(dir_drift)) / M_PI * 180.0 < 5) &&
            Eshower < 1200*units::MeV && max_dQ_dx < 2.4)
            flag_bad = true;
    } else {
        // Non-single-shower path: must have some valid track or be forward-going
        if (num_valid_tracks == 0 && angle_beam > 60 && n_vtx_segs <= 3)
            flag_bad = true;
        // 7017_969_48490: low-energy near-horizontal low-activity event
        if (Eshower < 200*units::MeV && angle_vertical < 10 && angle_drift < 5 &&
            max_length < 5*units::cm && num_valid_tracks <= 1)
            flag_bad = true;
    }

    ti.spt_flag_single_shower  = flag_single_shower;
    ti.spt_energy              = Eshower / units::MeV;
    ti.spt_shower_main_length  = shower_main_length / units::cm;
    ti.spt_shower_total_length = shower_total_length / units::cm;
    ti.spt_angle_beam          = angle_beam;
    ti.spt_angle_vertical      = angle_vertical;
    ti.spt_max_dQ_dx           = max_dQ_dx;
    ti.spt_angle_beam_1        = dir_shower1.angle(dir_beam) / M_PI * 180.0;
    ti.spt_angle_drift         = angle_drift;
    ti.spt_angle_drift_1       = std::fabs(M_PI/2.0 - dir_shower1.angle(dir_drift)) / M_PI * 180.0;
    ti.spt_num_valid_tracks    = num_valid_tracks;
    ti.spt_n_vtx_segs          = n_vtx_segs;
    ti.spt_max_length          = max_length / units::cm;
    ti.spt_flag                = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// multiple_showers
//
// Checks whether showers at the main vertex, or elsewhere in the event,
// carry enough energy to suggest that the "max_shower" candidate is actually
// a secondary photon in a multi-shower interaction rather than a primary
// nue electron.
//
// Prototype: NeutrinoID_nue_functions.h, WCPPID::NeutrinoID::multiple_showers()
// Fills: ti.mgo_*
//
// Note: pi0_identification and bad_reconstruction are called here as pure
// predicates; they receive a local dummy TaggerInfo to avoid corrupting ti.
// ===========================================================================
static bool multiple_showers(NuEContext& ctx, ShowerPtr max_shower,
                             double max_energy, TaggerInfo& ti) {
    bool flag_bad = false;

    // ------------------------------------------------------------------
    // First loop: showers connected to main_vertex with electron PDG and
    // connection type ≤ 1 (direct attachment).
    // ------------------------------------------------------------------
    double E_total     = 0;
    int    nshowers    = 0;
    double E_max_energy = 0;

    auto mv_it = ctx.map_vertex_to_shower.find(ctx.main_vertex);
    if (mv_it != ctx.map_vertex_to_shower.end()) {
        for (ShowerPtr shower : mv_it->second) {
            SegmentPtr sg = shower->start_segment();
            if (!sg || !sg->has_particle_info() || sg->particle_info()->pdg() != 11) continue;
            if (shower == max_shower) continue;
            auto [vtx, conn_type] = shower->get_start_vertex_and_type();
            if (conn_type > 1) continue;

            double E_shower = (shower->get_kine_best() != 0)
                              ? shower->get_kine_best() : shower->get_kine_charge();

            TaggerInfo dummy_ti{};
            bool flag_pi0 = pi0_identification(ctx, vtx, sg, shower, 15*units::MeV, dummy_ti);
            if (flag_pi0) continue;
            if (shower->get_total_length(sg->cluster()) < shower->get_total_length() * 0.1 &&
                shower->get_total_length(sg->cluster()) < 10*units::cm) continue;
            // 7010_532_26643
            if (segment_track_length(sg) > 80*units::cm) continue;

            if (E_shower > E_max_energy) E_max_energy = E_shower;
            if (E_shower > 50*units::MeV) {
                E_total += E_shower;
                ++nshowers;
            }
        }
    }

    if (E_max_energy > 0.6 * max_energy ||
        (E_max_energy > 0.45 * max_energy && max_energy - E_max_energy < 150*units::MeV))
        flag_bad = true;

    if ((E_total > 0.6 * max_energy ||
         (max_energy < 400*units::MeV && nshowers >= 2 && E_total > 0.3 * max_energy)) &&
        !flag_bad)
        flag_bad = true;

    // ------------------------------------------------------------------
    // Second loop: all showers in the event (not just at main_vertex),
    // excluding those in the main cluster.
    // ------------------------------------------------------------------
    double total_other_energy   = 0;
    double total_other_energy_1 = 0;
    int    total_num_showers    = 0;
    double E_max_energy_1       = 0;
    double E_max_energy_2       = 0;

    for (ShowerPtr shower : ctx.showers) {
        SegmentPtr sg = shower->start_segment();
        if (!sg || !sg->has_particle_info() || sg->particle_info()->pdg() != 11) continue;
        if (sg->cluster() == ctx.main_vertex->cluster()) continue;

        auto [vtx, conn_type] = shower->get_start_vertex_and_type();
        double E_shower = (shower->get_kine_best() != 0)
                          ? shower->get_kine_best() : shower->get_kine_charge();

        TaggerInfo dummy_ti{};
        if (bad_reconstruction(ctx, shower, dummy_ti)) continue;
        bool flag_pi0 = pi0_identification(ctx, vtx, sg, shower, 15*units::MeV, dummy_ti);

        if (flag_pi0) {
            // 7003_1682_84132: pi0 partner may share max_shower
            auto pio_it = ctx.map_shower_pio_id.find(shower);
            if (pio_it != ctx.map_shower_pio_id.end()) {
                auto& tmp_showers = ctx.map_pio_id_showers[pio_it->second];
                if (std::find(tmp_showers.begin(), tmp_showers.end(), max_shower) !=
                    tmp_showers.end()) {
                    if (E_shower > E_max_energy_1) E_max_energy_1 = E_shower;
                }
            }
        }

        if (flag_pi0) continue;

        if (conn_type <= 3) {
            total_other_energy += E_shower;
            if (vtx != ctx.main_vertex) total_other_energy_1 += E_shower;
            if (E_shower > 50*units::MeV) ++total_num_showers;
        }
        if (conn_type > 2) continue;
        if (E_shower > E_max_energy_2) E_max_energy_2 = E_shower;
    }

    if (E_max_energy_1 > max_energy * 0.75) flag_bad = true;

    // 7014_241_12058
    if (E_max_energy_2 > max_energy * 1.2 && max_energy < 250*units::MeV) flag_bad = true;

    if (!flag_bad &&
        ((max_energy < 250*units::MeV &&
          (total_other_energy - max_energy > 200*units::MeV ||
           (total_other_energy - max_energy > 60*units::MeV && total_num_showers >= 2))) ||
         (max_energy > 800*units::MeV && total_other_energy_1 > max_energy)))
        flag_bad = true;

    ti.mgo_energy              = max_energy / units::MeV;
    ti.mgo_max_energy          = E_max_energy / units::MeV;
    ti.mgo_total_energy        = E_total / units::MeV;
    ti.mgo_n_showers           = nshowers;
    ti.mgo_max_energy_1        = E_max_energy_1 / units::MeV;
    ti.mgo_max_energy_2        = E_max_energy_2 / units::MeV;
    ti.mgo_total_other_energy  = total_other_energy / units::MeV;
    ti.mgo_n_total_showers     = total_num_showers;
    ti.mgo_total_other_energy_1 = total_other_energy_1 / units::MeV;
    ti.mgo_flag                = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// other_showers
//
// Evaluates whether other showers in the event (outside the main cluster)
// are energetically or geometrically inconsistent with the max_shower being
// a primary nue electron.  Considers both direct (conn_type==1) and indirect
// (conn_type==2, within 72 cm) showers.
//
// Prototype: NeutrinoID_nue_functions.h, WCPPID::NeutrinoID::other_showers()
// Fills: ti.mgt_*
//
// Note: pi0_identification called as predicate; receives a local dummy TaggerInfo.
// ===========================================================================
static bool other_showers(NuEContext& ctx, ShowerPtr shower,
                          bool flag_single_shower, TaggerInfo& ti) {
    bool flag_bad = false;

    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    VertexPtr  vertex = shower->get_start_vertex_and_type().first;
    SegmentPtr sg     = shower->start_segment();
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));

    // ------------------------------------------------------------------
    // Quick survey: off-main-cluster showers for flag_single_shower path.
    // ------------------------------------------------------------------
    double total_other_energy = 0;
    double max_energy         = 0;

    for (ShowerPtr shower1 : ctx.showers) {
        SegmentPtr sg1 = shower1->start_segment();
        if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11) continue;
        if (sg1->cluster() == ctx.main_vertex->cluster()) continue;
        auto [vtx1, conn1] = shower1->get_start_vertex_and_type();
        double E_shower1 = (shower1->get_kine_best() != 0)
                           ? shower1->get_kine_best() : shower1->get_kine_charge();
        if (conn1 <= 3) total_other_energy += E_shower1;
        if (conn1 > 2) continue;
        if (E_shower1 > max_energy) max_energy = E_shower1;
    }

    if (flag_single_shower && max_energy > Eshower) flag_bad = true;
    if (flag_single_shower && Eshower < 150*units::MeV &&
        total_other_energy > 0.27 * Eshower) flag_bad = true;

    // ------------------------------------------------------------------
    // Detailed survey: classify competing showers as direct or indirect.
    // ------------------------------------------------------------------
    double E_direct_max_energy   = 0, E_direct_total_energy   = 0;
    double E_indirect_max_energy = 0, E_indirect_total_energy = 0;
    int    n_direct_showers      = 0;
    bool   flag_direct_max_pi0   = false;
    bool   flag_indirect_max_pi0 = false;
    double max_energy_1          = 0;

    for (ShowerPtr shower1 : ctx.showers) {
        SegmentPtr sg1 = shower1->start_segment();
        if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11) continue;
        if (shower1 == shower) continue;

        auto [vtx1, conn1] = shower1->get_start_vertex_and_type();
        double E_shower1 = (shower1->get_kine_best() != 0)
                           ? shower1->get_kine_best() : shower1->get_kine_charge();

        TaggerInfo dummy_ti{};
        bool flag_pi0 = pi0_identification(ctx, vtx1, sg1, shower1, 15*units::MeV, dummy_ti);

        if (flag_pi0) {
            // 7003_1682_84132: pi0 partner may share current shower
            auto pio_it = ctx.map_shower_pio_id.find(shower1);
            if (pio_it != ctx.map_shower_pio_id.end()) {
                auto& tmp_showers = ctx.map_pio_id_showers[pio_it->second];
                if (std::find(tmp_showers.begin(), tmp_showers.end(), shower) !=
                    tmp_showers.end()) {
                    if (E_shower1 > max_energy_1) max_energy_1 = E_shower1;
                }
            }
        }

        if (conn1 == 1) {
            // Direct: shower1 attached directly to some vertex
            // 7006_387_19382: skip if main cluster fraction is tiny
            if (shower1->get_total_length(sg1->cluster()) <
                    shower1->get_total_length() * 0.1 &&
                shower1->get_total_length(sg1->cluster()) < 10*units::cm) continue;
            // 7021_282_14130: skip long stem
            if (segment_track_length(sg1) > 80*units::cm) continue;
            // 6090_89_4498: skip pi0 if main cluster fraction tiny
            if (flag_pi0 &&
                shower1->get_total_length(sg1->cluster()) < 0.12 * shower1->get_total_length())
                continue;

            E_direct_total_energy += E_shower1;
            if (E_shower1 > E_direct_max_energy && vtx1 == ctx.main_vertex) {
                E_direct_max_energy = E_shower1;
                flag_direct_max_pi0 = ctx.map_shower_pio_id.count(shower1) > 0;
            }
            if (E_shower1 > 80*units::MeV) ++n_direct_showers;

        } else if (conn1 == 2) {
            // Indirect: shower1 connected via one intermediate vertex
            // Skip isolated non-shower-like single long segment (7021 pattern)
            if (shower1->get_num_segments() <= 2) {
                double sg1_len = segment_track_length(sg1);
                if (!sg1->flags_any(SegmentFlags::kShowerTrajectory) &&
                    !sg1->flags_any(SegmentFlags::kShowerTopology) &&
                    sg1_len > 45*units::cm &&
                    sg1_len > 0.95 * shower1->get_total_length()) continue;
            }
            double dis = ray_length(Ray{vertex_point, shower1->get_start_point()});
            // 6090_89_4498: too far away to be indirect partner
            if (dis > 72*units::cm) continue;
            double factor = (dis > 48*units::cm) ? 0.75 : 1.0;

            if (!flag_pi0) {
                E_indirect_total_energy += E_shower1;
                if (E_shower1 * factor > E_indirect_max_energy) {
                    E_indirect_max_energy = E_shower1 * factor;
                    flag_indirect_max_pi0 = ctx.map_shower_pio_id.count(shower1) > 0;
                }
            }
            if (E_shower1 > 80*units::MeV) { /* n_indirect_showers — not filled, keep for logic */ }
        }
    }

    if (max_energy_1 > Eshower * 0.75) flag_bad = true;
    if (E_indirect_max_energy > Eshower + 350*units::MeV ||
        E_direct_max_energy   > Eshower) flag_bad = true;
    if (Eshower < 1000*units::MeV && n_direct_showers > 0 &&
        E_direct_max_energy > 0.33 * Eshower) flag_bad = true;
    if (Eshower >= 1000*units::MeV && n_direct_showers > 0 &&
        E_direct_max_energy > 0.33 * Eshower &&
        E_direct_total_energy > 900*units::MeV) flag_bad = true;

    // 6748_57_2867 + 7004_1604_80229
    if (flag_indirect_max_pi0) {
        if (Eshower < 800*units::MeV &&
            E_indirect_total_energy - E_indirect_max_energy > Eshower &&
            E_indirect_max_energy > 0.5 * Eshower) flag_bad = true;
    } else {
        if (Eshower < 800*units::MeV &&
            E_indirect_total_energy > Eshower * 0.6 &&
            E_indirect_max_energy   > 0.5 * Eshower) flag_bad = true;
    }

    ti.mgt_flag_single_shower      = flag_single_shower;
    ti.mgt_max_energy              = max_energy / units::MeV;
    ti.mgt_energy                  = Eshower / units::MeV;
    ti.mgt_total_other_energy      = total_other_energy / units::MeV;
    ti.mgt_max_energy_1            = max_energy_1 / units::MeV;
    ti.mgt_e_indirect_max_energy   = E_indirect_max_energy / units::MeV;
    ti.mgt_e_direct_max_energy     = E_direct_max_energy / units::MeV;
    ti.mgt_n_direct_showers        = n_direct_showers;
    ti.mgt_e_direct_total_energy   = E_direct_total_energy / units::MeV;
    ti.mgt_e_indirect_total_energy = E_indirect_total_energy / units::MeV;
    ti.mgt_flag_indirect_max_pio   = flag_indirect_max_pi0;
    ti.mgt_flag                    = !flag_bad;

    return flag_bad;
}

// ===========================================================================
// vertex_inside_shower
//
// Detects two failure modes where vertex topology suggests the "shower" is
// actually a kink or nuclear interaction rather than an EM shower:
//
//   flag_bad1 (vis_1): A non-shower segment at the vertex is nearly
//     anti-parallel to the shower but has similar length — looks like a
//     broken track rather than an e/m shower.
//
//   flag_bad2 (vis_2): A segment at the vertex is nearly collinear with the
//     shower direction (nearly back-to-back after min over two reference
//     directions), combined with weak direction determination or high dQ/dx —
//     indicative of a crossing track or broken muon.
//
// Prototype: NeutrinoID_nue_functions.h, WCPPID::NeutrinoID::vertex_inside_shower()
// Fills: ti.vis_1_*, ti.vis_2_*, ti.vis_flag
//
// Porting note: prototype line 456 assigns `max_sg = sg` instead of
// `max_sg = sg1` — this is a prototype bug that is faithfully reproduced
// here to preserve BDT input values.
// ===========================================================================
static bool vertex_inside_shower(NuEContext& ctx, ShowerPtr shower, TaggerInfo& ti) {
    bool flag_bad1 = false;
    bool flag_bad2 = false;

    Vector dir_drift(1, 0, 0);
    Vector dir_beam(0, 0, 1);

    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    VertexPtr  vertex = shower->get_start_vertex_and_type().first;
    SegmentPtr sg     = shower->start_segment();
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));

    // ------------------------------------------------------------------
    // Block 1: check for segments nearly anti-parallel to the shower.
    // ------------------------------------------------------------------
    {
        Vector dir1 = shower_cal_dir_3vector(*shower, vertex_point, 30*units::cm);

        double     max_angle  = 0;
        SegmentPtr max_sg     = nullptr;   // NOTE: set to sg (not sg1) — see porting note
        int        num_good_tracks = 0;

        if (vertex && vertex->descriptor_valid()) {
            auto vd = vertex->get_descriptor();
            for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
                SegmentPtr sg1 = ctx.graph[*eit].segment;
                if (!sg1 || sg1 == sg) continue;
                Vector dir2    = segment_cal_dir_3vector(sg1, vertex_point, 15*units::cm);
                double angle   = dir2.angle(dir1) / M_PI * 180.0;
                if (!seg_is_shower(sg1) && !sg1->dir_weak()) ++num_good_tracks;
                if (angle > max_angle && segment_track_length(sg1) > 1.0*units::cm) {
                    max_angle = angle;
                    max_sg = sg;    // prototype assigns sg (start seg), not sg1
                }
            }
        }

        // Competing showers: max angle from other showers relative to shower direction
        double max_shower_angle = 0;
        Vector dir_long = shower_cal_dir_3vector(*shower, vertex_point, 100*units::cm);
        for (ShowerPtr shower1 : ctx.showers) {
            if (shower1 == shower) continue;
            if (shower1->get_start_vertex_and_type().second > 2) continue;
            double energy = (shower1->get_kine_best() != 0)
                            ? shower1->get_kine_best() : shower1->get_kine_charge();
            Point  sp1  = shower1->get_start_point();
            Vector dir2 = sp1 - vertex_point;
            Vector dir3 = shower_cal_dir_3vector(*shower1, sp1, 100*units::cm);
            if (energy > 30*units::MeV && energy > 0.2 * Eshower &&
                dir2.angle(dir3) / M_PI * 180.0 < 20) {
                double angle1 = dir_long.angle(dir3) / M_PI * 180.0;
                if (max_shower_angle < angle1) max_shower_angle = angle1;
            }
        }

        if (max_sg) {
            // tmp_length1 = max_sg length; max_sg = sg (start seg) per prototype
            double tmp_length1 = segment_track_length(max_sg);
            double tmp_length2 = segment_track_length(sg);
            size_t n_vtx_segs  = vertex && vertex->descriptor_valid()
                                 ? boost::out_degree(vertex->get_descriptor(), ctx.graph) : 0;

            if (n_vtx_segs >= 3 && Eshower < 500*units::MeV && num_good_tracks == 0 &&
                ((max_angle > 150 &&
                  (tmp_length1 < 15*units::cm || tmp_length2 < 15*units::cm) &&
                  std::max(tmp_length1, tmp_length2) < 25*units::cm) ||
                 ((max_angle > 170 || (max_shower_angle > 170 && max_angle > 120)) &&
                  (tmp_length1 < 25*units::cm || tmp_length2 < 25*units::cm) &&
                  std::max(tmp_length1, tmp_length2) < 35*units::cm)))
                flag_bad1 = true;
            else if (n_vtx_segs == 2 && Eshower < 500*units::MeV && num_good_tracks == 0 &&
                     max_angle > 150 &&
                     (sg->has_particle_info() && sg->particle_info()->pdg() == 13) &&
                     (tmp_length1 < 35*units::cm || tmp_length2 < 35*units::cm))
                flag_bad1 = true;

            ti.vis_1_filled          = 1.0f;
            ti.vis_1_n_vtx_segs      = n_vtx_segs;
            ti.vis_1_energy          = Eshower / units::MeV;
            ti.vis_1_num_good_tracks = num_good_tracks;
            ti.vis_1_max_angle       = max_angle;
            ti.vis_1_max_shower_angle = max_shower_angle;
            ti.vis_1_tmp_length1     = tmp_length1 / units::cm;
            ti.vis_1_tmp_length2     = tmp_length2 / units::cm;
            ti.vis_1_particle_type   = sg->has_particle_info()
                                       ? sg->particle_info()->pdg() : 0;
            ti.vis_1_flag            = !flag_bad1;
        }
    }

    // ------------------------------------------------------------------
    // Block 2: check for nearly-collinear segments (broken track topology).
    // ------------------------------------------------------------------
    {
        size_t n_vtx_segs = vertex && vertex->descriptor_valid()
                            ? boost::out_degree(vertex->get_descriptor(), ctx.graph) : 0;

        if (n_vtx_segs > 1) {
            // Shower direction (same logic as single_shower)
            Vector dir_shower;
            if (segment_track_length(sg) > 12*units::cm) {
                dir_shower = segment_cal_dir_3vector(sg, vertex_point, 15*units::cm);
            } else {
                dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
            }
            if (std::fabs(dir_shower.angle(dir_drift) / M_PI * 180.0 - 90.0) < 10.0 ||
                Eshower > 800*units::MeV)
                dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 25*units::cm);
            dir_shower = dir_shower.norm();

            // Stem direction (short range, for collinearity test)
            Vector dir2_sg = segment_cal_dir_3vector(sg, vertex_point, 6*units::cm);

            double max_angle = 0, max_angle1 = 0, max_length = 0, max_medium_dQ_dx = 0;
            int    max_weak_track = 0;
            double min_angle = 180, min_angle1 = 0, min_length = 0, min_medium_dQ_dx = 0;
            int    min_weak_track = 0;

            if (vertex->descriptor_valid()) {
                auto vd = vertex->get_descriptor();
                for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
                    SegmentPtr sg1 = ctx.graph[*eit].segment;
                    if (!sg1 || sg1 == sg) continue;
                    Vector dir1 = segment_cal_dir_3vector(sg1, vertex_point, 15*units::cm);

                    // Minimum of "how back-to-back vs shower" and "how back-to-back vs stem"
                    double angle = std::min(
                        180.0 - dir1.angle(dir_shower) / M_PI * 180.0,
                        180.0 - dir1.angle(dir2_sg)   / M_PI * 180.0);

                    double norm_dQ_dx = segment_median_dQ_dx(sg1) / (43e3 / units::cm);
                    double length     = segment_track_length(sg1);
                    bool   is_weak    = sg1->dir_weak();

                    if (angle > max_angle) {
                        max_angle        = angle;
                        max_angle1       = std::fabs(M_PI/2.0 - dir1.angle(dir_drift)) / M_PI * 180.0;
                        max_weak_track   = is_weak ? 1 : 0;
                        max_length       = length;
                        max_medium_dQ_dx = norm_dQ_dx;
                    }
                    if (angle < min_angle) {
                        min_angle        = angle;
                        min_angle1       = std::fabs(M_PI/2.0 - dir1.angle(dir_drift)) / M_PI * 180.0;
                        min_weak_track   = is_weak ? 1 : 0;
                        min_length       = length;
                        min_medium_dQ_dx = norm_dQ_dx;
                    }
                }
            }

            double iso_angle1   = std::fabs(M_PI/2.0 - dir_drift.angle(dir_shower)) / M_PI * 180.0;
            double angle_beam   = dir_beam.angle(dir_shower) / M_PI * 180.0;

            // 6090_85_4300
            if (n_vtx_segs == 2 &&
                ((min_angle < 25 && min_weak_track == 1) || min_angle < 20) &&
                angle_beam > 50)
                flag_bad2 = true;

            if (n_vtx_segs == 2 && min_angle < 70 && min_angle1 < 10 &&
                iso_angle1 < 10 && (iso_angle1 + min_angle1) < 15 &&
                min_medium_dQ_dx > 1.5 && min_medium_dQ_dx < 2.2) {
                flag_bad2 = true;
                // 7001_100_5003: short segment at large angle with long start seg is OK
                if (min_length < 4*units::cm && min_angle > 45 &&
                    segment_track_length(sg) > 30*units::cm)
                    flag_bad2 = false;
            }

            // 7003_1740_87003
            if (n_vtx_segs == 3 &&
                ((min_angle < 15 && min_medium_dQ_dx < 2.1) ||
                 (min_angle < 17.5 && min_length < 5.0*units::cm && min_medium_dQ_dx < 2.5)) &&
                ((min_weak_track == 1 && max_angle > 120) ||
                 (min_length < 6*units::cm && max_angle > 135 &&
                  min_angle < 12.5 && max_weak_track == 1))) {
                // 7004_8_428: exception for long non-weak max segment
                if (max_length > 40*units::cm && max_weak_track == 0) {
                    // no flag
                } else {
                    flag_bad2 = true;
                }
            }

            if (n_vtx_segs == 3 && min_angle < 5 && min_medium_dQ_dx < 2.1 &&
                min_length < 10*units::cm && max_angle > 90 && max_weak_track == 1)
                flag_bad2 = true;

            if (n_vtx_segs == 3 && min_angle < 35 && min_angle1 < 10 &&
                iso_angle1 < 10 && (iso_angle1 + min_angle1) < 15 &&
                min_medium_dQ_dx < 2.1 && min_weak_track == 1 && max_angle > 120)
                flag_bad2 = true;

            ti.vis_2_filled          = 1.0f;
            ti.vis_2_n_vtx_segs      = n_vtx_segs;
            ti.vis_2_min_angle       = min_angle;
            ti.vis_2_min_weak_track  = min_weak_track;
            ti.vis_2_angle_beam      = angle_beam;
            ti.vis_2_min_angle1      = min_angle1;
            ti.vis_2_iso_angle1      = iso_angle1;
            ti.vis_2_min_medium_dQ_dx = min_medium_dQ_dx;
            ti.vis_2_min_length      = min_length / units::cm;
            ti.vis_2_sg_length       = segment_track_length(sg) / units::cm;
            ti.vis_2_max_angle       = max_angle;
            ti.vis_2_max_weak_track  = max_weak_track;
            ti.vis_2_flag            = !flag_bad2;
        }
    }

    ti.vis_flag = !(flag_bad1 || flag_bad2);
    return flag_bad1 || flag_bad2;
}

// ===========================================================================
// broken_muon_id
//
// Follows the shower stem forward to check if the shower is actually a
// broken (gap-crossing) muon track.  Walks from the start segment outward,
// greedily chaining nearly-collinear segments (first within the shower, then
// across cluster gaps), and checks whether the resulting track length and
// straightness are more consistent with a muon than an EM shower.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::broken_muon_id()
// Fills: ti.brm_*
//
// Translation notes:
//   shower->get_map_seg_vtxs() → shower->fill_sets() to get all shower segs
//   shower->get_map_vtx_segs() → same fill_sets, then iterate out_edges
//     filtered to shower-internal segments
//   find_vertices(seg)         → find_vertices(ctx.graph, seg)
//   find_other_vertex(seg, v)  → find_other_vertex(ctx.graph, seg, v)
//   sg->get_direct_length()    → segment_track_direct_length(sg)  [defaults to full seg]
//   tmp_ids (set<cluster_id>)  → tmp_clusters (set<Facade::Cluster*>)
//   add_length                 → dropped (only in debug print, not in fills or cuts)
// ===========================================================================
static bool broken_muon_id(NuEContext& ctx, ShowerPtr shower, TaggerInfo& ti) {
    bool flag_bad = false;

    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    VertexPtr  vertex = shower->get_start_vertex_and_type().first;
    SegmentPtr sg     = shower->start_segment();
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vertex));

    // Pre-fill shower internal segment/vertex sets (replaces get_map_seg_vtxs / get_map_vtx_segs)
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    // -----------------------------------------------------------------
    // Walk the muon track: follow nearly-collinear segments forward from
    // the shower start, chaining through connected and nearby segments.
    // -----------------------------------------------------------------
    Vector dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);

    std::set<SegmentPtr> muon_segments;
    SegmentPtr curr_seg = sg;
    VertexPtr  curr_vtx = find_other_vertex(ctx.graph, curr_seg, vertex);
    muon_segments.insert(curr_seg);

    bool flag_continue = true;
    while (flag_continue) {
        flag_continue = false;

        Point  curr_vtx_pt = vtx_fit_pt(curr_vtx);
        Vector dir1 = segment_cal_dir_3vector(curr_seg, curr_vtx_pt, 15*units::cm);

        // Step A: look for a shower-internal connected segment that continues
        //         nearly collinearly (back-to-back, within 15°).
        if (curr_vtx && curr_vtx->descriptor_valid()) {
            auto vd = curr_vtx->get_descriptor();
            for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
                SegmentPtr sg1 = ctx.graph[*eit].segment;
                if (!sg1 || !shower_segs.count(sg1)) continue;
                if (muon_segments.count(sg1)) continue;
                Vector dir2 = segment_cal_dir_3vector(sg1, curr_vtx_pt, 15*units::cm);
                // back-to-back: 180° - angle < 15° AND length > 6cm
                if (180.0 - dir1.angle(dir2) / M_PI * 180.0 < 15 &&
                    segment_track_length(sg1) > 6*units::cm) {
                    flag_continue  = true;
                    curr_seg = sg1;
                    curr_vtx = find_other_vertex(ctx.graph, sg1, curr_vtx);
                    break;
                }
            }
        }

        // Step B: if no connected continuation found, look for a nearby segment
        //         in a different cluster (gap crossing).
        double min_dis  = 1e9;
        SegmentPtr min_seg = nullptr;
        VertexPtr  min_vtx_found = nullptr;

        if (!flag_continue) {
            for (SegmentPtr sg1 : shower_segs) {
                // Skip segments whose cluster is already in the muon track
                bool skip = false;
                for (SegmentPtr mseg : muon_segments) {
                    if (sg1->cluster() == mseg->cluster()) { skip = true; break; }
                }
                if (skip) continue;

                const auto& fits1 = sg1->fits();
                if (fits1.empty()) continue;
                Point front1 = fits1.front().point;
                Point back1  = fits1.back().point;

                double dis1 = ray_length(Ray{curr_vtx_pt, front1});
                double dis2 = ray_length(Ray{curr_vtx_pt, back1});

                // Choose the nearer endpoint
                Point  near_pt = (dis1 < dis2) ? front1 : back1;
                double near_dis = std::min(dis1, dis2);
                Vector dir2 = near_pt - curr_vtx_pt;
                Vector dir3 = segment_cal_dir_3vector(sg1, near_pt, 15*units::cm);

                double angle1 = 180.0 - dir1.angle(dir2) / M_PI * 180.0;
                double angle2 = dir2.angle(dir3)          / M_PI * 180.0;
                double angle3 = 180.0 - dir1.angle(dir3)  / M_PI * 180.0;

                bool close_collinear = ((std::min(angle1, angle2) < 10 &&
                                         angle1 + angle2 < 25) ||
                                        (angle3 < 15 && near_dis < 5*units::cm)) &&
                                       near_dis < 25*units::cm;
                bool strict_collinear = (std::min(angle1, angle2) < 5 &&
                                         angle1 + angle2 < 15) ||
                                        (angle3 < 10 && near_dis < 5*units::cm);
                bool far_collinear   = std::min(angle1, angle2) < 15 && angle3 < 30 &&
                                       near_dis > 30*units::cm &&
                                       segment_track_length(sg1) > 25*units::cm &&
                                       near_dis < 60*units::cm;

                bool passes = close_collinear ||
                              (strict_collinear && near_dis < 30*units::cm) ||
                              far_collinear;

                if (passes && near_dis < min_dis) {
                    min_dis = near_dis;
                    min_seg = sg1;
                    // Pick the farther endpoint as the new forward vertex
                    auto pair_vtxs = find_vertices(ctx.graph, min_seg);
                    double d3 = ray_length(Ray{curr_vtx_pt, vtx_fit_pt(pair_vtxs.first)});
                    double d4 = ray_length(Ray{curr_vtx_pt, vtx_fit_pt(pair_vtxs.second)});
                    min_vtx_found = (d4 > d3) ? pair_vtxs.second : pair_vtxs.first;
                }
            }

            if (min_seg) {
                flag_continue = true;
                curr_seg = min_seg;
                curr_vtx = min_vtx_found;
            }
        }

        if (flag_continue) {
            muon_segments.insert(curr_seg);
            // add_length (gap distance) intentionally dropped — only used in debug print
        }
    } // while

    // -----------------------------------------------------------------
    // Accumulate track properties over all muon segments.
    // -----------------------------------------------------------------
    double acc_length        = 0;
    double acc_direct_length = 0;
    std::set<Facade::Cluster*> tmp_clusters;

    for (SegmentPtr mseg : muon_segments) {
        acc_length        += segment_track_length(mseg);
        acc_direct_length += segment_track_direct_length(mseg);
        tmp_clusters.insert(mseg->cluster());
    }

    auto muon_range_fn = ctx.particle_data->get_range_function("muon");
    double Ep = muon_range_fn->scalar_function(acc_length / units::cm) * units::MeV;

    // Connected length: total length of shower-internal segments in muon clusters.
    // 7020_348_17421
    double connected_length = 0;
    for (SegmentPtr sg1 : shower_segs) {
        if (tmp_clusters.count(sg1->cluster()))
            connected_length += segment_track_length(sg1);
    }

    // 7022_42_2123: add segments in the main cluster that are nearly parallel
    //               to the shower direction but not already in muon_segments.
    {
        Vector dir_sg = segment_cal_dir_3vector(sg, vertex_point, 15*units::cm);
        for (SegmentPtr sg1 : shower_segs) {
            if (muon_segments.count(sg1)) continue;
            if (sg1->cluster() != sg->cluster()) continue;
            auto pair_vtxs = find_vertices(ctx.graph, sg1);
            Point pt1 = vtx_fit_pt(pair_vtxs.first);
            Point pt2 = vtx_fit_pt(pair_vtxs.second);
            Vector d1 = segment_cal_dir_3vector(sg1, pt1, 15*units::cm);
            Vector d2 = segment_cal_dir_3vector(sg1, pt2, 15*units::cm);
            double a1 = std::min(d1.angle(dir_sg) / M_PI * 180.0,
                                 180.0 - d1.angle(dir_sg) / M_PI * 180.0);
            double a2 = std::min(d2.angle(dir_sg) / M_PI * 180.0,
                                 180.0 - d2.angle(dir_sg) / M_PI * 180.0);
            if (std::min(a1, a2) < 10) muon_segments.insert(sg1);
        }
    }

    int num_muon_main = 0;
    for (SegmentPtr mseg : muon_segments) {
        if (mseg->cluster() == sg->cluster()) ++num_muon_main;
    }

    // Primary cut: multi-cluster muon track with sufficient straightness/length
    // at low shower energy.
    if (muon_segments.size() > 1 &&
        (Ep > Eshower * 0.55 ||
         acc_length > 0.65 * shower->get_total_length() ||
         connected_length > 0.95 * shower->get_total_length()) &&
        tmp_clusters.size() > 1 &&
        acc_direct_length > 0.94 * acc_length &&
        Eshower < 350*units::MeV) {
        // 7004_989_49482: cut only if shower is simple and muon dominates
        if (shower->get_num_main_segments() <= 3 &&
            shower->get_num_main_segments() - num_muon_main < 2 &&
            (shower->get_num_segments() < (int)muon_segments.size() + 6 ||
             acc_length > connected_length * 0.9 ||
             acc_length > 0.8 * shower->get_total_length()))
            flag_bad = true; // 6640_173_8673
    }

    ti.brm_n_mu_segs            = muon_segments.size();
    ti.brm_Ep                   = Ep / units::MeV;
    ti.brm_energy               = Eshower / units::MeV;
    ti.brm_acc_length           = acc_length / units::cm;
    ti.brm_shower_total_length  = shower->get_total_length() / units::cm;
    ti.brm_connected_length     = connected_length / units::cm;
    ti.brm_n_size               = tmp_clusters.size();
    ti.brm_acc_direct_length    = acc_direct_length / units::cm;
    ti.brm_n_shower_main_segs   = shower->get_num_segments();
    ti.brm_n_mu_main            = num_muon_main;
    ti.brm_flag                 = !flag_bad;

    return flag_bad;
}
