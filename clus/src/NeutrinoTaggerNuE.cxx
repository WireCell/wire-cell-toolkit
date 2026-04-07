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

// ===========================================================================
// mip_quality
//
// Checks two failure modes:
//   flag_overlap: the first few fit points of the shower stem overlap (in 2D
//     wire space) with another shower-internal segment — likely a crossing track.
//   flag_split: two electron-like showers at the vertex with no tracks — the
//     event looks like one shower was split into two by a reconstruction gap.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::mip_quality()
// Fills: ti.mip_quality_*
//
// Translation notes:
//   sg->get_closest_2d_dis(pt)       → segment_get_closest_2d_distances(sg,pt,ctx.apa,ctx.face)
//   map_vertex_segments[main_vertex] → boost::out_edges on ctx.main_vertex, check seg==sg1
//   n_protons                        → computed but never used in fills/cuts; dropped
// ===========================================================================
static bool mip_quality(NuEContext& ctx,
                        VertexPtr vertex, SegmentPtr sg, ShowerPtr shower,
                        TaggerInfo& ti) {
    bool flag_overlap = false;
    bool flag_split   = false;

    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    // ------------------------------------------------------------------
    // Overlap check: first ≤3 fit points near vertex vs shower-internal segs
    // ------------------------------------------------------------------
    {
        const auto& sg_fits = sg->fits();
        bool vertex_at_front = !sg_fits.empty() &&
            (ray_length(Ray{vtx_fit_pt(vertex), sg_fits.front().point}) <=
             ray_length(Ray{vtx_fit_pt(vertex), sg_fits.back().point}));

        // Collect indices of the ≤3 test fit points from the vertex end
        std::vector<int> test_indices;
        int n_fits = (int)sg_fits.size();
        if (vertex_at_front) {
            for (int i = 0; i < std::min(n_fits, 3); ++i) test_indices.push_back(i);
        } else {
            for (int i = n_fits - 1; i >= std::max(0, n_fits - 3); --i)
                test_indices.push_back(i);
        }

        // Far-end connectivity (used for exception: don't flag shared vertex)
        VertexPtr other_vtx = find_other_vertex(ctx.graph, sg, vertex);
        size_t nconnected = other_vtx && other_vtx->descriptor_valid()
                            ? boost::out_degree(other_vtx->get_descriptor(), ctx.graph) : 0;

        // Shower-internal segments for 2D overlap check
        IndexedSegmentSet shower_segs;
        IndexedVertexSet  shower_vtxs_tmp;
        shower->fill_sets(shower_vtxs_tmp, shower_segs, /*flag_exclude_start_segment=*/false);

        for (size_t k = 0; k < test_indices.size(); ++k) {
            int   orig_idx = test_indices[k];
            Point pt       = sg_fits[orig_idx].point;

            double min_u = 1e9, min_v = 1e9, min_w = 1e9;
            for (SegmentPtr sg1 : shower_segs) {
                if (sg1 == sg) continue;
                auto [u, v, w] = segment_get_closest_2d_distances(sg1, pt, ctx.apa, ctx.face);
                if (u < min_u) min_u = u;
                if (v < min_v) min_v = v;
                if (w < min_w) min_w = w;
            }

            if (min_u < 0.3*units::cm && min_v < 0.3*units::cm && min_w < 0.3*units::cm) {
                // Exception: skip if this is the vertex endpoint (shared vertex = 0,0,0 is OK)
                bool is_vertex_end = (k == 0);
                // Exception: skip if this is the other endpoint and it's a shared vertex
                bool is_other_end = vertex_at_front
                                    ? (orig_idx == n_fits - 1)
                                    : (orig_idx == 0);
                if (is_vertex_end && min_u == 0 && min_v == 0 && min_w == 0) {
                    // 7017_617_30888: don't flag the shared vertex endpoint
                } else if (is_other_end && min_u == 0 && min_v == 0 && min_w == 0 &&
                           nconnected == 2) {
                    // shared other-end vertex
                } else {
                    flag_overlap = true;
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Split check: two electron showers at vertex with no tracks
    // ------------------------------------------------------------------
    int n_showers = 0;
    int n_tracks  = 0;
    std::set<ShowerPtr> connected_showers;
    std::set<ShowerPtr> tmp_pi0_showers;

    {
        auto mv_it = ctx.map_vertex_to_shower.find(vertex);
        if (mv_it != ctx.map_vertex_to_shower.end()) {
            for (ShowerPtr shower1 : mv_it->second) {
                SegmentPtr sg1 = shower1->start_segment();
                if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11)
                    continue;
                // Check if sg1 is connected at main_vertex
                bool at_main = false;
                if (ctx.main_vertex && ctx.main_vertex->descriptor_valid()) {
                    for (auto [eit,eend] = boost::out_edges(
                             ctx.main_vertex->get_descriptor(), ctx.graph);
                         eit != eend; ++eit) {
                        if (ctx.graph[*eit].segment == sg1) { at_main = true; break; }
                    }
                }
                if (at_main) {
                    ++n_showers;
                    connected_showers.insert(shower1);
                }
                if (ctx.pi0_showers.count(shower1)) tmp_pi0_showers.insert(shower1);
            }
        }
        // Count non-shower tracks at vertex
        if (vertex && vertex->descriptor_valid()) {
            for (auto [eit,eend] = boost::out_edges(vertex->get_descriptor(), ctx.graph);
                 eit != eend; ++eit) {
                SegmentPtr sg1 = ctx.graph[*eit].segment;
                if (sg1 && !seg_is_shower(sg1)) ++n_tracks;
            }
        }
    }

    bool   flag_inside_pi0    = false;
    double shortest_length    = 1e9;
    double shortest_acc_length = 0;
    double shortest_angle     = 0;
    bool   flag_proton        = false;

    if (n_showers == 2 && n_tracks == 0) {
        flag_split = true;

        Point v_pt = vtx_fit_pt(vertex);
        Vector dir1 = shower_cal_dir_3vector(*shower, v_pt, 6*units::cm);

        for (ShowerPtr shower1 : connected_showers) {
            if (shower1 == shower) continue;
            SegmentPtr sg1  = shower1->start_segment();
            double norm_dQ  = segment_median_dQ_dx(sg1) / (43e3 / units::cm);
            double length1  = segment_track_length(sg1);
            double dQ_cut   = 0.8866 + 0.9533 * std::pow(18*units::cm / length1, 0.4234);
            if (norm_dQ > dQ_cut) { flag_split = false; flag_proton = true; }
            if (tmp_pi0_showers.count(shower1)) flag_inside_pi0 = true;

            Vector dir2 = shower_cal_dir_3vector(*shower1, v_pt, 6*units::cm);
            if (length1 < shortest_length) {
                shortest_length    = length1;
                shortest_angle     = dir1.angle(dir2) / M_PI * 180.0;
                if (std::isnan(shortest_angle)) shortest_angle = 0;
                shortest_acc_length = shower1->get_total_length(sg1->cluster());
            }
        }
        // 7004_365_18300
        if (!flag_inside_pi0 && !tmp_pi0_showers.empty()) flag_split = false;
        // 7010_1076_53830
        if (((shortest_angle > 45 && shortest_length > 20*units::cm) ||
             (shortest_angle > 35 && shortest_acc_length > 40*units::cm)) &&
            shortest_length < 1e9)
            flag_split = false;
    }

    bool flag_bad = ((Eshower < 800*units::MeV) && flag_overlap) ||
                    ((Eshower < 500*units::MeV) && flag_split);

    ti.mip_quality_flag          = !flag_bad;
    ti.mip_quality_energy        = Eshower / units::MeV;
    ti.mip_quality_overlap       = flag_overlap;
    ti.mip_quality_n_showers     = n_showers;
    ti.mip_quality_n_tracks      = n_tracks;
    ti.mip_quality_flag_inside_pi0 = flag_inside_pi0;
    ti.mip_quality_n_pi0_showers = tmp_pi0_showers.size();
    ti.mip_quality_shortest_length = (shortest_length > 10*units::m)
                                     ? 1000.0f
                                     : shortest_length / units::cm;
    ti.mip_quality_acc_length    = shortest_acc_length / units::cm;
    ti.mip_quality_shortest_angle = shortest_angle;
    ti.mip_quality_flag_proton   = flag_proton;
    ti.mip_quality_filled        = 1;

    return flag_bad;
}

// ===========================================================================
// high_energy_overlapping
//
// Checks for two types of overlap between the shower stem and other segments
// at the vertex (connection type == 1 only):
//
//   flag_overlap1: the minimum angle between the stem and any other segment
//     is small, with no valid tracks — likely an overlapping cluster.
//   flag_overlap2: consecutive fit points near the vertex are within 0.6 cm
//     of the most collinear segment, and that segment has high dQ/dx in the
//     overlap region — likely a collinear hadronic track.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::high_energy_overlapping()
// Fills: ti.hol_1_*, ti.hol_2_*, ti.hol_flag
//
// Translation notes:
//   (*it)->get_closest_point(pt).first → segment_get_closest_point(sg1, pt).first
//   min_sg->get_medium_dQ_dx(n1,n2)   → segment_median_dQ_dx(min_sg, n1, n2)
//   flag_start (vertex at front)       → geometric proximity check
// ===========================================================================
static bool high_energy_overlapping(NuEContext& ctx, ShowerPtr shower, TaggerInfo& ti) {
    bool flag_overlap1 = false;
    bool flag_overlap2 = false;

    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    VertexPtr  vtx = shower->get_start_vertex_and_type().first;
    SegmentPtr sg  = shower->start_segment();
    int conn_type  = shower->get_start_vertex_and_type().second;

    // Stem dQ/dx (normalized, first ≤3 fit points)
    auto vec_dQ_dx = shower->get_stem_dQ_dx(vtx, sg, 20);
    double max_dQ_dx = 0;
    for (size_t i = 0; i < vec_dQ_dx.size(); ++i) {
        if (vec_dQ_dx[i] > max_dQ_dx) max_dQ_dx = vec_dQ_dx[i];
        if (i == 2) break;
    }

    const auto& sg_fits = sg->fits();
    bool flag_start = !sg_fits.empty() &&
        (ray_length(Ray{vtx_fit_pt(vtx), sg_fits.front().point}) <=
         ray_length(Ray{vtx_fit_pt(vtx), sg_fits.back().point}));
    Point vtx_point = flag_start ? sg_fits.front().point : sg_fits.back().point;

    if (conn_type == 1 && vtx && vtx->descriptor_valid()) {
        // --------------------------------------------------------------
        // flag_overlap1: n_valid_tracks and min_angle at vertex
        // 7012_1195_59764 + 7017_1158_57929
        // --------------------------------------------------------------
        Vector dir1 = segment_cal_dir_3vector(sg, vtx_point, 15*units::cm);
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
                Point dir2_pt = vtx_point;
                Vector dir2 = segment_cal_dir_3vector(sg1, dir2_pt, 5*units::cm);
                if (dir2.magnitude() == 0) { flag_all_showers = false; continue; }
                double angle = dir1.angle(dir2) / M_PI * 180.0;
                if (angle < min_angle) { min_angle = angle; min_length = segment_track_length(sg1); }
            } else {
                flag_all_showers = false;
            }
            double norm_dQ = segment_median_dQ_dx(sg1) / (43e3 / units::cm);
            bool   is_proton = sg1->has_particle_info() && sg1->particle_info()->pdg() == 2212;
            if ((!sg1->dir_weak() || is_proton || segment_track_length(sg1) > 20*units::cm) &&
                !seg_is_shower(sg1))
                ++n_valid_tracks;
            else if (norm_dQ > 2.0 && segment_track_length(sg1) > 1.8*units::cm)
                ++n_valid_tracks; // 7010_20_1012
        }

        // Count other electron showers at this vertex
        int num_showers = 0;
        auto mv_it = ctx.map_vertex_to_shower.find(vtx);
        if (mv_it != ctx.map_vertex_to_shower.end()) {
            for (ShowerPtr shower1 : mv_it->second) {
                if (shower1 == shower) continue;
                SegmentPtr sg1 = shower1->start_segment();
                if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11)
                    continue;
                auto [vtx1, conn1] = shower1->get_start_vertex_and_type();
                double E1 = (shower1->get_kine_best() != 0)
                            ? shower1->get_kine_best() : shower1->get_kine_charge();
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

        ti.hol_1_n_valid_tracks  = n_valid_tracks;
        ti.hol_1_min_angle       = min_angle;
        ti.hol_1_energy          = Eshower / units::MeV;
        ti.hol_1_flag_all_shower = flag_all_showers;
        ti.hol_1_min_length      = min_length / units::cm;
        ti.hol_1_flag            = !flag_overlap1;

        // --------------------------------------------------------------
        // flag_overlap2: count consecutive close fit points and get dQ/dx
        // --------------------------------------------------------------
        {
            Vector dir1_8 = segment_cal_dir_3vector(sg, vtx_point, 8*units::cm);
            double min_ang2  = 180;
            SegmentPtr min_sg = nullptr;

            for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
                SegmentPtr sg1 = ctx.graph[*eit].segment;
                if (!sg1 || sg1 == sg) continue;
                Point dir2_pt = vtx_point;
                Vector dir2 = segment_cal_dir_3vector(sg1, dir2_pt, 5*units::cm);
                if (dir2.magnitude() == 0) continue;
                double ang = dir1_8.angle(dir2) / M_PI * 180.0;
                if (ang < min_ang2) { min_ang2 = ang; min_sg = sg1; }
            }

            // Count consecutive fit points within 0.6 cm of any other vertex segment
            int ncount = 0;
            const auto& pts = sg_fits;
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
            if (flag_start) iterate_pts(pts.begin(),  pts.end());
            else            iterate_pts(pts.rbegin(), pts.rend());

            double medium_dQ_dx = 0;
            if (min_sg) {
                const auto& min_fits = min_sg->fits();
                int n_min = (int)min_fits.size();
                bool min_front_near =
                    (ray_length(Ray{vtx_point, min_fits.front().point}) <=
                     ray_length(Ray{vtx_point, min_fits.back().point}));
                if (min_front_near)
                    medium_dQ_dx = segment_median_dQ_dx(min_sg, 0, ncount)
                                   / (43e3 / units::cm);
                else
                    medium_dQ_dx = segment_median_dQ_dx(min_sg, n_min-1-ncount, n_min-1)
                                   / (43e3 / units::cm);
            }

            if (min_ang2 < 15 && medium_dQ_dx > 0.95 && ncount > 5  && Eshower < 1500*units::MeV)
                flag_overlap2 = true;
            if (min_ang2 < 7.5 && medium_dQ_dx > 0.8  && ncount > 8  && Eshower < 1500*units::MeV)
                flag_overlap2 = true;
            if (min_ang2 < 5   && ncount > 12 && medium_dQ_dx > 0.5  && Eshower < 1500*units::MeV)
                flag_overlap2 = true;

            ti.hol_2_min_angle    = min_ang2;
            ti.hol_2_medium_dQ_dx = medium_dQ_dx;
            ti.hol_2_ncount       = ncount;
            ti.hol_2_energy       = Eshower / units::MeV;
            ti.hol_2_flag         = !flag_overlap2;
        }
    }

    bool flag_overlap = flag_overlap1 || flag_overlap2;
    ti.hol_flag = !flag_overlap;
    return flag_overlap;
}

// ===========================================================================
// low_energy_overlapping
//
// Three complementary checks for low-energy shower overlap:
//   flag_overlap_1: two shower-internal segments at a main-cluster vertex form
//     a small opening angle — likely a split track not a real shower fork.
//   flag_overlap_2: a short collinear muon-type segment at the vertex is nearly
//     parallel to the stem — a crossing muon stub.
//   flag_overlap_3: backward-going shower with no valid tracks, or isolated
//     shower cluster with too many outward-pointing hits.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::low_energy_overlapping()
// Fills: ti.lol_1_v_*, ti.lol_2_v_*, ti.lol_3_*, ti.lol_flag
//
// Translation notes:
//   shower->get_map_vtx_segs() / get_map_seg_vtxs() → fill_sets()
//   shower interior fit points for n_out → sg1->fits() interior (i=1..size-2)
//   map_vertex_segments[vtx].size()     → boost::out_degree on vtx
// ===========================================================================
static bool low_energy_overlapping(NuEContext& ctx, ShowerPtr shower, TaggerInfo& ti) {
    bool flag_overlap_1_save = false;
    bool flag_overlap_2_save = false;
    bool flag_overlap_3      = false;

    Vector dir_beam(0, 0, 1);
    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    VertexPtr  vtx      = shower->get_start_vertex_and_type().first;
    SegmentPtr sg       = shower->start_segment();
    int        conn_type = shower->get_start_vertex_and_type().second;
    Point vtx_point = seg_endpoint_near(sg, vtx_fit_pt(vtx));

    // Pre-fill shower internal sets
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    // Number of main-cluster shower segments
    int nseg = 0;
    for (SegmentPtr sg1 : shower_segs) {
        if (sg1->cluster() == sg->cluster()) ++nseg;
    }

    // Stem direction for angle and n_out computations
    Point vtx_pt_copy = vtx_point;
    Vector dir1_stem = segment_cal_dir_3vector(sg, vtx_pt_copy, 5*units::cm);
    double angle_beam = dir1_stem.angle(dir_beam) / M_PI * 180.0;

    // n_valid_tracks and min_angle at shower start vertex
    int    n_valid_tracks = 0;
    double min_angle_vtx  = 180;
    size_t n_vtx_segs_global = vtx && vtx->descriptor_valid()
                               ? boost::out_degree(vtx->get_descriptor(), ctx.graph) : 0;

    if (vtx && vtx->descriptor_valid()) {
        auto vd = vtx->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            Point dir2_pt = vtx_point;
            Vector dir2 = segment_cal_dir_3vector(sg1, dir2_pt, 5*units::cm);
            double tmp_angle = dir2.angle(dir1_stem) / M_PI * 180.0;
            if (tmp_angle < min_angle_vtx) min_angle_vtx = tmp_angle;
            bool is_proton = sg1->has_particle_info() && sg1->particle_info()->pdg() == 2212;
            if ((!sg1->dir_weak() || is_proton || segment_track_length(sg1) > 20*units::cm) &&
                !seg_is_shower(sg1))
                ++n_valid_tracks;
        }
    }

    // n_out / n_sum: shower hits outside the shower direction cone (15° half-angle)
    int n_sum = 0, n_out = 0;
    {
        Vector dir1_15 = segment_cal_dir_3vector(sg, vtx_pt_copy, 15*units::cm);

        // Shower-internal vertices
        for (VertexPtr vtx1 : shower_vtxs) {
            Vector dir2 = vtx_fit_pt(vtx1) - vtx_point;
            ++n_sum;
            if (dir1_15.angle(dir2) / M_PI * 180.0 > 15) ++n_out;
        }
        // Interior fit points of all shower-internal segments
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
    // flag_overlap_1: shower-internal vertices with 2 connected shower segs
    //   forming a small opening angle.
    // ------------------------------------------------------------------
    for (VertexPtr vtx1 : shower_vtxs) {
        if (!vtx1->cluster() || vtx1->cluster() != sg->cluster()) continue;
        if (!vtx1->descriptor_valid()) continue;

        // Collect shower-internal segments at vtx1
        std::vector<SegmentPtr> vtx_ss;
        for (auto [eit, eend] = boost::out_edges(vtx1->get_descriptor(), ctx.graph);
             eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (sg1 && shower_segs.count(sg1)) vtx_ss.push_back(sg1);
        }
        if (vtx_ss.empty()) continue;

        Point vp = vtx_fit_pt(vtx1);
        // Use first and last of vtx_ss (mirrors prototype's begin/rbegin of ordered set)
        Point d1p = vp, d2p = vp;
        Vector dv1 = segment_cal_dir_3vector(vtx_ss.front(), d1p, 5*units::cm);
        Vector dv2 = segment_cal_dir_3vector(vtx_ss.back(),  d2p, 5*units::cm);
        double open_angle = dv1.angle(dv2) / M_PI * 180.0;

        bool flag_ov1 = false;
        if (vtx_ss.size() == 2 &&
            open_angle < 36 && nseg == 2 && Eshower < 150*units::MeV &&
            n_vtx_segs_global == 1)
            flag_ov1 = true;

        ti.lol_1_v_energy.push_back(Eshower / units::MeV);
        ti.lol_1_v_vtx_n_segs.push_back(n_vtx_segs_global);
        ti.lol_1_v_nseg.push_back(nseg);
        ti.lol_1_v_angle.push_back(open_angle);
        ti.lol_1_v_flag.push_back(!flag_ov1);

        if (flag_ov1) flag_overlap_1_save = true;
    }

    // ------------------------------------------------------------------
    // flag_overlap_2: short collinear muon/weak segment at vertex
    // ------------------------------------------------------------------
    if (vtx && vtx->descriptor_valid()) {
        auto vd = vtx->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            bool flag_ov2 = false;
            Point dir2_pt = vtx_point;
            Vector dir2 = segment_cal_dir_3vector(sg1, dir2_pt, 5*units::cm);
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

            ti.lol_2_v_flag.push_back(!flag_ov2);
            ti.lol_2_v_length.push_back(len1 / units::cm);
            ti.lol_2_v_angle.push_back(ang2);
            ti.lol_2_v_type.push_back(sg1->has_particle_info()
                                       ? sg1->particle_info()->pdg() : 0);
            ti.lol_2_v_vtx_n_segs.push_back(n_vtx_segs_global);
            ti.lol_2_v_energy.push_back(Eshower / units::MeV);
            ti.lol_2_v_shower_main_length.push_back(main_len / units::cm);
            ti.lol_2_v_flag_dir_weak.push_back(sg1->dir_weak());

            if (flag_ov2) flag_overlap_2_save = true;
        }
    }

    // ------------------------------------------------------------------
    // flag_overlap_3: backward/isolated low-energy shower
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

    ti.lol_3_flag               = !flag_overlap_3;
    ti.lol_3_angle_beam         = angle_beam;
    ti.lol_3_min_angle          = min_angle_vtx;
    ti.lol_3_n_valid_tracks     = n_valid_tracks;
    ti.lol_3_vtx_n_segs         = n_vtx_segs_global;
    ti.lol_3_energy             = Eshower / units::MeV;
    ti.lol_3_shower_main_length = main_length / units::cm;
    ti.lol_3_n_sum              = n_sum;
    ti.lol_3_n_out              = n_out;

    bool flag_overlap = flag_overlap_1_save || flag_overlap_2_save || flag_overlap_3;
    ti.lol_flag = !flag_overlap;
    return flag_overlap;
}

// ===========================================================================
// bad_reconstruction_1
//
// Checks whether the shower stem direction is inconsistent with the shower's
// PCA axis (suggesting the shower is actually a mis-classified track).
// Nearly identical to stem_direction but with a different cut structure:
//   - Gated on flag_single_shower || num_valid_tracks == 0
//   - Extra unconditional cuts at 40° misalignment (7012) and for long
//     shower-trajectory stems (7012)
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::bad_reconstruction_1()
// Fills: ti.br2_*
// ===========================================================================
static bool bad_reconstruction_1(NuEContext& ctx, ShowerPtr shower,
                                  bool flag_single_shower, int num_valid_tracks,
                                  TaggerInfo& ti) {
    bool flag_bad = false;

    Vector dir_drift(1, 0, 0);
    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    SegmentPtr sg     = shower->start_segment();
    VertexPtr  vertex = shower->get_start_vertex_and_type().first;

    // PCA axis of shower points
    std::vector<Point> tmp_pts;
    shower->fill_point_vector(tmp_pts, /*flag_main=*/true);
    Vector dir1;
    if (!tmp_pts.empty()) dir1 = ctx.self.calc_PCA_main_axis(tmp_pts).second;

    // Determine which end of sg is at the main vertex
    const auto& sg_fits    = sg->fits();
    bool front_is_mv = !sg_fits.empty() &&
        (ray_length(Ray{vtx_fit_pt(ctx.main_vertex), sg_fits.front().point}) <=
         ray_length(Ray{vtx_fit_pt(ctx.main_vertex), sg_fits.back().point}));
    Point vertex_point = front_is_mv ? sg_fits.front().point : sg_fits.back().point;

    Vector dir_shower = shower_cal_dir_3vector(*shower, vertex_point, 100*units::cm);

    double angle  = 0;
    double angle1 = 0;
    double angle2 = std::fabs(dir1.angle(dir_drift) / M_PI * 180.0 - 90.0);
    double angle3 = 0;

    if (front_is_mv) {
        Point p = sg_fits.front().point;
        Vector dir2 = segment_cal_dir_3vector(sg, p, 5*units::cm);
        Vector dir3 = shower_cal_dir_3vector(*shower, p, 30*units::cm);
        angle  = dir1.angle(dir2) / M_PI * 180.0;
        if (angle > 90) angle = 180.0 - angle;
        angle1 = std::fabs(dir3.angle(dir_drift) / M_PI * 180.0 - 90.0);
        angle3 = dir_shower.angle(dir2) / M_PI * 180.0;
    } else {
        Point p = sg_fits.back().point;
        Vector dir2 = segment_cal_dir_3vector(sg, p, 5*units::cm);
        Vector dir3 = shower_cal_dir_3vector(*shower, p, 30*units::cm);
        angle  = dir1.angle(dir2) / M_PI * 180.0;
        if (angle > 90) angle = 180.0 - angle;
        angle1 = std::fabs(dir3.angle(dir_drift) / M_PI * 180.0 - 90.0);
        angle3 = dir_shower.angle(dir2) / M_PI * 180.0;
    }

    // Max angle at the far end of the stem (other_vertex)
    double max_angle = 0;
    VertexPtr other_vtx = find_other_vertex(ctx.graph, sg, vertex);
    if (other_vtx && other_vtx->descriptor_valid()) {
        Point ovp = vtx_fit_pt(other_vtx);
        Vector dir_1 = segment_cal_dir_3vector(sg, ovp, 10*units::cm);
        auto ov_vd = other_vtx->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(ov_vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            Point dir2_pt = ovp;
            Vector dir_2 = segment_cal_dir_3vector(sg1, dir2_pt, 10*units::cm);
            double a2 = dir_1.angle(dir_2) / M_PI * 180.0;
            if (a2 > max_angle) max_angle = a2;
        }
    }

    // Energy-dependent PCA misalignment cut
    if (flag_single_shower || num_valid_tracks == 0) {
        if (Eshower > 1000*units::MeV) {
            // No cut at very high energy
        } else if (Eshower > 500*units::MeV) {
            if ((angle1 > 10 || angle2 > 10) && angle > 30) {
                if (angle3 > 3) flag_bad = true;
            }
        } else {
            if ((angle > 25 && shower->get_num_main_segments() > 1 || angle > 30) &&
                (angle1 > 7.5 || angle2 > 7.5))
                flag_bad = true;
        }
    }

    // 7012_922_46106: unconditional large misalignment cut
    if (angle > 40 && (angle1 > 7.5 || angle2 > 7.5) && max_angle < 100)
        flag_bad = true;
    // 7012_785_39252: long shower-trajectory stem with misalignment
    if (angle > 20 && (angle1 > 7.5 || angle2 > 7.5) &&
        segment_track_length(sg) > 21*units::cm && Eshower < 600*units::MeV &&
        sg->flags_any(SegmentFlags::kShowerTrajectory))
        flag_bad = true;

    ti.br2_flag                = !flag_bad;
    ti.br2_flag_single_shower  = flag_single_shower;
    ti.br2_num_valid_tracks    = num_valid_tracks;
    ti.br2_energy              = Eshower / units::MeV;
    ti.br2_angle1              = angle1;
    ti.br2_angle2              = angle2;
    ti.br2_angle               = angle;
    ti.br2_angle3              = angle3;
    ti.br2_n_shower_main_segs  = shower->get_num_main_segments();
    ti.br2_max_angle           = max_angle;
    ti.br2_sg_length           = segment_track_length(sg) / units::cm;
    ti.br2_flag_sg_trajectory  = sg->flags_any(SegmentFlags::kShowerTrajectory);

    return flag_bad;
}

// ===========================================================================
// shower_to_wall
//
// Checks whether the shower points backward toward the detector wall (or
// another cluster/shower) rather than away from it, as expected for a
// genuine nue CC electron shower.  Four sub-flags:
//
//   flag_bad1: backward distance to wall < threshold (stw_1)
//   flag_bad2: another non-electron shower lies in the backward direction (stw_2)
//   flag_bad3: another cluster vertex lies in the backward direction (stw_3)
//   flag_bad4: another forward shower's end is within 3cm of the wall (stw_4)
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::shower_to_wall()
// Fills: ti.stw_1_*, ti.stw_2_v_*, ti.stw_3_v_*, ti.stw_4_v_*, ti.stw_flag
//
// Translation notes:
//   fid->inside_fiducial_volume(p, offset_x, &tol) → fiducial_utils->inside_fiducial_volume(p, tol)
//   shower->get_particle_type()                     → shower->get_particle_type()  (direct)
//   shower->get_end_point()                         → shower->get_end_point()      (direct)
//   stw_3 global vertex loop                        → graph_nodes(ctx.graph)
//   vtx1->get_wcpt().{x,y,z}                        → vtx1->wcpt().point.{x,y,z}
// ===========================================================================
static bool shower_to_wall(NuEContext& ctx, ShowerPtr shower,
                           double shower_energy, bool flag_single_shower,
                           TaggerInfo& ti) {
    bool flag_bad1      = false;
    bool flag_bad2_save = false;
    bool flag_bad3_save = false;
    bool flag_bad4_save = false;

    SegmentPtr sg     = shower->start_segment();
    VertexPtr  vertex = shower->get_start_vertex_and_type().first;

    // Tolerance: shrink all fiducial boundaries by 1.5 cm for wall-proximity walk
    const std::vector<double> stm_tol_vec(5, -1.5*units::cm);

    // Vertex end of the start segment + near-vertex dQ/dx (normalized)
    const auto& sg_fits = sg->fits();
    bool vertex_at_front = !sg_fits.empty() &&
        (ray_length(Ray{vtx_fit_pt(vertex), sg_fits.front().point}) <=
         ray_length(Ray{vtx_fit_pt(vertex), sg_fits.back().point}));
    Point vertex_point = vertex_at_front ? sg_fits.front().point : sg_fits.back().point;
    int n_fits = (int)sg_fits.size();
    double medium_dQ_dx = vertex_at_front
        ? segment_median_dQ_dx(sg, 0, 6)            / (43e3 / units::cm)
        : segment_median_dQ_dx(sg, n_fits-1-6, n_fits-1) / (43e3 / units::cm);

    // Stem dQ/dx (normalized, first ≤3 fit points)
    auto vec_dQ_dx = shower->get_stem_dQ_dx(vertex, sg, 20);
    double max_dQ_dx = 0;
    for (size_t i = 0; i < vec_dQ_dx.size(); ++i) {
        if (vec_dQ_dx[i] > max_dQ_dx) max_dQ_dx = vec_dQ_dx[i];
        if (i == 2) break;
    }

    // Backward directions from vertex (two ranges for stw_2 and stw_3)
    Vector shower_dir_15 = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
    Vector dir  = shower_dir_15.norm() * (-1.0);      // primary backward unit vector (15 cm)
    Vector dir2 = shower_cal_dir_3vector(*shower, vertex_point, 6*units::cm).norm() * (-1.0); // secondary (6 cm)

    // Walk backward from vertex_point to find distance to wall
    FiducialUtilsPtr fiducial_utils;
    if (ctx.main_cluster && ctx.main_cluster->grouping())
        fiducial_utils = ctx.main_cluster->grouping()->get_fiducialutils();

    double dis = 0;
    if (fiducial_utils) {
        const double step = 1*units::cm;
        Point test_p = vertex_point + dir * step;
        while (fiducial_utils->inside_fiducial_volume(test_p, stm_tol_vec)) {
            test_p = test_p + dir * step;
        }
        dis = ray_length(Ray{vertex_point, test_p});
    }

    // ------------------------------------------------------------------
    // flag_bad1: backward distance to wall
    // ------------------------------------------------------------------
    if (shower_energy < 300*units::MeV && dis < 15*units::cm && max_dQ_dx < 2.6 &&
        flag_single_shower)
        flag_bad1 = true;

    // Count other electron showers and pi0 showers at this vertex
    int n_other_shower = 0, n_pi0 = 0;
    {
        auto mv_it = ctx.map_vertex_to_shower.find(vertex);
        if (mv_it != ctx.map_vertex_to_shower.end()) {
            for (ShowerPtr shower1 : mv_it->second) {
                if (shower1 == shower) continue;
                if (shower1->get_start_vertex_and_type().second > 2) continue;
                if (shower1->get_particle_type() != 11) continue;
                double E1 = (shower1->get_kine_best() != 0)
                            ? shower1->get_kine_best() : shower1->get_kine_charge();
                if (E1 > 60*units::MeV) ++n_other_shower;
                if (E1 > 25*units::MeV &&
                    shower1->get_start_vertex_and_type().first == vertex &&
                    ctx.pi0_showers.count(shower1))
                    ++n_pi0;
            }
        }
    }

    // 7023_669_33467
    if ((n_pi0 < 2 || shower_energy > 1000*units::MeV) &&
        dis < 5*units::cm && flag_single_shower)
        flag_bad1 = true;

    // Count valid tracks at vertex (for non-single-shower path)
    int num_valid_tracks = 0;
    if (vertex && vertex->descriptor_valid()) {
        for (auto [eit, eend] = boost::out_edges(vertex->get_descriptor(), ctx.graph);
             eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || sg1 == sg) continue;
            if (!seg_is_shower(sg1) &&
                (!sg1->dir_weak() || segment_track_length(sg1) > 5*units::cm))
                ++num_valid_tracks;
        }
    }
    if (num_valid_tracks == 0 && dis < 3*units::cm && !flag_single_shower)
        flag_bad1 = true;

    ti.stw_1_energy             = shower_energy / units::MeV;
    ti.stw_1_dis                = dis / units::cm;
    ti.stw_1_dQ_dx              = max_dQ_dx;
    ti.stw_1_flag_single_shower = flag_single_shower;
    ti.stw_1_n_pi0              = n_pi0;
    ti.stw_1_num_valid_tracks   = num_valid_tracks;
    ti.stw_1_flag               = !flag_bad1;

    // ------------------------------------------------------------------
    // flag_bad2: other non-electron showers in backward direction
    // ------------------------------------------------------------------
    for (ShowerPtr shower1 : ctx.showers) {
        if (shower1 == shower) continue;
        if (shower1->get_total_length() == 0) continue;
        if (shower1->get_particle_type() == 11) continue; // only non-electron showers

        Point sp1 = shower1->get_start_point();
        Vector dir1 = sp1 - vertex_point;
        double dir1_mag = dir1.magnitude();

        double min_ang = std::min(dir1.angle(dir)  / M_PI * 180.0,
                                  dir1.angle(dir2) / M_PI * 180.0);
        bool flag_bad2 = false;
        if ((medium_dQ_dx > 1.3 || shower_energy < 300*units::MeV) &&
            min_ang < 15 && dir1_mag < 40*units::cm &&
            max_dQ_dx < 3.0 && flag_single_shower)
            flag_bad2 = true;

        ti.stw_2_v_medium_dQ_dx.push_back(medium_dQ_dx);
        ti.stw_2_v_energy.push_back(shower_energy / units::MeV);
        ti.stw_2_v_angle.push_back(min_ang);
        ti.stw_2_v_dir_length.push_back(dir1_mag / units::cm);
        ti.stw_2_v_max_dQ_dx.push_back(max_dQ_dx);
        ti.stw_2_v_flag.push_back(!flag_bad2);

        if (flag_bad2) flag_bad2_save = true;
    }

    // ------------------------------------------------------------------
    // flag_bad3: vertices of other clusters in backward direction
    // ------------------------------------------------------------------
    for (const auto& vd : graph_nodes(ctx.graph)) {
        VertexPtr vtx1 = ctx.graph[vd].vertex;
        if (!vtx1) continue;
        if (vtx1->cluster() == vertex->cluster()) continue;
        if (!vtx1->descriptor_valid()) continue;
        size_t n_conn = boost::out_degree(vd, ctx.graph);
        if (n_conn == 1) continue;

        // Use wcpt (wire cell point) as prototype does, not fit point
        Vector dir1 = vtx1->wcpt().point - vertex_point;
        double min_ang = std::min(dir1.angle(dir)  / M_PI * 180.0,
                                  dir1.angle(dir2) / M_PI * 180.0);
        bool flag_bad3 = false;
        if (min_ang < 15 && dir1.magnitude() < 40*units::cm &&
            (shower_energy < 300*units::MeV || medium_dQ_dx > 1.3) &&
            flag_single_shower)
            flag_bad3 = true;

        ti.stw_3_v_angle.push_back(min_ang);
        ti.stw_3_v_dir_length.push_back(dir1.magnitude() / units::cm);
        ti.stw_3_v_energy.push_back(shower_energy / units::MeV);
        ti.stw_3_v_medium_dQ_dx.push_back(medium_dQ_dx);
        ti.stw_3_v_flag.push_back(!flag_bad3);

        if (flag_bad3) flag_bad3_save = true;
    }

    // ------------------------------------------------------------------
    // flag_bad4: other vertex-connected showers whose end is near the wall
    // 7018_885_44275
    // ------------------------------------------------------------------
    {
        auto mv_it = ctx.map_vertex_to_shower.find(vertex);
        if (mv_it != ctx.map_vertex_to_shower.end()) {
            for (ShowerPtr shower1 : mv_it->second) {
                if (shower1 == shower) continue;
                if (shower1->get_start_vertex_and_type().second > 2) continue;

                Vector dir1 = shower_cal_dir_3vector(*shower1, shower1->get_start_point(),
                                                     15*units::cm);
                if (dir1.magnitude() == 0) continue;
                dir1 = dir1.norm();

                double dis1 = 0;
                if (dir1.angle(dir) / M_PI * 180.0 < 30 && fiducial_utils) {
                    const double step = 1*units::cm;
                    Point end_pt = shower1->get_end_point();
                    Point test_p = end_pt;
                    while (fiducial_utils->inside_fiducial_volume(test_p, stm_tol_vec)) {
                        test_p = test_p + dir1 * step;
                    }
                    dis1 = ray_length(Ray{end_pt, test_p});
                    if (dis1 < 3*units::cm && shower_energy < 500*units::MeV &&
                        flag_single_shower)
                        flag_bad4_save = true;
                }

                ti.stw_4_v_angle.push_back(dir1.angle(dir) / M_PI * 180.0);
                ti.stw_4_v_dis.push_back(dis1 / units::cm);
                ti.stw_4_v_energy.push_back(shower_energy / units::MeV);
                ti.stw_4_v_flag.push_back(!flag_bad4_save);
            }
        }
    }

    bool flag_bad = flag_bad1 || flag_bad2_save || flag_bad3_save || flag_bad4_save;
    ti.stw_flag = !flag_bad;
    return flag_bad;
}

// ===========================================================================
// single_shower_pio_tagger
//
// Identifies events where the main shower is likely a pi0 photon rather than
// a nue CC electron via two sub-checks:
//
//   flag_bad1: an indirect (conn_type==2) electron shower starts near the
//     main-vertex fit point and points away from it — likely a second pi0 γ.
//
//   flag_bad2: the segment at the far end of the shower (furthest point along
//     shower direction) has dQ/dx consistent with a track stub rather than EM,
//     while the shower start has low dQ/dx — suggests the "shower" is really
//     a broken track with a hadronic tail.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::single_shower_pio_tagger()
// Fills: ti.sig_1_v_*, ti.sig_2_v_*, ti.sig_flag
//
// Translation notes:
//   sg->get_dQ_vec() / get_dx_vec() → sg->fits()[i].dQ / fits()[i].dx
//   map_vtx_segs[max_vtx]           → boost::out_edges filtered to shower_segs
//   dir1.Dot(dir)                   → dir1.dot(dir)
// ===========================================================================
static bool single_shower_pio_tagger(NuEContext& ctx, ShowerPtr shower,
                                     bool flag_single_shower, TaggerInfo& ti) {
    bool flag_bad1_save = false;
    bool flag_bad2      = false;

    Vector dir_beam(0, 0, 1);
    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    VertexPtr  vtx = shower->get_start_vertex_and_type().first;
    SegmentPtr sg  = shower->start_segment();
    Point vertex_point = seg_endpoint_near(sg, vtx_fit_pt(vtx));

    Vector shower_dir = shower_cal_dir_3vector(*shower, vertex_point, 15*units::cm);
    double shower_angle = shower_dir.angle(dir_beam) / M_PI * 180.0;

    // ------------------------------------------------------------------
    // flag_bad1: indirect electron shower pointing away from vertex
    // ------------------------------------------------------------------
    auto mv_it = ctx.map_vertex_to_shower.find(vtx);
    if (mv_it != ctx.map_vertex_to_shower.end()) {
        for (ShowerPtr shower1 : mv_it->second) {
            SegmentPtr sg1 = shower1->start_segment();
            if (!sg1 || !sg1->has_particle_info() || sg1->particle_info()->pdg() != 11)
                continue;
            if (shower1 == shower) continue;
            if (shower1->get_start_vertex_and_type().second > 2) continue;

            double E1 = (shower1->get_kine_best() != 0)
                        ? shower1->get_kine_best() : shower1->get_kine_charge();
            bool flag_bad1 = false;

            if (shower1->get_start_vertex_and_type().second == 2) {
                Point  sp1  = shower1->get_start_point();
                Vector dir1 = sp1 - vtx_fit_pt(vtx);
                Vector dir2 = shower_cal_dir_3vector(*shower1, sp1, 15*units::cm);
                double ang  = dir1.angle(dir2) / M_PI * 180.0;
                if (ang < 30 && flag_single_shower &&
                    Eshower < 250*units::MeV && E1 > 60*units::MeV)
                    flag_bad1 = true;

                ti.sig_1_v_angle.push_back(ang);
                ti.sig_1_v_flag_single_shower.push_back(flag_single_shower);
                ti.sig_1_v_energy.push_back(Eshower / units::MeV);
                ti.sig_1_v_energy_1.push_back(E1 / units::MeV);
                ti.sig_1_v_flag.push_back(!flag_bad1);
                if (flag_bad1) flag_bad1_save = true;
            }
        }
    }

    // ------------------------------------------------------------------
    // flag_bad2: far-end segment has MIP-like dQ/dx while start is clean
    // ------------------------------------------------------------------
    // Pre-fill shower internal sets
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    // Find the shower vertex in the same cluster as vtx that maximises
    // the projection onto the shower direction (farthest forward vertex)
    double    max_dis = 0;
    VertexPtr max_vtx = nullptr;
    for (VertexPtr tmp_vtx : shower_vtxs) {
        if (tmp_vtx->cluster() != vtx->cluster()) continue;
        Vector dir1 = vtx_fit_pt(tmp_vtx) - vertex_point;
        double dis  = dir1.dot(shower_dir);
        if (dis > max_dis) { max_dis = dis; max_vtx = tmp_vtx; }
    }

    // Find the most anti-aligned shower-internal segment at max_vtx
    SegmentPtr max_sg = nullptr;
    double     max_angle_sg = 0;
    if (max_vtx && max_vtx->descriptor_valid()) {
        auto vd = max_vtx->get_descriptor();
        for (auto [eit, eend] = boost::out_edges(vd, ctx.graph); eit != eend; ++eit) {
            SegmentPtr sg1 = ctx.graph[*eit].segment;
            if (!sg1 || !shower_segs.count(sg1)) continue;
            Point max_vtx_pt = vtx_fit_pt(max_vtx);
            Vector dir1 = segment_cal_dir_3vector(sg1, max_vtx_pt, 15*units::cm);
            double ang = dir1.angle(shower_dir) / M_PI * 180.0;
            if (ang > max_angle_sg) { max_angle_sg = ang; max_sg = sg1; }
        }
    }

    if (max_vtx && max_sg) {
        // dQ/dx at far end of max_sg (from max_vtx end, normalized)
        const auto& max_fits = max_sg->fits();
        int n_max = (int)max_fits.size();
        bool max_front_near =
            !max_fits.empty() &&
            (ray_length(Ray{vtx_fit_pt(max_vtx), max_fits.front().point}) <=
             ray_length(Ray{vtx_fit_pt(max_vtx), max_fits.back().point}));
        double medium_dQ_dx = max_front_near
            ? segment_median_dQ_dx(max_sg, 0,       std::min(6, n_max-1)) / (43e3/units::cm)
            : segment_median_dQ_dx(max_sg, std::max(0, n_max-7), n_max-1) / (43e3/units::cm);

        // dQ/dx at shower start: max of first 3 per-fit-point normalized values
        double start_dQ_dx = 0;
        const auto& sg_fits = sg->fits();
        bool sg_front_near =
            !sg_fits.empty() &&
            (ray_length(Ray{vtx_fit_pt(vtx), sg_fits.front().point}) <=
             ray_length(Ray{vtx_fit_pt(vtx), sg_fits.back().point}));
        int ncount = 0;
        auto accumulate_start = [&](auto begin_it, auto end_it) {
            for (auto it = begin_it; it != end_it; ++it) {
                if (ncount > 2) break;
                double dqdx = it->dQ / (it->dx + 1e-9) / (43e3 / units::cm);
                if (dqdx > start_dQ_dx) start_dQ_dx = dqdx;
                ++ncount;
            }
        };
        if (sg_front_near) accumulate_start(sg_fits.begin(), sg_fits.end());
        else               accumulate_start(sg_fits.rbegin(), sg_fits.rend());

        // Cuts
        if ((Eshower < 250*units::MeV ||
             (Eshower < 500*units::MeV  && shower_angle > 120) ||
             (Eshower >= 500*units::MeV && shower_angle > 150)) && flag_single_shower) {
            // 7020_1108_55428
            if ((medium_dQ_dx > 1.6 && shower_angle > 60 && start_dQ_dx < 3.6) ||
                (medium_dQ_dx > 1.6 && shower_angle < 60 && start_dQ_dx < 2.5))
                flag_bad2 = true;
        }
        if (Eshower < 800*units::MeV && flag_single_shower) {
            if ((medium_dQ_dx > 2.0 && shower_angle > 60) ||
                (medium_dQ_dx > 2.0 && start_dQ_dx < 2.5))
                flag_bad2 = true;
        }

        ti.sig_2_v_energy.push_back(Eshower / units::MeV);
        ti.sig_2_v_shower_angle.push_back(shower_angle);
        ti.sig_2_v_flag_single_shower.push_back(flag_single_shower);
        ti.sig_2_v_medium_dQ_dx.push_back(medium_dQ_dx);
        ti.sig_2_v_start_dQ_dx.push_back(start_dQ_dx);
        ti.sig_2_v_flag.push_back(!flag_bad2);
    }

    bool flag_bad = flag_bad1_save || flag_bad2;
    ti.sig_flag = !flag_bad;
    return flag_bad;
}

// ===========================================================================
// gap_identification
//
// Checks for signal gaps in the shower stem near the interaction vertex.
// For each consecutive pair of fit points within 2.4 cm of the vertex,
// generates 3 sub-points and queries the charge point cloud for hits in
// each wire plane (U, V, W).  A sub-point is "bad" if fewer than 3 planes
// provide signal, dead-channel coverage, or prolongation exemption.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::gap_identification()
// Fills: ti.gap_*
//
// Translation notes:
//   main_cluster->check_direction(dir)        → inlined using grouping->wire_angles
//   ct_point_cloud->get_closest_points(p,r,pl) → grouping->get_closest_points(p,r,apa,face,pl)
//   ct_point_cloud->get_closest_dead_chs(p,pl,0)→grouping->get_closest_dead_chs(p,1,apa,face,pl)
//   apa/face obtained from ctx.dv->contained_by(test_p)
// ===========================================================================
static std::pair<bool, int> gap_identification(NuEContext& ctx,
                                               VertexPtr vertex, SegmentPtr sg,
                                               bool flag_single_shower,
                                               int num_valid_tracks,
                                               double E_shower,
                                               TaggerInfo& ti) {
    bool flag_gap = false;

    const auto& sg_fits = sg->fits();
    int n_fits = (int)sg_fits.size();
    bool vertex_at_front = n_fits > 0 &&
        (ray_length(Ray{vtx_fit_pt(vertex), sg_fits.front().point}) <=
         ray_length(Ray{vtx_fit_pt(vertex), sg_fits.back().point}));
    Point vertex_point = vertex_at_front ? sg_fits.front().point : sg_fits.back().point;

    // ------------------------------------------------------------------
    // Compute check_direction flags (inline of Graphs::check_direction)
    // ------------------------------------------------------------------
    bool flag_prolong_u = false, flag_prolong_v = false, flag_prolong_w = false;
    bool flag_parallel  = false;
    auto grouping = ctx.main_cluster ? ctx.main_cluster->grouping() : nullptr;
    if (grouping) {
        // Find fit point closest to 3cm from vertex (for direction estimate)
        Point closest_p = vertex_point;
        double min_d = 1e9;
        for (const auto& fit : sg_fits) {
            double d = std::fabs(ray_length(Ray{vertex_point, fit.point}) - 3*units::cm);
            if (d < min_d) { min_d = d; closest_p = fit.point; }
        }
        Vector dir_cd = closest_p - vertex_point;

        auto [angle_u, angle_v, angle_w] = grouping->wire_angles(ctx.apa, ctx.face);
        int drift_dirx = grouping->get_drift_dir().at(ctx.apa).at(ctx.face);
        Vector drift_abs(std::fabs(drift_dirx), 0, 0);
        Vector U_dir(0, std::cos(angle_u), std::sin(angle_u));
        Vector V_dir(0, std::cos(angle_v), std::sin(angle_v));
        Vector W_dir(0, std::cos(angle_w), std::sin(angle_w));
        Vector tempV1(0, dir_cd.y(), dir_cd.z());
        const double cut1 = 12.5 / 180.0 * M_PI;
        const double cut2 = 10.0 / 180.0 * M_PI;

        auto check_prolong = [&](const Vector& wire_dir) -> bool {
            double a  = tempV1.angle(wire_dir);
            double yz = std::sqrt(dir_cd.y()*dir_cd.y() + dir_cd.z()*dir_cd.z());
            Vector v5(std::fabs(dir_cd.x()), yz * std::sin(a), 0);
            return v5.angle(drift_abs) < cut1;
        };
        flag_prolong_u = check_prolong(U_dir);
        flag_prolong_v = check_prolong(V_dir);
        flag_prolong_w = check_prolong(W_dir);
        flag_parallel  = std::fabs(dir_cd.angle(drift_abs) - M_PI/2.0) < cut2;
    }

    // ------------------------------------------------------------------
    // Count bad sub-points along stem (within 2.4 cm of vertex)
    // ------------------------------------------------------------------
    int n_points = 0, n_bad = 0;

    auto query_sub_points = [&](int i_start, int i_end, int i_step) {
        for (int i = i_start; i != i_end; i += i_step) {
            int i_next = i + i_step;
            if (i_next < 0 || i_next >= n_fits) break;
            for (int j = 0; j < 3; ++j) {
                Point test_p = sg_fits[i].point +
                    (sg_fits[i_next].point - sg_fits[i].point) * (j / 3.0);

                int num_connect = 0, num_bad_ch = 0, num_spec = 0;
                if (grouping) {
                    auto wpid = ctx.dv->contained_by(test_p);
                    int q_apa = wpid.apa(), q_face = wpid.face();
                    if (q_apa >= 0) {
                        // U plane
                        if (!grouping->get_closest_points(test_p, 0.2*units::cm, q_apa, q_face, 0).empty())
                            ++num_connect;
                        else if (grouping->get_closest_dead_chs(test_p, 1, q_apa, q_face, 0))
                            ++num_bad_ch;
                        else if (flag_prolong_u) ++num_spec;
                        // V plane
                        if (!grouping->get_closest_points(test_p, 0.2*units::cm, q_apa, q_face, 1).empty())
                            ++num_connect;
                        else if (grouping->get_closest_dead_chs(test_p, 1, q_apa, q_face, 1))
                            ++num_bad_ch;
                        else if (flag_prolong_v) ++num_spec;
                        // W plane
                        if (!grouping->get_closest_points(test_p, 0.2*units::cm, q_apa, q_face, 2).empty())
                            ++num_connect;
                        else if (grouping->get_closest_dead_chs(test_p, 1, q_apa, q_face, 2))
                            ++num_bad_ch;
                        else if (flag_prolong_w) ++num_spec;
                    }
                }
                if (num_connect + num_bad_ch + num_spec == 3) {
                    if (n_bad == n_points && n_bad <= 5) n_bad = 0;
                } else {
                    ++n_bad;
                }
                ++n_points;
            }
            // Break once we're beyond 2.4cm from vertex
            Point far_pt = sg_fits[i_next].point;
            if (ray_length(Ray{vertex_point, far_pt}) > 2.4*units::cm) break;
        }
    };

    if (vertex_at_front)
        query_sub_points(0,       n_fits - 1, +1);
    else
        query_sub_points(n_fits-1, 0,          -1);

    // ------------------------------------------------------------------
    // Apply gap cuts (mirrors prototype exactly)
    // ------------------------------------------------------------------
    if (E_shower > 900*units::MeV) {
        if (!flag_single_shower && !flag_parallel) {
            if (E_shower > 1200*units::MeV) { if (n_bad > 2./3 * n_points) flag_gap = true; }
            else                             { if (n_bad > 1./3 * n_points) flag_gap = true; }
        }
        if (flag_parallel && !flag_single_shower) {
            if (n_bad > 0.5 * n_points) flag_gap = true;
        }
    } else if (E_shower > 150*units::MeV) {
        if (!flag_single_shower) {
            if (flag_parallel) { if (n_bad > 4) flag_gap = true; }
            else               { if (n_bad > 1) flag_gap = true; }
        } else {
            if (n_bad > 2) flag_gap = true;
        }
    } else {
        if (!flag_single_shower) {
            if (flag_parallel) { if (n_bad > 3) flag_gap = true; }
            else               { if (n_bad > 1) flag_gap = true; }
        } else {
            if (n_bad > 2) flag_gap = true;
        }
    }
    // 7021_521_26090
    if (n_bad >= 6 && E_shower < 1000*units::MeV) flag_gap = true;
    if (E_shower <= 900*units::MeV && n_bad > 1)  flag_gap = true;

    ti.gap_flag              = !flag_gap;
    ti.gap_flag_prolong_u    = flag_prolong_u;
    ti.gap_flag_prolong_v    = flag_prolong_v;
    ti.gap_flag_prolong_w    = flag_prolong_w;
    ti.gap_flag_parallel     = flag_parallel;
    ti.gap_n_points          = n_points;
    ti.gap_n_bad             = n_bad;
    ti.gap_energy            = E_shower / units::MeV;
    ti.gap_num_valid_tracks  = num_valid_tracks;
    ti.gap_flag_single_shower = flag_single_shower;
    ti.gap_filled            = 1;

    return {flag_gap, n_bad};
}

// ===========================================================================
// bad_reconstruction_3
//
// Checks two failure modes for "bad reconstruction" of the main cluster:
//
//   flag_bad1: the shower's main-cluster fraction is small and there are
//     other-cluster shower segments far from the main-cluster tip — the
//     shower likely crosses between clusters.
//
//   flag_bad2: the angular distribution of shower hits and vertices relative
//     to the shower direction is too narrow — consistent with a track rather
//     than an EM shower fan.
//
// Prototype: NeutrinoID_nue_tagger.h, WCPPID::NeutrinoID::bad_reconstruction_3()
// Fills: ti.br4_1_*, ti.br4_2_*, ti.br4_flag
// ===========================================================================
static bool bad_reconstruction_3(NuEContext& ctx,
                                  VertexPtr vertex, ShowerPtr shower,
                                  TaggerInfo& ti) {
    bool flag_bad1 = false, flag_bad2 = false;

    Vector drift_dir(1, 0, 0);
    double Eshower = (shower->get_kine_best() != 0)
                     ? shower->get_kine_best() : shower->get_kine_charge();

    double main_length  = vertex->cluster() ? shower->get_total_length(vertex->cluster()) : 0;
    double total_length = shower->get_total_length();

    // Pre-fill shower internal sets
    IndexedSegmentSet shower_segs;
    IndexedVertexSet  shower_vtxs;
    shower->fill_sets(shower_vtxs, shower_segs, /*flag_exclude_start_segment=*/false);

    // ------------------------------------------------------------------
    // br4_1: main-cluster fraction vs distance to closest off-cluster seg
    // ------------------------------------------------------------------
    // Find the shower vertex in the main cluster farthest from vertex
    double max_dis_v = 0;
    Point  max_p     = vtx_fit_pt(vertex);
    for (VertexPtr vtx1 : shower_vtxs) {
        if (vtx1->cluster() != vertex->cluster()) continue;
        double d = ray_length(Ray{vtx_fit_pt(vertex), vtx_fit_pt(vtx1)});
        if (d > max_dis_v) { max_dis_v = d; max_p = vtx_fit_pt(vtx1); }
    }

    // Find the off-cluster shower segment closest to max_p (skipping short segs)
    double min_dis  = 1e9;
    for (SegmentPtr sg1 : shower_segs) {
        if (sg1->cluster() == vertex->cluster()) continue;
        if (segment_track_length(sg1) < 6*units::cm) continue;
        double d = segment_get_closest_point(sg1, max_p).first;
        if (d < min_dis) min_dis = d;
    }

    // Accumulate close off-cluster segments (within min_dis)
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

    // 7006_489_24469: use refined min_dis if cluster has significant content
    SegmentPtr start_sg = shower->start_segment();
    if (acc_close_length > 10*units::cm ||
        (num_close >= 3 && acc_close_length > 4.5*units::cm) ||
        (start_sg && start_sg->flags_any(SegmentFlags::kAvoidMuonCheck)))
        min_dis = min_dis1;

    if (min_dis < 1e7) {
        if (main_length < 0.40 * total_length && min_dis > 40*units::cm) flag_bad1 = true;
        if (main_length < 0.25 * total_length && min_dis > 33*units::cm) flag_bad1 = true;
        if (main_length < 0.16 * total_length && min_dis > 23*units::cm) flag_bad1 = true;
        if (main_length < 0.10 * total_length && min_dis > 18*units::cm) flag_bad1 = true;
        if (main_length < 0.05 * total_length && min_dis >  8*units::cm) flag_bad1 = true;
        if (main_length < 8*units::cm && main_length < 0.1*total_length &&
            (min_dis > 8*units::cm && Eshower < 300*units::MeV || min_dis > 14*units::cm))
            flag_bad1 = true;

        // Exceptions
        // 7003_1226_61332
        if (flag_bad1 && start_sg && start_sg->flags_any(SegmentFlags::kAvoidMuonCheck) &&
            main_length > 12*units::cm && main_length > 0.1*total_length &&
            min_dis < 40*units::cm) flag_bad1 = false;
        // 7010_405_20296: start vertex with only 1 global edge
        VertexPtr start_vtx = shower->get_start_vertex_and_type().first;
        size_t n_vtx_segs = start_vtx && start_vtx->descriptor_valid()
                            ? boost::out_degree(start_vtx->get_descriptor(), ctx.graph) : 0;
        if (n_vtx_segs == 1 &&
            ((main_length > 20*units::cm && min_dis < 40*units::cm && main_length > 0.1*total_length) ||
             (main_length > 15*units::cm && min_dis < 32*units::cm && main_length > 0.15*total_length)))
            flag_bad1 = false;
        // 7049_1874_93741
        if (flag_bad1 && main_length > 30*units::cm && shower->get_num_main_segments() >= 4)
            flag_bad1 = false;
    }

    ti.br4_1_shower_main_length       = main_length / units::cm;
    ti.br4_1_shower_total_length      = total_length / units::cm;
    ti.br4_1_min_dis                  = min_dis / units::cm;
    ti.br4_1_energy                   = Eshower / units::MeV;
    ti.br4_1_flag_avoid_muon_check    = (start_sg && start_sg->flags_any(SegmentFlags::kAvoidMuonCheck));
    ti.br4_1_n_vtx_segs               = (start_sg ? [&](){
        VertexPtr sv = shower->get_start_vertex_and_type().first;
        return sv && sv->descriptor_valid() ? (int)boost::out_degree(sv->get_descriptor(), ctx.graph) : 0;
    }() : 0);
    ti.br4_1_n_main_segs              = shower->get_num_main_segments();
    ti.br4_1_flag                     = !flag_bad1;

    // ------------------------------------------------------------------
    // br4_2: angular distribution of shower hits relative to shower dir
    // ------------------------------------------------------------------
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
            for (size_t i = 1; i + 1 < fits1.size(); ++i) {
                count_dir(fits1[i].point - vp, is_main);
            }
            // Endpoint vertices
            auto [v1, v2] = find_vertices(ctx.graph, sg1);
            for (VertexPtr ep : {v1, v2}) {
                if (ep) count_dir(vtx_fit_pt(ep) - vp, is_main);
            }
        }

        if (ncount_45 < 0.7*ncount || ncount_25 < 0.6*ncount ||
            (ncount_25 < 0.8*ncount && ncount_15 < 0.3*ncount) ||
            (ncount_15 < 0.35*ncount && ncount_25 > 0.9*ncount && Eshower < 1000*units::MeV))
            flag_bad2 = true;

        double iso_angle  = std::fabs(dir.angle(drift_dir)    / M_PI * 180.0 - 90.0);
        double iso_angle1 = std::fabs(dir_sg.angle(drift_dir) / M_PI * 180.0 - 90.0);
        double sg_dir_angle = dir_sg.angle(dir) / M_PI * 180.0;

        bool cut_b = (ncount1_15 < 0.35*ncount1 && iso_angle > 15 &&
                      (ncount1_25 < 0.95*ncount1 && Eshower < 1000*units::MeV || Eshower >= 1000*units::MeV)) ||
                     (ncount1_15 < 0.2*ncount1 && ncount1_25 < 0.45*ncount1 && Eshower < 600*units::MeV) ||
                     (sg_dir_angle > 25 && std::max(iso_angle, iso_angle1) > 8 &&
                      (ncount1_15 < 0.8*ncount1 && Eshower < 1000*units::MeV || Eshower >= 1000*units::MeV)) ||
                     (sg_dir_angle > 20 && std::max(iso_angle, iso_angle1) > 5 &&
                      ncount1_15 < 0.5*ncount1);
        if (cut_b) flag_bad2 = true;

        if (ncount > 0) {
            ti.br4_2_ratio_45  = ncount_45  / (ncount + 1e-9);
            ti.br4_2_ratio_35  = ncount_35  / (ncount + 1e-9);
            ti.br4_2_ratio_25  = ncount_25  / (ncount + 1e-9);
            ti.br4_2_ratio_15  = ncount_15  / (ncount + 1e-9);
        } else {
            ti.br4_2_ratio_45 = ti.br4_2_ratio_35 = ti.br4_2_ratio_25 = ti.br4_2_ratio_15 = 1;
        }
        ti.br4_2_energy = Eshower / units::MeV;
        if (ncount1 > 0) {
            ti.br4_2_ratio1_45 = ncount1_45 / (ncount1 + 1e-9);
            ti.br4_2_ratio1_35 = ncount1_35 / (ncount1 + 1e-9);
            ti.br4_2_ratio1_25 = ncount1_25 / (ncount1 + 1e-9);
            ti.br4_2_ratio1_15 = ncount1_15 / (ncount1 + 1e-9);
        } else {
            ti.br4_2_ratio1_45 = ti.br4_2_ratio1_35 = ti.br4_2_ratio1_25 = ti.br4_2_ratio1_15 = 1;
        }
        ti.br4_2_iso_angle  = iso_angle;
        ti.br4_2_iso_angle1 = iso_angle1;
        ti.br4_2_angle      = sg_dir_angle;
        ti.br4_2_flag       = !flag_bad2;
    }

    bool flag_bad = flag_bad1 || flag_bad2;
    ti.br4_flag = !flag_bad;
    return flag_bad;
}
