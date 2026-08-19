#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/GraphTools.h"  // doc pr/93 r4: mir() for the orphan-track pass
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;
using namespace WireCell;
using WireCell::GraphTools::mir;  // doc pr/93 r4

// init_tagger_info: reset a TaggerInfo struct to its default values.
// In the toolkit we rely on C++ default-member-initializers on the struct
// itself (see NeutrinoTaggerInfo.h), so value-initializing the struct is
// sufficient — no 1200-line assignment list needed.
void PatternAlgorithms::init_tagger_info(TaggerInfo& ti)
{
    ti = TaggerInfo{};
}

// fill_kine_tree: reconstruct per-particle and total neutrino kinetic energy
// starting from main_vertex, traversing the PR graph, and collecting all
// connected showers and track segments.
//
// Translation notes from prototype (WCPPID::NeutrinoID::fill_kine_tree):
//   map_vertex_segments[vtx]              -> boost::out_edges(vtx->get_descriptor(), graph) + graph[*ei].segment
//   find_other_vertex(seg, vtx)           -> find_other_vertex(graph, seg, vtx)
//   shower->get_start_segment()           -> shower->start_segment()
//   shower->get_start_segment()->get_particle_type() -> shower->get_particle_type()
//   shower->get_start_segment()->get_particle_mass() -> shower->start_segment()->particle_info()->mass()
//   shower->get_start_segment()->get_length()        -> segment_track_length(shower->start_segment())
//   seg->get_particle_type()              -> seg->particle_info()->pdg()
//   seg->get_kine_best()                  -> seg->particle_info()->kinetic_energy()
//   seg->get_particle_mass()              -> seg->particle_info()->mass()
//   cal_kine_charge(seg)                  -> cal_kine_charge(seg, graph, track_fitter, dv)
//   seg->cal_kine_dQdx()                  -> segment_cal_kine_dQdx(seg, recomb_model)
//   seg->cal_kine_range()                 -> cal_kine_range(segment_track_length(seg), pdg, particle_data)
//   shower->get_start_vertex()            -> shower->get_start_vertex_and_type()
//   kine_pio_*                            -> pio_kine struct fields
//   SCE correction                        -> geom_helper->get_corrected_point(...) (or raw point if null)
KineInfo PatternAlgorithms::fill_kine_tree(
    VertexPtr main_vertex,
    IndexedShowerSet& showers,
    const Pi0KineFeatures& pio_kine,
    Graph& graph,
    TrackFitting& track_fitter,
    IDetectorVolumes::pointer dv,
    WireCell::IClusGeomHelper::pointer geom_helper,
    const Clus::ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    const IndexedShowerSet& pi0_showers,
    std::set<int>* dropped_satellites)
{
    KineInfo ktree{};

    // -------------------------------------------------------------------------
    // Neutrino vertex position with optional SCE correction
    // -------------------------------------------------------------------------
    Point nu_vtx = main_vertex->fit().point;

    if (geom_helper) {
        WirePlaneId wpid = dv->contained_by(nu_vtx);
        int apa  = wpid.apa();
        int face = wpid.face();
        Point corr = geom_helper->get_corrected_point(nu_vtx, IClusGeomHelper::SCE, apa, face);
        ktree.kine_nu_x_corr = static_cast<float>(corr.x() / units::cm);
        ktree.kine_nu_y_corr = static_cast<float>(corr.y() / units::cm);
        ktree.kine_nu_z_corr = static_cast<float>(corr.z() / units::cm);
    }
    else {
        // TODO: SCE correction requires a valid geom_helper; using raw vertex position for now.
        // doc pr/35 §10.5 (F4 = P5): the prototype SCE-corrects unconditionally
        // (kine.h:3-9); without a geom_helper the _corr branches are raw and
        // nothing downstream can tell.  Owner decision 2026-08-04: keep the raw
        // vertex on SBND, but say so at runtime.
        SPDLOG_LOGGER_WARN(s_log, "fill_kine_tree: no geom_helper -- kine_nu_*_corr are "
                                  "the RAW fitted vertex despite the _corr name");
        ktree.kine_nu_x_corr = static_cast<float>(nu_vtx.x() / units::cm);
        ktree.kine_nu_y_corr = static_cast<float>(nu_vtx.y() / units::cm);
        ktree.kine_nu_z_corr = static_cast<float>(nu_vtx.z() / units::cm);
    }

    // -------------------------------------------------------------------------
    // Build a map from a shower's start-segment to the shower itself, and
    // collect all vertices/segments already owned by showers.
    // -------------------------------------------------------------------------
    const double ave_binding_energy = 8.6 * units::MeV;

    IndexedVertexSet  used_vertices;
    IndexedSegmentSet used_segments;
    IndexedShowerSet  used_showers;

    // Mark all shower-internal vertices and segments as used.
    for (const ShowerPtr& shower : showers) {
        shower->fill_sets(used_vertices, used_segments, /*flag_exclude_start_segment=*/false);
    }

    // Map from start-segment -> shower, for fast lookup during graph traversal.
    std::map<SegmentPtr, ShowerPtr, SegmentIndexCmp> map_sg_shower;
    for (const ShowerPtr& shower : showers) {
        SegmentPtr start_sg = shower->start_segment();
        map_sg_shower[start_sg] = shower;
    }

    // -------------------------------------------------------------------------
    // Helper: the shower PDG this tree publishes.
    // doc pr/35 §10.2 (F1 = P1 + P8): the prototype reads the LIVE start-segment
    // PDG at every fill_kine_tree site (kine.h:53 :67 :175 :187) where the
    // toolkit reads Shower's cached particle_type, whose refresh path is
    // incomplete (P1).  kine_shower_pdg_live selects the live read; the
    // cached-vs-live PAIR is logged unconditionally (pr/32 F3 counter
    // precedent) so the mismatch population is measurable without an A/B --
    // the 11-vs-211 and cache-stuck-at-0 classes must stay distinguishable.
    // -------------------------------------------------------------------------
    auto shower_kine_pdg = [&](const ShowerPtr& shower) -> int {
        const int  cached   = shower->get_particle_type();
        const bool has_live = shower->start_segment() && shower->start_segment()->has_particle_info();
        const int  live     = has_live ? shower->start_segment()->particle_info()->pdg() : cached;
        if (cached != live) {
            SPDLOG_LOGGER_INFO(s_log, "fill_kine_tree: kine_shower_pdg cached={} live={} shower_id={}",
                               cached, live, shower->get_shower_id());
        }
        return (m_kine_charge.shower_pdg_live && has_live) ? live : cached;
    };

    // -------------------------------------------------------------------------
    // Helper: push one shower's kinematics into ktree vectors.
    // kine_energy_included is pushed by the caller (value differs by context).
    // -------------------------------------------------------------------------
    auto push_shower_kine = [&](const ShowerPtr& shower) {
        double kine_best   = shower->get_kine_best();
        double kine_charge = shower->get_kine_charge();
        double kine_range  = shower->get_kine_range();
        const int pdg      = shower_kine_pdg(shower);

        ktree.kine_energy_particle.push_back(static_cast<float>(kine_best / units::MeV));
        ktree.kine_particle_type.push_back(pdg);

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2); // charge
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1); // range
        else
            ktree.kine_energy_info.push_back(0); // dQdx

        // Add rest-mass correction for non-electrons/positrons
        if (pdg != 11) {
            SegmentPtr start_sg = shower->start_segment();
            if (start_sg && start_sg->particle_info()) {
                ktree.kine_reco_add_energy += static_cast<float>(
                    start_sg->particle_info()->mass() / units::MeV);
            }
        }
    };

    // -------------------------------------------------------------------------
    // Helper: push one track segment's kinematics into ktree vectors.
    // Returns the segment's PDG code.
    // -------------------------------------------------------------------------
    auto push_segment_kine = [&](SegmentPtr seg, int include_flag) -> int {
        int    pdg        = 0;
        double mass       = 0;
        double kine_best  = 0;

        if (seg->particle_info()) {
            pdg       = seg->particle_info()->pdg();
            mass      = seg->particle_info()->mass();
            kine_best = seg->particle_info()->kinetic_energy();
        }
        double kine_charge = cal_kine_charge(seg, graph, track_fitter, dv);
        double kine_range  = cal_kine_range(segment_track_length(seg), pdg, particle_data);

        ktree.kine_energy_particle.push_back(static_cast<float>(kine_best / units::MeV));
        ktree.kine_particle_type.push_back(pdg);

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2);
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1);
        else
            ktree.kine_energy_info.push_back(0);

        ktree.kine_energy_included.push_back(include_flag);

        if (pdg == 2212) { // proton: add binding energy
            ktree.kine_reco_add_energy += static_cast<float>(ave_binding_energy / units::MeV);
        }
        else if (pdg != 11) { // not electron: add rest mass
            ktree.kine_reco_add_energy += static_cast<float>(mass / units::MeV);
        }
        return pdg;
    };

    // -------------------------------------------------------------------------
    // First pass: segments directly connected to main_vertex
    // -------------------------------------------------------------------------
    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;

    const auto ei_begin_edges = sorted_out_edges(main_vertex->get_descriptor(), graph);
    for (auto ei : ei_begin_edges) {
        SegmentPtr seg = graph[ei].segment;

        auto it = map_sg_shower.find(seg);
        if (it != map_sg_shower.end()) {
            // This segment is a shower start-segment.
            push_shower_kine(it->second);
            ktree.kine_energy_included.push_back(1);
            used_showers.insert(it->second);
        }
        else {
            // Track segment.
            used_segments.insert(seg);
            VertexPtr other_vtx = find_other_vertex(graph, seg, main_vertex);
            segments_to_be_examined.emplace_back(other_vtx, seg);
            push_segment_kine(seg, 1);
        }
    }
    used_vertices.insert(main_vertex);

    // -------------------------------------------------------------------------
    // BFS traversal of the remaining track graph
    // -------------------------------------------------------------------------
    while (!segments_to_be_examined.empty()) {
        std::vector<std::pair<VertexPtr, SegmentPtr>> temp_segments;

        for (auto& [curr_vtx, prev_sg] : segments_to_be_examined) {
            if (used_vertices.count(curr_vtx)) continue;

            bool flag_reduce = false;
            int  prev_pdg = 0;
            if (prev_sg->particle_info()) prev_pdg = prev_sg->particle_info()->pdg();

            const auto ei2_begin_edges = sorted_out_edges(curr_vtx->get_descriptor(), graph);
            for (auto ei2 : ei2_begin_edges) {
                SegmentPtr curr_sg = graph[ei2].segment;
                if (curr_sg == prev_sg) continue;

                int curr_pdg = 0;
                if (curr_sg->particle_info()) curr_pdg = curr_sg->particle_info()->pdg();

                // Detect particle continuation (same type, or muon<->pion flip)
                if (curr_pdg == prev_pdg ||
                    (prev_pdg == 211 && curr_pdg == 13) ||
                    (prev_pdg == 13  && curr_pdg == 211))
                    flag_reduce = true;

                auto it2 = map_sg_shower.find(curr_sg);
                if (it2 == map_sg_shower.end()) {
                    // Track segment
                    if (used_segments.count(curr_sg)) continue;
                    used_segments.insert(curr_sg);

                    push_segment_kine(curr_sg, 1);

                    VertexPtr other_vtx = find_other_vertex(graph, curr_sg, curr_vtx);
                    if (!used_vertices.count(other_vtx))
                        temp_segments.emplace_back(other_vtx, curr_sg);
                }
                else {
                    // Shower
                    const ShowerPtr& shower = it2->second;
                    if (!used_showers.count(shower)) {
                        push_shower_kine(shower);
                        ktree.kine_energy_included.push_back(1);
                        used_showers.insert(shower);
                    }
                }
            }
            used_vertices.insert(curr_vtx);

            // If we detected a particle continuation, undo the rest-mass/binding-energy
            // added for prev_sg (it was already counted earlier in the chain).
            if (flag_reduce && prev_sg->particle_info()) {
                if (prev_pdg == 2212) {
                    ktree.kine_reco_add_energy -= static_cast<float>(ave_binding_energy / units::MeV);
                }
                else if (prev_pdg != 11) {
                    ktree.kine_reco_add_energy -= static_cast<float>(
                        prev_sg->particle_info()->mass() / units::MeV);
                }
            }
        }

        segments_to_be_examined = std::move(temp_segments);
    }

    // -------------------------------------------------------------------------
    // doc sbnd_xin/docs/pr/92 -- stray-satellite drop decision.
    //
    // The leftover pass below admits every BFS-unreached conn-2/3 shower
    // with no direction/distance check (the prototype has the identical
    // hole, NeutrinoID_kine.h:209-255), so overclustered cosmics and
    // second neutrinos are summed into kine_reco_Enu.  Decide here, on the
    // BFS-unreached candidates only: anything the BFS consumed is
    // genuinely graph-connected to the main vertex and never a candidate.
    // See NeutrinoPatternBase.h (pr/92 block) for the arms and the target
    // events.  WCT_KINE_SAT_PROBE prints per-candidate metrics for the
    // full population and NEVER drops (threshold tuning on
    // at-decision-time numbers; the probe line reports the would-verdict).
    // -------------------------------------------------------------------------
    const bool sat_probe = (std::getenv("WCT_KINE_SAT_PROBE") != nullptr);
    IndexedShowerSet sat_drop_set;
    if (m_kine_drop_stray_satellites || sat_probe) {
        const auto* main_cluster = main_vertex->cluster();
        const Point mv_pt = main_vertex->fit().point;
        auto angle_deg = [](const Vector& a, const Vector& b) -> double {
            const double ma = a.magnitude(), mb = b.magnitude();
            if (ma <= 0 || mb <= 0) return 0.0;
            double c = a.dot(b) / (ma * mb);
            c = std::max(-1.0, std::min(1.0, c));
            return std::acos(c) / M_PI * 180.0;
        };
        int n_sat_cand = 0;
        for (const ShowerPtr& shower : showers) {
            if (used_showers.count(shower)) continue;
            auto [svtx, conn] = shower->get_start_vertex_and_type();
            if (conn != 2 && conn != 3) continue;
            SegmentPtr start_sg = shower->start_segment();
            const auto* start_cl = start_sg ? start_sg->cluster() : nullptr;
            if (!start_cl || !main_cluster ||
                start_cl->get_cluster_id() == main_cluster->get_cluster_id()) continue;
            ++n_sat_cand;

            const Point sp = shower->get_start_point();
            const Point sv_pt = svtx ? (svtx->fit().valid() ? svtx->fit().point
                                                            : svtx->wcpt().point)
                                     : mv_pt;
            // Fresh axis from the start point into the shower body: the
            // STORED init_dir for conn 2/3 is exactly the vertex->start
            // chord (PRShower.cxx) and would always read 0 deg here.
            const Vector axis = shower_cal_dir_3vector(*shower, sp, m_kine_sat_axis_dis_cut);
            const bool in_main = svtx && (svtx == main_vertex ||
                (svtx->cluster() && svtx->cluster()->get_cluster_id() ==
                                        main_cluster->get_cluster_id()));
            const double d_sv   = (sp - sv_pt).magnitude();
            const double ang_sv = (d_sv > 0) ? angle_deg(axis, sp - sv_pt) : 0.0;
            const double d_mv   = (sp - mv_pt).magnitude();
            const double ang_mv = (d_mv > 0) ? angle_deg(axis, sp - mv_pt) : 0.0;
            const bool is_pi0   = pi0_showers.count(shower);
            const bool straight_cont =
                shower_start_is_track_continuation(graph, *shower, m_kine_sat_cont_kink);
            // pr/92 round 2 (owner retune): topology split.  A TRACK-like
            // satellite (straight-long start segment with little branching,
            // or a collinear continuation of an out-of-shower track) with a
            // bad direction is very likely overclustering; an EM-shower-like
            // satellite (branched, stubby trunk) is usually a genuinely
            // detached legit shower (NCpi0-like), so it is dropped only when
            // it is FAR from the main vertex AND direction-inconsistent --
            // the second-neutrino signature (389538: 169-250 cm), never the
            // nearby-fragment one (259542: 18-75 cm).  The EM angle is
            // folded (sign-insensitive): EM axis signs flip easily, and an
            // anti-aligned axis is still collinear with the vertex line
            // (52672's 453 MeV at 171 deg is a KEEP).
            const bool track_like = straight_cont ||
                (shower->get_num_segments() <= m_kine_sat_track_max_nseg &&
                 segment_is_straight_long_track(start_sg));
            const double ang_mv_fold = std::min(ang_mv, 180.0 - ang_mv);

            const char* verdict = "keep";
            do {
                if (is_pi0) break;
                if (shower->get_kine_best() <= m_kine_sat_min_energy) break;
                if (axis.magnitude() == 0) break;   // no axis: fail-safe keep
                if (d_sv < m_kine_sat_prox_max && in_main) break;
                if (track_like) {
                    if (ang_sv > m_kine_sat_angle_bad)          { verdict = "drop:A"; break; }
                    if ((d_sv > m_kine_sat_far_dis || !in_main) &&
                        ang_mv >= m_kine_sat_angle_main)        { verdict = "drop:B"; break; }
                    if (straight_cont)                          { verdict = "drop:C"; break; }
                }
                else if (d_mv > m_kine_sat_em_far_dis &&
                         ang_mv_fold >= m_kine_sat_angle_main)  { verdict = "drop:E"; break; }
            } while (false);

            if (sat_probe) {
                // One pre-built line per candidate (log lines tear mid-word
                // when interleaved -- write the TSV row atomically).
                std::ostringstream os;
                os << "KINE_SAT_PROBE\t" << shower->get_shower_id()
                   << '\t' << shower->get_kine_best() / units::MeV
                   << '\t' << conn
                   << '\t' << shower->get_particle_type()
                   << '\t' << shower->get_num_segments()
                   << '\t' << d_sv / units::cm
                   << '\t' << ang_sv
                   << '\t' << (in_main ? 1 : 0)
                   << '\t' << d_mv / units::cm
                   << '\t' << ang_mv
                   << '\t' << (straight_cont ? 1 : 0)
                   << '\t' << (is_pi0 ? 1 : 0)
                   << '\t' << (track_like ? 1 : 0)
                   << '\t' << verdict << '\n';
                std::cout << os.str() << std::flush;
            }
            else if (verdict[0] == 'd') {
                sat_drop_set.insert(shower);
                if (dropped_satellites) dropped_satellites->insert(shower->get_shower_id());
                SPDLOG_LOGGER_DEBUG(s_log,
                    "kine_sat DROP id={} arm={} E={:.1f} conn={} d_sv={:.1f} ang_sv={:.1f} in_main={} d_mv={:.1f} ang_mv={:.1f} cont={} track={}",
                    shower->get_shower_id(), verdict,
                    shower->get_kine_best() / units::MeV, conn,
                    d_sv / units::cm, ang_sv, in_main, d_mv / units::cm, ang_mv,
                    straight_cont, track_like);
            }
        }
        if (m_kine_drop_stray_satellites && !sat_probe) {
            SPDLOG_LOGGER_INFO(s_log, "kine_sat census: candidates={} dropped={}",
                               n_sat_cand, sat_drop_set.size());
        }
    }

    // -------------------------------------------------------------------------
    // Remaining showers not yet attached to the traversal above
    // (e.g. secondary showers with start vertex type <= 3)
    // -------------------------------------------------------------------------
    for (const ShowerPtr& shower : showers) {
        if (used_showers.count(shower)) continue;
        // doc pr/92: skip the whole iteration -- all four parallel kine
        // vectors and the proton binding-energy add stay consistent, and
        // the Enu sum below shrinks automatically.
        if (sat_drop_set.count(shower)) continue;

        auto [start_vtx, vtx_type] = shower->get_start_vertex_and_type();
        if (vtx_type > 3) continue;

        double kine_best   = shower->get_kine_best();
        double kine_charge = shower->get_kine_charge();
        double kine_range  = shower->get_kine_range();
        const int pdg      = shower_kine_pdg(shower);

        ktree.kine_energy_particle.push_back(static_cast<float>(kine_best / units::MeV));
        ktree.kine_particle_type.push_back(pdg);

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2);
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1);
        else
            ktree.kine_energy_info.push_back(0);

        ktree.kine_energy_included.push_back(vtx_type != 3 ? 1 : vtx_type);

        // Binding energy correction for proton showers with length > 5 cm
        if (pdg == 2212) {
            SegmentPtr start_sg = shower->start_segment();
            if (start_sg && segment_track_length(start_sg) > 5.0 * units::cm) {
                ktree.kine_reco_add_energy += static_cast<float>(ave_binding_energy / units::MeV);
            }
        }
        // electrons: no rest-mass correction; other non-electrons not present in remaining showers
        // since push_shower_kine would have already handled them in the BFS phase

        used_showers.insert(shower);
    }

    // -------------------------------------------------------------------------
    // doc sbnd_xin/docs/pr/93 round 4 (kine_count_orphan_tracks): count
    // confident straight-long main-cluster track segments that neither the
    // main-vertex BFS nor any shower claimed.  Such segments exist because
    // shower_cone_absorb_guard frees a graph-disconnected track from shower
    // membership (SBND 18255-315167: a 150.7cm pdg-2212 score-0.101 proton,
    // ~595 MeV KE, silently absent from kine_reco_Enu).  Shares
    // segment_orphan_confident_track with the PF-side pf_orphan_confident_
    // track knob so the two outputs describe the same particle set.
    // C++ default false => no pass => byte-identical.
    // -------------------------------------------------------------------------
    if (m_kine_count_orphan_tracks) {
        const auto* main_cl = main_vertex->cluster();
        // Deterministic candidate order: collect then sort by graph edge
        // index -- never pointer order.
        std::vector<std::pair<size_t, SegmentPtr>> orphan_cands;
        for (auto edesc : mir(boost::edges(graph))) {
            SegmentPtr seg = graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            if (used_segments.count(seg)) continue;
            if (map_sg_shower.count(seg)) continue;   // shower starts stay shower-owned
            const auto* cl = seg->cluster();
            if (!main_cl || !cl || cl->get_cluster_id() != main_cl->get_cluster_id()) continue;
            if (!segment_orphan_confident_track(seg, m_kine_orphan_track_min)) continue;
            orphan_cands.emplace_back(graph[edesc].index, seg);
        }
        std::sort(orphan_cands.begin(), orphan_cands.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        size_t n_orphan = 0;
        for (const auto& [eidx, seg] : orphan_cands) {
            used_segments.insert(seg);
            const int pdg = push_segment_kine(seg, 1);
            ++n_orphan;
            SPDLOG_LOGGER_INFO(s_log,
                "kine_count_orphan_tracks: COUNT seg idx={} cluster={} pdg={} ke_mev={:.2f} len_cm={:.1f}",
                eidx, seg->cluster()->get_cluster_id(), pdg,
                (seg->particle_info() ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0),
                segment_track_length(seg) / units::cm);
        }
        if (n_orphan) {
            SPDLOG_LOGGER_INFO(s_log, "kine_count_orphan_tracks: {} orphan track(s) added", n_orphan);
        }
    }

    // -------------------------------------------------------------------------
    // Total reconstructed neutrino energy
    // -------------------------------------------------------------------------
    ktree.kine_reco_Enu = 0;
    for (float e : ktree.kine_energy_particle) {
        ktree.kine_reco_Enu += e;
    }
    ktree.kine_reco_Enu += ktree.kine_reco_add_energy;

    // -------------------------------------------------------------------------
    // Pi0 kinematics (from pio_kine struct, angles converted to degrees)
    // -------------------------------------------------------------------------
    ktree.kine_pio_mass    = static_cast<float>(pio_kine.mass    / units::MeV);
    ktree.kine_pio_flag    = pio_kine.flag;
    ktree.kine_pio_vtx_dis = static_cast<float>(pio_kine.vtx_dis / units::cm);

    ktree.kine_pio_energy_1 = static_cast<float>(pio_kine.energy_1 / units::MeV);
    ktree.kine_pio_theta_1  = static_cast<float>(pio_kine.theta_1  / M_PI * 180.0);
    ktree.kine_pio_phi_1    = static_cast<float>(pio_kine.phi_1    / M_PI * 180.0);
    ktree.kine_pio_dis_1    = static_cast<float>(pio_kine.dis_1    / units::cm);

    ktree.kine_pio_energy_2 = static_cast<float>(pio_kine.energy_2 / units::MeV);
    ktree.kine_pio_theta_2  = static_cast<float>(pio_kine.theta_2  / M_PI * 180.0);
    ktree.kine_pio_phi_2    = static_cast<float>(pio_kine.phi_2    / M_PI * 180.0);
    ktree.kine_pio_dis_2    = static_cast<float>(pio_kine.dis_2    / units::cm);

    ktree.kine_pio_angle = static_cast<float>(pio_kine.angle / M_PI * 180.0);

    return ktree;
}
