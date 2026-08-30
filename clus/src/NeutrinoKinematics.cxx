#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/GraphTools.h"  // doc pr/93 r4: mir() for the orphan-track pass
#include <cmath>
#include <cstdlib>
#include <cstdio>   // doc pr/130 probe: std::fprintf
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

    // doc pr/101 (K2, kine_mass_rules): the paper's rest term.  mu/pi (and
    // any other non-nucleon, non-electron type) add the rest mass, nucleons
    // add the binding energy, electrons add nothing.  Returned in MeV.
    // Off => every add site keeps its legacy branch verbatim.
    auto rest_term_rules = [&](int pdg, double mass) -> double {
        const int a = std::abs(pdg);
        if (a == 11) return 0.0;
        if (a == 2212 || a == 2112) return ave_binding_energy / units::MeV;
        if (a == 13 || a == 211 || a == 321) return mass / units::MeV;
        if (a == 0 || a == 22 || a == 111) {
            SPDLOG_LOGGER_INFO(s_log, "kine_mass_rules: unhandled pdg={} -> no rest term", pdg);
            return 0.0;
        }
        return mass / units::MeV;
    };
    // Legacy value computed in parallel when the rules are on (census only).
    float add_energy_legacy = 0;
    int   n_2212_shower_graph = 0, n_leftover_nonem = 0, n_mainvtx_guard_skip = 0;

    IndexedVertexSet  used_vertices;
    IndexedSegmentSet used_segments;
    IndexedShowerSet  used_showers;

    // doc 85 sec 9 (owner request 2026-08-30): exactly the segments whose
    // kinetic energy reaches kine_reco_Enu -- recorded inside the two push
    // helpers so it cannot drift from what they actually summed.  NOT
    // used_segments, which is pre-seeded below with every shower's members
    // including showers this walk never counts.
    IndexedSegmentSet counted_segs;

    // Mark all shower-internal vertices and segments as used.
    for (const ShowerPtr& shower : showers) {
        shower->fill_sets(used_vertices, used_segments, /*flag_exclude_start_segment=*/false);
    }

    // doc sbnd_xin/docs/pr/128 (kine_count_conn4_near): conn-4 showers whose
    // material is the candidate's OWN -- closest approach to the main cluster
    // within m_kine_conn4_near_gap.  conn-4 means "cluster >80 cm from the
    // candidate" and dropping it from Enu is right for 490 of the 514 conn-4
    // showers measured (2815 of 3514 MeV, all >=50 cm away: the far-away
    // over-clustered activity the owner ruled must NOT be counted).  But
    // pr/74 conn3_unreachable, pr/123 pass4_prune and pr/124 pass4_prune2
    // stamp conn-4 on material they just rescued or shed -- 481.5 MeV in 4
    // events, e.g. SBND 18255-105074's two pdg-13 showers at 215.1 + 162.0
    // MeV inside the main cluster.  That energy is counted in NO object: the
    // prune removed it from its parent, whose kinematics were then recomputed.
    // Empty set when the knob is off => byte-identical.
    IndexedShowerSet conn4_keep_showers;
    if (m_kine_count_conn4_near && main_vertex) {
        const auto* main_cl = main_vertex->cluster();
        for (const ShowerPtr& shower : showers) {
            if (shower->get_start_vertex_and_type().second != 4) continue;
            if (!main_cl) continue;
            IndexedVertexSet sv; IndexedSegmentSet ss;
            shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
            double gap = -1.0;
            for (const auto& seg : ss) {
                for (const auto& f : seg->fits()) {
                    const double d = main_cl->get_closest_dis(f.point);
                    if (d >= 0 && (gap < 0 || d < gap)) gap = d;
                }
            }
            if (gap < 0 || gap > m_kine_conn4_near_gap) continue;
            conn4_keep_showers.insert(shower);
            SPDLOG_LOGGER_INFO(s_log,
                "kine_count_conn4_near: COUNT shower_id={} pdg={} ke_mev={:.2f} "
                "len_cm={:.1f} gap_cm={:.2f}",
                shower->get_shower_id(), shower->get_particle_type(),
                shower->get_kine_best() / units::MeV,
                shower->get_total_length() / units::cm, gap / units::cm);
        }
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

        // doc 85 sec 9: the shower's kine_best covers every member segment.
        {
            IndexedVertexSet vtmp;
            shower->fill_sets(vtmp, counted_segs, /*flag_exclude_start_segment=*/false);
        }

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2); // charge
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1); // range
        else
            ktree.kine_energy_info.push_back(0); // dQdx

        // Add rest-mass correction for non-electrons/positrons
        if (m_kine_charge.mass_rules) {
            SegmentPtr start_sg = shower->start_segment();
            const double mass = (start_sg && start_sg->particle_info()) ? start_sg->particle_info()->mass() : 0.0;
            if (pdg != 11 && start_sg && start_sg->particle_info())
                add_energy_legacy += static_cast<float>(mass / units::MeV);
            if (std::abs(pdg) == 2212) ++n_2212_shower_graph;
            ktree.kine_reco_add_energy += static_cast<float>(rest_term_rules(pdg, mass));
        }
        else if (pdg != 11) {
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
        counted_segs.insert(seg);   // doc 85 sec 9

        if (std::fabs(kine_best - kine_charge) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(2);
        else if (std::fabs(kine_best - kine_range) < 0.001 * kine_best)
            ktree.kine_energy_info.push_back(1);
        else
            ktree.kine_energy_info.push_back(0);

        ktree.kine_energy_included.push_back(include_flag);

        if (m_kine_charge.mass_rules) {
            if (pdg == 2212) add_energy_legacy += static_cast<float>(ave_binding_energy / units::MeV);
            else if (pdg != 11) add_energy_legacy += static_cast<float>(mass / units::MeV);
            ktree.kine_reco_add_energy += static_cast<float>(rest_term_rules(pdg, mass));
        }
        else if (pdg == 2212) { // proton: add binding energy
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
            // doc pr/101 (K5): a shower MEMBER attached to the main vertex
            // is already inside its shower's energy -- skip it the way the
            // BFS below does.  Off => legacy double count (prototype parity).
            if (m_kine_charge.mainvtx_used_guard && used_segments.count(seg)) {
                ++n_mainvtx_guard_skip;
                SPDLOG_LOGGER_INFO(s_log, "kine_mainvtx_guard: skipped shower-member seg idx={} pdg={} ke_mev={:.1f}",
                                   seg->get_graph_index(),
                                   (seg->particle_info() ? seg->particle_info()->pdg() : 0),
                                   (seg->particle_info() ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0));
                continue;
            }
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
                if (m_kine_charge.mass_rules) {
                    const double mass = prev_sg->particle_info()->mass();
                    if (prev_pdg == 2212) add_energy_legacy -= static_cast<float>(ave_binding_energy / units::MeV);
                    else if (prev_pdg != 11) add_energy_legacy -= static_cast<float>(mass / units::MeV);
                    ktree.kine_reco_add_energy -= static_cast<float>(rest_term_rules(prev_pdg, mass));
                }
                else if (prev_pdg == 2212) {
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
        // doc pr/128: a conn-4 shower kept by kine_count_conn4_near is the
        // candidate's own material and is counted like a conn-3 shower.
        // Empty keep set when the knob is off => legacy skip.
        if (vtx_type > 3 && !conn4_keep_showers.count(shower)) continue;

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
        if (m_kine_charge.mass_rules) {
            // doc pr/101 (K2): leftover showers are detached (conn-2/3)
            // objects.  Nucleon-typed ones get the binding energy behind the
            // legacy 5 cm start-segment gate (2212 as before, 2112 added);
            // mu/pi-typed ones stay MASSLESS as in legacy: the numu50 census
            // showed every such object (12 events, 30-119 cm conn-2 pieces
            // typed 13) sitting in an event whose muon mass is already
            // counted on the attached track -- a second piece of the same
            // muon, so a second mass would be the P4 double count.
            SegmentPtr start_sg = shower->start_segment();
            const bool pass_gate = start_sg && segment_track_length(start_sg) > 5.0 * units::cm;
            if (pdg == 2212 && pass_gate)
                add_energy_legacy += static_cast<float>(ave_binding_energy / units::MeV);
            if (pdg != 11) ++n_leftover_nonem;
            if (pass_gate && (std::abs(pdg) == 2212 || std::abs(pdg) == 2112))
                ktree.kine_reco_add_energy += static_cast<float>(ave_binding_energy / units::MeV);
        }
        else if (pdg == 2212) {
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
    // doc sbnd_xin/docs/pr/130 -- BYTE-NEUTRAL PROBE, no decision change.
    // kine_count_guard_freed (pr/123 r2, SBND ON) has exactly one predicate,
    // the kPass4GuardFreed flag -- no cluster, length, distance or
    // continuation term.  doc pr/129 found it counting 268.70 MeV of an
    // owner-adjudicated cosmic into kine_reco_Enu (SBND 18255-393505: the
    // segment chains end-to-end into one 295.0 cm object that never comes
    // within 68.9 cm of the neutrino vertex).  This tape reports, for every
    // candidate the pool is about to count, the four terms of the shipped
    // continuation test (doc pr/128) -- the question being whether that test
    // separates 393505 from 171572/94392 without any new threshold.
    //
    // Deliberately a SEPARATE block, not a `|| probe` widening of the
    // production block below (CLAUDE.md sec 2: the shipped path stays
    // byte-for-byte untouched).  The reference cloud is frozen at pool entry
    // so same-pool candidates cannot qualify each other.
    // Absent env var => nothing runs => byte-identical.
    // -------------------------------------------------------------------------
    const bool gf_probe = (std::getenv("WCT_KINE_GUARDFREED_PROBE") != nullptr);
    if (gf_probe && main_vertex) {
        std::vector<std::pair<size_t, SegmentPtr>> probe_cands;
        for (auto edesc : mir(boost::edges(graph))) {
            SegmentPtr seg = graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            if (used_segments.count(seg)) continue;
            if (map_sg_shower.count(seg)) continue;
            if (!seg->flags_any(SegmentFlags::kPass4GuardFreed)) continue;
            probe_cands.emplace_back(graph[edesc].index, seg);
        }
        std::sort(probe_cands.begin(), probe_cands.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        if (!probe_cands.empty()) {
            // Same reference cloud pr/128's pools test against: what this
            // tree actually counts (conn-4 showers only when kept), plus
            // every segment the traversal claimed.
            std::vector<SegmentPtr> ref_segs;
            {
                IndexedSegmentSet counted;
                for (const ShowerPtr& shower : showers) {
                    if (shower->get_start_vertex_and_type().second > 3 &&
                        !conn4_keep_showers.count(shower)) continue;
                    IndexedVertexSet sv; IndexedSegmentSet ss;
                    shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
                    counted.insert(ss.begin(), ss.end());
                }
                for (const auto& seg : used_segments) counted.insert(seg);
                for (const auto& seg : counted)
                    if (seg && !seg->fits().empty()) ref_segs.push_back(seg);
            }
            const auto* main_cl = main_vertex->cluster();
            const Point mv_pt = main_vertex->fit().point;
            for (const auto& [eidx, seg] : probe_cands) {
                const auto g = segment_continuation_geometry(seg, ref_segs);
                // Descriptive only (doc pr/129 ruled out a geometric bound):
                // closest approach to the main vertex, and the candidate's
                // own end-to-end extent.
                double d_mainvtx = -1.0, extent = 0.0;
                const auto& cf = seg->fits();
                if (!cf.empty()) {
                    d_mainvtx = 1e9;
                    for (const auto& f : cf)
                        d_mainvtx = std::min(d_mainvtx, (f.point - mv_pt).magnitude());
                    extent = (cf.back().point - cf.front().point).magnitude();
                }
                std::fprintf(stderr,
                    "PR130_GUARDFREED idx=%zu cluster=%d main_cluster=%d pdg=%d score=%.3f "
                    "ke_mev=%.2f len_cm=%.1f extent_cm=%.1f d_mainvtx_cm=%.2f "
                    "gap_cm=%.2f cand_end_cm=%.2f ref_end_cm=%.2f kink_deg=%.2f "
                    "ref_cluster=%d ref_idx=%d ref_len_cm=%.1f ref_pdg=%d "
                    "ref_ke_mev=%.2f ref_is_shower_start=%d ref_in_used=%d n_ref=%zu\n",
                    eidx,
                    seg->cluster() ? seg->cluster()->get_cluster_id() : -1,
                    main_cl ? main_cl->get_cluster_id() : -1,
                    (seg->particle_info() ? seg->particle_info()->pdg() : 0),
                    seg->particle_score(),
                    (seg->particle_info() ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0),
                    segment_track_length(seg) / units::cm,
                    extent / units::cm,
                    d_mainvtx < 0 ? -1.0 : d_mainvtx / units::cm,
                    g.gap < 0 ? -1.0 : g.gap / units::cm,
                    g.cand_end_dis > 1e8 ? -1.0 : g.cand_end_dis / units::cm,
                    g.ref_end_dis  > 1e8 ? -1.0 : g.ref_end_dis  / units::cm,
                    g.angle_deg,
                    (g.ref && g.ref->cluster()) ? g.ref->cluster()->get_cluster_id() : -1,
                    g.ref ? static_cast<int>(g.ref->get_graph_index()) : -1,
                    g.ref ? segment_track_length(g.ref) / units::cm : -1.0,
                    (g.ref && g.ref->particle_info()) ? g.ref->particle_info()->pdg() : 0,
                    (g.ref && g.ref->particle_info())
                        ? g.ref->particle_info()->kinetic_energy() / units::MeV : 0.0,
                    (g.ref && map_sg_shower.count(g.ref)) ? 1 : 0,
                    (g.ref && used_segments.count(g.ref)) ? 1 : 0,
                    ref_segs.size());
            }
        }
    }

    // -------------------------------------------------------------------------
    // doc sbnd_xin/docs/pr/123 round 2 (kine_count_guard_freed): count a
    // track the pass4 long-track guard declined (SegmentFlags::
    // kPass4GuardFreed) that neither the BFS nor any shower claimed.  The
    // kine twin of pf_orphan_guard_freed, so PF and kine describe the same
    // particle set (the pr/93 round-4 principle).  Outside the
    // kine_count_orphan_tracks scope on both axes: the freed track is
    // cross-cluster (main-cluster filter) and sentinel-scored (confident
    // test) -- SBND 18255-171572's 125cm muon, ~390 MeV silently absent
    // from kine_reco_Enu.  The flag is the predicate; the 120-segment
    // any-cluster unclaimed population stays uncounted.
    // C++ default false => no pass => byte-identical.
    // -------------------------------------------------------------------------
    if (m_kine_count_guard_freed) {
        std::vector<std::pair<size_t, SegmentPtr>> freed_cands;
        for (auto edesc : mir(boost::edges(graph))) {
            SegmentPtr seg = graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            if (used_segments.count(seg)) continue;
            if (map_sg_shower.count(seg)) continue;
            if (!seg->flags_any(SegmentFlags::kPass4GuardFreed)) continue;
            freed_cands.emplace_back(graph[edesc].index, seg);
        }
        std::sort(freed_cands.begin(), freed_cands.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        size_t n_freed = 0;
        for (const auto& [eidx, seg] : freed_cands) {
            // -----------------------------------------------------------------
            // doc sbnd_xin/docs/pr/129 -- THE POINTING TEST.  Owner, 2026-08-29:
            // "the key difference is the direction, if the direction of the
            //  track is point to the main vertex, it is more likely to be part
            //  of neutrino.  For overclustering they are generally not point to
            //  neutrino vertex."
            //
            // This pool is the ONLY one with no geometric admission test at
            // all -- its predicate is the kPass4GuardFreed flag and nothing
            // else -- which is how it came to count 268.70 MeV of an
            // owner-adjudicated cosmic into kine_reco_Enu (SBND 18255-393505).
            // The test asks the daughter question: does this track's line aim
            // back at the neutrino vertex, and does the vertex lie behind it?
            //
            // m_kine_guard_freed_impact == 0 => no test => every candidate
            // counted exactly as before => byte-identical.
            // -----------------------------------------------------------------
            if (m_kine_guard_freed_impact > 0 && main_vertex) {
                const auto p = segment_vertex_pointing(seg, main_vertex->fit().point);
                const bool aims = p.d_vtx >= 0 &&
                                  p.impact   <= m_kine_guard_freed_impact &&
                                  p.miss_deg <= m_kine_guard_freed_miss_deg;
                SPDLOG_LOGGER_INFO(s_log,
                    "kine_guard_freed_impact: seg idx={} cluster={} ke_mev={:.2f} "
                    "d_vtx_cm={:.2f} impact_cm={:.2f} miss_deg={:.1f} -> {}",
                    eidx, seg->cluster() ? seg->cluster()->get_cluster_id() : -1,
                    (seg->particle_info() ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0),
                    p.d_vtx < 0 ? -1.0 : p.d_vtx / units::cm,
                    p.impact > 1e8 ? -1.0 : p.impact / units::cm,
                    p.miss_deg, aims ? "COUNT" : "SKIP");
                if (!aims) continue;
            }
            used_segments.insert(seg);
            const int pdg = push_segment_kine(seg, 1);
            ++n_freed;
            SPDLOG_LOGGER_INFO(s_log,
                "kine_count_guard_freed: COUNT seg idx={} cluster={} pdg={} ke_mev={:.2f} len_cm={:.1f}",
                eidx, seg->cluster()->get_cluster_id(), pdg,
                (seg->particle_info() ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0),
                segment_track_length(seg) / units::cm);
        }
        if (n_freed) {
            SPDLOG_LOGGER_INFO(s_log, "kine_count_guard_freed: {} guard-freed track(s) added", n_freed);
        }
    }

    // -------------------------------------------------------------------------
    // doc sbnd_xin/docs/pr/128 (kine_count_near_cross_cluster): kine twin of
    // pf_orphan_near_cross_cluster.  An unclaimed CROSS-CLUSTER track segment
    // whose fit points TOUCH the counted candidate.  kine_count_orphan_tracks
    // (pr/93 r4) is main-cluster-only by an explicit cluster-id test, and
    // kine_count_guard_freed (pr/123 r2) needs the kPass4GuardFreed flag, so
    // this class -- 29 objects in 18 of 239 SBND events, 21 of them within
    // 10 cm of displayed content and many touching it at 0.00 cm -- is
    // counted by neither.  Predicate is segment_near_candidate_track, NOT
    // segment_orphan_confident_track: 10 of the 19 near-and-long objects
    // carry the score-100 sentinel stamp the latter rejects (doc pr/128 §3).
    //
    // Double counting: the candidate is not in used_segments and not a shower
    // start, and it lives in a different cluster, hence on different blobs
    // than any counted shower's charge integration -- so its charge is in no
    // counted object.  §7 of the doc verifies that per fired event rather
    // than resting on the argument.
    // C++ default false => no pass => byte-identical.
    // -------------------------------------------------------------------------
    if (m_kine_count_near_cross_cluster && main_vertex) {
        const auto* main_cl = main_vertex->cluster();
        // Reference cloud = what this tree actually counts: every counted
        // shower's segments (conn-4 ones only when kept) plus every segment
        // the traversal claimed.  Using raw used_segments would let a
        // candidate qualify by touching an uncounted conn-4 cosmic.
        std::vector<SegmentPtr> ref_segs;
        {
            IndexedSegmentSet counted;
            for (const ShowerPtr& shower : showers) {
                if (shower->get_start_vertex_and_type().second > 3 &&
                    !conn4_keep_showers.count(shower)) continue;
                IndexedVertexSet sv; IndexedSegmentSet ss;
                shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
                counted.insert(ss.begin(), ss.end());
            }
            for (const auto& seg : used_segments) counted.insert(seg);
            for (const auto& seg : counted)
                if (seg && !seg->fits().empty()) ref_segs.push_back(seg);
        }
        std::vector<std::pair<size_t, SegmentPtr>> near_cands;
        std::map<SegmentPtr, double, SegmentIndexCmp> near_gap;
        if (!ref_segs.empty()) {
            for (auto edesc : mir(boost::edges(graph))) {
                SegmentPtr seg = graph[edesc].segment;
                if (!seg || !seg->descriptor_valid()) continue;
                if (used_segments.count(seg)) continue;
                if (map_sg_shower.count(seg)) continue;   // shower starts stay shower-owned
                const auto* cl = seg->cluster();
                if (!main_cl || !cl) continue;
                if (cl->get_cluster_id() == main_cl->get_cluster_id()) continue;  // pr/93 r4 owns this
                if (!segment_near_candidate_track(seg, m_kine_near_min_len)) continue;
                // Continuation, not mere proximity -- see
                // segment_continuation_geometry's header (SBND 18255-72786:
                // gap alone admitted two cosmics and +1151 MeV).
                const auto g = segment_continuation_geometry(seg, ref_segs);
                if (g.gap < 0 || g.gap > m_kine_near_gap) continue;
                if (g.cand_end_dis > m_kine_near_end_tol) continue;
                if (g.ref_end_dis  > m_kine_near_end_tol) continue;
                if (g.angle_deg    > m_kine_near_kink_deg) continue;
                near_gap[seg] = g.gap;
                near_cands.emplace_back(graph[edesc].index, seg);
            }
        }
        std::sort(near_cands.begin(), near_cands.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        size_t n_near = 0;
        for (const auto& [eidx, seg] : near_cands) {
            used_segments.insert(seg);
            const int pdg = push_segment_kine(seg, 1);
            ++n_near;
            SPDLOG_LOGGER_INFO(s_log,
                "kine_count_near_cross_cluster: COUNT seg idx={} cluster={} pdg={} score={:.3f} "
                "ke_mev={:.2f} len_cm={:.1f} gap_cm={:.2f}",
                eidx, seg->cluster()->get_cluster_id(), pdg, seg->particle_score(),
                (seg->particle_info() ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0),
                segment_track_length(seg) / units::cm, near_gap.at(seg) / units::cm);
        }
        if (n_near) {
            SPDLOG_LOGGER_INFO(s_log,
                "kine_count_near_cross_cluster: {} near cross-cluster track(s) added", n_near);
        }
    }

    if (m_kine_charge.mass_rules || m_kine_charge.mainvtx_used_guard) {
        SPDLOG_LOGGER_INFO(s_log,
            "kine_mass_census: add_legacy={:.1f} add_rules={:.1f} n_2212_showers_graph={} n_leftover_nonEM={} n_mainvtx_guard_skip={}",
            add_energy_legacy, ktree.kine_reco_add_energy, n_2212_shower_graph, n_leftover_nonem, n_mainvtx_guard_skip);
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
    // doc 85 sec 9 -- the excluded-energy census (owner request 2026-08-30).
    //
    // Everything above is finished; this block only READS.  It appends to no
    // kine_* vector, writes only the five new KineInfo scalars, and never
    // touches kine_reco_Enu -- so the Enu of every event is bit-identical to
    // what the same binary would produce without this block.
    //
    // What it answers: kine_reco_Enu sums the objects the main-vertex walk
    // reached; how much reconstructed kinetic energy does the SAME PR graph
    // hold that the walk never summed?  That is the quantity docs pr/128 and
    // pr/131 could only measure behind an env-gated stderr probe.  Here it
    // becomes a stored variable, per candidate, always.
    //
    // Split (doc pr/131's own distinction): same cluster as the main vertex
    // ("main" -- material the candidate arguably owns, the population pr/65's
    // audit line can see) vs a different cluster ("other" -- the
    // cross-cluster activity named nowhere).  Cluster identity is compared by
    // id, never by pointer order (CLAUDE.md determinism rule); the graph edge
    // walk is the same mir(boost::edges(...)) the pr/93 r4 orphan pass uses.
    // -------------------------------------------------------------------------
    {
        const auto* main_cl = main_vertex ? main_vertex->cluster() : nullptr;
        const int main_cid = main_cl ? main_cl->get_cluster_id() : -1;

        double exc_main = 0, exc_other = 0;
        int n_exc = 0, n_nocl = 0, n_graph = 0, n_graph_main = 0;
        for (auto edesc : mir(boost::edges(graph))) {
            const SegmentPtr& seg = graph[edesc].segment;
            if (!seg) continue;
            ++n_graph;
            const auto* cl = seg->cluster();
            const bool same_cl = (cl && main_cid >= 0 && cl->get_cluster_id() == main_cid);
            if (same_cl) ++n_graph_main;
            if (counted_segs.count(seg)) continue;
            if (!seg->has_particle_info() || !seg->particle_info()) continue;
            const double ke = seg->particle_info()->kinetic_energy() / units::MeV;
            if (!(ke > 0)) continue;
            ++n_exc;
            if (!cl) ++n_nocl;
            if (same_cl) exc_main += ke;
            else exc_other += ke;
        }
        ktree.kine_energy_excluded_main  = static_cast<float>(exc_main);
        ktree.kine_energy_excluded_other = static_cast<float>(exc_other);
        ktree.kine_energy_excluded       = static_cast<float>(exc_main + exc_other);
        ktree.kine_n_excluded            = n_exc;

        // The rows the tree itself flags as not-included.  These ARE inside
        // kine_reco_Enu (the sum above is unconditional, in the toolkit and in
        // the prototype alike) -- reported so the two quantities cannot be
        // confused for one another.
        double flagged = 0;
        const size_t nrow = std::min(ktree.kine_energy_particle.size(),
                                     ktree.kine_energy_included.size());
        for (size_t i = 0; i < nrow; ++i) {
            if (ktree.kine_energy_included[i] != 1) flagged += ktree.kine_energy_particle[i];
        }
        ktree.kine_energy_flagged = static_cast<float>(flagged);

        // The three census counters after n_excluded are what makes a zero
        // readable: main=0 with n_graph_main > n_counted means the walk really
        // did consume every main-cluster segment, while main=0 with
        // n_graph_main=0 would mean the cluster-id compare never matched.
        SPDLOG_LOGGER_INFO(s_log,
            "kine_excluded_census: Enu={:.1f} excluded={:.1f} (main={:.1f} other={:.1f}) "
            "n_excluded={} flagged_inside_Enu={:.1f} main_cid={} n_graph_seg={} "
            "n_graph_main={} n_counted={} n_excluded_nocluster={}",
            ktree.kine_reco_Enu, ktree.kine_energy_excluded,
            ktree.kine_energy_excluded_main, ktree.kine_energy_excluded_other,
            ktree.kine_n_excluded, ktree.kine_energy_flagged,
            main_cid, n_graph, n_graph_main, (int)counted_segs.size(), n_nocl);
    }

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
