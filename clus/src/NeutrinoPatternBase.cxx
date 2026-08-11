#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellUtil/Logging.h"

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_map>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

// Named logger for this file.  At runtime set its level independently:
//   Log::set_level("debug", "clus.NeutrinoPattern");   // this file only
//   Log::set_level("debug", "clus");                    // whole clus subsystem
// At build time it is compiled away entirely when --with-spdlog-active-level > debug.
static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// Sentinel returned by proto_extend_point when no steiner point can be found.
static constexpr size_t INVALID_STEINER_INDEX = std::numeric_limits<size_t>::max();


std::vector<VertexPtr> PatternAlgorithms::find_cluster_vertices(Graph& graph, const Facade::Cluster& cluster)
{
    std::vector<VertexPtr> result;

    // Iterate in insertion order for deterministic results
    for (const auto& vd : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd].vertex;
        if (vtx && vtx->cluster() && vtx->cluster() == &cluster) {
            result.push_back(vtx);
        }
    }

    return result;
}

std::vector<SegmentPtr> PatternAlgorithms::find_cluster_segments(Graph& graph, const Facade::Cluster& cluster)
{
    std::vector<SegmentPtr> result;

    // Iterate in insertion order for deterministic results
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr seg = graph[ed].segment;
        if (seg && seg->cluster() && seg->cluster() == &cluster) {
            result.push_back(seg);
        }
    }

    return result;
}

bool PatternAlgorithms::seg_dir_weak(SegmentPtr seg) const
{
    return m_dir_weak_use_score ? segment_is_dir_weak(seg) : seg->dir_weak();
}

bool PatternAlgorithms::clean_up_graph(Graph& graph, const Facade::Cluster& cluster)
{
    bool modified = false;
    
    // First, find and remove all segments associated with this cluster
    std::vector<SegmentPtr> segments_to_remove = find_cluster_segments(graph, cluster);
    for (auto seg : segments_to_remove) {
        if (remove_segment(graph, seg)) {
            modified = true;
        }
    }
    
    // Then, find and remove all vertices associated with this cluster
    // Note: vertices that are still connected to other segments won't be removed
    // until their segments are removed first
    std::vector<VertexPtr> vertices_to_remove = find_cluster_vertices(graph, cluster);
    for (auto vtx : vertices_to_remove) {
        if (remove_vertex(graph, vtx)) {
            modified = true;
        }
    }
    
    return modified;
}

std::vector<Facade::geo_point_t> PatternAlgorithms::do_rough_path(const Facade::Cluster& cluster,Facade::geo_point_t& first_point, Facade::geo_point_t& last_point){
        // A present-but-empty steiner cloud is legal (SBND MC evt 11 --
        // Dataset::size() counts arrays, not points): the kNN then returns an
        // empty result and [0] below reads past the end.  Return the empty
        // path early; every caller already treats an empty path as "no rough
        // path" (the !has_pc branch below returns the same).
        if (!cluster.has_pc("steiner_pc") ||
            cluster.get_pc("steiner_pc").size_major() == 0) {
            return {};
        }

        // Find closest indices in the steiner point cloud
        auto first_knn_results = cluster.kd_steiner_knn(1, first_point, "steiner_pc");
        auto last_knn_results = cluster.kd_steiner_knn(1, last_point, "steiner_pc");
        if (first_knn_results.empty() || last_knn_results.empty()) {
            return {};
        }

        auto first_index = first_knn_results[0].first;  // Get the index from the first result
        auto last_index = last_knn_results[0].first;   // Get the index from the first result
 
        // 4. Use Steiner graph to find the shortest path
        //
        // doc pr/55 (2026-08-09): the owner asked, after doc pr/53's
        // relaxed_strict_img membership-graph fix, whether the PR/fitting
        // chain actually runs on that tightened graph. It does not: this is
        // the ONLY graph the trajectory router touches, and it is
        // "steiner_graph" -- a Steiner reduction of "ctpc_ref_pid"
        // (CreateSteinerGraph.cxx), which carries the uncapped MST
        // connect_graph()/connect_graph_with_reference() bridges and has no
        // relationship to protect_bundle's graph_name knob at all. Sentinel
        // log only, no behavior change.
        const auto& steiner_graph = cluster.find_graph("steiner_graph");
        SPDLOG_LOGGER_DEBUG(s_log,
            "pr55 do_rough_path: cluster {} flavor=steiner_graph n_vertices={} n_edges={} "
            "first=({:.1f},{:.1f},{:.1f}) last=({:.1f},{:.1f},{:.1f})",
            cluster.ident(), boost::num_vertices(steiner_graph), boost::num_edges(steiner_graph),
            first_point.x(), first_point.y(), first_point.z(),
            last_point.x(), last_point.y(), last_point.z());
        const std::vector<size_t>& path_indices =
            cluster.graph_algorithms("steiner_graph").shortest_path(first_index, last_index);

        std::vector<Facade::geo_point_t> path_points;
        if (!cluster.has_pc("steiner_pc")) return path_points;
        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

        for (size_t idx : path_indices) {
            path_points.emplace_back(x_coords[idx], y_coords[idx], z_coords[idx]);
        }
        return path_points;
}

void PatternAlgorithms::set_default_shower_particle_info(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex) {
    // Mirrors prototype ProtoSegment::get_particle_type() which always returns 11 for
    // any shower segment (flag_shower_trajectory || flag_shower_topology).
    // Any segment flagged as a shower but missing particle_info (e.g. because it was
    // newly classified as kShowerTrajectory after determine_direction ran) gets PDG=11.
    //
    // doc sbnd_xin/docs/pr/40 round 2 F5: an electron cannot father a proton.
    // This function's DEFAULT-fill role (below) is the single choke point
    // where a flag_shower segment with no PID gets defaulted to electron --
    // but it is NOT the only writer of pdg 11: determine_direction's
    // kShowerTopology branch (NeutrinoTrackShowerSep.cxx, prototype
    // determine_dir_shower_topology parity) writes type+mass UNCONDITIONALLY
    // at STAGE 3, before determine_main_vertex has run -- there is no
    // main_vertex yet at that call site to test against, so the check cannot
    // live there.  SBND evt 256587 seg 11079 (a 29 cm segment emanating from
    // the neutrino vertex whose far end abuts a PID'd, charge-confirmed
    // proton, segment 11080) is exactly this case: it already carries pdg 11
    // by the time this function -- called post-main_vertex, from
    // examine_direction -- ever sees it.  So the guard below runs as an
    // OVERRIDE on every existing pdg-11 shower segment, not only a DEFAULT
    // on segments still missing particle_info.  Designed divergence, not a
    // port-fidelity fix (the prototype has no such veto); see
    // porting_dictionary.md.  false = legacy = byte-identical.
    // ordered_edges: stable edge-index order (boost::edges is pointer order).
    static const bool dbg = std::getenv("WCT_PROTON_DAUGHTER_DEBUG") != nullptr;
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        if (!sg->flags_any(SegmentFlags::kShowerTrajectory) &&
            !sg->flags_any(SegmentFlags::kShowerTopology)) continue;

        if (dbg) {
            std::fprintf(stderr, "SET_DEFAULT_SHOWER_DEBUG clus=%d idx=%zu has_pi=%d pdg=%d m_shower_proton_daughter_pion=%d main_vertex=%p\n",
                         cluster.get_cluster_id(), sg->get_graph_index(), (int)sg->has_particle_info(),
                         sg->has_particle_info() ? sg->particle_info()->pdg() : 0,
                         (int)m_shower_proton_daughter_pion, (void*)main_vertex.get());
        }

        if (sg->has_particle_info()) {
            // OVERRIDE path: an earlier stage already assigned a pdg.  Only
            // touch it if that pdg is electron (11) -- or, under F7 below,
            // an already-relabelled pion (211) that still needs its shower
            // flags cleared -- and the proton-daughter rule fires; never
            // override a track PID.
            //
            // doc pr/40 round 4 F7: examine_direction can run more than
            // once, and the pdg()==211 re-entry case is exactly that second
            // pass -- widening the test to 11||211 makes the branch
            // idempotent (re-deriving the same 4-momentum for an
            // already-211 segment is a no-op) instead of silently skipping
            // the flag-clear on re-entry.
            const int cur_pdg = sg->particle_info()->pdg();
            const bool reentrant_pion = m_shower_proton_daughter_pion_dissolve && cur_pdg == 211;
            if (m_shower_proton_daughter_pion && (cur_pdg == 11 || reentrant_pion) &&
                segment_has_proton_daughter(graph, sg, main_vertex, m_mip_dqdx_median)) {
                const int pdg_code = 211;
                auto four_momentum = segment_cal_4mom(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx_median);
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    pdg_code, particle_data->get_particle_mass(pdg_code),
                    particle_data->pdg_to_name(pdg_code), four_momentum);
                sg->particle_info(pinfo);
                sg->particle_score(100.0);
                // F7: a pion is a track, not a shower -- see
                // m_shower_proton_daughter_pion_dissolve's docstring.
                if (m_shower_proton_daughter_pion_dissolve) {
                    sg->unset_flags(SegmentFlags::kShowerTrajectory);
                    sg->unset_flags(SegmentFlags::kShowerTopology);
                }
            }
            continue;
        }

        // DEFAULT-fill path: no particle_info at all yet.
        int pdg_code = 11;
        if (m_shower_proton_daughter_pion &&
            segment_has_proton_daughter(graph, sg, main_vertex, m_mip_dqdx_median)) {
            pdg_code = 211;
        }

        auto four_momentum = segment_cal_4mom(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx_median);
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            pdg_code,
            particle_data->get_particle_mass(pdg_code),
            particle_data->pdg_to_name(pdg_code),
            four_momentum
        );
        sg->particle_info(pinfo);
        sg->particle_score(100.0);
        // F7 (doc pr/40 round 4): same flag-clear on the default-fill path,
        // in case a fresh segment is relabelled pion on its first pass.
        if (m_shower_proton_daughter_pion_dissolve && pdg_code == 211) {
            sg->unset_flags(SegmentFlags::kShowerTrajectory);
            sg->unset_flags(SegmentFlags::kShowerTopology);
        }
    }
}

void PatternAlgorithms::override_muon_multi_proton_pion(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex) {
    // doc sbnd_xin/docs/pr/40 round 4 F8 -- see m_muon_multi_proton_pion's
    // docstring in NeutrinoPatternBase.h.  false = legacy = no-op.
    if (!m_muon_multi_proton_pion) return;

    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        if (sg->flags_any(SegmentFlags::kShowerTrajectory) ||
            sg->flags_any(SegmentFlags::kShowerTopology)) continue;
        if (!sg->has_particle_info() || std::abs(sg->particle_info()->pdg()) != 13) continue;

        if (!segment_at_multi_proton_vertex(graph, sg, main_vertex, m_mip_dqdx_median, /*min_protons=*/2)) continue;

        const int pdg_code = 211;
        auto four_momentum = segment_cal_4mom(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx_median);
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            pdg_code, particle_data->get_particle_mass(pdg_code),
            particle_data->pdg_to_name(pdg_code), four_momentum);
        sg->particle_info(pinfo);
        sg->particle_score(100.0);
    }
}

void PatternAlgorithms::override_michel_stem_muon(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex) {
    // doc sbnd_xin/docs/pr/40 round 6 F14 -- see m_michel_stem_muon_rescue's
    // docstring in NeutrinoPatternBase.h.  Widens the toolkit's own Michel
    // rescue ("a stopped proton cannot produce a Michel electron",
    // NeutrinoVertexFinder.cxx examine_direction's weak-direction branch):
    // that rescue only runs for seg_dir_weak segments and only when the
    // stopping vertex has EXACTLY two segments; a stopping muon whose Bragg
    // rise wins the proton template with a confident direction, ending at a
    // vertex that also carries stub fragments (SBND evt 54341 seg 18005,
    // stopping vertex degree 4), is out of its reach on both counts.
    // false = legacy = no-op.
    if (!m_michel_stem_muon_rescue) return;
    if (!main_vertex || !main_vertex->descriptor_valid()) return;

    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        if (sg->flags_any(SegmentFlags::kShowerTrajectory) ||
            sg->flags_any(SegmentFlags::kShowerTopology)) continue;
        if (!sg->has_particle_info() || !sg->particle_info()) continue;
        if (sg->particle_info()->pdg() != 2212) continue;
        // A genuine stopping proton is short; the rescue only claims a stem
        // that is long and straight (same helper as F10/F11/F12).
        if (!segment_is_straight_long_track(sg)) continue;

        // The stem must emanate from the main vertex; the Michel topology
        // test runs on the OTHER (stopping) end.
        auto vpair = find_vertices(graph, sg);
        VertexPtr far_vtx = nullptr;
        if (vpair.first == main_vertex) far_vtx = vpair.second;
        else if (vpair.second == main_vertex) far_vtx = vpair.first;
        else continue;
        if (!far_vtx || !far_vtx->descriptor_valid()) continue;

        // Same sibling test as the existing Michel branch (kShowerTrajectory
        // flag or electron PDG), degree restriction relaxed from ==2 to >=1
        // shower-like sibling among any number of stub fragments.
        bool flag_michel = false;
        for (auto edesc : sorted_out_edges(far_vtx->get_descriptor(), graph)) {
            SegmentPtr other_sg = graph[edesc].segment;
            if (!other_sg || other_sg == sg) continue;
            if (other_sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                (other_sg->has_particle_info() && other_sg->particle_info() &&
                 std::abs(other_sg->particle_info()->pdg()) == 11)) {
                flag_michel = true;
                break;
            }
        }
        if (!flag_michel) continue;

        // Relabel muon, recompute the 4-momentum -- same conventions as the
        // existing rescue's recompute (m_mip_dqdx scale, no score write).
        const int pdg_code = 13;
        auto four_momentum = segment_cal_4mom(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx);
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            pdg_code, particle_data->get_particle_mass(pdg_code),
            particle_data->pdg_to_name(pdg_code), four_momentum);
        sg->particle_info(pinfo);
    }
}

// doc sbnd_xin/docs/pr/43 round 2 K3 -- docstring in NeutrinoPatternBase.h.
// Late particle-info/flag reconciliation: runs after shower_clustering_with_nv
// and before the taggers so every downstream consumer (tagger features, kine,
// Bee PF tree, PR display) sees one consistent labeling.  false = no-op.
void PatternAlgorithms::reconcile_particle_flags(Graph& graph, VertexPtr main_vertex,
    IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower,
    ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower,
    ShowerIntMap& map_shower_pio_id,
    const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    if (!m_pid_flag_reconcile) return;
    if (!main_vertex || !main_vertex->descriptor_valid()) return;

    const char* dbg_env = std::getenv("WCT_PID_WRITE_DEBUG");
    const bool dbg = dbg_env && dbg_env[0] && dbg_env[0] != '0';

    Facade::Cluster* main_clu = nullptr;
    for (auto ed : sorted_out_edges(main_vertex->get_descriptor(), graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (sg) { main_clu = sg->cluster(); break; }
    }

    auto vertex_degree = [&graph](VertexPtr v) -> int {
        if (!v || !v->descriptor_valid()) return 0;
        int n = 0;
        for (auto ed : sorted_out_edges(v->get_descriptor(), graph)) if (graph[ed].segment) ++n;
        return n;
    };

    // A single-segment wrapper Shower whose segment is now a confirmed track
    // renders "shower  X MeV" from its stale cache; dissolve it so the
    // segment shows as a live track node.  Long-muon pseudo-showers (cached
    // type +-13) and pi0-paired showers are exempt; multi-segment showers
    // are structure, not wrappers, and are left alone.
    auto dissolve_wrapper = [&](SegmentPtr seg) {
        auto it = map_segment_in_shower.find(seg);
        if (it == map_segment_in_shower.end()) return;
        ShowerPtr sh = it->second;
        if (!sh) { map_segment_in_shower.erase(it); return; }
        if (std::abs(sh->get_particle_type()) == 13) return;
        if (sh->edges().size() > 1) return;
        auto pio_it = map_shower_pio_id.find(sh);
        if (pio_it != map_shower_pio_id.end() && pio_it->second >= 0) return;
        map_segment_in_shower.erase(it);
        showers.erase(sh);
        for (auto mit = map_vertex_in_shower.begin(); mit != map_vertex_in_shower.end();) {
            if (mit->second == sh) mit = map_vertex_in_shower.erase(mit);
            else ++mit;
        }
        for (auto& vs : map_vertex_to_shower) vs.second.erase(sh);
        if (dbg) std::fprintf(stderr, "PID_RECONCILE dissolve wrapper shower of seg (clus=%d idx=%zu)\n",
                              seg->cluster() ? seg->cluster()->get_cluster_id() : -1, seg->get_graph_index());
    };

    // ---- Rule 1: forced-electron terminal rescue behind a main-vertex
    // proton chain (evt 57661: proton 18003 -> stub 18005 -> seg 18007
    // forced pdg 11 / sentinel score 100 by the unconditional branches of
    // segment_determine_shower_direction_trajectory, never actually PID'd).
    for (auto e_it : sorted_out_edges(main_vertex->get_descriptor(), graph)) {
        SegmentPtr psg = graph[e_it].segment;
        if (!psg || !psg->has_particle_info() || !psg->particle_info()) continue;
        if (psg->particle_info()->pdg() != 2212) continue;

        auto pvs = find_vertices(graph, psg);
        VertexPtr cur_v = (pvs.first == main_vertex) ? pvs.second
                        : (pvs.second == main_vertex ? pvs.first : nullptr);
        SegmentPtr prev = psg;
        std::vector<SegmentPtr> stubs;
        for (int hop = 0; hop < 3 && cur_v && cur_v->descriptor_valid(); ++hop) {
            SegmentPtr next = nullptr;
            int nnbr = 0;
            for (auto ed : sorted_out_edges(cur_v->get_descriptor(), graph)) {
                SegmentPtr nbr = graph[ed].segment;
                if (!nbr || nbr == prev) continue;
                ++nnbr;
                next = nbr;
            }
            if (nnbr != 1 || !next) break;   // branching vertex or dead end
            const bool has_pi = next->has_particle_info() && next->particle_info();
            const int pdg = has_pi ? next->particle_info()->pdg() : 0;
            if (has_pi && pdg == 11 && next->particle_score() == 100.0) {
                // Terminal candidate: re-run ordinary track PID and trust a
                // confident non-electron conclusion.  Fresh ParticleInfo is
                // written by segment_determine_dir_track itself; on an
                // electron/no conclusion restore the prior state untouched.
                auto nvs = find_vertices(graph, next);
                const int start_n = vertex_degree(nvs.first);
                const int end_n = vertex_degree(nvs.second);
                auto saved_info = next->particle_info();
                const double saved_score = next->particle_score();
                segment_determine_dir_track(next, start_n, end_n, particle_data, recomb_model, m_mip_dqdx);
                const int new_pdg = (next->has_particle_info() && next->particle_info())
                                        ? next->particle_info()->pdg() : 0;
                if (std::abs(new_pdg) == 13 || std::abs(new_pdg) == 211 || new_pdg == 2212) {
                    next->unset_flags(SegmentFlags::kShowerTrajectory);
                    next->unset_flags(SegmentFlags::kShowerTopology);
                    dissolve_wrapper(next);
                    for (auto stub : stubs) {
                        if (!stub->has_particle_info() || !stub->particle_info() ||
                            std::abs(stub->particle_info()->pdg()) != 13) continue;
                        if (segment_track_length(stub) > 15*units::cm) continue;
                        auto fm = segment_cal_4mom(stub, 211, particle_data, recomb_model, m_mip_dqdx);
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(211,
                            particle_data->get_particle_mass(211), particle_data->pdg_to_name(211), fm);
                        stub->particle_info(pinfo);
                    }
                    if (dbg) std::fprintf(stderr, "PID_RECONCILE terminal rescue seg (clus=%d idx=%zu) 11 -> %d, %zu stub(s) -> 211\n",
                                          next->cluster() ? next->cluster()->get_cluster_id() : -1,
                                          next->get_graph_index(), new_pdg, stubs.size());
                }
                else {
                    next->particle_info(saved_info);
                    next->particle_score(saved_score);
                }
                break;
            }
            const bool is_shower_nbr = next->flags_any(SegmentFlags::kShowerTrajectory) ||
                                       next->flags_any(SegmentFlags::kShowerTopology) ||
                                       (has_pi && std::abs(pdg) == 11);
            if (is_shower_nbr || (has_pi && pdg == 2212)) break;
            stubs.push_back(next);
            auto nvs = find_vertices(graph, next);
            VertexPtr nv = (nvs.first == cur_v) ? nvs.second
                         : (nvs.second == cur_v ? nvs.first : nullptr);
            prev = next;
            cur_v = nv;
        }
    }

    // ---- Rule 2: consistency guard, main cluster only.  A segment whose
    // final pdg is a confirmed track type must not keep stale shower flags
    // (the "track still displayed as Shower" class), nor a stale
    // single-segment wrapper Shower.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != main_clu) continue;
        if (!sg->has_particle_info() || !sg->particle_info()) continue;
        const int apdg = std::abs(sg->particle_info()->pdg());
        if (apdg != 13 && apdg != 211 && apdg != 2212) continue;
        if (sg->flags_any(SegmentFlags::kShowerTrajectory) ||
            sg->flags_any(SegmentFlags::kShowerTopology)) {
            sg->unset_flags(SegmentFlags::kShowerTrajectory);
            sg->unset_flags(SegmentFlags::kShowerTopology);
            if (dbg) std::fprintf(stderr, "PID_RECONCILE clear stale shower flags (clus=%d idx=%zu pdg=%d)\n",
                                  sg->cluster() ? sg->cluster()->get_cluster_id() : -1,
                                  sg->get_graph_index(), sg->particle_info()->pdg());
        }
        dissolve_wrapper(sg);
    }
}

std::vector<Facade::geo_point_t> PatternAlgorithms::do_rough_path_reg_pc(const Facade::Cluster& cluster, Facade::geo_point_t& first_point, Facade::geo_point_t& last_point,  std::string graph_name){
    // Find closest indices in the regular point cloud using kd_knn
    auto first_knn_results = cluster.kd_knn(1, first_point);
    auto last_knn_results = cluster.kd_knn(1, last_point);
    // Empty kNN (zero-point cluster) -> [0] reads past the end; return the
    // empty path instead, same contract as do_rough_path above.
    if (first_knn_results.empty() || last_knn_results.empty()) {
        return {};
    }

    auto first_index = first_knn_results[0].first;  // Get the index from the first result
    auto last_index = last_knn_results[0].first;   // Get the index from the first result
    
    // Use the specified graph to find the shortest path
    const std::vector<size_t>& path_indices = 
        cluster.graph_algorithms(graph_name).shortest_path(first_index, last_index);
    
    // Convert indices to points using the regular point cloud
    std::vector<Facade::geo_point_t> path_points;
    const auto& points = cluster.points();  // Returns array of coordinate arrays [x_coords, y_coords, z_coords]
    const auto& x_coords = points[0];
    const auto& y_coords = points[1];
    const auto& z_coords = points[2];
    
    for (size_t idx : path_indices) {
        path_points.emplace_back(x_coords[idx], y_coords[idx], z_coords[idx]);
    }
    
    return path_points;
}


SegmentPtr PatternAlgorithms::create_segment_for_cluster(WireCell::Clus::Facade::Cluster& cluster, IDetectorVolumes::pointer dv, const std::vector<Facade::geo_point_t>& path_points, int dir){
     // Step 3: Prepare segment data
    std::vector<PR::WCPoint> wcpoints;
    // const auto transform = m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
    // Step 4: Create segment connecting the vertices
    auto segment = PR::make_segment();
    
    // create and associate Dynamic Point Cloud
    for (const auto& point : path_points) {
        PR::WCPoint wcp;
        wcp.point = point; 
        wcpoints.push_back(wcp);
    }

    // Step 5: Configure the segment
    segment->wcpts(wcpoints).cluster(&cluster).dirsign(dir); // direction: +1, 0, or -1
            
    // auto& wcpts = segment->wcpts();
    // for (size_t i=0;i!=path_points.size(); i++){
    //     std::cout << "A: " << i << " " << path_points.at(i) << " " << wcpts.at(i).point << std::endl;
    // }
    create_segment_point_cloud(segment, path_points, dv, "main");

    return segment;
}

SegmentPtr PatternAlgorithms::create_segment_from_vertices(Graph& graph, Facade::Cluster& cluster, VertexPtr v1, VertexPtr v2, IDetectorVolumes::pointer dv){
     // Create Segment using the vertices to derive a path 
    auto path_points = do_rough_path(cluster, v1->wcpt().point, v2->wcpt().point);
    
    // Check if path has enough points (similar to WCPPID check)
    if (path_points.size() <= 1) {
        return nullptr;
    }
    
    auto seg = create_segment_for_cluster(cluster, dv, path_points);
    WireCell::Clus::PR::add_segment(graph, seg, v1, v2);
    return seg;
}




bool PatternAlgorithms::find_iso_first_segment_endpoints(const Facade::Cluster& cluster,
                                                         Facade::geo_point_t& p1,
                                                         Facade::geo_point_t& p2) const
{
    // Isochronous first-segment endpoints (doc sbnd_xin/docs/pr/24 round 2).
    //
    // The legacy boundary metric (Facade_Cluster.cxx calculate_boundary_metric,
    // a faithful port of prototype PR3DCluster_path.h:530-536) scores a
    // candidate pair by |dx|/one_tick plus the count of live wires spanned in
    // U+V+W -- the wire-delta terms carry *0.0 coefficients.  On an
    // isochronous cluster (charge confined to a narrow drift slab, imaged as
    // a filled 2-D sheet) |dx| ~ 0, so the winner is the widest wire-footprint
    // pair: two corners on the sheet's edges, not the tip of its axis.  The
    // cheapest Steiner-graph paths between such corners are the sheet's
    // boundary edges, which is the observed two-edge-track fan.
    //
    // Here instead: endpoints = the extreme points along the principal axis of
    // the charge-qualified points, restricted to a tube around the axis line so
    // that on a sheet the pick is the axial tip at mid-width (the narrow
    // corner) rather than an edge corner, then snapped to a Steiner terminal
    // exactly like the legacy path.
    //
    // ROUND 3 (doc pr/24 sec 15).  Round 2 took the QUANTILE-TRIMMED axial
    // extreme and then the band point nearest the end-band centroid.  Both
    // steps pull the endpoint inward -- 8.4 / 15.1 cm on the 7200-point track
    // of SBND mcp1k evt 284794 -- and since nothing downstream of
    // init_first_segment is modified, find_other_segments correctly claimed the
    // uncovered tip as a separate segment: a 260.1 + 15.7 cm pair joined at a
    // 0.9 degree "vertex" in the middle of a straight track (evt 59899: an 8.0
    // degree junction where the legacy trajectory bends 1.6-2.8 degrees).  The
    // quantile trim now serves only the GATE measurements, where it belongs;
    // the endpoints are the true untrimmed tube extremes.  The aspect gate
    // keeps 1-D track-like clusters out of this branch altogether.

    const double length = cluster.get_length();
    if (length < m_iso_endpoint_min_length) return false;

    const int npts = cluster.npoints();
    if (npts < 10) return false;

    // Charge-qualified default-scope (T0-corrected) points: same exclusion
    // and blob-charge cut as the legacy boundary search
    // (Facade_Cluster.cxx get_two_boundary_wcps).  Keyed lookups only -- the
    // blob-charge memo is never iterated, so pointer keys are determinism-safe.
    std::unordered_map<const Facade::Blob*, double> blob_total_charge;
    std::vector<int> qual;
    qual.reserve(npts);
    for (int i = 0; i < npts; ++i) {
        if (cluster.is_point_excluded(i)) continue;
        const auto* blob = cluster.blob_with_point(i);
        auto bc_it = blob_total_charge.find(blob);
        if (bc_it == blob_total_charge.end()) {
            bc_it = blob_total_charge.emplace(blob, blob->estimate_total_charge()).first;
        }
        if (bc_it->second < 1500) continue;
        qual.push_back(i);
    }
    if (qual.size() < 10) return false;

    // Quantile-trimmed drift-x extent in the corrected frame.  A raw min-max
    // or blob-centre extent over-reads on this topology (multi-(apa,face)
    // transform inflation and stray-point tails, doc pr/24 sec 4.1): the
    // trimmed extent is the honest slab thickness.
    std::vector<double> xs;
    xs.reserve(qual.size());
    for (int i : qual) xs.push_back(cluster.point3d(i).x());
    std::sort(xs.begin(), xs.end());
    const size_t nq = xs.size();
    const size_t klo = static_cast<size_t>(m_iso_endpoint_xext_quantile * (nq - 1));
    const double xext = xs[nq - 1 - klo] - xs[klo];
    if (xext >= m_iso_endpoint_max_xext) return false;
    if (xext >= m_iso_endpoint_xext_frac * length) return false;

    // Principal axis of the qualified points: covariance + power iteration
    // (same numerical pattern as the local-PCA refinement below), seeded with
    // the cluster's global PCA axis 0 so the sign is stable.
    double cx = 0, cy = 0, cz = 0;
    for (int i : qual) {
        const auto p = cluster.point3d(i);
        cx += p.x();  cy += p.y();  cz += p.z();
    }
    const double nn = static_cast<double>(qual.size());
    cx /= nn;  cy /= nn;  cz /= nn;

    double Sxx=0, Sxy=0, Sxz=0, Syy=0, Syz=0, Szz=0;
    for (int i : qual) {
        const auto p = cluster.point3d(i);
        const double dx = p.x() - cx, dy = p.y() - cy, dz = p.z() - cz;
        Sxx += dx*dx;  Sxy += dx*dy;  Sxz += dx*dz;
        Syy += dy*dy;  Syz += dy*dz;  Szz += dz*dz;
    }

    // Power iteration on a 3x3 symmetric matrix, seeded with a global PCA axis
    // so the sign (and the pick order downstream) is stable.
    auto power_axis = [](const double S[6], const Facade::geo_point_t& seed,
                         double& ax, double& ay, double& az) -> bool {
        // S = {Sxx, Sxy, Sxz, Syy, Syz, Szz}
        ax = seed.x();  ay = seed.y();  az = seed.z();
        const double smag = std::sqrt(ax*ax + ay*ay + az*az);
        if (smag < 1e-12) return false;
        ax /= smag;  ay /= smag;  az /= smag;
        for (int iter = 0; iter < 10; ++iter) {
            const double wx = S[0]*ax + S[1]*ay + S[2]*az;
            const double wy = S[1]*ax + S[3]*ay + S[4]*az;
            const double wz = S[2]*ax + S[4]*ay + S[5]*az;
            const double mag = std::sqrt(wx*wx + wy*wy + wz*wz);
            if (mag < 1e-12) break;
            ax = wx/mag;  ay = wy/mag;  az = wz/mag;
        }
        if (ax*seed.x() + ay*seed.y() + az*seed.z() < 0) { ax = -ax; ay = -ay; az = -az; }
        return true;
    };

    const double S0[6] = {Sxx, Sxy, Sxz, Syy, Syz, Szz};
    const auto& seed = cluster.get_pca().axis.at(0);
    double vx = 0, vy = 0, vz = 0;
    if (!power_axis(S0, seed, vx, vy, vz)) return false;

    // Second axis (in-sheet transverse) for the aspect gate: deflate with the
    // Rayleigh quotient lambda0 = v0^T S v0 -- the power iteration normalises
    // the eigenvalue away, so it has to be recovered explicitly.
    double ux = 0, uy = 0, uz = 0;
    {
        const double Sv_x = Sxx*vx + Sxy*vy + Sxz*vz;
        const double Sv_y = Sxy*vx + Syy*vy + Syz*vz;
        const double Sv_z = Sxz*vx + Syz*vy + Szz*vz;
        const double lam0 = vx*Sv_x + vy*Sv_y + vz*Sv_z;
        const double S1[6] = {Sxx - lam0*vx*vx, Sxy - lam0*vx*vy, Sxz - lam0*vx*vz,
                              Syy - lam0*vy*vy, Syz - lam0*vy*vz, Szz - lam0*vz*vz};
        const auto& seed1 = cluster.get_pca().axis.at(1);
        if (!power_axis(S1, seed1, ux, uy, uz)) return false;
    }

    // Axis projections and transverse distances from the axis line through the
    // centroid.  The quantile-trimmed axial extremes s_lo/s_hi serve the GATE
    // only (degeneracy guard + aspect denominator); the endpoints below are
    // untrimmed.
    std::vector<double> ss, ts, dperp;
    ss.reserve(qual.size());
    ts.reserve(qual.size());
    dperp.reserve(qual.size());
    for (int i : qual) {
        const auto p = cluster.point3d(i);
        const double dx = p.x()-cx, dy = p.y()-cy, dz = p.z()-cz;
        const double s = dx*vx + dy*vy + dz*vz;
        ss.push_back(s);
        ts.push_back(dx*ux + dy*uy + dz*uz);
        const double r2 = dx*dx + dy*dy + dz*dz - s*s;
        dperp.push_back(std::sqrt(r2 > 0 ? r2 : 0.0));
    }
    std::vector<double> ss_sorted(ss);
    std::sort(ss_sorted.begin(), ss_sorted.end());
    const double s_lo = ss_sorted[klo];
    const double s_hi = ss_sorted[nq - 1 - klo];
    const double band = 3 * units::cm;
    if (s_hi - s_lo < 3 * band) return false;  // degenerate axis

    // Sheet-aspect gate: trimmed transverse extent over trimmed axial extent.
    // A track-like (1-D) cluster is thin in BOTH transverse directions and is
    // handed back to the legacy path here, so this branch can only act on
    // genuine 2-D sheets (doc pr/24 sec 15).
    std::vector<double> ts_sorted(ts);
    std::sort(ts_sorted.begin(), ts_sorted.end());
    const double ext0 = s_hi - s_lo;
    const double ext1 = ts_sorted[nq - 1 - klo] - ts_sorted[klo];
    const double aspect = ext1 / ext0;
    if (aspect < m_iso_endpoint_min_aspect) {
        SPDLOG_LOGGER_DEBUG(s_log,
            "iso endpoint: rejected by aspect (L={:.1f} cm, xext={:.1f} cm, aspect={:.3f} < {:.3f})",
            length/units::cm, xext/units::cm, aspect, m_iso_endpoint_min_aspect);
        return false;
    }

    // Endpoint = axial extreme FIRST, lateral centring SECOND.
    //
    // A first draft filtered to a tube of radius m_iso_endpoint_tube_radius
    // around the axis line and took the extreme of that subset.  Measured on
    // the 38 gated clusters of the doc pr/24 sec 15 set, that re-introduced the
    // very inward bias this round removes: a straight axis line through the
    // centroid leaves a long or curved object near its tips, so the tube held
    // only 383/2097 points on SBND mcp1k evt 282204 and its extreme sat 28.6 cm
    // INSIDE the round-2 endpoint (evt 284794: 20.8 cm).  The failure mode is
    // curvature over a long lever arm, not width.  So the axial extreme is now
    // taken over ALL qualified points -- the endpoint can never move inward,
    // which is what makes "cannot leave a stub for find_other_segments" true by
    // construction -- and the tube radius survives only as a diagnostic on how
    // laterally central the chosen point is.
    //
    // Lateral centring still matters, and is what keeps the sheet-corner
    // degeneracy from coming back: among the points within `band` of the axial
    // extreme, the one nearest the axis line is taken.  On a sheet's wide end
    // that is mid-width rather than a corner; on the narrow end it is the apex.
    auto support_at = [&](double s_j) -> ptrdiff_t {
        return std::upper_bound(ss_sorted.begin(), ss_sorted.end(), s_j + band)
             - std::lower_bound(ss_sorted.begin(), ss_sorted.end(), s_j - band);
    };

    auto pick_end_point = [&](int sign, double& s_raw, double& walk_in,
                              double& d_perp_pick) -> int {
        std::vector<size_t> order(qual.size());
        for (size_t k = 0; k < qual.size(); ++k) order[k] = k;
        std::sort(order.begin(), order.end(),
                  [&](size_t a, size_t b) {
                      const double sa = sign * ss[a], sb = sign * ss[b];
                      if (sa != sb) return sa > sb;   // most extreme first
                      return a < b;                   // deterministic tie-break
                  });
        s_raw = ss[order.front()];
        // Walk inward only past isolated spikes: the accepted axial position
        // must carry at least three qualified points (itself included) within
        // `band`.  No quantile discard; walk_in is logged and is ~0 on healthy
        // objects.
        double s_ext = 0;
        bool found = false;
        for (size_t j = 0; j < order.size(); ++j) {
            const double s_j = ss[order[j]];
            if (support_at(s_j) < 3) continue;
            s_ext = s_j;
            found = true;
            break;
        }
        if (!found) return -1;
        walk_in = std::abs(s_raw - s_ext);
        // Laterally most central point of the end band (strict < keeps the
        // lowest qualifying index on ties -- deterministic).
        double best_d = std::numeric_limits<double>::max();
        int best = -1;
        for (size_t k = 0; k < qual.size(); ++k) {
            if (std::abs(ss[k] - s_ext) > band) continue;
            if (dperp[k] < best_d) { best_d = dperp[k]; best = qual[k]; }
        }
        d_perp_pick = best_d;
        return best;
    };

    double s_rawA = 0, s_rawB = 0, walk_inA = 0, walk_inB = 0, dpA = 0, dpB = 0;
    const int iA = pick_end_point(-1, s_rawA, walk_inA, dpA);
    const int iB = pick_end_point(+1, s_rawB, walk_inB, dpB);
    if (iA < 0 || iB < 0 || iA == iB) return false;

    const auto pA = cluster.point3d(iA);
    const auto pB = cluster.point3d(iB);

    // Snap to Steiner terminals, same contract as
    // get_two_boundary_steiner_graph_idx: endpoints must be terminals
    // (original data points), not auxiliary Steiner nodes.
    if (!cluster.has_pc("steiner_pc")) return false;
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto xa = steiner_pc.get(coords.at(0));
    const auto ya = steiner_pc.get(coords.at(1));
    const auto za = steiner_pc.get(coords.at(2));
    const auto fa = steiner_pc.get("flag_steiner_terminal");
    if (!xa || !ya || !za || !fa) return false;
    const auto& sx = xa->elements<double>();
    const auto& sy = ya->elements<double>();
    const auto& sz = za->elements<double>();
    const auto& sterm = fa->elements<int>();
    const size_t nst = sx.size();
    if (nst < 2) return false;

    auto nearest_terminal = [&](const Facade::geo_point_t& target) -> int {
        double best_d2 = std::numeric_limits<double>::max();
        int best = -1;
        for (size_t i = 0; i < nst; ++i) {
            if (!sterm[i]) continue;
            const double dx = sx[i]-target.x(), dy = sy[i]-target.y(), dz = sz[i]-target.z();
            const double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < best_d2) { best_d2 = d2; best = static_cast<int>(i); }
        }
        return best;
    };

    const int tA = nearest_terminal(pA);
    const int tB = nearest_terminal(pB);
    if (tA < 0 || tB < 0 || tA == tB) return false;

    p1 = Facade::geo_point_t(sx[tA], sy[tA], sz[tA]);
    p2 = Facade::geo_point_t(sx[tB], sy[tB], sz[tB]);

    // `gain` = how much further out each endpoint now reaches than round 2's
    // quantile-trimmed extreme would have; `walk` = the spike-guard walk-in
    // (expected ~0; a persistent walk > 1-2 cm would mean the guard is
    // re-introducing the inward bias this round removes).
    SPDLOG_LOGGER_DEBUG(s_log,
        "iso endpoint: fired (L={:.1f} cm, xext={:.1f} cm, aspect={:.3f}, n={}) "
        "A=({:.1f},{:.1f},{:.1f}) B=({:.1f},{:.1f},{:.1f}) gain={:.1f}/{:.1f} cm "
        "walk={:.2f}/{:.2f} cm dperp={:.1f}/{:.1f} cm intube={}/{}",
        length/units::cm, xext/units::cm, aspect, qual.size(),
        p1.x()/units::cm, p1.y()/units::cm, p1.z()/units::cm,
        p2.x()/units::cm, p2.y()/units::cm, p2.z()/units::cm,
        (s_lo - s_rawA)/units::cm, (s_rawB - s_hi)/units::cm,
        walk_inA/units::cm, walk_inB/units::cm,
        dpA/units::cm, dpB/units::cm,
        dpA <= m_iso_endpoint_tube_radius, dpB <= m_iso_endpoint_tube_radius);
    return true;
}

SegmentPtr PatternAlgorithms::init_first_segment(Graph& graph, Facade::Cluster& cluster, Facade::Cluster* main_cluster,TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_back_search)
{
    using IFS_Clock = std::chrono::steady_clock;
    using IFS_MS = std::chrono::duration<double, std::milli>;
    auto t0 = IFS_Clock::now();
    // Require a NON-EMPTY steiner_pc BEFORE the boundary call below, which
    // throws on an empty cloud (zero-length arrays; SBND MC evt 11).
    if (!cluster.has_pc("steiner_pc") || cluster.get_pc("steiner_pc").size_major() == 0) return nullptr;
    Facade::geo_point_t first_pt, second_pt;

    // Isochronous sheet endpoints (m_iso_endpoint, C++ default OFF; doc
    // sbnd_xin/docs/pr/24 round 2).  When the gate fires, the whole legacy
    // block below -- boundary metric AND local-PCA refinement, both of which
    // degenerate on a filled 2-D sheet -- is bypassed.
    bool iso_endpoints_used = false;
    if (m_iso_endpoint) {
        iso_endpoints_used = find_iso_first_segment_endpoints(cluster, first_pt, second_pt);
    }
    if (!iso_endpoints_used) {
    // Get two boundary points from the cluster
    auto boundary_indices = cluster.get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc");
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();

    // Add the two boundary points as additional extreme point groups
    Facade::geo_point_t boundary_point_first(x_coords[boundary_indices.first], 
                                y_coords[boundary_indices.first],  
                                z_coords[boundary_indices.first]);
    Facade::geo_point_t boundary_point_second(x_coords[boundary_indices.second], 
                                y_coords[boundary_indices.second], 
                                z_coords[boundary_indices.second]);

    // // hack first and second point
    // boundary_point_first = Facade::geo_point_t(2285.68, 866.2, 6527.5);
    // boundary_point_second = Facade::geo_point_t(2296.69, 872.262, 5818);
    // // end hack
                                
    first_pt = boundary_point_first;
    second_pt = boundary_point_second;

    // Local PCA refinement of each boundary endpoint.
    //
    // The global boundary search returns the extreme Steiner point along a fixed
    // projection axis.  When the cluster curves near its end, the true tip may lie
    // "around a corner" and be missed.
    //
    // Fix: collect all Steiner points within r_local of the candidate endpoint,
    // compute their local PCA direction (first principal component via power
    // iteration), then find whichever of those points lies furthest from the local
    // centroid along that direction.  Using all neighbours equally via PCA avoids
    // the bias that arises from using the candidate endpoint itself to define the
    // direction.  The power iteration is seeded with fallback_dir (the global
    // outward direction) so the sign is consistent and convergence is fast.
    // doc pr/30 §11, P2.  m_first_seg_local_pca defaults TRUE -- unlike the
    // other pr/30 knobs this behaviour is already production, so the DEFAULT
    // is what ships today and OFF is the new option.  The point of the knob is
    // that an unconditional, un-knobbed departure from the prototype in the
    // single most load-bearing function of the stage could not be measured at
    // all; the prototype takes get_two_boundary_wcps(2) and uses it as-is
    // (NeutrinoID_proto_vertex.h:426).
    if (m_first_seg_local_pca) {
        const double r_local = 10.0 * units::cm;

        Facade::geo_vector_t global_dir(
            boundary_point_second.x() - boundary_point_first.x(),
            boundary_point_second.y() - boundary_point_first.y(),
            boundary_point_second.z() - boundary_point_first.z());
        const double global_mag = global_dir.magnitude();

        if (global_mag > 0) {
            global_dir = global_dir * (1.0 / global_mag);  // unit vector first→second

            auto refine_endpoint = [&](const Facade::geo_point_t& pt,
                                       const Facade::geo_vector_t& fallback_dir) -> Facade::geo_point_t {
                // Gather local Steiner neighbourhood and filter to terminal nodes.
                // Terminal nodes are the original cluster data points; intermediate
                // Steiner nodes are auxiliary connectivity points that would pull the
                // centroid and covariance off the true track axis.
                // Fall back to all Steiner points if there are too few terminals.
                auto all_nbrs = cluster.kd_steiner_radius(r_local, pt, "steiner_pc");

                std::vector<size_t> term_idx;
                term_idx.reserve(all_nbrs.size());
                for (auto& [idx, d2] : all_nbrs) {
                    if (flag_steiner_terminal[idx]) term_idx.push_back(idx);
                }
                // Use terminals if we have enough; otherwise fall back to all neighbours
                const bool use_terminals = (term_idx.size() >= 3);
                const std::vector<size_t>* use_idx = &term_idx;
                std::vector<size_t> all_idx;
                if (!use_terminals) {
                    if (all_nbrs.size() < 3) return pt;
                    all_idx.reserve(all_nbrs.size());
                    for (auto& [idx, d2] : all_nbrs) all_idx.push_back(idx);
                    use_idx = &all_idx;
                }

                // Centroid (from terminals, or all points as fallback)
                double cx = 0.0, cy = 0.0, cz = 0.0;
                for (size_t idx : *use_idx) {
                    cx += x_coords[idx]; cy += y_coords[idx]; cz += z_coords[idx];
                }
                const double nn = static_cast<double>(use_idx->size());
                cx /= nn;  cy /= nn;  cz /= nn;

                // Covariance matrix (unnormalised — scaling cancels in power iteration)
                double Sxx=0, Sxy=0, Sxz=0, Syy=0, Syz=0, Szz=0;
                for (size_t idx : *use_idx) {
                    double dx = x_coords[idx] - cx;
                    double dy = y_coords[idx] - cy;
                    double dz = z_coords[idx] - cz;
                    Sxx += dx*dx;  Sxy += dx*dy;  Sxz += dx*dz;
                    Syy += dy*dy;  Syz += dy*dz;  Szz += dz*dz;
                }

                // Power iteration: converges to the first principal component.
                // Seed with fallback_dir so the sign is correct from the start.
                double vx = fallback_dir.x(), vy = fallback_dir.y(), vz = fallback_dir.z();
                for (int iter = 0; iter < 10; ++iter) {
                    double wx = Sxx*vx + Sxy*vy + Sxz*vz;
                    double wy = Sxy*vx + Syy*vy + Syz*vz;
                    double wz = Sxz*vx + Syz*vy + Szz*vz;
                    double mag = std::sqrt(wx*wx + wy*wy + wz*wz);
                    if (mag < 1e-12) break;
                    vx = wx/mag;  vy = wy/mag;  vz = wz/mag;
                }

                // Ensure final direction agrees with fallback (flip if anti-parallel)
                if (vx*fallback_dir.x() + vy*fallback_dir.y() + vz*fallback_dir.z() < 0) {
                    vx = -vx;  vy = -vy;  vz = -vz;
                }

                // std::cout << "Refine endpoint: " << pt << " centroid(" << cx << "," << cy << "," << cz << ") local_dir(" << vx << "," << vy << "," << vz << ") fallback_dir(" << fallback_dir.x() << "," << fallback_dir.y() << "," << fallback_dir.z() << ") use_terminals=" << use_terminals << std::endl;

                // Find the terminal neighbour with the greatest projection from centroid
                // along local_dir.  Search only among terminals (or all points as fallback)
                // so the returned endpoint is a physical cluster point.
                // Initialise with pt so we only move if a genuinely further point exists.
                Facade::geo_point_t best = pt;
                double best_proj = (pt.x()-cx)*vx + (pt.y()-cy)*vy + (pt.z()-cz)*vz;
                for (size_t idx : *use_idx) {
                    Facade::geo_point_t nb(x_coords[idx], y_coords[idx], z_coords[idx]);
                    double proj = (nb.x()-cx)*vx + (nb.y()-cy)*vy + (nb.z()-cz)*vz;
                    if (proj > best_proj) { best_proj = proj; best = nb; }
                }
                return best;
            };

            Facade::geo_vector_t neg_global_dir(
                -global_dir.x(), -global_dir.y(), -global_dir.z());
            first_pt  = refine_endpoint(boundary_point_first,  neg_global_dir);
            second_pt = refine_endpoint(boundary_point_second, global_dir);

            // doc pr/30 §11, P2 instrumentation.  Unconditional and log-only:
            // this runs in the knob-ON (= production) arm too, which is the
            // arm that answers "how much does the refinement actually move
            // the seed of the whole stage".
            for (const auto& mv : {std::make_pair(boundary_point_first,  first_pt),
                                   std::make_pair(boundary_point_second, second_pt)}) {
                const double d = ray_length(Ray{mv.first, mv.second});
                g_port_audit.pca_refine_calls.fetch_add(1, std::memory_order_relaxed);
                if (d > 0) {
                    const uint64_t um = static_cast<uint64_t>(d / units::um);
                    g_port_audit.pca_refine_moved.fetch_add(1, std::memory_order_relaxed);
                    g_port_audit.pca_move_um_sum.fetch_add(um, std::memory_order_relaxed);
                    uint64_t prev = g_port_audit.pca_move_um_max.load(std::memory_order_relaxed);
                    while (um > prev &&
                           !g_port_audit.pca_move_um_max.compare_exchange_weak(prev, um,
                                                                std::memory_order_relaxed)) {}
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr30 P2 local-PCA endpoint moved {:.3f} cm: ({:.2f},{:.2f},{:.2f}) -> "
                        "({:.2f},{:.2f},{:.2f}) cluster {}",
                        d/units::cm,
                        mv.first.x()/units::cm, mv.first.y()/units::cm, mv.first.z()/units::cm,
                        mv.second.x()/units::cm, mv.second.y()/units::cm, mv.second.z()/units::cm,
                        cluster.get_cluster_id());
                }
            }

            // std::cout << boundary_point_first << " -> " << first_pt << std::endl;
            // std::cout << boundary_point_second << " -> " << second_pt << std::endl;
        }
    }
    } // !iso_endpoints_used: legacy boundary endpoints + local-PCA refinement

    // Determine the starting point based on whether this is the main cluster or not
    const bool is_main_flag = cluster.get_flag(Facade::Flags::main_cluster);
    const bool is_main_ptr  = (&cluster == main_cluster);
    if (is_main_flag != is_main_ptr) {
        SPDLOG_LOGGER_WARN(s_log, "init_first_segment: main_cluster flag ({}) disagrees with pointer comparison ({}); "
                           "using pointer comparison as authoritative", is_main_flag, is_main_ptr);
    }
    if (is_main_ptr) {
        // Main cluster: start from downstream (or upstream if flag_back_search)
        if (flag_back_search) {
            // Start from high z (upstream/backward)
            if (first_pt.z() < second_pt.z()) {
                std::swap(first_pt, second_pt);
            }
        } else {
            // Start from low z (downstream/forward)
            if (first_pt.z() > second_pt.z()) {
                std::swap(first_pt, second_pt);
            }
        }
        SPDLOG_LOGGER_TRACE(s_log, "init_first_segment: main cluster, flag_back_search={} -> start({:.2f},{:.2f},{:.2f})",
                            flag_back_search, first_pt.x(), first_pt.y(), first_pt.z());
    } else if (main_cluster) {
        // Non-main cluster: start from the point closest to main cluster
        // Find closest distances to main cluster's Steiner point cloud
        auto knn1 = main_cluster->kd_steiner_knn(1, first_pt, "steiner_pc");
        auto knn2 = main_cluster->kd_steiner_knn(1, second_pt, "steiner_pc");
        
        if (!knn1.empty() && !knn2.empty()) {
            double dis1 = std::sqrt(knn1[0].second);
            double dis2 = std::sqrt(knn2[0].second);
            
            // Start from the point closer to main cluster
            if (dis2 < dis1) {
                std::swap(first_pt, second_pt);
            }
            SPDLOG_LOGGER_TRACE(s_log, "init_first_segment: non-main cluster, dis_A={:.2f} dis_B={:.2f} -> start({:.2f},{:.2f},{:.2f})",
                                dis1, dis2, first_pt.x(), first_pt.y(), first_pt.z());
        }
    } else {
        // main_cluster is nullptr and this is not the main cluster:
        // ordering relative to main is undefined; apply ascending-z as a deterministic fallback
        SPDLOG_LOGGER_TRACE(s_log, "init_first_segment: main_cluster is nullptr for a non-main cluster; "
                            "falling back to ascending-z order");
        if (first_pt.z() > second_pt.z()) {
            std::swap(first_pt, second_pt);
        }
    }
    
    SPDLOG_LOGGER_TRACE(s_log, "init_first_segment: raw boundary pts  A({:.2f},{:.2f},{:.2f})  B({:.2f},{:.2f},{:.2f})  flag_back_search={}",
                        first_pt.x(), first_pt.y(), first_pt.z(),
                        second_pt.x(), second_pt.y(), second_pt.z(),
                        flag_back_search);

    // Create vertices for the endpoints
    VertexPtr v1 = make_vertex(graph);
    v1->wcpt().point = first_pt;
    v1->cluster(&cluster);
    VertexPtr v2 = make_vertex(graph);
    v2->wcpt().point = second_pt;
    v2->cluster(&cluster);

    auto seg = create_segment_from_vertices(graph, cluster, v1, v2, dv);
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "init_first_segment timing: do shortest path took {} ms", IFS_MS(IFS_Clock::now() - t0).count());
    t0 = IFS_Clock::now();

    if (!seg) {
        remove_vertex(graph, v1);
        remove_vertex(graph, v2);
        return nullptr;
    }
    SPDLOG_LOGGER_TRACE(s_log, "init_first_segment: Dijkstra path  npts={}", seg->wcpts().size());
    // for (size_t i = 0; i < seg->wcpts().size(); ++i) {
    //     SPDLOG_LOGGER_TRACE(s_log, "  [{}] ({:.2f},{:.2f},{:.2f})", i, 
    //                         seg->wcpts()[i].point.x(), seg->wcpts()[i].point.y(), seg->wcpts()[i].point.z());
    // }


    // perform fitting ...
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "init_first_segment timing: create segment and prepare data took {} ms", IFS_MS(IFS_Clock::now() - t0).count());
    t0 = IFS_Clock::now();
    track_fitter.add_segment(seg);
    track_fitter.do_single_tracking(seg, true, true, false, false, &cluster);
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "init_first_segment timing: do single_track fitting took {} ms", IFS_MS(IFS_Clock::now() - t0).count());
    t0 = IFS_Clock::now();

    const auto& fine_path = track_fitter.get_fine_tracking_path();
    const auto& dQ_vec = track_fitter.get_dQ();
    const auto& dx_vec = track_fitter.get_dx();
    const auto& pu_vec = track_fitter.get_pu();
    const auto& pv_vec = track_fitter.get_pv();
    const auto& pw_vec = track_fitter.get_pw();
    const auto& pt_vec = track_fitter.get_pt();
    const auto& chi2_vec = track_fitter.get_reduced_chi2();
    SPDLOG_LOGGER_TRACE(s_log, "init_first_segment: fitted path     npts={}", fine_path.size());
    // for (size_t i = 0; i < fine_path.size(); ++i) {
    //     SPDLOG_LOGGER_TRACE(s_log, "  [{}] ({:.2f},{:.2f},{:.2f})", i, 
    //                         fine_path[i].first.x(), fine_path[i].first.y(), fine_path[i].first.z());
    // }

    if (fine_path.size()>1) {
        v1->fit().point = fine_path.front().first;
        if (!dQ_vec.empty()) v1->fit().dQ = dQ_vec.front();
        if (!dx_vec.empty()) v1->fit().dx = dx_vec.front();
        if (!pu_vec.empty()) v1->fit().pu = pu_vec.front();
        if (!pv_vec.empty()) v1->fit().pv = pv_vec.front();
        if (!pw_vec.empty()) v1->fit().pw = pw_vec.front();
        if (!pt_vec.empty()) v1->fit().pt = pt_vec.front();
        if (!chi2_vec.empty()) v1->fit().reduced_chi2 = chi2_vec.front();
        
        v2->fit().point = fine_path.back().first;
        if (!dQ_vec.empty()) v2->fit().dQ = dQ_vec.back();
        if (!dx_vec.empty()) v2->fit().dx = dx_vec.back();
        if (!pu_vec.empty()) v2->fit().pu = pu_vec.back();
        if (!pv_vec.empty()) v2->fit().pv = pv_vec.back();
        if (!pw_vec.empty()) v2->fit().pw = pw_vec.back();
        if (!pt_vec.empty()) v2->fit().pt = pt_vec.back();
        if (!chi2_vec.empty()) v2->fit().reduced_chi2 = chi2_vec.back();

        // Set fit information for segment
        // std::vector<Fit> fits;
        // for (size_t i = 0; i < fine_path.size(); ++i) {
        //     Fit fit;
        //     fit.point = fine_path[i].first;
        //     if (i < dQ_vec.size()) fit.dQ = dQ_vec[i];
        //     if (i < dx_vec.size()) fit.dx = dx_vec[i];
        //     if (i < pu_vec.size()) fit.pu = pu_vec[i];
        //     if (i < pv_vec.size()) fit.pv = pv_vec[i];
        //     if (i < pw_vec.size()) fit.pw = pw_vec[i];
        //     if (i < pt_vec.size()) fit.pt = pt_vec[i];
        //     if (i < chi2_vec.size()) fit.reduced_chi2 = chi2_vec[i];
        //     fits.push_back(fit);
        // }
        // seg->fits(fits);
    }else{
        // Tracking failed, clean up
        remove_segment(graph, seg);
        remove_vertex(graph, v1);
        remove_vertex(graph, v2);
        return nullptr;
    }
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "init_first_segment timing: after fit assignment took {} ms", IFS_MS(IFS_Clock::now() - t0).count());
    
    return seg;
}


std::pair<Facade::geo_point_t,  size_t> PatternAlgorithms::proto_extend_point(const Facade::Cluster& cluster, Facade::geo_point_t& p, Facade::geo_vector_t& dir, Facade::geo_vector_t& dir_other, bool flag_continue, std::vector<Facade::geo_point_t>* walk_history){
    const double step_dis = 1.0 * units::cm;

    if (!cluster.has_pc("steiner_pc")) return {p, 0};

    // Get steiner point cloud data
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& steiner_x = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& steiner_y = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& steiner_z = steiner_pc.get(coords.at(2))->elements<double>();
    
    // Find closest point in steiner point cloud
    auto curr_knn_results = cluster.kd_steiner_knn(1, p, "steiner_pc");
    if (curr_knn_results.empty()) {
        return std::make_pair(p, INVALID_STEINER_INDEX); // No steiner point found
    }

    size_t curr_index = curr_knn_results[0].first;
    Facade::geo_point_t curr_wcp(steiner_x[curr_index], steiner_y[curr_index], steiner_z[curr_index]);
    Facade::geo_point_t next_wcp = curr_wcp;

    // (saved_start_wcp and saved_dir removed: set but never used)
    if (walk_history) walk_history->push_back(curr_wcp);

    // Forward search.  Cycle backstop: each pass moves to a *different*
    // point (mag != 0 is required) but nothing forbids a cycle on looping
    // point-cloud topology, and this walk had no bound (its relative
    // Facade_Cluster.cxx:343 caps at 400).  NOT 400 here: this walk can
    // legitimately follow a track for its full length in ~1 cm steps, and
    // the same cap on the TaggerCheckSTM copy of this walk was measured to
    // truncate real crawls.  4000 is beyond any SBND geometry (~6.4 m
    // diagonal = ~640 steps) while a genuine cycle still dies in
    // milliseconds.
    int walk_counter = 0;
    while(flag_continue && walk_counter < 4000){
        flag_continue = false;
        ++walk_counter;

        for (int i = 0; i != 3; i++){
            Facade::geo_point_t test_p(
                curr_wcp.x() + dir.x() * step_dis * (i + 1),
                curr_wcp.y() + dir.y() * step_dis * (i + 1),
                curr_wcp.z() + dir.z() * step_dis * (i + 1)
            );

            // Try steiner point cloud first
            auto next_knn_steiner = cluster.kd_steiner_knn(1, test_p, "steiner_pc");
            if (!next_knn_steiner.empty()) {
                size_t next_index = next_knn_steiner[0].first;
                next_wcp = Facade::geo_point_t(steiner_x[next_index], steiner_y[next_index], steiner_z[next_index]);
                Facade::geo_vector_t dir2(
                    next_wcp.x() - curr_wcp.x(),
                    next_wcp.y() - curr_wcp.y(),
                    next_wcp.z() - curr_wcp.z()
                );

                double mag2 = dir2.magnitude();
                if (mag2 != 0) {
                    double angle = std::acos(dir2.dot(dir) / mag2) / 3.1415926 * 180.0;
                    if (angle < 25.0) {
                        flag_continue = true;
                        curr_wcp = next_wcp;
                        curr_index = next_index;
                        dir = dir2.norm() + dir * 5;
                        dir = dir / dir.magnitude();
                        if (walk_history) walk_history->push_back(curr_wcp);
                        break;
                    }
                }
            }

            // Try regular point cloud
            auto closest_result = cluster.get_closest_wcpoint(test_p);
            // size_t regular_index = closest_result.first;
            next_wcp = closest_result.second;

            Facade::geo_vector_t dir1(
                next_wcp.x() - curr_wcp.x(),
                next_wcp.y() - curr_wcp.y(),
                next_wcp.z() - curr_wcp.z()
            );

            double mag1 = dir1.magnitude();
            if (mag1 != 0) {
                double angle = std::acos(dir1.dot(dir) / mag1) / 3.1415926 * 180.0;
                if (angle < 17.5) {
                    auto updated_knn = cluster.kd_steiner_knn(1, next_wcp, "steiner_pc");
                    if (!updated_knn.empty()) {
                        flag_continue = true;
                        curr_wcp = next_wcp;
                        curr_index = updated_knn[0].first;
                        dir = dir1.norm() + dir * 5;  // both terms are now dimensionless
                        dir = dir / dir.magnitude();
                        if (walk_history) walk_history->push_back(curr_wcp);
                        break;
                    }
                    // Secondary steiner KNN failed — do not advance, continue loop
                }
            }
        }
    }
    
    if (walk_counter >= 4000) {
        SPDLOG_LOGGER_WARN(s_log, "proto_extend_point: walk did not settle after 4000 passes; stopping to avoid an unbounded loop (cluster {})", cluster.get_cluster_id());
    }

    // Ensure we return the steiner point cloud position
    Facade::geo_point_t test_p(curr_wcp.x(), curr_wcp.y(), curr_wcp.z());
    auto final_knn = cluster.kd_steiner_knn(1, test_p, "steiner_pc");
    if (!final_knn.empty()) {
        curr_index = final_knn[0].first;
        curr_wcp = Facade::geo_point_t(steiner_x[curr_index], steiner_y[curr_index], steiner_z[curr_index]);
    }
    
    // Return: point, (cloud_type=2 for steiner, point_index)
    return std::make_pair(curr_wcp,  curr_index);
}

bool PatternAlgorithms::proto_break_tracks(const Facade::Cluster& cluster, const Facade::geo_point_t& first_wcp, Facade::geo_point_t& curr_wcp, const Facade::geo_point_t& last_wcp, std::list<Facade::geo_point_t>& wcps_list1, std::list<Facade::geo_point_t>& wcps_list2, bool flag_pass_check){
    
    // Calculate distances
    double dis1 = std::sqrt(std::pow(curr_wcp.x() - first_wcp.x(), 2) + 
                            std::pow(curr_wcp.y() - first_wcp.y(), 2) + 
                            std::pow(curr_wcp.z() - first_wcp.z(), 2));
    double dis2 = std::sqrt(std::pow(curr_wcp.x() - last_wcp.x(), 2) + 
                            std::pow(curr_wcp.y() - last_wcp.y(), 2) + 
                            std::pow(curr_wcp.z() - last_wcp.z(), 2));
    
    // Check if distances are sufficient or if we should pass the check
    if ((dis1 > 1.0 * units::cm && dis2 > 1.0 * units::cm) || flag_pass_check) {
        // Find shortest path from first_wcp to curr_wcp using steiner graph
        Facade::geo_point_t first_point = first_wcp;
        Facade::geo_point_t curr_point = curr_wcp;
        auto path1 = do_rough_path(cluster, first_point, curr_point);
        // Convert vector to list
        wcps_list1.clear();
        for (const auto& pt : path1) {
            wcps_list1.push_back(pt);
        }
        
        // Find shortest path from curr_wcp to last_wcp using steiner graph
        Facade::geo_point_t curr_point2 = curr_wcp;
        Facade::geo_point_t last_point = last_wcp;
        auto path2 = do_rough_path(cluster, curr_point2, last_point);
        // Convert vector to list
        wcps_list2.clear();
        for (const auto& pt : path2) {
            wcps_list2.push_back(pt);
        }
        
        // Remove overlapping points at the junction
        // Count how many points overlap from the end of list1 and beginning of list2
        int count = 0;
        if (!wcps_list1.empty() && !wcps_list2.empty()) {
            // Compare points from the end of list1 with the beginning of list2
            // Use reverse iterator for list1 and forward iterator for list2
            auto it1 = wcps_list1.rbegin();  // Start from the back of list1
            auto it2 = wcps_list2.begin();   // Start from the front of list2
            
            while (it1 != wcps_list1.rend() && it2 != wcps_list2.end()) {
                // Check if points are the same (within tolerance)
                double dx = it1->x() - it2->x();
                double dy = it1->y() - it2->y();
                double dz = it1->z() - it2->z();
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist < 0.01 * units::cm) {  // same point
                    count++;
                    ++it1;
                    ++it2;
                } else {
                    break;  // no more overlapping points
                }
            }
            
            // Remove overlapping points (keep one copy at the junction)
            for (int i = 0; i < count; i++) {
                if (i + 1 != count) {  // Keep the last overlapping point
                    if (!wcps_list1.empty()) wcps_list1.pop_back();
                    if (!wcps_list2.empty()) wcps_list2.pop_front();
                }
            }
        }
        
        // Snap curr_wcp to the steiner-graph-resolved junction point
        if (!wcps_list1.empty()) {
            curr_wcp = wcps_list1.back();
        }
        
        // Check if we have valid paths
        if (wcps_list1.size() <= 1 || wcps_list2.size() <= 1) {
            return false;
        }
        
        return true;
    } else {
        return false;
    }
}

bool PatternAlgorithms::replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr old_vertex, VertexPtr new_vertex, IDetectorVolumes::pointer dv){
    // Get the cluster from the old segment
    auto cluster = seg->cluster();
    if (!cluster) {
        return false;
    }
    
    // Get the other vertex connected to this segment (the one we'll keep)
    VertexPtr other_vertex = find_other_vertex(graph, seg, old_vertex);
    if (!other_vertex) {
        return false;
    }
    
    // Create new segment with the path points
    SegmentPtr new_seg = create_segment_from_vertices(graph, *cluster, other_vertex, new_vertex, dv);
    if (!new_seg) {
        return false;
    }
    
    // Remove the old segment (this will disconnect it from the graph)
    remove_segment(graph, seg);
    
    // Remove the old vertex if it no longer has any connected segments
    if (old_vertex->descriptor_valid()) {
        auto vd = old_vertex->get_descriptor();
        if (boost::degree(vd, graph) == 0) {
            remove_vertex(graph, old_vertex);
        }
    }
        
    // Update the output parameter
    seg = new_seg;
    
    return true;
}



bool PatternAlgorithms::replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr& vtx, std::list<Facade::geo_point_t>& path_point_list, Facade::geo_point_t& break_point, IDetectorVolumes::pointer dv) {
    // // Check that the vertex is only connected to one segment
    // if (!vtx->descriptor_valid()) {
    //     return false;
    // }
    // auto vd = vtx->get_descriptor();
    // if (boost::degree(vd, graph) != 1) {
    //     return false;  // Vertex is connected to more than one segment, cannot replace
    // }
    
    // Get the cluster from the old segment
    auto cluster = seg->cluster();
    if (!cluster) {
        return false;
    }
    
    // Get the other vertex connected to this segment (the one we'll keep)
    VertexPtr other_vertex = find_other_vertex(graph, seg, vtx);
    if (!other_vertex) {
        return false;
    }
    
    // Create new vertex at the break point
    VertexPtr new_vtx = make_vertex(graph);
    new_vtx->wcpt().point = break_point;
    new_vtx->cluster(cluster);
    
    // Convert list to vector for create_segment_for_cluster
    std::vector<Facade::geo_point_t> path_points;
    for (const auto& pt : path_point_list) {
        path_points.push_back(pt);
    }
    
    // Check if path has enough points
    if (path_points.size() <= 1) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Create new segment with the path points
    SegmentPtr new_seg = create_segment_for_cluster(*cluster, dv, path_points, seg->dirsign());
    if (!new_seg) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Remove the old segment (this will disconnect it from the graph)
    remove_segment(graph, seg);

    // Remove the old vertex, if it no longer has any connected segments
    if (vtx->descriptor_valid()) {
        auto vd = vtx->get_descriptor();
        if (boost::degree(vd, graph) == 0) remove_vertex(graph, vtx);
    }
    
    // Add the new segment connecting other_vertex and new_vtx
    add_segment(graph, new_seg, other_vertex, new_vtx);
    
    // Update the output parameters
    seg = new_seg;
    vtx = new_vtx;
    
    return true;
}

 bool PatternAlgorithms::break_segment_into_two(Graph& graph, VertexPtr vtx1, SegmentPtr seg, VertexPtr vtx2, std::list<Facade::geo_point_t>& path_point_list1, Facade::geo_point_t& break_point, std::list<Facade::geo_point_t>& path_point_list2, IDetectorVolumes::pointer dv, SegmentPtr& out_seg2){
    // Get the cluster from the old segment
    auto cluster = seg->cluster();
    if (!cluster) {
        return false;
    }
    
    // Verify that vtx1 and vtx2 are the endpoints of seg
    auto [v1, v2] = find_vertices(graph, seg);
    if ((v1 != vtx1 || v2 != vtx2) && (v1 != vtx2 || v2 != vtx1)) {
        return false;  // The provided vertices don't match the segment endpoints
    }
    
    // Create new vertex at the break point
    VertexPtr new_vtx = make_vertex(graph);
    new_vtx->wcpt().point = break_point;
    new_vtx->cluster(cluster);
    
    // Convert lists to vectors for create_segment_for_cluster
    std::vector<Facade::geo_point_t> path_points1;
    for (const auto& pt : path_point_list1) {
        path_points1.push_back(pt);
    }
    
    std::vector<Facade::geo_point_t> path_points2;
    for (const auto& pt : path_point_list2) {
        path_points2.push_back(pt);
    }
    
    // Check if paths have enough points
    if (path_points1.size() <= 1 || path_points2.size() <= 1) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Create first new segment with path_points1
    SegmentPtr new_seg1 = create_segment_for_cluster(*cluster, dv, path_points1, seg->dirsign());
    if (!new_seg1) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Create second new segment with path_points2
    SegmentPtr new_seg2 = create_segment_for_cluster(*cluster, dv, path_points2, seg->dirsign());
    if (!new_seg2) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Remove the old segment
    remove_segment(graph, seg);
    
    // Add the first new segment connecting vtx1 and new_vtx
    add_segment(graph, new_seg1, vtx1, new_vtx);
    
    // Add the second new segment connecting new_vtx and vtx2
    add_segment(graph, new_seg2, new_vtx, vtx2);
    
    out_seg2 = new_seg2;
    return true;
 }

 bool PatternAlgorithms::break_segments(Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, std::vector<SegmentPtr>& remaining_segments, float dis_cut) {
    bool flag_modified = false;
    int count = 0;
    std::set<size_t> saved_break_wcp_indices;
    using BS_Clock = std::chrono::steady_clock;
    using BS_MS = std::chrono::duration<double, std::milli>;
    BS_MS t_segment_search_kink{0}, t_proto_extend_point{0}, t_proto_break_tracks{0};
    BS_MS t_replace_segment_and_vertex{0}, t_break_segment_into_two{0}, t_do_multi_tracking{0};
    
    // Walk-halfway guard state: track dir1 from the previous break iteration.
    // Persists across outer-while pops so that break #2 can compare against break #1's dir1, etc.
    Facade::geo_vector_t dir1_prev(0, 0, 0);
    bool has_dir1_prev = false;

    // Capture the cluster being worked on before remaining_segments is
    // consumed: the post-process merge below must NOT pick a cluster off an
    // arbitrary (pointer-ordered) edge of the SHARED multi-cluster graph.
    Facade::Cluster* entry_cluster =
        remaining_segments.empty() ? nullptr : remaining_segments.back()->cluster();

    while(!remaining_segments.empty() && count < 2) {
        SegmentPtr curr_sg = remaining_segments.back();
        auto cluster = curr_sg->cluster();
        remaining_segments.pop_back();
        
        // Get the two vertices of this segment
        auto [start_v, end_v] = find_vertices(graph, curr_sg);
        if (!start_v || !end_v) {
            continue;
        }
        
        // Check if vertices match the segment endpoints
        const auto& wcpts = curr_sg->wcpts();
        if (wcpts.size() < 2) continue;
        
        auto front_pt = wcpts.front().point;
        auto back_pt = wcpts.back().point;
        
        // Determine which vertex is start and which is end based on point positions
        double dis_sv_front = ray_length(Ray{start_v->wcpt().point, front_pt});
        double dis_sv_back = ray_length(Ray{start_v->wcpt().point, back_pt});
        
        if (dis_sv_front > dis_sv_back) {
            std::swap(start_v, end_v);
        }
        
        // Initialize the start test point
        Facade::geo_point_t break_wcp = start_v->wcpt().point;
        const auto& point_vec = curr_sg->fits();  // use fit points, matching prototype's get_point_vec()
        Facade::geo_point_t test_start_p = point_vec.front().point;
        
        if (dis_cut > 0) {
            for (size_t i = 0; i < point_vec.size(); ++i) {
                double dis = ray_length(Ray{point_vec[i].point, point_vec.front().point});
                if (dis > dis_cut) {
                    test_start_p = point_vec[i].point;
                    break;
                }
            }
        }
        
        // Search for kinks and extend the break point.
        // kink_pass_counter: the visited-set progress guard below only applies
        // when break_idx is a valid steiner index; the INVALID_STEINER_INDEX
        // outcome (empty steiner kNN, or the no-steiner-pc {p, 0}-style return)
        // bypasses it, and the loop then re-runs segment_search_kink from the
        // same test_start_p with identical inputs -- a deterministic fixed
        // point, i.e. an infinite loop (same silent-no-op shape as the fixed
        // class-A hang, doc pr/11 sec 6.1).  Exact stationarity is detected
        // below; the pass cap is a generous backstop for a valid-index cycle
        // through the seen-before re-search branch.
        int kink_pass_counter = 0;
        bool have_prev_kink_pass = false;
        Facade::geo_point_t prev_pass_test_start_p, prev_pass_break_wcp;
        size_t prev_pass_break_idx = INVALID_STEINER_INDEX;
        // doc pr/48, 59335 fix (b): which accept criterion produced the
        // CURRENT break_wcp (-1 none, 0 A0 wide-cathode, 1-4 C1-C4).
        // accept_crit_cur is written by every segment_search_kink call;
        // accept_crit_break is copied from it only when that call actually
        // advanced break_wcp, so a final no-kink pass (dir1==0, crit -1)
        // cannot erase the criterion of the break that stands.  Only read
        // when m_kink_break_protect; the nullptr out-param when the knob is
        // off is byte-identical.
        int accept_crit_cur = -1;
        int accept_crit_break = -1;
        int* accept_crit_ptr = m_kink_break_protect ? &accept_crit_cur : nullptr;
        while(ray_length(Ray{start_v->wcpt().point, break_wcp}) <= 1.0 * units::cm &&
              ray_length(Ray{end_v->wcpt().point, break_wcp}) > 1.0 * units::cm) {
            if (++kink_pass_counter > 1000) {
                SPDLOG_LOGGER_WARN(s_log, "break_segments: kink loop passed 1000 iterations without leaving the start vertex (cluster {}); stopping to avoid an unbounded loop", cluster->get_cluster_id());
                break;
            }
            
            auto t_op = BS_Clock::now();
            auto kink_tuple = segment_search_kink(curr_sg, test_start_p, "fit", m_mip_dqdx_median, m_cathode_x, m_cathode_kink_xcut, m_cathode_wide_kink_angle, m_cathode_wide_kink_skirt, m_cathode_wide_kink_baseline, m_kink_walk_dqdx_stop, accept_crit_ptr, m_kink_dqdx_hot_ratio);
            t_segment_search_kink += BS_MS(BS_Clock::now() - t_op);
            auto& [kink_point, dir1, dir2, flag_continue] = kink_tuple;

            if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "break_segments 0: cluster={} test_start_p=({:.2f},{:.2f},{:.2f}) kink=({:.2f},{:.2f},{:.2f}) dir1=({:.2f},{:.2f},{:.2f}) dir2=({:.2f},{:.2f},{:.2f}) flag_continue={}",
                cluster->get_cluster_id(),
                test_start_p.x(), test_start_p.y(), test_start_p.z(),
                kink_point.x(), kink_point.y(), kink_point.z(),
                dir1.x(), dir1.y(), dir1.z(),
                dir2.x(), dir2.y(), dir2.z(),
                flag_continue);


            if (dir1.magnitude() != 0) {
                // Find the extreme point
                Facade::geo_vector_t dir1_geo(dir1.x(), dir1.y(), dir1.z());
                Facade::geo_vector_t dir2_geo(dir2.x(), dir2.y(), dir2.z());
                Facade::geo_point_t kink_geo(kink_point.x(), kink_point.y(), kink_point.z());

                t_op = BS_Clock::now();
                std::vector<Facade::geo_point_t> walk_hist;
                auto [break_pt, break_idx] = proto_extend_point(*cluster, kink_geo, dir1_geo, dir2_geo, flag_continue, &walk_hist);
                t_proto_extend_point += BS_MS(BS_Clock::now() - t_op);

                // Walk-halfway guard: if this walk's dir1 is anti-parallel to the previous
                // iteration's dir1, the walker followed a near-endpoint wiggle rather than the
                // segment's macro continuation. Replace break_pt with the walk-history midpoint
                // — an actual accepted steiner/wcp vertex — to shorten the overshoot.
                if (has_dir1_prev && walk_hist.size() >= 3) {
                    double dot_rev = dir1_geo.dot(dir1_prev);
                    if (dot_rev < -0.5) {
                        Facade::geo_point_t half_pt = walk_hist[walk_hist.size() / 2];
                        if (m_perf) SPDLOG_LOGGER_TRACE(s_log,
                            "break_segments walk_halfway: cluster={} dot(dir1,dir1_prev)={:.3f} "
                            "walk_size={} replacing break_pt ({:.2f},{:.2f},{:.2f}) with halfway ({:.2f},{:.2f},{:.2f})",
                            cluster->get_cluster_id(), dot_rev, walk_hist.size(),
                            break_pt.x(), break_pt.y(), break_pt.z(),
                            half_pt.x(), half_pt.y(), half_pt.z());
                        break_pt = half_pt;
                        auto half_knn = cluster->kd_steiner_knn(1, half_pt, "steiner_pc");
                        if (!half_knn.empty()) break_idx = half_knn[0].first;
                        else break_idx = INVALID_STEINER_INDEX;
                    }
                }
                // Save this iteration's dir1 for the next iteration's check.
                dir1_prev = dir1_geo;
                has_dir1_prev = true;

                break_wcp = break_pt;
                accept_crit_break = accept_crit_cur;   // this search's accept stands (pr/48 fix b)

                if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "break_segments 1: cluster={} kink=({:.2f},{:.2f},{:.2f}) break_pt=({:.2f},{:.2f},{:.2f}) break_idx={} flag_continue={}",
                    cluster->get_cluster_id(),
                    kink_geo.x(), kink_geo.y(), kink_geo.z(),
                    break_pt.x(), break_pt.y(), break_pt.z(),
                    break_idx, flag_continue);

                
                // Check if we've seen this break point before
                if (break_idx != INVALID_STEINER_INDEX &&
                    saved_break_wcp_indices.find(break_idx) != saved_break_wcp_indices.end()) {
                    test_start_p = kink_geo;
                    t_op = BS_Clock::now();
                    kink_tuple = segment_search_kink(curr_sg, test_start_p, "fit", m_mip_dqdx_median, m_cathode_x, m_cathode_kink_xcut, m_cathode_wide_kink_angle, m_cathode_wide_kink_skirt, m_cathode_wide_kink_baseline, m_kink_walk_dqdx_stop, accept_crit_ptr, m_kink_dqdx_hot_ratio);
                    t_segment_search_kink += BS_MS(BS_Clock::now() - t_op);
                    auto& [kink_point2, dir1_2, dir2_2, flag_continue2] = kink_tuple;
                    Facade::geo_vector_t dir1_geo2(dir1_2.x(), dir1_2.y(), dir1_2.z());
                    Facade::geo_vector_t dir2_geo2(dir2_2.x(), dir2_2.y(), dir2_2.z());
                    Facade::geo_point_t kink_geo2(kink_point2.x(), kink_point2.y(), kink_point2.z());
                    t_op = BS_Clock::now();
                    auto [break_pt2, break_idx2] = proto_extend_point(*cluster, kink_geo2, dir1_geo2, dir2_geo2, flag_continue2);
                    t_proto_extend_point += BS_MS(BS_Clock::now() - t_op);
                    break_wcp = break_pt2;
                    break_idx = break_idx2;
                    accept_crit_break = accept_crit_cur;   // re-search's accept stands (pr/48 fix b)
                } else {
                    if (break_idx != INVALID_STEINER_INDEX) {
                        saved_break_wcp_indices.insert(break_idx);
                    }
                }
                
                if (ray_length(Ray{start_v->wcpt().point, break_wcp}) <= 1.0 * units::cm &&
                    ray_length(Ray{end_v->wcpt().point, break_wcp}) > 1.0 * units::cm) {
                    test_start_p = kink_geo;
                }

                // Exact stationarity detector.  A pass is a function of
                // (test_start_p, dir1_prev, saved_break_wcp_indices) only --
                // segment_search_kink and proto_extend_point are pure in their
                // arguments.  If two consecutive passes see the same
                // test_start_p AND produce the same (break_wcp, break_idx),
                // then the next pass's inputs are also identical (dir1_prev is
                // at its fixed point and the visited set gained nothing new
                // that changed the outcome): the loop is in a period-1 cycle
                // and will never terminate.  This is the same
                // silent-no-op-turned-infinite-loop shape as the fixed class-A
                // hang (doc pr/11 sec 6.1); the INVALID_STEINER_INDEX outcome
                // (empty steiner kNN) bypasses the visited-set guard entirely
                // and is the easiest way to enter the cycle.
                const bool loop_would_continue =
                    ray_length(Ray{start_v->wcpt().point, break_wcp}) <= 1.0 * units::cm &&
                    ray_length(Ray{end_v->wcpt().point, break_wcp}) > 1.0 * units::cm;
                if (loop_would_continue && have_prev_kink_pass &&
                    prev_pass_test_start_p == test_start_p &&
                    prev_pass_break_wcp == break_wcp &&
                    prev_pass_break_idx == break_idx) {
                    SPDLOG_LOGGER_WARN(s_log, "break_segments: kink pass repeated exactly (break_idx={} valid={}) without leaving the start vertex (cluster {}); stopping to avoid an unbounded loop", break_idx == INVALID_STEINER_INDEX ? -1 : (long long)break_idx, break_idx != INVALID_STEINER_INDEX, cluster->get_cluster_id());
                    break;
                }
                prev_pass_test_start_p = test_start_p;
                prev_pass_break_wcp = break_wcp;
                prev_pass_break_idx = break_idx;
                have_prev_kink_pass = true;
            } else {
                break;
            }
        }
        
        // Check if we should break the segment
        if (ray_length(Ray{start_v->wcpt().point, break_wcp}) > 1.0 * units::cm) {
            std::list<Facade::geo_point_t> wcps_list1;
            std::list<Facade::geo_point_t> wcps_list2;
            
            bool flag_break;
            bool flag_pass_check = false;
            
            // Check if end vertex is close to break point and has only one connection
            if (ray_length(Ray{end_v->wcpt().point, break_wcp}) < 1.0 * units::cm) {
                auto vd = end_v->get_descriptor();
                if (boost::degree(vd, graph) == 1) {
                    flag_pass_check = true;
                }
            }
            
            { auto t_op_pbt = BS_Clock::now();
            flag_break = proto_break_tracks(*cluster, start_v->wcpt().point, break_wcp, 
                                           end_v->wcpt().point, wcps_list1, wcps_list2, flag_pass_check);
            t_proto_break_tracks += BS_MS(BS_Clock::now() - t_op_pbt); }
            
            if (flag_break) {
                if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "break_segments: cluster={} flag_break={} list1_front=({:.2f},{:.2f},{:.2f}) list1_back=({:.2f},{:.2f},{:.2f}) list2_front=({:.2f},{:.2f},{:.2f}) list2_back=({:.2f},{:.2f},{:.2f})",
                    cluster->get_cluster_id(), flag_break,
                    wcps_list1.front().x(), wcps_list1.front().y(), wcps_list1.front().z(),
                    wcps_list1.back().x(), wcps_list1.back().y(), wcps_list1.back().z(),
                    wcps_list2.front().x(), wcps_list2.front().y(), wcps_list2.front().z(),
                    wcps_list2.back().x(), wcps_list2.back().y(), wcps_list2.back().z());

                // Check geometry constraints
                Facade::geo_vector_t tv1 = end_v->wcpt().point - start_v->wcpt().point;
                Facade::geo_vector_t tv2 = end_v->wcpt().point - break_wcp;

                double min_dis = 1e9;
                for (const auto& wcp : wcps_list1) {
                    double dis = ray_length(Ray{wcp, end_v->wcpt().point});
                    if (dis < min_dis) min_dis = dis;
                }

                double angle = std::acos(tv1.dot(tv2) / (tv1.magnitude() * tv2.magnitude())) / 3.1415926 * 180.0;

                // Kink angle at break_wcp as seen in the steiner graph:
                // angle between the approach direction (second-to-last → last of list1)
                // and the departure direction (first → second of list2).
                // This is a physically meaningful angle that segment_search_kink already
                // required to be large (≥22-45°), so it is not zero for any valid kink.
                // A 15° floor guards against degenerate configurations where the steiner
                // path happens to land at a collinear point.
                double kink_angle_at_break = 0.0;
                if (wcps_list1.size() >= 2 && wcps_list2.size() >= 2) {
                    auto it1b = wcps_list1.end(); --it1b;   // last
                    auto it1a = it1b; --it1a;               // second to last
                    auto it2a = wcps_list2.begin();          // first
                    auto it2b = it2a; ++it2b;               // second
                    Facade::geo_vector_t approach = *it1b - *it1a;
                    Facade::geo_vector_t depart   = *it2b - *it2a;
                    double mag_ap = approach.magnitude(), mag_dep = depart.magnitude();
                    if (mag_ap > 0 && mag_dep > 0) {
                        double cos_kink = approach.dot(depart) / (mag_ap * mag_dep);
                        cos_kink = std::max(-1.0, std::min(1.0, cos_kink));
                        kink_angle_at_break = std::acos(cos_kink) / 3.1415926 * 180.0;
                    }
                }

                auto end_vd = end_v->get_descriptor();
                bool end_is_terminus = (boost::degree(end_vd, graph) == 1);

                // A degree-1 end vertex with a short tail and a real kink angle at the
                // break point is a track-endpoint fold-back artifact.  Absorb it via
                // replace_segment_and_vertex (trim the segment to break_wcp, remove the
                // stub) rather than materialising a tiny separate segment.
                // - min_dis < 2 cm:  the tail is too short to represent real physics
                // - kink_angle_at_break > 30°: the break falls at a genuine direction
                //   change in the steiner graph, not a collinear snap artefact
                // - angle < 45°: the stub (tv2=end-break) must be roughly aligned with
                //   the whole segment direction (tv1=end-start).  A fold-back artifact
                //   at the track endpoint has tv1≈tv2 (small angle, ~28° in the test
                //   event), whereas a real Michel electron or physical secondary
                //   diverges from the parent track axis (~79° for a 12.94 MeV Michel).
                //   This clause prevents the absorption of real short-track secondaries.
                // The original 1.5 cm+120° condition (for degree-!=1 junctions) is kept.
                bool use_replace = (end_is_terminus && min_dis / units::cm < 2.0
                                                    && kink_angle_at_break > 30.0
                                                    && angle < 45.0) ||
                                   (!end_is_terminus && min_dis / units::cm < 1.5 && angle > 120);

                // Check if we should replace end vertex instead of breaking
                if (use_replace) {
                    if (end_is_terminus) {
                        // Replace segment and end vertex
                        SegmentPtr new_seg = curr_sg;
                        VertexPtr new_vtx = end_v;
                        auto t_op_rsv = BS_Clock::now();
                        bool rsv_ok = replace_segment_and_vertex(graph, new_seg, new_vtx, wcps_list1, break_wcp, dv);
                        t_replace_segment_and_vertex += BS_MS(BS_Clock::now() - t_op_rsv);
                        if (rsv_ok) {
                            flag_modified = true;
                            // Perform tracking
                            // track_fitter.add_graph(&graph); added already
                            auto t_op_mt = BS_Clock::now();
                            // flag_exclusion stays hard false here, NOT m_fit_exclusion
                            // (doc pr/30 §11 P1): break_segments is one of only two
                            // places the PROTOTYPE also turns exclusion off --
                            // NeutrinoID_proto_vertex.h:722/751, "fit dQ/dx here, do
                            // not exclude others".  This is the one site pair where
                            // the port already matches; knobbing it would break the
                            // parity it is meant to restore.
                            track_fitter.do_multi_tracking(true, true, false, false, false, cluster);
                            t_do_multi_tracking += BS_MS(BS_Clock::now() - t_op_mt);
                        }
                    }
                    // degree!=1 + (min_dis<1.5 && angle>120): original "do nothing" preserved
                } else {
                    // Break segment into two
                    SegmentPtr out_seg2 = nullptr;
                    auto t_op_bst = BS_Clock::now();
                    bool bst_ok = break_segment_into_two(graph, start_v, curr_sg, end_v, wcps_list1, break_wcp, wcps_list2, dv, out_seg2);
                    t_break_segment_into_two += BS_MS(BS_Clock::now() - t_op_bst);
                    if (bst_ok) {
                        flag_modified = true;
                        // doc pr/48, 59335 fix (b): a break born from a
                        // dQ/dx-assisted accept (C4) or the wide-baseline
                        // cathode accept (A0) whose break point sits on
                        // BRAGG-HOT charge (5-point mean dQ/dx >
                        // kink_dqdx_hot_ratio x the median threshold; 59335's
                        // proton stub reads 2.5-6x) is high-confidence --
                        // protect its new vertex from examine_vertices_4's
                        // unconditional < 2 cm stub-absorption floor.
                        // Without the hot gate the protect fired on ~100/1000
                        // events (every C4 stub break); with it the footprint
                        // is the genuinely hot handful (doc pr/48 sec 9.6).
                        // The new vertex is the out_seg2 endpoint sitting
                        // exactly at break_wcp (break_segment_into_two sets
                        // its wcpt to break_point verbatim).
                        // Scope: only a stub SHORT enough for
                        // examine_vertices_4's < 2 cm absorption floor needs
                        // protecting -- a longer arm survives EV4 anyway,
                        // and flagging it would needlessly shield its vertex
                        // from every other examiner pass (46/1000 collateral
                        // movers before this scope, doc pr/48 sec 9.6).
                        bool protect_hot = false;
                        if (m_kink_break_protect &&
                            (accept_crit_break == 0 || accept_crit_break == 4) && out_seg2 &&
                            ray_length(Ray{break_wcp, end_v->wcpt().point}) < 2.0 * units::cm) {
                            const auto& pfits = out_seg2->fits().empty() ? curr_sg->fits() : out_seg2->fits();
                            double sum_dQ = 0, sum_dx = 0;
                            if (!pfits.empty()) {
                                size_t kmin = 0;
                                double dmin = 1e18;
                                for (size_t i = 0; i < pfits.size(); i++) {
                                    const double d = ray_length(Ray{pfits[i].point, break_wcp});
                                    if (d < dmin) { dmin = d; kmin = i; }
                                }
                                for (int j = -2; j <= 2; j++) {
                                    const int idx = static_cast<int>(kmin) + j;
                                    if (idx >= 0 && idx < static_cast<int>(pfits.size())) {
                                        sum_dQ += pfits[idx].dQ;
                                        sum_dx += pfits[idx].dx;
                                    }
                                }
                            }
                            protect_hot = sum_dx > 0 &&
                                (sum_dQ / sum_dx) > m_mip_dqdx_median * m_kink_dqdx_hot_ratio;
                        }
                        if (protect_hot) {
                            auto [nv1, nv2] = find_vertices(graph, out_seg2);
                            for (VertexPtr nv : {nv1, nv2}) {
                                if (nv && ray_length(Ray{nv->wcpt().point, break_wcp}) < 0.01 * units::cm) {
                                    nv->set_flags(VertexFlags::kProtectedBreak);
                                    SPDLOG_LOGGER_DEBUG(s_log,
                                        "break_segments: protected kink-break vertex at ({:.2f},{:.2f},{:.2f})cm crit={} (cluster {})",
                                        break_wcp.x()/units::cm, break_wcp.y()/units::cm,
                                        break_wcp.z()/units::cm, accept_crit_break,
                                        cluster->get_cluster_id());
                                }
                            }
                        }
                        // Perform tracking
                        // track_fitter.add_graph(&graph); added already
                        auto t_op_mt = BS_Clock::now();
                        // Hard false, as above: break_segments matches prototype
                        // NeutrinoID_proto_vertex.h:751 (doc pr/30 §11 P1).
                        track_fitter.do_multi_tracking(true, true, false, false, false, cluster);
                        t_do_multi_tracking += BS_MS(BS_Clock::now() - t_op_mt);
                        if (out_seg2) {
                            remaining_segments.push_back(out_seg2);
                        }
                    }
                }
            }

        }
    }

    // Post-process: merge any vertices that ended up at the same position
    // (within 0.1 cm).  The oscillating-break pattern leaves duplicate vertex
    // objects at identical wcpt locations; merging them collapses the zigzag
    // into the correct linear topology and triggers a refit.
    {
        // remaining_segments is always empty here; use the cluster captured at
        // entry (the shared graph holds many clusters, and "first edge" is
        // pointer-ordered — it used to pick a different cluster run-to-run).
        if (entry_cluster) {
            merge_nearby_vertices(graph, *entry_cluster, track_fitter, dv);
        }
    }

   

    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "break_segments timing: segment_search_kink={:.3f}ms proto_extend_point={:.3f}ms proto_break_tracks={:.3f}ms replace_segment_and_vertex={:.3f}ms break_segment_into_two={:.3f}ms do_multi_tracking={:.3f}ms",
        t_segment_search_kink.count(), t_proto_extend_point.count(), t_proto_break_tracks.count(),
        t_replace_segment_and_vertex.count(), t_break_segment_into_two.count(), t_do_multi_tracking.count());


    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "break_segments: cluster={} break tracks -- # of Vertices: {}; # of Segments: {}",
        11, boost::num_vertices(graph), boost::num_edges(graph));

        
    return flag_modified;
}


bool PatternAlgorithms::merge_nearby_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv)
{
    // Build an ordered list of vertices belonging to this cluster.
    std::vector<VertexPtr> all_vertices;
    for (const auto& vd : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd].vertex;
        if (vtx && vtx->cluster() == &cluster) {
            all_vertices.push_back(vtx);
        }
    }

    bool any_modified = false;

    // Pass 1: merge co-located vertices (within 0.1 cm).
    // After merging vtx_j into vtx_i, vtx_i retains its position, so we
    // continue scanning the tail of the list against the same vtx_i without
    // a full restart.  This reduces the cost from O(N² × merges) to O(N²).
    for (size_t i = 0; i < all_vertices.size(); ++i) {
        VertexPtr vtx1 = all_vertices[i];
        if (!vtx1 || !vtx1->descriptor_valid()) continue;

        for (size_t j = i + 1; j < all_vertices.size(); ) {
            VertexPtr vtx2 = all_vertices[j];
            if (!vtx2 || !vtx2->descriptor_valid()) { ++j; continue; }

            if (ray_length(Ray{vtx1->wcpt().point, vtx2->wcpt().point}) >= 0.1 * units::cm) {
                ++j; continue;
            }

            // Collect any segment directly connecting vtx1↔vtx2 (would become
            // a self-loop after the merge and must be removed beforehand).
            std::vector<SegmentPtr> direct_segs;
            auto v2d = vtx2->get_descriptor();
            for (auto eit : sorted_out_edges(v2d, graph)) {
                SegmentPtr sg = graph[eit].segment;
                if (!sg) continue;
                auto [sv1, sv2] = find_vertices(graph, sg);
                if ((sv1 == vtx1 || sv2 == vtx1) && (sv1 == vtx2 || sv2 == vtx2))
                    direct_segs.push_back(sg);
            }

            if (merge_vertex_into_another(graph, vtx2, vtx1, dv)) {
                for (auto sg : direct_segs) remove_segment(graph, sg);
                all_vertices.erase(all_vertices.begin() + j);
                any_modified = true;
                // Do NOT increment j: the element now at position j must also
                // be checked against vtx1 (which kept its position).
            } else {
                ++j;
            }
        }
    }

    // Second pass: remove segments completely covered by another segment.
    // Since the graph is linear (every vertex has at most 2 edges), this
    // means all fit points of the candidate segment lie within 0.3 cm of
    // some other segment.  When that happens the two endpoints of the
    // covered segment are effectively redundant — merge them.
    bool flag_covered = true;
    while (flag_covered) {
        flag_covered = false;

        // Collect all edges in a stable vector so we can iterate safely.
        // Stable edge-index order throughout: the loop below removes/merges
        // the FIRST covered segment, so iteration order is a graph mutation
        // decision and must not be pointer order.
        std::vector<SegmentPtr> all_segments;
        for (const auto& vd : ordered_nodes(graph)) {
            for (auto edesc : sorted_out_edges(vd, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (sg && sg->cluster() == &cluster) {
                    // Only add once (undirected graph exposes both directions)
                    if (std::find(all_segments.begin(), all_segments.end(), sg) == all_segments.end()) {
                        all_segments.push_back(sg);
                    }
                }
            }
        }

        for (SegmentPtr seg : all_segments) {
            const auto& fits = seg->fits();
            if (fits.empty()) continue;

            auto [vtx_a_pre, vtx_b_pre] = find_vertices(graph, seg);
            if (!vtx_a_pre || !vtx_b_pre) continue;
            if (!vtx_a_pre->descriptor_valid() || !vtx_b_pre->descriptor_valid()) continue;

            // Gather adjacent segments (at most 2 in a linear graph).
            // Use a fixed-size stack array — no heap allocation, no std::find.
            SegmentPtr neighbors[4];
            int n_neighbors = 0;
            for (VertexPtr end_vtx : {vtx_a_pre, vtx_b_pre}) {
                auto vd = end_vtx->get_descriptor();
                // Stable edge-index order: the first covered neighbor wins below.
                for (auto edesc : sorted_out_edges(vd, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (!sg || sg == seg) continue;
                    bool dup = false;
                    for (int k = 0; k < n_neighbors; ++k) if (neighbors[k] == sg) { dup = true; break; }
                    if (!dup && n_neighbors < 4) neighbors[n_neighbors++] = sg;
                }
            }

            // Check only adjacent segments
            for (int ni = 0; ni < n_neighbors; ++ni) {
                SegmentPtr seg_other = neighbors[ni];
                const auto& fits_other = seg_other->fits();
                if (fits_other.empty()) continue;

                // Test whether ALL fit points of seg lie within 0.3 cm of seg_other
                bool all_covered = true;
                for (const auto& fp : fits) {
                    auto [dis, _pt] = segment_get_closest_point(seg_other, fp.point, "fit");
                    if (dis >= 0.3 * units::cm) {
                        all_covered = false;
                        break;
                    }
                }

                if (!all_covered) continue;

                // seg is completely covered by seg_other.
                // The endpoint shared between seg and seg_other is the junction
                // vertex that should be eliminated; the non-shared endpoint lies
                // along seg_other's path and must be kept.
                // e.g.  vtx_p --(segP)-- vtx_a --(seg)-- vtx_b --(segQ)--
                //   seg covered by segP → vtx_a is shared → merge vtx_a into vtx_b
                auto [so_v1, so_v2] = find_vertices(graph, seg_other);
                VertexPtr vtx_to_remove = (vtx_a_pre == so_v1 || vtx_a_pre == so_v2) ? vtx_a_pre : vtx_b_pre;
                VertexPtr vtx_to_keep   = (vtx_to_remove == vtx_a_pre) ? vtx_b_pre : vtx_a_pre;

                // seg will become a self-loop after merge; remove it first.
                remove_segment(graph, seg);
                if (merge_vertex_into_another(graph, vtx_to_remove, vtx_to_keep, dv)) {
                    all_vertices.erase(std::remove(all_vertices.begin(), all_vertices.end(), vtx_to_remove), all_vertices.end());
                    flag_covered = true;
                    any_modified = true;
                }
                break;  // restart outer loop
            }
            if (flag_covered) break;
        }
    }

    if (any_modified) {
        // Hard false, deliberately NOT knobbed (doc pr/30 §11 P1).
        // merge_nearby_vertices is toolkit-only (doc pr/30 P6): the prototype's
        // nearest analogue, clean_up_maps_vertices_segments, has both call sites
        // commented out.  There is therefore no prototype answer for what
        // flag_exclusion should be here, and inventing one is exactly what
        // CLAUDE.md §5 rule 4 forbids.  It is also called from the tail of
        // break_segments, whose regime is exclusion-off.
        track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
    }

    return any_modified;
}

bool PatternAlgorithms::merge_two_segments_into_one(Graph& graph, SegmentPtr& seg1, VertexPtr& vtx, SegmentPtr& seg2, IDetectorVolumes::pointer dv){
    // Get cluster from seg1 (should be same as seg2)
    auto cluster = seg1->cluster();
    if (!cluster || cluster != seg2->cluster()) {
        return false;
    }
    
    // Get the other vertices (not vtx) from seg1 and seg2
    VertexPtr vtx1 = find_other_vertex(graph, seg1, vtx);
    VertexPtr vtx2 = find_other_vertex(graph, seg2, vtx);
    
    if (!vtx1 || !vtx2) {
        return false;
    }
    
    // Create new segment from vtx1 to vtx2
    SegmentPtr new_seg = create_segment_from_vertices(graph, *cluster, vtx1, vtx2, dv);
    if (!new_seg) {
        return false;
    }
    
    // Delete old segments
    remove_segment(graph, seg1);
    remove_segment(graph, seg2);

    // Delete the middle vertex (only if now isolated)
    if (vtx->descriptor_valid()) {
        auto vd = vtx->get_descriptor();
        if (boost::degree(vd, graph) == 0) {
            remove_vertex(graph, vtx);
        }
    }
    
    // Update output parameter
    seg1 = new_seg;
    
    return true;
}

bool PatternAlgorithms::merge_vertex_into_another(Graph& graph, VertexPtr& vtx_from, VertexPtr& vtx_to, IDetectorVolumes::pointer dv){
    if (!vtx_from || !vtx_to) {
        return false;
    }
    
    // Step 1: Check if there's a segment between vtx_from and vtx_to, and delete it
    SegmentPtr seg_between = find_segment(graph, vtx_from, vtx_to);
    if (seg_between) {
        remove_segment(graph, seg_between);
    }
    
    // Step 2 & 3: Check if vertices are at the same position
    double distance = ray_length(Ray{vtx_from->wcpt().point, vtx_to->wcpt().point});
    bool same_position = (distance < 0.1 * units::cm);
    
    // Get all segments connected to vtx_from
    if (!vtx_from->descriptor_valid()) {
        return false;
    }
    
    auto vd_from = vtx_from->get_descriptor();
    std::vector<std::pair<SegmentPtr, VertexPtr>> segments_to_reconnect;

    // Collect all segments and their other endpoints.  Stable edge-index
    // order: the reconnect loop below re-adds/recreates segments, assigning
    // NEW graph indices in this order — pointer order would permute segment
    // ids run-to-run.
    for (auto edesc : sorted_out_edges(vd_from, graph)) {
        SegmentPtr seg = graph[edesc].segment;
        if (seg) {
            VertexPtr other_vtx = find_other_vertex(graph, seg, vtx_from);
            if (other_vtx) {
                segments_to_reconnect.push_back(std::make_pair(seg, other_vtx));
            }
        }
    }
    if (std::getenv("WCT_DET_DEBUG") && std::getenv("WCT_DET_DEBUG")[0] == '2') {
        fprintf(stderr, "WCT_DETM merge from=(%.4f,%.4f,%.4f) to=(%.4f,%.4f,%.4f) nrec=%zu:",
                vtx_from->wcpt().point.x(), vtx_from->wcpt().point.y(), vtx_from->wcpt().point.z(),
                vtx_to->wcpt().point.x(), vtx_to->wcpt().point.y(), vtx_to->wcpt().point.z(),
                segments_to_reconnect.size());
        for (auto& [s, ov] : segments_to_reconnect)
            fprintf(stderr, " oldidx=%zu->(%.2f,%.2f,%.2f)", s->get_graph_index(),
                    ov->wcpt().point.x(), ov->wcpt().point.y(), ov->wcpt().point.z());
        fprintf(stderr, "\n");
    }
    
    // Process each segment
    for (auto& [old_seg, other_vtx] : segments_to_reconnect) {
        auto cluster = old_seg->cluster();
        if (!cluster) continue;
        
        SegmentPtr new_seg;
        
        if (same_position) {
            // Step 2: Vertices are at same position - just reconnect the segment
            // Extract the path from the old segment
            remove_segment(graph, old_seg);
            add_segment(graph, old_seg, vtx_to, other_vtx);
        } else {
            // Step 3: Vertices are not at same position - recalculate the path
            // Remove old segment first
            remove_segment(graph, old_seg);
            // Create new segment from vtx_to to other_vtx
            new_seg = create_segment_from_vertices(graph, *cluster, vtx_to, other_vtx, dv);
            // create_segment_from_vertices already adds the segment to the graph
        }
    }
    
    // Step 4: Delete vtx_from
    remove_vertex(graph, vtx_from);
    
    return true;
}

Facade::geo_vector_t PatternAlgorithms::vertex_get_dir(VertexPtr& vertex, Graph& graph, double dis_cut){
    if (!vertex || !vertex->cluster()) {
        return Facade::geo_vector_t(0, 0, 0);
    }
    
    Facade::geo_point_t center(0, 0, 0);
    int ncount = 0;
    
    // Get vertex position (use fit if available, otherwise wcpt)
    Facade::geo_point_t vtx_point = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
    
    // Loop through all segments in the graph.  Stable edge/node-index order:
    // 'center' FP-accumulates continuous coordinates, so pointer-order
    // iteration would shift the averaged direction run-to-run.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;

        // Skip if segment doesn't belong to same cluster
        if (!sg || sg->cluster() != vertex->cluster()) continue;

        // Get points from segment (skip first and last)
        const auto& fits = sg->fits();
        if (fits.size() > 2) {
            for (size_t i = 1; i + 1 < fits.size(); i++) {
                double dis = std::sqrt(std::pow(fits[i].point.x() - vtx_point.x(), 2) +
                                      std::pow(fits[i].point.y() - vtx_point.y(), 2) +
                                      std::pow(fits[i].point.z() - vtx_point.z(), 2));
                if (dis < dis_cut) {
                    center = center + Facade::geo_vector_t(fits[i].point.x(), fits[i].point.y(), fits[i].point.z());
                    ncount++;
                }
            }
        }
    }

    // Loop through all vertices in the graph
    for (const auto& vd : ordered_nodes(graph)) {
        VertexPtr other_vtx = graph[vd].vertex;
        
        // Skip if vertex doesn't belong to same cluster
        if (!other_vtx || other_vtx->cluster() != vertex->cluster()) continue;
        
        // Get other vertex position
        Facade::geo_point_t other_vtx_point = other_vtx->fit().valid() ? other_vtx->fit().point : other_vtx->wcpt().point;
        
        double dis = std::sqrt(std::pow(other_vtx_point.x() - vtx_point.x(), 2) +
                              std::pow(other_vtx_point.y() - vtx_point.y(), 2) +
                              std::pow(other_vtx_point.z() - vtx_point.z(), 2));
        if (dis < dis_cut) {
            center = center + Facade::geo_vector_t(other_vtx_point.x(), other_vtx_point.y(), other_vtx_point.z());
            ncount++;
        }
    }
    
    if (ncount == 0) {
        return Facade::geo_vector_t(0, 0, 0);
    }
    
    // Calculate average center
    center = Facade::geo_vector_t(center.x() / ncount, center.y() / ncount, center.z() / ncount);
    
    // Calculate direction from vertex to center
    Facade::geo_vector_t dir(center.x() - vtx_point.x(),
                            center.y() - vtx_point.y(),
                            center.z() - vtx_point.z());
    
    // Normalize to unit vector
    double mag = dir.magnitude();
    if (mag > 0) {
        dir = Facade::geo_vector_t(dir.x() / mag, dir.y() / mag, dir.z() / mag);
    }
    
    return dir;
}
Facade::geo_vector_t PatternAlgorithms::vertex_segment_get_dir(VertexPtr& vertex, SegmentPtr& segment, Graph& graph, double dis_cut){
    // Return zero vector if inputs are invalid
    if (!vertex || !segment) {
        return Facade::geo_vector_t(0, 0, 0);
    }
    
    // Check if this segment is connected to this vertex
    if (!vertex->descriptor_valid()) {
        return Facade::geo_vector_t(0, 0, 0);
    }
    
    bool segment_connected = false;
    auto vd = vertex->get_descriptor();
    const auto edge_range = sorted_out_edges(vd, graph);
    for (auto eit : edge_range) {
        if (graph[eit].segment == segment) {
            segment_connected = true;
            break;
        }
    }
    
    if (!segment_connected) {
        return Facade::geo_vector_t(0, 0, 0);
    }
    
    // Get vertex position (use fit if available, otherwise wcpt)
    Facade::geo_point_t vtx_point = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
    
    // Get fit points from segment (matching prototype's cal_dir_3vector which uses fit_pt_vec)
    const auto& pts = segment->fits();

    // Find the point on the segment whose distance from vertex is closest to dis_cut
    double min_dis = 1e9;
    Facade::geo_point_t min_point = vtx_point;

    for (size_t i = 0; i < pts.size(); i++) {
        double tmp_dis = std::sqrt(std::pow(pts[i].point.x() - vtx_point.x(), 2) +
                                   std::pow(pts[i].point.y() - vtx_point.y(), 2) +
                                   std::pow(pts[i].point.z() - vtx_point.z(), 2));
        if (std::fabs(tmp_dis - dis_cut) < min_dis) {
            min_dis = std::fabs(tmp_dis - dis_cut);
            min_point = pts[i].point;
        }
    }
    
    // Calculate direction from vertex to the found point
    Facade::geo_vector_t dir(min_point.x() - vtx_point.x(),
                            min_point.y() - vtx_point.y(),
                            min_point.z() - vtx_point.z());
    
    // Normalize to unit vector
    double mag = dir.magnitude();
    if (mag > 0) {
        dir = Facade::geo_vector_t(dir.x() / mag, dir.y() / mag, dir.z() / mag);
    }
    
    return dir;
}

bool PatternAlgorithms::find_proto_vertex(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_break_track, int nrounds_find_other_tracks, bool flag_back_search, const Clus::ParticleDataSet::pointer& particle_data){
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    // Check if steiner point cloud exists and has enough points
    if (!cluster.has_pc("steiner_pc")) return false;

    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    if (steiner_pc.size() < 2) return false;

    // Initialize first segment
    Facade::Cluster* main_cluster_ptr = cluster.get_flag(Facade::Flags::main_cluster) ? &cluster : nullptr;
    SegmentPtr sg1 = init_first_segment(graph, cluster, main_cluster_ptr, track_fitter, dv, flag_back_search);
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: init_first_segment took {} ms", MS(Clock::now() - t0).count());

    if (!sg1) return false;
    // std::cout << "Fits: " << sg1->fits().size() << std::endl;

    // Store initial pair of vertices for main cluster
    std::pair<VertexPtr, VertexPtr> main_cluster_initial_pair_vertices{nullptr, nullptr};
    bool is_main_cluster = cluster.get_flag(Facade::Flags::main_cluster);

    if (is_main_cluster) {
        main_cluster_initial_pair_vertices = find_vertices(graph, sg1);
    }

    // Check if segment has more than one point
    const auto& wcpts = sg1->wcpts();
    if (wcpts.size() <= 1) {
        return false;
    }

    
 
  

    // Break tracks and examine structure
    if (flag_break_track) {
        t0 = Clock::now();
        std::vector<SegmentPtr> remaining_segments;
        remaining_segments.push_back(sg1);
        break_segments(graph, track_fitter, dv, remaining_segments);
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: break_segments took {} ms", MS(Clock::now() - t0).count());

        //  // Debug: print all segments and vertices
        // {
        //     // Print all segments
        //     auto [ebegin, eend] = boost::edges(graph);
        //     for (auto eit = ebegin; eit != eend; ++eit) {
        //         SegmentPtr sg = graph[*eit].segment;
        //         if (!sg) continue;

        //         double length = segment_track_length(sg) / units::cm;
        //         auto [v1, v2] = find_vertices(graph, sg);
        //         WireCell::Point p1_wcp = v1 ? v1->wcpt().point : WireCell::Point(0, 0, 0);
        //         WireCell::Point p2_wcp = v2 ? v2->wcpt().point : WireCell::Point(0, 0, 0);
        //         bool p1_fit_valid = (v1 && v1->fit().valid());
        //         bool p2_fit_valid = (v2 && v2->fit().valid());
        //         WireCell::Point p1_fit = p1_fit_valid ? v1->fit().point : WireCell::Point(0, 0, 0);
        //         WireCell::Point p2_fit = p2_fit_valid ? v2->fit().point : WireCell::Point(0, 0, 0);

        //         std::cout << "[Segment] length=" << length << " cm"
        //             << "  v1_wcps=(" << p1_wcp.x() << "," << p1_wcp.y() << "," << p1_wcp.z() << ")"
        //             << "  v1_fit="
        //             << (p1_fit_valid ? "(" + std::to_string(p1_fit.x()) + "," + std::to_string(p1_fit.y()) + "," + std::to_string(p1_fit.z()) + ")" : std::string("invalid"))
        //             << "  v2_wcps=(" << p2_wcp.x() << "," << p2_wcp.y() << "," << p2_wcp.z() << ")"
        //             << "  v2_fit="
        //             << (p2_fit_valid ? "(" + std::to_string(p2_fit.x()) + "," + std::to_string(p2_fit.y()) + "," + std::to_string(p2_fit.z()) + ")" : std::string("invalid"))
        //                 << std::endl;
        //     }

        //     // Print all vertices
        //     auto [vbegin, vend] = boost::vertices(graph);
        //     for (auto vit = vbegin; vit != vend; ++vit) {
        //         VertexPtr vtx = graph[*vit].vertex;
        //         if (!vtx) continue;

        //         WireCell::Point wcp_pt = vtx->wcpt().point;
        //         bool fit_valid = vtx->fit().valid();
        //         WireCell::Point fit_pt = fit_valid ? vtx->fit().point : WireCell::Point(0, 0, 0);
        //         int degree = static_cast<int>(boost::degree(*vit, graph));

        //         std::cout << "[Vertex] wcps=(" << wcp_pt.x() << "," << wcp_pt.y() << "," << wcp_pt.z() << ")"
        //             << "  fit="
        //             << (fit_valid ? "(" + std::to_string(fit_pt.x()) + "," + std::to_string(fit_pt.y()) + "," + std::to_string(fit_pt.z()) + ")" : std::string("invalid"))
        //                 << "  num_segments=" << degree
        //                 << std::endl;
        //     }
        // }

        t0 = Clock::now();
        examine_structure(graph, cluster, track_fitter, dv);
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: examine_structure took {} ms", MS(Clock::now() - t0).count());
    } else {
        t0 = Clock::now();
        track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: do_multi_tracking (no break) took {} ms", MS(Clock::now() - t0).count());
    }


    // Find other segments
    for (int i = 0; i < nrounds_find_other_tracks; i++) {
        t0 = Clock::now();
        find_other_segments(graph, cluster, track_fitter, dv, flag_break_track);
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: find_other_segments round {} took {} ms", i, MS(Clock::now() - t0).count());
    }

  

    // For main cluster, merge tracks if angles are consistent
    if (is_main_cluster) {
        t0 = Clock::now();
        if (examine_structure_3(graph, cluster, track_fitter, dv)) {
            track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        }
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: examine_structure_3 took {} ms", MS(Clock::now() - t0).count());
    }

    // Examine the vertices
    t0 = Clock::now();
    examine_vertices(graph, cluster, track_fitter, dv);
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: examine_vertices took {} ms", MS(Clock::now() - t0).count());

    // doc pr/48: two-end residual-range back-to-back break.  Placed AFTER
    // examine_vertices so the new vertex is not visible to ES2/ES3/EV1-4 on
    // this pass; examine_partial_identical_segments skips degree <= 2 and
    // examine_vertices_3 only touches the two initial termini, so the break
    // survives to the final do_multi_tracking below, which re-fits the two
    // new arms.  Downstream re-runs (improve_vertex, examine_structure_final)
    // honor VertexFlags::kProtectedBreak.  No-op unless m_two_end_break and
    // particle_data are both set.
    if (m_two_end_break && particle_data) {
        t0 = Clock::now();
        break_two_end_dqdx(graph, cluster, dv, particle_data);
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: break_two_end_dqdx took {} ms", MS(Clock::now() - t0).count());
    }

    // Examine partial identical segments
    t0 = Clock::now();
    examine_partial_identical_segments(graph, cluster, track_fitter, dv);
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: examine_partial_identical_segments took {} ms", MS(Clock::now() - t0).count());

    // Examine the two initial points for main cluster
    if (is_main_cluster && main_cluster_initial_pair_vertices.first) {
        t0 = Clock::now();
        examine_vertices_3(graph, cluster, main_cluster_initial_pair_vertices, track_fitter, dv);
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: examine_vertices_3 took {} ms", MS(Clock::now() - t0).count());
    }

    // Final multi-tracking
    t0 = Clock::now();
    track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: final do_multi_tracking took {} ms", MS(Clock::now() - t0).count());

    // Verify that at least one segment for this cluster survived all the merging/cleanup.
    // If all segments were removed (e.g. Type II merge on the only segment), return false
    // so the caller can fall back to init_point_segment.
    bool has_segment = false;
    {
        for (const auto& ed : ordered_edges(graph)) {
            auto seg = graph[ed].segment;
            if (seg && seg->cluster() == &cluster) { has_segment = true; break; }
        }
    }
    if (!has_segment) {
        SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex: no segments remain for cluster after processing, returning false");
        if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: TOTAL took {} ms", MS(Clock::now() - t_total).count());
        return false;
    }

    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "find_proto_vertex timing: TOTAL took {} ms", MS(Clock::now() - t_total).count());

    return true;
}


bool PatternAlgorithms::break_two_end_dqdx(Graph& graph, Facade::Cluster& cluster, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data)
{
    // doc sbnd_xin/docs/pr/48 sec 6.  Every gate here is deliberately
    // conservative: any missing prerequisite silently declines (no throw, no
    // partial modification).
    if (!m_two_end_break || !particle_data) return false;
    if (!cluster.get_flag(Facade::Flags::main_cluster)) return false;

    // Topology gate: exactly one segment of this cluster longer than the stub
    // floor ("no non-stub prong" -- short attached stubs like 57903's 1.9 cm /
    // 57485's <= 3.6 cm are common even in genuine back-to-back events, so a
    // strict single-segment test would wrongly exclude them).
    SegmentPtr cand = nullptr;
    int n_long = 0;
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        if (segment_track_length(sg, 0) > m_teb_stub_max) {
            n_long++;
            cand = sg;
        }
    }
    if (n_long != 1 || !cand) return false;
    const auto& fits = cand->fits();
    if (fits.size() < 3) return false;

    // Containment gate, primary: two stopping ends is physically impossible
    // for a through-going or exiting track, so both fitted endpoints must be
    // inside the fiducial volume.  A missing FiducialUtils (ill-formed
    // ensemble) never fires -- conservative.
    auto grouping = cluster.grouping();
    if (!grouping) return false;
    auto fiducial_utils = grouping->get_fiducialutils();
    if (!fiducial_utils) return false;
    if (!fiducial_utils->inside_fiducial_volume(fits.front().point) ||
        !fiducial_utils->inside_fiducial_volume(fits.back().point)) return false;

    TwoEndBreakOptions opt;
    opt.mip_dqdx        = m_mip_dqdx;
    opt.mip_dqdx_median = m_mip_dqdx_median;
    opt.min_len         = m_teb_min_len;
    opt.min_arm         = m_teb_min_arm;
    opt.min_arm_pts     = m_teb_min_arm_pts;
    opt.accept_range    = m_teb_accept_range;
    opt.rise_r1         = m_teb_rise_r1;
    opt.rise_r2         = m_teb_rise_r2;
    opt.abs_end_min     = m_teb_abs_end_min;
    opt.dip_floor       = m_teb_dip_floor;
    opt.score_cap_r1    = m_teb_score_cap_r1;
    opt.score_cap_r2    = m_teb_score_cap_r2;
    opt.turn_angle      = m_teb_turn_angle;
    opt.turn_baseline   = m_teb_turn_baseline;
    opt.turn_skirt      = m_teb_turn_skirt;

    auto res = segment_two_end_break_scan(cand, particle_data, opt);
    SPDLOG_LOGGER_DEBUG(s_log,
        "break_two_end_dqdx: cluster {} seg len {:.1f}cm k*={} (dip {} turn {}) arms {:.1f}/{:.1f}cm "
        "J={:.3f} s15=({:.3f}{},{:.3f}{}) rise=({:.2f},{:.2f}) absmed=({:.2f},{:.2f})xMIP "
        "turn={:.1f}deg routes=({},{}) found={}",
        cluster.get_cluster_id(), segment_track_length(cand, 0)/units::cm, res.break_idx,
        res.idx_dip, res.idx_turn,
        res.arm_a_len/units::cm, res.arm_b_len/units::cm, res.joint_score,
        res.sA, res.flagA ? "F" : "f", res.sB, res.flagB ? "F" : "f",
        res.ratio_lo, res.ratio_hi,
        res.absmed_lo/m_mip_dqdx_median, res.absmed_hi/m_mip_dqdx_median,
        res.turn_deg, res.route1, res.route2, res.found);
    for (const auto& a : res.attempts) {
        SPDLOG_LOGGER_DEBUG(s_log,
            "break_two_end_dqdx:   cand idx={} m3={:.2f}xMIP sA={:.3f}{} sB={:.3f}{} accepted={}",
            a.idx, a.m3/m_mip_dqdx_median, a.sA, a.fA ? "F" : "f", a.sB, a.fB ? "F" : "f", a.accepted);
    }
    if (!res.found) return false;
    if (res.break_idx <= 0 || res.break_idx + 1 >= static_cast<int>(fits.size())) return false;

    const WireCell::Point break_pt = fits[res.break_idx].point;
    auto [ok, segs, vtx] = break_segment(graph, cand, break_pt, particle_data, m_recomb_model, dv);
    if (!ok || !vtx) return false;
    // break_segment does not associate the new vertex with a cluster; a
    // null-cluster vertex is invisible to determine_main_vertex's candidate
    // loops (they filter on vtx->cluster() == &cluster), which would defeat
    // the entire purpose of the break.  (Deliberately set HERE, not inside
    // break_segment -- its other caller keeps its current behavior.)
    vtx->cluster(&cluster);
    // Mark both arms: their travel direction away from the junction is
    // established by the accept itself (each arm's Bragg is at its outer
    // end).  determine_direction reconstructs the outward direction from
    // this flag + the kProtectedBreak endpoint and lets it stand over a
    // WEAK KS recompute -- the weak coin-flip landing "into the junction"
    // is exactly what made the scorer keep the old terminus on
    // 51513/56211/57485 while 57903's lucky flips selected the junction.
    // (A dirsign stamp cannot carry this: separate_track_shower's
    // shower_topo_reset zeroes dirsign on every segment.)
    auto [sg1, sg2] = segs;
    if (sg1) sg1->set_flags(SegmentFlags::kTwoEndBreakArm);
    if (sg2) sg2->set_flags(SegmentFlags::kTwoEndBreakArm);
    // Protect the new vertex from the downstream merge/absorb passes
    // (improve_vertex's examine_vertices re-run, examine_structure_final) --
    // a straight class-A junction is exactly the geometry ES2/ESF1 are
    // designed to re-merge.
    vtx->set_flags(VertexFlags::kProtectedBreak);
    SPDLOG_LOGGER_DEBUG(s_log,
        "break_two_end_dqdx: BROKE cluster {} at fit idx {} ({:.2f},{:.2f},{:.2f})cm route {}",
        cluster.get_cluster_id(), res.break_idx,
        break_pt.x()/units::cm, break_pt.y()/units::cm, break_pt.z()/units::cm,
        res.route1 ? 1 : 2);
    return true;
}

void PatternAlgorithms::init_point_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv) {
    // Get two boundary points from the cluster (using regular point cloud)
    auto boundary_wcps = cluster.get_two_boundary_wcps(false);

    // Ensure "relaxed_pid" graph is built and cached before do_rough_path_reg_pc,
    // which calls graph_algorithms() without dv/pcts and would throw if not cached.
    cluster.graph_algorithms("relaxed_pid", dv, track_fitter.get_pc_transforms());

    // Find shortest path using the regular point cloud with "relaxed_pid" graph
    auto path_points = do_rough_path_reg_pc(cluster, boundary_wcps.first, boundary_wcps.second, "relaxed_pid");
    
    // Check if path has enough points
    if (path_points.size() <= 1) {
        return;
    }
    
    // Create vertices for the endpoints
    VertexPtr v1 = make_vertex(graph);
    v1->wcpt().point = boundary_wcps.first;
    v1->cluster(&cluster);
    
    VertexPtr v2 = make_vertex(graph);
    v2->wcpt().point = boundary_wcps.second;
    v2->cluster(&cluster);
    
    // Create segment with the path points
    auto sg1 = create_segment_for_cluster(cluster, dv, path_points);
    if (!sg1) {
        remove_vertex(graph, v1);
        remove_vertex(graph, v2);
        return;
    }
    
    // Add segment to graph connecting the two vertices
    add_segment(graph, sg1, v1, v2);
    
    // Perform multi-tracking to fit the segment
    track_fitter.add_segment(sg1);
    track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
}

void PatternAlgorithms::transfer_info_from_segment_to_cluster(Graph& graph, Facade::Cluster& cluster,  const std::string& cloud_name){
    // Get the number of points in the cluster
    const size_t npoints = cluster.npoints();
    
    // Initialize arrays for segment ID and shower flag (-1 means no segment assigned)
    std::vector<int> point_segment_id(npoints, -1);
    std::vector<int> point_flag_shower(npoints, 0);
    
    // Iterate through all edges (segments) in the graph, in stable edge-index order.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr seg = graph[ed].segment;

        // Skip if segment is null or doesn't belong to this cluster
        if (!seg || seg->cluster() != &cluster) continue;

        // Get the edge index as the segment ID
        const auto& edge_bundle = graph[ed];
        int segment_id = static_cast<int>(edge_bundle.index);
        
        seg->set_id(segment_id);

        // Check if segment is a shower (either shower trajectory or shower topology)
        bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) || 
                        seg->flags_any(SegmentFlags::kShowerTopology);
        
        // Get the local-to-global index mapping for this segment's point cloud
        if (seg->has_global_indices(cloud_name)) {
            const auto& global_indices = seg->global_indices(cloud_name);
            
            // Map each local point to the global cluster point
            for (size_t local_idx = 0; local_idx < global_indices.size(); ++local_idx) {
                size_t global_idx = global_indices[local_idx];
                
                // Validate global index
                if (global_idx < npoints) {
                    point_segment_id[global_idx] = segment_id;
                    point_flag_shower[global_idx] = is_shower ? 1 : 0;
                }
            }
        }
    }
    
    // Invalidate cache before updating arrays
    cluster.invalidate_segment_data();
    
    // Add the arrays to the cluster's default point cloud as named arrays
    auto& local_pcs = cluster.local_pcs();
    
    // Get or create the default point cloud
    // The default point cloud should already exist, but we need to add arrays to it
    // We'll add to the root node's local_pcs
    auto& default_pc = local_pcs["3d"];  // The default 3D point cloud
    
    // Erase old arrays if they exist (allows re-adding)
    default_pc.erase("point_segment_id");
    default_pc.erase("point_flag_shower");
    
    // Add the arrays
    using namespace WireCell::PointCloud;
    default_pc.add("point_segment_id", Array(point_segment_id));
    default_pc.add("point_flag_shower", Array(point_flag_shower));
}


void PatternAlgorithms::print_segs_info(Graph& graph, Facade::Cluster& cluster, VertexPtr vertex){
    // Iterate through all segments in the graph (stable order for
    // reproducible diagnostic output)
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        
        // Skip if segment is null or doesn't belong to this cluster
        if (!sg || sg->cluster() != &cluster) continue;
        
        // If a specific vertex is provided, check if segment is connected to it
        if (vertex != nullptr) {
            auto [v1, v2] = find_vertices(graph, sg);
            if (v1 != vertex && v2 != vertex) continue;
        }
        
        // Determine if segment is "in" or "out" relative to the specific vertex
        int in_vertex = 0; // 0: no direction, -1: in, 1: out
        
        if (vertex != nullptr) {
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            WireCell::Point vtx_point = vertex->wcpt().point;
            auto front_pt = wcpts.front().point;
            auto back_pt = wcpts.back().point;
            
            double dis_front = ray_length(Ray{vtx_point, front_pt});
            double dis_back = ray_length(Ray{vtx_point, back_pt});
            
            bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
            
            // Check if segment points into or out of the vertex
            if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                in_vertex = -1; // pointing into vertex
            } else if ((flag_start && sg->dirsign() == 1) || (!flag_start && sg->dirsign() == -1)) {
                in_vertex = 1; // pointing out of vertex
            }
        }
        
        // Get segment properties
        int seg_id = sg->get_graph_index();
        double length = segment_track_length(sg) / units::cm;
        int flag_dir = sg->dirsign();
        int particle_type = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;
        double particle_mass = sg->has_particle_info() ? sg->particle_info()->mass() / units::MeV : 0;
        double kinetic_energy = sg->has_particle_info() ? (sg->particle_info()->energy() - sg->particle_info()->mass()) / units::MeV : 0;
        bool is_dir_weak = seg_dir_weak(sg);
        
        // Determine segment type and print
        if (sg->flags_any(SegmentFlags::kShowerTopology)) {
            std::cout << "print_segs_info: " << seg_id << " " << length << " S_topo "
                      << flag_dir << " " << particle_type << " " << particle_mass << " "
                      << kinetic_energy << " " << is_dir_weak << " " << in_vertex << "\n";
        } else if (sg->flags_any(SegmentFlags::kShowerTrajectory)) {
            std::cout << "print_segs_info: " << seg_id << " " << length << " S_traj "
                      << flag_dir << " " << particle_type << " " << particle_mass << " "
                      << kinetic_energy << " " << is_dir_weak << " " << in_vertex << "\n";
        } else {
            std::cout << "print_segs_info: " << seg_id << " " << length << " Track  "
                      << flag_dir << " " << particle_type << " " << particle_mass << " "
                      << kinetic_energy << " " << is_dir_weak << " " << in_vertex << "\n";
        }
    }
}


std::pair<Facade::geo_point_t, Facade::geo_vector_t> PatternAlgorithms::calc_PCA_main_axis(std::vector<Facade::geo_point_t>& points){
    Facade::geo_point_t center(0, 0, 0);
    int nsum = 0;
    
    // Calculate the center point
    for (size_t i = 0; i < points.size(); i++) {
        center += points[i];
        nsum++;
    }
    
    Facade::geo_vector_t PCA_main_axis(0, 0, 0);
    
    // Need at least 3 points for PCA
    if (nsum < 3) {
        center.set(0, 0, 0);
        return std::make_pair(center, PCA_main_axis);
    }
    
    center /= nsum;
    
    // Build covariance matrix
    Eigen::MatrixXd cov_matrix(3, 3);
    
    for (int i = 0; i != 3; i++) {
        for (int j = i; j != 3; j++) {
            cov_matrix(i, j) = 0;
            for (size_t k = 0; k < points.size(); k++) {
                cov_matrix(i, j) += (points[k][i] - center[i]) * (points[k][j] - center[j]);
            }
        }
    }
    
    // Fill in symmetric part
    cov_matrix(1, 0) = cov_matrix(0, 1);
    cov_matrix(2, 0) = cov_matrix(0, 2);
    cov_matrix(2, 1) = cov_matrix(1, 2);
    
    // Compute eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(cov_matrix);
    auto eigen_vectors = eigenSolver.eigenvectors();
    
    // Get the principal component (largest eigenvalue, which is index 2 in ascending order)
    int i = 2;  // Eigen returns eigenvalues in ascending order, so index 2 is the largest
    double norm = sqrt(eigen_vectors(0, i) * eigen_vectors(0, i) + 
                      eigen_vectors(1, i) * eigen_vectors(1, i) + 
                      eigen_vectors(2, i) * eigen_vectors(2, i));
    
    PCA_main_axis.set(eigen_vectors(0, i) / norm, 
                     eigen_vectors(1, i) / norm, 
                     eigen_vectors(2, i) / norm);
    
    return std::make_pair(center, PCA_main_axis);
}


Facade::geo_vector_t PatternAlgorithms::calc_dir_cluster(Graph& graph, Facade::Cluster& cluster, const Facade::geo_point_t& orig_p, double dis_cut){
    Facade::geo_point_t ave_p(0, 0, 0);
    int num = 0;
    const double dis_cut_sq = dis_cut * dis_cut;

    // Iterate through all segments in the graph that belong to this cluster.
    // Stable edge/node-index order: ave_p FP-accumulates continuous
    // coordinates, so pointer-order iteration would shift the direction.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;

        // Use fitted points (mirrors prototype's fit_pt_vec / get_point_vec()),
        // skipping first and last points as in the original.
        const auto& fits = sg->fits();
        for (size_t i = 1; i + 1 < fits.size(); i++) {
            const WireCell::Point& pt = fits[i].point;
            double dx = pt.x() - orig_p.x();
            double dy = pt.y() - orig_p.y();
            double dz = pt.z() - orig_p.z();
            if (dx*dx + dy*dy + dz*dz < dis_cut_sq) {
                ave_p.set(ave_p.x() + pt.x(), ave_p.y() + pt.y(), ave_p.z() + pt.z());
                num++;
            }
        }
    }

    // Iterate through all vertices in the graph that belong to this cluster
    for (const auto& vd : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;

        // Prefer fit point if available, otherwise fall back to wcpt
        const WireCell::Point& vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;

        double dx = vtx_pt.x() - orig_p.x();
        double dy = vtx_pt.y() - orig_p.y();
        double dz = vtx_pt.z() - orig_p.z();
        if (dx*dx + dy*dy + dz*dz < dis_cut_sq) {
            ave_p.set(ave_p.x() + vtx_pt.x(), ave_p.y() + vtx_pt.y(), ave_p.z() + vtx_pt.z());
            num++;
        }
    }

    // Calculate direction vector
    Facade::geo_vector_t dir(0, 0, 0);

    if (num > 0) {
        ave_p.set(ave_p.x() / num, ave_p.y() / num, ave_p.z() / num);
        dir.set(ave_p.x() - orig_p.x(), ave_p.y() - orig_p.y(), ave_p.z() - orig_p.z());
        double magnitude = std::sqrt(dir.x() * dir.x() + dir.y() * dir.y() + dir.z() * dir.z());
        if (magnitude > 0) {
            dir.set(dir.x() / magnitude, dir.y() / magnitude, dir.z() / magnitude);
        }
    }

    return dir;
}


   Facade::Cluster* PatternAlgorithms::swap_main_cluster(Facade::Cluster& new_main_cluster, Facade::Cluster& old_main_cluster, std::vector<Facade::Cluster*>& other_clusters){
       // doc pr/59: this function has no log line at any call site (DL rerank,
       // check_switch_main_cluster[_2], the two_end_break protected-cluster
       // path) -- a swap here is otherwise invisible except by cross-
       // referencing print_segs_info's graph indices against a cluster's
       // segment set by hand.  Env-gated (WCT_PR59_ASSOC_CENSUS unset => no
       // log line, no behavior change) to match the clustering_points_segments
       // sentinels in PRSegmentFunctions.cxx and TaggerCheckNeutrino.cxx.
       static const bool pr59_assoc_census = std::getenv("WCT_PR59_ASSOC_CENSUS") != nullptr;
       if (pr59_assoc_census) {
           SPDLOG_LOGGER_DEBUG(s_log,
               "pr59 assoc-census: swap_main_cluster {} -> {}",
               old_main_cluster.get_cluster_id(), new_main_cluster.get_cluster_id());
       }

       // Remove main_cluster flag from old main cluster (set to 0 to unset)
       old_main_cluster.set_flag(Facade::Flags::main_cluster, 0);
       
       // Add old main cluster to other_clusters
       other_clusters.push_back(&old_main_cluster);
       
       // Set main_cluster flag on new main cluster
       new_main_cluster.set_flag(Facade::Flags::main_cluster);
       
       // Remove new_main_cluster from other_clusters
       auto it = std::find_if(other_clusters.begin(), other_clusters.end(), 
                             [&new_main_cluster](Facade::Cluster* c) {
                                 return c == &new_main_cluster;
                             });
       
       if (it != other_clusters.end()) {
           other_clusters.erase(it);
       }
       
       return &new_main_cluster;
    }

    void PatternAlgorithms::examine_main_vertices(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters){
        if (!main_cluster) return;
        
        // Calculate cluster length cut
        double main_cluster_length = main_cluster->get_length();
        double cluster_length_cut = std::min(main_cluster_length * 0.6, 6.0 * units::cm);
        
        // First pass: remove short clusters without strong tracks
        std::vector<Facade::Cluster*> clusters_to_be_removed;
        
        for (auto& [cluster, vertex] : map_cluster_main_vertices) {
            if (!cluster || !vertex) continue;
            
            double length = cluster->get_length();
            
            if (length < cluster_length_cut) {
                bool flag_removed = true;
                
                // Check segments connected to this vertex
                if (vertex->descriptor_valid()) {
                    auto vd = vertex->get_descriptor();
                    const auto edge_range = sorted_out_edges(vd, graph);
                    
                    for (auto e_it : edge_range) {
                        SegmentPtr seg = graph[e_it].segment;
                        if (!seg) continue;
                        
                        bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                        seg->flags_any(SegmentFlags::kShowerTopology) ||
                                        (seg->has_particle_info() && seg->particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
                        int dirsign = seg->dirsign();
                        bool is_dir_weak = seg_dir_weak(seg);
                        double median_dqdx = segment_median_dQ_dx(seg) / m_mip_dqdx_median;
                        
                        // Keep if: not shower AND has direction AND (strong direction OR high dQ/dx)
                        if (!is_shower && dirsign != 0 && (!is_dir_weak || median_dqdx > 2.0)) {
                            flag_removed = false;
                            break;
                        }
                    }
                }
                
                // Additional check for very short clusters
                if (!flag_removed) {
                    // CreateSteinerGraph legitimately leaves a cluster without steiner
                    // products when create_steiner_tree finds <2 terminals ("only 0
                    // steiner terminal(s) found (need >=2), returning empty graph" ->
                    // "skipping transfer", clus/src/CreateSteinerGraph.cxx:226/241),
                    // and kd_steiner_knn then *throws* rather than returning empty
                    // (Facade_Cluster.cxx:3821).  The `if (!knn.empty())` blocks here
                    // and below already treat "no steiner answer" as "skip the test",
                    // so honour that intent instead of aborting the job: 67/1000
                    // SBND MCP2025C events died exactly here.
                    if (length < 5 * units::cm && main_cluster->has_pc("steiner_pc")) {
                        WireCell::Point vtx_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
                        auto knn = main_cluster->kd_steiner_knn(1, vtx_pt, "steiner_pc");
                        if (!knn.empty()) {
                            double closest_dis = std::sqrt(knn[0].second);
                            if (closest_dis > 100 * units::cm) {
                                flag_removed = true;
                            }
                        }
                    }
                }
                
                if (flag_removed) {
                    clusters_to_be_removed.push_back(cluster);
                }
            } else if (main_cluster->has_pc("steiner_pc")) {   // see the steiner guard above
                // For longer clusters, check if very far from main cluster
                WireCell::Point vtx_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
                auto knn = main_cluster->kd_steiner_knn(1, vtx_pt, "steiner_pc");
                if (!knn.empty()) {
                    double closest_dis = std::sqrt(knn[0].second);
                    if (closest_dis > 200 * units::cm) {
                        clusters_to_be_removed.push_back(cluster);
                    }
                }
            }
        }
        
        // Remove flagged clusters
        for (auto cluster : clusters_to_be_removed) {
            map_cluster_main_vertices.erase(cluster);
        }
        clusters_to_be_removed.clear();
        
        // Second pass: additional cuts if main cluster has a main vertex
        if (map_cluster_main_vertices.find(main_cluster) != map_cluster_main_vertices.end()) {
            // Check which clusters have only showers
            std::map<int, bool> map_cluster_id_shower;
            for (auto& [cluster, vertex] : map_cluster_main_vertices) {
                if (cluster) {
                    map_cluster_id_shower[cluster->get_cluster_id()] = true;
                }
            }
            
            // Check all segments to see if any cluster has non-shower segments
            for (const auto& ed : ordered_edges(graph)) {
                SegmentPtr sg = graph[ed].segment;
                if (!sg) continue;
                
                int cluster_id = sg->cluster() ? sg->cluster()->get_cluster_id() : -1;
                if (map_cluster_id_shower.find(cluster_id) == map_cluster_id_shower.end()) continue;
                
                bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                sg->flags_any(SegmentFlags::kShowerTopology) ||
                                (sg->has_particle_info() && sg->particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
                if (!is_shower) {
                    map_cluster_id_shower[cluster_id] = false;
                }
            }
            
            bool flag_main_vertex_all_showers = map_cluster_id_shower[main_cluster->get_cluster_id()];
            
            if (flag_main_vertex_all_showers) {
                // Calculate direction from main cluster's main vertex
                VertexPtr main_vtx = map_cluster_main_vertices[main_cluster];
                WireCell::Point main_vtx_pt = main_vtx->fit().valid() ? main_vtx->fit().point : main_vtx->wcpt().point;
                Facade::geo_vector_t dir_main = calc_dir_cluster(graph, *main_cluster, main_vtx_pt, 15 * units::cm);
                
                // Sort by cluster_id for deterministic evaluation order
                std::vector<std::pair<Facade::Cluster*, VertexPtr>> sorted_candidates(
                    map_cluster_main_vertices.begin(), map_cluster_main_vertices.end());
                std::sort(sorted_candidates.begin(), sorted_candidates.end(),
                          [](const auto& a, const auto& b) {
                              int aid = a.first ? a.first->get_cluster_id() : -1;
                              int bid = b.first ? b.first->get_cluster_id() : -1;
                              return aid < bid;
                          });

                for (auto& [cluster, vertex] : sorted_candidates) {
                    if (cluster == main_cluster || !vertex) continue;

                    WireCell::Point vtx_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
                    
                    // Get closest distance to main cluster
                    auto knn = main_cluster->kd_steiner_knn(1, vtx_pt, "steiner_pc");
                    double closest_dis = knn.empty() ? 1e9 : std::sqrt(knn[0].second);
                    
                    // Calculate direction vector from main vertex to this vertex
                    Facade::geo_vector_t dir1(vtx_pt.x() - main_vtx_pt.x(),
                                             vtx_pt.y() - main_vtx_pt.y(),
                                             vtx_pt.z() - main_vtx_pt.z());
                    
                    double angle = 0;
                    if (dir_main.magnitude() > 0 && dir1.magnitude() > 0) {
                        double cos_angle = std::clamp(dir_main.dot(dir1) / (dir_main.magnitude() * dir1.magnitude()), -1.0, 1.0);
                        angle = std::acos(cos_angle) / M_PI * 180.0;
                    }
                    
                    double cluster_length = cluster->get_length();
                    
                    if (angle < 10) {
                        // Cluster in same direction as main - check if small and close
                        if ((cluster_length < 15 * units::cm && closest_dis < 40 * units::cm) ||
                            (cluster_length < 7 * units::cm && closest_dis < 60 * units::cm)) {
                            clusters_to_be_removed.push_back(cluster);
                        }
                    } else if (angle > 160) {
                        // Cluster in opposite direction - check for main cluster swap
                        if (map_cluster_id_shower[cluster->get_cluster_id()] && 
                            cluster_length > 10 * units::cm && 
                            cluster_length > 0.5 * main_cluster_length) {
                            
                            // Calculate direction from this cluster
                            Facade::geo_vector_t dir2 = calc_dir_cluster(graph, *cluster, vtx_pt, 15 * units::cm);
                            
                            // Get closest distance between point clouds
                            auto closest_result = main_cluster->get_closest_points(*cluster);
                            double closest_dis_pc = std::get<2>(closest_result);
                            
                            double angle2 = 0;
                            if (dir2.magnitude() > 0 && dir_main.magnitude() > 0) {
                                double cos_angle2 = std::clamp(dir2.dot(dir_main) / (dir2.magnitude() * dir_main.magnitude()), -1.0, 1.0);
                                angle2 = std::acos(cos_angle2) / M_PI * 180.0;
                            }
                            
                            if (closest_dis_pc < 10 * units::cm && angle2 < 25) {
                                // Swap main cluster
                                main_cluster = swap_main_cluster(*cluster, *main_cluster, other_clusters);
                            }
                        }
                    }
                }
            }
            
            // Remove flagged clusters
            for (auto cluster : clusters_to_be_removed) {
                map_cluster_main_vertices.erase(cluster);
            }
            clusters_to_be_removed.clear();
        }
    }