// MCS muon momentum driver -- see the header and doc 80 secs 4/7/8.

#include "WireCellClus/MuonMCSDriver.h"

#include "WireCellClus/Facade_Cluster.h"  // Segment::cluster() -> get_cluster_id()
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellMcs/MuonMCS.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;
using WireCell::GraphTools::mir;

namespace {

    using PR::Graph;
    using PR::SegmentPtr;
    using PR::VertexPtr;

    // Deterministic segment collection: (graph edge index, segment), sorted
    // by index -- never pointer order (M4 discipline).
    std::vector<std::pair<size_t, SegmentPtr>> ordered_segments(Graph& graph)
    {
        std::vector<std::pair<size_t, SegmentPtr>> out;
        for (auto edesc : mir(boost::edges(graph))) {
            SegmentPtr seg = graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            out.emplace_back(graph[edesc].index, seg);
        }
        std::sort(out.begin(), out.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        return out;
    }

    double dist2(const Point& a, const Point& b)
    {
        const double dx = a.x() - b.x(), dy = a.y() - b.y(), dz = a.z() - b.z();
        return dx * dx + dy * dy + dz * dz;
    }

}  // namespace

void PR::mcs_fill_kine(KineInfo& kine, Graph& graph,
                       const IndexedSegmentSet& segments_in_long_muon,
                       const VertexPtr& main_vertex,
                       bool beam_gate_active,
                       const MuonMCSConfig& cfg,
                       WireCell::Log::logptr_t log)
{
    if (!cfg.enable) return;

    // Correctness gate (doc 80 sec 7.4): the bundle candidates are already
    // spill-coincident when the beam-window gate is active; without an active
    // gate there is no window to certify coincidence against, so skip.
    if (cfg.beam_window_only && !beam_gate_active) {
        SPDLOG_LOGGER_DEBUG(log, "mcs: beam window gate inactive, skipping (mcs_beam_window_only)");
        return;
    }

    // ---- muon selection (doc 80 sec 4.2) ----
    std::vector<SegmentPtr> muon_segments;
    SegmentPtr id_segment = nullptr;  // carrier of the kine_mcs_segment_id join key
    bool muon_typed_exists = false;

    auto all = ordered_segments(graph);
    for (const auto& [idx, seg] : all) {
        if (seg->particle_info() && std::abs(seg->particle_info()->pdg()) == 13) {
            muon_typed_exists = true;
            break;
        }
    }

    if (cfg.muon_source == "pf_muon") {
        // ubreco WireCellMCS_module.cc:197-207 reproduced on the PR graph,
        // with the owner-decided robustness changes: abs(pdg) == 13 (a
        // stamped -13 must not silently select nothing) and ranking by
        // kinetic energy (monotone-equivalent to ubreco's total energy at
        // fixed mass).  Ties break to the smaller graph index (first in the
        // ordered scan).
        double best_ke = -1;
        for (const auto& [idx, seg] : all) {
            const auto& pi = seg->particle_info();
            if (!pi || std::abs(pi->pdg()) != 13) continue;
            if (pi->kinetic_energy() > best_ke) {
                best_ke = pi->kinetic_energy();
                id_segment = seg;
            }
        }
        if (id_segment) muon_segments.push_back(id_segment);
    }
    else if (cfg.muon_source == "long_muon") {
        // the examine_direction chain (NeutrinoVertexFinder.cxx:1860-1919),
        // function-local to visit() -- the reason this driver is a call site.
        // IndexedSegmentSet iterates in segment-index order (deterministic).
        double best_len = -1;
        for (const auto& seg : segments_in_long_muon) {
            muon_segments.push_back(seg);
            const double len = segment_track_length(seg);
            if (len > best_len) {
                best_len = len;
                id_segment = seg;  // join key: the longest chain member
            }
        }
    }
    else if (cfg.muon_source == "longest_segment") {
        // PID-independent control arm: longest segment by track length.
        double best_len = -1;
        for (const auto& [idx, seg] : all) {
            const double len = segment_track_length(seg);
            if (len > best_len) {
                best_len = len;
                id_segment = seg;
            }
        }
        if (id_segment) muon_segments.push_back(id_segment);
    }
    else {
        SPDLOG_LOGGER_WARN(log, "mcs: unknown muon_source '{}', skipping", cfg.muon_source);
        return;
    }

    if (muon_segments.empty()) {
        if (muon_typed_exists) {
            // owner decision 6: prove the selection rule agrees with the
            // stamped PIDs rather than assume it
            SPDLOG_LOGGER_WARN(log, "mcs: a muon-typed segment exists but none was selected (source={})",
                               cfg.muon_source);
        }
        return;
    }

    double total_length_cm = 0;
    for (const auto& seg : muon_segments) { total_length_cm += segment_track_length(seg) / units::cm; }
    if (total_length_cm < cfg.muon_min_length_cm) {
        SPDLOG_LOGGER_DEBUG(log, "mcs: selected muon too short ({:.1f} < {} cm)", total_length_cm,
                            cfg.muon_min_length_cm);
        return;
    }

    // ---- endpoints: the extreme vertices' fit().point (owner decision 5) ----
    // For a chain, the extremes are the vertices used by exactly one chain
    // segment.  The endpoint nearer the main vertex becomes the muon START
    // (energy decreases along the track from there), mirroring ubreco's PF
    // Position()/EndPosition() -- which, like these, were NOT members of the
    // cloud handed to MCS (trim resolves them to the nearest cloud point and
    // re-inserts the true endpoint, upstream mcs.cxx:353-354).
    std::vector<std::pair<VertexPtr, int>> vertex_use;  // (vertex, #segments using it)
    for (const auto& seg : muon_segments) {
        auto [v1, v2] = find_vertices(graph, seg);
        for (const auto& v : { v1, v2 }) {
            if (!v) continue;
            bool found = false;
            for (auto& [vv, n] : vertex_use) {
                if (vv == v) {
                    n++;
                    found = true;
                    break;
                }
            }
            if (!found) vertex_use.emplace_back(v, 1);
        }
    }
    std::vector<VertexPtr> extremes;
    for (const auto& [v, n] : vertex_use) {
        if (n == 1) extremes.push_back(v);
    }
    if (extremes.size() != 2) {
        // loop or degenerate chain: fall back to the most separated pair
        double best = -1;
        VertexPtr a = nullptr, b = nullptr;
        for (size_t i = 0; i < vertex_use.size(); i++) {
            for (size_t j = i + 1; j < vertex_use.size(); j++) {
                double d = dist2(vertex_use[i].first->fit().point, vertex_use[j].first->fit().point);
                if (d > best) {
                    best = d;
                    a = vertex_use[i].first;
                    b = vertex_use[j].first;
                }
            }
        }
        if (!a || !b) {
            SPDLOG_LOGGER_WARN(log, "mcs: cannot determine muon endpoints, skipping");
            return;
        }
        extremes = { a, b };
    }
    VertexPtr vstart = extremes[0], vend = extremes[1];
    if (main_vertex) {
        if (dist2(extremes[1]->fit().point, main_vertex->fit().point) <
            dist2(extremes[0]->fit().point, main_vertex->fit().point)) {
            std::swap(vstart, vend);
        }
    }

    // ---- point harvest, internal units -> cm on the boundary (sec 4.1) ----
    std::vector<std::vector<double>> points;
    auto harvest = [&points](const SegmentPtr& seg) {
        for (const auto& f : seg->fits()) {
            points.push_back({ f.point.x() / units::cm, f.point.y() / units::cm,
                               f.point.z() / units::cm });
        }
    };
    if (cfg.point_source == "whole_event") {
        // ubreco-literal whole-event cloud: every segment's fit points plus
        // every vertex fit point.  Validation only (sec 7.3): trim is
        // O(N^2 log N) and this is 50x the default cloud.
        for (const auto& [idx, seg] : all) harvest(seg);
        std::vector<std::pair<size_t, VertexPtr>> nodes;
        for (auto ndesc : mir(boost::vertices(graph))) {
            const auto& v = graph[ndesc].vertex;
            if (!v) continue;
            nodes.emplace_back(graph[ndesc].index, v);
        }
        std::sort(nodes.begin(), nodes.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        for (const auto& [idx, v] : nodes) {
            points.push_back({ v->fit().point.x() / units::cm, v->fit().point.y() / units::cm,
                               v->fit().point.z() / units::cm });
        }
    }
    else {  // muon_segments (default, sec 7.3)
        for (const auto& seg : muon_segments) harvest(seg);
    }

    if ((int)points.size() > cfg.max_points) {
        SPDLOG_LOGGER_WARN(log, "mcs: cloud has {} points > mcs_max_points {}, skipping",
                           points.size(), cfg.max_points);
        return;
    }
    if (points.size() < 2) {
        SPDLOG_LOGGER_DEBUG(log, "mcs: only {} cloud points, skipping", points.size());
        return;
    }

    std::vector<double> start = { vstart->fit().point.x() / units::cm, vstart->fit().point.y() / units::cm,
                                  vstart->fit().point.z() / units::cm };
    std::vector<double> end = { vend->fit().point.x() / units::cm, vend->fit().point.y() / units::cm,
                                vend->fit().point.z() / units::cm };

    // ---- run (bug guards default FIXED, owner decision 4) ----
    Mcs::McsOptions opt;
    opt.cathode_x = cfg.cathode_x_cm;
    opt.cathode_xcut = cfg.cathode_xcut_cm;
    Mcs::MuonMCS mcs(opt);
    Mcs::McsResult res = mcs.run(start, end, points);

    // ---- stamp KineInfo scalars (sec 8; KE in MeV, matching
    //      kine_energy_particle's units) ----
    kine.kine_mcs_energy = (res.ke_MCS >= 0) ? static_cast<float>(res.ke_MCS) : -1.0f;
    kine.kine_mcs_ambiguity = (res.ambiguity_MCS >= 0) ? static_cast<float>(res.ambiguity_MCS) : -1.0f;
    kine.kine_mcs_tracklen = (res.mu_tracklen >= 0) ? static_cast<float>(res.mu_tracklen) : -1.0f;
    kine.kine_mcs_range_energy = (res.ke_tracklen >= 0) ? static_cast<float>(res.ke_tracklen) : -1.0f;
    const auto* cl = id_segment ? id_segment->cluster() : nullptr;
    kine.kine_mcs_segment_id = cl ? cl->get_cluster_id() * 1000 + static_cast<int>(id_segment->get_graph_index())
                                  : -1;

    SPDLOG_LOGGER_INFO(log,
                       "mcs: source={} nseg={} npoints={} len={:.1f}cm seg_id={} -> "
                       "ke_MCS={:.1f} MeV amb={:.3g} tracklen={:.1f}cm ke_range={:.1f} MeV "
                       "(nsegs14={} bad_path={} cathode_drop={}/{})",
                       cfg.muon_source, muon_segments.size(), points.size(), total_length_cm,
                       kine.kine_mcs_segment_id, kine.kine_mcs_energy, kine.kine_mcs_ambiguity,
                       kine.kine_mcs_tracklen, kine.kine_mcs_range_energy, res.nsegs, res.bad_path,
                       res.counters.cathode_seg_dropped, res.counters.cathode_angle_masked);
}
