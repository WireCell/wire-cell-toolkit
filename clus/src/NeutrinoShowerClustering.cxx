#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellUtil/Logging.h"
#include <chrono>
#include <Eigen/Dense>
#include <unordered_map>
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

// doc sbnd_xin/docs/pr/84 round 3: WCT_SHOWER_CREATE_DEBUG attributes every
// Shower to the pass that built it.  Two PR::Showers can be built on the SAME
// start segment (SBND 169626/174752/347129/394532), which makes the PF tree
// emit two nodes with one jsTree id and double-counts the shower in
// kine_reco_Enu -- neither is visible in map_segment_in_shower, whose writes
// are last-writer-wins.  Log/stderr only: no effect on emitted bytes.
static inline void pr84_probe_shower(const ShowerPtr& sh, const char* site) {
    static const bool dbg = std::getenv("WCT_SHOWER_CREATE_DEBUG") != nullptr;
    if (!dbg || !sh) return;
    auto seg = sh->start_segment();
    const int sid = seg ? (seg->id() >= 0 ? seg->id() : (int) seg->get_graph_index()) : -1;
    const auto* cl = seg ? seg->cluster() : nullptr;
    auto [svtx, conn] = sh->get_start_vertex_and_type();
    std::fprintf(stderr,
                 "SHOWER_CREATE_DEBUG site=%s shower_id=%d start_seg=%d conn=%d "
                 "start_vtx_gidx=%zu nseg=%d\n",
                 site, sh->get_shower_id(),
                 cl ? cl->get_cluster_id() * 1000 + sid : sid, conn,
                 svtx ? svtx->get_graph_index() : (size_t) -1,
                 sh->get_num_segments());
}

// doc sbnd_xin/docs/pr/91 round 1 -- display ids, identical to the encodings
// PrDisplayDump.cxx (pf_node_id / vertex_display_id) and the Bee PF writer use,
// so a probe line joins straight onto a calib dump row or an mc.json node.
static inline int pr91_seg_display_id(const SegmentPtr& s)
{
    if (!s) return -1;
    int sid = s->id();
    if (sid < 0) sid = static_cast<int>(s->get_graph_index());
    const auto* c = s->cluster();
    return c ? c->get_cluster_id() * 1000 + sid : sid;
}
static inline int pr91_vtx_display_id(const VertexPtr& v)
{
    if (!v) return -1;
    const auto* c = v->cluster();
    return (c ? c->get_cluster_id() : 0) * 1000 + static_cast<int>(v->get_graph_index());
}

// doc sbnd_xin/docs/pr/93 round 3: WCT_SHOWER_ABSORB_DEBUG -- call-site tags
// binding each complete_structure_with_start_segment walk (whose per-segment
// ADD/EXCLUDE lines live in PRShower.cxx) to the pass that invoked it, plus
// direct add_segment/add_shower membership writes.  A "site=" line precedes
// its walk_begin line by shower_start_seg.  stderr only, no effect on
// emitted bytes.
static inline bool pr93_absorb_dbg()
{
    static const bool dbg = std::getenv("WCT_SHOWER_ABSORB_DEBUG") != nullptr;
    return dbg;
}
static inline void pr93_probe_absorb_site(const char* site, const ShowerPtr& sh, bool guard)
{
    if (!pr93_absorb_dbg() || !sh) return;
    std::fprintf(stderr, "SHOWER_ABSORB site=%s shower_start_seg=%d guard=%d\n",
                 site, pr91_seg_display_id(sh->start_segment()), (int)guard);
}
static inline void pr93_probe_absorb_direct(const char* site, const ShowerPtr& sh, const SegmentPtr& sg)
{
    if (!pr93_absorb_dbg() || !sh) return;
    std::fprintf(stderr, "SHOWER_ABSORB DIRECT site=%s shower_start_seg=%d seg=%d pdg=%d\n",
                 site, pr91_seg_display_id(sh->start_segment()), pr91_seg_display_id(sg),
                 sg && sg->has_particle_info() && sg->particle_info() ? sg->particle_info()->pdg() : 0);
}
static inline void pr93_probe_absorb_splice(const char* site, const ShowerPtr& into, const ShowerPtr& from)
{
    if (!pr93_absorb_dbg() || !into || !from) return;
    std::fprintf(stderr, "SHOWER_ABSORB SPLICE site=%s into_start_seg=%d from_start_seg=%d from_nseg=%d\n",
                 site, pr91_seg_display_id(into->start_segment()),
                 pr91_seg_display_id(from->start_segment()), from->get_num_segments());
}

// doc sbnd_xin/docs/pr/91 round 1: WCT_SHOWER_MERGE_DEBUG prints, at every
// shower-merge decision site, the candidate pair and EVERY quantity in the
// condition together with the verdict -- fired or not -- plus a one-line reason
// whenever a whole pass is skipped.  The four structural gates that leave one
// physical EM shower split across several PR::Showers (G-A ordering, G-B
// examine_merge_showers' conn-1<-conn-2-only 10 deg test, G-C the 3 cm conn-3
// hatch in examine_shower_1, G-D in_other_clusters' downstream-only absorb)
// show up here as SKIPs rather than as failed comparisons, and telling those
// two apart is the whole point.  Log/stderr only: no effect on emitted bytes.
static inline bool pr91_merge_dbg()
{
    static const bool dbg = std::getenv("WCT_SHOWER_MERGE_DEBUG") != nullptr;
    return dbg;
}

// doc sbnd_xin/docs/pr/91 round 1: WCT_SHOWER_CONTENT_DEBUG dumps the exact
// membership of every PR::Shower.  This is the ONLY non-lossy source: the calib
// dump's `segment.shower_id` join stores one shower per segment
// (PrDisplayDump.cxx:432), so when two showers overlap -- SBND 347129 shower 4
// owns segment 51011, which is shower 2's own START segment; SBND 394532
// shower 4 owns shower 2's start segment 31010 -- the overlapped shower comes
// back with an empty member list and looks empty rather than nested.
//
// The owner's question is "what is contained in this EM shower", so each member
// carries its own length and charge plus a running fraction of the shower's
// kine_charge, not just an id.  Log/stderr only: no effect on emitted bytes.
static void pr91_probe_shower_content(const ShowerPtr& sh, Graph& graph, int idx)
{
    static const bool dbg = std::getenv("WCT_SHOWER_CONTENT_DEBUG") != nullptr;
    if (!dbg || !sh) return;

    auto [svtx, conn] = sh->get_start_vertex_and_type();
    const WireCell::Point sp = sh->get_start_point();
    const WireCell::Point ep = sh->get_end_point();
    const WireCell::Vector d15  = shower_cal_dir_3vector(*sh, sp, 15 * WireCell::units::cm);
    const WireCell::Vector d100 = shower_cal_dir_3vector(*sh, sp, 100 * WireCell::units::cm);

    // Two passes: the first totals charge so the second can print a fraction.
    double tot_dQ = 0, tot_len = 0;
    std::set<int> clusters;
    auto edges = ordered_edges(*sh, graph);
    for (auto edesc : edges) {
        SegmentPtr seg = graph[edesc].segment;
        if (!seg || !seg->descriptor_valid()) continue;
        if (seg->cluster()) clusters.insert(seg->cluster()->get_cluster_id());
        tot_len += segment_track_length(seg);
        for (const auto& f : seg->fits()) tot_dQ += f.dQ;
    }
    const double kq = sh->get_kine_charge();

    std::fprintf(stderr,
                 "SHOWER_CONTENT shower_id=%d node_id=%d conn=%d pdg=%d nseg=%d ncls=%zu "
                 "kine_best=%.3f kine_charge=%.3f kine_range=%.3f len=%.2f sum_len=%.2f "
                 "start_vtx=%d start=(%.3f,%.3f,%.3f) end=(%.3f,%.3f,%.3f) "
                 "dir15=(%.3f,%.3f,%.3f) dir100=(%.3f,%.3f,%.3f) idx=%d\n",
                 sh->get_shower_id(), pr91_seg_display_id(sh->start_segment()), conn,
                 sh->get_particle_type(), sh->get_num_segments(), clusters.size(),
                 sh->get_kine_best()/WireCell::units::MeV, kq/WireCell::units::MeV,
                 sh->get_kine_range()/WireCell::units::MeV,
                 sh->get_total_length()/WireCell::units::cm, tot_len/WireCell::units::cm,
                 pr91_vtx_display_id(svtx),
                 sp.x()/WireCell::units::cm, sp.y()/WireCell::units::cm, sp.z()/WireCell::units::cm,
                 ep.x()/WireCell::units::cm, ep.y()/WireCell::units::cm, ep.z()/WireCell::units::cm,
                 d15.x(), d15.y(), d15.z(), d100.x(), d100.y(), d100.z(), idx);

    for (auto edesc : edges) {
        SegmentPtr seg = graph[edesc].segment;
        if (!seg || !seg->descriptor_valid()) continue;
        double sdQ = 0, sdx = 0;
        for (const auto& f : seg->fits()) { sdQ += f.dQ; sdx += f.dx; }
        auto [v0, v1] = PR::find_vertices(graph, seg);
        const WireCell::Point p0 = v0 ? v0->fit().point : WireCell::Point(0,0,0);
        const WireCell::Point p1 = v1 ? v1->fit().point : WireCell::Point(0,0,0);
        std::fprintf(stderr,
                     "SHOWER_CONTENT   shower_id=%d seg=%d cluster=%d len=%.3f pdg=%d "
                     "dQ=%.1f dQdx=%.1f Efrac=%.4f E_est=%.3f flags=%c%c%c dirsign=%d "
                     "v0=%d(%.3f,%.3f,%.3f) v1=%d(%.3f,%.3f,%.3f)\n",
                     sh->get_shower_id(), pr91_seg_display_id(seg),
                     seg->cluster() ? seg->cluster()->get_cluster_id() : -1,
                     segment_track_length(seg)/WireCell::units::cm,
                     seg->has_particle_info() && seg->particle_info()
                         ? seg->particle_info()->pdg() : 0,
                     sdQ, sdx > 0 ? sdQ/(sdx/WireCell::units::cm) : 0.0,
                     tot_dQ > 0 ? sdQ/tot_dQ : 0.0,
                     tot_dQ > 0 ? sdQ/tot_dQ*kq/WireCell::units::MeV : 0.0,
                     seg->flags_any(SegmentFlags::kShowerTrajectory) ? 'J' : '-',
                     seg->flags_any(SegmentFlags::kShowerTopology)   ? 'T' : '-',
                     seg->flags_any(SegmentFlags::kAvoidMuonCheck)   ? 'A' : '-',
                     seg->dirsign(),
                     pr91_vtx_display_id(v0), p0.x()/WireCell::units::cm, p0.y()/WireCell::units::cm, p0.z()/WireCell::units::cm,
                     pr91_vtx_display_id(v1), p1.x()/WireCell::units::cm, p1.y()/WireCell::units::cm, p1.z()/WireCell::units::cm);
    }

    // Vertices in the view that no member segment touches -- the F1 suspects:
    // get_end_point() is a farthest-VERTEX search over exactly this set
    // (PRShower.cxx:1191, :1201), so a vertex here that belongs to another
    // shower's charge is how an end point lands outside the shower.
    std::set<node_descriptor> touched;
    for (auto edesc : edges) {
        touched.insert(boost::source(edesc, graph));
        touched.insert(boost::target(edesc, graph));
    }
    for (auto vdesc : ordered_nodes(*sh, graph)) {
        if (touched.count(vdesc)) continue;
        VertexPtr v = graph[vdesc].vertex;
        if (!v) continue;
        const WireCell::Point vp = v->fit().point;
        std::fprintf(stderr,
                     "SHOWER_CONTENT   shower_id=%d ORPHAN_VTX=%d cluster=%d (%.3f,%.3f,%.3f) "
                     "dis_from_start=%.3f is_start_vtx=%d\n",
                     sh->get_shower_id(), pr91_vtx_display_id(v),
                     v->cluster() ? v->cluster()->get_cluster_id() : -1,
                     vp.x()/WireCell::units::cm, vp.y()/WireCell::units::cm, vp.z()/WireCell::units::cm,
                     (sp - vp).magnitude()/WireCell::units::cm, v == svtx ? 1 : 0);
    }
}

namespace {
    struct cluster_point_info {
        Facade::Cluster* cluster;
        double min_angle;
        double min_dis;
        VertexPtr min_vertex;
        WireCell::Point min_point;
    };
    
    bool sortbydis(const cluster_point_info &a, const cluster_point_info &b) {
        if (a.min_dis != b.min_dis) return a.min_dis < b.min_dis;
        return a.cluster->get_cluster_id() < b.cluster->get_cluster_id();
    }
}

void PatternAlgorithms::update_shower_maps(IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters){
    // Clear all maps
    map_vertex_to_shower.clear();
    map_vertex_in_shower.clear();
    map_segment_in_shower.clear();
    used_shower_clusters.clear();
    
    // Iterate through all showers
    for (auto shower : showers) {
        // Map start vertex to shower
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        if (start_vtx) {
            map_vertex_to_shower[start_vtx].insert(shower);
        }
        
        // Fill maps using TrajectoryView - iterate through all vertices and segments in the shower
        TrajectoryView& traj = shower->fill_maps();
        
        // Fill map_vertex_in_shower with all vertices in this shower except start_vertex.
        // Matches prototype: set_start_vertex() never adds to map_vtx_segs, so fill_maps()
        // never includes start_vertex in map_vertex_in_shower.
        for (auto vdesc : traj.nodes()) {
            auto vtx = traj.view_graph()[vdesc].vertex;
            if (vtx && vtx != start_vtx) {
                map_vertex_in_shower[vtx] = shower;
            }
        }
        
        // Fill map_segment_in_shower with all segments in this shower
        for (auto edesc : traj.edges()) {
            auto seg = traj.view_graph()[edesc].segment;
            if (seg) {
                map_segment_in_shower[seg] = shower;
            }
        }
    }
    
    // Collect all cluster IDs from segments in the map
    for (auto it = map_segment_in_shower.begin(); it != map_segment_in_shower.end(); it++) {
        auto seg = it->first;
        if (seg && seg->cluster()) {
            used_shower_clusters.insert(seg->cluster());
        }
    }
}

// doc sbnd_xin/docs/pr/84 round 3 -- docstring at m_shower_dedup_start_seg
// (NeutrinoPatternBase.h).  One shower per start segment: the group collapses
// onto its most-directly-connected member, which absorbs the others.
//
// Grouping is keyed by the start segment's GRAPH INDEX -- the segment's true
// identity -- not by the display id `cluster_id*1000 + seg_id`, which is what
// the PF tree emits and is only the *symptom* (two distinct segments could in
// principle share a display id; pf_unique_node_ids covers that residue).
// Iteration order is the key order of a std::map<size_t,...> plus shower_id
// within a group, so no pointer is ever iterated (pr/34 sec 6).
int PatternAlgorithms::merge_showers_sharing_start_segment(IndexedShowerSet& showers)
{
    std::map<size_t, std::vector<ShowerPtr>> by_start_seg;
    for (const auto& sh : showers) {
        if (!sh) continue;
        auto seg = sh->start_segment();
        if (!seg || !seg->descriptor_valid()) continue;
        // conn 4 ("not clearly connected") is skipped by both the PF tree
        // (MultiAlgBlobClustering.cxx) and the kine tree, so it neither
        // duplicates a node nor double-counts energy: leave it alone.
        if (sh->get_start_vertex_and_type().second == 4) continue;
        by_start_seg[seg->get_graph_index()].push_back(sh);
    }

    int n_absorbed = 0;
    for (auto& [gidx, group] : by_start_seg) {
        if (group.size() < 2) continue;
        std::sort(group.begin(), group.end(),
                  [](const ShowerPtr& a, const ShowerPtr& b) {
                      const int ca = a->get_start_vertex_and_type().second;
                      const int cb = b->get_start_vertex_and_type().second;
                      if (ca != cb) return ca < cb;                 // most direct wins
                      const int na = a->get_num_segments(), nb = b->get_num_segments();
                      if (na != nb) return na > nb;                 // then the fuller view
                      const double ka = a->get_kine_best(), kb = b->get_kine_best();
                      if (ka != kb) return ka > kb;
                      return a->get_shower_id() < b->get_shower_id();  // stable
                  });
        auto keep = group.front();
        const int keep_conn = keep->get_start_vertex_and_type().second;
        for (size_t i = 1; i < group.size(); ++i) {
            auto drop = group[i];
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr84 shower_dedup: start_seg_gidx={} keep sid={} conn={} nseg={} ke={:.1f}MeV "
                "<- drop sid={} conn={} nseg={} ke={:.1f}MeV",
                gidx, keep->get_shower_id(), keep_conn, keep->get_num_segments(),
                keep->get_kine_best()/units::MeV,
                drop->get_shower_id(), drop->get_start_vertex_and_type().second,
                drop->get_num_segments(), drop->get_kine_best()/units::MeV);
            pr93_probe_absorb_splice("dedup_keep_drop", keep, drop);
            keep->add_shower(*drop);
            showers.erase(drop);
            ++n_absorbed;
        }
        // Force the production kinematics pass to recompute this shower --
        // and ONLY this shower: every other member of `showers` still carries
        // flag_kinematics and is skipped there.
        keep->set_flag_kinematics(false);
    }
    return n_absorbed;
}

// doc sbnd_xin/docs/pr/74 round 2 K4 -- docstring at m_shower_stem_backfill
// (NeutrinoPatternBase.h).  Runs between shower_clustering_in_other_clusters
// and the second kinematics pass, so both BFS-created main-cluster showers
// (90055's) and conn-3/4 attached showers (469665's) exist and absorbed stem
// charge is counted.  Only called when the knob is on => byte-identical off.
void PatternAlgorithms::stem_backfill(Graph& graph, VertexPtr main_vertex,
    IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower,
    ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower,
    ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon,
    const Clus::ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model)
{
    if (!main_vertex || !main_vertex->descriptor_valid()) return;

    // BFS parent map over the graph from main_vertex.  sorted_out_edges +
    // first-visit-wins makes the parent choice deterministic; the
    // pointer-keyed containers are only tested/inserted, never iterated.
    std::map<VertexPtr, std::pair<SegmentPtr, VertexPtr>> came_from;
    std::set<VertexPtr> visited{main_vertex};
    std::vector<VertexPtr> frontier{main_vertex};
    while (!frontier.empty()) {
        std::vector<VertexPtr> next_frontier;
        for (auto v : frontier) {
            if (!v || !v->descriptor_valid()) continue;
            for (auto e : sorted_out_edges(v->get_descriptor(), graph)) {
                SegmentPtr sg = graph[e].segment;
                if (!sg) continue;
                VertexPtr ov = find_other_vertex(graph, sg, v);
                if (!ov || visited.count(ov)) continue;
                visited.insert(ov);
                came_from.emplace(ov, std::make_pair(sg, v));
                next_frontier.push_back(ov);
            }
        }
        frontier = std::move(next_frontier);
    }

    // Deterministic shower order by start-segment graph index.
    std::vector<ShowerPtr> sorted_showers(showers.begin(), showers.end());
    std::sort(sorted_showers.begin(), sorted_showers.end(),
              [&](const ShowerPtr& a, const ShowerPtr& b) {
                  auto sa = a->start_segment(), sb = b->start_segment();
                  if (!sa && !sb) return false;
                  if (!sa) return true;
                  if (!sb) return false;
                  return graph[sa->get_descriptor()].index < graph[sb->get_descriptor()].index;
              });

    bool any_absorbed = false;
    for (auto shower : sorted_showers) {
        auto start_seg = shower->start_segment();
        if (!start_seg || !start_seg->has_particle_info() || !start_seg->particle_info() ||
            start_seg->particle_info()->pdg() != 11) continue;   // EM showers only
        if (shower->get_total_length() < m_stem_backfill_min_shower_len) continue;
        auto [attach_vtx, conn_type] = shower->get_start_vertex_and_type();
        VertexPtr cur = attach_vtx;
        SegmentPtr outer_stem = nullptr;   // round 3 Q1: outermost stem absorbed for THIS shower
        while (cur && cur != main_vertex && came_from.count(cur)) {
            auto [stem, prev] = came_from.at(cur);
            if (!stem || map_segment_in_shower.count(stem)) break;
            if (segments_in_long_muon.count(stem)) break;
            // doc pr/74 round 4: a stem the Michel-stem guard just separated
            // OUT of this very shower must not be absorbed straight back in.
            // Gated on that knob, so this is byte-identical whenever round 4
            // is off, whatever stem_backfill is doing.
            if (m_shower_traj_michel_stem &&
                stem->flags_any(SegmentFlags::kMuonStemGuard)) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr74r4 stem_backfill: shower(start gidx={}) chain gidx={} blocked: michel-stem guard muon",
                    start_seg->get_graph_index(), stem->get_graph_index());
                break;
            }
            // Junction guard: a PF orphan anchors via vtx_incoming_seg at a
            // vertex of its anchor segment; absorbing `stem` removes that
            // anchor and the orphan vanishes from the tree, audit-only
            // (measured: 268067 stranded a 595 MeV proton branch, 285567 two
            // protons totalling 442 MeV, 56982 a track fragment).  Both
            // endpoints of `stem` must therefore be clean: every incident
            // segment is shower-claimed, or the stem itself, or the next
            // chain candidate (which stays a PF track unless its own,
            // later, guard check passes).  The main vertex is exempt -- its
            // residents are PF roots and anchor directly.
            SegmentPtr stem_seg = stem;   // plain locals: structured bindings
            int conn_type_c = conn_type;  // cannot be lambda-captured pre-C++20
            SegmentPtr next_stem = (prev && came_from.count(prev)) ? came_from.at(prev).first : nullptr;
            auto junction_ok = [&](VertexPtr v, SegmentPtr allow2) -> bool {
                if (!v || !v->descriptor_valid() || v == main_vertex) return true;
                for (auto e : sorted_out_edges(v->get_descriptor(), graph)) {
                    SegmentPtr side = graph[e].segment;
                    if (!side || side == stem_seg || side == allow2) continue;
                    if (!map_segment_in_shower.count(side)) {
                        SPDLOG_LOGGER_DEBUG(s_log,
                            "pr74 stem_backfill: shower(start gidx={} conn={}) chain gidx={} blocked: non-shower sibling gidx={} at junction",
                            start_seg->get_graph_index(), conn_type_c, stem_seg->get_graph_index(),
                            side->get_graph_index());
                        return false;
                    }
                }
                return true;
            };
            if (!junction_ok(cur, nullptr) || !junction_ok(prev, next_stem)) break;
            const double len = segment_track_length(stem);
            const double med = segment_median_dQ_dx(stem);
            const double ratio = (m_mip_dqdx_median > 0 && med > 0) ? med / m_mip_dqdx_median : 0.0;
            // mip_lo: an EM trunk carries at least MIP-level charge.  A
            // sub-MIP stub is charge-poor debris whose absorption buys
            // nothing and can strand later-evicted shower fragments that
            // anchored on it (285567: 0.66x stub, two protons stranded).
            const bool ok = len < m_stem_backfill_max_len &&
                            ratio >= m_stem_backfill_mip_lo &&
                            ratio < m_stem_backfill_mip_hi;
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr74 stem_backfill: shower(start gidx={} conn={}) chain gidx={} len {:.1f}cm dqdx {:.2f}x -> {}",
                start_seg->get_graph_index(), conn_type, stem->get_graph_index(),
                len/units::cm, ratio, ok ? "absorb" : "stop");
            if (!ok) break;
            pr93_probe_absorb_direct("stem_backfill", shower, stem);
            shower->add_segment(stem, true);
            map_segment_in_shower[stem] = shower;   // claim immediately; full rebuild below
            any_absorbed = true;
            outer_stem = stem_seg;
            cur = prev;
        }

        // doc pr/74 round 3 Q1 -- RE-SEAT THE SHOWER ONTO THE ABSORBED STEM.
        //
        // Round 2 added the stem to the shower's MEMBERSHIP only and left
        // start_segment/start_vertex where they were.  Membership alone fixes
        // the paint layer but breaks the particle-flow tree: fill_bee_pf_tree
        // walks track-only segments from the main vertex
        // (MultiAlgBlobClustering.cxx:1204-1254) and an absorbed stem is no
        // longer a track, so the vertex it used to reach drops out of
        // vtx_incoming_seg.  Every object anchored there -- the shower itself
        // and every sibling shower / pseudo-gamma hanging off that vertex --
        // lost its parent and rendered as a top-level PF root far from the
        // neutrino vertex (measured on the round-2 production arm: 90055
        // 0 -> 7 dangling roots incl. the 2020 MeV shower itself, 469665
        // 0 -> 3, 138009 0 -> 1).  That is precisely the complaint doc pr/74
        // opened with for 18255-142421 -- "painted EM shower, missing from the
        // particle flow" -- re-created by the fix for it.
        //
        // Re-seating restores the anchor at the OTHER end of the stem chain
        // (`cur`, the main-vertex side of the last absorbed stem) using the
        // same (start_vertex, conn_type=1) + start_segment pair the BFS shower
        // builder itself uses (:240-241).  The shower then renders from the
        // stem base -- where the owner's hand scan says the shower starts --
        // and the vertex-set propagation at MultiAlgBlobClustering.cxx:1287
        // re-parents the stranded siblings underneath it.
        //
        // The start segment's own PDG stays whatever the tracker said; the
        // shower's type comes from update_particle_type()'s majority vote in
        // the kinematics pass that follows (doc pr/44; PRShower.cxx:842), which
        // relabels a track-typed start segment 13 -> 11 whenever shower length
        // dominates -- always true for a substantial shower with a <30 cm stem.
        if (outer_stem) {
            // conn_type 1 takes the shower's start_point from the start
            // segment's own fit endpoints, selected by dirsign
            // (PRShower.cxx:1141-1147) -- and a stem the tracker never
            // directioned has dirsign()==0, which assigns NEITHER endpoint and
            // silently leaves start_point at its previous value (measured:
            // 90055 kept (127,24,216), 469665 kept the EM blob at
            // (30,124,356), i.e. the very stale anchors this re-seat exists to
            // move).  Point the stem away from `cur` so start_point lands on
            // the stem base.  This is also the physically right sign: the
            // vertex nearer the neutrino is where the particle started.
            const WireCell::Point cur_pt =
                cur->fit().valid() ? cur->fit().point : cur->wcpt().point;
            const auto& stem_fits = outer_stem->fits();
            if (!stem_fits.empty()) {
                const double d_front = (stem_fits.front().point - cur_pt).magnitude();
                const double d_back  = (stem_fits.back().point  - cur_pt).magnitude();
                outer_stem->dirsign(d_front <= d_back ? 1 : -1);
            }
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr74 stem_backfill: re-seat shower(start gidx={} conn={}) -> start gidx={} conn=1 dirsign={}",
                start_seg->get_graph_index(), conn_type, outer_stem->get_graph_index(),
                outer_stem->dirsign());
            shower->set_start_vertex(cur, 1);
            shower->set_start_segment(outer_stem);
            // The majority vote relabels a track-typed start segment 13 -> 11
            // whenever shower length dominates (PRShower.cxx:842) -- it must
            // run BEFORE the kinematics recompute below, which copies the
            // start segment's PDG onto the shower verbatim
            // (PRShower.cxx:1121) and would otherwise turn a 2 GeV electron
            // shower into a muon.
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx,
                                         main_vertex, m_shower_proton_daughter_pion,
                                         m_mip_dqdx_median, m_shower_vote_track_pid_counts,
                                         m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
            // calc_kine_2 skips any shower whose kinematics flag is already
            // set (NeutrinoEnergyReco.cxx:327), which calc_kine_1 set for
            // every shower that existed then.  Round 2 therefore recomputed
            // NOTHING after absorbing: start_point/end_point/init_dir kept
            // their pre-absorption values and the absorbed stem's charge was
            // never counted (90055 stayed 2020 MeV, 469665 stayed 322 MeV
            // across the round-2 knob flip).  Clearing the flag is what makes
            // the re-seat and the charge accounting actually happen.
            shower->set_flag_kinematics(false);
        }
    }
    if (any_absorbed) {
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                           map_vertex_to_shower, used_shower_clusters);
    }
}

void PatternAlgorithms::shower_clustering_with_nv_in_main_cluster(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    if (!main_vertex || !main_vertex->descriptor_valid()) return;

    // BFS from main_vertex: find shower-flagged segments anywhere in the segment tree.
    // Segment ordering is guaranteed by sorted_out_edges() which sorts by stable graph index,
    // eliminating any pointer-address-dependent randomness.
    // std::set is required here because complete_structure_with_start_segment() takes set&.
    IndexedSegmentSet used_segments;
    std::vector<ShowerPtr> new_showers;

    // Seed BFS: all segments incident to main_vertex in stable index order.
    std::vector<std::pair<SegmentPtr, VertexPtr>> segments_to_examine;
    for (auto e : sorted_out_edges(main_vertex->get_descriptor(), graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg) continue;
        VertexPtr other_vtx = find_other_vertex(graph, seg, main_vertex);
        segments_to_examine.push_back({seg, other_vtx});
    }

    while (!segments_to_examine.empty()) {
        std::vector<std::pair<SegmentPtr, VertexPtr>> temp_segments;
        for (auto& [curr_sg, daughter_vtx] : segments_to_examine) {
            // Insert returns false if already present; skip without a second lookup.
            if (!used_segments.insert(curr_sg).second) continue;

            // get_flag_shower() = kShowerTrajectory || kShowerTopology || abs(pdg)==11
            bool is_shower_seg = curr_sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                 curr_sg->flags_any(SegmentFlags::kShowerTopology) ||
                                 (curr_sg->has_particle_info() && curr_sg->particle_info() &&
                                  std::abs(curr_sg->particle_info()->pdg()) == 11);
            bool in_long_muon = segments_in_long_muon.count(curr_sg) > 0;

            if (is_shower_seg || in_long_muon) {
                // Parent vertex is the other end from daughter_vtx
                VertexPtr parent_vtx = find_other_vertex(graph, curr_sg, daughter_vtx);
                ShowerPtr shower = std::make_shared<Shower>(graph);
                shower->set_start_vertex(parent_vtx, 1);
                shower->set_start_segment(curr_sg);

                // For long muon segments, record muon particle type on the shower
                if (curr_sg->has_particle_info() && curr_sg->particle_info() &&
                    std::abs(curr_sg->particle_info()->pdg()) == 13) {
                    shower->set_particle_type(curr_sg->particle_info()->pdg());
                    SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_with_nv_in_main_cluster: Main-cluster long muon {} : {}", new_showers.size(), curr_sg->particle_info()->pdg());
                } else {
                    SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_with_nv_in_main_cluster: Main-cluster shower {}", new_showers.size());
                }
                pr84_probe_shower(shower, "nv_in_main_cluster");
                new_showers.push_back(shower);
                // BFS does not descend into shower sub-tree
            } else {
                // Track-like segment: keep descending from daughter_vtx in stable index order.
                if (daughter_vtx && daughter_vtx->descriptor_valid()) {
                    for (auto e : sorted_out_edges(daughter_vtx->get_descriptor(), graph)) {
                        SegmentPtr next_sg = graph[e].segment;
                        if (next_sg && !used_segments.count(next_sg)) {
                            VertexPtr other_vtx = find_other_vertex(graph, next_sg, daughter_vtx);
                            temp_segments.push_back({next_sg, other_vtx});
                        }
                    }
                }
            }
        }
        segments_to_examine = std::move(temp_segments);
    }

    // Complete shower structure for all newly created showers.
    // used_segments (populated during BFS) prevents overlapping segment claims.
    for (auto shower : new_showers) {
        pr93_probe_absorb_site("in_main_cluster", shower, m_shower_absorb_track_guard);
        shower->complete_structure_with_start_segment(nv_bridge_seed(used_segments), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
        // Enforce electron type on start segment:
        //  - update_particle_type() handles multi-segment showers via majority vote
        //  - explicit PDG=0 guard catches single-segment showers skipped by update_particle_type()
        //    (long-muon start segments retain PDG=13 and are handled by the post-pass below)
        //
        // doc sbnd_xin/docs/pr/44: the parenthetical above is FALSE for a
        // MULTI-segment long-muon pseudo-shower -- update_particle_type's
        // majority vote counts every non-proton member (muons included) as
        // shower_length, so a pure muon chain ALWAYS trips the
        // `shower_length > track_length` branch and the start segment is
        // relabelled 13 -> 11 (PRShower.cxx:842).  The whole
        // update_particle_type call is a toolkit-only addition (18f09178,
        // 2026-03-31); the prototype completes the structure and goes
        // straight to the deliberate long-muon -> EM reclass loop below
        // (NeutrinoID_shower_clustering.h:1709-1717) -- a long-muon shower's
        // start segment is never re-typed here.  SBND 18255 evt 142421: the
        // ~143 cm collinear MIP chain 7023->7024->7018 lost seg 7024 (13->11)
        // to this vote, which then seeded a fake "e- 163 MeV" merged into the
        // pi0.  When the knob is on, a shower whose cached particle_type was
        // recorded 13 at the seeding above keeps its muon start segment
        // (prototype parity); EM showers (cached type 0) vote exactly as
        // before.  false = legacy = byte-identical.
        const bool keep_muon_type = m_shower_long_muon_keep_type &&
                                    std::abs(shower->get_particle_type()) == 13;
        if (!keep_muon_type)
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
        // PDG=0 guard: defensive fixup for shower-flagged start segments that
        // arrived without any ParticleInfo set.  Independent of
        // update_particle_type and still needed.
        auto start_seg = shower->start_segment();
        if (start_seg && start_seg->has_particle_info() && start_seg->particle_info() &&
            start_seg->particle_info()->pdg() == 0) {
            auto four_momentum = segment_cal_4mom(start_seg, 11, particle_data, recomb_model, m_mip_dqdx);
            start_seg->particle_info(std::make_shared<Aux::ParticleInfo>(
                11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11),
                four_momentum));
        }
        showers.insert(shower);
    }

    // Check if any long muon shower should be reclassified as an EM shower.
    // Condition: the shower has more non-muon segments than muon segments (by count and length).
    // Use index-stable sets so iteration order is deterministic across runs.
    IndexedSegmentSet tmp_segments;
    IndexedVertexSet  tmp_vertices;
    for (auto shower : showers) {
        // doc pr/33 F2: the prototype reads the START SEGMENT's PDG with an
        // EXACT muon test here (NeutrinoID_shower_clustering.h:1716,
        // get_start_segment()->get_particle_type()!=13); the port reads the
        // shower's cached type through std::abs.  Both readings feed the
        // audit counter; prototype parity needs BOTH knobs on (either alone
        // is neither tree's behavior).
        {
            const int type_shower = shower->get_particle_type();
            int type_startseg = 0;
            if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                type_startseg = shower->start_segment()->particle_info()->pdg();
            }
            g_pr33_audit.f2_calls[0]++;
            if ((std::abs(type_shower) != 13) != (type_startseg != 13)) g_pr33_audit.f2_disagree[0]++;
            const int type_used = m_shower_pdg_from_start_segment ? type_startseg : type_shower;
            if (m_shower_pdg_exact_muon_test ? (type_used != 13) : (std::abs(type_used) != 13)) continue;
        }
        if (!shower->start_segment()) continue;

        double n_muons = 0, length_muons = 0;
        double n_others = 0, length_others = 0;
        double max_muon_length = 0;
        SegmentPtr max_sg = nullptr;

        // Single pass: collect segments and gather statistics simultaneously
        // to avoid iterating shower->edges() twice.
        std::vector<SegmentPtr> shower_segs;
        // ordered_edges, not shower->edges(): edges() is an unordered_set hashed
        // on heap addresses (doc pr/28 §15.2).  length_muons/length_others are FP
        // accumulations and they feed a BRANCH --
        // `length_others > 0.33 * length_muons` below -- so the walk order is not
        // confined to rounding: it can decide whether this muon is retyped to an
        // EM shower.
        for (auto edesc : ordered_edges(*shower, graph)) {
            auto sg1 = shower->view_graph()[edesc].segment;
            if (!sg1) continue;
            shower_segs.push_back(sg1);
            double length = segment_track_length(sg1);
            if (segments_in_long_muon.count(sg1)) {
                n_muons++;
                length_muons += length;
                if (length > max_muon_length ||
                    (length == max_muon_length && max_sg && sg1->id() < max_sg->id())) {
                    max_muon_length = length;
                    max_sg = sg1;
                }
            } else {
                n_others++;
                length_others += length;
            }
        }

        if (n_others >= 2 * n_muons && length_others > 0.33 * length_muons &&
            n_muons > 0 && max_muon_length < 60 * units::cm) {
            SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_with_nv_in_main_cluster: Long muon converted to EM shower");
            for (auto sg1 : shower_segs) {
                if (sg1 == max_sg) sg1->set_flags(SegmentFlags::kAvoidMuonCheck);
                if (sg1->has_particle_info() && sg1->particle_info()) {
                    sg1->particle_info()->set_pdg(11);
                    pr40_probe_setpdg(sg1, 11, "NeutrinoShowerClustering.cxx:long_muon_to_EM");
                    sg1->particle_info()->set_mass(particle_data->get_particle_mass(11));
                }
                tmp_segments.insert(sg1);
                auto [vtx1, vtx2] = find_vertices(graph, sg1);
                if (vtx1) tmp_vertices.insert(vtx1);
                if (vtx2) tmp_vertices.insert(vtx2);
            }
        }
    }

    // Remove reclassified segments/vertices from the long muon containers
    for (auto seg : tmp_segments) segments_in_long_muon.erase(seg);
    for (auto vtx : tmp_vertices) vertices_in_long_muon.erase(vtx);

    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
}

void PatternAlgorithms::shower_clustering_connecting_to_main_vertex(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters){
    if (!main_vertex) return;

    // Build map_vertex_segments from graph (ordered for determinism)
    std::map<VertexPtr, std::vector<SegmentPtr>> map_vertex_segments;
    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg) continue;

        auto source_vdesc = boost::source(e, graph);
        auto target_vdesc = boost::target(e, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;

        if (v1) map_vertex_segments[v1].push_back(seg);
        if (v2) map_vertex_segments[v2].push_back(seg);
    }

    // Step 1: Collect segments from showers starting at main_vertex.
    // Sort by start-segment ID for deterministic order (showers is a ptr-address-ordered set).
    std::vector<ShowerPtr> main_vtx_showers;
    for (auto& shower : showers) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        if (start_vtx == main_vertex) main_vtx_showers.push_back(shower);
    }
    std::sort(main_vtx_showers.begin(), main_vtx_showers.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
        auto sa = a->start_segment();
        auto sb = b->start_segment();
        if (sa && sb) return sa->id() < sb->id();
        return sa < sb;
    });

    IndexedSegmentSet used_segments;
    for (auto& shower : main_vtx_showers) {
        for (auto edesc : shower->edges()) {
            auto seg = shower->view_graph()[edesc].segment;
            if (seg) used_segments.insert(seg);
        }
        // If there's only 1 segment and it's short, clear used_segments
        if (used_segments.size() == 1 && segment_track_length(*used_segments.begin()) < 8 * units::cm) {
            used_segments.clear();
        }
    }

    // Step 2: If used_segments is empty, try to create new showers
    if (used_segments.empty()) {
        std::set<ShowerPtr> del_showers;
        std::map<ShowerPtr, SegmentPtr> map_shower_max_sg;
        ShowerPtr max_shower = nullptr;
        double max_length = 0;

        // Get segments connected to main_vertex
        auto& main_vtx_segments = map_vertex_segments[main_vertex];

        for (auto sg : main_vtx_segments) {
            // Skip segments already claimed by a shower before the expensive BFS traversal.
            if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) continue;

            // Calculate total number of daughter segments for this segment.
            // doc pr/33 F1: the prototype calls calculate_num_daughter_tracks
            // here (NeutrinoID_shower_clustering.h:140, flag=true => count
            // everything with length > 0); the port calls _showers
            // (shower-flagged only).  Both computed unconditionally so the
            // audit counter bounds the reach of the restored callee.
            auto pair_result_legacy = calculate_num_daughter_showers(graph, main_vertex, sg);
            auto pair_result_proto  = calculate_num_daughter_tracks(graph, main_vertex, sg, true, 0);
            g_pr33_audit.f1_mv_calls++;
            if (pair_result_legacy.first != pair_result_proto.first) g_pr33_audit.f1_mv_differ++;
            auto pair_result = m_daughter_count_proto_main_vertex ? pair_result_proto : pair_result_legacy;

            // Get segment properties
            double medium_dQ_dx = segment_median_dQ_dx(sg);
            double medium_dQ_dx_1 = medium_dQ_dx / m_mip_dqdx_median;

            // Get particle type if available
            int particle_type = 0;
            if (sg->has_particle_info() && sg->particle_info()) {
                particle_type = sg->particle_info()->pdg();
            }

            // F1 audit: would the skip verdict below change with the other
            // callee's count?
            {
                auto skip_with = [&](int nd) {
                    return (particle_type == 11) ||
                           (particle_type == 2212 && ((medium_dQ_dx_1 > 1.45 && nd <= 3) || medium_dQ_dx_1 > 2.7)) ||
                           (particle_type == 211 && medium_dQ_dx_1 > 2.0);
                };
                if (skip_with(pair_result_legacy.first) != skip_with(pair_result_proto.first))
                    g_pr33_audit.f1_mv_gate_flip++;
            }

            // Skip segments with certain particle types and high dQ/dx
            if ((particle_type == 11) ||
                (particle_type == 2212 && ((medium_dQ_dx_1 > 1.45 && pair_result.first <= 3) || medium_dQ_dx_1 > 2.7)) ||
                (particle_type == 211 && medium_dQ_dx_1 > 2.0)) {
                continue;
            }

            // doc sbnd_xin/docs/pr/40 round 5 F10: also skip a candidate whose
            // OWN geometry is long and straight (segment_is_straight_long_
            // track) -- none of the criteria above inspect straightness or
            // dQ/dx of a not-yet-PID'd sg (particle_type==0 skips every
            // branch above).  SBND evt 54341 seg 18005 (21.3 cm, 0.99
            // direct/arc ratio, particle_type still 0 at this point) is
            // exactly this gap -- see m_shower_connect_main_vertex_straight_
            // guard's comment in NeutrinoPatternBase.h.  false = legacy =
            // byte-identical.
            if (m_shower_connect_main_vertex_straight_guard && segment_is_straight_long_track(sg)) {
                continue;
            }

            // doc sbnd_xin/docs/pr/40 round 6 F13: also skip a pion whose far
            // vertex carries a charge-confirmed proton daughter -- "an
            // electron cannot father a proton" (round-3 protected_pion, same
            // predicate as Shower::update_particle_type's guard, PRShower.cxx
            // and the F5 relabel condition itself, NeutrinoPatternBase.cxx).
            // SBND evt 55715 seg 15005 (pi+, 6.1 cm -- UNDER the 10 cm floor,
            // so the F10 branch above cannot save it): once F11 declassifies
            // its daughter 15007, this function selects 15005 as the EM
            // candidate and force-sets it to pdg 11 at the accept site below
            // (probe tag "connecting_to_main_vertex"), reverting the pion
            // against proton daughter 15006.  The legacy pdg==211 branch two
            // blocks up only skips on HIGH dQ/dx (>2.0x MIP), which a real
            // pion fails.  false = legacy = byte-identical.
            if (m_shower_connect_protected_pion_guard && particle_type == 211 &&
                segment_has_proton_daughter(graph, sg, main_vertex, m_mip_dqdx_median)) {
                continue;
            }

            // doc pr/93 round 4 (straight_cont_cross_cluster): a stem the
            // demotion pass claimed as the head of a cross-cluster muon
            // chain must not seed a conn-1 shower here (measured on
            // 18264-137238: this pass re-captured the demoted 14.6cm stem).
            // Empty set when the knob is off => byte-identical.
            if (m_straight_cont_cross_cluster && m_sccc_shield_segs.count(sg)) {
                continue;
            }

            // Create a new shower
            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(main_vertex, 1);
            shower->set_start_segment(sg);
            pr93_probe_absorb_site("connecting_to_main_vertex", shower, m_shower_absorb_track_guard);
            shower->complete_structure_with_start_segment(nv_bridge_seed(used_segments), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
            pr84_probe_shower(shower, "connecting_to_main_vertex");

            // Single pass over shower edges: accumulate segment stats, vertex counts,
            // and flag_good_track together to avoid iterating edges twice.
            int n_tracks = 0;
            double total_length = 0;
            double max_seg_length = 0;
            SegmentPtr max_sg = nullptr;
            bool flag_good_track = false;
            std::map<VertexPtr, int> vtx_segment_count;

            // ordered_edges: total_length is an FP accumulation feeding the
            // length/track cuts below, so the walk order can move a cut
            // (doc pr/28 §15.2).
            for (auto edesc : ordered_edges(*shower, graph)) {
                auto sg1 = shower->view_graph()[edesc].segment;
                if (!sg1) continue;

                double length = segment_track_length(sg1);
                double medium_dQ_dx_norm = segment_median_dQ_dx(sg1) / m_mip_dqdx_median;

                n_tracks++;
                total_length += length;

                // Track max segment; use segment ID as tie-breaker for determinism
                if (length > max_seg_length || (length == max_seg_length && max_sg && sg1->id() < max_sg->id())) {
                    max_seg_length = length;
                    max_sg = sg1;
                }

                // Accumulate per-vertex edge counts (for n_two_vtx / n_multi_vtx below)
                auto seg_edesc = sg1->get_descriptor();
                VertexPtr sv = graph[boost::source(seg_edesc, graph)].vertex;
                VertexPtr tv = graph[boost::target(seg_edesc, graph)].vertex;
                if (sv) vtx_segment_count[sv]++;
                if (tv) vtx_segment_count[tv]++;

                // flag_good_track: only evaluate while still false (early skip once set)
                if (!flag_good_track &&
                    !seg_dir_weak(sg1) &&
                    (length > 3.6 * units::cm || (length > 2.4 * units::cm && medium_dQ_dx_norm > 2.5))) {

                    VertexPtr v1 = graph[boost::source(seg_edesc, graph)].vertex;
                    VertexPtr v2 = graph[boost::target(seg_edesc, graph)].vertex;

                    VertexPtr end_vertex = nullptr;
                    if (sg1->dirsign() == 1) {
                        // Forward: end is at back of fits
                        if (!sg1->fits().empty() && v1 && v2) {
                            const auto& fits = sg1->fits();
                            double dist1 = (fits.back().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.back().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    } else {
                        // Reversed or unknown direction: end is at front of fits
                        if (!sg1->fits().empty() && v1 && v2) {
                            const auto& fits = sg1->fits();
                            double dist1 = (fits.front().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.front().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    }

                    if (end_vertex && map_vertex_segments.find(end_vertex) != map_vertex_segments.end()) {
                        if (map_vertex_segments[end_vertex].size() > 1) {
                            bool flag_non_ele = false;
                            for (auto sg2 : map_vertex_segments[end_vertex]) {
                                if (sg2 == sg1) continue;
                                // Matches prototype get_flag_shower(): shower topology/trajectory OR pdg==11
                                bool sg2_is_shower = sg2->flags_any(SegmentFlags::kShowerTrajectory) ||
                                                     sg2->flags_any(SegmentFlags::kShowerTopology) ||
                                                     (sg2->has_particle_info() && sg2->particle_info() &&
                                                      std::abs(sg2->particle_info()->pdg()) == 11);
                                if (!sg2_is_shower) {
                                    flag_non_ele = true;
                                    break;  // no need to check remaining segments
                                }
                            }
                            if (!flag_non_ele && map_vertex_segments[end_vertex].size() <= 3)
                                flag_good_track = true;
                        } else {
                            flag_good_track = true;
                        }
                    }
                }
            }

            // Tally vertex connectivity for topology cuts
            int n_multi_vtx = 0;
            int n_two_vtx = 0;
            for (auto& [vtx, count] : vtx_segment_count) {
                if (count == 2) n_two_vtx++;
                else if (count > 2) n_multi_vtx++;
            }

            // Apply selection criteria
            if (!flag_good_track && n_multi_vtx > 0 && max_seg_length < 65 * units::cm &&
                ((total_length < n_tracks * 27 * units::cm && total_length < 85 * units::cm) ||
                 (total_length < n_tracks * 18 * units::cm && total_length < 95 * units::cm)) &&
                n_two_vtx < 3) {

                map_shower_max_sg[shower] = max_sg;
                double shower_len = shower->get_total_length();
                if (shower_len > max_length ||
                    (shower_len == max_length && max_shower &&
                     shower->start_segment() && max_shower->start_segment() &&
                     shower->start_segment()->id() < max_shower->start_segment()->id())) {
                    max_shower = shower;
                    max_length = shower_len;
                }
            }
        }

        // Process selected showers
        for (auto& [shower, max_sg] : map_shower_max_sg) {
            if (shower == max_shower) {
                // Convert to EM shower (particle type 11 = electron)
                SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_connecting_to_main_vertex: Convert EM shower {}", shower->start_segment()->id());
                if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                    // doc sbnd_xin/docs/pr/40 round 9 (round 7 candidate 2c,
                    // D1 re-target): F10 vetoes only at seed time on the
                    // seed's OWN geometry; a short anchor collinear with a
                    // long straight sibling passes F10 and reaches this
                    // accept-time write.  Decline the pdg write only --
                    // kAvoidMuonCheck, shower insertion and conflict
                    // deletion below are structure, not PID, and stay
                    // untouched.  Knob off => byte-identical.
                    if (m_shower_connect_start_seg_straight_guard &&
                        segment_is_straight_long_track_or_continuation(graph, shower->start_segment(), m_sfv_kink_max)) {
                        SPDLOG_LOGGER_DEBUG(s_log,
                            "pr40r9 start_seg_straight_guard: decline e- write seg={}",
                            shower->start_segment()->id());
                    }
                    else {
                        shower->start_segment()->particle_info()->set_pdg(11);
                        pr40_probe_setpdg(shower->start_segment(), 11, "NeutrinoShowerClustering.cxx:connecting_to_main_vertex");
                    }
                }

                // Set avoid muon check flag on max segment
                if (max_sg) {
                    max_sg->set_flags(SegmentFlags::kAvoidMuonCheck);
                }

                // Add shower to collection
                showers.insert(shower);

                // Collect all vertices in the new shower for conflict detection
                std::set<VertexPtr> shower_vertices;
                for (auto vdesc : shower->nodes()) {
                    auto vtx = shower->view_graph()[vdesc].vertex;
                    if (vtx) shower_vertices.insert(vtx);
                }

                for (auto& shower1 : showers) {
                    if (shower == shower1) continue;
                    auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
                    if (conn_type1 == 1 && start_vtx1 && shower_vertices.count(start_vtx1)) {
                        del_showers.insert(shower1);
                    }
                }
            }
        }

        // Delete conflicting showers
        for (auto& shower1 : del_showers) {
            auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
            if (start_vtx1 != main_vertex) {
                showers.erase(shower1);
            }
        }

        // Update shower maps
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
    }
}

void PatternAlgorithms::shower_clustering_with_nv_from_main_cluster(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters){
    if (!main_vertex || !main_cluster) return;
    
    // Build map_segment_vertices from graph (ordered for determinism)
    std::map<SegmentPtr, std::set<VertexPtr>> map_segment_vertices;
    std::vector<SegmentPtr> seg_order;
    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg) continue;

        auto source_vdesc = boost::source(e, graph);
        auto target_vdesc = boost::target(e, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;

        if (map_segment_vertices.find(seg) == map_segment_vertices.end()) {
            seg_order.push_back(seg);
        }
        if (v1) map_segment_vertices[seg].insert(v1);
        if (v2) map_segment_vertices[seg].insert(v2);
    }
    
    std::map<ShowerPtr, WireCell::Vector> map_shower_dir;
    std::map<ShowerPtr, double> map_shower_angle_offset;
    WireCell::Vector drift_dir(1, 0, 0);
    SegmentPtr max_length_segment = nullptr;
    
    // Step 1: Find the maximum length segment in main_cluster
    {
        double max_length = 0;
        for (auto seg : seg_order) {
            if (seg->cluster() != main_cluster) continue;
            double length = segment_track_length(seg);
            if (length > max_length && length > 6 * units::cm) {
                max_length = length;
                max_length_segment = seg;
            }
        }
    }
    
    // Step 2: Build map_shower_dir for showers in main_cluster
    for (auto seg : seg_order) {
        if (seg->cluster() != main_cluster) continue;
        if (map_segment_in_shower.find(seg) == map_segment_in_shower.end()) continue;
        
        ShowerPtr shower = map_segment_in_shower[seg];
        
        // Skip long muons.  doc pr/33 F2 (the inverted site): the prototype
        // reads the SHOWER's cached type here
        // (NeutrinoID_shower_clustering.h:1800,
        // fabs(shower->get_particle_type())==13); the port reads the start
        // segment's PDG.  Both readings feed the audit counter.
        int particle_type = 0;
        if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
            particle_type = shower->start_segment()->particle_info()->pdg();
        }
        {
            const bool skip_legacy = std::abs(particle_type) == 13;
            const bool skip_proto  = std::abs(shower->get_particle_type()) == 13;
            g_pr33_audit.f2_calls[1]++;
            if (skip_legacy != skip_proto) g_pr33_audit.f2_disagree[1]++;
            if (m_shower_pdg_from_shower_type ? skip_proto : skip_legacy) continue;
        }
        
        double total_length = shower->get_total_length();
        
        if (seg == shower->start_segment()) {
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            WireCell::Point start_point = start_vtx ? (start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point) : WireCell::Point(0, 0, 0);
            
            bool is_shower_topology = seg->flags_any(SegmentFlags::kShowerTopology);
            double medium_dQ_dx = segment_median_dQ_dx(seg, 0, 100); // matches prototype get_medium_dQ_dx(0,100)
            double seg_length = segment_track_length(seg);
            
            if (is_shower_topology || shower->get_num_segments() > 2 || 
                (medium_dQ_dx > m_mip_dqdx_median * 1.5 && seg_length > 0)) {
                if (seg_length > 10 * units::cm) {
                    WireCell::Vector dir_shower = segment_cal_dir_3vector(seg, start_point, 15 * units::cm);
                    map_shower_dir[shower] = dir_shower;
                } else {
                    WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                    map_shower_dir[shower] = dir_shower;
                }
            } else if (shower->get_num_segments() <= 2) {
                if (total_length > 30 * units::cm) {
                    if (seg_length > 10 * units::cm) {
                        WireCell::Vector dir_shower = segment_cal_dir_3vector(seg, start_point, 15 * units::cm);
                        map_shower_dir[shower] = dir_shower;
                    } else {
                        WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                        map_shower_dir[shower] = dir_shower;
                    }
                }
            }
            
            // Very large shower
            if (total_length > 100 * units::cm) {
                WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 60 * units::cm);
                map_shower_dir[shower] = dir_shower;
            }
            
            // Max length segment or segment connected to main_vertex
            if (seg == max_length_segment && map_shower_dir.find(shower) == map_shower_dir.end()) {
                WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                map_shower_dir[shower] = dir_shower;
            } else if (map_shower_dir.find(shower) == map_shower_dir.end() && 
                      map_segment_vertices[seg].find(main_vertex) != map_segment_vertices[seg].end()) {
                if (seg_length > 5 * units::cm) {
                    WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                    map_shower_dir[shower] = dir_shower;
                }
            }
            
            // Check if parallel to drift direction
            if (map_shower_dir.find(shower) != map_shower_dir.end()) {
                map_shower_angle_offset[shower] = 0;
                double angle_to_drift = std::abs(map_shower_dir[shower].dot(drift_dir) / (map_shower_dir[shower].magnitude() * drift_dir.magnitude()));
                angle_to_drift = std::acos(std::clamp(angle_to_drift, -1.0, 1.0)) / M_PI * 180.0;
                if (std::abs(angle_to_drift - 90) < 5) {
                    map_shower_dir[shower] = shower_cal_dir_3vector(*shower, start_point, 50 * units::cm);
                    map_shower_angle_offset[shower] = 5;
                }
            }
        }
    }
    
    // Step 3: If no shower directions found, try to add segments based on closest distance
    if (map_shower_dir.empty()) {
        // Sort showers by start_segment ID for deterministic iteration order
        std::vector<ShowerPtr> showers_sorted(showers.begin(), showers.end());
        std::sort(showers_sorted.begin(), showers_sorted.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
            auto sa = a->start_segment(); auto sb = b->start_segment();
            if (sa && sb) return sa->id() < sb->id();
            return (sa ? 1 : 0) < (sb ? 1 : 0);
        });

        std::map<ShowerPtr, double> map_shower_length;
        for (auto shower : showers_sorted) {
            map_shower_length[shower] = shower->get_total_length();
        }

        bool flag_continue = true;
        while (flag_continue) {
            flag_continue = false;
            for (auto seg1 : seg_order) {
                if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
                if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;

                double min_dis = 1e9;
                ShowerPtr min_shower = nullptr;
                double seg1_length = segment_track_length(seg1);  // compute once outside inner loop

                for (auto shower : showers_sorted) {
                    int particle_type = 0;
                    if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                        particle_type = shower->start_segment()->particle_info()->pdg();
                    }
                    if (particle_type == 13) continue;

                    if (seg1_length > 0.75 * map_shower_length[shower]) continue;

                    double dis = shower_get_closest_dis(*shower, seg1);
                    if (dis < min_dis) {
                        min_dis = dis;
                        min_shower = shower;
                    }
                }

                if (min_shower && min_dis < 3.5 * units::cm) {
                    pr93_probe_absorb_direct("pass3_proximity", min_shower, seg1);
                    min_shower->add_segment(seg1, true);
                    map_shower_length[min_shower] = min_shower->get_total_length();
                    flag_continue = true;
                }
            }
            // Termination guard.  Shower::add_segment() -> TrajectoryView::add_segment()
            // silently no-ops when the segment's descriptor is invalid
            // (clus/src/PRTrajectoryView.cxx:52 returns false, and the Shower
            // wrapper discards it), so seg1 never enters map_segment_in_shower,
            // the next pass re-selects it, and flag_continue is set forever.
            // map_segment_in_shower is rebuilt below purely from each shower's
            // m_edges (Shower::fill_maps() is the identity), and showers only
            // gain segments here, so its size grows by exactly one per genuine
            // add: no growth across a pass == no progress == stop.
            const size_t n_seg_mapped_before = map_segment_in_shower.size();
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
            if (flag_continue && map_segment_in_shower.size() <= n_seg_mapped_before) {
                SPDLOG_LOGGER_WARN(s_log,
                    "shower_clustering_with_nv_from_main_cluster: a pass requested >=1 add_segment "
                    "but registered none ({} segments mapped, unchanged); stopping to avoid an "
                    "unbounded loop", map_segment_in_shower.size());
                flag_continue = false;
            }
        }
    }

    if (map_shower_dir.empty()) return;

    // Step 4: Precompute shower direction info sorted by start_segment ID.
    // Avoids pointer-ordered map iteration and redundant per-segment lookups.
    struct ShowerDirInfo {
        ShowerPtr shower;
        WireCell::Vector dir;
        WireCell::Point start_point;
        double angle_offset;
    };
    std::vector<ShowerDirInfo> shower_dir_info;
    shower_dir_info.reserve(map_shower_dir.size());
    for (auto& [shower, dir] : map_shower_dir) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        if (!start_vtx) continue;
        WireCell::Point sp = start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point;
        auto it = map_shower_angle_offset.find(shower);
        double ao = (it != map_shower_angle_offset.end()) ? it->second : 0;
        shower_dir_info.push_back({shower, dir, sp, ao});
    }
    std::sort(shower_dir_info.begin(), shower_dir_info.end(), [](const ShowerDirInfo& a, const ShowerDirInfo& b) {
        auto sa = a.shower->start_segment(); auto sb = b.shower->start_segment();
        if (sa && sb) return sa->id() < sb->id();
        return (sa ? 1 : 0) < (sb ? 1 : 0);
    });

    // Examine other segments and add to showers based on angle and distance
    for (auto seg1 : seg_order) {
        if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
        if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;

        double min_dis = 1e9;
        ShowerPtr min_shower = nullptr;

        for (auto& info : shower_dir_info) {
            const auto& dir = info.dir;
            const auto& start_point = info.start_point;
            double angle_offset = info.angle_offset;

            // Get closest point on segment to shower start vertex
            auto [dist, closest_pt] = segment_get_closest_point(seg1, start_point);

            // Vector from shower start to closest point
            WireCell::Vector v1(closest_pt.x() - start_point.x(),
                               closest_pt.y() - start_point.y(),
                               closest_pt.z() - start_point.z());

            double angle = std::acos(std::clamp(dir.dot(v1) / (dir.magnitude() * v1.magnitude()), -1.0, 1.0));
            angle = angle / M_PI * 180.0;

            // Check angle and distance criteria
            if ((angle < 25.0 + angle_offset && dist < 80 * units::cm) ||
                (angle < 12.5 + angle_offset * 8 / 5 && dist < 130 * units::cm) ||
                (angle < 5 + angle_offset * 2 && dist < 200 * units::cm)) {

                double dis = std::pow(dist * std::cos(angle * M_PI / 180.0), 2) / std::pow(40 * units::cm, 2) +
                            std::pow(dist * std::sin(angle * M_PI / 180.0), 2) / std::pow(5 * units::cm, 2);

                if (dis < min_dis) {
                    min_dis = dis;
                    min_shower = info.shower;
                }
            }
        }

        if (min_shower) {
            // doc sbnd_xin/docs/pr/93 Cause D (SBND 18255-315167): this
            // direction-cone absorber has no PID or straightness check, and
            // with shower_absorb_unreachable_main ON (pr/65, SBND
            // production) a graph-unreachable MAIN-cluster segment is
            // eligible here -- 315167's 150.7cm score-0.10 proton (seg
            // 8001) was cone-absorbed into a 15.7cm EM stub's shower,
            // inflating it to "e- 1046.7 MeV".  When on, decline the
            // absorb for a confidently-PID'd non-electron straight-long
            // track -- the SAME predicate as the flood-fill's
            // guard_excludes (PRShower.cxx, pr/40 F12).  The declined
            // segment stays unclaimed for later passes.  C++ default
            // false => byte-identical.
            const bool cone_guard_excludes = m_shower_cone_absorb_guard &&
                seg1->has_particle_info() && seg1->particle_info() &&
                seg1->particle_info()->pdg() != 0 &&
                std::abs(seg1->particle_info()->pdg()) != 11 &&
                segment_track_length(seg1) > m_shower_pid_guard_min_len &&
                segment_is_straight_long_track(seg1);
            if (cone_guard_excludes) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr93 cone_absorb_guard: decline absorb seg={} pdg={} len={:.1f}cm",
                    pr91_seg_display_id(seg1), seg1->particle_info()->pdg(),
                    segment_track_length(seg1)/units::cm);
            }
            else {
                pr93_probe_absorb_direct("pass3_cone", min_shower, seg1);
                min_shower->add_segment(seg1, true);
            }
        }
    }
    
    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);

    // After the sweep, some clusters may be partially claimed: one segment added to a
    // shower (which puts the cluster into used_shower_clusters) but sibling segments in
    // the same cluster were skipped because they failed the angle/distance criteria.
    // Downstream functions (shower_clustering_in_other_clusters) skip entire clusters
    // that are already in used_shower_clusters, so those orphaned segments would never
    // be assigned to any shower.  Fix: add them to the same shower as their sibling.
    {
        std::map<Facade::Cluster*, ShowerPtr, ClusterPtrCmp> cluster_to_shower;
        for (auto& [seg, shower] : map_segment_in_shower) {
            if (seg && seg->cluster()) {
                cluster_to_shower.emplace(seg->cluster(), shower);
            }
        }
        bool changed = false;
        for (auto seg1 : seg_order) {
            if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
            if (map_segment_in_shower.count(seg1)) continue;
            auto it = cluster_to_shower.find(seg1->cluster());
            if (it == cluster_to_shower.end()) continue;
            // doc sbnd_xin/docs/pr/93 Cause D, second seat (SBND
            // 18255-315167): after the cone absorber declines a segment,
            // this sibling backfill would re-adopt it (its cluster is
            // "partially claimed" by the shower's own stub members) AND
            // force-relabel it pdg 11 below -- the same knob and predicate
            // decline both.  C++ default false => byte-identical.
            if (m_shower_cone_absorb_guard &&
                seg1->has_particle_info() && seg1->particle_info() &&
                seg1->particle_info()->pdg() != 0 &&
                std::abs(seg1->particle_info()->pdg()) != 11 &&
                segment_track_length(seg1) > m_shower_pid_guard_min_len &&
                segment_is_straight_long_track(seg1)) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr93 cone_absorb_guard: decline sibling backfill seg={} pdg={} len={:.1f}cm",
                    pr91_seg_display_id(seg1), seg1->particle_info()->pdg(),
                    segment_track_length(seg1)/units::cm);
                continue;
            }
            pr93_probe_absorb_direct("pass3_cluster_map", it->second, seg1);
            it->second->add_segment(seg1, true);
            if (seg1->has_particle_info() && seg1->particle_info()) {
                seg1->particle_info()->set_pdg(11);
                pr40_probe_setpdg(seg1, 11, "NeutrinoShowerClustering.cxx:orphan_sibling_adopt");
                seg1->particle_info()->set_mass(0.511 * units::MeV);
            }
            changed = true;
        }
        if (changed) {
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
        }
    }
}

// doc sbnd_xin/docs/pr/40 round 9 B2 -- attempt the cross-cluster track
// bridge for one accepted Step-5 candidate of
// shower_clustering_with_nv_from_vertices.  Returns true iff the bridge was
// built (or already existed), in which case the caller must NOT create the
// conn-2 shower for this cluster.  Preconditions: the directional match has
// already accepted (angle gates), sg1 sits in `cluster` != main_cluster,
// `vertex` is a main-cluster vertex, `point` is on sg1.
//
// Owner directive (round 8): a long straight track pointing at the
// determined vertex across a small SP gap is one muon -- modify the GRAPH,
// after the neutrino vertex is determined.  The gap has no charge, so the
// bridge is a straight 2-point zero-charge segment (no do_rough_path
// routing -- that Dijkstra is single-cluster and would silently return an
// in-cluster detour; no refit -- do_multi_tracking is nonlocal and has no
// charge to fit).  sg1 keeps its own pdg: the F10 G2 lesson says ordinary
// track PID owns that decision once the electron overwrite is out of the
// way.
// doc sbnd_xin/docs/pr/93 round 4 -- the connect/register tail of
// nv_bridge_track (pr/40 round 9 B2), factored out verbatim so the sccc
// bridge replay (straight_cont_cross_cluster + sccc_bridge_body) can reuse
// it: find-or-create the synthetic main-cluster zero-charge bridge segment
// (a REAL graph edge -- the load-bearing piece for BOTH the PF and kine
// BFSes), shield the bridge + the rescued cluster's own segments from every
// shower flood-fill/absorber (the pr/54 lesson), and record the cluster id
// for shower_clustering_in_other_clusters and the PF-side
// pf_track_bridged_clusters gate (transported via TrackFitting).  Returns
// the bridge segment, or nullptr on failure.  Byte-identical refactor: the
// only caller-visible difference from the inlined original is that
// half1/half2 shielding stays at the nv_bridge_track call site.
SegmentPtr PatternAlgorithms::nv_bridge_connect(Graph& graph, Facade::Cluster* main_cluster,
                                                Facade::Cluster* cluster, VertexPtr vertex,
                                                VertexPtr far_vtx,
                                                const std::vector<SegmentPtr>& cluster_segs,
                                                TrackFitting& track_fitter,
                                                IDetectorVolumes::pointer dv)
{
    if (!far_vtx || far_vtx == vertex) return nullptr;

    SegmentPtr bridge = find_segment(graph, vertex, far_vtx);
    if (!bridge) {
        // Stamped main_cluster: same_cluster(bridge) holds in
        // fill_bee_pf_tree and the absorber guards treat it as claimed.
        std::vector<Facade::geo_point_t> pp{vertex->wcpt().point, far_vtx->wcpt().point};
        bridge = create_segment_for_cluster(*main_cluster, dv, pp, 0);
        if (!bridge) return nullptr;
        // Two synthetic zero-charge fits: downstream consumers (nusel
        // taggers) call fits().front()/back() unguarded, and an empty-fits
        // segment adjacent to the main vertex would be a novel state.
        WireCell::Point p0 = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
        WireCell::Point p1 = far_vtx->fit().valid() ? far_vtx->fit().point : far_vtx->wcpt().point;
        const double glen = (p1 - p0).magnitude();
        PR::Fit f0, f1;
        f0.point = p0; f0.index = 0; f0.dQ = 0; f0.dx = glen / 2;
        f1.point = p1; f1.index = 1; f1.dQ = 0; f1.dx = glen / 2;
        bridge->fits({f0, f1});
        add_segment(graph, bridge, vertex, far_vtx);
    }

    // Shield the bridge and the rescued cluster's own segments from every
    // shower flood-fill/absorber (the pr/54 lesson), and record the cluster
    // id both for shower_clustering_in_other_clusters and for the PF-side
    // pf_track_bridged_clusters gate (transported via TrackFitting).
    m_nv_bridge_shield_segs.insert(bridge);
    for (const auto& s : cluster_segs) m_nv_bridge_shield_segs.insert(s);
    m_nv_bridge_cluster_ids.insert(cluster->get_cluster_id());
    track_fitter.add_bridged_cluster_id(cluster->get_cluster_id());

    return bridge;
}

bool PatternAlgorithms::nv_bridge_track(Graph& graph, Facade::Cluster* main_cluster,
                                        Facade::Cluster* cluster, SegmentPtr sg1, VertexPtr vertex,
                                        const WireCell::Point& point,
                                        const std::vector<SegmentPtr>& cluster_segs,
                                        TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
                                        const Clus::ParticleDataSet::pointer& particle_data,
                                        const IRecombinationModel::pointer& recomb_model)
{
    // Straightness on the PRE-break candidate (or its collinear
    // continuation -- 286906's 8.7 cm anchor + 127 cm body).
    if (!segment_is_straight_long_track_or_continuation(graph, sg1, m_sfv_kink_max)) return false;

    // Exact steiner-cloud closest approach; -1 = unmeasurable = no bridge.
    const double gap = cluster_steiner_gap(*cluster, *main_cluster);
    if (gap < 0 || gap >= m_shower_nv_bridge_max_gap) return false;

    // Far-side vertex: break at an interior point (the same break the
    // legacy shower path performs), else the nearest endpoint vertex.
    VertexPtr far_vtx = nullptr;
    SegmentPtr half1 = nullptr, half2 = nullptr;
    const auto& fits = sg1->fits();
    const bool at_end = !fits.empty() &&
        ((fits.front().point - point).magnitude() < 0.01 * units::cm ||
         (fits.back().point - point).magnitude() < 0.01 * units::cm);
    if (!fits.empty() && !at_end) {
        auto [success, seg_pair, new_vtx] = break_segment(graph, sg1, point, particle_data,
                                                          recomb_model, dv, 1e9 * units::cm,
                                                          m_break_seg_orient);
        if (success && new_vtx) {
            far_vtx = new_vtx;
            half1 = seg_pair.first;
            half2 = seg_pair.second;
        }
    }
    if (!far_vtx) {
        auto [v1, v2] = find_vertices(graph, sg1);
        auto vpt = [](VertexPtr v) { return v->fit().valid() ? v->fit().point : v->wcpt().point; };
        if (v1 && v2) {
            far_vtx = ((vpt(v1) - point).magnitude() <= (vpt(v2) - point).magnitude()) ? v1 : v2;
        }
        else {
            far_vtx = v1 ? v1 : v2;
        }
    }
    SegmentPtr bridge = nv_bridge_connect(graph, main_cluster, cluster, vertex, far_vtx,
                                          cluster_segs, track_fitter, dv);
    if (!bridge) return false;
    if (half1) m_nv_bridge_shield_segs.insert(half1);
    if (half2) m_nv_bridge_shield_segs.insert(half2);

    SPDLOG_LOGGER_DEBUG(s_log,
        "pr40r9 nv_bridge: cluster {} -> main {} gap={:.2f}cm sg1={} len={:.2f}cm broke={}",
        cluster->get_cluster_id(), main_cluster->get_cluster_id(), gap / units::cm,
        sg1->id(), segment_track_length(sg1) / units::cm, half1 != nullptr);
    return true;
}

void PatternAlgorithms::shower_clustering_with_nv_from_vertices(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    if (!main_vertex || !main_cluster) return;
    
    // Build map_cluster_segments, map_segment_cluster, and seg_order.
    // seg_order holds unique segments in deterministic graph-index order for
    // all iteration loops, eliminating pointer-address-based ordering.
    std::map<Facade::Cluster*, std::vector<SegmentPtr>> map_cluster_segments;
    std::map<SegmentPtr, Facade::Cluster*> map_segment_cluster;
    std::vector<SegmentPtr> seg_order;

    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg || !seg->cluster()) continue;

        if (map_segment_cluster.emplace(seg, seg->cluster()).second) {
            map_cluster_segments[seg->cluster()].push_back(seg);
            seg_order.push_back(seg);
        }
    }
    
    // Step 1: Build map_cluster_center_point
    std::map<Facade::Cluster*, std::pair<WireCell::Point, double>> map_cluster_center_point;
    
    for (auto cluster : other_clusters) {
        auto it1 = map_cluster_segments.find(cluster);
        if (it1 == map_cluster_segments.end()) continue;
        
        double acc_length = 0;
        double acc_length1 = 0;
        WireCell::Point p(0, 0, 0);
        int np = 0;
        
        for (auto seg : it1->second) {
            if (map_segment_in_shower.find(seg) != map_segment_in_shower.end()) continue;
            
            bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) || seg->flags_any(SegmentFlags::kShowerTopology);
            int particle_type = 0;
            if (seg->has_particle_info() && seg->particle_info()) {
                particle_type = seg->particle_info()->pdg();
            }
            // doc pr/33 F4: the prototype's get_flag_shower() carries a third
            // disjunct this mirror dropped (ProtoSegment.cxx:1305-1312,
            // fabs(particle_type)==11).  It feeds the acc_length read below
            // AND the acc_length1 read (a pdg -11 segment lands in
            // acc_length1 without it), i.e. one gate at the center-point cut.
            if (!is_shower && std::abs(particle_type) == 11) {
                g_pr33_audit.f4_flip++;
                if (m_shower_flag_pdg_electron) is_shower = true;
            }

            if (is_shower || particle_type == 0 ||
                ((std::abs(particle_type) == 13 || std::abs(particle_type) == 211) && seg_dir_weak(seg))) {
                double length = segment_track_length(seg);
                acc_length += length;
                
                const auto& fits = seg->fits();
                for (const auto& fit : fits) {
                    p.set(p.x() + fit.point.x(), p.y() + fit.point.y(), p.z() + fit.point.z());
                    np++;
                }
            }
            
            if (particle_type != 11 && !is_shower) {
                acc_length1 += segment_track_length(seg);
            }
        }
        
        if ((acc_length > 1.0 * units::cm && acc_length >= acc_length1) || acc_length > 10 * units::cm) {
            if (np > 0) {
                p.set(p.x() / np, p.y() / np, p.z() / np);
            }
            map_cluster_center_point[cluster] = std::make_pair(p, acc_length);
        }
    }
    
    // Step 2: List main cluster vertices
    std::vector<VertexPtr> main_cluster_vertices;
    for (auto v : ordered_nodes(graph)) {
        VertexPtr vtx = graph[v].vertex;
        if (!vtx || !vtx->cluster() || vtx->cluster() != main_cluster) continue;

        if (vtx != main_vertex) {
            if (vertices_in_long_muon.find(vtx) != vertices_in_long_muon.end()) continue;
            if (map_vertex_in_shower.find(vtx) != map_vertex_in_shower.end()) continue;
        }
        main_cluster_vertices.push_back(vtx);
    }
    
    // Step 3: Analyze each cluster against main cluster vertices
    std::map<Facade::Cluster*, cluster_point_info> map_cluster_pi;
    
    // Get wpid_params for point cloud creation
    auto* grouping = main_cluster->grouping();
    if (!grouping) return;
    const auto& wpids = grouping->wpids();
    std::map<WirePlaneId, std::tuple<Facade::geo_point_t, double, double, double>> wpid_params;
    std::map<WirePlaneId, std::pair<Facade::geo_point_t, double>> wpid_U_dir, wpid_V_dir, wpid_W_dir;
    std::set<int> apas;
    Facade::compute_wireplane_params(wpids, dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
    
    for (auto cluster : other_clusters) {
        // doc pr/93 round 4 (straight_cont_cross_cluster): a cluster the
        // sccc replay already bridged as a muon-body chain is a TRACK
        // continuation, not an EM-shower candidate -- without this skip,
        // Path C broke the bridged cluster's entry stub off as a conn-2
        // electron (18264-137238).  Scoped to replay-bridged ids only:
        // nv_bridge_track's own Step-5 population keeps legacy behavior.
        // Empty set when the knob is off => byte-identical.
        if (m_straight_cont_cross_cluster &&
            m_sccc_bridged_cluster_ids.count(cluster->get_cluster_id())) continue;
        auto cpi = map_cluster_center_point.find(cluster);
        if (cpi == map_cluster_center_point.end()) continue;
        auto& center_pair = cpi->second;
        WireCell::Point center_p = center_pair.first;
        
        // Create point cloud from cluster segments
        std::vector<std::pair<Facade::geo_point_t, WirePlaneId>> point_plane_pairs;
        for (auto seg : map_cluster_segments[cluster]) {
            for (const auto& fit : seg->fits()) {
                WirePlaneId wpid = dv->contained_by(fit.point);
                point_plane_pairs.emplace_back(fit.point, wpid);
            }
        }
        auto dpc_points = Facade::make_points_direct(cluster, dv, wpid_params, point_plane_pairs, false);
        auto pcloud = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
        pcloud->add_points(dpc_points);
        
        cluster_point_info min_pi;
        min_pi.cluster = cluster;
        min_pi.min_angle = 90;
        min_pi.min_dis = 1e9;
        min_pi.min_vertex = nullptr;
        
        cluster_point_info main_pi;
        main_pi.cluster = cluster;
        main_pi.min_vertex = main_vertex;
        
        std::vector<double> query(3);
        for (auto vtx : main_cluster_vertices) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;

            // Get closest point using KD-tree
            auto& kd3d = pcloud->kd3d();
            query[0] = vtx_pt.x(); query[1] = vtx_pt.y(); query[2] = vtx_pt.z();
            auto results = kd3d.knn(1, query);
            
            if (results.empty()) continue;
            
            const size_t idx = results[0].first;
            const double dis = std::sqrt(results[0].second);  // KD-tree returns squared distance
            const WireCell::Point closest_pt = pcloud->point3d(idx);
            
            WireCell::Vector v1(closest_pt.x() - vtx_pt.x(),
                               closest_pt.y() - vtx_pt.y(),
                               closest_pt.z() - vtx_pt.z());
            WireCell::Vector v2(center_p.x() - closest_pt.x(),
                               center_p.y() - closest_pt.y(),
                               center_p.z() - closest_pt.z());
            
            WireCell::Point near_center = pcloud->get_center_point_radius(closest_pt, 2 * units::cm);
            WireCell::Vector v3(near_center.x() - closest_pt.x(),
                               near_center.y() - closest_pt.y(),
                               near_center.z() - closest_pt.z());
            
            double angle = std::acos(std::clamp(v1.dot(v2) / (v1.magnitude() * v2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            double angle3 = std::acos(std::clamp(v1.dot(v3) / (v1.magnitude() * v3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            
            if (angle < 30 || (dis < 5 * units::cm && angle < 45)) {
                angle = std::min(angle, angle3);
            }
            
            if (angle < 7.5) {
                if (dis * std::sin(angle / 180.0 * M_PI) < min_pi.min_dis * std::sin(min_pi.min_angle / 180.0 * M_PI) && angle < 90) {
                    min_pi.min_angle = angle;
                    min_pi.min_dis = dis;
                    min_pi.min_vertex = vtx;
                    min_pi.min_point = closest_pt;
                }
            } else {
                if (angle < min_pi.min_angle) {
                    min_pi.min_angle = angle;
                    min_pi.min_dis = dis;
                    min_pi.min_vertex = vtx;
                    min_pi.min_point = closest_pt;
                }
            }
            
            if (vtx == main_vertex) {
                main_pi.min_angle = angle;
                main_pi.min_dis = dis;
                main_pi.min_point = closest_pt;
            }
        }
        
        if (!min_pi.min_vertex) {
            min_pi.min_angle = main_pi.min_angle;
            min_pi.min_vertex = main_vertex;
            min_pi.min_point = main_pi.min_point;
            min_pi.min_dis = main_pi.min_dis;
        }
        
        WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
        WireCell::Point min_vtx_pt = min_pi.min_vertex->fit().valid() ? min_pi.min_vertex->fit().point : min_pi.min_vertex->wcpt().point;
        double vtx_dis = (main_vtx_pt - min_vtx_pt).magnitude();
        
        if (main_pi.min_angle < min_pi.min_angle + 3 && main_pi.min_dis < min_pi.min_dis * 1.2 &&
            (min_pi.min_angle > 0.9 * main_pi.min_angle || vtx_dis < 1.5 * units::cm)) {
            map_cluster_pi[cluster] = main_pi;
        } else {
            map_cluster_pi[cluster] = min_pi;
        }
    }
    
    // Step 4: Sort by distance (iterate other_clusters for deterministic pre-sort order)
    std::vector<cluster_point_info> vec_pi;
    for (auto cluster : other_clusters) {
        auto it = map_cluster_pi.find(cluster);
        if (it != map_cluster_pi.end()) {
            vec_pi.push_back(it->second);
        }
    }
    std::sort(vec_pi.begin(), vec_pi.end(), sortbydis);
    
    ClusterVertexMap map_cluster_associated_vertex;
    for (const auto& pi : vec_pi) {
        if (pi.min_angle < 10) {
            map_cluster_associated_vertex[pi.cluster] = pi.min_vertex;
        }
    }
    
    // Step 5: Process each cluster
    for (const auto& pi : vec_pi) {
        Facade::Cluster* cluster = pi.cluster;
        VertexPtr vertex = pi.min_vertex;
        WireCell::Point point = pi.min_point;
        SegmentPtr sg1 = nullptr;
        double angle = pi.min_angle;
        
        if (angle > 50 && pi.min_dis > 6 * units::cm) continue;
        if (angle > 60) continue;
        
        // Find segment at the point
        for (auto seg : map_cluster_segments[cluster]) {
            if (map_segment_in_shower.find(seg) != map_segment_in_shower.end()) continue;
            auto [dis, closest_pt] = segment_get_closest_point(seg, point);
            if (dis < 0.01 * units::cm) {
                sg1 = seg;
                break;
            }
        }
        
        if (!sg1) continue;

        // doc sbnd_xin/docs/pr/40 round 9 B2: a straight long track in
        // another cluster whose steiner-cloud gap to the main cluster is
        // below m_shower_nv_bridge_max_gap is a broken-off track (an SP
        // hole), not a shower.  Bridge the graph instead of fabricating a
        // conn-2 electron.  Knob off => predicate never evaluated =>
        // byte-identical legacy path below.
        if (m_shower_nv_bridge_track && sg1->cluster() != main_cluster &&
            nv_bridge_track(graph, main_cluster, cluster, sg1, vertex, point,
                            map_cluster_segments[cluster], track_fitter, dv,
                            particle_data, recomb_model)) {
            continue;
        }

        // Create new shower
        ShowerPtr shower = std::make_shared<Shower>(graph);
        shower->set_start_vertex(vertex, 2);
        showers.insert(shower);
        
        // Check if point is at segment endpoint
        const auto& fits = sg1->fits();
        if (!fits.empty()) {
            double dist_front = (fits.front().point - point).magnitude();
            double dist_back = (fits.back().point - point).magnitude();
            
            if (dist_front < 0.01 * units::cm || dist_back < 0.01 * units::cm) {
                shower->set_start_segment(sg1, true);
            } else {
                // Break segment at point
                auto [success, seg_pair, new_vtx] = break_segment(graph, sg1, point, particle_data, recomb_model, dv,
                                                                  1e9*units::cm, m_break_seg_orient);
                
                if (!success || !new_vtx) {
                    shower->set_start_segment(sg1, true);
                } else {
                    // Determine which segment to use based on direction to vertex
                    WireCell::Point vtx_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
                    WireCell::Vector v3(point.x() - vtx_pt.x(), point.y() - vtx_pt.y(), point.z() - vtx_pt.z());

                    WireCell::Vector v1 = segment_cal_dir_3vector(seg_pair.first, point, 5 * units::cm);
                    WireCell::Vector v2 = segment_cal_dir_3vector(seg_pair.second, point, 5 * units::cm);

                    double angle1 = std::acos(std::clamp(v1.dot(v3) / (v1.magnitude() * v3.magnitude()), -1.0, 1.0));
                    double angle2 = std::acos(std::clamp(v2.dot(v3) / (v2.magnitude() * v3.magnitude()), -1.0, 1.0));

                    if (angle1 < angle2) {
                        shower->set_start_segment(seg_pair.first, true);
                    } else {
                        shower->set_start_segment(seg_pair.second, true);
                    }
                }
            }
        } else {
            shower->set_start_segment(sg1);
        }
        
        pr84_probe_shower(shower, "nv_from_vertices_break");

        // doc pr/40 round 9: set when the Part-A guard below declines the
        // electron overwrite; consulted at the update_particle_type call
        // for this shower (which would otherwise redo the 13->11 via the
        // majority vote, PRShower.cxx:981-1009).
        bool sfv_guard_fired = false;

        // Set direction based on vertex proximity
        auto start_seg = shower->start_segment();
        if (start_seg) {
            const auto& seg_fits = start_seg->fits();
            if (!seg_fits.empty()) {
                WireCell::Point vtx_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
                double dis1 = (vtx_pt - seg_fits.front().point).magnitude();
                double dis2 = (vtx_pt - seg_fits.back().point).magnitude();
                
                if (dis1 < dis2) {
                    start_seg->dirsign(1);
                } else {
                    start_seg->dirsign(-1);
                }
            }
            
            // Set particle type to electron if needed
            int pdg = 0;
            if (start_seg->has_particle_info() && start_seg->particle_info()) {
                pdg = start_seg->particle_info()->pdg();
            }
            if (pdg == 0 || std::abs(pdg) == 13) {
                // doc sbnd_xin/docs/pr/40 round 9 (round 8 Part A): decline
                // the unconditional track->e- overwrite when the anchor (or
                // its collinear continuation across the shared vertex --
                // PATH C hands this block a broken sub-10cm half; 286906:
                // 8.68cm anchor at 4.9deg to the 126.89cm body) is a
                // straight long track.  Knob off => predicate never
                // evaluated => byte-identical.  The shower itself is kept:
                // with a track-typed start segment append_pseudo_shower
                // renders the 2112 neutron carrier, the owner-desired
                // "hadron -> neutron" display for a non-EM conn-2 object.
                if (m_shower_connect_from_vertices_straight_guard &&
                    segment_is_straight_long_track_or_continuation(graph, start_seg, m_sfv_kink_max)) {
                    sfv_guard_fired = true;
                    m_sfv_declined_anchors.insert(start_seg);
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr40r9 sfv_straight_guard: decline e- write seg={} clus={} pdg={}",
                        start_seg->id(),
                        start_seg->cluster() ? start_seg->cluster()->get_cluster_id() : -1,
                        pdg);
                }
                else {
                auto four_momentum = segment_cal_4mom(start_seg, 11, particle_data, recomb_model, m_mip_dqdx);

                // Create ParticleInfo for electron
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    11,                                          // electron PDG
                    particle_data->get_particle_mass(11),       // electron mass
                    particle_data->pdg_to_name(11),             // "electron"
                    four_momentum                                // 4-momentum
                );

                // Store particle info in start_segment
                start_seg->particle_info(pinfo);
                }
            }
        }
        
        // Complete shower structure
        IndexedSegmentSet used_segments;
        pr93_probe_absorb_site("from_vertices", shower, m_shower_absorb_track_guard);
        shower->complete_structure_with_start_segment(nv_bridge_seed(used_segments), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
        
        // Calculate shower direction
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        WireCell::Point start_pt = start_vtx ? (start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point) : point;
        
        WireCell::Vector dir_shower = segment_cal_dir_3vector(shower->start_segment(), point, 15 * units::cm);
        WireCell::Vector dir_main(point.x() - start_pt.x(), point.y() - start_pt.y(), point.z() - start_pt.z());
        
        if (std::acos(std::clamp(dir_shower.dot(dir_main) / (dir_shower.magnitude() * dir_main.magnitude()), -1.0, 1.0)) / M_PI * 180.0 > 30) {
            auto [_, test_p] = shower_get_closest_point(*shower, start_pt);
            dir_shower = shower_cal_dir_3vector(*shower, test_p, 30 * units::cm);
        }
        if (dir_shower.magnitude() < 0.001) dir_shower = dir_main;
        
        // Add segments from other clusters.
        // Cache the shower start front point to avoid re-evaluating per segment.
        const WireCell::Point shower_start_front = shower->start_segment()->fits().front().point;
        for (auto seg1 : seg_order) {
            if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
            if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
            if (seg1->cluster() == shower->start_segment()->cluster()) continue;

            auto it1 = map_cluster_associated_vertex.find(seg1->cluster());

            auto [pair_dis, pair_point] = segment_get_closest_point(seg1, start_pt);
            WireCell::Vector v1(pair_point.x() - start_pt.x(), pair_point.y() - start_pt.y(), pair_point.z() - start_pt.z());
            WireCell::Vector v2(pair_point.x() - point.x(), pair_point.y() - point.y(), pair_point.z() - point.z());

            double angle_v1 = std::acos(std::clamp(dir_shower.dot(v1) / (dir_shower.magnitude() * v1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            double angle_v2 = std::acos(std::clamp(dir_shower.dot(v2) / (dir_shower.magnitude() * v2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;

            // Filter early before the two expensive KD-tree calls below
            if (angle_v2 > 30) continue;

            double tmp_shower_dis = segment_get_closest_point(seg1, shower_start_front).first;
            double close_shower_dis = shower_get_closest_dis(*shower, seg1);

            if ((angle_v1 < 25 && (pair_dis < 80 * units::cm || close_shower_dis < 25 * units::cm)) ||
                (angle_v2 < 25 && (tmp_shower_dis < 40 * units::cm || close_shower_dis < 25 * units::cm)) ||
                (angle_v1 < 12.5 && (pair_dis < 120 * units::cm || close_shower_dis < 40 * units::cm)) ||
                (angle_v2 < 12.5 && (tmp_shower_dis < 80 * units::cm || close_shower_dis < 40 * units::cm))) {

                if (it1 != map_cluster_associated_vertex.end() && seg1->cluster() != shower->start_segment()->cluster()) {
                    if (it1->second != vertex) {
                        double dis1 = shower_get_dis(*shower, seg1);
                        if (dis1 > 25 * units::cm && dis1 > pair_dis * 0.4) {
                            continue;
                        }
                    }
                }
                pr93_probe_absorb_direct("pass4_angle", shower, seg1);
                shower->add_segment(seg1, true);
            }
        }

        // Update particle type
        // doc pr/40 round 9 (D3): skipped when the Part-A guard declined the
        // electron write above -- the vote counts a pdg-13 start segment
        // toward shower_length and would silently redo the 13->11.
        if (!sfv_guard_fired) {
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
        }

        bool tmp_flag = (shower->start_vertex() == main_vertex);
        SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_with_nv_from_vertices: Separated shower: {} {} {} {} {}",
            shower->start_segment()->cluster()->get_cluster_id() * 1000 + shower->start_segment()->id(),
            (shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info() ? shower->start_segment()->particle_info()->pdg() : 0),
            shower->get_num_segments(), tmp_flag, pi.min_dis / units::cm);

        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);

        // Iteratively add segments based on distance.
        // Sort showers once per iteration for deterministic ordering.
        {
            std::vector<ShowerPtr> sorted_showers(showers.begin(), showers.end());
            std::sort(sorted_showers.begin(), sorted_showers.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
                auto* sa = a->start_segment().get();
                auto* sb = b->start_segment().get();
                if (!sa || !sb) return sa < sb;
                int cid_a = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
                int cid_b = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
                if (cid_a != cid_b) return cid_a < cid_b;
                return sa->id() < sb->id();
            });

            std::map<ShowerPtr, double> map_shower_length;
            std::map<ShowerPtr, WireCell::Vector> map_shower_dir;

            for (auto shower1 : sorted_showers) {
                map_shower_length[shower1] = shower1->get_total_length();
                auto [start_vtx1, _] = shower1->get_start_vertex_and_type();
                WireCell::Point start_pt1 = start_vtx1 ? (start_vtx1->fit().valid() ? start_vtx1->fit().point : start_vtx1->wcpt().point) : WireCell::Point(0, 0, 0);
                auto [__, test_p] = shower_get_closest_point(*shower1, start_pt1);
                map_shower_dir[shower1] = shower_cal_dir_3vector(*shower1, test_p, 30 * units::cm);
            }

            bool flag_continue = true;
            while (flag_continue) {
                flag_continue = false;
                for (auto seg1 : seg_order) {
                    if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
                    if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;

                    double min_dis = 1e9;
                    ShowerPtr min_shower = nullptr;

                    for (auto shower1 : sorted_showers) {
                        if (segment_track_length(seg1) > 0.75 * map_shower_length[shower1]) continue;

                        auto [start_vtx1, _] = shower1->get_start_vertex_and_type();
                        WireCell::Point start_point1 = start_vtx1 ? (start_vtx1->fit().valid() ? start_vtx1->fit().point : start_vtx1->wcpt().point) : WireCell::Point(0, 0, 0);
                        auto [__, test_p] = segment_get_closest_point(seg1, start_point1);
                        auto [___, test_p1] = shower_get_closest_point(*shower1, start_point1);

                        WireCell::Vector tmp_dir(test_p.x() - test_p1.x(), test_p.y() - test_p1.y(), test_p.z() - test_p1.z());
                        double angle = std::acos(std::clamp(tmp_dir.dot(map_shower_dir[shower1]) / (tmp_dir.magnitude() * map_shower_dir[shower1].magnitude()), -1.0, 1.0)) / M_PI * 180.0;

                        double dis = shower_get_closest_dis(*shower1, seg1);
                        if (dis < min_dis && angle < 45) {
                            min_dis = dis;
                            min_shower = shower1;
                        }
                    }

                    if (min_shower && min_dis < 3.5 * units::cm) {
                        pr93_probe_absorb_direct("pass4_proximity", min_shower, seg1);
                        min_shower->add_segment(seg1, true);
                        map_shower_length[min_shower] = min_shower->get_total_length();
                        flag_continue = true;
                    }
                }
                // Termination guard -- see the identical guard in
                // shower_clustering_with_nv_in_main_cluster above.  This is the
                // site that hung for 8h17m on SBND MCP2025C evt 352365 (run
                // 18255/1): 100% CPU, byte-flat RSS, no I/O, gdb caught it
                // looping through segment_get_closest_point below.
                const size_t n_seg_mapped_before = map_segment_in_shower.size();
                update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
                if (flag_continue && map_segment_in_shower.size() <= n_seg_mapped_before) {
                    SPDLOG_LOGGER_WARN(s_log,
                        "shower_clustering_with_nv_from_vertices: a pass requested >=1 add_segment "
                        "but registered none ({} segments mapped, unchanged); stopping to avoid an "
                        "unbounded loop", map_segment_in_shower.size());
                    flag_continue = false;
                }
            }
        }
    }
    
    SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_with_nv_from_vertices: With separated-cluster shower: {}", showers.size());
}

void PatternAlgorithms::examine_merge_showers(IndexedShowerSet& showers, VertexPtr main_vertex, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex || map_vertex_to_shower.find(main_vertex) == map_vertex_to_shower.end()) {
        return;
    }

    auto& main_vertex_showers = map_vertex_to_shower[main_vertex];

    // Sort showers by start-segment graph index for deterministic iteration order,
    // avoiding dependence on pointer addresses (same pattern as ordered_edges/ordered_nodes).
    auto seg_index = [&](const ShowerPtr& s) -> size_t {
        auto seg = s->start_segment();
        return (seg && seg->descriptor_valid()) ? graph[seg->get_descriptor()].index
                                                : std::numeric_limits<size_t>::max();
    };
    std::vector<ShowerPtr> sorted_showers(main_vertex_showers.begin(), main_vertex_showers.end());
    std::sort(sorted_showers.begin(), sorted_showers.end(),
              [&](const ShowerPtr& a, const ShowerPtr& b) { return seg_index(a) < seg_index(b); });

    // Pre-classify showers and pre-compute type-2 directions once (avoids recomputing dir2
    // redundantly inside the inner loop for each type-1 shower that encounters it).
    std::vector<ShowerPtr> type1_showers, type2_showers;
    std::unordered_map<Shower*, WireCell::Vector> type2_dirs;
    for (auto& shower : sorted_showers) {
        // doc pr/33 F2: the prototype reads the start segment's PDG here
        // (NeutrinoID_shower_clustering.h:387/:394, exact ==13); the port
        // reads the shower's cached type.  Both feed the audit counter.
        {
            int type_startseg = 0;
            if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                type_startseg = shower->start_segment()->particle_info()->pdg();
            }
            const bool skip_legacy = shower->get_particle_type() == 13;
            const bool skip_proto  = type_startseg == 13;
            g_pr33_audit.f2_calls[2]++;
            if (skip_legacy != skip_proto) g_pr33_audit.f2_disagree[2]++;
            if (m_shower_pdg_from_start_segment ? skip_proto : skip_legacy) continue;
        }
        auto [sv, stype] = shower->get_start_vertex_and_type();
        if (stype == 1) {
            type1_showers.push_back(shower);
        } else if (stype == 2) {
            type2_dirs[shower.get()] = shower_cal_dir_3vector(*shower, shower->get_start_point(), 100 * units::cm);
            type2_showers.push_back(shower);
        }
    }

    if (type1_showers.empty() || type2_showers.empty()) {
        // doc pr/91 G-B: with no conn-1 shower on the main vertex, or no conn-2
        // one, this pass merges nothing at all whatever the geometry.
        if (pr91_merge_dbg())
            std::fprintf(stderr, "SHOWER_MERGE tag=merge_showers SKIP_PASS reason=%s "
                                 "n_type1=%zu n_type2=%zu\n",
                         type1_showers.empty() ? "no_conn1_at_main_vertex" : "no_conn2_at_main_vertex",
                         type1_showers.size(), type2_showers.size());
        return;
    }

    // Phase 1: collect merges in deterministic order.
    // 'claimed' prevents a type-2 shower from being absorbed by more than one type-1 shower.
    std::unordered_set<ShowerPtr> claimed;
    std::vector<std::pair<ShowerPtr, std::vector<ShowerPtr>>> merge_plan;

    for (auto& shower1 : type1_showers) {
        WireCell::Vector dir1 = shower_cal_dir_3vector(*shower1, shower1->get_start_point(), 100 * units::cm);
        std::vector<ShowerPtr> to_merge;
        for (auto& shower2 : type2_showers) {
            if (claimed.count(shower2)) continue;
            const WireCell::Vector& dir2 = type2_dirs.at(shower2.get());
            double cos_angle = dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude());
            double angle_deg = std::acos(std::max(-1.0, std::min(1.0, cos_angle))) * 180.0 / M_PI;
            const bool pr91_fires = angle_deg < 10.0;
            if (pr91_merge_dbg())
                std::fprintf(stderr, "SHOWER_MERGE tag=merge_showers keep_sid=%d keep_node=%d "
                                     "cand_sid=%d cand_node=%d angle100=%.3f cut=10.0 verdict=%s\n",
                             shower1->get_shower_id(), pr91_seg_display_id(shower1->start_segment()),
                             shower2->get_shower_id(), pr91_seg_display_id(shower2->start_segment()),
                             angle_deg, pr91_fires ? "MERGE" : "angle100>=10");
            if (pr91_fires) {
                to_merge.push_back(shower2);
                claimed.insert(shower2);
            }
        }
        if (!to_merge.empty()) {
            merge_plan.emplace_back(shower1, std::move(to_merge));
        }
    }

    // Phase 2: merge all collected shower2s into each shower1, then compute kinematics once.
    for (auto& [shower1, to_merge] : merge_plan) {
        for (auto& shower2 : to_merge) {
            pr93_probe_absorb_splice("examine_merge_showers", shower1, shower2);
            shower1->add_shower(*shower2);
        }
        shower1->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
        shower1->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
        double kine_charge = cal_kine_charge(shower1, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
        shower1->set_kine_charge(kine_charge);
        shower1->set_flag_kinematics(true);
    }

    for (auto& shower : claimed) {
        showers.erase(shower);
    }
    if (!claimed.empty()) {
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                          map_vertex_to_shower, used_shower_clusters);
    }
}


void PatternAlgorithms::shower_clustering_in_other_clusters(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_save){
    
    if (!main_vertex || !main_cluster) return;

    // Build map_vertex_segments and map_segment_vertices (ordered for determinism)
    std::map<VertexPtr, std::vector<SegmentPtr>> map_vertex_segments;
    std::map<SegmentPtr, std::set<VertexPtr>> map_segment_vertices;
    std::vector<SegmentPtr> seg_order;

    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg) continue;

        auto source_vdesc = boost::source(e, graph);
        auto target_vdesc = boost::target(e, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;


        if (map_segment_vertices.find(seg) == map_segment_vertices.end()) {
            seg_order.push_back(seg);
        }
        if (v1) {
            map_vertex_segments[v1].push_back(seg);
            map_segment_vertices[seg].insert(v1);
        }
        if (v2) {
            map_vertex_segments[v2].push_back(seg);
            map_segment_vertices[seg].insert(v2);
        }
    }
    
    // Build map_cluster_length
    std::map<Facade::Cluster*, double> map_cluster_length;
    for (auto cluster : other_clusters) {
        map_cluster_length[cluster] = cluster->get_length();
    }
    
    // Collect vertices in main cluster as well as existing showers (in deterministic order)
    std::vector<VertexPtr> vertices;
    for (auto v : ordered_nodes(graph)) {
        VertexPtr vtx = graph[v].vertex;
        if (!vtx) continue;
        if ((vtx->cluster() && vtx->cluster()->get_cluster_id() == main_cluster->get_cluster_id()) ||
            map_vertex_in_shower.find(vtx) != map_vertex_in_shower.end()) {
            vertices.push_back(vtx);
        }
    }
    
    // Process clusters in map_cluster_main_vertices
    for (auto& [cluster, vertex] : map_cluster_main_vertices) {
        if (used_shower_clusters.find(cluster) != used_shower_clusters.end()) continue;
        if (m_nv_bridge_cluster_ids.count(cluster->get_cluster_id())) continue;  // doc pr/40 round 9 B2: a bridged cluster is main-track material, not shower material (set empty when bridge off)
        if (map_cluster_length[cluster] < 4 * units::cm) continue;
        
        double min_dis = 1e9;
        VertexPtr min_vertex = nullptr;
        double main_dis = 1e9;
        
        for (auto vtx : vertices) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            
            // Get closest distance using cluster's get_closest_dis method
            double dis = cluster->get_closest_dis(vtx_pt);
            
            if (dis < min_dis) {
                min_dis = dis;
                min_vertex = vtx;
            }
            
            if (vtx == main_vertex) {
                main_dis = dis;
            }
        }
        
        if (min_dis > 0.8 * main_dis) {
            min_dis = main_dis;
            min_vertex = main_vertex;
        }
        
        // Find a shower segment starting at vertex.
        // Matches prototype's get_flag_shower(): kShowerTrajectory || kShowerTopology || abs(pdg)==11.
        SegmentPtr sg = nullptr;
        if (map_vertex_segments.find(vertex) != map_vertex_segments.end()) {
            for (auto seg : map_vertex_segments[vertex]) {
                bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                 seg->flags_any(SegmentFlags::kShowerTopology) ||
                                 (seg->has_particle_info() && seg->particle_info() &&
                                  std::abs(seg->particle_info()->pdg()) == 11);
                if (is_shower) {
                    sg = seg;
                    break;
                }
            }
        }
        
        int connection_type = 3;
        
        if (sg) {
            // Ensure the segment has a "fit" dpcloud before adding to the shower.
            // Segments from other clusters may have fits() but no dpcloud("fit") if
            // they were not individually track-fitted (create_segment_fit_point_cloud
            // was never called for them). Build it now so set_start_segment can copy it.
            if (!sg->dpcloud("fit") && !sg->fits().empty()) {
                create_segment_fit_point_cloud(sg, dv, "fit");
            }


            // Create new shower
            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(min_vertex, connection_type);
            shower->set_start_segment(sg);

            // Record the cluster's own vertex as the shower start point.
            // This is used by the merge-check loop below (and by subsequent
            // iterations of this loop) via get_start_point(). Without this call
            // get_start_point() returns {0,0,0} and the merge angles are wrong.
            WireCell::Point vertex_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
            shower->set_start_point(vertex_pt);
            const auto& fits = sg->fits();
            if (!fits.empty()) {
                double dis1 = (vertex_pt - fits.front().point).magnitude();
                double dis2 = (vertex_pt - fits.back().point).magnitude();
                
                if (dis1 < dis2) {
                    sg->dirsign(1);
                } else {
                    sg->dirsign(-1);
                }
            }
            
            // Complete shower structure
            IndexedSegmentSet used_segments;
            pr93_probe_absorb_site("in_other_clusters_A", shower, m_shower_absorb_track_guard);
            shower->complete_structure_with_start_segment(nv_bridge_seed(used_segments), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
            pr84_probe_shower(shower, "in_other_clusters_A");

            // Calculate shower direction
            WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, vertex_pt, 15 * units::cm);
            
            // Cluster with the rest - add segments based on angle and distance
            for (auto seg1 : seg_order) {
                if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
                if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
                if (seg1->cluster() == shower->start_segment()->cluster()) continue;
                
                // Find the closest point
                auto [pair_dis, pair_point] = segment_get_closest_point(seg1, vertex_pt);
                WireCell::Vector v1(pair_point.x() - vertex_pt.x(), 
                                   pair_point.y() - vertex_pt.y(), 
                                   pair_point.z() - vertex_pt.z());
                
                double angle = std::acos(std::clamp(dir_shower.dot(v1) / (dir_shower.magnitude() * v1.magnitude()), -1.0, 1.0));
                angle = angle / M_PI * 180.0;
                
                const bool pr91_fires = (angle < 25 && pair_dis < 80 * units::cm) ||
                                        (angle < 12.5 && pair_dis < 120 * units::cm);
                if (pr91_merge_dbg())
                    std::fprintf(stderr, "SHOWER_MERGE tag=other_cl_seg new_sid=%d new_node=%d "
                                         "cand_seg=%d cand_cluster=%d angle=%.3f dis=%.3f verdict=%s\n",
                                 shower->get_shower_id(), pr91_seg_display_id(shower->start_segment()),
                                 pr91_seg_display_id(seg1),
                                 seg1->cluster() ? seg1->cluster()->get_cluster_id() : -1,
                                 angle, pair_dis / units::cm,
                                 pr91_fires ? "ABSORB" : "outside_cone");
                if (pr91_fires) {
                    pr93_probe_absorb_direct("in_other_clusters_seg_cone", shower, seg1);
                    shower->add_segment(seg1, true);
                }
            }

            // Force PDG=0 or short+weak-muon start segment to electron before majority-vote.
            // Mirrors the pattern in sub-pass 2 (and shower_clustering_with_nv_from_vertices).
            {
                int particle_type_sp1 = 0;
                if (sg->has_particle_info() && sg->particle_info())
                    particle_type_sp1 = sg->particle_info()->pdg();
                if (particle_type_sp1 == 0 ||
                    (std::abs(particle_type_sp1) == 13 &&
                     segment_track_length(sg) < 40 * units::cm && seg_dir_weak(sg))) {
                    auto four_momentum = segment_cal_4mom(sg, 11, particle_data, recomb_model, m_mip_dqdx);
                    sg->particle_info(std::make_shared<Aux::ParticleInfo>(
                        11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11),
                        four_momentum));
                }
            }

            // Update particle type
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);

            // Check with other showers and merge if needed
            std::vector<ShowerPtr> showers_to_be_removed;
            for (auto shower1 : showers) {
                WireCell::Point start_pt1 = shower1->get_start_point();
                WireCell::Vector dir_shower1 = shower_cal_dir_3vector(*shower1, start_pt1, 15 * units::cm);
                WireCell::Vector dir2(start_pt1.x() - vertex_pt.x(),
                                     start_pt1.y() - vertex_pt.y(),
                                     start_pt1.z() - vertex_pt.z());
                
                double angle = std::acos(std::clamp(dir_shower.dot(dir_shower1) / (dir_shower.magnitude() * dir_shower1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                double angle1 = std::acos(std::clamp(dir_shower.dot(dir2) / (dir_shower.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                
                // doc pr/91 G-D: `dir_shower` is the NEW (downstream) shower's
                // axis and `dir2` runs from its vertex TO the candidate, so an
                // existing shower sitting UPSTREAM of this one fails angle1 by
                // construction -- the absorb only ever points downstream.
                const bool pr91_fires =
                    (angle < 25 && angle1 < 15 && dir2.magnitude() < 80 * units::cm) ||
                    (angle < 12.5 && angle1 < 7.5 && dir2.magnitude() < 120 * units::cm);
                if (pr91_merge_dbg())
                    std::fprintf(stderr, "SHOWER_MERGE tag=other_cl_absorb new_sid=%d new_node=%d "
                                         "cand_sid=%d cand_node=%d cand_conn=%d angle=%.3f angle1=%.3f "
                                         "dis=%.3f verdict=%s\n",
                                 shower->get_shower_id(), pr91_seg_display_id(shower->start_segment()),
                                 shower1->get_shower_id(), pr91_seg_display_id(shower1->start_segment()),
                                 shower1->get_start_vertex_and_type().second,
                                 angle, angle1, dir2.magnitude() / units::cm,
                                 pr91_fires ? "ABSORB" : "outside_cone");
                if (pr91_fires) {
                    pr93_probe_absorb_splice("in_other_clusters_cone", shower, shower1);
                    shower->add_shower(*shower1);
                    showers_to_be_removed.push_back(shower1);
                }
            }
            
            for (auto shower_to_remove : showers_to_be_removed) {
                showers.erase(shower_to_remove);
            }

            // Post-merge majority-vote and kinematics (prototype lines 1555-1556)
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
            shower->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);

            showers.insert(shower);
        }
    }
    
    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                      map_vertex_to_shower, used_shower_clusters);
    
    // Process remaining other_clusters not in map_cluster_main_vertices
    for (auto cluster : other_clusters) {
        if (used_shower_clusters.find(cluster) != used_shower_clusters.end()) continue;
        if (m_nv_bridge_cluster_ids.count(cluster->get_cluster_id())) continue;  // doc pr/40 round 9 B2: a bridged cluster is main-track material, not shower material (set empty when bridge off)
        
        // std::cout << "shower_clustering_in_other_clusters: Processing cluster " << cluster->get_cluster_id() << std::endl;

        double min_dis = 1e9;
        VertexPtr min_vertex = nullptr;
        double main_dis = 1e9;
        
        for (auto vtx : vertices) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            
            // Get closest distance using cluster's get_closest_dis method
            double dis = cluster->get_closest_dis(vtx_pt);
            
            if (dis < min_dis) {
                min_dis = dis;
                min_vertex = vtx;
            }
            
            if (vtx == main_vertex) {
                main_dis = dis;
            }
        }
        
        if (min_dis > 0.8 * main_dis) {
            min_dis = main_dis;
            min_vertex = main_vertex;
        }
        
        // Find a segment from this cluster in deterministic graph-index order
        SegmentPtr sg = nullptr;
        for (auto seg : seg_order) {
            if (seg->cluster() != cluster) continue;
            sg = seg;
            break;
        }
        
        int connection_type = 3;
        if (min_dis > 80 * units::cm) {
            connection_type = 4;
        }

        if (!flag_save) connection_type = 4;
        
        if (sg) {
            // Ensure the segment has a "fit" dpcloud before adding to the shower.
            if (!sg->dpcloud("fit") && !sg->fits().empty()) {
                create_segment_fit_point_cloud(sg, dv, "fit");
            }

            // std::cout << "shower_clustering_in_other_clusters: Found shower segment " << sg->id() << " for cluster " << cluster->get_cluster_id() << " " << sg->fits().size() << std::endl;


            // Create new shower
            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(min_vertex, connection_type);
            shower->set_start_segment(sg);

            // Set direction if not already set
            if (sg->dirsign() == 0) {
                // find_vertices returns (front_vtx, back_vtx) ordered by proximity
                // to sg->wcpts().front().point — correct geometric order, independent
                // of the arbitrary boost::source/target ordering on an undirected graph.
                auto [front_vtx, back_vtx] = find_vertices(graph, sg);

                if (front_vtx && back_vtx) {
                    if (map_vertex_segments[front_vtx].size() == 1 && map_vertex_segments[back_vtx].size() > 1) {
                        sg->dirsign(1);
                    } else if (map_vertex_segments[front_vtx].size() > 1 && map_vertex_segments[back_vtx].size() == 1) {
                        sg->dirsign(-1);
                    } else {
                        // Examine vertices based on distance to main_vertex
                        WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
                        const auto& fits = sg->fits();
                        if (!fits.empty()) {
                            double dis1 = (main_vtx_pt - fits.front().point).magnitude();
                            double dis2 = (main_vtx_pt - fits.back().point).magnitude();

                            if (dis1 < dis2) {
                                sg->dirsign(1);
                            } else {
                                sg->dirsign(-1);
                            }
                        }
                    }
                }
            }
            
            // Set particle type to electron if needed
            int particle_type = 0;
            if (sg->has_particle_info() && sg->particle_info()) {
                particle_type = sg->particle_info()->pdg();
            }
            
            if (particle_type == 0 || 
                (std::abs(particle_type) == 13 && segment_track_length(sg) < 40 * units::cm && seg_dir_weak(sg))) {
                auto four_momentum = segment_cal_4mom(sg, 11, particle_data, recomb_model, m_mip_dqdx);
                
                // Create ParticleInfo for electron
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    11,                                          // electron PDG
                    particle_data->get_particle_mass(11),       // electron mass
                    particle_data->pdg_to_name(11),             // "electron"
                    four_momentum                                // 4-momentum
                );
                
                sg->particle_info(pinfo);
            }
            
            // Complete shower structure
            IndexedSegmentSet used_segments;
            pr93_probe_absorb_site("in_other_clusters_B", shower, m_shower_absorb_track_guard);
            shower->complete_structure_with_start_segment(nv_bridge_seed(used_segments), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
            // Majority-vote correction for multi-segment showers whose start segment
            // has an unexpected PDG not covered by the explicit force-to-11 above.
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
            pr84_probe_shower(shower, "in_other_clusters_B");
            showers.insert(shower);
        }
    }

    // doc sbnd_xin/docs/pr/74 round 2 K5 = doc pr/65's deferred rung 2 --
    // docstring at m_shower_conn3_unreachable (NeutrinoPatternBase.h).
    // Extend the leftover-cluster branch just above to graph-unreachable,
    // unclaimed main-cluster segments (18255-142421 seg 7013): same
    // nearest-candidate-vertex anchor, same <80 cm conn-3/conn-4 split, same
    // component claim via complete_structure_with_start_segment.  false =
    // no pass = byte-identical.
    if (m_shower_conn3_unreachable && flag_save && main_vertex->descriptor_valid()) {
        const auto unreachable = unreachable_segments(graph, main_vertex);
        // doc pr/74 round 3 Q2 (18255-142421 seg 7013): the anchor search below
        // must only consider vertices the main vertex can actually reach.  The
        // promoted segment's OWN endpoints are in `vertices` too, and they sit
        // at distance 0 from it, so the unrestricted search always picked one
        // of them: the pseudo-gamma collapsed to zero length at the component's
        // far end and, because conn-3 derives start_point from the anchor
        // (PRShower.cxx:1140) and end_point from the farthest vertex, the
        // reconstructed e- ran BACKWARDS -- from the far end toward the
        // neutrino vertex instead of outward from it.
        const auto reachable_vtxs = reachable_vertices(graph, main_vertex);
        IndexedSegmentSet claimed_k5;
        for (auto seg : seg_order) {
            if (!seg || seg->cluster() != main_cluster) continue;
            if (!unreachable.count(seg)) continue;
            if (map_segment_in_shower.count(seg) || claimed_k5.count(seg)) continue;
            if (segment_track_length(seg) < m_conn3_unreachable_min_len) continue;

            // Nearest candidate vertex -- same preference rule as above.
            double min_dis = 1e9;
            VertexPtr min_vertex = nullptr;
            double main_dis = 1e9;
            for (auto vtx : vertices) {
                if (!reachable_vtxs.count(vtx)) continue;   // round 3 Q2: never anchor to your own component
                WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                double dis = segment_get_closest_point(seg, vtx_pt).first;
                if (dis < min_dis) {
                    min_dis = dis;
                    min_vertex = vtx;
                }
                if (vtx == main_vertex) main_dis = dis;
            }
            if (!min_vertex) continue;
            if (min_dis > 0.8 * main_dis) {
                min_dis = main_dis;
                min_vertex = main_vertex;
            }
            int connection_type = (min_dis > 80 * units::cm) ? 4 : 3;

            if (!seg->dpcloud("fit") && !seg->fits().empty()) {
                create_segment_fit_point_cloud(seg, dv, "fit");
            }

            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(min_vertex, connection_type);
            shower->set_start_segment(seg);

            // Force-to-electron rule of the leftover-cluster branch, plus the
            // low-confidence proton case this class exhibits (7013: pdg 2212
            // at particle_score 0.27 on a 2.32x-MIP EM chunk).  A confident
            // PID keeps its label but still becomes PF-visible.
            int particle_type = 0;
            if (seg->has_particle_info() && seg->particle_info()) {
                particle_type = seg->particle_info()->pdg();
            }
            if (particle_type == 0 ||
                (std::abs(particle_type) == 13 && segment_track_length(seg) < 40 * units::cm && seg_dir_weak(seg)) ||
                (particle_type == 2212 && seg->particle_score() < 0.3)) {
                // doc sbnd_xin/docs/pr/93 Cause B, third seat (SBND
                // 18255-315167): the pdg==2212 arm above was tuned for a
                // SHORT low-confidence EM chunk (142421: 0.27-score on a
                // 2.32x-MIP chunk); after the Cause-D guard frees 315167's
                // 150.7cm score-0.10 proton, this arm re-stamped it e-.
                // Same knob + floor as the other two acceptance seats.
                // C++ default false => byte-identical.
                if (m_shower_accept_pid_guard &&
                    segment_track_length(seg) > m_shower_pid_guard_min_len &&
                    segment_confident_nonelectron_pid(seg)) {
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr93 accept_pid_guard: decline e- write seg={} pdg={} score={:.3f} (conn3_unreachable_promote)",
                        pr91_seg_display_id(seg), particle_type, seg->particle_score());
                }
                else {
                    auto four_momentum = segment_cal_4mom(seg, 11, particle_data, recomb_model, m_mip_dqdx);
                    seg->particle_info(std::make_shared<Aux::ParticleInfo>(
                        11, particle_data->get_particle_mass(11), particle_data->pdg_to_name(11),
                        four_momentum));
                }
            }

            pr93_probe_absorb_site("conn3_unreachable", shower, m_shower_absorb_track_guard);
            shower->complete_structure_with_start_segment(nv_bridge_seed(claimed_k5), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr74 conn3_unreachable: promote gidx={} len {:.1f}cm conn={} anchor_dis {:.1f}cm",
                seg->get_graph_index(), segment_track_length(seg)/units::cm,
                connection_type, min_dis/units::cm);
            pr84_probe_shower(shower, "conn3_unreachable");
            showers.insert(shower);
        }
    }

    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                      map_vertex_to_shower, used_shower_clusters);
}


void PatternAlgorithms::examine_shower_1(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex) return;

    // Build map_vertex_segments (ordered for determinism)
    std::map<VertexPtr, std::vector<SegmentPtr>> map_vertex_segments;
    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg) continue;

        auto source_vdesc = boost::source(e, graph);
        auto target_vdesc = boost::target(e, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;

        if (v1) map_vertex_segments[v1].push_back(seg);
        if (v2) map_vertex_segments[v2].push_back(seg);
    }
    
    // Check if there is already a large EM shower connecting to main_vertex
    bool flag_skip = false;
    auto it = map_vertex_to_shower.find(main_vertex);
    if (it != map_vertex_to_shower.end()) {
        for (auto shower : it->second) {
            if (shower->start_segment() && shower->start_segment()->has_particle_info() && 
                shower->start_segment()->particle_info()->pdg() != 11) continue;
            
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            if (conn_type != 1) continue;
            
            double energy = 0;
            if (shower->get_kine_best() != 0) {
                energy = shower->get_kine_best();
            } else {
                energy = shower->get_kine_charge();
            }
            if (energy > 80 * units::MeV) flag_skip = true;
        }
    }
    
    bool flag_added = false;

    // doc pr/91 G-C: flag_skip only suppresses the first half; the conn-3 merge
    // in the second half still runs, so record which half is live.
    if (pr91_merge_dbg())
        std::fprintf(stderr, "SHOWER_MERGE tag=ex_shower1_merge ENTER flag_skip=%d nshowers=%zu\n",
                     flag_skip ? 1 : 0, showers.size());

    if (!flag_skip) {
        IndexedShowerSet used_showers;
        std::map<SegmentPtr, std::set<ShowerPtr>, SegmentIndexCmp> map_segment_showers;
        ShowerSegmentMap map_segment_new_shower;
        std::vector<SegmentPtr> seg_order;  // preserves outer-loop order for deterministic processing
        IndexedSegmentSet used_segments;
        IndexedShowerSet del_showers;

        WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;

        // Loop over segments at main_vertex
        auto it_mv = map_vertex_segments.find(main_vertex);
        if (it_mv != map_vertex_segments.end()) {
            for (auto sg : it_mv->second) {
                // Skip strong direction segments or certain particle types
                if (!seg_dir_weak(sg)) continue;
                
                int particle_type = 0;
                if (sg->has_particle_info() && sg->particle_info()) {
                    particle_type = sg->particle_info()->pdg();
                }
                
                double medium_dQ_dx = segment_median_dQ_dx(sg);
                if (particle_type == 2212 && medium_dQ_dx / m_mip_dqdx_median > 1.6) continue;
                if (particle_type == 11) continue;
                
                // Form a new shower
                ShowerPtr shower1 = std::make_shared<Shower>(graph);
                shower1->set_start_vertex(main_vertex, 1);
                shower1->set_start_segment(sg);
                pr93_probe_absorb_site("examine_shower_1_tmp", shower1, m_shower_absorb_track_guard);
                shower1->complete_structure_with_start_segment(nv_bridge_seed(used_segments), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
                pr84_probe_shower(shower1, "examine_shower_1_tmp");


                WireCell::Vector dir1 = segment_cal_dir_3vector(sg, main_vtx_pt, 15 * units::cm);

                // Build shower1_vertices once — shower1 does not change in the inner loop.
                TrajectoryView& traj1 = shower1->fill_maps();
                std::set<VertexPtr> shower1_vertices;
                for (auto vdesc : traj1.nodes()) {
                    auto vtx = traj1.view_graph()[vdesc].vertex;
                    if (vtx) shower1_vertices.insert(vtx);
                }

                // Check against existing showers.
                // Order of guards: cheapest first to avoid the two expensive KD-tree calls
                // (shower_get_closest_dis, shower_get_closest_point) for the majority of showers.
                for (auto shower : showers) {
                    // 1. Already claimed by an earlier segment — skip before any other work.
                    if (used_showers.count(shower)) continue;

                    // 2. Cheap scalar/pointer checks.
                    if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                        shower->start_segment()->particle_info()->pdg() != 11) continue;

                    auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                    if (start_vtx == main_vertex && conn_type == 1) continue;

                    // 3. First expensive KD-tree call — only for showers that pass cheap filters.
                    double min_dis = shower_get_closest_dis(*shower, sg);
                    if (conn_type > 2 && min_dis > 3 * units::cm) continue;

                    double energy = shower->get_kine_charge();

                    if (conn_type == 1 && energy > 80 * units::MeV) {
                        if (shower1_vertices.find(start_vtx) == shower1_vertices.end()) continue;
                        else {
                            WireCell::Vector dir3 = shower_cal_dir_3vector(*shower, shower->get_start_point(), 15 * units::cm);
                            double angle = std::acos(std::clamp(dir3.dot(dir1) / (dir3.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                            if (angle > 30) continue;
                        }
                    } else {
                        // Skip a conn_type=1 shower whose start_vertex is already inside shower1.
                        // Such a shower is a downstream decay product (e.g. a Michel electron that
                        // starts at the far-end junction of sg) — its presence cannot be used as
                        // evidence that sg itself is an EM shower.  High-energy (>80 MeV) downstream
                        // showers are handled by the if-branch above and are unaffected.
                        if (conn_type == 1 && shower1_vertices.count(start_vtx)) continue;

                        WireCell::Vector dir3 = shower_cal_dir_3vector(*shower, shower->get_start_point(), 15 * units::cm);

                        // Find closest vertex in shower1 to shower start point
                        WireCell::Point min_point;
                        double min_vtx_dis = 1e9;
                        for (auto vtx3 : shower1_vertices) {
                            WireCell::Point vtx3_pt = vtx3->fit().valid() ? vtx3->fit().point : vtx3->wcpt().point;
                            double dis = (vtx3_pt - shower->get_start_point()).magnitude();
                            if (dis < min_vtx_dis) {
                                min_vtx_dis = dis;
                                min_point = vtx3_pt;
                            }
                        }

                        WireCell::Vector dir4(shower->get_start_point().x() - min_point.x(),
                                             shower->get_start_point().y() - min_point.y(),
                                             shower->get_start_point().z() - min_point.z());

                        double angle3 = std::acos(std::clamp(dir3.dot(dir1) / (dir3.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double angle4 = std::acos(std::clamp(dir4.dot(dir1) / (dir4.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double tmp_angle = std::min(angle3, angle4);

                        if (energy > 25 * units::MeV && tmp_angle > 40) continue;
                    }

                    // 4. Second expensive KD-tree call — final angle test.
                    auto [shower_dis, shower_point] = shower_get_closest_point(*shower, main_vtx_pt);
                    WireCell::Vector dir2(shower_point.x() - main_vtx_pt.x(),
                                         shower_point.y() - main_vtx_pt.y(),
                                         shower_point.z() - main_vtx_pt.z());

                    double angle = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;

                    if ((angle < 15 && min_dis < 36 * units::cm) ||
                        (angle < 10 && min_dis < 46 * units::cm) ||
                        (angle < 7.5)) {
                        map_segment_showers[sg].insert(shower);
                        used_showers.insert(shower);
                    }
                }

                if (map_segment_showers.count(sg)) {
                    map_segment_new_shower[sg] = shower1;
                    seg_order.push_back(sg);  // record in outer-loop (graph-index) order
                }
            }
        }
        
        // Process segments in the same order they were discovered (graph-index order),
        // not in pointer order as a plain map iteration would give.
        for (auto sg : seg_order) {
            auto& associated_showers = map_segment_showers[sg];
            ShowerPtr shower1 = map_segment_new_shower[sg];
            int num_showers = associated_showers.size();
            
            double max_energy = 0;
            double total_energy = 0;
            for (auto shower : associated_showers) {
                double energy = shower->get_kine_charge();
                if (energy > max_energy) max_energy = energy;
                total_energy += energy;
            }
            
            // Analyze shower1 structure
            TrajectoryView& traj = shower1->fill_maps();
            std::set<VertexPtr> shower1_vertices;
            for (auto vdesc : traj.nodes()) {
                auto vtx = traj.view_graph()[vdesc].vertex;
                if (vtx) shower1_vertices.insert(vtx);
            }
            
            double max_length = 0;
            SegmentPtr max_sg = nullptr;
            int n_tracks = 0;
            int n_showers = 0;
            double total_length = 0;
            bool flag_good_track = false;
            
            // ordered_edges: total_length is an FP accumulation and it feeds the
            // hard cuts below (`total_length < 70cm`, `< n_tracks * 36cm`, ...),
            // so the walk order can move a cut (doc pr/28 §15.2).
            for (auto edesc : ordered_edges(traj, graph)) {
                auto sg1 = traj.view_graph()[edesc].segment;
                if (!sg1) continue;

                double length = segment_track_length(sg1);
                double medium_dQ_dx = segment_median_dQ_dx(sg1) / m_mip_dqdx_median;
                
                if (!seg_dir_weak(sg1)) {
                    // Find end vertex
                    auto seg_edesc = sg1->get_descriptor();
                    auto source_vdesc = boost::source(seg_edesc, graph);
                    auto target_vdesc = boost::target(seg_edesc, graph);
                    VertexPtr v1 = graph[source_vdesc].vertex;
                    VertexPtr v2 = graph[target_vdesc].vertex;
                    
                    VertexPtr end_vertex = nullptr;
                    if (sg1->dirsign() == 1) {
                        if (!sg1->fits().empty() && v1 && v2) {
                            auto& fits = sg1->fits();
                            double dist1 = (fits.back().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.back().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    } else if (sg1->dirsign() == -1) {
                        if (!sg1->fits().empty() && v1 && v2) {
                            auto& fits = sg1->fits();
                            double dist1 = (fits.front().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.front().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    }
                    
                    if (end_vertex && map_vertex_segments.find(end_vertex) != map_vertex_segments.end()) {
                        if (map_vertex_segments[end_vertex].size() > 1) {
                            bool flag_non_ele = false;
                            for (auto sg2 : map_vertex_segments[end_vertex]) {
                                if (sg2 == sg1) continue;
                                // Matches prototype get_flag_shower(): shower topology/trajectory OR pdg==11
                                bool sg2_is_shower = sg2->flags_any(SegmentFlags::kShowerTrajectory) ||
                                                     sg2->flags_any(SegmentFlags::kShowerTopology) ||
                                                     (sg2->has_particle_info() && sg2->particle_info() &&
                                                      std::abs(sg2->particle_info()->pdg()) == 11);
                                if (!sg2_is_shower) flag_non_ele = true;
                            }
                            if (!flag_non_ele && map_vertex_segments[end_vertex].size() <= 3) {
                                flag_good_track = true;
                            }
                        } else {
                            flag_good_track = true;
                        }
                    }
                }
                
                bool is_shower = sg1->flags_any(SegmentFlags::kShowerTrajectory) || 
                                sg1->flags_any(SegmentFlags::kShowerTopology);
                if (is_shower) n_showers++;
                
                n_tracks++;
                total_length += length;
                // Tie-broken on id(), matching the same max-segment idiom at
                // :188 and :327 in this file.  This one was the odd site out:
                // with a bare `<`, two equal-length segments were separated by
                // the walk order alone (doc pr/28 §15.8).
                if (max_length < length ||
                    (max_length == length && max_sg && sg1->id() < max_sg->id())) {
                    max_length = length;
                    max_sg = sg1;
                }
                (void)medium_dQ_dx; // To avoid unused variable warning
            }
            (void)flag_good_track;
            (void)n_showers;           

            // Check if should skip
            bool flag_skip_segment = false;
            for (auto shower_check : showers) {
                auto [start_vtx_check, conn_type_check] = shower_check->get_start_vertex_and_type();
                if (conn_type_check == 1 && associated_showers.find(shower_check) == associated_showers.end()) {
                    if (shower1_vertices.find(start_vtx_check) != shower1_vertices.end() && 
                        start_vtx_check != main_vertex && 
                        shower_check->get_kine_charge() > 60 * units::MeV) {
                        flag_skip_segment = true;
                    }
                }
            }
            
            // Decide whether to create this shower
            if (total_length < 70 * units::cm && 
                ((n_tracks == 1 && total_length < 60 * units::cm) ||
                 (n_tracks == 1 && total_length < 65 * units::cm && num_showers > 3 && total_energy > 150 * units::MeV) ||
                 total_length < n_tracks * 36 * units::cm) &&
                (total_energy > 50 * units::MeV || total_energy / units::MeV > total_length / units::cm * 0.75) &&
                !flag_skip_segment) {
                
                // Set particle type to electron
                if (shower1->start_segment() && shower1->start_segment()->has_particle_info() &&
                    shower1->start_segment()->particle_info()) {
                    // doc sbnd_xin/docs/pr/93 Cause B (SBND 18255-348471: a
                    // pdg=2212 score=0.23 proton force-relabelled e- here).
                    // Decline the pdg write ONLY -- kAvoidMuonCheck and the
                    // update_particle_type call below are untouched.  The
                    // m_shower_pid_guard_min_len floor (50cm) is required:
                    // real nueCC electron stems of 22-47cm carry equally
                    // confident (0.11-0.64) proton/muon template scores
                    // (the template competition never considers electron),
                    // and 9/48 nueCC48 events regressed without it.
                    // C++ default false => byte-identical.
                    if (m_shower_accept_pid_guard &&
                        segment_track_length(shower1->start_segment()) > m_shower_pid_guard_min_len &&
                        segment_confident_nonelectron_pid(shower1->start_segment())) {
                        SPDLOG_LOGGER_DEBUG(s_log,
                            "pr93 accept_pid_guard: decline e- write seg={} pdg={} score={:.3f} (new_shower_accepted)",
                            shower1->start_segment()->id(),
                            shower1->start_segment()->particle_info()->pdg(),
                            shower1->start_segment()->particle_score());
                    }
                    else {
                        shower1->start_segment()->particle_info()->set_pdg(11);
                        pr40_probe_setpdg(shower1->start_segment(), 11, "NeutrinoShowerClustering.cxx:new_shower_accepted");
                        shower1->start_segment()->particle_info()->set_mass(
                            particle_data->get_particle_mass(11));
                    }
                }
                shower1->start_segment()->set_flags(SegmentFlags::kAvoidMuonCheck);
                shower1->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);

                // Merge associated showers
                for (auto shower : associated_showers) {
                    del_showers.insert(shower);
                    pr93_probe_absorb_splice("examine_shower_1_assoc", shower1, shower);
                    shower1->add_shower(*shower);
                }

                shower1->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
                double kine_charge = cal_kine_charge(shower1, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
                shower1->set_kine_charge(kine_charge);
                shower1->set_flag_kinematics(true);

                // Dedup: remove any pre-existing shower at (main_vertex, conn_type=1) whose
                // start_segment is the same as this new shower1's start_segment.  Such a shower
                // is a stale single-segment wrapper left from shower_clustering_with_nv_in_main_cluster;
                // shower1 supersedes it.  Without this both objects survive in the showers set and
                // fill_bee_pf_tree emits duplicate nodes with identical cluster*1000+seg_id display ids.
                for (auto sh_old : showers) {
                    if (sh_old == shower1) continue;
                    if (sh_old->start_segment() != sg) continue;
                    auto [svtx_old, ctype_old] = sh_old->get_start_vertex_and_type();
                    if (svtx_old == main_vertex && ctype_old == 1) {
                        del_showers.insert(sh_old);
                    }
                }

                showers.insert(shower1);
                SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_in_other_clusters: Create a new low-energy shower: {} MeV", kine_charge / units::MeV);
                flag_added = true;
            }
        }
        
        // Remove deleted showers
        for (auto shower : del_showers) {
            showers.erase(shower);
        }
        
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                          map_vertex_to_shower, used_shower_clusters);
    }
    
    // Second part: merge existing showers if nothing was added
    // doc pr/91 G-C: this is the ONLY pass that can pull a conn-3 shower into
    // the main conn-1 shower, and it does not run at all when the first half
    // created a trial shower.
    if (pr91_merge_dbg() && flag_added)
        std::fprintf(stderr, "SHOWER_MERGE tag=ex_shower1_merge SKIP_PASS reason=flag_added\n");
    if (!flag_added) {
        // Index-ordered (not pointer-ordered): acc_energy below FP-sums over
        // the associated set, so iteration order must be content-stable.
        std::map<ShowerPtr, std::set<ShowerPtr, ShowerIndexCmp>, ShowerIndexCmp> map_shower_showers;
        WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
        
        auto it = map_vertex_to_shower.find(main_vertex);
        if (it != map_vertex_to_shower.end()) {
            for (auto shower : it->second) {
                if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                    shower->start_segment()->particle_info()->pdg() != 11) continue;
                
                auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                if (conn_type != 1) continue;
                
                WireCell::Vector dir1 = shower_cal_dir_3vector(*shower, main_vtx_pt, 15 * units::cm);

                // Build shower_vertices once — shower does not change in the inner loop.
                TrajectoryView& traj_shower = shower->fill_maps();
                std::set<VertexPtr> shower_vertices;
                for (auto vdesc : traj_shower.nodes()) {
                    auto vtx = traj_shower.view_graph()[vdesc].vertex;
                    if (vtx) shower_vertices.insert(vtx);
                }

                for (auto shower1 : showers) {
                    // doc pr/91 G-C -- one lambda so every reject prints the same
                    // fields; behaviour is untouched (each call is still a plain
                    // `continue` on the same condition, logging only).
                    auto pr91_rej = [&](const char* why, double md) {
                        if (!pr91_merge_dbg()) return;
                        std::fprintf(stderr, "SHOWER_MERGE tag=ex_shower1_merge parent_sid=%d "
                                             "parent_node=%d cand_sid=%d cand_node=%d cand_conn=%d "
                                             "cand_len=%.3f cand_ke=%.3f min_dis=%.3f verdict=%s\n",
                                     shower->get_shower_id(),
                                     pr91_seg_display_id(shower->start_segment()),
                                     shower1->get_shower_id(),
                                     pr91_seg_display_id(shower1->start_segment()),
                                     shower1->get_start_vertex_and_type().second,
                                     shower1->get_total_length() / units::cm,
                                     shower1->get_kine_charge() / units::MeV,
                                     md, why);
                    };
                    // Cheap guards before the expensive KD-tree call.
                    if (shower1->start_segment() && shower1->start_segment()->has_particle_info() &&
                        shower1->start_segment()->particle_info()->pdg() != 11) { pr91_rej("pdg_not_11", -1); continue; }

                    auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
                    if (start_vtx1 == main_vertex && conn_type1 == 1) { pr91_rej("is_the_parent_class", -1); continue; }
                    if (conn_type1 == 1 && shower_vertices.find(start_vtx1) == shower_vertices.end()) { pr91_rej("conn1_vtx_outside_parent", -1); continue; }
                    if (shower1->get_total_length() < 3 * units::cm) { pr91_rej("len_lt_3cm", -1); continue; }

                    // Expensive KD-tree call — only after cheap filters pass.
                    double min_dis = shower_get_closest_dis(*shower1, shower->start_segment());
                    // G-C: for every conn-3/4 shower this 3 cm gate to the PARENT'S
                    // START SEGMENT (not the parent's whole charge) is the only door.
                    if (conn_type1 > 2 && min_dis > 3 * units::cm) { pr91_rej("conn_gt2_min_dis_gt_3cm", min_dis / units::cm); continue; }

                    double energy = shower1->get_kine_charge();
                    
                    auto [shower_dis, shower_point] = shower_get_closest_point(*shower1, main_vtx_pt);
                    WireCell::Vector dir2(shower_point.x() - main_vtx_pt.x(),
                                         shower_point.y() - main_vtx_pt.y(),
                                         shower_point.z() - main_vtx_pt.z());
                    
                    double angle = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    WireCell::Vector dir3 = shower_cal_dir_3vector(*shower1, shower1->get_start_point(), 30 * units::cm);
                    double angle1 = std::acos(std::clamp(dir2.dot(dir3) / (dir2.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    const bool pr91_fires = angle < 15 && angle1 < 15 && min_dis < 28 * units::cm;
                    if (pr91_merge_dbg())
                        std::fprintf(stderr, "SHOWER_MERGE tag=ex_shower1_merge parent_sid=%d "
                                             "parent_node=%d cand_sid=%d cand_node=%d cand_conn=%d "
                                             "cand_len=%.3f cand_ke=%.3f min_dis=%.3f angle=%.3f "
                                             "angle1=%.3f verdict=%s\n",
                                     shower->get_shower_id(),
                                     pr91_seg_display_id(shower->start_segment()),
                                     shower1->get_shower_id(),
                                     pr91_seg_display_id(shower1->start_segment()),
                                     conn_type1, shower1->get_total_length() / units::cm,
                                     energy / units::MeV, min_dis / units::cm, angle, angle1,
                                     pr91_fires ? "CANDIDATE" : "geometry_fail");
                    if (pr91_fires) {
                        map_shower_showers[shower].insert(shower1);
                    }
                    (void)energy; // To avoid unused variable warning
                }
            }
        }
        
        // Find the shower combination with maximum energy
        std::set<ShowerPtr> del_showers;
        ShowerPtr max_shower = nullptr;
        double max_energy = 0;
        
        for (auto& [shower, associated_showers] : map_shower_showers) {
            double acc_energy = 0;
            if (shower->get_kine_best() != 0) {
                acc_energy += shower->get_kine_best();
            } else {
                acc_energy += shower->get_kine_charge();
            }
            
            for (auto shower1 : associated_showers) {
                if (shower1->get_kine_best() != 0) {
                    acc_energy += shower1->get_kine_best();
                } else {
                    acc_energy += shower1->get_kine_charge();
                }
            }
            
            if (acc_energy > max_energy ||
                (acc_energy == max_energy && max_shower &&
                 shower->start_segment() && max_shower->start_segment() &&
                 shower->start_segment()->id() < max_shower->start_segment()->id())) {
                max_energy = acc_energy;
                max_shower = shower;
            }
        }

        if (max_shower) {
            for (auto shower1 : map_shower_showers[max_shower]) {
                pr93_probe_absorb_splice("examine_showers_max", max_shower, shower1);
                max_shower->add_shower(*shower1);
                del_showers.insert(shower1);
            }
            
            max_shower->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
            max_shower->start_segment()->set_flags(SegmentFlags::kAvoidMuonCheck);
            double kine_charge = cal_kine_charge(max_shower, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
            max_shower->set_kine_charge(kine_charge);
            max_shower->set_flag_kinematics(true);
        }
        
        // Remove deleted showers
        for (auto shower : del_showers) {
            showers.erase(shower);
        }
        
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                          map_vertex_to_shower, used_shower_clusters);
    }
}


void PatternAlgorithms::examine_showers(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    
    if (!main_vertex) return;

    // Build map_vertex_segments
    std::map<VertexPtr, std::vector<SegmentPtr>> map_vertex_segments;
    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg) continue;
        auto source_vdesc = boost::source(e, graph);
        auto target_vdesc = boost::target(e, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;
        if (v1) map_vertex_segments[v1].push_back(seg);
        if (v2) map_vertex_segments[v2].push_back(seg);
    }

    WireCell::Point  main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
    WireCell::Vector drift_dir(1, 0, 0);

    // ---- Shower comparator by start-segment graph index (deterministic ordering) ----
    auto shower_cmp = [&](const ShowerPtr& a, const ShowerPtr& b) -> bool {
        auto sa = a->start_segment(), sb = b->start_segment();
        if (!sa && !sb) return false;
        if (!sa) return true;
        if (!sb) return false;
        return graph[sa->get_descriptor()].index < graph[sb->get_descriptor()].index;
    };

    // ---- Per-shower cached data (EM electrons only) ----
    // Pre-computing dir_internal, dir_from_main, composition, etc. avoids
    // O(N_seg * N_shower) calls to shower_cal_dir_3vector / fill_maps.
    struct EMShowerCache {
        ShowerPtr shower;
        int conn_type;
        VertexPtr start_vtx;
        double Eshower;               // kine_best (if non-zero) else kine_charge — for Case I
        WireCell::Vector dir_internal;   // shower_cal_dir_3vector at 100 cm
        WireCell::Vector dir_from_main;  // start_pt - main_vtx_pt (for Cases II/III)
        double angle_drift_internal;     // angle(dir_internal, drift_dir) in degrees
        double angle_fm_vs_int;          // angle(dir_from_main, dir_internal) in degrees (for Case III)
        bool composition_ok;             // passes track-fraction composition check (for Case I)
    };

    auto build_em_cache = [&](const IndexedShowerSet& shower_set) -> std::vector<EMShowerCache> {
        std::vector<ShowerPtr> sorted(shower_set.begin(), shower_set.end());
        std::sort(sorted.begin(), sorted.end(), shower_cmp);
        std::vector<EMShowerCache> result;
        result.reserve(sorted.size());
        for (const auto& sh : sorted) {
            if (!sh->start_segment() || !sh->start_segment()->has_particle_info() ||
                sh->start_segment()->particle_info()->pdg() != 11) continue;
            auto [sv, ct] = sh->get_start_vertex_and_type();
            double E      = sh->get_kine_best() != 0 ? sh->get_kine_best() : sh->get_kine_charge();
            WireCell::Point  sp    = sh->get_start_point();
            WireCell::Vector dir_i = shower_cal_dir_3vector(*sh, sp, 100 * units::cm);
            WireCell::Vector dir_fm(sp.x() - main_vtx_pt.x(), sp.y() - main_vtx_pt.y(), sp.z() - main_vtx_pt.z());
            double ang_drift = std::acos(std::clamp(dir_i.dot(drift_dir) / (dir_i.magnitude() * drift_dir.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            double ang_fm_i  = std::acos(std::clamp(dir_fm.dot(dir_i)   / (dir_fm.magnitude() * dir_i.magnitude()),    -1.0, 1.0)) / M_PI * 180.0;
            // Composition check: track-fraction within shower's start cluster
            TrajectoryView& traj = sh->fill_maps();
            double ttl = 0, ttrk = 0;
            // ordered_edges: ttl/ttrk are FP accumulations feeding the
            // `ttrk > 0.25 * ttl` composition test just below (doc pr/28 §15.2).
            for (auto edesc : ordered_edges(traj, graph)) {
                auto sg1 = traj.view_graph()[edesc].segment;
                if (!sg1 || sg1->cluster() != sh->start_segment()->cluster()) continue;
                double len = segment_track_length(sg1);
                ttl += len;
                if (!seg_dir_weak(sg1)) ttrk += len;
            }
            bool comp_ok = !(ttrk > 3 * units::cm && ttrk > 0.25 * ttl);
            result.push_back({sh, ct, sv, E, dir_i, dir_fm, ang_drift, ang_fm_i, comp_ok});
        }
        return result;
    };

    // Pre-cache Case II (showers at main_vertex) and Case III (all showers) — built once
    std::vector<EMShowerCache> cache_main_vtx;
    if (map_vertex_to_shower.count(main_vertex))
        cache_main_vtx = build_em_cache(map_vertex_to_shower[main_vertex]);
    std::vector<EMShowerCache> cache_all_showers = build_em_cache(showers);

    // Per-vertex lazy cache for Case I (populated on first access per vertex)
    std::map<VertexPtr, std::vector<EMShowerCache>> vtx_cache;

    // SegmentIndexCmp for deterministic merge-processing order
    std::map<SegmentPtr, ShowerPtr, SegmentIndexCmp> map_merge_seg_shower;
    std::set<ShowerPtr> del_showers;

    // ---- Segment loop ----
    if (map_vertex_segments.count(main_vertex)) {
        for (auto sg : map_vertex_segments[main_vertex]) {
            if (map_segment_in_shower.count(sg)) {
                // doc pr/33 F2: the prototype's muon test is EXACT here
                // (NeutrinoID_em_shower.h:10, get_particle_type()!=13), so a
                // pdg -13 segment is skipped there and processed by the abs()
                // port.  The !has_particle_info() term is NOT a divergence
                // (no info => type 0 => continue in both trees).
                const bool skip_legacy = !sg->has_particle_info() || std::abs(sg->particle_info()->pdg()) != 13;
                const bool skip_proto  = !sg->has_particle_info() || sg->particle_info()->pdg() != 13;
                g_pr33_audit.f2_calls[3]++;
                if (skip_legacy != skip_proto) g_pr33_audit.f2_disagree[3]++;
                if (m_shower_pdg_exact_muon_test ? skip_proto : skip_legacy) continue;
            }

            double sg_length = segment_track_length(sg);
            if ((sg_length > 45 * units::cm && !seg_dir_weak(sg)) || sg_length > 55 * units::cm) continue;

            // doc pr/93 round 4 (straight_cont_cross_cluster): a stem the
            // demotion pass claimed as the head of a cross-cluster muon
            // chain must not be re-adopted (and re-stamped pdg 11 with the
            // score-100 sentinel, which shower_accept_pid_guard cannot
            // decline) by this retarget.  Empty set when the knob is off =>
            // byte-identical.
            if (m_straight_cont_cross_cluster && m_sccc_shield_segs.count(sg)) continue;

            auto seg_edesc = sg->get_descriptor();
            VertexPtr v1   = graph[boost::source(seg_edesc, graph)].vertex;
            VertexPtr v2   = graph[boost::target(seg_edesc, graph)].vertex;
            VertexPtr vtx  = (v1 == main_vertex) ? v2 : v1;
            if (!vtx) continue;

            // doc pr/33 F1: the prototype calls calculate_num_daughter_tracks
            // here (NeutrinoID_em_shower.h:17, flag=false => track-only
            // lengths, length > 0); the port's _showers(...,false) sums every
            // daughter segment.  Both computed unconditionally for the counter.
            auto daughter_legacy = calculate_num_daughter_showers(graph, main_vertex, sg, false);
            auto daughter_proto  = calculate_num_daughter_tracks(graph, main_vertex, sg, false, 0);
            g_pr33_audit.f1_ex_calls++;
            if (daughter_legacy.second != daughter_proto.second) g_pr33_audit.f1_ex_differ++;
            double daughter_length = m_daughter_count_proto_examine_showers ? daughter_proto.second
                                                                            : daughter_legacy.second;

            // Pre-compute segment directions shared across Cases I, II, III
            WireCell::Point  vtx_pt           = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            WireCell::Vector dir_seg_from_vtx  = segment_cal_dir_3vector(sg, vtx_pt,      15 * units::cm);
            WireCell::Vector dir_seg_from_main = segment_cal_dir_3vector(sg, main_vtx_pt, 15 * units::cm);

            bool flag_checked = false;

            // ---- Case I: showers at the other vertex ----
            if (map_vertex_to_shower.count(vtx)) {
                // Lazy-populate vtx cache
                if (!vtx_cache.count(vtx))
                    vtx_cache[vtx] = build_em_cache(map_vertex_to_shower[vtx]);
                const auto& vtx_em = vtx_cache[vtx];

                // First pass: any high-energy directly connected shower? (composition not required)
                bool flag_tmp_connected = false;
                for (const auto& c : vtx_em) {
                    if (c.conn_type == 1 && c.Eshower > 60 * units::MeV) { flag_tmp_connected = true; break; }
                }

                // Second pass: evaluate merge conditions
                for (const auto& c : vtx_em) {
                    if (!c.composition_ok) continue;

                    if (c.conn_type == 1 && c.Eshower > 100 * units::MeV) flag_checked = true;

                    double angle1     = std::acos(std::clamp(dir_seg_from_vtx.dot(c.dir_internal)  / (dir_seg_from_vtx.magnitude()  * c.dir_internal.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle2     = std::acos(std::clamp(dir_seg_from_main.dot(c.dir_internal) / (dir_seg_from_main.magnitude() * c.dir_internal.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double tmp_angle  = std::min(180.0 - angle1, angle2);

                    // Lazy closest-distance: only computed when conn_type==2 conditions are reached
                    double closest_dis = -1.0;
                    auto get_closest_dis = [&]() -> double {
                        if (closest_dis < 0) {
                            auto [cd, cp] = segment_get_closest_point(sg, c.shower->get_start_point());
                            closest_dis = cd;
                        }
                        return closest_dis;
                    };

                    if (c.conn_type == 2) {
                        if ((!seg_dir_weak(sg) && sg_length > 3 * units::cm) || flag_tmp_connected) {
                            if (get_closest_dis() < 8 * units::cm && c.Eshower > 75 * units::MeV && tmp_angle < 6) {
                                // Allow this condition
                            } else {
                                continue;
                            }
                        }
                        if (get_closest_dis() > 20 * units::cm && c.Eshower < 150 * units::MeV && tmp_angle > 2.5) continue;
                    }

                    if ((c.Eshower > 800 * units::MeV && tmp_angle < 30) ||
                        (c.Eshower > 150 * units::MeV && tmp_angle < 10) ||
                        (c.Eshower > 150 * units::MeV && tmp_angle < 18 && c.conn_type == 1 && seg_dir_weak(sg)) ||
                        (c.Eshower > 100 * units::MeV && tmp_angle < 10 && sg_length < 25 * units::cm) ||
                        (c.Eshower > 250 * units::MeV && tmp_angle < 15) ||
                        (c.Eshower > 360 * units::MeV && tmp_angle < 25) ||
                        (c.Eshower > 100 * units::MeV && c.Eshower <= 150 * units::MeV && tmp_angle < 15 &&
                         sg_length < 25 * units::cm && flag_checked) ||
                        (c.Eshower > 60 * units::MeV && c.conn_type == 2 && seg_dir_weak(sg) &&
                         ((tmp_angle < 15   && get_closest_dis() < 18 * units::cm) ||
                          (tmp_angle < 17.5 && get_closest_dis() <  6 * units::cm)) &&
                         sg_length < 15 * units::cm) ||
                        (c.Eshower > 60 * units::MeV && c.conn_type == 2 && tmp_angle < 7.5 &&
                         get_closest_dis() < 8 * units::cm && sg_length < 20 * units::cm)) {
                        map_merge_seg_shower[sg] = c.shower;
                        continue;
                    }
                }
            }

            if (map_merge_seg_shower.count(sg)) continue;
            if (flag_checked) continue;
            if (!seg_dir_weak(sg) && sg_length > 6 * units::cm && daughter_length < 40 * units::cm) continue;

            // ---- Case II: showers at main_vertex (indirect connections only) ----
            {
                double angle_drift1 = std::acos(std::clamp(dir_seg_from_main.dot(drift_dir) / (dir_seg_from_main.magnitude() * drift_dir.magnitude()), -1.0, 1.0)) / M_PI * 180.0;

                for (const auto& c : cache_main_vtx) {
                    if (c.conn_type == 1) continue;

                    double angle_dir2 = std::acos(std::clamp(dir_seg_from_main.dot(c.dir_from_main) / (dir_seg_from_main.magnitude() * c.dir_from_main.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle_dir3 = std::acos(std::clamp(dir_seg_from_main.dot(c.dir_internal)  / (dir_seg_from_main.magnitude() * c.dir_internal.magnitude()),  -1.0, 1.0)) / M_PI * 180.0;

                    auto [min_dis, closest_pt] = segment_get_closest_point(sg, c.shower->get_start_point());

                    if (((c.shower->get_kine_charge() > 80 * units::MeV && angle_dir2 < 10) ||
                         (c.shower->get_kine_charge() > 50 * units::MeV && angle_dir2 < 3) ||
                         (c.shower->get_kine_charge() > 80 * units::MeV && angle_dir3 < 6 && angle_dir2 < 17.5) ||
                         (c.shower->get_kine_charge() > 80 * units::MeV && angle_dir3 < 6 &&
                          std::fabs(90 - angle_drift1) < 10 && std::fabs(90 - c.angle_drift_internal) < 10 && angle_dir2 < 30)) &&
                        (sg_length > 5 * units::cm || (sg_length > 3 * units::cm && min_dis < 2.0 * units::cm))) {
                        map_merge_seg_shower[sg] = c.shower;
                        continue;
                    }
                }
            }

            // ---- Case III: other showers not at vtx or main_vertex ----
            {
                for (const auto& c : cache_all_showers) {
                    if (c.start_vtx == vtx || c.start_vtx == main_vertex) continue;

                    if (c.conn_type <= 2) {
                        double angle_dir2 = std::acos(std::clamp(dir_seg_from_main.dot(c.dir_from_main) / (dir_seg_from_main.magnitude() * c.dir_from_main.magnitude()), -1.0, 1.0)) / M_PI * 180.0;

                        if ((((c.shower->get_kine_charge() > 80 * units::MeV && angle_dir2 < 15) ||
                              (c.shower->get_kine_charge() > 50 * units::MeV && angle_dir2 < 5)) &&
                             c.angle_fm_vs_int < 15 && (angle_dir2 + c.angle_fm_vs_int) < 25 && sg_length > 5 * units::cm)) {
                            map_merge_seg_shower[sg] = c.shower;
                            continue;
                        }
                    }
                }
            }
        }
    }

    // Mark showers for deletion if their start segment should be merged
    for (auto shower : showers) {
        SegmentPtr sg1 = shower->start_segment();
        if (sg1 && map_merge_seg_shower.count(sg1)) {
            sg1->set_flags(SegmentFlags::kAvoidMuonCheck);
            del_showers.insert(shower);
        }
    }

    // Delete marked showers first
    if (!del_showers.empty()) {
        for (auto shower1 : del_showers) showers.erase(shower1);
        del_showers.clear();
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                          map_vertex_to_shower, used_shower_clusters);
    }

    // Sort showers once for deterministic sub-shower absorption in the merge loop
    std::vector<ShowerPtr> showers_sorted_merge(showers.begin(), showers.end());
    std::sort(showers_sorted_merge.begin(), showers_sorted_merge.end(), shower_cmp);

    // Perform the merging; track updated showers in stable index order for post-merge pass
    std::vector<ShowerPtr> updated_showers_vec;
    std::set<ShowerPtr>    updated_showers_set;

    for (auto& [sg, shower] : map_merge_seg_shower) {
        SPDLOG_LOGGER_TRACE(s_log, "examine_showers: EM shower modification: {} -> {}", shower->start_segment()->id(), sg->id());
        if (!updated_showers_set.count(shower)) {
            updated_showers_vec.push_back(shower);
            updated_showers_set.insert(shower);
        }

        auto [pair_vertex, pair_conn_type] = shower->get_start_vertex_and_type();

        pr93_probe_absorb_direct("examine_showers_retarget_seed", shower, sg);
        shower->add_segment(sg);
        shower->set_start_vertex(main_vertex, 1);
        shower->set_start_segment(sg);
        shower->set_start_point(main_vtx_pt);
        IndexedSegmentSet tmp_used_segments;
        pr93_probe_absorb_site("examine_showers_retarget", shower, m_shower_absorb_track_guard);
        shower->complete_structure_with_start_segment(nv_bridge_seed(tmp_used_segments), "fit", "associate_points", m_shower_absorb_track_guard, m_shower_walk_visited_parity);  // doc pr/40 round 9 B2: shield pre-seed, no-op when bridge off
        if (pair_conn_type != 1) {
            if (segment_track_length(sg) > 44 * units::cm || seg_dir_weak(sg))
                sg->set_flags(SegmentFlags::kAvoidMuonCheck);
        } else {
            if (shower->get_num_main_segments() >= 3)
                sg->set_flags(SegmentFlags::kAvoidMuonCheck);
        }

        // Absorb other showers whose start vertex is now inside this shower
        TrajectoryView& traj = shower->fill_maps();
        std::set<VertexPtr> shower_vertices;
        for (auto vdesc : traj.nodes()) {
            auto vtx = traj.view_graph()[vdesc].vertex;
            if (vtx) shower_vertices.insert(vtx);
        }
        for (auto shower1 : showers_sorted_merge) {
            if (shower == shower1) continue;
            auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
            if (conn_type1 == 1 && start_vtx1 != main_vertex && shower_vertices.count(start_vtx1)) {
                pr93_probe_absorb_splice("examine_showers_vtxcontain", shower, shower1);
                shower->add_shower(*shower1);
                del_showers.insert(shower1);
            }
        }

        if (sg->has_particle_info() && sg->particle_info()) {
            // doc sbnd_xin/docs/pr/93 Cause B.  Same guard as
            // new_shower_accepted above (incl. the 50cm floor -- 8 of the 9
            // nueCC48 regressions without it were THIS site relabelling a
            // real electron's 22-47cm confidently-mis-scored stem); pdg
            // write only.  C++ default false => byte-identical.
            if (m_shower_accept_pid_guard &&
                segment_track_length(sg) > m_shower_pid_guard_min_len &&
                segment_confident_nonelectron_pid(sg)) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr93 accept_pid_guard: decline e- write seg={} pdg={} score={:.3f} (merged_shower_start_segment)",
                    sg->id(), sg->particle_info()->pdg(), sg->particle_score());
            }
            else {
                sg->particle_info()->set_pdg(11);
                pr40_probe_setpdg(sg, 11, "NeutrinoShowerClustering.cxx:merged_shower_start_segment");
            }
        }
        shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
        shower->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
        shower->set_kine_charge(cal_kine_charge(shower, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv));
        shower->set_flag_kinematics(true);
    }

    // Delete merged showers (not at main_vertex)
    for (auto shower1 : del_showers) {
        auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
        if (start_vtx1 != main_vertex) showers.erase(shower1);
    }
    del_showers.clear();

    // Post-merge pass: check updated showers against remaining showers (sorted for determinism)
    std::vector<ShowerPtr> showers_sorted_post(showers.begin(), showers.end());
    std::sort(showers_sorted_post.begin(), showers_sorted_post.end(), shower_cmp);

    for (auto shower : updated_showers_vec) {
        WireCell::Point  shower_start = shower->get_start_point();
        WireCell::Vector dir1 = shower_cal_dir_3vector(*shower, shower_start, 25 * units::cm);

        for (auto shower1 : showers_sorted_post) {
            if (updated_showers_set.count(shower1)) continue;
            // Fixed: correct particle-type filter (was inverted — admitted unclassified showers)
            if (!shower1->start_segment() || !shower1->start_segment()->has_particle_info() ||
                shower1->start_segment()->particle_info()->pdg() != 11) continue;

            auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
            if (conn_type1 != 2) continue;
            if (del_showers.count(shower1)) continue;

            auto [shower_vtx, shower_conn] = shower->get_start_vertex_and_type();
            WireCell::Point shower_vtx_pt = shower_vtx->fit().valid() ? shower_vtx->fit().point : shower_vtx->wcpt().point;

            WireCell::Vector dir2(shower1->get_start_point().x() - shower_vtx_pt.x(),
                                 shower1->get_start_point().y() - shower_vtx_pt.y(),
                                 shower1->get_start_point().z() - shower_vtx_pt.z());
            WireCell::Vector dir3 = shower_cal_dir_3vector(*shower1, shower1->get_start_point(), 25 * units::cm);

            double angle_dir2 = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            double angle_dir3 = std::acos(std::clamp(dir1.dot(dir3) / (dir1.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;

            if (angle_dir2 < 10 && angle_dir3 < 20) {
                pr93_probe_absorb_splice("examine_showers_angle", shower, shower1);
                shower->add_shower(*shower1);
                shower->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
                shower->set_kine_charge(cal_kine_charge(shower, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv));
                shower->set_flag_kinematics(true);
                del_showers.insert(shower1);
            }
        }
    }

    // Final deletion
    for (auto shower1 : del_showers) showers.erase(shower1);

    if (!map_merge_seg_shower.empty()) {
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                          map_vertex_to_shower, used_shower_clusters);
    }

    examine_shower_1(graph, main_vertex, showers, main_cluster, other_clusters, map_cluster_main_vertices,
                    map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters,
                    track_fitter, dv, particle_data, recomb_model);
}


void PatternAlgorithms::id_pi0_with_vertex(int& acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex) return;

    // -- Deterministic comparators (graph-index order, not pointer address) --
    auto shower_idx = [&](const ShowerPtr& s) -> size_t {
        auto seg = s->start_segment();
        return (seg && seg->descriptor_valid())
            ? graph[seg->get_descriptor()].index
            : std::numeric_limits<size_t>::max();
    };
    auto shower_cmp = [&](const ShowerPtr& a, const ShowerPtr& b) {
        return shower_idx(a) < shower_idx(b);
    };
    // Orders shower pairs (a1,a2) vs (b1,b2) by graph index of each element.
    auto shower_pair_cmp = [&](const std::pair<ShowerPtr,ShowerPtr>& a,
                                const std::pair<ShowerPtr,ShowerPtr>& b) {
        size_t ia1 = shower_idx(a.first),  ia2 = shower_idx(a.second);
        size_t ib1 = shower_idx(b.first),  ib2 = shower_idx(b.second);
        if (ia1 != ib1) return ia1 < ib1;
        return ia2 < ib2;
    };

    WireCell::Point main_vtx_pt = main_vertex->fit().valid()
        ? main_vertex->fit().point : main_vertex->wcpt().point;

    // -- Per-shower caches: avoid repeated calls inside loops ---------------
    // Covers all showers reachable through map_vertex_to_shower.
    std::map<ShowerPtr, std::pair<VertexPtr,int>> svc;  // {start_vertex, conn_type}
    std::map<ShowerPtr, double>                   kqc;  // kine_charge
    for (auto& [vtx, shower_set] : map_vertex_to_shower) {
        for (auto sh : shower_set) {
            if (!svc.count(sh)) {
                svc[sh] = sh->get_start_vertex_and_type();
                kqc[sh] = sh->get_kine_charge();
            }
        }
    }
    // Fallback to live call if a shower somehow escaped the cache.
    auto get_svc = [&](ShowerPtr sh) -> std::pair<VertexPtr,int> {
        auto it = svc.find(sh);
        return it != svc.end() ? it->second : sh->get_start_vertex_and_type();
    };
    auto get_kq = [&](ShowerPtr sh) -> double {
        auto it = kqc.find(sh);
        return it != kqc.end() ? it->second : sh->get_kine_charge();
    };

    // -- Build map_vertex_segments (ordered_edges gives deterministic order) --
    std::map<VertexPtr, std::vector<SegmentPtr>> map_vertex_segments;
    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (!seg) continue;
        VertexPtr v1 = graph[boost::source(e, graph)].vertex;
        VertexPtr v2 = graph[boost::target(e, graph)].vertex;
        if (v1) map_vertex_segments[v1].push_back(seg);
        if (v2) map_vertex_segments[v2].push_back(seg);
    }

    // Sort map_vertex_to_shower keys once; used for all loops below.
    std::vector<VertexPtr> vtx_sorted;
    vtx_sorted.reserve(map_vertex_to_shower.size());
    for (auto& [v, _] : map_vertex_to_shower) vtx_sorted.push_back(v);
    std::sort(vtx_sorted.begin(), vtx_sorted.end(), VertexIndexCmp{});

    // -- Disconnected showers (conn_type==2, non-muon), deterministic order --
    // Using set<..., shower_cmp> gives automatic dedup and graph-index ordering.
    std::set<ShowerPtr, decltype(shower_cmp)> disconnected_showers(shower_cmp);
    std::map<ShowerPtr, WireCell::Vector> map_shower_dir;  // lookup-only; pointer key OK

    // -- Candidate vertices in deterministic graph-index order ---------------
    IndexedVertexSet candidate_vertices;
    candidate_vertices.insert(main_vertex);

    // Single pass over map_vertex_to_shower (merged from the original two passes).
    for (auto vtx : vtx_sorted) {
        // candidate_vertices logic (was the second pass)
        bool flag_add = true;
        auto it_in = map_vertex_in_shower.find(vtx);
        if (it_in != map_vertex_in_shower.end()) {
            flag_add = false;
            auto [sv, ct] = get_svc(it_in->second);
            if (vtx == sv) flag_add = true;
        }
        if (flag_add) candidate_vertices.insert(vtx);

        // disconnected_showers / map_shower_dir logic (was the first pass)
        std::vector<ShowerPtr> sorted_showers(map_vertex_to_shower.at(vtx).begin(),
                                              map_vertex_to_shower.at(vtx).end());
        std::sort(sorted_showers.begin(), sorted_showers.end(), shower_cmp);
        for (auto shower : sorted_showers) {
            auto [start_vtx, conn_type] = get_svc(shower);
            if (conn_type == 2 && std::abs(shower->get_particle_type()) != 13) {
                disconnected_showers.insert(shower);
                if (!map_shower_dir.count(shower))
                    map_shower_dir[shower] = shower_cal_dir_3vector(
                        *shower, shower->get_start_point(), 15 * units::cm);
            } else if (conn_type == 1) {
                if (!map_shower_dir.count(shower)) {
                    WireCell::Point sp = start_vtx->fit().valid()
                        ? start_vtx->fit().point : start_vtx->wcpt().point;
                    map_shower_dir[shower] = shower_cal_dir_3vector(*shower, sp, 15 * units::cm);
                }
            }
        }
    }

    // -- Map shower pairs → masses and candidate vertices -------------------
    // Custom comparator eliminates pointer-address ordering of pair keys.
    std::map<std::pair<ShowerPtr,ShowerPtr>,
             std::vector<std::pair<double,VertexPtr>>,
             decltype(shower_pair_cmp)> map_shower_pair_mass_vertex(shower_pair_cmp);

    for (auto cand_vtx : candidate_vertices) {
        std::vector<ShowerPtr> tmp_showers;
        std::map<ShowerPtr, WireCell::Vector> local_dirs;

        WireCell::Point vtx_pt = cand_vtx->fit().valid()
            ? cand_vtx->fit().point : cand_vtx->wcpt().point;

        // Add directly connected showers (conn_type==1, non-muon), sorted.
        auto it2 = map_vertex_to_shower.find(cand_vtx);
        if (it2 != map_vertex_to_shower.end()) {
            std::vector<ShowerPtr> sv(it2->second.begin(), it2->second.end());
            std::sort(sv.begin(), sv.end(), shower_cmp);
            for (auto shower : sv) {
                auto [start_vtx, conn_type] = get_svc(shower);
                if (conn_type == 1 && std::abs(shower->get_particle_type()) != 13) {
                    tmp_showers.push_back(shower);
                    local_dirs[shower] = shower->get_init_dir();
                }
            }
        }

        // Add disconnected showers within 30° — set iteration is already deterministic.
        for (auto shower : disconnected_showers) {
            WireCell::Vector dir1 = map_shower_dir[shower];
            WireCell::Vector dir2(shower->get_start_point().x() - vtx_pt.x(),
                                  shower->get_start_point().y() - vtx_pt.y(),
                                  shower->get_start_point().z() - vtx_pt.z());
            auto [start_vtx, conn_type] = get_svc(shower);
            if (start_vtx == cand_vtx) {
                tmp_showers.push_back(shower);
                local_dirs[shower] = shower->get_init_dir();
            } else {
                double angle = std::acos(std::clamp(
                    dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                if (angle < 30) {
                    tmp_showers.push_back(shower);
                    local_dirs[shower] = dir2;
                }
            }
        }

        // Sort so (i < j) pairs are in graph-index order regardless of insertion order.
        std::sort(tmp_showers.begin(), tmp_showers.end(), shower_cmp);

        // Compute pi0 mass for each (i < j) pair.
        // Early-skip if both conn_type==1 avoids the acos+sqrt for ineligible pairs.
        for (size_t i = 0; i < tmp_showers.size(); i++) {
            ShowerPtr sh1 = tmp_showers[i];
            auto [sv1, ct1] = get_svc(sh1);
            WireCell::Vector dir1 = local_dirs[sh1];
            double kq1 = get_kq(sh1);

            for (size_t j = i + 1; j < tmp_showers.size(); j++) {
                ShowerPtr sh2 = tmp_showers[j];
                auto [sv2, ct2] = get_svc(sh2);
                if (ct1 == 1 && ct2 == 1) continue;  // ineligible — skip before expensive ops

                WireCell::Vector dir2 = local_dirs[sh2];
                double angle    = std::acos(std::clamp(
                    dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0));
                double mass_pio = std::sqrt(4 * kq1 * get_kq(sh2) * std::pow(std::sin(angle / 2.0), 2));
                map_shower_pair_mass_vertex[{sh1, sh2}].push_back({mass_pio, cand_vtx});
            }
        }
    }

    // -- Compute BDT features: highest-energy pair with a valid decay vertex --
    {
        double max_energy = 0;
        for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex) {
            ShowerPtr sh1 = shower_pair.first;
            ShowerPtr sh2 = shower_pair.second;
            double energy_1 = get_kq(sh1);
            double energy_2 = get_kq(sh2);
            if (energy_1 + energy_2 <= max_energy) continue;

            // Hoist start-vertex lookup out of the mass_vtx_vec inner loop.
            auto [sv1, ct1] = get_svc(sh1);
            auto [sv2, ct2] = get_svc(sh2);

            double best_vtx_dis = 1000 * units::cm;
            for (auto& [mass_pio, cand_vtx] : mass_vtx_vec) {
                if (cand_vtx != sv1 && cand_vtx != sv2) continue;

                WireCell::Point vtx_pt = cand_vtx->fit().valid()
                    ? cand_vtx->fit().point : cand_vtx->wcpt().point;
                double temp_dis = (vtx_pt - main_vtx_pt).magnitude();
                if (temp_dis >= best_vtx_dis) continue;

                best_vtx_dis = temp_dis;
                max_energy   = energy_1 + energy_2;

                WireCell::Point sp1 = sh1->get_start_point();
                WireCell::Point sp2 = sh2->get_start_point();
                double dis_1 = (vtx_pt - sp1).magnitude();
                double dis_2 = (vtx_pt - sp2).magnitude();

                // Prototype uses each shower's own start vertex for direction
                // when close (not the candidate pi0 vertex).
                WireCell::Point sv1_pt = sv1->fit().valid() ? sv1->fit().point : sv1->wcpt().point;
                WireCell::Point sv2_pt = sv2->fit().valid() ? sv2->fit().point : sv2->wcpt().point;
                WireCell::Vector dir1 = (dis_1 < 3 * units::cm)
                    ? shower_cal_dir_3vector(*sh1, sv1_pt, 15 * units::cm)
                    : WireCell::Vector(sp1.x()-vtx_pt.x(), sp1.y()-vtx_pt.y(), sp1.z()-vtx_pt.z());
                WireCell::Vector dir2 = (dis_2 < 3 * units::cm)
                    ? shower_cal_dir_3vector(*sh2, sv2_pt, 15 * units::cm)
                    : WireCell::Vector(sp2.x()-vtx_pt.x(), sp2.y()-vtx_pt.y(), sp2.z()-vtx_pt.z());

                pio_kine.flag     = 1;
                pio_kine.mass     = mass_pio;
                pio_kine.vtx_dis  = temp_dis;
                pio_kine.energy_1 = energy_1;
                pio_kine.energy_2 = energy_2;
                pio_kine.dis_1    = dis_1;
                pio_kine.dis_2    = dis_2;
                pio_kine.theta_1  = std::acos(std::clamp(dir1.z() / dir1.magnitude(), -1.0, 1.0));
                pio_kine.phi_1    = std::atan2(dir1.y(), dir1.x());
                pio_kine.theta_2  = std::acos(std::clamp(dir2.z() / dir2.magnitude(), -1.0, 1.0));
                pio_kine.phi_2    = std::atan2(dir2.y(), dir2.x());
                pio_kine.angle    = std::acos(std::clamp(
                    dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0));
            }
        }
    }

    // -- Find pi0 candidates iteratively ------------------------------------
    while (!map_shower_pair_mass_vertex.empty()) {
        double mass_diff = 1e9;
        double mass_save = 0;
        ShowerPtr shower_1 = nullptr;
        ShowerPtr shower_2 = nullptr;
        const double mass_offset = 10 * units::MeV;
        VertexPtr vtx = nullptr;
        double mass_penalty = 0;

        // Hoist start-vertex lookup out of the mass_vtx_vec inner loop.
        for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex) {
            auto [sv1, ct1] = get_svc(shower_pair.first);
            auto [sv2, ct2] = get_svc(shower_pair.second);
            double tmp_penalty = (ct1 == 2 && ct2 == 2) ? 6 * units::MeV : 0.0;

            for (auto& [mass, candidate_vtx] : mass_vtx_vec) {
                double delta = mass - 135 * units::MeV + mass_offset;
                if (delta >= 35 * units::MeV || delta <= -25 * units::MeV) continue;
                if (std::abs(delta) - tmp_penalty < std::abs(mass_diff) - mass_penalty) {
                    mass_diff    = delta;
                    mass_penalty = tmp_penalty;
                    mass_save    = mass;
                    shower_1     = shower_pair.first;
                    shower_2     = shower_pair.second;
                    vtx          = candidate_vtx;
                }
            }
        }

        if (mass_diff >= 35 * units::MeV || mass_diff <= -25 * units::MeV) break;

        pi0_showers.insert(shower_1);
        pi0_showers.insert(shower_2);

        int pio_id = acc_segment_id++;
        g_pr33_audit.f3_pi0_with_vertex++;
        map_shower_pio_id[shower_1]  = pio_id;
        map_shower_pio_id[shower_2]  = pio_id;
        map_pio_id_mass[pio_id]      = {mass_save, 1};
        map_pio_id_showers[pio_id].push_back(shower_1);
        map_pio_id_showers[pio_id].push_back(shower_2);

        auto [sv1, ct1] = get_svc(shower_1);
        if (sv1 != vtx) {
            shower_1->set_start_vertex(vtx, 2);
            shower_1->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
            svc[shower_1] = {vtx, 2};  // keep cache consistent
        }
        auto [sv2, ct2] = get_svc(shower_2);
        if (sv2 != vtx) {
            shower_2->set_start_vertex(vtx, 2);
            shower_2->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
            svc[shower_2] = {vtx, 2};
        }

        SPDLOG_LOGGER_TRACE(s_log, "examine_showers: Pi0 found with mass: {} MeV with {} MeV + {} MeV",
            mass_save / units::MeV, shower_1->get_kine_charge() / units::MeV, shower_2->get_kine_charge() / units::MeV);

        // Remove all pairs that involve either used shower.
        std::vector<std::pair<ShowerPtr,ShowerPtr>> to_remove;
        for (auto& [sp, _] : map_shower_pair_mass_vertex)
            if (sp.first == shower_1 || sp.first == shower_2 ||
                sp.second == shower_1 || sp.second == shower_2)
                to_remove.push_back(sp);
        for (auto& p : to_remove) map_shower_pair_mass_vertex.erase(p);
    }

    // -- Reclassify incoming muons/unknowns at pi0 decay vertices as pions --
    IndexedVertexSet pi0_vertices;
    for (auto shower : pi0_showers) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        pi0_vertices.insert(start_vtx);
    }

    for (auto pi0_vtx : pi0_vertices) {
        if (!map_vertex_segments.count(pi0_vtx)) continue;
        WireCell::Point vtx_pt = pi0_vtx->fit().valid()
            ? pi0_vtx->fit().point : pi0_vtx->wcpt().point;

        for (auto sg : map_vertex_segments[pi0_vtx]) {
            bool flag_start = false;
            if (!sg->fits().empty()) {
                double dist_front = (sg->fits().front().point - vtx_pt).magnitude();
                double dist_back  = (sg->fits().back().point  - vtx_pt).magnitude();
                flag_start = (dist_front < dist_back);
            }
            int dirsign_val = sg->dirsign();
            bool is_incoming = (flag_start && dirsign_val == -1) || (!flag_start && dirsign_val == 1);
            if (is_incoming && sg->has_particle_info()) {
                int pdg = sg->particle_info()->pdg();
                if (std::abs(pdg) == 13 || pdg == 0) {
                    sg->particle_info()->set_pdg(211);
                    sg->particle_info()->set_mass(139.57 * units::MeV);
                    // Recalculate 4-momentum under pion hypothesis, mirroring prototype:
                    // if (sg->get_particle_4mom(3)>0) sg->cal_4mom();
                    if (sg->particle_info()->kinetic_energy() > 0) {
                        auto four_momentum = segment_cal_4mom(sg, 211, particle_data, recomb_model, m_mip_dqdx);
                        sg->particle_info()->set_four_momentum(four_momentum);
                    }
                }
            }
        }
    }
}


void PatternAlgorithms::id_pi0_without_vertex(int& acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex) return;

    // Build main_vertex_segs by direct graph neighbor query: O(degree·log degree)
    // instead of scanning all edges via ordered_edges: O(E·log E).
    // Only main_vertex's segments are needed in this function.
    std::vector<SegmentPtr> main_vertex_segs;
    {
        auto vd = main_vertex->get_descriptor();
        for (auto eit : sorted_out_edges(vd, graph)) {
            SegmentPtr seg = graph[eit].segment;
            if (seg) main_vertex_segs.push_back(seg);
        }
        std::sort(main_vertex_segs.begin(), main_vertex_segs.end(), [&graph](const SegmentPtr& a, const SegmentPtr& b) {
            return graph[a->get_descriptor()].index < graph[b->get_descriptor()].index;
        });
    }

    // Check main vertex conditions
    if (main_vertex_segs.size() > 2) return;

    if (!main_vertex_segs.empty()) {
        auto first_seg = main_vertex_segs.front();
        auto last_seg = main_vertex_segs.back();
        
        if ((map_segment_in_shower.find(first_seg) == map_segment_in_shower.end() &&
             map_segment_in_shower.find(last_seg) == map_segment_in_shower.end()) ||
            segments_in_long_muon.find(first_seg) != segments_in_long_muon.end() ||
            segments_in_long_muon.find(last_seg) != segments_in_long_muon.end()) {
            return;
        }
    }
    
    // Stable shower comparator: order by start_segment graph index to eliminate
    // pointer-address-based non-determinism in all shower containers below.
    auto shower_less = [this, &graph](const ShowerPtr& a, const ShowerPtr& b) {
        auto sa = a->start_segment();
        auto sb = b->start_segment();
        if (sa && sb) {
            size_t ia = graph[sa->get_descriptor()].index;
            size_t ib = graph[sb->get_descriptor()].index;
            if (ia != ib) return ia < ib;
        }
        else if (!sa &&  sb) return true;
        else if ( sa && !sb) return false;
        // doc pr/33 F5: unconditional reachability counter for the
        // same-index fallback (0 on a manifest = never reached there);
        // knob-on orders by the stable per-run shower id (PRShower.cxx:45)
        // instead of heap addresses.  House-rule fix, prototype n/a.
        g_pr33_audit.f5_fallback_hits++;
        if (m_shower_less_id_tiebreak) return a->get_shower_id() < b->get_shower_id();
        return a.get() < b.get(); // same-index fallback: stable within a run
    };

    // Build good_showers set with stable ordering
    std::set<ShowerPtr, decltype(shower_less)> good_showers(shower_less);
    {
        auto it = map_vertex_to_shower.find(main_vertex);
        if (it != map_vertex_to_shower.end()) {
            for (auto shower : it->second) {
                if (pi0_showers.find(shower) != pi0_showers.end()) return;

                auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                if (conn_type == 1) {
                    good_showers.insert(shower);
                }
            }
        }

        if (good_showers.size() > 1) {
            ShowerPtr max_shower = nullptr;
            double max_energy = 0;
            for (auto shower : good_showers) { // iterates in shower_less order → deterministic tie-break
                double energy = shower->get_kine_charge();
                if (energy > max_energy) {
                    max_energy = energy;
                    max_shower = shower;
                }
            }
            good_showers.clear();
            if (max_shower) good_showers.insert(max_shower);
        }
    }
    
    // Check if we have exactly 2 segments at main vertex
    if (main_vertex_segs.size() == 2) {
        bool flag_return = true;
        int num_showers = 0;

        for (auto sg : main_vertex_segs) {
            if (map_segment_in_shower.find(sg) == map_segment_in_shower.end()) {
                double sg_length = segment_track_length(sg);
                if (sg_length < 1.2 * units::cm && (sg->dirsign() == 0 || seg_dir_weak(sg))) {
                    flag_return = false;
                }
            } else {
                num_showers++;
            }
        }

        if (flag_return && static_cast<size_t>(num_showers) == main_vertex_segs.size()) {
            flag_return = false;
        }
        if (flag_return) return;
    }
    
    // Build map of showers to rays (lines), ordered by shower_less for deterministic iteration
    std::map<ShowerPtr, WireCell::Ray, decltype(shower_less)> map_shower_ray(shower_less);
    WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
    
    // Add showers from main vertex
    auto it_main = map_vertex_to_shower.find(main_vertex);
    if (it_main != map_vertex_to_shower.end()) {
        for (auto shower : it_main->second) {
            // doc pr/33 F2: prototype reads the start segment's PDG here
            // (NeutrinoID_shower_clustering.h:497, exact ==13).
            {
                int type_startseg = 0;
                if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                    type_startseg = shower->start_segment()->particle_info()->pdg();
                }
                const bool skip_legacy = shower->get_particle_type() == 13;
                const bool skip_proto  = type_startseg == 13;
                g_pr33_audit.f2_calls[4]++;
                if (skip_legacy != skip_proto) g_pr33_audit.f2_disagree[4]++;
                if (m_shower_pdg_from_start_segment ? skip_proto : skip_legacy) continue;
            }
            if (shower->get_total_length() < 3 * units::cm) continue;
            if (pi0_showers.find(shower) != pi0_showers.end()) continue;

            WireCell::Point test_p = shower->get_start_point();
            WireCell::Vector dir = shower_cal_dir_3vector(*shower, test_p, 15 * units::cm);
            WireCell::Point p2(test_p.x() + dir.x(), test_p.y() + dir.y(), test_p.z() + dir.z());
            map_shower_ray[shower] = WireCell::Ray(test_p, p2);
        }
    }
    
    // Add showers from other vertices
    for (auto& [vtx, shower_set] : map_vertex_to_shower) {
        if (vtx == main_vertex) continue;
        
        for (auto shower : shower_set) {
            // doc pr/33 F2: prototype reads the start segment's PDG here
            // (NeutrinoID_shower_clustering.h:511, exact ==13).
            {
                int type_startseg = 0;
                if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                    type_startseg = shower->start_segment()->particle_info()->pdg();
                }
                const bool skip_legacy = shower->get_particle_type() == 13;
                const bool skip_proto  = type_startseg == 13;
                g_pr33_audit.f2_calls[5]++;
                if (skip_legacy != skip_proto) g_pr33_audit.f2_disagree[5]++;
                if (m_shower_pdg_from_start_segment ? skip_proto : skip_legacy) continue;
            }
            if (shower->get_total_length() < 3 * units::cm) continue;
            if (pi0_showers.find(shower) != pi0_showers.end()) continue;

            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            if (conn_type != 3) continue;
            
            {
                bool is_shower_seg = shower->start_segment()->flags_any(SegmentFlags::kShowerTrajectory) ||
                                     shower->start_segment()->flags_any(SegmentFlags::kShowerTopology);
                if (!is_shower_seg) {
                    // Also accept segments identified as electron by PID (prototype's get_flag_shower_dQdx)
                    int seg_pdg = 0;
                    if (shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info())
                        seg_pdg = shower->start_segment()->particle_info()->pdg();
                    if (std::abs(seg_pdg) != 11) continue;
                }
            }
            
            auto [closest_dis, test_p] = shower_get_closest_point(*shower, main_vtx_pt);
            WireCell::Vector dir = shower_cal_dir_3vector(*shower, test_p, 15 * units::cm);
            WireCell::Point p2(test_p.x() + dir.x(), test_p.y() + dir.y(), test_p.z() + dir.z());
            map_shower_ray[shower] = WireCell::Ray(test_p, p2);
        }
    }
    
    if (map_shower_ray.size() > 1) {
        // Calculate pi0 masses for shower pairs.
        // Use shower_less-based comparator on pair key to guarantee deterministic iteration order.
        auto pair_less = [&shower_less](const std::pair<ShowerPtr,ShowerPtr>& a,
                                        const std::pair<ShowerPtr,ShowerPtr>& b) {
            if (shower_less(a.first, b.first)) return true;
            if (shower_less(b.first, a.first)) return false;
            return shower_less(a.second, b.second);
        };
        std::map<std::pair<ShowerPtr, ShowerPtr>, std::pair<double, WireCell::Point>,
                 decltype(pair_less)> map_shower_pair_mass_point(pair_less);

        for (auto it = map_shower_ray.begin(); it != map_shower_ray.end(); ++it) {
            ShowerPtr shower_1 = it->first;
            WireCell::Ray ray1 = it->second;
            double length_1 = shower_1->get_total_length();

            // Start from next(it) — shower_1 == shower_2 is impossible, no skip needed
            for (auto it1 = std::next(it); it1 != map_shower_ray.end(); ++it1) {
                ShowerPtr shower_2 = it1->first;

                WireCell::Ray ray2 = it1->second;
                if (ray1.first == ray2.first) continue;
                
                double length_2 = shower_2->get_total_length();
                
                WireCell::Vector dir1 = ray_vector(ray1);
                WireCell::Vector dir2 = ray_vector(ray2);
                double angle_between = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0));
                if (angle_between == 0) continue;
                
                auto [p1_closest, p2_closest] = ray_closest_points(ray1, ray2);
                WireCell::Point center((p1_closest.x() + p2_closest.x()) / 2.0,
                                      (p1_closest.y() + p2_closest.y()) / 2.0,
                                      (p1_closest.z() + p2_closest.z()) / 2.0);
                
                if (length_1 > 15 * units::cm && length_2 > 15 * units::cm) {
                    WireCell::Vector dir_to_shower1(ray1.first.x() - center.x(),
                                                   ray1.first.y() - center.y(),
                                                   ray1.first.z() - center.z());
                    WireCell::Vector dir_to_shower2(ray2.first.x() - center.x(),
                                                   ray2.first.y() - center.y(),
                                                   ray2.first.z() - center.z());
                    
                    if (dir_to_shower1.magnitude() < 3 * units::cm) dir_to_shower1 = dir1;
                    if (dir_to_shower2.magnitude() < 3 * units::cm) dir_to_shower2 = dir2;
                    
                    double angle1 = std::acos(std::clamp(dir_to_shower1.dot(dir1) / (dir_to_shower1.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle2 = std::acos(std::clamp(dir_to_shower2.dot(dir2) / (dir_to_shower2.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    if (angle1 > 25 || angle2 > 25) continue;
                    
                    double angle = std::acos(std::clamp(dir_to_shower1.dot(dir_to_shower2) / (dir_to_shower1.magnitude() * dir_to_shower2.magnitude()), -1.0, 1.0));
                    double mass_pio = std::sqrt(4 * shower_1->get_kine_charge() * shower_2->get_kine_charge() * 
                                               std::pow(std::sin(angle / 2.0), 2));
                    map_shower_pair_mass_point[std::make_pair(shower_1, shower_2)] = std::make_pair(mass_pio, center);
                    
                } else if (length_1 > 15 * units::cm || length_2 > 15 * units::cm) {
                    WireCell::Vector dir_to_c1, dir_to_c2;

                    if (length_1 > length_2) {
                        // center already holds (p1_closest + p2_closest)/2 from above
                        auto [dis2, test_p] = shower_get_closest_point(*shower_2, center);
                        WireCell::Vector dir3 = shower_cal_dir_3vector(*shower_2, test_p, 15 * units::cm);
                        WireCell::Point p3(test_p.x() + dir3.x(), test_p.y() + dir3.y(), test_p.z() + dir3.z());
                        WireCell::Ray ray3(test_p, p3);
                        
                        auto [new_p1, new_p2] = ray_closest_points(ray1, ray3);
                        center = new_p1;
                        
                        dir_to_c1 = WireCell::Vector(ray1.first.x() - center.x(),
                                                     ray1.first.y() - center.y(),
                                                     ray1.first.z() - center.z());
                        dir_to_c2 = WireCell::Vector(test_p.x() - center.x(),
                                                     test_p.y() - center.y(),
                                                     test_p.z() - center.z());
                        
                        double angle1 = std::acos(std::clamp(dir_to_c1.dot(dir1) / (dir_to_c1.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double angle2 = std::acos(std::clamp(dir_to_c2.dot(dir3) / (dir_to_c2.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        if (angle1 > 25 || angle2 > 25) continue;
                        
                    } else {
                        // center already holds (p1_closest + p2_closest)/2 from above
                        auto [dis1, test_p] = shower_get_closest_point(*shower_1, center);
                        WireCell::Vector dir3 = shower_cal_dir_3vector(*shower_1, test_p, 15 * units::cm);
                        WireCell::Point p3(test_p.x() + dir3.x(), test_p.y() + dir3.y(), test_p.z() + dir3.z());
                        WireCell::Ray ray3(test_p, p3);
                        
                        auto [new_p1, new_p2] = ray_closest_points(ray3, ray2);
                        center = new_p2;
                        
                        dir_to_c2 = WireCell::Vector(ray2.first.x() - center.x(),
                                                     ray2.first.y() - center.y(),
                                                     ray2.first.z() - center.z());
                        dir_to_c1 = WireCell::Vector(test_p.x() - center.x(),
                                                     test_p.y() - center.y(),
                                                     test_p.z() - center.z());
                        
                        double angle1 = std::acos(std::clamp(dir_to_c1.dot(dir3) / (dir_to_c1.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double angle2 = std::acos(std::clamp(dir_to_c2.dot(dir2) / (dir_to_c2.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        if (angle1 > 25 || angle2 > 25) continue;
                    }
                    
                    double angle = std::acos(std::clamp(dir_to_c1.dot(dir_to_c2) / (dir_to_c1.magnitude() * dir_to_c2.magnitude()), -1.0, 1.0));
                    double mass_pio = std::sqrt(4 * shower_1->get_kine_charge() * shower_2->get_kine_charge() * 
                                               std::pow(std::sin(angle / 2.0), 2));
                    map_shower_pair_mass_point[std::make_pair(shower_1, shower_2)] = std::make_pair(mass_pio, center);
                    
                } else {
                    break; // both showers short: exit inner loop (prototype line 614)
                }
            }
        }
        
        // Fill BDT features for the highest-energy pi0 candidate (without-vertex case).
        // Reads pio_kine.energy_1 + energy_2 set by id_pi0_with_vertex as the starting threshold.
        {
            double max_energy = pio_kine.energy_1 + pio_kine.energy_2;
            for (auto& [shower_pair, mass_point] : map_shower_pair_mass_point) {
                ShowerPtr sh1 = shower_pair.first;
                ShowerPtr sh2 = shower_pair.second;
                double energy_1 = sh1->get_kine_charge();
                double energy_2 = sh2->get_kine_charge();
                if (energy_1 + energy_2 < max_energy) continue;
                max_energy = energy_1 + energy_2;

                WireCell::Point vtx_pt = mass_point.second;
                WireCell::Point sp1 = sh1->get_start_point();
                WireCell::Point sp2 = sh2->get_start_point();
                WireCell::Vector dir1(sp1.x() - vtx_pt.x(), sp1.y() - vtx_pt.y(), sp1.z() - vtx_pt.z());
                WireCell::Vector dir2(sp2.x() - vtx_pt.x(), sp2.y() - vtx_pt.y(), sp2.z() - vtx_pt.z());

                pio_kine.flag     = 2;
                pio_kine.mass     = mass_point.first;
                pio_kine.vtx_dis  = (vtx_pt - main_vtx_pt).magnitude();
                pio_kine.energy_1 = energy_1;
                pio_kine.energy_2 = energy_2;
                pio_kine.dis_1    = dir1.magnitude();
                pio_kine.dis_2    = dir2.magnitude();
                double mag1 = dir1.magnitude(), mag2 = dir2.magnitude();
                pio_kine.theta_1  = (mag1 > 0)
                    ? std::acos(std::clamp(dir1.z() / mag1, -1.0, 1.0))
                    : 0.0;
                pio_kine.theta_2  = (mag2 > 0)
                    ? std::acos(std::clamp(dir2.z() / mag2, -1.0, 1.0))
                    : 0.0;
                pio_kine.phi_1    = std::atan2(dir1.y(), dir1.x());
                pio_kine.phi_2    = std::atan2(dir2.y(), dir2.x());
                pio_kine.angle = (mag1 > 0 && mag2 > 0)
                    ? std::acos(std::clamp(dir1.dot(dir2) / (mag1 * mag2), -1.0, 1.0))
                    : 0.0;
            }
        }

        // Find best pi0 candidate
        double mass_diff = 1e9;
        double mass_save = 0;
        ShowerPtr shower_1 = nullptr;
        ShowerPtr shower_2 = nullptr;
        double mass_offset = 10 * units::MeV;
        WireCell::Point vtx_point;

        for (auto& [shower_pair, mass_point] : map_shower_pair_mass_point) {
            if (std::abs(mass_point.first - 135 * units::MeV + mass_offset) < mass_diff) {
                if (good_showers.find(shower_pair.first) == good_showers.end() && 
                    good_showers.find(shower_pair.second) == good_showers.end()) continue;
                
                shower_1 = shower_pair.first;
                shower_2 = shower_pair.second;
                mass_diff = std::abs(mass_point.first - 135 * units::MeV + mass_offset);
                mass_save = mass_point.first;
                vtx_point = mass_point.second;
            }
        }
        
        // If found good pi0, update everything
        if (mass_diff < 60 * units::MeV && shower_1 && shower_2) {
            pi0_showers.insert(shower_1);
            pi0_showers.insert(shower_2);
            
            int pio_id = acc_segment_id;
            acc_segment_id++;
            g_pr33_audit.f3_pi0_without_vertex++;

            map_shower_pio_id[shower_1] = pio_id;
            map_shower_pio_id[shower_2] = pio_id;
            map_pio_id_mass[pio_id] = std::make_pair(mass_save, 2);
            map_pio_id_showers[pio_id].push_back(shower_1);
            map_pio_id_showers[pio_id].push_back(shower_2);
            
            // Update main vertex position (hack) - set to reconstructed pi0 decay point
            main_vertex->fit().point = vtx_point;
            main_vertex->fit().dQ = 0;
            
            // Add other segments from main_vertex to showers
            auto [start_vtx_1, conn_type_1] = shower_1->get_start_vertex_and_type();
            if (start_vtx_1 == main_vertex && conn_type_1 == 1) {
                for (auto sg : main_vertex_segs) {
                    if (sg == shower_1->start_segment()) continue;
                    pr93_probe_absorb_direct("pi0_shower1", shower_1, sg);
                    shower_1->add_segment(sg, true);
                }
            }

            auto [start_vtx_2, conn_type_2] = shower_2->get_start_vertex_and_type();
            if (start_vtx_2 == main_vertex && conn_type_2 == 1) {
                for (auto sg : main_vertex_segs) {
                    if (sg == shower_2->start_segment()) continue;
                    pr93_probe_absorb_direct("pi0_shower2", shower_2, sg);
                    shower_2->add_segment(sg, true);
                }
            }
            
            shower_1->set_start_vertex(main_vertex, 2);
            shower_1->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
            
            shower_2->set_start_vertex(main_vertex, 2);
            shower_2->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
            
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                              map_vertex_to_shower, used_shower_clusters);
            
            SPDLOG_LOGGER_TRACE(s_log, "examine_showers: Pi0 (displaced vertex) found with mass: {} MeV with {} MeV + {} MeV",
                mass_save / units::MeV, shower_1->get_kine_charge() / units::MeV, shower_2->get_kine_charge() / units::MeV);
        }
    }
}


void PatternAlgorithms::shower_clustering_with_nv(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    // Diagnostic: print dirsign for all main-cluster segments at entry,
    // to verify examine_direction() ran correctly before this call.
    if (main_cluster) {
        for (auto e : ordered_edges(graph)) {
            SegmentPtr seg = graph[e].segment;
            if (!seg || seg->cluster() != main_cluster) continue;
            if (seg->fits().size() >= 2) {  // only segments with real fit data
                SPDLOG_LOGGER_TRACE(s_log,
                    "shower_clustering_with_nv entry: seg id={} nfits={} nwcpts={} dirsign={} dir_weak={}"
                    " shower_topo={} shower_traj={} pdg={}",
                    seg->id(), seg->fits().size(), seg->wcpts().size(),
                    seg->dirsign(), seg_dir_weak(seg) ? 1 : 0,
                    seg->flags_any(SegmentFlags::kShowerTopology) ? 1 : 0,
                    seg->flags_any(SegmentFlags::kShowerTrajectory) ? 1 : 0,
                    seg->has_particle_info() ? seg->particle_info()->pdg() : 0);
            }
        }
    }

    // doc sbnd_xin/docs/pr/65 round 3: main-cluster segments unreachable from
    // main_vertex, computed ONCE here -- PR-graph topology is frozen for the
    // whole span of this function (no graph mutators in this file; the last
    // mutation, main_vertex_graph_audit, runs before this call).  Recomputed
    // unconditionally so no state leaks between calls; empty when the knob is
    // off, which keeps every relaxed guard below byte-identical to legacy.
    // Such segments exist only when other_seg_keep_isolated (doc pr/54) kept
    // a residual segment as a disconnected component of the main cluster's
    // graph; the prototype cannot reach this state (attach-or-discard).
    // doc pr/40 round 9: the B2 bridge (nv_bridge_track, from_vertices) DOES
    // add graph edges after this point, but only main-cluster-vertex ->
    // other-cluster edges -- it cannot un-island a main-cluster segment, so
    // the unreachable set computed below stays valid for the whole call.
    m_absorb_unreachable_main_segs.clear();
    m_sfv_declined_anchors.clear();
    m_nv_bridge_shield_segs.clear();
    m_nv_bridge_cluster_ids.clear();
    track_fitter.clear_bridged_cluster_ids();
    m_sccc_bridged_cluster_ids.clear();
    if (m_shower_absorb_unreachable_main && main_vertex && main_cluster) {
        for (const auto& seg : unreachable_segments(graph, main_vertex)) {
            if (seg && seg->cluster() == main_cluster) {
                m_absorb_unreachable_main_segs.insert(seg);
            }
        }
        if (!m_absorb_unreachable_main_segs.empty()) {
            std::string ids;
            for (const auto& seg : m_absorb_unreachable_main_segs) {
                if (!ids.empty()) ids += ",";
                ids += std::to_string(seg->id());
            }
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr65 absorb_unreachable_main: {} graph-unreachable main-cluster segment(s) offered to shower absorbers: [{}]",
                m_absorb_unreachable_main_segs.size(), ids);
        }
    }

    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = t_total;
    MS t_in_main_cluster{0};
    MS t_connecting_to_main_vertex{0};
    MS t_from_main_cluster{0};
    MS t_from_vertices{0};
    MS t_calc_kine_1{0};
    MS t_examine_merge{0};
    MS t_in_other_clusters{0};
    MS t_calc_kine_2{0};
    MS t_examine_showers{0};
    MS t_id_pi0_with_vertex{0};
    MS t_id_pi0_without_vertex{0};
    
    // // Debug helper: check if cluster 933 has been added to used_shower_clusters
    // auto check_used_shower_cluster_933 = [&used_shower_clusters](const char* label) {
    //     bool found = false;
    //     for (auto* c : used_shower_clusters) {
    //         if (c->get_cluster_id() == 933) { found = true; break; }
    //     }
    //     std::cout << "[shower_clustering_with_nv] after " << label
    //               << ": cluster 933 in used_shower_clusters=" << (found ? "YES" : "NO")
    //               << " (size=" << used_shower_clusters.size() << ")" << std::endl;
    // };

    // Connect to the main cluster
    shower_clustering_with_nv_in_main_cluster(graph, main_vertex, showers,
                                              map_vertex_in_shower, map_segment_in_shower,
                                              map_vertex_to_shower, used_shower_clusters,
                                              vertices_in_long_muon, segments_in_long_muon,
                                              particle_data, recomb_model);
    t_in_main_cluster = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("shower_clustering_with_nv_in_main_cluster");
    
    // Examine things connecting to the main vertex
    shower_clustering_connecting_to_main_vertex(graph, main_vertex, showers,
                                                map_vertex_in_shower, map_segment_in_shower,
                                                map_vertex_to_shower, used_shower_clusters);
    t_connecting_to_main_vertex = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("shower_clustering_connecting_to_main_vertex");
    
    // Shower clustering from main cluster
    shower_clustering_with_nv_from_main_cluster(graph, main_vertex, main_cluster, showers,
                                                map_vertex_in_shower, map_segment_in_shower,
                                                map_vertex_to_shower, used_shower_clusters);
    t_from_main_cluster = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("shower_clustering_with_nv_from_main_cluster");
    
    // doc sbnd_xin/docs/pr/93 round 4 (sccc_bridge_body): replay the bridge
    // requests recorded by demote_cross_cluster_straight_stems.  The pass
    // itself runs BEFORE this function (between examine_direction and the
    // shower seeding), but the three bridge-state containers are cleared at
    // this function's entry -- so the pre-pass only RECORDS
    // {stem tip vertex, continuation cluster, continuation vertex} and the
    // actual nv_bridge_connect (real graph edge + shields + bridged-cluster
    // ids) happens HERE, at the same point in the pass sequence where
    // nv_bridge_track itself fires (inside from_vertices): passes 1-3
    // (in_main_cluster / connecting_to_main_vertex / from_main_cluster) see
    // the PRE-bridge graph exactly as legacy.  Replaying any earlier
    // exposes the bridged cluster to pass 2, which re-forms a conn-1
    // shower on the demoted stem (measured on 18264-137238: stem + entry
    // stubs re-labelled "e- 116 MeV" and the 81cm muon body orphaned).
    // The pr/65 unreachable-set invariant also holds (a bridge edge cannot
    // un-island a main-cluster segment).  Knob off => empty vector =>
    // byte-identical.
    if (!m_sccc_bridge_requests.empty() && main_cluster) {
        for (const auto& req : m_sccc_bridge_requests) {   // recorded order (deterministic)
            if (!req.cluster || !req.main_vtx || !req.far_vtx) continue;
            // Deterministic cluster_segs: the continuation cluster's own
            // segments in stable graph-edge-index order.
            std::vector<std::pair<size_t, SegmentPtr>> csegs;
            for (auto [eit, eend] = boost::edges(graph); eit != eend; ++eit) {
                SegmentPtr s = graph[*eit].segment;
                if (!s || s->cluster() != req.cluster) continue;
                csegs.emplace_back(graph[*eit].index, s);
            }
            std::sort(csegs.begin(), csegs.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            std::vector<SegmentPtr> cluster_segs;
            cluster_segs.reserve(csegs.size());
            for (const auto& [i, s] : csegs) cluster_segs.push_back(s);
            SegmentPtr bridge = nv_bridge_connect(graph, main_cluster, req.cluster,
                                                  req.main_vtx, req.far_vtx,
                                                  cluster_segs, track_fitter, dv);
            if (bridge) m_sccc_bridged_cluster_ids.insert(req.cluster->get_cluster_id());
            SPDLOG_LOGGER_DEBUG(s_log,
                "sccc bridge: cluster {} -> main {} main_vtx_idx={} far_vtx_idx={} bridge={}",
                req.cluster->get_cluster_id(), main_cluster->get_cluster_id(),
                req.main_vtx->get_graph_index(), req.far_vtx->get_graph_index(),
                bridge ? "OK" : "FAILED");
        }
        m_sccc_bridge_requests.clear();
    }

    // Shower clustering from vertices
    shower_clustering_with_nv_from_vertices(graph, main_vertex, main_cluster, other_clusters, showers,
                                           map_vertex_in_shower, map_segment_in_shower,
                                           map_vertex_to_shower, used_shower_clusters,
                                           vertices_in_long_muon, segments_in_long_muon,
                                           track_fitter, dv, particle_data, recomb_model);
    t_from_vertices = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("shower_clustering_with_nv_from_vertices");

    // Collect 2D charge maps ONCE here — after track fitting (shower_clustering_with_nv_from_vertices)
    // has populated the underlying charge data, but before any kinematics calculations.
    // All cal_kine_charge call sites below reuse m_charge_* to avoid repeated O(N_hits) collection.
    collect_charge_maps(track_fitter);
    SPDLOG_LOGGER_TRACE(s_log,
        "shower_clustering_with_nv: collected charge maps U={} V={} W={} hits, {} shower(s) before calc_kine_1",
        m_charge_2d_u.size(), m_charge_2d_v.size(), m_charge_2d_w.size(), showers.size());

    // Calculate shower kinematics
    calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                graph, track_fitter, dv, particle_data, recomb_model);
    t_calc_kine_1 = MS(Clock::now() - t0); t0 = Clock::now();
    
    // Examine and merge showers
    examine_merge_showers(showers, main_vertex, map_vertex_in_shower, map_segment_in_shower,
                         map_vertex_to_shower, used_shower_clusters,
                         vertices_in_long_muon, segments_in_long_muon,
                         graph, track_fitter, dv, particle_data, recomb_model);
    t_examine_merge = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("examine_merge_showers");
    
    // Check remaining clusters
    shower_clustering_in_other_clusters(graph, main_vertex, showers, main_cluster, other_clusters,
                                       map_cluster_main_vertices, map_vertex_in_shower,
                                       map_segment_in_shower, map_vertex_to_shower,
                                       used_shower_clusters, track_fitter, dv,
                                       particle_data, recomb_model, true);
    t_in_other_clusters = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("shower_clustering_in_other_clusters");

    // doc sbnd_xin/docs/pr/74 round 2 K4 -- absorb the walked-past track stem
    // between the main vertex and each substantial EM shower's attach vertex.
    // After in_other_clusters (so conn-3/4 showers exist) and before the
    // second kinematics pass (so absorbed charge is counted).  false = no
    // call = byte-identical.
    if (m_shower_stem_backfill) {
        stem_backfill(graph, main_vertex, showers, map_vertex_in_shower,
                      map_segment_in_shower, map_vertex_to_shower,
                      used_shower_clusters, segments_in_long_muon,
                      particle_data, recomb_model);
        t0 = Clock::now();
    }

    // Calculate shower kinematics again
    SPDLOG_LOGGER_TRACE(s_log,
        "shower_clustering_with_nv: {} shower(s) before calc_kine_2", showers.size());
    calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                graph, track_fitter, dv, particle_data, recomb_model);
    t_calc_kine_2 = MS(Clock::now() - t0); t0 = Clock::now();
    
    // Examine shower trunk and add to shower
    examine_showers(graph, main_vertex, showers, main_cluster, other_clusters,
                   map_cluster_main_vertices, map_vertex_in_shower, map_segment_in_shower,
                   map_vertex_to_shower, used_shower_clusters,
                   track_fitter, dv, particle_data, recomb_model);
    t_examine_showers = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("examine_showers");
    
    // doc sbnd_xin/docs/pr/84 round 3: collapse showers that share a start
    // segment.  Here and not earlier on purpose -- this is after EVERY pass
    // that can create or retarget a shower (examine_showers is the last), and
    // before the pi0 finders, so the pi0 pairing, fill_kine_tree and the Bee
    // PF writer all see the same de-duplicated set.  Knob off => no pass.
    if (m_shower_dedup_start_seg) {
        const int n_absorbed = merge_showers_sharing_start_segment(showers);
        if (n_absorbed > 0) {
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                               map_vertex_to_shower, used_shower_clusters);
            calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                        graph, track_fitter, dv, particle_data, recomb_model);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr84 shower_dedup: absorbed {} shower(s); {} remain", n_absorbed, showers.size());
        }
    }

    // doc sbnd_xin/docs/pr/93 round 4 (shower_detach_track_stem, SBND
    // 18255-348471 + 18255-292643): a TRACK-HEADED shower -- one whose final
    // start segment is a non-shower-flagged, non-zero, non-+-11-pdg track --
    // is really "track + daughter EM shower" glued into one object.  Peel the
    // main-cluster track prefix back into the track pool and re-root the
    // shower at the prefix's far vertex (conn 2), so the PF tree renders
    // "proton/pi+ track -> pseudo-gamma -> EM shower", kine costs the track
    // by range (not shower charge + rest mass), and the pi0 pairing sees a
    // real EM shower.  Placed after EVERY pass that creates/retargets/absorbs
    // showers and before the pi0 finders / kinematics consumers, so
    // everything downstream recomputes naturally.  This is a STRUCTURE fix
    // on the final shower set, deliberately NOT keyed on the pr/93 round-3
    // confident-score/50cm family: 292643's 22.9cm pi+ stem carries the
    // score-100 sentinel and sits under the floor.  Long-muon pseudo-showers
    // (cached type +-13, the discriminant every consumer already routes on)
    // are exempt.  C++ default false => no pass => byte-identical.
    if (m_shower_detach_track_stem) {
        int n_detached = 0;
        for (auto& shower : showers) {                       // IndexedShowerSet order
            VertexPtr sv = shower->start_vertex();
            SegmentPtr ss = shower->start_segment();
            if (!sv || !ss || !ss->descriptor_valid()) continue;
            if (std::abs(shower->get_particle_type()) == 13) continue;   // long-muon pseudo-shower
            if (segments_in_long_muon.count(ss)) continue;               // belt and braces
            if (!ss->has_particle_info() || !ss->particle_info()) continue;
            const int head_pdg = ss->particle_info()->pdg();
            if (head_pdg == 0 || std::abs(head_pdg) == 11) continue;
            if (ss->flags_any(SegmentFlags::kShowerTrajectory) ||
                ss->flags_any(SegmentFlags::kShowerTopology)) continue;

            // Walk the main-cluster track prefix from the start vertex:
            // follow the unique un-shower-flagged non-+-11 same-cluster
            // member continuation; stop where an EM member attaches or the
            // chain branches/ends.
            auto track_like = [&](SegmentPtr sg) -> bool {
                return sg && sg->descriptor_valid() && shower->has_edge(sg->get_descriptor())
                    && !sg->flags_any(SegmentFlags::kShowerTrajectory)
                    && !sg->flags_any(SegmentFlags::kShowerTopology)
                    && sg->has_particle_info() && sg->particle_info()
                    && sg->particle_info()->pdg() != 0
                    && std::abs(sg->particle_info()->pdg()) != 11
                    && sg->cluster() == ss->cluster();
            };
            std::vector<SegmentPtr> prefix;
            VertexPtr reroot = nullptr;
            VertexPtr cur_vtx = sv;
            SegmentPtr cur_seg = ss;
            bool walk_ok = true;
            for (int hop = 0; hop < 8; ++hop) {              // cycle guard
                if (!track_like(cur_seg)) break;
                prefix.push_back(cur_seg);
                VertexPtr far = find_other_vertex(graph, cur_seg, cur_vtx);
                if (!far || !far->descriptor_valid()) { walk_ok = false; break; }
                std::vector<SegmentPtr> cont;
                bool em_here = false;
                for (auto e : sorted_out_edges(far->get_descriptor(), graph)) {
                    SegmentPtr nx = graph[e].segment;
                    if (!nx || nx == cur_seg) continue;
                    if (!shower->has_edge(nx->get_descriptor())) continue;  // members only
                    if (track_like(nx)) cont.push_back(nx);
                    else em_here = true;
                }
                reroot = far;
                if (em_here || cont.size() != 1) break;
                cur_vtx = far;
                cur_seg = cont.front();
            }
            if (!walk_ok || prefix.empty() || !reroot) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr93 detach_track_stem: skip shower_id={} head_pdg={} (walk_ok={} prefix={})",
                    shower->get_shower_id(), head_pdg, walk_ok, prefix.size());
                continue;
            }
            const int peeled = shower->detach_track_prefix(prefix, reroot);
            if (!peeled) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr93 detach_track_stem: refuse shower_id={} head_pdg={} prefix={} nseg={}",
                    shower->get_shower_id(), head_pdg, prefix.size(), shower->get_num_segments());
                continue;
            }
            std::string peeled_ids;
            for (const auto& sg : prefix) {
                if (!peeled_ids.empty()) peeled_ids += ",";
                peeled_ids += std::to_string(sg->id());
            }
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr93 detach_track_stem: shower_id={} peel {} track seg(s) [{}] head_pdg={} "
                "reroot_vtx_idx={} conn=2 n_remain={}",
                shower->get_shower_id(), peeled, peeled_ids, head_pdg,
                graph[reroot->get_descriptor()].index, shower->get_num_segments());
            // Mirror the examine_showers merge-loop refresh: vote + kinematics
            // + charge over the REMAINING (EM) members only.
            shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex,
                                         m_shower_proton_daughter_pion, m_mip_dqdx_median,
                                         m_shower_vote_track_pid_counts, m_shower_accept_pid_guard,
                                         m_shower_pid_guard_min_len);
            shower->calculate_kinematics(particle_data, recomb_model,
                                         m_shower_endpoint_exclude_start_vertex,
                                         m_shower_endpoint_skip_orphan_vtx);
            shower->set_kine_charge(cal_kine_charge(shower, m_charge_2d_u, m_charge_2d_v,
                                                    m_charge_2d_w, m_map_apa_ch_plane_wires,
                                                    track_fitter, dv));
            shower->set_flag_kinematics(true);
            ++n_detached;
        }
        if (n_detached) {
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                               map_vertex_to_shower, used_shower_clusters);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr93 detach_track_stem: {} shower(s) split; maps rebuilt", n_detached);
        }
    }

    // Identify pi0 with vertex.
    // doc pr/33 F3: both finders get a reference to the same local copy, so
    // nothing propagates past this function either way (the caller's
    // variable is separately seeded by reference into ssm_tagger -- §10.10
    // amendment 1).  Knob-off restores the copy between the two finders, so
    // each seeds from the same base = the legacy by-value behavior.
    int pi0_acc = acc_segment_id;
    id_pi0_with_vertex(pi0_acc, pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass,
                      map_pio_id_saved_pair, pio_kine, graph, main_vertex, showers, main_cluster,
                      other_clusters, map_cluster_main_vertices, map_vertex_in_shower,
                      map_segment_in_shower, map_vertex_to_shower, used_shower_clusters,
                      track_fitter, dv, particle_data, recomb_model);
    t_id_pi0_with_vertex = MS(Clock::now() - t0); t0 = Clock::now();
    // check_used_shower_cluster_933("id_pi0_with_vertex");

    // Identify pi0 without vertex (displaced vertex)
    if (!m_pi0_id_shared_allocator) pi0_acc = acc_segment_id;
    id_pi0_without_vertex(pi0_acc, pi0_showers, map_shower_pio_id, map_pio_id_showers,
                         map_pio_id_mass, map_pio_id_saved_pair, pio_kine, graph, main_vertex, showers,
                         main_cluster, other_clusters, map_cluster_main_vertices,
                         map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower,
                         used_shower_clusters, segments_in_long_muon, track_fitter, dv, particle_data, recomb_model);
    t_id_pi0_without_vertex = MS(Clock::now() - t0);
    // check_used_shower_cluster_933("id_pi0_without_vertex");

    if (m_perf) {
        SPDLOG_LOGGER_TRACE(s_log,
            "shower_clustering_with_nv timing: "
            "in_main_cluster={:.3f}ms connecting_to_main_vertex={:.3f}ms "
            "from_main_cluster={:.3f}ms from_vertices={:.3f}ms "
            "calc_kine_1={:.3f}ms examine_merge={:.3f}ms "
            "in_other_clusters={:.3f}ms calc_kine_2={:.3f}ms "
            "examine_showers={:.3f}ms id_pi0_with_vertex={:.3f}ms "
            "id_pi0_without_vertex={:.3f}ms",
            t_in_main_cluster.count(), t_connecting_to_main_vertex.count(),
            t_from_main_cluster.count(), t_from_vertices.count(),
            t_calc_kine_1.count(), t_examine_merge.count(),
            t_in_other_clusters.count(), t_calc_kine_2.count(),
            t_examine_showers.count(), t_id_pi0_with_vertex.count(),
            t_id_pi0_without_vertex.count());
        SPDLOG_LOGGER_TRACE(s_log,
            "shower_clustering_with_nv timing: TOTAL={:.3f}ms",
            MS(Clock::now() - t_total).count());
    }

    // Debug summary: showers and pi0s
    {
        SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_with_nv: {} shower(s)", showers.size());
        int idx = 0;
        int pr91_idx = 0;
        for (auto& shower : showers) {
            if (!shower) continue;
            // Count unique clusters via map_segment_in_shower (already up-to-date)
            std::set<Facade::Cluster*> shower_clusters;
            for (auto& [seg, sh] : map_segment_in_shower) {
                if (sh == shower && seg && seg->cluster())
                    shower_clusters.insert(seg->cluster());
            }
            // Use data.start_point if set; fall back to start_vertex position.
            WireCell::Point sp = shower->get_start_point();
            // if (sp.x() == 0 && sp.y() == 0 && sp.z() == 0) {
            //     auto vtx = shower->start_vertex();
            //     if (vtx) sp = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            // }
            WireCell::Vector dir = shower->get_init_dir();
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            auto start_seg = shower->start_segment();
            int start_seg_dirsign = start_seg ? start_seg->dirsign() : -99;
            WireCell::Point vtx_pt(0, 0, 0);
            if (start_vtx) vtx_pt = start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point;
            SPDLOG_LOGGER_TRACE(s_log,
                "shower_clustering_with_nv:   shower[{}] pdg={} flag_shower={} conn={}"
                " nseg={} ncls={} kine_charge={:.1f}MeV kine_best={:.1f}MeV"
                " start=({:.1f},{:.1f},{:.1f})cm vtx=({:.1f},{:.1f},{:.1f})cm"
                " dir=({:.3f},{:.3f},{:.3f}) seg_dirsign={}",
                idx++,
                shower->get_particle_type(), shower->get_flag_shower() ? 1 : 0, conn_type,
                shower->get_num_segments(), shower_clusters.size(),
                shower->get_kine_charge() / units::MeV, shower->get_kine_best() / units::MeV,
                sp.x() / units::cm, sp.y() / units::cm, sp.z() / units::cm,
                vtx_pt.x() / units::cm, vtx_pt.y() / units::cm, vtx_pt.z() / units::cm,
                dir.x(), dir.y(), dir.z(), start_seg_dirsign);

            // doc pr/91 round 1 -- full membership + orphan-vertex dump, env-gated.
            // Own counter: `idx++` above lives inside SPDLOG_LOGGER_TRACE's
            // argument list and is compiled out with the macro.
            pr91_probe_shower_content(shower, graph, pr91_idx++);
        }

        SPDLOG_LOGGER_TRACE(s_log, "shower_clustering_with_nv: {} pi0(s)", map_pio_id_showers.size());
        for (auto& [pio_id, pi0_shower_vec] : map_pio_id_showers) {
            auto mass_it = map_pio_id_mass.find(pio_id);
            double mass = mass_it != map_pio_id_mass.end() ? mass_it->second.first : 0.0;
            int    flag = mass_it != map_pio_id_mass.end() ? mass_it->second.second : 0;
            SPDLOG_LOGGER_TRACE(s_log,
                "shower_clustering_with_nv:   pi0[id={}] flag={} mass={:.1f}MeV",
                pio_id, flag, mass / units::MeV);
            for (size_t i = 0; i < pi0_shower_vec.size(); ++i) {
                auto& s = pi0_shower_vec[i];
                if (!s) continue;
                WireCell::Point  sp  = s->get_start_point();
                WireCell::Vector dir = s->get_init_dir();
                SPDLOG_LOGGER_TRACE(s_log,
                    "shower_clustering_with_nv:     pi0[id={}] shower[{}] pdg={}"
                    " nseg={} kine_charge={:.1f}MeV"
                    " start=({:.1f},{:.1f},{:.1f})cm dir=({:.3f},{:.3f},{:.3f})",
                    pio_id, i,
                    s->get_particle_type(), s->get_num_segments(),
                    s->get_kine_charge() / units::MeV,
                    sp.x() / units::cm, sp.y() / units::cm, sp.z() / units::cm,
                    dir.x(), dir.y(), dir.z());
            }
        }
    }
}