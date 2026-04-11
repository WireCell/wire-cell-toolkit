// NeutrinoTaggerNuMu.cxx
//
// Ports of (NeutrinoID_numu_tagger.h):
//   WCPPID::NeutrinoID::numu_tagger()
//   WCPPID::NeutrinoID::count_daughters(ProtoSegment*)
//   WCPPID::NeutrinoID::count_daughters(WCShower*)
//
// Namespace/class: WireCell::Clus::PR::PatternAlgorithms
//
// Translation conventions (see neutrino_id_function_map.md):
//   map_vertex_segments[main_vertex]           → boost::out_edges(main_vertex->get_descriptor(), graph)
//   map_segment_vertices (all segs)            → boost::edges(graph) filtered by cluster id
//   sg->get_particle_type()                    → sg->particle_info()->pdg()
//   sg->get_length()                           → segment_track_length(sg)
//   sg->get_direct_length()                    → segment_track_direct_length(sg)
//   sg->get_medium_dQ_dx()                     → segment_median_dQ_dx(sg)
//   sg->get_flag_shower()                      → seg_is_shower(sg)
//   sg->get_flag_shower_topology()             → sg->flags_any(SegmentFlags::kShowerTopology)
//   sg->get_flag_avoid_muon_check()            → sg->flags_any(SegmentFlags::kAvoidMuonCheck)
//   find_vertices(sg)                          → find_vertices(graph, sg)
//   find_cont_muon_segment(sg, vtx)            → find_cont_muon_segment(graph, sg, vtx)
//   shower->get_start_segment()                → shower->start_segment()
//   shower->get_total_track_length()           → shower->get_total_track_length()
//   shower->get_total_length()                 → shower->get_total_length()
//   neutrino_type |= 1UL<<2 / 1UL<<3          → not tracked; caller uses returned flag_numu_cc
//
// NOT ported here (require TMVA, no toolkit dependency yet):
//   cal_numu_bdts(), cal_numu_bdts_xgboost()
//   cal_cosmict_2_4_bdt(), ..._3_5_bdt(), ..._6_bdt(), ..._7_bdt(), ..._8_bdt(), ..._10_bdt()
//   cal_numu_1_bdt(), cal_numu_2_bdt(), cal_numu_3_bdt()

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"
#include <cmath>
#include <set>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;
using namespace WireCell;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// ---------------------------------------------------------------------------
// Local helpers (mirrors NeutrinoTaggerCosmic.cxx helpers)
// ---------------------------------------------------------------------------

static inline Point vtx_fit_pt(VertexPtr v) {
    if (!v) return Point{};
    return v->fit().valid() ? v->fit().point : v->wcpt().point;
}

static inline bool seg_is_shower(SegmentPtr seg) {
    return seg->flags_any(SegmentFlags::kShowerTrajectory) ||
           seg->flags_any(SegmentFlags::kShowerTopology)   ||
           (seg->has_particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
}

// Returns true if sg is directly attached to main_vertex in graph.
// Prototype: map_vertex_segments[main_vertex].find(sg) != end()
static bool seg_at_main_vertex(SegmentPtr sg, VertexPtr main_vertex, const Graph& graph) {
    if (!main_vertex || !main_vertex->descriptor_valid()) return false;
    for (auto [eit, end] = boost::out_edges(main_vertex->get_descriptor(), graph); eit != end; ++eit) {
        if (graph[*eit].segment == sg) return true;
    }
    return false;
}

// ===========================================================================
// count_daughters (segment overload)
//
// Counts track daughters (n_daughter_tracks) and all daughters (n_daughter_all)
// at the far end of a muon segment, i.e. the end NOT connected to main_vertex.
//
// BFS is done from the vertex closer to main_vertex, going through sg
// (calculate_num_daughter_tracks counts what is reachable BEYOND sg from that vertex).
// The segment itself is subtracted from the counts (it's the muon, not a daughter).
//
// Prototype: WCPPID::NeutrinoID::count_daughters(ProtoSegment*) in NeutrinoID_numu_tagger.h.
// ===========================================================================
std::pair<int, int> PatternAlgorithms::count_daughters(
    Graph& graph, SegmentPtr max_muon, VertexPtr main_vertex)
{
    if (!max_muon || !main_vertex) return {0, 0};

    auto [v1, v2] = find_vertices(graph, max_muon);
    if (!v1 || !v2) return {0, 0};

    Point mv_pt = vtx_fit_pt(main_vertex);
    double dis1  = ray_length(Ray{vtx_fit_pt(v1), mv_pt});
    double dis2  = ray_length(Ray{vtx_fit_pt(v2), mv_pt});

    // Use the vertex closer to main_vertex as the BFS start.
    // calculate_num_daughter_tracks(graph, vtx, sg, ...) counts everything reachable
    // through sg from vtx — i.e., segments at the FAR end of sg.
    VertexPtr close_vtx = (dis1 < dis2) ? v1 : v2;

    int n_daughter_tracks = calculate_num_daughter_tracks(graph, close_vtx, max_muon, false, 3*units::cm).first;
    int n_daughter_all    = calculate_num_daughter_tracks(graph, close_vtx, max_muon, true,  3*units::cm).first;

    // Subtract the muon segment itself: it is a track (unless it is a shower)
    if (!seg_is_shower(max_muon)) n_daughter_tracks -= 1;
    n_daughter_all -= 1;

    return {n_daughter_tracks, n_daughter_all};
}

// ===========================================================================
// count_daughters (shower overload)
//
// Same semantics as the segment overload, but for a long-muon shower chain.
// Finds the last segment of the chain (via get_last_segment_vertex_long_muon),
// then counts daughters at its far end.
//
// Prototype: WCPPID::NeutrinoID::count_daughters(WCShower*) in NeutrinoID_numu_tagger.h.
// segments_in_long_muon_plain: plain std::set (no comparator), as required by
//   get_last_segment_vertex_long_muon. Convert from IndexedSegmentSet at call site.
// ===========================================================================
std::pair<int, int> PatternAlgorithms::count_daughters(
    Graph& graph,
    ShowerPtr max_long_muon,
    VertexPtr main_vertex,
    std::set<SegmentPtr>& segments_in_long_muon_plain)
{
    if (!max_long_muon || !main_vertex) return {0, 0};

    auto [last_sg, other_vtx] = max_long_muon->get_last_segment_vertex_long_muon(segments_in_long_muon_plain);
    if (!last_sg) return {0, 0};

    auto [v1, v2] = find_vertices(graph, last_sg);
    if (!v1 || !v2) return {0, 0};

    Point mv_pt = vtx_fit_pt(main_vertex);
    double dis1  = ray_length(Ray{vtx_fit_pt(v1), mv_pt});
    double dis2  = ray_length(Ray{vtx_fit_pt(v2), mv_pt});
    VertexPtr close_vtx = (dis1 < dis2) ? v1 : v2;

    int n_daughter_tracks = calculate_num_daughter_tracks(graph, close_vtx, last_sg, false, 3*units::cm).first;
    int n_daughter_all    = calculate_num_daughter_tracks(graph, close_vtx, last_sg, true,  3*units::cm).first;

    if (!seg_is_shower(last_sg)) n_daughter_tracks -= 1;
    n_daughter_all -= 1;

    return {n_daughter_tracks, n_daughter_all};
}

// ===========================================================================
// numu_tagger
//
// Identifies numu CC events via three complementary checks:
//   Flag 1 (numu_cc_1): a muon-like segment directly attached to main_vertex
//   Flag 2 (numu_cc_2): a long-muon shower (> 18 cm track length) from main cluster
//   Flag 3 (numu_cc_3): any muon-like track in the main cluster (not at main_vertex),
//                       longer than 25 cm (or 30 cm if no other tracks)
//
// Returns {flag_long_muon, max_muon_length}.
//   flag_long_muon = true if (max_muon_length > 100 cm || max_length_all > 120 cm) && numu_cc
//
// Prototype: WCPPID::NeutrinoID::numu_tagger() in NeutrinoID_numu_tagger.h.
// Prototype also sets neutrino_type bits — not done here; caller uses the returned bool.
// ===========================================================================
std::pair<bool, double> PatternAlgorithms::numu_tagger(
    Graph&            graph,
    VertexPtr         main_vertex,
    IndexedShowerSet& showers,
    IndexedSegmentSet& segments_in_long_muon,
    Facade::Cluster*  main_cluster,
    TaggerInfo&       ti)
{
    bool flag_long_muon    = false;
    bool flag_numu_cc      = false;
    bool flag_numu_cc_1_save = false;
    bool flag_numu_cc_2_save = false;
    bool flag_numu_cc_3    = false;

    // Muon length cut: segment shorter than dis_cut is not useful as a muon.
    const double dis_cut = 5 * units::cm;

    double max_muon_length = 0;
    SegmentPtr max_muon    = nullptr;
    ShowerPtr  max_long_muon = nullptr;

    const int main_cl_id = main_cluster ? main_cluster->get_cluster_id() : -1;

    // Convert IndexedSegmentSet → plain set for get_last_segment_vertex_long_muon
    std::set<SegmentPtr> segments_in_long_muon_plain(
        segments_in_long_muon.begin(), segments_in_long_muon.end());

    // -----------------------------------------------------------------------
    // Flag 1: direct muon check — segments at main_vertex.
    // Prototype: NeutrinoID_numu_tagger.h lines 27-66.
    // -----------------------------------------------------------------------
    if (main_vertex && main_vertex->descriptor_valid()) {
        auto vd = main_vertex->get_descriptor();
        for (auto [eit, eit_end] = boost::out_edges(vd, graph); eit != eit_end; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;

            bool flag_numu_cc_1 = false;
            double length        = segment_track_length(sg);
            double direct_length = segment_track_direct_length(sg);
            double medium_dQ_dx  = segment_median_dQ_dx(sg);
            double dQ_dx_cut     = 0.8866 + 0.9533 * std::pow(18*units::cm / length, 0.4234);
            int    pdg           = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;

            auto [n_daughter_tracks, n_daughter_all] = count_daughters(graph, sg, main_vertex);

            if (std::abs(pdg) == 13 &&
                length > dis_cut &&
                medium_dQ_dx < dQ_dx_cut * 43e3/units::cm &&
                (length > 40*units::cm || (length <= 40*units::cm && direct_length > 0.925 * length)) &&
                !(n_daughter_tracks > 1 || n_daughter_all - n_daughter_tracks > 2)) {
                flag_numu_cc_1 = true;
                if (length > max_muon_length) {
                    max_muon_length = length;
                    max_muon        = sg;
                }
            }

            ti.numu_cc_flag_1.push_back(static_cast<float>(flag_numu_cc_1));
            ti.numu_cc_1_particle_type.push_back(static_cast<float>(pdg));
            ti.numu_cc_1_length.push_back(static_cast<float>(length / units::cm));
            ti.numu_cc_1_direct_length.push_back(static_cast<float>(direct_length / units::cm));
            ti.numu_cc_1_medium_dQ_dx.push_back(static_cast<float>(medium_dQ_dx / (43e3/units::cm)));
            ti.numu_cc_1_dQ_dx_cut.push_back(static_cast<float>(dQ_dx_cut));
            ti.numu_cc_1_n_daughter_tracks.push_back(static_cast<float>(n_daughter_tracks));
            ti.numu_cc_1_n_daughter_all.push_back(static_cast<float>(n_daughter_all));

            if (flag_numu_cc_1) flag_numu_cc_1_save = true;
        }
    }

    // -----------------------------------------------------------------------
    // Flag 2: long-muon shower check.
    // Checks showers whose start segment is a muon in the main cluster.
    // Prototype: NeutrinoID_numu_tagger.h lines 69-104.
    // -----------------------------------------------------------------------
    for (ShowerPtr shower : showers) {
        bool flag_numu_cc_2 = false;
        SegmentPtr start_sg = shower->start_segment();
        if (!start_sg) continue;

        int sg_pdg   = start_sg->has_particle_info() ? start_sg->particle_info()->pdg() : 0;
        int sg_cl_id = start_sg->cluster() ? start_sg->cluster()->get_cluster_id() : -2;

        if (std::abs(sg_pdg) == 13 && main_cl_id == sg_cl_id) {
            double length       = shower->get_total_track_length();
            double total_length = shower->get_total_length();

            auto [n_daughter_tracks, n_daughter_all] =
                count_daughters(graph, shower, main_vertex, segments_in_long_muon_plain);

            if (length > 18*units::cm &&
                !(n_daughter_tracks > 1 || n_daughter_all - n_daughter_tracks > 2))
                flag_numu_cc_2 = true;

            if (length > max_muon_length) {
                max_muon_length = length;
                max_long_muon   = shower;
                max_muon        = nullptr;
            }

            ti.numu_cc_flag_2.push_back(static_cast<float>(flag_numu_cc_2));
            ti.numu_cc_2_length.push_back(static_cast<float>(length / units::cm));
            ti.numu_cc_2_total_length.push_back(static_cast<float>(total_length / units::cm));
            ti.numu_cc_2_n_daughter_tracks.push_back(static_cast<float>(n_daughter_tracks));
            ti.numu_cc_2_n_daughter_all.push_back(static_cast<float>(n_daughter_all));

            if (flag_numu_cc_2) flag_numu_cc_2_save = true;
        }
    }

    // -----------------------------------------------------------------------
    // Flag 3: any muon-like track in the main cluster, not directly at main_vertex.
    // Also accumulates acc_track_length and max_length_all over ALL non-shower
    // segments in the main cluster (including those at main_vertex).
    // Prototype: NeutrinoID_numu_tagger.h lines 108-209.
    // -----------------------------------------------------------------------
    double max_length      = 0;   // longest muon candidate length
    double acc_track_length = 0;  // total length of all tracks in main cluster
    double max_length_all   = 0;  // longest non-shower segment in main cluster
    SegmentPtr tmp_max_muon = nullptr;

    for (auto [eit, eit_end] = boost::edges(graph); eit != eit_end; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || !sg->cluster()) continue;
        if (sg->cluster()->get_cluster_id() != main_cl_id) continue;

        double length = segment_track_length(sg);

        // Update accumulators for ALL non-shower segments in main cluster
        // (done before the main_vertex skip below, matching prototype ordering).
        if (!seg_is_shower(sg) && max_length_all < length) max_length_all = length;
        if (!seg_is_shower(sg)) acc_track_length += length;

        // Segments directly connected to main_vertex are handled in flag 1 — skip here.
        if (seg_at_main_vertex(sg, main_vertex, graph)) continue;

        double direct_length = segment_track_direct_length(sg);
        double medium_dQ_dx  = segment_median_dQ_dx(sg);
        double dQ_dx_cut     = 0.8866 + 0.9533 * std::pow(18*units::cm / length, 0.4234);

        if (sg->flags_any(SegmentFlags::kAvoidMuonCheck)) continue;

        bool is_shower_topo = sg->flags_any(SegmentFlags::kShowerTopology);

        // Prototype condition: (shower && !topo) || !shower || length > 50cm
        // = !is_shower_topo || length > 50cm (since !shower implies !topo in this context)
        if (!is_shower_topo || length > 50*units::cm) {
            if (medium_dQ_dx < dQ_dx_cut * 1.05 * 43e3/units::cm &&
                medium_dQ_dx > 0.75 * 43e3/units::cm &&
                (direct_length > 0.925 * length || length > 50*units::cm)) {

                if (length > max_length) {
                    max_length   = length;
                    tmp_max_muon = sg;
                }

                // For long segments (> 50 cm), try extending through find_cont_muon_segment
                // to estimate the full muon length for max_muon_length tracking.
                // Prototype: NeutrinoID_numu_tagger.h lines 148-169.
                if (length > 50*units::cm) {
                    double tmp_length = length;
                    if (!is_shower_topo) {
                        auto [pv1, pv2] = find_vertices(graph, sg);
                        if (pv1) {
                            SegmentPtr ext_sg1 = find_cont_muon_segment(graph, sg, pv1).first;
                            if (ext_sg1) tmp_length += segment_track_length(ext_sg1);
                        }
                        if (pv2) {
                            SegmentPtr ext_sg2 = find_cont_muon_segment(graph, sg, pv2).first;
                            if (ext_sg2) tmp_length += segment_track_length(ext_sg2);
                        }
                    }
                    if (tmp_length > max_muon_length) {
                        max_muon_length = tmp_length;
                        max_muon        = sg;
                        max_long_muon   = nullptr;
                    }
                }
            }
        }
    }

    // Evaluate flag 3 using the longest non-main-vertex muon candidate.
    auto [n3_tracks, n3_all] = count_daughters(graph, tmp_max_muon, main_vertex);
    int tmp_particle_type = 0;

    if (((max_length > 25*units::cm && acc_track_length > 0) ||
         (max_length > 30*units::cm && acc_track_length == 0)) &&
        (max_length_all - max_muon_length <= 100*units::cm) &&
        !(n3_tracks > 1 || n3_all - n3_tracks > 2)) {

        tmp_particle_type = (tmp_max_muon && tmp_max_muon->has_particle_info())
                            ? tmp_max_muon->particle_info()->pdg() : 0;

        if (tmp_particle_type != 211) {
            if (max_length > max_muon_length) {
                max_muon_length = max_length;
                max_muon        = tmp_max_muon;
                max_long_muon   = nullptr;
            }
            flag_numu_cc_3 = true;
        }
    }

    ti.numu_cc_flag_3          = static_cast<float>(flag_numu_cc_3);
    ti.numu_cc_3_particle_type = static_cast<float>(tmp_particle_type);
    ti.numu_cc_3_max_length    = static_cast<float>(max_length / units::cm);
    ti.numu_cc_3_acc_track_length = static_cast<float>(acc_track_length / units::cm);
    ti.numu_cc_3_max_length_all= static_cast<float>(max_length_all / units::cm);
    ti.numu_cc_3_max_muon_length = static_cast<float>(max_muon_length / units::cm);
    ti.numu_cc_3_n_daughter_tracks = static_cast<float>(n3_tracks);
    ti.numu_cc_3_n_daughter_all    = static_cast<float>(n3_all);

    // -----------------------------------------------------------------------
    // Final numu CC decision.
    // Prototype: NeutrinoID_numu_tagger.h lines 240-260.
    // neutrino_type bit-setting is omitted; caller uses the returned bool.
    // -----------------------------------------------------------------------
    flag_numu_cc = flag_numu_cc_1_save || flag_numu_cc_2_save || flag_numu_cc_3;
    ti.numu_cc_flag = static_cast<float>(flag_numu_cc);

    if ((max_muon_length > 100*units::cm || max_length_all > 120*units::cm) && flag_numu_cc)
        flag_long_muon = true;

    return {flag_long_muon, max_muon_length};
}
