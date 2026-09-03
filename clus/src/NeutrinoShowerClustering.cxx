#include "WireCellClus/ExaminerPassBudget.h"
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
// doc sbnd_xin/docs/pr/130 item 4 part 10 -- WCT_SHOWER_XCLUS_DEBUG.
// Part 9 proved the shower WALK never reaches the charge the owner's scan
// wants (0 of 78 contention, 96.2% of it in clusters the shower does not hold
// at all), so the whole of the approved under-clustering is CROSS-CLUSTER
// acquisition.  The cross-cluster absorbers tape admissions only, so a segment
// that shower A wanted and shower B took leaves no record of A ever having
// looked -- exactly the gap the BLOCKED probe closed one level down.  This
// tape records the REJECTIONS at the pass-4 direct-cone acceptance, with the
// quantities the disjunction tests, plus every rival's arbitration metric.
// It answers: was A a candidate that lost, or never a candidate at all?
// stderr only; separate env var so prep_em_scan.py's ABSORB stream is untouched.
static inline bool pr130_xclus_dbg()
{
    static const bool dbg = std::getenv("WCT_SHOWER_XCLUS_DEBUG") != nullptr;
    return dbg;
}
// doc sbnd_xin/docs/pr/132 -- WCT_PI0_PAIR_DEBUG, the pi0-finder attribution
// tape.  The pr/126 census says WHICH hand pi0 the finders miss but not WHY:
// the pair pool, the 30-degree association verdicts, the mass-window/ranking
// verdicts and the without-vertex pre-gate that killed the event are all
// invisible in the outputs.  This probe prints each of those decisions.
// Shower ids use the calib-dump "id" encoding (pf_node_id ==
// pr91_seg_display_id of the start segment) so lines join to labels and
// census rows directly.  stderr only, no effect on emitted bytes.
static inline bool pr132_pi0_dbg()
{
    static const bool dbg = std::getenv("WCT_PI0_PAIR_DEBUG") != nullptr;
    return dbg;
}
// doc pr/132 round 10: the opening-angle census tape.  For every recorded
// path-1 pair, prints the LEGACY pair mass/angle (local_dirs: init_dir for
// attached, start-vtx ray for associated) next to the centroid-ray variant
// (vertex -> associate-cloud centroid for both members).  Round 7 measured
// the below-window true-pair population as angle-compressed (103798:
// charges 0.96/1.02 of label, mass 74) -- this tape measures whether the
// centroid ray decompresses it.  stderr only, byte-neutral.
static inline bool pr132_pi0_angle_dbg()
{
    static const bool dbg = std::getenv("WCT_PI0_ANGLE_DEBUG") != nullptr;
    return dbg;
}
static inline int pr132_pi0_shid(const ShowerPtr& s)
{
    return s ? pr91_seg_display_id(s->start_segment()) : -1;
}
// doc sbnd_xin/docs/pr/132 round 3 -- WCT_PI0_SUBSTRUCT_DEBUG, the over-merge
// substructure tape (132 doc sec 9.9 idea 2).  For every shower above 300 MeV
// after the two pi0 finders, PCA-decompose the associate-points cloud and
// report the two-cluster split along each of the two leading axes:
// count-fraction energies, the angle between the two half-centroid directions
// from the shower start, and the implied two-gamma mass.  Per-point charge is
// not available here, so energies are point-count-weighted -- this is a
// counting probe (how many over-merged events are recoverable at all), not a
// measurement.  stderr only, no effect on emitted bytes.
static inline bool pr132_pi0_substruct_dbg()
{
    static const bool dbg = std::getenv("WCT_PI0_SUBSTRUCT_DEBUG") != nullptr;
    return dbg;
}
// doc sbnd_xin/docs/pr/132 round 4 -- WCT_PI0_NCVTX_DEBUG, two tapes in one
// env (132 doc sec 10.9 fronts 1+2; owner scan verdicts sec 10.8):
//   PI0_NCVTX -- the NC vertex-in-shower recognizer.  Per shower above
//     50 MeV: distance main-vertex -> nearest associate-cloud point, distance
//     main-vertex -> shower start, and f_back = the fraction of associate
//     points BEHIND the main vertex along the vertex->centroid axis.  A
//     gamma truly starting at the vertex has f_back ~ 0; a vertex mis-seated
//     inside one arm (owner: 76346, 116962) has charge on both sides.
//   PI0_START -- the accepted-shower start-derivation audit.  Per pi0-paired
//     shower: distance start-vertex -> current start point, -> nearest
//     "fit"-cloud point (the conn-2 derivation, PRShower calculate_kinematics)
//     and -> nearest "associate_points" cloud point, plus the gap between the
//     two cloud answers.  169626 measured 13.4 cm of start bias (owner scan);
//     this tape tells whether the associate cloud closes it.
// stderr only, no effect on emitted bytes.
static inline bool pr132_pi0_ncvtx_dbg()
{
    static const bool dbg = std::getenv("WCT_PI0_NCVTX_DEBUG") != nullptr;
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

// doc sbnd_xin/docs/pr/121 round 1: the examine_shower_1 accept-branch dedup
// erases ANY pre-existing (main_vertex, conn-1) shower sharing the accepted
// shower1's start segment, with no re-homing.  Written for stale
// single-segment wrappers, its predicate never checks the segment count --
// SBND 17394-348471: it erased a 13-member 352.6 MeV shower (retargeted onto
// proton seg 12052 by examine_showers_retarget_seed), orphaning 12 EM
// segments from PF output (doc pr/115 sec 17.7).  Prints every dedup
// candidate with its size so the census can count multi-segment erasures.
// stderr only, no effect on emitted bytes.
static inline void pr121_probe_ex1_dedup(const ShowerPtr& keep, const ShowerPtr& old,
                                         int old_ctype, int old_svtx_main, int erase)
{
    if (!pr93_absorb_dbg() || !keep || !old) return;
    std::fprintf(stderr,
                 "SHOWER_ABSORB EX1_DEDUP into_start_seg=%d old_shower_id=%d old_nseg=%d "
                 "old_ctype=%d old_svtx_main=%d old_kine_mev=%.1f erase=%d\n",
                 pr91_seg_display_id(keep->start_segment()), old->get_shower_id(),
                 old->get_num_segments(), old_ctype, old_svtx_main,
                 old->get_kine_charge() / WireCell::units::MeV, erase);
}

// doc sbnd_xin/docs/pr/122 round 1: the in_main_cluster seeder accepts a
// segment as a shower root on the flag disjunction alone -- no length,
// straightness, or dQ/dx test (SBND 18255-54332: a straight 32.3 cm track
// carrying a kShowerTopology mis-flag seeded a 2.89e6 q_extra shower).  One
// line per accepted seed with the features a guard could use; the feature
// calls run only under the env gate.  stderr only, no effect on emitted bytes.
static inline void pr122_probe_seed(const SegmentPtr& sg, bool f_traj, bool f_topo,
                                    bool f_pdg11, bool in_long_muon, double mip_dqdx_median)
{
    if (!pr93_absorb_dbg() || !sg) return;
    const double len = segment_track_length(sg);
    const double med = segment_median_dQ_dx(sg);
    std::fprintf(stderr,
                 "SHOWER_SEED site=in_main_cluster seg=%d gidx=%zu pdg=%d traj=%d topo=%d pdg11=%d "
                 "long_muon=%d len_cm=%.2f med_dqdx_mip=%.3f straight=%d\n",
                 pr91_seg_display_id(sg), sg->get_graph_index(),
                 sg->has_particle_info() && sg->particle_info() ? sg->particle_info()->pdg() : 0,
                 (int)f_traj, (int)f_topo, (int)f_pdg11, (int)in_long_muon,
                 len / WireCell::units::cm, mip_dqdx_median > 0 ? med / mip_dqdx_median : -1.0,
                 (int)segment_is_straight_long_track(sg));
}

// doc sbnd_xin/docs/pr/123 round 1: the pass4_angle acceptance in
// shower_clustering_with_nv_from_vertices has no term for the shower's own
// extent or contiguity -- a stub 80 cm from the start vertex is absorbed even
// when detached from the body (doc pr/115 sec 16.5: 48% of wrongly-held
// marks; sec 17.6: 4.80e7 q_extra untouched by every shipped knob).  Owner
// over-reach definition 2026-08-28: a far member counts as over-reach when
// detached from the contiguous body (gap) OR track-like beyond the body.
// One line per ACCEPTED segment with every quantity in the acceptance
// disjunction plus the features that definition needs; all geometry is
// already computed unconditionally at the call site, only the dQ/dx median
// runs under the env gate.  stderr only, no effect on emitted bytes.
// Probe v2: body_dis (to the CURRENT body) is chain-laundered -- once a first
// far stub is absorbed, later stubs sit close to it and print small gaps
// (SBND 278420: first stub body_dis 23.9, the six followers 5-16).  snap_dis
// is the gap to the shower's membership AT LOOP ENTRY (chain-immune), the
// guard-usable form of the owner rule's "detached from the contiguous body".
static inline void pr123_probe_pass4_geom(const SegmentPtr& sg, const ShowerPtr& cur,
                                          const ShowerPtr& owner, double pair_dis,
                                          double front_dis, double body_dis,
                                          double snap_dis,
                                          double angle_v1, double angle_v2,
                                          double mip_dqdx_median)
{
    if (!pr93_absorb_dbg() || !sg || !cur || !owner) return;
    const double cm = WireCell::units::cm;
    int tier = 0;
    if      (angle_v1 < 25   && (pair_dis < 80 * cm  || body_dis < 25 * cm)) tier = 1;
    else if (angle_v2 < 25   && (front_dis < 40 * cm || body_dis < 25 * cm)) tier = 2;
    else if (angle_v1 < 12.5 && (pair_dis < 120 * cm || body_dis < 40 * cm)) tier = 3;
    else if (angle_v2 < 12.5 && (front_dis < 80 * cm || body_dis < 40 * cm)) tier = 4;
    const double med = segment_median_dQ_dx(sg);
    std::fprintf(stderr,
                 "SHOWER_ABSORB PASS4_GEOM seg=%d pdg=%d len_cm=%.2f med_dqdx_mip=%.3f "
                 "cur=%d cur_nseg=%d cur_len_cm=%.1f owner=%d divert=%d "
                 "pair_dis_cm=%.2f front_dis_cm=%.2f body_dis_cm=%.2f snap_dis_cm=%.2f "
                 "angle_v1=%.2f angle_v2=%.2f tier=%d\n",
                 pr91_seg_display_id(sg),
                 sg->has_particle_info() && sg->particle_info() ? sg->particle_info()->pdg() : 0,
                 segment_track_length(sg) / cm,
                 mip_dqdx_median > 0 ? med / mip_dqdx_median : -1.0,
                 pr91_seg_display_id(cur->start_segment()), cur->get_num_segments(),
                 cur->get_total_length() / cm,
                 pr91_seg_display_id(owner->start_segment()), (int)(owner != cur),
                 pair_dis / cm, front_dis / cm, body_dis / cm, snap_dis / cm,
                 angle_v1, angle_v2, tier);
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
            // doc sbnd_xin/docs/pr/120 -- byte-neutral admission census.  The
            // hand-scan OUT marks admitted here sit at 146-148 deg to the
            // FINAL shower axis; this prints the scan-equivalent angle (shower
            // dir15/dir60 at the CURRENT pre-re-seat start vs start->closest
            // stem point) for EVERY chain candidate, accepted or not, so the
            // census can measure whether admission-time geometry already
            // separates.  Env-gated stderr only: no effect on emitted bytes.
            if (pr93_absorb_dbg()) {
                const WireCell::Point sp0 = shower->get_start_point();
                const WireCell::Vector s15 = shower_cal_dir_3vector(*shower, sp0, 15 * units::cm);
                const WireCell::Vector s60 = shower_cal_dir_3vector(*shower, sp0, 60 * units::cm);
                auto [sdist, scp] = segment_get_closest_point(stem, sp0);
                const WireCell::Vector sv(scp.x() - sp0.x(), scp.y() - sp0.y(), scp.z() - sp0.z());
                auto p120_ang = [](const WireCell::Vector& ax, const WireCell::Vector& v) {
                    if (ax.magnitude() < 0.001 || v.magnitude() < 0.001) return -1.0;
                    return std::acos(std::clamp(ax.dot(v) / (ax.magnitude() * v.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                };
                // doc sbnd_xin/docs/pr/130 part 5 -- VERTEX-RELATIVE geometry.
                // Part 4 measured every feature above (pdg, len, dQ/dx ratio,
                // ang15, ang60, dist) as interleaved between the absorb-wanted
                // and decline-ok classes over the guard's COMPLETE 8-candidate
                // firing set, so no admission-time cut separates them.  Every
                // one of those is local to the shower/stem pair; none
                // references the neutrino vertex.  The physics that is missing:
                // a stem that is the shower's true PARENT lies BETWEEN the
                // vertex and the shower, so it is closer to the vertex than the
                // shower start (toward > 0) and sits on nearly the same ray out
                // of the vertex (vang small); a separate hadronic prong sits AT
                // the vertex and leaves along a different ray (vang large).
                // Same env gate, stderr only -- no effect on emitted bytes.
                double dvtx_start = -1.0, dvtx_stem = -1.0, vang = -1.0, toward = 0.0;
                if (main_vertex && main_vertex->descriptor_valid()) {
                    const WireCell::Point vp = main_vertex->fit().point;
                    const WireCell::Vector a(sp0.x() - vp.x(), sp0.y() - vp.y(), sp0.z() - vp.z());
                    dvtx_start = a.magnitude();
                    WireCell::Point best;
                    for (const auto& f : stem->fits()) {
                        const WireCell::Vector d(f.point.x() - vp.x(), f.point.y() - vp.y(),
                                                 f.point.z() - vp.z());
                        const double dm = d.magnitude();
                        if (dvtx_stem < 0 || dm < dvtx_stem) { dvtx_stem = dm; best = f.point; }
                    }
                    if (dvtx_stem >= 0) {
                        // >0 => absorbing this stem pulls the shower start
                        // TOWARD the neutrino vertex.
                        toward = dvtx_start - dvtx_stem;
                        const WireCell::Vector b(best.x() - vp.x(), best.y() - vp.y(),
                                                 best.z() - vp.z());
                        vang = p120_ang(a, b);
                    }
                }
                std::fprintf(stderr,
                             "SHOWER_ABSORB P120_STEM shower_start_seg=%d conn=%d seg=%d pdg=%d "
                             "len_cm=%.2f ratio=%.2f ok=%d ang15=%.2f ang60=%.2f dist_cm=%.2f "
                             "dvtx_start_cm=%.2f dvtx_stem_cm=%.2f toward_cm=%.2f vang=%.2f\n",
                             pr91_seg_display_id(start_seg), conn_type, pr91_seg_display_id(stem),
                             stem->has_particle_info() && stem->particle_info()
                                 ? stem->particle_info()->pdg() : 0,
                             len / units::cm, ratio, (int) ok,
                             p120_ang(s15, sv), p120_ang(s60, sv), sdist / units::cm,
                             dvtx_start < 0 ? -1.0 : dvtx_start / units::cm,
                             dvtx_stem  < 0 ? -1.0 : dvtx_stem  / units::cm,
                             toward / units::cm, vang);
            }
            if (!ok) break;
            // doc sbnd_xin/docs/pr/120 -- backward-stem guard.  Measured over
            // the 98-event emscan manifest: stem_backfill accepted 5 absorbs;
            // the 3 with a degenerate scan-equivalent angle (conn-1 showers
            // whose start already sits on the chain, dist=0) are trunk
            // extensions in "good"-note events, while BOTH measurable-angle
            // absorbs develop backward (~150 deg: the shower points away from
            // the chain) and both are scanner-condemned over-clustering
            // (evt47212 seg 2103 OUT-marked; evt281567 seg 95128 named in the
            // scan note).  Decline the absorb -- and stop the chain -- when
            // the angle is measurable and beyond the cut.  C++ default false
            // => byte-identical.
            if (m_stem_backfill_back_guard) {
                const WireCell::Point gsp = shower->get_start_point();
                const WireCell::Vector g15 = shower_cal_dir_3vector(*shower, gsp, 15 * units::cm);
                auto [gdist, gcp] = segment_get_closest_point(stem, gsp);
                const WireCell::Vector gv(gcp.x() - gsp.x(), gcp.y() - gsp.y(), gcp.z() - gsp.z());
                double gang = -1.0;
                if (g15.magnitude() > 0.001 && gv.magnitude() > 0.001) {
                    gang = std::acos(std::clamp(g15.dot(gv) / (g15.magnitude() * gv.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                }
                (void) gdist;
                // doc sbnd_xin/docs/pr/130 item B (stem_backfill_back_dvtx):
                // the owner ruled the guard WRONG on 2 of the 8 candidates it
                // can reach -- 292643 (declines an absorb that should happen,
                // -234.0 MeV) and 179369 (the decline leaves a spurious pi0,
                // +376.0 MeV).  The separator is the distance from the SHOWER
                // START to the neutrino vertex: a shower that starts far from
                // the vertex is detached, so a backward stem is plausibly its
                // parent; one that starts at the vertex has a backward stem
                // that is a separate prong.
                //
                // The 239-event P120_STEM census at the 2026-08-29 production
                // point proves this population is CLOSED: `if (!ok) break;`
                // above means the guard only ever sees chains the MIP window
                // accepted, and there are exactly 9 such rows in 8 events over
                // both manifests -- precisely the set the owner adjudicated.
                // The 10 further rows meeting the angle condition all carry
                // ok=0 and never reach here.  So this knob cannot touch an
                // unadjudicated candidate: absorb-wanted are 46.84 and 88.11
                // cm, every decline-ok is <= 44.34 cm.  0 = off = legacy.
                if (m_stem_backfill_back_dvtx > 0 && main_vertex
                    && main_vertex->descriptor_valid()) {
                    const WireCell::Point mvp = main_vertex->fit().point;
                    const WireCell::Vector dv(gsp.x() - mvp.x(), gsp.y() - mvp.y(),
                                              gsp.z() - mvp.z());
                    if (dv.magnitude() > m_stem_backfill_back_dvtx) {
                        if (pr93_absorb_dbg()) {
                            std::fprintf(stderr,
                                "SHOWER_ABSORB P130_BACK_DVTX shower_start_seg=%d seg=%d "
                                "dvtx_start_cm=%.2f > %.2f -- decline suppressed\n",
                                pr91_seg_display_id(start_seg), pr91_seg_display_id(stem),
                                dv.magnitude() / units::cm,
                                m_stem_backfill_back_dvtx / units::cm);
                        }
                        SPDLOG_LOGGER_DEBUG(s_log,
                            "pr130 stem_backfill_back_dvtx: suppress decline seg={} "
                            "dvtx_start={:.1f}cm > {:.1f}cm",
                            pr91_seg_display_id(stem), dv.magnitude() / units::cm,
                            m_stem_backfill_back_dvtx / units::cm);
                        gang = -1.0;   // fall through to the normal absorb path
                    }
                }
                if (gang >= 0 && gang > m_stem_backfill_back_ang) {
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr120 stem_backfill_back_guard: decline shower(start gidx={} conn={}) "
                        "chain gidx={} ang15={:.1f} > {:.1f}",
                        start_seg->get_graph_index(), conn_type, stem->get_graph_index(),
                        gang, m_stem_backfill_back_ang);
                    break;
                }
            }
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
    if (!main_vertex || !main_vertex->descriptor_valid()) {
        // doc pr/122 round 1 (case c): answers "was seeding ever attempted".
        if (pr93_absorb_dbg())
            std::fprintf(stderr, "SHOWER_SEED site=in_main_cluster ABORT no_main_vertex=%d\n",
                         main_vertex ? 0 : 1);
        return;
    }

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
            // (doc pr/122 round 1: disjuncts named so the seed probe can
            // report which one admitted the root; no semantic change.)
            const bool f_traj = curr_sg->flags_any(SegmentFlags::kShowerTrajectory);
            const bool f_topo = curr_sg->flags_any(SegmentFlags::kShowerTopology);
            const bool f_pdg11 = (curr_sg->has_particle_info() && curr_sg->particle_info() &&
                                  std::abs(curr_sg->particle_info()->pdg()) == 11);
            bool is_shower_seg = f_traj || f_topo || f_pdg11;
            bool in_long_muon = segments_in_long_muon.count(curr_sg) > 0;

            if (is_shower_seg || in_long_muon) {
                pr122_probe_seed(curr_sg, f_traj, f_topo, f_pdg11, in_long_muon, m_mip_dqdx_median);
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
            // pion fails.  doc 77 round 1 (2026-08-24):
            // shower_connect_protected_pion_guard removed -- measured dead,
            // never flipped (pr/40 sec 1459).  See
            // sbnd_xin/docs/77_knob-ledger.tsv.

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
        ExaminerPassCounter epb_shower_clustering_with_nv_from_main_cluster("shower_clustering_with_nv_from_main_cluster", main_cluster ? main_cluster->get_cluster_id() : -1);   // doc pdvd/26 round 2
        while (flag_continue) {
            if (epb_shower_clustering_with_nv_from_main_cluster.exceeded()) break;
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
        const ShowerDirInfo* min_info = nullptr;   // doc pr/120: census fields for the winner

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
                    min_info = &info;
                }
            }
        }

        // doc sbnd_xin/docs/pr/120 -- byte-neutral admission census for the
        // winning (segment, shower) pair: the site's own angle/dist (recomputed
        // from min_info, identical arithmetic to the loop) plus the
        // scan-equivalent angle (shower_cal_dir_3vector at the same start, 15
        // and 60 cm) so the census can measure how the two frames diverge --
        // the OUT marks admitted here sit at 113-134 deg to the FINAL axis yet
        // passed the <=30 deg cone above.  Env-gated stderr only.
        if (min_shower && min_info && pr93_absorb_dbg()) {
            auto [pdist, pcp] = segment_get_closest_point(seg1, min_info->start_point);
            const WireCell::Vector pv(pcp.x() - min_info->start_point.x(),
                                      pcp.y() - min_info->start_point.y(),
                                      pcp.z() - min_info->start_point.z());
            auto p120_ang = [](const WireCell::Vector& ax, const WireCell::Vector& v) {
                if (ax.magnitude() < 0.001 || v.magnitude() < 0.001) return -1.0;
                return std::acos(std::clamp(ax.dot(v) / (ax.magnitude() * v.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            };
            const WireCell::Vector c15 = shower_cal_dir_3vector(*min_shower, min_info->start_point, 15 * units::cm);
            const WireCell::Vector c60 = shower_cal_dir_3vector(*min_shower, min_info->start_point, 60 * units::cm);
            std::fprintf(stderr,
                         "SHOWER_ABSORB P120_P3CONE seg=%d pdg=%d len_cm=%.2f shower_start_seg=%d "
                         "site_ang=%.2f dist_cm=%.2f ang15=%.2f ang60=%.2f ao=%.1f\n",
                         pr91_seg_display_id(seg1),
                         seg1->has_particle_info() && seg1->particle_info()
                             ? seg1->particle_info()->pdg() : 0,
                         segment_track_length(seg1) / units::cm,
                         pr91_seg_display_id(min_shower->start_segment()),
                         p120_ang(min_info->dir, pv), pdist / units::cm,
                         p120_ang(c15, pv), p120_ang(c60, pv), min_info->angle_offset);
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
            // doc sbnd_xin/docs/pr/124 front C (shower_pass3_cone_guard_len):
            // the post-flip census over both manifests measured the 9
            // surviving labeled-OUT cone absorbs as ALL track-pdg
            // (13/211/2212: 52693 34cm mu at 3.8deg/114cm, 94392 29.8+46.8cm
            // mu, 175896 5.6+17.6cm p, ...), while labeled-IN track-pdg
            // members top out at 14.1cm (76346).  Decline the cone absorb
            // for a track-pdg segment longer than the knob; the pass4 guard
            // (50cm) does not cover this earlier site.  The declined seg is
            // flagged kPass4GuardFreed so the (SBND-ON) PF/kine orphan
            // passes claim it if nothing else does (the pr/123 r2 lesson).
            // 0 = off = legacy.
            bool p3_track_guard = false;
            if (m_shower_pass3_cone_guard_len > 0 &&
                segment_track_length(seg1) > m_shower_pass3_cone_guard_len &&
                seg1->has_particle_info() && seg1->particle_info()) {
                const int p3_pdg = std::abs(seg1->particle_info()->pdg());
                p3_track_guard = (p3_pdg == 13 || p3_pdg == 211 || p3_pdg == 2212);
            }
            if (cone_guard_excludes) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr93 cone_absorb_guard: decline absorb seg={} pdg={} len={:.1f}cm",
                    pr91_seg_display_id(seg1), seg1->particle_info()->pdg(),
                    segment_track_length(seg1)/units::cm);
            }
            else if (p3_track_guard) {
                seg1->set_flags(SegmentFlags::kPass4GuardFreed);
                if (pr93_absorb_dbg()) {
                    std::fprintf(stderr,
                        "SHOWER_ABSORB PASS3_CONE_GUARD seg=%d pdg=%d len_cm=%.1f declined=1\n",
                        pr91_seg_display_id(seg1), seg1->particle_info()->pdg(),
                        segment_track_length(seg1) / units::cm);
                }
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr124 pass3_cone_guard: decline seg={} pdg={} len={:.1f}cm",
                    pr91_seg_display_id(seg1), seg1->particle_info()->pdg(),
                    segment_track_length(seg1) / units::cm);
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
            // doc sbnd_xin/docs/pr/130 item 1b (shower_pass3_backfill_guard_len):
            // pr/93's predicate was given a second seat here (the block
            // above, SBND 315167) but pr/124's shower_pass3_cone_guard_len
            // was not -- so a track-pdg segment the pass-3 cone absorber
            // declined ~40 lines above is re-adopted into the SAME shower
            // here, and force-relabelled pdg 11 just below.  That defeats
            // pr/124's own stated intent: its decline sets kPass4GuardFreed
            // precisely so the PF/kine orphan passes can claim the segment,
            // not so this backfill can put it straight back.  Census over
            // both standard manifests: 6 segments in 5 events (137238
            // seg 145067, 175896 seg 66041, 176502 seg 109123, 396222 segs
            // 9080/9098, 415278 seg 24072), every one of them
            // pass3_cone_guard -> pass3_cluster_map.  Same predicate as
            // pr/124's; the flag is already set, so this only has to stop
            // the re-adoption.  0 = off = legacy => byte-identical.
            if (m_shower_pass3_backfill_guard_len > 0 &&
                segment_track_length(seg1) > m_shower_pass3_backfill_guard_len &&
                seg1->has_particle_info() && seg1->particle_info()) {
                const int bf_pdg = std::abs(seg1->particle_info()->pdg());
                if (bf_pdg == 13 || bf_pdg == 211 || bf_pdg == 2212) {
                    seg1->set_flags(SegmentFlags::kPass4GuardFreed);
                    if (pr93_absorb_dbg()) {
                        std::fprintf(stderr,
                            "SHOWER_ABSORB PASS3_BACKFILL_GUARD seg=%d pdg=%d len_cm=%.1f declined=1\n",
                            pr91_seg_display_id(seg1), seg1->particle_info()->pdg(),
                            segment_track_length(seg1) / units::cm);
                    }
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr130 pass3_backfill_guard: decline seg={} pdg={} len={:.1f}cm",
                        pr91_seg_display_id(seg1), seg1->particle_info()->pdg(),
                        segment_track_length(seg1) / units::cm);
                    continue;
                }
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
    if (m_shower_pass4_best_owner) {
        SPDLOG_LOGGER_DEBUG(s_log, "pr117 pass4_best_owner on: direct-cone accepts arbitrated over all showers");
    }
    
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
        // doc pr/97 D1: min_angle/min_dis/min_point are left INDETERMINATE by
        // the legacy path (and by the prototype) and are only written below if
        // main_vertex happens to be one of main_cluster_vertices.  When it is
        // not, the branch after the loop reads stale stack bytes.  Sentinels
        // make that case deterministically prefer min_pi.  Statement never
        // runs when the knob is off => legacy path bit-for-bit.
        if (m_shower_nv_main_pi_init) {
            main_pi.min_angle = 1e9;
            main_pi.min_dis = 1e9;
            main_pi.min_point.set(0, 0, 0);
        }
        
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

        // doc sbnd_xin/docs/pr/117 round 1 (shower_pass4_best_owner): rival
        // list for the direct-cone arbitration below -- every OTHER existing
        // shower (pass 1/2/3 and earlier pass-4 ones; this pass runs per
        // candidate cluster nearest-first, so earlier winners are visible
        // here), with its start point and 30cm local direction, in
        // cluster-id/segment-id order (the pass-4 proximity comparator) for
        // determinism.  Near-zero directions (1-segment stubs) are skipped so
        // they cannot steal.  Also built under the WCT_SHOWER_ABSORB_DEBUG
        // probe alone (knob off) to MEASURE the would-divert rate; the
        // divert itself only happens under the knob.  Knob off + no probe =>
        // list stays empty => the accept below is the legacy greedy add =>
        // byte-identical.
        struct Pass4Rival { ShowerPtr shower; WireCell::Point start_pt; WireCell::Vector dir; };
        std::vector<Pass4Rival> pass4_rivals;
        if (m_shower_pass4_best_owner || pr93_absorb_dbg()) {
            std::vector<ShowerPtr> rival_order(showers.begin(), showers.end());
            std::sort(rival_order.begin(), rival_order.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
                auto* sa = a->start_segment().get();
                auto* sb = b->start_segment().get();
                if (!sa || !sb) return sa < sb;
                int cid_a = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
                int cid_b = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
                if (cid_a != cid_b) return cid_a < cid_b;
                return sa->id() < sb->id();
            });
            for (auto& s1 : rival_order) {
                if (s1 == shower) continue;
                auto [rvtx, rtype] = s1->get_start_vertex_and_type();
                WireCell::Point sp1 = rvtx ? (rvtx->fit().valid() ? rvtx->fit().point : rvtx->wcpt().point) : WireCell::Point(0, 0, 0);
                auto [rd_, rtp] = shower_get_closest_point(*s1, sp1);
                WireCell::Vector rdir = shower_cal_dir_3vector(*s1, rtp, 30 * units::cm);
                if (rdir.magnitude() < 0.001) continue;
                pass4_rivals.push_back({s1, sp1, rdir});
            }
        }

        // Add segments from other clusters.
        // Cache the shower start front point to avoid re-evaluating per segment.
        const WireCell::Point shower_start_front = shower->start_segment()->fits().front().point;

        // doc pr/123: pass-entry membership snapshot for the probe's
        // chain-immune gap (env-gated, read-only, empty when probes off).
        std::vector<SegmentPtr> pr123_snapshot;
        if (pr93_absorb_dbg()) {
            for (auto edesc : shower->edges()) {
                auto sg_ = shower->view_graph()[edesc].segment;
                if (sg_) pr123_snapshot.push_back(sg_);
            }
        }
        for (auto seg1 : seg_order) {
            if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
            // doc pr/130 part 10: the decisive skip.  A segment already owned by
            // ANOTHER shower is dropped here BEFORE any geometry is computed, so the
            // owning shower's claim is invisible to every later shower -- the same
            // silent-skip shape the BLOCKED probe found in the walk one level down.
            if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) {
                if (pr130_xclus_dbg()) {
                    const auto& owner_sh = map_segment_in_shower.find(seg1)->second;
                    std::fprintf(stderr,
                        "SHOWER_XCLUS OWNED site=pass4_angle shower=%d seg=%d owner=%d\n",
                        pr91_seg_display_id(shower->start_segment()),
                        pr91_seg_display_id(seg1),
                        owner_sh ? pr91_seg_display_id(owner_sh->start_segment()) : -1);
                }
                continue;
            }
            if (seg1->cluster() == shower->start_segment()->cluster()) continue;

            auto it1 = map_cluster_associated_vertex.find(seg1->cluster());

            auto [pair_dis, pair_point] = segment_get_closest_point(seg1, start_pt);
            WireCell::Vector v1(pair_point.x() - start_pt.x(), pair_point.y() - start_pt.y(), pair_point.z() - start_pt.z());
            WireCell::Vector v2(pair_point.x() - point.x(), pair_point.y() - point.y(), pair_point.z() - point.z());

            double angle_v1 = std::acos(std::clamp(dir_shower.dot(v1) / (dir_shower.magnitude() * v1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            double angle_v2 = std::acos(std::clamp(dir_shower.dot(v2) / (dir_shower.magnitude() * v2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;

            // Filter early before the two expensive KD-tree calls below
            if (angle_v2 > 30) {
                // doc sbnd_xin/docs/pr/136 round 2
                // (shower_pass4_prefilter_v1_escape): this cheap filter is
                // STRICTER than the acceptance disjunction it guards -- clauses
                // 1 and 3 below test angle_v1 against pair_dis and never look at
                // angle_v2.  The escape re-tests exactly those two clauses, and
                // only on the halves whose quantities are already in hand
                // (angle_v1, pair_dis); the close_shower_dis halves would need
                // the KD-tree call this filter exists to avoid, so they stay out
                // (doc pr/136 sec 10.3 -- that is also why the measured ceiling
                // is a floor).  Knob off => the disjunction is never evaluated
                // and every candidate takes the legacy `continue` => the pass is
                // byte-identical, tape included.
                const bool p4_v1_escape =
                    m_shower_pass4_prefilter_v1_escape &&
                    (m_shower_pass4_prefilter_v1_max_v2 <= 0 ||
                     angle_v2 < m_shower_pass4_prefilter_v1_max_v2) &&
                    (m_shower_pass4_prefilter_v1_max_dis <= 0 ||
                     pair_dis < m_shower_pass4_prefilter_v1_max_dis) &&
                    ((angle_v1 < 25 && pair_dis < 80 * units::cm) ||
                     (angle_v1 < 12.5 && pair_dis < 120 * units::cm));
                if (!p4_v1_escape) {
                    // doc pr/130 part 10: this shower looked at this segment and
                    // dropped it on the cheap angle gate.  Gate tested first, so
                    // the branch costs one static bool when the env is unset.
                    if (pr130_xclus_dbg()) {
                        std::fprintf(stderr,
                            "SHOWER_XCLUS REJECT site=pass4_angle_early shower=%d seg=%d angle_v1=%.1f angle_v2=%.1f pair_dis_cm=%.1f\n",
                            pr91_seg_display_id(shower->start_segment()),
                            pr91_seg_display_id(seg1), angle_v1, angle_v2,
                            pair_dis / units::cm);
                    }
                    continue;
                }
                // The escape fired: fall through to the full disjunction.  Taped
                // so the ON arm's admissions are attributable segment by segment.
                if (pr130_xclus_dbg()) {
                    std::fprintf(stderr,
                        "SHOWER_XCLUS ESCAPE site=pass4_prefilter_v1 shower=%d seg=%d angle_v1=%.1f angle_v2=%.1f pair_dis_cm=%.1f\n",
                        pr91_seg_display_id(shower->start_segment()),
                        pr91_seg_display_id(seg1), angle_v1, angle_v2,
                        pair_dis / units::cm);
                }
            }

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
                // doc sbnd_xin/docs/pr/117 round 1 (shower_pass4_best_owner):
                // the cone has accepted seg1 -- arbitrate WHICH shower owns
                // it.  Current shower's metric is the better of its two
                // anchors (start vertex / closest-approach point, matching
                // the acceptance disjunction above); each rival competes
                // only if it passes the pass-3 cone disjunction from its own
                // start point, on the same 40cm x 5cm ellipsoid metric.
                // Empty rival list (knob off, no probe) => owner stays the
                // current shower and the tag stays pass4_angle => the two
                // legacy lines below are reached unchanged.
                ShowerPtr p4_owner = shower;
                if (!pass4_rivals.empty()) {
                    auto p4_ellip = [](double dist, double angle_deg) {
                        const double a = angle_deg * M_PI / 180.0;
                        return std::pow(dist * std::cos(a), 2) / std::pow(40 * units::cm, 2) +
                               std::pow(dist * std::sin(a), 2) / std::pow(5 * units::cm, 2);
                    };
                    const double metric_cur = std::min(p4_ellip(pair_dis, angle_v1),
                                                       p4_ellip(v2.magnitude(), angle_v2));
                    double best_metric = metric_cur;
                    ShowerPtr best_shower = shower;
                    for (const auto& r : pass4_rivals) {
                        auto [rdist, rcp] = segment_get_closest_point(seg1, r.start_pt);
                        WireCell::Vector rv1(rcp.x() - r.start_pt.x(), rcp.y() - r.start_pt.y(), rcp.z() - r.start_pt.z());
                        double rangle = std::acos(std::clamp(r.dir.dot(rv1) / (r.dir.magnitude() * rv1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        if (!((rangle < 25.0 && rdist < 80 * units::cm) ||
                              (rangle < 12.5 && rdist < 130 * units::cm) ||
                              (rangle < 5.0 && rdist < 200 * units::cm))) continue;
                        const double m = p4_ellip(rdist, rangle);
                        // doc pr/130 part 10: every rival's own metric, so a
                        // shower that competed and lost is distinguishable from
                        // one that was never in the list at all.
                        if (pr130_xclus_dbg()) {
                            std::fprintf(stderr,
                                "SHOWER_XCLUS RIVAL seg=%d cur=%d rival=%d rival_metric=%.4f rival_dist_cm=%.1f rival_angle=%.1f\n",
                                pr91_seg_display_id(seg1),
                                pr91_seg_display_id(shower->start_segment()),
                                pr91_seg_display_id(r.shower->start_segment()),
                                m, rdist / units::cm, rangle);
                        }
                        if (m < best_metric) { best_metric = m; best_shower = r.shower; }
                    }
                    if (pr93_absorb_dbg()) {
                        std::fprintf(stderr,
                            "SHOWER_ABSORB PASS4_OWNER seg=%d cur=%d cur_metric=%.4f best=%d best_metric=%.4f divert=%d applied=%d\n",
                            pr91_seg_display_id(seg1),
                            pr91_seg_display_id(shower->start_segment()),
                            metric_cur,
                            pr91_seg_display_id(best_shower->start_segment()),
                            best_metric,
                            (int)(best_shower != shower),
                            (int)m_shower_pass4_best_owner);
                    }
                    if (m_shower_pass4_best_owner) p4_owner = best_shower;
                }
                // doc sbnd_xin/docs/pr/123 round 1 (shower_pass4_track_guard_len):
                // the owner's track-like prong.  A long MIP-flat / track-pdg
                // segment riding the cone axis is a muon or hadron track, not
                // EM growth (SBND 171572: 125 cm mu- absorbed at 3.9 deg;
                // 393505: 108 cm mu-); the census shows zero labeled-good
                // absorbs above 50 cm.  Decline the absorb; 0 = off = legacy.
                if (m_shower_pass4_track_guard_len > 0) {
                    const double guard_len = segment_track_length(seg1);
                    if (guard_len > m_shower_pass4_track_guard_len) {
                        const int guard_pdg = seg1->has_particle_info() && seg1->particle_info()
                            ? std::abs(seg1->particle_info()->pdg()) : 0;
                        bool guard_trk = (guard_pdg == 13 || guard_pdg == 211 || guard_pdg == 2212);
                        if (!guard_trk && m_mip_dqdx_median > 0) {
                            guard_trk = segment_median_dQ_dx(seg1) < 1.3 * m_mip_dqdx_median;
                        }
                        if (guard_trk) {
                            // doc pr/123 round 2: mark the freed track so
                            // the (default-OFF) PF/kine orphan passes can
                            // claim it if nothing else does -- SBND
                            // 18255-171572's muon vanished from the PF tree
                            // (cross-cluster + score-100 sentinel, outside
                            // the pr/93 orphan machinery's scope).  Inert
                            // bit unless those knobs are on.
                            seg1->set_flags(SegmentFlags::kPass4GuardFreed);
                            if (pr93_absorb_dbg()) {
                                std::fprintf(stderr,
                                    "SHOWER_ABSORB PASS4_TRACK_GUARD seg=%d pdg=%d len_cm=%.1f declined=1\n",
                                    pr91_seg_display_id(seg1), guard_pdg,
                                    guard_len / units::cm);
                            }
                            SPDLOG_LOGGER_DEBUG(s_log,
                                "pr123 pass4_track_guard: decline seg={} pdg={} len={:.1f}cm",
                                seg1->id(), guard_pdg, guard_len / units::cm);
                            continue;
                        }
                    }
                }
                // doc pr/123: chain-immune gap -- min distance from seg1 to
                // the pass-entry membership (env-gated; -1 when unavailable).
                double pr123_snap_dis = -1.0;
                if (pr93_absorb_dbg()) {
                    for (const auto& sa : pr123_snapshot) {
                        for (const auto& f : sa->fits()) {
                            const double d_ = segment_get_closest_point(seg1, f.point).first;
                            if (d_ >= 0 && (pr123_snap_dis < 0 || d_ < pr123_snap_dis)) pr123_snap_dis = d_;
                        }
                    }
                }
                pr123_probe_pass4_geom(seg1, shower, p4_owner, pair_dis, tmp_shower_dis,
                                       close_shower_dis, pr123_snap_dis,
                                       angle_v1, angle_v2, m_mip_dqdx_median);
                pr93_probe_absorb_direct(p4_owner == shower ? "pass4_angle" : "pass4_angle_divert", p4_owner, seg1);
                p4_owner->add_segment(seg1, true);
            }
            // doc pr/130 part 10: the acceptance disjunction refused it.  All
            // five quantities are already computed, so the tape is free.
            else if (pr130_xclus_dbg()) {
                std::fprintf(stderr,
                    "SHOWER_XCLUS REJECT site=pass4_angle_cone shower=%d seg=%d angle_v1=%.1f angle_v2=%.1f pair_dis_cm=%.1f tmp_dis_cm=%.1f close_dis_cm=%.1f\n",
                    pr91_seg_display_id(shower->start_segment()),
                    pr91_seg_display_id(seg1), angle_v1, angle_v2,
                    pair_dis / units::cm, tmp_shower_dis / units::cm,
                    close_shower_dis / units::cm);
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
            ExaminerPassCounter epb_shower_clustering_with_nv_from_vertices("shower_clustering_with_nv_from_vertices", main_cluster ? main_cluster->get_cluster_id() : -1);   // doc pdvd/26 round 2
            while (flag_continue) {
                if (epb_shower_clustering_with_nv_from_vertices.exceeded()) break;
                flag_continue = false;
                for (auto seg1 : seg_order) {
                    if ((seg1->cluster() == main_cluster && !m_absorb_unreachable_main_segs.count(seg1)) || m_nv_bridge_shield_segs.count(seg1)) continue;  // doc pr/65 round 3: guard means "claimed by the main_vertex graph walk", so graph-unreachable main-cluster segments stay eligible (set empty when knob off => legacy); doc pr/40 round 9 B2: bridged-cluster segments stay un-absorbable (shield set empty when bridge off)
                    if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;

                    // doc sbnd_xin/docs/pr/130 item 1b (shower_pass4_prox_guard_len):
                    // pass4_proximity is the only pass-4 absorber with no
                    // track guard.  Its sibling pass4_angle (~50 lines above)
                    // declines a long track-pdg / MIP-flat segment via
                    // pr/123's shower_pass4_track_guard_len; this loop tests
                    // only distance (3.5 cm) and angle (45 deg), and the
                    // 0.75 x shower_length ceiling does not bite on a large
                    // shower.  The 239-event census found 4 segments that
                    // pr/93's cone_absorb_guard had ALREADY declined being
                    // re-admitted here: SBND 100222 seg 14003 (110 cm pdg-13,
                    // 34.5% of the labelled over-clustering charge on the
                    // 141-set) and 176502 segs 20008 / 20013 / 109141.  Same
                    // predicate and same kPass4GuardFreed release as pr/123's
                    // guard, so a declined track stays claimable by the PF /
                    // kine orphan passes.  0 = off = legacy => byte-identical.
                    if (m_shower_pass4_prox_guard_len > 0) {
                        const double prox_len = segment_track_length(seg1);
                        if (prox_len > m_shower_pass4_prox_guard_len) {
                            const int prox_pdg = seg1->has_particle_info() && seg1->particle_info()
                                ? std::abs(seg1->particle_info()->pdg()) : 0;
                            bool prox_trk = (prox_pdg == 13 || prox_pdg == 211 || prox_pdg == 2212);
                            if (!prox_trk && m_mip_dqdx_median > 0) {
                                prox_trk = segment_median_dQ_dx(seg1) < 1.3 * m_mip_dqdx_median;
                            }
                            if (prox_trk) {
                                seg1->set_flags(SegmentFlags::kPass4GuardFreed);
                                if (pr93_absorb_dbg()) {
                                    std::fprintf(stderr,
                                        "SHOWER_ABSORB PASS4_PROX_GUARD seg=%d pdg=%d len_cm=%.1f declined=1\n",
                                        pr91_seg_display_id(seg1), prox_pdg, prox_len / units::cm);
                                }
                                SPDLOG_LOGGER_DEBUG(s_log,
                                    "pr130 pass4_prox_guard: decline seg={} pdg={} len={:.1f}cm",
                                    pr91_seg_display_id(seg1), prox_pdg, prox_len / units::cm);
                                continue;
                            }
                        }
                    }

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


// doc sbnd_xin/docs/pr/118 -- charge-continuity connector between two showers.
// The pr/117 measurement showed blind gap cannot separate an under-clustered
// stub from an over-clustered one (winner/loser gap distributions overlap
// completely at 2.0-5.9 cm); the discriminating feature is whether CHARGE is
// continuous along the connector.  Endpoints come from an exact deterministic
// argmin over the fragment's member fit points against the absorber's fit
// cloud (shower_get_closest_dis is a seeded ping-pong approximation -- fine
// for the frozen legacy gates, wrong for walk endpoints).  The walk recipe is
// forked from NeutrinoGraphAudit.cxx straight_steiner_chain (itself from
// NeutrinoStructureExaminer.cxx:344-439): 0.6 cm steps, per sample
// dv->contained_by then transform->backward on the fragment cluster's t0 then
// is_good_point(0.2 cm, 0, 0); a sample outside every TPC is bad.  No early
// exit -- callers want the full distribution.  Shared by the byte-neutral
// pr/118 probe and (later) the knob path, so probe and knob measure ONE
// implementation.
namespace {
    struct Pr118Connector {
        double gap_exact = -1;             // exact frag-fit-point -> absorber-cloud distance
        WireCell::Point p_frag, p_abs;     // walk endpoints (fragment side, absorber side)
        SegmentPtr frag_seg;               // fragment member holding p_frag
        int nstep = 0, ngood = 0, nbad = 0, badrun = 0;
        double cont_frac = -1, qmed = -1, qfrac = -1;
        bool walked = false;
    };
}

static bool pr118_connector_endpoints(Shower& absorber, Shower& fragment, Graph& graph, Pr118Connector& out)
{
    bool found = false;
    for (auto e : ordered_edges(fragment, graph)) {
        SegmentPtr sb = graph[e].segment;
        if (!sb) continue;
        for (const auto& fit : sb->fits()) {
            if (!fit.valid()) continue;
            auto [d, p] = shower_get_closest_point(absorber, fit.point);
            if (d < 0) continue;           // absorber cloud missing/empty
            if (!found || d < out.gap_exact) {
                found = true;
                out.gap_exact = d;
                out.p_frag = fit.point;
                out.p_abs = p;
                out.frag_seg = sb;
            }
        }
    }
    return found;
}

static void pr118_connector_walk(Pr118Connector& c, TrackFitting& track_fitter, WireCell::IDetectorVolumes::pointer dv)
{
    c.walked = false;
    if (!c.frag_seg) return;
    auto* cluster = c.frag_seg->cluster();
    if (!cluster) return;
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster->get_scope_transform(cluster->get_default_scope()));
    auto grouping = cluster->grouping();
    if (!transform || !grouping) return;
    const double t0 = cluster->get_cluster_t0();

    const double step = 0.6 * WireCell::units::cm;
    const double dx = c.p_abs.x() - c.p_frag.x();
    const double dy = c.p_abs.y() - c.p_frag.y();
    const double dz = c.p_abs.z() - c.p_frag.z();
    const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    const int ncount = (int) std::round(dist / step);
    int run = 0;
    std::vector<double> qs;
    for (int i = 1; i < ncount; ++i) {
        WireCell::Point tp(c.p_frag.x() + dx / ncount * i,
                           c.p_frag.y() + dy / ncount * i,
                           c.p_frag.z() + dz / ncount * i);
        ++c.nstep;
        auto wpid = dv->contained_by(tp);
        bool good = false;
        if (wpid.face() != -1 && wpid.apa() != -1) {
            auto p_raw = transform->backward(tp, t0, wpid.face(), wpid.apa());
            good = grouping->is_good_point(p_raw, wpid.apa(), wpid.face(), 0.2 * WireCell::units::cm, 0, 0);
            qs.push_back(grouping->get_ave_3d_charge(p_raw, wpid.apa(), wpid.face()));
        }
        if (good) { ++c.ngood; run = 0; }
        else { ++c.nbad; c.badrun = std::max(c.badrun, ++run); }
    }
    // Zero interior samples (endpoints within one step) = touching showers:
    // continuous by construction.
    c.cont_frac = c.nstep > 0 ? (double) c.ngood / c.nstep : 1.0;
    if (!qs.empty()) {
        std::sort(qs.begin(), qs.end());
        c.qmed = qs[qs.size() / 2];
        int npos = 0;
        for (double q : qs) if (q > 0) ++npos;
        c.qfrac = (double) npos / qs.size();
    }
    c.walked = true;
}

// doc sbnd_xin/docs/pr/118 -- byte-neutral EM-EM shower-pair census under
// WCT_SHOWER_MERGE_DEBUG: for every fragment/absorber pair within reach
// (45 cm prefilter, gap_exact <= 30 cm -- covers the pr/115 orphan-stub body
// distance IQR upper quartile 29.7 cm) print gap (production recipe),
// gap_exact, the continuity walk, line charge, junction dQ/dx and the
// gamma-gamma guard inputs (mv1/mv2), with NO min-len floor and NO guard cut
// so offline analysis can apply any candidate predicate.  This is the
// measurement that pins the Phase-B continuity thresholds.
static void pr118_probe_continuity_pairs(const std::vector<ShowerPtr>& sorted_showers,
                                         std::map<ShowerPtr, double>& map_shower_length,
                                         VertexPtr main_vertex, Graph& graph,
                                         TrackFitting& track_fitter, WireCell::IDetectorVolumes::pointer dv)
{
    for (auto& shower2 : sorted_showers) {                       // fragment candidate
        if (std::abs(shower2->get_particle_type()) != 11) continue;
        const double len2 = map_shower_length[shower2];
        auto [sv2, t2] = shower2->get_start_vertex_and_type();
        const WireCell::Point p2 = shower2->get_start_point();
        for (auto& shower1 : sorted_showers) {                   // absorber candidate
            if (shower1 == shower2) continue;
            if (std::abs(shower1->get_particle_type()) != 11) continue;
            const double len1 = map_shower_length[shower1];
            if (!(len2 < len1)) continue;                        // fragment strictly shorter
            // cheap prefilter before the per-fit-point scan (a -1 "no cloud"
            // return passes here and is dropped by the endpoint argmin)
            if (shower_get_closest_point(*shower1, p2).first > 45 * WireCell::units::cm) continue;
            auto [sv1, t1] = shower1->get_start_vertex_and_type();

            Pr118Connector c;
            if (!pr118_connector_endpoints(*shower1, *shower2, graph, c)) continue;
            if (c.gap_exact > 30 * WireCell::units::cm) continue;

            // production gap recipe (merge_shower_fragments), for comparison
            double gap = 1e9;
            for (auto e : ordered_edges(*shower2, graph)) {
                SegmentPtr sb = graph[e].segment;
                if (!sb) continue;
                gap = std::min(gap, shower_get_closest_dis(*shower1, sb));
            }
            // local-pivot axis-folded angle (merge_shower_fragments recipe),
            // computed for every probed pair; -1 = unmeasurable direction
            double angle_fold = -1;
            {
                auto [d1_, p1] = shower_get_closest_point(*shower1, p2);
                WireCell::Vector dir1 = shower_cal_dir_3vector(*shower1, p1, 30 * WireCell::units::cm);
                WireCell::Vector dir2 = shower_cal_dir_3vector(*shower2, p2, 30 * WireCell::units::cm);
                if (dir1.magnitude() > 0.001 && dir2.magnitude() > 0.001) {
                    const double cf = std::abs(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()));
                    angle_fold = std::acos(std::clamp(cf, 0.0, 1.0)) / M_PI * 180.0;
                }
                (void) d1_;
            }
            pr118_connector_walk(c, track_fitter, dv);

            const double dqdx_frag = segment_median_dQ_dx(c.frag_seg);
            double dqdx_abs = -1;
            {   // absorber member nearest the junction, deterministic argmin
                double best = -1;
                for (auto e : ordered_edges(*shower1, graph)) {
                    SegmentPtr sa = graph[e].segment;
                    if (!sa) continue;
                    const double d = segment_get_closest_point(sa, c.p_frag).first;
                    if (d < 0) continue;
                    if (best < 0 || d < best) { best = d; dqdx_abs = segment_median_dQ_dx(sa); }
                }
            }
            // Absorber-AXIS geometry (pr/118 second measurement pass): the
            // local fold and the charge walk both measured non-separating;
            // the axis cone is the remaining candidate discriminator.  Angle
            // and distance of the fragment (start point p2 and junction
            // p_frag) seen from the absorber START along its 15 cm and
            // 100 cm axes, unfolded 0-180 deg; -1 when the axis direction is
            // degenerate.
            double ax15_ang = -1, ax100_ang = -1, ax_d = -1;
            double jx15_ang = -1, jx100_ang = -1, jx_d = -1;
            {
                const WireCell::Point s1 = shower1->get_start_point();
                WireCell::Vector a15 = shower_cal_dir_3vector(*shower1, s1, 30 * WireCell::units::cm);
                WireCell::Vector a100 = shower_cal_dir_3vector(*shower1, s1, 100 * WireCell::units::cm);
                auto ax_angle = [](const WireCell::Vector& ax, const WireCell::Vector& v) {
                    if (ax.magnitude() < 0.001 || v.magnitude() < 0.001) return -1.0;
                    return std::acos(std::clamp(ax.dot(v) / (ax.magnitude() * v.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                };
                WireCell::Vector v2(p2.x() - s1.x(), p2.y() - s1.y(), p2.z() - s1.z());
                WireCell::Vector vj(c.p_frag.x() - s1.x(), c.p_frag.y() - s1.y(), c.p_frag.z() - s1.z());
                ax_d = v2.magnitude();
                jx_d = vj.magnitude();
                ax15_ang = ax_angle(a15, v2);
                ax100_ang = ax_angle(a100, v2);
                jx15_ang = ax_angle(a15, vj);
                jx100_ang = ax_angle(a100, vj);
            }
            const auto* frag_cl = c.frag_seg->cluster();
            const auto* abs_cl = shower1->start_segment() ? shower1->start_segment()->cluster() : nullptr;
            std::fprintf(stderr, "SHOWER_MERGE tag=cont_probe keep_sid=%d keep_node=%d "
                                 "cand_sid=%d cand_node=%d conn1=%d conn2=%d mv1=%d mv2=%d "
                                 "len1=%.1fcm len2=%.1fcm gap=%.2fcm gap_exact=%.2fcm "
                                 "angle_fold=%.2f nstep=%d ngood=%d nbad=%d badrun=%d "
                                 "cont_frac=%.3f qmed=%.1f qfrac=%.3f dqdx_frag=%.1f "
                                 "dqdx_abs=%.1f t0_frag=%.1f t0_abs=%.1f walked=%d "
                                 "ax15_ang=%.2f ax100_ang=%.2f ax_d=%.2fcm "
                                 "jx15_ang=%.2f jx100_ang=%.2f jx_d=%.2fcm\n",
                         shower1->get_shower_id(), pr91_seg_display_id(shower1->start_segment()),
                         shower2->get_shower_id(), pr91_seg_display_id(shower2->start_segment()),
                         t1, t2,
                         (int) (sv1 == main_vertex && t1 <= 2),
                         (int) (sv2 == main_vertex && t2 <= 2),
                         len1 / WireCell::units::cm, len2 / WireCell::units::cm,
                         gap / WireCell::units::cm, c.gap_exact / WireCell::units::cm,
                         angle_fold, c.nstep, c.ngood, c.nbad, c.badrun,
                         c.cont_frac, c.qmed, c.qfrac, dqdx_frag, dqdx_abs,
                         frag_cl ? frag_cl->get_cluster_t0() / WireCell::units::us : -1.0,
                         abs_cl ? abs_cl->get_cluster_t0() / WireCell::units::us : -1.0,
                         (int) c.walked,
                         ax15_ang, ax100_ang, ax_d / WireCell::units::cm,
                         jx15_ang, jx100_ang, jx_d / WireCell::units::cm);
        }
    }
}

// doc sbnd_xin/docs/pr/119 -- byte-neutral shower-membership census under
// WCT_SHOWER_EXPEL_DEBUG.  The 29 hand-scan OUT marks (segments wrongly HELD
// by a shower) sit in 6 events, and on every one of them the wanted and
// unwanted sides live in DIFFERENT imaging clusters (pr/115 sec 16.4: "what
// those 6 need is a guard on the cross-cluster absorb"); 5 of the 6 have the
// shower's own root on the wrong side, so the unit of measurement is not the
// member but the (imaging-cluster x view-connected-component) GROUP, anchored
// at the component that contains the start segment.  The probe partitions
// every shower that way and prints per-group geometry/charge/junction
// features so the offline census (scripts/pr119_expel_census.py) can test
// the expel predicate before any behavior exists.  Log/stderr only: no
// effect on emitted bytes.
static inline bool pr119_expel_dbg()
{
    static const bool dbg = std::getenv("WCT_SHOWER_EXPEL_DEBUG") != nullptr;
    return dbg;
}

namespace {
    struct Pr119Group {
        int cluster_id = -1;
        bool anchor = false;               // holds the shower's start segment
        std::vector<SegmentPtr> segs;      // ascending segment graph index
    };
}

// Partition a shower's members into (imaging-cluster, view-connected
// component) groups: anchor first (the start segment's component restricted
// to the start segment's cluster), then remaining components seeded in
// ascending segment graph-index order, each restricted to its own cluster.
// Component membership is traversal-order independent; adjacency runs over
// stable vertex indices, never pointer order.  Shared by the Phase-A probe
// and the (future) Phase-B expel pass so the measurement is the mechanism.
static std::vector<Pr119Group> pr119_partition(Shower& sh, Graph& graph)
{
    std::vector<Pr119Group> out;
    std::vector<SegmentPtr> members;
    for (auto edesc : ordered_edges(sh, graph)) {
        SegmentPtr seg = graph[edesc].segment;
        if (!seg || !seg->descriptor_valid()) continue;
        members.push_back(seg);
    }
    std::sort(members.begin(), members.end(), SegmentIndexCmp{});
    if (members.empty()) return out;

    // vertex stable index -> member positions (the adjacency relation)
    std::map<size_t, std::vector<size_t>> vtx_to_members;
    std::vector<std::array<long long, 2>> member_vtx(members.size(), {-1, -1});
    for (size_t i = 0; i < members.size(); ++i) {
        auto [va, vb] = PR::find_vertices(graph, members[i]);
        int k = 0;
        for (VertexPtr v : {va, vb}) {
            if (!v || !v->descriptor_valid()) continue;
            const size_t idx = graph[v->get_descriptor()].index;
            vtx_to_members[idx].push_back(i);
            member_vtx[i][k++] = (long long) idx;
        }
    }
    auto cluster_of = [](const SegmentPtr& s) {
        return s->cluster() ? s->cluster()->get_cluster_id() : -1;
    };
    std::vector<char> visited(members.size(), 0);
    auto flood = [&](size_t seed, int cid, Pr119Group& g) {
        std::vector<size_t> stack{seed};
        visited[seed] = 1;
        while (!stack.empty()) {
            const size_t cur = stack.back();
            stack.pop_back();
            g.segs.push_back(members[cur]);
            for (long long vi : member_vtx[cur]) {
                if (vi < 0) continue;
                for (size_t nb : vtx_to_members[(size_t) vi]) {
                    if (visited[nb]) continue;
                    if (cluster_of(members[nb]) != cid) continue;
                    visited[nb] = 1;
                    stack.push_back(nb);
                }
            }
        }
        std::sort(g.segs.begin(), g.segs.end(), SegmentIndexCmp{});
    };

    SegmentPtr root = sh.start_segment();
    size_t root_pos = members.size();
    for (size_t i = 0; i < members.size(); ++i) {
        if (members[i] == root) { root_pos = i; break; }
    }
    if (root_pos < members.size()) {
        Pr119Group g;
        g.cluster_id = cluster_of(root);
        g.anchor = true;
        flood(root_pos, g.cluster_id, g);
        out.push_back(std::move(g));
    }
    for (size_t i = 0; i < members.size(); ++i) {
        if (visited[i]) continue;
        Pr119Group g;
        g.cluster_id = cluster_of(members[i]);
        flood(i, g.cluster_id, g);
        out.push_back(std::move(g));
    }
    return out;
}

// Exact group->anchor junction: min distance over both sides' valid fit
// points (deterministic double argmin, strict <).  Fills the Pr118Connector
// endpoint fields so pr118_connector_walk can run unchanged on the result.
static bool pr119_group_endpoints(const std::vector<SegmentPtr>& grp,
                                  const std::vector<SegmentPtr>& anchor,
                                  Pr118Connector& out)
{
    bool found = false;
    for (const auto& sb : grp) {
        for (const auto& fb : sb->fits()) {
            if (!fb.valid()) continue;
            for (const auto& sa : anchor) {
                for (const auto& fa : sa->fits()) {
                    if (!fa.valid()) continue;
                    const double d = (fb.point - fa.point).magnitude();
                    if (!found || d < out.gap_exact) {
                        found = true;
                        out.gap_exact = d;
                        out.p_frag = fb.point;
                        out.p_abs = fa.point;
                        out.frag_seg = sb;
                    }
                }
            }
        }
    }
    return found;
}

static void pr119_probe_expel_groups(IndexedShowerSet& showers, Graph& graph,
                                     VertexPtr main_vertex,
                                     TrackFitting& track_fitter,
                                     WireCell::IDetectorVolumes::pointer dv,
                                     double mip_dqdx_median)
{
    if (!pr119_expel_dbg()) return;
    const WireCell::Point mvp = (main_vertex && main_vertex->fit().valid())
        ? main_vertex->fit().point
        : (main_vertex ? main_vertex->wcpt().point : WireCell::Point(0, 0, 0));

    for (auto& sh : showers) {                       // IndexedShowerSet order
        if (!sh) continue;
        auto [svtx, conn] = sh->get_start_vertex_and_type();
        auto groups = pr119_partition(*sh, graph);
        if (groups.empty()) continue;

        double tot_dQ = 0, tot_len = 0;
        std::set<int> clusters;
        bool sh_touch_mv = false;
        for (const auto& g : groups) {
            clusters.insert(g.cluster_id);
            for (const auto& seg : g.segs) {
                tot_len += segment_track_length(seg);
                for (const auto& f : seg->fits()) tot_dQ += f.dQ;
                auto [va, vb] = PR::find_vertices(graph, seg);
                if (va == main_vertex || vb == main_vertex) sh_touch_mv = true;
            }
        }
        std::fprintf(stderr,
                     "EXPEL_SHOWER shower_id=%d node_id=%d conn=%d pdg=%d nseg=%d ncls=%zu ngrp=%zu "
                     "root_cluster=%d kine_charge=%.3f tot_len=%.2fcm tot_dQ=%.1f touches_main_vtx=%d\n",
                     sh->get_shower_id(), pr91_seg_display_id(sh->start_segment()), conn,
                     sh->get_particle_type(), sh->get_num_segments(), clusters.size(), groups.size(),
                     groups.front().anchor ? groups.front().cluster_id : -1,
                     sh->get_kine_charge() / WireCell::units::MeV,
                     tot_len / WireCell::units::cm, tot_dQ, (int) sh_touch_mv);
        if (groups.size() < 2) continue;             // single-group showers: header only

        const bool have_anchor = groups.front().anchor;
        const std::vector<SegmentPtr>& anchor_segs = groups.front().segs;
        const WireCell::Point sp = sh->get_start_point();
        const WireCell::Vector a30 = shower_cal_dir_3vector(*sh, sp, 30 * WireCell::units::cm);
        const WireCell::Vector a100 = shower_cal_dir_3vector(*sh, sp, 100 * WireCell::units::cm);
        auto ax_angle = [](const WireCell::Vector& ax, const WireCell::Vector& v) {
            if (ax.magnitude() < 0.001 || v.magnitude() < 0.001) return -1.0;
            return std::acos(std::clamp(ax.dot(v) / (ax.magnitude() * v.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
        };

        for (size_t gi = 0; gi < groups.size(); ++gi) {
            const auto& g = groups[gi];
            double g_dQ = 0, g_len = 0, g_maxlen = 0;
            int n_track_pid = 0;
            bool g_touch_mv = false;
            std::vector<double> dqdxs;
            std::set<size_t> g_vtx_idx;
            for (const auto& seg : g.segs) {
                const double sl = segment_track_length(seg);
                g_len += sl;
                g_maxlen = std::max(g_maxlen, sl);
                for (const auto& f : seg->fits()) {
                    g_dQ += f.dQ;
                    if (f.valid() && f.dx > 0) dqdxs.push_back(f.dQ / f.dx);
                }
                if (segment_confident_nonelectron_pid(seg)) ++n_track_pid;
                auto [va, vb] = PR::find_vertices(graph, seg);
                for (VertexPtr v : {va, vb}) {
                    if (!v || !v->descriptor_valid()) continue;
                    g_vtx_idx.insert(graph[v->get_descriptor()].index);
                    if (v == main_vertex) g_touch_mv = true;
                }
            }
            double med_dqdx_mip = -1;
            if (!dqdxs.empty() && mip_dqdx_median > 0) {
                std::sort(dqdxs.begin(), dqdxs.end());
                med_dqdx_mip = dqdxs[dqdxs.size() / 2] / mip_dqdx_median;
            }
            // how many showers hold at least one member of this group (>1 =
            // 47212-style double ownership: expel but do not spawn)
            int nsh_holding = 0;
            for (auto& other : showers) {
                if (!other) continue;
                bool holds = false;
                for (const auto& seg : g.segs) {
                    if (seg->descriptor_valid() && other->has_edge(seg->get_descriptor())) { holds = true; break; }
                }
                if (holds) ++nsh_holding;
            }
            // links to the rest of the shower: distinct vertices shared with
            // non-group members (0 = the group floats free of the view)
            std::set<size_t> other_vtx_idx;
            for (size_t oj = 0; oj < groups.size(); ++oj) {
                if (oj == gi) continue;
                for (const auto& seg : groups[oj].segs) {
                    auto [va, vb] = PR::find_vertices(graph, seg);
                    for (VertexPtr v : {va, vb}) {
                        if (v && v->descriptor_valid()) other_vtx_idx.insert(graph[v->get_descriptor()].index);
                    }
                }
            }
            int nlinks = 0;
            for (size_t idx : g_vtx_idx) {
                if (other_vtx_idx.count(idx)) ++nlinks;
            }

            // junction to the ANCHOR group + continuity walk + axis geometry
            Pr118Connector c;
            double ax30_ang = -1, ax100_ang = -1, dis_start = -1, grp_ang = -1;
            bool have_jx = false;
            if (!g.anchor && have_anchor) {
                have_jx = pr119_group_endpoints(g.segs, anchor_segs, c);
                if (have_jx) {
                    pr118_connector_walk(c, track_fitter, dv);
                    const WireCell::Vector vj(c.p_frag.x() - sp.x(), c.p_frag.y() - sp.y(),
                                              c.p_frag.z() - sp.z());
                    dis_start = vj.magnitude();
                    ax30_ang = ax_angle(a30, vj);
                    ax100_ang = ax_angle(a100, vj);
                    // group development direction: junction -> farthest group
                    // fit point, vs the shower axis
                    double far_d = -1;
                    WireCell::Point far_p = c.p_frag;
                    for (const auto& seg : g.segs) {
                        for (const auto& f : seg->fits()) {
                            if (!f.valid()) continue;
                            const double d = (f.point - c.p_frag).magnitude();
                            if (d > far_d) { far_d = d; far_p = f.point; }
                        }
                    }
                    const WireCell::Vector vg(far_p.x() - c.p_frag.x(), far_p.y() - c.p_frag.y(),
                                              far_p.z() - c.p_frag.z());
                    grp_ang = ax_angle(a30, vg);
                }
            }
            // distance to the main vertex (min over group fit points)
            double dis_main = -1;
            for (const auto& seg : g.segs) {
                for (const auto& f : seg->fits()) {
                    if (!f.valid()) continue;
                    const double d = (f.point - mvp).magnitude();
                    if (dis_main < 0 || d < dis_main) dis_main = d;
                }
            }
            // spawn-anchoring dry run (the conn3_unreachable recipe): nearest
            // non-group vertex of the shower view (or the main vertex) to the
            // group's longest member, 0.8*main_dis preference, <=80cm conn 3.
            int anchor_vtx = -1, conn_new = -1;
            double anchor_dis = -1;
            if (!g.anchor) {
                SegmentPtr seed = nullptr;
                for (const auto& seg : g.segs) {
                    if (!seed || segment_track_length(seg) > segment_track_length(seed)) seed = seg;
                }
                double min_dis = 1e9, main_dis = 1e9;
                VertexPtr min_vertex = nullptr;
                for (auto vdesc : ordered_nodes(*sh, graph)) {
                    VertexPtr v = graph[vdesc].vertex;
                    if (!v || !v->descriptor_valid()) continue;
                    if (g_vtx_idx.count(graph[v->get_descriptor()].index)) continue;  // never anchor to your own component
                    const WireCell::Point vp = v->fit().valid() ? v->fit().point : v->wcpt().point;
                    const double d = segment_get_closest_point(seed, vp).first;
                    if (d < 0) continue;
                    if (d < min_dis) { min_dis = d; min_vertex = v; }
                    if (v == main_vertex) main_dis = d;
                }
                if (main_vertex && main_vertex->descriptor_valid()
                    && !g_vtx_idx.count(graph[main_vertex->get_descriptor()].index)) {
                    const double d = segment_get_closest_point(seed, mvp).first;
                    if (d >= 0 && d < main_dis) main_dis = d;
                    if (d >= 0 && d < min_dis) { min_dis = d; min_vertex = main_vertex; }
                }
                if (min_vertex && min_dis > 0.8 * main_dis && main_dis < 1e8) {
                    min_dis = main_dis;
                    min_vertex = main_vertex;
                }
                if (min_vertex) {
                    anchor_vtx = pr91_vtx_display_id(min_vertex);
                    anchor_dis = min_dis;
                    conn_new = (min_dis > 80 * WireCell::units::cm) ? 4 : 3;
                }
            }

            std::fprintf(stderr,
                         "EXPEL_GROUP shower_id=%d grp=%zu cluster=%d anchor=%d nseg=%zu len=%.2fcm "
                         "dQ=%.1f qfrac=%.4f med_dqdx_mip=%.3f max_seglen=%.2fcm n_track_pid=%d "
                         "nsh_holding=%d touches_main_vtx=%d nlinks=%d gap_exact=%.2fcm "
                         "jx=(%.3f,%.3f,%.3f) nstep=%d cont_frac=%.3f conn_qmed=%.1f conn_qfrac=%.3f "
                         "walked=%d ax30_ang=%.2f ax100_ang=%.2f grp_ang=%.2f dis_start=%.2fcm "
                         "dis_main=%.2fcm anchor_vtx=%d anchor_dis=%.2fcm conn_new=%d t0_grp=%.1f\n",
                         sh->get_shower_id(), gi, g.cluster_id, (int) g.anchor, g.segs.size(),
                         g_len / WireCell::units::cm, g_dQ, tot_dQ > 0 ? g_dQ / tot_dQ : -1.0,
                         med_dqdx_mip, g_maxlen / WireCell::units::cm, n_track_pid,
                         nsh_holding, (int) g_touch_mv, nlinks,
                         have_jx ? c.gap_exact / WireCell::units::cm : -1.0,
                         have_jx ? c.p_frag.x() / WireCell::units::cm : 0.0,
                         have_jx ? c.p_frag.y() / WireCell::units::cm : 0.0,
                         have_jx ? c.p_frag.z() / WireCell::units::cm : 0.0,
                         c.nstep, c.cont_frac, c.qmed, c.qfrac, (int) c.walked,
                         ax30_ang, ax100_ang, grp_ang,
                         dis_start >= 0 ? dis_start / WireCell::units::cm : -1.0,
                         dis_main >= 0 ? dis_main / WireCell::units::cm : -1.0,
                         anchor_vtx, anchor_dis >= 0 ? anchor_dis / WireCell::units::cm : -1.0,
                         conn_new,
                         (g.segs.front()->cluster())
                             ? g.segs.front()->cluster()->get_cluster_t0() / WireCell::units::us : -1.0);

            for (const auto& seg : g.segs) {
                double sdQ = 0, sdx = 0;
                for (const auto& f : seg->fits()) { sdQ += f.dQ; sdx += f.dx; }
                auto [v0, v1] = PR::find_vertices(graph, seg);
                std::fprintf(stderr,
                             "EXPEL_MEMBER shower_id=%d seg=%d grp=%zu cluster=%d len=%.3fcm dQ=%.1f "
                             "dQdx=%.1f pdg=%d dirsign=%d v0=%d v1=%d\n",
                             sh->get_shower_id(), pr91_seg_display_id(seg), gi, g.cluster_id,
                             segment_track_length(seg) / WireCell::units::cm, sdQ,
                             sdx > 0 ? sdQ / (sdx / WireCell::units::cm) : 0.0,
                             seg->has_particle_info() && seg->particle_info()
                                 ? seg->particle_info()->pdg() : 0,
                             seg->dirsign(), pr91_vtx_display_id(v0), pr91_vtx_display_id(v1));
            }
        }
    }
}

// doc sbnd_xin/docs/pr/117 round 1 -- shower_merge_relax.  Rationale at the
// m_shower_merge_relax member block (NeutrinoPatternBase.h).  A late fragment
// consolidation pass, forked in SHAPE from examine_merge_showers' two-phase
// plan/execute structure (that production pass stays byte-untouched, M10) but
// with proximity + local-pivot geometry instead of the main-vertex 10-deg
// test.  Called only when the knob is on.
void PatternAlgorithms::merge_shower_fragments(Graph& graph, IndexedShowerSet& showers, VertexPtr main_vertex, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    if (showers.size() < 2) return;

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

    std::map<ShowerPtr, double> map_shower_length;   // lookup only, never iterated
    for (auto& s : sorted_showers) map_shower_length[s] = s->get_total_length();

    // doc sbnd_xin/docs/pr/118 -- byte-neutral pair census under
    // WCT_SHOWER_MERGE_DEBUG only; the production loops below are untouched.
    if (pr91_merge_dbg())
        pr118_probe_continuity_pairs(sorted_showers, map_shower_length, main_vertex, graph,
                                     track_fitter, dv);

    // Phase 1: plan merges in deterministic order, FRAGMENT-FIRST -- each
    // fragment picks its best (min-gap) absorber among the strictly-bigger
    // showers that pass the gates, rather than the first bigger shower in
    // iteration order claiming it.  'claimed' marks fragments already
    // promised; a shower claimed as a fragment cannot itself absorb.
    std::unordered_set<Shower*> claimed;   // membership-tested only, never iterated
    std::map<ShowerPtr, std::vector<ShowerPtr>> plan_by_absorber;  // lookup + ordered exec below

    for (auto& shower2 : sorted_showers) {                       // fragment candidate
        // EM <-> EM only.  Measured (doc pr/117 sec 7): without this guard
        // the two spurious directional merges in the 98-event set each had a
        // proton-typed side (37112: pdg-2212 fragment into an e-; 389538:
        // e- fragment into a pdg-2212 absorber on a GOOD event).  It also
        // keeps this pass out of the track pool entirely.
        if (std::abs(shower2->get_particle_type()) != 11) continue;
        if (plan_by_absorber.count(shower2)) continue;           // already absorbing: no chains
        const double len2 = map_shower_length[shower2];
        // Fragment length FLOOR.  A fragment too short to carry a 30cm
        // direction (the pr/115 'orphan' class measured: 40 of 41 are pdg-11
        // one-segment showers, median 1.7 cm) is not a candidate at all:
        // merging such stubs on gap alone measured net-NEGATIVE on the
        // marked set (blind proximity cannot tell an under-clustered stub
        // from an over-clustered one), and letting them into the angle test
        // fires on a random ~3% of noise directions (doc pr/117 sec 7).
        // pr/118: with the continuity knob on, short fragments DO enter the
        // loop -- for the continuity-only path; they never take the legacy
        // gap+angle path (zeroed below).
        if (len2 < m_shower_merge_relax_min_len && !m_shower_merge_relax_continuity) continue;
        auto [sv2, t2] = shower2->get_start_vertex_and_type();
        const WireCell::Point p2 = shower2->get_start_point();

        double best_gap = 1e9;
        ShowerPtr best_absorber = nullptr;
        for (auto& shower1 : sorted_showers) {                   // absorber candidate
            if (shower1 == shower2 || claimed.count(shower1.get())) continue;
            if (std::abs(shower1->get_particle_type()) != 11) continue;   // EM <-> EM only, see above
            const double len1 = map_shower_length[shower1];
            if (!(len2 < len1)) continue;                        // absorber is the strictly bigger object
            auto [sv1, t1] = shower1->get_start_vertex_and_type();
            // HARD gamma-gamma guard (not a knob): two showers both rooted
            // at the main vertex with conn type 1/2 are exactly the pi0
            // gamma-pair topology -- merging them is the legacy pass's
            // 10-deg call, never this pass's.
            if (sv1 == main_vertex && sv2 == main_vertex && t1 <= 2 && t2 <= 2) continue;

            // Body gap: min over the fragment's members of the closest
            // distance to the absorber's fit cloud.
            double gap = 1e9;
            for (auto e : ordered_edges(*shower2, graph)) {
                SegmentPtr sb = graph[e].segment;
                if (!sb) continue;
                gap = std::min(gap, shower_get_closest_dis(*shower1, sb));
            }
            // Local-pivot directions at the meeting point (NOT the main
            // vertex): fragments of one shower meet away from the vertex.
            // Axis-folded agreement -- a fragment CONTINUING the absorber
            // reads anti-parallel here (both dirs point into their own
            // body), so |cos| is the right test.
            double angle_fold = 180.0;
            if (gap < m_shower_merge_relax_dis) {
                auto [d1_, p1] = shower_get_closest_point(*shower1, p2);
                WireCell::Vector dir1 = shower_cal_dir_3vector(*shower1, p1, 30 * units::cm);
                WireCell::Vector dir2 = shower_cal_dir_3vector(*shower2, p2, 30 * units::cm);
                if (dir1.magnitude() > 0.001 && dir2.magnitude() > 0.001) {
                    const double c = std::abs(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()));
                    angle_fold = std::acos(std::clamp(c, 0.0, 1.0)) / M_PI * 180.0;
                }
            }
            bool fires = gap < m_shower_merge_relax_dis && angle_fold < m_shower_merge_relax_angle;
            // doc sbnd_xin/docs/pr/118 -- two-tier axis+charge admission path
            // (rationale and the census behind every threshold at the
            // m_shower_merge_relax_continuity member block).  T1 "touching
            // aligned" (any length): gap_exact <= t1_gap, axis < cont_axis,
            // local fold < t1_fold.  T2 "bright aligned stub" (below the
            // min_len floor only): gap_exact <= cont_gap, axis < cont_axis,
            // junction within cont_dmax, connector charged on every sample
            // with median line charge > cont_qmed.  Knob off => no new
            // computation, byte-identical.
            bool cont_admit = false;
            int cont_tier = 0;
            double gap_exact = -1;
            if (m_shower_merge_relax_continuity) {
                if (len2 < m_shower_merge_relax_min_len) fires = false;   // short fragments: tier paths only
                const double tier_gap_ceiling = std::max(m_shower_merge_relax_cont_gap,
                                                         m_shower_merge_relax_cont_t1_gap);
                if (!fires && shower_get_closest_point(*shower1, p2).first <= 45 * units::cm) {
                    Pr118Connector c;
                    if (pr118_connector_endpoints(*shower1, *shower2, graph, c) &&
                        c.gap_exact <= tier_gap_ceiling) {
                        gap_exact = c.gap_exact;
                        // absorber-axis angle at the junction: min over the
                        // 30 cm and 100 cm start directions (census: true
                        // merges sit at 2-6 deg; false neighbours isotropic)
                        double axang = 1e9;
                        double jx_d = -1;
                        {
                            const WireCell::Point s1 = shower1->get_start_point();
                            WireCell::Vector a15 = shower_cal_dir_3vector(*shower1, s1, 30 * units::cm);
                            WireCell::Vector a100 = shower_cal_dir_3vector(*shower1, s1, 100 * units::cm);
                            WireCell::Vector vj(c.p_frag.x() - s1.x(), c.p_frag.y() - s1.y(),
                                                c.p_frag.z() - s1.z());
                            jx_d = vj.magnitude();
                            if (jx_d > 0.001) {
                                for (const auto& ax : {a15, a100}) {
                                    if (ax.magnitude() < 0.001) continue;
                                    axang = std::min(axang,
                                        std::acos(std::clamp(ax.dot(vj) / (ax.magnitude() * jx_d), -1.0, 1.0)) / M_PI * 180.0);
                                }
                            }
                        }
                        if (axang < m_shower_merge_relax_cont_axis) {
                            if (c.gap_exact <= m_shower_merge_relax_cont_t1_gap) {
                                // T1: touching + aligned + loose local fold
                                double af = angle_fold;
                                if (!(gap < m_shower_merge_relax_dis)) {   // legacy fold not computed above
                                    auto [d1c, p1c] = shower_get_closest_point(*shower1, p2);
                                    WireCell::Vector dc1 = shower_cal_dir_3vector(*shower1, p1c, 30 * units::cm);
                                    WireCell::Vector dc2 = shower_cal_dir_3vector(*shower2, p2, 30 * units::cm);
                                    af = 180.0;
                                    if (dc1.magnitude() > 0.001 && dc2.magnitude() > 0.001) {
                                        const double cf = std::abs(dc1.dot(dc2) / (dc1.magnitude() * dc2.magnitude()));
                                        af = std::acos(std::clamp(cf, 0.0, 1.0)) / M_PI * 180.0;
                                    }
                                    (void) d1c;
                                }
                                if (af < m_shower_merge_relax_cont_t1_fold) {
                                    cont_admit = true;
                                    cont_tier = 1;
                                }
                            }
                            if (!cont_admit && len2 < m_shower_merge_relax_min_len &&
                                c.gap_exact <= m_shower_merge_relax_cont_gap &&
                                jx_d < m_shower_merge_relax_cont_dmax) {
                                // T2: bright aligned stub -- the connector must
                                // carry charge on every sample (a touching pair
                                // has no samples and is T1's case, not T2's)
                                pr118_connector_walk(c, track_fitter, dv);
                                if (c.walked && c.nstep > 0 &&
                                    c.qfrac >= m_shower_merge_relax_cont_frac &&
                                    c.qmed > m_shower_merge_relax_cont_qmed) {
                                    cont_admit = true;
                                    cont_tier = 2;
                                }
                            }
                        }
                    }
                }
            }
            if (pr91_merge_dbg())
                std::fprintf(stderr, "SHOWER_MERGE tag=merge_relax keep_sid=%d keep_node=%d "
                                     "cand_sid=%d cand_node=%d gap=%.2fcm angle_fold=%.2f len2=%.1fcm "
                                     "cut_dis=%.2f cut_angle=%.2f gap_exact=%.2fcm tier=%d verdict=%s\n",
                             shower1->get_shower_id(), pr91_seg_display_id(shower1->start_segment()),
                             shower2->get_shower_id(), pr91_seg_display_id(shower2->start_segment()),
                             gap / units::cm, angle_fold, len2 / units::cm,
                             m_shower_merge_relax_dis / units::cm, m_shower_merge_relax_angle,
                             gap_exact >= 0 ? gap_exact / units::cm : -1.0, cont_tier,
                             fires ? "MERGE"
                                   : cont_admit ? "MERGE_CONT"
                                   : (m_shower_merge_relax_continuity && gap_exact >= 0) ? "cont_fail"
                                   : (gap < m_shower_merge_relax_dis ? "angle_fold_fail" : "gap_fail"));
            // Rank by the admitting metric (gap for the legacy path,
            // gap_exact for the continuity path; the smaller when both).
            double rank = 1e9;
            if (fires) rank = std::min(rank, gap);
            if (cont_admit) rank = std::min(rank, gap_exact);
            if ((fires || cont_admit) && rank < best_gap) {
                best_gap = rank;
                best_absorber = shower1;
            }
        }
        if (best_absorber) {
            plan_by_absorber[best_absorber].push_back(shower2);
            claimed.insert(shower2.get());
        }
    }
    // Rebuild the ordered execution plan (absorbers in sorted order).
    std::vector<std::pair<ShowerPtr, std::vector<ShowerPtr>>> merge_plan;
    for (auto& shower1 : sorted_showers) {
        auto it = plan_by_absorber.find(shower1);
        if (it != plan_by_absorber.end()) merge_plan.emplace_back(shower1, std::move(it->second));
    }

    // Phase 2: execute -- the examine_merge_showers bookkeeping, verbatim.
    for (auto& [shower1, to_merge] : merge_plan) {
        for (auto& shower2 : to_merge) {
            pr93_probe_absorb_splice("merge_relax", shower1, shower2);
            shower1->add_shower(*shower2);
        }
        shower1->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
        shower1->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
        shower1->set_kine_charge(cal_kine_charge(shower1, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv));
        shower1->set_flag_kinematics(true);
    }
    int n_merged = 0;
    for (auto& [shower1, to_merge] : merge_plan) {               // deterministic erase order
        for (auto& shower2 : to_merge) {
            showers.erase(shower2);
            ++n_merged;
        }
    }
    if (n_merged) {
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                           map_vertex_to_shower, used_shower_clusters);
        SPDLOG_LOGGER_DEBUG(s_log,
            "pr117 merge_relax: {} fragment shower(s) absorbed; {} remain", n_merged, showers.size());
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
                // doc pr/130 part 10: the decisive skip.  A segment already owned by
                // ANOTHER shower is dropped here BEFORE any geometry is computed, so the
                // owning shower's claim is invisible to every later shower -- the same
                // silent-skip shape the BLOCKED probe found in the walk one level down.
                if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) {
                    if (pr130_xclus_dbg()) {
                        const auto& owner_sh = map_segment_in_shower.find(seg1)->second;
                        std::fprintf(stderr,
                            "SHOWER_XCLUS OWNED site=in_other_clusters shower=%d seg=%d owner=%d\n",
                            pr91_seg_display_id(shower->start_segment()),
                            pr91_seg_display_id(seg1),
                            owner_sh ? pr91_seg_display_id(owner_sh->start_segment()) : -1);
                    }
                    continue;
                }
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
                bool pr121_rehomed = false;
                for (auto sh_old : showers) {
                    if (sh_old == shower1) continue;
                    if (sh_old->start_segment() != sg) continue;
                    auto [svtx_old, ctype_old] = sh_old->get_start_vertex_and_type();
                    const bool erase_old = (svtx_old == main_vertex && ctype_old == 1);
                    pr121_probe_ex1_dedup(shower1, sh_old, ctype_old,
                                          svtx_old == main_vertex ? 1 : 0, erase_old ? 1 : 0);
                    if (erase_old) {
                        // doc pr/121 r1: a MULTI-segment dedup victim carries
                        // real membership (348471: 13 members, 352.6 MeV --
                        // see pr121_probe_ex1_dedup's comment); dropping it
                        // orphans every member.  Re-home into shower1 first,
                        // the same shape as the pr/84 start-seg dedup
                        // (keep->add_shower(*drop)).  false = legacy drop =
                        // byte-identical.
                        if (m_shower_ex1_dedup_rehome && sh_old->get_num_segments() > 1) {
                            SPDLOG_LOGGER_DEBUG(s_log,
                                "pr121 ex1_dedup_rehome: shower1 seg={} absorbs dedup victim shower_id={} nseg={} kine={:.1f} MeV",
                                pr91_seg_display_id(shower1->start_segment()),
                                sh_old->get_shower_id(), sh_old->get_num_segments(),
                                sh_old->get_kine_charge() / WireCell::units::MeV);
                            pr93_probe_absorb_splice("ex1_dedup_rehome", shower1, sh_old);
                            shower1->add_shower(*sh_old);
                            pr121_rehomed = true;
                        }
                        del_showers.insert(sh_old);
                    }
                }
                // doc pr/121 r1: the kinematics above predate the re-home;
                // recompute once so the absorbed membership is counted.
                if (pr121_rehomed) {
                    shower1->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
                    double kine_charge2 = cal_kine_charge(shower1, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
                    shower1->set_kine_charge(kine_charge2);
                    shower1->set_flag_kinematics(true);
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
                    // doc sbnd_xin/docs/pr/118 -- pr/91 P2: the same gate measured
                    // against the parent's WHOLE fit body (min over its ordered
                    // members; includes the start segment, so body_dis <=
                    // start_dis always -- strictly admissive, the 3 cm threshold
                    // unchanged).  The APPLY knob (shower_ex1_conn3_body_dis) was
                    // retired in doc 77 round 3 -- measured ZERO yield, 1 admit in
                    // 98 events and it failed the angles (pr/118 sec 4a).  What
                    // survives is the PROBE: body_dis is computed and printed under
                    // pr91_merge_dbg() with the would-be downstream angles, so the
                    // admit set and its angle-gate survival stay measurable (the
                    // 394532 lesson).  scripts/pr118_probe_census.py parses it.
                    if (pr91_merge_dbg()) {
                        double body_dis = min_dis;
                        for (auto edesc : ordered_edges(*shower, graph)) {
                            SegmentPtr mseg = graph[edesc].segment;
                            if (!mseg) continue;
                            body_dis = std::min(body_dis, shower_get_closest_dis(*shower1, mseg));
                        }
                        if (pr91_merge_dbg()) {
                        // read-only copy of the downstream recipes just below
                        auto [w_dis, w_pt] = shower_get_closest_point(*shower1, main_vtx_pt);
                        WireCell::Vector w_d2(w_pt.x() - main_vtx_pt.x(),
                                              w_pt.y() - main_vtx_pt.y(),
                                              w_pt.z() - main_vtx_pt.z());
                        double w_angle = std::acos(std::clamp(dir1.dot(w_d2) / (dir1.magnitude() * w_d2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        WireCell::Vector w_d3 = shower_cal_dir_3vector(*shower1, shower1->get_start_point(), 30 * units::cm);
                        double w_angle1 = std::acos(std::clamp(w_d2.dot(w_d3) / (w_d2.magnitude() * w_d3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        const bool legacy_gate_fail = conn_type1 > 2 && min_dis > 3 * units::cm;
                        const bool body_gate_fail   = conn_type1 > 2 && body_dis > 3 * units::cm;
                        const bool angles_pass      = w_angle < 15 && w_angle1 < 15 && body_dis < 28 * units::cm;
                        std::fprintf(stderr, "SHOWER_MERGE tag=ex_shower1_p2dis parent_sid=%d "
                                             "parent_node=%d cand_sid=%d cand_node=%d cand_conn=%d "
                                             "cand_len=%.3f start_dis=%.3f body_dis=%.3f angle=%.3f "
                                             "angle1=%.3f legacy_gate=%s body_gate=%s angles_pass=%d\n",
                                     shower->get_shower_id(),
                                     pr91_seg_display_id(shower->start_segment()),
                                     shower1->get_shower_id(),
                                     pr91_seg_display_id(shower1->start_segment()),
                                     conn_type1, shower1->get_total_length() / units::cm,
                                     min_dis / units::cm, body_dis / units::cm, w_angle, w_angle1,
                                     legacy_gate_fail ? "FAIL" : "PASS",
                                     body_gate_fail ? "FAIL" : "PASS", (int) angles_pass);
                        (void) w_dis;
                        }
                    }
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


// doc pr/132 round 3 -- body of the WCT_PI0_SUBSTRUCT_DEBUG tape (docstring at
// pr132_pi0_substruct_dbg above).  Deterministic: fixed-seed power iteration
// with one deflation; IndexedShowerSet iteration is id-ordered.
static void pr132_pi0_substruct_probe(const IndexedShowerSet& showers,
                                      const IndexedShowerSet& pi0_showers)
{
    for (auto shower : showers) {
        if (!shower) continue;
        const double kq = shower->get_kine_charge();
        if (kq < 300 * WireCell::units::MeV) continue;
        auto pc = shower->get_pcloud("associate_points");
        const size_t n = pc ? pc->npoints() : 0;
        if (n < 20) {
            std::fprintf(stderr, "PI0_SUBSTRUCT sh=%d E=%.1f npts=%zu skip=too_few_points\n",
                         pr132_pi0_shid(shower), kq / WireCell::units::MeV, n);
            continue;
        }
        // Centroid.
        double cx = 0, cy = 0, cz = 0;
        for (size_t i = 0; i < n; ++i) {
            auto p = pc->point3d(i);
            cx += p.x(); cy += p.y(); cz += p.z();
        }
        cx /= n; cy /= n; cz /= n;
        // Covariance (symmetric, upper triangle).
        double sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0, szz = 0;
        for (size_t i = 0; i < n; ++i) {
            auto p = pc->point3d(i);
            const double dx = p.x() - cx, dy = p.y() - cy, dz = p.z() - cz;
            sxx += dx * dx; sxy += dx * dy; sxz += dx * dz;
            syy += dy * dy; syz += dy * dz; szz += dz * dz;
        }
        auto matvec = [&](const double v[3], double o[3]) {
            o[0] = sxx * v[0] + sxy * v[1] + sxz * v[2];
            o[1] = sxy * v[0] + syy * v[1] + syz * v[2];
            o[2] = sxz * v[0] + syz * v[1] + szz * v[2];
        };
        // Two leading axes by power iteration + deflation against axis 1.
        double axes[2][3];
        double evals[2] = {0, 0};
        for (int ax = 0; ax < 2; ++ax) {
            double v[3] = {1.0, 0.7, 0.4};  // fixed seed
            for (int it = 0; it < 64; ++it) {
                if (ax == 1) {  // deflate: remove the axis-0 component
                    const double d = v[0]*axes[0][0] + v[1]*axes[0][1] + v[2]*axes[0][2];
                    v[0] -= d * axes[0][0]; v[1] -= d * axes[0][1]; v[2] -= d * axes[0][2];
                }
                double o[3];
                matvec(v, o);
                const double m = std::sqrt(o[0]*o[0] + o[1]*o[1] + o[2]*o[2]);
                if (m <= 0) break;
                v[0] = o[0] / m; v[1] = o[1] / m; v[2] = o[2] / m;
            }
            if (ax == 1) {  // final re-orthogonalization
                const double d = v[0]*axes[0][0] + v[1]*axes[0][1] + v[2]*axes[0][2];
                v[0] -= d * axes[0][0]; v[1] -= d * axes[0][1]; v[2] -= d * axes[0][2];
                const double m = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
                if (m > 0) { v[0] /= m; v[1] /= m; v[2] /= m; }
            }
            axes[ax][0] = v[0]; axes[ax][1] = v[1]; axes[ax][2] = v[2];
            double o[3];
            matvec(v, o);
            evals[ax] = (o[0]*v[0] + o[1]*v[1] + o[2]*v[2]) / n;
        }
        const WireCell::Point sp = shower->get_start_point();
        const int inpi0 = pi0_showers.count(shower) ? 1 : 0;
        for (int ax = 0; ax < 2; ++ax) {
            double c1x = 0, c1y = 0, c1z = 0, c2x = 0, c2y = 0, c2z = 0;
            size_t n1 = 0, n2 = 0;
            for (size_t i = 0; i < n; ++i) {
                auto p = pc->point3d(i);
                const double d = (p.x()-cx)*axes[ax][0] + (p.y()-cy)*axes[ax][1] + (p.z()-cz)*axes[ax][2];
                if (d >= 0) { c1x += p.x(); c1y += p.y(); c1z += p.z(); ++n1; }
                else        { c2x += p.x(); c2y += p.y(); c2z += p.z(); ++n2; }
            }
            if (n1 == 0 || n2 == 0) {
                std::fprintf(stderr, "PI0_SUBSTRUCT sh=%d E=%.1f npts=%zu inpi0=%d ax=%d skip=degenerate_split\n",
                             pr132_pi0_shid(shower), kq / WireCell::units::MeV, n, inpi0, ax);
                continue;
            }
            WireCell::Vector d1(c1x/n1 - sp.x(), c1y/n1 - sp.y(), c1z/n1 - sp.z());
            WireCell::Vector d2(c2x/n2 - sp.x(), c2y/n2 - sp.y(), c2z/n2 - sp.z());
            const double m1 = d1.magnitude(), m2 = d2.magnitude();
            if (m1 <= 0 || m2 <= 0) continue;
            const double angle = std::acos(std::clamp(d1.dot(d2) / (m1 * m2), -1.0, 1.0));
            const double e1 = kq * (double)n1 / (double)n;
            const double e2 = kq * (double)n2 / (double)n;
            const double mass = std::sqrt(4 * e1 * e2 * std::pow(std::sin(angle / 2.0), 2));
            std::fprintf(stderr,
                         "PI0_SUBSTRUCT sh=%d E=%.1f npts=%zu inpi0=%d ax=%d ev01=%.3f m=%.1f E1=%.1f E2=%.1f ang=%.1f n1=%zu n2=%zu\n",
                         pr132_pi0_shid(shower), kq / WireCell::units::MeV, n, inpi0, ax,
                         evals[0] > 0 ? evals[1] / evals[0] : -1.0,
                         mass / WireCell::units::MeV, e1 / WireCell::units::MeV, e2 / WireCell::units::MeV,
                         angle / M_PI * 180.0, n1, n2);
        }
    }
}

// doc pr/132 round 5 (K16) -- EM collinear-fragment merge, the build-time
// answer to the fragmentation class the pairing-time K12 could not fix (132
// doc sec 10.3) and the existing seats never see: the round-5 gamma ledger
// (pr132_gamma_ledger.py) measures the UNDER+SIBS fragments at 2.6-8.5 deg
// off the host gamma's axis but 25-108 cm downstream (specimens SBND
// 18255-76346 g1, 342199 g1+g2, 105946 g1, 54332 g1) -- far beyond the
// pr/118 axis-merge and pr/125 satellite-absorb reach (and above the 10 MeV
// satellite cap).  For each EM host (leading-first), absorb every remaining
// DETACHED (conn_type != 1 -- the primary-electron guard) EM shower whose
// start lies within the axis cone (angle < m_em_collinear_merge_deg from
// the host's 15-cm direction at its start) and within
// m_em_collinear_merge_dis of the host start; a fragment long enough to
// carry its own direction (> 10 cm) must also CONTINUE the host axis
// (< 30 deg, signed -- the pr/128 lesson: proximity alone admits cosmics).
// Hosts below m_em_collinear_merge_min_host never absorb.  Vertices are
// never touched (the owner's round-5 scope rule).  deg 0 (default) = off =
// byte-identical.
// ===================================================================== pr/138
// doc sbnd_xin/docs/pr/138 Phase B -- the EM SHOWER SPLITTER.
//
// Every other shower pass in this chain merges.  This one cuts, and it runs
// LAST (after shower_dedup_start_seg, pass4_prune_detached, pass4_prune_gap2,
// samevtx_absorb, satellite_absorb, em_collinear_merge and em_start_backext,
// all SBND production ON) because the owner's architecture is "merge them
// together, then separate cleanly".  It runs BEFORE the pi0 finders on purpose:
// a gamma pair that was over-clustered into one shower can only be PAIRED into
// a pi0 after it has been cut apart.
//
// TRIGGER AND KERNEL ARE THE SAME OBJECT, which is the whole borrow.  ATLAS
// topo-cluster splitting, CMS particle flow and GARLIC all seed on local
// maxima and let the SEED COUNT be the multiplicity decision; doc pr/137 sec 10
// found that shape in the literature and doc pr/138 sec A5.3 measured it here:
// `valley_best` -- the charge dip between the two best angular maxima -- scored
// AUC 0.930 and purity 0.857 at 50% efficiency against 164 owner hand labels,
// with the owner's own three stated factors (direction, size, distance) behind
// it in his order.
//
// The numbers are Phase A measurements, and NONE of them is free to tune here:
//   valley <= 0.95   sec A5.4; a fit-half threshold scan says 0.95 IS the knee
//                    (0.90 -> 0.875/0.913, 0.95 -> 0.917/0.917, 0.99 ->
//                    0.917/0.815) and the holdout has been opened once.
//   frac  >= 0.03    the trigger's own pair filter (split_model.propose()).
//   w_single(r) = 3.575 + 0.0283*r cm  the in-situ single-shower width null,
//                    FITTED on 346 SINGLE showers by pr137_null_model.py.  It
//                    is not a PDG number: the LAr Moliere radius is quoted in
//                    the doc for scale and is a threshold nowhere.  sigma(r) =
//                    w_single(r)/r radians, clipped to [2, 60] degrees, so a
//                    compact late-converting gamma is not smoothed away by the
//                    kernel that fits the near wide one (owner factors 2+3).
//   sep_scale 1.6, max_seeds 4, 25 arc samples, cap 4000 points
//                    frozen at the Phase A operating point -- moving any of
//                    them moves n_seed, which moves the trigger, which
//                    invalidates sec A5.4's fire list.
//
// Reference implementation and the source of every constant:
//   sbnd_xin/scripts/pr137_lib.py:391-455   profile_sigma_fn / angular_maxima
//   sbnd_xin/split_display/split_model.py:106-205  propose()
//   sbnd_xin/scripts/pr138_kernel_k.py      the sec B3 kernel comparison
//
// stderr tape WCT_SHOWER_SPLIT_DEBUG, byte-neutral, gate tested first so the
// cost with the env unset is one static bool.
static inline bool pr138_split_dbg()
{
    static const bool dbg = std::getenv("WCT_SHOWER_SPLIT_DEBUG") != nullptr;
    return dbg;
}

namespace {

// The angular point cloud of one shower, seen from the reference vertex.
struct Pr138Cloud {
    std::vector<WireCell::Vector> u;    // unit ray, vertex -> point
    std::vector<double> w;              // charge weight, clamped at 0
    std::vector<double> sig;            // per-point kernel bandwidth, radians
    double wsum{0};
};

struct Pr138Maxima {
    std::vector<WireCell::Vector> dirs;             // seed directions, density order
    std::vector<double> dens;                       // density at each seed
    std::vector<std::vector<double>> valley;        // pairwise valley depth
    std::vector<double> frac;                       // charge share nearest each seed
};

// pr137_lib.py:391 profile_sigma_fn -- w_single(r)/r, clipped.
inline double pr138_sigma(double r_cm)
{
    static const double lo = 2.0 * M_PI / 180.0, hi = 60.0 * M_PI / 180.0;
    const double r = std::max(r_cm, 1e-6);
    const double s = (3.575 + 0.0283 * r) / r;
    return std::clamp(s, lo, hi);
}

// pr137_lib.py:129 rays + :136 qwt.  dQ in the dump can be NEGATIVE
// (noise-subtracted); an unclamped weighted statistic then goes negative, so
// every weight in this module goes through the same clamp, and an all-zero
// cloud falls back to unit weights exactly as qwt() does.
Pr138Cloud pr138_cloud(const std::vector<SegmentPtr>& segs, const WireCell::Point& v,
                       size_t max_pts = 4000)
{
    Pr138Cloud C;
    std::vector<const PR::Fit*> fp;
    for (const auto& sg : segs) {
        if (!sg) continue;
        for (const auto& f : sg->fits()) if (f.valid()) fp.push_back(&f);
    }
    // Deterministic stride decimation above the cap.  The Phase A population
    // never reaches it (max 1414 points over 172 scanned objects), and the
    // offline reference draws a seeded random subsample there instead -- so a
    // fire when this line trips is NOT comparable to the offline number and
    // says so on the tape.
    const size_t stride = fp.size() > max_pts ? (fp.size() + max_pts - 1) / max_pts : 1;
    for (size_t i = 0; i < fp.size(); i += stride) {
        const WireCell::Vector d = fp[i]->point - v;
        double r = d.magnitude();
        if (r <= 0) r = 1e-9;
        C.u.push_back(d / r);
        C.w.push_back(std::max(fp[i]->dQ, 0.0));
        C.sig.push_back(pr138_sigma(r / WireCell::units::cm));
    }
    for (double x : C.w) C.wsum += x;
    if (C.wsum <= 0) { for (auto& x : C.w) x = 1.0; C.wsum = (double) C.w.size(); }
    return C;
}

// Charge-weighted angular density at one direction.  The bandwidth belongs to
// the SOURCE point, not the probe direction -- pr137_lib.py:445 broadcasts
// sig[None,:] over the second index and this must match it exactly.
inline double pr138_density(const Pr138Cloud& C, const WireCell::Vector& d)
{
    double s = 0;
    for (size_t j = 0; j < C.u.size(); ++j) {
        const double a = std::acos(std::clamp(d.dot(C.u[j]), -1.0, 1.0));
        const double t = a / C.sig[j];
        s += std::exp(-0.5 * t * t) * C.w[j];
    }
    return s;
}

// pr137_lib.py:405 angular_maxima.  Local maxima of the density, separated by
// more than the in-situ single-shower width, with the pairwise VALLEY between
// them: a genuine second core has a charge dip on the great-circle arc to the
// first; a bright patch inside ONE shower does not.
Pr138Maxima pr138_angular_maxima(const Pr138Cloud& C, double sep_scale = 1.6,
                                 size_t max_seeds = 4, int narc = 25)
{
    Pr138Maxima M;
    const size_t n = C.u.size();
    if (n < 8) return M;
    std::vector<double> dens(n, 0.0);
    for (size_t i = 0; i < n; ++i) dens[i] = pr138_density(C, C.u[i]);
    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = i;
    // Stable by construction: index breaks a density tie, so the seed set can
    // never depend on sort implementation or on pointer order (CLAUDE.md M4).
    std::stable_sort(order.begin(), order.end(), [&dens](size_t a, size_t b) {
        if (dens[a] != dens[b]) return dens[a] > dens[b];
        return a < b;
    });
    std::vector<size_t> seeds;
    for (size_t oi = 0; oi < order.size() && seeds.size() < max_seeds; ++oi) {
        const size_t i = order[oi];
        bool close = false;
        for (size_t j : seeds) {
            const double a = std::acos(std::clamp(C.u[i].dot(C.u[j]), -1.0, 1.0));
            if (a < sep_scale * std::max(C.sig[i], C.sig[j])) { close = true; break; }
        }
        if (!close) seeds.push_back(i);
    }
    const size_t k = seeds.size();
    M.dirs.resize(k); M.dens.resize(k); M.frac.assign(k, 0.0);
    M.valley.assign(k, std::vector<double>(k, 1.0));
    for (size_t i = 0; i < k; ++i) { M.dirs[i] = C.u[seeds[i]]; M.dens[i] = dens[seeds[i]]; }
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = i + 1; j < k; ++j) {
            double mn = -1.0;
            for (int t = 0; t < narc; ++t) {
                const double tt = (double) t / (double) (narc - 1);
                WireCell::Vector arc = M.dirs[i] * (1.0 - tt) + M.dirs[j] * tt;
                const double nn = arc.magnitude();
                arc = (nn > 0) ? arc / nn : M.dirs[i];
                const double d = pr138_density(C, arc);
                if (mn < 0 || d < mn) mn = d;
            }
            const double denom = std::max(std::min(M.dens[i], M.dens[j]), 1e-12);
            M.valley[i][j] = M.valley[j][i] = mn / denom;
        }
    }
    if (k) {
        for (size_t p = 0; p < n; ++p) {
            size_t best = 0; double bs = -2.0;
            for (size_t i = 0; i < k; ++i) {
                const double s = C.u[p].dot(M.dirs[i]);
                if (s > bs) { bs = s; best = i; }
            }
            M.frac[best] += C.w[p];
        }
        double f = 0; for (double x : M.frac) f += x;
        if (f > 0) for (auto& x : M.frac) x /= f;
    }
    return M;
}

// The ACCEPT test and the multiplicity decision, in one greedy walk -- doc
// pr/138 sec B2 + B3.  The brightest maximum is always a seed.  A later one
// joins iff it carries `min_frac` of the charge AND a charge valley at most
// `max_valley` separates it from a seed already accepted, i.e. it is a distinct
// core rather than a bright patch.  Firing is |accepted| >= 2, so the trigger
// IS the seed count, which is the shape sec 10 borrowed.
//
// Measured against the 164 owner-labelled EM objects (pr138_kernel_k.py):
// efficiency 0.791, purity 0.791.  The shipped Phase A rule -- best PAIR by
// minimum valley -- scores 0.791/0.773 and differs on exactly ONE object of
// 164 (evt281567 node99193, owner KEEP, which this rule correctly declines).
std::vector<size_t> pr138_accept(const Pr138Maxima& M, double max_valley,
                                 double min_frac, size_t cap)
{
    const size_t k = M.dirs.size();
    if (k < 1 || cap < 2) return {};
    std::vector<size_t> acc{0};
    for (size_t s = 1; s < k && acc.size() < cap; ++s) {
        if (M.frac[s] < min_frac) continue;
        double best = 2.0;
        for (size_t a : acc) best = std::min(best, M.valley[s][a]);
        if (best <= max_valley) acc.push_back(s);
    }
    if (acc.size() < 2) return {};
    return acc;
}

}  // namespace

// Single-linkage components of a member list at `gap` -- the BUNDLE, the
// spatially connected unit the owner's hand labels move wholesale.  Same idiom
// and same gap semantics as shower_pass4_prune_detached (:8607-8656); the
// distance is the exact minimum between the two segments' fit clouds, which is
// what the offline reference (scipy cKDTree) computes.
static std::vector<std::vector<size_t>> pr138_bundles(const std::vector<SegmentPtr>& mem,
                                                      double gap)
{
    const size_t n = mem.size();
    std::vector<size_t> parent(n);
    for (size_t i = 0; i < n; ++i) parent[i] = i;
    std::function<size_t(size_t)> find = [&](size_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (find(i) == find(j)) continue;
            double best = -1.0;
            for (const auto& fa : mem[i]->fits()) {
                if (!fa.valid()) continue;
                for (const auto& fb : mem[j]->fits()) {
                    if (!fb.valid()) continue;
                    const double d = (fa.point - fb.point).magnitude();
                    if (best < 0 || d < best) best = d;
                    if (best < gap) break;
                }
                if (best >= 0 && best < gap) break;
            }
            if (best >= 0 && best < gap) { size_t pi = find(i), pj = find(j); if (pi != pj) parent[pi] = pj; }
        }
    }
    std::map<size_t, std::vector<size_t>> comp;
    for (size_t i = 0; i < n; ++i) comp[find(i)].push_back(i);
    std::vector<std::vector<size_t>> out;
    for (auto& [root, idxs] : comp) out.push_back(idxs);
    // Descending charge, min-member-index tie-break: bundle 0 is the main body
    // and the order never depends on a pointer.
    std::vector<double> bq(out.size(), 0.0);
    for (size_t b = 0; b < out.size(); ++b)
        for (size_t i : out[b])
            for (const auto& f : mem[i]->fits()) if (f.valid()) bq[b] += std::max(f.dQ, 0.0);
    std::vector<size_t> ord(out.size());
    for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
        if (bq[a] != bq[b]) return bq[a] > bq[b];
        return out[a].front() < out[b].front();
    });
    std::vector<std::vector<size_t>> sorted;
    for (size_t i : ord) sorted.push_back(out[i]);
    return sorted;
}

// Assign one unit (a bundle, or a single segment) to an accepted seed.
// Returns (seed index into `acc`, that seed's share of the unit's charge).
//
// TWO MODES, and which one is used is decided by k, not by a free parameter.
//   k == 2   the unit's charge-weighted CENTROID RAY picks the seed.  This is
//            what split_model.propose() does and it is what the owner
//            corrected; conditioned on the trigger firing it reproduces his
//            two-way boundary with median charge agreement 1.000 and mean
//            0.974 over 27 objects (sec B3).  Nothing may touch it.
//   k >= 3   a bundle can now straddle two parts, and a centroid ray sends the
//            whole of it to one.  So the seed holding the PLURALITY of the
//            unit's own per-point charge wins, and a bundle whose winner holds
//            less than `snap` of it is cut at the SEGMENT level -- the finest
//            cut detach_member_set can make (PRShower.cxx:640-700).  Measured:
//            k>=3 mean agreement 0.573 -> 0.772.
static std::pair<size_t, double> pr138_place(const std::vector<SegmentPtr>& mem,
                                             const std::vector<size_t>& unit,
                                             const WireCell::Point& v,
                                             const std::vector<WireCell::Vector>& dirs,
                                             bool plurality)
{
    const size_t k = dirs.size();
    if (!plurality) {
        WireCell::Vector c(0, 0, 0);
        double wsum = 0;
        for (size_t i : unit)
            for (const auto& f : mem[i]->fits()) {
                if (!f.valid()) continue;
                const double w = std::max(f.dQ, 0.0);
                c = c + f.point * w; wsum += w;
            }
        if (wsum <= 0) {   // qwt() falls back to unit weights on an all-zero cloud
            c = WireCell::Vector(0, 0, 0); wsum = 0;
            for (size_t i : unit)
                for (const auto& f : mem[i]->fits()) {
                    if (!f.valid()) continue;
                    c = c + f.point; wsum += 1.0;
                }
        }
        if (wsum <= 0) return {k, 0.0};
        WireCell::Vector ray = c / wsum - v;
        const double n = ray.magnitude();
        if (n <= 0) return {k, 0.0};
        ray = ray / n;
        size_t best = 0; double bs = -2.0;
        for (size_t i = 0; i < k; ++i) { const double s = ray.dot(dirs[i]); if (s > bs) { bs = s; best = i; } }
        return {best, 1.0};
    }
    std::vector<double> share(k, 0.0);
    for (size_t i : unit)
        for (const auto& f : mem[i]->fits()) {
            if (!f.valid()) continue;
            WireCell::Vector r = f.point - v;
            const double n = r.magnitude();
            if (n <= 0) continue;
            r = r / n;
            size_t best = 0; double bs = -2.0;
            for (size_t j = 0; j < k; ++j) { const double s = r.dot(dirs[j]); if (s > bs) { bs = s; best = j; } }
            share[best] += std::max(f.dQ, 0.0);
        }
    double tot = 0; for (double x : share) tot += x;
    if (tot <= 0) return {k, 0.0};
    size_t best = 0;
    for (size_t j = 1; j < k; ++j) if (share[j] > share[best]) best = j;
    return {best, share[best] / tot};
}

// doc sbnd_xin/docs/pr/139 P1.2 -- the impact parameter b: how far the object's
// OWN charge-weighted principal axis misses the reference vertex.  pr/138 sec
// 4.2b killed void_frac (AUC 0.146, BACKWARDS) because everything measured
// ALONG the vertex ray is confounded -- seen from a wrong origin one shower
// genuinely looks like two well-separated ones.  b is perpendicular, and on the
// eight census movers every gain sits below 11 cm and every loss above 13.
//
// Transcribed from sbnd_xin/scripts/pr139_pointing.py:38-46, which builds
// X = (pts - c) * sqrt(w) and takes the first right singular vector.  The
// weighted covariance below has the same eigenvectors (X^T X), so the two
// agree to rounding; the tape prints b so that identity is CHECKED on the arm
// rather than asserted (the pr/138 sec B1 discipline).  Weights follow
// pr137_lib.qwt(): dQ clamped at zero, unit weights if the clamp empties it.
// NO decimation here -- the offline reference uses every point.
inline double pr139_impact(const std::vector<SegmentPtr>& segs, const WireCell::Point& v)
{
    std::vector<const PR::Fit*> fp;
    for (const auto& sg : segs) {
        if (!sg) continue;
        for (const auto& f : sg->fits()) if (f.valid()) fp.push_back(&f);
    }
    if (fp.size() < 8) return -1.0;
    std::vector<double> w(fp.size());
    double wsum = 0;
    for (size_t i = 0; i < fp.size(); ++i) { w[i] = std::max(fp[i]->dQ, 0.0); wsum += w[i]; }
    if (wsum <= 0) { for (auto& x : w) x = 1.0; wsum = (double) w.size(); }
    WireCell::Vector c(0, 0, 0);
    for (size_t i = 0; i < fp.size(); ++i) c = c + fp[i]->point * w[i];
    c = c / wsum;
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (size_t i = 0; i < fp.size(); ++i) {
        const WireCell::Vector d = fp[i]->point - c;
        const double dv[3] = {d.x(), d.y(), d.z()};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) cov(a, b) += w[i] * dv[a] * dv[b];
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
    if (es.info() != Eigen::Success) return -1.0;
    const Eigen::Vector3d e = es.eigenvectors().col(2);   // ascending eigenvalues
    WireCell::Vector ax(e(0), e(1), e(2));
    const double an = ax.magnitude();
    if (an <= 0) return -1.0;
    ax = ax / an;
    const WireCell::Vector d = c - v;
    return (d - ax * d.dot(ax)).magnitude();
}

/// doc sbnd_xin/docs/pr/138 Phase B -- the EM shower splitter.  Full docstring
/// and provenance at the pr138_split_dbg() block above.  No-op when
/// m_shower_split is false, which is the default and makes the pass
/// byte-identical.
void PatternAlgorithms::shower_split(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers,
    ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower,
    VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters,
    IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon,
    TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
    const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    const bool dbg = pr138_split_dbg();
    if (!m_shower_split && !dbg) return;
    if (!main_vertex) return;
    // The reference vertex is the event's NEUTRINO main vertex for every
    // candidate -- the same choice pr137_lib.build_population() makes, so the
    // offline features and these are measured from the same point.  NOTE the
    // pi0 finders re-seat this vertex at an accepted two-photon decay point
    // (:7886 and :6241, doc sec A1.4) and they run AFTER this pass, so on a
    // pi0 event the calib dump's main_vertex is NOT the point used here.  The
    // tape prints the vertex it used so the comparison is controlled rather
    // than explained away.
    const WireCell::Point v = main_vertex->fit().valid() ? main_vertex->fit().point
                                                         : main_vertex->wcpt().point;

    // Deterministic candidate order: cluster id, then segment id.  Never
    // iterate the shower set itself (CLAUDE.md determinism rule).
    std::vector<ShowerPtr> split_order(showers.begin(), showers.end());
    std::sort(split_order.begin(), split_order.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
        auto* sa = a->start_segment().get();
        auto* sb = b->start_segment().get();
        if (!sa || !sb) return sa < sb;
        const int ca = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
        const int cb = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
        if (ca != cb) return ca < cb;
        return sa->id() < sb->id();
    });

    const size_t cap = (size_t) std::max(2, m_shower_split_max_parts);
    int n_fired = 0, n_peeled = 0, n_shared_refused = 0, n_vetoed = 0, n_shed = 0;
    std::vector<ShowerPtr> em_retype;   // doc pr/139 sec 25
    int n_em_retyped = 0;
    std::vector<ShowerPtr> daughters;
    std::vector<ShowerPtr> dau_parent;   // aligned with `daughters` (doc pr/139 P1.4)

    // doc pr/139 P1.1 -- how many showers own each segment.  WCT showers MAY
    // share members, and a peel of shared segments makes a DUPLICATE object:
    // kine_charge is computed per shower from its own point clouds against the
    // 2D charge maps with no cross-shower dedup, so the charge lands on one of
    // the two arbitrarily.  Both failure modes are in the shipped arms
    // (work-pr138r2-c90on-*): evt281485's daughter gets 0.00 MeV while its
    // co-owner gains the whole 38.9, and evt165157's daughter is seeded on a
    // segment that already ROOTS another shower.  7 of 50 fired parents hold
    // shared members.  The dump's segment->shower field is single-valued and
    // cannot see this; the count below can.
    //
    // Counted into an unordered_map keyed by the segment's GRAPH INDEX, not by
    // pointer -- addition is order-free, and nothing here iterates the map.
    std::unordered_map<size_t, int> seg_nshower;
    if (m_shower_split_skip_shared) {
        for (const auto& sh : split_order) {
            if (!sh) continue;
            std::unordered_set<size_t> seen;
            for (auto edesc : ordered_edges(*sh, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg || !sg->descriptor_valid()) continue;
                const size_t idx = graph[sg->get_descriptor()].index;
                if (seen.insert(idx).second) seg_nshower[idx] += 1;
            }
        }
    }

    for (auto& shower : split_order) {
        SegmentPtr ss = shower ? shower->start_segment() : nullptr;
        if (!ss || !ss->descriptor_valid()) continue;
        // Member list in stable view order, unique -- transcribed from
        // shower_pass4_prune_detached (:8610-8618).
        std::vector<SegmentPtr> mem;
        {
            std::unordered_set<size_t> seen;
            for (auto edesc : ordered_edges(*shower, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg || !sg->descriptor_valid()) continue;
                if (seen.insert(graph[sg->get_descriptor()].index).second) mem.push_back(sg);
            }
        }
        if ((int) mem.size() < m_shower_split_min_nseg) continue;
        double Q = 0; size_t npts = 0;
        for (const auto& sg : mem)
            for (const auto& f : sg->fits()) if (f.valid()) { Q += f.dQ; ++npts; }
        if (Q < m_shower_split_min_charge) continue;

        const Pr138Cloud C = pr138_cloud(mem, v);
        // doc pr/139 sec 17 -- the seed cap is the REAL k>=3 constraint, not
        // `max_parts`.  Measured on the 2026-09-01 labels: n_seed sits at the
        // hardcoded 4 on 76 % of fired candidates and on ALL FOUR objects the
        // owner cut into k >= 3, whose k values are 3, 3, 5 and 7 -- a 4-seed
        // finder cannot express 7.  Default 4 keeps the shipped behaviour
        // byte for byte.
        const Pr138Maxima M = pr138_angular_maxima(
            C, 1.6, (size_t) std::max(2, m_shower_split_max_seeds));
        const std::vector<size_t> acc = pr138_accept(M, m_shower_split_max_valley,
                                                     m_shower_split_min_frac, cap);
        // valley_best / angle_best for the tape: the same pair statistic sec
        // A5.3 ranked, over the seeds carrying at least min_frac.
        double vbest = 1.0, abest = -1.0;
        for (size_t i = 0; i < M.dirs.size(); ++i)
            for (size_t j = i + 1; j < M.dirs.size(); ++j) {
                if (std::min(M.frac[i], M.frac[j]) < m_shower_split_min_frac) continue;
                if (M.valley[i][j] < vbest) {
                    vbest = M.valley[i][j];
                    abest = std::acos(std::clamp(M.dirs[i].dot(M.dirs[j]), -1.0, 1.0)) * 180.0 / M_PI;
                }
            }

        std::vector<size_t> part(mem.size(), 0);
        size_t nparts = 1;
        if (!acc.empty()) {
            std::vector<WireCell::Vector> dirs;
            for (size_t i : acc) dirs.push_back(M.dirs[i]);
            const bool plural = dirs.size() >= 3;
            for (const auto& b : pr138_bundles(mem, m_shower_split_bundle_gap)) {
                auto [g, dom] = pr138_place(mem, b, v, dirs, plural);
                if (g >= dirs.size()) continue;
                if (plural && dom < m_shower_split_snap && b.size() > 1) {
                    for (size_t i : b) {
                        auto [gi, di] = pr138_place(mem, {i}, v, dirs, true);
                        (void) di;
                        if (gi < dirs.size()) part[i] = gi;
                    }
                    continue;
                }
                for (size_t i : b) part[i] = g;
            }
            std::set<size_t> used(part.begin(), part.end());
            nparts = used.size();
        }
        // The HONESTY CHECK, transcribed from split_model.propose():205 -- the
        // seeds can be well separated and every bundle still fall to one of
        // them (evt389538 node19021: valley_best 0.091 yet one group), because
        // a minority lobe that is not its own connected component has no unit
        // to carry it.  Firing with one part is not a split.
        const bool fired = (!acc.empty() && nparts >= 2);
        if (fired) ++n_fired;

        // doc pr/139 P1.2 -- the impact-parameter veto.  Computed whenever the
        // tape is on (so the C++ value can be checked against
        // pr139_pointing.py on the same arm) or whenever the bound is armed.
        // m_shower_split_max_impact == 0 is the shipped behaviour: no bound.
        // A cloud too small for pr139_impact() returns -1 and is NEVER vetoed
        // -- absence of the feature must not silently suppress a fire.
        double bimp = -1.0;
        if (dbg || m_shower_split_max_impact > 0) bimp = pr139_impact(mem, v);
        const bool vetoed = fired && m_shower_split_max_impact > 0 && bimp >= 0
                            && bimp > m_shower_split_max_impact;
        if (vetoed) ++n_vetoed;

        if (dbg) {
            double vgap = -1.0;
            for (const auto& sg : mem)
                for (const auto& f : sg->fits())
                    if (f.valid()) {
                        const double d = (f.point - v).magnitude();
                        if (vgap < 0 || d < vgap) vgap = d;
                    }
            std::fprintf(stderr,
                "SHOWER_SPLIT cand shower=%d pdg=%d nseg=%zu npts=%zu Q=%.4g "
                "n_seed=%zu valley_best=%.4f angle_best=%.2f nacc=%zu nparts=%zu "
                "fired=%d decim=%d vtx=%.2f,%.2f,%.2f vgap_cm=%.2f "
                "vchi2=%.3f vdQ=%.4g vfit=%d b_cm=%.2f veto=%d\n",
                pr91_seg_display_id(ss), shower->get_particle_type(), mem.size(), npts, Q,
                M.dirs.size(), vbest, abest, acc.size(), nparts, fired ? 1 : 0,
                npts > 4000 ? 1 : 0, v.x() / units::cm, v.y() / units::cm, v.z() / units::cm,
                vgap >= 0 ? vgap / units::cm : -1.0,
                main_vertex->fit().reduced_chi2, main_vertex->fit().dQ,
                main_vertex->fit().valid() ? 1 : 0,
                bimp >= 0 ? bimp / units::cm : -1.0, vetoed ? 1 : 0);
            if (fired)
                for (size_t g = 0; g < acc.size(); ++g) {
                    std::string segs; double qg = 0; size_t ng = 0;
                    for (size_t i = 0; i < mem.size(); ++i) {
                        if (part[i] != g) continue;
                        ++ng;
                        for (const auto& f : mem[i]->fits()) if (f.valid()) qg += f.dQ;
                        if (!segs.empty()) segs += ",";
                        segs += std::to_string(pr91_seg_display_id(mem[i]));
                    }
                    if (!ng) continue;
                    std::fprintf(stderr, "SHOWER_SPLIT part shower=%d part=%zu nseg=%zu q=%.4g segs=%s\n",
                                 pr91_seg_display_id(ss), g, ng, qg, segs.c_str());
                }
        }

        if (!m_shower_split || !fired || vetoed) continue;
        // Long-muon pseudo-showers are never cut -- same exclusion as
        // shower_pass4_prune_detached (:8607).  The TAPE above deliberately
        // still records them: 8 of the 172 scanned objects are track-typed
        // (sec A1.6) and dropping them from the tape would shrink the
        // offline comparison population.
        if (std::abs(shower->get_particle_type()) == 13) continue;

        const auto ss_it = std::find(mem.begin(), mem.end(), ss);
        if (ss_it == mem.end()) continue;
        const size_t keep_part = part[std::distance(mem.begin(), ss_it)];
        const WireCell::Point sv_pt = shower->start_vertex()
            ? (shower->start_vertex()->fit().valid() ? shower->start_vertex()->fit().point
                                                     : shower->start_vertex()->wcpt().point)
            : WireCell::Point(0, 0, 0);
        for (size_t g = 0; g < acc.size(); ++g) {
            if (g == keep_part) continue;
            std::vector<SegmentPtr> comp;
            for (size_t i = 0; i < mem.size(); ++i) if (part[i] == g) comp.push_back(mem[i]);
            if (comp.empty()) continue;
            // doc pr/139 P1.1 -- refuse a component holding a segment that a
            // shower OTHER than this parent also owns.  A refusal, not a
            // repair: making a shared peel work means dividing 2D charge
            // between duplicate objects, which is the kine_charge dedup
            // problem (kine_charge_owned_scan, knob-gated, runs later).
            if (m_shower_split_skip_shared) {
                // Count rather than short-circuit: doc pr/139 §15 needs to know
                // WHETHER THE WHOLE COMPONENT IS CO-OWNED, not merely that some
                // member is, and the tape has to carry the split so the shed
                // rule below is priced instead of assumed.
                size_t nshared = 0;
                double q_shared = 0, q_comp = 0;
                for (const auto& sg : comp) {
                    if (!sg) continue;
                    double qs = 0;
                    for (const auto& f : sg->fits()) if (f.valid()) qs += f.dQ;
                    q_comp += qs;
                    if (!sg->descriptor_valid()) continue;
                    auto it = seg_nshower.find(graph[sg->get_descriptor()].index);
                    if (it != seg_nshower.end() && it->second > 1) { ++nshared; q_shared += qs; }
                }
                if (nshared) {
                    // doc pr/139 §15 -- the SHED.  When EVERY member of the
                    // refused component is also owned by another shower, the
                    // owner's cut can be honoured without building anything:
                    // detach the component from this parent and create NO
                    // daughter, because the charge already has a home in the
                    // co-owner.  That is the whole fix for evt281485, where all
                    // four members of part 1 are held by shower 91112 as well,
                    // so the peel that P1.1 refuses would have produced a
                    // kine_charge = 0 duplicate (§2) and the "drop the shared
                    // members and peel the remainder" alternative has NOTHING
                    // left to peel.
                    //
                    // The partial case is deliberately still refused.  Measured
                    // on the pr136 sidecars: evt165157 part 0 is 2 shared of 7
                    // and evt350354 part 1 is 1 of 12, and peeling their
                    // remainders would MAKE A CUT ON 165157 -- an object the
                    // owner labels KEEP.  A rule that fires there buys one
                    // confirmed cut and pays a false fire, so the shed is
                    // restricted to the unambiguous all-shared case.
                    const bool shed = m_shower_split_shed_shared && nshared == comp.size();
                    if (shed) {
                        const int nrem = shower->detach_member_set(comp);
                        if (nrem) {
                            ++n_shed;
                            if (dbg)
                                std::fprintf(stderr,
                                    "SHOWER_SPLIT shed shower=%d part=%zu nseg=%zu q=%.4g\n",
                                    pr91_seg_display_id(ss), g, comp.size(), q_comp);
                            SPDLOG_LOGGER_DEBUG(s_log,
                                "pr139 shower_split: shed all-shared component shower_id={} part={} nseg={}",
                                shower->get_shower_id(), g, comp.size());
                            continue;
                        }
                        // detach refused (the component holds the start
                        // segment): fall through to the ordinary refusal.
                    }
                    ++n_shared_refused;
                    if (dbg)
                        std::fprintf(stderr,
                            "SHOWER_SPLIT shared shower=%d part=%zu nseg=%zu nshared=%zu "
                            "q_excl_frac=%.3f\n",
                            pr91_seg_display_id(ss), g, comp.size(), nshared,
                            q_comp > 0 ? (q_comp - q_shared) / q_comp : 0.0);
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr139 shower_split: shared-member refusal shower_id={} part={} nseg={} nshared={}",
                        shower->get_shower_id(), g, comp.size(), nshared);
                    continue;
                }
            }
            const int removed = shower->detach_member_set(comp);
            if (!removed) {
                SPDLOG_LOGGER_DEBUG(s_log, "pr138 shower_split: refuse shower_id={} part={} nseg={}",
                                    shower->get_shower_id(), g, comp.size());
                continue;
            }
            // The daughter's start segment is the member NEAREST THE REFERENCE
            // VERTEX, not the one nearest the parent body the way
            // pass4_prune_detached picks it.  A prune re-seeds a fragment that
            // fell off; a split re-seeds a SHOWER, and the pi0 finders run
            // after this pass and read get_start_point() and get_init_dir() --
            // seeding at the downstream end would point the daughter back at
            // its sibling and poison exactly the pi0 mass this round is for.
            //
            // doc pr/139 P1.3 -- and that choice takes NO account of PDG, which
            // is the defect this knob fixes.  SBND runs
            // shower_pdg_from_start_segment, and Shower::calculate_kinematics
            // copies the start segment's pdg verbatim into data.particle_type
            // (PRShower.cxx:1617, :1738), so a daughter whose nearest member
            // happens to be muon-typed becomes a muon shower: 11 of 50
            // daughters in the shipped arm, four of them with an EM-MAJORITY
            // segment composition (evt84229 at 476 MeV, evt269774 at 182).
            // The harm is double.  A mu-typed shower is invisible to the pi0
            // finder except through K20, AND its energy is wrong: kine_charge
            // divides by recom*fudge chosen on get_flag_shower()
            // (NeutrinoEnergyReco.cxx:337-342), which the muon branch clears,
            // so at the SBND values 0.58*0.86 vs 0.87*0.95 the mu-typed energy
            // is LOW BY A FACTOR 1.657.
            // Knob on: prefer the nearest EM-typed member, fall back to the
            // legacy nearest-overall when the daughter has none.  Knob off:
            // the single pass below IS the legacy loop, byte for byte.
            SegmentPtr root_sg = comp.front();
            double root_d = -1.0;
            for (int rpass = 0; rpass < (m_shower_split_em_start ? 2 : 1); ++rpass) {
                const bool em_only = (m_shower_split_em_start && rpass == 0);
                for (const auto& sg : comp) {
                    if (em_only && !(sg->has_particle_info() && sg->particle_info()
                                     && std::abs(sg->particle_info()->pdg()) == 11)) continue;
                    double d = -1.0;
                    for (const auto& f : sg->fits())
                        if (f.valid()) { const double t = (f.point - v).magnitude(); if (d < 0 || t < d) d = t; }
                    if (d >= 0 && (root_d < 0 || d < root_d)) { root_d = d; root_sg = sg; }
                }
                if (em_only && root_d >= 0) break;   // an EM member won the seed
            }
            const double root_sv_dis = segment_get_closest_point(root_sg, sv_pt).first;
            const int conn = (root_sv_dis >= 0 && root_sv_dis <= 80 * units::cm) ? 3 : 4;
            ShowerPtr ns = std::make_shared<Shower>(graph);
            if (shower->start_vertex()) ns->set_start_vertex(shower->start_vertex(), conn);
            ns->set_start_segment(root_sg);
            for (const auto& sg : comp) if (sg != root_sg) ns->add_segment(sg, true);
            daughters.push_back(ns);
            dau_parent.push_back(shower);   // doc pr/139 P1.4: never re-home into the parent
            // doc pr/139 sec 25 -- the owner's rule, 2026-09-03: "the key issue
            // is track vs. shower separation.  If there is a track like MIP
            // muon is a descent candidate.  But if it is an EM shower, electron
            // would be the right choice."
            //
            // The separator cannot make that call on a short object and does not
            // fail loudly -- it DEFAULTS to track.  Measured on this manifest:
            // flag_shower is set on 0.1% of segments below 0.8 cm and 0.6%
            // between 0.8 and 1.2 cm, against 14.7% at 1.2-2.0 cm and 46.7% at
            // 2-4 cm.  So below ~1.2 cm "track-like" is an absence of evidence.
            //
            // A daughter that (a) is shorter than the knob, (b) holds no
            // EM-typed member of its own, and (c) was peeled off an EM-typed
            // parent, is a fragment of an electron shower by construction.  Type
            // it 11 rather than letting a 2-point stub be reported as a muon.
            // A LONGER daughter keeps whatever the separator says, which is the
            // owner's "a track like MIP muon is a decent candidate" half.
            if (m_shower_split_em_type_max_len > 0
                && std::abs(shower->get_particle_type()) == 11) {
                double dlen = 0;
                bool has_em = false;
                for (const auto& sg : comp) {
                    if (!sg) continue;
                    dlen += segment_track_length(sg);
                    if (sg->has_particle_info() && sg->particle_info()
                        && std::abs(sg->particle_info()->pdg()) == 11) has_em = true;
                }
                if (!has_em && dlen < m_shower_split_em_type_max_len) {
                    em_retype.push_back(ns);
                    if (dbg)
                        std::fprintf(stderr,
                            "SHOWER_SPLIT emtype shower=%d part=%zu new_start=%d nseg=%zu "
                            "len_cm=%.2f parent_pdg=%d\n",
                            pr91_seg_display_id(ss), g, pr91_seg_display_id(root_sg),
                            comp.size(), dlen / units::cm, shower->get_particle_type());
                }
            }
            ++n_peeled;
            if (dbg) {
                // The forward-seeding check the choice above exists to satisfy:
                // does the daughter's own charge sit DOWNSTREAM of its start?
                WireCell::Vector c(0, 0, 0); double wq = 0;
                for (const auto& sg : comp)
                    for (const auto& f : sg->fits())
                        if (f.valid()) { const double w = std::max(f.dQ, 0.0); c = c + f.point * w; wq += w; }
                double fwd = 0;
                if (wq > 0 && root_d >= 0) {
                    WireCell::Vector body = c / wq - v;
                    WireCell::Vector st(0, 0, 0);
                    for (const auto& f : root_sg->fits())
                        if (f.valid()) { st = f.point - v; break; }
                    if (body.magnitude() > 0 && st.magnitude() > 0)
                        fwd = body.norm().dot(st.norm());
                }
                std::fprintf(stderr,
                    "SHOWER_SPLIT peel shower=%d part=%zu new_start=%d nseg=%zu "
                    "conn=%d root_vtx_cm=%.2f fwd=%.3f\n",
                    pr91_seg_display_id(ss), g, pr91_seg_display_id(root_sg), comp.size(),
                    conn, root_d >= 0 ? root_d / units::cm : -1.0, fwd);
            }
        }
    }

    // A shed removes members from a parent without producing a daughter, so it
    // still needs the map rebuild and the kinematics refresh below -- the
    // parent's charge and PDG are now different.  `daughters` is empty in that
    // case and every loop over it no-ops.
    if (n_peeled || n_shed) {
        for (auto& ns : daughters) showers.insert(ns);
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                           map_vertex_to_shower, used_shower_clusters);
        // The splitter owns its refresh: the only recompute after this point is
        // recompute_shower_kine_charge_final, which is knob-gated and no-ops in
        // the production configuration.
        calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                    graph, track_fitter, dv, particle_data, recomb_model);

        // doc pr/139 sec 25 -- applied HERE, not at the peel, because
        // calculate_kinematics copies the start segment's pdg verbatim into
        // data.particle_type (PRShower.cxx:1617, :1738 -- the defect doc pr/122
        // named) and would overwrite an earlier write.  This fixes the REPORTED
        // type, which is what NeutrinoKinematics reads for the rest-mass term
        // and what the pi0 finders read.  It does NOT re-derive kine_charge,
        // which keeps the value computed under the start segment's typing; on a
        // sub-1.2 cm object that is a 1-3 MeV difference and is left alone
        // rather than re-running the whole refresh.
        for (auto& ns : em_retype) {
            if (!ns) continue;
            const int was = ns->get_particle_type();
            ns->set_particle_type(11);
            ns->set_flags(ShowerFlags::kShower);
            ++n_em_retyped;
            if (dbg)
                std::fprintf(stderr, "SHOWER_SPLIT emtype APPLIED shower_id=%d %d -> 11\n",
                             ns->get_shower_id(), was);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr139 shower_split: em-retype daughter shower_id={} pdg {} -> 11",
                ns->get_shower_id(), was);
        }

        // doc pr/139 P1.4 -- re-home the orphan daughter.  The owner's scan
        // says a cut that leaves an orphan has not finished the job, five
        // times ("keep, but this should be part of the earlier EM shower
        // cluster"), and again on 2026-08-31: "We should re-sit them into the
        // nearby showers."
        //
        // A FORK of pr125 satellite_absorb (:9864 in the shipped tree), not a
        // shared helper -- CLAUDE.md M10: that pass is production ON and stays
        // byte-for-byte untouched.  Three constraints satellite_absorb's own
        // predicates do not carry:
        //   (1) never re-home into the PARENT -- that silently undoes the
        //       split.  The owner says the EARLIER shower, i.e. a third object.
        //   (2) host must be EM-typed and out-charge the daughter (the
        //       satellite_absorb host rule).
        //   (3) the sibling-merge trap: if the nearest EM shower is the pi0
        //       PARTNER, this merges the two gammas and destroys the pair the
        //       split just enabled.  The census is the referee; a census below
        //       35 on this arm means checking the moved events for exactly
        //       that before anything else.
        // Knob off => the block never runs => byte-identical.
        if (m_shower_split_rehome && !daughters.empty()) {
            auto rh_points = [&](const ShowerPtr& sh, size_t cap) {
                std::vector<WireCell::Point> out;
                std::vector<const PR::Fit*> fp;
                std::unordered_set<size_t> seen;
                for (auto edesc : ordered_edges(*sh, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (!sg || !sg->descriptor_valid()) continue;
                    if (!seen.insert(graph[sg->get_descriptor()].index).second) continue;
                    for (const auto& f : sg->fits()) if (f.valid()) fp.push_back(&f);
                }
                const size_t stride = fp.size() > cap ? (fp.size() + cap - 1) / cap : 1;
                for (size_t i = 0; i < fp.size(); i += stride) out.push_back(fp[i]->point);
                return out;
            };
            // Deterministic host order -- the same (cluster id, segment id)
            // key the candidate loop above uses.
            std::vector<ShowerPtr> host_order(showers.begin(), showers.end());
            std::sort(host_order.begin(), host_order.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
                auto* sa = a->start_segment().get();
                auto* sb = b->start_segment().get();
                if (!sa || !sb) return sa < sb;
                const int ca = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
                const int cb = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
                if (ca != cb) return ca < cb;
                return sa->id() < sb->id();
            });
            std::set<Shower*> rh_is_dau;          // membership only
            for (const auto& d : daughters) if (d) rh_is_dau.insert(d.get());
            // Host clouds once, in host_order, before anything merges.
            std::vector<std::vector<WireCell::Point>> host_pts(host_order.size());
            std::vector<char> host_ok(host_order.size(), 0);
            for (size_t hi = 0; hi < host_order.size(); ++hi) {
                const ShowerPtr& h = host_order[hi];
                if (!h || rh_is_dau.count(h.get())) continue;
                if (std::abs(h->get_particle_type()) != 11) continue;
                host_pts[hi] = rh_points(h, 400);
                host_ok[hi] = !host_pts[hi].empty();
            }
            std::vector<std::pair<ShowerPtr, ShowerPtr>> rh_plan;   // (host, daughter)
            std::set<Shower*> rh_claimed;                            // membership only
            for (size_t di = 0; di < daughters.size(); ++di) {
                const ShowerPtr& dau = daughters[di];
                if (!dau || !showers.count(dau)) continue;
                const std::vector<WireCell::Point> dpts = rh_points(dau, 400);
                if (dpts.empty()) continue;
                ShowerPtr best; double best_gap = -1.0;
                for (size_t hi = 0; hi < host_order.size(); ++hi) {
                    if (!host_ok[hi]) continue;
                    const ShowerPtr& h = host_order[hi];
                    if (h == dau || h == dau_parent[di]) continue;
                    if (rh_claimed.count(h.get())) continue;
                    if (!(h->get_kine_charge() > dau->get_kine_charge())) continue;
                    double g = -1.0;
                    for (const auto& a : dpts)
                        for (const auto& b : host_pts[hi]) {
                            const double t = (a - b).magnitude();
                            if (g < 0 || t < g) g = t;
                        }
                    if (g >= 0 && g <= m_shower_split_rehome_gap && (best_gap < 0 || g < best_gap)) {
                        best_gap = g; best = h;
                    }
                }
                if (dbg)
                    std::fprintf(stderr,
                        "SHOWER_SPLIT rehome dau=%d parent=%d ke_mev=%.2f host=%d gap_cm=%.2f verdict=%s\n",
                        pr91_seg_display_id(dau->start_segment()),
                        dau_parent[di] ? pr91_seg_display_id(dau_parent[di]->start_segment()) : -1,
                        dau->get_kine_charge() / units::MeV,
                        best ? pr91_seg_display_id(best->start_segment()) : -1,
                        best_gap >= 0 ? best_gap / units::cm : -1.0,
                        best ? "REHOME" : "orphan");
                if (best) { rh_plan.emplace_back(best, dau); rh_claimed.insert(dau.get()); }
            }
            int n_rehomed = 0;
            std::set<Shower*> rh_recalc;                             // membership only
            for (auto& hd : rh_plan) {
                pr93_probe_absorb_splice("split_rehome", hd.first, hd.second);
                hd.first->add_shower(*hd.second);
                showers.erase(hd.second);
                rh_recalc.insert(hd.first.get());
                ++n_rehomed;
            }
            if (n_rehomed) {
                for (const auto& h : host_order) {
                    if (!h || !rh_recalc.count(h.get())) continue;
                    // Deliberately NO update_particle_type -- the same minimal
                    // absorb pr125 satellite_absorb settled on.
                    h->calculate_kinematics(particle_data, recomb_model,
                                            m_shower_endpoint_exclude_start_vertex,
                                            m_shower_endpoint_skip_orphan_vtx);
                    h->set_kine_charge(cal_kine_charge(h, m_charge_2d_u, m_charge_2d_v,
                                                       m_charge_2d_w, m_map_apa_ch_plane_wires,
                                                       track_fitter, dv));
                    h->set_flag_kinematics(true);
                }
                update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                                   map_vertex_to_shower, used_shower_clusters);
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr139 shower_split: {} daughter(s) re-homed; {} shower(s) now",
                    n_rehomed, showers.size());
            }
        }

        SPDLOG_LOGGER_DEBUG(s_log,
            "pr138 shower_split: {} candidate(s) fired, {} daughter(s) peeled; {} shower(s) now"
            " (pr139: {} shared-member refusal(s), {} impact veto(es), {} all-shared shed(s),"
            " {} em-retyped daughter(s))",
            n_fired, n_peeled, showers.size(), n_shared_refused, n_vetoed, n_shed, n_em_retyped);
    }
}

void PatternAlgorithms::em_collinear_merge(IndexedShowerSet& showers,
    ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower,
    VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters,
    TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
    const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    if (m_em_collinear_merge_deg <= 0) return;
    std::vector<ShowerPtr> pool;
    for (auto sh : showers) {
        if (!sh) continue;
        if (std::abs(sh->get_particle_type()) != 11) continue;
        pool.push_back(sh);
    }
    // Leading-first; deterministic tie-break by the stable per-run shower id.
    std::stable_sort(pool.begin(), pool.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
        const double ka = a->get_kine_charge(), kb = b->get_kine_charge();
        if (ka != kb) return ka > kb;
        return a->get_shower_id() < b->get_shower_id();
    });
    std::set<ShowerPtr> gone;  // lookup-only, never iterated
    bool did = false;
    for (size_t i = 0; i < pool.size(); ++i) {
        ShowerPtr host = pool[i];
        if (gone.count(host)) continue;
        if (host->get_kine_charge() < m_em_collinear_merge_min_host) continue;
        if (host->get_total_length() < 3 * units::cm) continue;
        const WireCell::Point hs = host->get_start_point();
        const WireCell::Vector axis = shower_cal_dir_3vector(*host, hs, 15 * units::cm);
        const double am = axis.magnitude();
        if (am <= 0) continue;
        std::vector<ShowerPtr> taken;
        for (size_t j = i + 1; j < pool.size(); ++j) {
            ShowerPtr frag = pool[j];
            if (gone.count(frag)) continue;
            auto [fsv, fct] = frag->get_start_vertex_and_type();
            if (fct == 1) continue;
            const WireCell::Point fs = frag->get_start_point();
            const WireCell::Vector v(fs.x() - hs.x(), fs.y() - hs.y(), fs.z() - hs.z());
            const double vm = v.magnitude();
            if (vm <= 0.5 * units::cm || vm > m_em_collinear_merge_dis) continue;
            const double ang = std::acos(std::clamp(v.dot(axis) / (vm * am), -1.0, 1.0)) / M_PI * 180.0;
            if (ang >= m_em_collinear_merge_deg) continue;
            if (frag->get_total_length() > 10 * units::cm) {
                const WireCell::Vector fax = shower_cal_dir_3vector(*frag, fs, 15 * units::cm);
                const double fm = fax.magnitude();
                if (fm > 0) {
                    const double a2 = std::acos(std::clamp(fax.dot(axis) / (fm * am), -1.0, 1.0)) / M_PI * 180.0;
                    if (a2 >= 30) {
                        if (pr93_absorb_dbg())
                            std::fprintf(stderr, "SHOWER_ABSORB EMCOL reject=continuation host=%d frag=%d ang=%.1f a2=%.1f\n",
                                         pr132_pi0_shid(host), pr132_pi0_shid(frag), ang, a2);
                        continue;
                    }
                }
            }
            if (pr93_absorb_dbg())
                std::fprintf(stderr, "SHOWER_ABSORB EMCOL take host=%d Eh=%.1f frag=%d Ef=%.1f ang=%.1f dis=%.1f fct=%d\n",
                             pr132_pi0_shid(host), host->get_kine_charge() / units::MeV,
                             pr132_pi0_shid(frag), frag->get_kine_charge() / units::MeV,
                             ang, vm / units::cm, fct);
            taken.push_back(frag);
        }
        for (auto& frag : taken) {
            pr93_probe_absorb_splice("em_collinear_merge", host, frag);
            host->add_shower(*frag);
            showers.erase(frag);
            gone.insert(frag);
            did = true;
        }
        if (!taken.empty()) {
            host->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
            double kine_charge = cal_kine_charge(host, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
            host->set_kine_charge(kine_charge);
            host->set_flag_kinematics(true);
        }
    }
    if (did)
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                          map_vertex_to_shower, used_shower_clusters);
}

// doc pr/132 round 6 (K17) -- EM shower start back-extension.  The round-5
// root-cause finding (132 doc sec 12.3): the reco EM start is often DEEP in
// the true gamma -- the round-6 start ledger measures 5 hand gammas with the
// reco start 9-51 cm DOWNSTREAM of the scanner's em_start_correction
// (105946 +11.9, 169626 +51.4, 342199 +9.3, 347129 +11.7, 76346 +30.2 cm)
// -- and the missing upstream charge, where it was reconstructed at all,
// lives in DETACHED EM fragments that sit BESIDE/BEHIND the reco start
// (41031 at 158.8 deg off the forward axis), unreachable by any forward
// cone.  This seat back-projects from each EM host's start along -axis and
// absorbs detached EM fragments whose endpoints fall inside the tube
// (perpendicular distance < m_em_backext_perp, upstream projection in
// (0, m_em_backext_len]); a conn-2/3 host's start is then re-seated to the
// upstream-most absorbed point (conn-1 starts stay vertex-anchored, and
// vertices are never touched -- the owner's scope rule).  Hadronic/track
// fragments (pdg != 11) never qualify: an incoming cosmic or track cannot
// be eaten.  perp 0 (default) = off = byte-identical.
void PatternAlgorithms::em_start_backext(IndexedShowerSet& showers,
    ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower,
    VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters,
    TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
    const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    if (m_em_backext_perp <= 0) return;
    std::vector<ShowerPtr> pool;
    for (auto sh : showers) {
        if (!sh) continue;
        if (std::abs(sh->get_particle_type()) != 11) continue;
        pool.push_back(sh);
    }
    std::stable_sort(pool.begin(), pool.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
        const double ka = a->get_kine_charge(), kb = b->get_kine_charge();
        if (ka != kb) return ka > kb;
        return a->get_shower_id() < b->get_shower_id();
    });
    std::set<ShowerPtr> gone;  // lookup-only
    bool did = false;
    for (size_t i = 0; i < pool.size(); ++i) {
        ShowerPtr host = pool[i];
        if (gone.count(host)) continue;
        if (host->get_total_length() < 3 * units::cm) continue;
        const WireCell::Point hs = host->get_start_point();
        const WireCell::Vector axis = shower_cal_dir_3vector(*host, hs, 15 * units::cm);
        const double am = axis.magnitude();
        if (am <= 0) continue;
        // upstream = -axis
        const WireCell::Vector up(-axis.x() / am, -axis.y() / am, -axis.z() / am);
        std::vector<std::pair<double, ShowerPtr>> taken;  // (upstream reach, frag)
        for (size_t j = 0; j < pool.size(); ++j) {
            if (j == i) continue;
            ShowerPtr frag = pool[j];
            if (gone.count(frag) || frag == host) continue;
            auto [fsv, fct] = frag->get_start_vertex_and_type();
            if (fct == 1) continue;
            double best_reach = -1;
            bool ok = false;
            for (const WireCell::Point& fp : {frag->get_start_point(), frag->get_end_point()}) {
                const WireCell::Vector v(fp.x() - hs.x(), fp.y() - hs.y(), fp.z() - hs.z());
                const double proj = v.dot(up);
                if (proj <= 0 || proj > m_em_backext_len) continue;
                const double perp2 = v.magnitude() * v.magnitude() - proj * proj;
                if (perp2 > m_em_backext_perp * m_em_backext_perp) continue;
                ok = true;
                if (proj > best_reach) best_reach = proj;
            }
            if (!ok) continue;
            // v2 CONTINUATION guard (round-6 r6be v1 postmortem): the tube
            // passes through the pi0 decay-vertex convergence region, where
            // the PARTNER gamma also starts -- v1 swallowed true gamma-2s
            // wholesale (census 31 -> 21 exact).  A fragment long enough to
            // carry its own direction must CONTINUE the host axis; a partner
            // gamma diverges by the opening angle and is vetoed.
            if (frag->get_total_length() > 10 * units::cm) {
                const WireCell::Vector fax = shower_cal_dir_3vector(*frag, frag->get_start_point(), 15 * units::cm);
                const double fm = fax.magnitude();
                if (fm > 0) {
                    const double a2 = std::acos(std::clamp(std::abs(fax.dot(axis)) / (fm * am), 0.0, 1.0)) / M_PI * 180.0;
                    if (a2 >= 30) {
                        if (pr93_absorb_dbg())
                            std::fprintf(stderr, "SHOWER_ABSORB EMBACK reject=continuation host=%d frag=%d a2=%.1f\n",
                                         pr132_pi0_shid(host), pr132_pi0_shid(frag), a2);
                        continue;
                    }
                }
            }
            if (pr93_absorb_dbg())
                std::fprintf(stderr, "SHOWER_ABSORB EMBACK take host=%d Eh=%.1f frag=%d Ef=%.1f reach=%.1f fct=%d\n",
                             pr132_pi0_shid(host), host->get_kine_charge() / units::MeV,
                             pr132_pi0_shid(frag), frag->get_kine_charge() / units::MeV,
                             best_reach / units::cm, fct);
            taken.push_back({best_reach, frag});
        }
        if (taken.empty()) continue;
        double far_reach = 0;
        WireCell::Point new_start = hs;
        for (auto& [reach, frag] : taken) {
            // the upstream-most endpoint of the absorbed set becomes the new start
            for (const WireCell::Point& fp : {frag->get_start_point(), frag->get_end_point()}) {
                const WireCell::Vector v(fp.x() - hs.x(), fp.y() - hs.y(), fp.z() - hs.z());
                const double proj = v.dot(up);
                const double perp2 = v.magnitude() * v.magnitude() - proj * proj;
                if (proj > far_reach && perp2 <= m_em_backext_perp * m_em_backext_perp) {
                    far_reach = proj;
                    new_start = fp;
                }
            }
            pr93_probe_absorb_splice("em_start_backext", host, frag);
            host->add_shower(*frag);
            showers.erase(frag);
            gone.insert(frag);
            did = true;
        }
        host->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
        double kine_charge = cal_kine_charge(host, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
        host->set_kine_charge(kine_charge);
        host->set_flag_kinematics(true);
        auto [hsv, hct] = host->get_start_vertex_and_type();
        if (hct != 1 && far_reach > 0) {
            if (pr93_absorb_dbg())
                std::fprintf(stderr, "SHOWER_ABSORB EMBACK reseat host=%d reach=%.1f\n",
                             pr132_pi0_shid(host), far_reach / units::cm);
            host->set_start_point(new_start);
        }
    }
    if (did)
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                          map_vertex_to_shower, used_shower_clusters);
}

// doc pr/132 round 4 -- body of the WCT_PI0_NCVTX_DEBUG tapes (docstring at
// pr132_pi0_ncvtx_dbg above).  Deterministic: IndexedShowerSet iteration is
// id-ordered; kd-tree queries are read-only.
static void pr132_pi0_ncvtx_probe(const IndexedShowerSet& showers,
                                  const IndexedShowerSet& pi0_showers,
                                  VertexPtr main_vertex,
                                  ShowerVertexMap& map_vertex_in_shower)
{
    if (!main_vertex) return;
    const WireCell::Point mvp = main_vertex->fit().valid()
        ? main_vertex->fit().point : main_vertex->wcpt().point;
    {
        auto it = map_vertex_in_shower.find(main_vertex);
        std::fprintf(stderr, "PI0_NCVTX mainvtx=%d in_shower_map=%d insh=%d\n",
                     pr91_vtx_display_id(main_vertex),
                     it != map_vertex_in_shower.end() ? 1 : 0,
                     it != map_vertex_in_shower.end() ? pr132_pi0_shid(it->second) : -1);
    }
    for (auto shower : showers) {
        if (!shower) continue;
        const double kq = shower->get_kine_charge();
        if (kq < 50 * WireCell::units::MeV) continue;
        auto pc = shower->get_pcloud("associate_points");
        const size_t n = pc ? pc->npoints() : 0;
        auto [sv, ct] = shower->get_start_vertex_and_type();
        const double d_start = (shower->get_start_point() - mvp).magnitude();
        double d_assoc = -1;
        double f_back = -1;
        if (n >= 10) {
            auto [dd, pp] = shower_get_closest_point(*shower, mvp, "associate_points");
            d_assoc = dd;
            double cx = 0, cy = 0, cz = 0;
            for (size_t i = 0; i < n; ++i) {
                auto p = pc->point3d(i);
                cx += p.x(); cy += p.y(); cz += p.z();
            }
            WireCell::Vector axis(cx / n - mvp.x(), cy / n - mvp.y(), cz / n - mvp.z());
            const double am = axis.magnitude();
            if (am > 0) {
                size_t nb = 0;
                for (size_t i = 0; i < n; ++i) {
                    auto p = pc->point3d(i);
                    const double d = (p.x()-mvp.x())*axis.x() + (p.y()-mvp.y())*axis.y() + (p.z()-mvp.z())*axis.z();
                    if (d < 0) ++nb;
                }
                f_back = (double)nb / (double)n;
            }
        }
        const int pdg = shower->get_particle_type();
        std::fprintf(stderr, "PI0_NCVTX sh=%d E=%.1f ct=%d pdg=%d d_assoc=%.1f d_start=%.1f f_back=%.3f n=%zu inpi0=%d\n",
                     pr132_pi0_shid(shower), kq / WireCell::units::MeV, ct, pdg,
                     d_assoc / WireCell::units::cm, d_start / WireCell::units::cm,
                     f_back, n, pi0_showers.count(shower) ? 1 : 0);
    }
    for (auto shower : pi0_showers) {
        if (!shower) continue;
        auto [sv, ct] = shower->get_start_vertex_and_type();
        if (!sv) continue;
        const WireCell::Point vp = sv->fit().valid() ? sv->fit().point : sv->wcpt().point;
        const double d_vtx_start = (shower->get_start_point() - vp).magnitude();
        auto [d_fit, p_fit]     = shower_get_closest_point(*shower, vp, "fit");
        auto [d_assoc, p_assoc] = shower_get_closest_point(*shower, vp, "associate_points");
        const double gap = (d_fit >= 0 && d_assoc >= 0) ? (p_fit - p_assoc).magnitude() : -1;
        std::fprintf(stderr, "PI0_START sh=%d E=%.1f ct=%d d_vtx_start=%.1f d_vtx_fit=%.1f d_vtx_assoc=%.1f gap=%.1f\n",
                     pr132_pi0_shid(shower), shower->get_kine_charge() / WireCell::units::MeV, ct,
                     d_vtx_start / WireCell::units::cm, d_fit / WireCell::units::cm,
                     d_assoc / WireCell::units::cm, gap / WireCell::units::cm);
    }
}

// doc pr/132 round 9 (K19): the gamma-ray back-projection NC vertex
// proposer.  Round 8 measured the without-vertex finder structurally unable
// to reach the owner's NC wrong-vertex class (doc sec 15): it runs AFTER
// path 1 consumes the gammas, and its start-point rays are corrupted by the
// vertex defect itself.  This runs BEFORE path 1, fires ONLY on the owner's
// reco-checkable NC signature -- the main vertex sits INSIDE a shower, not
// at its start (map_vertex_in_shower, the round-3 verdicts' "nu vertex
// inside one arm of an EM shower") -- and proposes the pi0 decay point as
// the closest-approach midpoint of two shower AXES back-projected from
// their starts (the em_geom pi0_backproject construction of the hand-label
// pass).  Acceptance needs: miss distance < knob, both crossings BEHIND the
// starts (within 120 cm), the pair mass at the proposed point inside the
// (100,160) window (K1 offset honored), shift 1-100 cm.  On accept the
// P2 precedent applies verbatim: move main_vertex->fit().point, register
// the pair (pio_type 2 = the displaced-vertex family the BDTs know), seat
// both starts on the moved vertex (conn 2), refresh the maps; path 1 then
// skips the paired showers (knob-gated skips in its pools).  0 (default)
// = off = byte-identical.
// doc pr/133 (K20): a shower-topology object typed muon (pdg +-13) can be a
// mistyped pi0 gamma -- the r9-verdict typing census found 3 of the 132 hand
// gammas in this class (specimens: SBND 18255-348691 sh 51080, 79.8 MeV,
// 5.4 cm from the owner's click; 166870 g2, 38.6 MeV), and unlike the
// A5/211 tag (K7) no readmit path exists.  Knob-on, such an object enters
// the pi0 pools iff it carries a real shower flag (topology/trajectory --
// for pdg 13 the get_flag_shower() pdg-11 clause is moot) and its start
// segment is not part of the long muon (the pseudo-shower wrap).  Accepted
// members are re-stamped EM by the K7/K8 machinery.  Knob off (default) =>
// returns false => every call site keeps the legacy exclusion =>
// byte-identical.
bool PatternAlgorithms::pi0_mu_shower_admit(const ShowerPtr& shower, const IndexedSegmentSet& segments_in_long_muon) const
{
    if (!m_pi0_admit_mu_showers) return false;
    if (!shower) return false;
    auto ss = shower->start_segment();
    if (ss && segments_in_long_muon.find(ss) != segments_in_long_muon.end()) {
        if (pr132_pi0_dbg())
            std::fprintf(stderr, "PI0_PAIR K20 mu-reject sh=%d why=longmu\n", pr132_pi0_shid(shower));
        return false;
    }
    // v2: the target class carries NO shower flag -- get_flag_shower() being
    // false is exactly why it got typed 13 (dbg2 tape: 348691 sh 51080 and
    // 166870 sh 85045 both refuse why=notflag).  Fall back to the file's own
    // shower-ish-muon idiom (three sites, e.g. the :3922 center-point cut):
    // short AND direction-weak.  A genuine traveling muon is long or
    // dir-strong; a wrapped long muon is caught above.
    // doc 77 round 4: pr/141 M2 exposed this 40 cm bound as a knob and
    // measured it INERT (seg_dir_weak, not length, is what refuses the
    // mu-typed candidates), so it is back to the file's own literal.
    const bool shower_ish = shower->get_flag_shower() ||
        (shower->get_total_length() < 40 * units::cm && ss && seg_dir_weak(ss));
    if (!shower_ish) {
        if (pr132_pi0_dbg())
            std::fprintf(stderr, "PI0_PAIR K20 mu-reject sh=%d why=trackish len=%.1f\n",
                         pr132_pi0_shid(shower), shower->get_total_length() / units::cm);
        return false;
    }
    if (pr132_pi0_dbg())
        std::fprintf(stderr, "PI0_PAIR K20 mu-admit sh=%d E=%.1f\n",
                     pr132_pi0_shid(shower), shower->get_kine_charge() / units::MeV);
    return true;
}

void PatternAlgorithms::id_pi0_backproject_vertex(int& acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    if (m_pi0_bp_miss <= 0) return;  // knob off = byte-identical
    if (!main_vertex) return;
    // The NC signature gate (owner diagnostic, reco-checkable).  v2 -- the
    // r9 smoke measured in_shower_map=0 on 180801/259542: the wrong vertex
    // is seated as the shower's START vertex, which map_vertex_in_shower
    // excludes by construction.  Use the P2-style segment test instead: a
    // pure-shower multi-prong vertex (>= 2 segments, ALL in showers, none
    // in a long muon).  CC track vertices and 1-prong primary electrons
    // are excluded.
    std::set<ShowerPtr> nc_sig_hosts;   // doc pr/133 K21; membership tests only (never iterated -- pointer-set rule)
    {
        std::vector<SegmentPtr> mv_segs;
        auto vd = main_vertex->get_descriptor();
        for (auto eit : sorted_out_edges(vd, graph)) {
            SegmentPtr seg = graph[eit].segment;
            if (seg) mv_segs.push_back(seg);
        }
        if (m_pi0_nc_sig_angle > 0) {
            // doc pr/133 K21 -- gate v4, the owner signature (2026-08-30):
            // "a) we are inside an EM shower b) there is at least another
            // major object [that] looks like [an EM] shower, and the
            // direction is not aligned ... we can assume the 2 gamma
            // (1 pi0) case."  ALL main-vertex prongs must sit inside
            // showers (the 406125 lesson: a p+mu track vertex passed v3),
            // none in the long muon; the hosting shower(s) are recorded and
            // the pairing below is restricted to host x non-aligned other.
            if (mv_segs.empty()) return;
            for (auto sg : mv_segs) {
                auto hit = map_segment_in_shower.find(sg);
                if (hit == map_segment_in_shower.end()) return;   // a track prong => not inside an EM shower
                if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()) return;
                nc_sig_hosts.insert(hit->second);
            }
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp gate=v4 vtx=%d n_host=%zu\n",
                             pr91_vtx_display_id(main_vertex), nc_sig_hosts.size());
        }
        else {
        if (mv_segs.size() < 2) return;
        // v3 (r9 smoke round 2): 180801/259542 carry a non-shower stub
        // prong, so "ALL in showers" over-excludes.  Use exactly the P2
        // gate2 semantics extended to any prong count: >= 1 main-vertex
        // segment in a shower, none in a long muon.  The tight acceptance
        // (miss cap, back-projection behind both starts, the (100,160)
        // window -- far tighter than P2's |m-125|<60) carries the safety;
        // the rounds-1-4 ADVERSE P2 masses (85.2/75.5) are window-excluded.
        bool any_shower = false;
        for (auto sg : mv_segs) {
            if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) any_shower = true;
            if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()) return;
        }
        if (!any_shower) return;
        }
    }

    const WireCell::Point vtx_pt = main_vertex->fit().valid()
        ? main_vertex->fit().point : main_vertex->wcpt().point;

    // Candidate gammas: EM, >= 3 cm, >= 20 MeV (the K13 floor constant --
    // a tiny partner must not relocate the vertex).  IndexedShowerSet
    // iterates in graph-index order: deterministic.
    std::vector<ShowerPtr> cands;
    std::map<ShowerPtr, WireCell::Vector> cdir;   // lookup-only
    std::map<ShowerPtr, WireCell::Vector> edir;   // v2: end-anchored local dir, lookup-only
    // doc pr/133 K21 v2 (owner 2026-08-30 "go ahead"): inside the signature
    // the partner floor drops to the knob (76346's TRUE gamma2 is 5.0 MeV,
    // below the 20 MeV floor -- the proposer paired a wrong partner).
    // 0 (default) or v3 mode = the legacy 20 MeV floor, byte-identical.
    const double cand_floor = (m_pi0_nc_sig_angle > 0 && m_pi0_nc_floor > 0)
        ? m_pi0_nc_floor : 20 * units::MeV;
    for (auto sh : showers) {
        // doc pr/133 K20: the mu-admit escape (default off = legacy skip).
        if (std::abs(sh->get_particle_type()) == 13 &&
            !pi0_mu_shower_admit(sh, segments_in_long_muon)) continue;
        if (sh->get_total_length() < 3 * units::cm) continue;
        if (sh->get_kine_charge() < cand_floor) continue;
        WireCell::Vector d = shower_cal_dir_3vector(*sh, sh->get_start_point(), 15 * units::cm);
        if (d.magnitude() <= 0) continue;
        cands.push_back(sh);
        cdir[sh] = d * (1.0 / d.magnitude());
        // doc pr/133 K21 v2: the 76346 lesson -- the wrong NC vertex is
        // SEATED as a gamma's start, so the reco orientation can be
        // inverted (76346: the reco END is the owner's truth click).  Build
        // the end-anchored ray too; the pairing tries both orientations.
        if (m_pi0_nc_sig_angle > 0) {
            WireCell::Vector de = shower_cal_dir_3vector(*sh, sh->get_end_point(), 15 * units::cm);
            if (de.magnitude() > 0) edir[sh] = de * (1.0 / de.magnitude());
        }
    }
    if (cands.size() < 2) return;
    // doc pr/133 K21 v2: the owner gate (b) made explicit -- "at least
    // another major object [that] looks like [an EM] shower".  With the
    // floor lowered, majorness (>= 20 MeV) must still hold for at least
    // one non-host candidate.
    if (m_pi0_nc_sig_angle > 0) {
        bool major_other = false;
        for (auto sh : cands)
            if (!nc_sig_hosts.count(sh) && sh->get_kine_charge() >= 20 * units::MeV) major_other = true;
        if (!major_other) return;
    }
    if (pr132_pi0_dbg())
        std::fprintf(stderr, "PI0_PAIR P0 bp gate vtx=%d n_cand=%zu\n",
                     pr91_vtx_display_id(main_vertex), cands.size());

    const double mass_offset = m_pi0_mass_offset;  // K1
    double best_abs = 1e9, best_mass = 0;
    double best_miss = 1e9;   // doc pr/133 K21 v2 ranking (v4 mode only)
    double best_ang = 0;      // doc pr/133 K21 v2.2; opening angle at accept, for the post-absorb mass refresh
    ShowerPtr b1 = nullptr, b2 = nullptr;
    WireCell::Point best_x;

    // doc pr/134 K22 (owner 2026-08-30: 116962 "still not good, the nu
    // vertex is not at the right place"): the six-fragment gamma problem.
    // Fragment-level rays cannot find the true crossing (each fragment's
    // local dir is a fragment artifact), and fragment-level charges put the
    // true crossing's mass out of the window while an in-window FALSE pair
    // (21072 x 55030, m=125.8) wins.  Complex-level pairing instead:
    //   host complex  = the nc_sig_hosts (co-seated at the wrong main),
    //   partner rays  = each non-host candidate's start/end anchor with the
    //                   local dir computed over the MERGED non-host cloud,
    //   masses        = merged charges (host sum x assigned-partner sum).
    // Under merged charges the 116962 false pair leaves the window
    // (125.8 -> 166.1) and the true crossing enters it (33.1 -> 146.4) --
    // the merge kills the wrong fire and enables the right one with the
    // SAME legacy window.  Offline final-sim (doc pr/134): 116962 fires at
    // (-70.1,-79.7,111.7) 44.3 cm upstream; 76346 reproduces its v2.2
    // vertex exactly; 21073/282979/175896/176502/275011/499423 all NO FIRE.
    // false (default) = the v2.1/v2.2 fragment-level path, byte-identical.
    bool mrg_mode = (m_pi0_nc_sig_angle > 0 && m_pi0_nc_frag_merge);
    std::vector<ShowerPtr> mrg_hosts, mrg_partners;   // graph-index order (never iterate the pointer set)
    std::vector<ShowerPtr> best_pmembers;             // partner-complex members of the winning pair
    if (mrg_mode) {
        for (auto sh : showers)
            if (nc_sig_hosts.count(sh)) mrg_hosts.push_back(sh);
        for (auto sh : cands)
            if (!nc_sig_hosts.count(sh)) mrg_partners.push_back(sh);
        if (mrg_hosts.empty() || mrg_partners.empty()) return;
        // The EM-host guard (full-arm k2223 v1 ADVERSE, SBND 499423): a
        // 211-typed attached prong passed the v4 gate as a "shower" host
        // and the merged partner charges lifted its false pair into the
        // window (vertex moved 4.7 cm off the hand click).  The owner
        // signature is "we are inside an EM shower": every host-complex
        // member must be EM-typed (pdg +-11; +-13 only through the K20
        // admit).  A non-EM host falls back to the fragment-level loop --
        // exactly the k21v5 behavior the owner accepted (499423: no fire).
        for (auto sh : mrg_hosts) {
            const int hpdg = std::abs(sh->get_particle_type());
            if (hpdg == 11) continue;
            if (hpdg == 13 && pi0_mu_shower_admit(sh, segments_in_long_muon)) continue;
            mrg_mode = false;
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp mrg refuse host=%d pdg=%d\n",
                             pr132_pi0_shid(sh), sh->get_particle_type());
            break;
        }
    }
    if (mrg_mode) {
        // ---- doc pr/134 K22 v2 (owner verdicts on the v1 Bee, 2026-08-30):
        // "cluster 58035 should be part of cluster 21072" and "the vertex
        // should be much upstream of both".  Three v2 rules, all derived
        // and verified offline before implementation (doc sec 9):
        //  (1) HOST-COLLINEAR PROMOTION: a partner candidate whose own PCA
        //      axis is within 25 deg of the host-complex PCA axis AND whose
        //      centroid sits within 15 deg of the host BACK-ray from the
        //      seat is the host gamma's upstream conversion stem -- it
        //      joins the host complex (116962: 58035 at 16.7/1.0 deg;
        //      67056 refused by the ray test at 60.7 deg).
        //  (2) PARTNER PCA RAYS: a multi-fragment partner gamma's direction
        //      only emerges from the merged cloud -- add the two
        //      PCA-axis rays (anchored at the extremal points) of the
        //      >= 35 MeV partners (sub-35 satellites would twist the axis;
        //      they stay assignment/PF material).
        //  (3) MERGE-MODE WINDOW FLOOR at m > 65 (delta > -60): the
        //      properly-grouped 116962 pair sits at m=79 -- the
        //      under-collection class (reco mass far below label, pr/132
        //      sec B) -- while every measured collateral crossing is
        //      <= 56.2 (282979 PCA) or >= 166 (fragment fires).
        auto cloud_points = [&graph](const std::vector<ShowerPtr>& members, std::vector<WireCell::Point>& out) {
            for (const auto& msh : members) {
                const auto& view = msh->view_graph();
                for (auto edesc : ordered_edges(*msh, graph)) {
                    SegmentPtr seg = view[edesc].segment;
                    if (!seg) continue;
                    for (const auto& fit : seg->fits()) out.push_back(fit.point);
                }
            }
        };
        auto points_pca = [](const std::vector<WireCell::Point>& pts, WireCell::Point& cen, WireCell::Vector& axis) -> bool {
            const size_t n = pts.size();
            if (n < 3) return false;
            double cx = 0, cy = 0, cz = 0;
            for (const auto& p : pts) { cx += p.x(); cy += p.y(); cz += p.z(); }
            cx /= n; cy /= n; cz /= n;
            cen = WireCell::Point(cx, cy, cz);
            double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
            for (const auto& p : pts) {
                const double dx = p.x()-cx, dy = p.y()-cy, dz = p.z()-cz;
                xx += dx*dx; xy += dx*dy; xz += dx*dz; yy += dy*dy; yz += dy*dz; zz += dz*dz;
            }
            // power iteration for the principal axis; deterministic seed =
            // the largest-variance coordinate axis
            double vx = (xx >= yy && xx >= zz) ? 1 : 0;
            double vy = (yy > xx && yy >= zz) ? 1 : 0;
            double vz = (vx == 0 && vy == 0) ? 1 : 0;
            for (int it = 0; it < 60; ++it) {
                const double nx = xx*vx + xy*vy + xz*vz;
                const double ny = xy*vx + yy*vy + yz*vz;
                const double nz = xz*vx + yz*vy + zz*vz;
                const double nm = std::sqrt(nx*nx + ny*ny + nz*nz);
                if (nm <= 0) return false;
                vx = nx/nm; vy = ny/nm; vz = nz/nm;
            }
            axis = WireCell::Vector(vx, vy, vz);
            return true;
        };
        {
            std::vector<WireCell::Point> hpts;
            cloud_points(mrg_hosts, hpts);
            WireCell::Point hcen; WireCell::Vector haxis;
            if (points_pca(hpts, hcen, haxis)) {
                // orient the host axis from the seat INTO the cloud
                const WireCell::Vector seat2cen(hcen.x()-vtx_pt.x(), hcen.y()-vtx_pt.y(), hcen.z()-vtx_pt.z());
                if (seat2cen.dot(haxis) < 0) haxis = haxis * -1.0;
                std::vector<ShowerPtr> keep;
                for (auto sh : mrg_partners) {
                    std::vector<WireCell::Point> fpts;
                    cloud_points({sh}, fpts);
                    WireCell::Point fcen; WireCell::Vector faxis;
                    bool promote = false;
                    if (points_pca(fpts, fcen, faxis)) {
                        const double a_ax = std::acos(std::clamp(std::abs(faxis.dot(haxis)), 0.0, 1.0)) * 180.0 / M_PI;
                        const WireCell::Vector disp(fcen.x()-vtx_pt.x(), fcen.y()-vtx_pt.y(), fcen.z()-vtx_pt.z());
                        const double dm = disp.magnitude();
                        const double a_ray = dm > 0
                            ? std::acos(std::clamp(disp.dot(haxis * -1.0) / dm, -1.0, 1.0)) * 180.0 / M_PI : 180.0;
                        promote = (a_ax < 25.0 && a_ray < 15.0);
                        if (promote && pr132_pi0_dbg())
                            std::fprintf(stderr, "PI0_PAIR P0 bp mrg promote sh=%d a_ax=%.1f a_ray=%.1f\n",
                                         pr132_pi0_shid(sh), a_ax, a_ray);
                    }
                    if (promote) mrg_hosts.push_back(sh);
                    else keep.push_back(sh);
                }
                mrg_partners.swap(keep);
                if (mrg_partners.empty()) return;   // whole event is one complex
            }
        }
        // Merged-cloud local dir: the multi-shower analogue of
        // shower_cal_dir_3vector (unweighted fit-point centroid within
        // 15 cm), accumulated across the member showers in vector order.
        auto mrg_dir = [&graph](const std::vector<ShowerPtr>& members, const WireCell::Point& p) -> WireCell::Vector {
            const double R2 = std::pow(15 * units::cm, 2);
            double sx = 0, sy = 0, sz = 0; int ncount = 0;
            for (const auto& msh : members) {
                const auto& view = msh->view_graph();
                for (auto edesc : ordered_edges(*msh, graph)) {
                    SegmentPtr seg = view[edesc].segment;
                    if (!seg) continue;
                    for (const auto& fit : seg->fits()) {
                        const double dx = fit.point.x() - p.x();
                        const double dy = fit.point.y() - p.y();
                        const double dz = fit.point.z() - p.z();
                        if (dx*dx + dy*dy + dz*dz < R2) { sx += fit.point.x(); sy += fit.point.y(); sz += fit.point.z(); ++ncount; }
                    }
                }
            }
            if (!ncount) return WireCell::Vector(0, 0, 0);
            WireCell::Vector dir(sx/ncount - p.x(), sy/ncount - p.y(), sz/ncount - p.z());
            const double mag = dir.magnitude();
            return mag > 0 ? dir * (1.0 / mag) : WireCell::Vector(0, 0, 0);
        };
        // Per-shower whole-cloud centroid (for the partner-charge assignment).
        auto sh_centroid = [&graph](const ShowerPtr& msh, WireCell::Point& out) -> bool {
            double sx = 0, sy = 0, sz = 0; int ncount = 0;
            const auto& view = msh->view_graph();
            for (auto edesc : ordered_edges(*msh, graph)) {
                SegmentPtr seg = view[edesc].segment;
                if (!seg) continue;
                for (const auto& fit : seg->fits()) { sx += fit.point.x(); sy += fit.point.y(); sz += fit.point.z(); ++ncount; }
            }
            if (!ncount) return false;
            out = WireCell::Point(sx/ncount, sy/ncount, sz/ncount);
            return true;
        };
        struct MrgRay { ShowerPtr sh; WireCell::Point a; WireCell::Vector d; };
        std::vector<MrgRay> hrays, prays;
        for (auto sh : mrg_hosts)
            for (int o = 0; o < 2; ++o) {
                const WireCell::Point a = o ? sh->get_end_point() : sh->get_start_point();
                const WireCell::Vector dd = mrg_dir(mrg_hosts, a);
                if (dd.magnitude() > 0) hrays.push_back({sh, a, dd});
            }
        for (auto sh : mrg_partners)
            for (int o = 0; o < 2; ++o) {
                const WireCell::Point a = o ? sh->get_end_point() : sh->get_start_point();
                const WireCell::Vector dd = mrg_dir(mrg_partners, a);
                if (dd.magnitude() > 0) prays.push_back({sh, a, dd});
            }
        // v2 rule (2): the two partner-cloud PCA rays (>= 35 MeV members).
        {
            std::vector<ShowerPtr> pmajor;
            ShowerPtr phead = nullptr;
            for (auto sh : mrg_partners)
                if (sh->get_kine_charge() >= 35 * units::MeV) {
                    pmajor.push_back(sh);
                    if (!phead || sh->get_kine_charge() > phead->get_kine_charge()) phead = sh;
                }
            std::vector<WireCell::Point> ppts;
            if (!pmajor.empty()) cloud_points(pmajor, ppts);
            WireCell::Point pcen; WireCell::Vector paxis;
            if (phead && points_pca(ppts, pcen, paxis)) {
                double tmin = 1e30, tmax = -1e30;
                WireCell::Point amin, amax;
                for (const auto& p : ppts) {
                    const double t = (p.x()-pcen.x())*paxis.x() + (p.y()-pcen.y())*paxis.y() + (p.z()-pcen.z())*paxis.z();
                    if (t < tmin) { tmin = t; amin = p; }
                    if (t > tmax) { tmax = t; amax = p; }
                }
                if (tmax > tmin) {
                    prays.push_back({phead, amin, paxis});
                    prays.push_back({phead, amax, paxis * -1.0});
                }
            }
        }
        double Eh = 0;
        ShowerPtr host_head = nullptr;
        for (auto sh : mrg_hosts) {
            Eh += sh->get_kine_charge();
            if (!host_head || sh->get_kine_charge() > host_head->get_kine_charge()) host_head = sh;
        }
        if (pr132_pi0_dbg())
            std::fprintf(stderr, "PI0_PAIR P0 bp mrg hosts=%zu partners=%zu Eh=%.1f head=%d\n",
                         mrg_hosts.size(), mrg_partners.size(), Eh / units::MeV,
                         pr132_pi0_shid(host_head));
        for (const auto& hr : hrays) {
            for (const auto& pr : prays) {
                const double axang = std::acos(std::clamp(hr.d.dot(pr.d), -1.0, 1.0)) * 180.0 / M_PI;
                if (axang < m_pi0_nc_sig_angle) continue;
                const WireCell::Vector w0v(hr.a.x()-pr.a.x(), hr.a.y()-pr.a.y(), hr.a.z()-pr.a.z());
                const double bb = hr.d.dot(pr.d);
                const double den = 1.0 - bb * bb;
                if (den < 1e-6) continue;
                const double dvv = hr.d.dot(w0v), evv = pr.d.dot(w0v);
                const double u1 = (bb * evv - dvv) / den;
                const double u2 = (evv - bb * dvv) / den;
                if (u1 > 1 * units::cm || u1 < -120 * units::cm) continue;
                if (u2 > 1 * units::cm || u2 < -120 * units::cm) continue;
                const WireCell::Point c1(hr.a.x()+u1*hr.d.x(), hr.a.y()+u1*hr.d.y(), hr.a.z()+u1*hr.d.z());
                const WireCell::Point c2(pr.a.x()+u2*pr.d.x(), pr.a.y()+u2*pr.d.y(), pr.a.z()+u2*pr.d.z());
                const double miss = (c1 - c2).magnitude();
                if (miss > m_pi0_bp_miss) continue;
                const WireCell::Point x((c1.x()+c2.x())/2, (c1.y()+c2.y())/2, (c1.z()+c2.z())/2);
                const double shift = (x - vtx_pt).magnitude();
                if (shift < 1 * units::cm || shift > 100 * units::cm) continue;
                const WireCell::Vector r1(hr.a.x()-x.x(), hr.a.y()-x.y(), hr.a.z()-x.z());
                const WireCell::Vector r2(pr.a.x()-x.x(), pr.a.y()-x.y(), pr.a.z()-x.z());
                if (r1.magnitude() < 1 * units::cm || r2.magnitude() < 1 * units::cm) continue;
                const double ang = std::acos(std::clamp(
                    r1.dot(r2) / (r1.magnitude() * r2.magnitude()), -1.0, 1.0));
                // Partner-charge assignment: the anchor fragment always
                // belongs to the partner (its charge is the conversion
                // stem); every other non-host candidate joins the side
                // whose vertex ray its cloud centroid sits closer to.
                double Ep = pr.sh->get_kine_charge();
                std::vector<ShowerPtr> pmem{pr.sh};
                for (auto osh : mrg_partners) {
                    if (osh == pr.sh) continue;
                    WireCell::Point cen;
                    if (!sh_centroid(osh, cen)) continue;
                    const WireCell::Vector vs(cen.x()-x.x(), cen.y()-x.y(), cen.z()-x.z());
                    const double vsm = vs.magnitude();
                    if (vsm <= 0) continue;
                    const double aH = std::acos(std::clamp(vs.dot(r1) / (vsm * r1.magnitude()), -1.0, 1.0));
                    const double aP = std::acos(std::clamp(vs.dot(r2) / (vsm * r2.magnitude()), -1.0, 1.0));
                    if (aP <= aH) { Ep += osh->get_kine_charge(); pmem.push_back(osh); }
                }
                const double m = std::sqrt(4 * Eh * Ep * std::pow(std::sin(ang / 2.0), 2));
                const double delta = m - 135 * units::MeV + mass_offset;
                if (pr132_pi0_dbg())
                    std::fprintf(stderr, "PI0_PAIR P0 bp try3 h=%d p=%d miss=%.1f u1=%.1f u2=%.1f shift=%.1f Ep=%.1f npm=%zu m=%.1f\n",
                                 pr132_pi0_shid(hr.sh), pr132_pi0_shid(pr.sh),
                                 miss / units::cm, u1 / units::cm, u2 / units::cm,
                                 shift / units::cm, Ep / units::MeV, pmem.size(), m / units::MeV);
                // v2 rule (3): merge-mode floor -60 (m > 65), upper edge legacy.
                if (delta >= 35 * units::MeV || delta <= -60 * units::MeV) continue;
                if (miss < best_miss - 1e-9 ||
                    (std::abs(miss - best_miss) <= 1e-9 && std::abs(delta) < best_abs)) {
                    best_miss = miss;
                    best_abs = std::abs(delta);
                    best_mass = m;
                    best_ang = ang;
                    b1 = host_head; b2 = pr.sh; best_x = x;
                    best_pmembers = pmem;
                }
            }
        }
    }
    else
    for (size_t i = 0; i < cands.size(); ++i) {
        for (size_t j = i + 1; j < cands.size(); ++j) {
            ShowerPtr s1 = cands[i], s2 = cands[j];
            const WireCell::Point p1 = s1->get_start_point(), p2 = s2->get_start_point();
            const WireCell::Vector d1 = cdir[s1], d2 = cdir[s2];
            // doc pr/133 K21: assume the 2-gamma (1 pi0) case -- one
            // member must be a vertex-hosting shower, and the two axes must
            // be misaligned by more than the knob (split-shower fragments
            // are near-collinear).  0 (default) = v3 behavior, unrestricted.
            //
            // v2 (owner 2026-08-30 "go ahead" after the r10-verdict scan):
            // end-agnostic rays -- both orientations of each gamma tried
            // (the 76346 inversion), ranked by geometric MISS (the
            // flat-mass finding: the window cannot validate a vertex), the
            // mass recorded but NOT gating.  Verified offline: end-flipped
            // gamma1 x the 5-MeV gamma2 reproduces the owner's hand vertex
            // to 0.5 cm (doc pr/133 sec 9).  v3 math below is untouched.
            if (m_pi0_nc_sig_angle > 0) {
                if (!nc_sig_hosts.count(s1) && !nc_sig_hosts.count(s2)) continue;
                // v2 co-seat veto (v3 of the veto): two candidates whose
                // START POINTS coincide are fragments of one complex
                // (116962: 21072 x 21073 x 21070 all start at the same
                // point on the wrong main; self-paired under flips at miss
                // 7.2) -- a pi0 pair must come from two distinct objects.
                // NOTE the first cut used seat-POINTER equality, which is
                // wrong: ct2/ct3 showers all seat at the main vertex object
                // (it killed the 76346 rescue 14059 x 49048, full-arm k21v2).
                if ((s1->get_start_point() - s2->get_start_point()).magnitude()
                    < 0.5 * units::cm) continue;
                for (int o1 = 0; o1 < 2; ++o1) {
                    if (o1 && !edir.count(s1)) continue;
                    const WireCell::Point a1 = o1 ? s1->get_end_point() : p1;
                    const WireCell::Vector dd1 = o1 ? edir[s1] : d1;
                    for (int o2 = 0; o2 < 2; ++o2) {
                        if (o2 && !edir.count(s2)) continue;
                        const WireCell::Point a2 = o2 ? s2->get_end_point() : p2;
                        const WireCell::Vector dd2 = o2 ? edir[s2] : d2;
                        const double axang = std::acos(std::clamp(dd1.dot(dd2), -1.0, 1.0)) * 180.0 / M_PI;
                        if (axang < m_pi0_nc_sig_angle) continue;
                        const WireCell::Vector w0v(a1.x()-a2.x(), a1.y()-a2.y(), a1.z()-a2.z());
                        const double bb = dd1.dot(dd2);
                        const double den = 1.0 - bb * bb;
                        if (den < 1e-6) continue;
                        const double dv = dd1.dot(w0v), ev = dd2.dot(w0v);
                        const double u1 = (bb * ev - dv) / den;
                        const double u2 = (ev - bb * dv) / den;
                        if (u1 > 1 * units::cm || u1 < -120 * units::cm) continue;
                        if (u2 > 1 * units::cm || u2 < -120 * units::cm) continue;
                        const WireCell::Point c1(a1.x()+u1*dd1.x(), a1.y()+u1*dd1.y(), a1.z()+u1*dd1.z());
                        const WireCell::Point c2(a2.x()+u2*dd2.x(), a2.y()+u2*dd2.y(), a2.z()+u2*dd2.z());
                        const double miss = (c1 - c2).magnitude();
                        if (miss > m_pi0_bp_miss) continue;
                        const WireCell::Point x((c1.x()+c2.x())/2, (c1.y()+c2.y())/2, (c1.z()+c2.z())/2);
                        const double shift = (x - vtx_pt).magnitude();
                        if (shift < 1 * units::cm || shift > 100 * units::cm) continue;
                        const WireCell::Vector r1(a1.x()-x.x(), a1.y()-x.y(), a1.z()-x.z());
                        const WireCell::Vector r2(a2.x()-x.x(), a2.y()-x.y(), a2.z()-x.z());
                        if (r1.magnitude() < 1 * units::cm || r2.magnitude() < 1 * units::cm) continue;
                        const double ang = std::acos(std::clamp(
                            r1.dot(r2) / (r1.magnitude() * r2.magnitude()), -1.0, 1.0));
                        const double kq1 = s1->get_kine_charge(), kq2 = s2->get_kine_charge();
                        const double m = std::sqrt(4 * kq1 * kq2 * std::pow(std::sin(ang / 2.0), 2));
                        const double delta = m - 135 * units::MeV + mass_offset;
                        if (pr132_pi0_dbg())
                            std::fprintf(stderr, "PI0_PAIR P0 bp try2 sh1=%d sh2=%d o1=%d o2=%d miss=%.1f u1=%.1f u2=%.1f shift=%.1f m=%.1f\n",
                                         pr132_pi0_shid(s1), pr132_pi0_shid(s2), o1, o2,
                                         miss / units::cm, u1 / units::cm, u2 / units::cm,
                                         shift / units::cm, m / units::MeV);
                        // v2.1: the window RESTORED as the sanity gate (the
                        // dbg7 tape: every collateral fire -- incl. both
                        // ADVERSE, 21073 off a dead-on click and 282979 --
                        // has a grossly non-pi0 mass 47.5-503.0, while both
                        // owner rescues sit in-window at 121.5/125.8).  The
                        // v1 failure was rays+ranking, never the window;
                        // ranking stays miss-first among in-window combos.
                        if (delta >= 35 * units::MeV || delta <= -25 * units::MeV) continue;
                        if (miss < best_miss - 1e-9 ||
                            (std::abs(miss - best_miss) <= 1e-9 && std::abs(delta) < best_abs)) {
                            best_miss = miss;
                            best_abs = std::abs(delta);
                            best_mass = m;
                            best_ang = ang;
                            b1 = s1; b2 = s2; best_x = x;
                        }
                    }
                }
                continue;   // v4 fully handled; the v3 block below stays v3-only
            }
            // Closest approach of the two axis lines.
            const WireCell::Vector w0(p1.x()-p2.x(), p1.y()-p2.y(), p1.z()-p2.z());
            const double b = d1.dot(d2);
            const double denom = 1.0 - b * b;
            if (denom < 1e-6) continue;                       // parallel axes
            const double d = d1.dot(w0), e = d2.dot(w0);
            const double t1 = (b * e - d) / denom;
            const double t2 = (e - b * d) / denom;
            // BACK-projection: the conversion point lies upstream of both
            // starts (1 cm forward slop for a mis-seated start).
            if (t1 > 1 * units::cm || t1 < -120 * units::cm) continue;
            if (t2 > 1 * units::cm || t2 < -120 * units::cm) continue;
            const WireCell::Point c1(p1.x()+t1*d1.x(), p1.y()+t1*d1.y(), p1.z()+t1*d1.z());
            const WireCell::Point c2(p2.x()+t2*d2.x(), p2.y()+t2*d2.y(), p2.z()+t2*d2.z());
            const double miss = (c1 - c2).magnitude();
            if (miss > m_pi0_bp_miss) continue;
            const WireCell::Point x((c1.x()+c2.x())/2, (c1.y()+c2.y())/2, (c1.z()+c2.z())/2);
            const double shift = (x - vtx_pt).magnitude();
            if (shift < 1 * units::cm || shift > 100 * units::cm) continue;
            const WireCell::Vector r1(p1.x()-x.x(), p1.y()-x.y(), p1.z()-x.z());
            const WireCell::Vector r2(p2.x()-x.x(), p2.y()-x.y(), p2.z()-x.z());
            if (r1.magnitude() < 1 * units::cm || r2.magnitude() < 1 * units::cm) continue;
            const double ang = std::acos(std::clamp(
                r1.dot(r2) / (r1.magnitude() * r2.magnitude()), -1.0, 1.0));
            const double kq1 = s1->get_kine_charge(), kq2 = s2->get_kine_charge();
            const double m = std::sqrt(4 * kq1 * kq2 * std::pow(std::sin(ang / 2.0), 2));
            const double delta = m - 135 * units::MeV + mass_offset;
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp try sh1=%d sh2=%d miss=%.1f t1=%.1f t2=%.1f shift=%.1f m=%.1f in=%d\n",
                             pr132_pi0_shid(s1), pr132_pi0_shid(s2), miss / units::cm,
                             t1 / units::cm, t2 / units::cm, shift / units::cm, m / units::MeV,
                             (delta < 35 * units::MeV && delta > -25 * units::MeV) ? 1 : 0);
            if (delta >= 35 * units::MeV || delta <= -25 * units::MeV) continue;
            if (std::abs(delta) < best_abs) {
                best_abs = std::abs(delta);
                best_mass = m;
                b1 = s1; b2 = s2; best_x = x;
            }
        }
    }
    if (!b1) return;

    if (pr132_pi0_dbg())
        std::fprintf(stderr, "PI0_PAIR P0 bp accept sh1=%d sh2=%d m=%.1f shift=%.1f\n",
                     pr132_pi0_shid(b1), pr132_pi0_shid(b2), best_mass / units::MeV,
                     (best_x - vtx_pt).magnitude() / units::cm);

    // The P2 acceptance, verbatim mechanics (vertex hack + registration).
    main_vertex->fit().point = best_x;
    main_vertex->fit().dQ = 0;
    pi0_showers.insert(b1);
    pi0_showers.insert(b2);
    const int pio_id = acc_segment_id++;
    map_shower_pio_id[b1] = pio_id;
    map_shower_pio_id[b2] = pio_id;
    map_pio_id_mass[pio_id] = {best_mass, 2};   // 2 = displaced-vertex family (BDT pio_type in-distribution)
    map_pio_id_showers[pio_id].push_back(b1);
    map_pio_id_showers[pio_id].push_back(b2);
    // doc pr/133 K20: a mu-typed member admitted by the knob is re-stamped
    // EM BEFORE its kinematics are recomputed (the K7/K8 precedent).
    if (m_pi0_admit_mu_showers) {
        for (auto& rsh : {b1, b2})
            if (std::abs(rsh->get_particle_type()) != 11)
                pi0_restamp_shower_em(rsh, particle_data, recomb_model);
    }
    // doc pr/134 K22: the accepted pair was formed at COMPLEX level -- make
    // the merge real (the K12/K18 accept-time precedent): the other host
    // members fold into the host head, the assigned partner members into
    // the partner anchor, then both gammas' charges and the registered mass
    // refresh from the merged objects at the accept-time opening angle.
    if (mrg_mode) {
        size_t n_mrg = 0;
        for (auto sh : mrg_hosts) {
            if (sh == b1) continue;
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp mrgabsorb sh=%d E=%.1f -> gamma=%d\n",
                             pr132_pi0_shid(sh), sh->get_kine_charge() / units::MeV,
                             pr132_pi0_shid(b1));
            pr93_probe_absorb_splice("nc_frag_merge", b1, sh);
            b1->add_shower(*sh);
            showers.erase(sh);
            ++n_mrg;
        }
        for (auto sh : best_pmembers) {
            if (sh == b2) continue;
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp mrgabsorb sh=%d E=%.1f -> gamma=%d\n",
                             pr132_pi0_shid(sh), sh->get_kine_charge() / units::MeV,
                             pr132_pi0_shid(b2));
            pr93_probe_absorb_splice("nc_frag_merge", b2, sh);
            b2->add_shower(*sh);
            showers.erase(sh);
            ++n_mrg;
        }
        if (n_mrg) {
            for (auto& g : {b1, b2}) {
                double kq = cal_kine_charge(g, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
                g->set_kine_charge(kq);
                g->set_flag_kinematics(true);
            }
            const double m_new = std::sqrt(4 * b1->get_kine_charge() * b2->get_kine_charge()
                                           * std::pow(std::sin(best_ang / 2.0), 2));
            map_pio_id_mass[pio_id] = {m_new, 2};
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp mrgmass n_absorb=%zu m %.1f -> %.1f\n",
                             n_mrg, best_mass / units::MeV, m_new / units::MeV);
        }
    }
    // doc pr/133 K21 v2.2 (owner 2026-08-30): "once the vertex is changed,
    // the electron connecting to the original vertex should be updated, or
    // included in the EM shower ... The other low-energy gammas can be
    // associated with the pi0 gammas."  Post-fire Particle-Flow update,
    // NC-signature mode only:
    //   (1) an object CO-STARTED with a pair gamma (< 1 cm: the old wrong
    //       vertex is seated INSIDE the shower, so its attached small
    //       "electron" shares the exact start point -- 76346's 14058,
    //       116962's 21070/21073) is the same complex: absorbed into that
    //       gamma regardless of energy;
    //   (2) any other LOW-ENERGY (< 35 MeV) shower whose start lies within
    //       the knob cone of a pair gamma's ray from the NEW vertex is a
    //       satellite of that gamma: absorbed too.
    // Long-muon members never touched.  The registered pio mass is then
    // refreshed from the post-absorb charges at the accept-time opening
    // angle (it may leave the window -- the window gated ACCEPTANCE; the
    // PF update is honest charge accounting).  0 (default) = off.
    if (m_pi0_nc_sig_angle > 0 && m_pi0_nc_pf_assoc > 0) {
        std::vector<ShowerPtr> pf_absorb;               // graph-index order
        std::map<ShowerPtr, ShowerPtr> pf_target;       // lookup-only (never iterated)
        WireCell::Vector g1(b1->get_start_point().x()-best_x.x(), b1->get_start_point().y()-best_x.y(), b1->get_start_point().z()-best_x.z());
        WireCell::Vector g2(b2->get_start_point().x()-best_x.x(), b2->get_start_point().y()-best_x.y(), b2->get_start_point().z()-best_x.z());
        const double g1m = g1.magnitude(), g2m = g2.magnitude();
        for (auto sh : showers) {
            if (sh == b1 || sh == b2) continue;
            auto ssl = sh->start_segment();
            if (ssl && segments_in_long_muon.find(ssl) != segments_in_long_muon.end()) continue;
            ShowerPtr tgt = nullptr;
            if ((sh->get_start_point() - b1->get_start_point()).magnitude() < 1 * units::cm) tgt = b1;
            else if ((sh->get_start_point() - b2->get_start_point()).magnitude() < 1 * units::cm) tgt = b2;
            else if (sh->get_kine_charge() < 35 * units::MeV) {
                WireCell::Vector vs(sh->get_start_point().x()-best_x.x(), sh->get_start_point().y()-best_x.y(), sh->get_start_point().z()-best_x.z());
                const double vsm = vs.magnitude();
                if (vsm > 1 * units::cm && g1m > 0 && g2m > 0) {
                    const double a1 = std::acos(std::clamp(vs.dot(g1) / (vsm * g1m), -1.0, 1.0)) * 180.0 / M_PI;
                    const double a2 = std::acos(std::clamp(vs.dot(g2) / (vsm * g2m), -1.0, 1.0)) * 180.0 / M_PI;
                    if (std::min(a1, a2) < m_pi0_nc_pf_assoc) tgt = (a1 <= a2) ? b1 : b2;
                }
            }
            if (tgt) { pf_absorb.push_back(sh); pf_target[sh] = tgt; }
        }
        for (auto sh : pf_absorb) {
            ShowerPtr tgt = pf_target[sh];
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp pfassoc sh=%d E=%.1f -> gamma=%d\n",
                             pr132_pi0_shid(sh), sh->get_kine_charge() / units::MeV,
                             pr132_pi0_shid(tgt));
            pr93_probe_absorb_splice("nc_pf_assoc", tgt, sh);
            tgt->add_shower(*sh);
            showers.erase(sh);
        }
        if (!pf_absorb.empty()) {
            for (auto& g : {b1, b2}) {
                double kq = cal_kine_charge(g, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
                g->set_kine_charge(kq);
                g->set_flag_kinematics(true);
            }
            const double m_new = std::sqrt(4 * b1->get_kine_charge() * b2->get_kine_charge()
                                           * std::pow(std::sin(best_ang / 2.0), 2));
            map_pio_id_mass[pio_id] = {m_new, 2};
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P0 bp pfmass n_absorb=%zu m %.1f -> %.1f\n",
                             pf_absorb.size(), best_mass / units::MeV, m_new / units::MeV);
        }
    }
    b1->set_start_vertex(main_vertex, 2);
    b1->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
    b2->set_start_vertex(main_vertex, 2);
    b2->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                       map_vertex_to_shower, used_shower_clusters);
}

void PatternAlgorithms::id_pi0_with_vertex(int& acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

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
    // doc 77 round 4: pr/141's M1 re-pricing of mu-typed pi0 candidates under
    // the SHOWER recombination+fudge (and M3's length floor on it) is retired
    // -- M1 was a REGRESSION, census 35 -> 34, because it moved the greedy
    // partner choice and broke 166870, the very event K20 was shipped to
    // rescue (doc pr/141 sec 9.1).  Every shower keeps the price its own type
    // gives it, as production always did.
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
    // doc pr/132 round 10: associate-cloud centroids, computed lazily and
    // ONLY under the angle tape (pointer-keyed map, lookup-only).
    std::map<ShowerPtr, WireCell::Point> ang_cent;
    auto get_centroid = [&](ShowerPtr sh) -> WireCell::Point {
        auto it = ang_cent.find(sh);
        if (it != ang_cent.end()) return it->second;
        WireCell::Point c(0, 0, 0);
        auto pc = sh->get_pcloud("associate_points");
        const size_t n = pc ? pc->npoints() : 0;
        if (n > 0) {
            double cx = 0, cy = 0, cz = 0;
            for (size_t i = 0; i < n; ++i) {
                auto p = pc->point3d(i);
                cx += p.x(); cy += p.y(); cz += p.z();
            }
            c = WireCell::Point(cx / n, cy / n, cz / n);
        }
        ang_cent[sh] = c;
        return c;
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
            // doc pr/132 round 9 (K19): proposer-paired showers stay paired.
            if (m_pi0_bp_miss > 0 && pi0_showers.count(shower)) continue;
            // doc pr/99 round 3 (A5): re-typed hadronic showers never enter
            // pi0 pairing; empty set (tag off) => no-op, byte-identical.
            // doc pr/132 round 2 (K7): the knob readmits them -- the A5 dQdx
            // track call is overruled by pairing evidence (accepted members
            // are re-stamped EM below).  false (default) = byte-identical.
            if (!m_pi0_readmit_retyped &&
                m_hadronic_retyped_shower_ids.count(shower->get_shower_id())) continue;
            auto [start_vtx, conn_type] = get_svc(shower);
            // doc pr/132 round 2 (K8): legacy disconnected pool is
            // conn_type==2 only; the knob also admits conn_type==3
            // (other-cluster) showers, which no path-1 pool ever saw
            // (specimen: SBND 18255-47212 g2).  false (default) = legacy.
            const bool ct_pool_ok = (conn_type == 2) ||
                                    (m_pi0_admit_type3 && conn_type == 3);
            // doc pr/133 K20: the mu-admit escape (default off = legacy).
            if (ct_pool_ok && (std::abs(shower->get_particle_type()) != 13 ||
                               pi0_mu_shower_admit(shower, segments_in_long_muon))) {
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

    if (pr132_pi0_dbg())
        std::fprintf(stderr, "PI0_PAIR P1 begin main_vtx=%d n_cand_vtx=%zu n_disc=%zu\n",
                     pr91_vtx_display_id(main_vertex), candidate_vertices.size(),
                     disconnected_showers.size());

    // -- Map shower pairs → masses and candidate vertices -------------------
    // Custom comparator eliminates pointer-address ordering of pair keys.
    std::map<std::pair<ShowerPtr,ShowerPtr>,
             std::vector<std::pair<double,VertexPtr>>,
             decltype(shower_pair_cmp)> map_shower_pair_mass_vertex(shower_pair_cmp);

    // doc pr/134 K24 round-2 fix: showers whose presence in the MAIN-vertex
    // pool exists ONLY because of the relaxed admission (they failed the
    // 30-deg dir test and the crumb bypass there).  Membership lookups only
    // -- never iterated (determinism rule).  Empty when the knob is off.
    std::set<ShowerPtr> main_relax_only;

    // doc pr/134 K24 round-2 fix (part 2, the ambiguity constraint): the
    // owner's rule is "if there is a direction ambiguity, the preference
    // should give to neutrino vertex" -- so the relaxed admission at main is
    // reserved for showers that DO legacy-admit at some other candidate
    // vertex (the ambiguity exists; 105946's 55063 admits at the proton-end
    // vertex).  A shower that admits NOWHERE is a pure relax artifact
    // (54332's 120071/128111: angle 50-176 deg at every vertex) and seeded
    // 19 unadjudicated pi0 groups across 239 events in the v4 arm.
    // Lookups only, never iterated.  Empty when the knob is off.
    std::set<ShowerPtr> main_amb_ok;
    if (m_pi0_prefer_main_vertex) {
        for (auto shower : disconnected_showers) {
            const WireCell::Vector dir1 = map_shower_dir[shower];
            auto [start_vtx, conn_type] = get_svc(shower);
            for (auto cand_vtx : candidate_vertices) {
                if (cand_vtx == main_vertex) continue;
                // Own-start admission counts as ambiguity (105946's 55063
                // seats at its own conversion vertex 48 cm away); the
                // ct1-partner rule at the pair-recording site is what keeps
                // pure artifacts (54332's 120071/128111) from pairing.
                if (start_vtx == cand_vtx) { main_amb_ok.insert(shower); break; }
                const WireCell::Point vtx_pt = cand_vtx->fit().valid()
                    ? cand_vtx->fit().point : cand_vtx->wcpt().point;
                const WireCell::Vector dir2(shower->get_start_point().x() - vtx_pt.x(),
                                            shower->get_start_point().y() - vtx_pt.y(),
                                            shower->get_start_point().z() - vtx_pt.z());
                const double angle = std::acos(std::clamp(
                    dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                const bool crumb_ok = (m_pi0_crumb_assoc_max > 0 &&
                                       get_kq(shower) < m_pi0_crumb_assoc_max);
                if (angle < m_pi0_assoc_angle_deg || crumb_ok) { main_amb_ok.insert(shower); break; }
            }
        }
    }

    // doc pr/132 round 3 (K12): virtual collinear merge (132 doc sec 9.9
    // idea 1).  A second gamma split into collinear detached fragments never
    // forms an in-window pair (specimen SBND 18255-54341: right topology,
    // pair mass 68.7 MeV -- half the charge is in sibling fragments).  When
    // the knob is on, at each candidate vertex the detached (conn_type!=1)
    // pool is greedily re-clustered: descending in kine_charge, each leading
    // fragment VIRTUALLY absorbs every remaining detached shower within the
    // vertex-ray cone (angle < knob, degrees).  The summed charge enters the
    // pair-mass computation only; absorbed fragments stop pairing on their
    // own at that vertex.  Attached (conn_type==1) showers neither host nor
    // get absorbed -- the primary-electron guard -- and the K3
    // attached-partner floor still tests the leading fragment's OWN charge
    // (deliberately stricter).  Shower objects are untouched unless the pair
    // is ACCEPTED; then the fragments are truly absorbed (P2 absorb
    // precedent, below).  0 (default) = off = byte-identical.  The maps are
    // pointer-keyed but lookup-only -- never iterated.
    std::map<VertexPtr, std::map<ShowerPtr, std::vector<ShowerPtr>>> cm_frags;
    std::map<VertexPtr, std::map<ShowerPtr, double>> cm_vkq;
    std::map<VertexPtr, std::set<ShowerPtr>> cm_absorbed;
    bool cm_did_absorb = false;
    // doc pr/134 K23: pairs accepted in this call, for the post-loop PF
    // satellite absorb.  Only filled knob-on (empty vector = byte-identical).
    std::vector<std::tuple<int, ShowerPtr, ShowerPtr, VertexPtr>> p1_pairs;

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
                // doc pr/132 round 9 (K19): a shower already paired by the
                // back-projection proposer stays paired.  Knob-off:
                // pi0_showers is empty at path-1 entry, unreachable.
                if (m_pi0_bp_miss > 0 && pi0_showers.count(shower)) continue;
                // doc pr/99 round 3 (A5): re-typed hadronic showers never enter
                // pi0 pairing; empty set (tag off) => no-op, byte-identical.
                // doc pr/132 round 2 (K7): knob readmits (see the pool above).
                if (!m_pi0_readmit_retyped &&
                    m_hadronic_retyped_shower_ids.count(shower->get_shower_id())) continue;
                auto [start_vtx, conn_type] = get_svc(shower);
                // doc pr/133 K20: the mu-admit escape (default off = legacy).
                if (conn_type == 1 && (std::abs(shower->get_particle_type()) != 13 ||
                                       pi0_mu_shower_admit(shower, segments_in_long_muon))) {
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
                // doc pr/132 round 2 (K9): a low-energy crumb's 15-cm PCA
                // direction is noise (specimen: SBND 18255-71872 g2, a
                // 23 MeV crumb at 73 deg) -- below the threshold the angle
                // test is skipped; the mass window, greedy selection and the
                // K3 attached-partner guard control the combinatorics.
                // 0 (default) = off = byte-identical.
                const bool crumb_ok = (m_pi0_crumb_assoc_max > 0 &&
                                       get_kq(shower) < m_pi0_crumb_assoc_max);
                if (pr132_pi0_dbg())
                    std::fprintf(stderr, "PI0_PAIR P1 assoc vtx=%d sh=%d E=%.1f angle=%.1f crumb=%d acc=%d\n",
                                 pr91_vtx_display_id(cand_vtx), pr132_pi0_shid(shower),
                                 get_kq(shower) / units::MeV, angle, crumb_ok ? 1 : 0,
                                 (angle < m_pi0_assoc_angle_deg || crumb_ok) ? 1 : 0);
                // doc pr/134 K24 (owner 2026-08-30: "the pi0 should be at
                // the neutrino vertex ... if there is a direction
                // ambiguity ... the preference should give to neutrino
                // vertex"): at the MAIN vertex only, a disconnected shower
                // enters the pool regardless of the 30-deg dir test.  The
                // shower's internal dir is deficit-biased (measured 50-75
                // deg off the flight line on 105946's 55063 -- it steered
                // the gamma to a proton-end vertex 48 cm away); the
                // conversion displacement (start - vertex) IS the gamma
                // direction, and that is the ray the recorded mass uses.
                // Round-2 fix: relax only for genuinely ambiguous showers
                // (legacy-admitted at some OTHER vertex -- main_amb_ok).
                const bool main_relax = (m_pi0_prefer_main_vertex && cand_vtx == main_vertex &&
                                         main_amb_ok.count(shower));
                const bool legacy_ok = (angle < m_pi0_assoc_angle_deg || crumb_ok);
                if (legacy_ok || main_relax) {  // doc pr/132 K2 + K9; pr/134 K24
                    tmp_showers.push_back(shower);
                    local_dirs[shower] = dir2;
                    // K24 round-2: remember relax-ONLY admissions -- a pair
                    // built on one ranks below any legacy pair at main
                    // (54332/506746: a relax-admitted fake partner displaced
                    // the true pair that had already won at main).
                    if (!legacy_ok) main_relax_only.insert(shower);
                }
            }
        }

        // Sort so (i < j) pairs are in graph-index order regardless of insertion order.
        std::sort(tmp_showers.begin(), tmp_showers.end(), shower_cmp);

        // doc pr/132 round 3 (K12): build this vertex's virtual merge
        // (rationale at the cm_* declarations above).
        if (m_pi0_collinear_merge_deg > 0) {
            std::vector<ShowerPtr> detached;
            for (auto sh : tmp_showers) {
                auto [sv, ct] = get_svc(sh);
                if (ct != 1) detached.push_back(sh);
            }
            // Leading-first; stable_sort on the shower_cmp-ordered list keeps
            // equal-charge ordering deterministic.
            std::stable_sort(detached.begin(), detached.end(),
                             [&](const ShowerPtr& a, const ShowerPtr& b)
                             { return get_kq(a) > get_kq(b); });
            auto& absorbed = cm_absorbed[cand_vtx];
            for (size_t i = 0; i < detached.size(); ++i) {
                ShowerPtr host = detached[i];
                if (absorbed.count(host)) continue;
                const WireCell::Vector hd = local_dirs[host];
                for (size_t j = i + 1; j < detached.size(); ++j) {
                    ShowerPtr frag = detached[j];
                    if (absorbed.count(frag)) continue;
                    const WireCell::Vector fd = local_dirs[frag];
                    const double a = std::acos(std::clamp(
                        hd.dot(fd) / (hd.magnitude() * fd.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    if (a >= m_pi0_collinear_merge_deg) continue;
                    if (!cm_vkq[cand_vtx].count(host)) cm_vkq[cand_vtx][host] = get_kq(host);
                    cm_vkq[cand_vtx][host] += get_kq(frag);
                    cm_frags[cand_vtx][host].push_back(frag);
                    absorbed.insert(frag);
                    if (pr132_pi0_dbg())
                        std::fprintf(stderr, "PI0_PAIR P1 cmerge vtx=%d host=%d frag=%d a=%.1f Eh=%.1f Ef=%.1f Ev=%.1f\n",
                                     pr91_vtx_display_id(cand_vtx), pr132_pi0_shid(host),
                                     pr132_pi0_shid(frag), a, get_kq(host) / units::MeV,
                                     get_kq(frag) / units::MeV, cm_vkq[cand_vtx][host] / units::MeV);
                }
            }
        }

        // Compute pi0 mass for each (i < j) pair.
        // Early-skip if both conn_type==1 avoids the acos+sqrt for ineligible pairs.
        for (size_t i = 0; i < tmp_showers.size(); i++) {
            ShowerPtr sh1 = tmp_showers[i];
            // doc pr/132 K12: an absorbed fragment pairs only via its host.
            if (m_pi0_collinear_merge_deg > 0 && cm_absorbed[cand_vtx].count(sh1)) continue;
            auto [sv1, ct1] = get_svc(sh1);
            WireCell::Vector dir1 = local_dirs[sh1];
            double kq1 = get_kq(sh1);
            if (m_pi0_collinear_merge_deg > 0 && cm_vkq[cand_vtx].count(sh1))
                kq1 = cm_vkq[cand_vtx][sh1];  // doc pr/132 K12: virtual sum

            for (size_t j = i + 1; j < tmp_showers.size(); j++) {
                ShowerPtr sh2 = tmp_showers[j];
                if (m_pi0_collinear_merge_deg > 0 && cm_absorbed[cand_vtx].count(sh2)) continue;
                auto [sv2, ct2] = get_svc(sh2);
                if (ct1 == 1 && ct2 == 1) continue;  // ineligible — skip before expensive ops

                // doc pr/134 K24 round-2 (the legacy-anchor rule): the
                // relaxed admission may COMPLETE a pair whose other member
                // is legacy-admitted at main (105946: 56056 passes the
                // 30-deg test at main, 18.5 deg, and anchors the relaxed
                // 55063), but two relax-only members cannot fabricate a
                // pair between themselves (54332: 120071 x 122091 m=121.6
                // and 16017 x 122091 m=109.8, every member angle-failed at
                // main, displaced the true pair at its true vertex).
                // Empty set when the knob is off = byte-identical.
                if (cand_vtx == main_vertex &&
                    main_relax_only.count(sh1) && main_relax_only.count(sh2)) continue;

                double kq2 = get_kq(sh2);
                if (m_pi0_collinear_merge_deg > 0 && cm_vkq[cand_vtx].count(sh2))
                    kq2 = cm_vkq[cand_vtx][sh2];  // doc pr/132 K12: virtual sum
                WireCell::Vector dir2 = local_dirs[sh2];
                double angle    = std::acos(std::clamp(
                    dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0));
                double mass_pio = std::sqrt(4 * kq1 * kq2 * std::pow(std::sin(angle / 2.0), 2));
                map_shower_pair_mass_vertex[{sh1, sh2}].push_back({mass_pio, cand_vtx});
                if (pr132_pi0_dbg())
                    std::fprintf(stderr, "PI0_PAIR P1 pair sh1=%d sh2=%d ct1=%d ct2=%d vtx=%d E1=%.1f E2=%.1f m=%.1f\n",
                                 pr132_pi0_shid(sh1), pr132_pi0_shid(sh2), ct1, ct2,
                                 pr91_vtx_display_id(cand_vtx), kq1 / units::MeV,
                                 kq2 / units::MeV, mass_pio / units::MeV);
                // doc pr/132 round 10: the opening-angle census tape
                // (docstring at pr132_pi0_angle_dbg).
                if (pr132_pi0_angle_dbg()) {
                    const WireCell::Point c1 = get_centroid(sh1);
                    const WireCell::Point c2 = get_centroid(sh2);
                    const WireCell::Vector r1(c1.x()-vtx_pt.x(), c1.y()-vtx_pt.y(), c1.z()-vtx_pt.z());
                    const WireCell::Vector r2(c2.x()-vtx_pt.x(), c2.y()-vtx_pt.y(), c2.z()-vtx_pt.z());
                    double a_cent = -1, m_cent = -1;
                    if (r1.magnitude() > 0.1 * units::cm && r2.magnitude() > 0.1 * units::cm) {
                        a_cent = std::acos(std::clamp(
                            r1.dot(r2) / (r1.magnitude() * r2.magnitude()), -1.0, 1.0));
                        m_cent = std::sqrt(4 * kq1 * kq2 * std::pow(std::sin(a_cent / 2.0), 2));
                    }
                    std::fprintf(stderr, "PI0_ANGLE vtx=%d sh1=%d sh2=%d ct1=%d ct2=%d E1=%.1f E2=%.1f m_start=%.1f a_start=%.2f m_cent=%.1f a_cent=%.2f\n",
                                 pr91_vtx_display_id(cand_vtx), pr132_pi0_shid(sh1), pr132_pi0_shid(sh2),
                                 ct1, ct2, kq1 / units::MeV, kq2 / units::MeV,
                                 mass_pio / units::MeV, angle / M_PI * 180.0,
                                 m_cent / units::MeV, a_cent >= 0 ? a_cent / M_PI * 180.0 : -1.0);
                }
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
        const double mass_offset = m_pi0_mass_offset;  // doc pr/132 K1; 10 MeV = legacy
        VertexPtr vtx = nullptr;
        double mass_penalty = 0;
        int best_tier = 0;   // doc pr/134 K24; knob-off everything is tier 0 = legacy comparator

        // Hoist start-vertex lookup out of the mass_vtx_vec inner loop.
        for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex) {
            auto [sv1, ct1] = get_svc(shower_pair.first);
            auto [sv2, ct2] = get_svc(shower_pair.second);
            // doc pr/132 K3: nueCC-fake guard.  A recorded pair has at most
            // one conn_type==1 member (ct1==1 && ct2==1 never recorded), and
            // such a pair exists only at that member's own start vertex -- so
            // sv==main_vertex identifies the "primary electron attached at
            // the nu vertex" topology.  Reject the pairing when the detached
            // partner is below the threshold.  Selection-loop only: the
            // T_KINE pio_kine scan above is untouched.  0 (default) = off =
            // byte-identical.
            if (m_pi0_attached_partner_min > 0) {
                const bool veto =
                    (ct1 == 1 && sv1 == main_vertex && get_kq(shower_pair.second) < m_pi0_attached_partner_min) ||
                    (ct2 == 1 && sv2 == main_vertex && get_kq(shower_pair.first)  < m_pi0_attached_partner_min);
                if (veto) {
                    if (pr132_pi0_dbg())
                        std::fprintf(stderr, "PI0_PAIR P1 veto sh1=%d sh2=%d E1=%.1f E2=%.1f\n",
                                     pr132_pi0_shid(shower_pair.first), pr132_pi0_shid(shower_pair.second),
                                     get_kq(shower_pair.first) / units::MeV,
                                     get_kq(shower_pair.second) / units::MeV);
                    continue;
                }
            }
            double tmp_penalty = (ct1 == 2 && ct2 == 2) ? 6 * units::MeV : 0.0;

            for (auto& [mass, candidate_vtx] : mass_vtx_vec) {
                double delta = mass - 135 * units::MeV + mass_offset;
                if (delta >= 35 * units::MeV || delta <= -25 * units::MeV) continue;
                // doc pr/134 K24 (round-2 tier fix): three-tier main-first
                // ranking among in-window candidates.  Tier 2 = at main with
                // BOTH members legacy-admitted there; tier 1 = at main but
                // riding the relaxed admission; tier 0 = non-main.  Higher
                // tier wins outright; the legacy |delta|-with-bonus key
                // decides within a tier.  Rationale: an in-window pair AT
                // the nu vertex beats any non-main pair (47212: same pair at
                // main m=133.7 lost to a pion-stub end m=132.8 by 0.9 MeV),
                // but a relax-admitted partner must NOT displace a pair that
                // already won at main under the legacy tests (54332: true
                // pair 27025x122091 displaced by relax-fake 120071; same on
                // 506746).  Knob off: tier stays 0 -- byte-identical.
                int tier = 0;
                if (m_pi0_prefer_main_vertex && candidate_vtx == main_vertex) {
                    const bool relaxed = main_relax_only.count(shower_pair.first) ||
                                         main_relax_only.count(shower_pair.second);
                    tier = relaxed ? 1 : 2;
                }
                bool better;
                if (tier != best_tier) better = (tier > best_tier);
                else better = (std::abs(delta) - tmp_penalty < std::abs(mass_diff) - mass_penalty);
                if (better) {
                    mass_diff    = delta;
                    mass_penalty = tmp_penalty;
                    mass_save    = mass;
                    shower_1     = shower_pair.first;
                    shower_2     = shower_pair.second;
                    vtx          = candidate_vtx;
                    best_tier    = tier;
                }
            }
        }

        if (mass_diff >= 35 * units::MeV || mass_diff <= -25 * units::MeV) {
            // doc pr/101 (L0): log-only.  Every remaining candidate pair
            // fell outside the (100,160) MeV window -- the incoming-track
            // pion stamp below will not fire for these vertices.
            double best_mass = -1, best_abs = 1e9;
            for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex)
                for (auto& [mass, candidate_vtx] : mass_vtx_vec) {
                    const double d = std::abs(mass - 135 * units::MeV + mass_offset);
                    if (d < best_abs) { best_abs = d; best_mass = mass; }
                }
            SPDLOG_LOGGER_DEBUG(s_log, "pi0 window reject: {} pair(s) left, best mass={:.1f} MeV delta={:.1f}",
                                map_shower_pair_mass_vertex.size(), best_mass / units::MeV,
                                (best_mass - 135 * units::MeV + mass_offset) / units::MeV);
            if (pr132_pi0_dbg())
                for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex) {
                    double bb = 1e9, bm = -1;
                    for (auto& [mass, candidate_vtx] : mass_vtx_vec) {
                        const double d = std::abs(mass - 135 * units::MeV + mass_offset);
                        if (d < bb) { bb = d; bm = mass; }
                    }
                    std::fprintf(stderr, "PI0_PAIR P1 winreject sh1=%d sh2=%d best_m=%.1f\n",
                                 pr132_pi0_shid(shower_pair.first), pr132_pi0_shid(shower_pair.second),
                                 bm / units::MeV);
                }

            // doc pr/132 round 7 (K18): acceptance-aware fragment merge.
            // Three geometric merges (K12 pairing cone, K16 forward cone,
            // K17 back tube) all measured net-negative: near a pi0 decay
            // vertex "own fragment" and "partner gamma" are locally
            // degenerate in charge geometry (132 doc sec 13.3).  The mass
            // constraint is the discriminator geometry lacks, so the merge
            // decision moves INTO the acceptance test: only after every
            // remaining pair fell outside the window, a BELOW-window pair
            // (a charge deficit only lowers mass, so absorbing charge can
            // only raise it back toward the window) may absorb ONE nearby
            // detached fragment, and only when the augmented mass lands
            // in-window.  Nothing in-window is ever destroyed: the rescue
            // runs only when no in-window candidate remains, and accepted
            // pairs already left the map.  The T_KINE pio_kine scan above is
            // untouched.  0 (default) = off = byte-identical.
            bool am_rescued = false;
            if (m_pi0_am_dis > 0) {
                double am_best_abs = 1e9, am_best_mass = 0;
                ShowerPtr am_s1 = nullptr, am_s2 = nullptr, am_host = nullptr, am_frag = nullptr;
                VertexPtr am_vtx = nullptr;
                for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex) {
                    auto [asv1, act1] = get_svc(shower_pair.first);
                    auto [asv2, act2] = get_svc(shower_pair.second);
                    // The K3 nueCC-fake guard holds in the rescue too.
                    if (m_pi0_attached_partner_min > 0) {
                        const bool veto =
                            (act1 == 1 && asv1 == main_vertex && get_kq(shower_pair.second) < m_pi0_attached_partner_min) ||
                            (act2 == 1 && asv2 == main_vertex && get_kq(shower_pair.first)  < m_pi0_attached_partner_min);
                        if (veto) continue;
                    }
                    for (auto& [mass, candidate_vtx] : mass_vtx_vec) {
                        if (mass - 135 * units::MeV + mass_offset > -25 * units::MeV) continue;
                        for (const auto& host : {shower_pair.first, shower_pair.second}) {
                            const double kqh = get_kq(host);
                            if (kqh <= 0) continue;
                            for (auto frag : disconnected_showers) {
                                if (frag == shower_pair.first || frag == shower_pair.second) continue;
                                if (pi0_showers.count(frag) || !showers.count(frag)) continue;
                                // v2 SUBORDINATION guard (round-7 r7am v1
                                // postmortem, SBND 18255-169626): a fragment
                                // never outweighs its host.  Without it the
                                // analytic boost sqrt((kqh+kqf)/kqh) lets a
                                // crumb pair fabricate any mass by absorbing
                                // a large shower (169626: a 1.3 MeV host
                                // "absorbed" the 511 MeV gamma, 5.5->111,
                                // and the accept blocked the without-vertex
                                // finder's true pair).  With it the boost
                                // caps at sqrt(2), so only pairs already at
                                // m >= 71 MeV can reach the window.
                                if (get_kq(frag) >= kqh) continue;
                                const double dis = (frag->get_start_point() - host->get_start_point()).magnitude();
                                if (dis > m_pi0_am_dis) continue;
                                const double m2 = mass * std::sqrt((kqh + get_kq(frag)) / kqh);
                                const double d2 = m2 - 135 * units::MeV + mass_offset;
                                if (pr132_pi0_dbg())
                                    std::fprintf(stderr, "PI0_PAIR P1 amtry sh1=%d sh2=%d host=%d frag=%d dis=%.1f m=%.1f m2=%.1f in=%d\n",
                                                 pr132_pi0_shid(shower_pair.first), pr132_pi0_shid(shower_pair.second),
                                                 pr132_pi0_shid(host), pr132_pi0_shid(frag), dis / units::cm,
                                                 mass / units::MeV, m2 / units::MeV,
                                                 (d2 < 35 * units::MeV && d2 > -25 * units::MeV) ? 1 : 0);
                                if (d2 >= 35 * units::MeV || d2 <= -25 * units::MeV) continue;
                                if (std::abs(d2) < am_best_abs) {
                                    am_best_abs  = std::abs(d2);
                                    am_best_mass = m2;
                                    am_s1 = shower_pair.first;
                                    am_s2 = shower_pair.second;
                                    am_host = host;
                                    am_frag = frag;
                                    am_vtx  = candidate_vtx;
                                }
                            }
                        }
                    }
                }
                if (am_frag) {
                    if (pr132_pi0_dbg())
                        std::fprintf(stderr, "PI0_PAIR P1 amrescue sh1=%d sh2=%d host=%d frag=%d m2=%.1f\n",
                                     pr132_pi0_shid(am_s1), pr132_pi0_shid(am_s2), pr132_pi0_shid(am_host),
                                     pr132_pi0_shid(am_frag), am_best_mass / units::MeV);
                    // Real absorb (the K12 accept-time precedent), then fall
                    // through to the shared accept path below.
                    pr93_probe_absorb_splice("pi0_accept_merge", am_host, am_frag);
                    am_host->add_shower(*am_frag);
                    showers.erase(am_frag);
                    am_host->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
                    double am_kq = cal_kine_charge(am_host, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
                    am_host->set_kine_charge(am_kq);
                    am_host->set_flag_kinematics(true);
                    kqc[am_host] = am_kq;  // keep the charge cache consistent
                    cm_did_absorb = true;  // update_shower_maps at loop exit
                    // The fragment's own recorded pairings are stale now.
                    std::vector<std::pair<ShowerPtr,ShowerPtr>> am_stale;
                    for (auto& [sp, _] : map_shower_pair_mass_vertex)
                        if (sp.first == am_frag || sp.second == am_frag)
                            am_stale.push_back(sp);
                    for (auto& p : am_stale) map_shower_pair_mass_vertex.erase(p);
                    shower_1  = am_s1;
                    shower_2  = am_s2;
                    vtx       = am_vtx;
                    mass_save = am_best_mass;
                    mass_diff = am_best_mass - 135 * units::MeV + mass_offset;
                    am_rescued = true;
                }
            }
            if (!am_rescued) break;
        }

        if (pr132_pi0_dbg())
            std::fprintf(stderr, "PI0_PAIR P1 accept sh1=%d sh2=%d vtx=%d m=%.1f\n",
                         pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2),
                         pr91_vtx_display_id(vtx), mass_save / units::MeV);
        pi0_showers.insert(shower_1);
        pi0_showers.insert(shower_2);

        // doc pr/132 round 2 (K7/K8/K9): an accepted pair member still typed
        // as a track (A5-retyped 211, PID 2212, ...) is re-stamped EM -- the
        // owner's direction: "the code can look at them, and update them to
        // an EM shower in the pi0 reconstruction code".  Gated on the rescue
        // family: all knobs off => no track-typed member can newly appear
        // here beyond legacy, and the stamp does not fire => byte-identical.
        if (m_pi0_readmit_retyped || m_pi0_admit_type3 || m_pi0_crumb_assoc_max > 0 || m_pi0_admit_mu_showers) {
            for (auto& rsh : {shower_1, shower_2})
                if (std::abs(rsh->get_particle_type()) != 11)
                    pi0_restamp_shower_em(rsh, particle_data, recomb_model);
        }

        // doc pr/132 round 3 (K12): the accepted pair's virtual merge becomes
        // real -- absorb the fragments into their host (the P2 accept-time
        // absorb precedent), refresh the host's kinematics/charge, and retire
        // the fragments' own pairings below.  Only reachable knob-on.
        std::set<ShowerPtr> cm_gone;  // lookup-only
        if (m_pi0_collinear_merge_deg > 0) {
            auto vit = cm_frags.find(vtx);
            if (vit != cm_frags.end()) {
                for (const auto& host : {shower_1, shower_2}) {
                    auto hit = vit->second.find(host);
                    if (hit == vit->second.end() || hit->second.empty()) continue;
                    bool did = false;
                    for (const auto& frag : hit->second) {
                        if (pi0_showers.count(frag)) continue;
                        if (pr132_pi0_dbg())
                            std::fprintf(stderr, "PI0_PAIR P1 cmabsorb host=%d frag=%d Ef=%.1f\n",
                                         pr132_pi0_shid(host), pr132_pi0_shid(frag),
                                         frag->get_kine_charge() / units::MeV);
                        pr93_probe_absorb_splice("pi0_collinear_merge", host, frag);
                        host->add_shower(*frag);
                        showers.erase(frag);
                        cm_gone.insert(frag);
                        did = true;
                    }
                    hit->second.clear();
                    if (did) {
                        host->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
                        double kine_charge = cal_kine_charge(host, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv);
                        host->set_kine_charge(kine_charge);
                        host->set_flag_kinematics(true);
                        cm_did_absorb = true;
                    }
                }
            }
        }

        int pio_id = acc_segment_id++;
        g_pr33_audit.f3_pi0_with_vertex++;
        map_shower_pio_id[shower_1]  = pio_id;
        map_shower_pio_id[shower_2]  = pio_id;
        map_pio_id_mass[pio_id]      = {mass_save, 1};
        map_pio_id_showers[pio_id].push_back(shower_1);
        map_pio_id_showers[pio_id].push_back(shower_2);
        if (m_pi0_pf_assoc > 0)   // doc pr/134 K23
            p1_pairs.emplace_back(pio_id, shower_1, shower_2, vtx);

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

        // doc 77 round 4: pr/132's K15 (re-seat an accepted conn-2 start on
        // the ASSOCIATE cloud) is retired.  It fired 22 times with the census
        // identical and no hand-pi0 effect -- 169626's 13.4 cm start bias is
        // charge the shower does not contain, not a wrong cloud choice
        // (pr/132 sec 11.3).  Starts stay on the fit cloud.

        SPDLOG_LOGGER_TRACE(s_log, "examine_showers: Pi0 found with mass: {} MeV with {} MeV + {} MeV",
            mass_save / units::MeV, shower_1->get_kine_charge() / units::MeV, shower_2->get_kine_charge() / units::MeV);

        // Remove all pairs that involve either used shower.
        std::vector<std::pair<ShowerPtr,ShowerPtr>> to_remove;
        for (auto& [sp, _] : map_shower_pair_mass_vertex)
            if (sp.first == shower_1 || sp.first == shower_2 ||
                sp.second == shower_1 || sp.second == shower_2)
                to_remove.push_back(sp);
        // doc pr/132 round 3 (K12): pairings that referenced a now-absorbed
        // fragment (recorded at ANY vertex) are stale too.
        if (!cm_gone.empty())
            for (auto& [sp, _] : map_shower_pair_mass_vertex)
                if (cm_gone.count(sp.first) || cm_gone.count(sp.second))
                    to_remove.push_back(sp);
        for (auto& p : to_remove) map_shower_pair_mass_vertex.erase(p);
    }

    // doc pr/134 K23: the with-vertex port of the v2.2 Particle-Flow update
    // (owner 2026-08-30, on the NC finder: "The other low-energy gammas can
    // be associated with the pi0 gammas").  After the whole selection loop,
    // each accepted pair gamma absorbs its DETACHED low-energy satellites:
    //   - E < 35 MeV (the v2.2 cap), never an attached (conn-1) primary,
    //   - start within (0.5, 120] cm of the gamma start and inside the knob
    //     cone of the gamma's own 15 cm axis (the gamma-ledger sibling
    //     criterion, doc pr/132), assigned to the closer-in-angle gamma,
    //   - never a pair member, never a long-muon member, and never a shower
    //     that still owns a recorded pi0 pairing of its own.
    // Charges then refresh and the registered mass rescales by the charge
    // ratio (the opening angle is unchanged).  0 (default) = off.
    if (m_pi0_pf_assoc > 0 && !p1_pairs.empty()) {
        for (auto& [pf_pid, pf_g1, pf_g2, pf_vtx] : p1_pairs) {
            const double kq1_old = pf_g1->get_kine_charge();
            const double kq2_old = pf_g2->get_kine_charge();
            const WireCell::Point pf_vp = pf_vtx->fit().valid()
                ? pf_vtx->fit().point : pf_vtx->wcpt().point;
            // Cone apex = the pair vertex; axis = the gamma FLIGHT direction
            // (vertex -> conversion start) when the conversion is displaced.
            // The reco shower's own local dir is deficit-biased -- on 105946
            // it points 50+ deg away from the label axis (the missing charge
            // IS the downstream tail), while vertex->start reproduces the
            // label axis to 5 deg.  Attached gammas (start at the vertex)
            // fall back to the 15 cm local dir.
            auto pf_axis = [&](const ShowerPtr& g) -> WireCell::Vector {
                WireCell::Vector ax = g->get_start_point() - pf_vp;
                if (ax.magnitude() > 2 * units::cm) return ax * (1.0 / ax.magnitude());
                ax = shower_cal_dir_3vector(*g, g->get_start_point(), 15 * units::cm);
                const double m = ax.magnitude();
                return m > 0 ? ax * (1.0 / m) : WireCell::Vector(0, 0, 0);
            };
            const WireCell::Vector ax1 = pf_axis(pf_g1);
            const WireCell::Vector ax2 = pf_axis(pf_g2);
            std::vector<ShowerPtr> pf_sat;             // graph-index order
            std::map<ShowerPtr, ShowerPtr> pf_tgt;     // lookup-only
            for (auto sh : showers) {
                if (sh == pf_g1 || sh == pf_g2) continue;
                if (pi0_showers.count(sh)) continue;
                if (sh->get_kine_charge() >= 35 * units::MeV) continue;
                auto [sat_sv, sat_ct] = sh->get_start_vertex_and_type();
                if (sat_ct == 1) continue;             // attached primaries stay
                auto ssl = sh->start_segment();
                if (ssl && segments_in_long_muon.find(ssl) != segments_in_long_muon.end()) continue;
                // No recorded-pairing guard: the selection loop above runs
                // until no acceptable pair remains, so every SURVIVING
                // recorded pairing is one the selector already rejected --
                // it is not a competing pi0 candidate (105946: all five
                // sub-35 MeV satellites own rejected pairings).
                double bang = 1e9; ShowerPtr tgt = nullptr;
                for (const auto& [g, ax] : {std::pair<ShowerPtr, WireCell::Vector>{pf_g1, ax1},
                                            std::pair<ShowerPtr, WireCell::Vector>{pf_g2, ax2}}) {
                    const double axm = ax.magnitude();
                    if (axm <= 0) continue;
                    const WireCell::Vector v = sh->get_start_point() - pf_vp;
                    const double vm = v.magnitude();
                    if (vm <= 0.5 * units::cm || vm > 120 * units::cm) continue;
                    const double ang = std::acos(std::clamp(v.dot(ax) / (vm * axm), -1.0, 1.0)) * 180.0 / M_PI;
                    if (ang < m_pi0_pf_assoc && ang < bang) { bang = ang; tgt = g; }
                }
                if (tgt) { pf_sat.push_back(sh); pf_tgt[sh] = tgt; }
            }
            if (pf_sat.empty()) continue;
            // Charges update ARITHMETICALLY (gamma + absorbed satellites),
            // never via a cal_kine_charge re-derivation: the full-arm k2223
            // v1 measured the refresh re-deriving stored kine_best charges
            // wholesale (47212 g1 221.9 -> 130.1 with no absorb into it;
            // 285567/314838/506746/71372 inflated 60-107 MeV beyond their
            // satellites) -- 8 ledger gammas left OK from the refresh alone.
            std::map<ShowerPtr, double> pf_add;   // lookup-only
            for (auto sh : pf_sat) {
                ShowerPtr tgt = pf_tgt[sh];
                if (pr132_pi0_dbg())
                    std::fprintf(stderr, "PI0_PAIR P1 pfassoc sh=%d E=%.1f -> gamma=%d\n",
                                 pr132_pi0_shid(sh), sh->get_kine_charge() / units::MeV,
                                 pr132_pi0_shid(tgt));
                pf_add[tgt] += sh->get_kine_charge();
                pr93_probe_absorb_splice("p1_pf_assoc", tgt, sh);
                tgt->add_shower(*sh);
                showers.erase(sh);
                kqc.erase(sh);
                svc.erase(sh);
            }
            for (auto& g : {pf_g1, pf_g2}) {
                auto ait = pf_add.find(g);
                if (ait == pf_add.end()) continue;
                const double kq = g->get_kine_charge() + ait->second;
                g->set_kine_charge(kq);
                kqc[g] = kq;
            }
            if (kq1_old > 0 && kq2_old > 0) {
                const double scale = std::sqrt((pf_g1->get_kine_charge() * pf_g2->get_kine_charge())
                                               / (kq1_old * kq2_old));
                auto mit = map_pio_id_mass.find(pf_pid);
                if (mit != map_pio_id_mass.end()) {
                    const double m_old = mit->second.first;
                    mit->second.first = m_old * scale;
                    if (pr132_pi0_dbg())
                        std::fprintf(stderr, "PI0_PAIR P1 pfmass pid=%d n_absorb=%zu m %.1f -> %.1f\n",
                                     pf_pid, pf_sat.size(), m_old / units::MeV,
                                     mit->second.first / units::MeV);
                }
            }
            cm_did_absorb = true;
        }
    }

    // doc pr/132 round 3 (K12) + round 7 (K18): fragments were physically
    // absorbed above -- refresh the membership maps (P2 does the same on
    // accept).  Knobs-off: cm_did_absorb is always false, no call,
    // byte-identical.
    if (cm_did_absorb)
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                          map_vertex_to_shower, used_shower_clusters);

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
                    SPDLOG_LOGGER_DEBUG(s_log, "pi0 incoming stamp: vtx idx={} seg idx={} pdg {} -> 211",
                                        pi0_vtx->get_graph_index(), sg->get_graph_index(), pdg);
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
    // doc pr/132 K5: 2 = legacy prong cap, byte-identical.
    if (main_vertex_segs.size() > static_cast<size_t>(m_pi0_nv_max_prongs)) {
        if (pr132_pi0_dbg())
            std::fprintf(stderr, "PI0_PAIR P2 return=gate1 nsegs=%zu\n", main_vertex_segs.size());
        return;
    }

    if (!main_vertex_segs.empty() && main_vertex_segs.size() <= 2) {
        auto first_seg = main_vertex_segs.front();
        auto last_seg = main_vertex_segs.back();
        
        if ((map_segment_in_shower.find(first_seg) == map_segment_in_shower.end() &&
             map_segment_in_shower.find(last_seg) == map_segment_in_shower.end()) ||
            segments_in_long_muon.find(first_seg) != segments_in_long_muon.end() ||
            segments_in_long_muon.find(last_seg) != segments_in_long_muon.end()) {
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P2 return=gate2 nsegs=%zu\n", main_vertex_segs.size());
            return;
        }
    }
    else if (main_vertex_segs.size() > 2) {
        // doc pr/132 K5 companion: reachable only when m_pi0_nv_max_prongs > 2
        // (GATE1 above returns otherwise), so the legacy path stays
        // byte-identical.  The front/back test generalizes to: at least one
        // main-vertex segment in a shower, none in a long muon.
        bool any_shower = false, any_muon = false;
        for (auto sg : main_vertex_segs) {
            if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) any_shower = true;
            if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()) any_muon = true;
        }
        if (!any_shower || any_muon) {
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P2 return=gate2 nsegs=%zu\n", main_vertex_segs.size());
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
                if (pi0_showers.find(shower) != pi0_showers.end()) {
                    // Path 2 is abandoned the moment ANY main-vertex shower
                    // is already pi0-paired.  doc 77 round 4: pr/132's K14,
                    // which skipped the paired shower and let path 2 run on
                    // the remainder, is retired -- the mechanism worked but
                    // its 239-event ledger was 0 rescues (pr/132 sec 11.2:
                    // path 1 consumes the gammas first, an ordering lock the
                    // retry cannot reach).
                    if (pr132_pi0_dbg()) {
                        int n_left = 0;
                        for (auto sh2 : it->second)
                            if (pi0_showers.find(sh2) == pi0_showers.end() &&
                                sh2->get_start_vertex_and_type().second == 1) n_left++;
                        std::fprintf(stderr, "PI0_PAIR P2 return=already_paired sh=%d unpaired_ct1_left=%d\n",
                                     pr132_pi0_shid(shower), n_left);
                    }
                    return;
                }

                // doc pr/99 round 3 (A5): re-typed hadronic showers never enter
                // pi0 pairing; empty set (tag off) => no-op, byte-identical.
                // doc pr/132 round 2 (K7): knob readmits -- the without-vertex
                // "one member attached at main" requirement then accepts a
                // retyped attached gamma (specimen: SBND 18255-169626, whose
                // true pair m=145.8 formed and died only on good=0).
                if (!m_pi0_readmit_retyped &&
                    m_hadronic_retyped_shower_ids.count(shower->get_shower_id())) continue;
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
        if (flag_return) {
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P2 return=gate3\n");
            return;
        }
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
                // doc pr/133 K20: the mu-admit escape (default off = legacy).
                if ((m_shower_pdg_from_start_segment ? skip_proto : skip_legacy) &&
                    !pi0_mu_shower_admit(shower, segments_in_long_muon)) continue;
            }
            if (shower->get_total_length() < 3 * units::cm) continue;
            if (pi0_showers.find(shower) != pi0_showers.end()) continue;

            WireCell::Point test_p = shower->get_start_point();
            WireCell::Vector dir = shower_cal_dir_3vector(*shower, test_p, 15 * units::cm);
            WireCell::Point p2(test_p.x() + dir.x(), test_p.y() + dir.y(), test_p.z() + dir.z());
            map_shower_ray[shower] = WireCell::Ray(test_p, p2);
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P2 ray src=main sh=%d E=%.1f len=%.1f\n",
                             pr132_pi0_shid(shower), shower->get_kine_charge() / units::MeV,
                             shower->get_total_length() / units::cm);
        }
    }
    
    // Add showers from other vertices
    for (auto& [vtx, shower_set] : map_vertex_to_shower) {
        if (vtx == main_vertex) continue;
        
        for (auto shower : shower_set) {
            // doc pr/99 round 3 (A5): re-typed hadronic showers never enter
            // pi0 pairing; empty set (tag off) => no-op, byte-identical.
            // doc pr/132 round 2 (K7): knob readmits (see id_pi0_with_vertex).
            if (!m_pi0_readmit_retyped &&
                m_hadronic_retyped_shower_ids.count(shower->get_shower_id())) continue;
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
                // doc pr/133 K20: the mu-admit escape (default off = legacy).
                if ((m_shower_pdg_from_start_segment ? skip_proto : skip_legacy) &&
                    !pi0_mu_shower_admit(shower, segments_in_long_muon)) continue;
            }
            if (shower->get_total_length() < 3 * units::cm) continue;
            if (pi0_showers.find(shower) != pi0_showers.end()) continue;

            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            // The pool is conn_type==3 only (prototype
            // NeutrinoID_shower_clustering.h:514).  doc 77 round 4: pr/132's
            // K4, which also admitted conn_type==2, is retired -- the P2
            // family measured ZERO NC rescues over two rounds and its only
            // new acceptances were ADVERSE vertex-draggers (pr/132 sec 6.3).
            const bool ct_ok = (conn_type == 3);
            if (!ct_ok) {
                if (pr132_pi0_dbg())
                    std::fprintf(stderr, "PI0_PAIR P2 rayreject src=other sh=%d ct=%d E=%.1f\n",
                                 pr132_pi0_shid(shower), conn_type,
                                 shower->get_kine_charge() / units::MeV);
                continue;
            }
            
            {
                bool is_shower_seg = shower->start_segment()->flags_any(SegmentFlags::kShowerTrajectory) ||
                                     shower->start_segment()->flags_any(SegmentFlags::kShowerTopology);
                if (!is_shower_seg) {
                    // Also accept segments identified as electron by PID (prototype's get_flag_shower_dQdx)
                    int seg_pdg = 0;
                    if (shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info())
                        seg_pdg = shower->start_segment()->particle_info()->pdg();
                    if (std::abs(seg_pdg) != 11) {
                        if (pr132_pi0_dbg())
                            std::fprintf(stderr, "PI0_PAIR P2 rayreject src=other sh=%d notshower pdg=%d\n",
                                         pr132_pi0_shid(shower), seg_pdg);
                        continue;
                    }
                }
            }
            
            auto [closest_dis, test_p] = shower_get_closest_point(*shower, main_vtx_pt);
            WireCell::Vector dir = shower_cal_dir_3vector(*shower, test_p, 15 * units::cm);
            WireCell::Point p2(test_p.x() + dir.x(), test_p.y() + dir.y(), test_p.z() + dir.z());
            map_shower_ray[shower] = WireCell::Ray(test_p, p2);
            if (pr132_pi0_dbg())
                std::fprintf(stderr, "PI0_PAIR P2 ray src=other sh=%d ct=%d E=%.1f len=%.1f\n",
                             pr132_pi0_shid(shower), conn_type,
                             shower->get_kine_charge() / units::MeV,
                             shower->get_total_length() / units::cm);
        }
    }
    
    if (pr132_pi0_dbg())
        std::fprintf(stderr, "PI0_PAIR P2 pool n_ray=%zu n_good=%zu\n",
                     map_shower_ray.size(), good_showers.size());

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
                    
                    if (angle1 > 25 || angle2 > 25) {
                        if (pr132_pi0_dbg())
                            std::fprintf(stderr, "PI0_PAIR P2 colreject sh1=%d sh2=%d a1=%.1f a2=%.1f\n",
                                         pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2), angle1, angle2);
                        continue;
                    }
                    
                    double angle = std::acos(std::clamp(dir_to_shower1.dot(dir_to_shower2) / (dir_to_shower1.magnitude() * dir_to_shower2.magnitude()), -1.0, 1.0));
                    double mass_pio = std::sqrt(4 * shower_1->get_kine_charge() * shower_2->get_kine_charge() * 
                                               std::pow(std::sin(angle / 2.0), 2));
                    map_shower_pair_mass_point[std::make_pair(shower_1, shower_2)] = std::make_pair(mass_pio, center);
                    if (pr132_pi0_dbg())
                        std::fprintf(stderr, "PI0_PAIR P2 pair br=LL sh1=%d sh2=%d m=%.1f\n",
                                     pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2),
                                     mass_pio / units::MeV);

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
                        if (angle1 > 25 || angle2 > 25) {
                        if (pr132_pi0_dbg())
                            std::fprintf(stderr, "PI0_PAIR P2 colreject sh1=%d sh2=%d a1=%.1f a2=%.1f\n",
                                         pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2), angle1, angle2);
                        continue;
                    }
                        
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
                        if (angle1 > 25 || angle2 > 25) {
                        if (pr132_pi0_dbg())
                            std::fprintf(stderr, "PI0_PAIR P2 colreject sh1=%d sh2=%d a1=%.1f a2=%.1f\n",
                                         pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2), angle1, angle2);
                        continue;
                    }
                    }
                    
                    double angle = std::acos(std::clamp(dir_to_c1.dot(dir_to_c2) / (dir_to_c1.magnitude() * dir_to_c2.magnitude()), -1.0, 1.0));
                    double mass_pio = std::sqrt(4 * shower_1->get_kine_charge() * shower_2->get_kine_charge() * 
                                               std::pow(std::sin(angle / 2.0), 2));
                    map_shower_pair_mass_point[std::make_pair(shower_1, shower_2)] = std::make_pair(mass_pio, center);
                    if (pr132_pi0_dbg())
                        std::fprintf(stderr, "PI0_PAIR P2 pair br=LS sh1=%d sh2=%d m=%.1f\n",
                                     pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2),
                                     mass_pio / units::MeV);

                } else {
                    if (pr132_pi0_dbg())
                        std::fprintf(stderr, "PI0_PAIR P2 bothshort sh1=%d sh2=%d\n",
                                     pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2));
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
        double mass_offset = m_pi0_mass_offset;  // doc pr/132 K1; 10 MeV = legacy
        WireCell::Point vtx_point;

        for (auto& [shower_pair, mass_point] : map_shower_pair_mass_point) {
            // doc pr/132 round 2 (K10): acceptance-quality gate.  A candidate
            // decay point far from the current main vertex is how the round-1
            // ADVERSE movers happened (122660 pulled the vertex 23 cm off the
            // truth click, 171143 5.8 cm) -- skip the pair so the next-best
            // eligible one can win instead.  0 (default) = off.
            if (m_pi0_nv_max_vtx_shift > 0 &&
                (mass_point.second - main_vtx_pt).magnitude() > m_pi0_nv_max_vtx_shift) {
                if (pr132_pi0_dbg())
                    std::fprintf(stderr, "PI0_PAIR P2 shiftrej sh1=%d sh2=%d m=%.1f shift=%.1f\n",
                                 pr132_pi0_shid(shower_pair.first), pr132_pi0_shid(shower_pair.second),
                                 mass_point.first / units::MeV,
                                 (mass_point.second - main_vtx_pt).magnitude() / units::cm);
                continue;
            }
            // doc pr/132 round 4 (K13): path-2 partner floor -- the K3
            // analog for the without-vertex finder.  A big shower plus a
            // tiny partner can land in-window and drag the main vertex with
            // it (round-2 sec 9.8: 116962's unresolved acceptance had a
            // 31 MeV partner); skip the pair when the smaller member is
            // below the floor.  Selection-loop only, T_KINE untouched.
            // 0 (default) = off = byte-identical.
            if (m_pi0_nv_partner_min > 0 &&
                std::min(shower_pair.first->get_kine_charge(),
                         shower_pair.second->get_kine_charge()) < m_pi0_nv_partner_min) {
                if (pr132_pi0_dbg())
                    std::fprintf(stderr, "PI0_PAIR P2 floorrej sh1=%d sh2=%d E1=%.1f E2=%.1f m=%.1f\n",
                                 pr132_pi0_shid(shower_pair.first), pr132_pi0_shid(shower_pair.second),
                                 shower_pair.first->get_kine_charge() / units::MeV,
                                 shower_pair.second->get_kine_charge() / units::MeV,
                                 mass_point.first / units::MeV);
                continue;
            }
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
        
        if (pr132_pi0_dbg()) {
            for (auto& [shower_pair, mass_point] : map_shower_pair_mass_point) {
                const bool good = good_showers.count(shower_pair.first) || good_showers.count(shower_pair.second);
                std::fprintf(stderr, "PI0_PAIR P2 pairsel sh1=%d sh2=%d m=%.1f good=%d\n",
                             pr132_pi0_shid(shower_pair.first), pr132_pi0_shid(shower_pair.second),
                             mass_point.first / units::MeV, good ? 1 : 0);
            }
            std::fprintf(stderr, "PI0_PAIR P2 best sh1=%d sh2=%d m=%.1f shift=%.1f acc=%d\n",
                         pr132_pi0_shid(shower_1), pr132_pi0_shid(shower_2),
                         mass_save / units::MeV,
                         shower_1 ? (vtx_point - main_vtx_pt).magnitude() / units::cm : -1.0,
                         (mass_diff < m_pi0_nv_mass_window && shower_1 && shower_2) ? 1 : 0);
        }

        // If found good pi0, update everything
        // doc pr/132 round 2 (K11): the acceptance-quality gate the data
        // supports -- the round-1 ADVERSE path-2 acceptances are low-mass
        // (122660 m=85.2, 171143 m=75.5) while the good ones sit near 135
        // (396222 m=133.5, 169626 m=138.9); a decay-point shift cap does NOT
        // separate them (shifts 23.0/5.8 vs 14.5/59.6 cm).  60 (default) =
        // legacy window = byte-identical.
        if (mass_diff < m_pi0_nv_mass_window && shower_1 && shower_2) {
            pi0_showers.insert(shower_1);
            pi0_showers.insert(shower_2);

            // doc pr/132 round 2 (K7/K8/K9): re-stamp track-typed accepted
            // members EM (see id_pi0_with_vertex) BEFORE the
            // calculate_kinematics calls below so kinematics reflect EM.
            if (m_pi0_readmit_retyped || m_pi0_admit_type3 || m_pi0_crumb_assoc_max > 0 || m_pi0_admit_mu_showers) {
                for (auto& rsh : {shower_1, shower_2})
                    if (std::abs(rsh->get_particle_type()) != 11)
                        pi0_restamp_shower_em(rsh, particle_data, recomb_model);
            }

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


// doc pr/132 round 2 (K7/K8/K9) -- reverse of the A5 hadronic retype (see
// shower_clustering_with_nv below) for a shower accepted into a pi0 pair:
// the pairing evidence (mass window + geometry) overrules the dQdx-based
// track call.  set_kine_best(0) restores the EM best-energy fall-through
// (PRShower.h get_kine_best returns kenergy_charge when kenergy_best==0),
// undoing apply_hadronic_dqdx_best's override; erasing the id from
// m_hadronic_retyped_shower_ids makes later A5-aware consumers see EM.
// Reached only when a rescue knob is on (call sites are gated).
void PatternAlgorithms::pi0_restamp_shower_em(const ShowerPtr& shower, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model)
{
    if (!shower) return;
    auto ss = shower->start_segment();
    if (ss && ss->has_particle_info() && ss->particle_info()) {
        ss->particle_info()->set_pdg(11);
        ss->particle_info()->set_mass(0.511 * units::MeV);
        if (ss->particle_info()->kinetic_energy() > 0) {
            auto four_momentum = segment_cal_4mom(ss, 11, particle_data, recomb_model, m_mip_dqdx);
            ss->particle_info()->set_four_momentum(four_momentum);
        }
    }
    shower->set_particle_type(11);
    shower->set_kine_best(0);
    m_hadronic_retyped_shower_ids.erase(shower->get_shower_id());
    SPDLOG_LOGGER_DEBUG(s_log, "pi0 restamp EM: shower id={} -> pdg 11 (doc pr/132 round 2 rescue)",
                        shower->get_shower_id());
    if (pr132_pi0_dbg())
        std::fprintf(stderr, "PI0_PAIR restamp sh=%d -> pdg 11\n", pr132_pi0_shid(shower));
}


void PatternAlgorithms::shower_clustering_with_nv(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    // doc pr/99 round 3 (A5): per-event state; nonempty only while
    // m_shower_hadronic_tag fired this event, so the clear is a no-op on the
    // legacy path (byte-identical).
    m_hadronic_retyped_shower_ids.clear();

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

    // doc sbnd_xin/docs/pr/117 round 1 (shower_merge_relax): late fragment
    // consolidation -- two PR::Showers that are one physical shower, rooted
    // away from the main vertex, have no legacy merge pass
    // (examine_merge_showers pairs only conn-1 x conn-2 AT the main vertex).
    // After examine_showers (the last pass that creates or retargets) and
    // before the dedup/detach/ghost family and the pi0 finders, so pairing
    // sees the consolidated set.  false = no call = byte-identical.
    if (m_shower_merge_relax) {
        merge_shower_fragments(graph, showers, main_vertex, map_vertex_in_shower,
                               map_segment_in_shower, map_vertex_to_shower,
                               used_shower_clusters, track_fitter, dv,
                               particle_data, recomb_model);
        t0 = Clock::now();
    }

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

    // doc sbnd_xin/docs/pr/99 round 2 -- shower_ghost_member_drop.  Design
    // block at the m_shower_ghost_* members (NeutrinoPatternBase.h).  The
    // op1-proj discriminator (NeutrinoGraphAudit.cxx op1-proj pass, doc
    // pr/83 r4) restated for shower MEMBERSHIP: a charge-starved shower
    // member whose fit points are 2D-shadowed in >= 2 of 3 wire views by a
    // healthy partner segment ANYWHERE in the graph (395148: the 23.4 cm
    // ghost member's shadow partner is the 65.7 cm track 2D-coincident
    // with it -- NOT a fellow member; measured 1.00/1.00/0.25).  Same
    // (apa,face) required.  Placed after shower dedup/detach and BEFORE
    // the pi0 finders so pairing never sees the ghost.  The overlap/charge
    // recipe is duplicated from op1-proj rather than shared (that
    // production pass stays byte-untouched).  C++ default false =>
    // no pass => byte-identical.
    if (m_shower_ghost_member_drop && m_mvga_dup_tol > 0) {
        auto grouping = main_cluster ? main_cluster->grouping() : nullptr;
        // Per-view 2D overlap of `fa` in `fb` (fraction of A's valid fit
        // points within m_mvga_dup_tol of B's in view coords
        // (x, cos(a)z - sin(a)y)); returns the 3 fractions sorted
        // descending, or empty on any inconsistency.  op1-proj recipe.
        auto view_overlaps = [&](const std::vector<Fit>& fa,
                                 const std::vector<Fit>& fb,
                                 const std::pair<int,int>& paf) {
            std::vector<double> out;
            if (!grouping) return out;
            const auto [au, av, aw] = grouping->wire_angles(paf.first, paf.second);
            for (double ang : {au, av, aw}) {
                const double c = std::cos(ang), s = std::sin(ang);
                int nin = 0, ntot = 0;
                for (const auto& pa : fa) {
                    if (!pa.valid()) continue;
                    ++ntot;
                    const double xa = pa.point.x();
                    const double wa = c*pa.point.z() - s*pa.point.y();
                    double best = 1e30;
                    for (const auto& pb : fb) {
                        if (!pb.valid()) continue;
                        const double dx = pb.point.x() - xa;
                        const double dw = c*pb.point.z() - s*pb.point.y() - wa;
                        best = std::min(best, dx*dx + dw*dw);
                    }
                    if (best < m_mvga_dup_tol*m_mvga_dup_tol) ++nin;
                }
                if (ntot == 0) return std::vector<double>{};
                out.push_back(static_cast<double>(nin)/ntot);
            }
            std::sort(out.rbegin(), out.rend());
            return out;
        };
        // Median per-point dQ/dx (full segment, NOT the op1-proj stem form
        // -- these members sit far from any vertex) + fraction of
        // non-positive charge points.  dQ <= 0 points stay IN the median --
        // a projective ghost's defining signature is negative fitted charge.
        struct GhostStats { double med_dqdx{0}; double frac_nonpos{0}; int n{0}; };
        auto member_stats = [](const SegmentPtr& sg) {
            GhostStats st;
            std::vector<double> dqdx;
            int nneg = 0;
            for (const auto& f : sg->fits()) {
                if (!f.valid() || f.dx <= 0) continue;
                dqdx.push_back(f.dQ / f.dx);
                if (f.dQ <= 0) ++nneg;
            }
            st.n = static_cast<int>(dqdx.size());
            if (st.n == 0) return st;
            std::sort(dqdx.begin(), dqdx.end());
            st.med_dqdx = dqdx[dqdx.size()/2];
            st.frac_nonpos = static_cast<double>(nneg) / st.n;
            return st;
        };

        int n_ghost_dropped = 0;
        for (auto& shower : showers) {                       // IndexedShowerSet order
            if (std::abs(shower->get_particle_type()) == 13) continue;  // long-muon pseudo-shower
            bool fired = true;
            for (int round = 0; round < 4 && fired; ++round) {  // >1 ghost per shower is rare; cap
                fired = false;
                std::vector<SegmentPtr> members;
                for (auto edesc : ordered_edges(*shower, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (sg && sg->descriptor_valid()) members.push_back(sg);
                }
                if (members.size() < 2) break;
                // Partner pool: every segment in the graph, stable index
                // order (the shadow partner is usually NOT a member).
                std::vector<SegmentPtr> partners;
                for (const auto& vd : ordered_nodes(graph)) {
                    for (auto edesc : sorted_out_edges(vd, graph)) {
                        SegmentPtr sg = graph[edesc].segment;
                        if (!sg || !sg->descriptor_valid()) continue;
                        if (std::find(partners.begin(), partners.end(), sg) != partners.end()) continue;
                        partners.push_back(sg);
                    }
                }
                std::sort(partners.begin(), partners.end(), SegmentIndexCmp{});
                for (size_t i = 0; i < members.size() && !fired; ++i) {
                    SegmentPtr cand = members[i];
                    if (cand == shower->start_segment()) continue;
                    if (segments_in_long_muon.count(cand)) continue;
                    if (segment_track_length(cand) < m_shower_ghost_min_len) continue;
                    const auto cst = member_stats(cand);
                    if (cst.n < 2) continue;
                    const double cratio = cst.med_dqdx / m_mip_dqdx_median;
                    const bool starved = (cratio <= m_shower_ghost_dqdx_ratio)
                                         || (cst.frac_nonpos > 0.5);
                    if (!starved) continue;
                    // Same-(apa,face) fit points only (op1-proj rule).
                    std::pair<int,int> paf{-1,-1};
                    bool paf_ok = true;
                    for (const auto& f : cand->fits()) {
                        if (!f.valid() || f.paf.first < 0) continue;
                        if (paf.first < 0) paf = f.paf;
                        else if (f.paf != paf) { paf_ok = false; break; }
                    }
                    if (!paf_ok || paf.first < 0) continue;
                    for (size_t j = 0; j < partners.size() && !fired; ++j) {
                        SegmentPtr part = partners[j];
                        if (part == cand) continue;
                        const auto pst = member_stats(part);
                        if (pst.n < 2) continue;
                        const double pratio = pst.med_dqdx / m_mip_dqdx_median;
                        if (pratio < 2*m_shower_ghost_dqdx_ratio) continue;  // partner must be healthy
                        if (pst.frac_nonpos > 0.5) continue;
                        bool part_paf_ok = true;
                        for (const auto& f : part->fits()) {
                            if (f.valid() && f.paf.first >= 0 && f.paf != paf) { part_paf_ok = false; break; }
                        }
                        if (!part_paf_ok) continue;
                        auto ov = view_overlaps(cand->fits(), part->fits(), paf);
                        if (ov.size() < 3) continue;
                        if (ov[1] > 0.25) {
                            SPDLOG_LOGGER_TRACE(s_log,
                                "pr99 ghost_member eval shower_id={} cand seg={} len={:.2f}cm "
                                "ratio={:.2f} fneg={:.2f} vs seg={} views={:.2f}/{:.2f}/{:.2f}",
                                shower->get_shower_id(), cand->id(),
                                segment_track_length(cand)/units::cm, cratio, cst.frac_nonpos,
                                part->id(), ov[0], ov[1], ov[2]);
                        }
                        if (ov[1] < m_shower_ghost_overlap_frac) continue;

                        // Verdict: ghost.  View removal first (leaf-only
                        // contract; refusal leaves everything untouched)...
                        auto [gv1, gv2] = find_vertices(graph, cand);
                        if (!shower->drop_ghost_member(cand)) {
                            SPDLOG_LOGGER_DEBUG(s_log,
                                "pr99 ghost_member refuse shower_id={} seg={} (leaf-only guard)",
                                shower->get_shower_id(), cand->id());
                            continue;
                        }
                        // ...then the graph deletion so the ghost leaves the
                        // PF/Bee display; degree-0 endpoints leave with it
                        // (never the main vertex, never a protected break --
                        // the mvga cleanup_vertex rule).
                        SPDLOG_LOGGER_DEBUG(s_log,
                            "pr99 ghost_member drop shower_id={} seg={} len={:.2f}cm "
                            "med_ratio={:.2f} fneg={:.2f} views={:.2f}/{:.2f}/{:.2f} "
                            "partner seg={} ratio={:.2f}",
                            shower->get_shower_id(), cand->id(),
                            segment_track_length(cand)/units::cm, cratio, cst.frac_nonpos,
                            ov[0], ov[1], ov[2], part->id(), pratio);
                        // doc pr/102 round 2 crash fix (SBND 18255-399998,
                        // rc=135 with other_seg_uncover_3d on): the graph
                        // removals below invalidate cand's/gv's descriptors,
                        // but drop_ghost_member cleaned only THIS shower's
                        // view -- any OTHER shower whose filter sets still
                        // hold those descriptors is left iterating freed
                        // storage (Shower::fill_sets UAF via ordered_nodes).
                        // Purge them from every view first, while the
                        // descriptors are still valid; TrajectoryView
                        // removal is a pure filter-set erase (no-op when the
                        // view never held the item), so a shower untouched
                        // by the ghost is byte-identically untouched here.
                        for (auto& other : showers) {
                            if (!other || other == shower) continue;
                            other->TrajectoryView::remove_segment(cand);
                        }
                        remove_segment(graph, cand);
                        for (VertexPtr gv : {gv1, gv2}) {
                            if (!gv || gv == main_vertex) continue;
                            if (gv->flags_any(VertexFlags::kProtectedBreak)) continue;
                            if (!gv->descriptor_valid()) continue;
                            if (boost::degree(gv->get_descriptor(), graph) == 0) {
                                for (auto& other : showers) {
                                    if (!other) continue;
                                    other->TrajectoryView::remove_vertex(gv);
                                }
                                remove_vertex(graph, gv);
                            }
                        }
                        // Mirror the detach refresh: vote + kinematics +
                        // charge over the remaining members.
                        shower->update_particle_type(particle_data, recomb_model, m_mip_dqdx,
                                                     main_vertex, m_shower_proton_daughter_pion,
                                                     m_mip_dqdx_median, m_shower_vote_track_pid_counts,
                                                     m_shower_accept_pid_guard,
                                                     m_shower_pid_guard_min_len);
                        shower->calculate_kinematics(particle_data, recomb_model,
                                                     m_shower_endpoint_exclude_start_vertex,
                                                     m_shower_endpoint_skip_orphan_vtx);
                        shower->set_kine_charge(cal_kine_charge(shower, m_charge_2d_u, m_charge_2d_v,
                                                                m_charge_2d_w, m_map_apa_ch_plane_wires,
                                                                track_fitter, dv));
                        shower->set_flag_kinematics(true);
                        ++n_ghost_dropped;
                        fired = true;
                    }
                }
            }
        }
        if (n_ghost_dropped) {
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                               map_vertex_to_shower, used_shower_clusters);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr99 ghost_member: {} member(s) dropped; maps rebuilt", n_ghost_dropped);
        }
    }

    // doc sbnd_xin/docs/pr/123 round 1 (shower_pass4_prune_detached): the
    // owner's over-reach line (decision 2026-08-28: detached-from-the-body
    // OR track-like-beyond-it), applied to the FINAL shower geometry.  The
    // per-absorb census measured that no absorb-time threshold separates --
    // legitimate EM growth chains far fragments exactly like over-reach does
    // (doc pr/123 sec 5) -- so the test runs here, on the finished object:
    // single-linkage components of the membership at m_shower_pass4_prune_gap
    // (owner: 40 cm); a component not containing the start segment is
    // detached and leaves the shower.  Disposition (owner decision): the
    // component RE-SEEDS as its own shower rooted at its member nearest the
    // kept body (conn 3 near / 4 far, the in_other_clusters typing), so no
    // segment is orphaned and the pi0 pairing sees the separated object.
    // Placed with the dedup/detach/ghost family: after every pass that
    // grows, merges, detaches or drops members, before the pr/119 census,
    // the hadronic tag, the final kine recompute and the pi0 finders.
    // C++ default false => no pass => byte-identical.
    if (m_shower_pass4_prune_detached && m_shower_pass4_prune_gap > 0) {
        std::vector<ShowerPtr> prune_order(showers.begin(), showers.end());
        std::sort(prune_order.begin(), prune_order.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
            auto* sa = a->start_segment().get();
            auto* sb = b->start_segment().get();
            if (!sa || !sb) return sa < sb;
            int cid_a = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
            int cid_b = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
            if (cid_a != cid_b) return cid_a < cid_b;
            return sa->id() < sb->id();
        });
        int n_comp_pruned = 0;
        std::vector<ShowerPtr> reseeded;
        for (auto& shower : prune_order) {
            SegmentPtr ss = shower->start_segment();
            if (!ss || !ss->descriptor_valid()) continue;
            if (std::abs(shower->get_particle_type()) == 13) continue;  // long-muon pseudo-shower
            // Member list in stable view order, unique.
            std::vector<SegmentPtr> mem;
            {
                std::unordered_set<size_t> seen;
                for (auto edesc : ordered_edges(*shower, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (!sg || !sg->descriptor_valid()) continue;
                    if (seen.insert(graph[sg->get_descriptor()].index).second) mem.push_back(sg);
                }
            }
            if (mem.size() < 2) continue;
            const size_t nmem = mem.size();
            size_t start_idx = nmem;
            for (size_t i = 0; i < nmem; ++i) if (mem[i] == ss) { start_idx = i; break; }
            if (start_idx == nmem) continue;
            // Pairwise min fit-cloud distances; a member with no usable fit
            // cloud is immune (forced into the kept component).
            auto pair_dis_fn = [&](const SegmentPtr& a, const SegmentPtr& b) -> double {
                double best = -1.0;
                for (const auto& f : a->fits()) {
                    if (!f.valid()) continue;
                    const double d = segment_get_closest_point(b, f.point).first;
                    if (d >= 0 && (best < 0 || d < best)) best = d;
                }
                return best;
            };
            std::vector<std::vector<double>> dmat(nmem, std::vector<double>(nmem, -1.0));
            for (size_t i = 0; i < nmem; ++i)
                for (size_t k = i + 1; k < nmem; ++k) {
                    double d = pair_dis_fn(mem[i], mem[k]);
                    if (d < 0) d = pair_dis_fn(mem[k], mem[i]);
                    dmat[i][k] = dmat[k][i] = d;
                }
            // Union-find at the gap; cloudless pairs (d < 0) never link.
            std::vector<size_t> parent(nmem);
            for (size_t i = 0; i < nmem; ++i) parent[i] = i;
            std::function<size_t(size_t)> find = [&](size_t x) {
                while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                return x;
            };
            for (size_t i = 0; i < nmem; ++i)
                for (size_t k = i + 1; k < nmem; ++k)
                    if (dmat[i][k] >= 0 && dmat[i][k] < m_shower_pass4_prune_gap) {
                        size_t pi_ = find(i), pk = find(k);
                        if (pi_ != pk) parent[pi_] = pk;
                    }
            // A member whose distance to EVERY other member is unknown is
            // immune: union it with the start component.
            for (size_t i = 0; i < nmem; ++i) {
                bool any = false;
                for (size_t k = 0; k < nmem; ++k) if (k != i && dmat[i][k] >= 0) { any = true; break; }
                if (!any) { size_t pi_ = find(i), pk = find(start_idx); if (pi_ != pk) parent[pi_] = pk; }
            }
            const size_t keep_root = find(start_idx);
            std::map<size_t, std::vector<size_t>> comps;   // root -> member indices (index order = stable)
            for (size_t i = 0; i < nmem; ++i) if (find(i) != keep_root) comps[find(i)].push_back(i);
            if (comps.empty()) continue;
            // Deterministic component order: by smallest member index.
            std::vector<std::vector<size_t>> comp_list;
            for (auto& [root_, idxs] : comps) comp_list.push_back(idxs);
            std::sort(comp_list.begin(), comp_list.end(),
                      [](const std::vector<size_t>& a, const std::vector<size_t>& b) { return a.front() < b.front(); });
            const WireCell::Point sv_pt = shower->start_vertex()
                ? (shower->start_vertex()->fit().valid() ? shower->start_vertex()->fit().point
                                                         : shower->start_vertex()->wcpt().point)
                : WireCell::Point(0, 0, 0);
            for (const auto& comp : comp_list) {
                std::vector<SegmentPtr> comp_segs;
                for (size_t i : comp) comp_segs.push_back(mem[i]);
                const int removed = shower->detach_member_set(comp_segs);
                if (!removed) {
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr123 pass4_prune: refuse shower_id={} comp_nseg={}",
                        shower->get_shower_id(), comp_segs.size());
                    continue;
                }
                // Root of the re-seed: component member nearest the KEPT
                // body; min-index tie-break (dmat is symmetric and stable).
                size_t root_i = comp.front();
                double root_d = -1.0;
                for (size_t i : comp) {
                    double d_keep = -1.0;
                    for (size_t k = 0; k < nmem; ++k) {
                        if (find(k) != keep_root) continue;
                        if (dmat[i][k] >= 0 && (d_keep < 0 || dmat[i][k] < d_keep)) d_keep = dmat[i][k];
                    }
                    if (d_keep >= 0 && (root_d < 0 || d_keep < root_d)) { root_d = d_keep; root_i = i; }
                }
                SegmentPtr root_sg = mem[root_i];
                const double root_sv_dis = segment_get_closest_point(root_sg, sv_pt).first;
                const int conn = (root_sv_dis >= 0 && root_sv_dis <= 80 * units::cm) ? 3 : 4;
                ShowerPtr ns = std::make_shared<Shower>(graph);
                if (shower->start_vertex()) ns->set_start_vertex(shower->start_vertex(), conn);
                ns->set_start_segment(root_sg);
                for (size_t i : comp) if (mem[i] != root_sg) ns->add_segment(mem[i], true);
                reseeded.push_back(ns);
                ++n_comp_pruned;
                if (pr93_absorb_dbg()) {
                    std::fprintf(stderr,
                        "SHOWER_ABSORB PASS4_PRUNE shower=%d comp_root=%d comp_nseg=%zu "
                        "gap_cm=%.1f root_body_dis_cm=%.2f conn=%d\n",
                        pr91_seg_display_id(ss), pr91_seg_display_id(root_sg), comp.size(),
                        m_shower_pass4_prune_gap / units::cm,
                        root_d >= 0 ? root_d / units::cm : -1.0, conn);
                }
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr123 pass4_prune: shower_id={} sheds {} detached seg(s), reseed root seg={} conn={}",
                    shower->get_shower_id(), removed, root_sg->id(), conn);
            }
        }
        if (n_comp_pruned) {
            for (auto& ns : reseeded) showers.insert(ns);
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                               map_vertex_to_shower, used_shower_clusters);
            calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                        graph, track_fitter, dv, particle_data, recomb_model);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr123 pass4_prune: {} detached component(s) re-seeded; {} shower(s) now",
                n_comp_pruned, showers.size());
        }
    }

    // doc sbnd_xin/docs/pr/124 front A (shower_pass4_prune_gap2): the 25-40
    // cm gap band.  The pr/123 prune's G=40 was deliberately conservative;
    // the offline component scan (pr124_gapband_scan.py, both label sets,
    // seg-level join) measured the band at 15 BAD / 2 COL components and a
    // ZERO-collateral qualifier pair on the FINAL body re-split at gap2:
    // a sub-component not holding the start segment is over-reach when its
    // charge centroid sits > prune2_ang degrees off the kept core (seen from
    // the start vertex) OR its pooled median dQ/dx exceeds prune2_mdqdx MIP
    // (heavily-ionizing hadronic blob, e.g. 406125's q_frac-0.949 proton at
    // 3.6x median-MIP).  Thresholds are in m_mip_dqdx_median multiples; the
    // scan's plateau-normalized 2.0 converts to 2.5 here (54657.7/43000).
    // Machinery forked from the tier-1 pass above (house fork-by-duplication;
    // tier 1 stays byte-identical).  Disposition identical: detach + re-seed
    // (conn 3/4).  Operating point measured with tier 1 ON at 40.  C++
    // default 0 => no pass => byte-identical.
    if (m_shower_pass4_prune_gap2 > 0) {
        std::vector<ShowerPtr> prune_order(showers.begin(), showers.end());
        std::sort(prune_order.begin(), prune_order.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
            auto* sa = a->start_segment().get();
            auto* sb = b->start_segment().get();
            if (!sa || !sb) return sa < sb;
            int cid_a = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
            int cid_b = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
            if (cid_a != cid_b) return cid_a < cid_b;
            return sa->id() < sb->id();
        });
        int n_comp_pruned = 0;
        std::vector<ShowerPtr> reseeded;
        for (auto& shower : prune_order) {
            SegmentPtr ss = shower->start_segment();
            if (!ss || !ss->descriptor_valid()) continue;
            if (std::abs(shower->get_particle_type()) == 13) continue;  // long-muon pseudo-shower
            std::vector<SegmentPtr> mem;
            {
                std::unordered_set<size_t> seen;
                for (auto edesc : ordered_edges(*shower, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (!sg || !sg->descriptor_valid()) continue;
                    if (seen.insert(graph[sg->get_descriptor()].index).second) mem.push_back(sg);
                }
            }
            if (mem.size() < 2) continue;
            const size_t nmem = mem.size();
            size_t start_idx = nmem;
            for (size_t i = 0; i < nmem; ++i) if (mem[i] == ss) { start_idx = i; break; }
            if (start_idx == nmem) continue;
            auto pair_dis_fn = [&](const SegmentPtr& a, const SegmentPtr& b) -> double {
                double best = -1.0;
                for (const auto& f : a->fits()) {
                    if (!f.valid()) continue;
                    const double d = segment_get_closest_point(b, f.point).first;
                    if (d >= 0 && (best < 0 || d < best)) best = d;
                }
                return best;
            };
            std::vector<std::vector<double>> dmat(nmem, std::vector<double>(nmem, -1.0));
            for (size_t i = 0; i < nmem; ++i)
                for (size_t k = i + 1; k < nmem; ++k) {
                    double d = pair_dis_fn(mem[i], mem[k]);
                    if (d < 0) d = pair_dis_fn(mem[k], mem[i]);
                    dmat[i][k] = dmat[k][i] = d;
                }
            std::vector<size_t> parent(nmem);
            for (size_t i = 0; i < nmem; ++i) parent[i] = i;
            std::function<size_t(size_t)> find = [&](size_t x) {
                while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                return x;
            };
            for (size_t i = 0; i < nmem; ++i)
                for (size_t k = i + 1; k < nmem; ++k)
                    if (dmat[i][k] >= 0 && dmat[i][k] < m_shower_pass4_prune_gap2) {
                        size_t pi_ = find(i), pk = find(k);
                        if (pi_ != pk) parent[pi_] = pk;
                    }
            for (size_t i = 0; i < nmem; ++i) {
                bool any = false;
                for (size_t k = 0; k < nmem; ++k) if (k != i && dmat[i][k] >= 0) { any = true; break; }
                if (!any) { size_t pi_ = find(i), pk = find(start_idx); if (pi_ != pk) parent[pi_] = pk; }
            }
            const size_t keep_root = find(start_idx);
            std::map<size_t, std::vector<size_t>> comps;
            for (size_t i = 0; i < nmem; ++i) if (find(i) != keep_root) comps[find(i)].push_back(i);
            if (comps.empty()) continue;
            std::vector<std::vector<size_t>> comp_list;
            for (auto& [root_, idxs] : comps) comp_list.push_back(idxs);
            std::sort(comp_list.begin(), comp_list.end(),
                      [](const std::vector<size_t>& a, const std::vector<size_t>& b) { return a.front() < b.front(); });
            const WireCell::Point sv_pt = shower->start_vertex()
                ? (shower->start_vertex()->fit().valid() ? shower->start_vertex()->fit().point
                                                         : shower->start_vertex()->wcpt().point)
                : WireCell::Point(0, 0, 0);
            // Angle apex: the start SEGMENT's own vertex (what the dump
            // writes as start_vertex_id and what the offline qualifier scan
            // measured) -- NOT the shower start vertex, which for conn-3/4
            // re-seeded showers is the far main vertex and compresses the
            // angle (168432: 59 deg local vs 24 deg from main).
            WireCell::Point apex_pt = sv_pt;
            {
                auto [p3sv, p3ev] = find_vertices(graph, ss);
                if (p3sv) apex_pt = p3sv->fit().valid() ? p3sv->fit().point : p3sv->wcpt().point;
            }
            // Charge centroid of a member-index set from valid fit points.
            auto centroid_fn = [&](const std::vector<size_t>& idxs, bool from_keep) -> std::pair<WireCell::Point, double> {
                WireCell::Point acc(0, 0, 0);
                double qsum = 0;
                for (size_t i = 0; i < nmem; ++i) {
                    const bool in_set = from_keep ? (find(i) == keep_root)
                                                  : (std::find(idxs.begin(), idxs.end(), i) != idxs.end());
                    if (!in_set) continue;
                    for (const auto& f : mem[i]->fits()) {
                        if (!f.valid()) continue;
                        const double q = std::abs(f.dQ);
                        acc += f.point * q;
                        qsum += q;
                    }
                }
                if (qsum > 0) acc = acc / qsum;
                return {acc, qsum};
            };
            const auto [kcen, kq] = centroid_fn({}, true);
            for (const auto& comp : comp_list) {
                const auto [ccen, cq] = centroid_fn(comp, false);
                double ang_deg = -1.0;
                if (kq > 0 && cq > 0) {
                    const WireCell::Vector vk = kcen - apex_pt;
                    const WireCell::Vector vc = ccen - apex_pt;
                    if (vk.magnitude() > 0.1 * units::cm && vc.magnitude() > 0.1 * units::cm) {
                        double c = vk.dot(vc) / (vk.magnitude() * vc.magnitude());
                        c = std::max(-1.0, std::min(1.0, c));
                        ang_deg = std::acos(c) / M_PI * 180.0;
                    }
                }
                double med_ratio = -1.0;
                if (m_mip_dqdx_median > 0) {
                    std::vector<double> dqdx;
                    for (size_t i : comp)
                        for (const auto& f : mem[i]->fits()) {
                            if (!f.valid() || f.dx <= 0) continue;
                            dqdx.push_back(std::abs(f.dQ) / f.dx);
                        }
                    if (!dqdx.empty()) {
                        std::sort(dqdx.begin(), dqdx.end());
                        med_ratio = dqdx[dqdx.size() / 2] / m_mip_dqdx_median;
                    }
                }
                const bool over = (ang_deg >= 0 && ang_deg > m_shower_pass4_prune2_ang) ||
                                  (med_ratio >= 0 && med_ratio > m_shower_pass4_prune2_mdqdx);
                if (pr93_absorb_dbg()) {
                    std::fprintf(stderr,
                        "SHOWER_ABSORB PASS4_PRUNE2 shower=%d comp_first=%d comp_nseg=%zu "
                        "ang_deg=%.1f med_mip=%.2f prune=%d\n",
                        pr91_seg_display_id(ss), pr91_seg_display_id(mem[comp.front()]),
                        comp.size(), ang_deg, med_ratio, (int)over);
                }
                if (!over) continue;
                std::vector<SegmentPtr> comp_segs;
                for (size_t i : comp) comp_segs.push_back(mem[i]);
                const int removed = shower->detach_member_set(comp_segs);
                if (!removed) {
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "pr124 pass4_prune2: refuse shower_id={} comp_nseg={}",
                        shower->get_shower_id(), comp_segs.size());
                    continue;
                }
                size_t root_i = comp.front();
                double root_d = -1.0;
                for (size_t i : comp) {
                    double d_keep = -1.0;
                    for (size_t k = 0; k < nmem; ++k) {
                        if (find(k) != keep_root) continue;
                        if (dmat[i][k] >= 0 && (d_keep < 0 || dmat[i][k] < d_keep)) d_keep = dmat[i][k];
                    }
                    if (d_keep >= 0 && (root_d < 0 || d_keep < root_d)) { root_d = d_keep; root_i = i; }
                }
                SegmentPtr root_sg = mem[root_i];
                const double root_sv_dis = segment_get_closest_point(root_sg, sv_pt).first;
                const int conn = (root_sv_dis >= 0 && root_sv_dis <= 80 * units::cm) ? 3 : 4;
                ShowerPtr ns = std::make_shared<Shower>(graph);
                if (shower->start_vertex()) ns->set_start_vertex(shower->start_vertex(), conn);
                ns->set_start_segment(root_sg);
                for (size_t i : comp) if (mem[i] != root_sg) ns->add_segment(mem[i], true);
                reseeded.push_back(ns);
                ++n_comp_pruned;
                SPDLOG_LOGGER_DEBUG(s_log,
                    "pr124 pass4_prune2: shower_id={} sheds {} band seg(s) (ang={:.1f} med_mip={:.2f}), reseed root seg={} conn={}",
                    shower->get_shower_id(), removed, ang_deg, med_ratio, root_sg->id(), conn);
            }
        }
        if (n_comp_pruned) {
            for (auto& ns : reseeded) showers.insert(ns);
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                               map_vertex_to_shower, used_shower_clusters);
            calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                        graph, track_fitter, dv, particle_data, recomb_model);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr124 pass4_prune2: {} band component(s) re-seeded; {} shower(s) now",
                n_comp_pruned, showers.size());
        }
    }

    // doc sbnd_xin/docs/pr/119 -- byte-neutral membership census at the
    // exact seat where a Phase-B expel pass would run: after every pass that
    // grows, merges, detaches or drops members, before the hadronic tag, the
    // final kine recompute and the pi0 finders.  Env-gated stderr only.
    pr119_probe_expel_groups(showers, graph, main_vertex, track_fitter, dv,
                             m_mip_dqdx_median);

    // doc pr/99 round 3 (A5, shower_hadronic_tag).  Design block at the
    // m_shower_hadronic_* members (NeutrinoPatternBase.h).  For every
    // claimed-EM conn-1 shower, walk the trajectory (fit-cloud points, arc
    // proxy s = |p - start_vtx|, the recipe of the validated python scan
    // sbnd_xin/scripts/analysis/pr99/pr99_transition_scan.py) and measure
    // (i) in-cylinder 3D imaged-point population growth over the first
    // scan_len (union of Cluster::kd_radius hits over {main, other}
    // clusters, per-cluster point-index dedup) and (ii) the terminal Bragg
    // rise of the member fits' median dQ/dx.  Hadronic verdict stamps the
    // START segment pdg 211 (the pi0 incoming-track stamp recipe) plus the
    // shower's cached type, and registers the shower id for the pi0-finder
    // guards.  Placed BEFORE the kine recompute and the pi0 finders.  A
    // DEBUG census line is emitted for every evaluated shower.  C++ default
    // false => no pass => byte-identical.
    if (m_shower_hadronic_tag && m_shower_hadronic_bin > 0) {
        std::vector<Facade::Cluster*> cluster_pool;
        if (main_cluster) cluster_pool.push_back(main_cluster);
        for (auto* oc : other_clusters) if (oc) cluster_pool.push_back(oc);

        // Event-level cloud of EVERY graph segment's fit points, with a
        // parallel row -> segment-graph-index map (append order preserved
        // by add_points), for the ownership test below.
        std::shared_ptr<Facade::DynamicPointCloud> all_fit;
        std::vector<size_t> row2segidx;
        for (auto e : ordered_edges(graph)) {
            SegmentPtr sg = graph[e].segment;
            if (!sg || !sg->descriptor_valid()) continue;
            auto dpc = sg->dpcloud("fit");
            if (!dpc || !dpc->npoints()) continue;
            if (!all_fit) {
                all_fit = std::make_shared<Facade::DynamicPointCloud>(dpc->get_wpid_params());
            } else {
                all_fit->merge_wpid_params(*dpc);
            }
            all_fit->add_points(*dpc);
            row2segidx.resize(all_fit->npoints(), graph[sg->get_descriptor()].index);
        }

        int n_retyped = 0;
        for (auto& shower : showers) {                       // IndexedShowerSet order
            if (!shower) continue;
            if (std::abs(shower->get_particle_type()) != 11) continue;
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            if (conn_type != 1 || !start_vtx) continue;
            SegmentPtr ss = shower->start_segment();
            if (!ss || !ss->descriptor_valid()) continue;
            if (segments_in_long_muon.count(ss)) continue;
            if (!ss->has_particle_info() || !ss->particle_info()) continue;

            auto fit_pc = shower->get_pcloud("fit");
            if (!fit_pc || fit_pc->npoints() < 2) continue;

            const WireCell::Point vp = start_vtx->fit().valid()
                ? start_vtx->fit().point : start_vtx->wcpt().point;

            const size_t npts = fit_pc->npoints();
            std::vector<double> s_arc(npts);
            double smax = 0;
            for (size_t i = 0; i < npts; ++i) {
                const WireCell::Point p = fit_pc->point3d(i);
                s_arc[i] = (p - vp).magnitude();
                if (s_arc[i] > smax) smax = s_arc[i];
            }
            if (smax < m_shower_hadronic_min_len) continue;  // too short to judge

            // Member set for the ownership test: this shower's segment
            // graph indices.
            std::set<size_t> member_idx;
            for (auto edesc : ordered_edges(*shower, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (sg && sg->descriptor_valid()) member_idx.insert(graph[sg->get_descriptor()].index);
            }

            // (i) growth: per-bin union population within r_cyl of the
            // bin's axis points, over the scan window only -- OWNERSHIP-
            // filtered: an imaged point counts only if its nearest fit
            // point (over ALL graph segments, all_fit cloud above) belongs
            // to a member segment.  Calibration (doc pr/99 r3): without the
            // filter, vertex-region activity from OTHER prongs inflates a
            // real electron's early bins and fakes a shrinking profile
            // (46363: 2365->335 raw but member-owned growth 3.1); with it,
            // all 36 nue-selected primaries read >= 2.3 while the misID'd
            // hadrons read <= 0.7.
            const double scan_len = std::min(smax, m_shower_hadronic_scan_len);
            const int nbins = std::max(1, static_cast<int>(std::ceil(scan_len / m_shower_hadronic_bin)));
            std::vector<double> bin_n(nbins, 0.0);
            for (size_t ci = 0; ci < cluster_pool.size(); ++ci) {
                const Facade::Cluster* cl = cluster_pool[ci];
                std::vector<std::set<size_t>> uniq(nbins);
                for (size_t i = 0; i < npts; ++i) {
                    if (s_arc[i] >= scan_len) continue;
                    const int b = std::min(nbins - 1, static_cast<int>(s_arc[i] / m_shower_hadronic_bin));
                    const WireCell::Point p = fit_pc->point3d(i);
                    for (const auto& [pidx, dist2] : cl->kd_radius(m_shower_hadronic_r_cyl, p)) {
                        uniq[b].insert(pidx);
                    }
                }
                const auto& cpts = cl->points();
                for (int b = 0; b < nbins; ++b) {
                    for (size_t pidx : uniq[b]) {
                        if (all_fit && !row2segidx.empty()) {
                            const WireCell::Point q(cpts[0][pidx], cpts[1][pidx], cpts[2][pidx]);
                            auto rr = all_fit->kd3d().knn(1, q);
                            if (rr.empty()) continue;
                            const size_t row = rr[0].first;
                            if (row >= row2segidx.size() || !member_idx.count(row2segidx[row])) continue;
                        }
                        bin_n[b] += 1.0;
                    }
                }
            }
            double n_early = 0, n_late = 0;
            for (int b = 0; b < nbins; ++b) {
                if ((b + 0.5) * m_shower_hadronic_bin < 0.5 * scan_len) n_early += bin_n[b];
                else n_late += bin_n[b];
            }
            // ends-ratio (the round-1/calibration definition): mean of the
            // last two populated bins over mean of the first two.
            double g_first = 0, g_last = 0;
            {
                int nf = 0, nl = 0;
                for (int b = 0; b < std::min(2, nbins); ++b)
                    if (bin_n[b] > 0) { g_first += bin_n[b]; ++nf; }
                for (int b = std::max(0, nbins - 2); b < nbins; ++b)
                    if (bin_n[b] > 0) { g_last += bin_n[b]; ++nl; }
                if (nf) g_first /= nf;
                if (nl) g_last /= nl;
            }
            const bool has_growth = (g_first > 0 && g_last > 0);
            const double growth = has_growth ? g_last / g_first : -1.0;

            // (ii) terminal Bragg: per-bin median member-fit dQ/dx over the
            // FULL trajectory; trunk = median of bins excluding the last
            // two, term = max of the last two.
            const int nbins_full = std::max(1, static_cast<int>(std::ceil(smax / m_shower_hadronic_bin)));
            std::vector<std::vector<double>> bin_dqdx(nbins_full);
            for (auto edesc : ordered_edges(*shower, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg || !sg->descriptor_valid()) continue;
                for (const auto& f : sg->fits()) {
                    if (!f.valid() || f.dx <= 0) continue;
                    const double sv = (f.point - vp).magnitude();
                    if (sv >= smax) continue;
                    const int b = std::min(nbins_full - 1, static_cast<int>(sv / m_shower_hadronic_bin));
                    bin_dqdx[b].push_back(f.dQ / f.dx);
                }
            }
            auto bin_med = [&](int b) -> double {
                auto& v = bin_dqdx[b];
                if (v.empty()) return -1.0;
                std::sort(v.begin(), v.end());
                return v[v.size() / 2];
            };
            std::vector<double> trunk_meds;
            double dqdx_term = -1.0;
            for (int b = 0; b < nbins_full; ++b) {
                const double m = bin_med(b);
                if (m < 0) continue;
                if (b < nbins_full - 2) trunk_meds.push_back(m);
                else if (m > dqdx_term) dqdx_term = m;
            }
            double dqdx_trunk = -1.0;
            if (!trunk_meds.empty()) {
                std::sort(trunk_meds.begin(), trunk_meds.end());
                dqdx_trunk = trunk_meds[trunk_meds.size() / 2];
            }
            const double bragg = (dqdx_trunk > 0 && dqdx_term > 0) ? dqdx_term / dqdx_trunk : -1.0;

            // Proton-like stem branch (395148: stem 3.2 x MIP, growth 0.87
            // -- above the main cut but hadronic).  Pair-conversion gamma
            // stems read ~2 x MIP, so the production value 3.0 keeps 1.5x
            // margin on gammas; high-stem REAL primaries are protected by
            // the growth ceiling (their measured growth >= 2.78).
            double stem_med = -1.0;
            if (m_shower_hadronic_stem_ratio > 0) {
                auto stem = shower->get_stem_dQ_dx(start_vtx, ss, 20, m_mip_dqdx_median);
                if (stem.size() >= 2) {
                    std::vector<double> head(stem.begin(),
                                             stem.begin() + std::min<size_t>(6, stem.size()));
                    std::sort(head.begin(), head.end());
                    stem_med = head[head.size() / 2];
                }
            }

            const bool verdict = has_growth
                && (growth < m_shower_hadronic_growth_max
                    || (bragg >= m_shower_hadronic_bragg_ratio
                        && growth < m_shower_hadronic_growth_bragg)
                    || (m_shower_hadronic_stem_ratio > 0
                        && stem_med >= m_shower_hadronic_stem_ratio
                        && growth < m_shower_hadronic_growth_bragg));

            SPDLOG_LOGGER_DEBUG(s_log,
                "A5 hadronic census: shower id={} pdg={} conn={} nseg={} smax={:.1f}cm "
                "growth={:.2f} n_early={:.0f} n_late={:.0f} dqdx_trunk={:.0f} dqdx_term={:.0f} "
                "bragg={:.2f} stem={:.2f} verdict={}",
                shower->get_shower_id(), shower->get_particle_type(), conn_type,
                shower->get_num_segments(), smax / units::cm,
                growth, n_early, n_late, dqdx_trunk, dqdx_term, bragg, stem_med, verdict ? 1 : 0);

            if (!verdict) continue;

            // Stamp: the pi0 incoming-track recipe (this file, id_pi0 stamp).
            ss->particle_info()->set_pdg(211);
            ss->particle_info()->set_mass(139.57 * units::MeV);
            if (ss->particle_info()->kinetic_energy() > 0) {
                auto four_momentum = segment_cal_4mom(ss, 211, particle_data, recomb_model, m_mip_dqdx);
                ss->particle_info()->set_four_momentum(four_momentum);
            }
            shower->set_particle_type(211);
            m_hadronic_retyped_shower_ids.insert(shower->get_shower_id());
            ++n_retyped;
            SPDLOG_LOGGER_DEBUG(s_log,
                "A5 hadronic retype: shower id={} start_seg len={:.1f}cm -> pdg 211 "
                "(growth={:.2f} bragg={:.2f})",
                shower->get_shower_id(), segment_track_length(ss) / units::cm, growth, bragg);
        }
        if (n_retyped) {
            SPDLOG_LOGGER_DEBUG(s_log, "A5 hadronic tag: {} shower(s) re-typed 211", n_retyped);
            // doc pr/101 (K3): the re-typed showers are now hadronic, so
            // their best energy follows the hadronic rule.  No-op when off.
            for (auto& shower : showers) {
                if (shower && m_hadronic_retyped_shower_ids.count(shower->get_shower_id()))
                    apply_hadronic_dqdx_best(shower);
            }
        }
    }

    // doc sbnd_xin/docs/pr/125 (owner 2026-08-29, SBND 18259-37112): a
    // track-typed shower that shares a NON-main start vertex with a bigger
    // EM shower at essentially zero cloud gap is the same physical shower
    // with a mis-typed half (37112: the pdg-2212 half of a gamma conversion,
    // stem dQ/dx 2.4-9 MIP, 1.28 cm from the 549 MeV gamma at shared vertex
    // 84104).  merge_shower_fragments cannot take it (EM<->EM only, by
    // design -- pr/117 sec 7); this dedicated pass absorbs it.  Placement
    // (measured twice): AFTER the pass4 prunes (before them, the tier-2
    // union-find re-detached what these passes folded in; after them, this
    // pass absorbs both 37112's proton stem, gap 1.28 cm, and the
    // prune-re-seeded comp 50031, gap 3.81 cm -- the owner-adjudicated "one
    // shower, everything") and AFTER shower_hadronic_tag (before it, the
    // tagger re-evaluated the ENLARGED showers and flipped 22 of them >15
    // MeV to pdg 211, including 69314's own primary electron -- production
    // tagging must see unmodified showers), before
    // recompute_shower_kine_charge_final + the pi0 finders.  Manifest scan
    // (docs/pr/pr125-pairs.tsv): the shared-non-main-vertex + gap<6cm +
    // track-typed gate fires on 37112 alone across both manifests; next
    // candidates sit at 17.4+ cm.  Knob off => no pass => byte-identical.
    if (m_shower_samevtx_track_absorb) {
        // Vertex identity by the dump's display id (cluster_id*1000 +
        // graph_index, PrDisplayDump::vertex_display_id) -- pointer identity
        // is unreliable across re-slicing (same graph vertex, different
        // object), which is what the offline scan compared.
        auto vtx_display_id = [](const VertexPtr& v) -> int {
            if (!v) return -1;
            const auto* cl = v->cluster();
            return (cl ? cl->get_cluster_id() : 0) * 1000 + static_cast<int>(v->get_graph_index());
        };
        const int main_vid = vtx_display_id(main_vertex);
        std::vector<ShowerPtr> sv_order(showers.begin(), showers.end());
        std::sort(sv_order.begin(), sv_order.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
            auto* sa = a->start_segment().get();
            auto* sb = b->start_segment().get();
            if (!sa || !sb) return sa < sb;
            int cid_a = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
            int cid_b = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
            if (cid_a != cid_b) return cid_a < cid_b;
            return sa->id() < sb->id();
        });
        auto collect_members = [&](const ShowerPtr& shw) {
            std::vector<SegmentPtr> mem;
            std::unordered_set<size_t> seen;
            for (auto edesc : ordered_edges(*shw, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg || !sg->descriptor_valid()) continue;
                if (seen.insert(graph[sg->get_descriptor()].index).second) mem.push_back(sg);
            }
            return mem;
        };
        auto cross_gap = [&](const std::vector<SegmentPtr>& ma,
                             const std::vector<SegmentPtr>& mb) -> double {
            double best = -1.0;
            for (const auto& a : ma) {
                for (const auto& f : a->fits()) {
                    if (!f.valid()) continue;
                    for (const auto& b : mb) {
                        const double d = segment_get_closest_point(b, f.point).first;
                        if (d >= 0 && (best < 0 || d < best)) best = d;
                    }
                }
            }
            return best;
        };
        std::vector<std::pair<ShowerPtr, ShowerPtr>> sv_plan;   // (absorber, fragment)
        for (auto& frag : sv_order) {
            const int fpdg = std::abs(frag->get_particle_type());
            if (!(fpdg == 13 || fpdg == 211 || fpdg == 2212)) continue;
            const int fvid = vtx_display_id(frag->start_vertex());
            if (fvid < 0 || fvid == main_vid) continue;
            const double flen = frag->get_total_length();
            if (m_shower_samevtx_absorb_max_len > 0 && flen > m_shower_samevtx_absorb_max_len) continue;
            // Floor: sub-5cm crumbs are satellite-class, not this pass's
            // business (first build absorbed 8 of them and crumb-triggered
            // host PID re-votes -- the measured 37112 gate is a SUBSTANTIAL
            // mis-typed fragment, 33.4 cm).
            if (m_shower_samevtx_absorb_min_len > 0 && flen < m_shower_samevtx_absorb_min_len) continue;
            auto fmem = collect_members(frag);
            if (fmem.empty()) continue;
            ShowerPtr best;
            double best_gap = -1.0;
            for (auto& host : sv_order) {
                if (host == frag) continue;
                if (std::abs(host->get_particle_type()) != 11) continue;
                if (vtx_display_id(host->start_vertex()) != fvid) continue;
                if (!(host->get_total_length() > flen)) continue;
                const double g = cross_gap(fmem, collect_members(host));
                if (g >= 0 && g < m_shower_samevtx_absorb_gap && (best_gap < 0 || g < best_gap)) {
                    best = host;
                    best_gap = g;
                }
            }
            if (pr91_merge_dbg())
                std::fprintf(stderr, "SHOWER_MERGE tag=samevtx_absorb frag_sid=%d frag_node=%d "
                                     "pdg=%d len_cm=%.1f host_sid=%d gap_cm=%.2f verdict=%s\n",
                             frag->get_shower_id(), pr91_seg_display_id(frag->start_segment()),
                             frag->get_particle_type(), flen / units::cm,
                             best ? best->get_shower_id() : -1,
                             best_gap >= 0 ? best_gap / units::cm : -1.0,
                             best ? "MERGE" : "no_host");
            if (best) sv_plan.emplace_back(best, frag);
        }
        int n_sv_merged = 0;
        for (auto& [host, frag] : sv_plan) {
            pr93_probe_absorb_splice("samevtx_absorb", host, frag);
            host->add_shower(*frag);
            host->update_particle_type(particle_data, recomb_model, m_mip_dqdx, main_vertex, m_shower_proton_daughter_pion, m_mip_dqdx_median, m_shower_vote_track_pid_counts, m_shower_accept_pid_guard, m_shower_pid_guard_min_len);
            host->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
            host->set_kine_charge(cal_kine_charge(host, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv));
            host->set_flag_kinematics(true);
            showers.erase(frag);
            ++n_sv_merged;
        }
        if (n_sv_merged) {
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                               map_vertex_to_shower, used_shower_clusters);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr125 samevtx_absorb: {} track-typed fragment shower(s) absorbed; {} remain",
                n_sv_merged, showers.size());
        }
    }

    // doc sbnd_xin/docs/pr/125 (owner 2026-08-29, SBND 18255-69314): the PF
    // list renders a fragmented EM shower as a "cascade of electrons" -- in
    // 69314, 38 pdg-11 entries of which 28 are below 5 MeV.  The crumbs are
    // separate PR::Showers seeded from graph-CONNECTED pdg-11 segments: their
    // start vertex is a vertex OF the big shower's own member chain.  Absorb
    // such satellites (EM-typed, conn 2/3, kine below the cap) into the
    // connected host.  Manifest scan (docs/pr/pr125-satellites.tsv +
    // vertex-connected variant): 31 IN-marked recoveries, 0 OUT-marked once
    // conn-4 (the disconnected class -- its nominal vertex link is an
    // artifact) is excluded.  Knob off => no pass => byte-identical.
    if (m_shower_satellite_absorb) {
        std::vector<ShowerPtr> sat_order(showers.begin(), showers.end());
        std::sort(sat_order.begin(), sat_order.end(), [](const ShowerPtr& a, const ShowerPtr& b) {
            auto* sa = a->start_segment().get();
            auto* sb = b->start_segment().get();
            if (!sa || !sb) return sa < sb;
            int cid_a = sa->cluster() ? sa->cluster()->get_cluster_id() : -1;
            int cid_b = sb->cluster() ? sb->cluster()->get_cluster_id() : -1;
            if (cid_a != cid_b) return cid_a < cid_b;
            return sa->id() < sb->id();
        });
        // Vertex identity by dump display id -- see the samevtx pass note.
        auto sat_vtx_display_id = [](const VertexPtr& v) -> int {
            if (!v) return -1;
            const auto* cl = v->cluster();
            return (cl ? cl->get_cluster_id() : 0) * 1000 + static_cast<int>(v->get_graph_index());
        };
        // Hosts: EM showers above the host floor.  Vertex-ID sets
        // snapshotted BEFORE any absorb (single pass, no chaining).
        std::vector<std::pair<ShowerPtr, std::set<int>>> sat_hosts;
        for (auto& host : sat_order) {
            if (std::abs(host->get_particle_type()) != 11) continue;
            if (host->get_kine_charge() < m_shower_satellite_absorb_host_mev) continue;
            std::set<int> vset;
            std::unordered_set<size_t> seen;
            for (auto edesc : ordered_edges(*host, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg || !sg->descriptor_valid()) continue;
                if (!seen.insert(graph[sg->get_descriptor()].index).second) continue;
                auto [v0, v1] = find_vertices(graph, sg);
                if (v0) vset.insert(sat_vtx_display_id(v0));
                if (v1) vset.insert(sat_vtx_display_id(v1));
            }
            if (!vset.empty()) sat_hosts.emplace_back(host, std::move(vset));
        }
        std::vector<std::pair<ShowerPtr, ShowerPtr>> sat_plan;   // (host, satellite)
        std::set<Shower*> sat_claimed;                           // membership only
        for (auto& sat : sat_order) {
            if (std::abs(sat->get_particle_type()) != 11) continue;
            const int conn = sat->get_start_vertex_and_type().second;
            if (conn != 2 && conn != 3) continue;
            if (sat->get_kine_charge() >= m_shower_satellite_absorb_max_mev) continue;
            // Attach vertex: the shower's stored start vertex, matched by
            // dump display id -- this is exactly what the dump writes as
            // start_vertex_id and what the offline scan joined on.  (First
            // build compared VertexPtr identity: conn-3 sats hold a
            // different object for the same graph vertex, 5/7 missed in
            // 69314; second build used find_vertices(start_seg), but the
            // attach vertex is NOT an endpoint of the sat's own segment.)
            const int sv_id = sat_vtx_display_id(sat->start_vertex());
            if (sv_id < 0) continue;
            ShowerPtr host;
            for (auto& [h, vset] : sat_hosts) {
                if (h == sat || sat_claimed.count(h.get())) continue;
                if (!(h->get_kine_charge() > sat->get_kine_charge())) continue;
                if (vset.count(sv_id)) { host = h; break; }
            }
            if (pr91_merge_dbg())
                std::fprintf(stderr, "SHOWER_MERGE tag=satellite_absorb sat_sid=%d sat_node=%d "
                                     "ke_mev=%.2f conn=%d host_sid=%d verdict=%s\n",
                             sat->get_shower_id(), pr91_seg_display_id(sat->start_segment()),
                             sat->get_kine_charge() / units::MeV, conn,
                             host ? host->get_shower_id() : -1,
                             host ? "MERGE" : "no_host");
            if (host) {
                sat_plan.emplace_back(host, sat);
                sat_claimed.insert(sat.get());
            }
        }
        int n_sat_merged = 0;
        std::set<Shower*> sat_recalc;                            // membership only
        for (auto& [host, sat] : sat_plan) {
            pr93_probe_absorb_splice("satellite_absorb", host, sat);
            host->add_shower(*sat);
            showers.erase(sat);
            sat_recalc.insert(host.get());
            ++n_sat_merged;
        }
        if (n_sat_merged) {
            for (auto& [host, vset] : sat_hosts) {
                if (!sat_recalc.count(host.get())) continue;
                // Deliberately NO update_particle_type here: a sub-10-MeV
                // crumb must not re-vote its host's PID (first build re-typed
                // 37112's 84070 e->pi via the vote; kinematics-only is the
                // minimal absorb).
                host->calculate_kinematics(particle_data, recomb_model, m_shower_endpoint_exclude_start_vertex, m_shower_endpoint_skip_orphan_vtx);
                host->set_kine_charge(cal_kine_charge(host, m_charge_2d_u, m_charge_2d_v, m_charge_2d_w, m_map_apa_ch_plane_wires, track_fitter, dv));
                host->set_flag_kinematics(true);
            }
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                               map_vertex_to_shower, used_shower_clusters);
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr125 satellite_absorb: {} satellite shower(s) absorbed; {} remain",
                n_sat_merged, showers.size());
        }
    }

    // doc pr/99 round 3 (C1/C1b, kine_charge_dedup / kine_charge_rebuild).
    // Final knob-gated recompute of every shower's kine_charge with
    // cross-shower 2D-cell ownership and/or member-true clouds -- placed
    // after ALL shower-structure passes (so every mid-pipeline gate saw
    // legacy values) and BEFORE the pi0 finders (which cache
    // get_kine_charge() at entry).  No-op when both knobs are off.
    // doc pr/132 round 5 (K16): EM collinear-fragment merge -- BEFORE the
    // final kine recompute and the pi0 finders, so merged charge feeds both.
    // No-op when the knob is off (default), byte-identical.
    em_collinear_merge(showers, map_vertex_in_shower, map_segment_in_shower,
                       map_vertex_to_shower, used_shower_clusters,
                       track_fitter, dv, particle_data, recomb_model);

    // doc pr/132 round 6 (K17): EM start back-extension -- after the forward
    // merge, before the final kine recompute and the finders.  No-op off.
    em_start_backext(showers, map_vertex_in_shower, map_segment_in_shower,
                     map_vertex_to_shower, used_shower_clusters,
                     track_fitter, dv, particle_data, recomb_model);

    // doc sbnd_xin/docs/pr/138 Phase B: the SPLITTER -- the last structure pass
    // and the only one that cuts.  After every merging pass (the owner's
    // "merge them together, then separate cleanly") and BEFORE the pi0
    // finders, so a gamma pair that was over-clustered into one shower can
    // still be paired.  No-op when the knob is off (default), byte-identical.
    shower_split(graph, main_vertex, showers, map_vertex_in_shower, map_segment_in_shower,
                 map_vertex_to_shower, used_shower_clusters, vertices_in_long_muon,
                 segments_in_long_muon, track_fitter, dv, particle_data, recomb_model);

    recompute_shower_kine_charge_final(showers, graph, track_fitter, dv);

    // Identify pi0 with vertex.
    // doc pr/33 F3: both finders get a reference to the same local copy, so
    // nothing propagates past this function either way (the caller's
    // variable is separately seeded by reference into ssm_tagger -- §10.10
    // amendment 1).  Knob-off restores the copy between the two finders, so
    // each seeds from the same base = the legacy by-value behavior.
    int pi0_acc = acc_segment_id;
    // doc pr/132 round 9 (K19): back-projection NC vertex proposer -- runs
    // BEFORE path 1 so a corrected NC vertex is set before the gammas are
    // consumed (the round-8 ordering lock).  No-op when the knob is off.
    id_pi0_backproject_vertex(pi0_acc, pi0_showers, map_shower_pio_id, map_pio_id_showers,
                              map_pio_id_mass, graph, main_vertex, showers,
                              map_vertex_in_shower, map_segment_in_shower,
                              map_vertex_to_shower, used_shower_clusters,
                              segments_in_long_muon, track_fitter, dv, particle_data, recomb_model);
    id_pi0_with_vertex(pi0_acc, pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass,
                      map_pio_id_saved_pair, pio_kine, graph, main_vertex, showers, main_cluster,
                      other_clusters, map_cluster_main_vertices, map_vertex_in_shower,
                      map_segment_in_shower, map_vertex_to_shower, used_shower_clusters,
                      segments_in_long_muon, track_fitter, dv, particle_data, recomb_model);
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

    // doc pr/132 round 3: over-merge substructure tape (byte-neutral, stderr
    // only; docstring at pr132_pi0_substruct_dbg).
    if (pr132_pi0_substruct_dbg())
        pr132_pi0_substruct_probe(showers, pi0_showers);

    // doc pr/132 round 4: NC vertex-in-shower recognizer + start-derivation
    // audit (byte-neutral, stderr only; docstring at pr132_pi0_ncvtx_dbg).
    if (pr132_pi0_ncvtx_dbg())
        pr132_pi0_ncvtx_probe(showers, pi0_showers, main_vertex, map_vertex_in_shower);

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