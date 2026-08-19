#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Logging.h"
#include <chrono>
#include <cstdlib>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// doc sbnd_xin/docs/pr/31 §12 (F1, was P1 + P3's 4-momentum half + P4).
//
// The prototype's reclassification idiom mutates two members and GUARDS the
// third (NeutrinoID_track_shower.h:372-374 and identically at 10 more sites):
//     sg->set_particle_type(11);
//     sg->set_particle_mass(mp.get_mass_electron());
//     if (sg->get_particle_4mom(3)>0) sg->cal_4mom();
// particle_4mom[3] = kine_energy + particle_mass (ProtoSegment.cxx:1437-1444),
// so ">0" means "previously computed": a segment that never had an energy
// keeps zeros.  The toolkit has no independent members -- reclassification
// constructs a whole new Aux::ParticleInfo -- which forces a 4-momentum
// decision the prototype never makes.  Fifteen sites made it three ways
// (11x unconditional recompute, 1x topology-branch recompute, 3x rest-mass
// (m,0,0,0) overwrite) and only the two guarded sites (improve_maps_one_in,
// improve_maps_shower_in_track_out's shower re-calc, which comments the
// guard) made it the prototype's way.
//
// preserve=false reproduces today's recompute-and-construct byte-for-byte
// (the shape-C sites keep their legacy rest-mass expression on their own
// else-branch).  preserve=true is the prototype's single rule at every
// site: recompute only where the prototype's guard passes
// (proto_recomputes && previously-had-energy), otherwise carry the
// existing 4-momentum forward verbatim, or zeros if there never was one.
// energy() reads m_four_momentum.e(), the exact analog of the prototype's
// get_particle_4mom(3).
//
// Why the preserve path cannot use the validating constructor: the carried
// states are exactly the ones Aux::ParticleInfo::validate_inputs forbids --
// all-zero (E < m) for never-computed, and an old on-shell 4-momentum
// carried across a mass change (energy-momentum relation violated).  The
// prototype holds both states as a matter of course (type and mass move,
// particle_4mom does not), and validate_inputs itself carries a
// commented-out "Allow zero 4-momentum as a placeholder" block for the
// first.  Rather than relax the shared aux validator, the preserve path
// constructs a legal placeholder and then writes the carried value through
// set_four_momentum() -- the class's own non-validating setter -- so KE
// lands at E - m (i.e. -m for never-computed), which is what any
// prototype consumer subtracting mass sees.
// doc sbnd_xin/docs/pr/40 round 2 F6: the !had branch below used to finish with
// pinfo->set_four_momentum(D4Vector(0,0,0,0)) -- an explicit, deliberate
// choice (see the block comment above) to make energy()==0 read exactly
// like the prototype's never-computed particle_4mom[3]==0.  But
// ParticleInfo::kinetic_energy() is e() - mass, so that same zero 4-vector
// makes kinetic_energy() == -mass: a NEGATIVE kinetic energy, wherever a
// caller (e.g. the Bee PF-tree writer, MultiAlgBlobClustering.cxx
// fill_bee_pf_tree) reads kinetic_energy() rather than energy() directly.
// Found while chasing pr/40 round 2 F4's zero-KE display defect; same failure
// shape, sharper (a negative MeV is worse than a zero one on a display no
// prototype consumer or Bee viewer ever actually reads as a raw
// particle_4mom[3] subtraction).  reclass_never_computed_ke_floor=true
// leaves the just-constructed (mass,0,0,0) placeholder in place instead
// of overwriting it to all-zero, i.e. never-computed reads KE==0, matching
// the had==true branch's own "keep it sane" spirit.  false = legacy
// set_four_momentum(0,0,0,0) = byte-identical.  Independent of `preserve`
// (this only touches the !preserve-recompute, !had sub-case) and of F4/F5.
static std::shared_ptr<WireCell::Aux::ParticleInfo> reclass_pinfo(
    const SegmentPtr& sg, int pdg_code,
    const WireCell::Clus::ParticleDataSet::pointer& particle_data,
    const WireCell::IRecombinationModel::pointer& recomb_model,
    double mip_scale, bool preserve, bool proto_recomputes,
    bool never_computed_ke_floor = false)
{
    const double mass = particle_data->get_particle_mass(pdg_code);
    const bool had = sg->has_particle_info() && sg->particle_info()->energy() > 0;
    if (!preserve || (proto_recomputes && had)) {
        auto four_momentum = segment_cal_4mom(sg, pdg_code, particle_data, recomb_model, mip_scale);
        return std::make_shared<WireCell::Aux::ParticleInfo>(
            pdg_code, mass, particle_data->pdg_to_name(pdg_code), four_momentum);
    }
    auto pinfo = std::make_shared<WireCell::Aux::ParticleInfo>(
        pdg_code, mass, particle_data->pdg_to_name(pdg_code),
        WireCell::D4Vector<double>(mass, 0, 0, 0));
    if (had) {
        pinfo->set_four_momentum(sg->particle_info()->four_momentum());
    } else if (!never_computed_ke_floor) {
        pinfo->set_four_momentum(WireCell::D4Vector<double>(0, 0, 0, 0));  // legacy: KE = -mass
    }
    // else: leave the (mass,0,0,0) placeholder constructed above -- KE == 0.
    return pinfo;
}

void PatternAlgorithms::clustering_points(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv, const std::string& cloud_name, double search_range, double scaling_2d){
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    // Collect all segments that belong to this cluster
    std::vector<SegmentPtr> segments;
    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (seg && seg->cluster() == &cluster) {
            segments.push_back(seg);
        }
    }
    // if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "clustering_points timing: collect segments ({}) took {} ms", segments.size(), MS(Clock::now() - t0).count());

    // Run clustering on the collected segments
    t0 = Clock::now();
    if (!segments.empty()) {
        clustering_points_segments(segments, dv, cloud_name, search_range, scaling_2d, m_assoc_reassign_orphans);
    }
    // if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "clustering_points timing: clustering_points_segments took {} ms", MS(Clock::now() - t0).count());

    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "clustering_points timing: TOTAL took {} ms", MS(Clock::now() - t_total).count());
}

// doc sbnd_xin/docs/pr/59 round 2: inert unless m_assoc_full_recluster (owner
// flip pending).  See the declaration in NeutrinoPatternBase.h for the full
// rationale; summary: clustering_points runs once per cluster, but a segment
// can be created afterward (examine_structure_final*/examine_vertices_1
// inside determine_main_vertex, confirmed for 18255-142421 seg 20 and
// 116944-71372 segs 19052/19053/136199) and never gets a chance to compete for
// points -- and if the cluster's main-ness is later swapped away
// (swap_main_cluster), it never gets a second chance either.  This helper is
// meant to be called wherever such a segment could just have been created.
size_t PatternAlgorithms::reassociate_cluster_orphans(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv) {
    if (!m_assoc_full_recluster) return 0;
    static const bool pr59r2_census = std::getenv("WCT_PR59_ASSOC_CENSUS") != nullptr;

    // Collect all segments that belong to this cluster, in stable edge-index
    // order (never boost::edges/pointer order -- CLAUDE.md determinism rule).
    std::vector<SegmentPtr> segments;
    for (auto e : ordered_edges(graph)) {
        SegmentPtr seg = graph[e].segment;
        if (seg && seg->cluster() == &cluster) {
            segments.push_back(seg);
        }
    }
    if (segments.empty()) return 0;

    // Any orphan (null or empty associate_points cloud) in the cluster?  If
    // not, this is a byte-identical no-op -- do not touch anything.
    std::map<SegmentPtr, size_t, SegmentIndexCmp> npts_before;
    bool any_orphan = false;
    for (auto seg : segments) {
        auto dpc = seg->dpcloud("associate_points");
        size_t n = dpc ? dpc->npoints() : 0;
        npts_before[seg] = n;
        if (n == 0) any_orphan = true;
    }
    if (!any_orphan) return 0;

    // Owner constraint: when this fires, delete the OLD associate_points for
    // the WHOLE cluster and establish new ones via a fresh full-cluster
    // competition -- never rescue just the orphan in isolation, which would
    // let it win points by default with no already-good sibling able to
    // contest for them (clustering_points_segments is a Voronoi + 2D
    // ghost-removal competition among exactly the segments handed to it).
    for (auto seg : segments) {
        if (seg->dpcloud("associate_points")) {
            seg->dpcloud("associate_points", nullptr);
        }
    }
    // doc pr/64 round 7: thread m_assoc_reassign_orphans through the pr/59
    // recluster path too -- this call is LIVE in SBND production
    // (assoc_full_recluster=true, wct-pr-perevt.jsonnet), so leaving it on
    // the trailing default would give rescued segments different
    // association rules than the main clustering_points pass just above.
    clustering_points_segments(segments, dv, "associate_points", 1.2*units::cm, 0.7, m_assoc_reassign_orphans);

    // Owner constraint: this must run BEFORE track/shower separation so the
    // classification pass can actually consume the new cloud -- but only for
    // the segments that were orphaned; re-running is_shower_topology/
    // is_shower_trajectory on an already-correctly-classified sibling is pure
    // blast radius (these are per-segment, non-competing tests, unlike
    // association).  Same two calls, same arguments, as separate_track_shower's
    // loop body just below in this file.
    size_t n_rescued = 0;
    for (auto seg : segments) {
        if (npts_before[seg] != 0) continue;  // was not orphaned
        auto dpc = seg->dpcloud("associate_points");
        size_t n_after = dpc ? dpc->npoints() : 0;
        if (n_after == 0) continue;  // still orphaned (e.g. lost the ghost-removal contest again)
        ++n_rescued;
        segment_is_shower_topology(seg, false, m_mip_dqdx_median, m_shower_topo_demote_len, m_shower_topo_reset, m_shower_topo_dqdx_guard);
        if (!seg->flags_any(SegmentFlags::kShowerTopology)) {
            segment_is_shower_trajectory(seg, 10*units::cm, m_mip_dqdx, m_shower_traj_straight_guard);
        }
    }

    if (pr59r2_census) {
        for (auto seg : segments) {
            auto dpc = seg->dpcloud("associate_points");
            size_t n_after = dpc ? dpc->npoints() : 0;
            size_t n_before = npts_before[seg];
            if (n_before == n_after) continue;
            const char* tag = (n_before == 0 && n_after > 0) ? "rescued"
                             : (n_after == 0) ? "lost"
                             : "moved";
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr59r2 recluster: cluster {} segment {} npts {} -> {} [{}]",
                cluster.get_cluster_id(), seg->get_graph_index(), n_before, n_after, tag);
        }
    }
    return n_rescued;
}

void PatternAlgorithms::separate_track_shower(Graph&graph, Facade::Cluster& cluster) {
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    MS t_topology{0}, t_trajectory{0};

    // Iterate through all edges (segments) in the graph, in stable edge-index
    // order (boost::edges on a setS graph is pointer order, which varies run
    // to run).
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr seg = graph[ed].segment;

        // Skip if segment is null or doesn't belong to this cluster
        if (!seg || seg->cluster() != &cluster) continue;

        // First check if segment is a shower topology
        auto t0 = Clock::now();
        segment_is_shower_topology(seg, false, m_mip_dqdx_median, m_shower_topo_demote_len, m_shower_topo_reset, m_shower_topo_dqdx_guard);
        t_topology += MS(Clock::now() - t0);

        // If not shower topology, check if it's a shower trajectory
        if (!seg->flags_any(SegmentFlags::kShowerTopology)) {
            t0 = Clock::now();
            segment_is_shower_trajectory(seg, 10*units::cm, m_mip_dqdx, m_shower_traj_straight_guard);
            t_trajectory += MS(Clock::now() - t0);
        }
    }

    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "separate_track_shower timing: shower_topology={:.3f}ms shower_trajectory={:.3f}ms TOTAL={:.3f}ms",
        t_topology.count(), t_trajectory.count(), MS(Clock::now() - t_total).count());
}

void PatternAlgorithms::determine_direction(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model) {
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    MS t_shower_traj{0}, t_shower_topo{0}, t_track{0};

    // Iterate through all edges (segments) in the graph, in stable edge-index order.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr seg = graph[ed].segment;

        // Skip if segment is null or doesn't belong to this cluster
        if (!seg || seg->cluster() != &cluster) continue;

        // Get the two vertices of this segment
        auto [start_v, end_v] = find_vertices(graph, seg);
        if (!start_v || !end_v) {
            SPDLOG_LOGGER_TRACE(s_log, "determine_direction: Error in finding vertices for a segment");
            continue;
        }

        // Check if vertices match the segment endpoints (start_v should be at front, end_v at back)
        const auto& wcpts = seg->wcpts();
        if (wcpts.size() < 2) continue;

        auto front_pt = wcpts.front().point;
        auto back_pt = wcpts.back().point;

        // Determine which vertex is start and which is end based on point positions
        double dis_sv_front = ray_length(Ray{start_v->wcpt().point, front_pt});
        double dis_sv_back = ray_length(Ray{start_v->wcpt().point, back_pt});

        if (dis_sv_front > dis_sv_back) {
            std::swap(start_v, end_v);
        }

        // Count number of segments connected to each vertex
        int start_n = 0, end_n = 0;
        if (start_v->descriptor_valid()) {
            start_n = boost::degree(start_v->get_descriptor(), graph);
        }
        if (end_v->descriptor_valid()) {
            end_n = boost::degree(end_v->get_descriptor(), graph);
        }

        // doc sbnd_xin/docs/pr/40: WCT_PID_TRACE_DEBUG, env-gated (same idiom
        // as WCT_PID_WRITE_DEBUG/WCT_SHOWER_TOPO_DEBUG), un-gates the
        // segment_determine_dir_track / segment_determine_shower_direction_trajectory
        // "Seg ... Track/S_traj dirsign dir_weak pdg mass KE score" TRACE line
        // (PRSegmentFunctions.cxx flag_print block).  Was hardcoded false
        // with no config or env path at all before doc pr/40's attribution
        // work needed it.
        static const bool s_pr40_pid_trace_debug = std::getenv("WCT_PID_TRACE_DEBUG") != nullptr;
        bool flag_print = s_pr40_pid_trace_debug;
        // if (seg->cluster() == main_cluster) flag_print = true;

        auto t0 = Clock::now();
        if (seg->flags_any(SegmentFlags::kShowerTrajectory)) {
            // Trajectory shower
            segment_determine_shower_direction_trajectory(seg, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, flag_print, track_pid_options());
            t_shower_traj += MS(Clock::now() - t0);
        } else if (seg->flags_any(SegmentFlags::kShowerTopology)) {
            // Topology shower: determine direction, then set electron particle info
            //
            // doc sbnd_xin/docs/pr/31 §11 (F2, was P2).  The prototype runs
            // determine_dir_shower_topology here (ProtoSegment.cxx:1677-1710),
            // which sets particle_type/particle_mass and does NOT touch
            // flag_dir; its determine_shower_direction() is called from one
            // place in the whole tree and that place is stage 4
            // (NeutrinoID_track_shower.h:1532,
            // compare_main_vertices_all_showers).  This call mutates only
            // dirsign -- dirsign(0) on entry, dirsign(flag_dir) at the end --
            // and its return is discarded, so skipping it leaves the direction
            // segment_is_shower_topology set, which is the prototype's state.
            // Default false => the call runs => byte-identical.
            if (!m_shower_topo_proto_dir) {
                segment_determine_shower_direction(seg, particle_data, recomb_model, "associate_points", m_mip_dqdx_median, 0.4*units::cm, m_mip_dqdx, m_dir_track_median_local);
            }
            {
                const int pdg_code = 11; // electron
                // doc pr/31 §12 F1 shape B: the prototype's
                // determine_dir_shower_topology writes type and mass only --
                // NO 4-momentum recompute here (proto_recomputes=false).
                auto pinfo = reclass_pinfo(seg, pdg_code, particle_data, recomb_model, m_mip_dqdx_median, m_reclass_preserve_4mom, false, m_reclass_never_computed_ke_floor);
                seg->particle_info(pinfo);
                seg->particle_score(100.0);
            }
            t_shower_topo += MS(Clock::now() - t0);
        } else {
            // Track
            // doc pr/48: an arm of a two-end dQ/dx break travels AWAY from
            // the break vertex -- the break's own two-arm stopping-template
            // accept established it (each arm's Bragg is at its outer end).
            // A dirsign stamp cannot carry this here (separate_track_shower's
            // segment_is_shower_topology zeroes dirsign on every segment
            // under shower_topo_reset), so the direction is reconstructed
            // from the graph: the arm's endpoint carrying
            // VertexFlags::kProtectedBreak IS the junction.  Let that stand
            // over a WEAK recompute -- the weak KS decision on a MIP-ish arm
            // is a coin flip, and landing "into the junction" defeats the
            // break's vertex at main-vertex scoring (51513/56211/57485 vs
            // 57903's lucky flips).  A STRONG recompute (!dir_weak) wins.
            // Inert unless m_two_end_break (the flag is only ever set by the
            // knob-gated pass).
            const bool teb_arm = m_two_end_break && seg->flags_any(SegmentFlags::kTwoEndBreakArm);
            int teb_outward = 0;
            if (teb_arm) {
                auto [v_front, v_back] = find_vertices(graph, seg);
                const bool fp = v_front && v_front->flags_any(VertexFlags::kProtectedBreak);
                const bool bp = v_back && v_back->flags_any(VertexFlags::kProtectedBreak);
                if (fp && !bp) teb_outward = 1;        // junction at front: travel front->back
                else if (bp && !fp) teb_outward = -1;  // junction at back:  travel back->front
            }
            segment_determine_dir_track(seg, start_n, end_n, particle_data, recomb_model, m_mip_dqdx_median, flag_print, track_pid_options());
            if (teb_arm) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "determine_direction: teb arm outward={} post_dirsign={} weak={}",
                    teb_outward, seg->dirsign(), seg_dir_weak(seg) ? 1 : 0);
            }
            if (teb_outward != 0 && seg->dirsign() != teb_outward && seg_dir_weak(seg)) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "determine_direction: two-end-break arm weak recompute dirsign {} -> outward {}",
                    seg->dirsign(), teb_outward);
                seg->dirsign(teb_outward);
            }
            t_track += MS(Clock::now() - t0);
        }

        {
            const char* seg_type = seg->flags_any(SegmentFlags::kShowerTrajectory) ? "S_traj"
                                 : seg->flags_any(SegmentFlags::kShowerTopology)   ? "S_topo"
                                 : "Track";
            double length = segment_track_length(seg);
            int    pdg    = seg->has_particle_info() ? seg->particle_info()->pdg()  : 0;
            SPDLOG_LOGGER_TRACE(s_log,
                "determine_direction: {} nfits={} nwcpts={} len={:.2f}cm dirsign={} dir_weak={}"
                " start_n={} end_n={} pdg={}",
                seg_type, seg->fits().size(), seg->wcpts().size(), length / units::cm,
                seg->dirsign(), seg_dir_weak(seg) ? 1 : 0,
                start_n, end_n, pdg);
        }
    }

    if (m_perf) SPDLOG_LOGGER_TRACE(s_log, "determine_direction timing: shower_traj={:.3f}ms shower_topo={:.3f}ms track={:.3f}ms TOTAL={:.3f}ms",
        t_shower_traj.count(), t_shower_topo.count(), t_track.count(), MS(Clock::now() - t_total).count());

   
}

std::pair<int, double> PatternAlgorithms::calculate_num_daughter_showers(Graph& graph, VertexPtr vertex, SegmentPtr segment, bool flag_count_shower) {
    int number_showers = 0;
    double acc_length = 0;
    
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;
    segments_to_be_examined.push_back(std::make_pair(vertex, segment));
    used_vertices.insert(vertex);
    
    while(segments_to_be_examined.size() > 0) {
        std::vector<std::pair<VertexPtr, SegmentPtr>> temp_segments;
        for (auto it = segments_to_be_examined.begin(); it != segments_to_be_examined.end(); it++) {
            VertexPtr prev_vtx = it->first;
            SegmentPtr current_sg = it->second;
            
            if (used_segments.find(current_sg) != used_segments.end()) continue; // looked at it before
            
            // Check if segment is a shower: trajectory flag, topology flag, or electron by dQ/dx
            // (matches prototype's get_flag_shower() = flag_shower_trajectory || flag_shower_topology || get_flag_shower_dQdx())
            bool is_shower = current_sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                             current_sg->flags_any(SegmentFlags::kShowerTopology) ||
                             (current_sg->has_particle_info() && std::abs(current_sg->particle_info()->pdg()) == 11);
            
            if (is_shower || (!flag_count_shower)) {
                number_showers++;
                acc_length += segment_track_length(current_sg);
            }
            used_segments.insert(current_sg);
            
            VertexPtr curr_vertex = find_other_vertex(graph, current_sg, prev_vtx);
            if (used_vertices.find(curr_vertex) != used_vertices.end()) continue;
            
            // Get all segments connected to curr_vertex (stable edge-index
            // order: acc_length is summed in BFS-frontier order).
            if (curr_vertex && curr_vertex->descriptor_valid()) {
                auto vd = curr_vertex->get_descriptor();
                for (auto edesc : sorted_out_edges(vd, graph)) {
                    SegmentPtr seg = graph[edesc].segment;
                    if (seg) {
                        temp_segments.push_back(std::make_pair(curr_vertex, seg));
                    }
                }
            }
            used_vertices.insert(curr_vertex);
        }
        segments_to_be_examined = temp_segments;
    }
    
    return std::make_pair(number_showers, acc_length);
}

// calculate_num_daughter_tracks: count tracks (non-shower segments) reachable from vertex
// via segment sg, skipping sg itself (used to find activity at the far end of a muon).
// Prototype: NeutrinoID::calculate_num_daughter_tracks in NeutrinoID_track_shower.h:724.
// flag_count_shower: if true also count shower segments; length_cut: only count if > cut.
std::pair<int, double> PatternAlgorithms::calculate_num_daughter_tracks(
    Graph& graph, VertexPtr vtx, SegmentPtr sg,
    bool flag_count_shower, double length_cut)
{
    int    number_tracks = 0;
    double acc_length    = 0;

    std::set<VertexPtr>  used_vertices;
    std::set<SegmentPtr> used_segments;

    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;
    segments_to_be_examined.push_back(std::make_pair(vtx, sg));
    used_vertices.insert(vtx);

    while (!segments_to_be_examined.empty()) {
        std::vector<std::pair<VertexPtr, SegmentPtr>> temp_segments;
        for (auto& [prev_vtx, current_sg] : segments_to_be_examined) {
            if (!used_segments.insert(current_sg).second) continue; // already seen

            bool is_shower = current_sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                             current_sg->flags_any(SegmentFlags::kShowerTopology) ||
                             (current_sg->has_particle_info() &&
                              std::abs(current_sg->particle_info()->pdg()) == 11);

            if (!is_shower || flag_count_shower) {
                double length = segment_track_length(current_sg);
                if (length > length_cut) {
                    acc_length += length;
                    number_tracks++;
                }
            }

            VertexPtr curr_vertex = find_other_vertex(graph, current_sg, prev_vtx);
            if (!curr_vertex || used_vertices.count(curr_vertex)) continue;
            used_vertices.insert(curr_vertex);

            if (curr_vertex->descriptor_valid()) {
                auto vd = curr_vertex->get_descriptor();
                // Stable edge-index order: acc_length sums in frontier order.
                for (auto edesc : sorted_out_edges(vd, graph)) {
                    SegmentPtr next_sg = graph[edesc].segment;
                    if (next_sg) temp_segments.push_back({curr_vertex, next_sg});
                }
            }
        }
        segments_to_be_examined = std::move(temp_segments);
    }

    return {number_tracks, acc_length};
}

// find_cont_muon_segment_nue: from vertex vtx, find a segment adjacent to sg that continues
// in roughly the same direction (opening angle < 12.5 deg) and has MIP-like dQ/dx.
// Returns {nullptr, nullptr} if no such continuation exists.
// Prototype: NeutrinoID::find_cont_muon_segment_nue in NeutrinoID_track_shower.h:2372.
std::pair<SegmentPtr, VertexPtr> PatternAlgorithms::find_cont_muon_segment_nue(
    Graph& graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx)
{
    SegmentPtr sg1  = nullptr;
    VertexPtr  vtx1 = nullptr;

    double sg_length  = segment_track_length(sg);
    WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;

    WireCell::Vector dir1 = segment_cal_dir_3vector(sg, vtx_pt, 15 * units::cm);
    // doc sbnd_xin/docs/pr/31 §12 (F5, was P6).  The prototype computes dir3
    // INSIDE the per-neighbour loop, always at 30 cm, whenever either length
    // qualifies (NeutrinoID_track_shower.h:2402-2408).  The toolkit's hoist is
    // correct -- dir3 does not depend on the loop variable -- but its fallback
    // to dir1 silently turns the reachable case "short reference segment, long
    // neighbour" (sg_length <= 30 cm < length) into a 15cm-vs-30cm comparison
    // where the prototype compares two 30 cm directions; angle1 feeds the
    // <12.5 deg muon-continuation test.  ON = unconditional 30 cm (still
    // hoisted).  OFF = today's conditional = byte-identical.
    WireCell::Vector dir3 = (m_cont_muon_dir3_30cm || sg_length > 30 * units::cm)
                                ? segment_cal_dir_3vector(sg, vtx_pt, 30 * units::cm)
                                : dir1;

    double max_length = 0;

    if (!vtx->descriptor_valid()) return {nullptr, nullptr};
    auto vd = vtx->get_descriptor();
    // Stable edge-index order: strict '>' argmax below keeps the first
    // candidate on a tie, so iteration order must not be pointer order.
    for (auto edesc : sorted_out_edges(vd, graph)) {
        SegmentPtr sg2 = graph[edesc].segment;
        if (!sg2 || sg2 == sg) continue;
        VertexPtr vtx2 = find_other_vertex(graph, sg2, vtx);

        double length = segment_track_length(sg2);
        double ratio  = segment_median_dQ_dx(sg2) / m_mip_dqdx_median;

        WireCell::Vector dir2 = segment_cal_dir_3vector(sg2, vtx_pt, 15 * units::cm);
        double angle = (M_PI - dir1.angle(dir2)) / M_PI * 180.0;

        double angle1 = angle;
        if (length > 30 * units::cm || sg_length > 30 * units::cm) {
            WireCell::Vector dir4 = segment_cal_dir_3vector(sg2, vtx_pt, 30 * units::cm);
            angle1 = (M_PI - dir3.angle(dir4)) / M_PI * 180.0;
        }

        bool angle_ok = (angle < 12.5 || angle1 < 12.5 ||
                         (sg_length < 6 * units::cm && (angle < 15 || angle1 < 15)));
        bool dqdx_ok  = (ratio < 1.3 || flag_ignore_dQ_dx);

        if (angle_ok && dqdx_ok) {
            double proj = length * std::cos(angle / 180.0 * M_PI);
            if (proj > max_length) {
                max_length = proj;
                sg1  = sg2;
                vtx1 = vtx2;
            }
        }
    }

    return {sg1, vtx1};
}

void PatternAlgorithms::examine_good_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data) {
    // Iterate through all edges (segments) in the graph, in stable edge-index order.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;

        // Skip if segment is null or doesn't belong to this cluster
        if (!sg || sg->cluster() != &cluster) continue;

        // Skip if segment is a shower (trajectory, topology, or electron by dQ/dx)
        // matches prototype get_flag_shower() = flag_shower_trajectory || flag_shower_topology || (particle_type==11)
        if (sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology) ||
            (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11)) continue;

        // Skip if no direction or weak direction
        if (sg->dirsign() == 0 || seg_dir_weak(sg)) continue;

        // Get the two vertices of this segment
        auto [vertex1, vertex2] = find_vertices(graph, sg);
        if (!vertex1 || !vertex2) continue;
        
        // Determine start and end vertices based on segment direction
        VertexPtr start_vertex = nullptr, end_vertex = nullptr;
        const auto& wcpts = sg->wcpts();
        if (wcpts.size() < 2) continue;
        
        auto front_pt = wcpts.front().point;
        auto back_pt = wcpts.back().point;
        
        if (sg->dirsign() == 1) {
            // Direction is forward (from front to back)
            double dis1_front = ray_length(Ray{vertex1->wcpt().point, front_pt});
            double dis1_back = ray_length(Ray{vertex1->wcpt().point, back_pt});
            if (dis1_front < dis1_back) {
                start_vertex = vertex1;
                end_vertex = vertex2;
            } else {
                start_vertex = vertex2;
                end_vertex = vertex1;
            }
        } else if (sg->dirsign() == -1) {
            // Direction is backward (from back to front)
            double dis1_front = ray_length(Ray{vertex1->wcpt().point, front_pt});
            double dis1_back = ray_length(Ray{vertex1->wcpt().point, back_pt});
            if (dis1_front < dis1_back) {
                start_vertex = vertex2;
                end_vertex = vertex1;
            } else {
                start_vertex = vertex1;
                end_vertex = vertex2;
            }
        }
        
        if (!start_vertex || !end_vertex) continue;
        
        // Calculate number of daughter showers
        auto result_pair = calculate_num_daughter_showers(graph, start_vertex, sg);
        int num_daughter_showers = result_pair.first;
        double length_daughter_showers = result_pair.second;
        
        // Calculate maximum angle between this segment and others at end_vertex
        double max_angle = 0;
        WireCell::Point end_pt = end_vertex->fit().valid() ? end_vertex->fit().point : end_vertex->wcpt().point;
        WireCell::Vector dir1 = segment_cal_dir_3vector(sg, end_pt, 15*units::cm);
        WireCell::Vector drift_dir(1, 0, 0);
        double min_para_angle = 1e9;
        
        // Get all segments connected to end_vertex
        if (end_vertex->descriptor_valid()) {
            auto vd = end_vertex->get_descriptor();
            const auto edge_range = sorted_out_edges(vd, graph);
            for (auto eit2 : edge_range) {
                SegmentPtr sg1 = graph[eit2].segment;
                if (!sg1 || sg1 == sg) continue;
                
                WireCell::Vector dir2 = segment_cal_dir_3vector(sg1, end_pt, 15*units::cm);
                double angle = std::acos(std::min(1.0, std::max(-1.0, dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude())))) / 3.1415926 * 180.0;
                if (angle > max_angle) max_angle = angle;
                
                angle = std::fabs(std::acos(std::min(1.0, std::max(-1.0, drift_dir.dot(dir2) / (drift_dir.magnitude() * dir2.magnitude())))) / 3.1415926 * 180.0 - 90.0);
                if (angle < min_para_angle) min_para_angle = angle;
            }
        }
        
        // Check if this track should be reclassified as an electron shower
        double drift_angle = std::fabs(std::acos(std::min(1.0, std::max(-1.0, drift_dir.dot(dir1) / (drift_dir.magnitude() * dir1.magnitude())))) / 3.1415926 * 180.0 - 90.0);
        double length = segment_track_length(sg);
        
        if ((num_daughter_showers >= 4 || (length_daughter_showers > 50*units::cm && num_daughter_showers >= 2)) &&
            (max_angle > 155 || (drift_angle < 15 && min_para_angle < 15 && min_para_angle + drift_angle < 25)) &&
            length < 15*units::cm) {
            
            // Reclassify as electron (PDG 11)
            // doc pr/31 §12 F1 shape C: the prototype writes type and mass only
            // (NeutrinoID_track_shower.h:257-261) -- the 4-momentum is untouched,
            // not zeroed to rest mass.
            if (m_reclass_preserve_4mom) {
                sg->particle_info(reclass_pinfo(sg, 11, particle_data, m_recomb_model, m_mip_dqdx, true, false, m_reclass_never_computed_ke_floor));
            }
            else {
                double em_mass = particle_data->get_particle_mass(11);
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    11,                                              // electron PDG
                    em_mass,                                         // electron mass
                    particle_data->pdg_to_name(11),                 // "e-"
                    WireCell::D4Vector<double>(em_mass, 0, 0, 0)    // at-rest 4-momentum
                );
                sg->particle_info(pinfo);
            }
            
            // Reset direction and mark as weak
            sg->dirsign(0);
            sg->dir_weak(true);
        }
        
        // Debug output (commented out)
        // std::cout << sg->get_id() << " " << sg->particle_type() << " " << num_daughter_showers << " " 
        //           << length/units::cm << " " << max_angle << " " << min_para_angle << " " << drift_angle << std::endl;
    }
}

void PatternAlgorithms::fix_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster){
    // Iterate through all vertices in the graph, in stable node-index order
    // (boost::vertices on a setS graph is pointer order).
    for (const auto& vd_it : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd_it].vertex;
        
        // Skip if vertex is null or doesn't belong to this cluster
        if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;
        
        // Check how many segments are connected to this vertex
        if (!vtx->descriptor_valid()) continue;
        auto vd = vtx->get_descriptor();
        if (boost::degree(vd, graph) <= 1) continue;
        
        int n_in = 0;
        int n_in_shower = 0;
        std::vector<SegmentPtr> in_tracks;
        
        // Get vertex position
        WireCell::Point vtx_point = vtx->wcpt().point;
        
        // Iterate through all segments connected to this vertex
        const auto edge_range = sorted_out_edges(vd, graph);
        for (auto eit : edge_range) {
            SegmentPtr sg = graph[eit].segment;
            if (!sg) continue;
            
            // Determine if this vertex is at the front or back of the segment
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            auto front_pt = wcpts.front().point;
            auto back_pt = wcpts.back().point;
            
            double dis_front = ray_length(Ray{vtx_point, front_pt});
            double dis_back = ray_length(Ray{vtx_point, back_pt});
            
            bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
            
            // Check if this segment is pointing "in" to the vertex
            // "in" means: (at front and direction is -1) OR (at back and direction is 1)
            if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                n_in++;

                // Check if it's a shower (trajectory, topology, or electron by dQ/dx)
                // matches prototype get_flag_shower() = flag_shower_trajectory || flag_shower_topology || (particle_type==11)
                if (sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology) ||
                    (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11)) {
                    n_in_shower++;
                } else {
                    in_tracks.push_back(sg);
                }
            }
        }
        
        // If there are multiple incoming tracks (not all showers), reset their directions
        if (n_in > 1 && n_in != n_in_shower) {
            for (auto it1 = in_tracks.begin(); it1 != in_tracks.end(); it1++) {
                (*it1)->dirsign(0);
                (*it1)->dir_weak(true);
            }
        }
    }
}

void PatternAlgorithms::fix_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster){
    // Iterate through all vertices in the graph.  Stable node-index order:
    // flipping a shower's dirsign at one vertex changes what its other
    // vertex sees, so vertex processing order affects the outcome.
    for (const auto& nd : ordered_nodes(graph)) {
        VertexPtr vtx = graph[nd].vertex;
        
        // Skip if vertex is null or doesn't belong to this cluster
        if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;
        
        // Check how many segments are connected to this vertex
        if (!vtx->descriptor_valid()) continue;
        auto vd = vtx->get_descriptor();
        if (boost::degree(vd, graph) <= 1) continue;
        
        std::vector<SegmentPtr> in_showers;
        bool flag_turn_shower_dir = false;
        
        // Get vertex position
        WireCell::Point vtx_point = vtx->wcpt().point;
        
        // Iterate through all segments connected to this vertex
        const auto edge_range = sorted_out_edges(vd, graph);
        for (auto eit : edge_range) {
            SegmentPtr sg = graph[eit].segment;
            if (!sg) continue;
            
            // Determine if this vertex is at the front or back of the segment
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            auto front_pt = wcpts.front().point;
            auto back_pt = wcpts.back().point;
            
            double dis_front = ray_length(Ray{vtx_point, front_pt});
            double dis_back = ray_length(Ray{vtx_point, back_pt});
            
            bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
            
            // Check if segment is a shower (matches prototype get_flag_shower())
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                             sg->flags_any(SegmentFlags::kShowerTopology) ||
                             (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
            
            // Check if this is an "incoming" segment (pointing into vertex)
            if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                if (is_shower) {
                    in_showers.push_back(sg);
                }
            }
            // Check if this is an "outgoing" segment (pointing away from vertex)
            else if ((flag_start && sg->dirsign() == 1) || (!flag_start && sg->dirsign() == -1)) {
                // If it's an outgoing non-shower track with strong direction
                if (!is_shower && !seg_dir_weak(sg)) {
                    flag_turn_shower_dir = true;
                }
            }
        }
        
        // If there's a strong outgoing track and incoming showers, flip shower directions
        if (flag_turn_shower_dir) {
            for (auto it1 = in_showers.begin(); it1 != in_showers.end(); it1++) {
                (*it1)->dirsign((*it1)->dirsign() * (-1));
                (*it1)->dir_weak(true);
            }
        }
    }
}

void PatternAlgorithms::improve_maps_one_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check){
    bool flag_update = true;
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    while(flag_update) {
        flag_update = false;
        
        // Iterate through all vertices in the graph
        for (auto vit : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vit].vertex;

            // Skip if vertex is null or doesn't belong to this cluster
            if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;

            // Check how many segments are connected to this vertex
            if (!vtx->descriptor_valid()) continue;
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) <= 1) continue;

            // Skip if already processed
            if (used_vertices.find(vtx) != used_vertices.end()) continue;

            int n_in = 0;
            std::vector<std::pair<SegmentPtr, bool>> map_sg_dir; // segment -> flag_start
            
            // Get vertex position
            WireCell::Point vtx_point = vtx->wcpt().point;
            
            // Iterate through all segments connected to this vertex
            const auto edge_range = sorted_out_edges(vd, graph);
            for (auto eit : edge_range) {
                SegmentPtr sg = graph[eit].segment;
                if (!sg) continue;
                
                // Skip if segment already processed
                if (used_segments.find(sg) != used_segments.end()) continue;
                
                // Determine if this vertex is at the front or back of the segment
                const auto& wcpts = sg->wcpts();
                if (wcpts.size() < 2) continue;
                
                auto front_pt = wcpts.front().point;
                auto back_pt = wcpts.back().point;
                
                double dis_front = ray_length(Ray{vtx_point, front_pt});
                double dis_back = ray_length(Ray{vtx_point, back_pt});
                
                bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
                
                // Check if this is an "incoming" segment (pointing into vertex)
                if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                    if (flag_strong_check) {
                        // Only count if direction is strong
                        if (!seg_dir_weak(sg)) n_in++;
                    } else {
                        n_in++;
                    }
                }
                
                // Collect segments with no or weak direction
                if (sg->dirsign() == 0 || seg_dir_weak(sg)) {
                    map_sg_dir.push_back({sg, flag_start});
                }
            }
            
            // If no segments to change direction, mark vertex as used
            if (map_sg_dir.size() == 0) {
                used_vertices.insert(vtx);
            }
            
            // If there are incoming segments, set all weak/no-direction segments to point out
            if (n_in > 0) {
                for (auto& [sg, flag_start] : map_sg_dir) {
                    
                    // Set direction to point away from vertex
                    if (flag_start) {
                        sg->dirsign(1);  // at front, point forward
                    } else {
                        sg->dirsign(-1); // at back, point backward
                    }
                    
                    // Recalculate 4-momentum if particle info exists
                    if (sg->has_particle_info()) {
                        int pdg_code = sg->particle_info()->pdg();
                        auto four_momentum = segment_cal_4mom(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx);
                        
                        // Update particle info with new 4-momentum
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(
                            pdg_code,
                            particle_data->get_particle_mass(pdg_code),
                            particle_data->pdg_to_name(pdg_code),
                            four_momentum
                        );
                        sg->particle_info(pinfo);
                    }
                    
                    sg->dir_weak(true);
                    used_segments.insert(sg);
                    flag_update = true;
                }
                used_vertices.insert(vtx);
            }
        }
    }
}

void PatternAlgorithms::improve_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check){
    bool flag_update = true;
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    while(flag_update) {
        flag_update = false;
        
        // Iterate through all vertices in the graph
        for (auto vit : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vit].vertex;

            // Skip if vertex is null or doesn't belong to this cluster
            if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;

            // Check how many segments are connected to this vertex
            if (!vtx->descriptor_valid()) continue;
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) <= 1) continue;

            // Skip if already processed
            if (used_vertices.find(vtx) != used_vertices.end()) continue;

            // int n_in = 0;
            int n_in_shower = 0;
            std::vector<SegmentPtr> out_tracks;
            std::vector<SegmentPtr> map_no_dir_segments; // segments with no direction
            
            // Get vertex position
            WireCell::Point vtx_point = vtx->wcpt().point;
            
            // Iterate through all segments connected to this vertex
            const auto edge_range = sorted_out_edges(vd, graph);
            for (auto eit : edge_range) {
                SegmentPtr sg = graph[eit].segment;
                if (!sg) continue;
                
                // Determine if this vertex is at the front or back of the segment
                const auto& wcpts = sg->wcpts();
                if (wcpts.size() < 2) continue;
                
                auto front_pt = wcpts.front().point;
                auto back_pt = wcpts.back().point;
                
                double dis_front = ray_length(Ray{vtx_point, front_pt});
                double dis_back = ray_length(Ray{vtx_point, back_pt});
                
                bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
                
                // matches prototype get_flag_shower() = flag_shower_trajectory || flag_shower_topology || (particle_type==11)
                bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                 sg->flags_any(SegmentFlags::kShowerTopology) ||
                                 (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);

                // Check if this is an "incoming" segment (pointing into vertex)
                if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                    if (is_shower) {
                        n_in_shower++;
                    }
                }
                // Check if this is an "outgoing" segment (pointing away from vertex)
                else if ((flag_start && sg->dirsign() == 1) || (!flag_start && sg->dirsign() == -1)) {
                    if (!is_shower) {
                        // Check if it's weak or has no particle type
                        bool no_particle_type = !sg->has_particle_info() || sg->particle_info()->pdg() == 0;
                        if (seg_dir_weak(sg) || (no_particle_type && !flag_strong_check)) {
                            out_tracks.push_back(sg);
                        }
                    }
                }
                // Segment with no direction
                else if (sg->dirsign() == 0) {
                    map_no_dir_segments.push_back(sg);
                }
            }
            
            // If there are incoming showers and outgoing tracks or no-direction segments
            if (n_in_shower > 0 && (out_tracks.size() > 0 || map_no_dir_segments.size() > 0)) {
                // Reclassify outgoing tracks as electrons
                for (auto it1 = out_tracks.begin(); it1 != out_tracks.end(); it1++) {
                    SegmentPtr sg1 = *it1;

                    // doc sbnd_xin/docs/pr/40 F2: spare a segment whose OWN
                    // median dQ/dx is decisively proton- or muon-like from
                    // this wholesale conversion.  false = legacy = every
                    // out_track (weak-direction or untyped) becomes electron
                    // unconditionally.
                    if (m_shower_reclass_dqdx_guard && segment_dqdx_spares_electron_reclass(sg1, m_mip_dqdx)) {
                        continue;
                    }

                    // Set as electron (PDG 11)
                    int pdg_code = 11;
                    auto pinfo = reclass_pinfo(sg1, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                    sg1->particle_info(pinfo);
                    sg1->dirsign(0);

                    flag_update = true;
                }
                
                // Process no-direction segments
                for (auto sg1 : map_no_dir_segments) {
                    if (used_segments.find(sg1) != used_segments.end()) continue;

                    // If it's not already a shower, reclassify as electron
                    // matches prototype: set particle_type=11 only if !get_flag_shower()
                    bool is_shower1 = sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                      sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                      (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11);
                    // doc sbnd_xin/docs/pr/40 F2 (same guard as the out_tracks
                    // loop above).
                    if (!is_shower1 && !(m_shower_reclass_dqdx_guard && segment_dqdx_spares_electron_reclass(sg1, m_mip_dqdx))) {
                        int pdg_code = 11;
                        auto pinfo = reclass_pinfo(sg1, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                        sg1->particle_info(pinfo);
                    }

                    // Prototype calls cal_4mom() for ALL segments here (including showers) if energy>0.
                    // Recalculate 4-momentum for showers that already have particle info with valid energy.
                    if (is_shower1 && sg1->has_particle_info() && sg1->particle_info()->energy() > 0) {
                        int pdg_code = sg1->particle_info()->pdg();
                        auto four_momentum = segment_cal_4mom(sg1, pdg_code, particle_data, recomb_model, m_mip_dqdx);
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(
                            pdg_code,
                            particle_data->get_particle_mass(pdg_code),
                            particle_data->pdg_to_name(pdg_code),
                            four_momentum
                        );
                        sg1->particle_info(pinfo);
                    }

                    sg1->dir_weak(true);
                    used_segments.insert(sg1);
                    flag_update = true;
                }
            }
            
            used_vertices.insert(vtx);
        }
    }
}

void PatternAlgorithms::improve_maps_no_dir_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    WireCell::Vector drift_dir(1, 0, 0);
    bool flag_update = true;
    
    while(flag_update) {
        flag_update = false;

        // Iterate through all edges (segments) in the graph.  Stable
        // edge-index order: reclassifying a segment mid-pass changes what
        // later segments sharing a vertex see, so order affects convergence.
        for (const auto& ed : ordered_edges(graph)) {
            SegmentPtr sg = graph[ed].segment;

            // Skip if segment is null or doesn't belong to this cluster
            if (!sg || !sg->cluster() || sg->cluster() != &cluster) continue;

            // Skip showers (trajectory, topology, or electron by dQ/dx)
            // matches prototype get_flag_shower() = flag_shower_trajectory || flag_shower_topology || (particle_type==11)
            if (sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology) ||
                (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11)) continue;

            double length = segment_track_length(sg);

            // Check if segment has no direction, weak direction, or is a proton
            int pdg = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;
            if (sg->dirsign() == 0 || seg_dir_weak(sg) || std::abs(pdg) == 2212) {
                
                auto two_vertices = find_vertices(graph, sg);
                if (!two_vertices.first || !two_vertices.second) continue;
                
                int nshowers[2] = {0, 0};
                int n_in[2] = {0, 0};
                int nmuons[2] = {0, 0};
                int nprotons[2] = {0, 0};
                
                // Get vertex descriptors
                if (!two_vertices.first->descriptor_valid() || !two_vertices.second->descriptor_valid()) continue;
                auto vd1 = two_vertices.first->get_descriptor();
                auto vd2 = two_vertices.second->get_descriptor();
                
                WireCell::Point vtx1_pt = two_vertices.first->wcpt().point;
                WireCell::Point vtx2_pt = two_vertices.second->wcpt().point;
                
                // Count segments at first vertex
                const auto edge_range1 = sorted_out_edges(vd1, graph);
                for (auto e_it : edge_range1) {
                    SegmentPtr sg1 = graph[e_it].segment;
                    if (!sg1) continue;
                    
                    const auto& wcpts = sg1->wcpts();
                    if (wcpts.size() < 2) continue;
                    
                    bool is_shower1 = sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                     sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                     (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11);
                    if (is_shower1) nshowers[0]++;
                    
                    auto front_pt = wcpts.front().point;
                    auto back_pt = wcpts.back().point;
                    double dis_front = ray_length(Ray{vtx1_pt, front_pt});
                    double dis_back = ray_length(Ray{vtx1_pt, back_pt});
                    bool flag_start = (dis_front < dis_back);
                    
                    if ((flag_start && sg1->dirsign() == -1) || (!flag_start && sg1->dirsign() == 1)) {
                        n_in[0]++;
                    }
                    
                    int pdg1 = sg1->has_particle_info() ? sg1->particle_info()->pdg() : 0;
                    if (std::abs(pdg1) == 13) nmuons[0]++;
                    if (std::abs(pdg1) == 2212) nprotons[0]++;
                }
                
                // Count segments at second vertex
                const auto edge_range2 = sorted_out_edges(vd2, graph);
                for (auto e_it : edge_range2) {
                    SegmentPtr sg1 = graph[e_it].segment;
                    if (!sg1) continue;
                    
                    const auto& wcpts = sg1->wcpts();
                    if (wcpts.size() < 2) continue;
                    
                    bool is_shower1 = sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                     sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                     (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11);
                    if (is_shower1) nshowers[1]++;
                    
                    auto front_pt = wcpts.front().point;
                    auto back_pt = wcpts.back().point;
                    double dis_front = ray_length(Ray{vtx2_pt, front_pt});
                    double dis_back = ray_length(Ray{vtx2_pt, back_pt});
                    bool flag_start = (dis_front < dis_back);
                    
                    if ((flag_start && sg1->dirsign() == -1) || (!flag_start && sg1->dirsign() == 1)) {
                        n_in[1]++;
                    }
                    
                    int pdg1 = sg1->has_particle_info() ? sg1->particle_info()->pdg() : 0;
                    if (std::abs(pdg1) == 13) nmuons[1]++;
                    if (std::abs(pdg1) == 2212) nprotons[1]++;
                }
                
                int nvtx1_segs = boost::degree(vd1, graph);
                int nvtx2_segs = boost::degree(vd2, graph);
                
                // Case A: Many showers and very short track
                if ((nshowers[0] + nshowers[1] > 2 && length < 5*units::cm) ||
                    (nshowers[0]+1 == nvtx1_segs && nshowers[1]+1 == nvtx2_segs &&
                     nshowers[0] > 0 && nshowers[1] > 0 && length < 5*units::cm)) {

                    int pdg_code = 11;
                    auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                    sg->particle_info(pinfo);
                    flag_update = true;
                }
                // Case C & D: First/second vertex all showers except current segment (proton)
                else if (nshowers[0]+1 == nvtx1_segs && nshowers[0] >= 2 && pdg == 2212) {
                    WireCell::Vector v1 = segment_cal_dir_3vector(sg, vtx1_pt, 5*units::cm);
                    double min_angle = 180;
                    
                    for (auto e_it : edge_range1) {
                        SegmentPtr sg2 = graph[e_it].segment;
                        if (!sg2 || sg2 == sg) continue;
                        WireCell::Vector v2 = segment_cal_dir_3vector(sg2, vtx1_pt, 5*units::cm);
                        double angle = std::abs(v1.angle(v2) / 3.14159265 * 180.0 - 180.0);
                        if (angle < min_angle) min_angle = angle;
                    }
                    
                    double dQ_dx_rms = segment_rms_dQ_dx(sg);
                    
                    if ((dQ_dx_rms > 1.0 * m_mip_dqdx_median && min_angle < 40) ||
                        (dQ_dx_rms > 0.75 * m_mip_dqdx_median && min_angle < 30) ||
                        (dQ_dx_rms > 0.4 * m_mip_dqdx_median && min_angle < 15)) {
                        
                        const auto& wcpts = sg->wcpts();
                        auto front_pt = wcpts.front().point;
                        auto back_pt = wcpts.back().point;
                        double dis_front = ray_length(Ray{vtx1_pt, front_pt});
                        double dis_back = ray_length(Ray{vtx1_pt, back_pt});
                        bool flag_start = (dis_front < dis_back);

                        if (flag_start)
                            sg->dirsign(-1);
                        else
                            sg->dirsign(1);

                        int pdg_code = 11;
                        auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                        sg->particle_info(pinfo);
                        flag_update = true;
                    }
                }
                else if (nshowers[1]+1 == nvtx2_segs && nshowers[1] >= 2 && pdg == 2212) {
                    WireCell::Vector v1 = segment_cal_dir_3vector(sg, vtx2_pt, 5*units::cm);
                    double min_angle = 180;
                    
                    for (auto e_it : edge_range2) {
                        SegmentPtr sg2 = graph[e_it].segment;
                        if (!sg2 || sg2 == sg) continue;
                        WireCell::Vector v2 = segment_cal_dir_3vector(sg2, vtx2_pt, 5*units::cm);
                        double angle = std::abs(v1.angle(v2) / 3.14159265 * 180.0 - 180.0);
                        if (angle < min_angle) min_angle = angle;
                    }
                    
                    double dQ_dx_rms = segment_rms_dQ_dx(sg);
                    
                    if ((dQ_dx_rms > 1.0 * m_mip_dqdx_median && min_angle < 40) ||
                        (dQ_dx_rms > 0.75 * m_mip_dqdx_median && min_angle < 30) ||
                        (dQ_dx_rms > 0.4 * m_mip_dqdx_median && min_angle < 15)) {
                        
                        const auto& wcpts = sg->wcpts();
                        auto front_pt = wcpts.front().point;
                        auto back_pt = wcpts.back().point;
                        double dis_front = ray_length(Ray{vtx2_pt, front_pt});
                        double dis_back = ray_length(Ray{vtx2_pt, back_pt});
                        bool flag_start = (dis_front < dis_back);

                        if (flag_start)
                            sg->dirsign(-1);
                        else
                            sg->dirsign(1);

                        int pdg_code = 11;
                        auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                        sg->particle_info(pinfo);
                        flag_update = true;
                    }
                }
                // Case E: Muon with specific topology conditions
                else if (std::abs(pdg) == 13 &&
                         ((nprotons[0] >= 0 && nmuons[0] >= 1 && nshowers[1]+1 == nvtx2_segs && nshowers[1] >= 2) ||
                          (nprotons[1] >= 0 && nmuons[1] >= 1 && nshowers[0]+1 == nvtx1_segs && nshowers[0] >= 2) ||
                          (((nprotons[0] >= 0 && nmuons[0] >= 1 && nshowers[1]+1 == nvtx2_segs && nshowers[1] >= 1) ||
                           (nprotons[1] >= 0 && nmuons[1] >= 1 && nshowers[0]+1 == nvtx1_segs && nshowers[0] >= 1)) &&
                          (sg->dirsign() == 0 || seg_dir_weak(sg))))) {
                    
                    double direct_length = segment_track_direct_length(sg);
                    
                    if ((direct_length < 34*units::cm && direct_length < 0.93 * length) ||
                        (length < 5*units::cm && ((nprotons[0] + nshowers[0] == 0 && nshowers[1] >= 2) ||
                                                   (nprotons[1] + nshowers[1] == 0 && nshowers[0] >= 2)))) {

                        // doc sbnd_xin/docs/pr/40 F2 (Case E: muon-topology
                        // demotion; same guard as the other two sites).  Only
                        // guards the CONVERSION, not entry to this branch, so
                        // a spared segment falls through to neither the
                        // conversion NOR the sibling Case F test below (both
                        // are the SAME if/else-if arm as the prototype's
                        // mutually-exclusive case selection).
                        if (!(m_shower_reclass_dqdx_guard && segment_dqdx_spares_electron_reclass(sg, m_mip_dqdx))) {
                            int pdg_code = 11;
                            auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                            sg->particle_info(pinfo);
                            flag_update = true;
                        }
                    }
                    // Case F: Check daughter showers
                    else if ((((nshowers[0]+nshowers[1] >= 2) && (nprotons[0]+nmuons[0]+nshowers[0] == 1 || nprotons[1]+nmuons[1]+nshowers[1] == 1)) ||
                              ((nshowers[0]+nshowers[1] >= 1) && (nprotons[0]+nmuons[0]+nshowers[0] > 1 || nprotons[1]+nmuons[1]+nshowers[1] > 1))) &&
                             length < 40*units::cm) {
                        
                        int num_s1 = 0, num_s2 = 0;
                        double length_s1 = 0, length_s2 = 0;
                        double max_angle1 = 0, max_angle2 = 0;
                        
                        WireCell::Vector dir1 = segment_cal_dir_3vector(sg, vtx1_pt, 15*units::cm);
                        for (auto e_it : edge_range1) {
                            SegmentPtr sg1 = graph[e_it].segment;
                            if (!sg1 || sg1 == sg) continue;
                            
                            WireCell::Vector dir2 = segment_cal_dir_3vector(sg1, vtx1_pt, 15*units::cm);
                            bool is_shower1 = sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                             sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                             (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11);
                            if (is_shower1) {
                                double angle = dir1.angle(dir2) / 3.14159265 * 180.0;
                                if (max_angle1 < angle) max_angle1 = angle;

                                auto pair_result = calculate_num_daughter_showers(graph, two_vertices.first, sg1);
                                num_s1 += pair_result.first;
                                length_s1 += pair_result.second;
                            }
                        }
                        
                        dir1 = segment_cal_dir_3vector(sg, vtx2_pt, 10*units::cm);
                        for (auto e_it : edge_range2) {
                            SegmentPtr sg1 = graph[e_it].segment;
                            if (!sg1 || sg1 == sg) continue;
                            
                            WireCell::Vector dir2 = segment_cal_dir_3vector(sg1, vtx2_pt, 15*units::cm);
                            bool is_shower1 = sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                             sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                             (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11);
                            if (is_shower1) {
                                double angle = dir1.angle(dir2) / 3.14159265 * 180.0;
                                if (max_angle2 < angle) max_angle2 = angle;

                                auto pair_result = calculate_num_daughter_showers(graph, two_vertices.second, sg1);
                                num_s2 += pair_result.first;
                                length_s2 += pair_result.second;
                            }
                        }
                        
                        if (((num_s1 >= 4 || (length_s1 > 50*units::cm && num_s1 >= 2)) && max_angle1 > 150) ||
                            ((num_s2 >= 4 || length_s2 > 50*units::cm) && max_angle2 > 150) ||
                            (length < 6*units::cm && ((num_s1 >= 4 && length_s1 > 20*units::cm) ||
                                                      (num_s2 >= 4 && length_s2 > 20*units::cm)))) {
                            
                            int pdg_code = 11;
                            auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                            sg->particle_info(pinfo);
                            flag_update = true;
                        }
                    }
                }
                // Case G: Muon with specific vertex connectivity
                else if (std::abs(pdg) == 13 && (sg->dirsign() == 0 || seg_dir_weak(sg)) &&
                         ((nmuons[0]+nprotons[0]+nshowers[0] == 1) || (nmuons[1]+nprotons[1]+nshowers[1] == 1)) &&
                         (nshowers[0] + nshowers[1] > 0 || segment_median_dQ_dx(sg) < 1.3*m_mip_dqdx_median)) {
                    
                    bool flag_change = false;
                    
                    if (nvtx1_segs == 2) {
                        SegmentPtr tmp_sg = nullptr;
                        for (auto e_it : edge_range1) {
                            SegmentPtr candidate = graph[e_it].segment;
                            if (candidate && candidate != sg) {
                                tmp_sg = candidate;
                                break;
                            }
                        }
                        if (tmp_sg) {
                            int tmp_pdg = tmp_sg->has_particle_info() ? tmp_sg->particle_info()->pdg() : 0;
                            if (tmp_pdg == 13 && segment_track_length(tmp_sg) > 4*length && length < 8*units::cm) {
                                flag_change = true;
                            }
                        }
                    } else if (nvtx2_segs == 2) {
                        SegmentPtr tmp_sg = nullptr;
                        for (auto e_it : edge_range2) {
                            SegmentPtr candidate = graph[e_it].segment;
                            if (candidate && candidate != sg) {
                                tmp_sg = candidate;
                                break;
                            }
                        }
                        if (tmp_sg) {
                            int tmp_pdg = tmp_sg->has_particle_info() ? tmp_sg->particle_info()->pdg() : 0;
                            if (tmp_pdg == 13 && segment_track_length(tmp_sg) > 4*length && length < 8*units::cm) {
                                flag_change = true;
                            }
                        }
                    }
                    
                    if (flag_change) {
                        int pdg_code = 11;
                        auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                        sg->particle_info(pinfo);
                        flag_update = true;
                    }
                }

                // Case B: Setting direction for segments between shower vertices
                if (((nshowers[0]+1 == nvtx1_segs) || nshowers[0] > 0) &&
                    ((nshowers[1]+1 == nvtx2_segs) || nshowers[1] > 0) &&
                    (nshowers[0] + nshowers[1] > 2) &&
                    ((nshowers[0]+1 == nvtx1_segs && nshowers[0] > 0) ||
                     (nshowers[1]+1 == nvtx2_segs && nshowers[1] > 0))) {
                    
                    if ((length < 25*units::cm && pdg != 11) || sg->dirsign() == 0) {
                        const auto& wcpts = sg->wcpts();
                        auto front_pt = wcpts.front().point;
                        auto back_pt = wcpts.back().point;
                        double dis_front = ray_length(Ray{vtx1_pt, front_pt});
                        double dis_back = ray_length(Ray{vtx1_pt, back_pt});
                        bool flag_start = (dis_front < dis_back);
                        
                        if (flag_start) {
                            if (nshowers[1] == 0) {
                                sg->dirsign(-1);
                            } else if (nshowers[0] == 0) {
                                sg->dirsign(1);
                            }
                        } else {
                            if (nshowers[1] == 0) {
                                sg->dirsign(1);
                            } else if (nshowers[0] == 0) {
                                sg->dirsign(-1);
                            }
                        }
                        sg->dir_weak(true);

                        int pdg_code = 11;
                        auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                        sg->particle_info(pinfo);
                        flag_update = true;
                    }
                }
                // Case H: No particle type, short length, high dQ/dx, has showers
                else if (pdg == 0 && length < 12*units::cm && 
                         (nshowers[0] + nshowers[1] > 0) && 
                         segment_median_dQ_dx(sg)/m_mip_dqdx_median > 1.2) {
                    
                    bool flag_change = false;
                    
                    auto pair_result1 = calculate_num_daughter_showers(graph, two_vertices.second, sg);
                    auto pair_result2 = calculate_num_daughter_showers(graph, two_vertices.first, sg);
                    
                    if (pair_result1.first > 2) {
                        WireCell::Vector v1 = segment_cal_dir_3vector(sg, vtx1_pt, 10*units::cm);
                        double min_angle = 180;
                        double para_angle = 90;
                        
                        for (auto e_it : edge_range1) {
                            SegmentPtr sg2 = graph[e_it].segment;
                            if (!sg2 || sg2 == sg) continue;
                            bool is_shower2 = sg2->flags_any(SegmentFlags::kShowerTrajectory) ||
                                             sg2->flags_any(SegmentFlags::kShowerTopology) ||
                                             (sg2->has_particle_info() && std::abs(sg2->particle_info()->pdg()) == 11);
                            if (!is_shower2) continue;

                            WireCell::Vector v2 = segment_cal_dir_3vector(sg2, vtx1_pt, 10*units::cm);
                            double angle = std::abs(v1.angle(v2) / 3.14159265 * 180.0 - 180.0);
                            if (angle < min_angle) {
                                min_angle = angle;
                                para_angle = std::abs(v2.angle(drift_dir) / 3.14159265 * 180.0 - 90);
                            }
                        }
                        
                        if (min_angle < 25 || 
                            (std::abs(v1.angle(drift_dir) / 3.14159265 * 180.0 - 90) < 10 && 
                             para_angle < 30 && min_angle < 45)) {
                            flag_change = true;
                        }
                    }
                    
                    if (!flag_change && pair_result2.first > 2) {
                        WireCell::Vector v1 = segment_cal_dir_3vector(sg, vtx2_pt, 10*units::cm);
                        double min_angle = 180;
                        double para_angle = 90;
                        
                        for (auto e_it : edge_range2) {
                            SegmentPtr sg2 = graph[e_it].segment;
                            if (!sg2 || sg2 == sg) continue;
                            bool is_shower2 = sg2->flags_any(SegmentFlags::kShowerTrajectory) ||
                                             sg2->flags_any(SegmentFlags::kShowerTopology) ||
                                             (sg2->has_particle_info() && std::abs(sg2->particle_info()->pdg()) == 11);
                            if (!is_shower2) continue;

                            WireCell::Vector v2 = segment_cal_dir_3vector(sg2, vtx2_pt, 10*units::cm);
                            double angle = std::abs(v1.angle(v2) / 3.14159265 * 180.0 - 180.0);
                            if (angle < min_angle) {
                                min_angle = angle;
                                para_angle = std::abs(v2.angle(drift_dir) / 3.14159265 * 180.0 - 90);
                            }
                        }
                        
                        if (min_angle < 25 || 
                            (std::abs(v1.angle(drift_dir) / 3.14159265 * 180.0 - 90) < 10 && 
                             para_angle < 10 && min_angle < 45)) {
                            flag_change = true;
                        }
                    }
                    
                    if (flag_change) {
                        int pdg_code = 11;
                        auto pinfo = reclass_pinfo(sg, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                        sg->particle_info(pinfo);
                        flag_update = true;
                    }
                }

            } // end if no direction or weak or proton
        } // loop over all segments
    } // while flag_update
}

void PatternAlgorithms::improve_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    bool flag_update = true;
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    while(flag_update) {
        flag_update = false;
        
        // Iterate through all vertices in the graph
        for (auto vit : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vit].vertex;

            // Skip if vertex is null or doesn't belong to this cluster
            if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;

            // Skip if vertex has only 1 segment
            if (!vtx->descriptor_valid()) continue;
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) <= 1) continue;
            
            // Skip if already processed
            if (used_vertices.find(vtx) != used_vertices.end()) continue;
            
            int n_in = 0;
            int n_in_shower = 0;
            std::vector<SegmentPtr> in_tracks;
            
            // Get vertex position
            WireCell::Point vtx_point = vtx->wcpt().point;
            
            // Iterate through all segments connected to this vertex
            const auto edge_range = sorted_out_edges(vd, graph);
            for (auto eit : edge_range) {
                SegmentPtr sg = graph[eit].segment;
                if (!sg) continue;
                
                // Determine if this vertex is at the front or back of the segment
                const auto& wcpts = sg->wcpts();
                if (wcpts.size() < 2) continue;
                
                auto front_pt = wcpts.front().point;
                auto back_pt = wcpts.back().point;
                
                double dis_front = ray_length(Ray{vtx_point, front_pt});
                double dis_back = ray_length(Ray{vtx_point, back_pt});
                
                bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
                
                // Check if this is an "incoming" segment (pointing into vertex)
                if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                    n_in++;
                    
                    // matches prototype get_flag_shower() = kShowerTrajectory || kShowerTopology || particle_type==11
                    bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                    sg->flags_any(SegmentFlags::kShowerTopology) ||
                                    (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
                    if (is_shower) {
                        n_in_shower++;
                    } else {
                        in_tracks.push_back(sg);
                    }
                }
            }
            
            // If there are multiple incoming segments and not all are showers
            if (n_in > 1 && n_in != n_in_shower) {
                // Reclassify all incoming tracks as electrons
                for (auto it1 = in_tracks.begin(); it1 != in_tracks.end(); it1++) {
                    SegmentPtr sg1 = *it1;
                    
                    int pdg_code = 11;
                    auto pinfo = reclass_pinfo(sg1, pdg_code, particle_data, recomb_model, m_mip_dqdx, m_reclass_preserve_4mom, true, m_reclass_never_computed_ke_floor);
                    sg1->particle_info(pinfo);
                    flag_update = true;
                }
            }
            
            used_vertices.insert(vtx);
        } // loop over all vertices
    } // while flag_update
}

void PatternAlgorithms::judge_no_dir_tracks_close_to_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, IDetectorVolumes::pointer dv){
    std::vector<SegmentPtr> shower_set;
    std::vector<SegmentPtr> no_dir_track_set;

    // Collect shower segments and no-direction track segments
    for (auto e : ordered_edges(graph)) {
        SegmentPtr sg = graph[e].segment;
        if (!sg || sg->cluster() != &cluster) continue;

        // matches prototype get_flag_shower() = kShowerTrajectory || kShowerTopology || particle_type==11
        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                        sg->flags_any(SegmentFlags::kShowerTopology) ||
                        (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);

        if (is_shower) {
            shower_set.push_back(sg);
        } else {
            if (sg->dirsign() == 0) {
                no_dir_track_set.push_back(sg);
            }
        }
    }

    // Process each no-direction track segment
    for (auto sg : no_dir_track_set) {
        bool flag_change = true;
        
        const auto& pts = sg->fits();//wcpts();
        
        // Check each point in the segment
        for (size_t i = 0; i < pts.size(); i++) {
            WireCell::Point test_p = pts.at(i).point;
            
            // Get apa and face for this point
            auto test_wpid = dv->contained_by(test_p);
            if (test_wpid.apa() == -1 || test_wpid.face() == -1) {
                flag_change = false;
                break;
            }
            
            int apa = test_wpid.apa();
            int face = test_wpid.face();
            
            double min_u_dis = 1e9;
            double min_v_dis = 1e9;
            double min_w_dis = 1e9;
            
            // Find minimum 2D distances to all shower segments
            for (auto it1 = shower_set.begin(); it1 != shower_set.end(); it1++) {
                auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(*it1, test_p, apa, face, "fit");
                
                if (dist_u < min_u_dis) min_u_dis = dist_u;
                if (dist_v < min_v_dis) min_v_dis = dist_v;
                if (dist_w < min_w_dis) min_w_dis = dist_w;
            }
            
            // If any distance exceeds threshold, don't reclassify
            if (min_u_dis > 0.6*units::cm || min_v_dis > 0.6*units::cm || min_w_dis > 0.6*units::cm) {
                flag_change = false;
                break;
            }
        }
        
        // Reclassify segment as electron if all points are close to showers
        if (flag_change) {
            int pdg_code = 11;
            // doc pr/31 §12 F1 shape C: prototype writes type and mass only
            // (NeutrinoID_track_shower.h:1238-1239), 4-momentum untouched.
            if (m_reclass_preserve_4mom) {
                sg->particle_info(reclass_pinfo(sg, pdg_code, particle_data, m_recomb_model, m_mip_dqdx, true, false, m_reclass_never_computed_ke_floor));
            }
            else {
                double em_mass = particle_data->get_particle_mass(pdg_code);
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    pdg_code,
                    em_mass,
                    particle_data->pdg_to_name(pdg_code),
                    WireCell::D4Vector<double>(em_mass, 0, 0, 0)
                );
                sg->particle_info(pinfo);
            }
        }
    }
}

bool PatternAlgorithms::examine_maps(Graph&graph, Facade::Cluster& cluster){
    bool flag_return = true;
    
    // Iterate through all vertices in the graph, in stable node-index order.
    for (const auto& vd_it : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd_it].vertex;

        // Skip if vertex is null or doesn't belong to this cluster
        if (!vtx || vtx->cluster() != &cluster) continue;

        // Skip vertices with only 1 segment
        if (!vtx->descriptor_valid()) continue;
        auto vd = vtx->get_descriptor();
        if (boost::degree(vd, graph) <= 1) continue;
        
        int n_in = 0;
        int n_in_shower = 0;
        int n_out_tracks = 0;
        
        // Get vertex position
        WireCell::Point vtx_point = vtx->wcpt().point;
        
        // Iterate through all segments connected to this vertex
        const auto edge_range = sorted_out_edges(vd, graph);
        for (auto eit : edge_range) {
            SegmentPtr sg = graph[eit].segment;
            if (!sg) continue;
            
            // Determine if this vertex is at the front or back of the segment
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            auto front_pt = wcpts.front().point;
            auto back_pt = wcpts.back().point;
            
            double dis_front = ray_length(Ray{vtx_point, front_pt});
            double dis_back = ray_length(Ray{vtx_point, back_pt});
            
            bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
            
            // matches prototype get_flag_shower() = kShowerTrajectory || kShowerTopology || particle_type==11
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                           sg->flags_any(SegmentFlags::kShowerTopology) ||
                           (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);

            // Check if this is an "incoming" segment (pointing into vertex)
            if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                n_in++;
                if (is_shower) {
                    n_in_shower++;
                }
            }

            // Check if this is an "outgoing" track (pointing away from vertex)
            if ((flag_start && sg->dirsign() == 1) || (!flag_start && sg->dirsign() == -1)) {
                if (!is_shower) {
                    n_out_tracks++;
                }
            }
        }
        
        // Check for violations
        if (n_in > 1 && n_in != n_in_shower) {
            SPDLOG_LOGGER_TRACE(s_log, "examine_maps: Wrong: Multiple ({}) particles into a vertex!", n_in);
            print_segs_info(graph, cluster, vtx);
            flag_return = false;
        }
        
        if (n_in_shower > 0 && n_out_tracks > 0) {
            SPDLOG_LOGGER_TRACE(s_log, "examine_maps: Wrong: {} showers in and {} tracks out!", n_in_shower, n_out_tracks);
            print_segs_info(graph, cluster, vtx);
            flag_return = false;
        }
    }
    
    return flag_return;
}

void PatternAlgorithms::examine_all_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data){
    int n_good_tracks = 0, n_tracks = 0, n_showers = 0;
    double length_good_tracks = 0, length_tracks = 0, length_showers = 0;
    double tracks_score = 0;
    SegmentPtr good_track = nullptr;
    
    double maximal_length = 0;
    SegmentPtr maximal_length_track = nullptr;
    
    // Count segments and their properties.  Stable edge-index order: the
    // length sums FP-accumulate and maximal_length_track is a tie-broken
    // argmax, both consumed by classification thresholds below.
    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &cluster) continue;

        double length = segment_track_length(sg);
        // matches prototype get_flag_shower() = kShowerTrajectory || kShowerTopology || particle_type==11
        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                        sg->flags_any(SegmentFlags::kShowerTopology) ||
                        (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);

        if (is_shower) {
            n_showers++;
            length_showers += length;
        } else {
            if (sg->dirsign() != 0 && !seg_dir_weak(sg)) {
                good_track = sg;
                n_good_tracks++;
                length_good_tracks += length;
            } else {
                n_tracks++;
                length_tracks += length;
                if (length > maximal_length) {
                    maximal_length = length;
                    maximal_length_track = sg;
                }
                if (sg->particle_score() != 100) tracks_score += sg->particle_score();
            }
        }
    }
    
    if (n_good_tracks + n_tracks + n_showers == 1) return;
    
    // If there is only one good track
    if (n_good_tracks == 1 && (length_good_tracks < 0.15 * (length_showers + length_tracks)) && length_good_tracks < 10*units::cm) {
        auto pair_vertices = find_vertices(graph, good_track);
        // doc sbnd_xin/docs/pr/31 §12 (F7, was P5; pr/30 F4's sibling).  The
        // two acceptance branches below are ASYMMETRIC -- the .second-side
        // branch has a relaxed 150-degree clause the .first-side branch lacks,
        // faithfully ported from the prototype -- but which physical vertex
        // lands in .first differs: prototype orders by vertex id
        // (NeutrinoID_proto_vertex.h:3227-3243), the toolkit by proximity to
        // the segment's first fit point (PRGraph.cxx find_vertices).  The
        // owner's positional clearance covers .first meaning "an end"; it does
        // not cover .first selecting a branch.  ON: order the pair by
        // get_graph_index() -- a creation-order counter, the SHAPE of the
        // prototype's get_id(), but the two trees create vertices in different
        // orders, so this restores A deterministic topological convention, not
        // provably the prototype's.  Site-local on purpose: reordering
        // find_vertices itself is pr/30 F4's decision and hits every caller.
        // OFF = today's proximity order = byte-identical.
        if (m_examine_showers_vertex_by_index && pair_vertices.first && pair_vertices.second
            && pair_vertices.first->get_graph_index() > pair_vertices.second->get_graph_index()) {
            std::swap(pair_vertices.first, pair_vertices.second);
        }

        int num_s1 = 0, num_s2 = 0;
        double length_s1 = 0, length_s2 = 0;
        
        auto pair_result1 = calculate_num_daughter_showers(graph, pair_vertices.first, good_track);
        auto pair_result2 = calculate_num_daughter_showers(graph, pair_vertices.second, good_track);
        num_s1 = pair_result1.first;
        length_s1 = pair_result1.second;
        num_s2 = pair_result2.first;
        length_s2 = pair_result2.second;
        
        if (num_s1 > 0 && length_s1 > length_good_tracks) {
            double max_angle = 0;
            WireCell::Point vtx2_pt = pair_vertices.second->fit().valid() ? pair_vertices.second->fit().point : pair_vertices.second->wcpt().point;
            WireCell::Vector dir1 = segment_cal_dir_3vector(good_track, vtx2_pt, 15*units::cm);
            
            if (pair_vertices.second->descriptor_valid()) {
                auto vd2 = pair_vertices.second->get_descriptor();
                const auto edge_range = sorted_out_edges(vd2, graph);
                for (auto e_it : edge_range) {
                    SegmentPtr sg1 = graph[e_it].segment;
                    if (!sg1 || sg1 == good_track) continue;
                    
                    WireCell::Vector dir2 = segment_cal_dir_3vector(sg1, vtx2_pt, 15*units::cm);
                    double angle = dir1.angle(dir2) / 3.14159265 * 180.0;
                    if (max_angle < angle) max_angle = angle;
                }
            }
            
            if (max_angle > 165 || (max_angle > 150 && length_good_tracks < 3.0*units::cm && length_good_tracks < 0.1 * length_showers)) {
                n_good_tracks = 0;
                length_tracks += length_good_tracks;
            }
        }
        
        if (num_s2 > 0 && length_s2 > length_good_tracks && n_good_tracks > 0) {
            double max_angle = 0;
            WireCell::Point vtx1_pt = pair_vertices.first->fit().valid() ? pair_vertices.first->fit().point : pair_vertices.first->wcpt().point;
            WireCell::Vector dir1 = segment_cal_dir_3vector(good_track, vtx1_pt, 15*units::cm);
            
            if (pair_vertices.first->descriptor_valid()) {
                auto vd1 = pair_vertices.first->get_descriptor();
                const auto edge_range = sorted_out_edges(vd1, graph);
                for (auto e_it : edge_range) {
                    SegmentPtr sg1 = graph[e_it].segment;
                    if (!sg1 || sg1 == good_track) continue;
                    
                    WireCell::Vector dir2 = segment_cal_dir_3vector(sg1, vtx1_pt, 15*units::cm);
                    double angle = dir1.angle(dir2) / 3.14159265 * 180.0;
                    if (max_angle < angle) max_angle = angle;
                }
            }
            
            if (max_angle > 165) {
                n_good_tracks = 0;
                length_tracks += length_good_tracks;
            }
        }
        
        // Check vertex connectivity and beam angle
        if (pair_vertices.first && pair_vertices.second) {
            int nvtx1_segs = 0, nvtx2_segs = 0;
            if (pair_vertices.first->descriptor_valid()) {
                nvtx1_segs = boost::degree(pair_vertices.first->get_descriptor(), graph);
            }
            if (pair_vertices.second->descriptor_valid()) {
                nvtx2_segs = boost::degree(pair_vertices.second->get_descriptor(), graph);
            }
            
            if (nvtx1_segs == 1 && nvtx2_segs > 1) {
                double max_length = 0;
                SegmentPtr max_segment = nullptr;
                
                if (pair_vertices.second->descriptor_valid()) {
                    auto vd2 = pair_vertices.second->get_descriptor();
                    // Stable edge-index order: tie-broken argmax below.
                    for (auto edesc : sorted_out_edges(vd2, graph)) {
                        SegmentPtr sg = graph[edesc].segment;
                        if (!sg) continue;

                        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                       sg->flags_any(SegmentFlags::kShowerTopology) ||
                                       (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
                        if (is_shower) {
                            double length = segment_track_length(sg);
                            if (length > max_length) {
                                max_length = length;
                                max_segment = sg;
                            }
                        }
                    }
                }

                if (max_segment != nullptr && max_length > 5*units::cm) {
                    WireCell::Point vtx2_pt = pair_vertices.second->fit().valid() ? pair_vertices.second->fit().point : pair_vertices.second->wcpt().point;
                    WireCell::Vector dir = segment_cal_dir_3vector(max_segment, vtx2_pt, 15*units::cm);
                    WireCell::Vector beam_dir(0, 0, 1);
                    if (beam_dir.angle(dir) / 3.14159265 * 180.0 > 90) {
                        n_good_tracks = 0;
                        length_tracks += length_good_tracks;
                    }
                }
            } else if (nvtx1_segs > 1 && nvtx2_segs == 1) {
                double max_length = 0;
                SegmentPtr max_segment = nullptr;
                
                if (pair_vertices.first->descriptor_valid()) {
                    auto vd1 = pair_vertices.first->get_descriptor();
                    // Stable edge-index order: tie-broken argmax below.
                    for (auto edesc : sorted_out_edges(vd1, graph)) {
                        SegmentPtr sg = graph[edesc].segment;
                        if (!sg) continue;

                        bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                       sg->flags_any(SegmentFlags::kShowerTopology) ||
                                       (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);
                        if (is_shower) {
                            double length = segment_track_length(sg);
                            if (length > max_length) {
                                max_length = length;
                                max_segment = sg;
                            }
                        }
                    }
                }

                if (max_segment != nullptr && max_length > 5*units::cm) {
                    WireCell::Point vtx1_pt = pair_vertices.first->fit().valid() ? pair_vertices.first->fit().point : pair_vertices.first->wcpt().point;
                    WireCell::Vector dir = segment_cal_dir_3vector(max_segment, vtx1_pt, 15*units::cm);
                    WireCell::Vector beam_dir(0, 0, 1);
                    if (beam_dir.angle(dir) / 3.14159265 * 180.0 > 90) {
                        n_good_tracks = 0;
                        length_tracks += length_good_tracks;
                    }
                }
            }
        }
    } else if (n_good_tracks == 0 && (n_tracks == 2 && length_tracks <= 35*units::cm)) {
        if (maximal_length_track != nullptr) {
            auto pair_vertices = find_vertices(graph, maximal_length_track);
            
            if (pair_vertices.first && pair_vertices.second) {
                int nvtx1_segs = 0, nvtx2_segs = 0;
                if (pair_vertices.first->descriptor_valid()) {
                    nvtx1_segs = boost::degree(pair_vertices.first->get_descriptor(), graph);
                }
                if (pair_vertices.second->descriptor_valid()) {
                    nvtx2_segs = boost::degree(pair_vertices.second->get_descriptor(), graph);
                }
                
                if (nvtx1_segs < nvtx2_segs) {
                    WireCell::Point vtx2_pt = pair_vertices.second->fit().valid() ? pair_vertices.second->fit().point : pair_vertices.second->wcpt().point;
                    WireCell::Vector dir = segment_cal_dir_3vector(maximal_length_track, vtx2_pt, 15*units::cm);
                    WireCell::Vector beam_dir(0, 0, 1);
                    if (beam_dir.angle(dir) / 3.14159265 * 180.0 > 100) {
                        n_tracks--;
                        length_tracks -= maximal_length;
                        n_showers++;
                        length_showers += maximal_length;
                    }
                } else if (nvtx1_segs > nvtx2_segs) {
                    WireCell::Point vtx1_pt = pair_vertices.first->fit().valid() ? pair_vertices.first->fit().point : pair_vertices.first->wcpt().point;
                    WireCell::Vector dir = segment_cal_dir_3vector(maximal_length_track, vtx1_pt, 15*units::cm);
                    WireCell::Vector beam_dir(0, 0, 1);
                    if (beam_dir.angle(dir) / 3.14159265 * 180.0 > 90) {
                        n_tracks--;
                        length_tracks -= maximal_length;
                        n_showers++;
                        length_showers += maximal_length;
                    }
                }
            }
        }
    }
    
    bool flag_change_showers = false;

    // Check main_cluster status
    bool is_main_cluster = cluster.get_flag(Facade::Flags::main_cluster);
    
    if (n_good_tracks == 0) {
        if (length_tracks < 1.0/3.0 * length_showers || (length_tracks < 2.0/3.0 * length_showers && n_tracks == 1)) {
            if ((length_showers + length_tracks) < 40*units::cm) {
                flag_change_showers = true;
            } else if (length_tracks < 0.18 * length_showers && ((length_showers + length_tracks) < 60*units::cm || length_tracks < 12*units::cm)) {
                flag_change_showers = true;
            } else if (length_tracks < 0.25 * length_showers && ((tracks_score == 0 && length_tracks < 30*units::cm) || length_tracks < 10*units::cm)) {
                flag_change_showers = true;
            } else if (n_tracks == 1 && tracks_score == 0 && length_tracks < 15*units::cm && length_tracks < 1.0/3.0 * length_showers) {
                flag_change_showers = true;
            }
        } else if ((length_tracks < 35*units::cm && length_tracks + length_showers < 50*units::cm && length_showers < 15*units::cm) && 
          (!is_main_cluster || 
                (is_main_cluster && 
                 (length_showers > 0.5*length_tracks ||
                  (length_showers > 0.3*length_tracks && n_showers >= 2) ||
                  (n_showers == 1 && n_tracks == 1 && length_showers > length_tracks * 0.3) ||
                  tracks_score == 0)))) {
                flag_change_showers = true;
                if (length_showers == 0 && n_tracks <= 2 && (is_main_cluster || length_tracks > 15*units::cm)) {
                    flag_change_showers = false;
                }
        } else if (length_tracks < 35*units::cm && length_tracks + length_showers < 50*units::cm && length_showers < 15*units::cm) {

            // matches prototype: set true then verify each non-shower touches a shower neighbor
            flag_change_showers = true;
            for (const auto& ed : ordered_edges(graph)) {
                SegmentPtr sg = graph[ed].segment;
                if (!sg || sg->cluster() != &cluster) continue;

                bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                               sg->flags_any(SegmentFlags::kShowerTopology) ||
                               (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);

                if (!is_shower) {
                    auto pair_vertices = find_vertices(graph, sg);
                    bool flag_shower = false;

                    if (pair_vertices.first && pair_vertices.first->descriptor_valid()) {
                        auto vd1 = pair_vertices.first->get_descriptor();
                        const auto edge_range = sorted_out_edges(vd1, graph);
                        for (auto e_it : edge_range) {
                            SegmentPtr sg1 = graph[e_it].segment;
                            if (sg1 && (sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                       sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                       (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11))) {
                                flag_shower = true;
                                break;
                            }
                        }
                    }

                    if (!flag_shower && pair_vertices.second && pair_vertices.second->descriptor_valid()) {
                        auto vd2 = pair_vertices.second->get_descriptor();
                        const auto edge_range = sorted_out_edges(vd2, graph);
                        for (auto e_it : edge_range) {
                            SegmentPtr sg1 = graph[e_it].segment;
                            if (sg1 && (sg1->flags_any(SegmentFlags::kShowerTrajectory) ||
                                       sg1->flags_any(SegmentFlags::kShowerTopology) ||
                                       (sg1->has_particle_info() && std::abs(sg1->particle_info()->pdg()) == 11))) {
                                flag_shower = true;
                                break;
                            }
                        }
                    }

                    if (!flag_shower) {
                        flag_change_showers = false;
                        break;
                    }
                }
            }
        }
    }
    
    if (flag_change_showers) {
        for (const auto& ed : ordered_edges(graph)) {
            SegmentPtr sg = graph[ed].segment;
            if (!sg || sg->cluster() != &cluster) continue;

            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) ||
                           sg->flags_any(SegmentFlags::kShowerTopology) ||
                           (sg->has_particle_info() && std::abs(sg->particle_info()->pdg()) == 11);

            // doc sbnd_xin/docs/pr/40 F2 (same guard as improve_maps_shower_
            // in_track_out).  false = legacy = every non-shower segment in a
            // shower-dominated cluster becomes electron unconditionally.
            //
            // doc sbnd_xin/docs/pr/40 round 10: shower_bragg_protect_start_
            // segment is an additive sibling spare -- see
            // segment_bragg_spares_electron_reclass's header comment.  false
            // = legacy = byte-identical.
            if (!is_shower &&
                ((m_shower_reclass_dqdx_guard && segment_dqdx_spares_electron_reclass(sg, m_mip_dqdx)) ||
                 (m_shower_bragg_protect_start_segment && segment_bragg_spares_electron_reclass(sg)))) {
                continue;
            }

            if (!is_shower) {
                int pdg_code = 11;
                // doc pr/31 §12 F1 shape C: prototype writes type and mass only
                // (NeutrinoID_track_shower.h:313-315), 4-momentum untouched.
                if (m_reclass_preserve_4mom) {
                    sg->particle_info(reclass_pinfo(sg, pdg_code, particle_data, m_recomb_model, m_mip_dqdx, true, false, m_reclass_never_computed_ke_floor));
                }
                else {
                    double electron_mass = particle_data->get_particle_mass(pdg_code);
                    auto pinfo = std::make_shared<Aux::ParticleInfo>(
                        pdg_code,
                        electron_mass,
                        particle_data->pdg_to_name(pdg_code),
                        WireCell::D4Vector<double>(electron_mass, 0, 0, 0)
                    );
                    sg->particle_info(pinfo);
                }
            }
        }
    }
}



void PatternAlgorithms::shower_determining_in_main_cluster(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, IDetectorVolumes::pointer dv){
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    // Examine good tracks first
    examine_good_tracks(graph, cluster, particle_data);
    MS t_examine_good_tracks(Clock::now() - t0); t0 = Clock::now();

    // If multiple tracks in, make them undetermined 
    fix_maps_multiple_tracks_in(graph, cluster);
    MS t_fix_multiple_tracks_in(Clock::now() - t0); t0 = Clock::now();

    // If one shower in and a good track out, reverse the shower
    fix_maps_shower_in_track_out(graph, cluster);
    MS t_fix_shower_in_track_out_1(Clock::now() - t0); t0 = Clock::now();

    // If there is one good track in, turn everything else to out
    improve_maps_one_in(graph, cluster, particle_data, recomb_model);
    MS t_improve_one_in(Clock::now() - t0); t0 = Clock::now();

    // If one shower in and a track out, change the track to shower
    improve_maps_shower_in_track_out(graph, cluster, particle_data, recomb_model);
    MS t_improve_shower_in_track_out_1(Clock::now() - t0); t0 = Clock::now();

    // Help to change tracks around shower to showers
    improve_maps_no_dir_tracks(graph, cluster, particle_data, recomb_model);
    MS t_improve_no_dir_tracks(Clock::now() - t0); t0 = Clock::now();

    // If one shower in and a track out, change the track to shower (no reverse flag)
    improve_maps_shower_in_track_out(graph, cluster, particle_data, recomb_model, false);
    MS t_improve_shower_in_track_out_2(Clock::now() - t0); t0 = Clock::now();

    // If multiple tracks in, change track to shower
    improve_maps_multiple_tracks_in(graph, cluster, particle_data, recomb_model);
    MS t_improve_multiple_tracks_in(Clock::now() - t0); t0 = Clock::now();

    // If one shower in and a good track out, reverse the shower
    fix_maps_shower_in_track_out(graph, cluster);
    MS t_fix_shower_in_track_out_2(Clock::now() - t0); t0 = Clock::now();

    // Judgement for no-direction tracks close to showers
    judge_no_dir_tracks_close_to_showers(graph, cluster, particle_data, dv);
    MS t_judge_no_dir_tracks(Clock::now() - t0); t0 = Clock::now();

    // Examine maps for physics violations
    examine_maps(graph, cluster);
    MS t_examine_maps(Clock::now() - t0); t0 = Clock::now();

    // Examine all showers comprehensively
    examine_all_showers(graph, cluster, particle_data);

    

    MS t_examine_all_showers(Clock::now() - t0);

    if (m_perf) {
        MS t_total_ms(Clock::now() - t_total);
        SPDLOG_LOGGER_TRACE(s_log,
            "shower_determining_in_main_cluster timing: "
            "examine_good_tracks={:.3f}ms fix_multiple_tracks_in={:.3f}ms "
            "fix_shower_in_track_out_1={:.3f}ms improve_one_in={:.3f}ms "
            "improve_shower_in_track_out_1={:.3f}ms improve_no_dir_tracks={:.3f}ms "
            "improve_shower_in_track_out_2={:.3f}ms improve_multiple_tracks_in={:.3f}ms "
            "fix_shower_in_track_out_2={:.3f}ms judge_no_dir_tracks={:.3f}ms "
            "examine_maps={:.3f}ms examine_all_showers={:.3f}ms ",
            t_examine_good_tracks.count(), t_fix_multiple_tracks_in.count(),
            t_fix_shower_in_track_out_1.count(), t_improve_one_in.count(),
            t_improve_shower_in_track_out_1.count(), t_improve_no_dir_tracks.count(),
            t_improve_shower_in_track_out_2.count(), t_improve_multiple_tracks_in.count(),
            t_fix_shower_in_track_out_2.count(), t_judge_no_dir_tracks.count(),
            t_examine_maps.count(), t_examine_all_showers.count());
        SPDLOG_LOGGER_TRACE(s_log, "shower_determining_in_main_cluster timing:  TOTAL={:.3f}ms", t_total_ms.count());

        // // Print final per-segment state, matching determine_direction format
        // auto [ebegin2, eend2] = boost::edges(graph);
        // for (auto eit = ebegin2; eit != eend2; ++eit) {
        //     SegmentPtr seg = graph[*eit].segment;
        //     if (!seg || seg->cluster() != &cluster) continue;

        //     std::string seg_type;
        //     if (seg->flags_any(SegmentFlags::kShowerTrajectory))
        //         seg_type = "Shower_traj";
        //     else if (seg->flags_any(SegmentFlags::kShowerTopology))
        //         seg_type = "Shower_topo";
        //     else if (seg->has_particle_info() && std::abs(seg->particle_info()->pdg()) == 11)
        //         seg_type = "Electron";
        //     else
        //         seg_type = "Track";

        //     double length = segment_track_length(seg, 0);
        //     int    pdg    = seg->has_particle_info() ? seg->particle_info()->pdg()  : 0;
        //     double mass   = seg->has_particle_info() ? seg->particle_info()->mass() / units::MeV : 0.0;
        //     double ke     = seg->has_particle_info() ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0;
        //     double score  = seg->particle_score();
        //     SPDLOG_LOGGER_TRACE(s_log,
        //         "shower_determining_in_main_cluster: {} len={:.2f}cm dir={} weak={} pdg={} mass={:.2f}MeV KE={:.2f}MeV score={:.3f}",
        //         seg_type, length / units::cm,
        //         seg->dirsign(), seg_dir_weak(seg) ? 1 : 0,
        //         pdg, mass, ke, score);
        // }
    }
}
