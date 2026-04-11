// NeutrinoTaggerSSM.cxx — Phase A: includes + local helpers
// Prototype: prototype_pid/src/NeutrinoID_ssm_tagger.h
// Translation conventions documented in plan (see neutrino_id_function_map.md)

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <set>
#include <map>
#include <array>
#include <algorithm>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;
using namespace WireCell;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// ---------------------------------------------------------------------------
// Local helpers
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

static inline double safe_acos(double x) {
    return std::acos(std::min(1.0, std::max(-1.0, x)));
}

// get_scores: returns {mu_fwd, p_fwd, e_fwd, mu_bck, p_bck, e_bck}
// Prototype: NeutrinoID::get_scores(ProtoSegment*)
static std::array<double,6> get_scores(SegmentPtr sg,
                                        const ParticleDataSet::pointer& particle_data)
{
    auto& fits = sg->fits();
    int size = (int)fits.size();
    if (size == 0) return {0,0,0,0,0,0};

    std::vector<double> L(size,0), dQ_dx(size,0), rL(size,0), rdQ_dx(size,0);
    double dis = 0;
    for (int i = 0; i < size; ++i) {
        double dq_dx_val = fits[i].dQ / (fits[i].dx / units::cm + 1e-9);
        dQ_dx[i] = dq_dx_val;
        rdQ_dx[size-1-i] = dq_dx_val;
        L[i] = dis;
        if (i+1 < size) {
            auto d = fits[i+1].point - fits[i].point;
            dis += std::sqrt(d.dot(d));
        }
    }
    for (int i = 0; i < size; ++i) rL[i] = L.back() - L[size-1-i];

    auto rf = do_track_comp(L,   dQ_dx,  15*units::cm, 0.0, particle_data);
    auto rb = do_track_comp(rL,  rdQ_dx, 15*units::cm, 0.0, particle_data);
    return { rf.at(1), rf.at(2), rf.at(3), rb.at(1), rb.at(2), rb.at(3) };
}

// Overload with break_point/dir — returns raw do_track_comp result (caller uses indices 1,2,3)
// Prototype: NeutrinoID::get_scores(ProtoSegment*, int break_point, int dir)
static std::vector<double> get_scores_bp(SegmentPtr sg, int break_point, int dir,
                                          const ParticleDataSet::pointer& particle_data)
{
    auto& fits = sg->fits();
    int full = (int)fits.size();
    int size  = (dir == -1) ? break_point : (full - break_point);
    int start = (dir == -1) ? 0           : break_point;
    if (size <= 0) return {0,0,0,0};

    std::vector<double> L(size,0), dQ_dx(size,0), rL(size,0), rdQ_dx(size,0);
    double dis = 0;
    for (int i = start; i < start+size; ++i) {
        int idx = i - start;
        double dq_dx_val = fits[i].dQ / (fits[i].dx / units::cm + 1e-9);
        dQ_dx[idx] = dq_dx_val;
        rdQ_dx[size-1-idx] = dq_dx_val;
        L[idx] = dis;
        if (i+1 < start+size) {
            auto d = fits[i+1].point - fits[i].point;
            dis += std::sqrt(d.dot(d));
        }
    }
    for (int i = 0; i < size; ++i) rL[i] = L.back() - L[size-1-i];

    if (dir == -1)
        return do_track_comp(rL, rdQ_dx, 15*units::cm, 0.0, particle_data);
    else
        return do_track_comp(L,  dQ_dx,  15*units::cm, 0.0, particle_data);
}

// get_containing_shower_info: {shower_start_seg_graph_index, ke, flag_shower}; id=-1 if none
static std::tuple<int,double,int> get_containing_shower_info(
    SegmentPtr seg, const ShowerSegmentMap& map_segment_in_shower)
{
    auto it = map_segment_in_shower.find(seg);
    if (it == map_segment_in_shower.end()) return {-1,-1,-1};
    ShowerPtr sh = it->second;
    SegmentPtr ss = sh->start_segment();
    int id = ss ? (int)ss->get_graph_index() : -1;
    double ke = sh->get_kine_best();
    if (ke == 0.0) ke = sh->get_kine_charge();
    int flag = sh->get_flag_shower() ? 1 : 0;
    return {id, ke, flag};
}

// fill_ssmsp: writes space points of one real segment into TaggerInfo ssmsp_* vectors.
// Prototype: NeutrinoID::fill_ssmsp(ProtoSegment*, int pdg, int mother, int dir)
static void fill_ssmsp(SegmentPtr sg, int pdg, int mother, int dir,
                       TaggerInfo& ti,
                       const ShowerSegmentMap& map_segment_in_shower,
                       const ParticleDataSet::pointer& particle_data)
{
    auto& fits = sg->fits();
    int npts = (int)fits.size();
    double length_cm = segment_track_length(sg) / units::cm;
    double residual_range = length_cm;

    auto [csh_id, csh_ke, csh_flag] = get_containing_shower_info(sg, map_segment_in_shower);
    int seg_id = (int)sg->get_graph_index();

    ti.ssmsp_Nsp.push_back(0);
    for (int point = 0; point < npts; ++point) {
        int tp = (dir == -1) ? (npts-point-1) : point;

        ti.ssmsp_id.push_back(seg_id);
        ti.ssmsp_pdg.push_back(pdg);
        ti.ssmsp_mother.push_back(mother);
        ti.ssmsp_x.push_back(fits[tp].point.x() / units::cm);
        ti.ssmsp_y.push_back(fits[tp].point.y() / units::cm);
        ti.ssmsp_z.push_back(fits[tp].point.z() / units::cm);
        ti.ssmsp_dQ.push_back(fits[tp].dQ);
        ti.ssmsp_dx.push_back(fits[tp].dx / units::cm);
        ti.ssmsp_KE.push_back(cal_kine_range(residual_range, std::abs(pdg), particle_data));
        ti.ssmsp_containing_shower_id.push_back(csh_id);
        ti.ssmsp_containing_shower_ke.push_back(csh_ke);
        ti.ssmsp_containing_shower_flag.push_back(csh_flag);

        if (point < npts-1) {
            int tp2 = (dir == -1) ? (tp-1) : (tp+1);
            auto d = fits[tp2].point - fits[tp].point;
            residual_range -= std::sqrt(d.dot(d)) / units::cm;
            if (residual_range < 0.0) residual_range = 0.0;
        }
        ti.ssmsp_Nsp_tot += 1;
        ti.ssmsp_Nsp.at(ti.ssmsp_Ntrack) += 1;
    }
    ti.ssmsp_Ntrack += 1;
}

// fill_ssmsp_pseudo overload 1: vertex→shower_start (no sg)
// Prototype: fill_ssmsp_psuedo(WCShower*, int mother, int acc_id)
static int fill_ssmsp_pseudo_1(ShowerPtr shower, int mother, int acc_id, TaggerInfo& ti)
{
    int pdg = (std::abs(shower->get_particle_type())==11 || std::abs(shower->get_particle_type())==22) ? 22 : 2112;
    int id = -(acc_id+1);
    ti.ssmsp_Nsp.push_back(2);

    // pt1: start vertex position
    auto vtx = shower->get_start_vertex_and_type().first;
    Point vp = vtx ? vtx_fit_pt(vtx) : Point(0,0,0);
    auto push_pt = [&](Point p) {
        ti.ssmsp_id.push_back(id); ti.ssmsp_pdg.push_back(pdg); ti.ssmsp_mother.push_back(mother);
        ti.ssmsp_containing_shower_id.push_back(-1); ti.ssmsp_containing_shower_ke.push_back(-1); ti.ssmsp_containing_shower_flag.push_back(-1);
        ti.ssmsp_x.push_back(p.x()/units::cm); ti.ssmsp_y.push_back(p.y()/units::cm); ti.ssmsp_z.push_back(p.z()/units::cm);
        ti.ssmsp_dQ.push_back(0); ti.ssmsp_dx.push_back(0); ti.ssmsp_KE.push_back(0);
        ti.ssmsp_Nsp_tot++;
    };
    push_pt(vp);

    // pt2: start of shower start segment
    SegmentPtr ss = shower->start_segment();
    Point sp(0,0,0);
    if (ss && !ss->fits().empty())
        sp = (ss->dirsign()==1) ? ss->fits().front().point : ss->fits().back().point;
    push_pt(sp);

    ti.ssmsp_Ntrack++;
    return id;
}

// fill_ssmsp_pseudo overload 2: shower_start→parent_sg
// Prototype: fill_ssmsp_psuedo(WCShower*, ProtoSegment* sg, int mother, int acc_id)
static int fill_ssmsp_pseudo_2(ShowerPtr shower, SegmentPtr sg, int mother, int acc_id, TaggerInfo& ti)
{
    int pdg = 22;
    int id = -(acc_id+1);
    SegmentPtr ss = shower->start_segment();
    int csh_id = ss ? (int)ss->get_graph_index() : -1;
    double csh_ke = shower->get_kine_best(); if (csh_ke==0.0) csh_ke = shower->get_kine_charge();
    int csh_flag = shower->get_flag_shower() ? 1 : 0;

    ti.ssmsp_Nsp.push_back(2);

    auto push_pt = [&](Point p, int cid, double cke, int cfl) {
        ti.ssmsp_id.push_back(id); ti.ssmsp_pdg.push_back(pdg); ti.ssmsp_mother.push_back(mother);
        ti.ssmsp_containing_shower_id.push_back(cid); ti.ssmsp_containing_shower_ke.push_back(cke); ti.ssmsp_containing_shower_flag.push_back(cfl);
        ti.ssmsp_x.push_back(p.x()/units::cm); ti.ssmsp_y.push_back(p.y()/units::cm); ti.ssmsp_z.push_back(p.z()/units::cm);
        ti.ssmsp_dQ.push_back(0); ti.ssmsp_dx.push_back(0); ti.ssmsp_KE.push_back(0);
        ti.ssmsp_Nsp_tot++;
    };

    // pt1: start of shower start segment
    Point sp1(0,0,0);
    if (ss && !ss->fits().empty())
        sp1 = (ss->dirsign()==1) ? ss->fits().front().point : ss->fits().back().point;
    push_pt(sp1, csh_id, csh_ke, csh_flag);

    // pt2: start of parent segment sg
    Point sp2(0,0,0);
    if (sg && !sg->fits().empty())
        sp2 = (sg->dirsign()==1) ? sg->fits().front().point : sg->fits().back().point;
    push_pt(sp2, csh_id, csh_ke, csh_flag);

    ti.ssmsp_Ntrack++;
    return id;
}

// fill_ssmsp_pseudo overload 3: mother_sg→daughter shower start
// Prototype: fill_ssmsp_psuedo(WCShower*, ProtoSegment* mother_sg, int acc_id)
static int fill_ssmsp_pseudo_3(ShowerPtr shower, SegmentPtr mother_sg, int acc_id, TaggerInfo& ti)
{
    int pdg = (std::abs(shower->get_particle_type())==11 || std::abs(shower->get_particle_type())==22) ? 22 : 2112;
    int id = -(acc_id+1);
    int mother = mother_sg ? (int)mother_sg->get_graph_index() : 0;

    ti.ssmsp_Nsp.push_back(2);

    auto push_pt = [&](Point p) {
        ti.ssmsp_id.push_back(id); ti.ssmsp_pdg.push_back(pdg); ti.ssmsp_mother.push_back(mother);
        ti.ssmsp_containing_shower_id.push_back(-1); ti.ssmsp_containing_shower_ke.push_back(-1); ti.ssmsp_containing_shower_flag.push_back(-1);
        ti.ssmsp_x.push_back(p.x()/units::cm); ti.ssmsp_y.push_back(p.y()/units::cm); ti.ssmsp_z.push_back(p.z()/units::cm);
        ti.ssmsp_dQ.push_back(0); ti.ssmsp_dx.push_back(0); ti.ssmsp_KE.push_back(0);
        ti.ssmsp_Nsp_tot++;
    };

    // pt1: end of mother_sg (or shower start vertex if null)
    Point sp1(0,0,0);
    if (mother_sg && !mother_sg->fits().empty())
        sp1 = (mother_sg->dirsign()==1) ? mother_sg->fits().back().point : mother_sg->fits().front().point;
    else {
        auto vtx = shower->get_start_vertex_and_type().first;
        if (vtx) sp1 = vtx_fit_pt(vtx);
    }
    push_pt(sp1);

    // pt2: start of daughter shower's start segment
    SegmentPtr ss = shower->start_segment();
    Point sp2(0,0,0);
    if (ss && !ss->fits().empty())
        sp2 = (ss->dirsign()==1) ? ss->fits().front().point : ss->fits().back().point;
    push_pt(sp2);

    ti.ssmsp_Ntrack++;
    return id;
}

// find_incoming_segment: first segment at vtx already in used_segments (the "parent" segment)
static SegmentPtr find_incoming_segment(const Graph& graph, VertexPtr vtx,
                                         const std::set<SegmentPtr>& used_segments)
{
    if (!vtx || !vtx->descriptor_valid()) return nullptr;
    for (auto [eit,end] = boost::out_edges(vtx->get_descriptor(), graph); eit != end; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (sg && used_segments.count(sg)) return sg;
    }
    return nullptr;
}

// fill_ssmsp_all: BFS from start_vtx + shower handling — shared by success/exit paths.
// In success path: ssm_sg != nullptr, ssm_main_vtx is the "neutrino vertex" side.
// In exit path:    ssm_sg == nullptr, start from main_vertex.
static void fill_ssmsp_all(
    Graph& graph, VertexPtr main_vertex,
    VertexPtr ssm_main_vtx, VertexPtr ssm_second_vtx,
    SegmentPtr ssm_sg, int dir_ssm,
    IndexedShowerSet& showers,
    ShowerVertexMap& map_vertex_in_shower,
    ShowerSegmentMap& map_segment_in_shower,
    int& acc_segment_id,
    const ParticleDataSet::pointer& particle_data,
    TaggerInfo& ti)
{
    std::set<SegmentPtr> used_segments;
    std::set<VertexPtr>  used_vertices;
    std::vector<std::pair<VertexPtr,SegmentPtr>> to_examine;
    int temp_acc = acc_segment_id;

    VertexPtr ref_vtx = ssm_sg ? ssm_main_vtx : main_vertex;

    if (ssm_sg) {
        // success: SSM first, then all at ssm_main_vtx as protons
        fill_ssmsp(ssm_sg, 13, 0, dir_ssm, ti, map_segment_in_shower, particle_data);
        used_segments.insert(ssm_sg);
        to_examine.push_back({ssm_second_vtx, ssm_sg});

        for (auto [eit,end] = boost::out_edges(ssm_main_vtx->get_descriptor(), graph); eit != end; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg || sg == ssm_sg) continue;
            fill_ssmsp(sg, 2212, 0, sg->dirsign(), ti, map_segment_in_shower, particle_data);
            used_segments.insert(sg);
            if (auto ov = find_other_vertex(graph, sg, ssm_main_vtx))
                to_examine.push_back({ov, sg});
        }
        used_vertices.insert(ssm_main_vtx);
    } else {
        // exit: start from main_vertex with nominal PDG
        for (auto [eit,end] = boost::out_edges(main_vertex->get_descriptor(), graph); eit != end; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;
            int pdg = sg->has_particle_info() ? sg->particle_info()->pdg() : 2212;
            fill_ssmsp(sg, pdg, 0, sg->dirsign(), ti, map_segment_in_shower, particle_data);
            used_segments.insert(sg);
            if (auto ov = find_other_vertex(graph, sg, main_vertex))
                to_examine.push_back({ov, sg});
        }
        used_vertices.insert(main_vertex);
    }

    // BFS daughters
    while (!to_examine.empty()) {
        std::vector<std::pair<VertexPtr,SegmentPtr>> next;
        for (auto& [curr_vtx, prev_sg] : to_examine) {
            if (!curr_vtx || used_vertices.count(curr_vtx)) continue;
            for (auto [eit,end] = boost::out_edges(curr_vtx->get_descriptor(), graph); eit != end; ++eit) {
                SegmentPtr curr_sg = graph[*eit].segment;
                if (!curr_sg || used_segments.count(curr_sg)) continue;
                used_segments.insert(curr_sg);
                int pdg = curr_sg->has_particle_info() ? curr_sg->particle_info()->pdg() : 2212;
                int mother = (int)prev_sg->get_graph_index();
                fill_ssmsp(curr_sg, pdg, mother, curr_sg->dirsign(), ti, map_segment_in_shower, particle_data);
                if (auto ov = find_other_vertex(graph, curr_sg, curr_vtx))
                    if (!used_vertices.count(ov))
                        next.push_back({ov, curr_sg});
            }
            used_vertices.insert(curr_vtx);
        }
        to_examine = std::move(next);
    }

    // showers not yet covered
    for (ShowerPtr shower : showers) {
        SegmentPtr curr_sg = shower->start_segment();
        if (!curr_sg) continue;
        auto [pair_vtx, conn_type] = shower->get_start_vertex_and_type();
        int temp_mother = 0;

        if (conn_type == 1) {
            if (used_segments.count(curr_sg)) continue;
            if (pair_vtx == ref_vtx) { used_segments.insert(curr_sg); continue; }
            SegmentPtr prev_sg = nullptr;
            if (map_vertex_in_shower.count(pair_vtx))
                prev_sg = map_vertex_in_shower.at(pair_vtx)->start_segment();
            else
                prev_sg = find_incoming_segment(graph, pair_vtx, used_segments);
            if (!prev_sg) continue;
            temp_mother = (int)prev_sg->get_graph_index();
            int pdg = shower->get_particle_type();
            fill_ssmsp(curr_sg, pdg, temp_mother, curr_sg->dirsign(), ti, map_segment_in_shower, particle_data);
        } else if (conn_type == 2 || conn_type == 3) {
            if (pair_vtx == ref_vtx) {
                temp_mother = fill_ssmsp_pseudo_1(shower, 0, temp_acc++, ti);
            } else {
                SegmentPtr prev_sg = nullptr;
                if (map_vertex_in_shower.count(pair_vtx))
                    prev_sg = map_vertex_in_shower.at(pair_vtx)->start_segment();
                else
                    prev_sg = find_incoming_segment(graph, pair_vtx, used_segments);
                if (!prev_sg)
                    temp_mother = fill_ssmsp_pseudo_1(shower, 0, temp_acc++, ti);
                else {
                    // success path: fill_ssmsp_pseudo_3 variant
                    temp_mother = fill_ssmsp_pseudo_3(shower, prev_sg, temp_acc++, ti);
                }
            }
            int pdg = curr_sg->has_particle_info() ? curr_sg->particle_info()->pdg() : 11;
            fill_ssmsp(curr_sg, pdg, temp_mother, curr_sg->dirsign(), ti, map_segment_in_shower, particle_data);
        } else {
            continue;
        }
        used_segments.insert(curr_sg);

        // other segs inside shower
        IndexedVertexSet sv; IndexedSegmentSet ss;
        shower->fill_sets(sv, ss, false);
        for (SegmentPtr shower_sg : ss) {
            if (used_segments.count(shower_sg)) continue;
            int pseudo_id = fill_ssmsp_pseudo_2(shower, shower_sg, temp_mother, temp_acc++, ti);
            int pdg = shower_sg->has_particle_info() ? shower_sg->particle_info()->pdg() : 11;
            fill_ssmsp(shower_sg, pdg, pseudo_id, shower_sg->dirsign(), ti, map_segment_in_shower, particle_data);
            used_segments.insert(shower_sg);
        }
    }
    acc_segment_id = temp_acc;
}

// ---------------------------------------------------------------------------
// ParticleBlock: top-2 tracks and top-2 showers at a vertex, with all their
// properties.  Used by ssm_tagger Phase C to replace the duplicate prim/daught
// loops.  fill_particle_block_at_vtx populates one block from a vertex.
// ---------------------------------------------------------------------------

struct ParticleSlot {
    double length        = -999;
    double direct_length = -999;
    int    pdg           = 0;
    double medium_dq_dx  = -999;
    double max_dev       = -999;
    double ke_cal        = -999;
    double ke_rng        = -999;
    double ke_rng_mu     = -999;
    double ke_rng_p      = -999;
    double ke_rng_e      = -999;
    double ke_best       = -999;   // showers: from shower obj; tracks: same as ke_rng
    Vector dir;
    double score_mu_fwd  = -999, score_p_fwd  = -999, score_e_fwd  = -999;
    double score_mu_bck  = -999, score_p_bck  = -999, score_e_bck  = -999;
    double add_daught_track_1  = -999, add_daught_all_1  = -999;
    double add_daught_track_5  = -999, add_daught_all_5  = -999;
    double add_daught_track_11 = -999, add_daught_all_11 = -999;
    SegmentPtr sg = nullptr;
};

struct ParticleBlock {
    int n_tracks_1=0, n_tracks_3=0, n_tracks_5=0, n_tracks_8=0, n_tracks_11=0;
    int n_all_1=0,    n_all_3=0,    n_all_5=0,    n_all_8=0,    n_all_11=0;
    ParticleSlot track1, track2;
    ParticleSlot shw1,   shw2;
};

// Iterate segments at vtx, skipping skip_sg; rank top-2 tracks and top-2
// showers by length.  Mirrors the two identical prim/daught loops in the
// prototype (NeutrinoID_ssm_tagger.h lines 794-1224).
static ParticleBlock fill_particle_block_at_vtx(
    PatternAlgorithms& self,
    Graph& graph,
    VertexPtr vtx,
    SegmentPtr skip_sg,
    ShowerSegmentMap& map_segment_in_shower,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model)
{
    ParticleBlock pb;
    if (!vtx || !vtx->descriptor_valid()) return pb;

    // daughter_counts: mirrors ssm_tagger's inner lambda (prototype pattern).
    auto dc = [&](SegmentPtr sg, double len_cm, bool is_shower) {
        double cut = len_cm * units::cm;
        int dt = self.calculate_num_daughter_tracks(graph, vtx, sg, false, cut).first;
        int da = self.calculate_num_daughter_tracks(graph, vtx, sg, true,  cut).first;
        if (!is_shower) { dt -= 1; da -= 1; }
        else            { da -= 1; }
        if (segment_track_length(sg)/units::cm < len_cm) {
            if (!is_shower) dt++;
            da++;
        }
        if (dt < 0) dt = 0;
        if (da < 0) da = 0;
        return std::make_pair(dt, da);
    };

    double len_t1=-1e9, len_t2=-1e9, len_s1=-1e9, len_s2=-1e9;

    for (auto [eit,end] = boost::out_edges(vtx->get_descriptor(), graph); eit != end; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg == skip_sg) continue;

        double sg_len    = segment_track_length(sg) / units::cm;
        double sg_dlen   = segment_track_direct_length(sg) / units::cm;
        int    sg_pdg    = std::abs(sg->has_particle_info() ? sg->particle_info()->pdg() : 0);
        double sg_med    = segment_median_dQ_dx(sg) / (43e3/units::cm);
        double sg_mdev   = segment_track_max_deviation(sg) / units::cm;
        double sg_ke_cal = segment_cal_kine_dQdx(sg, recomb_model);
        double sg_ke_rng = cal_kine_range(sg_len, sg_pdg ? sg_pdg : 13, particle_data);
        double sg_ke_rng_mu = cal_kine_range(sg_len, 13,   particle_data);
        double sg_ke_rng_p  = cal_kine_range(sg_len, 2212, particle_data);
        double sg_ke_rng_e  = cal_kine_range(sg_len, 11,   particle_data);
        Vector sg_dir = segment_cal_dir_3vector(sg);

        if (sg_len > 1)  { pb.n_all_1++;  if (!seg_is_shower(sg)) pb.n_tracks_1++; }
        if (sg_len > 3)  { pb.n_all_3++;  if (!seg_is_shower(sg)) pb.n_tracks_3++; }
        if (sg_len > 5)  { pb.n_all_5++;  if (!seg_is_shower(sg)) pb.n_tracks_5++; }
        if (sg_len > 8)  { pb.n_all_8++;  if (!seg_is_shower(sg)) pb.n_tracks_8++; }
        if (sg_len > 11) { pb.n_all_11++; if (!seg_is_shower(sg)) pb.n_tracks_11++; }

        bool is_shower = seg_is_shower(sg);
        double& len1 = is_shower ? len_s1 : len_t1;
        double& len2 = is_shower ? len_s2 : len_t2;
        ParticleSlot& slot1 = is_shower ? pb.shw1 : pb.track1;
        ParticleSlot& slot2 = is_shower ? pb.shw2 : pb.track2;

        if (sg_len < len2) continue;

        auto [dtr1,dal1]   = dc(sg, 1.0,  is_shower);
        auto [dtr5,dal5]   = dc(sg, 5.0,  is_shower);
        auto [dtr11,dal11] = dc(sg, 11.0, is_shower);

        double ke_best = sg_ke_rng;
        if (is_shower) {
            auto it_sh = map_segment_in_shower.find(sg);
            if (it_sh != map_segment_in_shower.end()) {
                ke_best = it_sh->second->get_kine_best();
                if (ke_best == 0.0) ke_best = it_sh->second->get_kine_charge();
            }
        }

        ParticleSlot s;
        s.length=sg_len; s.direct_length=sg_dlen; s.pdg=sg_pdg;
        s.medium_dq_dx=sg_med; s.max_dev=sg_mdev;
        s.ke_cal=sg_ke_cal; s.ke_rng=sg_ke_rng;
        s.ke_rng_mu=sg_ke_rng_mu; s.ke_rng_p=sg_ke_rng_p; s.ke_rng_e=sg_ke_rng_e;
        s.ke_best=ke_best; s.dir=sg_dir; s.sg=sg;
        s.add_daught_track_1=dtr1; s.add_daught_all_1=dal1;
        s.add_daught_track_5=dtr5; s.add_daught_all_5=dal5;
        s.add_daught_track_11=dtr11; s.add_daught_all_11=dal11;

        if (sg_len < len1) {
            len2 = sg_len; slot2 = s;
        } else {
            len2 = len1; slot2 = slot1;
            len1 = sg_len; slot1 = s;
        }
    }

    // Compute PID scores for all filled slots.
    auto fill_scores = [&](ParticleSlot& ps) {
        if (!ps.sg) return;
        auto sc = get_scores(ps.sg, particle_data);
        ps.score_mu_fwd=sc[0]; ps.score_p_fwd=sc[1]; ps.score_e_fwd=sc[2];
        ps.score_mu_bck=sc[3]; ps.score_p_bck=sc[4]; ps.score_e_bck=sc[5];
        if (ps.sg->dirsign() == -1) {
            std::swap(ps.score_mu_fwd, ps.score_mu_bck);
            std::swap(ps.score_p_fwd,  ps.score_p_bck);
            std::swap(ps.score_e_fwd,  ps.score_e_bck);
        }
    };
    fill_scores(pb.track1); fill_scores(pb.track2);
    fill_scores(pb.shw1);   fill_scores(pb.shw2);

    return pb;
}

// ===========================================================================
// ssm_tagger
//
// Identifies Short Straight Muon (SSM / KDAR-like) events.
// Returns true if an SSM is found; false otherwise (all ssm_* fields set to -999).
//
// Prototype: WCPPID::NeutrinoID::ssm_tagger() in NeutrinoID_ssm_tagger.h
// exit_ssm_tagger() merged as early-return lambda.
// print_ssm_tagger() omitted (pure logging).
// ===========================================================================
bool PatternAlgorithms::ssm_tagger(
    Graph& graph,
    VertexPtr main_vertex,
    IndexedShowerSet& showers,
    ShowerVertexMap& map_vertex_in_shower,
    ShowerSegmentMap& map_segment_in_shower,
    const Pi0KineFeatures& pio_kine,
    int flag_ssmsp,
    int& acc_segment_id,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    TaggerInfo& ti)
{
    // --- initialise count fields (always written) ---
    ti.ssmsp_Ntrack = 0;
    ti.ssmsp_Nsp_tot = 0;

    // Fixed reference directions (beam = z, drift = x, vertical = y)
    const Vector dir_beam(0,0,1);
    const Vector dir_vertical(0,1,0);
    const Vector target_dir(0.46,0.05,0.885);
    const Vector absorber_dir(0.33,0.75,-0.59);

    int  Nsm = 0;
    int  Nsm_wivtx = 0;

    // ------------------------------------------------------------------
    // Phase A: identify SSM candidates at main_vertex
    // Prototype lines 413-536
    // ------------------------------------------------------------------
    // map: segment → {has_vtx_activity, is_backwards, length_cm}
    std::map<SegmentPtr, std::tuple<bool,bool,double>, SegmentIndexCmp> all_ssm_sg;

    if (!main_vertex || !main_vertex->descriptor_valid()) {
        // nothing to do — fall through to exit path
    } else {
        for (auto [eit,end] = boost::out_edges(main_vertex->get_descriptor(), graph); eit != end; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;

            double sg_length        = segment_track_length(sg) / units::cm;
            double sg_direct_length = segment_track_direct_length(sg) / units::cm;
            int    sg_pdg           = sg->has_particle_info() ? sg->particle_info()->pdg() : 0;

            bool is_muon = (std::abs(sg_pdg) == 13);
            bool is_elec = (std::abs(sg_pdg) == 11);

            // length / PDG pre-cuts (prototype line 422)
            if (!((sg_length <= 46 && sg_length >= 1 && is_muon) ||
                  (sg_direct_length > 0.9*sg_length && sg_length <= 44 && sg_length >= 1 && is_elec)))
                continue;
            if (sg_length < 1) continue;
            if (segment_median_dQ_dx(sg) / (43e3/units::cm) < 0.95) continue;

            Nsm++;

            // compute d(dQ/dx) vector for this segment
            auto& fits = sg->fits();
            int nfits = (int)fits.size();
            std::vector<double> vec_d_dqdx;
            if (nfits > 1) {
                double last = fits[0].dQ / (fits[0].dx/units::cm);
                for (int i = 1; i < nfits; ++i) {
                    double cur = fits[i].dQ / (fits[i].dx/units::cm);
                    vec_d_dqdx.push_back(cur - last);
                    last = cur;
                }
            }

            // check first 5 and last 5 d(dQ/dx) for vertex activity
            bool va_fwd = false, va_bck = false;
            double max_fwd5=0, max_bck5=0, max_fwd3=0, max_bck3=0;

            int end4 = std::min(4, (int)vec_d_dqdx.size());
            int end2 = std::min(2, (int)vec_d_dqdx.size());
            for (int i = 0; i < end4; ++i) {
                double av = std::abs(vec_d_dqdx[i]);
                if (av > 0.7 && av > max_fwd5) { va_fwd = true; max_fwd5 = av; }
                if (av > 0.7 && av > max_fwd3 && i < end2) max_fwd3 = av;
            }
            int st4 = std::max(0, (int)vec_d_dqdx.size()-4);
            int st2 = std::max(0, (int)vec_d_dqdx.size()-2);
            for (int i = st4; i < (int)vec_d_dqdx.size(); ++i) {
                double av = std::abs(vec_d_dqdx[i]);
                if (av > 0.7 && av > max_bck5) { va_bck = true; max_bck5 = av; }
                if (av > 0.7 && av > max_bck3 && i >= st2) max_bck3 = av;
            }

            // resolve ambiguity (prototype lines 476-518)
            if (va_fwd && va_bck) {
                if (max_fwd5 < 1.0) va_fwd = false;
                if (max_bck5 < 1.0) va_bck = false;
                if (max_fwd5 == max_bck5) {
                    if (max_fwd3 == max_bck3) { va_fwd = va_bck = false; }
                    else if (max_fwd3 > max_bck3) { va_bck = false; }
                    else { va_fwd = false; }
                } else if (!va_fwd && !va_bck) {
                    if (max_fwd5 >= max_bck5) va_fwd = true; else va_bck = true;
                } else if (va_fwd && va_bck) {
                    if (max_fwd5 < 1.3) va_fwd = false;
                    if (max_bck5 < 1.3) va_bck = false;
                    if (!va_fwd && !va_bck) {
                        if (max_fwd5 >= max_bck5) va_fwd = true; else va_bck = true;
                    } else if (va_fwd && va_bck) {
                        if (max_fwd5 >= max_bck5) va_bck = false; else va_fwd = false;
                    }
                }
            }

            bool sg_flag_vtx = false, sg_flag_back = false;
            if (va_bck) {
                Nsm_wivtx++;
                sg_flag_vtx = true;
                if (sg->dirsign() == 1) sg_flag_back = true;
            } else if (va_fwd) {
                Nsm_wivtx++;
                sg_flag_vtx = true;
                if (sg->dirsign() == -1) sg_flag_back = true;
            }
            all_ssm_sg[sg] = {sg_flag_vtx, sg_flag_back, sg_length};
        }
    }

    ti.ssm_Nsm       = (float)Nsm;
    ti.ssm_Nsm_wivtx = (float)Nsm_wivtx;

    // -----------------------------------------------------------------------
    // early-exit lambda (exit_ssm_tagger) — sets all ssm_* to -999
    // Also fills ssmsp if flag_ssmsp > 0.
    // -----------------------------------------------------------------------
    auto exit_ssm = [&]() -> bool {
        ti.ssm_flag_st_kdar = 0;
        ti.ssm_dq_dx_fwd_1 = ti.ssm_dq_dx_fwd_2 = ti.ssm_dq_dx_fwd_3 = ti.ssm_dq_dx_fwd_4 = ti.ssm_dq_dx_fwd_5 = -999;
        ti.ssm_dq_dx_bck_1 = ti.ssm_dq_dx_bck_2 = ti.ssm_dq_dx_bck_3 = ti.ssm_dq_dx_bck_4 = ti.ssm_dq_dx_bck_5 = -999;
        ti.ssm_d_dq_dx_fwd_12 = ti.ssm_d_dq_dx_fwd_23 = ti.ssm_d_dq_dx_fwd_34 = ti.ssm_d_dq_dx_fwd_45 = -999;
        ti.ssm_d_dq_dx_bck_12 = ti.ssm_d_dq_dx_bck_23 = ti.ssm_d_dq_dx_bck_34 = ti.ssm_d_dq_dx_bck_45 = -999;
        ti.ssm_max_dq_dx_fwd_3 = ti.ssm_max_dq_dx_fwd_5 = ti.ssm_max_dq_dx_bck_3 = ti.ssm_max_dq_dx_bck_5 = -999;
        ti.ssm_max_d_dq_dx_fwd_3 = ti.ssm_max_d_dq_dx_fwd_5 = ti.ssm_max_d_dq_dx_bck_3 = ti.ssm_max_d_dq_dx_bck_5 = -999;
        ti.ssm_medium_dq_dx = ti.ssm_medium_dq_dx_bp = -999;
        ti.ssm_angle_to_z = ti.ssm_angle_to_target = ti.ssm_angle_to_absorber = ti.ssm_angle_to_vertical = -999;
        ti.ssm_x_dir = ti.ssm_y_dir = ti.ssm_z_dir = -999;
        ti.ssm_kine_energy = ti.ssm_kine_energy_reduced = -999;
        ti.ssm_vtx_activity = ti.ssm_pdg = ti.ssm_dQ_dx_cut = -999;
        ti.ssm_score_mu_fwd = ti.ssm_score_p_fwd = ti.ssm_score_e_fwd = -999;
        ti.ssm_score_mu_bck = ti.ssm_score_p_bck = ti.ssm_score_e_bck = -999;
        ti.ssm_score_mu_fwd_bp = ti.ssm_score_p_fwd_bp = ti.ssm_score_e_fwd_bp = -999;
        ti.ssm_length = ti.ssm_direct_length = ti.ssm_length_ratio = ti.ssm_max_dev = -999;
        ti.ssm_n_prim_tracks_1 = ti.ssm_n_prim_tracks_3 = ti.ssm_n_prim_tracks_5 = -999;
        ti.ssm_n_prim_tracks_8 = ti.ssm_n_prim_tracks_11 = -999;
        ti.ssm_n_all_tracks_1 = ti.ssm_n_all_tracks_3 = ti.ssm_n_all_tracks_5 = -999;
        ti.ssm_n_all_tracks_8 = ti.ssm_n_all_tracks_11 = -999;
        ti.ssm_n_daughter_tracks_1 = ti.ssm_n_daughter_tracks_3 = ti.ssm_n_daughter_tracks_5 = -999;
        ti.ssm_n_daughter_tracks_8 = ti.ssm_n_daughter_tracks_11 = -999;
        ti.ssm_n_daughter_all_1 = ti.ssm_n_daughter_all_3 = ti.ssm_n_daughter_all_5 = -999;
        ti.ssm_n_daughter_all_8 = ti.ssm_n_daughter_all_11 = -999;
        // prim/daught track/shw 1&2 and offvtx track/shw1 set below
        // (long repetitive block mirroring exit_ssm_tagger lines 2786-3081)
        auto set_track_block = [](float& pdg_, float& smu_f, float& sp_f, float& se_f,
                                   float& smu_b, float& sp_b, float& se_b,
                                   float& len, float& dlen, float& lrat, float& mdev,
                                   float& ke_rng, float& ke_rng_mu, float& ke_rng_p, float& ke_rng_e,
                                   float& ke_cal, float& med_dq, float& xd, float& yd, float& zd,
                                   float& dtr1, float& dal1, float& dtr5, float& dal5,
                                   float& dtr11, float& dal11) {
            pdg_=-999; smu_f=-999; sp_f=-999; se_f=-999;
            smu_b=-999; sp_b=-999; se_b=-999;
            len=-999; dlen=-999; lrat=-999; mdev=-999;
            ke_rng=-999; ke_rng_mu=-999; ke_rng_p=-999; ke_rng_e=-999;
            ke_cal=-999; med_dq=-999; xd=-999; yd=-999; zd=-999;
            dtr1=-999; dal1=-999; dtr5=-999; dal5=-999; dtr11=-999; dal11=-999;
        };
        auto set_shw_block = [](float& pdg_, float& smu_f, float& sp_f, float& se_f,
                                 float& smu_b, float& sp_b, float& se_b,
                                 float& len, float& dlen, float& lrat, float& mdev,
                                 float& ke_best, float& ke_rng, float& ke_rng_mu, float& ke_rng_p, float& ke_rng_e,
                                 float& ke_cal, float& med_dq, float& xd, float& yd, float& zd,
                                 float& dtr1, float& dal1, float& dtr5, float& dal5,
                                 float& dtr11, float& dal11) {
            pdg_=-999; smu_f=-999; sp_f=-999; se_f=-999;
            smu_b=-999; sp_b=-999; se_b=-999;
            len=-999; dlen=-999; lrat=-999; mdev=-999;
            ke_best=-999; ke_rng=-999; ke_rng_mu=-999; ke_rng_p=-999; ke_rng_e=-999;
            ke_cal=-999; med_dq=-999; xd=-999; yd=-999; zd=-999;
            dtr1=-999; dal1=-999; dtr5=-999; dal5=-999; dtr11=-999; dal11=-999;
        };

        set_track_block(ti.ssm_prim_track1_pdg,
            ti.ssm_prim_track1_score_mu_fwd, ti.ssm_prim_track1_score_p_fwd, ti.ssm_prim_track1_score_e_fwd,
            ti.ssm_prim_track1_score_mu_bck, ti.ssm_prim_track1_score_p_bck, ti.ssm_prim_track1_score_e_bck,
            ti.ssm_prim_track1_length, ti.ssm_prim_track1_direct_length, ti.ssm_prim_track1_length_ratio, ti.ssm_prim_track1_max_dev,
            ti.ssm_prim_track1_kine_energy_range, ti.ssm_prim_track1_kine_energy_range_mu, ti.ssm_prim_track1_kine_energy_range_p, ti.ssm_prim_track1_kine_energy_range_e,
            ti.ssm_prim_track1_kine_energy_cal, ti.ssm_prim_track1_medium_dq_dx,
            ti.ssm_prim_track1_x_dir, ti.ssm_prim_track1_y_dir, ti.ssm_prim_track1_z_dir,
            ti.ssm_prim_track1_add_daught_track_counts_1, ti.ssm_prim_track1_add_daught_all_counts_1,
            ti.ssm_prim_track1_add_daught_track_counts_5, ti.ssm_prim_track1_add_daught_all_counts_5,
            ti.ssm_prim_track1_add_daught_track_counts_11, ti.ssm_prim_track1_add_daught_all_counts_11);
        set_track_block(ti.ssm_prim_track2_pdg,
            ti.ssm_prim_track2_score_mu_fwd, ti.ssm_prim_track2_score_p_fwd, ti.ssm_prim_track2_score_e_fwd,
            ti.ssm_prim_track2_score_mu_bck, ti.ssm_prim_track2_score_p_bck, ti.ssm_prim_track2_score_e_bck,
            ti.ssm_prim_track2_length, ti.ssm_prim_track2_direct_length, ti.ssm_prim_track2_length_ratio, ti.ssm_prim_track2_max_dev,
            ti.ssm_prim_track2_kine_energy_range, ti.ssm_prim_track2_kine_energy_range_mu, ti.ssm_prim_track2_kine_energy_range_p, ti.ssm_prim_track2_kine_energy_range_e,
            ti.ssm_prim_track2_kine_energy_cal, ti.ssm_prim_track2_medium_dq_dx,
            ti.ssm_prim_track2_x_dir, ti.ssm_prim_track2_y_dir, ti.ssm_prim_track2_z_dir,
            ti.ssm_prim_track2_add_daught_track_counts_1, ti.ssm_prim_track2_add_daught_all_counts_1,
            ti.ssm_prim_track2_add_daught_track_counts_5, ti.ssm_prim_track2_add_daught_all_counts_5,
            ti.ssm_prim_track2_add_daught_track_counts_11, ti.ssm_prim_track2_add_daught_all_counts_11);
        set_track_block(ti.ssm_daught_track1_pdg,
            ti.ssm_daught_track1_score_mu_fwd, ti.ssm_daught_track1_score_p_fwd, ti.ssm_daught_track1_score_e_fwd,
            ti.ssm_daught_track1_score_mu_bck, ti.ssm_daught_track1_score_p_bck, ti.ssm_daught_track1_score_e_bck,
            ti.ssm_daught_track1_length, ti.ssm_daught_track1_direct_length, ti.ssm_daught_track1_length_ratio, ti.ssm_daught_track1_max_dev,
            ti.ssm_daught_track1_kine_energy_range, ti.ssm_daught_track1_kine_energy_range_mu, ti.ssm_daught_track1_kine_energy_range_p, ti.ssm_daught_track1_kine_energy_range_e,
            ti.ssm_daught_track1_kine_energy_cal, ti.ssm_daught_track1_medium_dq_dx,
            ti.ssm_daught_track1_x_dir, ti.ssm_daught_track1_y_dir, ti.ssm_daught_track1_z_dir,
            ti.ssm_daught_track1_add_daught_track_counts_1, ti.ssm_daught_track1_add_daught_all_counts_1,
            ti.ssm_daught_track1_add_daught_track_counts_5, ti.ssm_daught_track1_add_daught_all_counts_5,
            ti.ssm_daught_track1_add_daught_track_counts_11, ti.ssm_daught_track1_add_daught_all_counts_11);
        set_track_block(ti.ssm_daught_track2_pdg,
            ti.ssm_daught_track2_score_mu_fwd, ti.ssm_daught_track2_score_p_fwd, ti.ssm_daught_track2_score_e_fwd,
            ti.ssm_daught_track2_score_mu_bck, ti.ssm_daught_track2_score_p_bck, ti.ssm_daught_track2_score_e_bck,
            ti.ssm_daught_track2_length, ti.ssm_daught_track2_direct_length, ti.ssm_daught_track2_length_ratio, ti.ssm_daught_track2_max_dev,
            ti.ssm_daught_track2_kine_energy_range, ti.ssm_daught_track2_kine_energy_range_mu, ti.ssm_daught_track2_kine_energy_range_p, ti.ssm_daught_track2_kine_energy_range_e,
            ti.ssm_daught_track2_kine_energy_cal, ti.ssm_daught_track2_medium_dq_dx,
            ti.ssm_daught_track2_x_dir, ti.ssm_daught_track2_y_dir, ti.ssm_daught_track2_z_dir,
            ti.ssm_daught_track2_add_daught_track_counts_1, ti.ssm_daught_track2_add_daught_all_counts_1,
            ti.ssm_daught_track2_add_daught_track_counts_5, ti.ssm_daught_track2_add_daught_all_counts_5,
            ti.ssm_daught_track2_add_daught_track_counts_11, ti.ssm_daught_track2_add_daught_all_counts_11);
        set_shw_block(ti.ssm_prim_shw1_pdg,
            ti.ssm_prim_shw1_score_mu_fwd, ti.ssm_prim_shw1_score_p_fwd, ti.ssm_prim_shw1_score_e_fwd,
            ti.ssm_prim_shw1_score_mu_bck, ti.ssm_prim_shw1_score_p_bck, ti.ssm_prim_shw1_score_e_bck,
            ti.ssm_prim_shw1_length, ti.ssm_prim_shw1_direct_length, ti.ssm_prim_shw1_length_ratio, ti.ssm_prim_shw1_max_dev,
            ti.ssm_prim_shw1_kine_energy_best,
            ti.ssm_prim_shw1_kine_energy_range, ti.ssm_prim_shw1_kine_energy_range_mu, ti.ssm_prim_shw1_kine_energy_range_p, ti.ssm_prim_shw1_kine_energy_range_e,
            ti.ssm_prim_shw1_kine_energy_cal, ti.ssm_prim_shw1_medium_dq_dx,
            ti.ssm_prim_shw1_x_dir, ti.ssm_prim_shw1_y_dir, ti.ssm_prim_shw1_z_dir,
            ti.ssm_prim_shw1_add_daught_track_counts_1, ti.ssm_prim_shw1_add_daught_all_counts_1,
            ti.ssm_prim_shw1_add_daught_track_counts_5, ti.ssm_prim_shw1_add_daught_all_counts_5,
            ti.ssm_prim_shw1_add_daught_track_counts_11, ti.ssm_prim_shw1_add_daught_all_counts_11);
        set_shw_block(ti.ssm_prim_shw2_pdg,
            ti.ssm_prim_shw2_score_mu_fwd, ti.ssm_prim_shw2_score_p_fwd, ti.ssm_prim_shw2_score_e_fwd,
            ti.ssm_prim_shw2_score_mu_bck, ti.ssm_prim_shw2_score_p_bck, ti.ssm_prim_shw2_score_e_bck,
            ti.ssm_prim_shw2_length, ti.ssm_prim_shw2_direct_length, ti.ssm_prim_shw2_length_ratio, ti.ssm_prim_shw2_max_dev,
            ti.ssm_prim_shw2_kine_energy_best,
            ti.ssm_prim_shw2_kine_energy_range, ti.ssm_prim_shw2_kine_energy_range_mu, ti.ssm_prim_shw2_kine_energy_range_p, ti.ssm_prim_shw2_kine_energy_range_e,
            ti.ssm_prim_shw2_kine_energy_cal, ti.ssm_prim_shw2_medium_dq_dx,
            ti.ssm_prim_shw2_x_dir, ti.ssm_prim_shw2_y_dir, ti.ssm_prim_shw2_z_dir,
            ti.ssm_prim_shw2_add_daught_track_counts_1, ti.ssm_prim_shw2_add_daught_all_counts_1,
            ti.ssm_prim_shw2_add_daught_track_counts_5, ti.ssm_prim_shw2_add_daught_all_counts_5,
            ti.ssm_prim_shw2_add_daught_track_counts_11, ti.ssm_prim_shw2_add_daught_all_counts_11);
        set_shw_block(ti.ssm_daught_shw1_pdg,
            ti.ssm_daught_shw1_score_mu_fwd, ti.ssm_daught_shw1_score_p_fwd, ti.ssm_daught_shw1_score_e_fwd,
            ti.ssm_daught_shw1_score_mu_bck, ti.ssm_daught_shw1_score_p_bck, ti.ssm_daught_shw1_score_e_bck,
            ti.ssm_daught_shw1_length, ti.ssm_daught_shw1_direct_length, ti.ssm_daught_shw1_length_ratio, ti.ssm_daught_shw1_max_dev,
            ti.ssm_daught_shw1_kine_energy_best,
            ti.ssm_daught_shw1_kine_energy_range, ti.ssm_daught_shw1_kine_energy_range_mu, ti.ssm_daught_shw1_kine_energy_range_p, ti.ssm_daught_shw1_kine_energy_range_e,
            ti.ssm_daught_shw1_kine_energy_cal, ti.ssm_daught_shw1_medium_dq_dx,
            ti.ssm_daught_shw1_x_dir, ti.ssm_daught_shw1_y_dir, ti.ssm_daught_shw1_z_dir,
            ti.ssm_daught_shw1_add_daught_track_counts_1, ti.ssm_daught_shw1_add_daught_all_counts_1,
            ti.ssm_daught_shw1_add_daught_track_counts_5, ti.ssm_daught_shw1_add_daught_all_counts_5,
            ti.ssm_daught_shw1_add_daught_track_counts_11, ti.ssm_daught_shw1_add_daught_all_counts_11);
        set_shw_block(ti.ssm_daught_shw2_pdg,
            ti.ssm_daught_shw2_score_mu_fwd, ti.ssm_daught_shw2_score_p_fwd, ti.ssm_daught_shw2_score_e_fwd,
            ti.ssm_daught_shw2_score_mu_bck, ti.ssm_daught_shw2_score_p_bck, ti.ssm_daught_shw2_score_e_bck,
            ti.ssm_daught_shw2_length, ti.ssm_daught_shw2_direct_length, ti.ssm_daught_shw2_length_ratio, ti.ssm_daught_shw2_max_dev,
            ti.ssm_daught_shw2_kine_energy_best,
            ti.ssm_daught_shw2_kine_energy_range, ti.ssm_daught_shw2_kine_energy_range_mu, ti.ssm_daught_shw2_kine_energy_range_p, ti.ssm_daught_shw2_kine_energy_range_e,
            ti.ssm_daught_shw2_kine_energy_cal, ti.ssm_daught_shw2_medium_dq_dx,
            ti.ssm_daught_shw2_x_dir, ti.ssm_daught_shw2_y_dir, ti.ssm_daught_shw2_z_dir,
            ti.ssm_daught_shw2_add_daught_track_counts_1, ti.ssm_daught_shw2_add_daught_all_counts_1,
            ti.ssm_daught_shw2_add_daught_track_counts_5, ti.ssm_daught_shw2_add_daught_all_counts_5,
            ti.ssm_daught_shw2_add_daught_track_counts_11, ti.ssm_daught_shw2_add_daught_all_counts_11);

        // off-vertex
        ti.ssm_offvtx_length = ti.ssm_offvtx_energy = -999;
        ti.ssm_n_offvtx_tracks_1 = ti.ssm_n_offvtx_tracks_3 = ti.ssm_n_offvtx_tracks_5 = -999;
        ti.ssm_n_offvtx_tracks_8 = ti.ssm_n_offvtx_tracks_11 = -999;
        ti.ssm_n_offvtx_showers_1 = ti.ssm_n_offvtx_showers_3 = ti.ssm_n_offvtx_showers_5 = -999;
        ti.ssm_n_offvtx_showers_8 = ti.ssm_n_offvtx_showers_11 = -999;
        ti.ssm_offvtx_track1_pdg = ti.ssm_offvtx_track1_score_mu_fwd = ti.ssm_offvtx_track1_score_p_fwd = -999;
        ti.ssm_offvtx_track1_score_e_fwd = ti.ssm_offvtx_track1_score_mu_bck = ti.ssm_offvtx_track1_score_p_bck = -999;
        ti.ssm_offvtx_track1_score_e_bck = ti.ssm_offvtx_track1_length = ti.ssm_offvtx_track1_direct_length = -999;
        ti.ssm_offvtx_track1_max_dev = ti.ssm_offvtx_track1_kine_energy_range = -999;
        ti.ssm_offvtx_track1_kine_energy_range_mu = ti.ssm_offvtx_track1_kine_energy_range_p = ti.ssm_offvtx_track1_kine_energy_range_e = -999;
        ti.ssm_offvtx_track1_kine_energy_cal = ti.ssm_offvtx_track1_medium_dq_dx = -999;
        ti.ssm_offvtx_track1_x_dir = ti.ssm_offvtx_track1_y_dir = ti.ssm_offvtx_track1_z_dir = ti.ssm_offvtx_track1_dist_mainvtx = -999;
        ti.ssm_offvtx_shw1_pdg_offvtx = ti.ssm_offvtx_shw1_score_mu_fwd = ti.ssm_offvtx_shw1_score_p_fwd = -999;
        ti.ssm_offvtx_shw1_score_e_fwd = ti.ssm_offvtx_shw1_score_mu_bck = ti.ssm_offvtx_shw1_score_p_bck = -999;
        ti.ssm_offvtx_shw1_score_e_bck = ti.ssm_offvtx_shw1_length = ti.ssm_offvtx_shw1_direct_length = -999;
        ti.ssm_offvtx_shw1_max_dev = ti.ssm_offvtx_shw1_kine_energy_best = ti.ssm_offvtx_shw1_kine_energy_range = -999;
        ti.ssm_offvtx_shw1_kine_energy_range_mu = ti.ssm_offvtx_shw1_kine_energy_range_p = ti.ssm_offvtx_shw1_kine_energy_range_e = -999;
        ti.ssm_offvtx_shw1_kine_energy_cal = ti.ssm_offvtx_shw1_medium_dq_dx = -999;
        ti.ssm_offvtx_shw1_x_dir = ti.ssm_offvtx_shw1_y_dir = ti.ssm_offvtx_shw1_z_dir = ti.ssm_offvtx_shw1_dist_mainvtx = -999;

        // event-level angles
        ti.ssm_nu_angle_z = ti.ssm_nu_angle_target = ti.ssm_nu_angle_absorber = ti.ssm_nu_angle_vertical = -999;
        ti.ssm_con_nu_angle_z = ti.ssm_con_nu_angle_target = ti.ssm_con_nu_angle_absorber = ti.ssm_con_nu_angle_vertical = -999;
        ti.ssm_prim_nu_angle_z = ti.ssm_prim_nu_angle_target = ti.ssm_prim_nu_angle_absorber = ti.ssm_prim_nu_angle_vertical = -999;
        ti.ssm_track_angle_z = ti.ssm_track_angle_target = ti.ssm_track_angle_absorber = ti.ssm_track_angle_vertical = -999;
        ti.ssm_vtxX = ti.ssm_vtxY = ti.ssm_vtxZ = -999;

        if (flag_ssmsp > 0 && main_vertex && main_vertex->descriptor_valid()) {
            fill_ssmsp_all(graph, main_vertex,
                           main_vertex, nullptr, nullptr, 0,
                           showers, map_vertex_in_shower, map_segment_in_shower,
                           acc_segment_id, particle_data, ti);
        }
        return false;
    };

    if (Nsm == 0) return exit_ssm();

    // ------------------------------------------------------------------
    // Phase B: select best SSM, compute properties
    // Prototype lines 545-764
    // ------------------------------------------------------------------
    SegmentPtr ssm_sg = nullptr;
    double length = 0, direct_length = 0;
    bool backwards_muon = false;
    bool vtx_activity = false;

    for (auto& [sg, info] : all_ssm_sg) {
        bool has_vtx = std::get<0>(info);
        bool is_back = std::get<1>(info);
        double sg_len = std::get<2>(info);
        if (has_vtx && (sg_len > length || !vtx_activity)) {
            ssm_sg = sg; length = sg_len; backwards_muon = is_back; vtx_activity = has_vtx;
        } else if (sg_len > length && !vtx_activity) {
            ssm_sg = sg; length = sg_len; backwards_muon = is_back; vtx_activity = has_vtx;
        }
    }
    if (!ssm_sg) return exit_ssm();

    int pdg = std::abs(ssm_sg->has_particle_info() ? ssm_sg->particle_info()->pdg() : 0);
    direct_length = segment_track_direct_length(ssm_sg) / units::cm;

    int dir = ssm_sg->dirsign();
    if (backwards_muon) dir = -dir;

    // direction vectors
    Vector init_dir_10 = segment_cal_dir_3vector(ssm_sg, dir, 10, 0).norm();

    double x_dir = init_dir_10.x(), y_dir = init_dir_10.y(), z_dir = init_dir_10.z();

    double max_dev = segment_track_max_deviation(ssm_sg) / units::cm;

    // angles (using 10cm-direction for ti, consistent with prototype ssm_angle_to_z = angle_to_z_10)
    double angle_to_z_10        = safe_acos(init_dir_10.dot(dir_beam));
    double angle_to_target_10   = safe_acos(init_dir_10.dot(target_dir));
    double angle_to_absorber_10 = safe_acos(init_dir_10.dot(absorber_dir));
    double angle_to_vertical_10 = safe_acos(init_dir_10.dot(dir_vertical));

    // dQ/dx profile (prototype lines 603-718)
    auto& fits_ssm = ssm_sg->fits();
    int nfits_ssm = (int)fits_ssm.size();

    std::vector<double> vec_dqdx, vec_d_dqdx, vec_abs_d_dqdx;
    if (nfits_ssm > 0) {
        double last = fits_ssm[0].dQ / (fits_ssm[0].dx/units::cm);
        vec_dqdx.push_back(last);
        for (int i = 1; i < nfits_ssm; ++i) {
            double cur = fits_ssm[i].dQ / (fits_ssm[i].dx/units::cm);
            vec_dqdx.push_back(cur);
            vec_d_dqdx.push_back(cur - last);
            vec_abs_d_dqdx.push_back(std::abs(cur - last));
            last = cur;
        }
    }

    // break_point detection (prototype lines 621-662)
    int break_point_fwd = -1, break_point_bck = -1, break_point = 0;
    double reduced_muon_length = length;

    {
        int end4 = std::min(4, (int)vec_d_dqdx.size());
        for (int i = 0; i < end4; ++i)
            if (std::abs(vec_d_dqdx[i]) > 0.7) { break_point_fwd = i+1; }
    }
    if (dir == 1 && break_point_fwd >= 0) {
        for (int i = 0; i < break_point_fwd; ++i)
            reduced_muon_length -= std::abs(fits_ssm[i].dx/units::cm);
        break_point = break_point_fwd;
    }

    {
        int st4 = std::max(0, (int)vec_d_dqdx.size()-4);
        for (int i = st4; i < (int)vec_d_dqdx.size(); ++i)
            if (std::abs(vec_d_dqdx[i]) > 0.7) { if (break_point_bck<0) break_point_bck = i+1; }
    }
    if (dir == -1 && break_point_bck >= 0) {
        for (int i = break_point_bck; i < nfits_ssm; ++i)
            reduced_muon_length -= std::abs(fits_ssm[i].dx/units::cm);
        break_point = break_point_bck;
    }

    // point-by-point dQ/dx values (fwd and bck, prototype lines 664-717)
    auto get_dqdx = [&](int idx) -> double {
        if (idx < 0 || idx >= (int)vec_dqdx.size()) return 0.0;
        return vec_dqdx[idx];
    };
    auto get_d_dqdx = [&](int idx) -> double {
        if (idx < 0 || idx >= (int)vec_d_dqdx.size()) return 0.0;
        return vec_d_dqdx[idx];
    };

    double dq_dx_fwd_1 = get_dqdx(0), dq_dx_fwd_2 = get_dqdx(1);
    double dq_dx_fwd_3 = get_dqdx(2), dq_dx_fwd_4 = get_dqdx(3), dq_dx_fwd_5 = get_dqdx(4);
    int n = (int)vec_dqdx.size();
    double dq_dx_bck_1 = get_dqdx(n-1), dq_dx_bck_2 = get_dqdx(n-2);
    double dq_dx_bck_3 = get_dqdx(n-3), dq_dx_bck_4 = get_dqdx(n-4), dq_dx_bck_5 = get_dqdx(n-5);

    double d_dq_dx_fwd_12 = get_d_dqdx(0), d_dq_dx_fwd_23 = get_d_dqdx(1);
    double d_dq_dx_fwd_34 = get_d_dqdx(2), d_dq_dx_fwd_45 = get_d_dqdx(3);
    int nd = (int)vec_d_dqdx.size();
    double d_dq_dx_bck_12 = get_d_dqdx(nd-1), d_dq_dx_bck_23 = get_d_dqdx(nd-2);
    double d_dq_dx_bck_34 = get_d_dqdx(nd-3), d_dq_dx_bck_45 = get_d_dqdx(nd-4);

    // max dQ/dx and d(dQ/dx) in first/last 3 and 5 points
    auto max_in = [&](const std::vector<double>& v, int from, int to) -> double {
        if (v.empty() || from >= (int)v.size()) return 0.0;
        to = std::min(to, (int)v.size());
        return *std::max_element(v.begin()+from, v.begin()+to);
    };
    double max_dqdx_fwd3  = max_in(vec_dqdx, 0, 3);
    double max_dqdx_fwd5  = max_in(vec_dqdx, 0, 5);
    double max_dqdx_bck3  = max_in(vec_dqdx, std::max(0,n-3), n);
    double max_dqdx_bck5  = max_in(vec_dqdx, std::max(0,n-5), n);
    double max_d_fwd3 = max_in(vec_abs_d_dqdx, 0, 2);
    double max_d_fwd5 = max_in(vec_abs_d_dqdx, 0, 4);
    double max_d_bck3 = max_in(vec_abs_d_dqdx, std::max(0,nd-2), nd);
    double max_d_bck5 = max_in(vec_abs_d_dqdx, std::max(0,nd-4), nd);

    // swap fwd/bck if dir==-1 (prototype lines 704-718)
    if (dir == -1) {
        std::swap(dq_dx_fwd_1, dq_dx_bck_1); std::swap(dq_dx_fwd_2, dq_dx_bck_2);
        std::swap(dq_dx_fwd_3, dq_dx_bck_3); std::swap(dq_dx_fwd_4, dq_dx_bck_4);
        std::swap(dq_dx_fwd_5, dq_dx_bck_5);
        std::swap(d_dq_dx_fwd_12, d_dq_dx_bck_12); std::swap(d_dq_dx_fwd_23, d_dq_dx_bck_23);
        std::swap(d_dq_dx_fwd_34, d_dq_dx_bck_34); std::swap(d_dq_dx_fwd_45, d_dq_dx_bck_45);
        std::swap(max_dqdx_fwd3, max_dqdx_bck3); std::swap(max_dqdx_fwd5, max_dqdx_bck5);
        std::swap(max_d_fwd3, max_d_bck3); std::swap(max_d_fwd5, max_d_bck5);
    }

    // PID scores (prototype lines 720-756)
    auto sc = get_scores(ssm_sg, particle_data);
    double score_mu_fwd = sc[0], score_p_fwd = sc[1], score_e_fwd = sc[2];
    double score_mu_bck = sc[3], score_p_bck = sc[4], score_e_bck = sc[5];
    if (dir == -1) {
        std::swap(score_mu_fwd, score_mu_bck);
        std::swap(score_p_fwd,  score_p_bck);
        std::swap(score_e_fwd,  score_e_bck);
    }

    double score_mu_fwd_bp = score_mu_fwd, score_p_fwd_bp = score_p_fwd, score_e_fwd_bp = score_e_fwd;
    if (vtx_activity) {
        auto bp_sc = get_scores_bp(ssm_sg, break_point, dir, particle_data);
        score_mu_fwd_bp = bp_sc.at(1);
        score_p_fwd_bp  = bp_sc.at(2);
        score_e_fwd_bp  = bp_sc.at(3);
    }

    // catch degenerate break_point
    if ((break_point == nd && dir == 1 && vtx_activity) ||
        (break_point == 0  && dir == -1 && vtx_activity)) {
        reduced_muon_length = length;
        score_mu_fwd_bp = score_mu_fwd;
        score_p_fwd_bp  = score_p_fwd;
        score_e_fwd_bp  = score_e_fwd;
    }

    double dQ_dx_cut   = 0.8866 + 0.9533 * std::pow(18.0/length, 0.4234);
    double medium_dq_dx = segment_median_dQ_dx(ssm_sg) / (43e3/units::cm);
    double medium_dq_dx_bp = medium_dq_dx;
    if (dir == 1 && vtx_activity && break_point >= 0)
        medium_dq_dx_bp = segment_median_dQ_dx(ssm_sg, break_point, nfits_ssm) / (43e3/units::cm);
    else if (vtx_activity && break_point >= 0)
        medium_dq_dx_bp = segment_median_dQ_dx(ssm_sg, 0, break_point) / (43e3/units::cm);

    // Kinetic energy from range (prototype line 762: g_range->Eval(length) * MeV)
    double kine_energy         = cal_kine_range(length, 13, particle_data);
    double kine_energy_reduced = cal_kine_range(reduced_muon_length, 13, particle_data);

    // ------------------------------------------------------------------
    // Phase C: primary/daughter particle loops
    // Prototype lines 766-1224
    // Two structurally identical loops (one per vertex) are replaced by two
    // calls to fill_particle_block_at_vtx.  All per-particle data live in
    // pb_prim / pb_daught; Phase D swaps them via std::swap if backwards_muon.
    // ------------------------------------------------------------------
    auto pb_prim   = fill_particle_block_at_vtx(*this, graph, main_vertex, ssm_sg,
                         map_segment_in_shower, particle_data, recomb_model);
    VertexPtr second_vtx = find_other_vertex(graph, ssm_sg, main_vertex);
    auto pb_daught = fill_particle_block_at_vtx(*this, graph, second_vtx, ssm_sg,
                         map_segment_in_shower, particle_data, recomb_model);

    // PID scores for secondary particles (prototype lines 1226-1338)
    // Scores are now computed inside fill_particle_block_at_vtx and stored in
    // each ParticleSlot (score_mu_fwd / score_p_fwd / score_e_fwd / _bck).

    // ------------------------------------------------------------------
    // Phase D: backwards_muon swap + momentum + off-vertex loop
    // Prototype lines 1340-1749
    // ------------------------------------------------------------------

    // Swap prim<->daughter if muon is backwards (prototype lines 1340-1464).
    // ParticleBlock holds all fields, so a single std::swap replaces ~115 SWAP_PAIR lines.
    if (backwards_muon) std::swap(pb_prim, pb_daught);

    // momentum vectors and SSM vertex assignment (prototype lines 1466-1562)
    // mass from particle_data; momentum magnitude = sqrt(KE^2 + 2*KE*m)
    auto get_mass = [&](int abs_pdg) -> double {
        return particle_data->get_particle_mass(abs_pdg ? abs_pdg : 13);
    };
    auto make_mom = [&](Vector dir_v, double ke, int abs_pdg) -> Vector {
        double m = get_mass(abs_pdg); if (m < 0) m = 0;
        double pmag = std::sqrt(ke*ke + 2.0*ke*m);
        return dir_v * pmag;
    };

    Vector mom     = make_mom(init_dir_10, kine_energy, 13);
    Vector mom_prim_track1  = (pb_prim.track1.length > 0 && pb_prim.track1.sg)
        ? make_mom(pb_prim.track1.dir,   pb_prim.track1.ke_rng,  pb_prim.track1.pdg)  : Vector(0,0,0);
    Vector mom_prim_track2  = (pb_prim.track2.length > 0 && pb_prim.track2.sg)
        ? make_mom(pb_prim.track2.dir,   pb_prim.track2.ke_rng,  pb_prim.track2.pdg)  : Vector(0,0,0);
    Vector mom_daught_track1= (pb_daught.track1.length > 0 && pb_daught.track1.sg)
        ? make_mom(pb_daught.track1.dir, pb_daught.track1.ke_rng,pb_daught.track1.pdg): Vector(0,0,0);
    Vector mom_daught_track2= (pb_daught.track2.length > 0 && pb_daught.track2.sg)
        ? make_mom(pb_daught.track2.dir, pb_daught.track2.ke_rng,pb_daught.track2.pdg): Vector(0,0,0);
    Vector mom_prim_shw1    = (pb_prim.shw1.length > 0 && pb_prim.shw1.sg)
        ? make_mom(pb_prim.shw1.dir,     pb_prim.shw1.ke_rng,    pb_prim.shw1.pdg)    : Vector(0,0,0);
    Vector mom_prim_shw2    = (pb_prim.shw2.length > 0 && pb_prim.shw2.sg)
        ? make_mom(pb_prim.shw2.dir,     pb_prim.shw2.ke_rng,    pb_prim.shw2.pdg)    : Vector(0,0,0);
    Vector mom_daught_shw1  = (pb_daught.shw1.length > 0 && pb_daught.shw1.sg)
        ? make_mom(pb_daught.shw1.dir,   pb_daught.shw1.ke_best, pb_daught.shw1.pdg)  : Vector(0,0,0);
    Vector mom_daught_shw2  = (pb_daught.shw2.length > 0 && pb_daught.shw2.sg)
        ? make_mom(pb_daught.shw2.dir,   pb_daught.shw2.ke_best, pb_daught.shw2.pdg)  : Vector(0,0,0);

    // ssm vertex assignment (prototype line 1558-1562)
    VertexPtr ssm_main_vtx   = main_vertex;
    VertexPtr ssm_second_vtx = find_other_vertex(graph, ssm_sg, main_vertex);
    if (backwards_muon) std::swap(ssm_main_vtx, ssm_second_vtx);

    // neutrino direction angles (prototype lines 1706-1748)
    // also pi0 momentum (prototype lines 1707-1715)
    Vector mom_pi0(0,0,0);
    if (pio_kine.mass > 0) {
        double e1 = pio_kine.energy_1, t1 = pio_kine.theta_1/180.*M_PI, p1 = pio_kine.phi_1/180.*M_PI;
        double e2 = pio_kine.energy_2, t2 = pio_kine.theta_2/180.*M_PI, p2 = pio_kine.phi_2/180.*M_PI;
        Vector pv1(e1*std::sin(t1)*std::cos(p1), e1*std::sin(t1)*std::sin(p1), e1*std::cos(t1));
        Vector pv2(e2*std::sin(t2)*std::cos(p2), e2*std::sin(t2)*std::sin(p2), e2*std::cos(t2));
        mom_pi0 = pv1 + pv2;
    }

    auto compute_angles = [&](Vector v, double& az, double& at, double& aa, double& av) {
        double mag = v.magnitude();
        if (mag < 1e-9) { az=at=aa=av=-999; return; }
        Vector vn = v / mag;
        az = safe_acos(vn.dot(dir_beam));
        at = safe_acos(vn.dot(target_dir));
        aa = safe_acos(vn.dot(absorber_dir));
        av = safe_acos(vn.dot(dir_vertical));
    };

    double nu_angle_z, nu_angle_target, nu_angle_absorber, nu_angle_vertical;
    double con_nu_angle_z, con_nu_angle_target, con_nu_angle_absorber, con_nu_angle_vertical;
    double prim_nu_angle_z, prim_nu_angle_target, prim_nu_angle_absorber, prim_nu_angle_vertical;
    double track_angle_z, track_angle_target, track_angle_absorber, track_angle_vertical;

    Vector nu_all = mom + mom_prim_track1 + mom_prim_track2 + mom_prim_shw1 + mom_prim_shw2
                  + mom_daught_track1 + mom_daught_track2 + mom_daught_shw1 + mom_daught_shw2
                  + mom_pi0;  // off-vertex mom added later
    Vector con_nu = mom + mom_prim_track1 + mom_prim_track2 + mom_prim_shw1 + mom_prim_shw2
                  + mom_daught_track1 + mom_daught_track2 + mom_daught_shw1 + mom_daught_shw2;
    Vector prim_nu = mom + mom_prim_track1 + mom_prim_track2 + mom_prim_shw1 + mom_prim_shw2;
    Vector track_nu= mom + mom_prim_track1 + mom_prim_track2;

    // off-vertex loop (prototype lines 1564-1705)
    double off_vtx_length=0, off_vtx_energy=0;
    double n_offvtx_tracks_1=0,n_offvtx_tracks_3=0,n_offvtx_tracks_5=0,n_offvtx_tracks_8=0,n_offvtx_tracks_11=0;
    double n_offvtx_showers_1=0,n_offvtx_showers_3=0,n_offvtx_showers_5=0,n_offvtx_showers_8=0,n_offvtx_showers_11=0;

    double pdg_offvtx_track1=-999, score_mu_fwd_offvtx_track1=-999, score_p_fwd_offvtx_track1=-999, score_e_fwd_offvtx_track1=-999;
    double score_mu_bck_offvtx_track1=-999, score_p_bck_offvtx_track1=-999, score_e_bck_offvtx_track1=-999;
    double length_offvtx_track1=-999, direct_length_offvtx_track1=-999, max_dev_offvtx_track1=-999;
    double kine_energy_cal_offvtx_track1=-999, kine_energy_range_offvtx_track1=-999;
    double kine_energy_range_mu_offvtx_track1=-999, kine_energy_range_p_offvtx_track1=-999, kine_energy_range_e_offvtx_track1=-999;
    double medium_dq_dx_offvtx_track1=-999;
    double x_dir_offvtx_track1=-999, y_dir_offvtx_track1=-999, z_dir_offvtx_track1=-999;
    double dist_mainvtx_offvtx_track1=-999;
    double len_ovt1 = -1e9;

    double pdg_offvtx_shw1=-999, score_mu_fwd_offvtx_shw1=-999, score_p_fwd_offvtx_shw1=-999, score_e_fwd_offvtx_shw1=-999;
    double score_mu_bck_offvtx_shw1=-999, score_p_bck_offvtx_shw1=-999, score_e_bck_offvtx_shw1=-999;
    double length_offvtx_shw1=-999, direct_length_offvtx_shw1=-999, max_dev_offvtx_shw1=-999;
    double kine_energy_cal_offvtx_shw1=-999, kine_energy_range_offvtx_shw1=-999, kine_energy_best_offvtx_shw1=-999;
    double kine_energy_range_mu_offvtx_shw1=-999, kine_energy_range_p_offvtx_shw1=-999, kine_energy_range_e_offvtx_shw1=-999;
    double medium_dq_dx_offvtx_shw1=-999;
    double x_dir_offvtx_shw1=-999, y_dir_offvtx_shw1=-999, z_dir_offvtx_shw1=-999;
    double dist_mainvtx_offvtx_shw1=-999;
    double len_ovs1 = -1e9;
    Vector mom_offvtx_track1(0,0,0), mom_offvtx_shw1(0,0,0);

    Point main_vtx_pt = vtx_fit_pt(ssm_main_vtx);

    for (auto [eit,end] = boost::edges(graph); eit != end; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg) continue;

        // skip segments directly connected to ssm_main_vtx or ssm_second_vtx
        bool at_main = false, at_second = false;
        if (ssm_main_vtx && ssm_main_vtx->descriptor_valid()) {
            for (auto [e2,e2e] = boost::out_edges(ssm_main_vtx->get_descriptor(), graph); e2 != e2e; ++e2)
                if (graph[*e2].segment == sg) { at_main = true; break; }
        }
        if (!at_main && ssm_second_vtx && ssm_second_vtx->descriptor_valid()) {
            for (auto [e2,e2e] = boost::out_edges(ssm_second_vtx->get_descriptor(), graph); e2 != e2e; ++e2)
                if (graph[*e2].segment == sg) { at_second = true; break; }
        }
        if (at_main || at_second) continue;

        double sg_len  = segment_track_length(sg) / units::cm;
        int    sg_pdg  = std::abs(sg->has_particle_info() ? sg->particle_info()->pdg() : 0);
        double sg_dlen = segment_track_direct_length(sg) / units::cm;
        double sg_med  = segment_median_dQ_dx(sg) / (43e3/units::cm);
        double sg_mdev = segment_track_max_deviation(sg) / units::cm;
        double sg_ke_cal = segment_cal_kine_dQdx(sg, recomb_model);
        double sg_ke_rng = cal_kine_range(sg_len, sg_pdg ? sg_pdg : 13, particle_data);
        double sg_ke_rng_mu = cal_kine_range(sg_len, 13,   particle_data);
        double sg_ke_rng_p  = cal_kine_range(sg_len, 2212, particle_data);
        double sg_ke_rng_e  = cal_kine_range(sg_len, 11,   particle_data);
        Vector sg_dir = segment_cal_dir_3vector(sg);

        // distance from ssm_main_vtx to segment start
        Point seg_start = sg->fits().empty() ? Point(0,0,0) : sg->fits().front().point;
        double sep = (seg_start - main_vtx_pt).magnitude() / units::cm;
        if (sep > 80.0) continue;

        off_vtx_length += sg_len;
        off_vtx_energy += sg_ke_cal;

        bool is_shower_type = (sg_pdg == 22 || sg_pdg == 11);
        if (sg_len > 1)  { if (!is_shower_type) n_offvtx_tracks_1++;  else n_offvtx_showers_1++;  }
        if (sg_len > 3)  { if (!is_shower_type) n_offvtx_tracks_3++;  else n_offvtx_showers_3++;  }
        if (sg_len > 5)  { if (!is_shower_type) n_offvtx_tracks_5++;  else n_offvtx_showers_5++;  }
        if (sg_len > 8)  { if (!is_shower_type) n_offvtx_tracks_8++;  else n_offvtx_showers_8++;  }
        if (sg_len > 11) { if (!is_shower_type) n_offvtx_tracks_11++; else n_offvtx_showers_11++; }

        if (!is_shower_type && sg_len > len_ovt1) {
            len_ovt1 = sg_len;
            auto sc2 = get_scores(sg, particle_data);
            double smu_f=sc2[0], sp_f=sc2[1], se_f=sc2[2], smu_b=sc2[3], sp_b=sc2[4], se_b=sc2[5];
            if (sg->dirsign() == -1) { std::swap(smu_f,smu_b); std::swap(sp_f,sp_b); std::swap(se_f,se_b); }
            pdg_offvtx_track1=sg_pdg; score_mu_fwd_offvtx_track1=smu_f; score_p_fwd_offvtx_track1=sp_f; score_e_fwd_offvtx_track1=se_f;
            score_mu_bck_offvtx_track1=smu_b; score_p_bck_offvtx_track1=sp_b; score_e_bck_offvtx_track1=se_b;
            length_offvtx_track1=sg_len; direct_length_offvtx_track1=sg_dlen; max_dev_offvtx_track1=sg_mdev;
            kine_energy_cal_offvtx_track1=sg_ke_cal; kine_energy_range_offvtx_track1=sg_ke_rng;
            kine_energy_range_mu_offvtx_track1=sg_ke_rng_mu; kine_energy_range_p_offvtx_track1=sg_ke_rng_p; kine_energy_range_e_offvtx_track1=sg_ke_rng_e;
            medium_dq_dx_offvtx_track1=sg_med; dist_mainvtx_offvtx_track1=sep;
            double m = get_mass(sg_pdg ? sg_pdg : 13); if (m<0) m=0;
            double pmag = std::sqrt(sg_ke_rng*sg_ke_rng + 2.0*sg_ke_rng*m);
            mom_offvtx_track1 = sg_dir * pmag;
            x_dir_offvtx_track1=sg_dir.x(); y_dir_offvtx_track1=sg_dir.y(); z_dir_offvtx_track1=sg_dir.z();
        } else if (is_shower_type && sg_len > len_ovs1) {
            len_ovs1 = sg_len;
            auto sc2 = get_scores(sg, particle_data);
            double smu_f=sc2[0], sp_f=sc2[1], se_f=sc2[2], smu_b=sc2[3], sp_b=sc2[4], se_b=sc2[5];
            if (sg->dirsign() == -1) { std::swap(smu_f,smu_b); std::swap(sp_f,sp_b); std::swap(se_f,se_b); }
            pdg_offvtx_shw1=sg_pdg; score_mu_fwd_offvtx_shw1=smu_f; score_p_fwd_offvtx_shw1=sp_f; score_e_fwd_offvtx_shw1=se_f;
            score_mu_bck_offvtx_shw1=smu_b; score_p_bck_offvtx_shw1=sp_b; score_e_bck_offvtx_shw1=se_b;
            length_offvtx_shw1=sg_len; direct_length_offvtx_shw1=sg_dlen; max_dev_offvtx_shw1=sg_mdev;
            kine_energy_cal_offvtx_shw1=sg_ke_cal; kine_energy_range_offvtx_shw1=sg_ke_rng;
            kine_energy_best_offvtx_shw1=sg_ke_rng;
            auto it_sh = map_segment_in_shower.find(sg);
            if (it_sh != map_segment_in_shower.end()) {
                kine_energy_best_offvtx_shw1 = it_sh->second->get_kine_best();
                if (kine_energy_best_offvtx_shw1==0.0) kine_energy_best_offvtx_shw1 = it_sh->second->get_kine_charge();
            }
            kine_energy_range_mu_offvtx_shw1=sg_ke_rng_mu; kine_energy_range_p_offvtx_shw1=sg_ke_rng_p; kine_energy_range_e_offvtx_shw1=sg_ke_rng_e;
            medium_dq_dx_offvtx_shw1=sg_med; dist_mainvtx_offvtx_shw1=sep;
            double m = get_mass(sg_pdg ? sg_pdg : 11); if (m<0) m=0;
            double pmag = std::sqrt(kine_energy_best_offvtx_shw1*kine_energy_best_offvtx_shw1 + 2.0*kine_energy_best_offvtx_shw1*m);
            mom_offvtx_shw1 = sg_dir * pmag;
            x_dir_offvtx_shw1=sg_dir.x(); y_dir_offvtx_shw1=sg_dir.y(); z_dir_offvtx_shw1=sg_dir.z();
        }
    }

    // finalize neutrino angles (include off-vertex momentum)
    nu_all = nu_all + mom_offvtx_track1 + mom_offvtx_shw1;
    compute_angles(nu_all,  nu_angle_z,      nu_angle_target,      nu_angle_absorber,      nu_angle_vertical);
    compute_angles(con_nu,  con_nu_angle_z,  con_nu_angle_target,  con_nu_angle_absorber,  con_nu_angle_vertical);
    compute_angles(prim_nu, prim_nu_angle_z, prim_nu_angle_target, prim_nu_angle_absorber, prim_nu_angle_vertical);
    compute_angles(track_nu,track_angle_z,   track_angle_target,   track_angle_absorber,   track_angle_vertical);

    // ------------------------------------------------------------------
    // Phase E: SSMSP filling + flag_st_kdar + TaggerInfo assignment
    // Prototype lines 1750-2289
    // ------------------------------------------------------------------

    // fill ssmsp if flag_ssmsp >= 0 (success path)
    if (flag_ssmsp >= 0) {
        fill_ssmsp_all(graph, main_vertex,
                       ssm_main_vtx, ssm_second_vtx, ssm_sg, dir,
                       showers, map_vertex_in_shower, map_segment_in_shower,
                       acc_segment_id, particle_data, ti);
    }

    // flag_st_kdar: KDAR cut-based decision (prototype lines 1871-1882)
    bool flag_st_kdar = false;
    if (pb_prim.n_tracks_1==0 && pb_prim.n_all_3==0 && pb_daught.n_tracks_5==0 &&
        pb_daught.n_all_5<2 && Nsm_wivtx==1 &&
        !(pio_kine.mass > 70 && pio_kine.mass < 200))
        flag_st_kdar = true;

    // --- Assign TaggerInfo fields (prototype lines 1885-2284) ---
    ti.ssm_flag_st_kdar = flag_st_kdar ? 1.0f : 0.0f;

    // dQ/dx profile
    ti.ssm_dq_dx_fwd_1 = (float)dq_dx_fwd_1; ti.ssm_dq_dx_fwd_2 = (float)dq_dx_fwd_2;
    ti.ssm_dq_dx_fwd_3 = (float)dq_dx_fwd_3; ti.ssm_dq_dx_fwd_4 = (float)dq_dx_fwd_4; ti.ssm_dq_dx_fwd_5 = (float)dq_dx_fwd_5;
    ti.ssm_dq_dx_bck_1 = (float)dq_dx_bck_1; ti.ssm_dq_dx_bck_2 = (float)dq_dx_bck_2;
    ti.ssm_dq_dx_bck_3 = (float)dq_dx_bck_3; ti.ssm_dq_dx_bck_4 = (float)dq_dx_bck_4; ti.ssm_dq_dx_bck_5 = (float)dq_dx_bck_5;
    ti.ssm_d_dq_dx_fwd_12 = (float)d_dq_dx_fwd_12; ti.ssm_d_dq_dx_fwd_23 = (float)d_dq_dx_fwd_23;
    ti.ssm_d_dq_dx_fwd_34 = (float)d_dq_dx_fwd_34; ti.ssm_d_dq_dx_fwd_45 = (float)d_dq_dx_fwd_45;
    ti.ssm_d_dq_dx_bck_12 = (float)d_dq_dx_bck_12; ti.ssm_d_dq_dx_bck_23 = (float)d_dq_dx_bck_23;
    ti.ssm_d_dq_dx_bck_34 = (float)d_dq_dx_bck_34; ti.ssm_d_dq_dx_bck_45 = (float)d_dq_dx_bck_45;
    ti.ssm_max_dq_dx_fwd_3 = (float)max_dqdx_fwd3; ti.ssm_max_dq_dx_fwd_5 = (float)max_dqdx_fwd5;
    ti.ssm_max_dq_dx_bck_3 = (float)max_dqdx_bck3; ti.ssm_max_dq_dx_bck_5 = (float)max_dqdx_bck5;
    ti.ssm_max_d_dq_dx_fwd_3 = (float)max_d_fwd3; ti.ssm_max_d_dq_dx_fwd_5 = (float)max_d_fwd5;
    ti.ssm_max_d_dq_dx_bck_3 = (float)max_d_bck3; ti.ssm_max_d_dq_dx_bck_5 = (float)max_d_bck5;
    ti.ssm_medium_dq_dx = (float)medium_dq_dx; ti.ssm_medium_dq_dx_bp = (float)medium_dq_dx_bp;
    // angles (use 10cm direction, matching prototype ssm_angle_to_z = angle_to_z_10)
    ti.ssm_angle_to_z       = (float)angle_to_z_10;
    ti.ssm_angle_to_target  = (float)angle_to_target_10;
    ti.ssm_angle_to_absorber= (float)angle_to_absorber_10;
    ti.ssm_angle_to_vertical= (float)angle_to_vertical_10;
    ti.ssm_x_dir = (float)x_dir; ti.ssm_y_dir = (float)y_dir; ti.ssm_z_dir = (float)z_dir;
    ti.ssm_kine_energy         = (float)kine_energy;
    ti.ssm_kine_energy_reduced = (float)kine_energy_reduced;
    ti.ssm_vtx_activity = vtx_activity ? 1.0f : 0.0f;
    ti.ssm_pdg = (float)pdg;
    ti.ssm_dQ_dx_cut = (float)dQ_dx_cut;
    ti.ssm_score_mu_fwd = (float)score_mu_fwd; ti.ssm_score_p_fwd = (float)score_p_fwd; ti.ssm_score_e_fwd = (float)score_e_fwd;
    ti.ssm_score_mu_bck = (float)score_mu_bck; ti.ssm_score_p_bck = (float)score_p_bck; ti.ssm_score_e_bck = (float)score_e_bck;
    ti.ssm_score_mu_fwd_bp = (float)score_mu_fwd_bp; ti.ssm_score_p_fwd_bp = (float)score_p_fwd_bp; ti.ssm_score_e_fwd_bp = (float)score_e_fwd_bp;
    ti.ssm_length = (float)length; ti.ssm_direct_length = (float)direct_length;
    ti.ssm_length_ratio = (length > 0 && direct_length > 0) ? (float)(direct_length/length) : -999.0f;
    ti.ssm_max_dev = (float)max_dev;

    // particle counts
    ti.ssm_n_prim_tracks_1=(float)pb_prim.n_tracks_1; ti.ssm_n_prim_tracks_3=(float)pb_prim.n_tracks_3; ti.ssm_n_prim_tracks_5=(float)pb_prim.n_tracks_5;
    ti.ssm_n_prim_tracks_8=(float)pb_prim.n_tracks_8; ti.ssm_n_prim_tracks_11=(float)pb_prim.n_tracks_11;
    ti.ssm_n_all_tracks_1=(float)pb_prim.n_all_1; ti.ssm_n_all_tracks_3=(float)pb_prim.n_all_3; ti.ssm_n_all_tracks_5=(float)pb_prim.n_all_5;
    ti.ssm_n_all_tracks_8=(float)pb_prim.n_all_8; ti.ssm_n_all_tracks_11=(float)pb_prim.n_all_11;
    ti.ssm_n_daughter_tracks_1=(float)pb_daught.n_tracks_1; ti.ssm_n_daughter_tracks_3=(float)pb_daught.n_tracks_3; ti.ssm_n_daughter_tracks_5=(float)pb_daught.n_tracks_5;
    ti.ssm_n_daughter_tracks_8=(float)pb_daught.n_tracks_8; ti.ssm_n_daughter_tracks_11=(float)pb_daught.n_tracks_11;
    ti.ssm_n_daughter_all_1=(float)pb_daught.n_all_1; ti.ssm_n_daughter_all_3=(float)pb_daught.n_all_3; ti.ssm_n_daughter_all_5=(float)pb_daught.n_all_5;
    ti.ssm_n_daughter_all_8=(float)pb_daught.n_all_8; ti.ssm_n_daughter_all_11=(float)pb_daught.n_all_11;

    // prim track1
    auto lr = [](double dl, double l) -> float { return (dl>0&&l>0)?(float)(dl/l):-999.f; };
    ti.ssm_prim_track1_pdg=(float)pb_prim.track1.pdg;
    ti.ssm_prim_track1_score_mu_fwd=(float)pb_prim.track1.score_mu_fwd; ti.ssm_prim_track1_score_p_fwd=(float)pb_prim.track1.score_p_fwd; ti.ssm_prim_track1_score_e_fwd=(float)pb_prim.track1.score_e_fwd;
    ti.ssm_prim_track1_score_mu_bck=(float)pb_prim.track1.score_mu_bck; ti.ssm_prim_track1_score_p_bck=(float)pb_prim.track1.score_p_bck; ti.ssm_prim_track1_score_e_bck=(float)pb_prim.track1.score_e_bck;
    ti.ssm_prim_track1_length=(float)pb_prim.track1.length; ti.ssm_prim_track1_direct_length=(float)pb_prim.track1.direct_length;
    ti.ssm_prim_track1_length_ratio=lr(pb_prim.track1.direct_length,pb_prim.track1.length);
    ti.ssm_prim_track1_max_dev=(float)pb_prim.track1.max_dev;
    ti.ssm_prim_track1_kine_energy_range=(float)pb_prim.track1.ke_rng;
    ti.ssm_prim_track1_kine_energy_range_mu=(float)pb_prim.track1.ke_rng_mu;
    ti.ssm_prim_track1_kine_energy_range_p=(float)pb_prim.track1.ke_rng_p;
    ti.ssm_prim_track1_kine_energy_range_e=(float)pb_prim.track1.ke_rng_e;
    ti.ssm_prim_track1_kine_energy_cal=(float)pb_prim.track1.ke_cal;
    ti.ssm_prim_track1_medium_dq_dx=(float)pb_prim.track1.medium_dq_dx;
    ti.ssm_prim_track1_x_dir=(float)pb_prim.track1.dir.x(); ti.ssm_prim_track1_y_dir=(float)pb_prim.track1.dir.y(); ti.ssm_prim_track1_z_dir=(float)pb_prim.track1.dir.z();
    ti.ssm_prim_track1_add_daught_track_counts_1=(float)pb_prim.track1.add_daught_track_1; ti.ssm_prim_track1_add_daught_all_counts_1=(float)pb_prim.track1.add_daught_all_1;
    ti.ssm_prim_track1_add_daught_track_counts_5=(float)pb_prim.track1.add_daught_track_5; ti.ssm_prim_track1_add_daught_all_counts_5=(float)pb_prim.track1.add_daught_all_5;
    ti.ssm_prim_track1_add_daught_track_counts_11=(float)pb_prim.track1.add_daught_track_11; ti.ssm_prim_track1_add_daught_all_counts_11=(float)pb_prim.track1.add_daught_all_11;

    // prim track2
    ti.ssm_prim_track2_pdg=(float)pb_prim.track2.pdg;
    ti.ssm_prim_track2_score_mu_fwd=(float)pb_prim.track2.score_mu_fwd; ti.ssm_prim_track2_score_p_fwd=(float)pb_prim.track2.score_p_fwd; ti.ssm_prim_track2_score_e_fwd=(float)pb_prim.track2.score_e_fwd;
    ti.ssm_prim_track2_score_mu_bck=(float)pb_prim.track2.score_mu_bck; ti.ssm_prim_track2_score_p_bck=(float)pb_prim.track2.score_p_bck; ti.ssm_prim_track2_score_e_bck=(float)pb_prim.track2.score_e_bck;
    ti.ssm_prim_track2_length=(float)pb_prim.track2.length; ti.ssm_prim_track2_direct_length=(float)pb_prim.track2.direct_length;
    ti.ssm_prim_track2_length_ratio=lr(pb_prim.track2.direct_length,pb_prim.track2.length);
    ti.ssm_prim_track2_max_dev=(float)pb_prim.track2.max_dev;
    ti.ssm_prim_track2_kine_energy_range=(float)pb_prim.track2.ke_rng;
    ti.ssm_prim_track2_kine_energy_range_mu=(float)pb_prim.track2.ke_rng_mu;
    ti.ssm_prim_track2_kine_energy_range_p=(float)pb_prim.track2.ke_rng_p;
    ti.ssm_prim_track2_kine_energy_range_e=(float)pb_prim.track2.ke_rng_e;
    ti.ssm_prim_track2_kine_energy_cal=(float)pb_prim.track2.ke_cal;
    ti.ssm_prim_track2_medium_dq_dx=(float)pb_prim.track2.medium_dq_dx;
    ti.ssm_prim_track2_x_dir=(float)pb_prim.track2.dir.x(); ti.ssm_prim_track2_y_dir=(float)pb_prim.track2.dir.y(); ti.ssm_prim_track2_z_dir=(float)pb_prim.track2.dir.z();
    ti.ssm_prim_track2_add_daught_track_counts_1=(float)pb_prim.track2.add_daught_track_1; ti.ssm_prim_track2_add_daught_all_counts_1=(float)pb_prim.track2.add_daught_all_1;
    ti.ssm_prim_track2_add_daught_track_counts_5=(float)pb_prim.track2.add_daught_track_5; ti.ssm_prim_track2_add_daught_all_counts_5=(float)pb_prim.track2.add_daught_all_5;
    ti.ssm_prim_track2_add_daught_track_counts_11=(float)pb_prim.track2.add_daught_track_11; ti.ssm_prim_track2_add_daught_all_counts_11=(float)pb_prim.track2.add_daught_all_11;

    // daught track1
    ti.ssm_daught_track1_pdg=(float)pb_daught.track1.pdg;
    ti.ssm_daught_track1_score_mu_fwd=(float)pb_daught.track1.score_mu_fwd; ti.ssm_daught_track1_score_p_fwd=(float)pb_daught.track1.score_p_fwd; ti.ssm_daught_track1_score_e_fwd=(float)pb_daught.track1.score_e_fwd;
    ti.ssm_daught_track1_score_mu_bck=(float)pb_daught.track1.score_mu_bck; ti.ssm_daught_track1_score_p_bck=(float)pb_daught.track1.score_p_bck; ti.ssm_daught_track1_score_e_bck=(float)pb_daught.track1.score_e_bck;
    ti.ssm_daught_track1_length=(float)pb_daught.track1.length; ti.ssm_daught_track1_direct_length=(float)pb_daught.track1.direct_length;
    ti.ssm_daught_track1_length_ratio=lr(pb_daught.track1.direct_length,pb_daught.track1.length);
    ti.ssm_daught_track1_max_dev=(float)pb_daught.track1.max_dev;
    ti.ssm_daught_track1_kine_energy_range=(float)pb_daught.track1.ke_rng;
    ti.ssm_daught_track1_kine_energy_range_mu=(float)pb_daught.track1.ke_rng_mu;
    ti.ssm_daught_track1_kine_energy_range_p=(float)pb_daught.track1.ke_rng_p;
    ti.ssm_daught_track1_kine_energy_range_e=(float)pb_daught.track1.ke_rng_e;
    ti.ssm_daught_track1_kine_energy_cal=(float)pb_daught.track1.ke_cal;
    ti.ssm_daught_track1_medium_dq_dx=(float)pb_daught.track1.medium_dq_dx;
    ti.ssm_daught_track1_x_dir=(float)pb_daught.track1.dir.x(); ti.ssm_daught_track1_y_dir=(float)pb_daught.track1.dir.y(); ti.ssm_daught_track1_z_dir=(float)pb_daught.track1.dir.z();
    ti.ssm_daught_track1_add_daught_track_counts_1=(float)pb_daught.track1.add_daught_track_1; ti.ssm_daught_track1_add_daught_all_counts_1=(float)pb_daught.track1.add_daught_all_1;
    ti.ssm_daught_track1_add_daught_track_counts_5=(float)pb_daught.track1.add_daught_track_5; ti.ssm_daught_track1_add_daught_all_counts_5=(float)pb_daught.track1.add_daught_all_5;
    ti.ssm_daught_track1_add_daught_track_counts_11=(float)pb_daught.track1.add_daught_track_11; ti.ssm_daught_track1_add_daught_all_counts_11=(float)pb_daught.track1.add_daught_all_11;

    // daught track2
    ti.ssm_daught_track2_pdg=(float)pb_daught.track2.pdg;
    ti.ssm_daught_track2_score_mu_fwd=(float)pb_daught.track2.score_mu_fwd; ti.ssm_daught_track2_score_p_fwd=(float)pb_daught.track2.score_p_fwd; ti.ssm_daught_track2_score_e_fwd=(float)pb_daught.track2.score_e_fwd;
    ti.ssm_daught_track2_score_mu_bck=(float)pb_daught.track2.score_mu_bck; ti.ssm_daught_track2_score_p_bck=(float)pb_daught.track2.score_p_bck; ti.ssm_daught_track2_score_e_bck=(float)pb_daught.track2.score_e_bck;
    ti.ssm_daught_track2_length=(float)pb_daught.track2.length; ti.ssm_daught_track2_direct_length=(float)pb_daught.track2.direct_length;
    ti.ssm_daught_track2_length_ratio=lr(pb_daught.track2.direct_length,pb_daught.track2.length);
    ti.ssm_daught_track2_max_dev=(float)pb_daught.track2.max_dev;
    ti.ssm_daught_track2_kine_energy_range=(float)pb_daught.track2.ke_rng;
    ti.ssm_daught_track2_kine_energy_range_mu=(float)pb_daught.track2.ke_rng_mu;
    ti.ssm_daught_track2_kine_energy_range_p=(float)pb_daught.track2.ke_rng_p;
    ti.ssm_daught_track2_kine_energy_range_e=(float)pb_daught.track2.ke_rng_e;
    ti.ssm_daught_track2_kine_energy_cal=(float)pb_daught.track2.ke_cal;
    ti.ssm_daught_track2_medium_dq_dx=(float)pb_daught.track2.medium_dq_dx;
    ti.ssm_daught_track2_x_dir=(float)pb_daught.track2.dir.x(); ti.ssm_daught_track2_y_dir=(float)pb_daught.track2.dir.y(); ti.ssm_daught_track2_z_dir=(float)pb_daught.track2.dir.z();
    ti.ssm_daught_track2_add_daught_track_counts_1=(float)pb_daught.track2.add_daught_track_1; ti.ssm_daught_track2_add_daught_all_counts_1=(float)pb_daught.track2.add_daught_all_1;
    ti.ssm_daught_track2_add_daught_track_counts_5=(float)pb_daught.track2.add_daught_track_5; ti.ssm_daught_track2_add_daught_all_counts_5=(float)pb_daught.track2.add_daught_all_5;
    ti.ssm_daught_track2_add_daught_track_counts_11=(float)pb_daught.track2.add_daught_track_11; ti.ssm_daught_track2_add_daught_all_counts_11=(float)pb_daught.track2.add_daught_all_11;

    // prim shw1
    ti.ssm_prim_shw1_pdg=(float)pb_prim.shw1.pdg;
    ti.ssm_prim_shw1_score_mu_fwd=(float)pb_prim.shw1.score_mu_fwd; ti.ssm_prim_shw1_score_p_fwd=(float)pb_prim.shw1.score_p_fwd; ti.ssm_prim_shw1_score_e_fwd=(float)pb_prim.shw1.score_e_fwd;
    ti.ssm_prim_shw1_score_mu_bck=(float)pb_prim.shw1.score_mu_bck; ti.ssm_prim_shw1_score_p_bck=(float)pb_prim.shw1.score_p_bck; ti.ssm_prim_shw1_score_e_bck=(float)pb_prim.shw1.score_e_bck;
    ti.ssm_prim_shw1_length=(float)pb_prim.shw1.length; ti.ssm_prim_shw1_direct_length=(float)pb_prim.shw1.direct_length;
    ti.ssm_prim_shw1_length_ratio=lr(pb_prim.shw1.direct_length,pb_prim.shw1.length);
    ti.ssm_prim_shw1_max_dev=(float)pb_prim.shw1.max_dev;
    ti.ssm_prim_shw1_kine_energy_best=(float)pb_prim.shw1.ke_best;
    ti.ssm_prim_shw1_kine_energy_range=(float)pb_prim.shw1.ke_rng;
    ti.ssm_prim_shw1_kine_energy_range_mu=(float)pb_prim.shw1.ke_rng_mu;
    ti.ssm_prim_shw1_kine_energy_range_p=(float)pb_prim.shw1.ke_rng_p;
    ti.ssm_prim_shw1_kine_energy_range_e=(float)pb_prim.shw1.ke_rng_e;
    ti.ssm_prim_shw1_kine_energy_cal=(float)pb_prim.shw1.ke_cal;
    ti.ssm_prim_shw1_medium_dq_dx=(float)pb_prim.shw1.medium_dq_dx;
    ti.ssm_prim_shw1_x_dir=(float)pb_prim.shw1.dir.x(); ti.ssm_prim_shw1_y_dir=(float)pb_prim.shw1.dir.y(); ti.ssm_prim_shw1_z_dir=(float)pb_prim.shw1.dir.z();
    ti.ssm_prim_shw1_add_daught_track_counts_1=(float)pb_prim.shw1.add_daught_track_1; ti.ssm_prim_shw1_add_daught_all_counts_1=(float)pb_prim.shw1.add_daught_all_1;
    ti.ssm_prim_shw1_add_daught_track_counts_5=(float)pb_prim.shw1.add_daught_track_5; ti.ssm_prim_shw1_add_daught_all_counts_5=(float)pb_prim.shw1.add_daught_all_5;
    ti.ssm_prim_shw1_add_daught_track_counts_11=(float)pb_prim.shw1.add_daught_track_11; ti.ssm_prim_shw1_add_daught_all_counts_11=(float)pb_prim.shw1.add_daught_all_11;

    // prim shw2
    ti.ssm_prim_shw2_pdg=(float)pb_prim.shw2.pdg;
    ti.ssm_prim_shw2_score_mu_fwd=(float)pb_prim.shw2.score_mu_fwd; ti.ssm_prim_shw2_score_p_fwd=(float)pb_prim.shw2.score_p_fwd; ti.ssm_prim_shw2_score_e_fwd=(float)pb_prim.shw2.score_e_fwd;
    ti.ssm_prim_shw2_score_mu_bck=(float)pb_prim.shw2.score_mu_bck; ti.ssm_prim_shw2_score_p_bck=(float)pb_prim.shw2.score_p_bck; ti.ssm_prim_shw2_score_e_bck=(float)pb_prim.shw2.score_e_bck;
    ti.ssm_prim_shw2_length=(float)pb_prim.shw2.length; ti.ssm_prim_shw2_direct_length=(float)pb_prim.shw2.direct_length;
    ti.ssm_prim_shw2_length_ratio=lr(pb_prim.shw2.direct_length,pb_prim.shw2.length);
    ti.ssm_prim_shw2_max_dev=(float)pb_prim.shw2.max_dev;
    ti.ssm_prim_shw2_kine_energy_best=(float)pb_prim.shw2.ke_best;
    ti.ssm_prim_shw2_kine_energy_range=(float)pb_prim.shw2.ke_rng;
    ti.ssm_prim_shw2_kine_energy_range_mu=(float)pb_prim.shw2.ke_rng_mu;
    ti.ssm_prim_shw2_kine_energy_range_p=(float)pb_prim.shw2.ke_rng_p;
    ti.ssm_prim_shw2_kine_energy_range_e=(float)pb_prim.shw2.ke_rng_e;
    ti.ssm_prim_shw2_kine_energy_cal=(float)pb_prim.shw2.ke_cal;
    ti.ssm_prim_shw2_medium_dq_dx=(float)pb_prim.shw2.medium_dq_dx;
    ti.ssm_prim_shw2_x_dir=(float)pb_prim.shw2.dir.x(); ti.ssm_prim_shw2_y_dir=(float)pb_prim.shw2.dir.y(); ti.ssm_prim_shw2_z_dir=(float)pb_prim.shw2.dir.z();
    ti.ssm_prim_shw2_add_daught_track_counts_1=(float)pb_prim.shw2.add_daught_track_1; ti.ssm_prim_shw2_add_daught_all_counts_1=(float)pb_prim.shw2.add_daught_all_1;
    ti.ssm_prim_shw2_add_daught_track_counts_5=(float)pb_prim.shw2.add_daught_track_5; ti.ssm_prim_shw2_add_daught_all_counts_5=(float)pb_prim.shw2.add_daught_all_5;
    ti.ssm_prim_shw2_add_daught_track_counts_11=(float)pb_prim.shw2.add_daught_track_11; ti.ssm_prim_shw2_add_daught_all_counts_11=(float)pb_prim.shw2.add_daught_all_11;

    // daught shw1
    ti.ssm_daught_shw1_pdg=(float)pb_daught.shw1.pdg;
    ti.ssm_daught_shw1_score_mu_fwd=(float)pb_daught.shw1.score_mu_fwd; ti.ssm_daught_shw1_score_p_fwd=(float)pb_daught.shw1.score_p_fwd; ti.ssm_daught_shw1_score_e_fwd=(float)pb_daught.shw1.score_e_fwd;
    ti.ssm_daught_shw1_score_mu_bck=(float)pb_daught.shw1.score_mu_bck; ti.ssm_daught_shw1_score_p_bck=(float)pb_daught.shw1.score_p_bck; ti.ssm_daught_shw1_score_e_bck=(float)pb_daught.shw1.score_e_bck;
    ti.ssm_daught_shw1_length=(float)pb_daught.shw1.length; ti.ssm_daught_shw1_direct_length=(float)pb_daught.shw1.direct_length;
    ti.ssm_daught_shw1_length_ratio=lr(pb_daught.shw1.direct_length,pb_daught.shw1.length);
    ti.ssm_daught_shw1_max_dev=(float)pb_daught.shw1.max_dev;
    ti.ssm_daught_shw1_kine_energy_best=(float)pb_daught.shw1.ke_best;
    ti.ssm_daught_shw1_kine_energy_range=(float)pb_daught.shw1.ke_rng;
    ti.ssm_daught_shw1_kine_energy_range_mu=(float)pb_daught.shw1.ke_rng_mu;
    ti.ssm_daught_shw1_kine_energy_range_p=(float)pb_daught.shw1.ke_rng_p;
    ti.ssm_daught_shw1_kine_energy_range_e=(float)pb_daught.shw1.ke_rng_e;
    ti.ssm_daught_shw1_kine_energy_cal=(float)pb_daught.shw1.ke_cal;
    ti.ssm_daught_shw1_medium_dq_dx=(float)pb_daught.shw1.medium_dq_dx;
    ti.ssm_daught_shw1_x_dir=(float)pb_daught.shw1.dir.x(); ti.ssm_daught_shw1_y_dir=(float)pb_daught.shw1.dir.y(); ti.ssm_daught_shw1_z_dir=(float)pb_daught.shw1.dir.z();
    ti.ssm_daught_shw1_add_daught_track_counts_1=(float)pb_daught.shw1.add_daught_track_1; ti.ssm_daught_shw1_add_daught_all_counts_1=(float)pb_daught.shw1.add_daught_all_1;
    ti.ssm_daught_shw1_add_daught_track_counts_5=(float)pb_daught.shw1.add_daught_track_5; ti.ssm_daught_shw1_add_daught_all_counts_5=(float)pb_daught.shw1.add_daught_all_5;
    ti.ssm_daught_shw1_add_daught_track_counts_11=(float)pb_daught.shw1.add_daught_track_11; ti.ssm_daught_shw1_add_daught_all_counts_11=(float)pb_daught.shw1.add_daught_all_11;

    // daught shw2
    ti.ssm_daught_shw2_pdg=(float)pb_daught.shw2.pdg;
    ti.ssm_daught_shw2_score_mu_fwd=(float)pb_daught.shw2.score_mu_fwd; ti.ssm_daught_shw2_score_p_fwd=(float)pb_daught.shw2.score_p_fwd; ti.ssm_daught_shw2_score_e_fwd=(float)pb_daught.shw2.score_e_fwd;
    ti.ssm_daught_shw2_score_mu_bck=(float)pb_daught.shw2.score_mu_bck; ti.ssm_daught_shw2_score_p_bck=(float)pb_daught.shw2.score_p_bck; ti.ssm_daught_shw2_score_e_bck=(float)pb_daught.shw2.score_e_bck;
    ti.ssm_daught_shw2_length=(float)pb_daught.shw2.length; ti.ssm_daught_shw2_direct_length=(float)pb_daught.shw2.direct_length;
    ti.ssm_daught_shw2_length_ratio=lr(pb_daught.shw2.direct_length,pb_daught.shw2.length);
    ti.ssm_daught_shw2_max_dev=(float)pb_daught.shw2.max_dev;
    ti.ssm_daught_shw2_kine_energy_best=(float)pb_daught.shw2.ke_best;
    ti.ssm_daught_shw2_kine_energy_range=(float)pb_daught.shw2.ke_rng;
    ti.ssm_daught_shw2_kine_energy_range_mu=(float)pb_daught.shw2.ke_rng_mu;
    ti.ssm_daught_shw2_kine_energy_range_p=(float)pb_daught.shw2.ke_rng_p;
    ti.ssm_daught_shw2_kine_energy_range_e=(float)pb_daught.shw2.ke_rng_e;
    ti.ssm_daught_shw2_kine_energy_cal=(float)pb_daught.shw2.ke_cal;
    ti.ssm_daught_shw2_medium_dq_dx=(float)pb_daught.shw2.medium_dq_dx;
    ti.ssm_daught_shw2_x_dir=(float)pb_daught.shw2.dir.x(); ti.ssm_daught_shw2_y_dir=(float)pb_daught.shw2.dir.y(); ti.ssm_daught_shw2_z_dir=(float)pb_daught.shw2.dir.z();
    ti.ssm_daught_shw2_add_daught_track_counts_1=(float)pb_daught.shw2.add_daught_track_1; ti.ssm_daught_shw2_add_daught_all_counts_1=(float)pb_daught.shw2.add_daught_all_1;
    ti.ssm_daught_shw2_add_daught_track_counts_5=(float)pb_daught.shw2.add_daught_track_5; ti.ssm_daught_shw2_add_daught_all_counts_5=(float)pb_daught.shw2.add_daught_all_5;
    ti.ssm_daught_shw2_add_daught_track_counts_11=(float)pb_daught.shw2.add_daught_track_11; ti.ssm_daught_shw2_add_daught_all_counts_11=(float)pb_daught.shw2.add_daught_all_11;

    // event-level angles
    ti.ssm_nu_angle_z=(float)nu_angle_z; ti.ssm_nu_angle_target=(float)nu_angle_target;
    ti.ssm_nu_angle_absorber=(float)nu_angle_absorber; ti.ssm_nu_angle_vertical=(float)nu_angle_vertical;
    ti.ssm_con_nu_angle_z=(float)con_nu_angle_z; ti.ssm_con_nu_angle_target=(float)con_nu_angle_target;
    ti.ssm_con_nu_angle_absorber=(float)con_nu_angle_absorber; ti.ssm_con_nu_angle_vertical=(float)con_nu_angle_vertical;
    ti.ssm_prim_nu_angle_z=(float)prim_nu_angle_z; ti.ssm_prim_nu_angle_target=(float)prim_nu_angle_target;
    ti.ssm_prim_nu_angle_absorber=(float)prim_nu_angle_absorber; ti.ssm_prim_nu_angle_vertical=(float)prim_nu_angle_vertical;
    ti.ssm_track_angle_z=(float)track_angle_z; ti.ssm_track_angle_target=(float)track_angle_target;
    ti.ssm_track_angle_absorber=(float)track_angle_absorber; ti.ssm_track_angle_vertical=(float)track_angle_vertical;

    Point smvp = vtx_fit_pt(ssm_main_vtx);
    ti.ssm_vtxX = (float)(smvp.x()/units::cm);
    ti.ssm_vtxY = (float)(smvp.y()/units::cm);
    ti.ssm_vtxZ = (float)(smvp.z()/units::cm);

    // off-vertex
    ti.ssm_offvtx_length=(float)off_vtx_length; ti.ssm_offvtx_energy=(float)off_vtx_energy;
    ti.ssm_n_offvtx_tracks_1=(float)n_offvtx_tracks_1; ti.ssm_n_offvtx_tracks_3=(float)n_offvtx_tracks_3; ti.ssm_n_offvtx_tracks_5=(float)n_offvtx_tracks_5;
    ti.ssm_n_offvtx_tracks_8=(float)n_offvtx_tracks_8; ti.ssm_n_offvtx_tracks_11=(float)n_offvtx_tracks_11;
    ti.ssm_n_offvtx_showers_1=(float)n_offvtx_showers_1; ti.ssm_n_offvtx_showers_3=(float)n_offvtx_showers_3; ti.ssm_n_offvtx_showers_5=(float)n_offvtx_showers_5;
    ti.ssm_n_offvtx_showers_8=(float)n_offvtx_showers_8; ti.ssm_n_offvtx_showers_11=(float)n_offvtx_showers_11;
    ti.ssm_offvtx_track1_pdg=(float)pdg_offvtx_track1;
    ti.ssm_offvtx_track1_score_mu_fwd=(float)score_mu_fwd_offvtx_track1; ti.ssm_offvtx_track1_score_p_fwd=(float)score_p_fwd_offvtx_track1; ti.ssm_offvtx_track1_score_e_fwd=(float)score_e_fwd_offvtx_track1;
    ti.ssm_offvtx_track1_score_mu_bck=(float)score_mu_bck_offvtx_track1; ti.ssm_offvtx_track1_score_p_bck=(float)score_p_bck_offvtx_track1; ti.ssm_offvtx_track1_score_e_bck=(float)score_e_bck_offvtx_track1;
    ti.ssm_offvtx_track1_length=(float)length_offvtx_track1; ti.ssm_offvtx_track1_direct_length=(float)direct_length_offvtx_track1;
    ti.ssm_offvtx_track1_max_dev=(float)max_dev_offvtx_track1;
    ti.ssm_offvtx_track1_kine_energy_range=(float)kine_energy_range_offvtx_track1;
    ti.ssm_offvtx_track1_kine_energy_range_mu=(float)kine_energy_range_mu_offvtx_track1;
    ti.ssm_offvtx_track1_kine_energy_range_p=(float)kine_energy_range_p_offvtx_track1;
    ti.ssm_offvtx_track1_kine_energy_range_e=(float)kine_energy_range_e_offvtx_track1;
    ti.ssm_offvtx_track1_kine_energy_cal=(float)kine_energy_cal_offvtx_track1;
    ti.ssm_offvtx_track1_medium_dq_dx=(float)medium_dq_dx_offvtx_track1;
    ti.ssm_offvtx_track1_x_dir=(float)x_dir_offvtx_track1; ti.ssm_offvtx_track1_y_dir=(float)y_dir_offvtx_track1; ti.ssm_offvtx_track1_z_dir=(float)z_dir_offvtx_track1;
    ti.ssm_offvtx_track1_dist_mainvtx=(float)dist_mainvtx_offvtx_track1;
    ti.ssm_offvtx_shw1_pdg_offvtx=(float)pdg_offvtx_shw1;
    ti.ssm_offvtx_shw1_score_mu_fwd=(float)score_mu_fwd_offvtx_shw1; ti.ssm_offvtx_shw1_score_p_fwd=(float)score_p_fwd_offvtx_shw1; ti.ssm_offvtx_shw1_score_e_fwd=(float)score_e_fwd_offvtx_shw1;
    ti.ssm_offvtx_shw1_score_mu_bck=(float)score_mu_bck_offvtx_shw1; ti.ssm_offvtx_shw1_score_p_bck=(float)score_p_bck_offvtx_shw1; ti.ssm_offvtx_shw1_score_e_bck=(float)score_e_bck_offvtx_shw1;
    ti.ssm_offvtx_shw1_length=(float)length_offvtx_shw1; ti.ssm_offvtx_shw1_direct_length=(float)direct_length_offvtx_shw1;
    ti.ssm_offvtx_shw1_max_dev=(float)max_dev_offvtx_shw1;
    ti.ssm_offvtx_shw1_kine_energy_best=(float)kine_energy_best_offvtx_shw1;
    ti.ssm_offvtx_shw1_kine_energy_range=(float)kine_energy_range_offvtx_shw1;
    ti.ssm_offvtx_shw1_kine_energy_range_mu=(float)kine_energy_range_mu_offvtx_shw1;
    ti.ssm_offvtx_shw1_kine_energy_range_p=(float)kine_energy_range_p_offvtx_shw1;
    ti.ssm_offvtx_shw1_kine_energy_range_e=(float)kine_energy_range_e_offvtx_shw1;
    ti.ssm_offvtx_shw1_kine_energy_cal=(float)kine_energy_cal_offvtx_shw1;
    ti.ssm_offvtx_shw1_medium_dq_dx=(float)medium_dq_dx_offvtx_shw1;
    ti.ssm_offvtx_shw1_x_dir=(float)x_dir_offvtx_shw1; ti.ssm_offvtx_shw1_y_dir=(float)y_dir_offvtx_shw1; ti.ssm_offvtx_shw1_z_dir=(float)z_dir_offvtx_shw1;
    ti.ssm_offvtx_shw1_dist_mainvtx=(float)dist_mainvtx_offvtx_shw1;

    return true;
}
