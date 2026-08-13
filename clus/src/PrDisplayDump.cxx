// One-JSON-per-event dump feeding sbnd_xin/pr_display.  See the header for why
// this exists rather than the display reading mabc-pr.zip + tracking-pr.root,
// and sbnd_xin/docs/pr/26_pr-event-display.md for the schema.
//
// The PR-graph walk here is a fork-by-duplication of the two producers that
// already exist for these quantities (CLAUDE.md M10 -- both stay untouched):
//   - MultiAlgBlobClustering::fill_bee_points_from_pr_graph / _vertices  (Bee
//     track_fit / shower_track / vertices layers)
//   - Root::SbndPrMagnifyTrackingVisitor::write_proj_data / write_t_rec_data
//     (Magnify T_proj_data / T_rec_charge)
// Divergences from those two are deliberate and marked DELTA below.

#include "WireCellClus/PrDisplayDump.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRShower.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <map>
#include <set>
#include <tuple>

WIRECELL_FACTORY(PrDisplayDump, WireCell::Clus::PrDisplayDump,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;

// The node-id encoding fill_bee_pf_tree uses for a particle-flow node
// (MultiAlgBlobClustering.cxx, seg_display_id): cluster_id*1000 + segment id,
// with the graph index standing in when the segment id was never assigned.
// NeutrinoPatternBase.cxx sets seg->set_id(edge_bundle.index), so in practice
// the two agree -- but keep the same fallback so a segment that missed that
// pass still lands on the same number the Bee producer would print.
//
// This is the JOIN KEY between an mc.json shower node and the segments of that
// shower.  It must not drift from the Bee producer.
static int pf_node_id(const WireCell::Clus::PR::SegmentPtr& seg)
{
    if (!seg) return -1;
    int sid = seg->id();
    if (sid < 0) sid = static_cast<int>(seg->get_graph_index());
    const auto* cl = seg->cluster();
    return cl ? cl->get_cluster_id() * 1000 + sid : sid;
}

// Same encoding as the vertices[].id field below (cluster_id*1000 +
// graph_index): the id a segment's start/end vertex would carry in the
// "vertices" array, so the viewer can look a segment's endpoints up there
// without a second join key.
static int vertex_display_id(const WireCell::Clus::PR::VertexPtr& vtx)
{
    if (!vtx) return -1;
    const auto* cl = vtx->cluster();
    const int cluster_id = cl ? cl->get_cluster_id() : 0;
    return cluster_id * 1000 + static_cast<int>(vtx->get_graph_index());
}

Clus::PrDisplayDump::PrDisplayDump()
  : Aux::Logger("PrDisplayDump", "clus")
{
}

Clus::PrDisplayDump::~PrDisplayDump() {}

WireCell::Configuration Clus::PrDisplayDump::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = m_output_filename;
    cfg["grouping"] = m_grouping_name;
    cfg["anodes"] = Json::arrayValue;
    cfg["detector_volumes"] = "";
    cfg["runNo"] = m_runNo;
    cfg["subRunNo"] = m_subRunNo;
    cfg["eventNo"] = m_eventNo;
    cfg["dQdx_scale"] = m_dQdx_scale;
    cfg["dQdx_offset"] = m_dQdx_offset;
    cfg["nticks"] = m_nticks;
    cfg["proj_charge_min"] = m_proj_charge_min;
    cfg["pretty"] = m_pretty;
    cfg["particle_dataset"] = m_particle_dataset_name;
    cfg["mip_dqdx_median"] = m_mip_dqdx_median;
    cfg["mip_dqdx_flat"] = m_mip_dqdx_flat;
    cfg["dqdx_ref_grid_n"] = m_dqdx_ref_grid_n;
    cfg["dqdx_ref_grid_step"] = m_dqdx_ref_grid_step;
    return cfg;
}

void Clus::PrDisplayDump::configure(const WireCell::Configuration& cfg)
{
    m_output_filename = get<std::string>(cfg, "output_filename", m_output_filename);
    m_grouping_name = get<std::string>(cfg, "grouping", m_grouping_name);
    m_runNo = get<int>(cfg, "runNo", m_runNo);
    m_subRunNo = get<int>(cfg, "subRunNo", m_subRunNo);
    m_eventNo = get<int>(cfg, "eventNo", m_eventNo);
    m_dQdx_scale = get<double>(cfg, "dQdx_scale", m_dQdx_scale);
    m_dQdx_offset = get<double>(cfg, "dQdx_offset", m_dQdx_offset);
    m_nticks = get<int>(cfg, "nticks", m_nticks);
    m_proj_charge_min = get<double>(cfg, "proj_charge_min", m_proj_charge_min);
    m_pretty = get<bool>(cfg, "pretty", m_pretty);
    m_particle_dataset_name = get<std::string>(cfg, "particle_dataset", m_particle_dataset_name);
    m_mip_dqdx_median = get<double>(cfg, "mip_dqdx_median", m_mip_dqdx_median);
    m_mip_dqdx_flat = get<double>(cfg, "mip_dqdx_flat", m_mip_dqdx_flat);
    m_dqdx_ref_grid_n = get<int>(cfg, "dqdx_ref_grid_n", m_dqdx_ref_grid_n);
    m_dqdx_ref_grid_step = get<double>(cfg, "dqdx_ref_grid_step", m_dqdx_ref_grid_step);
    // doc sbnd_xin/docs/pr/45 -- mirror of MABC bee_points pseudo_shower_track_paint.
    m_pseudo_shower_track_paint = get<bool>(cfg, "pseudo_shower_track_paint", m_pseudo_shower_track_paint);

    for (auto anode_tn : cfg["anodes"]) {
        m_anodes.push_back(Factory::find_tn<IAnodePlane>(anode_tn.asString()));
    }
    if (cfg.isMember("detector_volumes") && !cfg["detector_volumes"].asString().empty()) {
        m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());
    }

    // Best-effort: a pipeline that runs pr_display without tagger_check_neutrino
    // or tagger_check_stm ahead of it never compiled a ParticleDataSet instance
    // in (sbnd/clus.jsonnet's tagger_uses gate) -- find_maybe_tn returns nullptr
    // rather than throwing, and dump_dqdx_ref() omits "dqdx_ref" in that case.
    if (!m_particle_dataset_name.empty()) {
        auto configurable = Factory::find_maybe_tn<IConfigurable>(m_particle_dataset_name);
        m_particle_data = std::dynamic_pointer_cast<Clus::ParticleDataSet>(configurable);
        if (!m_particle_data) {
            log->warn("ParticleDataSet '{}' not found -- dqdx_ref will be omitted from the dump",
                      m_particle_dataset_name);
        }
    }
}

// Same derivation as SbndPrMagnifyTrackingVisitor::chan_scheme().  DELTA: we
// record it in "meta" for the viewer's axis labels but store PER-APA wire
// indices in the data, so the display never has to undo the concatenation.
Clus::PrDisplayDump::ChanScheme Clus::PrDisplayDump::chan_scheme() const
{
    ChanScheme cs;
    for (const auto& anode : m_anodes) {
        for (const auto& face : anode->faces()) {
            if (!face) continue;
            auto planes = face->planes();
            for (size_t p = 0; p < planes.size() && p < 3; ++p) {
                const int n = static_cast<int>(planes[p]->channels().size());
                if (n > cs.nch[p]) cs.nch[p] = n;
            }
        }
    }
    const int napa = static_cast<int>(m_anodes.size());
    cs.base[0] = 0;
    cs.base[1] = napa * cs.nch[0];
    cs.base[2] = cs.base[1] + napa * cs.nch[1];
    return cs;
}

void Clus::PrDisplayDump::visit(Facade::Ensemble& ensemble) const
{
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        log->debug("no grouping '{}', nothing dumped", m_grouping_name);
        return;
    }
    auto& grouping = *groupings.at(0);

    if (!m_anodes.empty()) grouping.set_anodes(m_anodes);
    if (m_dv) grouping.set_detector_volumes(m_dv);

    auto tf = grouping.get_track_fitting();
    if (!tf) {
        log->warn("no TrackFitting in grouping '{}' -- is this stage after tagger_check_neutrino?",
                  m_grouping_name);
        return;
    }

    const auto cs = chan_scheme();

    Configuration top;
    top["meta"] = dump_meta(grouping, cs);

    Configuration graph = dump_graph(grouping);
    top["segments"] = graph["segments"];
    top["vertices"] = graph["vertices"];
    top["main_vertex"] = graph["main_vertex"];

    top["showers"] = dump_showers(grouping);
    top["kine"] = dump_kine(grouping);
    top["tagger"] = dump_tagger(grouping);

    top["track_shower"] = dump_track_shower(grouping);
    top["steiner"] = dump_steiner(grouping);
    top["proj"] = dump_proj(grouping);
    top["dead"] = dump_dead(grouping);
    top["dqdx_ref"] = dump_dqdx_ref();
    top["vertex_scoreboard"] = dump_vertex_scoreboard(grouping);

    Persist::dump(m_output_filename, top, m_pretty);

    log->debug("wrote {}: {} segment(s), {} vertex(es), {} shower(s), {} steiner cluster(s), {} proj plane(s)",
               m_output_filename, top["segments"].size(), top["vertices"].size(),
               top["showers"].size(), top["steiner"].size(), top["proj"].size());
}

// Five dQ/dx-vs-residual-range reference curves, sampled on a fixed cm grid
// from the SAME ParticleDataSet::get_dEdx_function() instances do_track_pid
// uses.  Emitted per event (even though the tables are static) so an old
// arm's JSON keeps its own templates rather than acquiring today's.
//
// UNIT TRAP -- get this backwards and it is a silent, plausible-looking 10x
// offset, not a crash.  Measured empirically (not just read off the code):
// the muon plateau the particle_dataset.jsonnet header documents is 54657.7
// e/cm at rr=59.5cm; `scalar_function(rr_cm)` ALONE reproduces that number,
// so the raw table is already e/cm and needs NO further scaling here.
//
// do_track_pid divides by units::cm (PRSegmentFunctions.cxx segment_do_track_pid
// / do_track_comp) for a DIFFERENT reason that does not apply here: its data
// side is `fits[i].dQ / fits[i].dx` with dx in INTERNAL length units (not
// pre-divided by cm), i.e. its own comparison scale is e/cm scaled down by
// units::cm=10 to match that internal-unit dx.  The dump's per-point dQ/dx
// (points[].dQ / points[].dx) uses dx ALREADY pre-divided by cm (fit_json,
// above), so it is already plain physical e/cm -- copying do_track_comp's
// division here would divide the template but not the data it is compared
// against, and land the template a spurious 10x below the measured points.
Configuration Clus::PrDisplayDump::dump_dqdx_ref() const
{
    if (!m_particle_data) return Json::nullValue;

    static const std::array<std::pair<const char*, const char*>, 5> particles{{
        {"muon", "muon"}, {"proton", "proton"}, {"pion", "pion"},
        {"kaon", "kaon"}, {"electron", "electron"},
    }};

    Configuration out;
    out["units"] = "e/cm";
    out["source"] = "ParticleDataSet";
    Configuration grid;
    grid["start"] = 0.0;
    grid["step"] = m_dqdx_ref_grid_step;
    grid["n"] = m_dqdx_ref_grid_n;
    out["grid"] = grid;

    for (const auto& [key, name] : particles) {
        auto fn = m_particle_data->get_dEdx_function(name);
        if (!fn) continue;
        Configuration vals = Json::arrayValue;
        for (int i = 0; i < m_dqdx_ref_grid_n; ++i) {
            const double rr_cm = i * m_dqdx_ref_grid_step;
            vals.append(fn->scalar_function(rr_cm));
        }
        out[key] = vals;
    }
    return out;
}

Configuration Clus::PrDisplayDump::dump_meta(Facade::Grouping& grouping, const ChanScheme& cs) const
{
    Configuration meta;
    meta["runNo"] = m_runNo;
    meta["subRunNo"] = m_subRunNo;
    meta["eventNo"] = m_eventNo;
    meta["nticks"] = m_nticks;
    meta["dQdx_scale"] = m_dQdx_scale;
    meta["dQdx_offset"] = m_dQdx_offset;
    // Positions are in cm and charges are RAW (no dQdx_scale/offset applied) --
    // unlike the Bee layers, which pre-bake the affine transform into q.  The
    // display re-applies it itself so the raw dQ stays available for dQ/dx.
    meta["length_unit"] = "cm";
    meta["charge_transform"] = "none";
    // Display-only MIP scales for the dQ/dx panel (sbnd_xin/docs/pr/42) -- see
    // the header comment on m_mip_dqdx_median for why these are NOT read from
    // the taggers' own members.
    meta["mip_dqdx_median"] = m_mip_dqdx_median;
    meta["mip_dqdx_flat"] = m_mip_dqdx_flat;

    for (int p = 0; p < 3; ++p) {
        meta["nch"].append(cs.nch[p]);
        meta["base"].append(cs.base[p]);
    }

    // Per (apa,face) ticks-per-slice: the divisor that turns PR::Fit::pt into
    // the same slice index the 2-D cells are keyed on.
    for (const auto& [apa, face_map] : grouping.get_nticks_per_slice()) {
        for (const auto& [face, n] : face_map) {
            Configuration row;
            row["apa"] = apa;
            row["face"] = face;
            row["nticks_per_slice"] = n;
            meta["nticks_per_slice"].append(row);
        }
    }

    // Which (apa,face) pairs actually carry clustered charge this event.
    for (const auto& wpid : grouping.wpids()) {
        Configuration row;
        row["apa"] = wpid.apa();
        row["face"] = wpid.face();
        meta["apa_faces"].append(row);
    }

    return meta;
}

// Segments as POLYLINES plus the PR-graph vertices.
//
// DELTA vs both existing producers: they flatten every fit point into one long
// row list (a ROOT tree / a Bee point array), which loses the segment grouping
// the display needs to draw a line.  Here each segment keeps its own ordered
// point list.  Residual range is computed exactly as
// SbndPrMagnifyTrackingVisitor::write_t_rec_data does.
Configuration Clus::PrDisplayDump::dump_graph(Facade::Grouping& grouping) const
{
    Configuration out;
    out["segments"] = Json::arrayValue;
    out["vertices"] = Json::arrayValue;
    out["main_vertex"] = Json::nullValue;

    auto pr_graph = grouping.get_pr_graph();
    if (!pr_graph) {
        log->debug("no PR graph: TaggerCheckNeutrino selected no candidate this event");
        return out;
    }

    const double cm = units::cm;

    // Shower membership, same construction as dump_track_shower(): a segment
    // absorbed into a shower from another cluster may carry none of the
    // per-segment shower flags, so the shower's own segment set is the only
    // authoritative answer.
    std::map<PR::SegmentPtr, PR::ShowerPtr, PR::SegmentIndexCmp> seg_to_shower;
    if (auto tf = grouping.get_track_fitting()) {
        for (const auto& shower : tf->get_showers()) {
            PR::IndexedVertexSet sv;
            PR::IndexedSegmentSet ss;
            shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
            for (const auto& seg : ss) seg_to_shower[seg] = shower;
        }
    }

    auto fit_json = [&cm](const PR::Fit& fit) {
        Configuration j;
        j["x"] = fit.point.x() / cm;
        j["y"] = fit.point.y() / cm;
        j["z"] = fit.point.z() / cm;
        j["dQ"] = fit.dQ;
        j["dx"] = fit.dx / cm;
        // FRACTIONAL wire coordinates -- doc pr/7 sec 1: truncating to an int
        // wire puts the drawn track a systematic half channel off the charge.
        j["pu"] = fit.pu;
        j["pv"] = fit.pv;
        j["pw"] = fit.pw;
        j["pt"] = fit.pt;
        j["apa"] = fit.paf.first;
        j["face"] = fit.paf.second;
        j["reduced_chi2"] = fit.reduced_chi2;
        return j;
    };

    // ordered_nodes/ordered_edges give a stable order; raw boost::vertices /
    // boost::edges on a setS graph iterate in POINTER order and reshuffle the
    // output run to run (CLAUDE.md determinism rule).
    for (auto node_desc : PR::ordered_nodes(*pr_graph)) {
        auto vertex = (*pr_graph)[node_desc].vertex;
        if (!vertex) continue;

        const int cluster_id = vertex->cluster() ? vertex->cluster()->get_cluster_id() : 0;
        const bool is_main = vertex->flags_any(PR::VertexFlags::kNeutrinoVertex);

        Configuration j;
        j["cluster_id"] = cluster_id;
        j["id"] = cluster_id * 1000 + static_cast<int>(vertex->get_graph_index());
        j["is_main"] = is_main;
        // doc sbnd_xin/docs/pr/32 §11 F4 (was P12): the prototype's
        // map_cluster_main_candidate_vertices, i.e. the per-cluster list of
        // vertices determine_main_vertex scored between before filtering.
        // Always false unless TaggerCheckNeutrino ran with
        // main_vertex_candidate_flag=true, so a dump taken without that knob
        // reads "no candidate information", not "no candidates".
        j["main_candidate"] = vertex->flags_any(PR::VertexFlags::kMainCandidate);
        j["degree"] = static_cast<int>(boost::out_degree(node_desc, *pr_graph));
        // How far improve_vertex/MyFCN moved this vertex off its seed point.
        // A main vertex sitting at fit_distance 0 did not move -- one of the
        // three gates in doc pr/27 sec 12 fired.
        j["fit_distance"] = vertex->fit_distance() / cm;
        j["fit"] = fit_json(vertex->fit());
        out["vertices"].append(j);

        if (is_main && out["main_vertex"].isNull()) {
            Configuration mv;
            mv["x"] = vertex->fit().point.x() / cm;
            mv["y"] = vertex->fit().point.y() / cm;
            mv["z"] = vertex->fit().point.z() / cm;
            mv["cluster_id"] = cluster_id;
            out["main_vertex"] = mv;
        }
    }

    for (auto edge_desc : PR::ordered_edges(*pr_graph)) {
        auto segment = (*pr_graph)[edge_desc].segment;
        if (!segment) continue;

        const auto& fits = segment->fits();
        if (fits.empty()) continue;

        const int cluster_id = segment->cluster() ? segment->cluster()->get_cluster_id() : 0;
        const bool is_shower = segment->flags_any(PR::SegmentFlags::kShowerTrajectory) ||
                               segment->flags_any(PR::SegmentFlags::kShowerTopology);

        Configuration j;
        j["cluster_id"] = cluster_id;
        j["id"] = cluster_id * 1000 + static_cast<int>(segment->get_graph_index());
        j["flag_shower"] = is_shower;
        j["particle_id"] = segment->has_particle_info() ? segment->particle_info()->pdg()
                                                        : (is_shower ? 1 : 4);
        j["dirsign"] = segment->dirsign();
        j["is_main_cluster"] = segment->cluster()
            ? segment->cluster()->get_flag(Facade::Flags::main_cluster) : false;
        // PID/direction quality, already computed by segment_determine_dir_track
        // -- plain getters, no PID function invoked here (see header comment).
        j["particle_score"] = segment->particle_score();
        j["dir_weak"] = segment->dir_weak();
        // The particle-flow join key: -1 for a segment in no shower, otherwise
        // the mc.json node id of the shower this segment belongs to.  Clicking
        // that node in the display highlights every segment carrying this id.
        {
            auto sit = seg_to_shower.find(segment);
            j["shower_id"] = (sit == seg_to_shower.end())
                ? -1 : pf_node_id(sit->second->start_segment());
        }

        // Residual range, verbatim from write_t_rec_data: accumulate arclength,
        // orient by dirsign, then blank (-1) an end that meets a branching
        // vertex since "range to the end" is meaningless there.
        std::vector<double> rr(fits.size(), 0);
        {
            std::vector<double> L(fits.size(), 0);
            double acc = 0;
            for (size_t i = 0; i + 1 < fits.size(); ++i) {
                acc += (fits[i + 1].point - fits[i].point).magnitude();
                L[i + 1] = acc;
            }
            if (segment->dirsign() == 1) {
                for (size_t i = 0; i < fits.size(); ++i) {
                    rr[fits.size() - 1 - i] = L.back() - L[fits.size() - 1 - i];
                }
            }
            else {
                rr = L;
            }
            auto [start_vtx, end_vtx] = PR::find_vertices(*pr_graph, segment);
            if (start_vtx && boost::out_degree(start_vtx->get_descriptor(), *pr_graph) > 1) {
                rr.front() = -1;
            }
            if (end_vtx && boost::out_degree(end_vtx->get_descriptor(), *pr_graph) > 1) {
                rr.back() = -1;
            }
            // Total fitted arc length + the vertices[] ids of this segment's
            // two ends -- lets the viewer orient a start-mode (shower stem)
            // plot from the shower's own start vertex without recomputing
            // find_vertices() itself.  Also fixes the PF table's L (cm)
            // column, blank for track rows today (only showers[].total_length
            // populates it).
            j["length"] = L.back() / cm;
            j["start_vertex_id"] = vertex_display_id(start_vtx);
            j["end_vertex_id"] = vertex_display_id(end_vtx);
        }

        for (size_t i = 0; i < fits.size(); ++i) {
            Configuration p = fit_json(fits[i]);
            p["rr"] = rr[i] / cm;
            j["points"].append(p);
        }
        out["segments"].append(j);
    }

    return out;
}

// One row per reconstructed shower, keyed by the SAME id the Bee particle-flow
// tree (mc.json) puts on its shower nodes.
//
// DELTA vs fill_bee_pf_tree: that producer bakes the energy into a display
// STRING ("e-  1148 MeV") and drops everything else.  The display needs the
// numbers, so they are emitted here and the tree itself is left to mc.json --
// there is exactly one particle-flow producer and this is not it.
Configuration Clus::PrDisplayDump::dump_showers(Facade::Grouping& grouping) const
{
    Configuration out = Json::arrayValue;

    auto tf = grouping.get_track_fitting();
    if (!tf) return out;

    const double cm = units::cm;
    const double MeV = units::MeV;

    const auto& map_shower_pio_id = tf->get_map_shower_pio_id();
    const auto& map_pio_id_mass = tf->get_map_pio_id_mass();

    auto point_json = [&cm](const WireCell::Point& p) {
        Configuration j;
        j["x"] = p.x() / cm;
        j["y"] = p.y() / cm;
        j["z"] = p.z() / cm;
        return j;
    };

    // IndexedShowerSet is ordered by shower index, not by pointer.
    for (const auto& shower : tf->get_showers()) {
        if (!shower) continue;

        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();

        Configuration j;
        j["id"] = pf_node_id(shower->start_segment());
        j["shower_id"] = shower->get_shower_id();
        j["particle_id"] = shower->get_particle_type();
        j["kine_best"] = shower->get_kine_best() / MeV;
        j["kine_range"] = shower->get_kine_range() / MeV;
        j["kine_dQdx"] = shower->get_kine_dQdx() / MeV;
        j["kine_charge"] = shower->get_kine_charge() / MeV;
        j["flag_kinematics"] = shower->get_flag_kinematics();
        j["start"] = point_json(shower->get_start_point());
        j["end"] = point_json(shower->get_end_point());
        // 1 direct, 2 indirect, 3 gap, 4 not clearly connected (mc.json drops
        // type 4 entirely -- such a shower has a row here and no PF node).
        j["start_connection_type"] = conn_type;
        j["start_vertex_id"] = start_vtx
            ? (start_vtx->cluster() ? start_vtx->cluster()->get_cluster_id() * 1000 : 0)
              + static_cast<int>(start_vtx->get_graph_index())
            : -1;
        j["num_segments"] = shower->get_num_segments();
        j["num_main_segments"] = shower->get_num_main_segments();
        j["total_length"] = shower->get_total_length() / cm;

        // The literal stem samples (<=20, in MIP units) the nue/single-photon
        // taggers cut on -- Shower::get_stem_dQ_dx is read-only (walks
        // m_full_graph/fits()/wcpts(), appends only to a local vector; see
        // header comment).  m_mip_dqdx_median is a plain e/cm display
        // constant, so pass its internal-unit form to match the function's
        // own 43000/units::cm default convention.
        j["stem_dqdx"] = Json::arrayValue;
        if (start_vtx) {
            auto stem = shower->get_stem_dQ_dx(start_vtx, shower->start_segment(), 20,
                                                m_mip_dqdx_median / units::cm);
            for (double v : stem) j["stem_dqdx"].append(v);
        }

        auto pit = map_shower_pio_id.find(shower);
        if (pit == map_shower_pio_id.end()) {
            j["pio_id"] = -1;
            j["pio_mass"] = -1.0;
        }
        else {
            j["pio_id"] = pit->second;
            auto mit = map_pio_id_mass.find(pit->second);
            j["pio_mass"] = (mit == map_pio_id_mass.end()) ? -1.0
                                                           : mit->second.first / MeV;
        }
        out.append(j);
    }

    log->debug("showers: {} row(s)", out.size());
    return out;
}

// KineInfo verbatim.  Already in MeV and cm at the source
// (NeutrinoKinematics.cxx divides by units::MeV / units::cm before storing), so
// nothing is rescaled here.
Configuration Clus::PrDisplayDump::dump_kine(Facade::Grouping& grouping) const
{
    Configuration out;

    auto tf = grouping.get_track_fitting();
    if (!tf) return out;
    const auto& k = tf->get_kine_info();

    out["kine_reco_Enu"] = k.kine_reco_Enu;
    out["kine_reco_add_energy"] = k.kine_reco_add_energy;
    out["kine_nu_x_corr"] = k.kine_nu_x_corr;
    out["kine_nu_y_corr"] = k.kine_nu_y_corr;
    out["kine_nu_z_corr"] = k.kine_nu_z_corr;

    out["kine_energy_particle"] = Json::arrayValue;
    out["kine_energy_info"] = Json::arrayValue;
    out["kine_particle_type"] = Json::arrayValue;
    out["kine_energy_included"] = Json::arrayValue;
    for (float v : k.kine_energy_particle) out["kine_energy_particle"].append(v);
    for (int v : k.kine_energy_info) out["kine_energy_info"].append(v);
    for (int v : k.kine_particle_type) out["kine_particle_type"].append(v);
    for (int v : k.kine_energy_included) out["kine_energy_included"].append(v);

    out["kine_pio_mass"] = k.kine_pio_mass;
    out["kine_pio_flag"] = k.kine_pio_flag;
    out["kine_pio_vtx_dis"] = k.kine_pio_vtx_dis;
    out["kine_pio_energy_1"] = k.kine_pio_energy_1;
    out["kine_pio_theta_1"] = k.kine_pio_theta_1;
    out["kine_pio_phi_1"] = k.kine_pio_phi_1;
    out["kine_pio_dis_1"] = k.kine_pio_dis_1;
    out["kine_pio_energy_2"] = k.kine_pio_energy_2;
    out["kine_pio_theta_2"] = k.kine_pio_theta_2;
    out["kine_pio_phi_2"] = k.kine_pio_phi_2;
    out["kine_pio_dis_2"] = k.kine_pio_dis_2;
    out["kine_pio_angle"] = k.kine_pio_angle;

    return out;
}

// The selection numbers: the BDT scores, the cosmic tagger's verdict and its
// per-test decomposition.  The sub-BDT scores come along so a score that looks
// wrong can be decomposed on the spot.
//
// This stage is appended AFTER numu_bdt_scorer / nue_bdt_scorer in the PR
// pipeline (run_pr_chain_batch.sh), which is what makes these non-zero -- both
// scorers write through TrackFitting::get_tagger_info_mutable().
//
// ONLY FIELDS THAT ARE ACTUALLY COMPUTED ARE EMITTED.  TaggerInfo carries the
// uBooNE ntuple schema, a large part of which belongs to the legacy TMVA BDT
// path (prototype cal_bdts() / cal_numu_bdts()).  Neither has a caller: the
// prototype selects the xgboost variants at NeutrinoID.cxx:277,282
// (flag_bdt == 1) and only those are ported here.  A dumped field that no code
// ever assigns reads as a physics answer of 0 on the display, which is worse
// than its absence -- so the dead slots are left out and named in doc pr/26.
// doc sbnd_xin/docs/pr/75 -- HOW the neutrino vertex was chosen, not just
// where it landed.  The two selectors (doc pr/52 sec 1) compare numbers that
// exist nowhere in any artifact today: compare_main_vertices' additive score
// per candidate, and the DL rerank's top-K voxel scores + seven composite
// terms + accept decision.  A hand scan without them can only say
// correct/incorrect; with them it says WHICH TERM was wrong, which is what
// re-tuning the acceptance layer needs.
//
// An EMPTY object means the `vertex_scoreboard` knob was off -- read it as
// "no scoreboard was taken", NEVER as "no candidates existed".  Same
// discipline as vertices[].main_candidate above.
//
// Rows are emitted sorted by vertex_id: the sources they were collected from
// (map_vertex_num, snap_map) are both VertexPtr-keyed, i.e. ADDRESS-ordered,
// so an unsorted dump would reshuffle run to run (CLAUDE.md determinism rule).
Configuration Clus::PrDisplayDump::dump_vertex_scoreboard(Facade::Grouping& grouping) const
{
    Configuration out;

    auto tf = grouping.get_track_fitting();
    if (!tf) return out;
    const auto& b = tf->get_vertex_scoreboard();
    if (!b.filled) return out;

    out["filled"] = true;
    out["route"] = b.route;
    out["weights_missing"] = b.weights_missing;
    out["dl_ran"] = b.dl_ran;
    out["dl_rerank"] = b.dl_rerank;
    out["dl_accepted"] = b.dl_accepted;
    out["dl_best_score"] = b.dl_best_score;
    out["dl_min_accept_score"] = b.dl_min_accept_score;
    out["dl_score_scale"] = b.dl_score_scale;
    out["dl_top_k"] = b.dl_top_k;
    out["final_vertex_id"] = b.final_vertex_id;
    out["final_x"] = b.final_x;
    out["final_y"] = b.final_y;
    out["final_z"] = b.final_z;

    out["voxels"] = Json::arrayValue;
    for (const auto& v : b.voxels) {
        Configuration j;
        j["rank"] = v.rank;
        j["x"] = v.x;
        j["y"] = v.y;
        j["z"] = v.z;
        j["dl_score"] = v.dl_score;
        out["voxels"].append(j);
    }

    std::vector<const PR::VertexScoreRow*> rows;
    rows.reserve(b.rows.size());
    for (const auto& r : b.rows) rows.push_back(&r);
    std::sort(rows.begin(), rows.end(),
              [](const PR::VertexScoreRow* a, const PR::VertexScoreRow* c) {
                  return a->vertex_id < c->vertex_id;
              });

    out["rows"] = Json::arrayValue;
    for (const auto* r : rows) {
        Configuration j;
        // Joins to vertices[].id -- same cluster_id*1000 + graph_index encoding.
        j["vertex_id"] = r->vertex_id;
        j["cluster_id"] = r->cluster_id;
        j["x"] = r->x;
        j["y"] = r->y;
        j["z"] = r->z;
        j["trad_scored"] = r->trad_scored;
        j["trad_score"] = r->trad_score;
        j["trad_winner"] = r->trad_winner;
        // dl_snapped false means the DL had NO OPINION on this vertex (no
        // voxel snapped to it), which is not a zero score.
        j["dl_snapped"] = r->dl_snapped;
        j["voxel_rank"] = r->voxel_rank;
        j["dl_score"] = r->dl_score;
        j["snap_dis"] = r->snap_dis;
        j["host_length"] = r->host_length;
        j["s_dl"] = r->s_dl;
        j["s_snap"] = r->s_snap;
        j["s_fwd_z"] = r->s_fwd_z;
        j["s_clen"] = r->s_clen;
        j["s_isol"] = r->s_isol;
        j["s_main"] = r->s_main;
        j["s_fv"] = r->s_fv;
        j["total"] = r->total;
        j["dl_winner"] = r->dl_winner;
        // All-zero terms with this set mean "removed before scoring", not
        // "scored zero".
        j["skipped_by_swap_guard"] = r->skipped_by_swap_guard;
        out["rows"].append(j);
    }

    return out;
}

Configuration Clus::PrDisplayDump::dump_tagger(Facade::Grouping& grouping) const
{
    Configuration out;

    auto tf = grouping.get_track_fitting();
    if (!tf) return out;
    const auto& t = tf->get_tagger_info();

    // SBND books the uBooNE-TRAINED weight XMLs (sbnd/clus.jsonnet), so these
    // scores carry availability and relative ranking only.  The caveat travels
    // WITH the data so a display cannot show the number without it.
    out["weights"] = "uboone-trained -- UNCALIBRATED on SBND (doc pr/2 gap G1)";

    // ---- the two live discriminants ---------------------------------------
    out["nue_score"] = t.nue_score;
    out["numu_score"] = t.numu_score;
    out["match_isFC"] = t.match_isFC;

    // ---- the cosmic answer -------------------------------------------------
    // cosmict_flag is the cosmic tagger's VERDICT: the OR of its ten tests
    // (NeutrinoTaggerCosmic.cxx:1342-1347).  It is the field to read when
    // asking "did the cosmic tagger fire on this event".
    //
    // cosmic_flag is NOT a second verdict and NOT the same polarity.  It is a
    // BDT input feature written only inside the flag-9 block, and it is exactly
    // !cosmict_flag_9: 0 means flag 9 fired, 1 means either the flag-9 block
    // never ran (cosmic_filled == 0) or it ran and a neutrino-like shower at
    // the main vertex rescued the event (cosmic_filled == 1).  Its in-class
    // default is 1, so 1 on its own proves nothing -- hence cosmic_filled is
    // emitted beside it.  (Prototype: NeutrinoID_cosmic_tagger.h:781-783.)
    out["cosmict_flag"] = t.cosmict_flag;
    out["cosmic_flag"] = t.cosmic_flag;
    out["cosmic_filled"] = t.cosmic_filled;

    // Per-test decomposition -- which of the ten fired.  This is what makes a
    // cosmic tag actionable when tuning; the verdict alone does not say why.
    out["cosmict_flag_1"] = t.cosmict_flag_1;
    out["cosmict_flag_2"] = t.cosmict_flag_2;
    out["cosmict_flag_3"] = t.cosmict_flag_3;
    out["cosmict_flag_4"] = t.cosmict_flag_4;
    out["cosmict_flag_5"] = t.cosmict_flag_5;
    out["cosmict_flag_6"] = t.cosmict_flag_6;
    out["cosmict_flag_7"] = t.cosmict_flag_7;
    out["cosmict_flag_8"] = t.cosmict_flag_8;
    out["cosmict_flag_9"] = t.cosmict_flag_9;

    // The "did this test even run" companions.  Tests 2-8 each fill their
    // feature block only when their topology precondition is met (a muon
    // candidate at the main vertex, a stopped muon with a secondary, ...), so
    // a 0 flag with a 0 _filled means NOT EVALUATED, while a 0 flag with a 1
    // _filled means evaluated and not fired.  Without these the display cannot
    // tell an inactive tagger from a quiet one.  Tests 1 and 9 have no _filled
    // of their own (1 always evaluates; 9 uses cosmic_filled above) and test
    // 10's is the length of its vector.
    out["cosmict_2_filled"] = t.cosmict_2_filled;
    out["cosmict_3_filled"] = t.cosmict_3_filled;
    out["cosmict_4_filled"] = t.cosmict_4_filled;
    out["cosmict_5_filled"] = t.cosmict_5_filled;
    out["cosmict_6_filled"] = t.cosmict_6_filled;
    out["cosmict_7_filled"] = t.cosmict_7_filled;
    out["cosmict_8_filled"] = t.cosmict_8_filled;

    // Test 10 is per-(vertex,segment), not per-event: cosmict_flag_10 is a
    // vector with one entry per near-front-face candidate, and the tagger ORs
    // it into cosmict_flag through a local that is never stored.  Emit BOTH the
    // vector and the derived OR, because an EMPTY vector ("no candidate was
    // ever examined") is a different statement from an all-false one and only
    // the vector distinguishes them.
    {
        Configuration arr(Json::arrayValue);
        bool any10 = false;
        for (float f : t.cosmict_flag_10) {
            arr.append(f);
            if (f != 0) any10 = true;
        }
        out["cosmict_flag_10"] = arr;
        out["cosmict_flag_10_any"] = any10 ? 1.0f : 0.0f;
    }

    // ---- numu sub-scores (the four the xgboost path fills) -----------------
    // cosmict_2_4 / 3_5 / 6 / 7 / 8_score are legacy cal_numu_bdts() slots and
    // are never assigned; cosmict_score is never assigned in the prototype
    // either (initialised to 0 at NeutrinoID.cxx:3365 and left there).
    out["cosmict_10_score"] = t.cosmict_10_score;
    out["numu_1_score"] = t.numu_1_score;
    out["numu_2_score"] = t.numu_2_score;
    out["numu_3_score"] = t.numu_3_score;

    // ---- nue sub-scores (the fifteen the xgboost path fills) ---------------
    // UbooneNueBDTScorer.cxx:1631-1645.  The other fifteen uBooNE nue slots
    // (mipid/gap/hol_lol/cme_anc/mgo_mgt/br1/br3/stemdir_br2/trimuon/br4_tro/
    // mipquality/pio_1/stw_spt/vis_1/vis_2) belong to the legacy cal_bdts()
    // path (flag_bdt == 2), which is not ported.
    out["br3_3_score"] = t.br3_3_score;
    out["br3_5_score"] = t.br3_5_score;
    out["br3_6_score"] = t.br3_6_score;
    out["pio_2_score"] = t.pio_2_score;
    out["stw_2_score"] = t.stw_2_score;
    out["stw_3_score"] = t.stw_3_score;
    out["stw_4_score"] = t.stw_4_score;
    out["sig_1_score"] = t.sig_1_score;
    out["sig_2_score"] = t.sig_2_score;
    out["lol_1_score"] = t.lol_1_score;
    out["lol_2_score"] = t.lol_2_score;
    out["tro_1_score"] = t.tro_1_score;
    out["tro_2_score"] = t.tro_2_score;
    out["tro_4_score"] = t.tro_4_score;
    out["tro_5_score"] = t.tro_5_score;

    // photon_flag is deliberately NOT emitted: TaggerCheckNeutrino.cxx:813
    // discards singlephoton_tagger()'s return value, where the prototype does
    // `if (flag_sp) tagger_info.photon_flag = true;` (NeutrinoID.cxx:271).  The
    // field is therefore always 0 here regardless of the tagger's answer.
    // Reported in doc pr/26 sec 7.5, not fixed in this change.

    log->debug("tagger: nue_score {:.3f} numu_score {:.3f} cosmict_flag {:.0f} "
               "(flags 1-9 {}{}{}{}{}{}{}{}{} 10 {}) cosmic_flag {:.0f} filled {:.0f} "
               "Enu {:.1f} MeV",
               t.nue_score, t.numu_score, t.cosmict_flag,
               (int)t.cosmict_flag_1, (int)t.cosmict_flag_2, (int)t.cosmict_flag_3,
               (int)t.cosmict_flag_4, (int)t.cosmict_flag_5, (int)t.cosmict_flag_6,
               (int)t.cosmict_flag_7, (int)t.cosmict_flag_8, (int)t.cosmict_flag_9,
               out["cosmict_flag_10_any"].asFloat() != 0 ? 1 : 0,
               t.cosmic_flag, t.cosmic_filled,
               tf->get_kine_info().kine_reco_Enu);
    return out;
}

// The associated 3-D points with a per-point track/shower answer -- the same
// content as the Bee "shower_track" layer, which is what the owner sees today.
// Shower membership (not the per-segment flags) is authoritative: a segment
// absorbed into a shower from another cluster may never have had
// kShowerTrajectory/kShowerTopology or pdg=11 set on it.
Configuration Clus::PrDisplayDump::dump_track_shower(Facade::Grouping& grouping) const
{
    Configuration out;
    out["x"] = Json::arrayValue;
    out["y"] = Json::arrayValue;
    out["z"] = Json::arrayValue;
    out["flag_shower"] = Json::arrayValue;
    out["cluster_id"] = Json::arrayValue;
    out["particle_id"] = Json::arrayValue;

    auto pr_graph = grouping.get_pr_graph();
    auto tf = grouping.get_track_fitting();
    if (!pr_graph || !tf) return out;

    std::map<PR::SegmentPtr, PR::ShowerPtr, PR::SegmentIndexCmp> seg_to_shower;
    for (const auto& shower : tf->get_showers()) {
        PR::IndexedVertexSet sv;
        PR::IndexedSegmentSet ss;
        shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
        for (const auto& seg : ss) seg_to_shower[seg] = shower;
    }

    const double cm = units::cm;
    size_t npts = 0;
    for (auto edge_desc : PR::ordered_edges(*pr_graph)) {
        auto segment = (*pr_graph)[edge_desc].segment;
        if (!segment) continue;

        auto dpc = segment->dpcloud("associate_points");
        if (!dpc) continue;

        auto shower_it = seg_to_shower.find(segment);
        bool is_shower = (shower_it != seg_to_shower.end()) ||
            segment->flags_any(PR::SegmentFlags::kShowerTrajectory) ||
            segment->flags_any(PR::SegmentFlags::kShowerTopology) ||
            (segment->has_particle_info() && std::abs(segment->particle_info()->pdg()) == 11);
        // Long-muon pseudo-shower (cached type +-13): the PF tree shows a
        // muon; paint as track for consistency (doc sbnd_xin/docs/pr/45,
        // mirrors MultiAlgBlobClustering's gated override).
        if (m_pseudo_shower_track_paint && shower_it != seg_to_shower.end() &&
            std::abs(shower_it->second->get_particle_type()) == 13) {
            is_shower = false;
        }

        const int cluster_id = segment->cluster() ? segment->cluster()->get_cluster_id() : 0;
        // Per-particle id, the prototype convention MABC's `particle_ids` knob
        // turns on: a segment inside a shower inherits the shower's start
        // segment, otherwise it is its own particle.
        const int particle_id = [&]() -> int {
            if (shower_it == seg_to_shower.end()) {
                int sid = segment->id();
                if (sid < 0) sid = static_cast<int>(segment->get_graph_index());
                return cluster_id * 1000 + sid;
            }
            auto start_seg = shower_it->second->start_segment();
            if (!start_seg) return cluster_id;
            int sid = start_seg->id();
            if (sid < 0) sid = static_cast<int>(start_seg->get_graph_index());
            const auto* cl = start_seg->cluster();
            return cl ? cl->get_cluster_id() * 1000 + sid : sid;
        }();

        const size_t n = dpc->npoints();
        for (size_t ip = 0; ip < n; ++ip) {
            const auto p = dpc->point3d(ip);
            out["x"].append(p.x() / cm);
            out["y"].append(p.y() / cm);
            out["z"].append(p.z() / cm);
            out["flag_shower"].append(is_shower ? 1 : 0);
            out["cluster_id"].append(cluster_id);
            out["particle_id"].append(particle_id);
        }
        npts += n;
    }

    log->debug("track_shower: {} associated point(s)", npts);
    return out;
}

// The Steiner skeleton, per cluster, with the terminal flag that is ALREADY a
// persisted column of steiner_pc (SteinerGrapher.cxx step 9) -- so the display
// gets both the skeleton and the terminal subset with no clustering-stage
// change and no pctree schema change.
Configuration Clus::PrDisplayDump::dump_steiner(Facade::Grouping& grouping) const
{
    Configuration out = Json::arrayValue;

    // Cluster-id order, not child order: the emitted order must not depend on
    // anything pointer-like.
    std::map<int, const Facade::Cluster*> by_id;
    for (const auto* cluster : grouping.children()) {
        if (!cluster) continue;
        by_id[cluster->get_cluster_id()] = cluster;
    }

    const double cm = units::cm;
    size_t nterm = 0, npts = 0;
    for (const auto& [cid, cluster] : by_id) {
        if (!cluster->has_pc("steiner_pc")) continue;
        const auto& pc = cluster->get_pc("steiner_pc");

        const auto& coords = cluster->get_default_scope().coords;
        const auto xa = pc.get(coords.at(0));
        const auto ya = pc.get(coords.at(1));
        const auto za = pc.get(coords.at(2));
        if (!xa || !ya || !za) continue;   // an empty steiner_pc has no arrays at all
        const auto fa = pc.get("flag_steiner_terminal");

        const auto& xs = xa->elements<double>();
        const auto& ys = ya->elements<double>();
        const auto& zs = za->elements<double>();

        Configuration j;
        j["cluster_id"] = cid;
        j["is_main_cluster"] = cluster->get_flag(Facade::Flags::main_cluster);
        for (size_t i = 0; i < xs.size(); ++i) {
            j["x"].append(xs[i] / cm);
            j["y"].append(ys[i] / cm);
            j["z"].append(zs[i] / cm);
        }
        if (fa) {
            const auto& fs = fa->elements<int>();
            for (size_t i = 0; i < fs.size(); ++i) {
                j["flag_terminal"].append(fs[i]);
                if (fs[i]) ++nterm;
            }
        }
        npts += xs.size();
        out.append(j);
    }

    log->debug("steiner: {} cluster(s), {} point(s), {} terminal(s)", out.size(), npts, nterm);
    return out;
}

// The 2-D fitted charge, grouped by (apa, face, plane) so the viewer's six
// panels are a direct read.
//
// DELTA vs write_proj_data, which groups by cluster id and emits GLOBAL channel
// numbers: the panels are per-TPC-per-plane, so grouping that way and keeping
// the per-APA wire index is what the display actually wants.  The
// determinism-critical part is kept: accumulate into an ORDERED map and take
// the cluster ids of a cell from a sorted vector, never from the iteration
// order of fc.clusters (a std::set<Cluster*>).
//
// KNOWN NON-DETERMINISM, INHERITED, NOT INTRODUCED HERE: `charge_pred` is
// run-dependent on cells claimed by more than one cluster -- 1379 of 13507
// cells (10.2%) changed between two runs of SBND evt 18255/388 under
// `setarch -R`, while wire/slice/charge/charge_err/cluster_id/flag were all
// bit-identical, as was every other section of this dump.
//
// Cause: TrackFitting::assemble_fitted_charge_2d (clus/src/TrackFitting.cxx:1139)
// merges the per-cluster snapshots last-writer-wins while iterating
// m_cluster_fitted_charge_2d, a std::map<Facade::Cluster*, ...> -- a
// pointer-keyed container, whose order is not reproducible.  charge and
// charge_err depend only on the readout so the overwrite is benign for them;
// pred_charge is per-cluster, so the winner varies.
//
// Blast radius is diagnostic-only: the merged map is read by this dumper, the
// three Magnify tracking writers and TaggerCheckSTM's stm_fit record.  No
// tagger verdict, Bee layer or pctree tensor depends on it, which is why no
// A/B gate has ever caught it.  Reported in sbnd_xin/docs/pr/26, NOT fixed
// here (fixing it belongs in TrackFitting, with its own gate).  Until then,
// do not read the display's per-cell measured-vs-predicted comparison as a
// stable number.
Configuration Clus::PrDisplayDump::dump_proj(Facade::Grouping& grouping) const
{
    Configuration out = Json::arrayValue;

    auto tf = grouping.get_track_fitting();
    if (!tf) return out;

    const auto& fitted = tf->get_fitted_charge_2d();
    if (fitted.empty()) {
        log->warn("fitted_charge_2d is empty -- the 2-D panels will be blank");
        return out;
    }

    auto nticks_map = grouping.get_nticks_per_slice();

    size_t ncell = 0;
    for (const auto& [afp, wt_map] : fitted) {          // std::map -> ordered
        const int apa = std::get<0>(afp);
        const int face = std::get<1>(afp);
        const int plane = std::get<2>(afp);

        int nticks_per_slice = 1;
        auto ait = nticks_map.find(apa);
        if (ait != nticks_map.end()) {
            auto fit_ = ait->second.find(face);
            if (fit_ != ait->second.end()) nticks_per_slice = fit_->second;
        }

        Configuration j;
        j["apa"] = apa;
        j["face"] = face;
        j["plane"] = plane;
        j["nticks_per_slice"] = nticks_per_slice;
        j["wire"] = Json::arrayValue;
        j["slice"] = Json::arrayValue;
        j["charge"] = Json::arrayValue;
        j["charge_err"] = Json::arrayValue;
        j["charge_pred"] = Json::arrayValue;
        j["flag"] = Json::arrayValue;
        j["cluster_id"] = Json::arrayValue;

        for (const auto& [wt, fc] : wt_map) {           // std::map -> ordered
            if (fc.charge <= m_proj_charge_min) continue;

            // Lowest cluster id owning the cell.  A cell can be shared; the
            // display only needs "which cluster does this belong to" for
            // highlighting, and taking the sorted-first id makes that answer
            // independent of the pointer-set iteration order.
            std::vector<int> cids;
            cids.reserve(fc.clusters.size());
            for (auto* cl : fc.clusters) {
                if (cl) cids.push_back(cl->get_cluster_id());
            }
            std::sort(cids.begin(), cids.end());
            const int cid = cids.empty() ? -1 : cids.front();

            j["wire"].append(wt.first);
            j["slice"].append(wt.second / nticks_per_slice);
            j["charge"].append(fc.charge);
            j["charge_err"].append(fc.charge_err);
            j["charge_pred"].append(fc.pred_charge);
            j["flag"].append(fc.flag);
            j["cluster_id"].append(cid);
            ++ncell;
        }
        out.append(j);
    }

    log->debug("proj: {} (apa,face,plane) group(s), {} cell(s)", out.size(), ncell);
    return out;
}

// Dead channel ranges, so the 2-D panels can shade what was never measured.
// Same clamping as SbndPrMagnifyTrackingVisitor: an unbounded dead x extent
// converts to +-1e9 ticks and makes the pads unreadable.
Configuration Clus::PrDisplayDump::dump_dead(Facade::Grouping& grouping) const
{
    Configuration out = Json::arrayValue;

    std::set<std::pair<int, int>> apa_face_set;
    for (const auto& wpid : grouping.wpids()) {
        apa_face_set.insert({wpid.apa(), wpid.face()});
    }

    auto nticks_map = grouping.get_nticks_per_slice();

    for (const auto& [apa, face] : apa_face_set) {
        // get_all_dead_chs speaks TICKS while the 2-D cells above are keyed on
        // SLICES.  The Magnify writer leaves that mismatch to its reader (its
        // T_bad_ch is ticks, its T_proj_data is slices); emit both here so the
        // display cannot get it wrong.
        int nticks_per_slice = 1;
        auto ait = nticks_map.find(apa);
        if (ait != nticks_map.end()) {
            auto fit_ = ait->second.find(face);
            if (fit_ != ait->second.end()) nticks_per_slice = fit_->second;
        }

        for (int pind = 0; pind < 3; ++pind) {
            try {
                for (const auto& [ch, tr] : grouping.get_all_dead_chs(apa, face, pind)) {
                    int t0 = std::max(0, std::min(m_nticks, tr.first));
                    int t1 = std::max(0, std::min(m_nticks, tr.second));
                    if (t1 <= t0) { t0 = 0; t1 = m_nticks; }
                    Configuration j;
                    j["apa"] = apa;
                    j["face"] = face;
                    j["plane"] = pind;
                    j["wire"] = ch;
                    j["t0"] = t0;                        // ticks
                    j["t1"] = t1;
                    j["s0"] = t0 / nticks_per_slice;     // slices, = proj.slice units
                    j["s1"] = t1 / nticks_per_slice;
                    out.append(j);
                }
            }
            catch (const std::exception& e) {
                log->warn("dead channels unavailable for apa={} face={} plane={}: {}",
                          apa, face, pind, e.what());
            }
        }
    }
    return out;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
