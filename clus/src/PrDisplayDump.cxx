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

    for (auto anode_tn : cfg["anodes"]) {
        m_anodes.push_back(Factory::find_tn<IAnodePlane>(anode_tn.asString()));
    }
    if (cfg.isMember("detector_volumes") && !cfg["detector_volumes"].asString().empty()) {
        m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());
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

    top["track_shower"] = dump_track_shower(grouping);
    top["steiner"] = dump_steiner(grouping);
    top["proj"] = dump_proj(grouping);
    top["dead"] = dump_dead(grouping);

    Persist::dump(m_output_filename, top, m_pretty);

    log->debug("wrote {}: {} segment(s), {} vertex(es), {} steiner cluster(s), {} proj plane(s)",
               m_output_filename, top["segments"].size(), top["vertices"].size(),
               top["steiner"].size(), top["proj"].size());
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
        j["degree"] = static_cast<int>(boost::out_degree(node_desc, *pr_graph));
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
        const bool is_shower = (shower_it != seg_to_shower.end()) ||
            segment->flags_any(PR::SegmentFlags::kShowerTrajectory) ||
            segment->flags_any(PR::SegmentFlags::kShowerTopology) ||
            (segment->has_particle_info() && std::abs(segment->particle_info()->pdg()) == 11);

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
