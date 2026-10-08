// CascadeDeghosting (wcfm doc 14).  See WireCellImg/CascadeDeghosting.h.

#include "WireCellImg/CascadeDeghosting.h"
#include "WireCellImg/CascadeGraph.h"
#include "WireCellImg/CellSteiner.h"
#include "WireCellImg/GeomClusteringUtil.h"

#include "WireCellAux/CascadeRun.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

WIRECELL_FACTORY(CascadeDeghosting, WireCell::Img::CascadeDeghosting,
                 WireCell::INamed,
                 WireCell::Img::IClusterFrameJoin, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;

Img::CascadeDeghosting::CascadeDeghosting()
  : Aux::Logger("CascadeDeghosting", "img")
{
}
Img::CascadeDeghosting::~CascadeDeghosting() {}

WireCell::Configuration Img::CascadeDeghosting::default_configuration() const
{
    Configuration cfg;
    cfg["levels"] = Json::arrayValue;
    cfg["cut_max_depth"] = m_cut_max_depth;
    cfg["cut_min_length"] = m_cut_min_length;
    cfg["cut_nudge"] = m_cut_nudge;
    cfg["guard"] = m_guard;
    cfg["max_level_nodes"] = m_max_level_nodes;
    cfg["repair"] = m_repair;
    cfg["repair_p_term"] = m_repair_p_term;
    cfg["repair_q_floor"] = m_repair_q_floor;
    cfg["repair_budget"] = m_repair_budget;
    cfg["policy"] = m_policy;
    cfg["charge_tag"] = m_charge_tag;
    cfg["charge_scale"] = m_charge_scale;
    cfg["slice_start_relative"] = m_slice_start_relative;
    cfg["uncer_cut"] = m_uncer_cut;
    cfg["ident_base"] = m_ident_base;
    cfg["dump_dir"] = m_dump_dir;
    cfg["nthreads"] = m_nthreads;
    cfg["iso_fallback"] = Json::nullValue;   // off; see the header for the members
    cfg["final_guard"] = m_final_guard;
    cfg["keep_slices"] = m_keep_slices;
    return cfg;
}

void Img::CascadeDeghosting::configure(const WireCell::Configuration& cfg)
{
    m_cut_max_depth = get(cfg, "cut_max_depth", m_cut_max_depth);
    m_cut_min_length = get(cfg, "cut_min_length", m_cut_min_length);
    m_cut_nudge = get(cfg, "cut_nudge", m_cut_nudge);
    m_guard = get(cfg, "guard", m_guard);
    m_max_level_nodes = get(cfg, "max_level_nodes", m_max_level_nodes);
    m_repair = get(cfg, "repair", m_repair);
    m_repair_p_term = get(cfg, "repair_p_term", m_repair_p_term);
    m_repair_q_floor = get(cfg, "repair_q_floor", m_repair_q_floor);
    m_repair_budget = get(cfg, "repair_budget", m_repair_budget);
    m_policy = get(cfg, "policy", m_policy);
    m_charge_tag = get(cfg, "charge_tag", m_charge_tag);
    m_charge_scale = get(cfg, "charge_scale", m_charge_scale);
    m_slice_start_relative = get(cfg, "slice_start_relative", m_slice_start_relative);
    m_uncer_cut = get(cfg, "uncer_cut", m_uncer_cut);
    m_ident_base = get(cfg, "ident_base", m_ident_base);
    m_dump_dir = get(cfg, "dump_dir", m_dump_dir);
    m_nthreads = std::max(1, get(cfg, "nthreads", m_nthreads));
    m_final_guard = get(cfg, "final_guard", m_final_guard);
    m_keep_slices = get(cfg, "keep_slices", m_keep_slices);
    if (m_final_guard || m_keep_slices) {
        log->debug("wcfm doc 22 knobs: final_guard={} keep_slices={}", m_final_guard, m_keep_slices);
    }
    const auto& jiso = cfg["iso_fallback"];
    m_iso = jiso.isObject();
    if (m_iso) {
        m_iso_nmin = get(jiso, "nmin", m_iso_nmin);
        m_iso_mmin = get(jiso, "mmin", m_iso_mmin);
        m_iso_amin = get(jiso, "amin", m_iso_amin);
        m_iso_t = get(jiso, "t_keep", m_iso_t);
        m_iso_amb_lo = get(jiso, "amb_lo", m_iso_amb_lo);
        m_iso_amb_hi = get(jiso, "amb_hi", m_iso_amb_hi);
        log->debug("iso_fallback on: nmin={} mmin={} amin={} t_keep={} amb=({}, {})", m_iso_nmin, m_iso_mmin, m_iso_amin,
                   m_iso_t, m_iso_amb_lo, m_iso_amb_hi);
    }

    m_levels.clear();
    for (const auto& jl : cfg["levels"]) {
        LevelCfg lc;
        lc.width = get<int>(jl, "width", 0);
        lc.forward_tn = get<std::string>(jl, "forward", "");
        lc.superwire = get<int>(jl, "superwire", 1);
        lc.threshold = get<double>(jl, "threshold", 0.0);
        if (lc.forward_tn.empty()) {
            THROW(ValueError() << errmsg{"CascadeDeghosting: every level needs a forward"});
        }
        if (m_levels.empty() != (lc.width <= 0)) {
            THROW(ValueError() << errmsg{"CascadeDeghosting: level 0 is uncut (width 0), every later level has a width"});
        }
        lc.forward = Factory::find_tn<ITensorForward>(lc.forward_tn);
        m_levels.push_back(lc);
    }
    if (m_levels.empty()) {
        THROW(ValueError() << errmsg{"CascadeDeghosting: no levels"});
    }
    if (m_repair_budget <= 0) {
        const double t = m_levels.back().threshold;
        m_repair_budget = Cascade::default_repair_budget(t);
    }
    std::string desc;
    for (const auto& lc : m_levels) {
        desc += String::format(" [width=%d k=%d thr=%.4f %s]", lc.width, lc.superwire, lc.threshold,
                               lc.forward_tn.c_str());
    }
    log->debug("levels:{} cut max_depth={} min_length={} guard={} max_level_nodes={} repair={} budget={:.4f} policy={} charge_scale={} nthreads={}", desc,
               m_cut_max_depth, m_cut_min_length, m_guard, m_max_level_nodes, m_repair, m_repair_budget, m_policy, m_charge_scale,
               m_nthreads);
}

namespace {
    // BlobClustering.cxx add_slice / add_blobs, duplicated (fork rule: the production file is untouched)
    void add_slice(cluster_indexed_graph_t& grind, const ISlice::pointer& islice)
    {
        if (grind.has(islice)) {
            return;
        }
        for (const auto& ichv : islice->activity()) {
            const IChannel::pointer ich = ichv.first;
            if (grind.has(ich)) {
                continue;
            }
            for (const auto& iwire : ich->wires()) {
                grind.edge(ich, iwire);
            }
        }
    }
    void add_blobs(cluster_indexed_graph_t& grind, const IBlob::vector& iblobs)
    {
        for (const auto& iblob : iblobs) {
            auto islice = iblob->slice();
            add_slice(grind, islice);
            grind.edge(islice, iblob);
            auto iface = iblob->face();
            auto wire_planes = iface->planes();
            const auto& shape = iblob->shape();
            for (const auto& strip : shape.strips()) {
                const int num_nonplane_layers = 2;
                int iplane = strip.layer - num_nonplane_layers;
                if (iplane < 0) {
                    continue;
                }
                const auto& wires = wire_planes[iplane]->wires();
                for (int wip = strip.bounds.first; wip < strip.bounds.second and wip < int(wires.size()); ++wip) {
                    grind.edge(iblob, wires[wip]);
                }
            }
        }
    }
}  // namespace

bool Img::CascadeDeghosting::operator()(const input_tuple_type& intup, output_pointer& out)
{
    out = nullptr;
    const auto& in = std::get<0>(intup);
    const auto& frame = std::get<1>(intup);
    if (!in) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    const auto& gr = in->graph();

    bool use_frame = false;
    if (frame) {
        const size_t ntr = m_charge_tag.empty() ? frame->traces()->size() : frame->tagged_traces(m_charge_tag).size();
        use_frame = ntr > 0;
        if (!use_frame) log->warn("call={} frame {} has no trace tagged \"{}\": charge from the slice activity", m_count,
                                  frame->ident(), m_charge_tag);
    }
    else {
        log->warn("call={} no frame: charge from the slice activity (not the training charge)", m_count);
    }
    auto sc = use_frame ? Cascade::make_slice_charge_frame(gr, frame, m_charge_tag, m_charge_scale, m_slice_start_relative)
                        : Cascade::make_slice_charge(gr, m_charge_scale, m_uncer_cut);
    std::vector<IBlob::pointer> cur;
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        if (gr[vtx].code() == 'b') cur.push_back(std::get<IBlob::pointer>(gr[vtx].ptr));
    }
    // ---- the level loop: Aux::Cascade::run_cascade (pdvd doc 130 phase 2; the body moved as it was)
    Cascade::RunParams rp;
    rp.cut.nudge = m_cut_nudge;
    rp.cut.max_depth = m_cut_max_depth;
    rp.cut.min_length = m_cut_min_length;
    rp.guard = m_guard;
    rp.max_level_nodes = m_max_level_nodes;
    rp.repair = m_repair;
    rp.steiner.p_term = m_repair_p_term;
    rp.steiner.q_floor = m_repair_q_floor;
    rp.steiner.budget = m_repair_budget;
    rp.policy = m_policy;
    rp.nthreads = m_nthreads;
    rp.iso = m_iso;
    rp.iso_params.nmin = (int) m_iso_nmin;
    rp.iso_params.mmin = m_iso_mmin;
    rp.iso_params.amin = m_iso_amin;
    rp.iso_params.t_keep = m_iso_t;
    rp.iso_params.amb_lo = m_iso_amb_lo;
    rp.iso_params.amb_hi = m_iso_amb_hi;
    rp.final_guard = m_final_guard;
    rp.dump_dir = m_dump_dir;
    rp.ident_base = m_ident_base;
    std::vector<Cascade::LevelSpec> levels;
    for (const auto& lc : m_levels) {
        levels.push_back({lc.width, lc.superwire, lc.threshold, lc.forward, lc.forward_tn});
    }
    const auto res = Cascade::run_cascade(cur, sc, levels, rp, in->ident(), log, m_count);
    const auto& kept = res.kept;
    const size_t nin = res.nin;
    const size_t nkeep_thr = res.nkeep_thr;
    const auto& srep = res.srep;

    // ---- the output cluster: kept blobs, renumbered, built as BlobClustering builds one
    IBlob::vector outblobs;
    outblobs.reserve(kept.size());
    int ident = m_ident_base;
    for (const auto& b : kept) {
        outblobs.push_back(std::make_shared<Aux::SimpleBlob>(ident++, b->value(), b->uncertainty(), b->shape(),
                                                             b->slice(), b->face()));
    }
    std::vector<IBlob::vector> per(sc.slice_of.size());
    for (const auto& b : outblobs) per[sc.index(b->slice())].push_back(b);
    IBlobSet::vector sets;
    for (size_t s = 0; s < per.size(); ++s) {
        if (per[s].empty()) continue;
        sets.push_back(std::make_shared<Aux::SimpleBlobSet>((int) s, sc.slice_of[s], per[s]));
    }
    cluster_indexed_graph_t grind;
    for (auto it = sets.begin(); it != sets.end(); ++it) {
        add_blobs(grind, (*it)->blobs());
        Img::geom_clustering(grind, it, sets.end(), m_policy);
    }
    if (m_keep_slices) {   // wcfm doc 22: every input slice node of a time whose blobs were all dropped is kept,
                           // with its activity (PointTreeBuilding's ctpc reads every s-node), in input vertex order
        size_t nkept_slices = 0;
        for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
            if (gr[vtx].code() != 's') continue;
            const auto& islice = std::get<ISlice::pointer>(gr[vtx].ptr);
            if (!per[sc.index(islice)].empty() || grind.has(islice)) continue;
            add_slice(grind, islice);
            grind.vertex(islice);
            ++nkept_slices;
        }
        log->debug("call={} cluster={} keep_slices: {} blob-less slice nodes kept", m_count, in->ident(), nkept_slices);
    }
    out = std::make_shared<Aux::SimpleCluster>(std::move(grind.graph()), in->ident());
    const double dt = std::chrono::duration<double>(clock::now() - t0).count();
    log->debug("call={} cluster={} blobs in={} kept at threshold={} weak dropped={} bridges={} added={} out={} "
               "vertices={} edges={} t={:.2f}s",
               m_count, in->ident(), nin, nkeep_thr, srep.nweak_cells, srep.nbridges, srep.nadded, outblobs.size(),
               boost::num_vertices(out->graph()), boost::num_edges(out->graph()), dt);
    ++m_count;
    return true;
}
