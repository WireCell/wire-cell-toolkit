// This reimplements the ideas from:
// https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ChargeSolving.cxx
// and
// https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/MatrixSolver.cxx


#include "WireCellImg/ChargeSolving.h"
#include "WireCellImg/CSGraph.h"

#include "WireCellAux/SimpleCluster.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

#include <iterator>

WIRECELL_FACTORY(ChargeSolving, WireCell::Img::ChargeSolving,
                 WireCell::INamed,
                 WireCell::IClusterFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Img::CS;

Img::ChargeSolving::ChargeSolving()
    : Aux::Logger("ChargeSolving", "img")
{
}

Img::ChargeSolving::~ChargeSolving() {}

WireCell::Configuration Img::ChargeSolving::default_configuration() const
{
    WireCell::Configuration cfg;
    cfg["meas_value_threshold"] = m_meas_thresh.value();
    cfg["meas_error_threshold"] = m_meas_thresh.uncertainty();
    cfg["blob_value_threshold"] = m_blob_thresh.value();
    cfg["blob_error_threshold"] = m_blob_thresh.uncertainty();

    cfg["lasso_tolerance"] = m_lasso_tolerance;
    cfg["lasso_minnorm"] = m_lasso_minnorm;
    for (const auto& one : m_weighting_strategies) {
        cfg["weighting_strategies"].append(one);
    }
    cfg["whiten"] = m_whiten;

    return cfg;
}



void blob_weight_uniform(const cluster_graph_t& /*cgraph*/, graph_t& csg)
{
    for (auto desc : vertex_range(csg)) {
        auto& vtx = csg[desc];
        if (vtx.kind != node_t::blob) {
            continue;
        }
        vtx.value.uncertainty(1.0);
    }
}

void blob_weight_simple(const cluster_graph_t& cgraph, graph_t& csg)
{
    for (auto desc : vertex_range(csg)) {
        auto& vtx = csg[desc];
        if (vtx.kind != node_t::blob) {
            continue;
        }
        // get neighbors in original cluster graph
        // count how many unique slices have connected blobs
        // add 1.0.  that's the simple weight.
        std::unordered_set<int> slice_idents;
        slice_idents.insert(csg[boost::graph_bundle].islice->ident());
        for (auto edge : mir(boost::out_edges(vtx.orig_desc, cgraph)))
        {
            vdesc_t ndesc = boost::target(edge, cgraph);
            const auto& nnode = cgraph[ndesc];
            if (nnode.code() == 'b') {
                slice_idents.insert(get<blob_t>(nnode.ptr)->slice()->ident());
            }
        }
        vtx.value.uncertainty((float) slice_idents.size());
    }
}
// fixme: implement distance.

// Weighting function lookup
using blob_weighting_f = std::function<void(const cluster_graph_t& cgraph, graph_t& csg)>;
using blob_weighting_lut = std::unordered_map<std::string, blob_weighting_f>;
static blob_weighting_lut gStrategies{
    {"uniform", blob_weight_uniform},
    {"simple", blob_weight_simple},
    // distance....
};

// this depends on gStrategies.
void Img::ChargeSolving::configure(const WireCell::Configuration& cfg)
{
    float mt_val = get(cfg,"meas_value_threshold", m_meas_thresh.value());
    float mt_err = get(cfg,"meas_error_threshold", m_meas_thresh.uncertainty());
    m_meas_thresh = value_t(mt_val, mt_err);

    float bt_val = get(cfg,"blob_value_threshold", m_blob_thresh.value());
    float bt_err = get(cfg,"blob_error_threshold", m_blob_thresh.uncertainty());
    m_blob_thresh = value_t(bt_val, bt_err);

    m_lasso_tolerance = get(cfg, "lasso_tolerance", m_lasso_tolerance);
    m_lasso_minnorm = get(cfg, "lasso_minnorm", m_lasso_minnorm);

    m_weighting_strategies =
        get<std::vector<std::string>>(cfg, "weighting_strategies",
                                      m_weighting_strategies);
    for (const auto& strat : m_weighting_strategies) {
        if (gStrategies.find(strat) == gStrategies.end()) {
            THROW(ValueError() << errmsg{"Unknown weighting strategy: " + strat});
        }
    }
    m_whiten = get<bool>(cfg, "whiten", m_whiten);
}


static void dump_cg(const cluster_graph_t& cg, Log::logptr_t& log)
{
    size_t mcount{0}, bcount{0};
    value_t bval;
    for (const auto& vtx : mir(boost::vertices(cg))) {
        const auto& node = cg[vtx];
        if (node.code() == 'b') {
            const auto iblob = get<blob_t>(node.ptr);
            bval += value_t(iblob->value(), iblob->uncertainty());
            ++bcount;
            continue;
        }
        if (node.code() == 'm') {
            const auto imeas = get<meas_t>(node.ptr);
            ++mcount;
        }
    }
    log->debug("cluster graph: vertices={} edges={} $blob={} bval={} #meas={}",
               boost::num_vertices(cg), boost::num_edges(cg),
               bcount, bval, mcount);
}

static void dump_sg(const graph_t& sg, Log::logptr_t& log)
{
    value_t bval, mval;
    const auto& gval = sg[boost::graph_bundle];

    int nblob{0}, nmeas{0};

    for (const auto& vtx : mir(boost::vertices(sg))) {
        const auto& node = sg[vtx];
        
        if (node.kind == node_t::blob) {
            bval += node.value;
            ++nblob;
            continue;
        }
        if (node.kind == node_t::meas) {
            mval += node.value;
            ++nmeas;
            continue;
        }
    }

    log->trace("cs graph: slice={} index={} vertices={} edges={} nblob={} bval={} nmeas={} mval={}",
               gval.islice->ident(), gval.index,
               boost::num_vertices(sg), boost::num_edges(sg),
               nblob, bval, nmeas, mval);
}

value_t total_value(const graph_t& sg, node_t::Kind kind)
{
    value_t ret;
    for (auto desc : vertex_range(sg)) {
        const auto& node = sg[desc];
        if (node.kind == kind) {
            ret += node.value;
        }
    }
    return ret;
}


bool Img::ChargeSolving::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    // Separate the big graph spanning the whole frame into connected
    // b-m subgraphs with all the info needed for solving each round.
    graph_vector_t sgs;
    const auto in_graph = in->graph();
    dump_cg(in_graph, log);
    unpack(in_graph, std::back_inserter(sgs), m_meas_thresh);

    // for (const auto& sg : sgs) {
    //     dump_sg(sg, log);
    // }

    const size_t nstrats = m_weighting_strategies.size();

    std::vector<float> blob_threshold(nstrats, m_blob_thresh.value());

    SolveParams sparams{Ress::Params{Ress::lasso}, 1000, m_whiten};
    for (size_t ind = 0; ind < nstrats; ++ind) {
        const auto& strategy = m_weighting_strategies[ind];
        log->debug("cluster: {} strategy={}",
                   in->ident(), strategy);

        auto& blob_weighter = gStrategies[strategy];

        std::transform(sgs.begin(), sgs.end(), sgs.begin(),
                       [&](graph_t& sg) {
                           //dump_sg(sg, log);
                           blob_weighter(in_graph, sg);
                           auto tmp_csg = solve(sg, sparams);
                           return prune(tmp_csg, blob_threshold[ind]);
                       });
    }
    for (const auto& sg : sgs) {
        dump_sg(sg, log);
    }
    log->debug("count={} cluster={} solved nsubclusters={} over nstrategies={}",
               m_count,
               in->ident(), sgs.size(), m_weighting_strategies.size());

    auto packed = repack(in_graph, sgs);

    log->debug("call={} cluster={} nvertices={} nedges={}",
               m_count,
               in->ident(),
               boost::num_vertices(packed),
               boost::num_edges(packed));

    out = std::make_shared<Aux::SimpleCluster>(packed, in->ident());
    ++m_count;
    return true;
}


