
#include "WireCellImg/GlobalGeomClustering.h"
#include "WireCellImg/CSGraph.h"
#include "WireCellImg/GeomClusteringUtil.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"

#include <iterator>
#include <chrono>
#include <fstream>

WIRECELL_FACTORY(GlobalGeomClustering, WireCell::Img::GlobalGeomClustering, WireCell::INamed, WireCell::IClusterFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;

Img::GlobalGeomClustering::GlobalGeomClustering()
  : Aux::Logger("GlobalGeomClustering", "img")
{
}

Img::GlobalGeomClustering::~GlobalGeomClustering() {}

WireCell::Configuration Img::GlobalGeomClustering::default_configuration() const
{
    WireCell::Configuration cfg;
    cfg["clustering_policy"] = m_clustering_policy;
    return cfg;
}

void Img::GlobalGeomClustering::configure(const WireCell::Configuration& cfg)
{
    log->debug("{}", cfg);
    m_clustering_policy = get<std::string>(cfg, "clustering_policy", m_clustering_policy);
}

namespace {
    void dump_cg(const cluster_graph_t& cg, Log::logptr_t& log)
    {
        size_t mcount{0}, bcount{0};
        CS::value_t bval;
        for (const auto& vtx : GraphTools::mir(boost::vertices(cg))) {
            const auto& node = cg[vtx];
            if (node.code() == 'b') {
                const auto iblob = get<cluster_node_t::blob_t>(node.ptr);
                bval += CS::value_t(iblob->value(), iblob->uncertainty());
                ++bcount;
                continue;
            }
            if (node.code() == 'm') {
                const auto imeas = get<cluster_node_t::meas_t>(node.ptr);
                ++mcount;
            }
        }
        log->debug("cluster graph: vertices={} edges={} #blob={} bval={} #meas={}", boost::num_vertices(cg),
                   boost::num_edges(cg), bcount, bval, mcount);
    }
}  // namespace

bool GlobalGeomClustering::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }
    const auto in_graph = in->graph();
    dump_cg(in_graph, log);

    /// DEBUGONLY:
    std::unordered_map<cluster_vertex_t, int> clusters;
    using Filtered =
        typename boost::filtered_graph<cluster_graph_t, boost::keep_all, std::function<bool(cluster_vertex_t)> >;
    {
        Filtered bcg(in_graph, {}, [&](auto vtx) { return in_graph[vtx].code() == 'b'; });
        auto nclust = boost::connected_components(bcg, boost::make_assoc_property_map(clusters));
        log->debug("in #clusters {}", nclust);
    }

    /// rm old b-b edges
    using EFiltered =
        typename boost::filtered_graph<cluster_graph_t, std::function<bool(cluster_edge_t)>, boost::keep_all>;
    EFiltered fg_no_bb(in_graph,
                       [&](auto edge) {
                           auto source = boost::source(edge, in_graph);
                           auto target = boost::target(edge, in_graph);
                           if (in_graph[source].code() == 'b' and in_graph[target].code() == 'b') {
                               return false;
                           }
                           return true;
                       },
                       {});

    /// make new edges
    cluster_graph_t cg_new_bb;
    boost::copy_graph(fg_no_bb, cg_new_bb);
    /// DEBUGONLY:
    log->debug("rm bb:");
    dump_cg(cg_new_bb, log);
    grouped_geom_clustering(cg_new_bb, m_clustering_policy);

    /// DEBUGONLY:
    {
        Filtered bcg(cg_new_bb, {}, [&](auto vtx) { return cg_new_bb[vtx].code() == 'b'; });
        auto nclust = boost::connected_components(bcg, boost::make_assoc_property_map(clusters));
        log->debug("out #clusters {}", nclust);
    }

    // debug info
    log->debug("in_graph:");
    dump_cg(in_graph, log);
    log->debug("out:");
    dump_cg(cg_new_bb, log);

    out = std::make_shared<Aux::SimpleCluster>(cg_new_bb, in->ident());
    return true;
}
