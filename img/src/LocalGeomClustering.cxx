
#include "WireCellImg/LocalGeomClustering.h"
#include "WireCellImg/CSGraph.h"
#include "WireCellImg/GeomClusteringUtil.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"

#include <iterator>
#include <chrono>
#include <fstream>

WIRECELL_FACTORY(LocalGeomClustering, WireCell::Img::LocalGeomClustering, WireCell::INamed, WireCell::IClusterFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;

Img::LocalGeomClustering::LocalGeomClustering()
  : Aux::Logger("LocalGeomClustering", "img")
{
}

Img::LocalGeomClustering::~LocalGeomClustering() {}

WireCell::Configuration Img::LocalGeomClustering::default_configuration() const
{
    WireCell::Configuration cfg;
    cfg["dryrun"] = m_dryrun;
    cfg["clustering_policy"] = m_clustering_policy;
    return cfg;
}

void Img::LocalGeomClustering::configure(const WireCell::Configuration& cfg)
{
    log->debug("{}", cfg);
    m_dryrun = get<bool>(cfg, "dryrun", m_dryrun);
    m_clustering_policy = get<std::string>(cfg, "clustering_policy", m_clustering_policy);
}

bool LocalGeomClustering::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }
    const auto in_graph = in->graph();
    log->debug("in_graph: {}", dumps(in_graph));

    /// Find b-b clusters using old edges
    std::unordered_map<cluster_vertex_t, int> clusters;
    std::unordered_set<int> clusterids;
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
    /// FIXME: guaranteed identical vdesc?
    boost::copy_graph(fg_no_bb, cg_new_bb);
    /// DEBUGONLY:
    log->debug("rm bb: {}", dumps(cg_new_bb));
    grouped_geom_clustering(cg_new_bb, m_clustering_policy, clusters);

    /// DEBUGONLY:
    {
        Filtered bcg(cg_new_bb, {}, [&](auto vtx) { return cg_new_bb[vtx].code() == 'b'; });
        auto nclust = boost::connected_components(bcg, boost::make_assoc_property_map(clusters));
        log->debug("out #clusters {}", nclust);
    }

    /// DEBUGONLY:
    log->debug("in: {}", dumps(in_graph));
    log->debug("out: {}", dumps(cg_new_bb));

    out = std::make_shared<Aux::SimpleCluster>(cg_new_bb, in->ident());
    if (m_dryrun) {
        out = std::make_shared<Aux::SimpleCluster>(in_graph, in->ident());
    }
    return true;
}
