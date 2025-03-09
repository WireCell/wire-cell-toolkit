#include "WireCellImg/ClusterScopeFilter.h"
#include "WireCellImg/CSGraph.h"
#include "WireCellImg/GeomClusteringUtil.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/TimeKeeper.h"

#include <iterator>
#include <chrono>
#include <fstream>

WIRECELL_FACTORY(ClusterScopeFilter, WireCell::Img::ClusterScopeFilter, WireCell::INamed, WireCell::IClusterFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;

Img::ClusterScopeFilter::ClusterScopeFilter()
  : Aux::Logger("ClusterScopeFilter", "img")
{
}

Img::ClusterScopeFilter::~ClusterScopeFilter() {}

WireCell::Configuration Img::ClusterScopeFilter::default_configuration() const
{
    WireCell::Configuration cfg;
    return cfg;
}

void Img::ClusterScopeFilter::configure(const WireCell::Configuration& cfg)
{
    m_face_index = get<int>(cfg, "face_index", m_face_index);

    log->debug("{}", cfg);
}

bool ClusterScopeFilter::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }

    TimeKeeper tk(fmt::format("ClusterScopeFilter"));

    const auto in_graph = in->graph();
    log->debug("in_graph: {}", dumps(in_graph));

    log->debug(tk(fmt::format("start delete some blobs")));
    using VFiltered =
        typename boost::filtered_graph<cluster_graph_t, boost::keep_all, std::function<bool(cluster_vertex_t)> >;
    VFiltered filtered_graph(in_graph, {}, [&](auto vtx) {
        if (in_graph[vtx].code() != 'b') return false;
        const auto iblob = std::get<cluster_node_t::blob_t>(in_graph[vtx].ptr);
        if (iblob->face()->which() == m_face_index) return true;
        return false;
    });

    WireCell::cluster_graph_t out_graph;
    boost::copy_graph(filtered_graph, out_graph);
    log->debug("out_graph: {}", dumps(out_graph));

    /// output
    log->debug("in_graph: {}", dumps(in_graph));
    log->debug("out_graph: {}", dumps(out_graph));

    out = std::make_shared<Aux::SimpleCluster>(out_graph, in->ident());
    return true;
}
