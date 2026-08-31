
#include "WireCellImg/TestClusterShadow.h"
#include "WireCellImg/CSGraph.h"
#include "WireCellImg/Projection2D.h"

#include "WireCellAux/ClusterShadow.h"
#include "WireCellAux/SimpleCluster.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Rectangles.h"
#include "WireCellUtil/GraphTools.h"

#include <iterator>
#include <chrono>
#include <fstream>

WIRECELL_FACTORY(TestClusterShadow, WireCell::Img::TestClusterShadow,
                 WireCell::INamed,
                 WireCell::IClusterFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Img;
using namespace WireCell::Img::Projection2D;

Img::TestClusterShadow::TestClusterShadow()
    : Aux::Logger("TestClusterShadow", "img")
{
}

Img::TestClusterShadow::~TestClusterShadow() {}

WireCell::Configuration Img::TestClusterShadow::default_configuration() const
{
    WireCell::Configuration cfg;

    return cfg;
}

void Img::TestClusterShadow::configure(const WireCell::Configuration& cfg)
{
}

static void dump_cg(const cluster_graph_t& cg, Log::logptr_t& log)
{
    using value_t = ISlice::value_t;
    size_t mcount{0}, bcount{0};
    value_t bval;
    for (const auto& vtx : GraphTools::mir(boost::vertices(cg))) {
        const auto& node = cg[vtx];
        if (node.code() == 'b') {
            const auto iblob = get<cluster_node_t::blob_t>(node.ptr);
            bval += value_t(iblob->value(), iblob->uncertainty());
            ++bcount;
            continue;
        }
        if (node.code() == 'm') {
            const auto imeas = get<cluster_node_t::meas_t>(node.ptr);
            ++mcount;
        }
    }
    log->debug("cluster graph: vertices={} edges={} $blob={} bval={} #meas={}",
               boost::num_vertices(cg), boost::num_edges(cg),
               bcount, bval, mcount);
}

bool Img::TestClusterShadow::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }
    const auto& in_graph = in->graph();
    dump_cg(in_graph, log);
    auto& out_graph = in_graph;
    BlobShadow::graph_t bsgraph = BlobShadow::shadow(out_graph, 'w'); // or 'c'
    ClusterShadow::blob_cluster_map_t clusters;
    auto cs_graph = ClusterShadow::shadow(out_graph, bsgraph, clusters);
    log->debug("clusters.size(): {} ", clusters.size());

    // make a cluster -> blob map
    std::unordered_map<ClusterShadow::vdesc_t, std::set<cluster_vertex_t> > c2b;
    for(auto p : clusters) {
        c2b[p.second].insert(p.first);
    }

    // Count edges
    std::map<std::string, int> counters = {
        {"cs_edges", 0},
        {"cs_vertices", 0},
        {"nblobs", 0},
        {"nblobs_filter", 0},
        {"nblobs_clusters_tmp", 0},
        {"nblobs_clusters", 0}
    };

    {
        using cvertex_t = typename boost::graph_traits<cluster_graph_t>::vertex_descriptor;
        using Filtered = typename boost::filtered_graph<cluster_graph_t, boost::keep_all,
                                                        std::function<bool(cvertex_t)> >;
        Filtered bcg(out_graph, {}, [&](auto vtx) { return out_graph[vtx].code() == 'b'; });
        for (const auto& bcg_vertex : GraphTools::mir(boost::vertices(bcg))) {
            if (out_graph[bcg_vertex].code() != 'b') {
                log->debug("{} is not 'b'",bcg_vertex);
            }
            counters["nblobs_filter"] += 1;
        }
        ClusterShadow::blob_cluster_map_t clusters_tmp;
        boost::connected_components(bcg, boost::make_assoc_property_map(clusters_tmp));
        counters["nblobs_clusters_tmp"] += clusters_tmp.size();
    }
    counters["nblobs"] += boost::num_vertices(bsgraph);
    counters["nblobs_clusters"] += clusters.size();
    counters["cs_vertices"] += boost::num_vertices(cs_graph);
    std::ofstream fout("TestClusterShadow.log");
    for (const auto& cs_edge : GraphTools::mir(boost::edges(cs_graph))) {
        counters["cs_edges"] += 1;
        auto tail = boost::source(cs_edge, cs_graph);
        auto head = boost::target(cs_edge, cs_graph);
        auto b_tail = c2b[tail];
        auto b_head = c2b[head];
        log->debug("tail: {}, head: {}", tail, head);
        std::stringstream ss;
        ss << tail << " {";
        for (auto& b : c2b[tail]) {
            ss << b << " ";
        }
        ss << "} " << cs_graph[cs_edge].wpid.layer() << "-> ";
        ss << head << " {";
        for (auto& b : c2b[head]) {
            ss << b << " ";
        }
        ss << "} \n";
        fout << ss.str();
    }
    fout.close();
    for (auto c : counters) {
        log->debug("{} : {} ", c.first, c.second);
    }

    dump_cg(out_graph, log);
    out = std::make_shared<Aux::SimpleCluster>(out_graph, in->ident());
    return true;
}


