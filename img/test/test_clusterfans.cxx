#include "WireCellImg/ClusterFanin.h"
#include "WireCellImg/ClusterFanout.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Testing.h"

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;
using WireCell::GraphTools::mir;

template<typename Fan>
void fan_cfg(Fan& fan, int multiplicity)
{
    auto cfg = fan.default_configuration();
    cfg["multiplicity"] = multiplicity;
    fan.configure(cfg);
}

template<typename Graph>
void dump(const Graph& graph, const std::string& name="")
{
    std::cerr << name << " vertices:";
    for (auto vtx : mir(boost::vertices(graph))) {
        std::cerr << " " << vtx;
    }
    std::cerr << "\n" << name << " edges:";
    for (auto edge : mir(boost::edges(graph))) {
        auto tvtx = boost::source(edge, graph);
        auto hvtx = boost::target(edge, graph);
        std::cerr << " (" << tvtx << "," << hvtx << ")";
    }
    std::cerr << "\n";
}

int main()
{
    const int multiplicity = 2;
    const int ident = 42;

    ClusterFanout cfo;
    ClusterFanin cfi;
    fan_cfg(cfo, multiplicity);
    fan_cfg(cfi, multiplicity);

    cluster_graph_t gorig;
    size_t num_vertices=0, num_edges=0;
    auto v1 = boost::add_vertex(gorig);
    ++num_vertices;
    auto v2 = boost::add_vertex(gorig);
    ++num_vertices;
    boost::add_edge(v1, v2, gorig);
    ++num_edges;

    dump(gorig, "gorig");

    auto corig = std::make_shared<SimpleCluster>(gorig, ident);

    // Fanout cluster
    ICluster::vector cmany;
    Assert(cfo(corig, cmany));
    Assert(cmany.size() == multiplicity);
    for (size_t ind=0; ind < (size_t)multiplicity; ++ind) {
        Assert(cmany[ind]->ident() == ident);
    }

    // Fanin cluster
    ICluster::pointer cout = nullptr;
    Assert(cfi(cmany, cout));
    Assert(cout);
    Assert(cout->ident() == ident);
    const auto& gout = cout->graph();
    dump(gout, "gout");
    Assert(boost::num_vertices(gout) == multiplicity*num_vertices);
    Assert(boost::num_edges(gout) == multiplicity*num_edges);

    // Fanout EOF
    Assert(cfo(nullptr, cmany));
    Assert(cmany.size() == multiplicity);
    for (size_t ind=0; ind < (size_t)multiplicity; ++ind) {
        Assert(cmany[ind] == nullptr);
    }

    // Fanin EOF
    Assert(cfi(cmany, cout));
    Assert(cout == nullptr);

    return 0;
}
