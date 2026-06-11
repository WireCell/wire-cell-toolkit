
#include "WireCellImg/DeadLiveMerging.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/SimpleCluster.h"

WIRECELL_FACTORY(DeadLiveMerging, WireCell::Img::DeadLiveMerging,
                 WireCell::INamed,
                 WireCell::IClusterFanin,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;

DeadLiveMerging::DeadLiveMerging()
    : Aux::Logger("DeadLiveMerging", "img")
{
}

void DeadLiveMerging::configure(const WireCell::Configuration& cfg)
{
    int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m <= 0) {
        raise<ValueError>("DeadLiveMerging multiplicity must be > 0");
    }
    m_multiplicity = m;

    m_tags.resize(m);

    // Tag entire input frame worth of traces in the output frame.
    auto jtags = cfg["tags"];
    for (int ind = 0; ind < m; ++ind) {
        m_tags[ind] = convert<std::string>(jtags[ind], "");
    }
}

WireCell::Configuration DeadLiveMerging::default_configuration() const
{
    Configuration cfg;
    return cfg;
}

std::vector<std::string> WireCell::Img::DeadLiveMerging::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

namespace {
}

bool DeadLiveMerging::operator()(const input_vector& in, output_pointer& out)
{
    out = nullptr;
    if (in.empty()) {
        return true;
    }

    // use reference_wrapper to avoid copying for this
    // std::vector<cluster_graph_t> in_cgraphs;
    std::vector<std::reference_wrapper<const cluster_graph_t>> in_cgraphs;

    for (const auto& inp : in) {
        if (!inp) {
            return true;
        }
        log->debug("Input: {}", Aux::dumps(inp->graph()));
        in_cgraphs.push_back(inp->graph());
    }

    auto merged_graph = Aux::merge_graphs(in_cgraphs);

    /// 1, tag dead/live blobs

    /// 2, make new edges between dead/live blobs

    // dummy implementation
    // out = in.front();
    // make SimpleCluster from merged graph
    /// FIXME: which ident to use?
    out = std::make_shared<Aux::SimpleCluster>(std::move(merged_graph), 0);
    log->debug("Output: {}", Aux::dumps(out->graph()));
    return true;
}