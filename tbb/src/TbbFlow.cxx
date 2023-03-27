#include "WireCellTbb/TbbFlow.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"

#include <string>

WIRECELL_FACTORY(TbbFlow, WireCellTbb::TbbFlow,
                 WireCell::IApplication, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCellTbb;

TbbFlow::TbbFlow()
    : Aux::Logger("TbbFlow", "tbb")
{
}

TbbFlow::~TbbFlow() {}

Configuration TbbFlow::default_configuration() const
{
    Configuration cfg;

    cfg["edges"] = Json::arrayValue;
    return cfg;
}

void TbbFlow::configure(const Configuration& cfg)
{
    std::string type, name, desc = get<std::string>(cfg, "dfp", "TbbDataFlowGraph");
    std::tie(type, name) = String::parse_pair(desc);
    m_dfp = Factory::lookup<IDataFlowGraph>(type, name);

    m_dfpgraph.configure(cfg["edges"]);

    log->debug("connecting graph with {} edges", cfg["edges"].size());

    for (auto thc : m_dfpgraph.connections()) {
        auto tail_tn = get<0>(thc);
        auto head_tn = get<1>(thc);
        auto conn = get<2>(thc);

        log->debug("connect: {}:{} ( {} -> {} ) {}:{}", tail_tn.type, tail_tn.name, conn.tail, conn.head,
                 head_tn.type, head_tn.name);

        INode::pointer tail_node = WireCell::Factory::lookup<INode>(tail_tn.type, tail_tn.name);
        INode::pointer head_node = WireCell::Factory::lookup<INode>(head_tn.type, head_tn.name);

        bool ok = m_dfp->connect(tail_node, head_node, conn.tail, conn.head);
        if (!ok) {
            THROW(ValueError() << errmsg{"edge connection error"});
        }
    }
}

void TbbFlow::execute()
{
    if (!m_dfp) {
        log->critical("TbbFlow: not configured");
        return;
    }

    log->debug("TbbFlow: run: ");
    m_dfp->run();
}
