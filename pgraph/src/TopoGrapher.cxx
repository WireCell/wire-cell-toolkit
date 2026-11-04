#include "WireCellPgraph/TopoGrapher.h"
#include "WireCellPgraph/Factory.h"
#include "WireCellIface/INode.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(TopoGrapher, WireCell::Pgraph::TopoGrapher, WireCell::IApplication, WireCell::IConfigurable)

using WireCell::get;
using namespace WireCell::Pgraph;

WireCell::Configuration TopoGrapher::default_configuration() const
{
    Configuration cfg;
    cfg["edges"] = Json::arrayValue;
    return cfg;
}

static std::pair<WireCell::INode::pointer, int> get_node(WireCell::Configuration jone)
{
    using namespace WireCell;
    std::string node = jone["node"].asString();

    // We should NOT be the one creating this component.
    auto nptr = WireCell::Factory::find_maybe_tn<INode>(node);
    if (!nptr) {
        raise<ValueError>("failed to get node \"%s\"", node);
    }

    int port = get(jone, "port", 0);
    return std::make_pair(nptr, port);
}

static std::string port_to_string(WireCell::Configuration node)
{
    std::stringstream ss;
    ss << node["node"].asString() << ":" << node["port"].asString();
    return ss.str();
}

static std::string edge_to_string(WireCell::Configuration tail,
                                  WireCell::Configuration head)
{
    return port_to_string(tail) + " -> " + port_to_string(head);
}

void TopoGrapher::configure(const WireCell::Configuration& cfg)
{
    m_verbosity = get(cfg, "verbosity", m_verbosity);
    m_graph.set_enable_em((m_verbosity == 2));

    Pgraph::Factory fac;
    log->debug("connecting: {} edges", cfg["edges"].size());
    for (auto jedge : cfg["edges"]) {
        auto tail = get_node(jedge["tail"]);
        auto head = get_node(jedge["head"]);

        SPDLOG_LOGGER_TRACE(log, "connecting: {}", jedge);

        bool ok = m_graph.connect(fac(tail.first), fac(head.first), tail.second, head.second);
        if (!ok) {
            log->critical("failed to connect edge: {}", jedge);
            raise<ValueError>("failed to connect edge %s", edge_to_string(jedge["tail"], jedge["head"]));
        }
    }
    if (!m_graph.connected()) {
        log->critical("graph not fully connected");
        raise<ValueError>("graph not fully connected");
    }
}

void TopoGrapher::execute()
{
    log->debug("executing graph (drain scheduler)");
    m_graph.execute_drain();
    log->debug("graph execution complete");
    if (m_verbosity) {
        m_graph.print_timers(m_verbosity == 2);
    }
}

TopoGrapher::TopoGrapher()
    : Aux::Logger("TopoGrapher", "pgraph")
{
}
TopoGrapher::~TopoGrapher() {}
