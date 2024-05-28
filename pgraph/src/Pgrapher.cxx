#include "WireCellPgraph/Pgrapher.h"
#include "WireCellPgraph/Factory.h"
#include "WireCellIface/INode.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(Pgrapher, WireCell::Pgraph::Pgrapher, WireCell::IApplication, WireCell::IConfigurable)

using WireCell::get;
using namespace WireCell::Pgraph;

WireCell::Configuration Pgrapher::default_configuration() const
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

static
std::string port_to_string(WireCell::Configuration node)
{
    std::stringstream ss;
    ss << node["node"].asString() << ":" << node["port"].asString();
    return ss.str();
}

static
std::string edge_to_string(WireCell::Configuration tail,
                           WireCell::Configuration head)
{
    return port_to_string(tail) + " -> " + port_to_string(head);
}

void Pgrapher::configure(const WireCell::Configuration& cfg)

{
    m_verbosity = get(cfg, "verbosity", m_verbosity);

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

void Pgrapher::execute()
{
    log->debug("executing graph");
    m_graph.execute();
    log->debug("graph execution complete");
    if (m_verbosity) {
        m_graph.print_timers(m_verbosity == 2);
    }
}

Pgrapher::Pgrapher()
    : Aux::Logger("Pgrapher", "pgraph")
{
}
Pgrapher::~Pgrapher() {}
