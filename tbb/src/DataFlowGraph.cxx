#include "WireCellTbb/DataFlowGraph.h"

#include "WireCellUtil/Type.h"
#include "WireCellUtil/NamedFactory.h"

#include <tbb/global_control.h>

#include <iostream>

WIRECELL_FACTORY(TbbDataFlowGraph, WireCellTbb::DataFlowGraph, WireCell::IDataFlowGraph, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;
using namespace WireCellTbb;

DataFlowGraph::DataFlowGraph(int max_threads)
    : WireCell::Aux::Logger("DataFlowGraph", "tbb")
    , m_graph()
    , m_factory(m_graph)
{
}

DataFlowGraph::~DataFlowGraph() {}

Configuration DataFlowGraph::default_configuration() const
{
    Configuration cfg;
    cfg["max_threads"] = 0;
    cfg["summary"] = m_summary;
    return cfg;
}

void DataFlowGraph::configure(const Configuration& cfg)
{
    if (! cfg["max_threads"].isNull()) {
        m_thread_limit = cfg["max_threads"].asInt();
    }
    m_summary = get(cfg, "summary", m_summary);
}

bool DataFlowGraph::connect(INode::pointer tail, INode::pointer head, size_t sport, size_t rport)
{
    using namespace WireCellTbb;

    const std::string tname = demangle(tail->signature());
    const std::string hname = demangle(head->signature());

    {                           // check WCT interface level info
        const auto& ttypes = tail->output_types();
        if (sport < 0 || ttypes.size() <= sport) {
            log->critical("bad tail port index: {} out of {} for {}", sport, ttypes.size(), tname);
            return false;
        }
        const auto& htypes = head->input_types();
        if (rport < 0 || htypes.size() <= rport) {
            log->critical("bad head port index: {} out of {} for {}", rport, htypes.size(), hname);
            return false;
        }

        if (ttypes[sport] != htypes[rport]) {
            log->critical("edge type mismatch: tail:{}[{}]={} head={}[{}]={}",
                          tname, sport, demangle(ttypes[sport]),
                          hname, rport, demangle(htypes[rport]));
            return false;
        }
    }

    Node mytail = m_factory(tail);
    if (!mytail) {
        log->critical("no tail node wrapper for {}", tname);
        return false;
    }

    Node myhead = m_factory(head);
    if (!myhead) {
        log->critical("no head node wrapper for {}", hname);
        return false;
    }

    auto sports = mytail->sender_ports();
    if (sport < 0 || sports.size() <= sport) {
        log->critical("bad sender port index: {} out of {} for {}", sport, sports.size(), tname);
        return false;
    }

    auto rports = myhead->receiver_ports();
    if (rport < 0 || rports.size() <= rport) {
        log->critical("bad receiver port index: {} out of {} for {}", rport, rports.size(), hname);
        return false;
    }

    sender_type* s = sports[sport];
    if (!s) {
        log->critical("no sender port {} for {}", sport, tname);
        return false;
    }

    receiver_type* r = rports[rport];
    if (!s) {
        log->critical("no receiver port {} for {}", rport, hname);
        return false;
    }

    make_edge(*s, *r);
    m_nodes.insert(mytail);
    m_nodes.insert(myhead);
    return true;
}

bool DataFlowGraph::run()
{
    for (auto it : m_factory.seen()) {
        //log->debug("Initialize node of type: {}", demangle(it.first->signature()));
        it.second->initialize();
    }

    std::unique_ptr<tbb::global_control> gc;
    if (m_thread_limit) {
        gc = std::make_unique<tbb::global_control>(
            tbb::global_control::max_allowed_parallelism,
            m_thread_limit);
    }
    m_graph.wait_for_all();

    if (m_summary) {
        std::vector<WireCellTbb::Node> nodes(m_nodes.begin(), m_nodes.end());
        std::sort(nodes.begin(), nodes.end(),
                  [](const WireCellTbb::Node& a, const WireCellTbb::Node& b) {
                      return a->info().runtime() > b->info().runtime();
                  });
        for (const auto& node : nodes) {
            std::stringstream ss;
            ss << node->info();
            if (m_summary < 2) {
                log->debug(ss.str());
            }
            else {
                log->info(ss.str());
            }
        }
    }

    return true;
}
