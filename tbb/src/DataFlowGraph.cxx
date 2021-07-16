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
  : m_graph()
  , m_factory(m_graph)
  , l(Log::logger("tbb"))
{
}

DataFlowGraph::~DataFlowGraph() {}

Configuration DataFlowGraph::default_configuration() const
{
    Configuration cfg;
    cfg["max_threads"] = 0;
    return cfg;
}

void DataFlowGraph::configure(const Configuration& cfg)
{
    if (! cfg["max_threads"].isNull()) {
        m_thread_limit = cfg["max_threads"].asInt();
    }
}

bool DataFlowGraph::connect(INode::pointer tail, INode::pointer head, size_t sport, size_t rport)
{
    using namespace WireCellTbb;

    Node mytail = m_factory(tail);
    if (!mytail) {
        l->critical("DFP: failed to get tail node wrapper for {}", demangle(tail->signature()));
        return false;
    }

    Node myhead = m_factory(head);
    if (!myhead) {
        l->critical("DFP: failed to get head node wrapper for {}", demangle(head->signature()));
        return false;
    }

    auto sports = mytail->sender_ports();
    if (sport < 0 || sports.size() <= sport) {
        l->critical("DFP: bad sender port number: {}", sport);
        return false;
    }

    auto rports = myhead->receiver_ports();
    if (rport < 0 || rports.size() <= rport) {
        l->critical("DFP: bad receiver port number: {}", rport);
        return false;
    }

    sender_type* s = sports[sport];
    if (!s) {
        l->critical("DFP: failed to get sender port {}", sport);
        return false;
    }

    receiver_type* r = rports[rport];
    if (!s) {
        l->critical("DFP: failed to get receiver port {}", rport);
        return false;
    }

    make_edge(*s, *r);
    return true;
}

bool DataFlowGraph::run()
{
    for (auto it : m_factory.seen()) {
        l->debug("Initialize node of type: {}", demangle(it.first->signature()));
        it.second->initialize();
    }

    std::unique_ptr<tbb::global_control> gc;
    if (m_thread_limit) {
        gc = std::make_unique<tbb::global_control>(
            tbb::global_control::max_allowed_parallelism,
            m_thread_limit);
    }
    m_graph.wait_for_all();

    return true;
}
