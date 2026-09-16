#include "WireCellTbb/DataFlowGraph.h"

#include "WireCellUtil/Type.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"

#include <tbb/global_control.h>
#include <tbb/task_arena.h>

#include <iostream>

WIRECELL_FACTORY(TbbDataFlowGraph, WireCellTbb::DataFlowGraph, WireCell::IDataFlowGraph, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;
using namespace WireCellTbb;

DataFlowGraph::DataFlowGraph(int max_threads)
    : WireCell::Aux::Logger("DataFlowGraph", "tbb")
    , m_thread_limit(max_threads)
{
    // The graph/factory are built lazily in ensure() so they can be constructed
    // inside a sized task_arena when a thread limit is set.
}

DataFlowGraph::~DataFlowGraph() {}

void DataFlowGraph::ensure()
{
    if (m_graph) {
        return;
    }
    auto make = [this]() {
        m_graph = std::make_unique<tbb::flow::graph>();
        m_factory = std::make_unique<WrapperFactory>(*m_graph);
    };
    if (m_thread_limit > 0) {
        // A sized arena bounds Wire-Cell node parallelism to m_thread_limit
        // worker slots.  Unlike tbb::global_control(max_allowed_parallelism),
        // this is arena-LOCAL: it does not throttle a host framework's own TBB
        // work (e.g. Phlex) and creates no excess workers that busy-spin.
        m_arena = std::make_unique<tbb::task_arena>(m_thread_limit);
        m_arena->execute(make);
    }
    else {
        // No WCT limit: build in the current arena so the graph inherits whatever
        // the host framework (or the default arena) provides.
        make();
    }
}

Configuration DataFlowGraph::default_configuration() const
{
    Configuration cfg;
    cfg["max_threads"] = 0;
    cfg["summary"] = m_summary;
    // If set, write per-node execution intervals as JSON to this file.
    cfg["timeline"] = m_timeline;
    return cfg;
}

void DataFlowGraph::configure(const Configuration& cfg)
{
    if (! cfg["max_threads"].isNull()) {
        m_thread_limit = cfg["max_threads"].asInt();
    }
    m_summary = get(cfg, "summary", m_summary);
    m_timeline = get<std::string>(cfg, "timeline", m_timeline);
    // Enabling interval collection is global to all NodeInfo, so set it here,
    // before the graph is connected/built.
    if (!m_timeline.empty()) {
        WireCellTbb::NodeInfo::set_collect_intervals(true);
    }
}

bool DataFlowGraph::connect(INode::pointer tail, INode::pointer head, size_t sport, size_t rport)
{
    ensure();
    // Create/wire the TBB nodes inside the (possibly sized) arena so they belong
    // to it.  See ensure() and spng-fra.12.
    bool result = false;
    in_arena([&]() { result = do_connect(tail, head, sport, rport); });
    return result;
}

bool DataFlowGraph::do_connect(INode::pointer tail, INode::pointer head, size_t sport, size_t rport)
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

    Node mytail = (*m_factory)(tail);
    if (!mytail) {
        log->critical("no tail node wrapper for {}", tname);
        return false;
    }

    Node myhead = (*m_factory)(head);
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
    ensure();

    for (auto it : m_factory->seen()) {
        //log->debug("Initialize node of type: {}", demangle(it.first->signature()));
        it.second->initialize();
    }

    // Run the graph in its arena.  When a thread limit is configured the graph
    // was built inside a sized task_arena (ensure()), which bounds Wire-Cell
    // node parallelism to that worker count.  This replaces the former
    // tbb::global_control(max_allowed_parallelism) approach, which kept a full
    // hardware-sized worker pool that busy-spun (sched_yield) while a long node
    // ran -- burning every core and polluting the std::clock-based per-node
    // core-sec -- and which also composed by minimum with a host framework's
    // own global_control (clamping e.g. Phlex to WCT's limit).  See spng-fra.12.
    in_arena([this]() { m_graph->wait_for_all(); });

    if (m_summary) {
        double coretot_s = 0;
        double walltot_ms = 0;

        std::vector<WireCellTbb::Node> nodes(m_nodes.begin(), m_nodes.end());
        std::sort(nodes.begin(), nodes.end(),
                  [](const WireCellTbb::Node& a, const WireCellTbb::Node& b) {
                      return a->info().runtime() > b->info().runtime();
                  });
        for (const auto& node : nodes) {
            const auto& info = node->info();

            coretot_s += info.coretime();
            walltot_ms += std::chrono::duration_cast<std::chrono::milliseconds>(info.runtime()).count();

            std::stringstream ss;
            ss << info;
            if (m_summary < 2) {
                log->debug(ss.str());
            }
            else {
                log->info(ss.str());
            }
        }
        if (m_summary < 2) {
            log->debug("totals: wall={:.3f} s, core={:.3f} s, (counts include any EOS)",
                       walltot_ms / 1000.0, coretot_s);
        }
        else {
            log->info("totals: wall={:.3f} s, core={:.3f} s, (counts include any EOS)",
                      walltot_ms / 1000.0, coretot_s);
        }

    }

    if (!m_timeline.empty()) {
        // Emit the per-node execution timeline as JSON: each node's [start,end]
        // wall-clock (CLOCK_REALTIME) intervals in seconds since epoch, plus the
        // wall/core sums and call count.  Aligns with an external memory sampler.
        Configuration jroot;
        jroot["clock"] = "CLOCK_REALTIME";
        jroot["unit"] = "seconds";
        Configuration jnodes(Json::arrayValue);
        for (const auto& node : m_nodes) {
            const auto& info = node->info();
            Configuration jn;
            jn["instance"] = info.instance_name();
            jn["class"] = WireCell::type(*info.inode());
            jn["calls"] = (Json::UInt64) info.calls();
            jn["wall_sum"] = info.runtime().count();
            jn["core_sum"] = info.coretime();
            Configuration jints(Json::arrayValue);
            for (const auto& iv : info.intervals()) {
                Configuration jiv(Json::arrayValue);
                jiv.append(iv.first);
                jiv.append(iv.second);
                jints.append(jiv);
            }
            jn["intervals"] = jints;
            jnodes.append(jn);
        }
        jroot["nodes"] = jnodes;
        WireCell::Persist::dump(m_timeline, jroot);
        log->debug("wrote node timeline ({} nodes) to {}", m_nodes.size(), m_timeline);
    }

    return true;
}
