#include "WireCellPgraph/Graph.h"
#include "WireCellUtil/Type.h"
#include "WireCellUtil/String.h"

#include <unordered_map>
#include <unordered_set>
#include <ctime>
#include <boost/algorithm/string.hpp>

using WireCell::demangle;
using WireCell::String::format;
using namespace WireCell::Pgraph;

Graph::Graph()
  : l(Log::logger("pgraph"))
  , l_timer(Log::logger("timer"))
{
}

void Graph::add_node(Node* node)
{
    if (!node) {
        raise<ValueError>("Pgraph::Graph given null node");
    }
    m_nodes[node->instance()] = node;
}

void Graph::set_enable_em(bool flag) { m_enable_em = flag; }

bool Graph::connect(Node* tail, Node* head, size_t tpind, size_t hpind)
{
    auto& tports = tail->output_ports();
    auto& hports = head->input_ports();

    if (tpind >= tports.size()) {
        l->critical("tail port multiplicity mismatch with {}: {} >= {}", tail->ident(), tpind, tports.size());
        raise<ValueError>("tail port signature mismatch");
    }

    if (hpind >= hports.size()) {
        l->critical("head port multiplicity mismatch with {}: {} >= {}", head->ident(), hpind, hports.size());
        raise<ValueError>("tail port signature mismatch");
    }

    Port& tport = tports[tpind];
    Port& hport = hports[hpind];

    if (tport.signature() != hport.signature()) {
        l->critical("port signature mismatch: {}[{}] -> {}[{}]",
                    tail->ident(), tport.signature(),
                    head->ident(), hport.signature());
        raise<ValueError>("port signature mismatch");
        return false;
    }

    m_edges.push_back(std::make_pair(tail, head));
    Edge edge = std::make_shared<Queue>();

    tport.plug(edge);
    hport.plug(edge);

    add_node(tail);
    add_node(head);

    m_edges_forward[tail].push_back(head);
    m_edges_backward[head].push_back(tail);

    // was at in TRACE macro, maybe deserves a bit more prominence 
    l->debug("connect {}:({}:{}) -> {}({}:{})", tail->ident(), demangle(tport.signature()), tpind,
             head->ident(), demangle(hport.signature()), hpind);

    return true;
}

std::vector<Node*> Graph::sort_kahn()
{
    std::unordered_map<Node*, int> nincoming;
    for (auto th : m_edges) {
        nincoming[th.first] += 0;  // make sure all nodes represented
        nincoming[th.second] += 1;
    }

    std::vector<Node*> ret;
    std::map<size_t, Node*> seeds;

    for (auto it : nincoming) {
        if (it.second == 0) {
            Node* node = it.first;
            seeds[node->instance()] = node;
        }
    }

    while (!seeds.empty()) {
        auto [instance, node] = *seeds.begin();
        seeds.erase(instance);
        ret.push_back(node);

        for (auto h : m_edges_forward[node]) {
            nincoming[h] -= 1;
            if (nincoming[h] == 0) {
                seeds[h->instance()] = h;
            }
        }
    }
    return ret;
}

int Graph::execute_upstream(Node* node)
{
    int count = 0;
    for (auto parent : m_edges_backward[node]) {
        bool ok = call_node(parent);
        if (ok) {
            ++count;
            continue;
        }
        count += execute_upstream(parent);
    }
    bool ok = call_node(node);
    if (ok) {
        ++count;
    }
    return count;
}

// this bool indicates exception, and is probably ignored
bool Graph::execute()
{
    auto nodes = sort_kahn();
    l->debug("executing with {} nodes", nodes.size());

    while (true) {
        int count = 0;
        bool did_something = false;

        for (auto nit = nodes.rbegin(); nit != nodes.rend(); ++nit, ++count) {
            Node* node = *nit;

            auto core = std::clock();
            auto wall = std::chrono::high_resolution_clock::now();

            bool ok = call_node(node);

            m_nodes_timer[node].core += (std::clock() - core) / (double) CLOCKS_PER_SEC;

            auto delta = std::chrono::high_resolution_clock::now() - wall;
            const double wall_ms = std::chrono::duration_cast<std::chrono::milliseconds>(delta).count();
            m_nodes_timer[node].wall += wall_ms / 1000.0;

            if (ok) {
                if (m_enable_em) {
                  m_em(format("called %d: %s", count, node->ident()));
                }
                SPDLOG_LOGGER_TRACE(l, "ran node {}: {}", count, node->ident());
                did_something = true;
                break;  // start again from bottom of graph
            }
        }

        if (!did_something) {
            return true;  // it's okay to do nothing.
        }
    }
    return true;  // shouldn't reach
}

bool Graph::call_node(Node* node)
{
    if (!node) {
        l->error("graph call: got nullptr node");
        return false;
    }
    bool ok = (*node)();
    // this can be very noisy but useful to uncomment to understand
    // the graph execution order.
    if (ok) {
        SPDLOG_LOGGER_TRACE(l, "graph call [{}] called: {}", ok, node->ident());
    }
    return ok;
}

bool Graph::connected()
{
    bool okay = true;
    for (auto& [instance, node] : m_nodes) {
        auto bad = node->disconnected_ports();
        if (bad.empty()) continue;
        okay = false;
        l->warn("disconnected node #{}: {}", instance, node->ident());
        for (const auto& p : bad) {
            l->warn("\tempty: {}", p.ident());
        }
        for (const auto& p : node->input_ports()) {
            l->warn("\tiport: {}", p.ident());
        }
        for (const auto& p : node->output_ports()) {
            l->warn("\toport: {}", p.ident());
        }
    }
    return okay;
}

void Graph::print_timers(bool include_execmon) const
{
    std::multimap<float, std::pair<Node*, TimerClocks>> m;
    for (auto it : m_nodes_timer) { // to sort by wall time
        m.emplace(it.second.wall, std::make_pair(it.first, it.second));
    }

    double total_wall = 0;
    double total_core = 0;
    for (auto it = m.rbegin(); it != m.rend(); ++it) {
        auto [node, tc] = it->second;

        l_timer->info("Timer: {:.3f} wall-sec, {:.3f} core-sec: ({}) \"{}\"",
                      tc.wall, tc.core, node->cpptype_name(), node->instance_name()); 
        total_wall += tc.wall;
        total_core += tc.core;
    }
    l_timer->info("Timer: Total {:.3f} wall-sec, {:.3f} core-sec", total_wall, total_core);

    if (include_execmon) {
        l_timer->debug("ExecMon:\n{}", m_em.summary());
    }
}
