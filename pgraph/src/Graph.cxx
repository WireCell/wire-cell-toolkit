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

void Graph::add_node(Node* node) { m_nodes.insert(node); }

bool Graph::connect(Node* tail, Node* head, size_t tpind, size_t hpind)
{
    Port& tport = tail->output_ports()[tpind];
    Port& hport = head->input_ports()[hpind];
    if (tport.signature() != hport.signature()) {
        l->critical("port signature mismatch: \"{}\" != \"{}\"", tport.signature(), hport.signature());
        THROW(ValueError() << errmsg{"port signature mismatch"});
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
    std::unordered_set<Node*> seeds;

    for (auto it : nincoming) {
        if (it.second == 0) {
            seeds.insert(it.first);
        }
    }

    while (!seeds.empty()) {
        Node* t = *seeds.begin();
        seeds.erase(t);
        ret.push_back(t);

        for (auto h : m_edges_forward[t]) {
            nincoming[h] -= 1;
            if (nincoming[h] == 0) {
                seeds.insert(h);
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

    for (Node* node : nodes) {
        m_nodes_timer[node] = 0.0;
    }

    while (true) {
        int count = 0;
        bool did_something = false;

        for (auto nit = nodes.rbegin(); nit != nodes.rend(); ++nit, ++count) {
            Node* node = *nit;

            auto start = std::clock();

            bool ok = call_node(node);

            m_nodes_timer[node] += (std::clock() - start) / (double) CLOCKS_PER_SEC;

            if (ok) {
                m_em(format("called %d: %s", count, node->ident()));
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
    for (auto n : m_nodes) {
        auto bad = n->disconnected_ports();
        if (bad.empty()) continue;
        okay = false;
        l->warn("disconnected node: {}", n->ident());
        for (const auto& p : bad) {
            l->warn("\tport: {}", p.ident());
        }
    }
    return okay;
}

void Graph::print_timers(bool include_execmon) const
{
    std::multimap<float, Node*> m;
    double total_time = 0;
    for (auto it : m_nodes_timer) {
        m.emplace(it.second, it.first);
    }
    std::vector<Node*> ordered;
    for (auto it = m.rbegin(); it != m.rend(); ++it) {
        ordered.push_back(it->second);
        std::string iden = it->second->ident();
        std::vector<std::string> tags;
        boost::split(tags, iden, [](char c) { return c == ' '; });
        l_timer->info("Timer: {} : {} sec", tags[2].substr(5), it->first);
        total_time += it->first;
    }
    l_timer->info("Timer: Total node execution : {} sec", total_time);

    if (include_execmon) {
        l_timer->debug("ExecMon:\n{}", m_em.summary());
    }
}
