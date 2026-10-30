#include "WireCellPgraph/Graph.h"
#include "WireCellUtil/Type.h"
#include "WireCellUtil/String.h"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <ctime>
#include <boost/algorithm/string.hpp>
#include "WireCellPgraph/TorchnvTools.h"

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

// Always fire the sink-most runnable node so no more than one datum sits
// on any edge, discovering the runnable node incrementally instead of
// rescanning the whole topological order after every firing.
//
// Each node's topological index is its fixed "rank" (sink-most = highest).
// A node's runnability can change ONLY when a neighbor fires: a child
// consuming empties a parent's output edge (clearing the parent's "don't
// run while output pending" guard), and a parent producing fills a child's
// input.  So after each successful firing we re-offer just that node plus
// its parents and children to a rank-ordered ready set, and pop the
// highest rank.  Rank is a total order, so this fires the same node
// sequence as a reverse scan of the topological order restarted after each
// firing (the original algorithm), but in O((V+E) log V) instead of O(V^2).
//
// The returned bool indicates exception, and is probably ignored.
bool Graph::execute()
{
    auto nodes = sort_kahn();
    const size_t nnodes = nodes.size();
    l->debug("executing with {} nodes", nnodes);

    std::unordered_map<Node*, size_t> rank;
    rank.reserve(nnodes);
    for (size_t ind = 0; ind < nnodes; ++ind) {
        rank[nodes[ind]] = ind;
    }

    std::set<size_t> ready;             // ranks of queued candidate nodes
    std::vector<char> queued(nnodes, 0);  // dedup: at most one entry per node

    auto enqueue = [&](Node* n) {
        auto it = rank.find(n);
        if (it == rank.end()) return;   // not in the sorted DAG (e.g. in a cycle)
        const size_t r = it->second;
        if (queued[r]) return;
        queued[r] = 1;
        ready.insert(r);
    };

    // Seed with sources (no input ports).  Every other runnability change
    // is triggered by a neighbor firing, so seeding sources is sufficient.
    for (size_t ind = 0; ind < nnodes; ++ind) {
        if (nodes[ind]->input_ports().empty()) {
            enqueue(nodes[ind]);
        }
    }

    while (!ready.empty()) {
        auto rit = std::prev(ready.end());  // highest rank == sink-most
        const size_t r = *rit;
        ready.erase(rit);
        queued[r] = 0;
        Node* node = nodes[r];

        NVTX_SCOPED_RANGE((node->instance_name()).c_str());
        auto core = std::clock();
        auto wall = std::chrono::high_resolution_clock::now();

        bool ok = call_node(node);

        m_nodes_timer[node].core += (std::clock() - core) / (double) CLOCKS_PER_SEC;
        auto delta = std::chrono::high_resolution_clock::now() - wall;
        const double wall_ms = std::chrono::duration_cast<std::chrono::milliseconds>(delta).count();
        m_nodes_timer[node].wall += wall_ms / 1000.0;

        if (!ok) {
            // Was not actually runnable; it will be re-offered when a
            // neighbor's firing next changes its input/output state.
            continue;
        }

        if (m_enable_em) {
            m_em(String::format("called: %s", node->ident()));
        }
        SPDLOG_LOGGER_TRACE(l, "ran node: {}", node->ident());

        // Only this node and its immediate neighbors can have become runnable.
        enqueue(node);
        for (auto h : m_edges_forward[node]) enqueue(h);
        for (auto t : m_edges_backward[node]) enqueue(t);
    }

    return true;
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

        l_timer->info("Timer: {:.3f} wall-sec, {:.3f} core-sec:  ({}) \"{}\"",
                      tc.wall, tc.core, node->cpptype_name(), node->instance_name());
        total_wall += tc.wall;
        total_core += tc.core;
    }
    l_timer->info("Timer: Total {:.3f} wall-sec, {:.3f} core-sec", total_wall, total_core);

    if (include_execmon) {
        l_timer->debug("ExecMon:\n{}", m_em.summary());
    }
}
