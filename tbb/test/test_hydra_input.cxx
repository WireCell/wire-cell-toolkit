// Exercise TBB's queuing and ability to hold inputs in the input queue.

#include "WireCellTbb/HydraCat.h"

#include <cassert>
#include <vector>
#include <iostream>
#include <sstream> // avoid threads chopping up i/o lines.
#include <chrono>
#include <thread>

using namespace WireCell;

// message types
using WireCellTbb::msg_t;
using WireCellTbb::tagged_msg_t;

// port types
using WireCellTbb::sender_type;
using WireCellTbb::receiver_type;

// node types
using WireCellTbb::src_node;
using sink_node = tbb::flow::function_node<tagged_msg_t>;

using WireCellTbb::build_hydra_input;
using WireCellTbb::receiver_port_vector;

// The "payload" for the messages will be a simple int giving the ID number of
// the source from which the message originated.

class SourceBody {
    const int id;
    const size_t sleepms;
    const size_t max;
    mutable size_t count{0};

  public:
    ~SourceBody() {}

    SourceBody(int id, size_t sleepms=100, size_t max=100)
        : id(id), sleepms(sleepms), max(max) {
        std::cerr << "source " << id << ", sleep=" << sleepms << "ms, max=" << max << "\n";
    }
    msg_t operator()(tbb::flow_control& fc) {
        std::this_thread::sleep_for(std::chrono::milliseconds(sleepms));
        if (count < max) {
            std::stringstream ss;
            ss << "-> " << id << ": " << count
               << " " << (void*) this << "\n";
            std::cerr << ss.str();
            boost::any aid = id;
            return std::make_pair(count++, aid);
        }
        fc.stop();
        return {};
    }
};

struct SinkBody {
    ~SinkBody() {}

    SinkBody() {
        std::cerr << "sink\n";
    }

    tbb::flow::continue_msg operator()(const tagged_msg_t& in)
    {
        auto& [p,m] = in;
        auto& [n,aid] = m;
        auto id = boost::any_cast<int>(aid);
        std::stringstream ss;
        ss << "<- " << id << ": " << n << " " << p << "\n";
        std::cerr << ss.str();
        return {};
    }
};

void test_input()
{
    constexpr size_t multiplicity = 2;

    std::cerr << "test_hydra: constructing with multiplicity " << multiplicity << "\n";

    tbb::flow::graph graph;

    // Build up the "front end" input side of the hydra by hand.
    std::vector<tbb::flow::graph_node*> hydra_internal_nodes;
    receiver_port_vector rx;
    tbb::flow::sender<tagged_msg_t>* fe = build_hydra_input<multiplicity>(graph, rx, hydra_internal_nodes);
    assert(rx.size() == multiplicity);

    // The sink here mimics the input side of the hydra backend output half.
    sink_node sink(graph, 1, SinkBody());
    tbb::flow::make_edge(*fe, sink);

    std::vector<size_t> sleeps = {100,200};
    std::vector<src_node> src_nodes;
    for (size_t ind=0; ind<multiplicity; ++ind) {
        src_nodes.emplace_back(graph, SourceBody(ind, sleeps[ind]));
        tbb::flow::make_edge(src_nodes[ind], *rx[ind]);
    }

    for (size_t ind=0; ind<multiplicity; ++ind) {
        src_nodes[ind].activate();
    }

    std::cerr << "test_hydra: running graph:\n";
    graph.wait_for_all();
    std::cerr << "test_hydra: done\n";
}


int main()
{
    test_input();
    return 0;
}
