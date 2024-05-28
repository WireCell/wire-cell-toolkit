// A minimally complete example of a tbb dfp

#include "tbb_mock.h"

#include "WireCellTbb/NodeWrapper.h"
#include "WireCellTbb/SinkCat.h"
#include "WireCellTbb/SourceCat.h"
#include "WireCellTbb/QueuedoutCat.h"

#include "WireCellUtil/Testing.h"

#include <tbb/flow_graph.h>

#include <string>
#include <deque>
#include <iostream>

using namespace std;

// fixme: this fakes the factory until we clean up nodes to allow
// empty c'tor and use configuration.
WireCell::INode::pointer get_node(const std::string& node_desc)
{
    using namespace WireCell;

    if (node_desc == "source") {  // note actual desc should be class or class:inst
        return INode::pointer(new WireCellTbb::MockDepoSource);
    }
    if (node_desc == "drift") {  // note actual desc should be class or class:inst
        return INode::pointer(new WireCellTbb::MockDrifter);
    }
    if (node_desc == "sink") {  // note actual desc should be class or class:inst
        return INode::pointer(new WireCellTbb::MockDepoSink);
    }
    return nullptr;
}

WireCellTbb::Node make_node(tbb::flow::graph& graph, const std::string& node_desc, const std::string& name)
{
    using namespace WireCell;
    using namespace WireCellTbb;

    INode::pointer wcnode = get_node(node_desc);
    if (!wcnode) {
        cerr << "Failed to get node for " << node_desc << endl;
        return nullptr;
    }
    auto inamed = std::dynamic_pointer_cast<INamed>(wcnode);
    if (inamed) {
        inamed->set_name(name);
    }

    cerr << "Getting node from category: " << wcnode->category() << endl;
    switch (wcnode->category()) {
    case INode::sourceNode:
        return Node(new SourceNodeWrapper(graph, wcnode));
    case INode::sinkNode:
        return Node(new SinkNodeWrapper(graph, wcnode));
    // case INode::functionNode:
    // 	return Node(new FunctionWrapper(graph, wcnode));
    case INode::queuedoutNode:
        return Node(new QueuedoutWrapper(graph, wcnode));
    default:
        return nullptr;
    }
    return nullptr;
}

bool connect(WireCellTbb::Node sender, WireCellTbb::Node receiver, size_t sport = 0, size_t rport = 0);
bool connect(WireCellTbb::Node sender, WireCellTbb::Node receiver, size_t sport, size_t rport)
{
    using namespace WireCellTbb;

    Assert(sender);
    Assert(receiver);
    auto sports = sender->sender_ports();
    auto rports = receiver->receiver_ports();

    Assert(sports.size() > sport);
    Assert(rports.size() > rport);

    sender_type* s = sports[sport];
    receiver_type* r = rports[rport];
    Assert(s);
    Assert(r);

    cerr << "Connecting " << s << " and " << r << endl;
    make_edge(*s, *r);
    return true;
}

int main()
{
    using namespace WireCellTbb;

    tbb::flow::graph graph;
    Node source = make_node(graph, "source", "source1");
    Assert(source);
    Node drift = make_node(graph, "drift", "drift1");
    Assert(drift);
    Node sink = make_node(graph, "sink", "sink1");
    Assert(sink);

    Assert(connect(source, drift));
    Assert(connect(drift, sink));

    // fixme: in general all nodes should be initialize()'d but so far only source nodes need it.
    source->initialize();

    graph.wait_for_all();

    auto source_info = source->info();
    std::cerr << source_info << std::endl;
    auto drift_info = drift->info();
    std::cerr << drift_info << std::endl;
    auto sink_info = sink->info();
    std::cerr << sink_info << std::endl;
    return 0;
}
