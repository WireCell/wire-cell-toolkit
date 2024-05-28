#include "WireCellTbb/WrapperFactory.h"
#include "WireCellTbb/SourceCat.h"
#include "WireCellTbb/SinkCat.h"
#include "WireCellTbb/JoinCat.h"
#include "WireCellTbb/FaninCat.h"
#include "WireCellTbb/FanoutCat.h"
#include "WireCellTbb/QueuedoutCat.h"
#include "WireCellTbb/FunctionCat.h"
#include "WireCellTbb/HydraCat.h"

using namespace WireCell;
using namespace WireCellTbb;

WrapperFactory::WrapperFactory(tbb::flow::graph& graph)
  : m_graph(graph)
{
    // gotta add one for each Wire Cell category
    bind_maker<SourceNodeWrapper>(INode::sourceNode);
    bind_maker<SinkNodeWrapper>(INode::sinkNode);
    bind_maker<QueuedoutWrapper>(INode::queuedoutNode);
    bind_maker<JoinWrapper>(INode::joinNode);
    bind_maker<FaninWrapper>(INode::faninNode);
    bind_maker<FanoutWrapper>(INode::fanoutNode);
    bind_maker<FunctionWrapper>(INode::functionNode);
    bind_maker<HydraWrapper>(INode::hydraNode);
    // fixme: add more ...
}

Node WrapperFactory::operator()(INode::pointer wcnode)
{
    auto nit = m_nodes.find(wcnode);
    if (nit != m_nodes.end()) {
        return nit->second;
    }

    auto mit = m_factory.find(wcnode->category());
    if (mit == m_factory.end()) {
        return nullptr;
    }
    auto maker = mit->second;

    Node node = (*maker)(m_graph, wcnode);
    m_nodes[wcnode] = node;
    return node;
}
