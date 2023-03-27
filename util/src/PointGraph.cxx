#include "WireCellUtil/PointGraph.h"

using namespace WireCell;

PointGraph::PointGraph()
    : m_mquery(m_nodes)
{
}
        
PointGraph::PointGraph(const PointCloud::Dataset& nodes)
    : m_nodes(nodes)
    , m_mquery(m_nodes)
{
}

PointGraph::PointGraph(const PointCloud::Dataset& nodes,
                       const PointCloud::Dataset& edges)
    : m_nodes(nodes)
    , m_edges(edges)
    , m_mquery(m_nodes)
{
}

PointGraph::PointGraph(const PointGraph& rhs)
    : m_nodes(rhs.m_nodes)
    , m_edges(rhs.m_edges)
    , m_mquery(m_nodes)
{
}
PointGraph& PointGraph::operator=(const PointGraph& rhs)
{
    m_nodes = rhs.m_nodes;
    m_edges = rhs.m_edges;
    m_mquery = KDTree::MultiQuery(m_nodes);
    return *this;
}
