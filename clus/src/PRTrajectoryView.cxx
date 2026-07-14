#include "WireCellClus/PRTrajectoryView.h"
#include "WireCellClus/PRGraph.h"

namespace WireCell::Clus::PR {

    bool TrajectoryViewNodePredicate::operator()(const node_descriptor& desc) const
    {
        return view.has_node(desc);
    }


    bool TrajectoryViewEdgePredicate::operator()(const edge_descriptor& desc) const
    {
        return view.has_edge(desc);
    }
    


    TrajectoryView::TrajectoryView(Graph& graph)
        : m_graph(graph, TrajectoryViewNodePredicate(*this), TrajectoryViewEdgePredicate(*this))
        , m_nodes(0)
        , m_edges(0, EdgeDescriptorHash(graph), EdgeDescriptorEqual(graph))
    {}

    TrajectoryView::~TrajectoryView() {}

    const TrajectoryView::view_graph_type& TrajectoryView::view_graph() const
    {
        return m_graph;
    }

    bool TrajectoryView::has_node(node_descriptor desc) const
    {
        return m_nodes.count(desc) > 0;
    }

    bool TrajectoryView::has_edge(edge_descriptor desc) const
    {
        return m_edges.count(desc) > 0;
    }

    bool TrajectoryView::add_vertex(VertexPtr vtx)
    {
        if (! vtx->descriptor_valid()) {
            return false;
        }
        m_nodes.insert(vtx->get_descriptor());
        return true;
    }


    bool TrajectoryView::add_segment(SegmentPtr seg)
    {
        if (! seg->descriptor_valid()) {
            return false;
        }
        m_edges.insert(seg->get_descriptor());
        return true;
    }


    bool TrajectoryView::remove_vertex(VertexPtr vtx)
    {
        return m_nodes.erase(vtx->get_descriptor()) == 1;
    }
    
    bool TrajectoryView::remove_segment(SegmentPtr seg)
    {
        return m_edges.erase(seg->get_descriptor()) == 1;
    }
}
