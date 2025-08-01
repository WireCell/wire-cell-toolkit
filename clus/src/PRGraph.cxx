#include "WireCellClus/PRGraph.h"

namespace WireCell::Clus::PR {

    bool add_vertex(Graph& g, VertexPtr vtx)
    {
        if (vtx->descriptor_valid()) {
            return false;
        }
        auto& gb = g[boost::graph_bundle];
        const size_t index = gb.num_node_indices;
        auto desc = boost::add_vertex(NodeBundle{vtx, index}, g);
        ++ gb.num_node_indices;
        vtx->set_descriptor(desc);
        return true;
    }

    bool remove_vertex(Graph& graph, VertexPtr vtx)
    {
        if (! vtx->descriptor_valid()) { return false; }
        auto desc = vtx->get_descriptor();
        boost::remove_vertex(desc, graph);
        vtx->invalidate_descriptor();
        return true;
    }

    bool add_segment(Graph& g, SegmentPtr seg, VertexPtr vtx1, VertexPtr vtx2)
    {
        bool changed = false;
        changed = add_vertex(g, vtx1) || changed;
        changed = add_vertex(g, vtx2) || changed;

        if (seg->descriptor_valid()) {
            return changed;
        }

        auto& gb = g[boost::graph_bundle];
        const size_t index = gb.num_edge_indices;
        auto [desc,added] = boost::add_edge(vtx1->get_descriptor(),
                                            vtx2->get_descriptor(),
                                            EdgeBundle{seg, index}, g);

        seg->set_descriptor(desc);

        // Edge was added
        if (added) {
            ++ gb.num_edge_indices;
            seg->set_descriptor(desc);
            return true;
        }

        // Edge already existed, assure its object is this one.
        g[desc].segment = seg;

        return changed;
    }

    bool remove_segment(Graph& graph, SegmentPtr seg)
    {
        if (! seg->descriptor_valid()) { return false; }
        auto desc = seg->get_descriptor();
        boost::remove_edge(desc, graph);
        seg->invalidate_descriptor();
        return true;
    }


    std::pair<VertexPtr, VertexPtr> find_endpoints(Graph& graph, SegmentPtr seg)
    {
        if (! seg->descriptor_valid()) { return std::pair<VertexPtr, VertexPtr>{}; }

        auto ed = seg->get_descriptor();

        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);

        auto [ed2,ingraph] = boost::edge(vd1, vd2, graph);
        if (!ingraph)  { return std::pair<VertexPtr, VertexPtr>{}; }

        return std::make_pair(graph[vd1].vertex, graph[vd2].vertex);
    }

}

    
