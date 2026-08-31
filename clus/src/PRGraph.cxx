#include "WireCellClus/PRGraph.h"

#include <cstdio>
#include <cstdlib>
#include <execinfo.h>

namespace WireCell::Clus::PR {

    // Determinism-debug (WCT_DET_DEBUG=2): log every segment creation with
    // its stable index and a mini backtrace so creation-order divergence
    // between two runs can be diffed and attributed to the calling function.
    static bool det_dbg_graph() {
        static const char* v = std::getenv("WCT_DET_DEBUG");
        static const bool on = (v && v[0] == '2');
        return on;
    }

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
        vtx->set_graph_index(index);
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
            seg->set_graph_index(index);
            if (det_dbg_graph()) {
                const auto& p1 = vtx1->wcpt().point;
                const auto& p2 = vtx2->wcpt().point;
                double wf[3] = {0,0,0}, wb[3] = {0,0,0};
                if (!seg->wcpts().empty()) {
                    const auto& f = seg->wcpts().front().point;
                    const auto& b = seg->wcpts().back().point;
                    wf[0]=f.x(); wf[1]=f.y(); wf[2]=f.z();
                    wb[0]=b.x(); wb[1]=b.y(); wb[2]=b.z();
                }
                fprintf(stderr, "WCT_DETA seg idx=%zu nw=%zu v1=(%.4f,%.4f,%.4f) v2=(%.4f,%.4f,%.4f) wf=(%.4f,%.4f,%.4f) wb=(%.4f,%.4f,%.4f)\n",
                        index, seg->wcpts().size(), p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(),
                        wf[0], wf[1], wf[2], wb[0], wb[1], wb[2]);
                void* bt[8];
                int nbt = backtrace(bt, 8);
                char** syms = backtrace_symbols(bt, nbt);
                for (int i = 1; i < nbt && i < 6; ++i)
                    fprintf(stderr, "WCT_DETA   bt[%d] %s\n", i, syms ? syms[i] : "?");
                free(syms);
            }
            return true;
        }

        // Edge already existed, assure its object is this one.
        // Inherit the existing edge's graph index so m_graph_index is not left at SIZE_MAX.
        g[desc].segment = seg;
        seg->set_graph_index(g[desc].index);

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


    std::pair<VertexPtr, VertexPtr> find_vertices(Graph& graph, SegmentPtr seg)
    {
        if (! seg->descriptor_valid()) { return std::pair<VertexPtr, VertexPtr>{}; }

        auto ed = seg->get_descriptor();

        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);

        auto vtx1 = graph[vd1].vertex;
        auto vtx2 = graph[vd2].vertex;

        // Use wcpts if available, fall back to fits, then return unordered pair.
        WireCell::Point ept;
        if (!seg->wcpts().empty()) {
            ept = seg->wcpts().front().point;
        } else if (!seg->fits().empty()) {
            ept = seg->fits().front().point;
        } else {
            return std::make_pair(vtx1, vtx2);
        }

        double d1 = ray_length(Ray{vtx1->wcpt().point, ept});
        double d2 = ray_length(Ray{vtx2->wcpt().point, ept});

        if (d1 < d2) {
            return std::make_pair(vtx1, vtx2);
        }
        return std::make_pair(vtx2, vtx1);        
    }

    VertexPtr find_other_vertex(Graph& graph, SegmentPtr seg, VertexPtr vertex)
    {
        if (! seg->descriptor_valid()) { return nullptr; }
        if (! vertex->descriptor_valid()) { return nullptr; }

        auto ed = seg->get_descriptor();
        auto vd = vertex->get_descriptor();

        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);

        auto [ed2,ingraph] = boost::edge(vd1, vd2, graph);
        if (!ingraph)  { return nullptr; }

        // Check if the given vertex is connected to this segment
        if (vd == vd1) {
            return graph[vd2].vertex;
        }
        else if (vd == vd2) {
            return graph[vd1].vertex;
        }

        // Given vertex is not connected to this segment
        return nullptr;
    }

    SegmentPtr find_segment(Graph& graph, VertexPtr vtx1, VertexPtr vtx2)
    {
        if (! vtx1->descriptor_valid()) { return nullptr; }
        if (! vtx2->descriptor_valid()) { return nullptr; }

        auto vd1 = vtx1->get_descriptor();
        auto vd2 = vtx2->get_descriptor();

        // Check if edge exists between the two vertices
        auto [ed, exists] = boost::edge(vd1, vd2, graph);
        if (!exists) { return nullptr; }

        // Return the segment associated with this edge
        return graph[ed].segment;
    }

}

    
