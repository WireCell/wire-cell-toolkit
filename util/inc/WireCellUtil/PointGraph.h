#ifndef WIRECELLUTIL_POINTGRAPH
#define WIRECELLUTIL_POINTGRAPH

#include "WireCellUtil/KDTree.h"

#include <boost/graph/adjacency_list.hpp>

namespace WireCell {

    /** A kd-tree queriable point cloud with graph-edge relationships
     * between points.
     *
     * This is largely a unifying container of these three things.
     */
    class PointGraph {
      public:

        ~PointGraph() = default;

        PointGraph();
        
        PointGraph(const PointGraph& rhs);
        PointGraph& operator=(const PointGraph& rhs);

        explicit PointGraph(const PointCloud::Dataset& nodes);

        PointGraph(const PointCloud::Dataset& nodes,
                   const PointCloud::Dataset& edges);

        KDTree::MultiQuery& mquery() { return m_mquery; }
        const KDTree::MultiQuery& mquery() const { return m_mquery; }

        PointCloud::Dataset& nodes() { return m_nodes; }
        const PointCloud::Dataset& nodes() const { return m_nodes; }

        PointCloud::Dataset& edges() { return m_edges; }
        const PointCloud::Dataset& edges() const { return m_edges; }

        /** Return a boost graph representation of the current state.
         *
         * The tails/heads may name the arrays which provide edge
         * endpoints.
         *
         * The IndexType should match the type of these arrays.
         * 
         * Neither vertex nor edge properties are set.
         */
        template <typename G, typename IndexType = int>
        G boost_graph(const std::string& tails = "tails",
                      const std::string& heads = "heads") const
        {
            using vdesc_t = typename boost::graph_traits<G>::vertex_descriptor;
            const size_t nvertices = m_nodes.size_major();
            G g(nvertices);
            if (nvertices == 0) {
                return g;
            }
            const size_t nedges = m_edges.size_major();
            if (nedges == 0) {
                return g;
            }
            const auto& t = m_edges.get(tails).elements<IndexType>();
            const auto& h = m_edges.get(heads).elements<IndexType>();
            for (size_t ind=0; ind<nedges; ++ind) {
                const vdesc_t tvtx = t[ind];
                const vdesc_t hvtx = h[ind];
                boost::add_edge(tvtx, hvtx, g);
            }
            return g;
        }


      private:

        PointCloud::Dataset m_nodes, m_edges;
        KDTree::MultiQuery m_mquery;

    };


}

#endif
