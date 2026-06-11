#include "WireCellIface/ICluster.h"
#include "WireCellUtil/Graph.h"
namespace WireCell::Aux {

    class SimpleCluster : public ICluster {
       public:
        SimpleCluster(const cluster_graph_t& g, int ident = 0)
          : m_ident(ident)
        {
            boost::copy_graph(g, m_graph);
        }
        // Move overload: a caller whose graph is dead after construction
        // avoids the full copy_graph + destruction of the source (one
        // whole-graph copy per pipeline stage).  Iteration-order safe:
        // the moved graph is the very object the copy would reproduce.
        // m_graph(std::move(g)) would silently call the copy ctor
        // (adjacency_list has no move members), hence move_graph().
        SimpleCluster(cluster_graph_t&& g, int ident = 0)
          : m_ident(ident)
        {
            WireCell::move_graph(m_graph, g);
        }
        virtual int ident() const { return m_ident; }
        virtual ~SimpleCluster() {}
        const cluster_graph_t& graph() const { return m_graph; }

        // Non-const access for creators
        cluster_graph_t& graph() { return m_graph; }

       private:
        int m_ident;
        cluster_graph_t m_graph;
    };
}  // namespace WireCell::Aux
