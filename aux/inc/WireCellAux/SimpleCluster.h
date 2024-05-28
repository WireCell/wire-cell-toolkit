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
