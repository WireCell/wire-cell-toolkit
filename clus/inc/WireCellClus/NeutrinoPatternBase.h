#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"

namespace WireCell::Clus::PR {
    class PatternAlgorithms{
        public:
        std::set<VertexPtr> find_vertices(Graph& graph, const Facade::Cluster& cluster);
        std::set<SegmentPtr> find_segments(Graph& graph, const Facade::Cluster& cluster);
        bool clean_up_graph(Graph& graph, const Facade::Cluster& cluster);
    };
}
