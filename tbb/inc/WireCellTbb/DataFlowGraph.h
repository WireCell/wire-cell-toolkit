#ifndef WIRECELLTBB_DATAFLOWGRAPH
#define WIRECELLTBB_DATAFLOWGRAPH

#include "WireCellIface/IDataFlowGraph.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellTbb/NodeWrapper.h"
#include "WireCellTbb/WrapperFactory.h"

#include <map>
#include <string>
#include <unordered_set>

namespace WireCellTbb {

    class DataFlowGraph : public WireCell::Aux::Logger,
                          public WireCell::IDataFlowGraph,
                          public WireCell::IConfigurable
    {
      public:
        DataFlowGraph(int max_threads = 0);
        virtual ~DataFlowGraph();

        /// Connect two nodes so that data runs from tail to head.
        /// Return false on error.
        virtual bool connect(WireCell::INode::pointer tail, WireCell::INode::pointer head, size_t tail_port = 0,
                             size_t head_port = 0);

        /// Run the graph, return false on error.
        virtual bool run();

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

      private:
        tbb::flow::graph m_graph;    // here lives the TBB graph
        WrapperFactory m_factory;

        int m_thread_limit{0};  // 0 means no limit

        // if 0, no summary logged, else log at level 1=debug, 2=info
        int m_summary{1};
        std::unordered_set<WireCellTbb::Node> m_nodes;
    };

}  // namespace WireCellTbb

#endif
