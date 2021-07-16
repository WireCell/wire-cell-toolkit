#ifndef WIRECELLTBB_DATAFLOWGRAPH
#define WIRECELLTBB_DATAFLOWGRAPH

#include "WireCellIface/IDataFlowGraph.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellTbb/NodeWrapper.h"
#include "WireCellTbb/WrapperFactory.h"

#include <map>
#include <string>

namespace WireCellTbb {

    class DataFlowGraph : public WireCell::IDataFlowGraph, public WireCell::IConfigurable {
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
        WireCell::Log::logptr_t l;
        int m_thread_limit{0};  // 0 means no limit
    };

}  // namespace WireCellTbb

#endif
