#ifndef WIRECELLTBB_DATAFLOWGRAPH
#define WIRECELLTBB_DATAFLOWGRAPH

#include "WireCellIface/IDataFlowGraph.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellTbb/NodeWrapper.h"
#include "WireCellTbb/WrapperFactory.h"

#include <tbb/task_arena.h>

#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

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
        // The TBB graph (and the wrapper factory that creates its nodes) are
        // built lazily so they can be constructed INSIDE a sized task_arena when
        // a thread limit is configured.  A tbb::flow::graph binds to the arena of
        // its construction, so limiting node parallelism requires building,
        // connecting, and running the graph within that arena -- not merely
        // wrapping wait_for_all().  See spng-fra.12 and test_dfp_thread_limit.
        std::unique_ptr<tbb::task_arena> m_arena;   // null => inherit host arena
        std::unique_ptr<tbb::flow::graph> m_graph;
        std::unique_ptr<WrapperFactory> m_factory;

        int m_thread_limit{0};  // 0 means no limit (inherit host/default arena)

        // if 0, no summary logged, else log at level 1=debug, 2=info
        int m_summary{1};
        std::unordered_set<WireCellTbb::Node> m_nodes;

        // Lazily build arena (if limited) + graph + factory, once.
        void ensure();

        // The body of connect(), run inside the sized arena so the TBB nodes it
        // creates belong to that arena.
        bool do_connect(WireCell::INode::pointer tail, WireCell::INode::pointer head,
                        size_t tail_port, size_t head_port);

        // Run f() inside the sized arena if one exists, else in the current arena.
        template <class F>
        void in_arena(F&& f)
        {
            if (m_arena) { m_arena->execute(std::forward<F>(f)); }
            else { f(); }
        }
    };

}  // namespace WireCellTbb

#endif
