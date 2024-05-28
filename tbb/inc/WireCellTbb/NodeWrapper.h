#ifndef WIRECELLTBB_NODEWRAPPER
#define WIRECELLTBB_NODEWRAPPER

#include "WireCellIface/INode.h"
#include "WireCellIface/INamed.h"
#include "WireCellUtil/TupleHelpers.h"

#include <tbb/flow_graph.h>
#include <boost/any.hpp>
#include <iostream>
#include <utility>              // make_index_sequence
#include <string>
#include <memory>
#include <chrono>
#include <map>

namespace WireCellTbb {

    // Message type passed through WCT nodes
    using wct_t = boost::any;

    // We combine WCT data with a sequence number to assure order.
    using seqno_t = size_t;

    // The message type used by all TBB flow_graph nodes
    using msg_t = std::pair<seqno_t, wct_t>;
    using msg_vector_t = std::vector<msg_t>;

    // Broken out sender/receiver types, vectors of their pointers
    typedef tbb::flow::sender<msg_t> sender_type;
    typedef tbb::flow::receiver<msg_t> receiver_type;

    typedef std::vector<sender_type*> sender_port_vector;
    typedef std::vector<receiver_type*> receiver_port_vector;

    // Each WCT node is serviced by a small subgraph of TBB nodes.
    using src_node = tbb::flow::input_node<msg_t>;
    using func_node = tbb::flow::function_node<msg_t, msg_t>;
    using mfunc_node = tbb::flow::multifunction_node<msg_t, std::tuple<msg_t>>;
    using seq_node = tbb::flow::sequencer_node<msg_t>;
    using sink_node = tbb::flow::function_node<msg_t>;

    // tuple type nodes include join_node, split_node and indexer_node

    // A helper to provide info about the node.
    class NodeInfo {
      public:

        void set(WireCell::INode::pointer wcnode) { m_inode = wcnode; }

        WireCell::INode::pointer inode() const { return m_inode; }

        // Return the instance name, only works if INode is also an
        // INamed.
        std::string instance_name() const {
            auto inamed = std::dynamic_pointer_cast<WireCell::INamed>(m_inode);
            if (inamed) {
                return inamed->get_name();
            }
            return "(unknown)";
        }

        // Start/stop the stop watch.
        void start() {
            m_clock = std::chrono::high_resolution_clock::now();
        }
        void stop() {
            duration_t delta = std::chrono::high_resolution_clock::now() - m_clock;
            m_runtime += delta;
            if (delta > m_maxrt) {
                m_maxrt = delta;
            }
            ++m_calls;
        }

        //using duration_t = std::chrono::high_resolution_clock::duration;
        using duration_t = std::chrono::duration<double>;

        // Total runtime for all executions.
        duration_t runtime() const {
            return m_runtime;
        };
        // Maximum runtime for any execution.
        duration_t max_runtime() const {
            return m_maxrt;
        };

        size_t calls() const {
            return m_calls;
        }

      private:
        WireCell::INode::pointer m_inode;
        duration_t m_runtime {0}, m_maxrt{0};
        using clock_t = std::chrono::high_resolution_clock;
        using time_point_t = clock_t::time_point;
        time_point_t m_clock;

        size_t m_calls{0};
    };
    std::ostream& operator<<(std::ostream& os, const NodeInfo& info);

    // A base facade which expose sender/receiver ports and provide
    // initialize hook.  There is one NodeWrapper for each node
    // category.
    class NodeWrapper {
       public:
        virtual ~NodeWrapper() {}

        virtual sender_port_vector sender_ports() { return sender_port_vector(); }
        virtual receiver_port_vector receiver_ports() { return receiver_port_vector(); }

        // call before running graph
        virtual void initialize() {}

        const NodeInfo& info() const { return m_info; };
      protected:
        NodeInfo m_info;
    };

    // expose the wrappers only as a shared pointer
    typedef std::shared_ptr<NodeWrapper> Node;


    // tuple helpers

    template <typename Tuple, std::size_t... Is>
    msg_vector_t as_msg_vector(const Tuple& tup, std::index_sequence<Is...>)
    {
        return {std::get<Is>(tup)...};
    }
    template <typename Tuple>
    msg_vector_t as_msg_vector(const Tuple& tup)
    {
        return as_msg_vector(tup, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
    }


    // join
    template <typename Tuple, std::size_t... Is>
    receiver_port_vector receiver_ports(tbb::flow::join_node<Tuple>& jn, std::index_sequence<Is...>)
    {
        return {dynamic_cast<receiver_type*>(&tbb::flow::input_port<Is>(jn))...};
    }
    /// Return receiver ports of a join node as a vector.
    template <typename Tuple>
    receiver_port_vector receiver_ports(tbb::flow::join_node<Tuple>& jn)
    {
        return receiver_ports(jn, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
    }


    // split
    template <typename Tuple, std::size_t... Is>
    sender_port_vector sender_ports(tbb::flow::split_node<Tuple>& sp, std::index_sequence<Is...>)
    {
        return {dynamic_cast<sender_type*>(&tbb::flow::output_port<Is>(sp))...};
    }
    /// Return sender ports of a split node as a vector.
    template <typename Tuple>
    sender_port_vector sender_ports(tbb::flow::split_node<Tuple>& sp)
    {
        return sender_ports(sp, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
    }



}  // namespace WireCellTbb

#endif
