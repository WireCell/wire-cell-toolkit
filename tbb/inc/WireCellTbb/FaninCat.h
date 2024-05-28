#ifndef WIRECELLTBB_FANINCAT
#define WIRECELLTBB_FANINCAT

#include "WireCellIface/IFaninNode.h"
#include "WireCellUtil/TupleHelpers.h"
#include "WireCellUtil/Testing.h"
#include "WireCellTbb/NodeWrapper.h"

#include <iostream> // fixme: for non-error handling

namespace WireCellTbb {

    // Body for a TBB join node.
    template <typename TupleType>
    class FaninBody {
        WireCell::IFaninNodeBase::pointer m_wcnode;
        NodeInfo& m_info;

        mutable seqno_t m_seqno{0};

      public:
        typedef typename WireCell::IFaninNodeBase::any_vector any_vector;
        typedef typename WireCell::tuple_helper<TupleType> helper_type;

        FaninBody(WireCell::INode::pointer wcnode, NodeInfo& info)
            : m_wcnode(std::dynamic_pointer_cast<WireCell::IFaninNodeBase>(wcnode))
            , m_info(info)
        {
            Assert(m_wcnode);
        }

        msg_t operator()(const TupleType& tup) const
        {
            auto msg_vec = as_msg_vector(tup);
            any_vector in;
            for (auto& msg : msg_vec) {
                in.push_back(msg.second);
            }
            wct_t out;
            m_info.start();
            bool ok = (*m_wcnode)(in, out);
            m_info.stop();
            if (!ok) {
                std::cerr << "TbbFlow: fanin node return false ignored\n";
            }
            return msg_t(m_seqno++, out);
        }
    };

    template <std::size_t N>
    receiver_port_vector build_faniner(tbb::flow::graph& graph,
                                       WireCell::INode::pointer wcnode,
                                       std::vector<tbb::flow::graph_node*>& nodes,
                                       NodeInfo& info)
    {
        typedef typename WireCell::type_repeater<N, msg_t>::type TupleType;

        // this node is fully TBB and joins N receiver ports into a tuple
        typedef tbb::flow::join_node<TupleType> tbb_join_node_type;
        tbb_join_node_type* jn = new tbb_join_node_type(graph);

        // this node takes user WC body and runs it after converting input tuple to vector
        typedef tbb::flow::function_node<TupleType, msg_t> joining_node;
        joining_node* fn = new joining_node(graph, 1 /*wcnode->concurrency()*/,
                                            FaninBody<TupleType>(wcnode, info));

        tbb::flow::make_edge(*jn, *fn);

        auto* sn = new seq_node(graph, [](const msg_t& m) {return m.first;});
        tbb::flow::make_edge(*fn, *sn);

        nodes.push_back(sn);    // first
        nodes.push_back(fn);
        nodes.push_back(jn);
        return receiver_ports(*jn);
    }

    // Wrap the TBB (compound) node
    class FaninWrapper : public NodeWrapper {
        std::vector<tbb::flow::graph_node*> m_nodes;
        receiver_port_vector m_receiver_ports;

      public:
        FaninWrapper(tbb::flow::graph& graph, WireCell::INode::pointer wcnode)
        {
            m_info.set(wcnode);
            int nin = wcnode->input_types().size();
            // an exhaustive switch to convert from run-time to compile-time types and enumerations.
            Assert(nin > 0 && nin <= 10);  // fixme: exception instead?
            if (1 == nin) m_receiver_ports = build_faniner<1>(graph, wcnode, m_nodes, m_info);
            if (2 == nin) m_receiver_ports = build_faniner<2>(graph, wcnode, m_nodes, m_info);
            if (3 == nin) m_receiver_ports = build_faniner<3>(graph, wcnode, m_nodes, m_info);
            if (4 == nin) m_receiver_ports = build_faniner<4>(graph, wcnode, m_nodes, m_info);
            if (5 == nin) m_receiver_ports = build_faniner<5>(graph, wcnode, m_nodes, m_info);
            if (6 == nin) m_receiver_ports = build_faniner<6>(graph, wcnode, m_nodes, m_info);
            if (7 == nin) m_receiver_ports = build_faniner<7>(graph, wcnode, m_nodes, m_info);
            if (8 == nin) m_receiver_ports = build_faniner<8>(graph, wcnode, m_nodes, m_info);
            if (9 == nin) m_receiver_ports = build_faniner<9>(graph, wcnode, m_nodes, m_info);
            if (10 == nin) m_receiver_ports = build_faniner<10>(graph, wcnode, m_nodes, m_info);
        }

        virtual receiver_port_vector receiver_ports() { return m_receiver_ports; }

        virtual sender_port_vector sender_ports()
        {
            auto ptr = dynamic_cast<sender_type*>(m_nodes[0]);
            return sender_port_vector{ptr};
        }
    };

}  // namespace WireCellTbb

#endif
