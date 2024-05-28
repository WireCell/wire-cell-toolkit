#ifndef WIRECELLTBB_JOINCAT
#define WIRECELLTBB_JOINCAT

#include "WireCellIface/IJoinNode.h"
#include "WireCellUtil/Testing.h"
#include "WireCellTbb/NodeWrapper.h"

#include <iostream> // fixme: for non-error handling

namespace WireCellTbb {

    // Body for a TBB join node.
    template <typename TupleType>
    class JoinBody {
        WireCell::IJoinNodeBase::pointer m_wcnode;
        NodeInfo& m_info;

        mutable seqno_t m_seqno{0};

       public:
        typedef typename WireCell::IJoinNodeBase::any_vector any_vector;
        typedef typename WireCell::tuple_helper<TupleType> helper_type;

        JoinBody(WireCell::INode::pointer wcnode, NodeInfo& info)
            : m_wcnode(std::dynamic_pointer_cast<WireCell::IJoinNodeBase>(wcnode))
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
                std::cerr << "TbbFlow: join node return false ignored\n";
            }
            return msg_t(m_seqno++, out);
        }
    };

    template <std::size_t N>
    receiver_port_vector build_joiner(tbb::flow::graph& graph,
                                      WireCell::INode::pointer wcnode,
                                      std::vector<tbb::flow::graph_node*>& nodes,
                                      NodeInfo& info)
    {
        typedef typename WireCell::type_repeater<N, msg_t>::type TupleType;

        // this node takes user WC body and runs it after converting input tuple to vector
        typedef tbb::flow::function_node<TupleType, msg_t> joining_node;
        auto* fn = new joining_node(graph, 1 /*wcnode->concurrency()*/,
                                    JoinBody<TupleType>(wcnode, info));

        // this node is fully TBB and joins N receiver ports into a tuple
        typedef tbb::flow::join_node<TupleType> tbb_join_node_type;
        auto* jn = new tbb_join_node_type(graph);

        tbb::flow::make_edge(*jn, *fn);

        auto* sn = new seq_node(graph, [](const msg_t& m) {return m.first;});
        tbb::flow::make_edge(*fn, *sn);

        nodes.push_back(sn);    // first
        nodes.push_back(jn);
        nodes.push_back(fn);
        return receiver_ports(*jn);
    }

    // Wrap the TBB (compound) node
    class JoinWrapper : public NodeWrapper {
        std::vector<tbb::flow::graph_node*> m_nodes;
        receiver_port_vector m_receiver_ports;

       public:
        JoinWrapper(tbb::flow::graph& graph, WireCell::INode::pointer wcnode)
        {
            m_info.set(wcnode);
            int nin = wcnode->input_types().size();
            // an exhaustive switch to convert from run-time to compile-time types and enumerations.
            Assert(nin > 0 && nin <= 3);  // fixme: exception instead?
            if (1 == nin) m_receiver_ports = build_joiner<1>(graph, wcnode, m_nodes, m_info);
            if (2 == nin) m_receiver_ports = build_joiner<2>(graph, wcnode, m_nodes, m_info);
            if (3 == nin) m_receiver_ports = build_joiner<3>(graph, wcnode, m_nodes, m_info);
        }

        virtual receiver_port_vector receiver_ports() {
            return m_receiver_ports;
        }

        virtual sender_port_vector sender_ports()
        {
            auto ptr = dynamic_cast<sender_type*>(m_nodes[0]);
            return sender_port_vector{ptr};
        }
    };

}  // namespace WireCellTbb

#endif
