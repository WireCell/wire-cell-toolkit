#ifndef WIRECELLTBB_FANOUTCAT
#define WIRECELLTBB_FANOUTCAT

#include "WireCellIface/IFanoutNode.h"
#include "WireCellUtil/TupleHelpers.h"
#include "WireCellUtil/Testing.h"
#include "WireCellTbb/NodeWrapper.h"

#include <iostream>

namespace WireCellTbb {

    // Body for a TBB split node.
    template <typename std::size_t N>
    class FanoutBody {
        WireCell::IFanoutNodeBase::pointer m_wcnode;
        NodeInfo& m_info;

        // This is maybe better provided as a std::atomic.  However,
        // the lack of copy constructor conflicts with the need for
        // the body to be copied.  For now, we are safe as we ignore
        // wcnode's desired concurrency and force the parent
        // function_node to have concurrency=1.
        mutable size_t m_seqno{0};

      public:
        typedef typename WireCell::IFanoutNodeBase::any_vector any_vector;
        typedef typename WireCell::type_repeater<N, msg_t>::type TupleType;

        FanoutBody(WireCell::INode::pointer wcnode, NodeInfo& info)
            : m_wcnode(std::dynamic_pointer_cast<WireCell::IFanoutNodeBase>(wcnode))
            , m_info(info)
        {
            Assert(m_wcnode);
        }

        TupleType operator()(msg_t in) const
        {
            any_vector anyvec;
            m_info.start();
            bool ok = (*m_wcnode)(in.second, anyvec);
            m_info.stop();
            if (!ok ) {
                std::cerr << "TbbFlow: fanout call fails\n";
            }
            msg_vector_t savec;
            const size_t seqno = m_seqno++;
            for (auto& a : anyvec) {
                savec.push_back(msg_t(seqno, a));
            }
            return WireCell::vectorToTuple<N>(savec);
        }
    };



    template <std::size_t N>
    sender_port_vector build_fanouter(tbb::flow::graph& graph,
                                      WireCell::INode::pointer wcnode,
                                      std::vector<tbb::flow::graph_node*>& nodes,
                                       NodeInfo& info)
    {
        using TupleType = typename WireCell::type_repeater<N, msg_t>::type;
        using tuple_func_node = tbb::flow::function_node<msg_t, TupleType>;

        // This node takes user WC body and runs it after converting input verctor to tuple.
        // We force-set concurency=1 so as not to upset to seqno in the body.
        auto* fn = new tuple_func_node(graph, 1 /*wcnode->concurrency()*/,
                                       FanoutBody<N>(wcnode, info));
        // Below requires first to be the WCT body caller
        nodes.push_back(fn);

        // this node is fully TBB and splits N sender ports into a tuple
        auto* sp = new tbb::flow::split_node<TupleType>(graph);
        nodes.push_back(sp);

        tbb::flow::make_edge(*fn, *sp);

        sender_port_vector outports;

        // The split_node can introduce out-of-order as it has
        // unlimited concurency.  After each sender_node output port
        // we add a sequencer based on the seqno.
        for (sender_type* one : sender_ports<TupleType>(*sp)) {

            auto* seq = new seq_node(graph, [](const msg_t& sa) { return sa.first; });
            nodes.push_back(seq);
            tbb::flow::make_edge(*one, *seq);

            outports.push_back(dynamic_cast<sender_type*>(seq));
        }
        return outports;
    }

    // Wrap the TBB (compound) node
    class FanoutWrapper : public NodeWrapper {
        std::vector<tbb::flow::graph_node*> m_nodes;
        sender_port_vector m_sender_ports;

      public:
        FanoutWrapper(tbb::flow::graph& graph, WireCell::INode::pointer wcnode)
        {
            m_info.set(wcnode);
            int nout = wcnode->output_types().size();
            // an exhaustive switch to convert from run-time to compile-time types and enumerations.
            Assert(nout > 0 && nout <= 10);  // fixme: exception instead?
            if (1 == nout) m_sender_ports = build_fanouter<1>(graph, wcnode, m_nodes, m_info);
            if (2 == nout) m_sender_ports = build_fanouter<2>(graph, wcnode, m_nodes, m_info);
            if (3 == nout) m_sender_ports = build_fanouter<3>(graph, wcnode, m_nodes, m_info);
            if (4 == nout) m_sender_ports = build_fanouter<4>(graph, wcnode, m_nodes, m_info);
            if (5 == nout) m_sender_ports = build_fanouter<5>(graph, wcnode, m_nodes, m_info);
            if (6 == nout) m_sender_ports = build_fanouter<6>(graph, wcnode, m_nodes, m_info);
            if (7 == nout) m_sender_ports = build_fanouter<7>(graph, wcnode, m_nodes, m_info);
            if (8 == nout) m_sender_ports = build_fanouter<8>(graph, wcnode, m_nodes, m_info);
            if (9 == nout) m_sender_ports = build_fanouter<9>(graph, wcnode, m_nodes, m_info);
            if (10 == nout) m_sender_ports = build_fanouter<10>(graph, wcnode, m_nodes, m_info);
        }

        virtual sender_port_vector sender_ports()
        {
            return m_sender_ports;
        }

        virtual receiver_port_vector receiver_ports()
        {
            auto ptr = dynamic_cast<receiver_type*>(m_nodes[0]);
            return receiver_port_vector{ptr};
        }
    };

}  // namespace WireCellTbb

#endif
