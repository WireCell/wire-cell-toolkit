#ifndef WIRECELLTBB_FUNCTIONCAT
#define WIRECELLTBB_FUNCTIONCAT

#include "WireCellIface/IFunctionNode.h"
#include "WireCellTbb/NodeWrapper.h"

#include <iostream>  // temporary, don't ignore the error code, chump!

namespace WireCellTbb {

    // Body for a TBB function node.
    class FunctionBody {
        WireCell::IFunctionNodeBase::pointer m_wcnode;
        NodeInfo& m_info;

        mutable seqno_t m_seqno{0};

      public:
        FunctionBody(WireCell::INode::pointer wcnode, NodeInfo& info)
            : m_wcnode(std::dynamic_pointer_cast<WireCell::IFunctionNodeBase>(wcnode))
            , m_info(info)
        {
        }
        virtual ~FunctionBody() {}

        msg_t operator()(const msg_t& in) const
        {
            wct_t out;
            m_info.start();
            bool ok = (*m_wcnode)(in.second, out);
            m_info.stop();
            if (!ok) {
                std::cerr << "TbbFlow: function node return false ignored\n";
            }
            return msg_t(m_seqno++, out);
        }
    };

    // Wrap the TBB (compound) node
    class FunctionWrapper : public NodeWrapper {
        tbb::flow::graph_node *m_fn, *m_sn;

      public:

        FunctionWrapper(tbb::flow::graph& graph, WireCell::INode::pointer wcnode)
        {
            m_info.set(wcnode);
            auto fn = new func_node(graph, 1 /*wcnode->concurrency()*/, FunctionBody(wcnode, m_info));
            auto sn = new seq_node(graph, [](const msg_t& m) {return m.first;});
            tbb::flow::make_edge(*fn, *sn);
            m_fn = fn;
            m_sn = sn;
        }

        virtual receiver_port_vector receiver_ports() { return {dynamic_cast<receiver_type*>(m_fn)}; }

        virtual sender_port_vector sender_ports() { return {dynamic_cast<sender_type*>(m_sn)}; }
    };

}  // namespace WireCellTbb

#endif
