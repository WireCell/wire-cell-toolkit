/** This adapts the queued out node category to a TBB multifunction node.
 */

#ifndef WIRECELLTBB_QUEUEDOUT
#define WIRECELLTBB_QUEUEDOUT

#include "WireCellTbb/NodeWrapper.h"
#include "WireCellIface/IQueuedoutNode.h"

#include <iostream>             // debug

namespace WireCellTbb {

    class QueuedoutBody {
        WireCell::IQueuedoutNodeBase::pointer m_wcnode;
        NodeInfo& m_info;

        seqno_t m_seqno{0};

      public:
        using mfunc_port = mfunc_node::output_ports_type;

        ~QueuedoutBody() {}
        QueuedoutBody(WireCell::INode::pointer wcnode, NodeInfo& info)
            : m_wcnode(std::dynamic_pointer_cast<WireCell::IQueuedoutNodeBase>(wcnode))
            , m_info(info)
        {
        }
        void operator()(const msg_t& in, mfunc_port& out)
        {
            WireCell::IQueuedoutNodeBase::queuedany outq;
            m_info.start();
            bool ok = (*m_wcnode)(in.second, outq);
            m_info.stop();
            if (!ok) {
                std::cerr << "TbbFlow: queuedout node return false ignored\n";
                return;
            }
            for (auto a : outq) {
                // does not block.
                bool accepted = std::get<0>(out).try_put(msg_t(m_seqno++, a));
                if (!accepted) {
                    std::cerr << "TbbFlow: unaccepted try_put " << m_wcnode->signature() << std::endl;
                }
            }
        }
    };

    class QueuedoutWrapper : public NodeWrapper {
        tbb::flow::graph_node *m_fn, *m_qn;

      public:

        QueuedoutWrapper(tbb::flow::graph& graph, WireCell::INode::pointer wcnode)
        {
            m_info.set(wcnode);
            auto fn = new mfunc_node(graph, wcnode->concurrency(), QueuedoutBody(wcnode, m_info));
            auto qn = new seq_node(graph, [](const msg_t& m) {return m.first;});
            tbb::flow::make_edge(*fn, *qn);
            m_fn = fn;
            m_qn = qn;
        }
        virtual ~QueuedoutWrapper()
        {
            /// we haven't been doing this religiously....
            // delete m_tbbnode;
            // m_tbbnode = nullptr;
        }

        virtual receiver_port_vector receiver_ports()
        {
            auto ptr = dynamic_cast<receiver_type*>(m_fn);
            return receiver_port_vector{ptr};
        }

        virtual sender_port_vector sender_ports()
        {
            auto ptr = dynamic_cast<sender_type*>(m_qn);
            return sender_port_vector{ptr};
        }

    };

}  // namespace WireCellTbb
#endif
