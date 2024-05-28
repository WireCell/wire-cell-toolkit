#ifndef WIRECELLTBB_SINKCAT
#define WIRECELLTBB_SINKCAT

#include "WireCellTbb/NodeWrapper.h"
#include "WireCellIface/ISinkNode.h"

#include <iostream> // fixme: for non-error handling

namespace WireCellTbb {

    // adapter to convert from WC sink node to TBB sink node body.
    class SinkBody {
        WireCell::ISinkNodeBase::pointer m_wcnode;
        NodeInfo& m_info;

      public:
        ~SinkBody() {}

        SinkBody(WireCell::INode::pointer wcnode, NodeInfo& info)
            : m_wcnode(std::dynamic_pointer_cast<WireCell::ISinkNodeBase>(wcnode))
            , m_info(info)
        {
              
        }
        tbb::flow::continue_msg operator()(const msg_t& in)
        {
            m_info.start();
            bool ok = (*m_wcnode)(in.second);
            m_info.stop();
            if (!ok) {
                std::cerr << "TbbFlow: sink node return false ignored\n";
            }
            return {};
        }
    };

    // implement facade to access ports for sink nodes
    class SinkNodeWrapper : public NodeWrapper {
        sink_node* m_tbbnode;

      public:
        SinkNodeWrapper(tbb::flow::graph& graph, WireCell::INode::pointer wcnode)
            : m_tbbnode(new sink_node(graph, 1 /*wcnode->concurrency()*/, SinkBody(wcnode, m_info)))
        {
            m_info.set(wcnode);
        }
        ~SinkNodeWrapper() {
            delete m_tbbnode;
        }
        virtual receiver_port_vector receiver_ports()
        {
            auto ptr = dynamic_cast<receiver_type*>(m_tbbnode);
            return receiver_port_vector{ptr};
        }
    };
}  // namespace WireCellTbb
#endif
