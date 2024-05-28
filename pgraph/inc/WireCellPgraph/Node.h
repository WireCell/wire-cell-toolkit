#ifndef WIRECELL_PGRAPH_NODE
#define WIRECELL_PGRAPH_NODE

#include "WireCellPgraph/Port.h"

namespace WireCell {
    namespace Pgraph {

        // A node in the DFP graph must inherit from Node.
        class Node {
           public:
            Node() {}  // constructures may wish to resize/populate m_ports.
            virtual ~Node() {}

            // Concrete Node must implement this to consume inputs
            // and/or produce outputs.
            virtual bool operator()() = 0;

            // Concrete node must return some instance identifier.
            virtual std::string ident() = 0;

            Port& iport(size_t ind = 0) { return port(Port::input, ind); }
            Port& oport(size_t ind = 0) { return port(Port::output, ind); }

            PortList& input_ports() { return m_ports[Port::input]; }
            PortList& output_ports() { return m_ports[Port::output]; }

            Port& port(Port::Type type, size_t ind = 0)
            {
                if (ind >= m_ports[type].size()) {
                    THROW(ValueError() << errmsg{"unknown port"});
                }
                return m_ports[type][ind];
            }
            Port& port(Port::Type type, const std::string& name)
            {
                for (size_t ind = 0; ind < m_ports[type].size(); ++ind) {
                    if (m_ports[type][ind].name() != name) {
                        continue;
                    }
                    return port(type, ind);
                }
                THROW(ValueError() << errmsg{"unknown port"});
            }

            bool connected()
            {
                for (auto& p : input_ports()) {
                    if (!p.edge()) {
                        return false;
                    }
                }
                for (auto& p : output_ports()) {
                    if (!p.edge()) {
                        return false;
                    }
                }
                return true;
            }

            std::vector<Port> disconnected_ports()
            {
                std::vector<Port> ret;
                for (auto& p : input_ports()) {
                    if (!p.edge()) {
                        ret.push_back(p);
                    }
                }
                for (auto& p : output_ports()) {
                    if (!p.edge()) {
                        ret.push_back(p);
                    }
                }
                return ret;
            }

           protected:
            // Concrete class should fill during construction
            PortList m_ports[Port::ntypes];
        };
    }  // namespace Pgraph
}  // namespace WireCell

#endif
