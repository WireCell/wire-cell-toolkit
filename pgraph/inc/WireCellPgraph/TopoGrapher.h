/** A TopoGrapher is a drop-in alternative to Pgrapher.

    It runs the exact same data-flow graph with the exact same execution
    policy -- always fire the sink-most runnable node so that no more than
    one datum is ever held on an edge (minimizing "in flight" data) -- and
    so fires the identical sequence of nodes and produces identical output.

    The only difference from Pgrapher is the scheduling algorithm.  Pgrapher
    rescans the full topological node order after every single firing, which
    costs O(V^2) and becomes pathological for graphs with many nodes upstream
    of the sinks (e.g. wide fan-out trees).  TopoGrapher instead maintains a
    rank-ordered ready set and updates it incrementally as nodes fire, giving
    O((V+E) log V).

    Configuration is identical to Pgrapher (an "edges" list and optional
    "verbosity").  See Pgrapher.h for the configuration format.  Pgrapher is
    retained unchanged; choose the engine via the `wire-cell -a` application
    name and the config component `type`.
 */

#ifndef WIRECELL_PGRAPH_TOPOGRAPHER
#define WIRECELL_PGRAPH_TOPOGRAPHER

#include "WireCellIface/IApplication.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellPgraph/Graph.h"

namespace WireCell::Pgraph {

    class TopoGrapher : public Aux::Logger,
                        public IApplication,
                        public IConfigurable
    {
      public:
        TopoGrapher();
        virtual ~TopoGrapher();

        // IApplication
        virtual void execute();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

      private:
        Graph m_graph;
        int m_verbosity{1};
    };

}  // namespace WireCell::Pgraph

#endif
