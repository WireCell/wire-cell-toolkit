/** A Pgrapher is an application that runs components via a
    user-configured directect acyclic graph following wire-cell data
    flow processing paradigm.

    The individual components must be configured on their own as
    usual.

    The graph is specified as a list of edges.  Each edge has a tail
    (source) and a head (destination) endpoint.  An endpoint is
    specified by a "typename" and an optional "port".  As usual, the
    typename is in the form "type:name" or just "type" and as
    typically constructed by the Jsonnet function wc.tn() from a
    previously defined configuration object.  The port is an index an
    integer default to 0.

    A configuration might look like:

    local cosmics = { type:"TrackDepos", name:"cosmics", ...};
    [ // main configuration sequence which includes config for
      // TrackDepos:cosmics, TrackDepos:beam,
      // DepoJoiner (which may not yet exist) but not
      // DumpDepos as it is not a configurable....,
    cosmics, beam, ...,
    {
    type:"Pgrapher",
    data:{
      edges:[
        {
          tail:{node:wc.tn(cosmics)},
          head:{node:wc.tn(joiner), port:0}
        },
        {
          tail:{node:wc.tn(beam)},
          head:{node:"DepoJoiner", port:1}
        },
        {
          tail:{node:"DepoJoiner"},
          head:{node:"DumpDepos"}
        },
      ],
    }}]

    Note the port number need really only be specified in the second
    edge in order to send the data to port #1.  port #0 is default.

    As when configuraing a component itself, the name need only be
    specified in an edge pair if not using the default (empty string).

    A "verbosity" sets how much info is printed about the graph execution.  0 is
    none, 1 is default and gives summary of time, 2 also includes ExecMon
    tracing.

 */

#ifndef WIRECELL_PGRAPH_PGRAPHER
#define WIRECELL_PGRAPH_PGRAPHER

#include "WireCellIface/IApplication.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellPgraph/Graph.h"

namespace WireCell::Pgraph {

    class Pgrapher : public Aux::Logger,
                     public IApplication,
                     public IConfigurable
    {
      public:
        Pgrapher();
        virtual ~Pgrapher();

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
