/** Blob shadow.

    A "blob shadow" measures the amount of an "overlap" shared by a
    pair blobs.

    There are two definitions for shadow overlap which depend on the
    basis for their calculation.  They are qualifed either as "wire"
    or "channel" shadows.  For detectors with a single wire (segment)
    per channel (eg, no wrapped wires) the two are degenerate.  Two
    blobs may have a channel shadow because they have channels in
    common but may not have any shared wire segments so would lack any
    wire shadow.  The contrary may not arise. A pair of blobs with a
    wire shadow must also have a channel shadow.

    A shadow is characterized by the WirePlaneId and by its extent in
    the pitch direction as bound by extreme channel or wire indices.

    For a channel shadow, the WirePlaneId::face() value is not material
    due to wire wrapping.  The channel shadow bounds are given by a pair of
    wire-attachment-numbers such as as returned by IChannel::index().

    For a wire shadow, the WirePlaneId::face() value is meaningful
    even in for detectors with wrapped wires.  The wire shadow bounds
    are represented by a pair of wire-in-plane number such as returned
    by IWire::index().
    
    The collection of shadows is represented by a "channel shadow
    graph" or a "wire shadow graph" with each node representing a blob
    and an edge representing their shadow.  The edge provides
    parameters representing WirePlaneId and shadow extent.  A shadow
    graph holds a graph-level property which is a "char" giving the
    code 'w' if it holds wire shadows or 'c' if it holds channel
    shadows.

*/
#ifndef WIRECELL_AUX_BLOBSHADOW
#define WIRECELL_AUX_BLOBSHADOW

#include "WireCellIface/ICluster.h"
#include "WireCellIface/WirePlaneId.h"

namespace WireCell::Aux::BlobShadow {

    struct Graph {
        // the shadow type code w=wire, c=channel
        char stype;
    };

    struct Node {
        // Vertex descriptor into original cluster graph.
        cluster_vertex_t desc;
    };

    struct Edge {
        // The bounds of the shadow in terms of
        // wire-atachment-numbers or wire-in-plane indices.
        // Reminder, "end" is one past the range so that "end-beg"
        // gives the width of the shadow in units of channels or
        // wires.
        int beg{0}, end{0};

        // The wire plane in which the shadow exists.  The
        // wpid.face() value is only relevant for a wire shadow.
        WirePlaneId wpid{0};
    };

    using graph_t = boost::adjacency_list<boost::setS, boost::vecS,
                                          boost::undirectedS,
                                          Node, Edge, Graph>;
    using vdesc_t = boost::graph_traits<graph_t>::vertex_descriptor;
    using edesc_t = boost::graph_traits<graph_t>::edge_descriptor;


    // Return a shadow graph for either a wire-shadow (code='w') or a
    // channel-shadow (code='c').
    graph_t shadow(const cluster_graph_t& cgraph, char code);
     
}

#endif
