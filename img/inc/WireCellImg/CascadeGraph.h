/** Forwarding header: the cascade level graph moved to WireCellAux/CascadeGraph.h (pdvd doc 130 phase 2).
    Img::Cascade keeps every name of Aux::Cascade through the using-directive (qualified lookup follows it), plus
    geom_pairs_graph, the test reference that needs Img::geom_clustering and therefore stays in this plugin.
 */
#ifndef WIRECELLIMG_CASCADEGRAPH
#define WIRECELLIMG_CASCADEGRAPH

#include "WireCellAux/CascadeGraph.h"

namespace WireCell::Img::Cascade {
    using namespace WireCell::Aux::Cascade;

    /// The pairs geom_pairs computes, read back from a cluster graph filled by Img::geom_clustering (reference).
    std::vector<std::array<int64_t, 2>> geom_pairs_graph(const std::vector<IBlob::pointer>& blobs,
                                                         const std::vector<int>& sidx,
                                                         const std::vector<ISlice::pointer>& slice_of,
                                                         const std::string& policy);
}

#endif
