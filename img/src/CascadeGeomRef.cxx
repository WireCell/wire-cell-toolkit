// The reference cross-slice pair builder of the cascade level graph (wcfm doc 14): the same pairs as
// Aux::Cascade::geom_pairs, read back from a cluster graph filled by Img::geom_clustering itself.  Test-only
// (doctest_cascade_deghosting); it stays in img because img owns geom_clustering (pdvd doc 130 phase 2 moved the
// rest of CascadeGraph.cxx to aux).

#include "WireCellImg/CascadeGraph.h"
#include "WireCellImg/GeomClusteringUtil.h"

#include "WireCellAux/SimpleBlob.h"

#include <algorithm>
#include <unordered_map>

using namespace WireCell;
using namespace WireCell::Img;

std::vector<std::array<int64_t, 2>> Cascade::geom_pairs_graph(const std::vector<IBlob::pointer>& blobs,
                                                              const std::vector<int>& sidx,
                                                              const std::vector<ISlice::pointer>& slice_of,
                                                              const std::string& policy)
{
    const size_t N = blobs.size();
    std::vector<int> set_slice;
    const auto per = Aux::Cascade::slice_sets(sidx, slice_of.size(), set_slice);
    IBlobSet::vector sets;
    for (size_t k = 0; k < per.size(); ++k) {
        IBlob::vector v;
        for (auto i : per[k]) v.push_back(blobs[i]);
        sets.push_back(std::make_shared<Aux::SimpleBlobSet>(set_slice[k], slice_of[set_slice[k]], v));
    }
    cluster_indexed_graph_t grind;
    for (auto it = sets.begin(); it != sets.end(); ++it) {
        Img::geom_clustering(grind, it, sets.end(), policy);
    }
    std::unordered_map<const IBlob*, int64_t> bidx;  // lookup only
    bidx.reserve(N);
    for (size_t i = 0; i < N; ++i) bidx[blobs[i].get()] = (int64_t) i;
    std::vector<std::array<int64_t, 2>> bb;
    const auto& g = grind.graph();
    for (auto e : boost::make_iterator_range(boost::edges(g))) {
        const auto& na = g[boost::source(e, g)];
        const auto& nb = g[boost::target(e, g)];
        if (na.code() != 'b' || nb.code() != 'b') continue;
        const int64_t a = bidx.at(std::get<IBlob::pointer>(na.ptr).get());
        const int64_t b = bidx.at(std::get<IBlob::pointer>(nb.ptr).get());
        if (a == b) continue;
        bb.push_back({std::min(a, b), std::max(a, b)});
    }
    std::sort(bb.begin(), bb.end());
    bb.erase(std::unique(bb.begin(), bb.end()), bb.end());
    return bb;
}

