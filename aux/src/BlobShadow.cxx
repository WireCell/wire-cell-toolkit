#include "WireCellAux/BlobShadow.h"
#include "WireCellIface/ICluster.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellUtil/Graph.h"

#include <unordered_map>

using namespace WireCell;
using namespace WireCell::GraphTools;
using namespace WireCell::Aux;
using namespace WireCell::Aux::BlobShadow;

// #include <iostream>             // debug
// template <typename Gr>
// void dump(typename boost::graph_traits<Gr>::vertex_descriptor vtx, const Gr& gr, const std::string& msg="")
// {
//     std::cerr << "vtx:" << vtx << " " << gr[vtx].code() << gr[vtx].ident() << " " << msg << "\n";
// }

namespace {
    // Walk b->w(->c) on the undirected cluster graph directly.  The
    // former directed copy (type_directed) only served to stop BFS from
    // crawling back up; every hop here already filters its target by
    // type code, which subsumes that.  Order-identical to the directed
    // walk: setS out-edge sets order by target descriptor in both
    // graphs, and type_directed preserved vertex descriptors.
    void connected_leaves(std::vector<cluster_vertex_t> & leaves,
    const cluster_graph_t & cgraph, const cluster_vertex_t & bvtx, char leaf_code) {
        // leaf_code can be 'w' or 'c'
        std::unordered_set<char> valid_codes = {'w', 'c'};
        if (valid_codes.find(leaf_code) == valid_codes.end()) {
            // TODO: make some noise?
            return;
        }
        for(auto bedge : mir(boost::out_edges(bvtx, cgraph))) {
            auto wvtx = boost::target(bedge, cgraph);
            if (cgraph[wvtx].code() != 'w') {
                continue;
            }
            if (leaf_code == 'w') {
                leaves.push_back(wvtx);
                continue;
            }
            // if not 'w', find 'c'
            for(auto wedge : mir(boost::out_edges(wvtx, cgraph))) {
                auto cvtx = boost::target(wedge, cgraph);
                if (cgraph[cvtx].code() != 'c') {
                    continue;
                }
                leaves.push_back(wvtx);
            }
        }
    }

    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator()(std::pair<T1, T2> const& pair) const
        {
            std::size_t h1 = std::hash<T1>()(pair.first);
            std::size_t h2 = std::hash<T2>()(pair.second);
            return h1 ^ h2;
        }
    };

    // Dedup index for (blob pair, layer) -> edge slot.  Replaces the
    // former edge_range() scan, which required multisetS out-edges and
    // was linear in the parallel-edge count.  The pair key is
    // order-normalized; shadows are undirected so direction never
    // mattered.
    using layer_edge_key_t = std::pair<uint64_t, int>;
    using layer_edge_map_t = std::unordered_map<layer_edge_key_t, size_t, pair_hash>;

    layer_edge_key_t layer_edge_key(BlobShadow::vdesc_t v1, BlobShadow::vdesc_t v2,
                                    WirePlaneLayer_t layer)
    {
        if (v2 < v1) { std::swap(v1, v2); }
        return {(static_cast<uint64_t>(v1) << 32) | static_cast<uint64_t>(v2),
                static_cast<int>(layer)};
    }
}  // namespace

BlobShadow::Shadows BlobShadow::shadow_list(const cluster_graph_t& cgraph, char leaf_code)
{
    BlobShadow::Shadows shadows; // will return
    shadows.stype = leaf_code;
    layer_edge_map_t layer_edges;

    // Loop over blobs, to load up output nodes, old->new map.
    std::unordered_map<cluster_vertex_t, BlobShadow::vdesc_t> c2bs;
    for (auto bvtx : mir(boost::vertices(cgraph))) {
        if (cgraph[bvtx].code() == 'b') {
            c2bs[bvtx] = shadows.nodes.size();
            shadows.nodes.push_back({bvtx});
        }
    }

    // Loop again, to pick up slices and make edges
    for (auto svtx : mir(boost::vertices(cgraph))) {
        if (cgraph[svtx].code() != 's') {
            continue;
        }

        // Keep track of every blob in a slice from whence we came to a leaf.
        std::unordered_map<cluster_vertex_t, std::vector<cluster_vertex_t>> leaf2blob;

        // Loop over the neighbor blobs of current slice.
        for (auto bedge : mir(boost::out_edges(svtx, cgraph))) {
            auto bvtx = boost::target(bedge, cgraph);

            std::vector<cluster_vertex_t> leaves;
            connected_leaves(leaves, cgraph, bvtx, leaf_code);

            for (auto lvtx : leaves) {
                leaf2blob[lvtx].push_back(bvtx);
            }
        }

        // Process leaves.
        for (const auto& [lvtx, bvtxs] : leaf2blob) {
            size_t nblobs = bvtxs.size();

            if (nblobs < 2) {
                continue;
            }

            // pair-wise blobs
            for (size_t ind1=0; ind1 < nblobs-1; ++ind1) {
                auto bvtx1 = bvtxs[ind1];

                auto bs_vtx1 = c2bs[bvtx1];

                for (size_t ind2=ind1+1; ind2 < nblobs; ++ind2) {
                    auto bvtx2 = bvtxs[ind2];
                    auto bs_vtx2 = c2bs[bvtx2];

                    // figure out WirePlaneId of this lvtx
                    WirePlaneId wpid(0);
                    int index{-1};
                    const auto& obj = cgraph[lvtx];
                    const char lcode = obj.code();
                    if (lcode == 'w') { // wire
                        auto iptr = get<IWire::pointer>(obj.ptr);
                        assert(iptr);
                        wpid = iptr->planeid();
                        index = iptr->index();
                    }
                    else if (lcode == 'c') { // channel
                        auto iptr = get<IChannel::pointer>(obj.ptr);
                        assert(iptr);
                        wpid = iptr->planeid();
                        index = iptr->index();
                    }
                    else {
                        continue;
                    }

                    // if a same-layer edge already exists, refresh its beg-end
                    const auto key = layer_edge_key(bs_vtx1, bs_vtx2, wpid.layer());
                    auto [lit, inserted] = layer_edges.try_emplace(key, shadows.edges.size());
                    if (!inserted) {
                        Edge& eobj = shadows.edges[lit->second].prop;
                        eobj.beg = std::min(eobj.beg, index);
                        eobj.end = std::max(eobj.end, index + 1);
                        continue;
                    }

                    // first encounter: a new edge
                    shadows.edges.push_back({bs_vtx1, bs_vtx2, Edge{index, index + 1, wpid}});
                }
            }
        }
    }

    return shadows;
}

BlobShadow::graph_t BlobShadow::shadow(const cluster_graph_t& cgraph, char leaf_code)
{
    const auto shadows = shadow_list(cgraph, leaf_code);

    BlobShadow::graph_t bsgraph;
    bsgraph[boost::graph_bundle].stype = shadows.stype;
    for (const auto& node : shadows.nodes) {
        boost::add_vertex(node, bsgraph);
    }
    // Adding in list order preserves the first-encounter edges() order.
    for (const auto& er : shadows.edges) {
        auto [edge, added] = boost::add_edge(er.v1, er.v2, er.prop, bsgraph);
        if (!added) {
            // this should not happen
            THROW(RuntimeError() << errmsg{"edge not added!"});
        }
    }
    return bsgraph;
}

