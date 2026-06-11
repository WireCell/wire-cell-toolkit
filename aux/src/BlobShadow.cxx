#include "WireCellAux/BlobShadow.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellIface/ICluster.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellUtil/Graph.h"
#include <boost/graph/breadth_first_search.hpp>

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



// A BFS visitor that records vertices of the given code.
template<typename Gr>
struct LeafVisitor : public boost::default_bfs_visitor
{
    using graph_t = Gr;
    using vertex_t = typename boost::graph_traits<Gr>::vertex_descriptor;
    using edge_t = typename boost::graph_traits<Gr>::edge_descriptor;
    
    char code;
    std::vector<vertex_t>& leaves;

    void examine_vertex(vertex_t v, const graph_t& g) const {
        if (g[v].code() == code) {
            leaves.push_back(v);
        }
    }
};

namespace {
    void connected_leaves(std::vector<cluster::directed::vertex_t> & leaves,
    const cluster::directed::graph_t & dgraph, const cluster::directed::vertex_t & bvtx, char leaf_code) {
        // leaf_code can be 'w' or 'c'
        std::unordered_set<char> valid_codes = {'w', 'c'};
        if (valid_codes.find(leaf_code) == valid_codes.end()) {
            // TODO: make some noise?
            return;
        }
        for(auto bedge : mir(boost::out_edges(bvtx, dgraph))) {
            auto wvtx = boost::target(bedge, dgraph);
            if (dgraph[wvtx].code() != 'w') {
                continue;
            }
            if (leaf_code == 'w') {
                leaves.push_back(wvtx);
                continue;
            }
            // if not 'w', find 'c'
            for(auto wedge : mir(boost::out_edges(wvtx, dgraph))) {
                auto cvtx = boost::target(wedge, dgraph);
                if (dgraph[cvtx].code() != 'c') {
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
    using dvertex_t = cluster::directed::vertex_t;
    // using dedge_t = cluster::directed::edge_t;

    BlobShadow::Shadows shadows; // will return
    shadows.stype = leaf_code;
    layer_edge_map_t layer_edges;

    // Convert to directed so BFS does not "crawl up" from leaves.
    auto dgraph = cluster::directed::type_directed(cgraph);

    // Loop over blobs, to load up output nodes, old->new map.
    std::unordered_map<dvertex_t, BlobShadow::vdesc_t> c2bs;
    for (auto bvtx : mir(boost::vertices(dgraph))) {
        if (dgraph[bvtx].code() == 'b') {
            c2bs[bvtx] = shadows.nodes.size();
            shadows.nodes.push_back({bvtx});
        }
    }

    // Loop again, to pick up slices and make edges
    for (auto svtx : mir(boost::vertices(dgraph))) {
        if (dgraph[svtx].code() != 's') {
            continue;
        }

        // Keep track of every blob in a slice from whence we came to a leaf.
        std::unordered_map<cluster_vertex_t, std::vector<cluster_vertex_t>> leaf2blob;

        // Loop over the neighbor blobs of current slice.
        for (auto bedge : mir(boost::out_edges(svtx, dgraph))) {
            auto bvtx = boost::target(bedge, dgraph);

            std::vector<dvertex_t> leaves;
            // LeafVisitor<cluster::directed::graph_t> leafvis{{}, leaf_code, leaves, verbose};
            // boost::breadth_first_search(dgraph, bvtx, boost::visitor(leafvis));
            connected_leaves(leaves, dgraph, bvtx, leaf_code);

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
                    const auto& obj = dgraph[lvtx];
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

