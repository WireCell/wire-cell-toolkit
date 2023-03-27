#include "WireCellAux/BlobShadow.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellIface/ICluster.h"
#include "WireCellUtil/GraphTools.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

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

BlobShadow::graph_t BlobShadow::shadow(const cluster_graph_t& cgraph, char leaf_code)
{
    using dvertex_t = cluster::directed::vertex_t;

    BlobShadow::graph_t bsgraph; // will return

    // Convert to directed so BFS does not "crawl up" from leaves.
    auto dgraph = cluster::directed::type_directed(cgraph);

    // Loop over blobs, to load up output graph, old->new map.
    std::unordered_map<dvertex_t, BlobShadow::vdesc_t> c2bs;
    for (auto bvtx : mir(boost::vertices(dgraph))) {
        if (dgraph[bvtx].code() == 'b') {
            c2bs[bvtx] = boost::add_vertex({bvtx}, bsgraph);
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
            LeafVisitor<cluster::directed::graph_t> leafvis{{}, leaf_code, leaves};
            boost::breadth_first_search(dgraph, bvtx, boost::visitor(leafvis));

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

                    auto [edge, added] = boost::add_edge(bs_vtx1, bs_vtx2, bsgraph);

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
                    Edge& eobj = bsgraph[edge];
                    eobj.wpid = wpid;
                    if (added) {    // first time
                        eobj.beg = index;
                        eobj.end = index+1;
                    }
                    else {
                        eobj.beg = std::min(eobj.beg, index);
                        eobj.end = std::max(eobj.end, index+1);
                    }
                }
            }
        }
    }
    return bsgraph;
}

