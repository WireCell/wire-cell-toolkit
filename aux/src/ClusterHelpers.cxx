// Very general helpers.
// See also ClusterHelper*.cxx

#include "WireCellAux/ClusterHelpers.h"

#include "WireCellIface/ISlice.h"
#include "WireCellIface/IFrame.h"

#include "WireCellUtil/GraphTools.h"

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/copy.hpp>

using namespace WireCell;
using WireCell::GraphTools::mir;
using slice_t = cluster_node_t::slice_t;

// maybe useful to export
ISlice::vector WireCell::Aux::find_slices(const cluster_graph_t& gr)
{
    ISlice::vector ret;
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        const auto& vobj = gr[vtx];
        if (vobj.ptr.index() != 4) {
            continue;
        }

        auto islice = std::get<slice_t>(vobj.ptr);
        ret.push_back(islice);
    }
    return ret;
}

IFrame::pointer WireCell::Aux::find_frame(const cluster_graph_t& gr)
{
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        const auto& vobj = gr[vtx];
        if (vobj.ptr.index() != 4) {
            continue;
        }

        auto islice = std::get<slice_t>(vobj.ptr);
        return islice->frame();
    }
    return nullptr;
}

// std::string WireCell::Aux::name(const ICluster& cluster)
// {
//     std::stringstream ss;
//     ss << "cluster_" << cluster.ident();
//     return ss.str();
// }

Aux::blobs_by_slice_t Aux::blobs_by_slice(const cluster_graph_t& gr)
{
    Aux::blobs_by_slice_t ret;
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        const auto& vobj = gr[vtx];
        if (vobj.code() != 's') {
            continue;
        }
        auto islice = std::get<cluster_node_t::slice_t>(vobj.ptr);
        IBlob::vector blobs;
        for (auto nnvtx : mir(boost::adjacent_vertices(vtx, gr))) {
            const auto& nnobj = gr[nnvtx];
            if (nnobj.code() == 'b') {
                auto iblob = std::get<cluster_node_t::blob_t>(nnobj.ptr);
                blobs.push_back(iblob);
            }
        }
        ret[islice] = blobs;
    }
    return ret;
}

Aux::code_counts_t Aux::counts_by_type(const cluster_graph_t& gr)
{
    Aux::code_counts_t ret;
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        const auto& vobj = gr[vtx];
        const char code = vobj.code();
        ret[code] += 1;
    }
    return ret;
}
Aux::cluster_vertices_type_map_t Aux::vertices_type_map(const cluster_graph_t& cgraph)
{
    Aux::cluster_vertices_type_map_t ret;
    for (auto vtx : boost::make_iterator_range(boost::vertices(cgraph))) {
        const auto& vobj = cgraph[vtx];
        const char code = vobj.code();
        ret[code].push_back(vtx);
    }
    return ret;
}

using Filtered = boost::filtered_graph<cluster_graph_t, boost::keep_all,
                                       std::function<bool(cluster_vertex_t)> >;

static
bool inletter(char c, const std::string& letters)
{
    for (char l : letters) {
        if (l == c) return true;
    }
    return false;
}

/// Return an sbw graph given an sbwc cluster graph.  The c-nodes are trimmed.
WireCell::cluster_graph_t Aux::extract_sbw(const cluster_graph_t& cgraph)
{
    Filtered fg(cgraph, {}, [&](cluster_vertex_t vtx) {
        return inletter(cgraph[vtx].code(), "sbw");
    });
    
    cluster_graph_t cg;
    boost::copy_graph(fg, cg);
    return cg;    
}


/// Return an sbc graph given an sbwc cluster graph.  The w-nodes
/// are removed and b-c edges made.
WireCell::cluster_graph_t Aux::extract_sbc(const cluster_graph_t& cgraph)
{
    cluster_graph_t cg = cgraph;

    const auto nvtxs = boost::num_vertices(cg);

    std::unordered_map<int, cluster_vertex_t> chid2vtx;
    for (auto vtx : mir(boost::vertices(cg))) {
        const auto& vobj = cg[vtx];
        const char code = vobj.code();
        if (code == 'w') {
            cluster_vertex_t cvtx{nvtxs};
            std::vector<cluster_vertex_t> bvtxs;
            for (auto nnvtx : mir(boost::adjacent_vertices(vtx, cg))) {
                const auto& nnobj = cg[nnvtx];
                const char nncode = nnobj.code();
                switch (nncode) {
                    case 'c':
                        cvtx = nnvtx;
                        break;
                    case 'b':
                        bvtxs.push_back(nnvtx);
                        break;
                }
            }
            if (cvtx < nvtxs) { // really should never fail
                for (auto bvtx : bvtxs) {
                    boost::add_edge(bvtx, cvtx, cg);
                }
            }
        }
    }
    Filtered fg(cg, {}, [&](cluster_vertex_t vtx) {
        return inletter(cg[vtx].code(), "sbc");
    });

    cluster_graph_t cg2;
    boost::copy_graph(fg, cg2);
    return cg2;
}
WireCell::cluster_graph_t Aux::extract_sbX(const cluster_graph_t& cgraph, char code)
{
    cluster_graph_t ret;
    if (code == 'w') {
        ret = extract_sbw(cgraph);
    }
    else if (code == 'c') {
        ret = extract_sbc(cgraph);
    }
    return ret;
}


Aux::cluster::directed::graph_t
Aux::cluster::directed::type_directed(const cluster_graph_t& cgraph,
                                      const std::string& order)
{
    std::map<char, size_t> code_index;
    for (size_t ind=0; ind<order.size(); ++ind) {
        code_index[order[ind]] = ind+1; // 0 is illegal
    }

    Aux::cluster::directed::graph_t dgraph;

    // Output retains input's vertex ordering and properties.
    for (auto vtx : mir(boost::vertices(cgraph))) {
        boost::add_vertex(cgraph[vtx], dgraph);
    }
    for (const auto& edge : mir(boost::edges(cgraph))) {
        auto t = boost::source(edge, cgraph);
        auto h = boost::target(edge, cgraph);

        size_t tnum = code_index[cgraph[t].code()];
        size_t hnum = code_index[cgraph[h].code()];
        if (tnum > 0 and hnum > 0 and tnum > hnum) {
            std::swap(t,h);
        }
        boost::add_edge(t,h,dgraph);
    }
    return dgraph;
}
