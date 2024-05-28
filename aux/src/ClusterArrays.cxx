/**
   Commentary, because this code is a bit messy.

   to_arrays converts a cluster graph to cluster array schema.  It
   delegates to all the init_<type>() functions.  init_slice() also
   handles channels.  The cluster array schema columns are partly
   encded in *_col constants.

   to_cluster is the inverse.  It requires a set of IAnodePlanes.  The
   function body is largely cribbed from ClusterHelpersLoader with the
   exception that here we exactly restore the vertex descriptors
   (something that needs fixing in ClusterHelpersLoader).

 */

#include "WireCellAux/ClusterArrays.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/GraphTools.h"

// for cluster reconstitution
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleWire.h"
#include "WireCellAux/SimpleChannel.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleMeasure.h"

#include "WireCellUtil/Graph.h"



using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Aux::ClusterArrays;
using namespace WireCell::GraphTools;

using channel_t = cluster_node_t::channel_t;
using wire_t = cluster_node_t::wire_t;
using blob_t = cluster_node_t::blob_t;
using slice_t = cluster_node_t::slice_t;
using meas_t = cluster_node_t::meas_t;

// Some hard-wired schema values
// See ClusterArray.org docs and keep in sync.
// For the common node array columns.
// Fixme: remove these when move to Dataset store
static const size_t desc_col=0; // all node and edge types
static const size_t ident_col=1; // all node types
static const size_t sigv_col=2;  // all but wires
static const size_t sigu_col=3;  // all but wires

// Edges do not have an integer descriptor, per se but we put an edge
// counter in the desc_col.
static const size_t tail_col=1;  // edge tail
static const size_t head_col=2;  // edge head


template<typename ArrayType>
void rezero(ArrayType& arr, const std::vector<size_t>& sv)
{
    // resize and zero
    arr.resize(sv);
    std::fill_n(arr.data(), arr.num_elements(), 0);
}


// CAVEAT: BIG FAT BODGE HACK
//
// APPEND new c-nodes to represent activity nodes.  The old c-nodes
// will have their IChannel::pointer NULLED.  While the old c-nodes
// represent physical channels the new c-nodes represent activity in a
// channel AND in a slice.  These new c-nodes are used to fill the
// activity array and together hold the information in the "activity
// map" of the ISlice instances.  Old c-nodes with nullptr shoudl be
// ignored but are left in order not to invalidate vertex descriptors.
//
cluster_graph_t ClusterArrays::bodge_channel_slice(cluster_graph_t graph)
{
    std::vector<cluster_vertex_t> old_cvtx;
    for (const auto& vtx : vertex_range(graph)) {
        if ('c' == graph[vtx].code()) {
            old_cvtx.push_back(vtx);
        }
    }

    // Get known wire nodes
    std::unordered_map<IWire::pointer, cluster_vertex_t> wire_nodes;
    for (auto wvtx : vertex_range(graph)) {
        auto wnode = graph[wvtx];
        if (node_code(wnode) == 'w') {
            wire_nodes[std::get<wire_t>(wnode.ptr)] = wvtx;
        }
    }
    
    // Walk each slice.
    for (const auto& svtx : vertex_range(graph)) {
        const auto& snode = graph[svtx];
        if (node_code(snode) != 's') { continue; }

        // Add a new cnode for each activity.  Make c-w and c-s edges.
        // Remember new cnodes to later make c-m edges.
        std::unordered_map<IChannel::pointer, cluster_vertex_t> active_cnodes;
        const auto& islice = std::get<slice_t>(snode.ptr);
        auto activity = islice->activity(); // modifiable copy
        for (const auto& [ich, val] : activity) {
            auto cvtx = boost::add_vertex(cluster_node_t(ich), graph);
            active_cnodes[ich] = cvtx;
            boost::add_edge(cvtx, svtx, graph);

            for (const auto& iwire : ich->wires()) {
                auto wit = wire_nodes.find(iwire);
                cluster_vertex_t wvtx;
                if (wit == wire_nodes.end()) {
                    wvtx = boost::add_vertex(cluster_node_t(iwire), graph);
                    wire_nodes[iwire] = wvtx;
                }
                else {
                    wvtx = wit->second;
                }
                boost::add_edge(cvtx, wvtx, graph);
            }
        }

        // walk s-b-m to make new c-m edge
        for (auto bvtx : mir(boost::adjacent_vertices(svtx, graph))) {
            const auto& bnode = graph[bvtx];
            if (node_code(bnode) != 'b') { continue; }
            for (auto mvtx : mir(boost::adjacent_vertices(bvtx, graph))) {
                const auto& mnode = graph[mvtx];
                if (node_code(mnode) != 'm') { continue; }
                for (auto cvtx : mir(boost::adjacent_vertices(mvtx, graph))) {
                    const auto& cnode = graph[cvtx]; // old
                    if (node_code(cnode) != 'a') { continue; }
                    auto active_cvtx = active_cnodes[std::get<channel_t>(cnode.ptr)];
                    boost::add_edge(active_cvtx, mvtx, graph);
                }
            }
        }
    }

    // Invalidate the original ("old") c-nodes.  These the code()
    // method for these nodes will now return 0 instead of 'c'.  Use
    // is_old_chan() to test.
    for (auto cvtx : old_cvtx) {
        graph[cvtx].ptr = (size_t)0;
    }

    return graph;
}

static bool is_old_chan(const cluster_node_t& node)
{
    // a nulled ptr makes code() return 0.
    if (! node.code()) {
        return true;        
    }
    return false;
}

///// used only in commented out debugging.  Comment it out to avoid -Werror=unused-function
// static
// void dump(const cluster_graph_t& graph)
// {
//     size_t nold_chan=0;
//     std::unordered_map<char, std::unordered_set<int>> code2idents;
//     for (const auto& vdesc : vertex_range(graph)) {
//         const auto& node = graph[vdesc];
//         if (is_old_chan(node)) {
//             ++nold_chan;
//             continue;
//         }
//         char code = node.code();
//         int ident = node.ident();
//         code2idents[code].insert(ident);
//     }
//     const std::map<std::string, size_t>& counts = count(graph);
//     std::cerr << "nodes: ";
//     for (const auto& [code,num] : counts) {
//         if (code.size() != 1) { continue; }
//         if (code[0] == 'c') {
//             std::cerr << code << ":" << num-nold_chan << "(" << nold_chan << ") " ;
//             continue;
//         }
//         std::cerr << code << ":" << num << " " ;
//     }
//     std::cerr << "edges: ";
//     for (const auto& [code,num] : counts) {
//         if (code.size() == 2) {
//             std::cerr << code << ":" << num << " ";
//         }
//     }
//     std::cerr << "\n";
// }

using node_row_t = node_array_t::array_view<1>::type;

// Local context for converting from graph to arrays
struct ToArrays {
    node_array_set_t& nas;
    edge_array_set_t& eas;

    ToArrays(node_array_set_t& nas, edge_array_set_t& eas) : nas(nas), eas(eas) {}

    struct store_address_t {
        node_code_t code;
        node_array_t::size_type index;
    };
    std::unordered_map<cluster_vertex_t, store_address_t> v2s;
    void store_address(cluster_vertex_t vtx, node_code_t code,
                       node_array_t::size_type index)
    {
        v2s[vtx] = store_address_t{code,index};
    }

    node_row_t node_row(cluster_vertex_t vtx) {
        const auto sa = v2s.at(vtx);
        typedef boost::multi_array_types::index_range range;
        return nas.at(sa.code)[boost::indices[sa.index][range()]];
    }
};

using chid2desc_t = std::unordered_map<int, cluster_vertex_t>;

// must process s+c together to spread activity onto channel
static void init_slice(const cluster_graph_t& graph,
                       ToArrays& ta,
                       cluster_vertex_t svtx,
                       const chid2desc_t& chid2desc)
{
    node_row_t srow = ta.node_row(svtx);
    const auto& snode = graph[svtx];

    const auto& islice = std::get<slice_t>(snode.ptr);
    const int sident = islice->ident();
    const auto& iframe = islice->frame();
    const int fident = iframe ? iframe->ident() : -1;
    const auto& activity = islice->activity();
        
    // The slices get a sum of all their activity.

    srow[desc_col] = svtx;
    srow[ident_col] = sident;

    for (const auto& [_, sig] : activity) {
        srow[sigv_col] += sig.value();
        srow[sigu_col] += sig.uncertainty();
    }

    {   // slice start and span, just after sigv/sigu.
        node_array_t::size_type col = sigu_col;
        srow[++col] = fident;
        srow[++col] = islice->start();
        srow[++col] = islice->span();
    }
    
    // This is bodged activity
    for (auto cvtx : mir(boost::adjacent_vertices(svtx, graph))) {
        const auto& cnode = graph[cvtx];
        if (node_code(cnode) != 'a') { continue; }
        if (is_old_chan(cnode)) { continue; }

        auto ichan = std::get<channel_t>(cnode.ptr);
        const int chid = ichan->ident();
        auto sigit = activity.find(ichan);
        if (sigit == activity.end()) {
            continue;
        }
        auto cval = sigit->second;
        auto crow = ta.node_row(cvtx);
        crow[desc_col] = chid2desc.at(chid);
        crow[ident_col] = chid;
        crow[sigv_col] = cval.value();
        crow[sigu_col] = cval.uncertainty();
        crow[sigu_col+1] = ichan->index();
        crow[sigu_col+2] = ichan->planeid().ident();
    }
}

static void init_measure(const cluster_graph_t& graph,
                         node_row_t mrow, cluster_vertex_t mvtx)
{
    const auto& mnode = graph[mvtx];

    auto imeas = std::get<meas_t>(mnode.ptr);
    auto mval = imeas->signal();
    int ident = imeas->ident();
    mrow[desc_col] = mvtx;
    mrow[ident_col] = ident;
    mrow[sigv_col] = mval.value();
    mrow[sigu_col] = mval.uncertainty();
    mrow[sigu_col+1] = imeas->planeid().ident();
}

static void init_wire(const cluster_graph_t& graph,
                      node_row_t wrow, cluster_vertex_t wvtx)
{
    const auto& wnode = graph[wvtx];

    const auto& iwire = get<wire_t>(wnode.ptr);
    const auto [p1,p2] = iwire->ray();

    int ident = iwire->ident();
    wrow[desc_col] = wvtx;
    wrow[ident_col] = ident;

    node_array_t::size_type col = ident_col;

    wrow[++col] = iwire->index();
    wrow[++col] = iwire->segment();
    wrow[++col] = iwire->channel();
    wrow[++col] = iwire->planeid().ident();

    wrow[++col] = p1.x();
    wrow[++col] = p1.y();
    wrow[++col] = p1.z();
    wrow[++col] = p2.x();
    wrow[++col] = p2.y();
    wrow[++col] = p2.z();
}

static void init_blob(const cluster_graph_t& graph,
                      node_row_t brow, cluster_vertex_t bvtx)
{
    const auto& bnode = graph[bvtx];

    const auto& iblob = get<blob_t>(bnode.ptr);
    const auto iface = iblob->face();
    const auto islice = iblob->slice();
    const auto& coords = iface->raygrid();

    const auto& shape = iblob->shape();
    const auto& strips = shape.strips();

    int ident = iblob->ident();
    brow[desc_col] = bvtx;
    brow[ident_col] = ident;
    brow[sigv_col] = iblob->value();
    brow[sigu_col] = iblob->uncertainty();

    node_array_t::size_type col = sigu_col + 1;

    // 3:faceid, 4:sident, 5:start, 6:span
    brow[col++] = WirePlaneId(kUnknownLayer, iface->which(), iface->anode()).ident();
    brow[col++] = islice->ident();
    brow[col++] = islice->start();
    brow[col++] = islice->span();

    // 3x2 wire bounds: [7:13]
    for (size_t sind = 0; sind < 3; ++sind) {
        size_t strip_index = 2 + sind; // skip "virtual" layers
        const auto& bounds = strips[strip_index].bounds;

        brow[col++] = bounds.first;  // these are 
        brow[col++] = bounds.second; // WIP numbers
    }

    // 13, [14:38]
    //// corners is an awkward array.  Each of the 3 pairs of planes
    //// has 4 corners so a total of 12 but only those corners
    //// contained by each of the three "strips" that make up the blob
    //// are relevant and given.  Best we can do, I think is to use
    //// 1+2*12 = 25 columns.  First to given the number of corners,
    //// than that number of (y,z) pairs and zero pad any remainder.
    //// Alternatively we may define a second array type of shape 3 x
    //// N with columns (rowid, y, z).

    const auto& corners = shape.corners();
    brow[col++] = corners.size();
    for (const auto& c : corners) {
        auto pos = coords.ray_crossing(c.first, c.second);
        brow[col++] = pos[1];
        brow[col++] = pos[2];
    }
}

void ClusterArrays::to_arrays(const cluster_graph_t& cin_graph,
                              node_array_set_t& nas, edge_array_set_t& eas)
{
    ToArrays ta(nas, eas);
    
    // Remember original channel vertex descriptor
    chid2desc_t chid2desc;
    for (const auto& vtx : vertex_range(cin_graph)) {
        const auto& node = cin_graph[vtx];
        auto nc = node.code();
        if (nc == 'c') {
            chid2desc[node.ident()] = vtx;
        }
    }

    // Change c-nodes from representing physical channels to
    // representing per-slice channels associated with activites.
    const cluster_graph_t graph = bodge_channel_slice(cin_graph);

    // {                           // debugging
    //     std::cerr << "pre-bodge:  ";
    //     dump(cin_graph);
    //     std::cerr << "post-bodge: ";
    //     dump(graph);
    // }

    // Get number of rows in per node type arrays.
    std::unordered_map<node_code_t, size_t> counts;
    for (const auto& vdesc : vertex_range(graph)) {
        const auto& node = graph[vdesc];
        if (is_old_chan(node)) {
            continue;
        }
        node_code_t code = node_code(node);
        ++counts[code];
    }

    // Number of columns in per node type arrays.  This hardwires
    // cluster array schema so one must carefully count columns listed
    // in the ClusterArray.org document.
    std::unordered_map<node_code_t, size_t> node_array_ncols = {
        // [desc, ident, value, uncertainty, index, wpid]
        {'a',1+3+2},            // 6
        // [desc, ident, wip, seg, ch, plane, [tail x,y,z], [head x,y,z]]
        {'w',1+11},             // 12
        // [desc,ident,sigv,sigu, faceid,sliceid,start,span, [3x2 wire bounds], nc,[12 corner y,z]]
        {'b',1+3+4+3*2+1+12*2}, // 39
        // [desc,ident,sigv,sigu, frameid, start, span]
        {'s',1+3+3},            // 7
        // [desc,ident,sigv,sigu, wpid]
        {'m',1+4},              // 5
    };

    // Allocate the node arrays.
    for (const auto& [code, count] : counts) {
        const std::vector<size_t> shape = {count, node_array_ncols[code]};
        rezero(nas[code], shape);
    }

    // Loop vertices, build vtx<->row lookups and set ident column.
    // Skip any nulled c-nodes.
    counts.clear();             // reuse to track nrows.
    for (const auto& vdesc : vertex_range(graph)) {
        const auto& node = graph[vdesc];
        if (is_old_chan(node)) {
            continue;
        }
        const auto code = node_code(node);

        node_array_t::size_type row = counts[code];
        ++counts[code];

        ta.store_address(vdesc, code, row);

        nas.at(code)[row][ident_col] = node.ident();
    }
    
    // One more time to process things that require graph traversal.
    for (const auto& vdesc : vertex_range(graph)) {
        const auto& node = graph[vdesc];
        if (is_old_chan(node)) {
            continue;
        }

        const auto code = node_code(node);

        if (code == 's') { // slices and channels
            init_slice(graph, ta, vdesc, chid2desc);
            continue;
        }
        if (code == 'm') {
            init_measure(graph, ta.node_row(vdesc), vdesc);
            continue;
        }
        if (code == 'w') {
            init_wire(graph, ta.node_row(vdesc), vdesc);
            continue;
        }
        if (code == 'b') {
            init_blob(graph, ta.node_row(vdesc), vdesc);
            continue;
        }
    }

    // And the edges.  First to collect the data and learn the sizes.
    std::unordered_map<edge_code_t, size_t> ecounts;
    for (const auto& edge : mir(boost::edges(graph))) {
        auto tvtx = boost::source(edge, graph);
        auto hvtx = boost::target(edge, graph);

        const auto& tnode = graph[tvtx];
        if (is_old_chan(tnode)) continue;
        const auto& hnode = graph[hvtx];
        if (is_old_chan(hnode)) continue;

        auto tsa = ta.v2s.at(tvtx);
        auto hsa = ta.v2s.at(hvtx);

        // Does not need this as edge_code() does it internally.
        if (hsa.code < tsa.code) {
            std::swap(tsa, hsa);
            std::swap(tvtx, hvtx);
        }

        auto ecode = edge_code(tsa.code, hsa.code);
        ++ecounts[ecode];
    }
    for (const auto& [ecode, esize] : ecounts) {
        // rezero(edge_array(ecode), {esize, 3});
        rezero(eas[ecode], {esize, 3});
    }

    // And again, this time with feeling!
    ecounts.clear();
    int edesc = 0;
    for (const auto& edge : mir(boost::edges(graph))) {
        auto tvtx = boost::source(edge, graph);
        auto hvtx = boost::target(edge, graph);

        const auto& tnode = graph[tvtx];
        if (is_old_chan(tnode)) continue;
        const auto& hnode = graph[hvtx];
        if (is_old_chan(hnode)) continue;

        auto tsa = ta.v2s.at(tvtx);
        auto hsa = ta.v2s.at(hvtx);

        // Do manual ordering to keep code/index correlated
        if (hsa.code < tsa.code) {
            std::swap(tsa, hsa);
            std::swap(tvtx, hvtx);
        }

        auto ecode = edge_code(tsa.code, hsa.code);
        edge_array_t::size_type erow = ecounts[ecode];
        ++ecounts[ecode];
        eas.at(ecode)[erow][desc_col] = edesc++;
        eas.at(ecode)[erow][tail_col] = tsa.index;
        eas.at(ecode)[erow][head_col] = hsa.index;
    }
}


cluster_graph_t ClusterArrays::to_cluster(const node_array_set_t& nas,
                                          const edge_array_set_t& eas,
                                          const anodes_t& anodes)
{
    std::unordered_map<int, IAnodeFace::pointer> known_faces;
    for (const auto& anode : anodes) {
        for (const auto& face : anode->faces()) {
            const WirePlaneId afid(kUnknownLayer, face->which(), face->anode());
            known_faces[afid.ident()] = face;
        }
    }
    auto find_face = [&](WirePlaneId wpid) -> IAnodeFace::pointer {
        WirePlaneId afid(kUnknownLayer, wpid.face(), wpid.apa());

        auto it = known_faces.find(afid.ident());
        if (it == known_faces.end()) {
            THROW(ValueError() << errmsg("missing face.  Provide cluster_graph anodes consistent with data"));
        }
        return it->second;
    };


    cluster_graph_t graph;
    typedef boost::multi_array_types::index_range range;

    // Collect unique vertices (channels will be dup'ed).
    {
        std::unordered_set<cluster_vertex_t> all_the_vertices;
        for (const auto& [nc, arr] : nas) {
            const size_t nrows = arr.shape()[0];
            for (size_t ind=0; ind<nrows; ++ind) {
                cluster_vertex_t vtx = arr[ind][desc_col];
                all_the_vertices.insert(vtx);
            }
        }
        const size_t nnodes = all_the_vertices.size();
        for (size_t ind=0; ind<nnodes; ++ind) {
            boost::add_vertex(graph);
        }
    }
    // Hook up all the edges and fill activities when we get "asedges".  
    for (const auto& [ec, earr] : eas) {
        if (ec == edge_code('a','s')) {
            // a-s edges are cluster array schema and not cluster
            // graph schema.  They will be used to fill slice activity
            // maps below.  But we need c-nodes made whole first.
            continue;
        }

        const std::string ecs = to_string(ec);
        const auto& tarr = nas.at(ecs[0]);
        const auto& harr = nas.at(ecs[1]);
        const size_t nedges = earr.shape()[0];

        for (size_t ind=0; ind<nedges; ++ind) {
            // ignore old edge ordering (desc_col).
            const auto& tind = earr[ind][tail_col];
            const auto& hind = earr[ind][head_col];
            cluster_vertex_t tvtx = tarr[tind][desc_col];
            cluster_vertex_t hvtx = harr[hind][desc_col];

            // When dealing with edges from a-nodes, we simply re-add
            // them for all activities in the same channel.  This
            // redundance is absorbed by cluster_graph_t using a setS
            // for edges.  In the end, we get only one edge.
            boost::add_edge(tvtx, hvtx, graph);
        }
    }

    // Graph structure is made, now need to populate the node data.

    {                           // wires
        const auto& arr = nas.at('w');
        const size_t nnodes = arr.shape()[0];
        for (size_t ind=0; ind<nnodes; ++ind) {
            const auto& wrow = arr[boost::indices[ind][range()]];
        
            const cluster_vertex_t wvtx = wrow[desc_col];
            node_array_t::size_type col = ident_col;
            const int ident = wrow[col++];
            const int index = wrow[col++];
            const int seg = wrow[col++];
            const int chid = wrow[col++];
            const WirePlaneId wpid(wrow[col++]);
            const double x1 = wrow[col++];
            const double y1 = wrow[col++];
            const double z1 = wrow[col++];
            const double x2 = wrow[col++];
            const double y2 = wrow[col++];
            const double z2 = wrow[col++];
            Ray ray(Point(x1,y1,z1), Point(x2,y2,z2));

            IWire::pointer iwire = std::make_shared<SimpleWire>(wpid, ident, index, chid, ray, seg);
            graph[wvtx].ptr = iwire;
        }
    }


    std::unordered_map<int, IChannel::pointer> chid2ich;
    {                           // activity->channels
        const auto& arr = nas.at('a');
        const size_t nnodes = arr.shape()[0];
        for (size_t ind=0; ind<nnodes; ++ind) {
            const auto& crow = arr[boost::indices[ind][range()]];
            const cluster_vertex_t cvtx = crow[desc_col];

            const int chid = crow[ident_col];
            // const float val = crow[sigv_col];
            // const float unc = crow[sigu_col];
            // cvtx2sig[cvtx] = {val,unc};

            // We have same cvtx for many activity in a given channel.
            // Restore the c-node only for the first time seen.
            if (graph[cvtx].code() == 'c') { // code is 0 if nullptr
                continue;
            }

            const int index = crow[sigu_col+1];
            // const int wpid = crow[sigu_col+2];

            // Check neighbors to get the just added wires.
            IWire::vector wires;
            for (auto nnvtx : mir(boost::adjacent_vertices(cvtx, graph))) {
                const auto& nnobj = graph[nnvtx];
                if (nnobj.code() == 'w') {
                    auto iwire = get<wire_t>(nnobj.ptr);
                    wires.push_back(iwire);
                }
            }

            IChannel::pointer ich = std::make_shared<SimpleChannel>(chid, index, wires);
            chid2ich[chid] = ich;
            graph[cvtx].ptr = ich;
        }
    }

    std::vector<ISlice::map_t> slice_activities;
    {                  // Collect slice activities by s-node row index
        const auto& edges = eas.at(edge_code('a','s'));
        const size_t nedges = edges.shape()[0];
        const auto& ans = nas.at('a');
        const size_t narows = ans.shape()[0];
        // const auto& sns = nas.at('s');
        slice_activities.resize(ans.size());
        for (size_t ind=0; ind<nedges; ++ind) {
            const size_t aind = edges[ind][tail_col];
            const size_t sind = edges[ind][head_col];
            if (aind < 0 or aind >= narows) {
                std::cerr << "edge to activity index " << aind << " out of total " << narows << "\n";
                THROW(ValueError() << errmsg{"cluster array data is corrupt"});
            }
            const int  chid = ans[aind][ident_col];
            const float val = ans[aind][sigv_col];
            const float unc = ans[aind][sigu_col];
            IChannel::pointer ich = chid2ich.at(chid);
            slice_activities[sind][ich] = {val,unc};
        }
    }

    std::unordered_map<int, std::shared_ptr<SimpleFrame>> made_frames;
    {                           // slices
        const auto& arr = nas.at('s');
        const size_t nnodes = arr.shape()[0];
        for (size_t ind=0; ind<nnodes; ++ind) {
            const auto& srow = arr[boost::indices[ind][range()]];

            const cluster_vertex_t svtx = srow[desc_col];
            const int ident = srow[ident_col];
            node_array_t::size_type col = sigu_col;
            const int frame_ident = srow[++col];
            const double start = srow[++col];
            const double span = srow[++col];
            
            std::shared_ptr<SimpleFrame> sframe;
            auto fit = made_frames.find(frame_ident);
            if (fit == made_frames.end()) {
                sframe = std::make_shared<SimpleFrame>(frame_ident, start /*, tick?*/);
                made_frames[frame_ident] = sframe;
            }
            else {
                sframe = fit->second;
                if (sframe->time() > start) {
                    sframe->set_time(start); // the best we can do
                }
            }

            const auto& activities = slice_activities.at(ind);

            ISlice::pointer islice = std::make_shared<SimpleSlice>(
                sframe, ident, start, span, activities);
            graph[svtx].ptr = islice;
        }
    }

    // An ICluster output by BlobClustering has no 'm' nodes.
    if (nas.find('m') != nas.end()) { // measure
        const auto& arr = nas.at('m');
        const size_t nnodes = arr.shape()[0];
        for (size_t ind=0; ind<nnodes; ++ind) {
            const auto& mrow = arr[boost::indices[ind][range()]];

            const cluster_vertex_t mvtx = mrow[desc_col];
            const int ident = mrow[ident_col];
            const float val = mrow[sigv_col];
            const float unc = mrow[sigu_col];
            const WirePlaneId wpid(mrow[sigu_col+1]);
            IMeasure::value_t sig{val,unc};

            IChannel::vector chans;
            for (auto nnvtx : mir(boost::adjacent_vertices(mvtx, graph))) {
                const auto& nnobj = graph[nnvtx];
                if (nnobj.code() == 'c') {
                    chans.push_back( get<channel_t>(nnobj.ptr) );
                }
            }

            IMeasure::pointer imeas = std::make_shared<SimpleMeasure>(ident, wpid, sig, chans);
            graph[mvtx].ptr = imeas;
        }        
    }

    {                           // blobs
        const auto& arr = nas.at('b');
        const size_t nnodes = arr.shape()[0];
        for (size_t ind=0; ind<nnodes; ++ind) {
            const auto& brow = arr[boost::indices[ind][range()]];

            const cluster_vertex_t bvtx = brow[desc_col];
            const int ident = brow[ident_col];
            const double val = brow[sigv_col];
            const double unc = brow[sigu_col];

            node_array_t::size_type col = sigu_col;
            const WirePlaneId face_wpid(brow[++col]);

            // find slice and fill view strips
            std::map<int, RayGrid::Strip> strips;
            strips[0] = RayGrid::Strip{0, {0,1}};
            strips[1] = RayGrid::Strip{1, {0,1}};
            ISlice::pointer islice = nullptr;
            for (auto nnvtx : mir(boost::adjacent_vertices(bvtx, graph))) {
                const auto& nnobj = graph[nnvtx];
                if (nnobj.code() == 's') {
                    islice = get<slice_t>(nnobj.ptr);
                    continue;
                }
                if (nnobj.code() != 'w') {
                    continue;
                }
                auto iwire = get<wire_t>(nnobj.ptr);
                auto wpid = iwire->planeid();
                const int wip = iwire->index();
                const int layer = 2 + wpid.index();
                auto sit = strips.find(layer);
                if (sit == strips.end()) {
                    strips[layer] = RayGrid::Strip{layer, {wip,wip+1}};
                    continue;
                }
                auto& strip = sit->second;
                strip.bounds.first = std::min(strip.bounds.first, wip);
                strip.bounds.second = std::max(strip.bounds.second, wip+1);
            }

            auto iface = find_face(face_wpid);
            if (!iface) {
                raise<ValueError>("failed to find face for wpid %d", face_wpid.ident());
            }

            const RayGrid::Coordinates& coords = iface->raygrid();
            RayGrid::Blob bshape;
            for (const auto& [layer, strip] : strips) {
                bshape.add(coords, strip);
            }

            IBlob::pointer iblob = std::make_shared<SimpleBlob>(ident, val, unc, bshape, islice, iface);
            graph[bvtx].ptr = iblob;
        }
    }

    return graph;
}

