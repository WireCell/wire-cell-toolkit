/** Provide some helper functions for working with ICluster.

 */

#ifndef WIRECELLAUX_CLUSTERJSONIFY
#define WIRECELLAUX_CLUSTERJSONIFY

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/ISlice.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/Graph.h"
#include <sstream>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include <string>

namespace WireCell::Aux {

    std::string dumps(const cluster_graph_t& cgraph);

    /// Return JSON representation of the cluster.
    Json::Value jsonify(const cluster_graph_t& cgraph);

    /// A loader of cluster JSON files.  See also ClusterArrays.
    ///
    class ClusterLoader {
    public:

        // A context consisting of anode planes is needed to resolve
        // some JSON to references to objects.  / ValueError is thrown
        // if faces are not uniquely identified.
        using anodes_t = std::vector<IAnodePlane::pointer>;
        ClusterLoader(const anodes_t& anodes);
        
        // Return cluster graph corresponding to data in jgraph.  If
        // any frames are loaded, provide them.  Otherwise a
        // referenced frame ID is resolved into an IFrame is
        // constructed which only holds ident.  ValueError is thrown
        // jgraph is bad.
        cluster_graph_t load(const Json::Value& jgraph,
                             const IFrame::vector& known_frames = {}) const;

        // Return known anode face by id or nullptr if id not found;
        IAnodeFace::pointer face(WirePlaneId faceid) const;

    private:
        using faces_t = std::unordered_map<int, IAnodeFace::pointer>;
        faces_t m_faces;
    };

    /// Return GraphViz dot representation of a cluster graph like
    /// graph.
    template <typename CGraph>
    std::string dotify(const CGraph& cgraph) {
        std::stringstream ss;
        using vertex_t = typename boost::graph_traits<CGraph>::vertex_descriptor;
        boost::write_graphviz(ss, cgraph, [&](std::ostream& out, vertex_t v) {
            const auto& dat = cgraph[v];
            out << "[label=\"" << dat.code() << dat.ident() << "\"]";
        });
        return ss.str() + "\n";
    }

    // Return a cluster graph like graph with any nodes that have type
    // codes missing gfrom the "codes" string excluded.
    template <typename CGraph>
    CGraph type_filter(const CGraph& cgraph, const std::string& codes) {
        using vertex_t = typename boost::graph_traits<CGraph>::vertex_descriptor;
        using Filtered = typename boost::filtered_graph<CGraph, boost::keep_all,
                                               std::function<bool(vertex_t)> >;
        std::unordered_set<char> keep(codes.begin(), codes.end());
        Filtered fg(cgraph, {}, [&](vertex_t vtx) {
            return keep.find(cgraph[vtx].code()) != keep.end();
        });
        CGraph gr;
        boost::copy_graph(fg, gr);
        return gr;
    }


    /// Return the slices in the cluster.
    ISlice::vector find_slices(const cluster_graph_t& cgraph);

    /// Return the frame of the first slice in the cluster.  Note, in
    /// principle, clusters can span frames.
    IFrame::pointer find_frame(const cluster_graph_t& cgraph);

    /// Return the blobs for each slice
    using blobs_by_slice_t = std::unordered_map<ISlice::pointer, IBlob::vector>;
    blobs_by_slice_t blobs_by_slice(const cluster_graph_t& cgraph);

    /// Return number of vertices of each type by type code
    using code_counts_t = std::unordered_map<char, size_t>;
    code_counts_t counts_by_type(const cluster_graph_t& cgraph);

    // Return the vertex descriptors for the given type code (in "cwsbm")
    // std::vector<cluster_vertex_t> vdesc_by_code(const cluster_graph_t& cgraph, char code);

    // Return a flat vector of vertices in cgraph by their type code.
    // Vectors have no particular ordering.
    using cluster_vertices_t = std::vector<cluster_vertex_t>;
    using cluster_vertices_type_map_t = std::unordered_map<char, cluster_vertices_t>;
    cluster_vertices_type_map_t vertices_type_map(const cluster_graph_t& cgraph);

    namespace cluster::directed {
        using graph_t  = boost::adjacency_list<boost::setS, boost::vecS, boost::directedS, cluster_node_t>;
        using vertex_t = boost::graph_traits<graph_t>::vertex_descriptor;
        using edge_t = boost::graph_traits<graph_t>::edge_descriptor;

        /// Return a directed version of the cgraph with direction
        /// determined by a vertex type code order.  The vertices and
        /// edges and their order within their sets are unchanged in
        /// the returned graph.  However, the order of an individual
        /// edge edpoints (which is "source" vs "target") is subject
        /// change based on the "order" string which provides cluster
        /// node type codes in ascending value.  Eg, with the default
        /// "order" string, an edge formed with an s-type and a b-type
        /// will produce an s-b edge.  Edges between like types (eg
        /// b-b edge) or edges with the type code of one or both
        /// endpoints that are not represented in the "order" string
        /// will be output in their original order.
        graph_t type_directed(const cluster_graph_t& cgraph,
                              const std::string& order = "sbmwc");
    
    }
    /// Return an s-b-w graph given an s-b-w-c cluster graph.  The
    /// c-nodes (and any m-nodes) are removed.
    cluster_graph_t extract_sbw(const cluster_graph_t& cgraph);

    /// Return an sbc graph given an sbwc cluster graph.  The w-nodes
    /// are removed and b-c edges made.  Any m-nodes are also removed.
    cluster_graph_t extract_sbc(const cluster_graph_t& cgraph);

    /// Simply return an sbw or sbc depending on if code is 'w' or 'c'
    /// or return an empty graph.
    cluster_graph_t extract_sbX(const cluster_graph_t& cgraph, char code);

    /// Return counts indexed by node and/or edge code
    std::map<std::string, size_t> count(const cluster_graph_t& cgraph, bool nodes=true, bool edges=true);

    /// function to merge  boost graph
    template <typename GraphType>
    GraphType merge_graphs(const std::vector<std::reference_wrapper<const GraphType>>& graphs)
    {
        GraphType merged_graph;
        // merge the graphs
        for (const auto& graphr : graphs) {
            const auto& graph = graphr.get();
            std::unordered_map<size_t, size_t> vertex_map;
            // add the vertices
            for (const auto& vtx : GraphTools::mir(boost::vertices(graph))) {
                const auto& node = graph[vtx];
                auto new_vtx = boost::add_vertex(node, merged_graph);
                vertex_map[vtx] = new_vtx;
            }
            // add the edges
            for (const auto& edg : GraphTools::mir(boost::edges(graph))) {
                const auto& src = boost::source(edg, graph);
                const auto& tgt = boost::target(edg, graph);
                const auto& edge = graph[edg];
                auto new_src = vertex_map[src];
                auto new_tgt = vertex_map[tgt];
                boost::add_edge(new_src, new_tgt, edge, merged_graph);
            }
        }
        return merged_graph;
    }
}

#endif

