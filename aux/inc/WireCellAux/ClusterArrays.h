#ifndef WIRECELLAUX_CLUSTERARRAYS
#define WIRECELLAUX_CLUSTERARRAYS

#include "ClusterHelpers.h"
#include "WireCellIface/ICluster.h"
#include <boost/multi_array.hpp>
#include <boost/range/counting_range.hpp>

#include <vector>
#include <unordered_map>

namespace WireCell::Aux {

    /**
       ClusterArrays produces array representations from ICluster data.

       See also the aux/docs/ClusterArrays.org document.
    */
    class ClusterArrays {
      public:

        ClusterArrays(const cluster_graph_t& graph);
        ClusterArrays();
        ~ClusterArrays();

        // Initialize.  Causes a clear().
        void init(const cluster_graph_t& graph);

        // Clear all stored data.  This will invalidate any arrays
        // previously accessed and held by const reference.
        void clear();

        // Return node type codes represented by the graph.
        using node_code_t = char;
        std::vector<node_code_t> node_codes() const;

        // Return edge type codes represented by the graph.
        using edge_code_t = int;
        std::vector<edge_code_t> edge_codes() const;

        // Return a properly formed edge type code from two node type
        // codes.  It assures order by code
        edge_code_t edge_code(node_code_t nc1, node_code_t nc2);
        std::string edge_code_str(edge_code_t ec);

        // Edge array holds edges of one type.  An type is simply the
        // ordered composition of node types of the endpoints.  Each
        // edge is represented as a pair of indices (2 columns) into
        // the pair of node type arrays.  We use int (instead of
        // size_t) to save storage space and because the use of
        // unique, positive "ident" values, also of type int,
        // effectively limits us to 2^31 nodes of any one type.
        using edge_array_t = boost::multi_array<int, 2>;

        /// Return the edge array for the given pair of node type
        /// codes.
        const edge_array_t& edge_array(edge_code_t ec) const;

        // Each node type produces a 2D array of doubles.
        using node_array_t = boost::multi_array<double, 2>;

        /// Return the node array of the given node type code.
        const node_array_t& node_array(node_code_t nc) const;

      private:
        using node_store_t = std::unordered_map<node_code_t, node_array_t>;
        using edge_store_t = std::unordered_map<edge_code_t, edge_array_t>;
        node_store_t m_na;
        edge_store_t m_ea;

        // We must keep a bidirectional mapping between graph and
        // array stores.
        struct store_address_t {
            node_code_t code;
            node_array_t::size_type index;
        };
        using v2s_t = std::unordered_map<cluster_vertex_t, store_address_t>;
        v2s_t m_v2s;

        using node_row_t = node_array_t::array_view<1>::type;
        node_row_t node_row(cluster_vertex_t vtx);
        
        // heavy lifters
        void init_blob(const cluster_graph_t& graph, cluster_vertex_t vtx);
        void init_wire(const cluster_graph_t& graph, cluster_vertex_t vtx);
        void init_signals(const cluster_graph_t& graph, cluster_vertex_t vtx);

    };
}

#endif
