#ifndef WIRECELLAUX_CLUSTERARRAYS
#define WIRECELLAUX_CLUSTERARRAYS

#include "ClusterHelpers.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/ITensorSet.h"
#include "WireCellUtil/MultiArray.h"
#include <boost/range/counting_range.hpp>

#include <vector>
#include <unordered_map>

namespace WireCell::Aux::ClusterArrays {

    // See aux/docs/ClusterArrays.org document for details.


    // Node object attributes are all coerced into double precision
    // floats.
    using node_element_t = double;
    using node_array_t = boost::multi_array<node_element_t, 2>;

    // Edges are represented as indices into arrays.
    using edge_element_t = int;
    using edge_array_t = boost::multi_array<edge_element_t, 2>;

    // Every node type is a single ASCII letter in "abmsw".  Note,
    // there is no 'c' nodes as c-node channels are converted into
    // a-node activities.
    using node_code_t = char;

    // Edge types are formed by tail and head node type codes in
    // alphabetical order with tail occupying the upper 8 bits.
    using edge_code_t = int16_t;

    // Build edge code from two node codes
    inline edge_code_t edge_code(node_code_t nc1, node_code_t nc2)
    {
        if (nc2 < nc1) {
            std::swap(nc1, nc2);
        }
        return ((nc1&0xff)<<8) | (nc2&0xff);
    }
    // Edge code as string
    inline std::string to_string(ClusterArrays::edge_code_t ec) {
        std::string ret;
        ret.push_back((char)( (ec&0xff00) >> 8 ));
        ret.push_back((char)( (ec&0xff) ));
        return ret;
    }


    // A set of node arrays indexed by their node type.
    using node_array_set_t = std::unordered_map<node_code_t, node_array_t>;
    // A set of edge arrays indexed by their ege type.
    using edge_array_set_t = std::unordered_map<edge_code_t, edge_array_t>;

    inline
    std::vector<node_code_t> node_codes(const node_array_set_t& nas) {
        std::vector<node_code_t> ret;
        for (const auto& [code, arr] : nas) {
            ret.push_back(code);
        }
        return ret;
    }
    inline
    std::vector<edge_code_t> edge_codes(const edge_array_set_t& eas) {
        std::vector<ClusterArrays::edge_code_t> ret;
        for (const auto& [code, arr] : eas) {
            ret.push_back(code);
        }
        return ret;
    }
        
    // Return the node's code.  This returns the code in the
    // cluster array schema.  Specifically, 'c' is translated to
    // 'a'.
    inline node_code_t node_code(const cluster_node_t& node) {
        char nc = node.code();
        if (nc == 'c') return 'a';
        return nc;
    }

    // Return new graph same as input but with new channel nodes
    // appended that represent channels from which activity in the
    // slice was collected.  Original channel nodes have their pointer
    // nulled but are otherwise LEFT IN PLACE in order not to
    // invalidate neither node nor edge descriptors.
    cluster_graph_t bodge_channel_slice(cluster_graph_t graph);


    // Fill array sets from cluster.
    void to_arrays(const cluster_graph_t& cgraph,
                   node_array_set_t& nas, edge_array_set_t& eas);

    // Return cluster given array sets.  The anodes must provide the
    // faces referenced by the blobs.
    using anodes_t = std::vector<IAnodePlane::pointer>;
    cluster_graph_t to_cluster(const node_array_set_t& nas,
                               const edge_array_set_t& eas,
                               const anodes_t& anodes);


    /// See TensorDM for conversion between ICluster and ITensor.

}

#endif
