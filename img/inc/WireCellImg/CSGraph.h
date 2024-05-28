/**
   CSGraph - charge solving graph

   This is used internally by ChargeSolving.

 */

#ifndef WIRECELLIMG_CSGRAPH
#define WIRECELLIMG_CSGRAPH

#include "WireCellIface/ICluster.h"
#include "WireCellUtil/IndexedSet.h"
#include "WireCellUtil/Ress.h"
#include "WireCellUtil/Eigen.h"

#include <iterator>

namespace WireCell::Img::CS {

    
    // A value holds (central,uncertainty) for measures and blobs.
    // For blobs, uncertainty is interpeted as the weight.
    using value_t = ISlice::value_t;

    // For matrix representation of the graphs.
    using double_vector_t = Eigen::VectorXd;
    using double_matrix_t = Eigen::MatrixXd;

    // Here, "node" implies a cluster graph vertex payload object.
    using channel_t = cluster_node_t::channel_t;
    using wire_t = cluster_node_t::wire_t;
    using blob_t = cluster_node_t::blob_t;
    using slice_t = cluster_node_t::slice_t;
    using meas_t = cluster_node_t::meas_t;

    // The type associated with a graph_t as a whole
    struct graphprop_t {
        // The input islice from which this graph was derived
        slice_t islice;

        // An arbitrary index counting the original connected subgraph
        // from the slice from which this graph derives
        size_t index{0};

        // The two chi2 components from a solution.
        double chi2_base{0};
        double chi2_l1{0};
    };

    // The vertex payload type associated to a charge solving graph
    // vertex.
    struct node_t {
        // Remember where we came from in the original cluster graph.
        cluster_vertex_t orig_desc;

        // Mark what kind of the vertex this is.  
        enum Kind { unknown=0, meas, blob };
        Kind kind;

        // A number by which a collection of vertices of the same kind may
        // be ordered
        int ordering;

        // Both m-kind and b-kind vertex has a value as
        // (central,uncertainty) pairs.  For blob type, the central value
        // may be used on input to the solver as a starting value and
        // holds the result of the solving on output.  The uncertainty is
        // interpreted as a weight and is untouched by solving.  For meas
        // type, the value is the sum of the per channel values.
        value_t value;
    };

    //
    // The grpah type.  Just b's and m's holding ready to use info
    //
    using graph_t = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                          node_t, boost::no_property, graphprop_t>;
    using vdesc_t = boost::graph_traits<graph_t>::vertex_descriptor;
    using edesc_t = boost::graph_traits<graph_t>::edge_descriptor;
    using vertex_iter_t = boost::graph_traits<graph_t>::vertex_iterator;
    
    // make easier iteration.
    // for (auto v : mir(boost::vertices(graph))) {
    //    for (auto edge : mir(boost::out_edges(v, graph)) { ... }
    // }
    template <typename It> boost::iterator_range<It> mir(std::pair<It, It> const& p) {
        return boost::make_iterator_range(p.first, p.second);
    }
    template <typename It> boost::iterator_range<It> mir(It b, It e) {
        return boost::make_iterator_range(b, e);
    }
    template <typename Gr> 
    boost::iterator_range<typename boost::graph_traits<Gr>::vertex_iterator> vertex_range(Gr& g) {
        return mir(boost::vertices(g));
    }

    // Used, eg when forming subgraphs.
    using graph_vector_t = std::vector<graph_t>;

    // Return selection of vertex descriptions of kind ordered by
    // .ordering and as an indexed set.
    using indexed_vdescs_t = IndexedSet<CS::vdesc_t>;
    indexed_vdescs_t select_ordered(const graph_t& csg,
                                    node_t::Kind kind);
                           
    // Break up input CS graph into "connected subgraphs".
    void connected_subgraphs(const graph_t& slice_graph,
                             std::back_insert_iterator<graph_vector_t> subgraphs_out);

    // Break a cluster_graph into a vector of connected subgraphs
    // represented by a CS graph_t.  Measurements with value less than
    // meas_thresh are ignored.
    void unpack(const cluster_graph_t& cgraph,
                std::back_insert_iterator<graph_vector_t> subgraphs,
                const value_t& meas_thresh);

    // Reform a cluster graph from a collection of CS graphs.
    cluster_graph_t repack(const cluster_graph_t& cgraph,
                           const std::vector<graph_t>& csgs);

    // Solve the graph, returning a new one holding solution in blob nodes
    struct SolveParams {
        // affect
        // 1, fitting model config
        // 2, how to process single blob graph
        enum Config{uboone=42,simple};
        Config config{uboone};
        double scale{1000};
        // If true, apply matrix transform of the cholesky
        // decomposition on the measure covariance matrix.  If False,
        // measurement uncertainties are not considered.
        bool whiten{true};
    };
    graph_t solve(const graph_t& csg, const SolveParams& params, const bool verbose=false);

    // Return a copy of the input retaining only blobs with signal
    // value above or equal to threshold.  Any blobless measures are
    // also removed.
    graph_t prune(const graph_t& csg, float threshold=0);
}

namespace std {
    template <>
    struct hash<WireCell::Img::CS::node_t> {
        std::size_t operator()(const WireCell::Img::CS::node_t& n) const {
            return n.ordering;
        }
    };
}
#endif
