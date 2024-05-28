/** An ICluster is predominantly a graph of typed nodes.
 *
 * Node types are shared pointers to types in a fixed set of types
 * related to WC 3D imaging.  An ICluster is more than a "cluster of
 * blobs" (aka more than a "geometric" or "trajectory" cluster).  Node
 * types with their code letters inlcude:
 *
 * - c :: IChannel
 * - w :: IWire
 * - b :: IBlob 
 * - m :: IMeasure
 * - s :: ISlice
 *
 * An ICluster is indirectly associated to one or more IFrames through
 * its ISlice nodes.  Thus it may "span" less than, more than or
 * exactly one IFrame.
 */

#ifndef WIRECELL_ICLUSTER
#define WIRECELL_ICLUSTER

#include "WireCellIface/IChannel.h"
#include "WireCellIface/IWire.h"
#include "WireCellIface/IBlob.h"
#include "WireCellIface/ISlice.h"
#include "WireCellIface/IMeasure.h"
#include "WireCellIface/WirePlaneId.h"

#include "WireCellUtil/IndexedGraph.h"

#include "WireCellUtil/Graph.h"



#include <variant>

namespace WireCell {

    struct cluster_node_t {
        using channel_t = IChannel::pointer;
        using wire_t = IWire::pointer;
        using blob_t = IBlob::pointer;
        using slice_t = ISlice::pointer;
        using meas_t = IMeasure::pointer;
        using ptr_t = std::variant<size_t, channel_t, wire_t, blob_t, slice_t, meas_t>;

        ptr_t ptr;

        cluster_node_t()
          : ptr()
        {
        }
        cluster_node_t(const ptr_t& p)
          : ptr(p)
        {
        }
        cluster_node_t(const channel_t& p)
          : ptr(p)
        {
        }
        cluster_node_t(const wire_t& p)
          : ptr(p)
        {
        }
        cluster_node_t(const blob_t& p)
          : ptr(p)
        {
        }
        cluster_node_t(const slice_t& p)
          : ptr(p)
        {
        }
        cluster_node_t(const meas_t& p)
          : ptr(p)
        {
        }

        cluster_node_t(const cluster_node_t& other)
          : ptr(other.ptr)
        {
        }

        // A node type is represented by a single character code.  The
        // codes are in a canonical order.  Note, the index to a code
        // is one less than the ptr.index().
        const static std::string known_codes; // "cwbsm"

        // Return the index of a known code in [0,4].  The index of a
        // code is one less than ptr.index() and an illegal code
        // returns index beyond the known_codes() string (ie 5).
        static size_t code_index(char code);

        // Helper: return the letter code for the type of the ptr or
        // \0 if code is type is undefined.
        inline char code() const
        {
            auto ind = ptr.index();
            if (ind == std::variant_npos) {
                return 0;
            }
            return known_codes[ind-1];
        }

        bool operator==(const cluster_node_t& other) const { return ptr == other.ptr; }
        cluster_node_t& operator=(const cluster_node_t& other)
        {
            ptr = other.ptr;
            return *this;
        }

        // Access the underlying IData ident number, eg to use to
        // produce an ordering.
        int ident() const
        {
            using ident_f = std::function<int()>;
            std::vector<ident_f> ofs {
                [](){return 0;}, 
                [&](){return std::get<channel_t>(ptr)->ident();},
                [&](){return std::get<wire_t>(ptr)->ident();},
                [&](){return std::get<blob_t>(ptr)->ident();},
                [&](){return std::get<slice_t>(ptr)->ident();},
                [&](){return std::get<meas_t>(ptr)->ident();},
            };
            const auto ind = ptr.index();
            if (ind == std::variant_npos) {
                return 0;
            }
            return ofs[ind]();
        }
    };                          // struct cluster_node_t

}  // namespace WireCell

namespace std {
    template <>
    struct hash<WireCell::cluster_node_t> {
        std::size_t operator()(const WireCell::cluster_node_t& n) const
        {
            size_t h = 0;
            switch (n.ptr.index()) {
            case 0:
                h = std::get<0>(n.ptr);
                break;
            case 1:
                h = (std::size_t) std::get<1>(n.ptr).get();
                break;
            case 2:
                h = (std::size_t) std::get<2>(n.ptr).get();
                break;
            case 3:
                h = (std::size_t) std::get<3>(n.ptr).get();
                break;
            case 4:
                h = (std::size_t) std::get<4>(n.ptr).get();
                break;
            case 5:
                h = (std::size_t) std::get<5>(n.ptr).get();
                break;
            }
            return h;
        }
    };
}  // namespace std

namespace WireCell {

    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, cluster_node_t> cluster_graph_t;
    typedef boost::graph_traits<cluster_graph_t>::vertex_descriptor cluster_vertex_t;
    typedef boost::graph_traits<cluster_graph_t>::edge_descriptor cluster_edge_t;
    typedef boost::graph_traits<cluster_graph_t>::vertex_iterator cluster_vertex_iter_t;

    // The actual ICluster interface.
    //
    // It is small and essentially delievers a cluster_grapht_.  All
    // the rest of the code defines that type with some helpers.
    class ICluster : public IData<ICluster> {
       public:
        virtual ~ICluster();

        /// Return an identifying number.
        virtual int ident() const = 0;

        // Access the graph.
        virtual const cluster_graph_t& graph() const = 0;
    };

    typedef IndexedGraph<cluster_node_t> cluster_indexed_graph_t;

    template <typename Type>
    inline std::vector<Type> oftype(const cluster_indexed_graph_t& g)
    {
        std::vector<Type> ret;
        for (const auto& v : boost::make_iterator_range(boost::vertices(g.graph()))) {
            const auto& vp = g.graph()[v];
            if (std::holds_alternative<Type>(vp.ptr)) {
                ret.push_back(std::get<Type>(vp.ptr));
            }
        }
        return ret;
    }

    // 
    template <typename Type>
    inline std::vector<Type> neighbors_oftype(const cluster_indexed_graph_t& g, const cluster_node_t& n)
    {
        std::vector<Type> ret;
        for (const auto& vp : g.neighbors(n)) {
            if (std::holds_alternative<Type>(vp.ptr)) {
                ret.push_back(std::get<Type>(vp.ptr));
            }
        }
        return ret;
    }

    // descriptor version of oftype
    inline std::vector<cluster_vertex_t> oftype(const WireCell::cluster_graph_t& cg, const char typecode)
    {
        std::vector<cluster_vertex_t> ret;
        for (const auto& vtx : boost::make_iterator_range(boost::vertices(cg))) {
            if (cg[vtx].code() != typecode) continue;
            ret.push_back(vtx);
        }
        return ret;
    }

    // descriptor version of neighbors/neighbors_oftype
    inline std::vector<cluster_vertex_t> neighbors(const WireCell::cluster_graph_t& cg, const cluster_vertex_t& vd)
    {
        std::vector<cluster_vertex_t> ret;
        for (auto edge : boost::make_iterator_range(boost::out_edges(vd, cg))) {
            cluster_vertex_t neigh = boost::target(edge, cg);
            ret.push_back(neigh);
        }
        return ret;
    }

    template <typename Type>
    inline std::vector<cluster_vertex_t> neighbors_oftype(const WireCell::cluster_graph_t& cg, const cluster_vertex_t& vd)
    {
        std::vector<cluster_vertex_t> ret;
        for (const auto& vp : neighbors(cg, vd)) {
            if (std::holds_alternative<Type>(cg[vp].ptr)) {
                ret.push_back(vp);
            }
        }
        return ret;
    }

}  // namespace WireCell

#endif
