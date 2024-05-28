#include "WireCellImg/BlobGrouping.h"

#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/SimpleMeasure.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/Graph.h"



WIRECELL_FACTORY(BlobGrouping, WireCell::Img::BlobGrouping,
                 WireCell::INamed,
                 WireCell::IClusterFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::GraphTools;


// blob/channel graph with original graph vertices on per-slice basis.
namespace bcdesc {
    struct nprop_t { cluster_vertex_t desc; };
    struct eprop_t { };

    using graph_t = boost::adjacency_list<boost::setS, 
                                          boost::vecS, 
                                          boost::undirectedS,
                                          nprop_t, eprop_t>;

    using vdesc_t = boost::graph_traits<graph_t>::vertex_descriptor;
    using edesc_t = boost::graph_traits<graph_t>::edge_descriptor;
}


static
void doit(cluster_graph_t& cgraph)
{
    // Count to set measure idents.
    int tot_meas = 0;

    for (const auto& svtx : vertex_range(cgraph)) {
        const auto& snode = cgraph[svtx];
        const char scode = snode.code();
        if (scode != 's') {
            continue;
        }

        const auto islice = std::get<ISlice::pointer>(snode.ptr);
        ISlice::map_t activity = islice->activity();

        // Recieve b and c reached from this s per plane
        std::vector<bcdesc::graph_t> bcs(3); // fixme: hard-code 3 planes
        std::vector< std::unordered_map<int, bcdesc::vdesc_t> > uniq_chans(3);

        for (auto bvtx : mir(boost::adjacent_vertices(svtx, cgraph))) {
            const auto& bnode = cgraph[bvtx];
            const char bcode = bnode.code();
            if (bcode != 'b') {
                continue;
            }

            const auto iblob = std::get<IBlob::pointer>(bnode.ptr);

            // add this blob to each plane graph
            std::vector<bcdesc::vdesc_t> bidx;
            for (auto& bc : bcs) {
                bidx.push_back(boost::add_vertex({ bvtx }, bc));
            }

            // wires
            for (auto wvtx : mir(boost::adjacent_vertices(bvtx, cgraph))) {
                const auto& wnode = cgraph[wvtx];
                const char wcode = wnode.code();
                if (wcode != 'w') {
                    continue;
                }
                
                const auto iwire = std::get<IWire::pointer>(wnode.ptr);

                // channels
                for (auto cvtx : mir(boost::adjacent_vertices(wvtx, cgraph))) {
                    const auto& cnode = cgraph[cvtx];
                    const char ccode = cnode.code();
                    if (ccode != 'c') {
                        continue;
                    }

                    auto ich = std::get<IChannel::pointer>(cnode.ptr);
                    const auto wpid = ich->planeid();
                    const auto pind = wpid.index();
                    const auto cident = ich->ident();

                    auto& bc = bcs[pind];

                    bcdesc::vdesc_t cidx;
                    auto& uniq_chan = uniq_chans[pind];
                    auto ucit = uniq_chan.find(cident);
                    if (ucit == uniq_chan.end()) {
                        cidx = boost::add_vertex({ cvtx }, bc);
                        uniq_chan[cident] = cidx;
                    }
                    else {
                        cidx = ucit->second;
                    }

                    boost::add_edge(bidx[pind], cidx, bc);
                } // channels
            } // wires
        } // blobs
        
        // now have a per-slice, per-plane graph of blob-channel vertices
        for (size_t iplane = 0; iplane < 3; ++iplane) {
            const auto& bc = bcs[iplane];

            size_t nverts = boost::num_vertices(bc);

            // every cc's int locally names a measure
            // std::unordered_map<bcdesc::vdesc_t, int> cc;
            // auto ncomponents = boost::connected_components(bc, boost::make_assoc_property_map(cc));
            std::vector<int> cc(nverts);
            auto ncomponents = boost::connected_components(bc, &cc[0]);
            
            std::vector<Aux::SimpleMeasure*> measures(ncomponents, nullptr);
            std::vector<cluster_vertex_t> mvtxs(ncomponents);

            for (size_t idx = 0; idx < nverts; ++idx) {
            // for (const auto& [idx, ncomp] : cc) {
                size_t ncomp = cc[idx];
                Aux::SimpleMeasure* sm=measures[ncomp];
                cluster_vertex_t mvtx=mvtxs[ncomp];
                if (!sm) {
                    sm = measures[ncomp] = new Aux::SimpleMeasure(tot_meas++);
                    mvtx = mvtxs[ncomp] = boost::add_vertex(IMeasure::pointer(sm), cgraph);
                }

                // edge to measure for either blob or channel
                const auto vtx = bc[idx].desc;
                boost::add_edge(vtx, mvtx, cgraph);

                // if channel, accumulate signal
                const auto& node = cgraph[vtx];
                if (node.code() == 'c') {
                    auto ich = std::get<IChannel::pointer>(node.ptr);
                    const auto sig = activity[ich];
                    const auto wpid = ich->planeid();
                    sm->sig += sig;
                    sm->wpid = wpid;
                }
            }
        } // planes
    } // slices
}


Img::BlobGrouping::BlobGrouping()
    : Aux::Logger("BlobGrouping", "img")
{
}

Img::BlobGrouping::~BlobGrouping() {}

void Img::BlobGrouping::configure(const WireCell::Configuration& cfg) {}

WireCell::Configuration Img::BlobGrouping::default_configuration() const
{
    WireCell::Configuration cfg;
    return cfg;
}

using Filtered = boost::filtered_graph<cluster_graph_t, boost::keep_all, boost::function<bool(cluster_vertex_t)> >;

bool Img::BlobGrouping::operator()(const input_pointer& in, output_pointer& out)
{
    if (!in) {
        log->debug("EOS");
        out = nullptr;
        ++m_count;
        return true;
    }

    const cluster_graph_t& ingr = in->graph();

    // Remove measures and copy to make a writable graph to receive
    // new measures.
    Filtered fg(ingr, {}, [&](cluster_vertex_t vtx) {
        return ingr[vtx].code() != 'm';
    });
    cluster_graph_t outgr;
    boost::copy_graph(fg, outgr);

    doit(outgr);

    log->debug("call={} cluster={} nvertices={} nedges={}",
               m_count,
               in->ident(),
               boost::num_vertices(outgr),
               boost::num_edges(outgr));

    out = std::make_shared<Aux::SimpleCluster>(outgr, in->ident());
    ++m_count;
    return true;


    // new way
    // [13:38:06.958] I [ timer  ] Timer: WireCell::Img::BlobGrouping : 8.557092 sec

    // old way, using IndexedGraph ...
    // [12:23:35.901] I [ timer  ] Timer: WireCell::Img::BlobGrouping : 40.23048 sec
}
