/** ClusterFileSink -> ClusterFileSource (numpy format) round trip of a cluster
 *  that has slices, channels and wires but NO blobs.
 *
 *  Such clusters are real: when a deghosting pass tags every blob of an anode
 *  (seen on an isochronous track in the DUNE FD-HD 1x2x6 workspace,
 *  wcp-porting-img/wcfm/docs/02, 2026-09-24) the sink writes the s/w/a node
 *  arrays and the a-s/a-w edges but no 'b' array, and the loader
 *  (Aux::ClusterArrays to_cluster) used to do nas.at('b') unconditionally and
 *  throw std::out_of_range ("unordered_map::at"), aborting the clustering job.
 *  The measure array was already guarded; the blob array must be too.
 */

#include "WireCellSio/ClusterFileSink.h"
#include "WireCellSio/ClusterFileSource.h"
#include "WireCellAux/SimpleChannel.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleWire.h"
#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cstdio>
#include <cstdlib>
#include <string>

using namespace WireCell;

namespace {
    std::string tmpname(const std::string& base)
    {
        const char* dir = std::getenv("TMPDIR");
        return std::string(dir ? dir : "/tmp") + "/" + base;
    }
}

TEST_CASE("clusterfile numpy round trip of a cluster with no blobs")
{
    // The source needs anodes to resolve faces of blobs; there are none here
    // but the configuration requires a non-empty list.
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    auto anode = anodes[0];

    // s-node with activity on a few channels, each channel with one wire (c-w edges)
    const int nchan = 4;
    cluster_graph_t g;
    ISlice::map_t activity;
    std::vector<cluster_vertex_t> cvs;
    for (int i = 0; i < nchan; ++i) {
        const int chid = 100 + i;
        const WirePlaneId wpid(kWlayer, 0, 0);
        Ray ray(Point(0, -1000, 3 * i), Point(0, 1000, 3 * i));
        auto iwire = std::make_shared<Aux::SimpleWire>(wpid, i, i, chid, ray, 0);
        auto ichan = std::make_shared<Aux::SimpleChannel>(chid, i, IWire::vector{iwire});
        activity[ichan] = ISlice::value_t(10.0 * (i + 1), 1.0);
        auto cv = boost::add_vertex(cluster_node_t(IChannel::pointer(ichan)), g);
        auto wv = boost::add_vertex(cluster_node_t(IWire::pointer(iwire)), g);
        boost::add_edge(cv, wv, g);
        cvs.push_back(cv);
    }
    auto frame = std::make_shared<Aux::SimpleFrame>(100, 0.0);
    auto islice = std::make_shared<Aux::SimpleSlice>(frame, 0, 0.0, 2.0 * units::us, activity);
    auto sv = boost::add_vertex(cluster_node_t(ISlice::pointer(islice)), g);
    for (auto cv : cvs) boost::add_edge(sv, cv, g);
    auto icl = std::make_shared<Aux::SimpleCluster>(g, 100);

    const std::string fname = tmpname("wct-clusterfile-noblobs.tar");
    std::remove(fname.c_str());
    {
        Sio::ClusterFileSink sink;
        auto cfg = sink.default_configuration();
        cfg["outname"] = fname;
        cfg["format"] = "numpy";
        sink.configure(cfg);
        REQUIRE(sink(icl));
        ICluster::pointer eos;
        REQUIRE(sink(eos));
    }

    Sio::ClusterFileSource src;
    auto cfg = src.default_configuration();
    cfg["inname"] = fname;
    cfg["anodes"][0] = "AnodePlane:" + std::to_string(anode->ident());
    src.configure(cfg);
    ICluster::pointer out;
    REQUIRE_NOTHROW(src(out));   // threw std::out_of_range before the fix
    REQUIRE(out);
    CHECK(out->ident() == 100);

    int nb = 0, ns = 0, nc = 0, nw = 0;
    for (const auto& v : boost::make_iterator_range(boost::vertices(out->graph()))) {
        switch (out->graph()[v].code()) {
        case 'b': ++nb; break;
        case 's': ++ns; break;
        case 'c': ++nc; break;
        case 'w': ++nw; break;
        default: break;
        }
    }
    CHECK(nb == 0);
    CHECK(ns == 1);
    CHECK(nc == nchan);
    CHECK(nw == nchan);
}
