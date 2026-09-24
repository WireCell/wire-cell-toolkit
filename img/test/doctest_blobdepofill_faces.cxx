/** BlobDepoFill must fill blobs on BOTH faces of an anode.
 *
 *  Its "depo is behind the face" test compared the raw x offset of a depo
 *  from the wire plane's Pimpos origin (the response plane) with zero.  That
 *  is only right for a face whose drift volume lies at +x.  On a -x face
 *  (dirx = -1) every depo in the volume has x below the response plane, so
 *  every depo was skipped and every blob of the face came out with zero true
 *  charge -- found on the DUNE FD-HD 1x2x6 geometry, whose APAs drift both
 *  ways (wcp-porting-img/wcfm/docs/02, 2026-09-24).  A second edge: a depo
 *  drifted exactly onto the response plane and round-tripped through a
 *  float32 depo file reads back a few um off, on either side, so the test
 *  also needs a tolerance.
 *
 *  The test tiles the same small activity on face 0 and on face 1 of one
 *  FD-HD APA, puts one depo a few um BEHIND each face's response plane at the
 *  centre of a blob, and requires both faces to receive the depo's charge.
 *  Before the fix face 1 received exactly 0.
 */

#include "WireCellImg/BlobDepoFill.h"
#include "WireCellImg/GridTiling.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>

using namespace WireCell;

namespace {

    IAnodePlane::pointer make_fdhd_anode(int ident)
    {
        Testing::load_plugins({"WireCellGen", "WireCellSigProc"});
        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = "dune10kt-1x2x6-wires-larsoft-v1.json.bz2";
            icfg->configure(cfg);
        }
        const std::string tn = String::format("AnodePlane:%d", ident);
        auto icfg = Factory::lookup_tn<IConfigurable>(tn);
        auto cfg = icfg->default_configuration();
        cfg["ident"] = ident;
        cfg["wire_schema"] = "WireSchemaFile";
        // both faces live, mirror images: wires at |x| ~ 30-40 mm, drift to +-3.63 m
        const double apa_plane = 39.5355 * units::mm, res_plane = 130.0155 * units::mm, cpa = 3629.1625 * units::mm;
        cfg["faces"][0]["anode"] = apa_plane;
        cfg["faces"][0]["response"] = res_plane;
        cfg["faces"][0]["cathode"] = cpa;
        cfg["faces"][1]["anode"] = -apa_plane;
        cfg["faces"][1]["response"] = -res_plane;
        cfg["faces"][1]["cathode"] = -cpa;
        icfg->configure(cfg);
        return Factory::find_tn<IAnodePlane>(tn);
    }

    ISlice::pointer make_slice(IAnodeFace::pointer face, double start, double span)
    {
        auto frame = std::make_shared<Aux::SimpleFrame>(100, 0.0);
        ISlice::map_t activity;
        for (const auto& plane : face->planes()) {
            const auto& chans = plane->channels();
            const size_t mid = chans.size() / 2;
            for (size_t i = mid; i < std::min(mid + 12, chans.size()); ++i) {
                activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
            }
        }
        return std::make_shared<Aux::SimpleSlice>(frame, 0, start, span, activity);
    }

    // Fill the blobs tiled on `face` from one depo of charge q placed `behind` (along
    // -dirx) the response plane at the centre of the first blob; return the summed fill.
    double filled_charge(IAnodePlane::pointer anode, int iface, double behind, double q)
    {
        auto face = anode->face(iface);
        REQUIRE(face);
        const double span = 2.0 * units::us;
        auto slice = make_slice(face, 0.0, span);

        Img::GridTiling gt;
        auto gcfg = gt.default_configuration();
        gcfg["anode"] = String::format("AnodePlane:%d", anode->ident());
        gcfg["face"] = iface;
        gt.configure(gcfg);
        IBlobSet::pointer bs;
        REQUIRE(gt(slice, bs));
        REQUIRE(bs);
        const auto blobs = bs->blobs();
        REQUIRE(blobs.size() > 0);

        // the depo: y,z at the centroid of blob 0's corners, x just behind the response plane
        const auto& coords = face->raygrid();
        Point cen(0, 0, 0);
        int ncorner = 0;
        for (const auto& [a, b] : blobs[0]->shape().corners()) {
            cen += coords.ray_crossing(a, b);
            ++ncorner;
        }
        REQUIRE(ncorner >= 3);
        cen = cen / ncorner;
        const double xresp = face->planes()[2]->pimpos()->origin()[0];
        const Point pos(xresp - face->dirx() * behind, cen.y(), cen.z());
        auto depo = std::make_shared<Aux::SimpleDepo>(0.5 * span, pos, q, nullptr, 0.1 * units::mm, 1.0 * units::mm);
        auto depos = std::make_shared<Aux::SimpleDepoSet>(0, IDepo::vector{depo});

        cluster_graph_t g;
        for (const auto& b : blobs) {
            boost::add_vertex(cluster_node_t(IBlob::pointer(b)), g);
        }
        auto icl = std::make_shared<Aux::SimpleCluster>(g, 100);

        Img::BlobDepoFill fill;
        auto fcfg = fill.default_configuration();
        fcfg["speed"] = 1.6 * units::mm / units::us;
        fcfg["time_offset"] = 0.0;
        fill.configure(fcfg);
        ICluster::pointer out;
        REQUIRE(fill(std::make_tuple(icl, depos), out));
        REQUIRE(out);

        double sum = 0;
        for (const auto& v : boost::make_iterator_range(boost::vertices(out->graph()))) {
            const auto& node = out->graph()[v];
            if (node.code() != 'b') continue;
            sum += std::get<IBlob::pointer>(node.ptr)->value();
        }
        return sum;
    }
}

TEST_CASE("blobdepofill fills both faces of an anode")
{
    auto anode = make_fdhd_anode(8);
    REQUIRE(anode);
    REQUIRE(anode->face(0)->dirx() == 1);
    REQUIRE(anode->face(1)->dirx() == -1);

    const double q = 5000.0;
    // 3 um behind the response plane: the float32 round-trip offset seen in practice
    const double f0 = filled_charge(anode, 0, 3e-6 * units::mm, q);
    const double f1 = filled_charge(anode, 1, 3e-6 * units::mm, q);
    CHECK(f0 > 0.5 * q);
    CHECK(f1 > 0.5 * q);                 // was exactly 0 before the fix
    CHECK(std::abs(f0 - f1) < 0.05 * q); // the two faces are mirror images

    // a depo well behind the face (between the wires and the response plane) is still skipped
    CHECK(filled_charge(anode, 0, 50 * units::mm, q) == 0.0);
    CHECK(filled_charge(anode, 1, 50 * units::mm, q) == 0.0);
}
