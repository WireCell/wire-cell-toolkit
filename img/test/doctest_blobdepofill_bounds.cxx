/** BlobDepoFill must give a blob no charge from a depo that lies outside the
 *  blob's extent along the primary wire.
 *
 *  The fill integrates each depo's transverse Gaussian along the primary (W)
 *  wire between the blob's bounds [wlo, whi] (from its U/V strips) clipped to
 *  the depo's +-nsigma window.  When the two do not overlap the clipped bounds
 *  are inverted, and gbounds() -- by design symmetric in its arguments --
 *  returned the Gaussian mass of the GAP between the window edge and the blob:
 *  0.5*erfc(nsigma/sqrt 2) ~ 1.3e-3 per (depo, wire) at nsigma 3.  Every blob
 *  elsewhere on the same W wire in the slice therefore received tens of
 *  electrons per depo -- "diffusion tails" in W-wire columns hundreds of cm
 *  from the track, 94 % of the "true" cells of the FD-HD 1x2x6 truth tiers
 *  (wcp-porting-img/wcfm/docs/05 sec 8, 2026-09-25).
 *
 *  The test tiles one slice with TWO separated U/V activity groups sharing
 *  the same W wires, puts one depo at the centre of a blob of the first group
 *  and requires every blob of the second group (> 50 cm away along the wire)
 *  to receive exactly 0.  Before the fix each of them received ~1e-3 of the
 *  depo's charge.
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

    // One slice: U and V active in a few separated groups, W active over a wide range, so
    // the tiling makes blobs at several places along the same W wires.
    ISlice::pointer make_slice(IAnodeFace::pointer face, double start, double span, size_t shift)
    {
        auto frame = std::make_shared<Aux::SimpleFrame>(100, 0.0);
        ISlice::map_t activity;
        const auto& planes = face->planes();
        for (size_t ip = 0; ip < planes.size(); ++ip) {
            const auto& chans = planes[ip]->channels();
            const size_t mid = chans.size() / 2;
            std::vector<std::pair<size_t, size_t>> ranges;
            if (ip < 2) {
                ranges = {{mid, mid + 8}, {mid + shift, mid + shift + 8}, {mid - shift, mid - shift + 8}};
            }
            else {
                ranges = {{mid - shift, mid + shift}};
            }
            for (const auto& [lo, hi] : ranges) {
                for (size_t i = lo; i < std::min(hi, chans.size()); ++i) {
                    activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
                }
            }
        }
        return std::make_shared<Aux::SimpleSlice>(frame, 0, start, span, activity);
    }

    Point centroid(const IAnodeFace::pointer& face, const IBlob::pointer& blob)
    {
        const auto& coords = face->raygrid();
        Point cen(0, 0, 0);
        int n = 0;
        for (const auto& [a, b] : blob->shape().corners()) {
            cen += coords.ray_crossing(a, b);
            ++n;
        }
        return n ? cen / n : cen;
    }
}

TEST_CASE("blobdepofill gives nothing to blobs outside the depo's extent along the wire")
{
    auto anode = make_fdhd_anode(8);
    REQUIRE(anode);
    auto face = anode->face(0);
    REQUIRE(face);
    const double span = 2.0 * units::us;
    auto slice = make_slice(face, 0.0, span, 200);

    Img::GridTiling gt;
    auto gcfg = gt.default_configuration();
    gcfg["anode"] = "AnodePlane:8";
    gcfg["face"] = 0;
    gt.configure(gcfg);
    IBlobSet::pointer bs;
    REQUIRE(gt(slice, bs));
    REQUIRE(bs);
    const auto blobs = bs->blobs();
    REQUIRE(blobs.size() >= 2);

    // the depo: at the centroid of blob 0, just in front of the response plane
    const Point cen0 = centroid(face, blobs[0]);
    const double xresp = face->planes()[2]->pimpos()->origin()[0];
    const Point pos(xresp + face->dirx() * 1.0 * units::mm, cen0.y(), cen0.z());
    const double q = 100000.0;
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
    fcfg["nsigma"] = 3.0;
    fill.configure(fcfg);
    ICluster::pointer out;
    REQUIRE(fill(std::make_tuple(icl, depos), out));
    REQUIRE(out);

    // W wires of the depo's blob
    const auto w0 = blobs[0]->shape().strips()[4].bounds;
    double near = 0, far = 0;
    int nfar = 0, nfar_sharedw = 0, nfar_nonzero = 0;
    for (const auto& v : boost::make_iterator_range(boost::vertices(out->graph()))) {
        const auto& node = out->graph()[v];
        if (node.code() != 'b') continue;
        auto b = std::get<IBlob::pointer>(node.ptr);
        const double d = (centroid(face, b) - cen0).magnitude();
        if (d > 50 * units::cm) {
            ++nfar;
            far += b->value();
            const auto w = b->shape().strips()[4].bounds;
            if (w.first < w0.second && w0.first < w.second) ++nfar_sharedw;
            if (b->value() != 0.0) ++nfar_nonzero;
        }
        else {
            near += b->value();
        }
    }
    REQUIRE(nfar > 0);
    REQUIRE(nfar_sharedw > 0);      // far blobs on the depo's own W wires: the leak's victims
    CHECK(near > 0.5 * q);          // the depo's blob got the charge
    CHECK(nfar_nonzero == 0);       // before the fix: every far blob on the depo's W wires > 0
    CHECK(far == 0.0);
}
