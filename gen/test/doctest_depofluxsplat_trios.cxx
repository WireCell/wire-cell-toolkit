/** Tests for DepoFluxSplat's optional depo-to-channel ("trio") output.

    A depo is seen by all three planes at once, so the channels it lights in
    each plane share a cause.  The frame sums every depo's contribution and so
    discards which depo went where; trio_file recovers that association.

    The anode here is generated rather than loaded, so these tests need no
    external data files.
 */

#include "WireCellGen/DepoFluxSplat.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDepoFramer.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/ITrace.h"

#include "WireCellUtil/DetectorWires.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/NumpyHelper.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/doctest.h"

#include <map>
#include <set>
#include <string>

using namespace WireCell;

namespace {

    const double the_tick = 0.5 * units::us;
    const double the_speed = 1.6 * units::mm / units::us;
    const double the_origin = 10 * units::cm;

    // A uboone-like three plane anode written to a temporary file, so the test
    // carries no dependence on installed data.  Returns the anode's type:name.
    std::string make_anode(const std::string& wires_file)
    {
        WireSchema::StoreDB storedb;
        for (int ipl = 0; ipl < 3; ++ipl) {
            auto& plane = WireSchema::get_append(storedb, ipl, 0, 0, 0);
            DetectorWires::plane(storedb, plane, DetectorWires::uboone[ipl]);
        }
        DetectorWires::flat_channels(storedb);
        WireSchema::Store store(std::make_shared<WireSchema::StoreDB>(storedb));
        WireSchema::dump(wires_file.c_str(), store);

        auto& pm = PluginManager::instance();
        pm.add("WireCellGen");

        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = wires_file;
            icfg->configure(cfg);
        }
        {
            auto icfg = Factory::lookup_tn<IConfigurable>("AnodePlane:0");
            auto cfg = icfg->default_configuration();
            cfg["ident"] = 0;
            cfg["wire_schema"] = "WireSchemaFile";
            cfg["faces"][0]["anode"] = 0.0;
            cfg["faces"][0]["response"] = the_origin;
            cfg["faces"][0]["cathode"] = 2.5604 * units::m;
            icfg->configure(cfg);
        }
        return "AnodePlane:0";
    }

    Configuration splat_config(const std::string& anode_tn,
                               const std::string& trio_file = "")
    {
        auto icfg = Factory::lookup<IConfigurable>("DepoFluxSplat");
        auto cfg = icfg->default_configuration();
        cfg["anode"] = anode_tn;
        cfg["drift_speed"] = the_speed;
        cfg["response_plane"] = the_origin;
        cfg["tick"] = the_tick;
        cfg["window_start"] = 0.0;
        cfg["window_duration"] = 400 * units::us;
        cfg["trio_file"] = trio_file;
        return cfg;
    }

    // Depos spread over the sensitive volume with enough extent that each
    // covers several wires and ticks -- the case the tiling exists for.
    IDepoSet::pointer make_depos(size_t n = 12)
    {
        IDepo::vector depos;
        for (size_t i = 0; i < n; ++i) {
            const double y = (-400.0 + 80.0 * i) * units::mm;
            const double z = (1000.0 + 250.0 * i) * units::mm;
            depos.push_back(std::make_shared<Aux::SimpleDepo>(
                                10 * units::us, Point(1 * units::m, y, z),
                                -1000.0, nullptr,
                                2.0 * units::mm,     // extent_long
                                2.5 * units::mm));   // extent_tran
        }
        return std::make_shared<Aux::SimpleDepoSet>(0, depos);
    }

    // (channel, tick) of every non-zero sample of a frame.
    std::set<std::pair<int,int>> frame_pixels(const IFrame::pointer& frame)
    {
        std::set<std::pair<int,int>> out;
        for (const auto& trace : *frame->traces()) {
            const int chid = trace->channel();
            const int tbin = trace->tbin();
            const auto& q = trace->charge();
            for (size_t i = 0; i < q.size(); ++i) {
                if (q[i] != 0.0) out.insert({chid, tbin + (int)i});
            }
        }
        return out;
    }

    struct Trios {
        Array::array_xxi index;
        Eigen::ArrayXXd value;
        size_t size() const { return (size_t)index.rows(); }
    };

    Trios load_trios(const std::string& fname, int ident = 0)
    {
        Trios t;
        Numpy::load2d(t.index, "trio_index_" + std::to_string(ident), fname);
        Numpy::load2d(t.value, "trio_value_" + std::to_string(ident), fname);
        return t;
    }
}

TEST_SUITE("depofluxsplat trios") {

TEST_CASE("trios land on pixels the frame actually lit")
{
    // This is the test that matters.  A trio asserts that three channels share
    // a cause; if one of them is not even in the frame this component produced,
    // the association describes something other than the truth it ships with.
    Persist::TempDir td;
    const std::string wires_file = (td.path() / "wires.json.bz2").string();
    const std::string trio_file = (td.path() / "trios.npz").string();
    const auto anode_tn = make_anode(wires_file);

    auto icfg = Factory::lookup<IConfigurable>("DepoFluxSplat");
    icfg->configure(splat_config(anode_tn, trio_file));
    auto splat = Factory::find<IDepoFramer>("DepoFluxSplat");

    IFrame::pointer frame;
    REQUIRE(splat->operator()(make_depos(), frame));
    REQUIRE(frame);

    auto trios = load_trios(trio_file, frame->ident());
    REQUIRE(trios.size() > 0);
    CHECK(trios.index.cols() == 6);
    CHECK(trios.value.rows() == trios.index.rows());

    const auto pixels = frame_pixels(frame);
    REQUIRE(!pixels.empty());

    size_t missing = 0;
    for (size_t irow = 0; irow < trios.size(); ++irow) {
        for (int ipl = 0; ipl < 3; ++ipl) {
            const int tbin = trios.index(irow, ipl);
            const int chid = trios.index(irow, 3 + ipl);
            if (!pixels.count({chid, tbin})) ++missing;
        }
    }
    CHECK(missing == 0);

    // Charge is positive-definite here (all depos same sign) and every trio
    // must carry some, or the row is noise.
    for (size_t irow = 0; irow < trios.size(); ++irow) {
        CHECK(trios.value(irow, 0) != 0.0);
    }
}

TEST_CASE("the three channels of a trio are distinct planes")
{
    Persist::TempDir td;
    const std::string wires_file = (td.path() / "wires.json.bz2").string();
    const std::string trio_file = (td.path() / "trios.npz").string();
    const auto anode_tn = make_anode(wires_file);

    auto icfg = Factory::lookup<IConfigurable>("DepoFluxSplat");
    icfg->configure(splat_config(anode_tn, trio_file));
    auto splat = Factory::find<IDepoFramer>("DepoFluxSplat");
    auto anode = Factory::find_tn<IAnodePlane>(anode_tn);

    IFrame::pointer frame;
    REQUIRE(splat->operator()(make_depos(), frame));
    auto trios = load_trios(trio_file, frame->ident());
    REQUIRE(trios.size() > 0);

    // Map channel -> plane index, so a trio can be checked to span U, V and W
    // rather than, say, naming three wires of one plane.
    std::map<int,int> chan2plane;
    for (auto face : anode->faces()) {
        for (auto plane : face->planes()) {
            for (auto wire : plane->wires()) {
                chan2plane[wire->channel()] = plane->planeid().index();
            }
        }
    }

    for (size_t irow = 0; irow < trios.size(); ++irow) {
        for (int ipl = 0; ipl < 3; ++ipl) {
            const int chid = trios.index(irow, 3 + ipl);
            REQUIRE(chan2plane.count(chid));
            CHECK(chan2plane[chid] == ipl);
        }
    }
}

TEST_CASE("off by default: no file and an unchanged frame")
{
    Persist::TempDir td;
    const std::string wires_file = (td.path() / "wires.json.bz2").string();
    const std::string trio_file = (td.path() / "unwanted.npz").string();
    const auto anode_tn = make_anode(wires_file);

    auto icfg = Factory::lookup<IConfigurable>("DepoFluxSplat");
    auto splat = Factory::find<IDepoFramer>("DepoFluxSplat");

    icfg->configure(splat_config(anode_tn, trio_file));
    IFrame::pointer with;
    REQUIRE(splat->operator()(make_depos(), with));

    icfg->configure(splat_config(anode_tn, ""));   // trio_file unset
    IFrame::pointer without;
    REQUIRE(splat->operator()(make_depos(), without));

    CHECK(!Persist::exists((td.path() / "never-written.npz").string()));

    // Asking for trios must not perturb the frame at all.
    CHECK(frame_pixels(with) == frame_pixels(without));
}

TEST_CASE("time_offsets shift each plane's tick independently")
{
    // The three channels of a trio need not be hot at the same tick, which is
    // why the row carries a tick per plane rather than one shared tick.
    Persist::TempDir td;
    const std::string wires_file = (td.path() / "wires.json.bz2").string();
    const auto anode_tn = make_anode(wires_file);
    auto icfg = Factory::lookup<IConfigurable>("DepoFluxSplat");
    auto splat = Factory::find<IDepoFramer>("DepoFluxSplat");

    const std::string plain_file = (td.path() / "plain.npz").string();
    icfg->configure(splat_config(anode_tn, plain_file));
    IFrame::pointer f0;
    REQUIRE(splat->operator()(make_depos(), f0));
    auto plain = load_trios(plain_file, f0->ident());

    const std::string offset_file = (td.path() / "offset.npz").string();
    auto cfg = splat_config(anode_tn, offset_file);
    cfg["time_offsets"][0] = 0.0;
    cfg["time_offsets"][1] = 2 * the_tick;
    cfg["time_offsets"][2] = 5 * the_tick;
    icfg->configure(cfg);
    IFrame::pointer f1;
    REQUIRE(splat->operator()(make_depos(), f1));
    auto offset = load_trios(offset_file, f1->ident());

    REQUIRE(plain.size() > 0);
    REQUIRE(plain.size() == offset.size());

    // Same channels, ticks shifted by exactly the configured offsets.
    for (size_t irow = 0; irow < plain.size(); ++irow) {
        for (int ipl = 0; ipl < 3; ++ipl) {
            CHECK(offset.index(irow, 3 + ipl) == plain.index(irow, 3 + ipl));
        }
        CHECK(offset.index(irow, 0) - plain.index(irow, 0) == 0);
        CHECK(offset.index(irow, 1) - plain.index(irow, 1) == 2);
        CHECK(offset.index(irow, 2) - plain.index(irow, 2) == 5);
    }
}

}
