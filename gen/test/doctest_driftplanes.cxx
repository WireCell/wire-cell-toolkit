// AnodePlane and Drifter may take their drift-axis planes (anode, response,
// cathode) from an IDriftPlanes component instead of from their own
// configuration so that both use one source, such as the geometry of a host
// application.

#include "WireCellGen/Drifter.h"

#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDriftPlanes.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

namespace {

    const double anode_x = 1 * units::m;
    const double response_x = 1.1 * units::m;
    const double cathode_x = 3 * units::m;

    Configuration face_planes()
    {
        Configuration jf;
        jf["anode"] = anode_x;
        jf["response"] = response_x;
        jf["cathode"] = cathode_x;
        return jf;
    }

    // Configure the "DriftPlanes" provider: the given anodes each with the
    // front face sensitive and the back face not.
    void configure_drift_planes(const std::vector<int>& idents)
    {
        auto icfg = Factory::lookup<IConfigurable>("DriftPlanes");
        auto cfg = icfg->default_configuration();
        for (int ident : idents) {
            Configuration jvol;
            jvol["wires"] = ident;
            jvol["faces"][0] = face_planes();
            jvol["faces"][1] = Json::nullValue;
            cfg["volumes"].append(jvol);
        }
        icfg->configure(cfg);
    }
}  // namespace

TEST_CASE("gen anodeplane takes faces from drift planes")
{
    auto anodes = Testing::anodes("pdsp");  // configures the "WireSchemaFile" wire schema
    REQUIRE(anodes.size() > 0);
    const int ident = anodes[0]->ident();
    configure_drift_planes({ident});

    auto icfg = Factory::lookup<IConfigurable>("AnodePlane", "driftplanes");
    auto cfg = icfg->default_configuration();
    cfg["ident"] = ident;
    cfg["wire_schema"] = "WireSchemaFile";

    SUBCASE("planes come from the provider")
    {
        cfg["drift_planes"] = "DriftPlanes";
        icfg->configure(cfg);
        auto anode = Factory::find<IAnodePlane>("AnodePlane", "driftplanes");

        auto front = anode->face(0);
        REQUIRE(front);
        const Ray bb = front->sensitive().bounds();
        CHECK(std::min(bb.first.x(), bb.second.x()) == doctest::Approx(anode_x));
        CHECK(std::max(bb.first.x(), bb.second.x()) == doctest::Approx(cathode_x));
        CHECK(front->dirx() == +1);
        for (auto plane : front->planes()) {
            CHECK(plane->pimpos()->origin().x() == doctest::Approx(response_x));
        }

        auto back = anode->face(1);
        REQUIRE(back);
        CHECK(back->sensitive().empty());
    }

    SUBCASE("inline faces and a provider are exclusive")
    {
        cfg["drift_planes"] = "DriftPlanes";
        cfg["faces"][0] = face_planes();
        CHECK_THROWS(icfg->configure(cfg));
    }
}

TEST_CASE("gen drifter takes xregions from drift planes")
{
    PluginManager::instance().add("WireCellGen");
    configure_drift_planes({0, 1});  // same planes: one xregion
    {
        auto icfg = Factory::lookup<IConfigurable>("Random");
        icfg->configure(icfg->default_configuration());
    }

    Gen::Drifter drifter;
    auto cfg = drifter.default_configuration();
    cfg["drift_planes"] = "DriftPlanes";

    SUBCASE("regions come from the provider")
    {
        drifter.configure(cfg);

        IDrifter::output_queue out;
        auto inside = std::make_shared<Aux::SimpleDepo>(0, Point(2 * units::m, 0, 0), -1000);
        auto beyond = std::make_shared<Aux::SimpleDepo>(0, Point(4 * units::m, 0, 0), -1000);
        drifter(inside, out);
        drifter(beyond, out);
        drifter(nullptr, out);

        std::vector<IDepo::pointer> drifted;
        for (const auto& d : out) {
            if (d) drifted.push_back(d);
        }
        REQUIRE(drifted.size() == 1);
        CHECK(drifted[0]->pos().x() == doctest::Approx(response_x));
        CHECK(drifted[0]->prior() == inside);
    }

    SUBCASE("inline xregions and a provider are exclusive")
    {
        cfg["xregions"][0] = face_planes();
        CHECK_THROWS(drifter.configure(cfg));
    }
}

TEST_CASE("gen driftplanes from configuration")
{
    PluginManager::instance().add("WireCellGen");
    configure_drift_planes({3, 7});
    auto dp = Factory::find<IDriftPlanes>("DriftPlanes");

    auto xr = dp->xregions();
    REQUIRE(xr.size() == 1);
    CHECK(xr[0]["response"].asDouble() == doctest::Approx(response_x));

    auto f = dp->faces(7);
    REQUIRE(f.size() == 2);
    CHECK(f[0]["cathode"].asDouble() == doctest::Approx(cathode_x));
    CHECK(f[1].isNull());

    CHECK_THROWS(dp->faces(5));
}
