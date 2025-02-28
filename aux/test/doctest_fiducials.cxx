#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Persist.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IConfigurable.h"
using namespace WireCell;
using spdlog::debug;

TEST_CASE("aux boxfiducial") {
    // PluginManager& pm = PluginManager::instance();
    // pm.add("WireCellAux");

    {
        auto icfg = Factory::lookup<IConfigurable>("BoxFiducial");
        auto cfg = icfg->default_configuration();
        cfg["bounds"]["tail"]["x"] = 0;
        cfg["bounds"]["tail"]["y"] = 0;
        cfg["bounds"]["tail"]["z"] = 0;
        cfg["bounds"]["head"]["x"] = 10;
        cfg["bounds"]["head"]["y"] = 10;
        cfg["bounds"]["head"]["z"] = 10;
        icfg->configure(cfg);
    }

    {
        auto fv = Factory::lookup_tn<IFiducial>("BoxFiducial");

        CHECK(fv->contained(Point(0,0,0)));
        CHECK(fv->contained(Point(10,10,10)));
        CHECK(! fv->contained(Point(10.1,10.1,10.1)));
        CHECK(! fv->contained(Point(-0.1,-0.1,-0.1)));
    }

}


TEST_CASE("aux polyfiducial") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellAux");

    const std::string jsonnet_config = R"(
[
    {                           // each slab
        local npts = n + 3,
        local ang = (2*3.1415)/npts,

        min: n,
        max: n+1,
        corners:[
            [std.cos(ang*i), std.sin(ang*i)]
            for i in std.range(0, npts-1)]
    } for n in std.range(0,5)
]
)";

    {
        auto icfg = Factory::lookup<IConfigurable>("PolyFiducial");
        auto cfg = icfg->default_configuration();
        cfg["slabs"] = Persist::loads(jsonnet_config);
        icfg->configure(cfg);
    }

    {
        auto fv = Factory::lookup_tn<IFiducial>("PolyFiducial");

        CHECK(! fv->contained(Point(-1,0,0))); // before on axis
        CHECK(  fv->contained(Point( 0,0,0))); // low edge on axis
        CHECK(  fv->contained(Point( 6,0,0))); // high edge on axis
        CHECK(! fv->contained(Point(6.1,0,0))); // after on axis

        // first slab is a triangle with a corner at (y=1,z=0)
        CHECK(  fv->contained(Point(0.5, 0.9999, 0)) );
        CHECK(! fv->contained(Point(0.5, 1.0001, 0)) );
        CHECK(! fv->contained(Point(0.5, 1.0, 0.1)) );

    }

}
