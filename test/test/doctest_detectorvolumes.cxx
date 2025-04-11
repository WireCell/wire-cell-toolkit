#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IConfigurable.h"

#include <chrono>

using namespace WireCell;
using spdlog::debug;

// largely from:
//
// wcsonnet apps/test/anode-dumper.jsonnet
//
/// I don't think these "faces" are even right but it doesn't matter for this test.
std::string json_config = R"(
[
    {
	"data" : 
	    {
		"filename" : "protodune-wires-larsoft-v4.json.bz2"
	    },
	"name" : "protodune-wires-larsoft-v4.json.bz2",
	"type" : "WireSchemaFile"
    },
    {
	"data" : 
	    {
		"faces" : 
		    [
			{
			    "anode" : -3574.5999999999999,
			    "cathode" : -2584.6000000000004,
			    "response" : -3484.5999999999999
			},
			{
			    "anode" : -3699.2000000000003,
			    "cathode" : -4689.1999999999998,
			    "response" : -3789.2000000000003
			}
		    ],
		"ident" : 0,
		"wire_schema" : "WireSchemaFile:protodune-wires-larsoft-v4.json.bz2"
	    },
	"name" : "0",
	"type" : "AnodePlane"
    },
    {
	"data" : 
	    {
		"faces" : 
		    [
			{
			    "anode" : 3699.2000000000003,
			    "cathode" : 4689.1999999999998,
			    "response" : 3789.2000000000003
			},
			{
			    "anode" : 3574.5999999999999,
			    "cathode" : 2584.6000000000004,
			    "response" : 3484.5999999999999
			}
		    ],
		"ident" : 1,
		"wire_schema" : "WireSchemaFile:protodune-wires-larsoft-v4.json.bz2"
	    },
	"name" : "1",
	"type" : "AnodePlane"
    },
    {
	"data" : 
	    {
		"faces" : 
		    [
			{
			    "anode" : -3574.5999999999999,
			    "cathode" : -2584.6000000000004,
			    "response" : -3484.5999999999999
			},
			{
			    "anode" : -3699.2000000000003,
			    "cathode" : -4689.1999999999998,
			    "response" : -3789.2000000000003
			}
		    ],
		"ident" : 2,
		"wire_schema" : "WireSchemaFile:protodune-wires-larsoft-v4.json.bz2"
	    },
	"name" : "2",
	"type" : "AnodePlane"
    },
    {
	"data" : 
	    {
		"faces" : 
		    [
			{
			    "anode" : 3699.2000000000003,
			    "cathode" : 4689.1999999999998,
			    "response" : 3789.2000000000003
			},
			{
			    "anode" : 3574.5999999999999,
			    "cathode" : 2584.6000000000004,
			    "response" : 3484.5999999999999
			}
		    ],
		"ident" : 3,
		"wire_schema" : "WireSchemaFile:protodune-wires-larsoft-v4.json.bz2"
	    },
	"name" : "3",
	"type" : "AnodePlane"
    },
    {
	"data" : 
	    {
		"faces" : 
		    [
			{
			    "anode" : -3574.5999999999999,
			    "cathode" : -2584.6000000000004,
			    "response" : -3484.5999999999999
			},
			{
			    "anode" : -3699.2000000000003,
			    "cathode" : -4689.1999999999998,
			    "response" : -3789.2000000000003
			}
		    ],
		"ident" : 4,
		"wire_schema" : "WireSchemaFile:protodune-wires-larsoft-v4.json.bz2"
	    },
	"name" : "4",
	"type" : "AnodePlane"
    },
    {
	"data" : 
	    {
		"faces" : 
		    [
			{
			    "anode" : 3699.2000000000003,
			    "cathode" : 4689.1999999999998,
			    "response" : 3789.2000000000003
			},
			{
			    "anode" : 3574.5999999999999,
			    "cathode" : 2584.6000000000004,
			    "response" : 3484.5999999999999
			}
		    ],
		"ident" : 5,
		"wire_schema" : "WireSchemaFile:protodune-wires-larsoft-v4.json.bz2"
	    },
	"name" : "5",
	"type" : "AnodePlane"
    },
    {
        "type": "DetectorVolumes",
        "name": "",
        "data": {
            "anodes": ["AnodePlane:0","AnodePlane:1","AnodePlane:2","AnodePlane:3","AnodePlane:4","AnodePlane:5"],
            "metadata": {"default":"default", 
                        "a0f0pV":"a0f0pV", "a0f0pU":"a0f0pU", "a0f0pW":"a0f0pW",
                        "a0f0pA":"a0f0pA", "overall":"overall"}

        }
    }
]
)";

TEST_CASE("detectorvolumes")
{
    PluginManager& pm = PluginManager::instance();
    // fixme, we need to move anodeplane, etc into aux
    pm.add("WireCellSigProc");
    pm.add("WireCellGen");
    pm.add("WireCellAux");
    
    auto cfg = Persist::loads(json_config);

    for (const auto& c : cfg) {
        auto type = get<std::string>(c, "type");
        auto name = get<std::string>(c, "name");
        auto iface = Factory::lookup<Interface>(type, name);
    }
    for (auto c : cfg) {
        auto type = get<std::string>(c, "type");
        auto name = get<std::string>(c, "name");
        auto cfgobj = Factory::find_maybe<IConfigurable>(type, name);
        if (!cfgobj) {
            continue;
        }
        Configuration cfg = cfgobj->default_configuration();
        cfg = update(cfg, c["data"]);
        cfgobj->configure(cfg);  // throws
    }

    auto dv = Factory::lookup_tn<IDetectorVolumes>("DetectorVolumes");
    REQUIRE(dv);
    auto fv = Factory::lookup_tn<IFiducial>("DetectorVolumes");
    REQUIRE(fv);

    // $ wirecell-util wires-info pdsp
    // anode:0 face:0 X=[-3594.16,-3584.63]mm Y=[76.10,6066.70]mm Z=[0.00,2306.73]mm
    //	0: x=-3584.63mm dx=9.5250mm n=1148 pitch=(4.6691 +/- 0.000006 [4.6684<4.6693], p0=4.6691) 
    //	1: x=-3589.39mm dx=4.7620mm n=1148 pitch=(4.6691 +/- 0.000005 [4.6684<4.6693], p0=4.6691) 
    //	2: x=-3594.16mm dx=0.0000mm n=480 pitch=(4.7920 +/- 0.000000 [4.7920<4.7920], p0=4.7920)
    CHECK(false == dv->contained_by(Point(0,0,0)).valid());
    CHECK(false == fv->contained(Point(0,0,0)));
    CHECK(false == dv->contained_by(Point(-3400*units::mm, 0, 0)).valid());
    CHECK(false == fv->contained(Point(-3400*units::mm, 0, 0)));
    auto wpid = dv->contained_by(Point(-3500*units::mm, 100*units::mm, 100*units::mm));
    CHECK(true == wpid.valid());
    CHECK(wpid.apa() == 0);
    CHECK(wpid.face() == 0);
    CHECK(fv->contained(Point(-3500*units::mm, 100*units::mm, 100*units::mm)));

//	std::cout << "haha " << fv->contained(Point(-3600*units::mm, 100*units::mm, 100*units::mm)) << " " << dv->is_in_overall_volume(Point(-3600*units::mm, 100*units::mm, 100*units::mm)) << std::endl;

   // inside the overall detector volume
   CHECK(false == fv->contained(Point(-3600*units::mm, 100*units::mm, 100*units::mm)) );	

    {
        const size_t ncalls = 100000;

        Point outside(-3400*units::mm, 0, 0);
        auto t0 = std::chrono::steady_clock::now();
        for (size_t count = 0; count < ncalls; ++count) {
            fv->contained(outside);
        }
        auto t1 = std::chrono::steady_clock::now();

        Point inside(-3500*units::mm, 100*units::mm, 100*units::mm);
        for (size_t count = 0; count < ncalls; ++count) {
            fv->contained(inside);
        }
        auto t2 = std::chrono::steady_clock::now();

        std::chrono::duration<double> dt1 = t1-t0;
        std::chrono::duration<double> dt2 = t2-t1;

        debug("n={}: out={}s {}Hz, in={}s {}Hz",
              ncalls,
              dt1.count(), ncalls/dt1.count(),
              dt2.count(), ncalls/dt2.count());

        // Original implementation with no compiler optimization:
        // n=100000: out=0.110162781s  907747.5994365102Hz, in=0.1092375s    915436.548804211Hz
        // n=100000: out=0.005630255s 17761184.884165995Hz, in=0.020494938s 4879253.599108228Hz
    }

    auto wpid_u = wpid.to_u();
    CHECK(wpid_u.index() == 0);
    CHECK(wpid_u.valid());
    // CHECK(wpid_u == true);

    auto wdir = dv->wire_direction(wpid_u);
    debug("wire direction: {}", wdir);
    CHECK(wdir.magnitude() == 1.0);
    auto pitch = dv->pitch_vector(wpid_u);
    debug("pitch vector: {}", pitch);
    CHECK(pitch.magnitude() > 1*units::mm);

    {                           // test metadata
        std::vector<WirePlaneId> wpids = {
            WirePlaneId(WirePlaneLayer_t::kUlayer, 0, 0),
            WirePlaneId(WirePlaneLayer_t::kVlayer, 0, 0),
            WirePlaneId(WirePlaneLayer_t::kWlayer, 0, 0),
            WirePlaneId(WirePlaneLayer_t::kAllLayers, 0, 0),
        };

        for (const auto& wpid : wpids) {
            auto md = dv->metadata(wpid);
            std::string have = md.asString();
            std::string want = wpid.name();
            debug("metadata: have:{} want:{}", have, want);
            CHECK(have == want);
        }
        {
            WirePlaneId wpid(kUnknownLayer, 1, 0); // not in configuration
            auto md = dv->metadata(wpid);
            std::string have = md.asString();
            std::string want = "default";
            debug("metadata: have:{} want:{}", have, want);
            CHECK(have == want);
        }
        {
            WirePlaneId wpid(0); // overall
            auto md = dv->metadata(wpid);
            std::string have = md.asString();
            std::string want = "overall";
            debug("metadata: have:{} want:{}", have, want);
            CHECK(have == want);
        }

    }

}

TEST_CASE("envfiducial") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellSigProc");
    pm.add("WireCellGen");
    pm.add("WireCellAux");
    
    // This test "cheats a little by reusing the same config as for testing
    // DetectorVolumes but change the type EnvFiducial.  The extra "metadata"
    // attribute that DetectorVolumes accepts is not honored and ignored by
    // EnvFiducial.
    auto cfg = Persist::loads(json_config);
    cfg[cfg.size()-1]["type"] = "EnvFiducial";

    for (const auto& c : cfg) {
        auto type = get<std::string>(c, "type");
        auto name = get<std::string>(c, "name");
        auto iface = Factory::lookup<Interface>(type, name);
    }
    for (auto c : cfg) {
        auto type = get<std::string>(c, "type");
        auto name = get<std::string>(c, "name");
        auto cfgobj = Factory::find_maybe<IConfigurable>(type, name);
        if (!cfgobj) {
            continue;
        }
        Configuration cfg = cfgobj->default_configuration();
        cfg = update(cfg, c["data"]);
        cfgobj->configure(cfg);  // throws
    }

    auto fv = Factory::lookup_tn<IFiducial>("EnvFiducial");
    REQUIRE(fv);

    // EnvFiducial bounds: [(-4689.2 76.1 3.3497) --> (4689.2 6060 6944.58)]
    {
        CHECK(! fv->contained(Point(0,0,0)));
        CHECK(! fv->contained(Point(10,10,10)));
        CHECK(  fv->contained(Point(0, 100, 10)));
        CHECK(! fv->contained(Point(-5000, 100, 10)));
        CHECK(! fv->contained(Point( 5000, 100, 10)));
    }

}
