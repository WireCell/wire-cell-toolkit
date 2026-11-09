// Test helper: the anode of APA 0 of the DUNE FD HD 1x2x6 geometry, built
// from its wires file found on WIRECELL_PATH, or nullptr if the file is not
// found.  Also provides a BoxRecombination.

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"

namespace trig_test {

    inline WireCell::IAnodePlane::pointer fdhd_anode()
    {
        using namespace WireCell;
        const std::string wires = "dune10kt-1x2x6-wires-larsoft-v1.json.bz2";
        if (Persist::resolve(wires).empty()) return nullptr;
        static IAnodePlane::pointer anode;
        if (anode) return anode;

        PluginManager::instance().add("WireCellGen");
        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile", "trigtest");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = wires;
            icfg->configure(cfg);
        }
        {
            auto icfg = Factory::lookup<IConfigurable>("AnodePlane", "trigtest");
            auto cfg = icfg->default_configuration();
            cfg["ident"] = 0;
            cfg["wire_schema"] = "WireSchemaFile:trigtest";
            // as cfg/pgrapher/experiment/dune10kt-1x2x6/params.jsonnet
            const double anode_x = 39.5355 * units::mm, response_x = 130.0155 * units::mm,
                         cathode_x = 3629.1625 * units::mm;
            cfg["faces"][0]["anode"] = anode_x;
            cfg["faces"][0]["response"] = response_x;
            cfg["faces"][0]["cathode"] = cathode_x;
            cfg["faces"][1]["anode"] = -anode_x;
            cfg["faces"][1]["response"] = -response_x;
            cfg["faces"][1]["cathode"] = -cathode_x;
            icfg->configure(cfg);
        }
        anode = Factory::find_tn<IAnodePlane>("AnodePlane:trigtest");
        return anode;
    }

    inline std::string box_recombination()
    {
        using namespace WireCell;
        PluginManager::instance().add("WireCellGen");
        auto icfg = Factory::lookup<IConfigurable>("BoxRecombination", "trigtest");
        icfg->configure(icfg->default_configuration());
        return "BoxRecombination:trigtest";
    }
}  // namespace trig_test
