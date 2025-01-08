// DO NO USE
//
// This code should not be considered for any user examples.  This
// reproduces some internal configuration stuff so that units tests
// need not specify a configuration file.
//

#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWireSchema.h"
#include "WireCellIface/IDFT.h"

#include <vector>
#include <string>

using namespace WireCell;
using namespace std;

std::vector<std::string> known_dets = {"uboone", "apa", "protodune-larsoft"};
void known_det(std::string maybe)
{
    for (auto det : known_dets) {
        if (maybe == det) {
            return;
        }
    }
    THROW(ValueError() << errmsg{String::format("unknown detector: %s", maybe)});
}
std::vector<std::string> anode_loader(std::string detector)
{
    known_det(detector);

    std::vector<std::string> ret;
    {
        int nanodes = 1;
        // Note: these files must be located via WIRECELL_PATH
        std::string ws_fname = "microboone-celltree-wires-v2.1.json.bz2";
        std::string fr_fname = "ub-10-half.json.bz2";
        if (detector == "uboone") {
            ws_fname = "microboone-celltree-wires-v2.1.json.bz2";
            fr_fname = "ub-10-half.json.bz2";
        }
        if (detector == "apa") {
            ws_fname = "apa-wires.json.bz2";
            fr_fname = "garfield-1d-3planes-21wires-6impacts-dune-v1.json.bz2";
        }
        if (detector == "protodune-larsoft") {
            ws_fname = "protodune-wires-larsoft-v1.json.bz2";
            fr_fname = "garfield-1d-3planes-21wires-6impacts-dune-v1.json.bz2";
            nanodes = 6;
        }

        PluginManager& pm = PluginManager::instance();
        pm.add("WireCellSigProc");
        pm.add("WireCellGen");
        pm.add("WireCellAux");

        const std::string fr_tn = "FieldResponse";
        const std::string ws_tn = "WireSchemaFile";

        {
            auto icfg = Factory::lookup<IConfigurable>(fr_tn);
            auto cfg = icfg->default_configuration();
            cfg["filename"] = fr_fname;
            icfg->configure(cfg);
        }
        {
            auto icfg = Factory::lookup<IConfigurable>(ws_tn);
            auto cfg = icfg->default_configuration();
            cfg["filename"] = ws_fname;
            icfg->configure(cfg);
        }
        {
            // If FftwDFT grows to be an IConfigurable, this needs to
            // change to suit.
            Factory::lookup<IDFT>("FftwDFT"); 
        }

        for (int ianode = 0; ianode < nanodes; ++ianode) {
            std::string tn = String::format("AnodePlane:%d", ianode);
            ret.push_back(tn);
            cerr << "Configuring: " << tn << "\n";
            auto icfg = Factory::lookup_tn<IConfigurable>(tn);
            auto cfg = icfg->default_configuration();
            cfg["ident"] = ianode;
            cfg["wire_schema"] = ws_tn;
            cfg["faces"][0]["anode"] = 0;
            cfg["faces"][0]["response"] = 10 * units::cm - 6 * units::mm;
            cfg["faces"][0]["cathode"] = 2.5604 * units::m;
            cerr << cfg << endl;
            icfg->configure(cfg);
        }
    }
    return ret;
}
