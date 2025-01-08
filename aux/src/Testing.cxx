#include "WireCellAux/Testing.h"

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/IWireSchema.h"

#include <algorithm>            // find

using namespace WireCell;
using namespace WireCell::Aux;

WireCell::IDFT::pointer Testing::get_dft()
{
    // load_plugins({"WireCellAux"});
    // we are in Aux, so no need to load!

    return Factory::lookup<IDFT>("FftwDFT"); // not configurable
}


WireCell::IRandom::pointer Testing::get_random(const std::string& name)
{
    using namespace WireCell;
    load_plugins({"WireCellGen"}); // fixme: should move Random to Aux.
    auto rngcfg = Factory::lookup<IConfigurable>("Random", name);
    auto cfg = rngcfg->default_configuration();
    cfg["generator"] = name;    // default or twister
    rngcfg->configure(cfg);
    return Factory::find<IRandom>("Random", name);
}


IAnodePlane::vector Testing::anodes(std::string detector)
{
    // Probably should move WSF and FR to aux....
    Testing::load_plugins({"WireCellGen", "WireCellSigProc"});

    auto kds = detectors();
    const auto& kd = kds[detector];
    if (kd.isNull()) {
        raise<ValueError>("unknown detector \"%s\"", detector);
    }        

    const std::string ws_tn = "WireSchemaFile";
    {
        auto icfg = Factory::lookup<IConfigurable>(ws_tn);
        auto cfg = icfg->default_configuration();
        auto wires = kd["wires"];
        // std::cerr << "using wires file: " << wires << "\n";
        cfg["filename"] = wires;
        icfg->configure(cfg);
    }

    const std::string fr_tn = "FieldResponse";
    {
        auto icfg = Factory::lookup<IConfigurable>(fr_tn);
        auto cfg = icfg->default_configuration();
        auto fields = kd["fields"];
        if (fields.isNull()) {
          raise<ValueError>("unknown fields file for detector \"%s\"", detector);
        }
        if (fields.isArray()) {
          fields = fields[0];
        }
        cfg["filename"] = fields;
        icfg->configure(cfg);
    }

    auto iwsf = Factory::lookup_tn<IWireSchema>(ws_tn);
    const auto& wstore = iwsf->wire_schema_store();
    const auto& wdets = wstore.detectors();
    const auto& wanodes = wstore.anodes();

    IAnodePlane::vector ret;
    for (const auto& det : wdets) {
        for (const auto& anode_index : det.anodes) {
            const auto& anode = wanodes[anode_index];
            const int ianode = anode.ident;

            // std::cerr << "anode index: " << anode_index << " anode ident: " << anode.ident << "\n";
            
            std::string tn = String::format("AnodePlane:%d", ianode);
            auto icfg = Factory::lookup_tn<IConfigurable>(tn);
            auto cfg = icfg->default_configuration();
            cfg["ident"] = ianode;
            cfg["wire_schema"] = ws_tn;

            // FIXME: this is NOT general and needs to be retrieved from detectors.jsonnet somehow!!!!
            cfg["faces"][0]["anode"] = 10 * units::cm - 6 * units::mm;
            cfg["faces"][0]["response"] = 10 * units::cm - 6 * units::mm;
            cfg["faces"][0]["cathode"] = 2.5604 * units::m;

            icfg->configure(cfg);
            // std::cerr << "Configured AnodePlane with: " << cfg << "\n";

            ret.push_back(Factory::find_tn<IAnodePlane>(tn));
        }
    }
    return ret;
}


