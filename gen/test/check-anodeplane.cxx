#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"


using namespace WireCell;

int main(int argc, char* argv[])
{
    std::string detname = "pdsp";
    if (argc > 1) {
        detname = argv[1];
    }



    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");

    {
        auto wfile = Persist::resolve(detname, "wires");
        if (wfile.empty()) {
            std::cerr << "No wires file for " << detname << "\n";
            return 1;
        }
        auto icfg = Factory::lookup_tn<IConfigurable>("WireSchemaFile");
        auto cfg = icfg->default_configuration();
        cfg["filename"] = wfile;
        icfg->configure(cfg);
    }
    {
        auto icfg = Factory::lookup_tn<IConfigurable>("AnodePlane");
        auto cfg = icfg->default_configuration();
        cfg["ident"] = 0;
        cfg["wire_schema"] ="WireSchemaFile";
        // bugs
        cfg["faces"][0]["anode"] = 0;
        cfg["faces"][0]["response"] = -10*units::cm;
        cfg["faces"][0]["cathode"] = -1*units::m;
        cfg["faces"][1]["anode"] = 0;
        cfg["faces"][1]["response"] = +10*units::cm;
        cfg["faces"][1]["cathode"] = +1*units::m;
        icfg->configure(cfg);
    }

    auto anode = Factory::find_tn<IAnodePlane>("AnodePlane");
    if (not anode) {
        return 1;
    }

    for (auto face : anode->faces()) {
        std::cerr << "face ident=" << face->ident() << " which=" << face->which() << " dirx=" << face->dirx() << "\n";
        for (auto plane: face->planes()) {
            const auto& all_chans = plane->channels();
            IChannel::vector end_chans = {all_chans.front(), all_chans.back()};
            const auto& chans = end_chans;
            const auto& wires = plane->wires();
            std::cerr << "plane ident=" << plane->ident() << " nchannels=" << chans.size() << " nwires=" << wires.size() << "\n";
            for (auto chan : chans) {
                auto wire0 = chan->wires()[0];
                std::cerr << "chan ident=" << chan->ident() << " index="<<chan->ident()
                          << " wire0id=" << wire0->ident() << " wire0ind=" << wire0->index()
                          << " wire0head=" << wire0->ray().second
                          << "\n";
            }
        }
    }
    return 0;
}
