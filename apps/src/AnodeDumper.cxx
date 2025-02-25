#include "WireCellApps/AnodeDumper.h"

#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ConfigManager.h"
#include "WireCellUtil/Logging.h"


WIRECELL_FACTORY(AnodeDumper, WireCellApps::AnodeDumper, WireCell::IApplication, WireCell::IConfigurable)

using spdlog::info;
using spdlog::warn;

using namespace WireCell;
using namespace WireCellApps;

AnodeDumper::AnodeDumper()
  : m_cfg(default_configuration())
{
}

AnodeDumper::~AnodeDumper() {}

void AnodeDumper::configure(const Configuration& config) { m_cfg = config; }

WireCell::Configuration AnodeDumper::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = "/dev/stdout";
    cfg["anodes"] = Json::arrayValue; // list of IAnodePlane type/names
    return cfg;
}

void AnodeDumper::execute()
{
    std::vector<std::string> anode_tns;
    for (auto jone : m_cfg["anodes"]) {
        anode_tns.push_back(jone.asString());
    }

    Configuration sum;

    int anode_index=0;
    for (const auto& anode_tn : anode_tns) {
        const auto& ianode = Factory::find_tn<IAnodePlane>(anode_tn);
        auto& janode = sum["anodes"][anode_index++];

        janode["ident"] = ianode->ident();
        janode["nfaces"] = ianode->faces().size();
        janode["nchannels"] = ianode->channels().size();

        int face_index = 0;
        for (const auto& iface : ianode->faces()) {
            auto& jface = janode["faces"][face_index++];
            jface["ident"] = iface->ident();
            jface["which"] = iface->which();
            jface["dirx"] = iface->dirx();
            jface["aid"] = iface->anode();
            jface["nplanes"] = iface->nplanes();

            int plane_index = 0;
            for (const auto& iplane : iface->planes()) {
                auto& jplane = jface["planes"][plane_index++];
                jplane["ident"] = iplane->ident();
                jplane["wpid"] = iplane->planeid().ident();
                jplane["nchannels"] = iplane->channels().size();
                jplane["nwires"] = iplane->wires().size();
            }
        }
    }

    Persist::dump(get<std::string>(m_cfg, "filename"), sum, true);
}
