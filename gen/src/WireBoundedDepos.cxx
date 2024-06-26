#include "WireCellGen/WireBoundedDepos.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IWirePlane.h"

#include <boost/range.hpp>


WIRECELL_FACTORY(WireBoundedDepos,
                 WireCell::Gen::WireBoundedDepos,
                 WireCell::INamed, WireCell::IDrifter, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Gen::WireBoundedDepos::WireBoundedDepos()
  : Aux::Logger("WireBoundedDepos", "gen")
{
}
Gen::WireBoundedDepos::~WireBoundedDepos()
{
}

WireCell::Configuration Gen::WireBoundedDepos::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = "";

    // List of wire regions
    // [{
    //   plane: <plane-number>,
    //   min: <min-wire-number>,
    //   max: <max-wire-number,
    // },...]
    cfg["wires"] = Json::arrayValue;

    cfg["mode"] = "accept";

    return cfg;
}

void Gen::WireBoundedDepos::configure(const WireCell::Configuration& cfg)
{
    // we just want the pimpos but best to keep a handle on the anode
    // since it's shared_ptr memory.
    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());

    for (auto face : m_anode->faces()) {
        if (face->planes().empty()) {
            continue;
        }
        for (auto plane : face->planes()) {
            const size_t iplane = plane->ident();
            if (m_pimpos.size() <= iplane) {
                m_pimpos.resize(iplane + 1);
            }
            m_pimpos[iplane] = plane->pimpos();
        }
        break;  // fixme:
    }
    m_accept = get<std::string>(cfg, "mode", "accept") == "accept";

    auto jregions = cfg["regions"];
    for (auto jregion : jregions) {
        wire_region_t wr;
        for (auto jtrio : jregion) {
            wr.push_back(wire_bounds_t(jtrio["plane"].asInt(), jtrio["min"].asInt(), jtrio["max"].asInt()));
        }
        m_regions.push_back(wr);
    }

    log->debug("accept={} nregions={} nplanes={}", m_accept, m_regions.size(), m_pimpos.size());
}

bool Gen::WireBoundedDepos::operator()(const input_pointer& depo, output_queue& outq)
{
    if (!depo) {  // EOS protocol
        outq.push_back(nullptr);
        return true;
    }

    // lazily calculate nearest wire
    const size_t nplanes = m_pimpos.size();
    std::vector<bool> already(nplanes, false);
    std::vector<int> closest_wire(nplanes, -1);

    for (const auto& region : m_regions) {
        bool inregion = true;
        for (const auto& trio : region) {
            const int iplane = get<0>(trio);
            const int imin = get<1>(trio);
            const int imax = get<2>(trio);

            if (!already[iplane]) {  // lazy
                const double pitch = m_pimpos[iplane]->distance(depo->pos(), 2);
                closest_wire[iplane] = m_pimpos[iplane]->region_binning().bin(pitch);
                already[iplane] = true;
            }
            const int iwire = closest_wire[iplane];

            if (iwire < imin or iwire > imax) {
                inregion = false;
                break;
            }
        }

        if (inregion) {
            if (m_accept) {
                outq.push_back(depo);
            }
            // accept or reject, we landed this depo in a configured
            // region.
            return true;
        }
    }
    if (!m_accept) {  // depo missed all rejections.
        outq.push_back(depo);
    }
    return true;
}

