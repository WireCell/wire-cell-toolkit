#include "WireCellGen/DepoTransformYZ.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

WIRECELL_FACTORY(DepoTransformYZ, WireCell::Gen::DepoTransformYZ,
                 WireCell::INamed,
                 WireCell::IDepoFramer, WireCell::IConfigurable)

using namespace WireCell;

Gen::DepoTransformYZ::DepoTransformYZ() {}
Gen::DepoTransformYZ::~DepoTransformYZ() {}

WireCell::Configuration Gen::DepoTransformYZ::default_configuration() const
{
    auto cfg = DepoTransform::default_configuration();
    // type:name of the IYZMap service (required)
    cfg["yzmap"] = "";
    // anode key in the map; empty means the "anode" type:name
    cfg["yzmap_key"] = "";
    return cfg;
}

void Gen::DepoTransformYZ::configure(const WireCell::Configuration& cfg)
{
    DepoTransform::configure(cfg);

    const std::string yzmap_tn = get<std::string>(cfg, "yzmap", "");
    if (yzmap_tn.empty()) {
        THROW(ValueError() << errmsg{"DepoTransformYZ: \"yzmap\" (an IYZMap type:name) is required"});
    }
    m_yzmap = Factory::find_tn<IYZMap>(yzmap_tn);

    const std::string anode_tn = get<std::string>(cfg, "anode", "");
    m_key = get<std::string>(cfg, "yzmap_key", "");
    if (m_key.empty()) {
        m_key = anode_tn;
    }

    // Fail at configure, not with silently zeroed charge: IYZMap
    // returns 0 for an anode or plane that the map lacks.
    auto anode = Factory::find_tn<IAnodePlane>(anode_tn);
    for (auto face : anode->faces()) {
        if (!face) continue;
        const auto center = face->sensitive().center();
        for (auto plane : face->planes()) {
            const int iplane = plane->planeid().index();
            const double val = m_yzmap->value(m_key, iplane, center.y(), center.z());
            if (val == 0.0) {
                THROW(ValueError() << errmsg{
                        "DepoTransformYZ: map \"" + yzmap_tn + "\" has no non-zero entry for key \""
                        + m_key + "\" plane " + std::to_string(iplane)});
            }
            log->debug("yzmap {} key \"{}\" face {} plane {}: value {} at (y,z)=({:.1f},{:.1f}) cm",
                       yzmap_tn, m_key, face->ident(), iplane, val,
                       center.y() / units::cm, center.z() / units::cm);
        }
    }
}

IDepo::pointer Gen::DepoTransformYZ::modify_depo(WirePlaneId wpid, IDepo::pointer depo)
{
    const auto& pos = depo->pos();
    const double scale = m_yzmap->value(m_key, wpid.index(), pos.y(), pos.z());
    if (scale == 1.0) {
        return depo;
    }
    return std::make_shared<Aux::SimpleDepo>(depo->time(), pos, depo->charge() * scale, depo,
                                              depo->extent_long(), depo->extent_tran(),
                                              depo->id(), depo->pdg(), depo->energy());
}
