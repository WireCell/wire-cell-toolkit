#include "WireCellGen/DriftPlanes.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(DriftPlanes, WireCell::Gen::DriftPlanes,
                 WireCell::IDriftPlanes, WireCell::IConfigurable)

using namespace WireCell;

Gen::DriftPlanes::DriftPlanes() {}
Gen::DriftPlanes::~DriftPlanes() {}

WireCell::Configuration Gen::DriftPlanes::default_configuration() const
{
    Configuration cfg;
    cfg["volumes"] = Json::arrayValue;
    return cfg;
}

void Gen::DriftPlanes::configure(const WireCell::Configuration& cfg)
{
    m_faces.clear();
    m_xregions = Json::arrayValue;
    for (const auto& jvol : cfg["volumes"]) {
        const int ident = jvol["wires"].asInt();
        const auto& jfaces = jvol["faces"];
        if (m_faces.count(ident)) {
            raise<ValueError>("DriftPlanes: anode %d given more than once", ident);
        }
        m_faces[ident] = jfaces;
        for (const auto& jface : jfaces) {
            if (jface.isNull()) continue;
            bool seen = false;
            for (const auto& jxr : m_xregions) {
                if (jxr == jface) {
                    seen = true;
                    break;
                }
            }
            if (!seen) m_xregions.append(jface);
        }
    }
}

Configuration Gen::DriftPlanes::faces(int anode_ident) const
{
    auto it = m_faces.find(anode_ident);
    if (it == m_faces.end()) {
        raise<KeyError>("DriftPlanes: no faces for anode %d", anode_ident);
    }
    return it->second;
}

Configuration Gen::DriftPlanes::xregions() const { return m_xregions; }
