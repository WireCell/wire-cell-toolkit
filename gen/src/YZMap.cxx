#include "WireCellGen/YZMap.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Exceptions.h"

#include <algorithm>
#include <cmath>

WIRECELL_FACTORY(YZMap, WireCell::Gen::YZMap,
                 WireCell::INamed,
                 WireCell::IYZMap, WireCell::IConfigurable)

using namespace WireCell;

Gen::YZMap::YZMap()
    : Aux::Logger("YZMap", "gen")
{
}

Gen::YZMap::~YZMap() {}

WireCell::Configuration Gen::YZMap::default_configuration() const
{
    Configuration cfg;
    cfg["filename"]   = "";
    cfg["bin_width"]  = m_bin_width;
    cfg["bin_height"] = m_bin_height;
    cfg["yoffset"]    = m_yoffset;
    cfg["zoffset"]    = m_zoffset;
    cfg["nbinsy"]     = m_nbinsy;
    cfg["nbinsz"]     = m_nbinsz;
    return cfg;
}

void Gen::YZMap::configure(const WireCell::Configuration& cfg)
{
    const std::string filename = cfg["filename"].asString();
    if (filename.empty()) {
        THROW(ValueError() << errmsg{"YZMap: filename is required"});
    }

    m_bin_width  = get<double>(cfg, "bin_width",  m_bin_width);
    m_bin_height = get<double>(cfg, "bin_height", m_bin_height);
    m_yoffset    = get<double>(cfg, "yoffset",    m_yoffset);
    m_zoffset    = get<double>(cfg, "zoffset",    m_zoffset);
    m_nbinsy     = get<int>   (cfg, "nbinsy",     m_nbinsy);
    m_nbinsz     = get<int>   (cfg, "nbinsz",     m_nbinsz);

    auto jmap = Persist::load(filename);

    m_map.clear();
    for (const auto& anode_name : jmap.getMemberNames()) {
        auto& plane_map = m_map[anode_name];
        const auto& jplanes = jmap[anode_name];
        for (const auto& plane_str : jplanes.getMemberNames()) {
            int plane = std::stoi(plane_str);
            auto& grid = plane_map[plane];
            const auto& jgrid = jplanes[plane_str];
            const int nbinsz = jgrid.size();
            grid.resize(nbinsz);
            for (int binz = 0; binz < nbinsz; ++binz) {
                const auto& jrow = jgrid[binz];
                const int nbinsy = jrow.size();
                grid[binz].resize(nbinsy);
                for (int biny = 0; biny < nbinsy; ++biny) {
                    grid[binz][biny] = jrow[biny].asDouble();
                }
            }
        }
    }

    log->debug("loaded {} anodes from {}", m_map.size(), filename);
}

double Gen::YZMap::value(const std::string& anode_name,
                         int plane,
                         double y, double z) const
{
    auto ait = m_map.find(anode_name);
    if (ait == m_map.end()) {
        return 0.0;
    }
    auto pit = ait->second.find(plane);
    if (pit == ait->second.end()) {
        return 0.0;
    }
    const auto& grid = pit->second;

    int biny = (int)std::floor((y + m_yoffset) / m_bin_height);
    int binz = (int)std::floor((z + m_zoffset) / m_bin_width);

    biny = std::clamp(biny, 0, m_nbinsy);
    binz = std::clamp(binz, 0, m_nbinsz);

    return grid[binz][biny];
}
