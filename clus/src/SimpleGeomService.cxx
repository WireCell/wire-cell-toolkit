#include "WireCellClus/SimpleGeomService.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

#include <string>

WIRECELL_FACTORY(SimpleGeomService, WireCell::Clus::SimpleGeomService,
                 WireCell::INamed,
                 WireCell::IGeomService, WireCell::IConfigurable)

using namespace WireCell;
using namespace std;
using WireCell::String::format;

Clus::SimpleGeomService::SimpleGeomService()
  : Aux::Logger("SimpleGeomService", "geom")
{
}

WireCell::Configuration Clus::SimpleGeomService::default_configuration() const
{
    Configuration cfg;
    return cfg;
}

void Clus::SimpleGeomService::configure(const WireCell::Configuration& cfg)
{
    m_tpcparams = cfg;

    for (const auto& apaface : m_tpcparams.getMemberNames()) {
        const Json::Value& tpcparam = m_tpcparams[apaface];
        FV fv;
        fv.m_FV_xmin = get<double>(tpcparam, "FV_xmin", fv.m_FV_xmin);
        fv.m_FV_xmax = get<double>(tpcparam, "FV_xmax", fv.m_FV_xmax);
        fv.m_FV_ymin = get<double>(tpcparam, "FV_ymin", fv.m_FV_ymin);
        fv.m_FV_ymax = get<double>(tpcparam, "FV_ymax", fv.m_FV_ymax);
        fv.m_FV_zmin = get<double>(tpcparam, "FV_zmin", fv.m_FV_zmin);
        fv.m_FV_zmax = get<double>(tpcparam, "FV_zmax", fv.m_FV_zmax);
        fv.m_FV_xmin_margin = get<double>(tpcparam, "FV_xmin_margin", fv.m_FV_xmin_margin);
        fv.m_FV_xmax_margin = get<double>(tpcparam, "FV_xmax_margin", fv.m_FV_xmax_margin);
        fv.m_FV_ymin_margin = get<double>(tpcparam, "FV_ymin_margin", fv.m_FV_ymin_margin);
        fv.m_FV_ymax_margin = get<double>(tpcparam, "FV_ymax_margin", fv.m_FV_ymax_margin);
        fv.m_FV_zmin_margin = get<double>(tpcparam, "FV_zmin_margin", fv.m_FV_zmin_margin);
        fv.m_FV_zmax_margin = get<double>(tpcparam, "FV_zmax_margin", fv.m_FV_zmax_margin);
    }
    
}

WireCell::Configuration Clus::SimpleGeomService::get_params(const int apa, const int face) const
{
    const std::string apa_face = format("a%df%d", apa, face);
    WireCell::Configuration cfg = m_tpcparams[apa_face];
    return cfg;
}

bool Clus::SimpleGeomService::is_in_FV(const WireCell::Point& point, const int apa, const int face) const
{
    const std::string apa_face = format("a%df%d", apa, face);
    if (m_FV_map.find(apa_face) == m_FV_map.end()) {
        raise<ValueError>("failed to find face for wpid %d", face);
    }
    const FV& fv = m_FV_map.at(apa_face);
    if (point.x() < fv.m_FV_xmin + fv.m_FV_xmin_margin || point.x() > fv.m_FV_xmax - fv.m_FV_xmax_margin) {
        return false;
    }
    if (point.y() < fv.m_FV_ymin + fv.m_FV_ymin_margin || point.y() > fv.m_FV_ymax - fv.m_FV_ymax_margin) {
        return false;
    }
    if (point.z() < fv.m_FV_zmin + fv.m_FV_zmin_margin || point.z() > fv.m_FV_zmax - fv.m_FV_zmax_margin) {
        return false;
    }
    return true;
}

bool Clus::SimpleGeomService::is_in_FV_dim(const WireCell::Point& point, const int dim, const double margin, const int apa, const int face) const
{
    const std::string apa_face = format("a%df%d", apa, face);
    if (m_FV_map.find(apa_face) == m_FV_map.end()) {
        raise<ValueError>("failed to find face for wpid %d", face);
    }
    const FV& fv = m_FV_map.at(apa_face);
    if (dim == 0) {
        if (point.x() < fv.m_FV_xmin + margin || point.x() > fv.m_FV_xmax - margin) {
            return false;
        }
    }
    if (dim == 1) {
        if (point.y() < fv.m_FV_ymin + margin || point.y() > fv.m_FV_ymax - margin) {
            return false;
        }
    }
    if (dim == 2) {
        if (point.z() < fv.m_FV_zmin + margin || point.z() > fv.m_FV_zmax - margin) {
            return false;
        }
    }
    return true;
}

WireCell::Point Clus::SimpleGeomService::get_corrected_point(const WireCell::Point& point, const WireCell::IGeomService::CorrectionType type, const int apa, const int face) const
{
    if (type != WireCell::IGeomService::CorrectionType::NONE) {
        raise<ValueError>("failed to find face for wpid %d", type);
    }
    return point;
}