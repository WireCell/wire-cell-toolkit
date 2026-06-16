#include "WireCellRoot/SCEFieldTH3.h"
#include "WireCellUtil/NamedFactory.h"

#include <TFile.h>
#include <TH3F.h>
#include <TAxis.h>

#include <stdexcept>

WIRECELL_FACTORY(SCEFieldTH3, WireCell::Root::SCEFieldTH3,
                 WireCell::ISCEField, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Root;

SCEFieldTH3::SCEFieldTH3()
    : m_th3_name_E("TrueBkwd_Displacement_X_E"),
      m_th3_name_W("TrueBkwd_Displacement_X_W"),
      m_th3_name_E_y("TrueBkwd_Displacement_Y_E"),
      m_th3_name_W_y("TrueBkwd_Displacement_Y_W"),
      m_th3_name_E_z("TrueBkwd_Displacement_Z_E"),
      m_th3_name_W_z("TrueBkwd_Displacement_Z_W"),
      l(Log::logger("wct"))
{
}

SCEFieldTH3::~SCEFieldTH3() = default;

Configuration SCEFieldTH3::default_configuration() const
{
    Configuration cfg;
    cfg["sce_map_file"]  = m_sce_map_file;
    cfg["th3_name_E"]    = m_th3_name_E;
    cfg["th3_name_W"]    = m_th3_name_W;
    cfg["th3_name_E_y"]  = m_th3_name_E_y;
    cfg["th3_name_W_y"]  = m_th3_name_W_y;
    cfg["th3_name_E_z"]  = m_th3_name_E_z;
    cfg["th3_name_W_z"]  = m_th3_name_W_z;
    cfg["axis_unit_mm"]  = m_axis_unit_mm;
    cfg["sign"]          = m_sign;
    return cfg;
}

void SCEFieldTH3::load_pair(const std::string& nameE, const std::string& nameW,
                            TH3F*& hE, TH3F*& hW, bool hard)
{
    hE = dynamic_cast<TH3F*>(m_file->Get(nameE.c_str()));
    hW = dynamic_cast<TH3F*>(m_file->Get(nameW.c_str()));
    if (!hE || !hW) {
        if (hard) {
            throw std::runtime_error("SCEFieldTH3: missing TH3F '" + nameE
                                     + "' or '" + nameW + "' in " + m_sce_map_file);
        }
        l->warn("SCEFieldTH3: optional TH3 '{}' / '{}' absent -- that component is a no-op",
                nameE, nameW);
        hE = nullptr; hW = nullptr;
        return;
    }
    hE->SetDirectory(nullptr);
    hW->SetDirectory(nullptr);
}

void SCEFieldTH3::configure(const Configuration& cfg)
{
    m_sce_map_file = get<std::string>(cfg, "sce_map_file", m_sce_map_file);
    m_th3_name_E   = get<std::string>(cfg, "th3_name_E",   m_th3_name_E);
    m_th3_name_W   = get<std::string>(cfg, "th3_name_W",   m_th3_name_W);
    m_th3_name_E_y = get<std::string>(cfg, "th3_name_E_y", m_th3_name_E_y);
    m_th3_name_W_y = get<std::string>(cfg, "th3_name_W_y", m_th3_name_W_y);
    m_th3_name_E_z = get<std::string>(cfg, "th3_name_E_z", m_th3_name_E_z);
    m_th3_name_W_z = get<std::string>(cfg, "th3_name_W_z", m_th3_name_W_z);
    m_axis_unit_mm = get<bool>       (cfg, "axis_unit_mm", m_axis_unit_mm);
    m_sign         = get<double>     (cfg, "sign",         m_sign);

    if (m_sce_map_file.empty()) {
        l->info("SCEFieldTH3: no sce_map_file -- field is a no-op");
        return;
    }

    l->info("SCEFieldTH3: opening SCE map {}", m_sce_map_file);
    m_file.reset(TFile::Open(m_sce_map_file.c_str(), "READ"));
    if (!m_file || m_file->IsZombie()) {
        throw std::runtime_error("SCEFieldTH3: cannot open " + m_sce_map_file);
    }

    load_pair(m_th3_name_E,   m_th3_name_W,   m_hE,   m_hW,   /*hard=*/true);
    load_pair(m_th3_name_E_y, m_th3_name_W_y, m_hE_y, m_hW_y, /*hard=*/false);
    load_pair(m_th3_name_E_z, m_th3_name_W_z, m_hE_z, m_hW_z, /*hard=*/false);

    l->info("SCEFieldTH3: loaded X[E='{}' W='{}'] Y[{}] Z[{}] axis_mm={} sign={}",
            m_th3_name_E, m_th3_name_W,
            (m_hE_y ? "ok" : "absent"), (m_hE_z ? "ok" : "absent"),
            m_axis_unit_mm, m_sign);
}

double SCEFieldTH3::clamp_axis(const TAxis* a, double v)
{
    const double lo  = a->GetBinCenter(1);
    const double hi  = a->GetBinCenter(a->GetNbins());
    const double pad = 1e-3 * a->GetBinWidth(1);
    if (v < lo + pad) return lo + pad;
    if (v > hi - pad) return hi - pad;
    return v;
}

double SCEFieldTH3::interp(TH3F* hE, TH3F* hW, int apa,
                           double x, double y, double z) const
{
    TH3F* h = (apa == 0) ? hE : (apa == 1) ? hW : nullptr;
    if (!h) return 0.0;

    // WCT input is mm; map axes are cm by default (mm if m_axis_unit_mm).
    const double s_in  = m_axis_unit_mm ? 1.0 : 0.1;   // mm -> axis unit
    const double s_out = m_axis_unit_mm ? 1.0 : 10.0;  // map value -> mm

    const double xq = clamp_axis(h->GetXaxis(), x * s_in);
    const double yq = clamp_axis(h->GetYaxis(), y * s_in);
    const double zq = clamp_axis(h->GetZaxis(), z * s_in);
    const double d_map = h->Interpolate(xq, yq, zq);

    return m_sign * d_map * s_out;  // signed displacement in mm
}

double SCEFieldTH3::displacement_x(int apa, double x, double y, double z) const
{
    return interp(m_hE, m_hW, apa, x, y, z);
}

double SCEFieldTH3::displacement_y(int apa, double x, double y, double z) const
{
    return interp(m_hE_y, m_hW_y, apa, x, y, z);
}

double SCEFieldTH3::displacement_z(int apa, double x, double y, double z) const
{
    return interp(m_hE_z, m_hW_z, apa, x, y, z);
}
