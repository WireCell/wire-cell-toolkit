#include "WireCellRoot/UbooneGeomHelper.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"

// ROOT
#include "TFile.h"
#include "TH3.h"

#include <string>

WIRECELL_FACTORY(UbooneGeomHelper, WireCell::Root::UbooneGeomHelper,
                 WireCell::INamed,
                 WireCell::IClusGeomHelper, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Root;
using WireCell::String::format;

UbooneGeomHelper::UbooneGeomHelper()
    : Aux::Logger("UbooneGeomHelper", "geom")
{
}

UbooneGeomHelper::~UbooneGeomHelper()
{
    delete m_h3_Dx;
    delete m_h3_Dy;
    delete m_h3_Dz;
}

WireCell::Configuration UbooneGeomHelper::default_configuration() const
{
    Configuration cfg;
    cfg["sce_offsets_file"] = "";  // path resolved via wc.resolve(); empty = no SCE correction
    return cfg;
}

void UbooneGeomHelper::configure(const WireCell::Configuration& cfg)
{
    // Store tpcparams (all keys except sce_offsets_file) for get_params / FV queries.
    m_tpcparams = cfg;
    m_tpcparams.removeMember("sce_offsets_file");

    for (const auto& apaface : m_tpcparams.getMemberNames()) {
        const Json::Value& tpcparam = m_tpcparams[apaface];
        FV fv;
        fv.m_FV_xmin        = get<double>(tpcparam, "FV_xmin",        fv.m_FV_xmin);
        fv.m_FV_xmax        = get<double>(tpcparam, "FV_xmax",        fv.m_FV_xmax);
        fv.m_FV_ymin        = get<double>(tpcparam, "FV_ymin",        fv.m_FV_ymin);
        fv.m_FV_ymax        = get<double>(tpcparam, "FV_ymax",        fv.m_FV_ymax);
        fv.m_FV_zmin        = get<double>(tpcparam, "FV_zmin",        fv.m_FV_zmin);
        fv.m_FV_zmax        = get<double>(tpcparam, "FV_zmax",        fv.m_FV_zmax);
        fv.m_FV_xmin_margin = get<double>(tpcparam, "FV_xmin_margin", fv.m_FV_xmin_margin);
        fv.m_FV_xmax_margin = get<double>(tpcparam, "FV_xmax_margin", fv.m_FV_xmax_margin);
        fv.m_FV_ymin_margin = get<double>(tpcparam, "FV_ymin_margin", fv.m_FV_ymin_margin);
        fv.m_FV_ymax_margin = get<double>(tpcparam, "FV_ymax_margin", fv.m_FV_ymax_margin);
        fv.m_FV_zmin_margin = get<double>(tpcparam, "FV_zmin_margin", fv.m_FV_zmin_margin);
        fv.m_FV_zmax_margin = get<double>(tpcparam, "FV_zmax_margin", fv.m_FV_zmax_margin);
        m_FV_map[apaface] = fv;
    }

    // Load SCE offset histograms if a file path is given.
    const std::string sce_raw = get<std::string>(cfg, "sce_offsets_file", "");
    if (!sce_raw.empty()) {
        std::string sce_path = Persist::resolve(sce_raw);
        if (sce_path.empty()) {
            raise<ValueError>("UbooneGeomHelper: sce_offsets_file not found: %s", sce_raw.c_str());
        }
        load_sce_offsets(sce_path);
        SPDLOG_LOGGER_INFO(log, "Loaded SCE offsets from {}", sce_path);
    }
}

// Load UBoone SCE position-correction histograms from a ROOT file.
// The file must contain TH3F histograms named "hDx", "hDy", "hDz"
// (e.g. SCEoffsets_dataDriven_combined_bkwd_Jan18.root).
// Content is copied into owned TH3D members so the file can be closed.
void UbooneGeomHelper::load_sce_offsets(const std::string& filename)
{
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        raise<ValueError>("UbooneGeomHelper: cannot open SCE offsets file: %s", filename.c_str());
    }

    auto copy_th3 = [&](const char* name) -> TH3D* {
        TH3F* h = dynamic_cast<TH3F*>(file->Get(name));
        if (!h) {
            raise<ValueError>("UbooneGeomHelper: histogram %s not found in %s",
                              name, filename.c_str());
        }
        int nx = h->GetNbinsX(), ny = h->GetNbinsY(), nz = h->GetNbinsZ();
        TH3D* out = new TH3D(
            (std::string("h3_") + name).c_str(), (std::string("h3_") + name).c_str(),
            nx, h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(nx),
            ny, h->GetYaxis()->GetBinLowEdge(1), h->GetYaxis()->GetBinUpEdge(ny),
            nz, h->GetZaxis()->GetBinLowEdge(1), h->GetZaxis()->GetBinUpEdge(nz));
        out->SetDirectory(nullptr);  // detach from any file/directory
        for (int i = 1; i <= nx; ++i)
            for (int j = 1; j <= ny; ++j)
                for (int k = 1; k <= nz; ++k)
                    out->SetBinContent(i, j, k, h->GetBinContent(i, j, k));
        return out;
    };

    delete m_h3_Dx; m_h3_Dx = copy_th3("hDx");
    delete m_h3_Dy; m_h3_Dy = copy_th3("hDy");
    delete m_h3_Dz; m_h3_Dz = copy_th3("hDz");

    file->Close();
    delete file;
}

WireCell::Configuration UbooneGeomHelper::get_params(const int apa, const int face) const
{
    const std::string apa_face = format("a%df%d", apa, face);
    if (!m_tpcparams.isMember(apa_face)) {
        raise<ValueError>("UbooneGeomHelper: no params for apa %d face %d", apa, face);
    }
    return m_tpcparams[apa_face];
}

bool UbooneGeomHelper::is_in_FV(const WireCell::Point& point, const int apa, const int face) const
{
    const std::string apa_face = format("a%df%d", apa, face);
    auto it = m_FV_map.find(apa_face);
    if (it == m_FV_map.end()) {
        raise<ValueError>("UbooneGeomHelper: no FV for apa %d face %d", apa, face);
    }
    const FV& fv = it->second;
    if (point.x() < fv.m_FV_xmin + fv.m_FV_xmin_margin || point.x() > fv.m_FV_xmax - fv.m_FV_xmax_margin) return false;
    if (point.y() < fv.m_FV_ymin + fv.m_FV_ymin_margin || point.y() > fv.m_FV_ymax - fv.m_FV_ymax_margin) return false;
    if (point.z() < fv.m_FV_zmin + fv.m_FV_zmin_margin || point.z() > fv.m_FV_zmax - fv.m_FV_zmax_margin) return false;
    return true;
}

bool UbooneGeomHelper::is_in_FV_dim(const WireCell::Point& point, const int dim, const double margin,
                                     const int apa, const int face) const
{
    const std::string apa_face = format("a%df%d", apa, face);
    auto it = m_FV_map.find(apa_face);
    if (it == m_FV_map.end()) {
        raise<ValueError>("UbooneGeomHelper: no FV for apa %d face %d", apa, face);
    }
    const FV& fv = it->second;
    if (dim == 0) return point.x() >= fv.m_FV_xmin + margin && point.x() <= fv.m_FV_xmax - margin;
    if (dim == 1) return point.y() >= fv.m_FV_ymin + margin && point.y() <= fv.m_FV_ymax - margin;
    if (dim == 2) return point.z() >= fv.m_FV_zmin + margin && point.z() <= fv.m_FV_zmax - margin;
    return true;
}

// Apply UBoone data-driven SCE position correction.
//
// Ported from WCP::TPCParams::func_pos_SCE_correction() in
// prototype_base/data/src/TPCParams.cxx.
//
// UBoone-specific constants (do not generalise to other detectors):
//   active volume: x=256 cm, y=232.5 cm, z=1037 cm
//   histogram axis scaling: scale_x=2.50/2.56, scale_y=2.50/2.33, scale_z=10.0/10.37
WireCell::Point UbooneGeomHelper::get_corrected_point(
    const WireCell::Point& point,
    const WireCell::IClusGeomHelper::CorrectionType type,
    const int apa, const int face) const
{
    if (type == IClusGeomHelper::NONE) {
        return point;
    }
    if (type != IClusGeomHelper::SCE) {
        raise<ValueError>("UbooneGeomHelper: unsupported CorrectionType %d", static_cast<int>(type));
    }
    if (!m_h3_Dx || !m_h3_Dy || !m_h3_Dz) {
        // SCE histograms not loaded (no sce_offsets_file configured) — return unchanged.
        SPDLOG_LOGGER_WARN(log, "get_corrected_point(SCE) called but no SCE offsets loaded; returning raw point");
        return point;
    }

    const double y_length = 232.5;  // cm
    const double scale_x  = 2.50 / 2.56;
    const double scale_y  = 2.50 / 2.33;
    const double scale_z  = 10.0 / 10.37;

    double px = point.x() / units::cm;
    double py = point.y() / units::cm;
    double pz = point.z() / units::cm;

    // Shift y origin to bottom of active volume
    py += y_length / 2.0;

    // Clamp to valid histogram range
    if (px < 0.001)    px = 0.001;
    if (px > 255.999)  px = 255.998;
    if (py < 0.001)    py = 0.001;
    if (py > 232.499)  py = 232.498;
    if (pz < 0.001)    pz = 0.001;
    if (pz > 1036.999) pz = 1036.998;

    // Map to histogram coordinate system (scaled metres)
    double cw_x = 2.50 - scale_x * (px / 100.0);
    double cw_y = scale_y * (py / 100.0);
    double cw_z = scale_z * (pz / 100.0);

    // Interpolate displacement corrections (cm)
    double corr_x = m_h3_Dx->Interpolate(cw_x, cw_y, cw_z);
    double corr_y = m_h3_Dy->Interpolate(cw_x, cw_y, cw_z);
    double corr_z = m_h3_Dz->Interpolate(cw_x, cw_y, cw_z);

    // Apply corrections, undo y shift
    double out_x = (px - corr_x / scale_x) * units::cm;
    double out_y = (py + corr_y / scale_y - y_length / 2.0) * units::cm;
    double out_z = (pz + corr_z / scale_z) * units::cm;

    return WireCell::Point(out_x, out_y, out_z);
}
