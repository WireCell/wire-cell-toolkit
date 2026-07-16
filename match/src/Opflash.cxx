#include "WireCellMatch/Opflash.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Match;

void Opflash::init(double time_, std::vector<double> pe, double threshold, int nchan,
                   const PEErr& pe_err, const std::vector<double>* pe_scale)
{
    flash_id    = 0;
    m_nchan     = nchan;
    m_threshold = threshold;
    type        = 0;
    low_time    = 0;
    high_time   = 0;
    time        = time_;

    PE = std::move(pe);
    PE.resize(nchan, 0);
    PE_err.resize(nchan, 1);

    // Optional per-channel measured-PE gain correction, applied before PE_err /
    // total / fired are derived so they all stay self-consistent with the
    // corrected measurement. nullptr (or short vector) => identity.
    if (pe_scale) {
        const int ns = static_cast<int>(pe_scale->size());
        for (int i = 0; i < nchan && i < ns; ++i) PE[i] *= (*pe_scale)[i];
    }

    total_PE = 0;
    for (int i = 0; i < nchan; ++i) {
        PE_err[i] = (PE[i] < pe_err.knee) ? pe_err.floor : pe_err.frac * PE[i];
        total_PE += PE[i];
        if (PE[i] > threshold) fired_channels.push_back(i);
    }
}

Opflash::Opflash(double time_, std::vector<double> pe, double threshold, int nchan,
                 const PEErr& pe_err, const std::vector<double>* pe_scale)
{
    init(time_, std::move(pe), threshold, nchan, pe_err, pe_scale);
}

Opflash::Opflash(const Clus::Facade::Flash& flash, double threshold, int nchan,
                 const PEErr& pe_err, const std::vector<double>* pe_scale)
{
    init(flash.time(), flash.pes(nchan), threshold, nchan, pe_err, pe_scale);
    flash_id = flash.ident();

    // Per-channel saturation flags, carried on the light-PC "error" field
    // (1 = the channel's hits overlap a DAPHNE rail in this flash; see
    // FlashTensorToOpticalPCs).  All-zero (and unread) unless the light
    // chain ran with flag_saturation and QLMatching enables
    // use_saturation_flag.
    const auto idents = flash.idents();
    const auto errors = flash.errors();
    for (size_t i = 0; i < idents.size() && i < errors.size(); ++i) {
        if (errors[i] > 0.5) {
            const int ch = idents[i];
            if (ch >= 0 && ch < nchan) {
                if (sat.empty()) sat.assign(nchan, 0);
                sat[ch] = 1;
            }
        }
    }

    // Per-channel readout-coverage fractions (sparse "flashcov" PC rows,
    // only channels with coverage < 1 of this flash's window; see
    // FlashTensorToOpticalPCs).  Empty (and get_cov == 1.0) unless the
    // light chain ran with emit_coverage; consumed by QLMatching
    // use_coverage_flag.
    const auto cov_idents = flash.cov_idents();
    const auto covs = flash.covs();
    for (size_t i = 0; i < cov_idents.size() && i < covs.size(); ++i) {
        const int ch = cov_idents[i];
        if (ch >= 0 && ch < nchan) {
            if (cov.empty()) cov.assign(nchan, 1.0f);
            cov[ch] = (float)covs[i];
        }
    }
}

void Opflash::inflate_nodata_err(double err, double cov_min)
{
    if (err <= 0 || cov.empty()) return;
    for (int i = 0; i < m_nchan; ++i) {
        if (get_cov(i) < cov_min) PE_err[i] = std::max(PE_err[i], err);
    }
}

Opflash::~Opflash() = default;

bool Opflash::get_fired(int ch) const
{
    return std::find(fired_channels.begin(), fired_channels.end(), ch) != fired_channels.end();
}
