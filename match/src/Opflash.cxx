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
}

Opflash::~Opflash() = default;

bool Opflash::get_fired(int ch) const
{
    return std::find(fired_channels.begin(), fired_channels.end(), ch) != fired_channels.end();
}
