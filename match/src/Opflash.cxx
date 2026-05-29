#include "WireCellMatch/Opflash.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Match;

void Opflash::init(double time_, std::vector<double> pe, double threshold, int nchan)
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

    total_PE = 0;
    for (int i = 0; i < nchan; ++i) {
        PE_err[i] = (PE[i] < 1) ? 0.3 : 0.3 * PE[i];
        total_PE += PE[i];
        if (PE[i] > threshold) fired_channels.push_back(i);
    }
}

Opflash::Opflash(double time_, std::vector<double> pe, double threshold, int nchan)
{
    init(time_, std::move(pe), threshold, nchan);
}

Opflash::Opflash(const Clus::Facade::Flash& flash, double threshold, int nchan)
{
    init(flash.time(), flash.pes(nchan), threshold, nchan);
    flash_id = flash.ident();
}

Opflash::~Opflash() = default;

bool Opflash::get_fired(int ch) const
{
    return std::find(fired_channels.begin(), fired_channels.end(), ch) != fired_channels.end();
}
