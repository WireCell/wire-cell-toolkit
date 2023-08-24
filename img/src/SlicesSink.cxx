#include "WireCellImg/SlicesSink.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

WIRECELL_FACTORY(SlicesSink, WireCell::Img::SlicesSink, WireCell::ISliceFrameSink, WireCell::IConfigurable)

using namespace WireCell;
using namespace std;

Img::SlicesSink::SlicesSink() {}

Img::SlicesSink::~SlicesSink() {}

WireCell::Configuration Img::SlicesSink::default_configuration() const
{
    Configuration cfg;

    return cfg;
}

void Img::SlicesSink::configure(const WireCell::Configuration& cfg) { m_cfg = cfg; }

bool Img::SlicesSink::operator()(const ISliceFrame::pointer& sf)
{
    if (!sf) {
        return true;
    }

    // auto slices = sf->slices();

    // for (auto slice : slices) {
    //     auto cvmap = slice->activity();
    //     double qtot = 0;
    //     for (const auto& cv : cvmap) {
    //         qtot += cv.second;
    //     }
    // }
    return true;
}
