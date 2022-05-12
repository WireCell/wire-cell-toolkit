#include "WireCellSio/NumpyFrameSaver.h"

#include "WireCellAux/FrameTools.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/NumpyHelper.h"

#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <sstream>

WIRECELL_FACTORY(NumpyFrameSaver, WireCell::Sio::NumpyFrameSaver, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Numpy::save2d;

Sio::NumpyFrameSaver::NumpyFrameSaver()
    : Aux::Logger("NumpyFrameSaver", "sio")
{
}

Sio::NumpyFrameSaver::~NumpyFrameSaver() {}

// WireCell::Configuration Sio::NumpyFrameSaver::default_configuration() const
// {
//     nljs_t nljs = m_cfg;
//     return nljs.get<Json::Value>();
// }

// void Sio::NumpyFrameSaver::configure(const WireCell::Configuration& config)
// {
//     nljs_t nljs = config;
//     m_cfg = nljs.get<config_t>();
//     if (m_cfg.frame_tags.empty()) {
//         // internally we use an empty tag to mean all tags
//         m_cfg.frame_tags.push_back("");
//     }
// }
void Sio::NumpyFrameSaver::configured()
{
    if (m_cfg.frame_tags.empty()) {
        // internally we use an empty tag to mean all tags
        m_cfg.frame_tags.push_back("");
    }
}

bool Sio::NumpyFrameSaver::operator()(const IFrame::pointer& inframe, IFrame::pointer& outframe)
{
    if (!inframe) {
        log->debug("EOS at call={}", m_save_count);
        outframe = nullptr;
        return true;
    }

    outframe = inframe;  // pass through actual frame

    const std::string mode = "a";

    // Eigen3 array is indexed as (irow, icol) or (ichan, itick)
    // one row is one channel, one column is a tick.
    // Numpy saves reversed dimensions: {ncols, nrows} aka {ntick, nchan} dimensions.

    std::stringstream ss;
    ss << "NumpyFrameSaver: see frame #" << inframe->ident() << " with " << inframe->traces()->size()
       << " traces with frame tags:";
    for (auto t : inframe->frame_tags()) {
        ss << " \"" << t << "\"";
    }
    ss << " and trace tags:";
    for (auto t : inframe->trace_tags()) {
        ss << " \"" << t << "\"";
    }
    ss << " looking for tags:";
    for (auto ft : m_cfg.frame_tags) {
        if (ft.empty()) {
            ss << " <all>";
        }
        else {
            ss << " \"" << ft << "\"";
        }
    }
    log->debug(ss.str());

    for (auto tag : m_cfg.frame_tags) {
        auto traces = Aux::tagged_traces(inframe, tag);
        log->debug("NumpyFrameSaver: save {} tagged as {}", traces.size(), tag);
        if (traces.empty()) {
            log->warn("NumpyFrameSaver: no traces for tag: \"{}\"", tag);
            continue;
        }
        auto channels = Aux::channels(traces);
        std::sort(channels.begin(), channels.end());
        auto chbeg = channels.begin();
        auto chend = std::unique(chbeg, channels.end());
        auto tbinmm = Aux::tbin_range(traces);

        // fixme: may want to give user some config over tbin range to save.
        const size_t ncols = tbinmm.second - tbinmm.first;
        const size_t nrows = std::distance(chbeg, chend);
        log->debug("NumpyFrameSaver: saving ncols={} nrows={}", ncols, nrows);

        Array::array_xxf arr = Array::array_xxf::Zero(nrows, ncols) + m_cfg.baseline;
        Aux::fill(arr, traces, channels.begin(), chend, tbinmm.first);
        arr = arr * m_cfg.scale + m_cfg.offset;

        {  // the 2D frame array
            const std::string aname = String::format("frame_%s_%d", tag.c_str(), m_save_count);
            if (m_cfg.digitize) {
                Array::array_xxs sarr = arr.cast<short>();
                save2d(sarr, aname, m_cfg.filename, mode);
            }
            else {
                save2d(arr, aname, m_cfg.filename, mode);
            }
            log->debug("NumpyFrameSaver: saved {} with {} channels {} ticks @t={} ms qtot={}", aname, nrows, ncols,
                     inframe->time() / units::ms, arr.sum());
        }

        {  // the channel array
            const std::string aname = String::format("channels_%s_%d", tag.c_str(), m_save_count);
            cnpy::npz_save(m_cfg.filename, aname, channels.data(), {nrows}, mode);
            
        }

        {  // the tick array
            const std::string aname = String::format("tickinfo_%s_%d", tag.c_str(), m_save_count);
            const std::vector<double> tickinfo{inframe->time(), inframe->tick(), (double) tbinmm.first};
            cnpy::npz_save(m_cfg.filename, aname, tickinfo.data(), {3}, mode);
        }
    }

    ++m_save_count;
    return true;
}
