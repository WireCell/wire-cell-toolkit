#include "WireCellSio/FrameFileSink.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Waveform.h"

#include "WireCellAux/FrameTools.h"


WIRECELL_FACTORY(FrameFileSink, WireCell::Sio::FrameFileSink,
                 WireCell::INamed,
                 WireCell::IFrameSink, WireCell::ITerminal,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Stream;
using namespace WireCell::Waveform;

Sio::FrameFileSink::FrameFileSink()
    : Aux::Logger("FrameFileSink", "io")
{
}

Sio::FrameFileSink::~FrameFileSink()
{
}

void Sio::FrameFileSink::finalize()
{
    log->debug("closing {} after {} calls",
               m_outname, m_count);
    m_out.pop();
}

WireCell::Configuration Sio::FrameFileSink::default_configuration() const
{
    Configuration cfg;
    // Output tar file name.
    cfg["outname"] = m_outname;

    // Select which traces to consider.
    cfg["tags"] = Json::arrayValue;

    // If digitize is true, then samples as 16 bit ints.  Otherwise
    // save as 32 bit floats.
    cfg["digitize"] = false;

    // This number is set to the waveform sample array before any
    // charge is added.
    cfg["baseline"] = 0.0;

    // This number will be multiplied to each waveform sample before
    // casting to dtype.
    cfg["scale"] = 1.0;

    // This number will be added to each scaled waveform sample before
    // casting to dtype.
    cfg["offset"] = 0.0;

    cfg["masks"] = true;

    // If the "dense" option is given the frame array will be extended
    // and padded prior to output.  The value of "dense" should be an
    // object with keys "chbeg" and "chend" which give half-inclusive
    // range of channel IDs (chend is +1 past the last) and "tbbeg"
    // and "tbend" which gives likewise is half-inclusive of "tbin"
    // range.  Like with non-dense saving, the "offset" gives the
    // padding value.

    return cfg;
}

void Sio::FrameFileSink::configure(const WireCell::Configuration& cfg)
{
    m_outname = get(cfg, "outname", m_outname);

    m_masks = get(cfg, "masks", m_masks);

    m_out.clear();
    output_filters(m_out, m_outname);
    if (m_out.size() < 1) {
        THROW(ValueError() << errmsg{"FrameFileSink: unsupported outname: " + m_outname});
    }

    m_tags.clear();
    for (auto jtag : cfg["tags"]) {
        m_tags.push_back(jtag.asString());
    }
    if (m_tags.empty()) {
        m_tags.push_back("*");
    }
    std::string comma="";
    std::string stags="";
    for (const auto& tag : m_tags) {
        stags += comma + tag;
        comma = ",";
    }

    m_digitize = get<bool>(cfg, "digitize", m_digitize);
    m_baseline = get(cfg, "baseline", m_baseline);
    m_scale = get(cfg, "scale", m_scale);
    m_offset = get(cfg, "offset", m_offset);

    if (cfg["dense"].isNull() or cfg["dense"].empty()) {
        log->debug("save {} with baseline={} scale={} offset={} digitize={} to {}",
                   stags, m_baseline, m_scale, m_offset, m_digitize, m_outname);
        return;
    }
    m_dense = true;
    auto dense = cfg["dense"];
    m_chbeg = dense["chbeg"].asInt();
    m_chend = dense["chend"].asInt();
    m_tbbeg = dense["tbbeg"].asInt();
    m_tbend = dense["tbend"].asInt();
    log->debug("save {} with baseline={} scale={} offset={} digitize={} ch=[{},{}] tbin=[{},{}] to {}",
               stags, m_baseline, m_scale, m_offset, m_digitize,
               m_chbeg, m_chend, m_tbbeg, m_tbend,
               m_outname);
}

void Sio::FrameFileSink::one_tag(const IFrame::pointer& frame,
                                 const std::string& tag)
{
    ITrace::vector traces;
    if (tag == "*") {
        // The designer of IFrame was a dumb dumb in chosing to return
        // a shared vector.  He is also to blame for writing this
        // comment.
        auto sv = frame->traces();
        traces.insert(traces.begin(), sv->begin(), sv->end());
    }
    else {
        traces = Aux::tagged_traces(frame, tag);
    }
    if (traces.empty()) {
        log->warn("call={} frame={} ntraces={} tag=\"{}\" zero traces",
                   m_count, frame->ident(),traces.size(), tag);
        return;
    }

    log->debug("call={} frame={} ntraces={} tag=\"{}\"",
               m_count, frame->ident(),traces.size(), tag);

    Aux::channel_list channels;
    std::pair<int, int> tbinmm;
    if (m_dense) {
        channels.resize(m_chend-m_chbeg, 0);
        std::iota(channels.begin(), channels.end(), m_chbeg);
        tbinmm = std::make_pair(m_tbbeg, m_tbend);
    }
    else {
        auto tmp = Aux::channels(traces);
        std::sort(tmp.begin(), tmp.end());
        auto chbeg = tmp.begin();
        auto chend = std::unique(chbeg, tmp.end());
        channels.insert(channels.begin(), chbeg, chend);
        tbinmm = Aux::tbin_range(traces);
    }
    

    const size_t ncols = tbinmm.second - tbinmm.first;
    const size_t nrows = std::distance(channels.begin(), channels.end());

    Array::array_xxf arr = Array::array_xxf::Zero(nrows, ncols) + m_baseline;
    Aux::fill(arr, traces, channels.begin(), channels.end(), tbinmm.first);
    arr = arr * m_scale + m_offset;

    {  // the 2D frame array
        const std::string aname = String::format("frame_%s_%d.npy", tag.c_str(), frame->ident());
        if (m_digitize) {
            Array::array_xxs sarr = arr.cast<short>();
            write(m_out, aname, sarr);
        }
        else {
            write(m_out, aname, arr);
        }
    }

    {  // the channel array
        const std::string aname = String::format("channels_%s_%d.npy", tag.c_str(), frame->ident());
        write(m_out, aname, channels);
    }

    {  // the tick array
        const std::string aname = String::format("tickinfo_%s_%d.npy", tag.c_str(), frame->ident());
        const std::vector<double> tickinfo{frame->time(), frame->tick(), (double) tbinmm.first};
        write(m_out, aname, tickinfo);
    }
    m_out.flush();

}

static
size_t cms_size(const ChannelMasks& cms)
{
    size_t size = 0;
    for (const auto& [chid, brl] : cms) {
        size += brl.size();
    }
    return size;
}

void Sio::FrameFileSink::masks(const IFrame::pointer& frame)
{
    // - cmm :: string -> ChannelMasks
    // - ChannelMasks :: int -> BinRangeList
    // - BinRangeList :: vector<BinRange>
    // - BinRange :: pair<int,int>

    // Flatten each ChannelMasks to array shape (N,3) with colums:
    // (chid, begin_time_bin, end_time_bin).  chid values will be
    // duplicated for each BinRange in a BinRangeList.
    constexpr size_t ncols = 3;

    auto cmm = frame->masks();
    if (cmm.empty()) {
        log->debug("no channel mask maps at call {}", m_count);
        return;
    }

    for (const auto& [name, cms] : cmm) {

        const size_t nrows = cms_size(cms);
        Eigen::Array<int, Eigen::Dynamic, ncols> arr(nrows, ncols);
        size_t irow = 0;
        for (const auto& [chid, brl] : cms) {
            for (const auto& [tbeg, tend] : brl) {
                arr(irow, 0) = chid;
                arr(irow, 1) = tbeg;
                arr(irow, 2) = tend;
                ++irow;
            }
        }

        const std::string aname = String::format("chanmask_%s_%d.npy", name.c_str(), frame->ident());
        write(m_out, aname, arr);

        log->debug("save {} with {} entries", aname, nrows);
    }
}


bool Sio::FrameFileSink::operator()(const IFrame::pointer& frame)
{
    if (! frame) { // eos
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    for (auto tag : m_tags) {
        one_tag(frame, tag);
    }
    if (m_masks) {
        masks(frame);
    }
    ++m_count;
    return true;
}
