#include "WireCellSio/FrameFileSource.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellAux/FrameTools.h"

WIRECELL_FACTORY(FrameFileSource, WireCell::Sio::FrameFileSource,
                 WireCell::INamed,
                 WireCell::IFrameSource, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Sio;
using namespace WireCell::Stream;
using namespace WireCell::String;
using namespace WireCell::Waveform;

FrameFileSource::FrameFileSource()
    : Aux::Logger("FrameFileSource", "io")
{
}

FrameFileSource::~FrameFileSource()
{
}

WireCell::Configuration FrameFileSource::default_configuration() const
{
    Configuration cfg;
    // Input tar stream
    cfg["inname"] = m_inname;

    // Set of traces to consider.
    cfg["tags"] = Json::arrayValue;

    cfg["frame_tags"] = Json::arrayValue;

    return cfg;
}

void FrameFileSource::configure(const WireCell::Configuration& cfg)
{
    m_inname = get(cfg, "inname", m_inname);

    m_in.clear();
    input_filters(m_in, m_inname);
    if (m_in.size() < 1) {
        THROW(ValueError() << errmsg{"FrameFileSource: unsupported inname: " + m_inname});
    }

    m_tags.clear();
    for (auto jtag : cfg["tags"]) {
        m_tags.push_back(jtag.asString());
    }
    if (m_tags.empty()) {
        m_tags.push_back("*");
    }

    m_frame_tags.clear();
    for (auto jtag : cfg["frame_tags"]) {
        m_frame_tags.push_back(jtag.asString());
    }
}

bool FrameFileSource::matches(const std::string& tag)
{
    for (const auto& maybe : m_tags) {
        if (maybe == "*") {
            return true;
        }
        if (maybe == tag) {
            return true;
        }
        log->debug("read unmatched tag: {} != ", tag, maybe);
    }
    return false;    
}

bool FrameFileSource::read()
{
    clear();

    custard::read(m_in, m_cur.fname, m_cur.fsize);
    m_cur.pig.read(m_in);

    // log->debug("read file \"{}\" from stream with {} bytes", m_cur.fname, m_cur.fsize);

    if (m_cur.fname.empty()) {
        // log->warn("read empty file name from stream");
        return false;
    }
    auto parts = split(m_cur.fname, "_");
    if (parts.size() != 3) {
        // log->warn("read parse file name failed got {} parts from {}", parts.size(), m_cur.fname);
        return false;
    }
    auto rparts = split(parts[2], ".");
    m_cur.ident = std::atoi(rparts[0].c_str());
    m_cur.okay = true;
    m_cur.type = parts[0];
    m_cur.tag = parts[1];
    m_cur.ext = rparts[1];
    return true;
}

IFrame::pointer FrameFileSource::load()
{
    int ident = -1;

    // Collect portions of the frame data;
    using trace_array_t = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    struct framelet_t { 
        std::vector<double> tickinfo;
        std::vector<int> channels;
        trace_array_t trace_array;
        std::string tag{""};
    };
    std::vector<framelet_t> framelets; // ordered
    std::map<std::string, size_t> framelet_indices;
    auto get_framelet = [&](const std::string& tag) -> framelet_t& {
        auto it = framelet_indices.find(tag);
        if (it == framelet_indices.end()) {
            const size_t oldsize = framelets.size();
            framelet_indices[tag] = oldsize;
            framelets.resize(oldsize+1);
            return framelets.back();
        }
        return framelets[it->second];
    };

    // fill this along the way
    ChannelMaskMap cmm;

    while (true) {

        // fixme: throw on errors but return nullptr on EOF.

        // Read next element if current cursor is empty
        if (m_cur.fsize == 0) {
            this->read();

            if (m_cur.fsize == 0) {
                // zero read means EOF
                break;
            }

            if (!m_in) {
                log->error("call={}, bad pig read with file={}", m_count, m_inname);
                return nullptr;
            }
            if (!m_cur.okay) {
                log->error("call={}, failed to parse npy file name {} in file={}",
                           m_count, m_cur.fname, m_inname);
                return nullptr;
            }
        }

        if (ident < 0) {
            // first time
            ident = m_cur.ident;
        }
        if (ident != m_cur.ident) {
            // we see next frame, done with this current
            break;       
        }

        if (m_cur.type == "chanmask") {
            Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> cmsarr;
            bool ok = pigenc::eigen::load(m_cur.pig, cmsarr);
            if (!ok) {
                log->error("call={}, cmm load failed tag=\"{}\" file={}",
                           m_count, m_cur.tag, m_inname);
                return nullptr;
            }
            const int ncols = cmsarr.cols();
            if (ncols != 3) {
                log->error("call={}, cmm wrong size nrows={}!=3 tag=\"{}\" file={}",
                           m_count, ncols, m_cur.tag, m_inname);
                return nullptr;
            }
            ChannelMasks cms;
            int nrows = cmsarr.rows();
            for (int irow=0; irow<nrows; ++irow) {
                auto row = cmsarr.row(irow);
                cms[row[0]].push_back(std::make_pair(row[1], row[2]));
            }
            cmm[m_cur.tag] = cms;

            clear();
            continue;
        }
    

        if (! matches(m_cur.tag)) {
            log->debug("call={}, skipping unmatched tag=\"{}\" file={}",
                       m_count, m_cur.tag, m_inname);
            clear();
            continue;
        }

        if (m_cur.type == "frame") {
            auto& framelet = get_framelet(m_cur.tag);
            framelet.tag = m_cur.tag;
            auto dtype = m_cur.pig.header().dtype();
            if (dtype == "<i2" or dtype == "i2") { // ADC short ints
                Array::array_xxs sarr;
                bool ok = pigenc::eigen::load(m_cur.pig, sarr);
                if (!ok) {
                    log->error("call={}, short load failed tag=\"{}\" file={}",
                               m_count, m_cur.tag, m_inname);
                    return nullptr;
                }
                framelet.trace_array = sarr.cast<float>();
            }
            else {
                bool ok = pigenc::eigen::load(m_cur.pig, framelet.trace_array);
                if (!ok) {
                    log->error("call={}, float load failed tag=\"{}\" file={}",
                               m_count, m_cur.tag, m_inname);
                    return nullptr;
                }
            }        
            clear();
            continue;
        }

        if (m_cur.type == "channels") {
            auto& framelet = get_framelet(m_cur.tag);
            bool ok = pigenc::stl::load(m_cur.pig, framelet.channels);
            if (!ok) {
                log->error("call={}, channel load failed tag=\"{}\" file={}",
                           m_count, m_cur.tag, m_inname);
                return nullptr;
            }
            clear();
            continue;
        }

        if (m_cur.type == "tickinfo") {
            auto& framelet = get_framelet(m_cur.tag);
            bool ok = pigenc::stl::load(m_cur.pig, framelet.tickinfo);
            if (!ok) {
                log->error("call={}, tickinfo load failed tag=\"{}\" file={}",
                           m_count, m_cur.tag, m_inname);
                return nullptr;
            }
            clear();
            continue;
        }

        log->warn("call={} skipping unsupported input: {}", m_count, m_cur.fname);
        clear();
        continue;

    } // loading loop

    if (framelets.empty()) {
        log->error("call={}, no frame arrays loaded ident={} file={}",
                   m_count, ident, m_inname);
        return nullptr;
    }

    // collect traces from framelets in order
    ITrace::vector all_traces;
    for (auto& framelet : framelets) {

        const size_t nrows = framelet.trace_array.rows();
        const size_t ncols = framelet.trace_array.cols();

        if (nrows != framelet.channels.size()) {
            log->error("call={}, mismatch in frame and channel array sizes ident={} file={}",
                       m_count, ident, m_inname);
            return nullptr;
        }

        const int tbin0 = (int)framelet.tickinfo[2];

        for (size_t irow=0; irow < nrows; ++irow) {
            int chid = framelet.channels[irow];
            auto row = framelet.trace_array.row(irow);
            ITrace::ChargeSequence charges(row.data(), row.data() + ncols);
            auto itrace = std::make_shared<Aux::SimpleTrace>(chid, tbin0, charges);
            all_traces.push_back(itrace);
        }
    }

    const double time = framelets[0].tickinfo[0];
    const double tick = framelets[0].tickinfo[1];

    auto sframe = std::make_shared<Aux::SimpleFrame>(ident, time, all_traces, tick);
    for (auto ftag : m_frame_tags) {
        sframe->tag_frame(ftag);
    }

    // Tag traces of each framelet
    size_t last_index=0;
    for (auto& framelet : framelets) {
        const size_t size = framelet.channels.size();
        std::vector<size_t> inds(size);
        std::iota(inds.begin(), inds.end(), last_index);
        last_index += size;
        sframe->tag_traces(m_cur.tag, inds);
    }        
    return sframe;
}

void FrameFileSource::clear()
{
    m_cur.pig.clear();
    m_cur = entry_t();
}

bool FrameFileSource::operator()(IFrame::pointer& frame)
{
    frame = nullptr;
    if (m_eos_sent) {
        return false;
    }
    frame = load();
    if (frame) {
        log->debug("load frame call={}", m_count ++);
    }
    else {
        log->debug("EOS at call={}", m_count ++);
        m_eos_sent = true;
    }
    
    return true;
}
