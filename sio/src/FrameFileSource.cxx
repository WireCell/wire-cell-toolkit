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
        raise<ValueError>("FrameFileSource: unsupported inname: %s", m_inname);
    }

    m_tags.clear();
    for (auto jtag : cfg["tags"]) {
        m_tags.push_back(jtag.asString());
    }

    m_frame_tags.clear();
    for (auto jtag : cfg["frame_tags"]) {
        m_frame_tags.push_back(jtag.asString());
    }
}

bool FrameFileSource::is_excluded(const std::string& tag)
{
    if (is_tagged(tag)) return false;
    if (is_untagged(tag)) return false;
    return true;
}

bool FrameFileSource::is_untagged(const std::string& tag)
{
    return tag.empty() or tag == "*";
}

bool FrameFileSource::is_tagged(const std::string& tag)
{
    if (is_untagged(tag)) {
        // avoid misinterpreting "" or "*" in m_tags
        return false;
    }

    for (const auto& maybe : m_tags) {
        if (maybe == tag) {
            return true;
        }
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
    log->debug("read type={} tag={} ident={} ext={} at call={}",
               m_cur.type, m_cur.tag, m_cur.ident, m_cur.ext, m_count);
    return true;
}

IFrame::pointer FrameFileSource::load()
{
    int ident = -1;

    // Collect portions of the frame data;
    using trace_array_t = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    struct framelet_t { 
        std::vector<double> tickinfo, summary;
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
            log->trace("call={}, make framelet for \"{}\" -> {}", m_count, tag, oldsize);
            return framelets.back();
        }
        log->trace("call={}, get framelet for \"{}\" -> {}", m_count, tag, it->second);
        return framelets[it->second];
    };

    // fill this along the way
    ChannelMaskMap cmm;

    while (true) {

        // Loop over arrays, falling through until we recognize the array type.

        // fixme: throw on errors but return nullptr on EOF.

        // Read next element if current cursor is empty
        if (m_cur.fsize == 0) {
            this->read();
            log->trace("call={}, read fsize={} ident={}, okay={}", m_count, m_cur.fsize, m_cur.ident, m_cur.okay);

            if (m_cur.fsize == 0) {
                // zero read means EOF
                log->trace("call={}, zero read.  EOS", m_count);
                break;
            }

            if (!m_in) {
                log->error("call={}, bad pig read with file={}", m_count, m_inname);
                THROW(IOError() << errmsg{"bad pig read with file " + m_inname});
            }
            if (!m_cur.okay) {
                log->error("call={}, failed to parse npy file name {} in file={}",
                           m_count, m_cur.fname, m_inname);
                THROW(IOError() << errmsg{"numpy parse error with file " + m_inname});
            }
        }

        if (ident < 0) {
            // first time
            ident = m_cur.ident;
        }
        if (ident != m_cur.ident) {
            // we see next frame, done with this current
            log->trace("call={}, found next frame: {} != {}", m_count, ident, m_cur.ident);
            break;       
        }

        log->trace("call={}, loading type \"{}\", tag \"{}\"",
                   m_count, m_cur.type, m_cur.tag);

        if (m_cur.type == "chanmask") {
            Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> cmsarr;
            bool ok = pigenc::eigen::load(m_cur.pig, cmsarr);
            if (!ok) {
                log->error("call={}, cmm load failed tag=\"{}\" file={}",
                           m_count, m_cur.tag, m_inname);
                THROW(IOError() << errmsg{"numpy parse error of chanmask with file " + m_inname});
            }
            const int ncols = cmsarr.cols();
            if (ncols != 3) {
                log->error("call={}, cmm wrong size nrows={}!=3 tag=\"{}\" file={}",
                           m_count, ncols, m_cur.tag, m_inname);
                THROW(IOError() << errmsg{"wrong size chanmask with file " + m_inname});
            }
            ChannelMasks cms;
            int nrows = cmsarr.rows();
            for (int irow=0; irow<nrows; ++irow) {
                auto row = cmsarr.row(irow);
                cms[row[0]].push_back(std::make_pair(row[1], row[2]));
            }
            cmm[m_cur.tag] = cms;
            log->debug("call={}, add chanmask of size {} for tag \"{}\"",
                       m_count, cms.size(), m_cur.tag);
            clear();
            continue;
        }
    
        if (is_excluded(m_cur.tag)) {
            log->debug("call={}, skipping tag=\"{}\" file={}",
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
                    THROW(IOError() << errmsg{"int frame read error with file " + m_inname});
                }
                framelet.trace_array = sarr.cast<float>();
            }
            else {
                bool ok = pigenc::eigen::load(m_cur.pig, framelet.trace_array);
                if (!ok) {
                    log->error("call={}, float load failed tag=\"{}\" file={}",
                               m_count, m_cur.tag, m_inname);
                    THROW(IOError() << errmsg{"float frame read error with file " + m_inname});
                }
            }        
            log->trace("call={}, load {} traces in frame with tag=\"{}\" have {}",
                       m_count, framelet.trace_array.rows(), m_cur.tag, framelets.size());
            clear();
            continue;
        }

        if (m_cur.type == "channels") {
            auto& framelet = get_framelet(m_cur.tag);
            bool ok = pigenc::stl::load(m_cur.pig, framelet.channels);
            if (!ok) {
                log->error("call={}, channel load failed tag=\"{}\" file={}",
                           m_count, m_cur.tag, m_inname);
                THROW(IOError() << errmsg{"channel load error with file " + m_inname});
            }
            log->trace("call={}, load {} channels with tag=\"{}\" have {}",
                       m_count, framelet.channels.size(), m_cur.tag, framelets.size());
            clear();
            continue;
        }

        if (m_cur.type == "tickinfo") {
            auto& framelet = get_framelet(m_cur.tag);
            bool ok = pigenc::stl::load(m_cur.pig, framelet.tickinfo);
            if (!ok) {
                log->error("call={}, tickinfo load failed tag=\"{}\" file={}",
                           m_count, m_cur.tag, m_inname);
                THROW(IOError() << errmsg{"tickinfo load error with file " + m_inname});
            }
            log->trace("call={}, load tickinfo with tag=\"{}\" have {}",
                       m_count, m_cur.tag, framelets.size());
            clear();
            continue;
        }

        if (m_cur.type == "summary") {
            auto& framelet = get_framelet(m_cur.tag);
            bool ok = pigenc::stl::load(m_cur.pig, framelet.summary);
            if (!ok) {
                log->error("call={}, summary load failed tag=\"{}\" file={}",
                           m_count, m_cur.tag, m_inname);
                THROW(IOError() << errmsg{"summary load error with file " + m_inname});
            }
            log->trace("call={}, load trace with tag=\"{}\" have {}",
                       m_count, m_cur.tag, framelets.size());
            clear();
            continue;
        }

        log->warn("call={} skipping unsupported input: {}",
                  m_count, m_cur.fname);
        clear();
        continue;

    } // loading loop

    if (framelets.empty()) {
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
            THROW(IOError() << errmsg{"frame channel mismatch with file " + m_inname});
        }

        const int tbin0 = (int)framelet.tickinfo[2];

        for (size_t irow=0; irow < nrows; ++irow) {
            int chid = framelet.channels[irow];
            auto row = framelet.trace_array.row(irow);
            ITrace::ChargeSequence charges(row.data(), row.data() + ncols);
            auto itrace = std::make_shared<Aux::SimpleTrace>(chid, tbin0, charges);
            all_traces.push_back(itrace);
        }
        log->trace("call={}, add {} traces from \"{}\" to total {}",
                   m_count, nrows, framelet.tag, all_traces.size());
    }

    const double time = framelets[0].tickinfo[0];
    const double tick = framelets[0].tickinfo[1];

    auto sframe = std::make_shared<Aux::SimpleFrame>(ident, time, all_traces, tick, cmm);
    for (auto ftag : m_frame_tags) {
        sframe->tag_frame(ftag);
    }

    // Tag traces of each framelet
    size_t last_index=0;
    for (auto& framelet : framelets) {
        const size_t size = framelet.channels.size();

        if (is_untagged(framelet.tag)) {
            log->trace("call={}, {} untagged with (\"{}\")", m_count, size, framelet.tag);
            last_index += size;
            continue;
        }

        std::vector<size_t> inds(size);
        std::iota(inds.begin(), inds.end(), last_index);
        last_index += size;
        if (framelet.summary.empty()) {
            sframe->tag_traces(framelet.tag, inds);
        }
        else {
            sframe->tag_traces(framelet.tag, inds, framelet.summary);
        }
        log->trace("call={}, tag {} traces with \"{}\"", m_count, size, framelet.tag);
    }        

    log->trace("call={}, loaded frame with {} tags, {} traces",
               m_count, framelets.size(), all_traces.size());

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

    frame = load();             // throws

    if (frame) {
        log->debug("call={} load frame: {}", m_count++, Aux::taginfo(frame));
    }
    else {
        log->debug("EOS at call={}", m_count++);
        m_eos_sent = true;
    }
    
    return true;
}
