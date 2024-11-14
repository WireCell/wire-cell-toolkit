#include "WireCellSigProc/OmnibusNoiseFilter.h"

#include "WireCellSigProc/Diagnostics.h"

#include "WireCellUtil/Response.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include "WireCellUtil/NamedFactory.h"
// #include "WireCellUtil/ExecMon.h" // debugging

#include "WireCellAux/FrameTools.h"

#include <unordered_map>

WIRECELL_FACTORY(OmnibusNoiseFilter,
                 WireCell::SigProc::OmnibusNoiseFilter,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

using namespace WireCell;

using namespace WireCell::SigProc;

OmnibusNoiseFilter::OmnibusNoiseFilter(std::string intag, std::string outtag)
  : Aux::Logger("OmnibusNoiseFilter", "sigproc")
  , m_nticks(0)
  , m_intag(intag)    // orig
  , m_outtag(outtag)  // raw
{
}
OmnibusNoiseFilter::~OmnibusNoiseFilter() {}

void OmnibusNoiseFilter::configure(const WireCell::Configuration& cfg)
{
    // std::cerr << "OmnibusNoiseFilter: configuring with:\n" << cfg << std::endl;
    m_nticks = (size_t) get(cfg, "nticks", (int) m_nticks);
    if (!cfg["nsamples"].isNull()) {
        log->warn("\"nsamples\" is an obsolete parameter, use \"nticks\"");
        // we don't throw here on assumption that nticks is still provided
    }

    auto jmm = cfg["maskmap"];
    for (auto name : jmm.getMemberNames()) {
        m_maskmap[name] = jmm[name].asString();
        //	std::cerr << name << " " << m_maskmap[name] << std::endl;
    }

    for (auto jf : cfg["channel_filters"]) {
        auto filt = Factory::find_tn<IChannelFilter>(jf.asString());
        log->debug("adding channel filter: {} \"{}\"", m_perchan.size(), jf.asString());
        m_perchan.push_back(filt);
    }
    for (auto jf : cfg["channel_status_filters"]) {
        auto filt = Factory::find_tn<IChannelFilter>(jf.asString());
        log->debug("adding channel status filter: {} \"{}\"", m_perchan_status.size(),
                   jf.asString());
        m_perchan_status.push_back(filt);
    }

    auto jcndb = cfg["noisedb"];
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(jcndb.asString());
    log->debug("using channel noise DB object: \"{}\"", jcndb.asString());
    log->warn("Use of CHND is deprecated and will be removed in the future.");

    // We should use multi-group channel filters, the interface below is maintained
    // to ensure backwards compatibility. See:
    // https://github.com/WireCell/wire-cell-toolkit/issues/327
    if ( !cfg["grouped_filters"].empty() ) {
        MGCF mgcf;
        mgcf.channelgroups = m_noisedb->coherent_channels();
        for (auto jf : cfg["grouped_filters"]) {
            auto filt = Factory::find_tn<IChannelFilter>(jf.asString());
            log->debug("adding grouped filter: {} \"{}\"", mgcf.filters.size(), jf.asString());
            mgcf.filters.push_back(filt);
        }
        m_multigroup_chanfilters.push_back(mgcf);
    }

    for (auto jmgcf: cfg["multigroup_chanfilters"]) {
        MGCF mgcf;
        auto jgroups = jmgcf["channelgroups"];
        for (auto jgroup : jgroups) {
            std::vector<int> channel_group;
            for (auto jch : jgroup) {
                channel_group.push_back(jch.asInt());
            }
            mgcf.channelgroups.push_back(channel_group);
        }
        for (auto jf: jmgcf["filters"]){
            auto filt = Factory::find_tn<IChannelFilter>(jf.asString());
            log->debug("adding grouped filter: {} \"{}\"", mgcf.filters.size(), jf.asString());
            mgcf.filters.push_back(filt);
        }
        m_multigroup_chanfilters.push_back(mgcf);
    }

    m_intag = get(cfg, "intraces", m_intag);
    m_outtag = get(cfg, "outtraces", m_outtag);
}

WireCell::Configuration OmnibusNoiseFilter::default_configuration() const
{
    Configuration cfg;
    cfg["nticks"] = (int) m_nticks;
    cfg["maskmap"]["chirp"] = "bad";
    cfg["maskmap"]["noisy"] = "bad";

    cfg["channel_filters"][0] = "mbOneChannelNoise";
    cfg["channel_status_filters"][0] = "mbOneChannelStatus";
    cfg["grouped_filters"][0] = "mbCoherentNoiseSub";

    // user must supply.  "OmniChannelNoiseDB" is a likely choice.
    // Avoid SimpleChannelNoiseDB.
    cfg["noisedb"] = "";

    // The tags for input and output traces
    cfg["intraces"] = m_intag;
    cfg["outtraces"] = m_outtag;
    return cfg;
}

bool OmnibusNoiseFilter::operator()(const input_pointer& inframe, output_pointer& outframe)
{
    if (!inframe) {  // eos
        outframe = nullptr;
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }
    log->debug("call={} input frame: {}", m_count, Aux::taginfo(inframe));

    auto traces = Aux::tagged_traces(inframe, m_intag);
    if (traces.empty()) {
        log->warn("no traces for tag \"{}\", sending empty frame at call={}",
                  m_intag, m_count);
        outframe = std::make_shared<Aux::SimpleFrame>(
            inframe->ident(), inframe->time(),
            std::make_shared<ITrace::vector>(), inframe->tick());
        ++m_count;
        return true;
    }

    if (! m_nticks) {
        // Warning: this implicitly assumes a dense frame (ie, all tbin=0 and all waveforms same size).
        // It also won't stop triggering a warning inside OneChannelNoise if there is a mismatch.
        m_nticks = traces.at(0)->charge().size();
    }

    // For now, just collect any and all masks and interpret them as "bad".
    Waveform::ChannelMaskMap input_cmm = inframe->masks();
    Waveform::ChannelMaskMap cmm;
    Waveform::merge(cmm, input_cmm, m_maskmap);

    // Get the ones from database and then merge
    std::vector<int> bad_channels = m_noisedb->bad_channels();
    {
        Waveform::BinRange bad_bins;
        bad_bins.first = 0;
        bad_bins.second = (int) m_nticks;
        Waveform::ChannelMasks temp;
        for (size_t i = 0; i < bad_channels.size(); i++) {
            temp[bad_channels.at(i)].push_back(bad_bins);
            // std::cout << temp.size() << " " << temp[bad_channels.at(i)].size() << std::endl;
        }
        Waveform::ChannelMaskMap temp_map;
        temp_map["bad"] = temp;
        Waveform::merge(cmm, temp_map, m_maskmap);
    }

    int nchanged_samples = 0;

    // Collect our working area indexed by channel.
    std::unordered_map<int, Aux::SimpleTrace*> bychan;
    for (auto trace : traces) {
        int ch = trace->channel();

        // make working area directly in simple trace to avoid memory fragmentation
        auto signal = new Aux::SimpleTrace(ch, 0, m_nticks);
        bychan[ch] = signal;

        // if good
        if (find(bad_channels.begin(), bad_channels.end(), ch) == bad_channels.end()) {
            auto const& charge = trace->charge();
            const size_t ncharges = charge.size();

            signal->charge().assign(charge.begin(), charge.begin() + std::min(m_nticks, ncharges));
            signal->charge().resize(m_nticks, 0.0);

            if (ncharges != m_nticks) {
                nchanged_samples += std::abs((int) m_nticks - (int) ncharges);
            }
        }

        int filt_count = 0;
        for (auto filter : m_perchan) {
            auto masks = filter->apply(ch, signal->charge());

            // fixme: probably should assure these masks do not lead to out-of-bounds...

            Waveform::merge(cmm, masks, m_maskmap);
            ++filt_count;
        }
    }
    traces.clear();  // done with our copy of vector of shared pointers

    if (nchanged_samples) {
        log->warn("warning, truncated or extended {} samples", nchanged_samples);
    }

    // int group_counter = 0;
    int nunknownchans = 0;
    for (const auto& mgcf : m_multigroup_chanfilters) {
        for (auto group : mgcf.channelgroups) {

            int flag = 1;

            IChannelFilter::channel_signals_t chgrp;
            for (auto ch : group) {  // fix me: check if we don't actually have this channel
                                    // std::cout << group_counter << " " << ch << " " << std::endl;
                if (bychan.find(ch) == bychan.end()) {
                    ++nunknownchans;
                    flag = 0;
                }
                else {
                    chgrp[ch] = bychan[ch]->charge();  // copy...
                }
            }
            if (flag == 0) continue;

            for (auto filter : mgcf.filters) {
                auto masks = filter->apply(chgrp);
                Waveform::merge(cmm, masks, m_maskmap);
            }

            for (auto cs : chgrp) {
                // cs.second; // copy
                bychan[cs.first]->charge().assign(cs.second.begin(), cs.second.end());
            }

        } // end of channel groups

    } // end of MGCF

    if (nunknownchans) {
        log->debug("{} unknown channels (probably the channel selector is in use)", nunknownchans);
    }

    // run status
    for (auto& it : bychan) {
        const int ch = it.first;
        IChannelFilter::signal_t& signal = it.second->charge();
        for (auto filter : m_perchan_status) {
            auto masks = filter->apply(ch, signal);

            Waveform::merge(cmm, masks, m_maskmap);
        }
    }

    ITrace::vector itraces;
    for (auto& cs : bychan) {  // fixme: that tbin though
        itraces.push_back(ITrace::pointer(cs.second));
    }

    bychan.clear();

    auto sframe = new Aux::SimpleFrame(inframe->ident(), inframe->time(), itraces, inframe->tick(), cmm);
    IFrame::trace_list_t indices(itraces.size());
    for (size_t ind = 0; ind < itraces.size(); ++ind) {
        indices[ind] = ind;
    }
    if (! m_outtag.empty()) {
        sframe->tag_traces(m_outtag, indices);
    }
    sframe->tag_frame("noisefilter"); // fixme: this is unnecessary and should be removed 
    outframe = IFrame::pointer(sframe);

    log->debug("call={} output frame: {}", m_count, Aux::taginfo(outframe));
    log->debug("call={}, frame={}, ntraces={}, nticks={} intag={} outtag={}",
               m_count,
               sframe->ident(), itraces.size(), m_nticks,
               m_intag, m_outtag);

    ++m_count;
    return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
