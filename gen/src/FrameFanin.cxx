#include "WireCellGen/FrameFanin.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellAux/SimpleFrame.h"

#include <iostream>
#include <map>

WIRECELL_FACTORY(FrameFanin, WireCell::Gen::FrameFanin, WireCell::IFrameFanin, WireCell::IConfigurable)

using namespace WireCell;

Gen::FrameFanin::FrameFanin(size_t multiplicity)
    : Aux::Logger("FrameFanin", "glue")
{
}
Gen::FrameFanin::~FrameFanin() {}

WireCell::Configuration Gen::FrameFanin::default_configuration() const
{
    Configuration cfg;

    // The number of input ports.  If set to zero the fanin will have DYNAMIC
    // MULTIPLICITY and simply accept what it is given.  
    cfg["multiplicity"] = (int) m_multiplicity;

    // A non-null entry in this array is taken as a string and used to
    // tag traces which arrive on the corresponding input port when
    // they are placed to the output frame.  This creates a tagged
    // subset of output traces corresponding to the entire set of
    // traces in the input frame from the associated port.
    //
    // When DYNAMIC MULTIPLICITY is used (multiplicity=0) the last tag in this
    // array is applied to the all ports numbered greater or equal to the size
    // of the tag array.
    cfg["tags"] = Json::arrayValue;

    // Tag rules are an array, one element per input port.  Each
    // element is an object keyed with "frame" or "trace".  Each of
    // their values are an object keyed by a regular expression
    // (regex) and with values that are a single tag or an array of
    // tags. See util/test_tagrules for examples.
    //
    // When DYNAMIC MULTIPLICITY (multiplicity=0) the last tag rule in this
    // array is applied to all ports numbered greater or equal to the size of
    // the tag rules array.
    cfg["tag_rules"] = Json::arrayValue;

    // See also comments at top of FrameFanin.h

    return cfg;
}
static bool has_traces(const IFrame::pointer& fr)
{
    auto traces = fr->traces();
    return traces && !traces->empty();
}

// The ruleset index for an input port: its own, or the last one given when
// there are fewer rulesets than ports (0, giving no tags, when none).
static size_t rule_index(const tagrules::Context& ft, size_t port, const std::string& name)
{
    const size_t n = ft.nrules(name);
    if (n == 0) {
        return 0;
    }
    return port < n ? port : n - 1;
}

std::string Gen::FrameFanin::tag(size_t port) const
{
    if (m_tags.empty()) return "";
    if (port < m_tags.size()) return m_tags.at(port);
    return m_tags.back();
}

void Gen::FrameFanin::configure(const WireCell::Configuration& cfg)
{
    int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m < 0) {
        log->critical("illegal multiplicity: {}", m);
        THROW(ValueError() << errmsg{"FrameFanin multiplicity must be positive"});
    }
    if (m == 0) {
        log->debug("FrameFanin using dynamic multiplicity");
    }

    m_multiplicity = m;

    // If tag is not empty string it will be used to tag the traces from the
    // corresponding input port.
    for (const auto& jtag : cfg["tags"]) {
        if (jtag.isNull()) {
            m_tags.push_back("");
        }
        else {
            m_tags.push_back(jtag.asString());
        }
    }

    // frame/trace tagging at the port level
    m_ft.configure(cfg["tag_rules"]);
}

std::vector<std::string> Gen::FrameFanin::input_types()
{
    if (!m_multiplicity) {
        return std::vector<std::string>();
    }
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool Gen::FrameFanin::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    size_t neos = 0;
    for (const auto& fr : invec) {
        if (!fr) {
            ++neos;
        }
    }

    if (neos) {
        log->debug("EOS at call={} with {}", m_count, neos);
        ++m_count;
        return true;
    }

    const size_t multiplicity = invec.size();
    if (m_multiplicity && multiplicity != m_multiplicity) {
        log->critical("input vector size={} my multiplicity={}", multiplicity, m_multiplicity);
        THROW(ValueError() << errmsg{"input vector size mismatch"});
    }

    std::vector<IFrame::trace_list_t> by_port(multiplicity);

    std::vector<std::tuple<tagrules::tag_t, IFrame::trace_list_t, IFrame::trace_summary_t> > stash;

    std::stringstream info;
    info << "call=" << m_count << " ";

    tagrules::tagset_t fouttags;
    ITrace::vector out_traces;
    IFrame::pointer one = nullptr;
    
    Waveform::ChannelMaskMap out_cmm;
    // creating an empty name map; only want combine the two maskmaps from each APA into a single maskmap, NOT merge masks
    std::map<std::string, std::string> empty_name_map;

    for (size_t iport = 0; iport < multiplicity; ++iport) {
        const size_t trace_offset = out_traces.size();

        const auto& fr = invec[iport];
        info << "input " << iport << ": " << Aux::taginfo(fr) << " ";

        // The output frame's ident, time and tick are those of the first
        // input frame with traces: an empty frame (eg DeposOrBust's for an
        // anode without depos) has no meaningful time or tick.
        if (!one || (!has_traces(one) && has_traces(fr))) {
            one = fr;
        }
        auto traces = fr->traces();
        auto masks = fr->masks();
        Waveform::merge(out_cmm, masks, empty_name_map); 

        {  // collect output frame tags by tranforming each input frame
           // tag based on user rules.
            auto fintags = fr->frame_tags();
            fintags.push_back("");
            auto fo = m_ft.transform(rule_index(m_ft, iport, "frame"), "frame", fintags);
            fouttags.insert(fo.begin(), fo.end());
        }


        // collect transformed trace tags, any trace summary and their
        // offset trace indeices. annoying_factor *= 10;
        for (auto inttag : fr->trace_tags()) {
            tagrules::tagset_t touttags = m_ft.transform(rule_index(m_ft, iport, "trace"), "trace", inttag);
            if (touttags.empty()) {
                continue;
            }

            const auto& summary = fr->trace_summary(inttag);
            IFrame::trace_list_t outtrinds;
            for (size_t trind : fr->tagged_traces(inttag)) {
                outtrinds.push_back(trind + trace_offset);
            }
            for (auto ot : touttags) {
                // need to stash them until after creating the output frame.
                stash.push_back(std::make_tuple(ot, outtrinds, summary));
            }
        };

        if (tag(iport).size()) {
            size_t tsize = 0;
            if (traces) {
                tsize = traces->size();
            }
            IFrame::trace_list_t tl(tsize);
            std::iota(tl.begin(), tl.end(), trace_offset);
            by_port[iport] = tl;
        }

        // Insert one intput frames traces.
        if (traces) {
            out_traces.insert(out_traces.end(), traces->begin(), traces->end());
        }
    }

    auto sf = new Aux::SimpleFrame(one->ident(), one->time(), out_traces, one->tick(), out_cmm);
    for (size_t iport = 0; iport < multiplicity; ++iport) {
        if (tag(iport).size()) {
            sf->tag_traces(tag(iport), by_port[iport]);
        }
    }

    for (auto ftag : fouttags) {
        sf->tag_frame(ftag);
    }

    // Several ports may give traces to the same output tag (eg each
    // anode's "gauss" into one detector-wide "gauss"): the tag covers the
    // traces of all of them.  A summary is kept only if every contribution
    // to the tag has one aligned with its traces.
    std::map<tagrules::tag_t, std::pair<IFrame::trace_list_t, IFrame::trace_summary_t>> merged;
    std::map<tagrules::tag_t, bool> summarized;
    for (const auto& ttt : stash) {
        const auto& ot = get<0>(ttt);
        const auto& inds = get<1>(ttt);
        const auto& summ = get<2>(ttt);
        auto& [minds, msumm] = merged[ot];
        auto sit = summarized.emplace(ot, true).first;
        sit->second = sit->second && summ.size() == inds.size();
        minds.insert(minds.end(), inds.begin(), inds.end());
        msumm.insert(msumm.end(), summ.begin(), summ.end());
    }
    for (auto& [ot, is] : merged) {
        if (summarized[ot] && !is.second.empty()) {
            sf->tag_traces(ot, is.first, is.second);
        }
        else {
            sf->tag_traces(ot, is.first);
        }
    }

    out = IFrame::pointer(sf);
    info << "output: " << Aux::taginfo(out);
    log->debug(info.str());

    ++m_count;
    return true;
}
