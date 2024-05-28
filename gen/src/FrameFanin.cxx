#include "WireCellGen/FrameFanin.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellAux/SimpleFrame.h"

#include <iostream>

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

        if (!one) {
            one = fr;
        }
        auto traces = fr->traces();
        auto masks = fr->masks();
        Waveform::merge(out_cmm, masks, empty_name_map); 

        {  // collect output frame tags by tranforming each input frame
           // tag based on user rules.
            auto fintags = fr->frame_tags();
            fintags.push_back("");
            const size_t n = m_ft.nrules("frame");
            auto fo = m_ft.transform(iport >= n ? n-1 : iport, "frame", fintags);
            fouttags.insert(fo.begin(), fo.end());
        }


        // collect transformed trace tags, any trace summary and their
        // offset trace indeices. annoying_factor *= 10;
        for (auto inttag : fr->trace_tags()) {
            const size_t n = m_ft.nrules("frame");
            tagrules::tagset_t touttags = m_ft.transform(iport >= n ? n-1 : iport, "trace", inttag);
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
            IFrame::trace_list_t tl(traces->size());
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

    for (auto ttt : stash) {
        sf->tag_traces(get<0>(ttt), get<1>(ttt), get<2>(ttt));
    }

    out = IFrame::pointer(sf);
    info << "output: " << Aux::taginfo(out);
    log->debug(info.str());

    ++m_count;
    return true;
}
