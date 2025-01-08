#include "WireCellGen/FrameFanout.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"

WIRECELL_FACTORY(FrameFanout, WireCell::Gen::FrameFanout,
                 WireCell::INamed,
                 WireCell::IFrameFanout, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::SimpleFrame;

Gen::FrameFanout::FrameFanout(size_t multiplicity)
    : Aux::Logger("FrameFanout", "glue")
    , m_multiplicity(multiplicity)
{
}
Gen::FrameFanout::~FrameFanout() {}

WireCell::Configuration Gen::FrameFanout::default_configuration() const
{
    Configuration cfg;
    // How many output ports
    cfg["multiplicity"] = (int) m_multiplicity;
    // If "tag_rules" is not given or empty then the frame is fanned
    // out in a trivial manner.
    // 
    // Tag rules are an array, one element per output port.  Each
    // element is an object keyed with "frame" or "trace".  Each of
    // their values are an object keyed by a regular expression
    // (regex) and with values that are a single tag or an array of
    // tags.
    //
    // The frame traces are forwarded in any case, tag filtering only
    // affects tags (and their trace indices in the case of trace
    // tags).
    //
    // Each frame tag found in the input is sent through the rules and
    // the resulting output tags are added to the set of frame tags of
    // the output.
    //
    // Trace tag work similarly.  Each trace tag is sent through the
    // rules and tags that the trace rule returns will be placed in
    // the output frame along with the trace index vector of the input
    // trace tag. 
    // 
    // See util/test_tagrules for examples.
    cfg["tag_rules"] = Json::arrayValue;
    return cfg;
}
void Gen::FrameFanout::configure(const WireCell::Configuration& cfg)
{
    int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m <= 0) {
        THROW(ValueError() << errmsg{"FrameFanout multiplicity must be positive"});
    }
    m_multiplicity = m;

    auto tr = cfg["tag_rules"];
    if (tr.isNull() or tr.empty()) {
        m_trivial = true;
        log->debug("trivial {}-way fanout", m_multiplicity);
        return;
    }
    m_trivial = false;
    m_ft.configure(tr);
    log->debug("rule based {}-way fanout", m_multiplicity);
}

std::vector<std::string> Gen::FrameFanout::output_types()
{
    const std::string tname = std::string(typeid(output_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool Gen::FrameFanout::operator()(const input_pointer& in, output_vector& outv)
{
    outv.resize(m_multiplicity);

    std::stringstream info;
    info << "call=" << m_count << ": ";
    ++m_count;

    if (!in) {  //  pass on EOS
        for (size_t ind = 0; ind < m_multiplicity; ++ind) {
            outv[ind] = in;
        }
        info << "see EOS";
        log->debug(info.str());
        return true;
    }

    if (m_trivial) {
        info << "input->output x"<<m_multiplicity<<": " << Aux::taginfo(in);
        log->debug(info.str());

        for (size_t ind = 0; ind < m_multiplicity; ++ind) {
            outv[ind] = in;
        }
        return true;
    }

    // O.w. we are rule driven.

    info << "input: " << Aux::taginfo(in) << " ";

    auto fintags = in->frame_tags();
    // Add empty tag to allow rules to create an output tag even if none exist.
    // No equivalent should be made for trace tags.
    fintags.push_back("");

    for (size_t ind = 0; ind < m_multiplicity; ++ind) {
        // Basic frame stays the same.
        auto sfout = new SimpleFrame(in->ident(), in->time(), *in->traces(), in->tick(), in->masks());

        // Transform any frame tags based on a per output port ruleset
        auto fouttags = m_ft.transform(ind, "frame", fintags);

        for (auto ftag : fouttags) {
            sfout->tag_frame(ftag);
        }

        for (auto inttag : in->trace_tags()) {
            tagrules::tagset_t touttags = m_ft.transform(ind, "trace", inttag);
            if (touttags.empty()) {
                continue;
            }
            const auto& traces = in->tagged_traces(inttag);
            const auto& summary = in->trace_summary(inttag);

            for (auto otag : touttags) {
                sfout->tag_traces(otag, traces, summary);
            }
        };

        auto out = IFrame::pointer(sfout);
        info << "output " << ind << ": " << Aux::taginfo(out) << " ";
        outv[ind] = out;
    }

    log->debug(info.str());

    return true;
}
