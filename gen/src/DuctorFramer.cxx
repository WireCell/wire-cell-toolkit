#include "WireCellGen/DuctorFramer.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(DuctorFramer,
                 WireCell::Gen::DuctorFramer,
                 WireCell::INamed,
                 WireCell::IDepoFramer,
                 WireCell::IConfigurable)

using namespace WireCell;

Gen::DuctorFramer::DuctorFramer()
    : Aux::Logger("DuctorFramer", "gen")
{
}

Gen::DuctorFramer::~DuctorFramer()
{
}

WireCell::Configuration Gen::DuctorFramer::default_configuration() const
{
    Configuration cfg;
    // The IDuctor that does the transformation
    cfg["ductor"] = Json::nullValue;
    // The IFrameFanin that merges frames in the case that the IDuctor returns
    // more than one over a given depo set.
    cfg["fanin"] = Json::nullValue;
    return cfg;
}

void Gen::DuctorFramer::configure(const WireCell::Configuration& cfg)
{
    auto ductor_tn = get<std::string>(cfg, "ductor", "");
    m_ductor = Factory::find_tn<IDuctor>(ductor_tn);
    auto fanin_tn = get<std::string>(cfg, "fanin", "");    
    m_fanin = Factory::find_tn<IFrameFanin>(fanin_tn);
}


bool Gen::DuctorFramer::operator()(const input_pointer& in, output_pointer& out)
{
    out=nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    IDepo::vector in_depos(in->depos()->begin(), in->depos()->end());
    in_depos.push_back(nullptr); // input EOS

    IFrameFanin::input_vector all_frames;

    for (auto idepo : in_depos) {
        IDuctor::output_queue more;        
        (*m_ductor)(idepo, more);
        all_frames.insert(all_frames.end(), more.begin(), more.end());
    }

    // The EOS should come through
    if (all_frames.back()) {
        raise<ValueError>("ductor did not pass through an EOS");
    }
    all_frames.pop_back();

    const IFrameFanin::input_vector& call_frames = all_frames;
    IFrameFanin::output_pointer op;
    bool ok = (*m_fanin)(call_frames, op);
    if (ok) {
        out = op;
    }
    else {
        log->warn("DuctorFramer fanin failed");
    }
    ++m_count;
    return true;
}

