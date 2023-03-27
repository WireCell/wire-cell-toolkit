#include "WireCellSigProc/NoiseModeler.h"

#include "WireCellAux/FrameTools.h"

#include "WireCellIface/ITraceRanker.h"
#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Interpolate.h"
#include "WireCellUtil/NamedFactory.h"

#include <vector>

WIRECELL_FACTORY(NoiseModeler, WireCell::SigProc::NoiseModeler,
                 WireCell::INamed,
                 WireCell::IFrameSink, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SigProc;

NoiseModeler::NoiseModeler()
    : Aux::Logger("NoiseModeler", "sigproc")
{
}

NoiseModeler::~NoiseModeler() 
{
}

WireCell::Configuration NoiseModeler::default_configuration() const
{
    Configuration cfg;

    cfg["dft"] = "FftwDFT";

    // type:name of a trace noise ranker
    cfg["isnoise"] = "NoiseRanker";
    // A rank from NoiseRanker less than this fails to be considered
    // noise.
    cfg["threshold"] = 0.9;

    // The "groups" provides an array of objects with attributes
    // groupID, an integer, and channels, an array of integers or it
    // provides the name of a file which provides this data.  If null,
    // a single group with ID 0 is assumed.  Any channels which are
    // not resolved to a group will be resolved to group ID 0.  
    cfg["groups"] = Json::nullValue;

    // Number of frequency samples to use in collecting the spectra.
    // If unset, the smallest power-of-2 greater than the number of
    // ticks in the first trace encountered will be used.
    cfg["nfft"] = Json::nullValue;
    // The number of samples of the mean half spectrum to output.  If
    // not set, smallest power of two larger than sqrt(nticks) is
    // used.
    cfg["nhalfout"] = Json::nullValue;

    // The name of file to write at terminate time.  It should have an
    // extension indicating uncompressed .json or compressed
    // .json.bz2.
    //
    // See aux/docs/noise.org for details of this format.
    cfg["outname"] = Json::nullValue;

    return cfg;
}


void NoiseModeler::configure(const WireCell::Configuration& cfg)
{
    int errors = 0;

    auto dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_nfft = get(cfg, "nfft", 0);

    m_chgrp.clear();
    auto jgroups = cfg["groups"];
    if (jgroups.isNull()) {
        log->debug("no \"groups\" parameter given, will produce single group output");
    }
    else {
        if (jgroups.isString()) {
            jgroups = Persist::load(jgroups.asString());
        }
        for (const auto& jg : jgroups) {
            const int gid = jg["groupID"].asInt();
            for (auto jch : jg["channels"]) {
                m_chgrp[jch.asInt()] = gid;
            }
        }
    }

    // Must be at least this high to be noise.
    auto threshold = get(cfg, "threshold", 0.9);
    auto isnoise_tn = get<std::string>(cfg, "isnoise", "NoiseRanker");
    if (isnoise_tn.empty()) {
        log->critical("configuration parameter \"isnoise\" is empty");
        ++errors;
    }
    auto isnoise = Factory::find_tn<ITraceRanker>(isnoise_tn);

    auto jo = cfg["outname"];
    if (jo.isNull()) {
        log->critical("configuration parameter \"outname\" is empty");
        ++errors;
    }
    m_outname = jo.asString();

    if (errors) {
        THROW(ValueError() << errmsg{"configuration error"});
    }

    m_isnoise = [=](const ITrace::pointer& trace) -> bool {
        const double rank = isnoise->rank(trace);
        return rank >= threshold;
    };        
        
}

NoiseModeler::Collector& NoiseModeler::nc(int chid)
{
    int group = 0;
    auto git = m_chgrp.find(chid);
    if (git != m_chgrp.end()) {
        group = git->second;
    }
     
    if (m_gcol.empty()) {
        std::set<int> ugroups;        
        for (auto [c,g] : m_chgrp) {
            ugroups.insert(g);
        }
        // ugroups.insert(0);      // why this?
        for (int g : ugroups) {
            m_gcol.emplace(g, Collector(m_dft, m_nfft));
        }
    }
    return m_gcol[group];
}


// we expect voltage-level frames
bool NoiseModeler::operator()(const input_pointer& frame)
{
    if (!frame) {
        log->debug("EOS at count={}", m_count);
        ++m_count;
        return true;
    }

    auto traces = Aux::untagged_traces(frame);
    if (traces.empty()) {
        log->debug("no traces in input frame {} at call={}",
                   frame->ident(), m_count);
        ++m_count;
        return true;
    }
    
    if (m_tick == 0) {
        m_tick = frame->tick();
    }

    int got=0;
    for (const auto& trace : traces) {
        const auto& vwave = trace->charge();

        // default set-on-first-seen pattern
        if (! m_nfft) {
            m_nfft = vwave.size();
            log->debug("set nfft={} based on first seen", m_nfft);
        }
        if (!m_nhalfout) {
            m_nhalfout = ceil(sqrt(m_nfft));
            log->debug("set nsave={} based on first seen", m_nhalfout);
        }

        if (! m_isnoise(trace)) {
            continue;
        }

        ++got;
        auto& col = nc(trace->channel());
        col.add(vwave.begin(), vwave.end());
    }

    log->debug("got {} noise traces out of {} at frame {} in {} groups at count={}",
               got, traces.size(), frame->ident(), m_gcol.size(), m_count);

    ++m_count;
    return true;
}


void NoiseModeler::finalize()
{
    const size_t nout = 2*m_nhalfout;
    Json::Value jout = Json::arrayValue;
    for (const auto& [gid, col] : m_gcol) {

        if (col.nwaves() == 0) {
            log->warn("collected no waves for group {}", gid);
            continue;
        }

        Json::Value jentry = Json::objectValue;

        // Set the group number.  ("groupID" is deprecated)
        jentry["group"] = gid;

        // This is the original sampling period and file needs it so
        // user can know Fnyquist, etc.
        jentry["period"] = m_tick;

        // Note, "nsamples" means at least three different things
        // depending on context!

        // In the data file, "nsamples" means the number of samples
        // (ticks) in the original time series waveform.  The
        // Collector truncates all waves to the size of the first wave
        // size and that's just what the file wants.
        const int nticks = col.nticks(); 
        jentry["nsamples"] = nticks;

        // But, the collected spectra are size nfft which may differ
        // from nticks.  
        auto camps = col.amplitude();

        // To interpolate between the collect size nfft and the
        // reqquested output size 2*nhalfout we need the corresponding
        // Rayleigh frequencies (bin size).
        const double Fcol = 1.0/(m_tick*m_nfft);
        const double Fout = 1.0/(m_tick*nout);

        // Interpolate full spectrum nfft -> nhalfout
        linterp<float> terp(camps.begin(), camps.end(), 0, Fcol);
        std::vector<float> amps(m_nhalfout);
        terp(amps.begin(), m_nhalfout, 0, Fout);

        // copy in
        Json::Value jamps = Json::arrayValue;
        Json::Value jfreq = Json::arrayValue;
        for (size_t ind=0; ind < m_nhalfout; ++ind) {
            jamps.append(amps[ind]);
            jfreq.append(ind * Fout);
        }
        jentry["freqs"] = jfreq;
        jentry["amps"] = jamps;
        jout.append(jentry);
        log->debug("group {} nwaves={} nsamples={} nfft={}, nsave={}",
                   gid, col.nwaves(), nticks, m_nfft, m_nhalfout);
    }
    Persist::dump(m_outname, jout, true);
}
