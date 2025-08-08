#include "WireCellGen/Fourdee.h"

#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ConfigManager.h"
#include "WireCellUtil/ExecMon.h"

#include "WireCellGen/GenPipeline.h"

WIRECELL_FACTORY(FourDee, WireCell::Gen::Fourdee, WireCell::IApplication, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Gen::Fourdee::Fourdee()
  : m_depos(nullptr)
  , m_drifter(nullptr)
  , m_ductor(nullptr)
  , m_dissonance(nullptr)
  , m_digitizer(nullptr)
  , m_output(nullptr)
{
}

Gen::Fourdee::~Fourdee() {}

WireCell::Configuration Gen::Fourdee::default_configuration() const
{
    Configuration cfg;

    // the 4 d's and proof the developer can not count:
    put(cfg, "DepoSource", "TrackDepos");
    put(cfg, "DepoFilter", "");
    put(cfg, "Drifter", "Drifter");
    put(cfg, "Ductor", "Ductor");
    put(cfg, "Dissonance", "SilentNoise");
    put(cfg, "Digitizer", "Digitizer");
    put(cfg, "Filter", "");
    put(cfg, "FrameSink", "DumpFrames");

    return cfg;
}

void Gen::Fourdee::configure(const Configuration& thecfg)
{
    std::string tn = "";
    Configuration cfg = thecfg;

    cerr << "Gen::Fourdee:configure:\n";

    tn = get<string>(cfg, "DepoSource");
    cerr << "\tDepoSource: " << tn << endl;
    m_depos = Factory::find_maybe_tn<IDepoSource>(tn);

    tn = get<string>(cfg, "DepoFilter");
    if (tn.empty()) {
        m_depofilter = nullptr;
        cerr << "\tDepoFilter: none\n";
    }
    else {
        cerr << "\tDepoFilter: " << tn << endl;
        m_depofilter = Factory::find_maybe_tn<IDepoFilter>(tn);
    }

    tn = get<string>(cfg, "Drifter");
    cerr << "\tDrifter: " << tn << endl;
    m_drifter = Factory::find_maybe_tn<IDrifter>(tn);

    tn = get<string>(cfg, "Ductor");
    cerr << "\tDuctor: " << tn << endl;
    m_ductor = Factory::find_maybe_tn<IDuctor>(tn);

    tn = get<string>(cfg, "Dissonance", "");
    if (tn.empty()) {  // noise is optional
        m_dissonance = nullptr;
        cerr << "\tDissonance: none\n";
    }
    else {
        m_dissonance = Factory::find_maybe_tn<IFrameSource>(tn);
        cerr << "\tDissonance: " << tn << endl;
    }

    tn = get<string>(cfg, "Digitizer", "");
    if (tn.empty()) {  // digitizer is optional, voltage saved w.o. it.
        m_digitizer = nullptr;
        cerr << "\tDigitizer: none\n";
    }
    else {
        m_digitizer = Factory::find_maybe_tn<IFrameFilter>(tn);
        cerr << "\tDigitizer: " << tn << endl;
    }

    tn = get<string>(cfg, "Filter", "");
    if (tn.empty()) {  // filter is optional
        m_filter = nullptr;
        cerr << "\tFilter: none\n";
    }
    else {
        m_filter = Factory::find_maybe_tn<IFrameFilter>(tn);
        cerr << "\tFilter: " << tn << endl;
    }

    tn = get<string>(cfg, "FrameSink", "");
    if (tn.empty()) {  // sink is optional
        m_output = nullptr;
        cerr << "\tSink: none\n";
    }
    else {
        m_output = Factory::find_maybe_tn<IFrameSink>(tn);
        cerr << "\tSink: " << tn << endl;
    }
}

void dump(const IFrame::pointer frame)
{
    if (!frame) {
        cerr << "Fourdee: dump: empty frame\n";
        return;
    }

    for (auto tag : frame->frame_tags()) {
        const auto& tlist = frame->tagged_traces(tag);
        cerr << "Fourdee: frame tag: " << tag << " with " << tlist.size() << " traces\n";
    }
    for (auto tag : frame->trace_tags()) {
        const auto& tlist = frame->tagged_traces(tag);
        cerr << "Fourdee: trace tag: " << tag << " with " << tlist.size() << " traces\n";
    }

    auto traces = frame->traces();
    const int ntraces = traces->size();

    if (ntraces <= 0) {
        cerr << "Fourdee: dump: no traces\n";
        return;
    }

    std::vector<int> tbins, tlens;
    for (auto trace : *traces) {
        const int tbin = trace->tbin();
        tbins.push_back(tbin);
        tlens.push_back(tbin + trace->charge().size());
    }

    int tmin = *(std::minmax_element(tbins.begin(), tbins.end()).first);
    int tmax = *(std::minmax_element(tlens.begin(), tlens.end()).second);

    cerr << "frame: #" << frame->ident() << " @" << frame->time() / units::ms << "ms with " << ntraces
         << " traces, tbins in: "
         << "[" << tmin << "," << tmax << "]" << endl;
}

template <typename DEPOS>
void dump(DEPOS& depos)
{
    if (depos.empty() or (depos.size() == 1 and !depos[0])) {
        std::cerr << "Fourdee::dump: empty depos set\n";
        return;
    }

    std::vector<double> t, x, y, z;
    double qtot = 0.0;
    double qorig = 0.0;

    for (auto depo : depos) {
        if (!depo) {
            cerr << "Gen::Fourdee: null depo" << endl;
            break;
        }
        auto prior = depo->prior();
        if (!prior) {
            cerr << "Gen::Fourdee: null prior depo" << endl;
        }
        else {
            qorig += prior->charge();
        }
        t.push_back(depo->time());
        x.push_back(depo->pos().x());
        y.push_back(depo->pos().y());
        z.push_back(depo->pos().z());
        qtot += depo->charge();
    }

    auto tmm = std::minmax_element(t.begin(), t.end());
    auto xmm = std::minmax_element(x.begin(), x.end());
    auto ymm = std::minmax_element(y.begin(), y.end());
    auto zmm = std::minmax_element(z.begin(), z.end());
    const int ndepos = depos.size();

    std::cerr << "Gen::FourDee: drifted " << ndepos << ", extent:\n"
              << "\tt in [ " << (*tmm.first) / units::us << "," << (*tmm.second) / units::us << "]us,\n"
              << "\tx in [" << (*xmm.first) / units::mm << "," << (*xmm.second) / units::mm << "]mm,\n"
              << "\ty in [" << (*ymm.first) / units::mm << "," << (*ymm.second) / units::mm << "]mm,\n"
              << "\tz in [" << (*zmm.first) / units::mm << "," << (*zmm.second) / units::mm << "]mm,\n"
              << "\tQtot=" << qtot / units::eplus << ", <Qtot>=" << qtot / ndepos / units::eplus << " "
              << " Qorig=" << qorig / units::eplus << ", <Qorig>=" << qorig / ndepos / units::eplus << " electrons"
              << std::endl;
}

// Implementation warning: this violates node and proc generality in
// order to coerce a join into a linear pipeline.
class NoiseAdderProc : public FilterProc {
   public:
    NoiseAdderProc(IFrameSource::pointer nn)
      : noise_node(nn)
    {
    }
    virtual ~NoiseAdderProc() {}

    virtual Pipe& input_pipe() { return iq; }
    virtual Pipe& output_pipe() { return oq; }

    virtual bool operator()()
    {
        if (iq.empty()) {
            return false;
        }

        const IFrame::pointer iframe = boost::any_cast<const IFrame::pointer>(iq.front());
        iq.pop();

        if (!iframe) {  // eos
            std::cerr << "NoiseAdderProc eos\n";
            boost::any out = iframe;
            oq.push(out);
            return true;
        }

        IFrame::pointer nframe;
        bool ok = (*noise_node)(nframe);
        if (!ok) return false;
        nframe = Aux::sum(IFrame::vector{iframe, nframe}, iframe->ident());
        boost::any anyout = nframe;
        oq.push(anyout);
        return true;
    }

   private:
    Pipe iq, oq;
    IFrameSource::pointer noise_node;
};

void Gen::Fourdee::execute()
{
    execute_new();
    // execute_old();
}

void Gen::Fourdee::execute_new()
{
    if (!m_depos) {
        cerr << "Fourdee: no depos, can't do much" << endl;
        return;
    }

    if (!m_ductor and (m_digitizer or m_dissonance or m_digitizer or m_filter or m_output)) {
        std::cerr << "Fourdee: a Ductor is required for subsequent pipeline stages\n";
        return;
    }

    Pipeline pipeline;
    auto source = new SourceNodeProc(m_depos);
    Proc* last_proc = source;

    if (m_drifter) {  // depo in, depo out
        cerr << "Pipeline adding #" << pipeline.size() << ": " << type(*m_drifter) << endl;
        auto proc = new QueuedNodeProc(m_drifter);
        last_proc = join(pipeline, last_proc, proc);
    }
    if (m_depofilter) {  // depo in, depo out
        cerr << "Pipeline adding #" << pipeline.size() << ": " << type(*m_depofilter) << endl;
        auto proc = new FunctionNodeProc(m_depofilter);
        last_proc = join(pipeline, last_proc, proc);
    }

    if (m_ductor) {  // depo in, zero or more frames out
        cerr << "Pipeline adding #" << pipeline.size() << ": " << type(*m_ductor) << endl;
        auto proc = new QueuedNodeProc(m_ductor);
        last_proc = join(pipeline, last_proc, proc);
    }

    if (m_dissonance) {  // frame in, frame out
        cerr << "Pipeline adding #" << pipeline.size() << ": " << type(*m_dissonance) << endl;
        auto proc = new NoiseAdderProc(m_dissonance);
        last_proc = join(pipeline, last_proc, proc);
    }

    if (m_digitizer) {  // frame in, frame out
        cerr << "Pipeline adding #" << pipeline.size() << ": " << type(*m_digitizer) << endl;
        auto proc = new FunctionNodeProc(m_digitizer);
        last_proc = join(pipeline, last_proc, proc);
    }

    if (m_filter) {  // frame in, frame out
        cerr << "Pipeline adding #" << pipeline.size() << ": " << type(*m_filter) << endl;
        auto proc = new FunctionNodeProc(m_filter);
        last_proc = join(pipeline, last_proc, proc);
    }

    if (m_output) {  // frame in, full stop.
        cerr << "Pipeline adding #" << pipeline.size() << ": " << type(*m_output) << endl;
        auto proc = new SinkNodeProc(m_output);
        last_proc = join(pipeline, last_proc, proc);
    }
    else {
        cerr << "Pipeline adding #" << pipeline.size() << ": sink\n";
        auto sink = new DropSinkProc;
        last_proc = join(pipeline, last_proc, sink);
    }

    // truly last proc needs to be added by hand.
    pipeline.push_back(last_proc);

    // A "drain end first" strategy gives attention to draining
    // queues the more toward the the end of the pipeline they are.
    // This keeps the overall pipeline empty and leads to "pulse" type
    // data movement.  It's stupid for multiprocessing but in a single
    // thread will keep memory usage low.
    size_t pipelen = pipeline.size();
    while (true) {
        bool did_something = false;
        for (size_t ind = pipelen - 1; ind > 0; --ind) {  // source has no input
            SinkProc* proc = dynamic_cast<SinkProc*>(pipeline[ind]);
            if (proc->input_pipe().empty()) {
                continue;
            }
            bool ok = (*proc)();
            if (!ok) {
                std::cerr << "Pipeline failed\n";
                return;
            }
            // std::cerr << "Executed process " << ind << std::endl;
            did_something = true;
            break;
        }
        if (!did_something) {
            bool ok = (*pipeline[0])();
            if (!ok) {
                std::cerr << "Source empty\n";
                return;
            }
            // std::cerr << "Executed source\n";
        }
        // otherwise, go through pipeline again
    }
}
void Gen::Fourdee::execute_old()
{
    if (!m_depos) cerr << "Fourdee: no depos" << endl;
    if (!m_drifter) cerr << "Fourdee: no drifter" << endl;
    if (!m_ductor) cerr << "Fourdee: no ductor" << endl;
    if (!m_dissonance) cerr << "Fourdee: no dissonance" << endl;
    if (!m_digitizer) cerr << "Fourdee: no digitizer" << endl;
    if (!m_filter) cerr << "Fourdee: no filter" << endl;
    if (!m_output) cerr << "Fourdee: no sink" << endl;

    // here we make a manual pipeline.  In a "real" app this might be
    // a DFP executed by TBB.
    int count = 0;
    int ndrifted = 0;
    ExecMon em;
    cerr << "Gen::Fourdee: starting\n";
    while (true) {
        ++count;

        IDepo::pointer depo;
        if (!(*m_depos)(depo)) {
            cerr << "DepoSource is empty\n";
        }
        if (!depo) {
            cerr << "Got null depo from source at loop iteration " << count << endl;
        }

        IDrifter::output_queue drifted;
        if (!(*m_drifter)(depo, drifted)) {
            cerr << "Stopping on " << type(*m_drifter) << endl;
            goto bail;
        }
        if (drifted.empty()) {
            continue;
        }
        if (m_depofilter) {
            IDrifter::output_queue fdrifted;
            for (auto drifted_depo : drifted) {
                IDepo::pointer fdepo;
                if ((*m_depofilter)(drifted_depo, fdepo)) {
                    fdrifted.push_back(fdepo);
                }
            }
            drifted = fdrifted;
        }

        ndrifted += drifted.size();
        cerr << "Gen::FourDee: seen " << ndrifted << " drifted\n";
        dump(drifted);

        for (auto drifted_depo : drifted) {
            IDuctor::output_queue frames;
            if (!(*m_ductor)(drifted_depo, frames)) {
                cerr << "Stopping on " << type(*m_ductor) << endl;
                goto bail;
            }
            if (frames.empty()) {
                continue;
            }
            cerr << "Gen::FourDee: got " << frames.size() << " frames\n";

            for (IFrameFilter::input_pointer voltframe : frames) {
                em("got frame");
                cerr << "voltframe: ";
                dump(voltframe);

                if (m_dissonance) {
                    cerr << "Gen::FourDee: including noise\n";
                    IFrame::pointer noise;
                    if (!(*m_dissonance)(noise)) {
                        cerr << "Stopping on " << type(*m_dissonance) << endl;
                        goto bail;
                    }
                    if (noise) {
                        cerr << "noiseframe: ";
                        dump(noise);
                        voltframe = Aux::sum(IFrame::vector{voltframe, noise}, voltframe->ident());
                        em("got noise");
                    }
                    else {
                        cerr << "Noise source is empty\n";
                    }
                }

                IFrameFilter::output_pointer adcframe;
                if (m_digitizer) {
                    if (!(*m_digitizer)(voltframe, adcframe)) {
                        cerr << "Stopping on " << type(*m_digitizer) << endl;
                        goto bail;
                    }
                    em("digitized");
                    cerr << "digiframe: ";
                    dump(adcframe);
                }
                else {
                    adcframe = voltframe;
                }

                if (m_filter) {
                    IFrameFilter::output_pointer filtframe;
                    if (!(*m_filter)(adcframe, filtframe)) {
                        cerr << "Stopping on " << type(*m_filter) << endl;
                        goto bail;
                    }
                    adcframe = filtframe;
                    em("filtered");
                }

                if (m_output) {
                    if (!(*m_output)(adcframe)) {
                        cerr << "Stopping on " << type(*m_output) << endl;
                        goto bail;
                    }
                    em("output");
                }
                cerr << "Gen::Fourdee: frame with " << adcframe->traces()->size() << " traces\n";
            }
        }
    }
bail:  // what's this weird syntax?  What is this, BASIC?
    cerr << em.summary() << endl;
}
