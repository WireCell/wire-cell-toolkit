#include "WireCellGen/TrackSegmentSampler.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <cmath>

WIRECELL_FACTORY(TrackSegmentSampler, WireCell::Gen::TrackSegmentSampler, WireCell::INamed,
                 WireCell::ITrackSegmentSampler, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Gen;

TrackSegmentSampler::TrackSegmentSampler()
  : Aux::Logger("TrackSegmentSampler", "gen")
  , m_step_size(1.0 * units::mm)
  , m_ionization("quanta")
  , m_recombination("BoxRecombination")
{
}
TrackSegmentSampler::~TrackSegmentSampler() {}

WireCell::Configuration TrackSegmentSampler::default_configuration() const
{
    Configuration cfg;
    cfg["ionization"] = m_ionization;
    cfg["recombination"] = m_recombination;
    cfg["step_size"] = m_step_size;
    return cfg;
}

void TrackSegmentSampler::configure(const WireCell::Configuration& cfg)
{
    m_step_size = get<double>(cfg, "step_size", m_step_size);
    if (m_step_size <= 0) {
        THROW(ValueError() << errmsg{"TrackSegmentSampler: step_size must be positive"});
    }

    m_ionization = get<std::string>(cfg, "ionization", m_ionization);
    m_recombination = get<std::string>(cfg, "recombination", m_recombination);

    // Select the ionization model ONCE; operator() has no mode dispatch.
    if (m_ionization == "quanta") {
        m_ionize = [](const ITrackSegment& seg) -> double {
            const double ne = seg.n_electrons();
            // ITrackSegment reports a negative n_electrons() to mean "not
            // provided" (producers set a -1 sentinel).  Treat a clearly
            // negative value as that error, but clamp small sub-count
            // negatives (eg floating-point noise around an exact-zero yield)
            // up to zero so a vanishing ionization count does not abort a job.
            if (ne < -0.5) {
                THROW(ValueError() << errmsg{"TrackSegmentSampler: segment lacks n_electrons "
                                             "required by ionization mode 'quanta'"});
            }
            return -std::max(0.0, ne);  // electrons are negative charge
        };
        log->debug("ionization mode 'quanta' (producer n_electrons, no recombination model)");
        return;
    }
    if (m_ionization == "recombination") {
        auto model = Factory::find_tn<IRecombinationModel>(m_recombination);
        const double step = m_step_size;
        m_ionize = [model, step](const ITrackSegment& seg) -> double {
            // dE/dX from the full charged path length; fall back to the
            // straight-line length, then to one step, for degenerate input.
            double dX = seg.track_length();
            if (dX <= 0) {
                dX = (seg.stop() - seg.start()).magnitude();
            }
            if (dX <= 0) {
                dX = step;
            }
            return (*model)(seg.energy(), dX);
        };
        log->debug("ionization mode 'recombination' via '{}'", m_recombination);
        return;
    }
    THROW(ValueError() << errmsg{"TrackSegmentSampler: unknown ionization mode '" + m_ionization +
                                 "' (expected 'recombination' or 'quanta')"});
}

bool TrackSegmentSampler::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }
    if (!m_ionize) {
        THROW(ValueError() << errmsg{"TrackSegmentSampler: not configured"});
    }

    IDepo::vector depos;
    auto segments = in->segments();
    if (segments) {
        for (const auto& seg : *segments) {
            const Point step_vec = seg->stop() - seg->start();
            const double dist = step_vec.magnitude();
            const int nsamples = std::max(1, static_cast<int>(std::ceil(dist / m_step_size)));

            const double charge = m_ionize(*seg) / nsamples;
            const double energy = seg->energy() / nsamples;
            const double dtime = seg->stop_time() - seg->start_time();

            for (int ind = 0; ind < nsamples; ++ind) {
                const double frac = (ind + 0.5) / nsamples;
                const Point pos = seg->start() + frac * step_vec;
                const double time = seg->start_time() + frac * dtime;
                depos.push_back(std::make_shared<Aux::SimpleDepo>(time, pos, charge, nullptr, 0.0,
                                                                  0.0, seg->id(), seg->pdg(),
                                                                  energy));
            }
        }
    }

    log->debug("call={} segments={} depos={}", m_count, segments ? segments->size() : 0,
               depos.size());
    ++m_count;
    out = std::make_shared<Aux::SimpleDepoSet>(in->ident(), depos);
    return true;
}
