#include "WireCellTrig/TCFinder.h"
#include "WireCellTrig/TrigConfig.h"

#include "WireCellAux/SimpleTriggerSets.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>

WIRECELL_FACTORY(TCFinder, WireCell::Trig::TCFinder, WireCell::INamed, WireCell::ITriggerCandidateFinder,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Trig;

TCFinder::TCFinder()
  : Aux::Logger("TCFinder", "trig")
  , m_period(500 * units::ns)
{
    m_clock.tick = 16 * units::ns;
    m_common.tc_type = Daq::TriggerCandidateData::Type::kSupernova;
}

TCFinder::~TCFinder() {}

WireCell::Configuration TCFinder::default_configuration() const
{
    Configuration cfg;
    cfg["algorithm"] = m_algorithm;
    cfg["config"] = Json::objectValue;
    cfg["store_inputs"] = m_store_inputs;
    cfg["lateness"] = m_lateness;
    cfg["prescale"] = 1;
    cfg["tc_type"] = (int) m_common.tc_type;
    cfg["period"] = m_period;
    cfg["clock_tick"] = m_clock.tick;
    cfg["epoch"] = m_clock.epoch;
    return cfg;
}

void TCFinder::configure(const WireCell::Configuration& cfg)
{
    m_algorithm = get(cfg, "algorithm", m_algorithm);
    if (m_algorithm != "Burst" && m_algorithm != "Supernova") {
        THROW(ValueError() << errmsg{"TCFinder: unknown algorithm: " + m_algorithm});
    }
    m_alg_cfg = cfg["config"];
    m_common = tcmaker_config_from_json(cfg);
    m_common.store_inputs = true;  // the node strips, keeping the count
    m_store_inputs = get(cfg, "store_inputs", m_store_inputs);
    m_lateness = get(cfg, "lateness", m_lateness);
    m_period = get(cfg, "period", m_period);
    m_clock.tick = get(cfg, "clock_tick", m_clock.tick);
    m_clock.epoch = get(cfg, "epoch", m_clock.epoch);
    make_core();  // also checks the algorithm configuration
    m_have_last = false;
    log->debug("algorithm={} store_inputs={} lateness={}", m_algorithm, m_store_inputs, m_lateness);
}

void TCFinder::make_core()
{
    m_burst.reset();
    m_daq.reset();
    if (m_algorithm == "Burst") {
        m_burst = std::make_unique<TCMakerBurst>(burst_config_from_json(m_alg_cfg));
        return;
    }
    m_daq = std::make_unique<TCMakerSupernova>(supernova_tc_config_from_json(m_alg_cfg));
    m_daq->set_common(m_common);
}

void TCFinder::feed(const TriggerActivity& ta, std::vector<TriggerCandidate>& tcs)
{
    if (m_burst) {
        (*m_burst)(ta, tcs);
        return;
    }
    std::vector<Daq::TriggerCandidate> dtcs;
    (*m_daq)(to_daq(ta, m_clock, m_period), dtcs);
    for (const auto& dtc : dtcs) tcs.push_back(from_daq(dtc, m_clock));
}

void TCFinder::flush(double until, std::vector<TriggerCandidate>& tcs)
{
    if (m_burst) {
        m_burst->flush(until, tcs);
        return;
    }
    const auto ticks = std::isinf(until) ? Daq::invalid_timestamp : m_clock.to_ticks(until);
    std::vector<Daq::TriggerCandidate> dtcs;
    m_daq->flush(ticks, dtcs);
    for (const auto& dtc : dtcs) tcs.push_back(from_daq(dtc, m_clock));
}

static void finish(std::vector<TriggerCandidate>& tcs, bool store_inputs)
{
    if (!store_inputs) {
        for (auto& tc : tcs) tc.inputs.clear();
    }
    std::stable_sort(tcs.begin(), tcs.end(),
                     [](const auto& a, const auto& b) { return a.time_start < b.time_start; });
}

bool TCFinder::operator()(const input_pointer& in, output_queue& outq)
{
    const double inf = std::numeric_limits<double>::infinity();
    if (!in) {  // end of stream
        if (m_have_last) {
            std::vector<TriggerCandidate> tcs;
            flush(inf, tcs);
            finish(tcs, m_store_inputs);
            log->debug("end of stream, {} TCs flushed", tcs.size());
            outq.push_back(std::make_shared<Aux::SimpleTriggerCandidateSet>(m_count++, m_detector, m_last_end,
                                                                            m_last_end, std::move(tcs)));
        }
        if (m_burst) {
            log->debug("burst: {} TAs seen, {} selected, {} late, {} dropped", m_burst->nseen(),
                       m_burst->nselected(), m_burst->nlate(), m_burst->ndropped());
        }
        make_core();
        m_have_last = false;
        outq.push_back(nullptr);
        return true;
    }

    std::vector<TriggerCandidate> tcs;
    if (m_have_last && std::abs(in->time_begin() - m_last_end) > 1e-3 * m_clock.tick) {
        log->debug("coverage gap before set {}, restarting", in->ident());
        flush(inf, tcs);
        make_core();
    }
    m_detector = in->detector();

    for (const auto& ta : in->activities()) feed(ta, tcs);

    const double until = in->time_end() - m_lateness;
    if (m_burst || until > m_clock.epoch) flush(until, tcs);

    finish(tcs, m_store_inputs);
    log->debug("set {}: {} TAs, {} TCs", in->ident(), in->activities().size(), tcs.size());
    outq.push_back(std::make_shared<Aux::SimpleTriggerCandidateSet>(m_count++, m_detector, in->time_begin(),
                                                                    in->time_end(), std::move(tcs)));
    m_last_end = in->time_end();
    m_have_last = true;
    return true;
}
