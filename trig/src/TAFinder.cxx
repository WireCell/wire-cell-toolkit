#include "WireCellTrig/TAFinder.h"
#include "WireCellTrig/TrigConfig.h"

#include "WireCellAux/SimpleTriggerSets.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>

WIRECELL_FACTORY(TAFinder, WireCell::Trig::TAFinder, WireCell::INamed, WireCell::ITriggerActivityFinder,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Trig;

TAFinder::TAFinder()
  : Aux::Logger("TAFinder", "trig")
  , m_period(500 * units::ns)
{
    m_clock.tick = 16 * units::ns;
}

TAFinder::~TAFinder() {}

WireCell::Configuration TAFinder::default_configuration() const
{
    Configuration cfg;
    cfg["algorithm"] = m_algorithm;
    cfg["config"] = Json::objectValue;
    cfg["planes"] = Json::arrayValue;
    cfg["prescale"] = 1;
    cfg["store_inputs"] = m_store_inputs;
    cfg["anode"] = "";
    cfg["period"] = m_period;
    cfg["clock_tick"] = m_clock.tick;
    cfg["epoch"] = m_clock.epoch;
    cfg["latency"] = m_latency;
    cfg["streaming"] = m_streaming;
    return cfg;
}

void TAFinder::configure(const WireCell::Configuration& cfg)
{
    m_algorithm = get(cfg, "algorithm", m_algorithm);
    if (m_algorithm != "Supernova" && m_algorithm != "SWIFT") {
        THROW(ValueError() << errmsg{"TAFinder: unknown algorithm: " + m_algorithm});
    }
    m_alg_cfg = cfg["config"];
    m_plane_cfgs.clear();
    for (const auto& one : cfg["planes"]) {
        Configuration merged = m_alg_cfg;
        update(merged, one);
        m_plane_cfgs.push_back(merged);
    }
    m_common = tamaker_config_from_json(cfg);
    m_common.store_inputs = true;  // the node strips, keeping the count
    m_store_inputs = get(cfg, "store_inputs", m_store_inputs);
    m_period = get(cfg, "period", m_period);
    m_clock.tick = get(cfg, "clock_tick", m_clock.tick);
    m_clock.epoch = get(cfg, "epoch", m_clock.epoch);
    m_latency = get(cfg, "latency", m_latency);
    m_streaming = get(cfg, "streaming", m_streaming);
    const std::string anode_tn = get<std::string>(cfg, "anode", "");
    m_anode = anode_tn.empty() ? nullptr : Factory::find_tn<IAnodePlane>(anode_tn);
    if (!m_plane_cfgs.empty() && !m_anode) {
        THROW(ValueError() << errmsg{"TAFinder: per-plane configuration requires an anode"});
    }
    m_cores.clear();
    m_have_last = false;
    log->debug("algorithm={} store_inputs={} anode={} latency={}", m_algorithm, m_store_inputs,
               anode_tn, m_latency);
}

TAFinder::Core& TAFinder::core(int channel)
{
    WirePlaneId wpid(kUnknownLayer);
    if (m_anode) {
        wpid = m_anode->resolve(channel);
        if (!wpid.valid()) {
            THROW(ValueError() << errmsg{"TAFinder: channel " + std::to_string(channel) + " not in anode"});
        }
    }
    const int key = wpid.ident();
    auto it = m_cores.find(key);
    if (it != m_cores.end()) return it->second;

    Configuration acfg = m_alg_cfg;
    if (!m_plane_cfgs.empty()) {
        const size_t index = wpid.index();
        if (index >= m_plane_cfgs.size()) {
            THROW(ValueError() << errmsg{"TAFinder: no configuration for plane index " + std::to_string(index)});
        }
        acfg = m_plane_cfgs[index];
    }
    Core one;
    one.wpid = wpid;
    if (m_algorithm == "Supernova") {
        one.maker = std::make_unique<TAMakerSupernova>(supernova_ta_config_from_json(acfg, m_period / m_clock.tick));
    }
    else {
        one.maker = std::make_unique<TAMakerSWIFT>(swift_config_from_json(acfg));
    }
    one.maker->set_common(m_common);
    return m_cores.emplace(key, std::move(one)).first->second;
}

void TAFinder::add_tas(const Core& core, std::vector<Daq::TriggerActivity>& daq_tas,
                       std::vector<TriggerActivity>& tas) const
{
    for (const auto& dta : daq_tas) {
        auto ta = from_daq(dta, m_clock);
        ta.detector = m_detector;
        ta.plane = core.wpid;
        if (!m_store_inputs) ta.inputs.clear();
        tas.push_back(std::move(ta));
    }
    daq_tas.clear();
}

void TAFinder::flush_cores(Daq::timestamp_t until, std::vector<TriggerActivity>& tas)
{
    std::vector<Daq::TriggerActivity> daq_tas;
    for (auto& [key, one] : m_cores) {
        one.maker->flush(until, daq_tas);
        add_tas(one, daq_tas, tas);
    }
}

static void sort_tas(std::vector<TriggerActivity>& tas)
{
    std::stable_sort(tas.begin(), tas.end(),
                     [](const auto& a, const auto& b) { return a.time_start < b.time_start; });
}

bool TAFinder::operator()(const input_pointer& in, output_queue& outq)
{
    if (!in) {  // end of stream
        if (m_have_last) {
            std::vector<TriggerActivity> tas;
            flush_cores(Daq::invalid_timestamp, tas);
            sort_tas(tas);
            log->debug("end of stream, {} TAs flushed", tas.size());
            outq.push_back(std::make_shared<Aux::SimpleTriggerActivitySet>(m_count++, m_detector, m_last_end,
                                                                           m_last_end, std::move(tas)));
        }
        m_cores.clear();
        m_have_last = false;
        outq.push_back(nullptr);
        return true;
    }

    std::vector<TriggerActivity> tas;
    if (m_have_last && std::abs(in->time_begin() - m_last_end) > 1e-3 * m_clock.tick) {
        log->debug("coverage gap before set {}, restarting", in->ident());
        flush_cores(Daq::invalid_timestamp, tas);
        m_cores.clear();
    }
    m_detector = in->detector();

    // The set is ordered by time; feed each TP to its plane's core.
    std::vector<Daq::TriggerActivity> daq_tas;
    for (const auto& tp : in->primitives()) {
        auto& one = core(tp.channel);
        one.maker->operator()(to_daq(tp, m_clock, m_period), daq_tas);
        add_tas(one, daq_tas, tas);
    }

    if (!m_streaming) {
        flush_cores(Daq::invalid_timestamp, tas);
        m_cores.clear();
    }
    else {
        const double until = in->time_end() - m_latency;
        if (until > m_clock.epoch) {
            flush_cores(m_clock.to_ticks(until), tas);
        }
    }

    sort_tas(tas);
    log->debug("set {}: {} TPs, {} TAs", in->ident(), in->primitives().size(), tas.size());
    outq.push_back(std::make_shared<Aux::SimpleTriggerActivitySet>(m_count++, m_detector, in->time_begin(),
                                                                   in->time_end(), std::move(tas)));
    m_last_end = in->time_end();
    m_have_last = m_streaming;
    return true;
}
