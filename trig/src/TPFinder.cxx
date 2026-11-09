#include "WireCellTrig/TPFinder.h"
#include "WireCellTrig/TrigConfig.h"

#include "WireCellAux/SimpleTriggerSets.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>

WIRECELL_FACTORY(TPFinder, WireCell::Trig::TPFinder, WireCell::INamed, WireCell::ITriggerPrimitiveFinder,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Trig;

TPFinder::TPFinder()
  : Aux::Logger("TPFinder", "trig")
  , m_period(500 * units::ns)
{
}

TPFinder::~TPFinder() {}

WireCell::Configuration TPFinder::default_configuration() const
{
    Configuration cfg;
    cfg["period"] = m_period;
    cfg["anode"] = "";
    cfg["detector"] = m_detector;
    cfg["tag"] = m_tag;
    cfg["streaming"] = m_streaming;
    cfg["tpg"] = tpg_config_to_json(TPGConfig{});
    cfg["planes"] = Json::arrayValue;
    return cfg;
}

void TPFinder::configure(const WireCell::Configuration& cfg)
{
    m_period = get(cfg, "period", m_period);
    m_tag = get(cfg, "tag", m_tag);
    m_streaming = get(cfg, "streaming", m_streaming);
    m_detector = get(cfg, "detector", m_detector);
    const std::string anode_tn = get<std::string>(cfg, "anode", "");
    if (anode_tn.empty()) {
        m_anode = nullptr;
    }
    else {
        m_anode = Factory::find_tn<IAnodePlane>(anode_tn);
        m_detector = m_anode->ident();
    }

    const auto& tpg = cfg["tpg"];
    m_settle = get(tpg, "settle", 0);
    m_plane_cfgs.clear();
    const auto& planes = cfg["planes"];
    if (planes.isArray() && planes.size() > 0) {
        for (const auto& one : planes) {
            Configuration merged = tpg;
            update(merged, one);
            m_plane_cfgs.push_back(tpg_config_from_json(merged));
        }
    }
    else {
        m_plane_cfgs.push_back(tpg_config_from_json(tpg));
    }
    if (m_plane_cfgs.size() > 1 && !m_anode) {
        THROW(ValueError() << errmsg{"TPFinder: per-plane configuration requires an anode"});
    }
    m_channels.clear();
    m_have_last = false;
    log->debug("period={} streaming={} detector={} with {} plane configurations", m_period, m_streaming,
               m_detector, m_plane_cfgs.size());
}

const TPGConfig& TPFinder::plane_config(int channel)
{
    if (m_plane_cfgs.size() == 1) return m_plane_cfgs[0];
    const auto wpid = m_anode->resolve(channel);
    if (!wpid.valid()) {
        THROW(ValueError() << errmsg{"TPFinder: channel " + std::to_string(channel) + " not in anode"});
    }
    const size_t index = wpid.index();
    if (index >= m_plane_cfgs.size()) {
        THROW(ValueError() << errmsg{"TPFinder: no configuration for plane index " + std::to_string(index)});
    }
    return m_plane_cfgs[index];
}

TriggerPrimitive TPFinder::to_tp(int channel, const ThresholdHit& hit) const
{
    TriggerPrimitive tp;
    tp.channel = channel;
    tp.time_start = hit.sample_start * m_period + m_resid;
    tp.samples_over_threshold = hit.samples_over_threshold;
    tp.samples_to_peak = hit.samples_to_peak;
    tp.adc_integral = hit.adc_integral;
    tp.adc_peak = hit.adc_peak;
    return tp;
}

void TPFinder::flush_all(ITriggerPrimitiveSet::primitives_t& tps)
{
    std::vector<ThresholdHit> hits;
    for (auto& [channel, chan] : m_channels) {
        hits.clear();
        chan.tpg.flush(chan.next, hits);
        for (const auto& hit : hits) tps.push_back(to_tp(channel, hit));
    }
    m_channels.clear();
}

static void sort_tps(ITriggerPrimitiveSet::primitives_t& tps)
{
    std::sort(tps.begin(), tps.end(), [](const auto& a, const auto& b) {
        if (a.time_start == b.time_start) return a.channel < b.channel;
        return a.time_start < b.time_start;
    });
}

static int16_t to_adc(float charge)
{
    const long val = std::lround(charge);
    return (int16_t) std::clamp<long>(val, std::numeric_limits<int16_t>::min(), std::numeric_limits<int16_t>::max());
}

bool TPFinder::operator()(const input_pointer& in, output_queue& outq)
{
    if (!in) {  // end of stream
        if (m_streaming && m_have_last) {
            ITriggerPrimitiveSet::primitives_t tps;
            flush_all(tps);
            sort_tps(tps);
            log->debug("end of stream, {} open TPs flushed", tps.size());
            outq.push_back(std::make_shared<Aux::SimpleTriggerPrimitiveSet>(m_count++, m_detector, m_last_end,
                                                                            m_last_end, std::move(tps)));
        }
        m_channels.clear();
        m_have_last = false;
        outq.push_back(nullptr);
        return true;
    }

    if (std::abs(in->tick() - m_period) > 1e-6 * m_period) {
        THROW(ValueError() << errmsg{"TPFinder: frame tick " + std::to_string(in->tick()) +
                                     " differs from configured period " + std::to_string(m_period)});
    }

    ITriggerPrimitiveSet::primitives_t tps;

    const double t0 = in->time();
    const int64_t k0 = std::llround(t0 / m_period);
    const double resid = t0 - k0 * m_period;
    if (m_streaming && m_have_last && std::abs(resid - m_resid) > 1e-3 * m_period) {
        log->warn("frame {} is off the sample grid of the stream, restarting", in->ident());
        flush_all(tps);
    }
    m_resid = resid;

    auto traces = in->traces();
    std::vector<size_t> indices;
    if (m_tag.empty()) {
        indices.resize(traces->size());
        for (size_t ind = 0; ind < indices.size(); ++ind) indices[ind] = ind;
    }
    else {
        indices = in->tagged_traces(m_tag);
    }

    int64_t cov_first = std::numeric_limits<int64_t>::max();
    int64_t cov_last = std::numeric_limits<int64_t>::min();
    std::vector<int16_t> samples;
    std::vector<ThresholdHit> hits;
    for (size_t ind : indices) {
        const auto& trace = traces->at(ind);
        const int channel = trace->channel();
        const auto& charge = trace->charge();
        const int64_t first = k0 + trace->tbin();
        const int64_t next = first + (int64_t) charge.size();
        cov_first = std::min(cov_first, first);
        cov_last = std::max(cov_last, next);

        samples.resize(charge.size());
        std::transform(charge.begin(), charge.end(), samples.begin(), to_adc);
        hits.clear();

        if (!m_streaming) {
            TPGChannel tpg(plane_config(channel));
            if (m_settle) tpg.settle(samples, m_settle);
            tpg(samples, first, hits);
            tpg.flush(next, hits);
        }
        else {
            auto it = m_channels.find(channel);
            if (it != m_channels.end() && it->second.next != first) {
                it->second.tpg.flush(it->second.next, hits);
                m_channels.erase(it);
                it = m_channels.end();
            }
            if (it == m_channels.end()) {
                it = m_channels.emplace(channel, Channel{TPGChannel(plane_config(channel)), first}).first;
                if (m_settle) it->second.tpg.settle(samples, m_settle);
            }
            it->second.tpg(samples, first, hits);
            it->second.next = next;
        }
        for (const auto& hit : hits) tps.push_back(to_tp(channel, hit));
    }

    double begin = t0, end = t0;
    if (cov_first <= cov_last) {
        begin = cov_first * m_period + resid;
        end = cov_last * m_period + resid;
    }
    sort_tps(tps);
    log->debug("frame {}: {} traces, {} TPs, coverage [{}, {}) us", in->ident(), indices.size(), tps.size(),
               begin / units::us, end / units::us);
    outq.push_back(
        std::make_shared<Aux::SimpleTriggerPrimitiveSet>(m_count++, m_detector, begin, end, std::move(tps)));
    m_last_end = end;
    m_have_last = true;
    return true;
}
