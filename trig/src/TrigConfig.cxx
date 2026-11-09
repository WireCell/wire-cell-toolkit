#include "WireCellTrig/TrigConfig.h"

#include "WireCellUtil/Exceptions.h"

using namespace WireCell;
using namespace WireCell::Trig;

namespace {
    // Read an unsigned integer through int64_t, keeping the default if absent.
    template <typename T>
    T get_uint(const Configuration& cfg, const std::string& key, T def)
    {
        if (!cfg.isMember(key) || cfg[key].isNull()) return def;
        const int64_t val = cfg[key].asInt64();
        if (val < 0) {
            THROW(ValueError() << errmsg{"trig config: negative value for unsigned key " + key});
        }
        return (T) val;
    }
}  // namespace

TPGConfig WireCell::Trig::tpg_config_from_json(const Configuration& cfg)
{
    TPGConfig tpg;
    tpg.accum_limit = get(cfg, "accum_limit", tpg.accum_limit);

    const std::string init = get<std::string>(cfg, "pedestal_init", "mode");
    if (init == "zero") tpg.pedestal_init = TPGConfig::PedestalInit::zero;
    else if (init == "mode") tpg.pedestal_init = TPGConfig::PedestalInit::mode;
    else if (init == "fixed") tpg.pedestal_init = TPGConfig::PedestalInit::fixed;
    else THROW(ValueError() << errmsg{"trig config: unknown pedestal_init: " + init});
    tpg.pedestal = (int16_t) get(cfg, "pedestal", (int) tpg.pedestal);

    const std::string filter = get<std::string>(cfg, "filter", "none");
    if (filter == "none") tpg.filter = TPGConfig::Filter::none;
    else if (filter == "running_sum") tpg.filter = TPGConfig::Filter::running_sum;
    else if (filter == "abs_running_sum") tpg.filter = TPGConfig::Filter::abs_running_sum;
    else THROW(ValueError() << errmsg{"trig config: unknown filter: " + filter});
    tpg.r_value = get(cfg, "r_value", tpg.r_value);
    tpg.scale_factor = get(cfg, "scale_factor", tpg.scale_factor);
    tpg.saturate_low = get(cfg, "saturate_low", tpg.saturate_low);
    tpg.filter_pedestal = get(cfg, "filter_pedestal", tpg.filter_pedestal);
    tpg.threshold = (int16_t) get(cfg, "threshold", (int) tpg.threshold);
    tpg.min_samples_over_threshold =
        get_uint(cfg, "min_samples_over_threshold", tpg.min_samples_over_threshold);
    return tpg;
}

Configuration WireCell::Trig::tpg_config_to_json(const TPGConfig& tpg)
{
    Configuration cfg;
    cfg["accum_limit"] = tpg.accum_limit;
    const char* inits[] = {"zero", "mode", "fixed"};
    cfg["pedestal_init"] = inits[(int) tpg.pedestal_init];
    cfg["pedestal"] = tpg.pedestal;
    const char* filters[] = {"none", "running_sum", "abs_running_sum"};
    cfg["filter"] = filters[(int) tpg.filter];
    cfg["r_value"] = tpg.r_value;
    cfg["scale_factor"] = tpg.scale_factor;
    cfg["saturate_low"] = tpg.saturate_low;
    cfg["filter_pedestal"] = tpg.filter_pedestal;
    cfg["threshold"] = tpg.threshold;
    cfg["min_samples_over_threshold"] = tpg.min_samples_over_threshold;
    cfg["settle"] = 0;
    return cfg;
}

TAMaker::Config WireCell::Trig::tamaker_config_from_json(const Configuration& cfg)
{
    TAMaker::Config com;
    com.prescale = get_uint(cfg, "prescale", com.prescale);
    com.max_samples_over_threshold = get_uint(cfg, "max_samples_over_threshold", com.max_samples_over_threshold);
    return com;
}

TCMaker::Config WireCell::Trig::tcmaker_config_from_json(const Configuration& cfg)
{
    TCMaker::Config com;
    com.prescale = get_uint(cfg, "prescale", com.prescale);
    com.tc_type = (Daq::TriggerCandidateData::Type) get(cfg, "tc_type", (int) com.tc_type);
    return com;
}

TAMakerSupernova::Config WireCell::Trig::supernova_ta_config_from_json(const Configuration& cfg,
                                                                       double default_ticks_per_sample)
{
    TAMakerSupernova::Config sn;
    sn.time_tolerance = get(cfg, "time_tolerance", (int64_t) sn.time_tolerance);
    sn.channel_tolerance = get(cfg, "channel_tolerance", (int) sn.channel_tolerance);
    sn.both = get(cfg, "both", sn.both);
    sn.ticks_per_sample = get(cfg, "ticks_per_sample", default_ticks_per_sample);
    return sn;
}

TCMakerSupernova::Config WireCell::Trig::supernova_tc_config_from_json(const Configuration& cfg)
{
    TCMakerSupernova::Config sn;
    sn.time_window = get(cfg, "time_window", (int64_t) sn.time_window);
    sn.threshold = get_uint(cfg, "threshold", sn.threshold);
    sn.hit_threshold = get_uint(cfg, "hit_threshold", sn.hit_threshold);
    return sn;
}

TAMakerSWIFT::Config WireCell::Trig::swift_config_from_json(const Configuration& cfg)
{
    TAMakerSWIFT::Config sw;
    sw.window_length = get_uint(cfg, "window_length", sw.window_length);
    sw.inspect_energy_threshold = get_uint(cfg, "inspect_energy_threshold_sadc", sw.inspect_energy_threshold);
    sw.accept_energy_threshold = get_uint(cfg, "accept_energy_threshold_sadc", sw.accept_energy_threshold);
    sw.min_adc_peak = get_uint(cfg, "min_adc_peak", sw.min_adc_peak);
    sw.min_samples_over_threshold = get_uint(cfg, "min_samples_over_threshold", sw.min_samples_over_threshold);
    sw.cm_per_tick = get(cfg, "cm_per_tick", sw.cm_per_tick);
    sw.wire_pitch = get(cfg, "wire_pitch", sw.wire_pitch);
    sw.min_samples = get(cfg, "min_samples", sw.min_samples);
    sw.epsilon = get(cfg, "epsilon", sw.epsilon);
    sw.cluster_energy_cut = get_uint(cfg, "cluster_energy_cut_sadc", sw.cluster_energy_cut);
    if (sw.window_length == 0) {
        THROW(ValueError() << errmsg{"trig config: SWIFT window_length must be positive"});
    }
    return sw;
}

TCMakerBurst::Config WireCell::Trig::burst_config_from_json(const Configuration& cfg)
{
    TCMakerBurst::Config bc;
    bc.window = get(cfg, "window", bc.window);
    bc.threshold = get_uint(cfg, "threshold", bc.threshold);
    bc.holdoff = get(cfg, "holdoff", bc.holdoff);
    bc.tc_type = get(cfg, "tc_type", bc.tc_type);
    bc.algorithm = get(cfg, "algorithm", bc.algorithm);
    auto& sel = bc.select;
    if (cfg.isMember("algorithms") && !cfg["algorithms"].isNull()) {
        sel.algorithms.clear();
        for (const auto& one : cfg["algorithms"]) sel.algorithms.push_back(one.asInt());
    }
    sel.min_inputs = get_uint(cfg, "min_inputs", sel.min_inputs);
    sel.max_inputs = get_uint(cfg, "max_inputs", sel.max_inputs);
    sel.min_adc_integral = get_uint(cfg, "min_adc_integral", sel.min_adc_integral);
    sel.max_adc_integral = get_uint(cfg, "max_adc_integral", sel.max_adc_integral);
    sel.max_channel_span = get(cfg, "max_channel_span", sel.max_channel_span);
    sel.max_time_span = get(cfg, "max_time_span", sel.max_time_span);
    if (bc.window <= 0) {
        THROW(ValueError() << errmsg{"trig config: burst window must be positive"});
    }
    if (bc.threshold == 0) {
        THROW(ValueError() << errmsg{"trig config: burst threshold must be positive"});
    }
    return bc;
}
