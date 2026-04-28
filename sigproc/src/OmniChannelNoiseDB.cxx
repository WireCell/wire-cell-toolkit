#include "WireCellSigProc/OmniChannelNoiseDB.h"
#include "WireCellAux/DftTools.h"
#include "WireCellUtil/Response.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Point.h"

#include <algorithm>
#include <cmath>

WIRECELL_FACTORY(OmniChannelNoiseDB, WireCell::SigProc::OmniChannelNoiseDB, WireCell::IChannelNoiseDatabase,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SigProc;
using WireCell::Aux::DftTools::fwd_r2c;

OmniChannelNoiseDB::OmniChannelNoiseDB()
  : m_tick(0.5 * units::us)
  , m_nsamples(9600)
  , m_rc_layers(2)
  , log(Log::logger("sigproc"))
{
}
OmniChannelNoiseDB::~OmniChannelNoiseDB()
{
    // for (auto it = m_db.begin(); it!= m_db.end(); it++){
    //    delete it->second;
    // }
    // m_db.clear();
}

void OmniChannelNoiseDB::cache_wire_lengths()
{
    if (m_wire_length_cached) return;
    for (int ch : m_anode->channels()) {
        double total = 0.0;
        for (auto wire : m_anode->wires(ch)) {
            total += ray_length(wire->ray());
        }
        m_wire_length_cm[ch] = total / units::cm;
    }
    m_wire_length_cached = true;
    log->debug("OmniChannelNoiseDB: wire-length cache built for {} channels", m_wire_length_cm.size());
}

OmniChannelNoiseDB::ChannelInfo::ChannelInfo()
  : chid(-1)
  , nominal_baseline(0.0)
  , gain_correction(1.0)
  , response_offset(0.0)
  , min_rms_cut(0.5)
  , max_rms_cut(10.0)
  , pad_window_front(0.0)
  , pad_window_back(0.0)
  , decon_limit(0.02)
  , decon_lf_cutoff(0.08)
  , adc_limit(0.0)
  , decon_limit1(0.08)
  , protection_factor(5.0)
  , min_adc_limit(50)
  , roi_min_max_ratio(0.8)
  , rcrc(nullptr)
  , config(nullptr)
  , noise(nullptr)
  , response(nullptr)
{
}

WireCell::Configuration OmniChannelNoiseDB::default_configuration() const
{
    Configuration cfg;
    // The assumed time-domain sample period used to make response spetra.
    cfg["tick"] = m_tick;
    // The number of *frequency domain* sample.  This is NOT anything
    // related to the length of a waveform to which this class may be
    // applied, except by accident or contrivance.
    cfg["nsamples"] = m_nsamples;
    cfg["anode"] = "AnodePlane";
    cfg["field_response"] = "FieldResponse";

    /// These must be provided
    cfg["groups"] = Json::arrayValue;
    cfg["channel_info"] = Json::arrayValue;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use

    return cfg;
}

/*
  Interpret and return a list of channels for JSON like:

  // just one channel
  channels: 42,

  or

  // explicit list of channels
  channels: [1,42,107],

  or

  // inclusive range of channels
  channels: { first: 0, last: 2400 },

  or

  // all channels in a wire plane
  channels: { wpid: wc.WirePlaneId(wc.kWlayer) },
*/
std::vector<int> OmniChannelNoiseDB::parse_channels(const Json::Value& jchannels)
{
    std::vector<int> ret;

    // single channel
    if (jchannels.isInt()) {
        ret.push_back(jchannels.asInt());
        return ret;
    }

    // array of explicit channels
    if (jchannels.isArray()) {
        const int nch = jchannels.size();
        ret.resize(nch);
        for (int ind = 0; ind < nch; ++ind) {
            ret[ind] = jchannels[ind].asInt();
        }
        return ret;
    }

    // else, assume an object

    // range
    if (jchannels.isMember("first") && jchannels.isMember("last")) {
        const int chf = jchannels["first"].asInt();
        const int chl = jchannels["last"].asInt();
        const int nch = chl - chf + 1;
        ret.resize(nch);
        for (int ind = 0; ind < nch; ++ind) {
            ret[ind] = chf + ind;
        }
        return ret;
    }

    // wire plane id
    if (jchannels.isMember("wpid")) {
        WirePlaneId wpid(jchannels["wpid"].asInt());
        for (auto ch : m_anode->channels()) {
            if (m_anode->resolve(ch) == wpid) {
                ret.push_back(ch);
            }
        }
        return ret;
    }

    return ret;
}

OmniChannelNoiseDB::shared_filter_t OmniChannelNoiseDB::make_filter(std::complex<float> defval)
{
    return std::make_shared<filter_t>(m_nsamples, defval);
}
OmniChannelNoiseDB::shared_filter_t OmniChannelNoiseDB::default_filter()
{
    // Not cached: must reflect the current m_nsamples (which may change via set_nsamples).
    return make_filter();
}

OmniChannelNoiseDB::shared_filter_t OmniChannelNoiseDB::parse_freqmasks(Json::Value jfm)
{
    if (jfm.isNull()) {
        return default_filter();
    }

    auto spectrum = make_filter(std::complex<float>(1, 0));
    for (auto jone : jfm) {
        double value = jone["value"].asDouble();
        int lo = std::max(jone["lobin"].asInt(), 0);
        int hi = std::min(jone["hibin"].asInt(), m_nsamples - 1);
        // std::cerr << "freqmasks: set [" << lo << "," << hi << "] to " << value << std::endl;
        for (int ind = lo; ind <= hi; ++ind) {  // inclusive
            spectrum->at(ind) = value;
        }
    }
    return spectrum;
}

OmniChannelNoiseDB::shared_filter_t OmniChannelNoiseDB::parse_rcrc(Json::Value jrcrc, int nrc)
{
    if (jrcrc.isNull()) {
        return default_filter();
    }
    const double rcrc_val = jrcrc.asDouble();
    const int key = int(round(1000 * rcrc_val / units::ms));
    auto it = m_rcrc_cache.find(key);
    if (it != m_rcrc_cache.end()) {
        return it->second;
    }

    Response::SimpleRC rcres(rcrc_val, m_tick);
    // auto signal = rcres.generate(WireCell::Binning(m_nsamples, 0, m_nsamples*m_tick));
    auto signal = rcres.generate(WireCell::Waveform::Domain(0, m_nsamples * m_tick), m_nsamples);

    Waveform::compseq_t spectrum = fwd_r2c(m_dft, signal);

    // get the square of it because there are two RC filters
    Waveform::compseq_t spectrum2 = spectrum;
    // Waveform::scale(spectrum2,spectrum);

    // std::cerr << "[wgu] parse_rcrc nrc= " << nrc << std::endl;
    nrc--;
    while (nrc > 0) {
        // std::cerr << "[wgu] more nrc= " << nrc << std::endl;
        Waveform::scale(spectrum2, spectrum);
        nrc--;
    }

    // std::cerr << "OmniChannelNoiseDB:: get rcrc as: " << rcrc_val
    //           << " sum=" << Waveform::sum(spectrum2)
    //           << std::endl;

    auto ret = std::make_shared<filter_t>(spectrum2);
    m_rcrc_cache[key] = ret;
    return ret;
}

double OmniChannelNoiseDB::parse_gain(Json::Value jreconfig)
{
    if (jreconfig.empty()) {
        return 1.0;
    }

    const double from_gain = jreconfig["from"]["gain"].asDouble();
    const double to_gain = jreconfig["to"]["gain"].asDouble();
    return to_gain / from_gain;
}

OmniChannelNoiseDB::shared_filter_t OmniChannelNoiseDB::parse_reconfig(Json::Value jreconfig)
{
    if (jreconfig.empty()) {
        return default_filter();
    }

    const double from_gain = jreconfig["from"]["gain"].asDouble();
    const double from_shaping = jreconfig["from"]["shaping"].asDouble();
    const double to_gain = jreconfig["to"]["gain"].asDouble();
    const double to_shaping = jreconfig["to"]["shaping"].asDouble();

    return get_reconfig(from_gain, from_shaping, to_gain, to_shaping);
}
OmniChannelNoiseDB::shared_filter_t OmniChannelNoiseDB::get_reconfig(double from_gain, double from_shaping,
                                                                     double to_gain, double to_shaping)
{
    // kind of evil.
    int key = int(round(10.0 * from_gain / (units::mV / units::fC))) << 24 |
              int(round(10.0 * from_shaping / units::us)) << 16 |
              int(round(10.0 * to_gain / (units::mV / units::fC))) << 8 | int(round(10.0 * to_shaping / units::us));

    // std::cerr << "KEY:" << key
    //           << " fg="<<from_gain/(units::mV/units::fC) << " mV/fC,"
    //           << " fs=" << from_shaping/units::us << " us,"
    //           << " tg="<<to_gain/(units::mV/units::fC) << " mV/fC,"
    //           << " ts="<<to_shaping/units::us << " us,"
    //           << " m_tick=" << m_tick/units::us << " us."
    //           << std::endl;

    auto it = m_reconfig_cache.find(key);
    if (it != m_reconfig_cache.end()) {
        return it->second;
    }

    Response::ColdElec from_ce(from_gain, from_shaping);
    Response::ColdElec to_ce(to_gain, to_shaping);
    // auto to_sig   =   to_ce.generate(WireCell::Binning(m_nsamples, 0, m_nsamples*m_tick));
    // auto from_sig = from_ce.generate(WireCell::Binning(m_nsamples, 0, m_nsamples*m_tick));
    auto to_sig = to_ce.generate(WireCell::Waveform::Domain(0, m_nsamples * m_tick), m_nsamples);
    auto from_sig = from_ce.generate(WireCell::Waveform::Domain(0, m_nsamples * m_tick), m_nsamples);

    auto to_filt = fwd_r2c(m_dft, to_sig);
    auto from_filt = fwd_r2c(m_dft, from_sig);

    // auto from_filt_sum = Waveform::sum(from_filt);
    // auto to_filt_sum   = Waveform::sum(to_filt);

    Waveform::shrink(to_filt, from_filt);  // divide
    auto filt = std::make_shared<filter_t>(to_filt);

    //    std::cerr << filt->at(0) << " " << filt->at(1) << std::endl;

    // std::cerr << "OmniChannelNoiseDB: "
    //           << " from_sig sum=" << Waveform::sum(from_sig)
    //           << " to_sig sum=" << Waveform::sum(to_sig)
    //           << " from_filt sum=" << from_filt_sum
    //           << " to_filt sum=" << to_filt_sum
    //           << " rat_filt sum=" << Waveform::sum(to_filt)
    //           << std::endl;

    m_reconfig_cache[key] = filt;
    return filt;
}
void OmniChannelNoiseDB::set_misconfigured(const std::vector<int>& channels, double from_gain, double from_shaping,
                                           double to_gain, double to_shaping, bool reset)
{
    if (reset) {
        auto def = default_filter();
        for (auto& it : m_db) {
            it.second.config = def;
        }
    }
    auto val = get_reconfig(from_gain, from_shaping, to_gain, to_shaping);
    for (int ch : channels) {
        // m_db.at(ch).config = val;
        // dbget(ch).config = val;
        get_ci(ch).config = val;
        m_miscfg_channels.push_back(ch);
    }
}

OmniChannelNoiseDB::shared_filter_t OmniChannelNoiseDB::parse_response(Json::Value jfilt)
{
    if (jfilt.isMember("wpid")) {
        WirePlaneId wpid(jfilt["wpid"].asInt());
        auto it = m_response_cache.find(wpid.ident());
        if (it != m_response_cache.end()) {
            return it->second;
        }
        auto const& fr = m_fr->field_response();
        auto fravg = Response::wire_region_average(fr);
        auto const& pr = fravg.planes[wpid.index()];

        // full length waveform
        std::vector<float> waveform(m_nsamples, 0.0);
        for (auto const& path : pr.paths) {
            auto const& current = path.current;
            const size_t nsamp = std::min(m_nsamples, (int) current.size());
            for (size_t ind = 0; ind < nsamp; ++ind) {
                waveform[ind] += current[ind];
            }
        }
        auto spectrum = fwd_r2c(m_dft, waveform);
        auto ret = std::make_shared<filter_t>(spectrum);
        m_response_cache[wpid.ident()] = ret;
        return ret;
    }

    if (jfilt.isMember("waveform") && jfilt.isMember("waveformid")) {
        int id = jfilt["waveformid"].asInt();
        auto it = m_waveform_cache.find(id);
        if (it != m_waveform_cache.end()) {
            return it->second;
        }

        auto jwave = jfilt["waveform"];
        const int nsamp = std::min(m_nsamples, (int) jwave.size());

        // Explicitly given waveform
        std::vector<float> waveform(m_nsamples, 0.0);
        for (int ind = 0; ind < nsamp; ++ind) {
            waveform[ind] = jwave[ind].asFloat();
        }

        auto spectrum = fwd_r2c(m_dft, waveform);
        auto ret = std::make_shared<filter_t>(spectrum);
        m_waveform_cache[id] = ret;
        return ret;
    }

    // this default return is special in that it's empty instead of
    // being a flat, unity spectrum.
    static shared_filter_t empty = std::make_shared<filter_t>();
    return empty;
}

OmniChannelNoiseDB::ChannelInfo& OmniChannelNoiseDB::get_ci(int chid)
{
    auto it = m_db.find(chid);
    if (it == m_db.end()) {
        THROW(KeyError() << errmsg{String::format("no db info for channel %d", chid)});
    }
    return it->second;
    // return m_db.at(chid);
}

// template <typename Type>
// void dump_cfg(const std::string& name, std::vector<int> chans, Type val)
// {
//     std::sort(chans.begin(), chans.end());
//     std::cerr << "OmniChannelNoiseDB: setting " << name << " to " << val << " on " << chans.size() << ":[" <<
//     chans.front() << "," << chans.back() << "]\n";
// }

void OmniChannelNoiseDB::update_channels(Json::Value cfg)
{
    auto chans = parse_channels(cfg["channels"]);

    if (cfg.isMember("nominal_baseline")) {
        double val = cfg["nominal_baseline"].asDouble();
        // dump_cfg("baseline", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).nominal_baseline = val;
            // dbget(ch).nominal_baseline = val;
            get_ci(ch).nominal_baseline = val;
        }
    }
    if (cfg.isMember("gain_correction")) {
        double val = cfg["gain_correction"].asDouble();
        // dump_cfg("gain", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).gain_correction = val;
            // dbget(ch).gain_correction = val;
            get_ci(ch).gain_correction = val;
        }
    }
    // fixme: why have two ways to set the same thing?
    {
        auto jfilt = cfg["reconfig"];
        if (!jfilt.isNull()) {
            auto val = parse_gain(jfilt);
            // dump_cfg("gain", chans, val);
            for (int ch : chans) {
                // m_db.at(ch).gain_correction = val;
                // dbget(ch).gain_correction = val;
                get_ci(ch).gain_correction = val;
            }
        }
    }

    if (cfg.isMember("response_offset")) {
        double val = cfg["response_offset"].asDouble();
        // dump_cfg("offset", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).response_offset = val;
            // dbget(ch).response_offset = val;
            get_ci(ch).response_offset = val;
        }
    }
    if (cfg.isMember("min_rms_cut")) {
        const auto& jv = cfg["min_rms_cut"];
        const bool linear = jv.isObject()
            && jv.get("type", "").asString() == "linear_in_wirelength";
        double scalar_val = linear ? 0.0 : jv.asDouble();
        double l0 = 0, v0 = 0, l1 = 1, v1 = 0;
        if (linear) {
            cache_wire_lengths();
            l0 = jv["l0"].asDouble(); v0 = jv["v0"].asDouble();
            l1 = jv["l1"].asDouble(); v1 = jv["v1"].asDouble();
        }
        for (int ch : chans) {
            if (linear) {
                auto it = m_wire_length_cm.find(ch);
                double L = (it == m_wire_length_cm.end()) ? l0 : it->second;
                double t = std::clamp((L - l0) / (l1 - l0), 0.0, 1.0);
                get_ci(ch).min_rms_cut = v0 + t * (v1 - v0);
            }
            else {
                get_ci(ch).min_rms_cut = scalar_val;
            }
        }
    }
    if (cfg.isMember("max_rms_cut")) {
        const auto& jv = cfg["max_rms_cut"];
        const bool linear = jv.isObject()
            && jv.get("type", "").asString() == "linear_in_wirelength";
        double scalar_val = linear ? 0.0 : jv.asDouble();
        double l0 = 0, v0 = 0, l1 = 1, v1 = 0;
        if (linear) {
            cache_wire_lengths();
            l0 = jv["l0"].asDouble(); v0 = jv["v0"].asDouble();
            l1 = jv["l1"].asDouble(); v1 = jv["v1"].asDouble();
        }
        for (int ch : chans) {
            if (linear) {
                auto it = m_wire_length_cm.find(ch);
                double L = (it == m_wire_length_cm.end()) ? l0 : it->second;
                double t = std::clamp((L - l0) / (l1 - l0), 0.0, 1.0);
                get_ci(ch).max_rms_cut = v0 + t * (v1 - v0);
            }
            else {
                get_ci(ch).max_rms_cut = scalar_val;
            }
        }
    }
    if (cfg.isMember("pad_window_front")) {
        int val = cfg["pad_window_front"].asDouble();
        // dump_cfg("padfront", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).pad_window_front = val;
            // dbget(ch).pad_window_front = val;
            get_ci(ch).pad_window_front = val;
        }
    }
    if (cfg.isMember("pad_window_back")) {
        int val = cfg["pad_window_back"].asDouble();
        // dump_cfg("padback", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).pad_window_back = val;
            // dbget(ch).pad_window_back = val;
            get_ci(ch).pad_window_back = val;
        }
    }

    if (cfg.isMember("decon_limit")) {
        float val = cfg["decon_limit"].asDouble();
        // dump_cfg("deconlimit", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).decon_limit = val;
            // dbget(ch).decon_limit = val;
            get_ci(ch).decon_limit = val;
        }
    }

    if (cfg.isMember("decon_lf_cutoff")) {
        float val = cfg["decon_lf_cutoff"].asDouble();
        // dump_cfg("deconlfcutoff", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).decon_limit = val;
            // dbget(ch).decon_limit = val;
            get_ci(ch).decon_lf_cutoff = val;
        }
    }

    if (cfg.isMember("decon_limit1")) {
        float val = cfg["decon_limit1"].asDouble();
        // dump_cfg("deconlimit1", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).decon_limit1 = val;
            // dbget(ch).decon_limit1 = val;
            get_ci(ch).decon_limit1 = val;
        }
    }
    if (cfg.isMember("adc_limit")) {
        float val = cfg["adc_limit"].asDouble();
        // dump_cfg("adclimit", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).adc_limit = val;
            // dbget(ch).adc_limit = val;
            get_ci(ch).adc_limit = val;
        }
    }
    if (cfg.isMember("protection_factor")) {
        float val = cfg["protection_factor"].asDouble();
        // dump_cfg("protectionfactor", chans, val);
        for (int ch : chans) {
            get_ci(ch).protection_factor = val;
        }
    }
    if (cfg.isMember("min_adc_limit")) {
        float val = cfg["min_adc_limit"].asDouble();
        // dump_cfg("minadclimit", chans, val);
        for (int ch : chans) {
            // m_db.at(ch).adc_limit = val;
            // dbget(ch).adc_limit = val;
            get_ci(ch).min_adc_limit = val;
        }
    }
    if (cfg.isMember("roi_min_max_ratio")) {
        float val = cfg["roi_min_max_ratio"].asDouble();
        // dump_cfg("roiminmaxratio", chans, val);
        for (int ch : chans) {
            get_ci(ch).roi_min_max_ratio = val;
        }
    }

    {
        auto jfilt = cfg["rcrc"];
        if (!jfilt.isNull()) {
            if (cfg.isMember("rc_layers")) {
                m_rc_layers = cfg["rc_layers"].asInt();
            }
            // std::cerr << "rc_layers = " << m_rc_layers << std::endl;
            auto val = parse_rcrc(jfilt, m_rc_layers);
            // dump_cfg("rcrc", chans, Waveform::sum(*val));
            for (int ch : chans) {
                // m_db.at(ch).rcrc = val;
                // dbget(ch).rcrc = val;
                get_ci(ch).rcrc = val;
            }
        }
    }
    {
        auto jfilt = cfg["reconfig"];
        if (!jfilt.isNull()) {
            auto val = parse_reconfig(jfilt);
            // dump_cfg("reconfig", chans, Waveform::sum(*val));
            for (int ch : chans) {
                // m_db.at(ch).config = val;
                // dbget(ch).config = val;
                get_ci(ch).config = val;
                // fill in misconfgured channels
                // std::cout <<" miscfg_channels fill: "<< ch <<"\n"
                if (!jfilt.empty()) m_miscfg_channels.push_back(ch);
            }
        }
    }
    {
        auto jfilt = cfg["freqmasks"];
        if (!jfilt.isNull()) {
            auto val = parse_freqmasks(jfilt);
            // dump_cfg("freqmasks", chans, Waveform::sum(*val));
            // std::cerr << jfilt << std::endl;
            for (int ch : chans) {
                // m_db.at(ch).noise = val;
                // dbget(ch).noise = val;
                get_ci(ch).noise = val;
            }
        }
    }
    {
        auto jfilt = cfg["response"];
        if (!jfilt.isNull()) {
            auto val = parse_response(jfilt);
            // dump_cfg("response", chans, Waveform::sum(*val));
            for (int ch : chans) {
                // m_db.at(ch).response = val;
                // dbget(ch).response = val;
                get_ci(ch).response = val;
            }
        }
    }
}

void OmniChannelNoiseDB::configure(const WireCell::Configuration& cfg)
{
    m_tick = get(cfg, "tick", m_tick);
    m_nsamples = get(cfg, "nsamples", m_nsamples);
    std::string anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);
    std::string fr_tn = get<std::string>(cfg, "field_response", "FieldResponse");
    m_fr = Factory::find_tn<IFieldResponse>(fr_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    // WARNING: this assumes channel numbers count from 0 with no gaps!
    // int nchans = m_anode->channels().size();
    // std::cerr << "noise database with " << nchans << " channels\n";
    // m_db.resize(nchans);

    // clear any previous config, and recover the memory
    // for (auto it = m_db.begin(); it!= m_db.end(); it++){
    //    delete it->second;
    // }
    // m_db.clear();
    for (auto ch : m_anode->channels()) {
        // m_db.insert(std::make_pair(ch, new ChannelInfo));
        m_db[ch] = ChannelInfo();
    }

    m_channel_groups.clear();
    auto jgroups = cfg["groups"];
    for (auto jgroup : jgroups) {
        std::vector<int> channel_group;
        for (auto jch : jgroup) {
            channel_group.push_back(jch.asInt());
        }
        m_channel_groups.push_back(channel_group);
    }
    m_bad_channels.clear();
    for (auto jch : cfg["bad"]) {
        m_bad_channels.push_back(jch.asInt());
    }
    std::sort(m_bad_channels.begin(), m_bad_channels.end());
    if (m_bad_channels.size()) {
        log->debug("OmniChannelNoiseDB: setting {}:[{},{}] bad channels", m_bad_channels.size(), m_bad_channels.front(),
                   m_bad_channels.back());
    }

    m_miscfg_channels.clear();
    m_cfg_channel_info = cfg["channel_info"];
    for (auto jci : m_cfg_channel_info) {
        update_channels(jci);
    }
}

void OmniChannelNoiseDB::set_nsamples(int n)
{
    if (n == m_nsamples) return;
    log->info("OmniChannelNoiseDB: nsamples {} -> {}, rebuilding spectra", m_nsamples, n);
    m_nsamples = n;
    m_rcrc_cache.clear();
    m_reconfig_cache.clear();
    m_waveform_cache.clear();
    m_response_cache.clear();
    for (auto& kv : m_db) {
        kv.second.rcrc = nullptr;
        kv.second.config = nullptr;
        kv.second.noise = nullptr;
        kv.second.response = nullptr;
    }
    for (auto jci : m_cfg_channel_info) {
        update_channels(jci);
    }
}

double OmniChannelNoiseDB::sample_time() const { return m_tick; }

double OmniChannelNoiseDB::nominal_baseline(int channel) const { return dbget(channel).nominal_baseline; }

double OmniChannelNoiseDB::gain_correction(int channel) const { return dbget(channel).gain_correction; }

double OmniChannelNoiseDB::response_offset(int channel) const { return dbget(channel).response_offset; }

double OmniChannelNoiseDB::min_rms_cut(int channel) const { return dbget(channel).min_rms_cut; }

double OmniChannelNoiseDB::max_rms_cut(int channel) const { return dbget(channel).max_rms_cut; }

int OmniChannelNoiseDB::pad_window_front(int channel) const { return dbget(channel).pad_window_front; }

int OmniChannelNoiseDB::pad_window_back(int channel) const { return dbget(channel).pad_window_back; }

float OmniChannelNoiseDB::coherent_nf_decon_limit(int channel) const { return dbget(channel).decon_limit; }

float OmniChannelNoiseDB::coherent_nf_decon_lf_cutoff(int channel) const { return dbget(channel).decon_lf_cutoff; }

float OmniChannelNoiseDB::coherent_nf_decon_limit1(int channel) const { return dbget(channel).decon_limit1; }

float OmniChannelNoiseDB::coherent_nf_adc_limit(int channel) const { return dbget(channel).adc_limit; }

float OmniChannelNoiseDB::coherent_nf_protection_factor(int channel) const { return dbget(channel).protection_factor; }

float OmniChannelNoiseDB::coherent_nf_min_adc_limit(int channel) const { return dbget(channel).min_adc_limit; }

float OmniChannelNoiseDB::coherent_nf_roi_min_max_ratio(int channel) const { return dbget(channel).roi_min_max_ratio; }

const IChannelNoiseDatabase::filter_t& OmniChannelNoiseDB::rcrc(int channel) const
{
    auto filt = dbget(channel).rcrc;
    if (filt) {
        return *filt;
    }
    static filter_t dummy;
    return dummy;
}

const IChannelNoiseDatabase::filter_t& OmniChannelNoiseDB::config(int channel) const
{
    auto filt = dbget(channel).config;
    if (filt) {
        return *filt;
    }
    static filter_t dummy;
    return dummy;
}

const IChannelNoiseDatabase::filter_t& OmniChannelNoiseDB::noise(int channel) const
{
    auto filt = dbget(channel).noise;
    if (filt) {
        return *filt;
    }
    static filter_t dummy;
    return dummy;
}

const IChannelNoiseDatabase::filter_t& OmniChannelNoiseDB::response(int channel) const
{
    auto filt = dbget(channel).response;
    if (filt) {
        return *filt;
    }
    static filter_t dummy;
    return dummy;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
