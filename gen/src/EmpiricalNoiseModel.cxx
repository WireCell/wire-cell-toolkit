
#include "WireCellGen/EmpiricalNoiseModel.h"

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Response.h"  // fixme: should remove direct dependency
#include "WireCellUtil/Waveform.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/FFTBestLength.h"

#include "WireCellAux/DftTools.h"

#include <iostream>  // debug

WIRECELL_FACTORY(EmpiricalNoiseModel,
                 WireCell::Gen::EmpiricalNoiseModel,
                 WireCell::INamed,
                 WireCell::IChannelSpectrum,
                 WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::DftTools::fwd_r2c;

static int fft_length(int nsamples)
{
    return fft_best_length(nsamples);
    // return nsamples;
}


Gen::EmpiricalNoiseModel::EmpiricalNoiseModel(const std::string& spectra_file, const int nsamples, const double period,
                                              const double wire_length_scale,
                                              // const double time_scale,
                                              // const double gain_scale,
                                              // const double freq_scale,
                                              const std::string anode_tn, const std::string chanstat_tn)
  : Aux::Logger("EmpiricalNoiseModel", "gen")
  , m_spectra_file(spectra_file)
  , m_nsamples(nsamples)
  , m_period(period)
  , m_wlres(wire_length_scale)
  // , m_tres(time_scale)
  // , m_gres(gain_scale)
  // , m_fres(freq_scale)
  , m_anode_tn(anode_tn)
  , m_chanstat_tn(chanstat_tn)
  , m_len_amp_cache(4)
{
    m_fft_length = fft_length(m_nsamples);
    gen_elec_resp_default();
}

Gen::EmpiricalNoiseModel::~EmpiricalNoiseModel() {}

void Gen::EmpiricalNoiseModel::gen_elec_resp_default()
{
    // Calculate the frequencies.  Despite the name of this method,
    // these frequencies are not limited to just electronics response.
    m_elec_resp_freq.resize(m_fft_length, 0);
    for (unsigned int i = 0; i != m_elec_resp_freq.size(); i++) {
        if (i <= m_elec_resp_freq.size() / 2.) {
            m_elec_resp_freq.at(i) =
                i / (m_elec_resp_freq.size() * 1.0) * 1. / m_period;  // the second half is useless ...
        }
        else {
            m_elec_resp_freq.at(i) = (m_elec_resp_freq.size() - i) / (m_elec_resp_freq.size() * 1.0) * 1. /
                                     m_period;  // the second half is useless ...
        }
    }
}

WireCell::Configuration Gen::EmpiricalNoiseModel::default_configuration() const
{
    Configuration cfg;

    /// The file holding the spectral data
    cfg["spectra_file"] = m_spectra_file;
    cfg["nsamples"] = m_nsamples;  // number of samples up to Nyquist frequency
    cfg["period"] = m_period;
    cfg["wire_length_scale"] = m_wlres;  //
    cfg["gain_scale"] = m_gres;
    cfg["shaping_scale"] = m_sres;
    // cfg["time_scale"] = m_tres;
    // cfg["gain_scale"] = m_gres;
    // cfg["freq_scale"] = m_fres;
    cfg["anode"] = m_anode_tn;  // name of IAnodePlane component
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use

    return cfg;
}

void Gen::EmpiricalNoiseModel::resample(NoiseSpectrum& spectrum) const
{
    log->debug("m_nsamples={} fft_length={} m_period={} spec: nsamples={} period={}",
               m_nsamples, m_fft_length, m_period, spectrum.nsamples, spectrum.period);

    if (spectrum.nsamples == m_fft_length && spectrum.period == m_period) {
        return;              // natural sampling is what is requested.
    }

    // Scale the amplitude.  Note, period term moved to *AddNoise.
    double scale = sqrt(m_fft_length / (spectrum.nsamples * 1.0));
    // double scale = sqrt(m_fft_length / (spectrum.nsamples * 1.0) * spectrum.period / (m_period * 1.0));

    spectrum.constant *= scale;
    for (unsigned int ind = 0; ind != spectrum.amps.size(); ind++) {
        spectrum.amps[ind] *= scale;
    }

    // interpolate ...

    amplitude_t temp_amplitudes(m_fft_length, 0);
    int count_low = 0;
    int count_high = 1;
    double mu = 0;
    for (int i = 0; i != m_fft_length; i++) {
        double frequency = m_elec_resp_freq.at(i);
        if (frequency <= spectrum.freqs[0]) {
            count_low = 0;
            count_high = 1;
            mu = 0;
        }
        else if (frequency >= spectrum.freqs.back()) {
            count_low = spectrum.freqs.size() - 2;
            count_high = spectrum.freqs.size() - 1;
            mu = 1;
        }
        else {
            for (unsigned int j = 0; j != spectrum.freqs.size(); j++) {
                if (frequency > spectrum.freqs.at(j)) {
                    count_low = j;
                    count_high = j + 1;
                }
                else {
                    break;
                }
            }
            mu = (frequency - spectrum.freqs.at(count_low)) /
                 (spectrum.freqs.at(count_high) - spectrum.freqs.at(count_low));
        }

        temp_amplitudes.at(i) = (1 - mu) * spectrum.amps[count_low] + mu * spectrum.amps[count_high];
    }

    spectrum.amps = temp_amplitudes;
    // spectrum.amps.push_back(spectrum.constant);

    return;
}

void Gen::EmpiricalNoiseModel::configure(const WireCell::Configuration& cfg)
{
    int errors = 0;

    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

    m_chanstat_tn = get(cfg, "chanstat", m_chanstat_tn);
    if (m_chanstat_tn != "") {  // allow for an empty channel status, no deviation from norm
        m_chanstat = Factory::find_tn<IChannelStatus>(m_chanstat_tn);
    }
    log->debug("using anode {}, chanstat {}", m_anode_tn, m_chanstat_tn);

    auto jspectra = cfg["spectra_file"]; // best if this were renamed "spectra"...
    if (jspectra.isNull()) {
        log->critical("required configuration parameter \"spectra_file\" is empty");
        ++errors;
    }
    if (jspectra.isString()) {
        m_spectra_file = jspectra.asString();
        jspectra = Persist::load(m_spectra_file);
    }

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_nsamples = get(cfg, "nsamples", m_nsamples);
    m_fft_length = fft_length(m_nsamples);
    m_period = get(cfg, "period", m_period);
    m_wlres = get(cfg, "wire_length_scale", m_wlres);
    // m_tres = get(cfg, "time_scale", m_tres);
    // m_gres = get(cfg, "gain_scale", m_gres);
    // m_fres = get(cfg, "freq_scale", m_fres);

    // needs m_period and m_fft_length
    gen_elec_resp_default();

    // Load spectral data.  Fixme: should break out this code
    // separate.
    m_spectral_data.clear();
    for (const auto& jentry : jspectra) {
        NoiseSpectrum* nsptr = new NoiseSpectrum();

        // optional for this model

        nsptr->plane = jentry["plane"].asInt();
        nsptr->gain = jentry["gain"].asFloat();        // *m_gres;
        nsptr->shaping = jentry["shaping"].asFloat();  // * m_tres;
        nsptr->wirelen = jentry["wirelen"].asFloat();  // * m_wlres;
        nsptr->constant = get(jentry, "const", 0.0);

        // required from all WCT noise spectra

        nsptr->nsamples = jentry["nsamples"].asInt();
        nsptr->period = jentry["period"].asFloat();    // * m_tres;

        auto jfreqs = jentry["freqs"];
        const int nfreqs = jfreqs.size();
        nsptr->freqs.resize(nfreqs, 0.0);
        for (int ind = 0; ind < nfreqs; ++ind) {
            nsptr->freqs[ind] = jfreqs[ind].asFloat();  // * m_fres;
        }
        auto jamps = jentry["amps"];
        const int namps = jamps.size();
        nsptr->amps.resize(namps, 0.0);
        for (int ind = 0; ind < namps; ++ind) {
            nsptr->amps[ind] = jamps[ind].asFloat();
        }

        resample(*nsptr);
        m_spectral_data[nsptr->plane].push_back(nsptr);  // assumes ordered by wire length!
        log->debug("nwanted={} plane={} ntold={} ngot={} ninput={}",
                   m_nsamples, nsptr->plane, nsptr->nsamples, nsptr->amps.size(), nfreqs);
    }

    if (errors) {
        THROW(ValueError() << errmsg{"EmpiricalNoiseModel: configuration errors"});
    }
}


Gen::EmpiricalNoiseModel::spectrum_data_t
Gen::EmpiricalNoiseModel::interpolate_wire_length(int iplane, double wire_length) const
{
    auto it = m_spectral_data.find(iplane);
    if (it == m_spectral_data.end()) {
        THROW(ValueError() << errmsg{"EmpiricalNoiseModel: unknown plane index"});
    }

    const auto& spectra = it->second;

    // Linearly interpolate spectral data to wire length.

    // Any wire lengths
    // which are out of bounds causes the nearest result to be
    // returned (flat extrapolation)
    const NoiseSpectrum* front = spectra.front();
    if (wire_length <= front->wirelen) {
        return spectrum_data_t{front->constant, front->amps};
    }
    const NoiseSpectrum* back = spectra.back();
    if (wire_length >= back->wirelen) {
        return spectrum_data_t{back->constant, back->amps};
    }

    const int nspectra = spectra.size();
    for (int ind = 1; ind < nspectra; ++ind) {
        const NoiseSpectrum* hi = spectra[ind];
        if (hi->wirelen < wire_length) {
            continue;
        }
        const NoiseSpectrum* lo = spectra[ind - 1];

        const double delta = hi->wirelen - lo->wirelen;
        const double dist = wire_length - lo->wirelen;
        const double mu = dist / delta;

        const int nsamp = lo->amps.size();
        amplitude_t amp(nsamp);
        for (int isamp = 0; isamp < nsamp; ++isamp) {
            amp[isamp] = lo->amps[isamp] * (1 - mu) + hi->amps[isamp] * mu;
        }
        double constant = lo->constant*(1-mu) + hi->constant*mu;
        return spectrum_data_t{constant, amp};
    }
    // should not get here
    THROW(ValueError() << errmsg{"EmpiricalNoiseModel: data corruption"});
}

const std::vector<float>& Gen::EmpiricalNoiseModel::freq() const
{
    // assume frequency is universal ...
    const int iplane = 0;
    auto it = m_spectral_data.find(iplane);
    const auto& spectra = it->second;
    return spectra.front()->freqs;
}

const double Gen::EmpiricalNoiseModel::shaping_time(int chid) const
{
    auto wpid = m_anode->resolve(chid);
    const int iplane = wpid.index();
    auto it = m_spectral_data.find(iplane);
    const auto& spectra = it->second;
    return spectra.front()->shaping;
}

const double Gen::EmpiricalNoiseModel::gain(int chid) const
{
    auto wpid = m_anode->resolve(chid);
    const int iplane = wpid.index();
    auto it = m_spectral_data.find(iplane);
    const auto& spectra = it->second;
    return spectra.front()->gain;
}

const Gen::EmpiricalNoiseModel::spectrum_data_t&
Gen::EmpiricalNoiseModel::get_spectrum_data(int iplane, int ilen) const
{
    if (iplane < 0 or iplane >= (int)m_len_amp_cache.size()) {
        THROW(ValueError() << errmsg{"EmpiricalNoiseModel: illegal plane index"});
    }

    auto& amp_cache = m_len_amp_cache[iplane];
    auto lenamp = amp_cache.find(ilen);
    if (lenamp == amp_cache.end()) {
        amp_cache[ilen] = interpolate_wire_length(iplane, ilen * m_wlres);
    }
    lenamp = amp_cache.find(ilen);
    return lenamp->second;
}

const IChannelSpectrum::amplitude_t&
Gen::EmpiricalNoiseModel::channel_spectrum(int chid) const
{
    auto packkeyit = m_channel_packkey_cache.find(chid);
    if (packkeyit != m_channel_packkey_cache.end()) {
        unsigned int packkey = packkeyit->second;
        auto cacit = m_channel_amp_cache.find(packkey);
        if (cacit != m_channel_amp_cache.end()) {
            return cacit->second;
        }
    }

    // get truncated wire length for cache
    int ilen = 0;
    auto chlen = m_chid_to_intlen.find(chid);
    if (chlen == m_chid_to_intlen.end()) {  // new channel
        auto wires = m_anode->wires(chid);  // sum up wire length
        double len = 0.0;
        for (auto wire : wires) {
            len += ray_length(wire->ray());
        }
        // cache every cm ...
        ilen = int(len / m_wlres);
        // there might be  aproblem with the wire length's unit
        // ilen = int(len);
        m_chid_to_intlen[chid] = ilen; // fixme: not thread safe
    }
    else {
        ilen = chlen->second;
    }

    // saved content
    auto wpid = m_anode->resolve(chid);
    const int iplane = wpid.index();

    // Must copy as next apply channel-specific changes.
    auto [constant, amp] = get_spectrum_data(iplane, ilen);
    const size_t namps = amp.size();

    // Convert from as-built to as-designed.
    const double db_gain = gain(chid);
    const double db_shaping = shaping_time(chid);
    double ch_gain = db_gain, ch_shaping = db_shaping;
    if (m_chanstat) {  // allow for deviation from nominal
        ch_gain = m_chanstat->preamp_gain(chid);
        ch_shaping = m_chanstat->preamp_shaping(chid);
    }

    // Convert gain.
    if (fabs(ch_gain - db_gain) > 0.01 * ch_gain) {
        Waveform::scale(amp, ch_gain / db_gain);
    }

    // Convert shaping
    if (fabs(ch_shaping - db_shaping) > 0.01 * ch_shaping) {
        // scale the amplitude by different shaping time ...
        int nconfig = ch_shaping / units::us / 0.1;
        auto resp1 = m_elec_resp_cache.find(nconfig);
        if (resp1 == m_elec_resp_cache.end()) {
            Response::ColdElec elec_resp(10, ch_shaping);  // default at 1 mV/fC
            auto sig = elec_resp.generate(WireCell::Waveform::Domain(0, namps * m_period), namps);
            auto filt = fwd_r2c(m_dft, sig);
            auto ele_resp_amp = Waveform::magnitude(filt);

            ele_resp_amp.resize(m_elec_resp_freq.size());
            m_elec_resp_cache[nconfig] = ele_resp_amp;
        }
        resp1 = m_elec_resp_cache.find(nconfig);

        nconfig = db_shaping / units::us / 0.1;
        auto resp2 = m_elec_resp_cache.find(nconfig);
        if (resp2 == m_elec_resp_cache.end()) {
            Response::ColdElec elec_resp(10, db_shaping);  // default at 1 mV/fC
            auto sig = elec_resp.generate(WireCell::Waveform::Domain(0, namps * m_period), namps);
            auto filt = fwd_r2c(m_dft, sig);
            auto ele_resp_amp = Waveform::magnitude(filt);

            ele_resp_amp.resize(m_elec_resp_freq.size());
            m_elec_resp_cache[nconfig] = ele_resp_amp;
        }
        resp2 = m_elec_resp_cache.find(nconfig);

        for (size_t ind = 0; ind < namps; ++ind) {
            amp[ind] *= resp1->second[ind] / resp2->second[ind];
        }
    }

    // add the constant terms ...
    const double constant_squared = constant*constant;
    for (auto& val : amp) {
        val = sqrt(val*val + constant_squared); // units still in mV
    }

    struct PWGS {
        unsigned int p;
        unsigned int w;
        unsigned int g;
        unsigned int s;
        unsigned int operator() () {
            return (
                ((p & 0b111) << 29) // 0 - 65535 length steps
                | ((w & 0b111111111111) << 16) // 0 - 65535 length steps
                | ((g & 0xFF) << 8) // 0 - 255 gain steps
                | (s & 0xFF)
                ); // 0 - 255 shaping steps
        } 
    } pack{(unsigned int)iplane, (unsigned int)ilen, (unsigned int)(ch_gain/m_gres), (unsigned int)(ch_shaping/m_sres)};

    unsigned int packkey = pack();
    m_channel_packkey_cache[chid] = packkey;
    // m_channel_amp_cache.clear();
    m_channel_amp_cache[packkey] = amp;
    return channel_spectrum(chid);
}
