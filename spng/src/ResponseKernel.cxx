#include "WireCellSpng/ResponseKernel.h"
#include "WireCellSpng/Util.h"
#include "WireCellSpng/DeconTools.h"
#include "WireCellSpng/Convo.h"
#include "WireCellSpng/TorchLMN.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include "boost/container_hash/hash.hpp"

WIRECELL_FACTORY(SPNGResponseKernel,
                 WireCell::SPNG::ResponseKernel,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

using namespace WireCell::HanaJsonCPP;                         
using torch::nn::functional::pad;
using torch::nn::functional::PadFuncOptions;

namespace WireCell::SPNG {

    ResponseKernel::ResponseKernel()
        : Logger("ResponseKernel", "spng")
        , m_cache(1)
    {}
    ResponseKernel::ResponseKernel(const ResponseKernelConfig& cfg)
        : Logger("ResponseKernel", "spng")
        , m_cfg(cfg)
        , m_cache(m_cfg.capacity)
    {
        configme();
    }

    WireCell::Configuration ResponseKernel::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = this->Logger::default_configuration();
        update(cfg, cfg2);
        cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    void ResponseKernel::configure(const WireCell::Configuration& config)
    {
        this->ContextBase::configure(config);
        this->Logger::configure(config);
        from_json(m_cfg, config);
        configme();
    }
    
    // convolve then downsample
    void ResponseKernel::configme_ctd()
    {
        log->debug("configuring FR*ER for convolve-then-downsample");
        auto ier = Factory::find_tn<IWaveform>(m_cfg.elec_response);
        auto ifr = Factory::find_tn<IFieldResponse>(m_cfg.field_response);

        const double er_period = ier->waveform_period();
        // Some FRs have FP round off problems in the period.
        const double fr_period_fine = round(ifr->field_response().period/units::ns)*units::ns;

        // The natural number of samples for ER.  WARNING: old config sets this
        // to nticks in ADC readout which is absurdly long.  FIXME: instead of
        // fighting history, add a dedicated "nerticks" config to this component.
        const size_t nerticks = ier->waveform_samples().size();
        const double nerticks_duration = nerticks * er_period;
        const int64_t nerticks_fine = nerticks * er_period/fr_period_fine;
        Binning er_bins(nerticks_fine, 0, nerticks_duration);
        log->debug("ER: {} ticks, {} us duration", nerticks_fine, nerticks_duration/units::us);

        std::vector<float> erv = ier->waveform_samples(er_bins);
        /// we hold on to er_fine so better not from_blob it!
        // auto er_fine = torch::from_blob(erv.data(), nerticks_fine, torch::kFloat32);
        auto er_fine = torch::tensor(erv, torch::kFloat32);
        m_raw_er = er_fine;

        // Still at FR sampling, eg 100ns
        m_raw_fr_full_fine = fr_tensor(ifr->field_response(), m_cfg.plane_index);
        auto fr_fine = fr_average_tensor(ifr->field_response(), m_cfg.plane_index);
        m_raw_fr_avg_fine = fr_fine;

        // Find size to assure linear convolution.
        const int64_t linear_size = linear_shape(fr_fine.size(1), er_fine.size(0));

        // pad is LAST FIRST
        auto fr_pad = pad(fr_fine, PadFuncOptions({0, linear_size - fr_fine.size(1), 0, 0}));
        auto er_pad = pad(er_fine, PadFuncOptions({0, linear_size - er_fine.size(0)}));
                          
        auto FR = torch::fft::fft(fr_pad, {}, 1);
        auto ER = torch::fft::fft(er_pad).unsqueeze(0);
        log->debug("FR: {}", to_string(FR));
        log->debug("ER: {}", to_string(ER));
        auto FRER = FR*ER;      // "natural spectrum"

        auto frer_fine = torch::real(torch::fft::ifft(FRER, {}, 1));
        m_raw_fr = frer_fine; // repurposed to hold pre-resampled.

        auto frer_coarse = LMN::resample_interval(frer_fine, fr_period_fine, er_period, 1);

        log->debug("nerticks={}, er_period={}, er_fine={}, fr_period={}, fr_fine={}, linear_size={}, fr_pad={}, er_pad={}, frer_fine={}, frer_coarse={}",
                   nerticks, er_period, to_string(er_fine),
                   fr_period_fine, to_string(fr_fine),
                   linear_size,
                   to_string(fr_pad), to_string(er_pad),
                   to_string(frer_fine), to_string(frer_coarse));

        m_response_waveform = frer_coarse;
    }

    // downsample then convolve
    void ResponseKernel::configme_dtc()
    {
        log->debug("configuring FR*ER for downsample-then-convolve");
        auto ier = Factory::find_tn<IWaveform>(m_cfg.elec_response);
        std::vector<float> erv = ier->waveform_samples();
        // We hold on to this so better not from_blob it!
        // auto er = torch::from_blob(erv.data(), {static_cast<int64_t>(erv.size())}, torch::kFloat32);
        auto er = torch::tensor(erv, torch::kFloat32);
        m_raw_er = er;

        auto ifr = Factory::find_tn<IFieldResponse>(m_cfg.field_response);
        // Still at FR sampling, eg 100ns
        m_raw_fr_full_fine = fr_tensor(ifr->field_response(), m_cfg.plane_index);
        auto fr_fine = fr_average_tensor(ifr->field_response(), m_cfg.plane_index);
        m_raw_fr_avg_fine = fr_fine;

        // We now must do two transformations prior to the FR*ER convolution.
        //
        // 1. Downsample FR from the fine FR sampling period to the coarse ER
        // sampling period.  This is done using LMN resampling which applies a
        // "perfect" low-pass filter.
        //
        // 2. Pad both FR and ER time dimension in INTERVAL space to assure
        // linear convolution.
        //
        // I'm sure there is a trick to minimize the FFT round trips, but for
        // now we do the dumb, straightforward thing.  It is a one-shot anyways.


        const double er_period = ier->waveform_period();
        // Some FRs have FP round off problems in the period.
        const double fr_period_fine = round(ifr->field_response().period/units::ns)*units::ns;

        // Resample under default interpolation interpretation.  We keep the
        // default to return it to interval space as we must pad for the FR*ER
        // convolution.
        auto fr = LMN::resample_interval(fr_fine, fr_period_fine, er_period, 1);
        m_raw_fr = fr;

        // Find size to assure linear convolution.
        const int64_t linear_size = linear_shape(fr.size(1), er.size(0));

        // pad is LAST FIRST
        auto fr_pad = pad(fr, PadFuncOptions({0, linear_size - fr.size(1), 0, 0}));
        auto er_pad = pad(er, PadFuncOptions({0, linear_size - er.size(0)}));
                          
        log->debug("linear_size={} er={} fr_fine={} fr={} fr_pad={} er_pad={}",
                   linear_size,
                   to_string(er), to_string(fr_fine), to_string(fr),
                   to_string(fr_pad), to_string(er_pad));

        auto FR = torch::fft::fft(fr_pad, {}, 1);
        auto ER = torch::fft::fft(er_pad).unsqueeze(0);
        log->debug("FR: {}", to_string(FR));
        log->debug("ER: {}", to_string(ER));
        auto FRER = FR*ER;      // "natural spectrum"

        auto frer = torch::fft::ifft(FRER, {}, 1);
        m_response_waveform = torch::real(frer);

    }
    
    void ResponseKernel::configme()
    {
        if (m_cfg.plane_index < 0 || m_cfg.plane_index>3) {
            raise<ValueError>("illegal plane ID: %d", m_cfg.plane_index);
        }
        m_cache.reset_capacity(m_cfg.capacity);

        if (false) {
            configme_dtc();
        }
        else {
            configme_ctd();
        }

        if (m_cfg.debug_filename.size()) {
            log->debug("writing debug file: {}", m_cfg.debug_filename);
            write_debug(m_cfg.debug_filename);
        }

    }

    void ResponseKernel::write_debug(const std::string& filename) const
    {
        using tensor_map = torch::Dict<std::string, torch::Tensor>;
        tensor_map to_save;
        to_save.insert("field_response_full_fine", m_raw_fr_full_fine);
        to_save.insert("field_response_avg_fine", m_raw_fr_avg_fine);
        to_save.insert("field_response", m_raw_fr);
        to_save.insert("elec_response", m_raw_er);

        to_save.insert("response_waveform", m_response_waveform);
        std::vector<int64_t> s = {100,2000};
        to_save.insert("padded_waveform", pad_waveform(s));
        to_save.insert("padded_spectrum", make_spectrum(s));
        auto data = torch::pickle_save(to_save);
        std::ofstream output_file(filename, std::ios::binary);
        output_file.write(data.data(), data.size());
    }

    size_t ResponseKernel::make_cache_key(const shape_t& shape) const
    {
        size_t h = 0;
        for (const auto& s : shape) {
            boost::hash_combine(h, s);
        }
        return h;
    }

    torch::Tensor ResponseKernel::make_spectrum(const std::vector<int64_t> & shape) const
    {
        return torch::fft::fft2(pad_waveform(shape));
    }

    torch::Tensor ResponseKernel::pad_waveform(const std::vector<int64_t> & shape) const
    {
        auto tmp = LMN::resize_middle(m_response_waveform, shape[0], 0);
        return LMN::resize(tmp, shape[1], 1);
    }


    torch::Tensor ResponseKernel::spectrum(const std::vector<int64_t> & shape) const
    {
        if (m_cfg.capacity == 0) {
            return make_spectrum(shape);
        }

        const size_t key = make_cache_key(shape);
        torch::Tensor ten = m_cache.get_or_compute(key, [&]() {
            return this->make_spectrum(shape);
        });
        return ten;
    }

    std::vector<int64_t> ResponseKernel::shape() const {
        return m_response_waveform.sizes().vec();
    }

    std::vector<int64_t> ResponseKernel::shifts() const {
        std::vector<int64_t> ret = {0,0};
        return ret;
    }

}
