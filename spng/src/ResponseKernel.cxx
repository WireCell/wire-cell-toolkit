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
    
    void ResponseKernel::configme()
    {
        if (m_cfg.plane_index < 0 || m_cfg.plane_index>3) {
            raise<ValueError>("illegal plane ID: %d", m_cfg.plane_index);
        }
        m_cache.reset_capacity(m_cfg.capacity);

        //
        // FIXME: ER gets configured with like 6000 ticks.
        //
        auto ier = Factory::find_tn<IWaveform>(m_cfg.elec_response);
        auto erv = ier->waveform_samples();
        auto er = torch::from_blob(erv.data(), {static_cast<int64_t>(erv.size())}, torch::kFloat32);

        auto ifr = Factory::find_tn<IFieldResponse>(m_cfg.field_response);
        // Still at FR sampling, eg 100ns
        auto fr_fine = fr_average_tensor(ifr, m_cfg.plane_index);

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

    size_t ResponseKernel::make_cache_key(const shape_t& shape) const {
        size_t h = 0;
        for (const auto& s : shape) {
            boost::hash_combine(h, s);
        }
        return h;
    }

    torch::Tensor ResponseKernel::make_spectrum(const std::vector<int64_t> & shape) const
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
