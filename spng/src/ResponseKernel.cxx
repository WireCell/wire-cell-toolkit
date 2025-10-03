#include "WireCellSpng/ResponseKernel.h"
#include "WireCellSpng/Util.h"
#include "WireCellSpng/DeconTools.h"
#include "WireCellSpng/Convo.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include "boost/container_hash/hash.hpp"

using namespace WireCell::HanaJsonCPP;                         
using torch::nn::functional::pad;
using torch::nn::functional::PadFuncOptions;

WIRECELL_FACTORY(SPNGResponseKernel,
                 WireCell::SPNG::ResponseKernel,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)


namespace WireCell::SPNG {

    ResponseKernel::ResponseKernel()
        : m_cache(1)
    {}
    ResponseKernel::ResponseKernel(const ResponseKernelConfig& cfg)
        : m_cfg(cfg)
        , m_cache(m_cfg.capacity)
    {
        configme();
    }

    WireCell::Configuration ResponseKernel::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    void ResponseKernel::configure(const WireCell::Configuration& config)
    {
        this->ContextBase::configure(config);
        from_json(m_cfg, config);
        configme();
    }
    
    void ResponseKernel::configme()
    {
        if (m_cfg.plane_id < 0 || m_cfg.plane_id>3) {
            raise<ValueError>("illegal plane ID: %d", m_cfg.plane_id);
        }
        m_cache.reset_capacity(m_cfg.capacity);

        auto ier = Factory::find_tn<IWaveform>(m_cfg.elec_response);
        auto erv = ier->waveform_samples();
        auto er = torch::from_blob(erv.data(), {static_cast<int64_t>(erv.size())}, torch::kFloat32);

        auto ifr = Factory::find_tn<IFieldResponse>(m_cfg.field_response);
        // Still at FR sampling, eg 100ns
        auto fr_fine = fr_average_tensor(ifr, m_cfg.plane_id);

        // We now must do two transformations.
        //
        // 1. Downsample FR from the fine FR sampling period to the coarse ER
        // sampling period.  This is done by removing N=1-fine/coarse FREQUENCY
        // bins near the Nyquist frequency.  That is, we wish to move the
        // Nyquist frequency from 1/fine to 1/coarse.
        //
        // 2. Pad both FR and ER time dimension in INTERVAL space to assure
        // linear convolution.
        //
        // I'm sure there is a trick to minimize the FFT round trips, but for
        // now we do the dumb straightforward thing.


        const int64_t er_nticks = er.size(0);
        const double er_period = ier->waveform_period();

        const int64_t fr_nticks_fine = fr_fine.size(1);
        const double fr_period_fine = ifr->field_response().period;

        // Find the half spectrum size of the downsampled FR.
        const int64_t fr_nticks_half = 0.5 * fr_nticks_fine * (fr_period_fine / er_period);

        // Chop it out 
        auto fr_spec = torch::fft::fft(fr_fine, {}, 1);
        fr_spec = fr_spec.narrow(1, fr_nticks_half, fr_nticks_fine - fr_nticks_half);
        auto fr = torch::real(torch::fft::ifft(fr_spec, {}, 1)); // coarse
        
        // This is the convolution target size.
        auto dft_nticks = linear_shape(fr.size(1), er.size(1));

        auto fr_pad = pad(fr, PadFuncOptions({0, 0, 0, dft_nticks - fr.size(1)}));
        auto er_pad = pad(er, PadFuncOptions({0, dft_nticks - fr.size(0)}));
                          
        auto FR = torch::fft::fft(fr_pad, {}, 1);
        auto ER = torch::fft::fft(er_pad).unsqueeze(0);
        auto FRER = FR*ER;      // "natural spectrum"

        // we got yet back to interval space again because we will need to pad
        // there to do the Fourier-space interpolation.
        m_response_waveform = torch::real(torch::fft::ifft(FRER, {}, 1));
    }


    size_t ResponseKernel::make_cache_key(const shape_t& shape) const {
        size_t h = 0;
        for (const auto& s : shape) {
            boost::hash_combine(h, s);
        }
        return h;
    }


    torch::Tensor ResponseKernel::spectrum(const std::vector<int64_t> & shape) const
    {
        torch::Tensor ten;
        raise<AssertionError>("not implemented");
        return ten;
    }

    std::vector<int64_t> ResponseKernel::shape() const {
        return vshape(m_natural_spectrum.sizes());
    }

}
