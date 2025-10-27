#include "WireCellSpng/ResponseKernel.h"
#include "WireCellSpng/Util.h"
#include "WireCellSpng/DeconTools.h"
#include "WireCellSpng/Convo.h"
#include "WireCellSpng/TorchLMN.h"

#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include "WireCellUtil/NamedFactory.h"

#include "boost/container_hash/hash.hpp"

WIRECELL_FACTORY(SPNGResponseKernel,
                 WireCell::SPNG::ResponseKernel,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

using namespace WireCell::HanaJsonCPP;                         

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
    void ResponseKernel::configme()
    {
        if (m_cfg.plane_index < 0 || m_cfg.plane_index>3) {
            raise<ValueError>("illegal plane ID: %d", m_cfg.plane_index);
        }
        m_cache.reset_capacity(m_cfg.capacity);

        if (m_cfg.scale == 0.0) {
            raise<ValueError>("bogus zero scale, check config");
        }

        /// for debugging we may save some tensors to file for checkout
        using tensor_map = torch::Dict<std::string, torch::Tensor>;
        tensor_map to_save;
        auto maybe_save = [&](torch::Tensor ten, std::string name) {
            if (m_cfg.debug_filename.size()) {
                to_save.insert(name, ten);
            }
        };


        log->debug("configuring FR*ER for convolve-then-downsample");
        auto ier = Factory::find_tn<IWaveform>(m_cfg.elec_response);
        auto ifr = Factory::find_tn<IFieldResponse>(m_cfg.field_response);

        // Some FRs have FP round off problems in their sample period.
        const double fr_period_fine = round(ifr->field_response().period/units::ns)*units::ns;

        // Sample ER at fr sampling period.
        Binning er_sampling(m_cfg.elec_duration/fr_period_fine, 0, m_cfg.elec_duration);
        std::vector<float> erv = ier->waveform_samples(er_sampling);
        auto er_fine = torch::tensor(erv, torch::kFloat32);
        maybe_save(er_fine, "fine_sampled_er");

        // Note, this full FR does not have symmetries expanded.
        maybe_save(fr_tensor(ifr->field_response(), m_cfg.plane_index), "fine_full_fr");

        // Average FR is symmetric
        auto fr_fine = fr_average_tensor(ifr->field_response(), m_cfg.plane_index);
        maybe_save(fr_fine, "fine_avg_fr");

        // Find size to assure linear convolution.
        const int64_t linear_size = linear_shape(fr_fine.size(1), er_fine.size(0));

        auto fr_pad = resize_tensor_tail(fr_fine, 1, linear_size);
        auto er_pad = resize_tensor_tail(er_fine, 0, linear_size);
                          
        auto FR = torch::fft::fft(fr_pad, {}, 1); // [current/electron]
        auto ER = torch::fft::fft(er_pad).unsqueeze(0); // [voltage/charge]
        auto FRER = FR*ER;      // The "natural spectrum".
        auto frer_fine = torch::real(torch::fft::ifft(FRER, {}, 1));
        maybe_save(frer_fine, "fine_frer");

        auto frer_coarse = LMN::resample_interval(frer_fine, fr_period_fine, m_cfg.period, 1);

        /// Convert to units of voltage/electron and apply user scale.
        frer_coarse = frer_coarse * (float)m_cfg.period * (float)m_cfg.scale;

        log->debug("er_period={}, scale={}, er_fine={}, fr_period={}, fr_fine={}, linear_size={}, fr_pad={}, er_pad={}, frer_fine={}, frer_coarse={}",
                   er_sampling.binsize(), m_cfg.scale, to_string(er_fine),
                   fr_period_fine, to_string(fr_fine),
                   linear_size,
                   to_string(fr_pad), to_string(er_pad),
                   to_string(frer_fine), to_string(frer_coarse));

        m_response_waveform = frer_coarse;

        if (m_cfg.debug_filename.size()) {
            to_save.insert("response_waveform", m_response_waveform);
            std::vector<int64_t> s = {100,2000};
            to_save.insert("padded_waveform", pad_waveform(s));
            to_save.insert("padded_spectrum", make_spectrum(s));
            auto data = torch::pickle_save(to_save);
            std::ofstream output_file(m_cfg.debug_filename, std::ios::binary);
            output_file.write(data.data(), data.size());
        }

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
        const int64_t caxis = 0;
        const int64_t cindex = middle_index(m_response_waveform.size(caxis));
        const int64_t csize = shape[caxis];
        auto tmp = resize_tensor(m_response_waveform, caxis, cindex, csize);

        const int64_t taxis = 1;
        const int64_t tsize = shape[taxis];
        return resize_tensor_tail(tmp, taxis, tsize);
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
