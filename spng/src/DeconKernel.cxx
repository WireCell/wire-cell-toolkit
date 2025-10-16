#include "WireCellSpng/DeconKernel.h"

#include "WireCellUtil/HanaJsonCPP.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include "boost/container_hash/hash.hpp"


WIRECELL_FACTORY(SPNGDeconKernel,
                 WireCell::SPNG::DeconKernel,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {

    DeconKernel::DeconKernel()
        : Logger("DeconKernel", "spng")
        , m_cache(1)
    {}

    DeconKernel::DeconKernel(const DeconKernelConfig& cfg)
        : Logger("DeconKernel", "spng")
        , m_cfg(cfg)
        , m_cache(m_cfg.capacity)
    {
        configme();
    }

    WireCell::Configuration DeconKernel::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    void DeconKernel::configure(const WireCell::Configuration& config)
    {
        this->ContextBase::configure(config);
        // use the hana serializer
        from_json(m_cfg, config);
        configme();
    }

    void DeconKernel::configme()
    {
        m_cache.reset_capacity(m_cfg.capacity);

        m_filter = Factory::find_tn<ITorchSpectrum>(m_cfg.filter);
        m_response = Factory::find_tn<ITorchSpectrum>(m_cfg.response);

        if (m_cfg.debug_filename.size()) {
            log->debug("writing debug file: {}", m_cfg.debug_filename);
            write_debug(m_cfg.debug_filename);
        }
    }

    void DeconKernel::write_debug(const std::string& filename) const
    {
        /// There is no "native" size for the filter since it is a sampled
        /// analytical function.  Ultimately it is requested to be the padded
        /// size required for M*F/R.  Here we will kludge that size without
        /// going overboard..
        std::vector<int64_t> shape = {100,1000};
        
        using tensor_map = torch::Dict<std::string, torch::Tensor>;
        tensor_map to_save;

        auto kern = spectrum(shape);
        to_save.insert("decon_kernel", kern);

        // Some quick and dirty impulse functions.
        auto meas = torch::zeros(shape);
        meas.index({0, 0}) = 1.0;
        meas.index({25, 999}) = 1.0;
        meas.index({50, 100}) = 1.0;
        meas.index({75, 500}) = 1.0;
        meas = torch::fft::fft2(meas);
        auto sig = torch::real(torch::fft::ifft2(meas*kern));
        to_save.insert("decon_signal", sig);

        auto data = torch::pickle_save(to_save);
        std::ofstream output_file(filename, std::ios::binary);
        output_file.write(data.data(), data.size());
    }

    size_t DeconKernel::make_cache_key(const shape_t& shape) const {
        size_t h = 0;
        for (const auto& s : shape) {
            boost::hash_combine(h, s);
        }
        return h;
    }

    torch::Tensor DeconKernel::make_spectrum(const std::vector<int64_t> & shape) const
    {
        auto num = m_filter->spectrum(shape);
        auto den = m_response->spectrum(shape);

        // return num / den;
        // Instead, apply implicit filter to avoid divide-by-zero NaN.
        auto zeros = (den == 0);
        return torch::where(zeros, torch::zeros_like(num), num / den);
    }

    torch::Tensor DeconKernel::spectrum(const std::vector<int64_t> & shape) const
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

    std::vector<int64_t> DeconKernel::shape() const {
        return m_response->shape();
    }

    std::vector<int64_t> DeconKernel::shifts() const {
        std::vector<int64_t> ret = {0,0};
        return ret;
    }
}
