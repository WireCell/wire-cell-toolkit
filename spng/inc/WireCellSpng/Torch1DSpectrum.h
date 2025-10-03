/** This component provides field response data as read in from a "WCT
 * field response" JSON file */

#ifndef WIRECELLSPNG_TORCH1DSPECTRUM
#define WIRECELLSPNG_TORCH1DSPECTRUM

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Units.h"
#include "WireCellIface/IFilterWaveform.h"
#include <boost/compute/detail/lru_cache.hpp>

namespace WireCell::SPNG {

    class Torch1DSpectrum : public Logger,
                            public ContextBase,
                            public ITorchSpectrum,
                            virtual public IConfigurable {
    public:
        // Create directly with the JSON data file or delay that
        // for configuration.
        Torch1DSpectrum();

        virtual ~Torch1DSpectrum();

        // ITorchSpectrum
        virtual torch::Tensor spectrum() const;
        virtual torch::Tensor spectrum(const std::vector<int64_t> & shape);
        virtual std::vector<int64_t> shape() const { return m_shape; }

        // IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        /// Get any shifts of the response
        virtual std::vector<int64_t> shifts() const {return {};};

    private:
        torch::Tensor m_total_spectrum;
        std::vector<std::shared_ptr<IFilterWaveform>> m_spectra{};
        std::vector<torch::Tensor> m_spectrum_tensors{};
        boost::compute::detail::lru_cache<std::vector<int64_t>, torch::Tensor> m_cache;

        int64_t m_default_length = 0;


        bool m_debug_force_cpu = false;

        //List holding the type & name of spectra i.e. {HfFilter, Wiener_tight_V}
        std::vector<std::string> m_spectra_tns{};

        mutable std::vector<int64_t> m_shape;
    };
 
}  // namespace WireCell::SPNG
#endif
