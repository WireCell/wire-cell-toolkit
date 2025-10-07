/**  Provide a monolithic kernel for filtered-response deconvolution.
 *
 */

#ifndef WIRECELL_SPNG_DECONKERNEL
#define WIRECELL_SPNG_DECONKERNEL

#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchSpectrum.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFilterWaveform.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellUtil/ThreadSafeCache.h"

namespace WireCell::SPNG {

    /// The configuration for DeconKernel which accepts a JSON object of the
    /// same schema.
    struct DeconKernelConfig {

        /// The type/name of the ITorchSpectrum providing the filter to
        /// (numerator of the kernel).  This typically combines a "wire" and
        /// "time" filters
        std::string filter{""};

        /// The type/name of the ITorchSpectrum providing the response
        /// (denominator of the kernel).  This is typically the FR*ER
        /// convolution.
        std::string response{""};

        /// Optional, the maximum number of spectra to hold in the LRU cache
        /// (not including the "natural" spectra).  By default it only caches 1
        /// providing last-cache.  Set larger for jobs that expect to see a
        /// variety of argument values passed to spectrum(shape).
        ///
        /// Note: as the DeconKernel is the combination of two other kernels, to
        /// reduce memory usage if you enable caching here you likely want to
        /// disable caching there or vice versa.  
        int capacity{1};
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::DeconKernelConfig,
                        filter,
                        response,
                        capacity);

namespace WireCell::SPNG {

    /** The DeconKernel provides ITorchSpectrum to deliver a single
     * CONVOLUTIONAL kernel to perform filtered-response DECONVOLUTION for an
     * anode wire plane.
     *
     * Note: currently DeconKernel assumes filter is a sampling of a
     * Fourier-domain function and thus a simple ratio is formed (in
     * Fourier-space) with the filter and response.  This is inherently linear
     * and not convolutional and thus no additional padding is considered.
     */
    struct DeconKernel : public ContextBase, 
                         public ITorchSpectrum, virtual public IConfigurable {

        DeconKernel();
        DeconKernel(const DeconKernelConfig& cfg);
        virtual ~DeconKernel() = default;

        /// ITorchSpectrum methods.

        /// Return the ratio of filter/response in the given shape.
        using shape_t = std::vector<int64_t>;
        virtual torch::Tensor spectrum(const shape_t & shape) const;

        /// Get the "natural" size of the kernel.  This returns the natural
        /// shape of the response.
        virtual std::vector<int64_t> shape() const;

        // deprecated
        virtual std::vector<int64_t> shifts() const;

        // IConfigurable - see DeconKernelConfig for configuration 
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
        // Non interface methods.  Unconditionally and without cache, produce
        // version of the natural spectrum at the given shape..
        torch::Tensor make_spectrum(const std::vector<int64_t> & shape) const;

        // return hash of a shape
        size_t make_cache_key(const shape_t& shape) const;

    private:

        DeconKernelConfig m_cfg;
        void configme();

        using tensor_cache_t = ThreadSafeCache<size_t, torch::Tensor>;
        mutable tensor_cache_t m_cache;

        ITorchSpectrum::pointer m_filter;
        ITorchSpectrum::pointer m_response;        

    };

}


#endif
