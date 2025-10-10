/** Provide a monolithic kernel for 2D detector response.
 *
 * This provides FR*ER and excludes long responses like RC.
 */

#ifndef WIRECELL_SPNG_RESPONSEKERNEL
#define WIRECELL_SPNG_RESPONSEKERNEL

#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ITorchSpectrum.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFilterWaveform.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellUtil/ThreadSafeCache.h"

namespace WireCell::SPNG {

    /// The configuration for ResponseKernel which accepts a JSON object of the
    /// same schema.
    struct ResponseKernelConfig {

        /// Required typename of an IFieldResponse.  This is assembled into a
        /// zero-sample-centered representation that has no artificial shift.
        /// See spng/docs/decon.org for details of what this means.
        std::string field_response{""};

        /// Required typename of an IWaveform giving the 1D electronics response.
        /// Sample 0 must represent t=0.
        std::string elec_response{""};

        /// Required plane index on which use from the larger field response.
        int plane_index{-1};

        /// Optional, the maximum number of spectra to hold in the LRU cache
        /// (not including the "natural" spectra).  By default it only caches 1
        /// providing last-cache.  Set larger for jobs that expect to see a
        /// variety of argument values passed to spectrum(shape).
        int capacity{1};

    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::ResponseKernelConfig,
                        field_response,
                        elec_response,
                        plane_index,
                        capacity);

namespace WireCell::SPNG {
    /** The ResponseKernel provides FR*ER spectrum
     *
     * Each spectrum tensor row spans the (time) frequency dimension and each
     * column spans the (channel) "periodicity" dimension.
     *
     * This component is thread safe and caches its tensors.
     */
    struct ResponseKernel : public ContextBase, public Logger,
                            public ITorchSpectrum, virtual public IConfigurable {

        ResponseKernel();
        ResponseKernel(const ResponseKernelConfig& cfg);
        virtual ~ResponseKernel() = default;

        /** The spectrum has no artificial shifts.
         *
         * It is formed by central padding the channel dimension and end padding
         * the time dimension in interval space prior to fft2d() to produce the
         * output.
         *
         * It is LRU-cached by shape.
         */
        virtual torch::Tensor spectrum(const shape_t & shape) const;
        virtual std::vector<int64_t> shape() const;

        /// Deprecated.
        virtual std::vector<int64_t> shifts() const;

        // IConfigurable - see ResponseKernelConfig for configuration 
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
        // Non interface methods.  Unconditionally and without cache, produce
        // version of the natural spectrum at the given shape..
        torch::Tensor make_spectrum(const std::vector<int64_t> & shape) const;

        // return hash of a shape
        size_t make_cache_key(const shape_t& shape) const;

    private:

        ResponseKernelConfig m_cfg;
        void configme();

        using tensor_cache_t = ThreadSafeCache<size_t, torch::Tensor>;
        mutable tensor_cache_t m_cache;

        // This is Fourier-space "natural" kernel, as calculated once in
        // configure() and for return by spectrum().
        torch::Tensor m_natural_spectrum;

        // This is natural FR*ER in interval space.  We keep this in order to
        // pad to match request in spectrum(shape).
        torch::Tensor m_response_waveform;

    };

}


#endif
