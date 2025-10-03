/**  Provide a monolithic kernel for filtered-response deconvolution.
 *
 */

#ifndef WIRECELL_SPNG_DECONKERNEL
#define WIRECELL_SPNG_DECONKERNEL

#include "WireCellSpng/Logger.h"
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

        /// Required typename of an IFieldResponse.  This is assembled into a
        /// zero-sample-centered representation that has no artificial shift.
        /// See spng/docs/decon.org for details of what this means.
        std::string field_response{""};

        /// Required typename of an IWaveform giving the 1D electronics response.
        /// Sample 0 must represent t=0.
        std::string elec_response{""};

        /// Required typename of a IFilterWaveform giving 1D filter function in
        /// time-Fourier (frequency) space.
        std::string time_filter{""};

        /// Required typename of a IFilterWaveform giving 1D filter function in
        /// channel-Fourier (periodicity) space.
        std::string channel_filter{""};
        
        /// Required plane ID on which use from the larger field response.
        int plane_id{-1};

        /// Optional, the maximum number of spectra to hold in the LRU cache
        /// (not including the "natural" spectra).  By default it only caches 1
        /// providing last-cache.  Set larger for jobs that expect to see a
        /// variety of argument values passed to spectrum(shape).
        int capacity{1};

    };

    /** The DeconKernel provides a single convolutional kernel to perform
     * filtered-response deconvolution for an anode wire plane.
     *
     * Each spectrum tensor row spans the (time) frequency dimension and each
     * column spans the (channel) "periodicity" dimension.
     *
     * The spectrum is formed as a ratio of filter/response in Fourier-space.
     *
     * The filter is the outer product of two 1D filters, one for each
     * dimension.  Filters are sampled analytic functions on the Fourier space
     * domain.
     *
     * The response is the convolution of 2D field and 1D electronics responses.
     *
     * FIXME: adhere to the fixed FIXME's in ITorchSpectrum when they are fixed.
     *
     * This component is thread safe and caches its tensors.
     *
     */
    struct DeconKernel : public ContextBase, public Logger,
                         public ITorchSpectrum, virtual public IConfigurable {

        DeconKernel();
        virtual ~DeconKernel() = default;

        /// ITorchSpectrum methods.

        /// Return the "natural" filter/response kernel.
        ///
        /// "Natural" here means the filter/response kernel as constructed from
        /// the configured filters and responses.  It is of minimal size
        /// required so that FR*ER is a linear convolution.
        ///
        /// The response is zero-sample-centered and so has no "artificial"
        /// shifts.  See spectrum(shape) for discussion of padding based shifts.
        virtual torch::Tensor spectrum() const;

        /// Return the kernel in a specific shape.
        ///
        /// Conceptually, the returned spectrum is an interpolation in Fourier
        /// space of the "natural" spectrum.  The exact nature of the
        /// interpolation is relevant the interpretation of the result.  It
        /// preserves the "no artificial shift" of the natural spectrum by
        /// interpolating in the following manner:
        ///
        /// 1. Central-pad the channel dimension of the interval-space
        /// representation of the natural spectrum.
        ///
        /// 2. End-pad the time dimension of the interval-space representation
        /// of the natural spectrum.
        ///
        /// 3. Apply fft2().  (FIXME: once we clarify how to represent full vs
        /// half spectra we may have rfft()+fft() optimization).
        ///
        using shape_t = std::vector<int64_t>;
        virtual torch::Tensor spectrum(const shape_t & shape) const;

        /// Get the "natural" size of the response.  FIXME: is this not
        /// redundant with spectrum().sizes()?
        virtual std::vector<int64_t> shape() const = 0;

        /// Get any shifts of the response.  FIXME: what is this supposed to
        /// mean?
        ///
        /// For now, it will return the size of the response denominator as this
        /// allows the user to torch::roll() the time dimension of a convolution
        /// with this spectrum so that earliest possible sample is in column
        /// zero, latest possible is in the last column.
        virtual std::vector<int64_t> shifts() const = 0;

        // IConfigurable - see DeconKernelConfig for configuration 
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
    private:

        DeconKernelConfig m_cfg;

        using tensor_cache_t = ThreadSafeCache<size_t, torch::Tensor>;
        tensor_cache_t m_cache;

        // return hash of a shape
        size_t make_cache_key(const shape_t& shape) const;

        // This is Fourier-space "natural" kernel, as calculated once in
        // configure() and for return by spectrum().
        torch::Tensor m_natural_spectrum;

        // This is natural FR*ER in interval space.  We keep this in order to
        // pad to match request in spectrum(shape).
        torch::Tensor m_response_waveform;

        // The filter for each dimension.  We keep these original (not tensor)
        // as we must sample them differently for the "natural" and custom
        // spectra shapes.
        IFilterWaveform::pointer m_time_filter;
        IFilterWaveform::pointer m_channel_filter;        


    };

}


#endif
