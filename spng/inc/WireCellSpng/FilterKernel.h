#ifndef WIRECELL_SPNG_FILTERKERNEL
#define WIRECELL_SPNG_FILTERKERNEL

#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchSpectrum.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFilterWaveform.h"
#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    /// The configuration for FilterKernel which accepts a JSON object of the
    /// same schema.
    struct FilterKernelConfig {

        /// Required string specifying the kind of filter.  Two are kinds of
        /// filters supported:
        ///
        /// - "highpass" :: pass high frequencies, attenuate low frequencies.
        ///                 (WCT mainline calls this "lf").
        ///
        /// The "highpass" filter is: 1 - exp(-(freq/scale)^power)
        ///
        /// - "lowpass" :: pass low frequencies, attenuate high frequencies.
        ///                (WCT mainline calls this "hf")
        /// The "lowpass" filter is: exp(-0.5 (freq/scale)^power).
        ///
        std::string kind{""};

        /// Required, give the interval-domain sampling period in
        /// interval-domain units.  This may be time units (eg 500ns) for time
        /// domain or unitless (eg 1.0) for wire/channel domain.  As usual, the
        /// sampling frequency is 1/period, Nyquist frequency is 0.5/period.
        double period=-1.0;

        /// Required, the "scale" parameter in units of frequency.
        ///
        /// The scale values for the various filters requires tuning for each
        /// detector and different values for induction and collection planes.
        double scale = -1.0;

        /// Optional, the "power" parameter, unitless.  Default is 2.0 for both
        /// kinds.
        double power = 2.0;

        /// Optional bool, default is true.  When true, zero out the "zero
        /// frequency bin" which effectively leaves the baseline unchanged.
        /// O.w. retain whatever value the underlying filter function provides.
        bool ignore_baseline = true;


    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::FilterKernelConfig,
                        kind,
                        period,
                        scale,
                        power,
                        ignore_baseline);


namespace WireCell::SPNG {
    /** The FilterKernel provides a 1D convolutional kernel intended to filter
     * other (de)convolutional parts.
     *
     * It returns 1D real-valued tensors.
     *
     * The filter spectra dimensions are (channel periodicity, temporal frequency).
     *
     * The 2D filters are constructed from two 1D filters combined via outer
     * product.  Each 1D filter is a sampled analytical function provided by
     * ITorchFilterWaveform.
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
    struct FilterKernel : public ContextBase, 
                         public ITorchSpectrum, virtual public IConfigurable {

        FilterKernel();
        FilterKernel(const FilterKernelConfig& cfg);
        virtual ~FilterKernel() = default;

        /// ITorchSpectrum methods.

        /// Return the filter spectrum.
        virtual torch::Tensor spectrum(const shape_t & shape) const;

        /// Filters are sampled from analytic functions which have no native
        /// size/shape.
        virtual shape_t shape() const { return shape_t{0,}; }

        // Deprecated.
        virtual shape_t shifts() const { return shape_t{0,}; }

        // IConfigurable - see FilterKernelConfig for configuration 
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
    private:

        FilterKernelConfig m_cfg;
        void configme();        // trigger application of m_cfg;

        using filter_function = std::function<float(float)>;
        filter_function m_func;


    };

}


#endif
