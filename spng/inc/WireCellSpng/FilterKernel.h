#ifndef WIRECELL_SPNG_FILTERKERNEL
#define WIRECELL_SPNG_FILTERKERNEL

#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ITorchSpectrum.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFilterWaveform.h"
#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    /// The configuration for FilterKernel which accepts a JSON object of the
    /// same schema.
    struct FilterKernelAxisConfig {

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
    struct FilterKernelConfig {

        /// The full N-D filter kernel is the outer product of a 1D filter
        /// kernels applied ALONG a dimension.  That is, for the 2D case,
        /// axis[0] describes the filter applied along rows at each column.
        ///
        /// Note: currently this is supported only up to 3D but support for
        /// higher dimensions can be added.
        std::vector<FilterKernelAxisConfig> axis{};

        /// If set, save some debug output to the named file.  File is produced
        /// by torch::pickle_save().
        std::string debug_filename = "";
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::FilterKernelAxisConfig,
                        kind,
                        period,
                        scale,
                        power,
                        ignore_baseline);
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::FilterKernelConfig, axis, debug_filename);


namespace WireCell::SPNG {
    /** The FilterKernel provides up to 3-D Fourier-space filter comprised of
     * outer-product of up to 3 1D filters.
     *
     * Higher dimensions can be supported with code changes. 
     *
     * It is a real-valued, Fourier-space sampling of analytical filters.  As
     * such, its shape should not incur any additional padding when used as part
     * of a linear convolution.
     *
     * This is typically used for decon in which case it should be configured
     * with two dimensions to produce shape:
     *
     *   (channel periodicity, time frequency) 
     * 
     */
    struct FilterKernel : public ContextBase, 
                          public Logger,
                          public ITorchSpectrum, virtual public IConfigurable {

        FilterKernel();
        FilterKernel(const FilterKernelConfig& cfg);
        virtual ~FilterKernel() = default;

        /// ITorchSpectrum methods.

        /// Return the filter spectrum.
        virtual torch::Tensor spectrum(const shape_t & shape) const;

        /// Filters are sampled from analytic functions which have no native
        /// size/shape.  This returns zeros.
        virtual shape_t shape() const { return shape_t(m_funcs.size(), 0); }

        // Deprecated.  Returns zeros.
        virtual shape_t shifts() const { return shape_t(m_funcs.size(), 0); }

        // IConfigurable - see FilterKernelConfig for configuration 
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // Non-API methods

        /// Get the 1D filter spectrum for the given axis.
        torch::Tensor spectrum_axis(const shape_t & shape, int64_t axis) const;

        // Write a debug output file.
        void write_debug(const std::string& filename) const;


    private:

        FilterKernelConfig m_cfg;
        void configme();        // trigger application of m_cfg;

        using filter_function = std::function<float(float)>;
        // Per-dimension filter functions
        std::vector<filter_function> m_funcs;


    };

}


#endif
