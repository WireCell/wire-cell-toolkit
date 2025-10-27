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
#include "WireCellUtil/Units.h"

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

        /// The duration of the electronics response sampling.  This should be
        /// chosen NOT LONGER than needed in order to sample where the ER is
        /// non-zero.  DO NOT set this to some ADC readout duration.  The ER
        /// will be sampled at the sampling period of the FR.  Note, IWaveform
        /// provides a default sampling and we ignore it.
        double elec_duration{20*units::us};

        /// The sample period for the kernel.  This must match the period of the
        /// target tensor to which this response will be applied.
        double period{500*units::ns};

        /// Required index of the plane to select from the full field response.
        int plane_index{-1};

        /// A multiplicative scale applied to the natural response tensor.
        ///
        /// The FR is originally in units of [current].  The ER is in units
        /// [voltage/charge] and its coarse sample period T is multiplied to the
        /// convolution of the two to give T(FR*ER) in units [voltage].  When
        /// the response is used to deconvolve an ADC waveform tensor, it is
        /// convenient to convert the response kernel from units of [voltage] to
        /// a unitless ADC [count].  This is done by giving a scale that is in
        /// units of [count/voltage].  In Jsonnet, and assuming the conventional
        /// "adc" object, this would be
        ///
        /// @code(jsonnet)
        /// adc.gain * (1 << adc.resolution) / (adc.fullscale[1] - adc.fullscale[0])
        /// @endcode
        ///
        /// Note, this leaves the kernel with a nominal baseline of 0 ADC count
        /// and thus expects the ADC waveform tensor to be baseline-subtracted.
        ///
        double scale{1.0};

        /// Specify how to pad to achieve requested spectrum shape.
        ///
        /// Current modes are:
        ///
        /// - "head" prepend padding to the front of the time dimension.  Relevant for decon.
        /// - "tail" append padding to the back of the time dimension.  Relevant for convo.
        ///
        /// This padding policy should be matched by whatever is (de)convolving
        /// with this kernel.  
        ///
        /// Currently, these options only control the time dimension.  Channel
        /// dimension is always tail-padded.
        std::string padding{"head"};

        /// Optional, the maximum number of spectra to hold in the LRU cache
        /// (not including the "natural" spectra).  By default it only caches 1
        /// providing last-cache.  Set larger for jobs that expect to see a
        /// variety of argument values passed to spectrum(shape).
        int capacity{1};

        /// If set, save some debug output to the named file.  File is produced
        /// by torch::pickle_save().
        std::string debug_filename = "";
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::ResponseKernelConfig,
                        field_response,
                        elec_response,
                        elec_duration,
                        period,
                        plane_index,
                        scale,
                        padding,
                        capacity,
                        debug_filename);

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
        // version of the response waveform at the given shape..
        torch::Tensor pad_waveform(const std::vector<int64_t> & shape) const;
        // Pad and apply fft2
        torch::Tensor make_spectrum(const std::vector<int64_t> & shape) const;

        // return hash of a shape
        size_t make_cache_key(const shape_t& shape) const;

    private:

        ResponseKernelConfig m_cfg;
        void configme();

        using tensor_cache_t = ThreadSafeCache<size_t, torch::Tensor>;
        mutable tensor_cache_t m_cache;

        // This is natural FR*ER in interval space.  We keep this in order to
        // pad to match request in spectrum(shape).
        torch::Tensor m_response_waveform;

    };

}


#endif
