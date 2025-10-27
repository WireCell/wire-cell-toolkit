#ifndef WIRECELL_SPNG_KERNELCONVOVLE
#define WIRECELL_SPNG_KERNELCONVOVLE

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/TorchLMN.h"
#include "WireCellSpng/ITorchTensorFilter.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellSpng/DFT.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"

namespace WireCell::SPNG {

    struct KernelConvolveAxisConfig {

        /// Subtract a baseline.  If true, the median value is found along this
        /// axis and subtracted.  Eg, if this axis is "time", each "channel" is
        /// median-subtracted.
        bool baseline = false;

        ///  Optional, (default false), if the convolution over this axis is
        ///  cyclic (true) or linear (false, default).  Note, if cyclic is true,
        ///  the dimension will not be padded at all, including for "faster FFT
        ///  size".  If false (linear) the dimension will be padded at least
        ///  large enough to assure linear convolution and potentially slightly
        ///  larger to hit a nearby "faster FFT size".
        bool cyclic = false;

        /// Specify how to pad this axis.
        ///
        /// Current modes are:
        ///
        /// - "head" prepend padding to the front of the time dimension.
        /// - "tail" append padding to the back of the time dimension (default).
        /// - "none" do NO padding.
        ///
        /// No padding can be reasonable if the kernel is a pure Fourier-space
        /// filter.
        ///
        /// Note, see roll_mode "decon" which is equivalent to "head' padding in
        /// time.
        std::string padding{"tail"};

        /// Optional, (default 0), crop the interval-space dimension.  The
        /// integer is interpreted in the following way:
        ///
        /// - 0 (default) does no cropping.
        ///
        /// - positive, retain this many low-side samples
        ///
        /// - -1, remove any high-side samples that were added to reach "faster FFT
        ///       size".  Note, if "faster" is false, this setting does not crop.
        ///
        /// - -2, remove all padding, leave dimension same size as input.
        ///
        /// Note: cropping a cyclic convolution likely is nonsensical.
        int crop = 0;

        /// Optional, apply a "roll" to the interval-space convolution result
        /// (before a crop).  This will cycle the tensor to move sample
        /// originally at index zero to index "roll".
        ///
        /// See section "Shifts" in the spng/docs/decon.org document for
        /// discussion on this.
        int roll = 0;

        /// In ADDITION to an explicit roll number, a set of canned roll modes
        /// can be applied.
        ///
        /// - "decon" :: This will roll the dimension by size of the kernel.  It
        /// can be be appropriate for the time dimension for SP decon.  When
        /// "tail" padding in time is used then a "decon" roll will cycle the
        /// "negative time" samples to the front of the time dimension.  This
        /// changes the time of the traces by the kernel size x sample period.
        /// This is equivalent to "head" padding in time.
        std::string roll_mode="";

    };

    struct KernelConvolveConfig {

        /// Required, the type/name of an ITorchSpectrum providing the kernel.
        std::string kernel{""};

        /// Optional.  A per-axis configuration.  If not given, defaults are
        /// used. See KernelConvolveAxisConfig for parameters that can be given.
        std::vector<KernelConvolveAxisConfig> axis{};

        /// Optional, default is false.  If true, seek a "faster fft size" for
        /// dimensions that are not cyclically convolved.  
        bool faster=false;

        /// Set the "tag" metadata attribute on the produced tensor.
        std::string tag="";

        /// Apply a format string to the metadata to produce a datapath for the
        /// output tensor.  If not provided (default) the produced tensor will
        /// retain the same datapath as the input tensor.  If provided, the
        /// input datapath is stored as the "derived_from" attribute.
        std::string datapath_format="";

        // fixme: tag and datapath_format should be uplifted to a base class
        // handling config and application.

        /// If set, save some debug output to the named file.  File is produced
        /// by torch::pickle_save().
        std::string debug_filename = "";
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::KernelConvolveAxisConfig,
                        baseline, cyclic, padding, crop, roll, roll_mode);
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::KernelConvolveConfig,
                        kernel, axis, faster, tag, datapath_format, debug_filename);

namespace WireCell::SPNG {
    /** Apply a 2D convolution with a kernel to input tensor.
     *
     * The kernel is configurable and fully determines the nature of the
     * convolution.  For the "WC signal processing deconvolution" configure this
     * node to use a DeconKernel.  See KernelConvolveConfig for this and other
     * configuration parameters.
     *
     * The shape of the input tensor may be either: (nbatches, nrows, ncols) or
     * (nrows, ncols) matching (nrows, ncols) of the strictly 2D kernel.
     */
    class KernelConvolve : public ContextBase,
                           public Logger,
                           public ITorchTensorFilter,
                           virtual public IConfigurable {

    public:
        
        KernelConvolve();
        KernelConvolve(const KernelConvolveConfig& cfg);
        virtual ~KernelConvolve() = default;

        /// ITorchTensorFilter
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        // IConfigurable - see KernelConvolveConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;


    private:

        void configme();        // called from configured constructor and configure()

        KernelConvolveConfig m_cfg;
        ITorchSpectrum::pointer m_kernel;
        FasterDftSize m_faster;

        std::vector<int> m_roll;
    };
}

#endif
