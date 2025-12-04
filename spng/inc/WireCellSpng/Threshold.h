#ifndef WIRECELL_SPNG_THRESHOLD
#define WIRECELL_SPNG_THRESHOLD

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchTensorFilter.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"

namespace WireCell::SPNG {

    struct ThresholdConfig {
        /// Nominal absolute threshold.
        ///
        /// This is in addition to an RMS-based threshold.
        ///
        /// Warning, this quantity is compared directly to the waveforms and so
        /// must be expressed in commensurate units.  WARNING: OSP sets this to
        /// 1.0 (no units).
        double nominal = 0.0;

        /// A multiple of the calcualted RMS to add to the nominal_threshold to
        /// form applied threshold.  If zero (default) RMS does not contribute
        /// to the threshold.  OSP defaults chooses 5.0 for collection and 3.0
        /// for induction.
        double rms_nsigma = 0.0;

        /// The dimension over which to calculate RMS.  Negative axis counts
        /// from end.  Default -1 will use RMS over time
        int rms_axis = -1;

        /// If RMS contributes to the threshold, a non-zero max value will
        /// exclude any larger samples.
        double rms_max_value = 0;

        /// If binary is true (default) the output tensor is of type bool and is
        /// value gives Boolean true for above threshold, false if at or below
        /// threshold.  If binary is false, values at or below threshold are set
        /// to zero and remaining values are left unchanged from input.
        bool binary = true;


        /// Note, OSP calculates RMS in a "creative" way.  We may want too add
        /// an option to invoke creativity.  For now, RMS is calculated in the
        /// simple torch::std().

    };

}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::ThresholdConfig, nominal, rms_nsigma, rms_axis, rms_max_value, binary);

namespace WireCell::SPNG {

    struct Threshold : public ContextBase,
                       public Logger,
                       public ITorchTensorFilter,
                       virtual public IConfigurable {
        Threshold();
        Threshold(const ThresholdConfig& cfg);
        virtual ~Threshold() = default;

        /// ITorchTensorFilter
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        // IConfigurable - see ThresholdConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // Non-API methods
        // Apply nominal and rms based threshold
        torch::Tensor rms_threshold(torch::Tensor tensor);
        // Apply just nominal.
        torch::Tensor nominal_threshold(torch::Tensor tensor);


    private:
        ThresholdConfig m_cfg;
    };

}
#endif
