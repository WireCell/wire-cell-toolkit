/** Per-channel trigger primitive generation (TPG) chain.

    A TPGChannel applies, sample by sample, the stages of the DAQ TPG
    (tpglibs) and of the dunetrigger TPAlgTPC* tools:

    1. pedestal subtraction by a frugal streaming median,
    2. an optional filter, the running sum or absolute running sum,
    3. an optional second frugal pedestal subtracted from the filter output
       (as dunetrigger's AbsRunningSum tool does),
    4. threshold hit finding.

    The state persists across calls to operator() so a waveform may be fed in
    pieces (eg successive frames of a stream).  It counts in samples and knows
    nothing about the sample period.  It uses no WCT interface types.
 */

#ifndef WIRECELLTRIG_TPGCHANNEL
#define WIRECELLTRIG_TPGCHANNEL

#include "WireCellTrig/TPGMath.h"

#include <cstdint>
#include <vector>

namespace WireCell::Trig {

    struct TPGConfig {
        enum class Filter { none, running_sum, abs_running_sum };
        enum class PedestalInit { zero, mode, fixed };

        /// Stage 1.  The DAQ starts the pedestal at zero; dunetrigger at the
        /// mode of the waveform it is given.
        int accum_limit{10};
        PedestalInit pedestal_init{PedestalInit::mode};
        int16_t pedestal{0};  // the value for PedestalInit::fixed

        /// Stage 2.
        Filter filter{Filter::none};
        int r_value{9};
        int scale_factor{5};
        bool saturate_low{false};  // dunetrigger running sum clamps low

        /// Stage 3.
        bool filter_pedestal{false};

        /// Stage 4.  A negative threshold disables TP finding (dunetrigger).
        int16_t threshold{0};
        uint16_t min_samples_over_threshold{0};
    };

    class TPGChannel {
       public:
        TPGChannel(const TPGConfig& cfg = TPGConfig{});

        /// Process samples, the first at absolute sample count "first".
        /// Hits are appended to "hits".  The pedestal is initialized from
        /// these samples on the first call after construction or reset() if
        /// configured to use the mode.
        void operator()(const std::vector<int16_t>& samples, int64_t first, std::vector<ThresholdHit>& hits);

        /// Run samples through stages 1-3 to settle the pedestals and filter
        /// without finding hits (dunetrigger "stabilisation_period").
        void settle(const std::vector<int16_t>& samples, size_t nsamples);

        /// Close any open hit as if the next sample, at absolute sample count
        /// "next", were under threshold.
        void flush(int64_t next, std::vector<ThresholdHit>& hits);

        /// Forget all state, eg on a gap in the stream.
        void reset();

        /// True until the first samples are seen after construction or reset.
        bool fresh() const { return m_fresh; }

        const TPGConfig& config() const { return m_cfg; }

       private:
        TPGConfig m_cfg;
        bool m_fresh{true};
        FrugalPedestal m_pedestal;
        RunningSum m_rs;
        AbsRunningSum m_ars;
        FrugalPedestal m_filter_pedestal;
        ThresholdHitter m_hitter;

        void initialize(const std::vector<int16_t>& samples);
        int16_t filter(int16_t sample);
    };

}  // namespace WireCell::Trig

#endif
