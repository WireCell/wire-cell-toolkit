/** Integer arithmetic of the DUNE DAQ trigger primitive generation (TPG).

    These are art-free, per-channel building blocks that reproduce the
    fixed-point arithmetic of the DAQ TPG (DUNE-DAQ tpglibs, "Naive"
    processors) and of the dunetrigger TPAlgTPC* tools that emulate it.  They
    keep their state in plain members so that a node may hold one instance per
    channel and carry it across successive frames when emulating a continuous
    stream.

    Everything here counts in samples.  Nothing depends on the sample period.

    Where tpglibs and dunetrigger differ, tpglibs is the default and the
    dunetrigger behavior is an option.  See trig/docs/dunetrigger-summary.org
    and trig/docs/timing.org.
 */

#ifndef WIRECELLTRIG_TPGMATH
#define WIRECELLTRIG_TPGMATH

#include <cstddef>
#include <cstdint>
#include <vector>

namespace WireCell::Trig {

    /// Approximate a/b with the reciprocal-multiply used by the AVX2 TPG
    /// (tpglibs _naive_div_int16).  The rounding is "unorthodox" and is
    /// reproduced exactly.
    int16_t avx2_divide(int16_t a, int16_t b);

    /// Return the most frequent value, the smallest one in case of a tie.
    /// Returns 0 for empty input.  This is how the dunetrigger tools initialize
    /// the pedestal of a channel from its whole waveform.  The DAQ starts from
    /// a pedestal of 0 and relies on the stream to converge.
    int16_t mode_pedestal(const std::vector<int16_t>& adcs);

    /// The "frugal streaming" median estimate of a pedestal.
    struct FrugalPedestal {
        int accum_limit{10};
        int16_t pedestal{0};
        int16_t accum{0};

        /// Update the estimate with one sample and return the new pedestal.
        int16_t operator()(int16_t sample);
    };

    /// The recursive running sum, RS = r*(RS/10) + s*(x/10), computed in 32
    /// bits and clamped at INT16_MAX as tpglibs does.  Below INT16_MIN tpglibs
    /// wraps while dunetrigger clamps; set saturate_low for the latter.
    struct RunningSum {
        int r_value{9};
        int scale_factor{5};
        bool saturate_low{false};
        int16_t sum{0};

        int16_t operator()(int16_t sample);
    };

    /// As RunningSum but of |x|.  Note, dunetrigger's AbsRunningSum tool also
    /// subtracts a second FrugalPedestal from the sum; compose that explicitly.
    struct AbsRunningSum : public RunningSum {
        AbsRunningSum() { r_value = 8; }
        int16_t operator()(int16_t sample);
    };

    /// A hit found by ThresholdHitter, in the terms of the current DAQ TP
    /// format.  sample_start counts samples in whatever absolute sample count
    /// the caller supplies.
    struct ThresholdHit {
        int64_t sample_start{0};             // first sample over threshold
        uint16_t samples_over_threshold{0};
        uint16_t samples_to_peak{0};         // from sample_start to the peak
        uint32_t adc_integral{0};
        uint16_t adc_peak{0};
    };

    /// Track a run of samples strictly over threshold and emit a hit when the
    /// run ends.  This is the hit logic shared by tpglibs and the dunetrigger
    /// tools.
    ///
    /// Unlike the dunetrigger tools, which drop a hit still open at the end of
    /// a waveform, the state persists across calls and flush() may close it.
    struct ThresholdHitter {
        int16_t threshold{0};

        /// Hits with fewer samples over threshold are dropped (tpglibs
        /// "sot_minima", configured per plane).  Zero keeps all.
        uint16_t min_samples_over_threshold{0};

        /// Feed one (pedestal subtracted, possibly filtered) sample at the
        /// absolute sample count.  Return true and fill hit if a hit just ended
        /// and passes the minimum length.
        bool operator()(int64_t sample_count, int16_t sample, ThresholdHit& hit);

        /// Close any open hit as if a sample under threshold arrived at
        /// sample_count.  Return true and fill hit if one passes.
        bool flush(int64_t sample_count, ThresholdHit& hit);

        /// True if a hit is in progress.
        bool open() const { return m_over; }

        /// Discard any hit in progress.
        void reset();

       private:
        bool m_over{false};
        uint16_t m_tover{0};
        uint16_t m_peak_time{0};
        uint16_t m_peak_adc{0};
        uint32_t m_charge{0};
    };

}  // namespace WireCell::Trig

#endif
