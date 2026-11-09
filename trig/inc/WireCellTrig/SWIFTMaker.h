/** The DAQ "SWIFT" TA maker.

    Ported from DUNE-DAQ triggeralgs develop (450867e)
    TriggerActivityMakerSWIFT.

    TPs (kept if their peak OR their length passes a minimum) are binned in
    fixed time windows aligned to multiples of window_length.  When a window
    is complete its summed adc_integral decides: at least accept_threshold
    gives a TA, at least inspect_threshold gives a TA only if the dominant
    DBSCAN cluster of its TPs in (channel*wire_pitch, time*cm_per_tick) has
    more than cluster_energy_cut, else nothing.  The TA spans the window and
    holds all kept TPs of the window.

    Deviations from triggeralgs:

    - min_samples_over_threshold defaults to 8 samples.  triggeralgs keeps a
      default of 256, meant as 8 samples of 32 ticks in the old TP format, but
      compares it with a sample count.
    - epsilon is a real number.  triggeralgs parses it as an integer.
    - flush(until) closes the window only once until reaches its end, then
      starts an empty window after it.  triggeralgs closes it regardless and
      keeps it open.
 */

#ifndef WIRECELLTRIG_SWIFTMAKER
#define WIRECELLTRIG_SWIFTMAKER

#include "WireCellTrig/Makers.h"

namespace WireCell::Trig {

    class TAMakerSWIFT : public TAMaker {
       public:
        struct Config {
            /// Window length, ticks.
            uint64_t window_length{32000};
            /// Window summed adc_integral thresholds.
            uint64_t inspect_energy_threshold{15000};
            uint64_t accept_energy_threshold{55000};
            /// A TP is kept if adc_peak or samples_over_threshold reaches its minimum.
            uint16_t min_adc_peak{80};
            uint16_t min_samples_over_threshold{8};
            /// Clustering: distance per tick and per channel, both cm.
            float cm_per_tick{0.016f * 0.16f};
            float wire_pitch{0.48f};
            /// DBSCAN: minimum neighbourhood size (including the point) and radius (cm).
            int min_samples{2};
            float epsilon{2.0f};
            /// The dominant cluster of an inspected window must exceed this.
            uint64_t cluster_energy_cut{22000};
        };

        TAMakerSWIFT();
        explicit TAMakerSWIFT(const Config& cfg);

        /// The summed adc_integral of the largest DBSCAN cluster, exposed for tests.
        uint64_t dominant_cluster_energy(const std::vector<input_type>& tps, Daq::timestamp_t window_start) const;

       protected:
        virtual bool preprocess(const input_type& tp);
        virtual void process(const input_type& tp, output_vector& out);
        virtual void do_flush(Daq::timestamp_t until, output_vector& out);

       private:
        Config m_cfg;
        bool m_initialised{false};
        Daq::timestamp_t m_window_start{0};
        uint64_t m_window_energy{0};
        Daq::TriggerActivity m_ta;

        void reset_window(Daq::timestamp_t window_start);
        void close_window(output_vector& out);
        void set_ta_attributes();
    };

}  // namespace WireCell::Trig

#endif
