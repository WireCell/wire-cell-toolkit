/** The DAQ "Supernova" TA and TC makers.

    Ported from DUNE-DAQ triggeralgs develop (450867e)
    TAMakerSupernovaAlgorithm and TCMakerSupernovaAlgorithm.  The defaults
    reproduce them, including their hardcoded values, now configurable.

    TA: greedy clustering of successive TPs.  A TP joins the current cluster
    if it is close to it in time OR in channel (both=false, as the DAQ) or in
    time AND channel (both=true).  Otherwise the cluster is emitted and the TP
    starts a new one.  Every cluster becomes a TA, even of one TP.

    TC: count the TAs with more than hit_threshold TPs within the last
    time_window; when more than threshold are counted, emit a TC of them and
    start counting again.

    Deviations from triggeralgs:

    - The TA "no cluster yet" state is a flag, not time_start == 0.
    - TA flush(until) emits the open cluster only if there is one and no TP
      at or after "until" could join it, then resets.  With either-or
      merging that is only at end of stream (until = invalid_timestamp); with
      both, once until reaches time_end + time_tolerance.  triggeralgs emits
      it always (an empty TA at start) and does not reset.
    - The TC time_start (time_window before the TA) is clamped at zero
      instead of wrapping around.
    - Sample counts are converted to ticks with ticks_per_sample (DAQ: 32).
 */

#ifndef WIRECELLTRIG_SUPERNOVAMAKERS
#define WIRECELLTRIG_SUPERNOVAMAKERS

#include "WireCellTrig/Makers.h"

namespace WireCell::Trig {

    class TAMakerSupernova : public TAMaker {
       public:
        struct Config {
            /// Maximum time difference between TP and cluster edges, ticks.
            Daq::timestamp_diff_t time_tolerance{250};
            /// Maximum channel difference beyond the cluster edges, plus one.
            Daq::channel_diff_t channel_tolerance{2};
            /// Require time AND channel consistency to join a cluster.
            bool both{false};
            double ticks_per_sample{32};
        };

        TAMakerSupernova();
        explicit TAMakerSupernova(const Config& cfg);

       protected:
        virtual void process(const input_type& tp, output_vector& out);
        virtual void do_flush(Daq::timestamp_t until, output_vector& out);

       private:
        Config m_cfg;
        bool m_active{false};
        Daq::TriggerActivity m_ta;

        Daq::timestamp_t time_end_of(const input_type& tp) const;
        Daq::timestamp_t time_peak_of(const input_type& tp) const;
        bool is_time_consistent(const input_type& tp) const;
        bool is_channel_consistent(const input_type& tp) const;
        void start(const input_type& tp);
    };

    class TCMakerSupernova : public TCMaker {
       public:
        struct Config {
            /// Sliding window in which TAs are counted, ticks.
            Daq::timestamp_diff_t time_window{500'000'000};
            /// Emit a TC when more than this many TAs are counted.
            uint16_t threshold{3};
            /// Count only TAs with more than this many TPs.
            uint16_t hit_threshold{2};
        };

        TCMakerSupernova();
        explicit TCMakerSupernova(const Config& cfg);

       protected:
        virtual void process(const input_type& ta, output_vector& out);

       private:
        Config m_cfg;
        std::vector<Daq::TriggerActivityData> m_activity;
    };

}  // namespace WireCell::Trig

#endif
