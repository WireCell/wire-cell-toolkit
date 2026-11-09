/** A TC maker that looks for a burst of TAs, as from a supernova.

    TAs passing a selection are counted in a sliding time window.  When the
    count reaches a threshold a TC is made of the TAs in the window.  This is
    a counting experiment: with selected background TAs at rate b the count
    in a window T is Poisson with mean b*T, which sets the fake TC rate for a
    given threshold.

    Unlike the ported DAQ algorithms this one is native to WCT: it works on
    the WCT trigger objects in WCT units (no DAQ clock).

    Time ordering.  The window slides over the TA "time_start".  TAs need
    not arrive in time order: they are buffered and processed in order up to
    the "until" given to flush(), which promises no TA with an earlier time
    will come.  A TA arriving after that promise (late) is still counted if it
    falls in the current window, else it is dropped and counted as such.

    Since the count only rises when a TA enters the window, checking at each
    TA finds every crossing.  The TC has

    - time_candidate :: the time of the TA that reached the threshold.
    - time_start :: the earliest TA time_start in the window.
    - time_end :: the latest TA time_end in the window.
    - detectors :: of the TAs in the window, sorted and unique.
    - ninputs, inputs :: the TAs in the window.

    After a TC no other is made until time_candidate + holdoff.  The window
    keeps sliding, it is not cleared.
 */

#ifndef WIRECELLTRIG_BURSTMAKER
#define WIRECELLTRIG_BURSTMAKER

#include "WireCellIface/ITriggerCandidateSet.h"

#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <vector>

namespace WireCell::Trig {

    class TCMakerBurst {
       public:
        struct Selection {
            /// TA algorithms to count; empty counts all.
            std::vector<int> algorithms;
            /// Bounds on the number of TPs of a TA, inclusive.
            size_t min_inputs{0};
            size_t max_inputs{std::numeric_limits<size_t>::max()};
            /// Bounds on the TA ADC integral, inclusive.
            uint64_t min_adc_integral{0};
            uint64_t max_adc_integral{std::numeric_limits<uint64_t>::max()};
            /// Upper bounds on the TA extent, inclusive.  Channel span is
            /// channel_end - channel_start, time span time_end - time_start.
            int max_channel_span{std::numeric_limits<int>::max()};
            double max_time_span{std::numeric_limits<double>::infinity()};

            bool operator()(const TriggerActivity& ta) const;
        };

        struct Config {
            /// Sliding window length.  A TA at t counts with those in
            /// (t - window, t].
            double window{1e9};  // 1 s in WCT units
            /// Make a TC when this many selected TAs are in the window.
            size_t threshold{10};
            /// After a TC, make no other before time_candidate + holdoff.  A
            /// negative value means the window length.
            double holdoff{-1};
            /// The TC type (DAQ TriggerCandidateData::Type, default kSupernova).
            int tc_type{3};
            /// The TC algorithm (default DAQ TriggerCandidateData kCustom).
            int algorithm{8};
            Selection select;
        };

        using output_vector = std::vector<TriggerCandidate>;

        TCMakerBurst();
        explicit TCMakerBurst(const Config& cfg);

        /// Feed one TA.  It is buffered, TCs come from flush().
        void operator()(const TriggerActivity& ta, output_vector& out);

        /// Process the buffered TAs before "until", which promises no TA with
        /// an earlier time_start will come.  Infinity processes all.
        void flush(double until, output_vector& out);

        /// Forget all state.
        void reset();

        // Counters, for monitoring.
        size_t nseen() const { return m_nseen; }
        size_t nselected() const { return m_nselected; }
        size_t nlate() const { return m_nlate; }
        size_t ndropped() const { return m_ndropped; }

       private:
        Config m_cfg;
        std::multimap<double, TriggerActivity> m_pending;  // by time_start
        std::deque<TriggerActivity> m_window;              // sorted by time_start
        double m_processed{-std::numeric_limits<double>::infinity()};  // time up to which processed
        double m_head{-std::numeric_limits<double>::infinity()};       // latest time in window
        double m_armed{-std::numeric_limits<double>::infinity()};      // no TC before this
        size_t m_nseen{0}, m_nselected{0}, m_nlate{0}, m_ndropped{0};

        void enter(const TriggerActivity& ta, output_vector& out);
        void late(const TriggerActivity& ta, output_vector& out);
        void check(double now, output_vector& out);
    };

}  // namespace WireCell::Trig

#endif
