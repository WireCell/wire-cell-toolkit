/** Find trigger candidates in TA sets.

    Algorithms:

    - "Burst" :: TCMakerBurst, the count of selected TAs in a sliding window
      reaching a threshold (see BurstMaker.h).  Native, in WCT units.

    - "Supernova" :: the DAQ triggeralgs TCMakerSupernova port, in DAQ units.
      TAs are converted with the DAQ clock (see DaqConvert.h) and the TCs
      made are converted back: their TAs lose the plane and the TPs.

    One core sees all TAs of the input, which may merge several detector
    units (see TASetFanin).  For each input set the TAs are fed, then the
    core is flushed up to the input coverage end less "lateness".  One TC set
    is output per input set, with the input coverage and detector.  A gap in
    coverage flushes and restarts the core.  On end of stream the core is
    flushed completely into a last, empty-coverage set.

    Configuration:

    - algorithm :: "Burst" (default) or "Supernova".
    - config :: the algorithm keys (see TrigConfig.h).  For "Burst" in WCT
      units, for "Supernova" in triggeralgs units.
    - store_inputs :: whether output TCs keep their TAs (default true).
    - lateness :: how far before the coverage end TAs may still come
      (default 0).  TA finders can emit a TA long after it started, so the
      upstream latency bounds this only for some of them.  "Burst" counts
      late TAs that still fall in its window and drops the others.
    - prescale, tc_type :: "Supernova" only, as triggeralgs.
    - clock_tick, epoch, period :: "Supernova" only, the DAQ clock (default
      16 ns, 0) and the sample period (500 ns).
 */

#ifndef WIRECELLTRIG_TCFINDER
#define WIRECELLTRIG_TCFINDER

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITriggerCandidateFinder.h"
#include "WireCellTrig/BurstMaker.h"
#include "WireCellTrig/DaqConvert.h"
#include "WireCellTrig/Makers.h"

#include <memory>

namespace WireCell::Trig {

    class TCFinder : public Aux::Logger, public ITriggerCandidateFinder, public IConfigurable {
       public:
        TCFinder();
        virtual ~TCFinder();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_queue& outq);

       private:
        std::string m_algorithm{"Burst"};
        Configuration m_alg_cfg;
        TCMaker::Config m_common;
        bool m_store_inputs{true};
        double m_lateness{0};
        double m_period;
        DaqClock m_clock;

        std::unique_ptr<TCMakerBurst> m_burst;
        std::unique_ptr<TCMaker> m_daq;

        double m_last_end{0};
        bool m_have_last{false};
        int m_detector{-1};
        int m_count{0};

        void make_core();
        void feed(const TriggerActivity& ta, std::vector<TriggerCandidate>& tcs);
        // until = infinity for all
        void flush(double until, std::vector<TriggerCandidate>& tcs);
    };

}  // namespace WireCell::Trig

#endif
