/** Find trigger activities in TP sets with a ported DAQ TA algorithm.

    One algorithm core runs per wire plane (with an anode) or one for all
    TPs (without).  TPs are converted to DAQ units and the TAs made are
    converted back, see DaqConvert.h.

    For each input set the TPs are fed in time order, then every core is
    flushed up to the input coverage end less "latency".  One TA set is output
    per input set, with the input coverage.  A gap in coverage flushes and
    restarts the cores.  On end of stream the cores are flushed completely
    into a last, empty-coverage set.

    With streaming=false each input set stands alone, as for independent
    simulated frames: after its TPs the cores are flushed completely and
    restarted, and end of stream adds no set.

    Configuration:

    - algorithm :: "Supernova" or "SWIFT".
    - config :: the algorithm keys (see TrigConfig.h), in triggeralgs units,
      and planes :: an optional array of per-plane objects (by
      WirePlaneId::index()) overriding them.
    - prescale, max_samples_over_threshold :: as triggeralgs.
    - store_inputs :: whether output TAs keep their TPs (default true).
    - anode :: optional IAnodePlane type:name.
    - period :: the sample period (default 500 ns), giving ticks_per_sample.
    - clock_tick :: the DAQ clock tick (default 16 ns).
    - epoch :: the time of DAQ tick zero (default 0).
    - latency :: how far before the coverage end TPs may still come
      (default 0, use the TP finder latency when streaming).
    - streaming :: whether state carries across sets (default true).
 */

#ifndef WIRECELLTRIG_TAFINDER
#define WIRECELLTRIG_TAFINDER

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITriggerActivityFinder.h"
#include "WireCellTrig/DaqConvert.h"
#include "WireCellTrig/Makers.h"

#include <map>
#include <memory>

namespace WireCell::Trig {

    class TAFinder : public Aux::Logger, public ITriggerActivityFinder, public IConfigurable {
       public:
        TAFinder();
        virtual ~TAFinder();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_queue& outq);

       private:
        std::string m_algorithm{"Supernova"};
        Configuration m_alg_cfg;
        std::vector<Configuration> m_plane_cfgs;
        TAMaker::Config m_common;
        bool m_store_inputs{true};
        IAnodePlane::pointer m_anode;
        double m_period;
        DaqClock m_clock;
        double m_latency{0};
        bool m_streaming{true};

        struct Core {
            WirePlaneId wpid{kUnknownLayer};
            std::unique_ptr<TAMaker> maker;
        };
        std::map<int, Core> m_cores;  // by WirePlaneId::ident(), or 0
        double m_last_end{0};
        bool m_have_last{false};
        int m_detector{-1};
        int m_count{0};

        Core& core(int channel);
        void flush_cores(Daq::timestamp_t until, std::vector<TriggerActivity>& tas);
        void add_tas(const Core& core, std::vector<Daq::TriggerActivity>& daq_tas,
                     std::vector<TriggerActivity>& tas) const;
    };

}  // namespace WireCell::Trig

#endif
