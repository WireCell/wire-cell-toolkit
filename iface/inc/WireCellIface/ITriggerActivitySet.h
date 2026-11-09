/** Trigger activities (TAs) and sets of them.

    A TA is a group of TPs from one detector unit, usually one plane, that a
    TA finder judged to be worth passing on.  Whether a TA keeps copies of its
    input TPs is a choice of its finder; ninputs counts them either way.

    Times are absolute in WCT units.  The set follows the stream contract
    described in ITriggerPrimitiveSet.h.  A TA may start before the coverage
    of the set it is in, by up to the latency of its finder.
 */

#ifndef WIRECELL_ITRIGGERACTIVITYSET
#define WIRECELL_ITRIGGERACTIVITYSET

#include "WireCellIface/ITriggerPrimitiveSet.h"
#include "WireCellIface/WirePlaneId.h"

namespace WireCell {

    struct TriggerActivity {
        double time_start{0}, time_end{0}, time_peak{0}, time_activity{0};
        int channel_start{-1}, channel_end{-1}, channel_peak{-1};
        uint64_t adc_integral{0};
        uint16_t adc_peak{0};
        int detector{-1};                  // anode ident
        WirePlaneId plane{kUnknownLayer};  // plane the TA was made on, if any
        int algorithm{0};                  // finder specific
        size_t ninputs{0};                 // number of TPs whether stored or not
        std::vector<TriggerPrimitive> inputs;
    };

    class ITriggerActivitySet : public IData<ITriggerActivitySet> {
       public:
        typedef std::vector<TriggerActivity> activities_t;

        virtual ~ITriggerActivitySet();

        virtual int ident() const = 0;
        /// The detector unit, or -1 if the set merges several.
        virtual int detector() const = 0;
        virtual double time_begin() const = 0;
        virtual double time_end() const = 0;
        /// The TAs, ordered by time_start.
        virtual const activities_t& activities() const = 0;
    };

}  // namespace WireCell

#endif
