/** Trigger candidates (TCs) and sets of them.

    A TC is a group of TAs, possibly from several detector units, that a TC
    finder judged to be a trigger candidate.  Whether a TC keeps copies of
    its input TAs is a choice of its finder; ninputs counts them either way.

    Times are absolute in WCT units.  The set follows the stream contract
    described in ITriggerPrimitiveSet.h.
 */

#ifndef WIRECELL_ITRIGGERCANDIDATESET
#define WIRECELL_ITRIGGERCANDIDATESET

#include "WireCellIface/ITriggerActivitySet.h"

namespace WireCell {

    struct TriggerCandidate {
        double time_start{0}, time_end{0}, time_candidate{0};
        std::vector<int> detectors;  // the detector units of the input TAs
        int type{0};                 // as DAQ TriggerCandidateData::Type
        int algorithm{0};            // finder specific
        size_t ninputs{0};           // number of TAs whether stored or not
        std::vector<TriggerActivity> inputs;
    };

    class ITriggerCandidateSet : public IData<ITriggerCandidateSet> {
       public:
        typedef std::vector<TriggerCandidate> candidates_t;

        virtual ~ITriggerCandidateSet();

        virtual int ident() const = 0;
        /// The detector unit, or -1 if the set merges several.
        virtual int detector() const = 0;
        virtual double time_begin() const = 0;
        virtual double time_end() const = 0;
        /// The TCs, ordered by time_start.
        virtual const candidates_t& candidates() const = 0;
    };

}  // namespace WireCell

#endif
