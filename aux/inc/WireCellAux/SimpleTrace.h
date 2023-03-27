#ifndef WIRECELLAUX_SIMPLETRACE
#define WIRECELLAUX_SIMPLETRACE

#include "WireCellIface/ITrace.h"

#include <map>
#include <vector>

namespace WireCell::Aux {

    /** This concrete trace is filled by time bin and charge.
     *
     * It provides the results of the filling such that the ChargeSequence
     * is trivially (exactly) zero suppressed but only at the ends.  Any
     * zeros bounded by non-zero charge are kept.
     */
    struct SimpleTrace : public ITrace {
        int m_chid;
        int m_tbin;
        ChargeSequence m_charge;

        SimpleTrace(int chid, int tbin, const ChargeSequence& charge);
        SimpleTrace(int chid, int tbin, size_t ncharges);

        // while someone has me as a concrete object they can modify.
        ChargeSequence& charge() { return m_charge; }

        // ITrace interface:

        virtual int channel() const;

        virtual int tbin() const;

        virtual const ChargeSequence& charge() const;
    };

}  // namespace WireCell::Aux

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
