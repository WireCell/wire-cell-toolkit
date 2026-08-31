#include "WireCellAux/SimpleTrace.h"

using namespace WireCell;

Aux::SimpleTrace::SimpleTrace(int chid, int tbin, const ChargeSequence& charge)
  : m_chid(chid)
  , m_tbin(tbin)
  , m_charge(charge)
{
}
Aux::SimpleTrace::SimpleTrace(int chid, int tbin, ChargeSequence&& charge)
  : m_chid(chid)
  , m_tbin(tbin)
  , m_charge(std::move(charge))
{
}
Aux::SimpleTrace::SimpleTrace(int chid, int tbin, size_t ncharges)
  : m_chid(chid)
  , m_tbin(tbin)
  , m_charge(ncharges, 0.0)
{
}

int Aux::SimpleTrace::channel() const { return m_chid; }

int Aux::SimpleTrace::tbin() const { return m_tbin; }

const ITrace::ChargeSequence& Aux::SimpleTrace::charge() const { return m_charge; }

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
