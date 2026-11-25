#include "WireCellEdep/Ionization.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(Ionization, WireCell::Edep::Ionization, WireCell::INamed, WireCell::ISimTruthSegments)

using namespace WireCell;
using namespace WireCell::Edep;

Ionization::Ionization()
  : Aux::Logger("Ionization", "edep")
{
}

Ionization::~Ionization() {}

bool Ionization::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {  // end of stream
        return true;
    }
    out = in->segments();  // strip the segments; re-emit the same object, no copy
    const size_t n = (out && out->segments()) ? out->segments()->size() : 0;
    log->debug("stripped {} segments from sim truth {}", n, in->ident());
    return true;
}
