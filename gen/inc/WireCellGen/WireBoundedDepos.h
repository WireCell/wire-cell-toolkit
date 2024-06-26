#ifndef WIRECELL_GEN_WIREBOUNDEDDEPOS
#define WIRECELL_GEN_WIREBOUNDEDDEPOS

#include "WireCellIface/IDrifter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/Logger.h"

#include <tuple>
#include <vector>

namespace WireCell::Gen {

    /**
       WireBoundedDepos outputs depos based on which wires they
       "land".  To "land" means to drift antiparallel to the
       X-axis.  The set of wire regions on which depos may or may
       not "land" is configured as an array of objects with three
       integer attributes:

       [{
       plane: <plane-number>,
       min: <min-wire-number>,
       max: <max-wire-number,
       },...]

       Note the range is inclusive of the max.

       Wire numbers must be counted along the positive pitch
       direction and starting at 0.

       The filter can operate in "accept" or "reject" mode.  In
       "accept" mode, a depo which "lands" on any configured wire
       region will be output.  In "reject" mode a depo will not be
       output if it explicitly lands on any configured wire
       region.  The first wire region landed will determine the
       fate of the depo.  Default is "accept" mode.

       If users desire to bound depos by the intersection of wires
       from multiple planes they may pipeline multiple
       WireBoundDepos serially, each configured to accept or
       reject wire regions defined for a given plane.


    */
    class WireBoundedDepos : public Aux::Logger, public IDrifter, public IConfigurable {
    public:
        WireBoundedDepos();
        virtual ~WireBoundedDepos();
        virtual bool operator()(const input_pointer& depo, output_queue& outq);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

    private:
        IAnodePlane::pointer m_anode;
        bool m_accept{true};
        std::vector<const Pimpos*> m_pimpos;

        typedef std::tuple<int, int, int> wire_bounds_t;
        typedef std::vector<wire_bounds_t> wire_region_t;
        std::vector<wire_region_t> m_regions;
    };

}  // namespace WireCell::Gen

#endif
