/**
   Drift planes given directly by configuration.

   The "volumes" parameter takes the same form as params.det.volumes in the
   pgrapher configuration: a list of objects, each with "wires" (the anode
   ident) and "faces" (as for AnodePlane).  The xregions are the distinct
   non-null faces in order of first appearance, as sim.drifter has made them.
 */

#ifndef WIRECELL_GEN_DRIFTPLANES
#define WIRECELL_GEN_DRIFTPLANES

#include "WireCellIface/IDriftPlanes.h"
#include "WireCellIface/IConfigurable.h"

#include <map>

namespace WireCell {
    namespace Gen {

        class DriftPlanes : public IDriftPlanes, public IConfigurable {
           public:
            DriftPlanes();
            virtual ~DriftPlanes();

            // IDriftPlanes
            virtual Configuration faces(int anode_ident) const;
            virtual Configuration xregions() const;

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
            std::map<int, Configuration> m_faces;
            Configuration m_xregions;
        };
    }  // namespace Gen
}  // namespace WireCell

#endif
