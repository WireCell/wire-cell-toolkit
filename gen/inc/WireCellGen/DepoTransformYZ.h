/** DepoTransform with a per-plane YZ charge scale.

    Same as DepoTransform, but before a depo enters the simulation of
    a wire plane, its charge is multiplied by the value of an IYZMap
    for that plane at the depo's (y, z).  Each plane gets its own
    factor, at no extra simulation cost.  Use it to put a data-driven
    YZ non-uniformity (one map per plane and anode) into the
    simulation.

    Configuration: everything DepoTransform takes, plus

    - yzmap: type:name of an IYZMap service (required).
    - yzmap_key: the anode key looked up in the map.  Default: the
      "anode" type:name (e.g. "AnodePlane:anode0"), the convention of
      the YZMap files.

    configure() checks that the map has an entry, with a non-zero value
    at the centre of each face, for every plane of the anode.
 */

#ifndef WIRECELLGEN_DEPOTRANSFORMYZ
#define WIRECELLGEN_DEPOTRANSFORMYZ

#include "WireCellGen/DepoTransform.h"
#include "WireCellIface/IYZMap.h"

#include <string>

namespace WireCell::Gen {

    class DepoTransformYZ : public DepoTransform {
      public:
        DepoTransformYZ();
        virtual ~DepoTransformYZ();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        /// Scale the depo charge by the map value of this plane.
        virtual IDepo::pointer modify_depo(WirePlaneId wpid, IDepo::pointer depo);

      private:
        IYZMap::pointer m_yzmap{nullptr};
        std::string m_key;
    };

}  // namespace WireCell::Gen

#endif
