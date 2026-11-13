#ifndef WIRECELL_IPHOTONHITSET
#define WIRECELL_IPHOTONHITSET

#include "WireCellIface/IPhotonHit.h"

namespace WireCell {

    /** An interface to the set of optical-photon hits of one event.
     */
    class IPhotonHitSet : public IData<IPhotonHitSet> {
       public:
        virtual ~IPhotonHitSet();

        /// An identifier unique to this set (typically the event number).
        virtual int ident() const = 0;

        /// The photon hits in this set.
        virtual IPhotonHit::shared_vector hits() const = 0;

        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
