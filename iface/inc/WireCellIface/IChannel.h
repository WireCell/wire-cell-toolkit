/**
   IChannel embodies static information about a single front end
   electronics channel to which some number of wire segments in an
   conductor feeds.
 */

#ifndef WIRECELLIFACE_ICHANNEL
#define WIRECELLIFACE_ICHANNEL

#include "WireCellIface/IWire.h"
#include "WireCellIface/IData.h"

namespace WireCell {
    class IChannel : public IData<IChannel> {
       public:
        virtual ~IChannel();

        /// Detector-dependent, globally unique channel ID number.
        /// Negative is illegal, not guaranteed consecutive over all
        /// given channels.  This same number may be called simply
        /// "channel" in other contexts.
        virtual int ident() const = 0;

        /// The channel index is the "Wire Attachment Number".  The
        /// WAN counts along points of attachment of the zero'
        /// zero-segment wires for a wire plane.  Note, for detectors
        /// with wrapped wires this index also wraps.  It counts in
        /// the same direction as the WIP number of IWire::index but
        /// starts at zero even with wrapped wire detectors.
        virtual int index() const = 0;

        /// Wire segments ordered in increasing distance from channel input.
        virtual const IWire::vector& wires() const = 0;

        /// The ID of the plane of wire zero.  This is just sugar.
        virtual WirePlaneId planeid() const;

        /// Provide a global order number for this channel.
        virtual size_t global_order() const;

    };
}  // namespace WireCell

#endif
