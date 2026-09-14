#ifndef WIRECELL_ISIMTRUTH
#define WIRECELL_ISIMTRUTH

#include "WireCellIface/IData.h"
#include "WireCellIface/IPrimaryVertexSet.h"
#include "WireCellIface/ITrajectorySet.h"
#include "WireCellIface/ITrackSegmentSet.h"
#include "WireCellIface/IPhotonHitSet.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {

    /** An interface to the full truth of one tracking-simulation event -- the
     *  WCT model of an edep-sim TG4Event.
     *
     *  It aggregates the event's constituent products: the primaries handed to
     *  the simulation, the resulting particle trajectory tree, the ionization
     *  track segments, and the optical-photon hits.  Each constituent is
     *  returned as an IData shared pointer to the SAME object a standalone
     *  producer would emit, so a downstream node that wants only one part (eg
     *  the segments for the ionization chain) can consume an ISimTruth and
     *  re-emit that single constituent as-is, with no copy.  Any constituent
     *  may be null if the producer did not provide it.
     */
    class ISimTruth : public IData<ISimTruth> {
       public:
        virtual ~ISimTruth();

        /// An identifier unique to this event.
        virtual int ident() const = 0;

        /// The primary vertices/particles handed to the simulation.
        virtual IPrimaryVertexSet::pointer primaries() const = 0;

        /// The particle trajectory tree.
        virtual ITrajectorySet::pointer trajectories() const = 0;

        /// The ionization track segments.
        virtual ITrackSegmentSet::pointer segments() const = 0;

        /// The optical-photon hits.
        virtual IPhotonHitSet::pointer photons() const = 0;

        /// Optional event-level attributes (eg "run", "subrun").
        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
