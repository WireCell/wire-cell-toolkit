#ifndef WIRECELL_IPRIMARYVERTEXSET
#define WIRECELL_IPRIMARYVERTEXSET

#include "WireCellIface/IData.h"
#include "WireCellIface/IPrimaryVertex.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {

    /** An interface to a set of primary vertices: the complete primary
     *  kinematics of one event to be handed to a tracking simulation.
     *
     *  This is the WCT analog of edep-sim's TG4PrimaryVertexContainer (a
     *  std::vector<TG4PrimaryVertex>) and the natural input product to a
     *  tracking node (eg the edep/ "Ionization" and "ParticleTracking"
     *  nodes).
     */
    class IPrimaryVertexSet : public IData<IPrimaryVertexSet> {
       public:
        virtual ~IPrimaryVertexSet();

        /// An identifier unique to this set (typically the event number).
        /// A tracking producer may also use it to seed the simulation.
        virtual int ident() const = 0;

        /// The primary vertices in this set.
        virtual IPrimaryVertex::shared_vector vertices() const = 0;

        /// Optional event-level attributes (eg "run", "subrun").  An empty
        /// (null) Configuration if the producer provides none.
        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
