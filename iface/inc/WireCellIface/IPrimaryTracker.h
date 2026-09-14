#ifndef WIRECELL_IPRIMARYTRACKER
#define WIRECELL_IPRIMARYTRACKER

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IPrimaryVertexSet.h"
#include "WireCellIface/ITrackSegmentSet.h"

namespace WireCell {

    /** PrimaryVertexSet -> TrackSegmentSet
     *
     *  A primary tracker runs (or otherwise realizes) a tracking simulation
     *  over a set of primary vertices and summarizes the resulting energy
     *  deposition as a set of track segments.  It is the primaries-in,
     *  segments-out counterpart to ITrackSegmentSampler (segments-in,
     *  depos-out), so the two compose into a full primaries -> depos chain.
     */
    class IPrimaryTracker : public IFunctionNode<IPrimaryVertexSet, ITrackSegmentSet> {
       public:
        virtual ~IPrimaryTracker();

        typedef std::shared_ptr<IPrimaryTracker> pointer;

        virtual std::string signature() { return typeid(IPrimaryTracker).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };

}  // namespace WireCell

#endif
