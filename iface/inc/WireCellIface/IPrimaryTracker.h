#ifndef WIRECELL_IPRIMARYTRACKER
#define WIRECELL_IPRIMARYTRACKER

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IPrimaryVertexSet.h"
#include "WireCellIface/ISimTruth.h"

namespace WireCell {

    /** PrimaryVertexSet -> SimTruth
     *
     *  A primary tracker runs a tracking simulation over a set of primary
     *  vertices and produces the event's full simulation truth (ISimTruth):
     *  optionally the primaries, the trajectory tree, the ionization track
     *  segments, and the optical-photon hits.  Which constituents are populated
     *  is a property of the producer's configuration.
     *
     *  Downstream "stripper" nodes (eg ISimTruthSegments) then pull a single
     *  constituent out of the ISimTruth for the rest of a chain.
     */
    class IPrimaryTracker : public IFunctionNode<IPrimaryVertexSet, ISimTruth> {
       public:
        virtual ~IPrimaryTracker();

        typedef std::shared_ptr<IPrimaryTracker> pointer;

        virtual std::string signature() { return typeid(IPrimaryTracker).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };

}  // namespace WireCell

#endif
