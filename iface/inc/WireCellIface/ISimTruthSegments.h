#ifndef WIRECELL_ISIMTRUTHSEGMENTS
#define WIRECELL_ISIMTRUTHSEGMENTS

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/ISimTruth.h"
#include "WireCellIface/ITrackSegmentSet.h"

namespace WireCell {

    /** SimTruth -> TrackSegmentSet
     *
     *  A "stripper" that pulls the ionization track segments out of an
     *  ISimTruth and emits them for the rest of a chain (eg the
     *  ITrackSegmentSampler -> depos path).  Because ISimTruth returns its
     *  constituents as shared pointers to the same objects a standalone
     *  producer would emit, this can be a pass-through with no copy.
     *
     *  (This is the first of an anticipated family of ISimTruth stripper
     *  interfaces -- others would pull the primaries, trajectory tree or photon
     *  hits.)
     */
    class ISimTruthSegments : public IFunctionNode<ISimTruth, ITrackSegmentSet> {
       public:
        virtual ~ISimTruthSegments();

        typedef std::shared_ptr<ISimTruthSegments> pointer;

        virtual std::string signature() { return typeid(ISimTruthSegments).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };

}  // namespace WireCell

#endif
