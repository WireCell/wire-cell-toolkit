#ifndef WIRECELL_ITRACKSEGMENTSAMPLER
#define WIRECELL_ITRACKSEGMENTSAMPLER

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/ITrackSegmentSet.h"
#include "WireCellIface/IDepoSet.h"

namespace WireCell {

    /** TrackSegmentSet -> DepoSet
     *
     *  A track segment sampler converts each track segment in a set
     *  into some number of point-like depositions, applying an
     *  ionization model to determine the number of electrons.
     */
    class ITrackSegmentSampler : public IFunctionNode<ITrackSegmentSet, IDepoSet> {
       public:
        virtual ~ITrackSegmentSampler();

        typedef std::shared_ptr<ITrackSegmentSampler> pointer;

        virtual std::string signature() { return typeid(ITrackSegmentSampler).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };

}  // namespace WireCell

#endif
