#ifndef WIRECELL_ITRACKSEGMENTSET
#define WIRECELL_ITRACKSEGMENTSET

#include "WireCellIface/ITrackSegment.h"

namespace WireCell {

    /** An interface to information about a set of track segments.
     */
    class ITrackSegmentSet : public IData<ITrackSegmentSet> {
       public:
        virtual ~ITrackSegmentSet();

        /// Return some identifier number that is unique to this set.
        virtual int ident() const = 0;

        /// Return the track segments in this set.
        virtual ITrackSegment::shared_vector segments() const = 0;
    };

}  // namespace WireCell

#endif
