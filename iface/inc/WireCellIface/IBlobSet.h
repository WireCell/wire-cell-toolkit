/** A blob set holds a collection of blobs spanning a SINGLE time slice.
 *
 * Undefined behavior will occur if the blob set spans multiple time slices.
 *
 * See also ICluster which allows more rich associations.
 */

#ifndef WIRECELL_IBLOBSET
#define WIRECELL_IBLOBSET

#include "WireCellIface/IData.h"
#include "WireCellIface/IBlob.h"
#include "WireCellIface/ISlice.h"

namespace WireCell {

    class IBlobSet : public IData<IBlobSet> {
       public:
        virtual ~IBlobSet();

        /// Return some identifier number that is unique to this set.
        virtual int ident() const = 0;

        /// A slice relevant to this set.  This may be given even if there are
        /// no blobs (which have their own pointer to a slice).  It may also be
        /// nullptr.  Blobs in the set may span multiple slices.
        virtual ISlice::pointer slice() const = 0;

        /// Return the blobs in this set.  There is no ordering
        /// requirement.
        virtual IBlob::vector blobs() const = 0;

        /// Return a vector of the underlying IBlob::shape() in order of
        /// blobs(), (this method is just some sugar).
        virtual RayGrid::blobs_t shapes() const;
    };
}  // namespace WireCell

#endif
