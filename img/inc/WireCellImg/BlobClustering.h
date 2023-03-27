/** Cluster blobs.

    This takes a stream of IBlobSets and mints a new ICluster on EOS
    or when the frame ident changes.

    An output ICluster will have the same ident number as the frame
    from which the input blobsets originated.

    The produced ICluster has b, w, c and s nodes with edges:

    - slice-blob
    - blob-wire
    - channel-wire
    - blob-blob

    The input IBlobSets need not be delivered in time order but the
    input stream must not mix IBlobSets between frame boundaries.

    The blob-blob edges will only be formed if their slices that have
    no large "gap" with the maximum gap given by the "spans"
    configuraiton parameter multiplied by the current slice being
    considered.

    Note, that input blob sets and thus their blobs may be held
    between calls (via their shared pointers).

 */

#ifndef WIRECELLIMG_BLOBCLUSTERING
#define WIRECELLIMG_BLOBCLUSTERING

#include "WireCellIface/IClustering.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IWire.h"
#include "WireCellIface/IChannel.h"
#include "WireCellIface/IBlob.h"
#include "WireCellIface/ISlice.h"
#include "WireCellIface/IBlobSet.h"

#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Img {

        class BlobClustering : public Aux::Logger, public IClustering, public IConfigurable {
           public:
            BlobClustering();
            virtual ~BlobClustering();

            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

            virtual bool operator()(const input_pointer& blobset, output_queue& clusters);

           private:
            // Define how many spans away from the current slice may
            // the next slice be and still be considered not to form a
            // "gap".  This is multiplier to the current slice so that
            // a value of 1.0 means that if the next slice does not
            // abut the current slice end time then it would indicate
            // a "gap"
            double m_spans{1.0};

            // Collect blob sets
            IBlobSet::vector m_cache;

            // flush graph to output queue
            void flush(output_queue& clusters);

            // Add the newbs to the graph.  Return true if a flush is needed (eg, because of a gap)
            bool graph_bs(const input_pointer& newbs);

            // Return true if newbs is from a new frame
            bool new_frame(const input_pointer& newbs) const;

            // Return the ident from the frame of the first cached blob or 0 if empty.
            int cur_ident() const;

            int m_count{0};

        };
    }  // namespace Img
}  // namespace WireCell

#endif
