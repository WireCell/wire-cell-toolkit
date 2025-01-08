/** Produce a PointCloud::Dataset representing a set of blobs.

    The interface places no restriction on the form of the
    representation.
 */

#ifndef WIRECELL_IBLOBSAMPLER
#define WIRECELL_IBLOBSAMPLER

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellIface/IBlob.h"

#include <tuple>

namespace WireCell {

    class IBlobSampler : public IComponent<IBlobSampler> {
      public:

        virtual ~IBlobSampler();
        
        // Return 0: a dataset of points that sample a blob.  The blob
        // index identifies the blob in the caller's context.
        // Return 1: Auxilliary data
        virtual std::tuple<PointCloud::Dataset, PointCloud::Dataset> sample_blob(const IBlob::pointer& blob,
            int blob_index=0) = 0;

    };

}

#endif


