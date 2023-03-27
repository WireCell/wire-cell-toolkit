/** Produce a PointCloud::Dataset representing a set of blobs.

    The interface places no restriction on the form of the
    representation.
 */

#ifndef WIRECELL_IBLOBSAMPLER
#define WIRECELL_IBLOBSAMPLER

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/PointCloud.h"
#include "WireCellIface/IBlob.h"

namespace WireCell {

    class IBlobSampler : public IComponent<IBlobSampler> {
      public:

        virtual ~IBlobSampler() = default;
        
        virtual PointCloud::Dataset sample_blobs(const IBlob::vector& blobs) = 0;

    };

}

#endif


