// An API to interface with Bee that deals in IData.

#ifndef WIRECELLAUX_BEE
#define WIRECELLAUX_BEE

#include "WireCellUtil/Bee.h"
#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/IBlobSampler.h"

namespace WireCell::Aux::Bee {

    using WireCell::Bee::Points;

    // Convert a blob set to a Bee Points object by applying the sampler.  The
    // sampler is expected to produce point cloud arrays named: "x", "y", "z"
    // for point positions and may produce "charge_val" for a charge measure.
    Points dump(IBlobSet::pointer bs, IBlobSampler::pointer sampler);
    Points dump(const IBlob::vector& blobs, IBlobSampler::pointer sampler);

}


#endif
