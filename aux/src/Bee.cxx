#include "WireCellAux/Bee.h"
#include "WireCellUtil/RaySamplers.h"

using namespace WireCell;

Aux::Bee::Object Aux::Bee::dump(IBlobSet::pointer bs, IBlobSampler::pointer sampler)
{
    return dump(bs->blobs(), sampler);
}

Aux::Bee::Object Aux::Bee::dump(const IBlob::vector& blobs, IBlobSampler::pointer sampler)
{
    Aux::Bee::Object bee;    

    for (const auto& iblob : blobs) {
        auto pc = sampler->sample_blob(iblob, iblob->ident());
        auto x = pc.get("x")->elements<double>();
        auto y = pc.get("y")->elements<double>();
        auto z = pc.get("z")->elements<double>();
        PointCloud::Array::span_t<double> q;
        auto qarr = pc.get("charge_val");
        if (qarr) {
            q = qarr->elements<double>();
        }

        const size_t npts = x.size();
        for (size_t ind=0; ind<npts; ++ind) {
            double charge = 0;
            if (q.size()) {
                charge = q[ind];
            }
            bee.append({x[ind], y[ind], z[ind]}, charge);
        }
    }
    return bee;
}