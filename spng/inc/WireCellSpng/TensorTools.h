#ifndef WIRECELL_SPNG_TORCHTOOLS
#define WIRECELL_SPNG_TORCHTOOLS

#include "WireCellSpng/ITorchTensorSet.h"


namespace WireCell::SPNG {

}


#include "WireCellUtil/Configuration.h"


namespace WireCell::SPNG {

    // FIXME: Remove.  This already exists Configuration.h as update().  It also
    // is not generic but defines specific data model policy.
  void metadata_passthrough(
    const WireCell::Configuration & metadata_in,
    WireCell::Configuration & metadata_out,
    const Json::Value & passing_values);

    // FIXME:  Remove, eventually.  Not generic, obsoleted by TDM
    std::vector<torch::IValue> from_itensor(const ITorchTensorSet::pointer& in, bool is_gpu = false);
    ITorchTensorSet::pointer to_itensor(const std::vector<torch::IValue>& inputs);
}

#endif


