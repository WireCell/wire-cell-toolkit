#ifndef WIRECELL_SPNG_ITORCHTENSOR
#define WIRECELL_SPNG_ITORCHTENSOR

#include <torch/torch.h>
#include "WireCellIface/IData.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {
class ITorchTensor : public IData<ITorchTensor> {
  public:
    virtual torch::Tensor tensor() const = 0;
    virtual Configuration metadata() const = 0;
    virtual std::string dtype() const = 0;
    virtual std::vector<int64_t> shape() const = 0;
    virtual torch::Device device() const = 0;
};
}

#endif