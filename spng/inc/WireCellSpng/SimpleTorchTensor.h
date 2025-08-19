#ifndef WIRECELL_SPNG_SIMPLETORCHTENSOR
#define WIRECELL_SPNG_SIMPLETORCHTENSOR 
#include "WireCellSpng/ITorchTensor.h"

namespace WireCell {
class SimpleTorchTensor: public ITorchTensor {
  public:
    SimpleTorchTensor(torch::Tensor tensor,
                const Configuration& md = Json::nullValue)
      : m_tensor(tensor), m_md(md) {
        m_tensor.requires_grad_(false); //Turn off the computational graph
      }

    virtual torch::Tensor tensor() const { return m_tensor.detach().clone(); }
    virtual Configuration metadata() const { return m_md; }
    virtual std::string dtype() const { return torch::toString(m_tensor.dtype()); }
    virtual std::vector<int64_t> shape() const { return m_tensor.sizes().vec(); }
    virtual torch::Device device() const { return m_tensor.device(); }

  private:

    torch::Tensor m_tensor;
    Configuration m_md;

};

}
#endif