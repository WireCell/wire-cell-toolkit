/** An interface to a "forward" operator on a Torch set */

#ifndef WIRECELL_ITorchFORWARD
#define WIRECELL_ITorchFORWARD

#include "WireCellUtil/IComponent.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell::SPNG{
  class ITorchForward : public IComponent<ITorchForward> {
    public:
      virtual ~ITorchForward(){};

      virtual torch::Tensor forward(const torch::Tensor& input) const = 0;
  };
}  // namespace WireCell

#endif  // WIRECELL_ITorchFORWARD
