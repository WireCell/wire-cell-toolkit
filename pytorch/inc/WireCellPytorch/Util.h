/** Utilities
 */

#ifndef WIRECELLPYTORCH_UTIL
#define WIRECELLPYTORCH_UTIL

#include "WireCellIface/ITensorSet.h"

#include "WireCellPytorch/Torch.h"  // One-stop header.

#include <string>

namespace WireCell {
    namespace Pytorch {

        std::string dump(const torch::Tensor &ten);

        ITensorSet::pointer to_itensor(const std::vector<torch::IValue> &inputs);
        std::vector<torch::IValue> from_itensor(const ITensorSet::pointer &inputs, const bool gpu = false);
        void write_torch_to_npy(const torch::Tensor &ten, const std::string &filename);
        std::string tensor_shape_string(const torch::Tensor& tensor);

    };  // namespace Pytorch
}  // namespace WireCell

#endif  // WIRECELLPYTORCH_UTIL
