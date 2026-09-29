/** Apply a TorchScript module to a whole tensor set: every tensor is one positional argument.

    TorchService (and Util.cxx's from_itensor / to_itensor) accept only 4-D float32 tensors and one output
    tensor.  A graph network needs integer edge lists, 1-D node arrays and several outputs, so this service:

    - passes the tensors of the input set, in order, as the positional arguments of the module's forward();
    - accepts any rank and the element types float, double, int32_t and int64_t ("f4", "f8", "i4", "i8");
      inputs are wrapped without a copy (torch::from_blob) on CPU, so the module must not modify its inputs;
    - returns a tensor, or a tuple or list of tensors, as an output set of the same element types, copied out.

    The output set carries the input set's ident and metadata.  TorchService is left untouched (wcfm doc 14).

    Configuration:
      model   (string) TorchScript file, resolved with WIRECELL_PATH
      device  (string) "cpu" (default), "gpu" or "gpuN"
      semaphore (string, optional) type:name of an ISemaphore; default "Semaphore:torch-<device>"
 */

#ifndef WIRECELLPYTORCH_TORCHTENSORSETSERVICE
#define WIRECELLPYTORCH_TORCHTENSORSETSERVICE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorForward.h"
#include "WireCellAux/Logger.h"
#include "WireCellPytorch/TorchContext.h"

#include "WireCellPytorch/Torch.h"  // One-stop header.

#include <vector>

namespace WireCell::Pytorch {

    /// The element type of an ITensor as a torch dtype; throws ValueError for types other than
    /// float, double, int32_t and int64_t.
    torch::Dtype torch_dtype_of(const ITensor::pointer& ten);

    /// Wrap every tensor of the set, in order, as an IValue on the device (no copy on CPU).
    std::vector<torch::IValue> tensorset_to_ivalues(const ITensorSet::pointer& in, const torch::Device& dev);

    /// A tensor, or a tuple / list of tensors, as a vector of ITensors (CPU copies, same element type).
    ITensor::vector ivalue_to_tensors(const torch::IValue& out);

    /// The service's forward on a given module (exposed for tests).
    ITensorSet::pointer forward_tensorset(torch::jit::script::Module& module, const ITensorSet::pointer& in,
                                          const torch::Device& dev);

    class TorchTensorSetService : public Aux::Logger, public ITensorForward, public IConfigurable {
      public:
        TorchTensorSetService();
        virtual ~TorchTensorSetService() {}

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        virtual ITensorSet::pointer forward(const ITensorSet::pointer& input) const;

      private:
        mutable torch::jit::script::Module m_module;
        TorchContext m_ctx;
        std::string m_model_path;
    };
}  // namespace WireCell::Pytorch

#endif  // WIRECELLPYTORCH_TORCHTENSORSETSERVICE
