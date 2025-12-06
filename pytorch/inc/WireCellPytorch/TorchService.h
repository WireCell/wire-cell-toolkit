/** Apply a torch script module to "forward" an input tensor. */

#ifndef WIRECELLPYTORCH_TORCHSERVICE
#define WIRECELLPYTORCH_TORCHSERVICE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorForward.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"
#include "WireCellPytorch/TorchContext.h"

#include "WireCellPytorch/Torch.h"  // One-stop header.

namespace WireCell::Pytorch {
    class TorchService : public Aux::Logger,
                         public ITensorForward,
                         public IConfigurable
    {
      public:
        TorchService();
        virtual ~TorchService() {}

        // IConfigurable interface
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // ITensorSetFilter interface.  This is thread-safe.
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& input) const;

      private:

        // for read-only access, claim is that .forward() is thread
        // safe.  However .forward() is not const so we must make this
        // mutable.
        mutable torch::jit::script::Module m_module;

        // Even though thread safe, we want to honor a per device
        // semaphore to give user chance ot limit us.
        TorchContext m_ctx;

    };
}  // namespace WireCell::Pytorch

#endif  // WIRECELLPYTORCH_TORCHSERVICE
