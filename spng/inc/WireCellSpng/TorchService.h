/** Apply a torch script module to "forward" an input tensor. */

#ifndef WIRECELL_SPNG_TORCHSERVICE
#define WIRECELL_SPNG_TORCHSERVICE

#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchForward.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"
#include "WireCellSpng/TorchContext.h"

#include "WireCellSpng/Torch.h"  // One-stop header.

namespace WireCell::SPNG {
    class TorchService : public Aux::Logger,
                         public ITorchForward,
                         public IConfigurable
    {
      public:
        TorchService();
        virtual ~TorchService() {}

        // IConfigurable interface
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // ITensorSetFilter interface.  This is thread-safe.
        virtual ITorchTensorSet::pointer forward(const ITorchTensorSet::pointer& input) const;

      private:

        // for read-only access, claim is that .forward() is thread
        // safe.  However .forward() is not const so we must make this
        // mutable.
        mutable torch::jit::script::Module m_module;
        torch::Device m_device{torch::kCPU}; // Default to CPU

        // Even though thread safe, we want to honor a per device
        // semaphore to give user chance ot limit us.

    };
}  // namespace WireCell::SPNG

#endif  // WIRECELLSPNG_TORCHSERVICE
