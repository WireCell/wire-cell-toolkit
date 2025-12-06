/** A wrapper for pytorch torchscript
 */

#ifndef WIRECELLPYTORCH_TORCHSCRIPT
#define WIRECELLPYTORCH_TORCHSCRIPT

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/Logger.h"

#include "WireCellPytorch/Torch.h"  // One-stop header.

namespace WireCell::Pytorch {
    class TorchScript : public Aux::Logger,
                        public ITensorSetFilter,
                        public IConfigurable
    {
      public:
        TorchScript();
        virtual ~TorchScript() {}

        // IConfigurable interface
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // ITensorSetFilter interface
        virtual bool operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out);

      private:

        bool m_gpu{true};
        double m_wait_time{500*units::ms};

        torch::jit::script::Module m_module;
    };
}

#endif  // WIRECELLPYTORCH_TORCHSCRIPT
