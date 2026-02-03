#include "WireCellSpng/Expand.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Torch.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGExpand,
                 WireCell::SPNG::Expand,
                 WireCell::SPNG::Expand::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);

using WireCell::HanaJsonCPP::from_json;
using WireCell::HanaJsonCPP::to_json;
using WireCell::String::has;

namespace WireCell::SPNG {

    Expand::Expand()
        : FanoutBase<ITorchTensor>("Expand")
    {
    }

    void Expand::fanout_separate(const input_pointer& in, output_vector& outv)
    {
        Configuration md = in->metadata();
        const int64_t nout = outv.size();
        auto inten = in->tensor();
        const int64_t nin = inten.size(m_config.dim);

        if (nout != nin) {
            raise<ValueError>("mismatch between dimension (%d) of size (%d) and configured (%d) fanout size",
                              m_config.dim, nin, nout);
        }

        auto outtens = m_op(inten);

        for (int64_t ind=0; ind<nout; ++ind) {
            outv[ind] = std::make_shared<SimpleTorchTensor>(outtens[ind], md);
        }
    }

    WireCell::Configuration Expand::default_configuration() const
    {
        auto cfg = FanoutBase<ITorchTensor>::default_configuration();
        update(cfg, HanaJsonCPP::to_json(m_config));
        return cfg;
    }

    void Expand::configure(const WireCell::Configuration& jconfig) {
        FanoutBase<ITorchTensor>::configure(jconfig);
        from_json(m_config, jconfig);

        // Functions requiring 'dim'

        if (m_config.operation == "unbind") {
            int captured_dim = m_config.dim;
            m_op = [captured_dim](const torch::Tensor& tensor) {
                return torch::unbind(tensor, captured_dim);
            };
        }

        // Functions ignoring 'dim'
        // none, yet....

        else {
            raise<ValueError>("unsupported reduce operation: %s", m_config.operation);
        }
    }


}
