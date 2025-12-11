#include "WireCellSpng/TensorForwardTS.h"
#include "WireCellSpng/TorchScript.h"
#include "WireCellSpng/Util.h"

#include "WireCellSpng/HanaConfigurable.h"

#include "WireCellUtil/Persist.h"
#include "WireCellUtil/NamedFactory.h"


WIRECELL_FACTORY(SPNGTensorForwardTS,
                 WireCell::SPNG::TensorForwardTS,
                 WireCell::SPNG::ITensorForward,
                 WireCell::IConfigurable,
                 WireCell::INamed)


using namespace WireCell::HanaJsonCPP;                         


namespace WireCell::SPNG {

    TensorForwardTS::TensorForwardTS()
        : Logger("TensorForwardTS", "spng")
    {
    }

    void TensorForwardTS::configure(const WireCell::Configuration& cfg)
    {
        WireCell::configure_bases<TensorForwardTS, ContextBase, Logger>(this, cfg);
        from_json(m_config, cfg);

        auto path = Persist::resolve(m_config.ts_filename);
        if (path.empty()) {
            raise<ValueError>("failed to find Torch Script model file %s",
                              m_config.ts_filename);
        }
        log->debug("Loading TorchScript model from: {} to device {}",
                   path, to_string(device()));

        m_module = torch::jit::load(path, device()); // may throw
    }

    WireCell::Configuration TensorForwardTS::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<TensorForwardTS, ContextBase, Logger>(this);
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

    torch::Tensor TensorForwardTS::forward(const torch::Tensor& input) const
    {
        /// Do NOT assert this in a service.  Do it in a node.
        // TorchSemaphore sem(context());

        std::vector<torch::IValue> value = { input };
        auto result = m_module.forward(value);
        return result.toTensor();
    }
}
