#include "WireCellSpng/TensorForwardTS.h"
#include "WireCellSpng/TorchScript.h"
#include "WireCellSpng/Util.h"

#include "WireCellSpng/HanaConfigurable.h"

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

        try {
            log->debug("Loading TorchScript model from: {} to device {}",
                       m_config.ts_filename, to_string(device()));
        
            m_module = torch::jit::load(m_config.ts_filename, device());
        }
        catch (const c10::Error& e) {
            log->critical("error loading model: \"{}\" to device \"{}\": {}",
                          m_config.ts_filename, to_string(device()), e.what());
            throw;
        }
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
        TorchSemaphore sem(context());

        std::vector<torch::IValue> value = { input };
        auto result = m_module.forward(value);
        return result.toTensor();
    }
}
