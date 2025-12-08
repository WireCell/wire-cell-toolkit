#include "WireCellSpng/TensorForwardTS.h"
#include "WireCellSpng/TorchScript.h"
#include "WireCellSpng/Util.h"
#include "WireCellSpng/SimpleTorchTensor.h"

#include "WireCellSpng/HanaConfigurable.h"

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

    ITorchTensor::pointer TensorForwardTS::forward(const ITorchTensor::pointer& input) const
    {
        if (! input) {
            return input;
        }
        TorchSemaphore sem(context());

        std::vector<torch::IValue> value = { to(input->tensor()) };
        auto result = m_module.forward(value);
        return std::make_shared<SimpleTorchTensor>(result.toTensor(),
                                                   input->metadata());        
    }
}
