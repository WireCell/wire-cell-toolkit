#include "WireCellSpng/TensorForward.h"
#include "WireCellSpng/HanaConfigurable.h"

using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {

    TensorForward::TensorForward()
        : Logger("TensorForward", "spng")
    {
    }

    bool TensorForward::operator()(const input_pointer& in, output_pointer& out)
    {
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        logit(in, "forwarding");
        out = m_forward->forward(in);
        logit(out, "forwarded");
        next_count();
        return true;
    }

    void TensorForward::configure(const WireCell::Configuration& config)
    {
        WireCell::configure_bases<TensorForward, ContextBase, Logger>(this, config);
        from_json(m_config, config);
        if (m_config.forward.empty()) {
            log->critical("no configured forward service");
            raise<ValueError>("TensorForward not given a forward service");
        }
    }

    WireCell::Configuration TensorForward::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<TensorForward, ContextBase, Logger>(this);
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

}

