#include "WireCellSpng/Resampler.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/SimpleTorchTensor.h"

namespace WireCell::SPNG {

    Resampler::Resampler()
        : Logger("Resampler", "spng") {}


    bool Resampler::operator()(const input_pointer& in, output_pointer& out) {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }
        if (m_config.ratio == 1.0) {
            out = in;
            logit(in, "pass through");
            next_count();
            return true;
        }

        auto md = in->metadata();
        auto tensor = in->tensor();

        out = std::make_shared<SimpleTorchTensor>(tensor, md);
        logit(out, "resampled");
        next_count();
        return true;
    }


    void Resampler::configure(const WireCell::Configuration& config) {
        WireCell::configure_bases<Resampler, ContextBase, Logger>(this, config);
        HanaJsonCPP::from_json(m_config, config);
        m_norm = LMN::norm(m_config.norm);
    }

    WireCell::Configuration Resampler::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<Resampler, ContextBase, Logger>(this);
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
