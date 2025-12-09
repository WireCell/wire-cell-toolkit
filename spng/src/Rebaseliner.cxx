#include "WireCellSpng/Rebaseliner.h"
#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/SimpleTorchTensor.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGRebaseliner,
                 WireCell::SPNG::Rebaseliner,
                 WireCell::ITorchTensorFilter,
                 WireCell::IConfigurable,
                 WireCell::INamed);


namespace WireCell::SPNG {

    Rebaseliner::Rebaseliner()
        : Logger("Rebaseliner", "spng") {}


    bool Rebaseliner::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        auto md = in->metadata();
        auto tensor = in->tensor();

        tensor = rebaseline(tensor, m_config.dim, m_config.threshold);

        out = std::make_shared<SimpleTorchTensor>(tensor, md);
        logit(out, "rebaselined");
        next_count();
        return true;
    }


    void Rebaseliner::configure(const WireCell::Configuration& config)
    {
        // log->debug("config: {}", config);
        WireCell::configure_bases<Rebaseliner, ContextBase, Logger>(this, config);
        HanaJsonCPP::from_json(m_config, config);
    }

    WireCell::Configuration Rebaseliner::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<Rebaseliner, ContextBase, Logger>(this);
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
