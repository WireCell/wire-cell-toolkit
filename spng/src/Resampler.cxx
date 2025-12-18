#include "WireCellSpng/Resampler.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/SimpleTorchTensor.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGResampler,
                 WireCell::SPNG::Resampler,
                 WireCell::ITorchTensorFilter,
                 WireCell::IConfigurable,
                 WireCell::INamed);


namespace WireCell::SPNG {

    Resampler::Resampler()
        : Logger("Resampler", "spng") {}


    bool Resampler::operator()(const input_pointer& in, output_pointer& out)
    {

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
        // Fixme: TDM MD handling still needs thought
        // md["datapath"] = md["datapath"].asString() + "/Resampler/" + get_name();

        auto tensor = in->tensor();

        tensor = resample_interval(tensor, 1.0, 1.0/m_config.ratio, m_config.dim, m_norm);

        out = std::make_shared<SimpleTorchTensor>(tensor, md);
        logit(out, "resampled");
        next_count();
        return true;
    }


    void Resampler::configure(const WireCell::Configuration& config)
    {
        // log->debug("config: {}", config);
        WireCell::configure_bases<Resampler, ContextBase, Logger>(this, config);
        HanaJsonCPP::from_json(m_config, config);

        m_norm = LMN::norm(m_config.norm);

        if (m_config.ratio > 1.0) {
            log->debug("will upsample by {}x (verbosity={})", m_config.ratio, verbosity());
        }
        else if (m_config.ratio < 1.0) {
            log->debug("will downsample by {}x (verbosity={})", 1.0/m_config.ratio, verbosity());
        }
        else {
            log->debug("not resampling, will pass through, I hope you know what you are doing");
        }
    }

    WireCell::Configuration Resampler::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<Resampler, ContextBase, Logger>(this);
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
