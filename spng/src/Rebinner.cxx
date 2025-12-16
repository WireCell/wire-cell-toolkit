#include "WireCellSpng/Rebinner.h"
#include "WireCellSpng/Rebin.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/SimpleTorchTensor.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGRebinner,
                 WireCell::SPNG::Rebinner,
                 WireCell::ITorchTensorFilter,
                 WireCell::IConfigurable,
                 WireCell::INamed);


namespace WireCell::SPNG {

    Rebinner::Rebinner()
        : Logger("Rebinner", "spng") {}


    bool Rebinner::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }
        if (std::abs(m_config.factor) == 1) {
            out = in;
            logit(in, "pass through");
            next_count();
            return true;
        }

        auto md = in->metadata();
        // Fixme: TDM MD handling still needs thought
        md["datapath"] = md["datapath"].asString() + "/Rebinner/" + get_name();
        auto tensor = in->tensor();

        tensor = m_rebin_func(tensor);

        out = std::make_shared<SimpleTorchTensor>(tensor, md);
        logit(out, "rebinned");
        next_count();
        return true;
    }


    void Rebinner::configure(const WireCell::Configuration& config)
    {
        // log->debug("config: {}", config);
        WireCell::configure_bases<Rebinner, ContextBase, Logger>(this, config);
        HanaJsonCPP::from_json(m_config, config);

        if (m_config.factor == 0) {
            raise<ValueError>("rebin factor of 0");
        }
        if (std::abs(m_config.factor) == 1) {
            log->debug("not resampling, will pass through, I hope you know what you are doing");
        }

        Rebin::Normalization mode = Rebin::parse_mode(m_config.norm);
        int64_t factor = m_config.factor;
        int64_t dim = m_config.dim;
        if (factor > 1) {
            m_rebin_func = [factor, dim, mode](const torch::Tensor& tensor) -> torch::Tensor {
                return Rebin::downsample(tensor, factor, dim, mode);
            };
        }
        else {
            factor = -factor;
            m_rebin_func = [factor, dim, mode](const torch::Tensor& tensor) -> torch::Tensor {
                return Rebin::upsample(tensor, factor, dim, mode);
            };
        }
    }

    WireCell::Configuration Rebinner::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<Rebinner, ContextBase, Logger>(this);
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
