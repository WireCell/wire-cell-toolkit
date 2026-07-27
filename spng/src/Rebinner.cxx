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
        : TensorFilter("Rebinner", "spng") {}

    ITorchTensor::pointer Rebinner::filter_tensor(const ITorchTensor::pointer& in)
    {
        if (std::abs(m_config.factor) == 1) {
            logit(in, "pass through");
            return in;
        }

        auto tensor = in->tensor();

        tensor = m_rebin_func(tensor);

        return std::make_shared<SimpleTorchTensor>(tensor, in->metadata());
    }


    void Rebinner::configure(const WireCell::Configuration& config)
    {
        this->TensorFilter::configure(config);
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
        auto cfg = this->TensorFilter::default_configuration();
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
