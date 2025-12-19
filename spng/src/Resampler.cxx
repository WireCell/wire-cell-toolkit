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
        : TensorFilter("Resampler", "spng") {}


    ITorchTensor::pointer Resampler::filter_tensor(const ITorchTensor::pointer& in)
    {
        if (m_config.ratio == 1.0) {
            logit(in, "pass through");
            return in;
        }

        auto tensor = in->tensor();

        tensor = resample_interval(tensor, 1.0, 1.0/m_config.ratio, m_config.dim, m_norm);

        return std::make_shared<SimpleTorchTensor>(tensor, in->metadata());
    }


    void Resampler::configure(const WireCell::Configuration& config)
    {
        this->TensorFilter::configure(config);
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
        auto cfg = this->TensorFilter::default_configuration();
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
