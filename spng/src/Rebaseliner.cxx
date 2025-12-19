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
        : TensorFilter("Rebaseliner", "spng") {}


    ITorchTensor::pointer Rebaseliner::filter_tensor(const ITorchTensor::pointer& in)
    {
        auto tensor = in->tensor();

        tensor = rebaseline(tensor, m_config.dim, m_config.threshold);

        return std::make_shared<SimpleTorchTensor>(tensor, in->metadata());
    }


    void Rebaseliner::configure(const WireCell::Configuration& config)
    {
        this->TensorFilter::configure(config);
        HanaJsonCPP::from_json(m_config, config);
    }

    WireCell::Configuration Rebaseliner::default_configuration() const
    {
        auto cfg = this->TensorFilter::default_configuration();
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
