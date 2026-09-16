#include "WireCellSpng/Rebaseliner.h"
#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Util.h"

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
        log->debug("Rebaselining Tensor {} with {} nonzero entries",
            to_string(tensor),
            tensor.count_nonzero().item<int64_t>()
        );
        tensor = rebaseline_zero(tensor,
                                 m_config.dim,
                                 m_config.consequtive_zeros,
                                 m_config.min_roi_size,
                                 m_config.shrink_size,
                                 m_config.remove_small,
                                 m_config.remove_negative);
        log->debug("Finished Rebaselining");
        return std::make_shared<SimpleTorchTensor>(tensor, in->metadata());
    }


    void Rebaseliner::configure(const WireCell::Configuration& config)
    {
        this->TensorFilter::configure(config);
        HanaJsonCPP::from_json(m_config, config);
        log->debug("Rebaselining with options: {}",
            m_config.as_str()
        );
    }

    WireCell::Configuration Rebaseliner::default_configuration() const
    {
        auto cfg = this->TensorFilter::default_configuration();
        auto cfg2 = HanaJsonCPP::to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }


}
