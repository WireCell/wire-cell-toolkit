#include "WireCellSpng/Threshold.h"

#include "WireCellSpng/Threshold.h"
#include "WireCellSpng/Convo.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/TdmTools.h"
#include "WireCellSpng/HanaConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/NamedFactory.h"

#include "boost/container_hash/hash.hpp"


WIRECELL_FACTORY(SPNGThreshold,
                 WireCell::SPNG::Threshold,
                 WireCell::ITorchTensorFilter,
                 WireCell::IConfigurable)

using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {

    Threshold::Threshold()
        : TensorFilter("Threshold", "spng")
    {}

    Threshold::Threshold(const ThresholdConfig& cfg)
        : TensorFilter("Threshold", "spng")
        , m_config(cfg)
    {
    }

    WireCell::Configuration Threshold::default_configuration() const
    {
        auto cfg = this->TensorFilter::default_configuration();
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

    void Threshold::configure(const WireCell::Configuration& config)
    {
        this->TensorFilter::configure(config);
        from_json(m_config, config);

        log->debug("nominal={}, RMS: nsigma={} max={} axis={} binary={}",
                   m_config.nominal,
                   m_config.rms_nsigma,
                   m_config.rms_max_value,
                   m_config.rms_axis,
                   m_config.binary);
    }

    torch::Tensor Threshold::rms_threshold(torch::Tensor tensor)
    {
        auto clamped = tensor.clone();
        // rms_nsigma, rms_axis, rms_max_value
        if (m_config.rms_max_value > 0) {
            auto big = torch::abs(tensor) > m_config.rms_max_value;
            clamped.masked_fill_(big, 0);
        }            
        // keepdims, unbiased
        auto thresh = m_config.rms_nsigma*torch::std(clamped, m_config.rms_axis, true, true) + m_config.nominal;
        log->debug("threshold: rms of rms:{} thresh={}",
                   torch::std(thresh, 0, true, true).item<float>(), to_string(thresh));

        if (m_config.binary) {
            return tensor > thresh;
        }
        // Note, overwrite input.
        tensor.masked_fill_(tensor <= thresh, 0);
        return tensor;        
    }

    torch::Tensor Threshold::nominal_threshold(torch::Tensor tensor)
    {
        if (m_config.binary) {
            return tensor > m_config.nominal;
        }
        tensor.masked_fill_(tensor <= m_config.nominal, 0);
        return tensor;
    }

    ITorchTensor::pointer Threshold::filter_tensor(const ITorchTensor::pointer& in)
    {
        auto tensor = to(in->tensor());
        if (m_config.rms_nsigma == 0.0) {
            tensor = nominal_threshold(tensor);
        }
        else {
            tensor = rms_threshold(tensor);
        }

        return std::make_shared<SimpleTorchTensor>(tensor, in->metadata());
    }

}

