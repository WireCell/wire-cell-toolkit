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
        : Logger("Threshold", "spng")
    {}

    Threshold::Threshold(const ThresholdConfig& cfg)
        : Logger("Threshold", "spng")
        , m_cfg(cfg)
    {
    }

    WireCell::Configuration Threshold::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<Threshold, ContextBase, Logger>(this);
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    void Threshold::configure(const WireCell::Configuration& config)
    {
        WireCell::configure_bases<Threshold, ContextBase, Logger>(this, config);
        from_json(m_cfg, config);
    }

    torch::Tensor Threshold::rms_threshold(torch::Tensor tensor)
    {
        // rms_nsigma, rms_axis, rms_max_value
        auto big = torch::abs(tensor) > m_cfg.rms_max_value;
        auto clamped = tensor.clone();
        clamped.masked_fill_(big, 0);
        // keepdims, unbiased
        auto thresh = m_cfg.rms_nsigma*torch::std(clamped, m_cfg.rms_axis, true, true) + m_cfg.nominal;

        if (m_cfg.binary) {
            return tensor > thresh;
        }
        // Note, overwrite input.
        tensor.masked_fill_(tensor <= thresh, 0);
        return tensor;        
    }

    torch::Tensor Threshold::nominal_threshold(torch::Tensor tensor)
    {
        if (m_cfg.binary) {
            return tensor > m_cfg.nominal;
        }
        tensor.masked_fill_(tensor <= m_cfg.nominal, 0);
        return tensor;
    }

    bool Threshold::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }
        logit(in, "input");

        TorchSemaphore sem(context());

        auto tensor = in->tensor();
        if (m_cfg.rms_nsigma == 0.0) {
            tensor = nominal_threshold(tensor);
        }
        else {
            tensor = rms_threshold(tensor);
        }
        out = std::make_shared<SimpleTorchTensor>(tensor, in->metadata());

        logit(out, "output");
        next_count();
        return true;
    }

}

