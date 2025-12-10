#include "WireCellSpng/TensorForward.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Util.h"


#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTensorForward,
                 WireCell::SPNG::TensorForward,
                 WireCell::ITorchTensorFilter,
                 WireCell::IConfigurable,
                 WireCell::INamed)


using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {

    TensorForward::TensorForward()
        : Logger("TensorForward", "spng")
    {
    }

    bool TensorForward::operator()(const input_pointer& in, output_pointer& out)
    {
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        logit(in, "forwarding");

        auto inten = in->tensor();
        auto md = in->metadata();

        const int ndims = inten.dim();
        const int64_t batch_dimension = m_config.batch_dimension;
        if (ndims == m_config.ndims) { // unbatched
            inten = inten.unsqueeze(batch_dimension);
        }

        log->debug("forwarding tensor: {}", to_string(inten));

        torch::Tensor result;

        std::vector<torch::Tensor> batches = {inten};
        if (m_config.nbatch) { // per-batch forwarding
            batches = torch::chunk(inten, batch_dimension, m_config.nbatch);
        }

        for (size_t bind = 0; bind<batches.size(); ++bind) {
            auto one = m_forward->forward(batches[bind]);
            batches.push_back(one);
            // fixme: add inner loop on chunks of time.
        }

        result = torch::cat(batches, batch_dimension);

        out = std::make_shared<SimpleTorchTensor>(result, md);

        logit(out, "forwarded");
        next_count();
        return true;
    }

    void TensorForward::configure(const WireCell::Configuration& config)
    {
        log->debug("configured with: {}", config);

        WireCell::configure_bases<TensorForward, ContextBase, Logger>(this, config);
        from_json(m_config, config);
        if (m_config.forward.empty()) {
            log->critical("no configured forward service in:{}", config);
            raise<ValueError>("TensorForward not given a forward service");
        }
        log->debug("using forward service: \"{}\"", m_config.forward);
        m_forward = Factory::find_tn<ITensorForward>(m_config.forward);
        
    }

    WireCell::Configuration TensorForward::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<TensorForward, ContextBase, Logger>(this);
        log->debug("default base config: {}", cfg);
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

}

