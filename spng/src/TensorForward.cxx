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
        TorchSemaphore sem(this->context());

        auto inten = to(in->tensor());
        auto md = in->metadata();
        // md["datapath"] = md["datapath"].asString() + "/TensorForward/" + get_name();

        const int ndims = inten.dim();
        const int64_t batch_dimension = m_config.batch_dimension;

        if (ndims == m_config.ndims) {
            // Input is unbatched.  We normalize to batched for here
            inten = inten.unsqueeze(batch_dimension);
        }

        // nbatch=0 really means "1 batch but no batch dimension" so make 0 into
        // 1 so we can uniformly call chunk().  We'll debatch below when nbatch
        // is really 0.
        int nbatch = m_config.nbatch;
        if (!nbatch) nbatch = 1; 
        auto batches = torch::chunk(inten, nbatch, batch_dimension);
        size_t nbatches = batches.size();

        std::vector<torch::Tensor> results;
        for (size_t bind = 0; bind<nbatches; ++bind) {
            auto the_batch = batches[bind];
            if (m_config.nbatch == 0) {
                the_batch = the_batch.squeeze(batch_dimension);
            }

            log->debug("forwarding {} {}/{}: input:{}",
                       m_config.nbatch ? "batch" : "unbatched",
                       bind, nbatches, to_string(the_batch));

            auto one = m_forward->forward(the_batch);
            log->debug("forwarded batch {}/{}: output:{}",
                       bind, nbatches, to_string(one));
            results.push_back(one);
            // fixme: add inner loop on chunks of time.
        }

        torch::Tensor result;

        if (nbatches > 1) {
            result = torch::cat(results, batch_dimension);
        }
        else {
            result = results[0];
        }

        // If we added a batch dim only to make the model happy, remove it from
        // the output.
        if (m_config.nbatch == 0) {
            result = result.squeeze(batch_dimension);
        }

        out = std::make_shared<SimpleTorchTensor>(result, md);

        logit(out, "forwarded");
        next_count();
        return true;
    }

    void TensorForward::configure(const WireCell::Configuration& config)
    {
        WireCell::configure_bases<TensorForward, ContextBase, Logger>(this, config);
        from_json(m_config, config);
        if (m_config.forward.empty()) {
            log->critical("no configured forward service in:{}", config);
            raise<ValueError>("TensorForward not given a forward service");
        }
        log->debug("using forward service: \"{}\" nbatch={}, device={}",
                   m_config.forward, m_config.nbatch, to_string(device()));
        m_forward = Factory::find_tn<ITensorForward>(m_config.forward);
        
    }

    WireCell::Configuration TensorForward::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<TensorForward, ContextBase, Logger>(this);
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

}

