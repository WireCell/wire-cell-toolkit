#include "WireCellSpng/TensorSetPickleSink.h"
#include "WireCellSpng/Torch.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTensorSetPickleSink,
                 WireCell::SPNG::TensorSetPickleSink,
                 WireCell::ITorchTensorSetSink,
                 WireCell::IConfigurable,
                 WireCell::INamed);


using WireCell::HanaJsonCPP::to_json;
using WireCell::HanaJsonCPP::from_json;

namespace WireCell::SPNG {

    TensorSetPickleSink::TensorSetPickleSink()
        : Logger("TensorSetPickleSink", "spng")
    {
    }

    void TensorSetPickleSink::configure(const WireCell::Configuration& jconfig)
    {
        this->Logger::configure(jconfig);
        HanaJsonCPP::from_json(m_config, jconfig);

        m_output = std::ofstream(m_config.filename, std::ios::binary);

    }

    WireCell::Configuration TensorSetPickleSink::default_configuration() const
    {
        auto cfg = this->Logger::default_configuration();
        auto cfg2 = to_json(m_config);
        return update(cfg, cfg2);
    }

    bool TensorSetPickleSink::operator()(const ITorchTensorSet::pointer& in)
    {
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }
        logit(in, "sinking");

        using tensor_map = torch::Dict<std::string, torch::Tensor>;
        tensor_map to_save;

        const auto& itensors = in->tensors();
        const size_t ntensors = itensors->size();

        for (size_t ind=0; ind<ntensors; ++ind) {

            // Give standard name that is unique even if tensor is not TDM compliant.
            std::string scount = std::to_string(get_count());
            std::string sind = std::to_string(ind);
            std::string name = "exec" + scount + "-tensor" + sind;

            // Embellish with TDM info
            auto iten = itensors->at(ind);
            auto md = iten->metadata();
            auto dt = get<std::string>(md, "datatype", "");
            auto dp = get<std::string>(md, "datapath", "");
            dp = String::replace(dp, "/", "_");
            if (! dt.empty()) {
                name += "-" + dt;
            }
            if (! dp.empty()) {
                name += "-" + dp;
            }
            to_save.insert(name, iten->tensor());

        }
        
        auto data = torch::pickle_save(to_save);
        m_output.write(data.data(), data.size());
        m_output.flush();

        next_count();
        return true;
    }


}
