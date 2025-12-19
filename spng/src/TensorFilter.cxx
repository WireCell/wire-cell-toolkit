#include "WireCellSpng/TensorFilter.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/TdmTools.h"
#include "WireCellSpng/HanaConfigurable.h"

using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {

    TensorFilter::TensorFilter(const std::string& type_name, const std::string& group_name)
        : Logger(type_name, group_name)
    {
    }

    void TensorFilter::maybe_save(torch::Tensor ten, std::string name) {
        if (m_config.debug_filename.empty()) {
            return;
        }
        //if (!batched) {
        //  ten = ten.squeeze(0);
        //}
        m_to_save.insert(name, ten);
    }

    bool TensorFilter::operator()(const input_pointer& in, output_pointer& out)
    {
        m_to_save.clear();

        out = nullptr;
        if (! in) {
            logit("EOS");
            next_count();
            return true;
        }

        logit(in, "input");

        TorchSemaphore sem(context());

        out = filter_tensor(in);

        auto md = TDM::derive_metadata(in->metadata(), out->metadata(), m_config.datapath_format, m_config.tag);
        out = std::make_shared<SimpleTorchTensor>(out->tensor(), md);

        if (m_config.debug_filename.size() && m_to_save.size()) {
            std::string filename = fmt::format(m_config.debug_filename, fmt::arg("ident", get_count()));
            auto data = torch::pickle_save(m_to_save);
            std::ofstream output_file(filename, std::ios::binary);
            output_file.write(data.data(), data.size());
        }

        logit(out, "output");
        next_count();
        return true;
    }
    

    void TensorFilter::configure(const WireCell::Configuration& config)
    {
        WireCell::configure_bases<TensorFilter, ContextBase, Logger>(this, config);
        from_json(m_config, config);
    }
    
    WireCell::Configuration TensorFilter::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<TensorFilter, ContextBase, Logger>(this);
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }
    

}
