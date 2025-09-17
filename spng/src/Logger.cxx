#include "WireCellSpng/Logger.h"
#include "WireCellSpng/Util.h"

#include <sstream>

namespace WireCell::SPNG {


    Logger::Logger(const std::string& logname, const std::string& pkgname)
        : Aux::Logger(logname, pkgname)
    {
    }

    void Logger::configure(const WireCell::Configuration& cfg)
    {
        m_verbose = get<int>(cfg, "verbose", m_verbose);
    }
    
    WireCell::Configuration Logger::default_configuration() const
    {
        Configuration cfg;
        cfg["verbose"] = m_verbose;
        return cfg;
    }
    
    void Logger::logit(const std::string& context) const
    {
        if (m_verbose == 0) return;

        log->debug("{}: call={}", context, m_count);
    }

    void Logger::logit(const ITorchTensorSet::pointer& ts, const std::string& context) const
    {
        if (m_verbose == 0) return;

        if (!ts) {
            log->debug("{}: call={}, null tensor set", context, m_count);
            return;
        }

        const size_t ntens = ts->tensors()->size();
        log->debug("{} call={}, {} tensors in set", context, m_count, ntens);

        if (m_verbose == 1) {
            return;
        }

        --m_verbose;
        for (const auto& iten : *ts->tensors()) {
            logit(iten, context);
        }
        ++m_verbose;
    }
    
    void Logger::logit(const ITorchTensor::pointer& ten, const std::string& context) const
    {
        if (m_verbose == 0) return;

        auto md = ten->metadata();
        std::string comma = "";
        std::stringstream ss;
        for (auto siz : ten->shape()) {
            ss << comma << siz << " ";
            comma = ", ";
        }


        std::string parent = get<std::string>(md, "parent", "");;
        if (parent.size()) {
            parent = "<--" + parent;
            
            log->debug("{} <{}> {} ({}) <{}> @{} {}",
                       context,
                       get<std::string>(md, "datatype", ""),
                       get<std::string>(md, "datapath", ""),
                       ss.str(),
                       ten->dtype(),
                       Torch::to_string(ten->device()),
                       parent);
        }

        if (m_verbose == 1) return;

        // fixme: some level 2 verbosity?

    }

    void Logger::logit(const TensorIndex& ti, const std::string& context) const
    {
        if (m_verbose == 0) return;
        
        log->debug("{} call={}, ident={}, {} tensors, {} parents (TDM)",
                   context, ti.ident(), m_count, ti.ntensors(), ti.nparents());

        if (m_verbose == 1) return;

        --m_verbose;
        for (const auto& node : ti.tree().depth()) {
            if (! node.value) { // skip empty root node
                continue;
            }
            logit(node.value, context);
        }
        ++m_verbose;

    }

}
