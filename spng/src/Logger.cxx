#include "WireCellSpng/Logger.h"
#include "WireCellSpng/Util.h"

#include "WireCellUtil/Units.h"

#include <sstream>

// Used in formatting the log line prefix.
static
std::string centered(std::string s, size_t target, char space = ' ')
{
    for (size_t now = s.size(); now < target; ++now) {
        if (now%2) {
            s.insert(s.end(), space);
        }
        else {
            s.insert(s.begin(), space);
        }
    }
    return s;
}

namespace WireCell::SPNG {


    Logger::Logger(const std::string& group_name)
        : log(Log::logger(group_name))

    {
        m_hana_config.group_name = group_name;
    }

    Logger::Logger(const std::string& log_name, const std::string& group_name)
        : log(Log::logger(group_name))
    {
    }

    std::string Logger::get_name() const
    {
        if (m_interface_name.empty()) {
            return m_hana_config.group_name+m_hana_config.type_name;
        }
        return m_interface_name;
    }
    
    void Logger::set_name(const std::string& name)
    {
        if (name == m_interface_name) {
            return;
        }

        m_interface_name = name;

        std::string log_name = m_hana_config.group_name + "/" + m_hana_config.type_name + "/" + m_interface_name;

        // Make allow unique sinks so we can set unique pattern.
        // Note, this spawns an SPDLOG thread....
        log = Log::logger(log_name, false);

        // Set the pattern
        std::stringstream ss;
        ss << "[%H:%M:%S.%03e] %L [" << centered(m_hana_config.group_name, 8) << "] <" << m_hana_config.type_name << ":" << m_interface_name << "> %v %@";

        log->set_pattern(ss.str());

        log->debug("log name: \"{}\"", log_name);
    }

    
    void Logger::logit(const std::string& context) const
    {
        if (m_verbosity == 0) return;

        log->debug("{}: call={}", context, m_count);
    }

    void Logger::logit(const ITorchTensor::vector& tens, const std::string& context) const
    {
        if (m_verbosity == 0) return;

        const size_t ntens = tens.size();
        log->debug("{} call={}, {} tensors in set", context, m_count, ntens);

        if (m_verbosity == 1) {
            return;
        }

        --m_verbosity;
        for (const auto& iten : tens) {
            logit(iten, context);
        }
        ++m_verbosity;
    }

    void Logger::logit(const ITorchTensorSet::pointer& ts, const std::string& context) const
    {
        if (m_verbosity == 0) return;

        if (!ts) {
            log->debug("{}: call={}, null tensor set", context, m_count);
            return;
        }

        logit(*ts->tensors(), context);
    }
    
    void Logger::logit(const ITorchTensor::pointer& ten, const std::string& context) const
    {
        if (m_verbosity == 0) return;

        auto md = ten->metadata();
        std::string comma = "";
        std::stringstream ss;
        for (auto siz : ten->shape()) {
            ss << comma << siz;
            comma = ", ";
        }

        auto datatype = get<std::string>(md, "datatype", "");

        std::string parent = get<std::string>(md, "parent", "");;
        if (parent.size()) {
            parent = "<-- " + parent;
        }            
        log->debug("{} <{}> {} ({}) <{}> @{} {}",
                   context,
                   datatype,
                   get<std::string>(md, "datapath", ""),
                   ss.str(),
                   ten->dtype(),
                   to_string(ten->device()),
                   parent);
        if (datatype == "traces") {
            log->debug("\ttbin={}, time={}ms period={}ns tag:\"{}\"",
                       get(md, "tbin", 0),
                       get(md, "time", 0.0)/units::ms,
                       get(md, "period", 0.0)/units::ns,
                       get<std::string>(md, "tag", ""));
        }

        if (m_verbosity == 1) return;

        // fixme: some level 2 verbosity?

    }

    void Logger::logit(const TensorIndex& ti, const std::string& context) const
    {
        if (m_verbosity == 0) return;
        
        log->debug("{} call={}, ident={}, {} tensors, {} parents (TDM)",
                   context, ti.ident(), m_count, ti.ntensors(), ti.nparents());

        if (m_verbosity == 1) return;

        --m_verbosity;
        for (const auto& node : ti.tree().depth()) {
            if (! node.value) { // skip empty root node
                continue;
            }
            logit(node.value, context);
        }
        ++m_verbosity;

    }

    void Logger::set_verbosity(int level)
    {
        m_hana_config.verbosity = level;
    }
    int Logger::verbosity() const
    {
        return m_hana_config.verbosity;
    }
    void Logger::configured()
    {
        // we need this mutable.
        m_verbosity = m_hana_config.verbosity;
        debug("logger with verbosity: {}", m_verbosity);

    }

}
