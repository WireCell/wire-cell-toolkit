#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/TensorIndex.h"

namespace WireCell::SPNG {

    FanBase::FanBase(const std::string& logname, const std::string& pkgname)
        : Aux::Logger(logname, pkgname)
    {
    }

    void FanBase::configure(const WireCell::Configuration& cfg)
    {
        m_multiplicity = get(cfg, "multiplicity", m_multiplicity);
        m_quiet = get<bool>(cfg, "quiet", m_quiet);

    }

    WireCell::Configuration FanBase::default_configuration() const
    {
        Configuration cfg;
        cfg["multiplicity"] = m_multiplicity;
        cfg["quiet"] = m_quiet;
        return cfg;
    }

    void FanBase::maybe_log(const ITorchTensorSet::pointer& ts, const std::string context) const
    {
        if (m_quiet) return;

        if (!ts) {
            log->debug("{}: call={}", context, m_count);
            return;
        }


        const size_t ntens = ts->tensors()->size();

        try {
            TensorIndex ti(ts); // maybe too expensive just for a line of log file...
            size_t nparents = ti.nparents();
            log->debug("{}: call={}, {} tensors, {} parents (TDM)", context, m_count, ntens, nparents);
        }
        catch (const ValueError& err) {
            log->debug("{}: call={}, {} tensors (not TDM)", context, m_count, ntens);
        }
    }


}
