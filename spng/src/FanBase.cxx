#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/TensorIndex.h"

namespace WireCell::SPNG {

    FanBase::FanBase(const std::string& logname, const std::string& pkgname)
        : Logger(logname, pkgname)
    {
    }

    void FanBase::configure(const WireCell::Configuration& cfg)
    {
        this->Logger::configure(cfg);
        m_multiplicity = get(cfg, "multiplicity", (int)m_multiplicity);
    }

    WireCell::Configuration FanBase::default_configuration() const
    {
        Configuration cfg = this->Logger::default_configuration();
        cfg["multiplicity"] = (int)m_multiplicity;
        return cfg;
    }
}
