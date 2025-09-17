#include "WireCellSpng/ContextBase.h"

#include <iostream> // debug

namespace WireCell::SPNG {

    WireCell::Configuration ContextBase::default_configuration() const
    {
        Configuration cfg;
        cfg["device"] = "cpu";
        // "semaphore" is not specified by default
        return cfg;
    }

    void ContextBase::configure(const WireCell::Configuration& cfg)
    {
        auto devname = get<std::string>(cfg, "device", "cpu");
        auto semname = get<std::string>(cfg, "semaphore", "");
        m_ctx.connect(devname, semname);
    }

}
