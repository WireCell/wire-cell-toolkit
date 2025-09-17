#ifndef WIRECELL_SPNG_FANBASE
#define WIRECELL_SPNG_FANBASE

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/Logger.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    class FanBase : public Logger {
    public:
        FanBase(const std::string& logname, const std::string& pkgnam="spng");
        virtual ~FanBase() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    protected:

        /// Configuration
        ///
        /// - multiplicity :: the fan-out multiplicity.
        size_t m_multiplicity{2};
    };
}
#endif
