#ifndef WIRECELL_SPNG_TENSORSETLOGGER
#define WIRECELL_SPNG_TENSORSETLOGGER

#include "WireCellAux/Logger.h"

namespace WireCell::SPNG {

    class TensorSetLogger : public Aux::Logger, public IConfigurable {
    public:

        TensorSetLogger(const std::string& logname, const std::string& pkgnam="spng");

        // IConfigurable

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    };
}

#endif
