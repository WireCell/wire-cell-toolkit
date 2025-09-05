#ifndef WIRECELL_SPNG_FANBASE
#define WIRECELL_SPNG_FANBASE

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    class FanBase : public Aux::Logger, public IConfigurable {
    public:
        FanBase(const std::string& logname, const std::string& pkgnam="spng");
        virtual ~FanBase() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    protected:

        /// Emit standard log line for the state of the tensor set after combination
        virtual void maybe_log(const ITorchTensorSet::pointer& ts, const std::string context="") const;

        /// Configuration
        ///
        /// - multiplicity :: the fan-out multiplicity.
        int m_multiplicity{2};

        /// Configuration:
        ///
        /// - quiet ::  Set true to not call any logging.  Default is false.
        bool m_quiet{false};


        // Track calls.  Subclass must increment property
        mutable size_t m_count{0};
    };
}
#endif
