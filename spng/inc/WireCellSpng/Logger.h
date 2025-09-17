#ifndef WIRECELL_SPNG_LOGGER
#define WIRECELL_SPNG_LOGGER

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/TensorIndex.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {

    /// This class provides an Aux::Logger with some "standard" log methods and
    /// configuration.
    class Logger : public Aux::Logger, public virtual IConfigurable {
    public:

        Logger(const std::string& logname, const std::string& pkgname="spng");
        virtual ~Logger() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    protected:

        /// Subclass call these to chirp about a tensor set in a standard way,
        /// and reactive to "verbose" configuration.
        virtual void logit(const std::string& context) const;
        virtual void logit(const ITorchTensorSet::pointer& ts, const std::string& context="") const;
        virtual void logit(const ITorchTensor::pointer& ten, const std::string& context="") const;
        virtual void logit(const TensorIndex& ti, const std::string& context="") const;

        /// Configuration:
        ///
        /// - verbose ::  integer.  0: no logging emitted, 1: one line (default), 2: multiple lines.
        // mutable as we -- and ++ for recursion.
        mutable int m_verbose{1};


        /// Track calls.  A subclass MUST increment this once per operator() call.
        mutable size_t m_count{0};
    };
}

#endif
