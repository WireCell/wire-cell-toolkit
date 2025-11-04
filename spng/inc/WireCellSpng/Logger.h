#ifndef WIRECELL_SPNG_LOGGER
#define WIRECELL_SPNG_LOGGER

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/TensorIndex.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellIface/INamed.h"
#include "WireCellUtil/Logging.h"

namespace WireCell::SPNG {

    /** @brief Logging configuration
     */
    struct LoggerConfig {
        /// The name for the underlying logger, should be set to component name.
        std::string type_name = "SPNG"; 
        /// the name for the logger group, should be set to package name, okay to leave as-is.
        std::string group_name = "spng"; 
        /// These names should NOT be instance names.  See set_name() and
        /// get_name() in the Logger class.  

        // The verbosity level.  See comments in Logger class.
        int verbosity = 0;
        
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::LoggerConfig, type_name, group_name, verbosity);

namespace WireCell::SPNG {
    /** @brief Logging functionality.

        This may be used as a base class to provide a "log" data member or the
        class may be used as a class member variable.
    */
    struct Logger : public HanaConfigurable<Logger, LoggerConfig>, virtual public INamed {

        /** Construct a logger with a default "gropu name" equal to the package name.

            User can later set the "type name".
         */
        Logger(const std::string& group_name="spng");
        /** Construct logger with fully qualified names

        */
        Logger(const std::string& type_name, const std::string& group_name="spng");
        virtual ~Logger() = default;

        // INamed
        /** Return a name for this instance.

            This will return the name set by set_name() or if none set, a
            default is constructed from the configured log and gropu names.
        */
        virtual std::string get_name() const;

        /** Set the instance name.

            If this class, or its subclass, is constructed via the NamedFactory,
            the instance name is set based on user configuration and generally
            should not be reset by the user.
        */
        virtual void set_name(const std::string& name);

        /// Subclass call these to chirp about a tensor set in a standard way,
        /// and reactive to "verbosity" configuration.
        virtual void logit(const std::string& context) const;
        virtual void logit(const ITorchTensorSet::pointer& ts, const std::string& context="") const;
        virtual void logit(const ITorchTensor::pointer& ten, const std::string& context="") const;
        virtual void logit(const ITorchTensor::vector& tens, const std::string& context="") const;
        virtual void logit(const TensorIndex& ti, const std::string& context="") const;

        /** @brief A numerical verbosity level controlling logit()'s

            0 means emit nothing.

            1+ means level of hierarchy iteration.

            Example of the iteration:, logging a tensor set at 1 only prints
            things about the set itself.  Logging at 2 will iterate the set and
            print a summary each tensor.  At 3, details of each tensor will be
            printed.
        */
        void set_verbosity(int level);
        int verbosity() const;

        // Respond to configuration change.
        virtual void configured();

        
        /** @brief The logger includes a "count" in log messages.

            The logger does not automatically advance the count.  The user
            should call next() to do that.

            This is intended to be used to executions of the DFP node that uses
            a logger.
        */
        int set_count(int count) { m_count = count; return m_count; };
        int next_count() { return ++m_count; }
        int get_count() const { return m_count; }

        /** @brief Replicate main log methods.

            This allows almost the same logging API if the Logger is used as a
            base class or as a class member variable eg named "log".  
        */
        template <typename... Args>
        void debug(Args&&... args) {
            log->debug(std::forward<Args>(args)...);
        }
        template <typename... Args>
        void info(Args&&... args) {
            log->info(std::forward<Args>(args)...);
        }
        template <typename... Args>
        void warn(Args&&... args) {
            log->warn(std::forward<Args>(args)...);
        }
        template <typename... Args>
        void critical(Args&&... args) {
            log->critical(std::forward<Args>(args)...);
        }

    protected:

        /// The actual log object can be accessed via subclass or use debug(), etc.
        WireCell::Log::logptr_t log;

        /// User should call next_count() to increment.
        mutable size_t m_count{0};

    private:

        /// See get_name().
        std::string m_interface_name = "";

        /// We want to mutate this to handle recursive calls.
        mutable int m_verbosity = 0;

    };
}

#endif
